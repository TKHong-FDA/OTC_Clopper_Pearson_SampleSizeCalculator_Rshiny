### Library
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(samplingbook)
library(DT)
library(dplyr)
library(ggplot2)


### Precision-based Sample Size Calculator
clopper.pearson.sample.size <- function(min = 20, max = 1500, stepsize = 5, thres, max_precision,
                                        alpha = 0.05) {
  p = (thres+max_precision)/100   # Using equation: d (precision) <= d* (maximum precision) = p (expected proportion) - p0 (threshold)
  n <- seq(min, max, by = stepsize)
  fin <- data.frame()
  stop <- 0
  
  for (i in n) {
    CPCI <- Sprop(m = round(i * p), n = i, level = (1 - alpha))$ci$cp
    act_proportion <- round(i * p)*100/i  # Observed proportion can be different from the expected proportion due to rounding
    lb <- CPCI[1]*100
    ub <- CPCI[2]*100
    lb_err <- (act_proportion - max_precision - lb)
    act_precision <- act_proportion - lb   # Actual precision (<= d*) : Observed proportion - Derived lower bound
    
    if ((lb_err < 0)&(lb > thres)) {
      df <- data.frame(Proportion = act_proportion,
                       Precision = act_precision,
                       Size = i)
      stop = stop + 1
    }
    if (stop >0) {break}
  }
  
  if (nrow(df)!=0) {
    colnames(df) <- c( "Observed Proportion (%)", "Precision (%)", "Sample Size")
    result = df
  } else {
    result = "Out of range"
  }
  
  return(result)
}

### Power-based Sample Size Calculator
powerbased_CP = function(p0, p1, alpha = 0.05, power, stepsize=5, max=1500) {
  ntotal <- seq(from=20, to=max, by=stepsize)  
  
  powers = c()
  lower_crit = c()
  upper_crit = c()
  actual_alpha = c()
  stop = 0
  for (n in ntotal) {
    # Determine rejection region using binomial test under H0
    x_vals <- 0:n
    cis = sapply(x_vals, function(x)
      binom.test(x, n, p=p0, alternative = 'two.sided')$conf.int)
    reject_rej = apply(cis, 2, function(x) {
      if (p0<x[1]|p0>x[2]) {return(T)} else {return(F)}
    })
    
    # Rejection region 
    reject_x <- x_vals[reject_rej]
    accept_x <- x_vals[!reject_rej]
    
    # Critical values
    lower_crit <- c(lower_crit, min(accept_x)-1)
    upper_crit <- c(upper_crit, max(accept_x)+1)
    
    # Actual alpha (<= prespecified alpha 0.05)
    actual_alpha <- c(actual_alpha, round(max(sapply(reject_x, function(x)
      binom.test(x, n, p=p0, alternative = 'two.sided')$p.value)),4))
    
    # Power: probability under p_alt of rejecting H0
    cal_power <- round(sum(dbinom(reject_x, n, p1)),3)
    powers = c(powers, cal_power)
    
    if (cal_power >= power) {stop = stop + 1}
    if (stop >0) {
      ntotal = ntotal[ntotal<=n]
      break
    }
  }
  ## This will not be shown in the result; however, they are mimicking PROC POWER onesamplefreq test=exact result ##
  minitab = data.frame(Ntotal = ntotal, lower = lower_crit, upper = upper_crit, alpha = actual_alpha, power = powers)
  colnames(minitab) = c("N Total", "Lower Crit Val", "Upper Crit Val", "Actual Alpha", "Power")
  ##################################################################################################################
  
  # Result
  result = data.frame(Power = powers*100, p0 = p0*100, p1 = p1*100, size = ntotal )
  result = result %>% filter((Power>=power*100)) %>% mutate(diff = abs(Power-power*100))
  if (nrow(result)!=0) {
    result = result %>% filter(size == min(size)) %>% dplyr::select(-diff)
    colnames(result) = c("Power (%)", "Null Hypothesis Proportion (H0)", 
                         "Alternative Hypothesis Proportion (H1)", #"Mean 95% Clopper-Pearson LB (%)", 
                         "Required sample size")
  } else {
    result = "Out of range"
  }
  
  
  return(result)
}

### Plot generator under Monte Carlo simulation
Visualization_power = function(alpha, p0, p1, n, iter = 10000) {
  
  ## Exact Calculation
  x_vals <- 0:n
  pvals <- sapply(x_vals, function(x)
    binom.test(x, n, p = p0, alternative = "two.sided")$p.value)
  
  # Rejection region (p < alpha)
  reject_x <- x_vals[pvals < 0.05]
  accept_x <- x_vals[pvals >= 0.05]
  lower_crit <- min(accept_x)-1
  upper_crit <- max(accept_x)+1
  actual_alpha <- round(max(pvals[pvals<0.05]),4)
  
  # Power: probability under p_alt of rejecting H0
  cal_power <- round(sum(dbinom(reject_x, n, p1)),3)
  
  ## Simulation
  set.seed(10903)
  
  # Simulation (Null hypothesis)
  x0 = rbinom(iter, n, p0)  # Simulated data 
  p0_lb = suppressWarnings(sapply(x0, function(x) Sprop(m = x, n = n, level = (1 - alpha))$ci$cp)[1,]) # Clopper-Pearson 95% LB
  
  # Simulation (Alternative hypothesis)
  x1 = rbinom(iter, n, p1)  # Simulated data 
  p1_lb = suppressWarnings(sapply(x1, function(x) Sprop(m = x, n = n, level = (1 - alpha))$ci$cp)[1,]) # Clopper-Pearson 95% LB
  
  # Calculate standard error for the alternative hypothesis
  lbs = c(p0_lb, p1_lb)
  
  # Create a data frame for plotting
  df <- data.frame(Hypothesis = c(rep("Null Hypothesis", length(p0_lb)),
                                  rep("Alternative Hypothesis", length(p1_lb))),
                   LB = lbs)
  
  # Compute alpha and power
  alpha_observed <- mean(p0_lb > p0)
  power_observed <- mean(p1_lb > p0)
  
  # Create the plot
  ggplot(df, aes(x = LB, fill = Hypothesis)) +
    geom_histogram(position = "identity", alpha = 0.7, bins = 50) +
    geom_vline(xintercept = p0, linetype = "dashed", color = "red") +
    labs(x = "Lower Bound of 95% Clopper-Pearson Confidence Interval", y = "Count", 
         title = "Statistical Power Analysis: Monte Carlo Simulation with 10,000 Iterations",
         subtitle = paste0("Analytical Power =", round(cal_power, 3)*100,"%, ", 
                           "Simulated Power = ", round(power_observed, 3)*100,"%")) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
}

# UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Sample Size Calculator: Clopper-Pearson Method"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h3("Overview", style = "margin-bottom: 10px;"),
        p("This application estimates the required sample size based on the two-sided Clopper-Pearson exact confidence interval (95% CI)."),
        p("The sample size search range is fixed between 20 and 1500, with a default increment of 5 per iteration."),
        p("Please provide the necessary inputs below to begin the calculation.")
      ),
      tabsetPanel(id = "method_type",
                  tabPanel("Power-Based",
                           numericInput("power", "Power (%)", value = NA, min = 0.1, max = 99.9, step = 0.1),
                           numericInput("p0", "Null Hypothesis Proportion (%)", value = NA, min = 0.1, max = 99.9, step = 0.1),
                           numericInput("p1", "Alternative Hypothesis Proportion (%)", value = NA, min = 0.1, max = 99.9, step = 0.1),
                           actionButton("run_power", "Run Power-Based Calculation")
                  ),
                  tabPanel("Precision-Based",
                           numericInput("threshold", "Target Threshold (%)", value = NA, min = 0.1, max = 99.9, step = 0.1),
                           numericInput("max_precision", "Maximum Precision (%)", value = NA, min = 0.1, max = 10, step = 0.1),
                           actionButton("run_precision", "Run Precision-Based Calculation")
                  )
      )
    ),
    mainPanel(
      uiOutput("error_message"),
      div(id = "loading_message", "Processing... Please wait.", 
          style = "display: none; font-weight: bold; color: red;"),
      
      withSpinner(dataTableOutput("table"), type = 6),
      
      # Add plot output
      plotOutput("power_plot")
    )
  )
)

server <- function(input, output) {
  
  # Create reactive values to store the error state and message
  error_state <- reactiveVal(FALSE)
  error_message <- reactiveVal("")
  
  # Create reactive values to store input parameters
  params <- reactiveValues(
    power = NULL,
    p0 = NULL,
    p1 = NULL,
    threshold = NULL,
    max_precision = NULL,
    method = NULL
  )
  
  # Create a reactive value to store the calculation result
  result <- reactiveVal(NULL)
  
  # Observe Power-Based action button
  observeEvent(input$run_power, {
    params$method <- "Power-Based"
    params$power <- input$power
    params$p0 <- input$p0
    params$p1 <- input$p1
    
    req(params$power, params$p0, params$p1)
    
    withProgress(message = "Calculating required sample size...", value = 0.5, {
      calc_result <- powerbased_CP(
        p0 = params$p0 / 100,
        p1 = params$p1 / 100,
        power = params$power / 100,
        stepsize = 5,
        max = 1500
      )
    })
    
    if (identical(calc_result, "Out of range")) {
      error_state(TRUE)
      error_message("The calculation resulted in a sample size outside the acceptable range
                    of 20 to 1,500 subjects. Please adjust your input parameters, particularly 
                    ensuring that the alternative hypothesis proportion is sufficiently larger
                    than the null hypothesis proportion.")
      result(NULL)
    } else {
      error_state(FALSE)
      result(calc_result)
    }
    
  })
  
  # Observe Precision-Based action button
  observeEvent(input$run_precision, {
    params$method <- "Precision-Based"
    params$threshold <- input$threshold
    params$max_precision <- input$max_precision
    
    req(params$threshold, params$max_precision)
    
    withProgress(message = "Calculating required sample size...", value = 0.5, {
      calc_result <- clopper.pearson.sample.size(
        thres = params$threshold,
        max_precision = params$max_precision,
        stepsize = 5
      )
    })
    
    if (identical(calc_result, "Out of range")) {
      error_state(TRUE)
      error_message("The calculation resulted in a sample size outside the acceptable range
                    of 20 to 1,500 subjects. Please adjust your input parameters, particularly 
                    ensuring that the precision is sufficiently large to yield a sample size within
                    the specified range.")
      result(NULL)
    } else {
      error_state(FALSE)
      result(calc_result)
    }
  })
  
  # Create an output for the error message
  output$error_message <- renderUI({
    if(error_state()) {
      div(style = "color: red; font-weight: bold; margin-bottom: 15px;", error_message())
    }
  })
  
  # Render the data table
  output$table <- renderDataTable({
    req(result())
    if (params$method == "Power-Based") {
      datatable(result(), options = list(pageLength = 10)) %>% 
        formatRound(c(1,2,3), digits = 1)
    } else {
      datatable(result(), options = list(pageLength = 10)) %>% 
        formatRound(c(1,2), digits = 1)
    }
  })
  
  # Add this reactive expression to generate the plot
  plot_data <- reactive({
    req(result())
    
    if (params$method == "Power-Based" && !is.null(result())) {
      Visualization_power(
        alpha = 0.05,
        p0 = params$p0 / 100,
        p1 = params$p1 / 100,
        n = result()$`Required sample size`
      )
    } else if (params$method == "Precision-Based" && !is.null(result()) ) {
      Visualization_power(
        alpha = 0.05,
        p0 = params$threshold / 100,
        p1 = result()$`Observed Proportion (%)` / 100,
        n = result()$`Sample Size`
      )
    }
  })
  
  # Render the plot
  output$power_plot <- renderPlot({
    req(plot_data())
    plot_data()
  })
  
}


# Launch App
shinyApp(ui = ui, server = server)

