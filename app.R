### Library
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(dplyr)
library(ggplot2)

### Precision-based Sample Size Calculator
clopper.pearson.sample.size <- function(min = 20, max = 1500, stepsize = 5, thres, exp_proportion, max_precision,
                                        alpha = 0.025) {
  p = exp_proportion/100  # Use the actual expected proportion entered by user
  n <- seq(min, max, by = stepsize)
  fin <- data.frame()
  stop <- 0
  res <- NA
  
  for (i in n) {
    lbd_CPCI <- binom.test(x = floor(i * p), n = i, p = p, alternative = 'greater', conf.level = (1-alpha))$conf.int[1]
    act_proportion <- floor(i * p)*100/i  # Observed proportion can be different from the expected proportion due to flooring
    lb <- lbd_CPCI*100
    lb_err <- (act_proportion - max_precision - lb)
    act_precision <- act_proportion - lb  # Actual precision : Observed proportion - Derived lower bound
    
    # Check both conditions: lower bound > threshold AND actual precision <= target precision
    if (lb > thres && act_precision <= max_precision) {
      res <- data.frame(Proportion = act_proportion,
                        Precision = act_precision,
                        Size = i)
      stop = stop + 1
    }
    if (stop >0) {break}
  }
  
  if (sum(!is.na(res))>0) {
    colnames(res) <- c( "Observed Proportion (%)", "Observed Precision (%)", "Sample Size")
    result = res
  } else {
    result = "Out of range"
  }
  
  return(result)
}

### Power-based Sample Size Calculator
powerbased_CP = function(p0, p1, alpha = 0.025, power, stepsize=5, max=1500) {
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
      binom.test(x, n, p=p0, alternative = 'greater', conf.level = (1-alpha))$conf.int)
    reject_rej = apply(cis, 2, function(x) {
      if (p0<x[1]) {return(T)} else {return(F)}
    })
    if (sum(reject_rej) >0 ) {
      # Rejection region 
      reject_x <- x_vals[reject_rej]
      accept_x <- x_vals[!reject_rej]
      
      # Critical values
      lower_crit <- c(lower_crit, NA) #min(accept_x)-1)
      upper_crit <- c(upper_crit, max(accept_x)+1)
      
      # Actual alpha (<= prespecified alpha 0.025)
      actual_alpha <- c(actual_alpha, round(max(sapply(reject_x, function(x)
        binom.test(x, n, p=p0, alternative = 'greater', conf.level = (1-alpha))$p.value)),4))
      
      # Power: probability under p_alt of rejecting H0
      cal_power <- round(sum(dbinom(reject_x, n, p1)),3)
      powers = c(powers, cal_power)
      
      if (cal_power >= power) {stop = stop + 1}
      if (stop >0) {
        ntotal = ntotal[ntotal<=n]
        break
      }  
    } else {
      # Critical values
      lower_crit <- c(lower_crit, NA)
      upper_crit <- c(upper_crit, NA)
      actual_alpha <- c(actual_alpha, 0)
      powers = c(powers, 0)
      
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
                         "Alternative Hypothesis Proportion (H1)", 
                         "Required sample size")
  } else {
    result = "Out of range"
  }
  
  
  return(result)
}

### Plot generator under Monte Carlo (MC) simulation
### Note: Power can be analyticially computed due to the nature of exact method;
###       however, we just provide MC simulation result for confirmation/visualization
Visualization_power = function(alpha, p0, p1, n, iter = 10000) {
  
  ## Exact Calculation
  x_vals <- 0:n
  pvals <- sapply(x_vals, function(x)
    binom.test(x, n, p=p0, alternative = 'greater', conf.level = (1-alpha))$p.value)
  
  # Rejection region (p < alpha)
  reject_x <- x_vals[pvals < alpha]
  
  if (length(reject_x) == 0) {
    # Handle case where no rejection region exists
    return(ggplot() + 
             geom_text(aes(x = 0.5, y = 0.5, label = "No rejection region found"), 
                       size = 6) +
             theme_void())
  } else {
    accept_x <- x_vals[pvals >= alpha]
    lower_crit <- min(accept_x)-1
    upper_crit <- max(accept_x)+1
    actual_alpha <- round(max(pvals[pvals<alpha]),4)
    
    # Power: probability under p_alt of rejecting H0
    cal_power <- round(sum(dbinom(reject_x, n, p1)),3)
    
    ## Simulation
    set.seed(10903)
    
    # Simulation (Null hypothesis)
    x0 = rbinom(iter, n, p0)  # Simulated data 
    p0_lb = suppressWarnings(sapply(x0, function(x) binom.test(x, n, p=p0, alternative = 'greater', conf.level = (1-alpha))$conf.int[1])) # Clopper-Pearson LB
    
    # Simulation (Alternative hypothesis)
    x1 = rbinom(iter, n, p1)  # Simulated data 
    p1_lb = suppressWarnings(sapply(x1, function(x) binom.test(x, n, p=p0, alternative = 'greater', conf.level = (1-alpha))$conf.int[1])) # Clopper-Pearson LB
    
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
      labs(x = "Lower Bound of 97.5% Clopper-Pearson One-sided Confidence Interval", y = "Count", 
           title = "Statistical Power Analysis: Monte Carlo Simulation with 10,000 Iterations",
           subtitle = paste0("Analytical Power = ", round(cal_power, 3)*100,"%, ", 
                             "Simulated Power = ", round(power_observed, 3)*100,"%")) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set2")
  }
}

# UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Sample Size Calculator: Clopper-Pearson Method"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h3("Overview", style = "margin-bottom: 10px;"),
        p("This application estimates the required sample size based on the one-sided Clopper-Pearson exact confidence interval (alpha = 0.025)."),
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
                           numericInput("exp_proportion", "Expected Proportion (%)", value = NA, min = 0.1, max = 99.9, step = 0.1),
                           div(style = "margin: 10px 0; padding: 8px; background-color: #f8f9fa; border-left: 4px solid #007bff;",
                               textOutput("threshold_display")),
                           numericInput("max_precision", "Target Precision (%)", value = NA, min = 0.1, max = 99.9, step = 0.1),
                           div(style = "margin: 5px 0; padding: 5px; background-color: #fff3cd; border-left: 3px solid #ffc107; font-size: 0.9em;",
                               "Note: Target Precision will auto-update when Expected Proportion or Threshold changes, but you can override it manually."),
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

server <- function(input, output, session) {
  
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
  
  # Target precision (modified to handle manual input)
  max_precision <- reactive({
    if (!is.na(input$threshold) && !is.na(input$exp_proportion)) {
      calculated_precision <- input$exp_proportion - input$threshold
      
      # If user hasn't manually set max_precision or if it's NA, use calculated value
      if (is.na(input$max_precision)) {
        return(calculated_precision)
      } else {
        return(input$max_precision)
      }
    } else {
      return(input$max_precision)
    }
  })
  
  # Observer to auto-update max_precision input when threshold or exp_proportion changes
  observe({
    if (!is.na(input$threshold) && !is.na(input$exp_proportion)) {
      calculated_precision <- input$exp_proportion - input$threshold
      
      # Only update if the current max_precision is NA or if it matches the previous calculation
      # This prevents overwriting user's manual input
      if (is.na(input$max_precision)) {
        updateNumericInput(session, "max_precision", value = calculated_precision)
      }
    }
  })
  
  # To display the expected proportion (updated)
  output$threshold_display <- renderText({
    if (!is.na(input$threshold) && !is.na(input$exp_proportion)) {
      calculated_precision <- input$exp_proportion - input$threshold
      current_precision <- max_precision()
      
      if (!is.null(current_precision)) {
        if (abs(current_precision - calculated_precision) < 0.001) {
          paste0("Target Precision: ", round(current_precision, 1), "% (Auto-calculated)")
        } else {
          paste0("Target Precision: ", round(current_precision, 1), "% (Manually set | Auto-calculated would be: ", round(calculated_precision, 1), "%)")
        }
      } else {
        "Target Precision: Please enter valid threshold and expected proportion values"
      }
    } else {
      "Target Precision: Please enter valid threshold and expected proportion values"
    }
  })
  
  # Observe Power-Based action button
  observeEvent(input$run_power, {
    params$method <- "Power-Based"
    params$power <- input$power
    params$p0 <- input$p0
    params$p1 <- input$p1
    
    req(params$power, params$p0, params$p1)
    
    # Error message case 1
    if (params$p1 <= params$p0) {
      error_state(TRUE)
      error_message("Alternative hypothesis proportion must be greater than null hypothesis proportion.")
      result(NULL)
      return()
    }
    
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
  
  # Observe Precision-Based action button (updated)
  observeEvent(input$run_precision, {
    params$method <- "Precision-Based"
    params$threshold <- input$threshold
    params$p1 <- input$exp_proportion
    params$max_precision <- max_precision()  # Use the reactive value
    
    req(params$threshold, params$p1, params$max_precision)
    
    # Error message case 1
    if (params$max_precision <= 0) {
      error_state(TRUE)
      error_message("Target precision must be greater than 0.")
      result(NULL)
      return()
    }
    
    # Additional error check: if manually set precision is larger than expected proportion
    if (params$max_precision >= params$p1) {
      error_state(TRUE)
      error_message("Target precision must be less than the expected proportion.")
      result(NULL)
      return()
    }
    
    withProgress(message = "Calculating required sample size...", value = 0.5, {
      calc_result <- clopper.pearson.sample.size(
        thres = params$threshold,
        exp_proportion = params$p1,  # Pass expected proportion separately
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
  
  # Reactive expression to generate the plot
  plot_data <- reactive({
    req(result())
    
    if (params$method == "Power-Based" && !is.null(result())) {
      Visualization_power(
        alpha = 0.025,
        p0 = params$p0 / 100,
        p1 = params$p1 / 100,
        n = result()$`Required sample size`
      )
    } else if (params$method == "Precision-Based" && !is.null(result()) ) {
      Visualization_power(
        alpha = 0.025,
        p0 = params$threshold / 100,
        p1 = params$p1 / 100,
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



