# Sample Size Calculator for OTC Application: Clopper-Pearson Exact Method


## Overview
An interactive R Shiny web application for calculating sample sizes using the two-sided Clopper-Pearson exact confidence interval method. This tool supports both power-based and precision-based sample size calculations with real-time visualization based on simulation.


## Features
- **Power-Based Calculation**: Determine required sample size based on desired statistical power, null hypothesis proportion, and alternative hypothesis proportion
- **Precision-Based Calculation**: Calculate sample size based on target threshold and maximum precision requirements
- **Interactive Visualization**: Monte Carlo simulation plots showing statistical power analysis with 10,000 iterations


## Methodology
The application implements the Clopper-Pearson exact method for constructing confidence intervals, which provides conservative coverage probabilities. Sample size calculations are performed within a range of 20 to 1,500 subjects with a default increment of 5 per iteration.


## Application Structure
### Core Functions
- `powerbased_CP()`: Power-based sample size calculation
- `clopper.pearson.sample.size()`: Precision-based sample size calculation
- `Visualization_power()`: Monte Carlo simulation for statistical power visualization

### Input Parameters
**Power-Based Method:**
- Power (%) - Desired statistical power (0.1-99.9%)
- Null Hypothesis Proportion (%) - H₀ proportion (0.1-99.9%)
- Alternative Hypothesis Proportion (%) - H₁ proportion (0.1-99.9%)

**Precision-Based Method:**
- Target Threshold (%) - Prespecified performance threshold (0.1-99.9%)
- Maximum Precision (%) - Maximum allowable precision (0.1-10%)

### Output
- **Results Table**: Sample size requirements with relevant statistical parameters
- **Visualization**: Statistical power visualization under the proposed sample size and the parameters


## Example Usage (https://taekwonhong.shinyapps.io/CPCI_App/)
1. Select either "Power-Based" or "Precision-Based" tab
2. Enter required parameters in the input fields
3. Click the corresponding "Run" button
4. Review results in the data table and visualization plot


## Contact
Taekwon Hong
Division of Biometrics VII, OB, OTS, CDER
U.S. Food and Drug Administration
Department of Health and Human Services
taekwon.hong@fda.hhs.gov
