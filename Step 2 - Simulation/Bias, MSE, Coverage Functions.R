###########################################################
################### Simulation Functions ##################
###########################################################
#true_value <- c(-3.3524072, 0.8329091, 1.0000000, 3.0000000, 3.0000000, 0.6931472)

calculate_metrics <- function(imp_results, true_value) {
  Parameters <- imp_results$Parameters
  parameter_estimates <- imp_results$Estimates
  total_variance <- imp_results$`Total Var.`
  lower_ci <- imp_results$`Lower.CI`
  upper_ci <- imp_results$`Upper.CI`
  imp_results$true_value <- true_value
  True_Value <- imp_results$true_value
  Estimates <- imp_results$Estimates
  
  bias_ind <- parameter_estimates - True_Value
  mse_ind <- (parameter_estimates - True_Value)^2
  coverage_check <- (lower_ci <= True_Value & upper_ci >= True_Value)
  
  return(tibble(
    Parameters = Parameters,
    Estimates = Estimates,
    Total_Variance = total_variance,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    Bias = bias_ind,
    MSE = mse_ind,
    Coverage = coverage_check
  ))
}
#calculate_metrics(test_results, true_value)
