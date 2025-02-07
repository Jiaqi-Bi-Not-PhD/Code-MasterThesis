####################################################################
####################### Rubin's Rule ###############################
############## Pooling Results & Variance Estimation ###############
####################################################################
Pooling <- function(model, imputed_data, robust) {
  theta <- lapply(model, function(x) x$estimates) 
  if(robust == FALSE) SE <- lapply(model, function(x) x$se)
  else SE <- lapply(model, function(x) x$se.robust)
  
  ## Var-Cov matrices
  if (robust == TRUE) varcov_matrices <- lapply(model, function(x) x$varcov.robust)
  else varcov_matrices <- lapply(model, function(x) x$varcov)
  
  theta_matrix <- do.call(rbind, theta)
  SE_matrix <- do.call(rbind, SE)
  
  k <- ncol(theta_matrix) # number of parameters
  m <- nrow(theta_matrix) # number of imputations
  n <- nrow(imputed_data[[1]])
  
  Var_W <- colMeans(SE_matrix^2, na.rm = TRUE)
  theta_means <- colMeans(theta_matrix, na.rm = TRUE) 
  names_theta <- names(theta_means)
  Var_B <- colSums((theta_matrix - matrix(theta_means, nrow = m, ncol = k, byrow = TRUE))^2, na.rm = TRUE) / (m - 1)
  Var_T <- Var_W + Var_B + Var_B/m
  
  Wald_pooled <- ( theta_means^2 )/Var_T
  r <- (Var_B + Var_B/m)/(Var_W) # r
  rho <- ( Var_B + Var_B/m )/(Var_T) # rho
  df_old <- (m-1)*( (1+1/r)^2 ) # Old df
  df_obs <- ( ((n-k)+1)/((n-k)+3) ) * (n-k) * (1-rho) # df obs
  df_adj <- (df_old * df_obs)/(df_old + df_obs) # adjusted df
  
  p_values <- 2 * pt(-abs(sqrt(Wald_pooled)), df_adj) # P Value
  
  SE_Pooled <- sqrt(Var_T)
  Lower_CI <- theta_means - qt((1-0.05/2), df_adj) * SE_Pooled
  Upper_CI <- theta_means + qt((1-0.05/2), df_adj) * SE_Pooled
  
  ## Pooling the varcov matrix
  pooled_within_varcov <- Reduce("+", varcov_matrices) / m
  
  ## Between-imputation variance for the variance-covariance matrix
  mean_theta <- colMeans(theta_matrix)
  between_varcov <- Reduce("+", lapply(1:m, function(i) {
    (theta_matrix[i, ] - mean_theta) %*% t(theta_matrix[i, ] - mean_theta)
  })) / (m - 1)
  
  ## Total pooled variance-covariance matrix
  pooled_varcov <- pooled_within_varcov + (1 + 1/m) * between_varcov
  
  return(list(tibble(Parameters = names_theta,
              Estimates = theta_means, 
              `Within Imp. Var.` = Var_W, 
              `Between Imp. Var.` = Var_B,
              `Total Var.` = Var_T, 
              `p values` = p_values,
              Lower.CI = Lower_CI, 
              Upper.CI = Upper_CI),
              `Pooled Var-Cov` = pooled_varcov))
}

