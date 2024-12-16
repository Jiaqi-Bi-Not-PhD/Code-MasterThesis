library(ggplot2)
library(MASS)  
library(tmvtnorm)


penf <- function(est, x, age, base.dist="Weibull", frailty.dist=NULL, agemin=18, cuts=NULL){
  nbase <- length(est) - length(x)
  
  base.est <- exp(est[1:nbase])
  if(base.dist == "lognormal") base.est[1] <- est[1]
  
  if(is.null(frailty.dist) || frailty.dist == "none"){
    xbeta <- sum(est[-c(1:nbase)] * x)
    H <- cumhaz(base.dist, age - agemin, base.est, cuts = cuts) * exp(xbeta)
    pen <- 1 - exp(-H)
  } else {
    k <- est[nbase + 1]
    xbeta <- sum(est[-c(1:(nbase + 1))] * x)
    H <- cumhaz(base.dist, age - agemin, base.est, cuts = cuts) * exp(xbeta)
    pen <- 1 - laplace(dist = frailty.dist, g = H, k = k)
  }
  
  return(pen)
}

generate_penetrance_ci <- function(est, vcov_matrix, frailty_dist, x, age_seq, n_simulations = 1000) {
  
  theta_replicates <- mvrnorm(n_simulations, mu = est, Sigma = vcov_matrix)
  
  penetrance_values <- matrix(NA, nrow = n_simulations, ncol = length(age_seq))
  
  for (i in 1:n_simulations) {
    penetrance_values[i, ] <- sapply(age_seq, function(age) penf(theta_replicates[i, ], x, age, frailty.dist = frailty_dist, agemin = 18))
  }
  
  valid_indices <- complete.cases(penetrance_values)
  penetrance_values <- penetrance_values[valid_indices, ]
  
  ci_lower <- apply(penetrance_values, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_upper <- apply(penetrance_values, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  return(list(ci_lower = ci_lower, ci_upper = ci_upper))
}

###################################################################################################################
##################### BRCA 1 Gamma CCA #######################
###################################################################################################################

vcov_matrix <- matrix(
  c(
    1.293566e-02, -4.272685e-05, -4.199538e-02, -1.312291e-03, -1.780518e-04,
    -4.272685e-05,  3.184451e-03,  7.033335e-03,  3.408487e-04, -1.716147e-02,
    -4.199538e-02,  7.033335e-03,  1.70360492e-01, 1.16226710e-02, -5.18932485e-02,
    -1.312291e-03,  3.408487e-04,  1.16226710e-02, 2.87436785e-02,  1.58956010e-03,
    -1.780518e-04, -1.716147e-02, -5.18932485e-02, 1.58956010e-03,  2.927867246e-01
  ),
  nrow = 5, byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa")
  )
)
est <- c(-4.0384, 1.1787, 1.0179, 0.2795, 0.9274)
est[5] <- exp(est[5])
frailty_dist <- "gamma"
#vcov_matrix <- BRCA1_lognormal_K_H0T$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_noKH0T_lognormal[[1]]$Estimates # Change the estimates
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Gamma Frailty Model - CCA)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Gamma CCA.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Log Normal CCA #######################
###################################################################################################################

vcov_matrix <- matrix(
  c(
    0.0137733677, -0.0007847323, -0.0448828760, -0.0011170178,  0.0126216789,
    -0.0007847323,  0.0033322024,  0.0066324015,  0.0004570801, -0.0177779774,
    -0.0448828763,  0.0066324015,  0.1719708320,  0.0115754581, -0.0421486000,
    -0.0011170178,  0.0004570801,  0.0115754580,  0.0284608208,  0.0002683010,
    0.0126216789, -0.0177779774, -0.0421486000,  0.0002683010,  0.2618058170
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa")
  )
)
est <- c(-4.0950737,1.1878431,1.0058789,0.2841896 ,0.8075770)
est[5] <- exp(est[5])
frailty_dist <- "lognormal"
#vcov_matrix <- BRCA1_lognormal_K_H0T$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_noKH0T_lognormal[[1]]$Estimates # Change the estimates
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Log Normal Frailty Model - CCA)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Log Normal CCA.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Gamma MI-LOGT #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA1_noKlogT_gamma$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma[[1]]$Estimates # Change the estimates
vcov_matrix <- matrix(
  c(
    0.0044305287,  0.0006626435, -0.0119907570, -0.0001639542,  0.0036745846,
    0.0006626435,  0.0017720878,  0.0034422570,  0.0004739397, -0.0087251476,
    -0.0119907575,  0.0034422567,  0.0624314070,  0.0089224728, -0.0411338313,
    -0.0001639542,  0.0004739397,  0.0089224730,  0.0280736935, -0.0009240991,
    0.0036745846, -0.0087251476, -0.0411338310, -0.0009240991,  0.2377505054
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
est <- c(-4.52, 1.17, 2.24, 0.288, 1.11)
#rownames(vcov_matrix) <- colnames(vcov_matrix) <- c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
#est <- c(-4.51, 0.931, 2.27, 0.513, 1.11)

est[5] <- exp(est[5])
frailty_dist <- "gamma"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Gamma Frailty Model - MI-LOGT)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Gamma MI LOGT.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Log Normal MI-LOGT #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA1_noKlogT_lognormal$`Pooled Var-Cov` # Change the covariance matrix
vcov_matrix <- matrix(
  c(
    0.0057054232,  0.0003024115, -0.0149381552,  0.0014274240,  0.0140691078,
    0.0003024115,  0.0018047543,  0.0034404395,  0.0004369887, -0.0092760568,
    -0.0149381552,  0.0034404395,  0.0658700800,  0.0068919102, -0.0402174460,
    0.0014274240,  0.0004369887,  0.0068919102,  0.0364333611,  0.0075892554,
    0.0140691078, -0.0092760568, -0.0402174460,  0.0075892554,  0.2370849280
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
est <- c(-4.56,1.18,2.22,0.337,1.18)
#est <- BRCA1_noKlogT_lognormal[[1]]$Estimates # Change the estimates
est[5] <- exp(est[5])
frailty_dist <- "lognormal"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Log Normal Frailty Model - MI-LOGT)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Log Normal MI LOGT.png")


data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Gamma MI-H0T #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA1_noKH0T_gamma$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_noKH0T_gamma[[1]]$Estimates # Change the estimates
est <- c(-4.67, 0.91, 2.22, 0.144, 1.11)
vcov_matrix <- matrix(
  c(
    0.0015119685,  0.0004660541, -0.0023043883, -0.0005262026,  0.0006230351,
    0.0004660541,  0.0008322693,  0.0013422866,  0.0006732604, -0.0017020301,
    -0.0023043883,  0.0013422866,  0.0179534310,  0.0141412393, -0.0043958069,
    -0.0005262026,  0.0006732604,  0.0141412393,  0.0333606367,  0.0012628642,
    0.0006230351, -0.0017020301, -0.0043958069,  0.0012628642,  0.0529569948
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)

est[5] <- exp(est[5])
frailty_dist <- "gamma"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Gamma Frailty Model - MI-H0T)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Gamma MI H0T.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Log Normal MI-H0T #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA1_noKH0T_lognormal$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_noKH0T_lognormal[[1]]$Estimates # Change the estimates
vcov_matrix <- matrix(
  c(
    0.0021302707,  0.0004758841, -0.0023187909,  0.0023178650,  0.0030060558,
    0.0004758841,  0.0008667426,  0.0011702750,  0.0008484437, -0.0018170584,
    -0.0023187909,  0.0011702750,  0.0127700090,  0.0032656665, -0.0045506741,
    0.0023178650,  0.0008484437,  0.0032656665,  0.0274482496, -0.0007117185,
    0.0030060558, -0.0018170584, -0.0045506741, -0.0007117185,  0.0448953596
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
rownames(vcov_matrix) <- colnames(vcov_matrix) <- c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
est <- c(-4.72, 0.92, 2.23, 0.338, 1.06)
est[5] <- exp(est[5])
frailty_dist <- "lognormal"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Log Normal Frailty Model - MI-H0T)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Log Normal MI H0T.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Gamma MI-K-LOGT #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA2_K_LOGT_gamma$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA2_K_LOGT_gamma[[1]]$Estimates # Change the estimates
vcov_matrix <- matrix(
  c(
    1.717539e-03,  8.971586e-05, -3.534895e-03, -1.114951e-02,  1.109067e-03,
    8.971586e-05,  1.227137e-03,  3.429787e-03,  1.220222e-02, -2.310984e-03,
    -3.534895e-03,  3.429787e-03,  0.025627868,   0.07512227,  -0.006606065,
    -1.114951e-02,  1.220222e-02,  0.075122268,   0.37178233,  -0.018080979,
    1.109067e-03, -2.310984e-03, -0.006606065,  -0.01808098,   0.050574562
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
rownames(vcov_matrix) <- colnames(vcov_matrix) <- c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
est <- c(-4.68, 0.916, 2.19, 0.186, 1.08)

est[5] <- exp(est[5])
frailty_dist <- "gamma"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Gamma MI K LOGT.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Log Normal MI-K-LOGT #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA2_K_LOGT_lognormal$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA2_K_LOGT_lognormal[[1]]$Estimates # Change the estimates
vcov_matrix <- matrix(
  c(
    0.0026617019,  0.0007958518, -0.0018831321,  0.0077113961,  0.0043064465,
    0.0007958518,  0.0009781821,  0.0013074096,  0.0039820152, -0.0012526757,
    -0.0018831321,  0.0013074096,  0.0186281350,  0.0305413570, -0.0056467750,
    0.0077113961,  0.0039820152,  0.0305413570,  0.1827159740,  0.0050528350,
    0.0043064465, -0.0012526757, -0.0056467750,  0.0050528350,  0.0462872930
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
rownames(vcov_matrix) <- colnames(vcov_matrix) <- c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
est <- c(-4.73, 0.924, 2.25, 0.394, 1.04)

est[5] <- exp(est[5])
frailty_dist <- "lognormal"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Log Normal Frailty Model - MI-K-LOGT)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Log Normal MI K LOGT.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Gamma MI-K-H0T #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA1_gamma_K_H0T$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_gamma_K_H0T[[1]]$Estimates # Change the estimates
vcov_matrix <- matrix(
  c(
    0.0023082874,  0.0007931672, -0.0018804293,  0.0094845536,  0.0003482212,
    0.0007931672,  0.0009251295,  0.0009800761,  0.0033033603, -0.0016451061,
    -0.0018804293,  0.0009800761,  0.0123647992,  0.0033401483, -0.0033017236,
    0.0094845536,  0.0033033603,  0.0033401483,  0.1105324592,  0.0003615726,
    0.0003482212, -0.0016451061, -0.0033017236,  0.0003615726,  0.0500374001
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
est <- c(-4.66, 0.912, 2.16, 0.116, 1.10)

est[5] <- exp(est[5])
frailty_dist <- "gamma"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Gamma Frailty Model - MI-K-H0T)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Gamma MI K H0T.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]

###################################################################################################################
##################### BRCA 1 Log Normal MI-K-H0T #######################
###################################################################################################################

#vcov_matrix <- CCA_model_lognormal$varcov
#est <- CCA_model_lognormal$estimates
#vcov_matrix <- BRCA1_lognormal_K_H0T$`Pooled Var-Cov` # Change the covariance matrix
#est <- BRCA1_lognormal_K_H0T[[1]]$Estimates # Change the estimates
vcov_matrix <- matrix(
  c(
    0.0020247021,  0.0005627869, -0.0017932137,  0.0067338527,  0.0024358692,
    0.0005627869,  0.0008850863,  0.0008990831,  0.0016800469, -0.0019272160,
    -0.0017932137,  0.0008990831,  0.0156947923,  0.0170759161, -0.0009442950,
    0.0067338527,  0.0016800469,  0.0170759161,  0.1562720550,  0.0061156510,
    0.0024358692, -0.0019272160, -0.0009442950,  0.0061156510,  0.0437705530
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
  )
)
rownames(vcov_matrix) <- colnames(vcov_matrix) <- c("log.lambda", "log.rho", "mgeneI", "PRS_I", "log.kappa")
est <- c(-4.73, 0.92, 2.16, 0.101, 1.01)

est[5] <- exp(est[5])
frailty_dist <- "lognormal"
#est <- c(-4.73, 0.92, 2.22, 0.23, exp(1.05))
J <- diag(5)
J[5,5] <- est[5]
transformed_vcov_matrix <- J %*% vcov_matrix %*% t(J)
transformed_vcov_matrix[5, 5] <- est[5]^2 * vcov_matrix[5, 5]
vcov_matrix <- transformed_vcov_matrix
#vcov_matrix <- BRCA1_noKlogT_gamma$varcov # Change the covariance matrix
#est <- BRCA1_noKlogT_gamma$estimates # Change the estimates
#est <- c(-4.04, 1.18, 1.02, 0.28, 0.92)
age <- seq(18, 70, by = 1)

data <- expand.grid(age = age, mgene = c(0, 1), PRS = c(-0.5, 0, 0.5))

data$ci_lower <- NA
data$ci_upper <- NA

# Define a function to process each row
process_row <- function(i, data, est, frailty_dist, vcov_matrix, age) {
  x <- c(data$mgene[i], data$PRS[i])
  ci_results <- generate_penetrance_ci(est, frailty_dist = frailty_dist, vcov_matrix, x, age, n_simulations = 1000)
  
  # Return the results as a list
  list(
    ci_lower = ci_results$ci_lower[which(age == data$age[i])],
    ci_upper = ci_results$ci_upper[which(age == data$age[i])]
  )
}

# Apply the function in parallel
results <- mclapply(1:nrow(data), function(i) process_row(i, data, est, frailty_dist, vcov_matrix, age), mc.cores = 6)

# Combine the results into the data frame
data$ci_lower <- sapply(results, function(res) res$ci_lower)
data$ci_upper <- sapply(results, function(res) res$ci_upper)


# Apply LOESS smoothing to CI lines
for (m in unique(data$mgene)) {
  for (p in unique(data$PRS)) {
    sub_data <- data[data$mgene == m & data$PRS == p, ]
    
    if(nrow(sub_data) > 1) {  
      loess_lower <- loess(ci_lower ~ age, data = sub_data, span = 0.2)
      loess_upper <- loess(ci_upper ~ age, data = sub_data, span = 0.2)
      
      data$ci_lower[data$mgene == m & data$PRS == p] <- predict(loess_lower)
      data$ci_upper[data$mgene == m & data$PRS == p] <- predict(loess_upper)
    }
  }
}

data$penetrance <- sapply(1:nrow(data), function(i) penf(est, c(data$mgene[i], data$PRS[i]), data$age[i], frailty.dist = frailty_dist, agemin = 18))

data$mgene <- factor(data$mgene, levels = c(0, 1), labels = c("Non Mutation Carrier", "Mutation Carrier"))
data$PRS <- factor(data$PRS, levels = c(-0.5, 0, 0.5), labels = c("PRS = -0.5", "PRS = 0", "PRS = 0.5"))

data <- data |>
  mutate(ci_lower = ifelse(ci_lower < 0, 0, ci_lower), 
         ci_upper = ifelse(ci_upper > 1, 1, ci_upper))

# Plotting
brca1_gammaCCA <- ggplot(data, aes(x = age, y = penetrance, color = PRS, fill = PRS)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, linetype = "88") +
  labs(title = "BRCA1 Penetrance Function Plot (Log Normal Frailty Model - MI-K-H0T)",
       x = "Age",
       y = "Penetrance",
       color = "PRS",
       fill = "PRS") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),    # Larger plot title
    axis.title = element_text(size = 20),                                # Larger axis titles
    axis.text = element_text(size = 20),                                 # Larger axis tick labels
    legend.title = element_text(size = 18),                              # Larger legend title
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  facet_wrap(~ mgene)

#brca1_gammaCCA + labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - MI-K-LOGT)")

print(brca1_gammaCCA)
show_and_save(brca1_gammaCCA, "~/Desktop/MSc Biostatistics Western/Jiaqi-Master-Thesis/Thesis Paper/figures/Penetrance plot for Data Analysis; Ver 3.0/BRCA1/BRCA1 Penetrance Plot; Log Normal MI K H0T.png")

data2 <- data |>
  mutate(ci_lower = round(ci_lower * 100, 2),
         ci_upper = round(ci_upper * 100, 2),
         penetrance = round(penetrance * 100,2) )

data2[data2$age == 40 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Non Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Non Mutation Carrier",]

data2[data2$age == 40 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 50 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 60 & data2$mgene=="Mutation Carrier",]
data2[data2$age == 70 & data2$mgene=="Mutation Carrier",]




