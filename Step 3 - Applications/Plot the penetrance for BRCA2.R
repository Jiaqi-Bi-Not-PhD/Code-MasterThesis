Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

library(ggplot2)
library(MASS)
library(tmvtnorm)
library(parallel)
library(dplyr)


source("FamEvent/R/cumhaz.R")
source("FamEvent/R/laplace.R")
source("FamEvent/R/dlaplace.R")
source("FamEvent/R/gh.R")

show_and_save <- function(plot, filename) {
  print(plot)
  ggsave(filename, plot = plot, width = 11, height = 7.36, units = "in")
}

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
    0.0143879461,  0.0005711908, -0.0571507215,  0.0010446755, -0.0026940757,
    0.0005711908,  0.0047664607,  0.0090410320,  0.0037287773, -0.0113490466,
    -0.0571507215,  0.0090410320,  0.2958790770,  0.0256793100, -0.0177389720,
    0.0010446755,  0.0037287773,  0.0256793100,  0.0879675550, -0.0138973060,
    -0.0026940757, -0.0113490466, -0.0177389720, -0.0138973060,  0.1197691830
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa")
  )
)
est <- c(-4.00226304 , 1.44845709 , 1.13622066 , 0.26012284, -0.01915905 )
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
  labs(title = "BRCA2 Penetrance Function Plot (Gamma Frailty Model - CCA)",
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

show_and_save(brca1_gammaCCA, "BRCA2 Penetrance Plot; Gamma CCA.png")



###################################################################################################################
##################### BRCA 1 Log Normal CCA #######################
###################################################################################################################

vcov_matrix <- matrix(
  c(
    0.0145881719, -0.0003384238, -0.0574628925, -0.0001497695,  0.0072946483,
    -0.0003384238,  0.0051530205,  0.0090221934,  0.0035293537, -0.0143796688,
    -0.0574628925,  0.0090221934,  0.2891708660,  0.0234289277, -0.0235353560,
    -0.0001497695,  0.0035293537,  0.0234289277,  0.0888205933, -0.0149836565,
    0.0072946483, -0.0143796688, -0.0235353560, -0.0149836565,  0.1202486570
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa"),
    c("log.lambda", "log.rho", "mgeneI", "PRS", "log.kappa")
  )
)
est <- c(-4.11120140 , 1.44777619 , 1.05652312 , 0.25058446, -0.03469157)
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
  labs(title = "BRCA2 Penetrance Function Plot (Log Normal Frailty Model - CCA)",
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
show_and_save(brca1_gammaCCA, "BRCA2 Penetrance Plot; Log Normal CCA.png")






