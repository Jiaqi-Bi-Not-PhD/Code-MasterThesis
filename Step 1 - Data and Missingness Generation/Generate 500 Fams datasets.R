library(parallel)
library(dplyr)
library(survival)
library(parallelly)
library(mice)
#library(FamEvent)
library(truncnorm)
library(MASS)
library(kinship2)

#source("Functions; Before Generating Data.R")
source("find_closest.R")
#source("gh_REVISED; May 27.R")
#source("hazards.R")
#source("hermite.R")
#source("laplace.R")
#source("loglik_frailty_REVISED; May 27.R")
#source("penmodel_REVISED; May 27.R")
source("Suppress Cat.R")
source("Delete Males.R")
source("Suppress Cat.R")
source("Delete Males.R")
source("FamEvent/R/cumhaz.R")
source("FamEvent/R/hazards.R")
source("FamEvent/R/gh.R")
source("FamEvent/R/penmodel.R")
source("FamEvent/R/loglik_frailty.R")
source("FamEvent/R/dlaplace.R")
source("FamEvent/R/laplace.R")

source("FamEvent/R/familyDesign.R")
source("FamEvent/R/fgeneZX.R")
source("FamEvent/R/Pgene.R")
source("FamEvent/R/surv.dist.R")
source("FamEvent/R/survp.dist.R")
source("FamEvent/R/inv.surv.R")
source("FamEvent/R/inv2.surv.R")
source("FamEvent/R/parents.g.R")
source("FamEvent/R/kids.g.R")

source("FamEvent/R/simfam.R")

source("familyStructure_REVISED; July 3 2024.R")
source("loglik_frailty; REVISED; July 3 2024.R")

#######################################################################################
## True thetas - penetrance function, compare with real data
true_value2 <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))

n_simulations <- 500
n_cores <- parallelly::availableCores()
#cl <- parallel::makeCluster(n_cores)

#######################################################################################
################################ Generate 1000 datasets ###############################
#######################################################################################
## 500, 200, 50 families => Gamma or LogNormal
Nfam_values <- 500
frailty_dist_values <- "gamma"

scenarios <- expand.grid(Nfam = Nfam_values, frailty_dist = frailty_dist_values)
scenarios_list <- split(scenarios, seq(nrow(scenarios)))

run_simulation <- function(params, error_index) {
  #library(survival)
  #library(survPen)
  
  Nfam <- params$Nfam
  frailty_dist <- params$frailty_dist
  
  repeat{
  results <- tryCatch({
  ## Data Generation
  famx <- simfam(N.fam = Nfam, design = "pop+", variation = "frailty", 
                             base.dist = "Weibull", frailty.dist = frailty_dist, interaction = FALSE,
                             add.x = TRUE, x.dist = "mvnormal", x.parms = c(0, 1), depend = exp(1.15), 
                             base.parms = c(exp(-4.71), exp(0.804)), vbeta = c(0, 2.13, 1.0), agemin = 18,
                 allelefreq = 0.0021) 
  #famx <- famx |> dplyr::filter(gender == 0)
    
  return(famx) 
  }, error = function(e) NULL )
  if (!is.null(results)) break
  }
}
#ch <- penmodel(survival::Surv(time, status) ~ mgene + newx, cluster = "famID", gvar = "mgene", 
#         design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
#         agemin = 18, data = famx,
#         parms = c(exp(-4.71), exp(0.804), 2.13, 1.0, exp(log(0.5))))
#summary(ch)

#kinship_mat <- with(famx, kinship2::kinship(id = indID, dadid = fatherID, momid = motherID,
#                                            sex = gender))
#kinship_mat_sparse <- Matrix::Matrix(kinship_mat, sparse = TRUE)
#kinship_mat_sparse <- 2 * kinship_mat_sparse
#kinship_mat <- as.matrix(kinship_mat_sparse)

#Iden_mat <- diag(nrow = nrow(famx)) # Identity matrix n*n
#Iden_mat_sparse <- Matrix::Matrix(Iden_mat, sparse = TRUE)

#model_test <- lme4qtl::relmatLmer(newx ~ log(time) + status + mgene + (1|indID), data = famx, relmat = list(indID = #kinship_mat))
#summary(model_test)
## Generate 500 complete datasets for 2*3 = 6 types of data
#n_cores <- parallelly::availableCores()   # 7 cores on my own computer

# Create a cluster

#clusterEvalQ(cl, {
#  library(survival)
#  library(survPen)
#})

#parallel::clusterExport(cl, c("penmodel", "run_simulation", "gh",
#                    "loglik_frailty", "cumhaz", "hazards", "Surv", "dlaplace",
#                    "laplace", "scenarios_list", "n_simulations"))  

## Apply the function in parallel for all scenarios
#results_parallel <- parallel::parLapply(cl, seq_along(scenarios_list), function(scenario_idx) {
#  scenario_params <- scenarios_list[[scenario_idx]]
#  lapply(1:n_simulations, function(sim_index) {
#    run_simulation(scenario_params, sim_index)
#  })
#})

#total_simulations <- length(scenarios_list) * n_simulations
#results_parallel <- parallel::parLapply(cl, 1:total_simulations, function(i) {
#  scenario_idx <- ((i - 1) %/% n_simulations) + 1
#  sim_index <- ((i - 1) %% n_simulations) + 1
#  scenario_params <- scenarios_list[[scenario_idx]]
#  run_simulation(scenario_params, sim_index)
#})

results_parallel <- vector("list", length(scenarios_list))

for (scenario_idx in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[scenario_idx]]
  
  results_parallel[[scenario_idx]] <- parallel::mclapply(1:n_simulations, function(sim_index) {
    run_simulation(scenario_params, sim_index)
  }, mc.cores = n_cores)
}



#stopCluster(cl)

for (i in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[i]]
  scenario_name <- paste("simulated_dataset_500_Complete", scenario_params$frailty_dist, scenario_params$Nfam, "fam", sep = "_")
  
  assign(scenario_name, results_parallel[[i]])
}

saveRDS(simulated_dataset_500_Complete_gamma_500_fam, file = "simulated_dataset_500_Complete_gamma_500_fam.RData")

#######################################################################################
################################ Generate 1000 datasets ends ##########################
#######################################################################################
Complete_Data_Gamma <- simulated_dataset_500_Complete_gamma_500_fam
Complete_Data_Gamma <- lapply(Complete_Data_Gamma,
                              delete_males)


true_value <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
true_value2 <- c(-4.71, 0.804, 2.13, 1, 1.15)
Parameters <- c("log.lambda", "log.rho", "mgene", "newx", "log.kappa")
True_df <- tibble::tibble(Parameters, true_value2)

simulation_penmodel <- function(i, data_list, frailty.dist) {
  data_sim <- data_list[[i]]
  tryCatch({
    Complete_model <- penmodel(survival::Surv(time, status) ~ mgene + newx, cluster = "famID", 
                               gvar = "mgene", 
                               design = "pop+", base.dist = "Weibull", frailty.dist = "lognormal", 
                               agemin = 18, 
                               data = data_sim,
                               parms = true_value)
    return(list(Estimates = Complete_model$estimates,
                SE = Complete_model$se))
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)
  })
}


results_Gamma <- parallel::mclapply(1:500, simulation_penmodel,
                                    data_list = Complete_Data_Gamma, 
                                    frailty.dist = "gamma", mc.cores = n_cores)

Gamma_est <- data.frame( do.call(rbind, lapply(results_Gamma, function(x) x$Estimates)) )
Gamma_est
Gamma_SE <- data.frame( do.call(rbind, lapply(results_Gamma, function(x) x$SE)) )
Gamma_SE

colnames(Gamma_SE) <- paste0(colnames(Gamma_SE), "_SE")
combined_data <- cbind(Gamma_est, Gamma_SE)
clean_data <- na.omit(combined_data)
n_cols <- ncol(Gamma_est)
Gamma_est <- clean_data[,1:n_cols]
Gamma_SE <- clean_data[,(n_cols + 1) : (2*n_cols)]
colnames(Gamma_SE) <- colnames(Gamma_est)

nsim_left_gamma <- nrow(Gamma_est)
nsim_left_gamma


print("Gamma 500 Simulations - Complete Data - Males Deleted - Gamma fit with Log Normal")

True_df |>
  rowwise() |>
  mutate(Parameters = Parameters, 
         Ave.Est = round( mean(Gamma_est[[Parameters]], na.rm = TRUE), 3 ),
         #CI_Length = mean(Gamma_est[[Parameters]] - 1.96 * Gamma_est[[Parameters]], na.rm = TRUE) - mean(Gamma_est[[Parameters]] + 1.96 * Gamma_est[[Parameters]], na.rm = TRUE),
         Bias = round( mean(Gamma_est[[Parameters]], na.rm = TRUE) - true_value2, 3 ),
         Model_SE = round( mean(Gamma_SE[[Parameters]], na.rm = TRUE), 3 ), 
         EmpSE = round( sd(Gamma_est[[Parameters]], na.rm = TRUE), 3 ), 
         RMSE = round( sqrt(mean((Gamma_est[[Parameters]] - true_value2)^2)), 3 ),
         Coverage = mean((Gamma_est[[Parameters]] - 1.96 * Gamma_SE[[Parameters]] <= true_value2) & (true_value2 <= Gamma_est[[Parameters]] + 1.96 * Gamma_SE[[Parameters]]), na.rm = TRUE) ) |>
  print(width = Inf)

#######################################################################################
################################ Generate 1000 datasets ###############################
#######################################################################################
## 500, 200, 50 families => Gamma or LogNormal
Nfam_values <- 500
frailty_dist_values <- "lognormal"

scenarios <- expand.grid(Nfam = Nfam_values, frailty_dist = frailty_dist_values)
scenarios_list <- split(scenarios, seq(nrow(scenarios)))

run_simulation <- function(params, error_index) {
  #library(survival)
  #library(survPen)
  
  Nfam <- params$Nfam
  frailty_dist <- params$frailty_dist
  
  repeat{
  results <- tryCatch({
  ## Data Generation
  famx <- simfam(N.fam = Nfam, design = "pop+", variation = "frailty", 
                 base.dist = "Weibull", frailty.dist = frailty_dist, interaction = FALSE,
                 add.x = TRUE, x.dist = "mvnormal", x.parms = c(0, 1), depend = exp(1.15), 
                 base.parms = c(exp(-4.71), exp(0.804)), vbeta = c(0, 2.13, 1), agemin = 18,
                 allelefreq = 0.0021) 
  #famx <- famx |> dplyr::filter(gender == 0)
  
  return(famx)
  }, error = function(e) NULL 
  )
  if (!is.null(results)) break
  }
}


## Generate 500 complete datasets for 2*3 = 6 types of data
#n_cores <- parallelly::availableCores()   # 7 cores on my own computer

# Create a cluster

#clusterEvalQ(cl, {
#  library(survival)
#  library(survPen)
#})

#parallel::clusterExport(cl, c("penmodel", "run_simulation", "gh",
#                    "loglik_frailty", "cumhaz", "hazards", "Surv", "dlaplace",
#                    "laplace", "scenarios_list", "n_simulations"))  

## Apply the function in parallel for all scenarios
#results_parallel <- parallel::parLapply(cl, seq_along(scenarios_list), function(scenario_idx) {
#  scenario_params <- scenarios_list[[scenario_idx]]
#  lapply(1:n_simulations, function(sim_index) {
#    run_simulation(scenario_params, sim_index)
#  })
#})

#total_simulations <- length(scenarios_list) * n_simulations
#results_parallel <- parallel::parLapply(cl, 1:total_simulations, function(i) {
#  scenario_idx <- ((i - 1) %/% n_simulations) + 1
#  sim_index <- ((i - 1) %% n_simulations) + 1
#  scenario_params <- scenarios_list[[scenario_idx]]
#  run_simulation(scenario_params, sim_index)
#})

results_parallel <- vector("list", length(scenarios_list))

for (scenario_idx in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[scenario_idx]]
  
  results_parallel[[scenario_idx]] <- parallel::mclapply(1:n_simulations, function(sim_index) {
    run_simulation(scenario_params, sim_index)
  }, mc.cores = n_cores)
}



#stopCluster(cl)

for (i in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[i]]
  scenario_name <- paste("simulated_dataset_500_Complete", scenario_params$frailty_dist, scenario_params$Nfam, "fam", sep = "_")
  
  assign(scenario_name, results_parallel[[i]])
}

saveRDS(simulated_dataset_500_Complete_lognormal_500_fam, file = "simulated_dataset_500_Complete_lognormal_500_fam.RData")

#######################################################################################
################################ Generate 1000 datasets ends ##########################
#######################################################################################
Complete_Data_LogNormal <- simulated_dataset_500_Complete_lognormal_500_fam
Complete_Data_LogNormal <- lapply(Complete_Data_LogNormal,
                                  delete_males)


true_value <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
true_value2 <- c(-4.71, 0.804, 2.13, 1, 1.15)
Parameters <- c("log.lambda", "log.rho", "mgene", "newx", "log.kappa")
True_df <- tibble::tibble(Parameters, true_value2)

simulation_penmodel <- function(i, data_list, frailty.dist) {
  data_sim <- data_list[[i]]
  tryCatch({
    Complete_model <- quiet(penmodel(survival::Surv(time, status) ~ mgene + newx, cluster = "famID", 
                                     gvar = "mgene", 
                                     design = "pop+", base.dist = "Weibull", frailty.dist = "gamma", 
                                     agemin = 18, 
                                     data = data_sim,
                                     parms = true_value))
    return(list(Estimates = Complete_model$estimates,
                SE = Complete_model$se))
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)
  })
}


results_LogNormal <- parallel::mclapply(1:500, simulation_penmodel,
                                        data_list = Complete_Data_LogNormal, 
                                        frailty.dist = "lognormal", mc.cores = n_cores)

LogNormal_est <- data.frame( do.call(rbind, lapply(results_LogNormal, function(x) x$Estimates)) )
LogNormal_est
LogNormal_SE <- data.frame( do.call(rbind, lapply(results_LogNormal, function(x) x$SE)) )
LogNormal_SE

colnames(LogNormal_SE) <- paste0(colnames(LogNormal_SE), "_SE")
combined_data <- cbind(LogNormal_est, LogNormal_SE)
clean_data <- na.omit(combined_data)
n_cols <- ncol(LogNormal_est)
LogNormal_est <- clean_data[,1:n_cols]
LogNormal_SE <- clean_data[,(n_cols + 1) : (2*n_cols)]
colnames(LogNormal_SE) <- colnames(LogNormal_est)

nsim_left_lognormal <- nrow(LogNormal_est)
nsim_left_lognormal


print("Log Normal 500 Simulations - Complete Data - Males deleted - Log Normal fit with Gamma")

True_df |>
  rowwise() |>
  mutate(Parameters = Parameters, 
         Ave.Est = round( mean(LogNormal_est[[Parameters]], na.rm = TRUE), 3 ),
         #CI_Length = mean(LogNormal_est[[Parameters]] + 1.96 * LogNormal_est[[Parameters]], na.rm = TRUE) - mean(LogNormal_est[[Parameters]] - 1.96 * LogNormal_est[[Parameters]], na.rm = TRUE),
         Bias = round( mean(LogNormal_est[[Parameters]], na.rm = TRUE) - true_value2, 3 ),
         Model_SE = round( mean(LogNormal_SE[[Parameters]], na.rm = TRUE), 3 ), 
         EmpSE = round( sd(LogNormal_est[[Parameters]], na.rm = TRUE), 3 ),
         RMSE = round( sqrt(mean((LogNormal_est[[Parameters]] - true_value2)^2, na.rm = TRUE)), 3 ),
         Coverage = mean((LogNormal_est[[Parameters]] - 1.96 * LogNormal_SE[[Parameters]] <= true_value2) & (true_value2 <= LogNormal_est[[Parameters]] + 1.96 * LogNormal_SE[[Parameters]]), na.rm = TRUE),
  ) |>
  print(width = Inf)

#######################################################################################
################################ Generate 5000 datasets ###############################
#######################################################################################
## 500, 200, 50 families => Gamma or LogNormal
Nfam_values <- 500
frailty_dist_values <- "lognormal"

scenarios <- expand.grid(Nfam = Nfam_values, frailty_dist = frailty_dist_values)
scenarios_list <- split(scenarios, seq(nrow(scenarios)))

run_simulation <- function(params, error_index) {
  #library(survival)
  #library(survPen)
  
  Nfam <- params$Nfam
  frailty_dist <- params$frailty_dist
  
  repeat{
  results <- tryCatch({
  ## Data Generation
  famx <- simfam(N.fam = Nfam, design = "pop+", variation = "kinship", 
                 base.dist = "Weibull", frailty.dist = frailty_dist, interaction = FALSE,
                 add.x = TRUE, x.dist = "mvnormal", x.parms = c(0, 1), depend = exp(1.15), 
                 base.parms = c(exp(-4.71), exp(0.804)), vbeta = c(0, 2.13, 1), agemin = 18,
                 allelefreq = 0.0021) 
  #famx <- famx |> dplyr::filter(gender == 0)
  
  return(famx)
  }, error = function(e) NULL 
  )
  if (!is.null(results)) break
  }
}


## Generate 500 complete datasets for 2*3 = 6 types of data
#n_cores <- parallelly::availableCores()   # 7 cores on my own computer

# Create a cluster

#clusterEvalQ(cl, {
#  library(survival)
#  library(survPen)
#})

#parallel::clusterExport(cl, c("penmodel", "run_simulation", "gh",
#                    "loglik_frailty", "cumhaz", "hazards", "Surv", "dlaplace",
#                    "laplace", "scenarios_list", "n_simulations"))  

## Apply the function in parallel for all scenarios
#results_parallel <- parallel::parLapply(cl, seq_along(scenarios_list), function(scenario_idx) {
#  scenario_params <- scenarios_list[[scenario_idx]]
#  lapply(1:n_simulations, function(sim_index) {
#    run_simulation(scenario_params, sim_index)
#  })
#})

#total_simulations <- length(scenarios_list) * n_simulations
#results_parallel <- parallel::parLapply(cl, 1:total_simulations, function(i) {
#  scenario_idx <- ((i - 1) %/% n_simulations) + 1
#  sim_index <- ((i - 1) %% n_simulations) + 1
#  scenario_params <- scenarios_list[[scenario_idx]]
#  run_simulation(scenario_params, sim_index)
#})

results_parallel <- vector("list", length(scenarios_list))

for (scenario_idx in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[scenario_idx]]
  
  results_parallel[[scenario_idx]] <- parallel::mclapply(1:n_simulations, function(sim_index) {
    run_simulation(scenario_params, sim_index)
  }, mc.cores = n_cores)
}



#stopCluster(cl)

for (i in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[i]]
  scenario_name <- paste("simulated_dataset_500_Complete", scenario_params$frailty_dist, scenario_params$Nfam, "fam", sep = "_")
  
  assign(scenario_name, results_parallel[[i]])
}

saveRDS(simulated_dataset_500_Complete_lognormal_500_fam, file = "simulated_dataset_500_Complete_kinship_500_fam.RData")

#######################################################################################
################################ Generate 5000 datasets ends ##########################
#######################################################################################
delete_males <- function(data) {
  data <- data |> dplyr::filter(gender == 0)
}

Complete_Data_LogNormal_kinship <- simulated_dataset_500_Complete_lognormal_500_fam
Complete_Data_LogNormal_kinship <- lapply(Complete_Data_LogNormal_kinship,
                                          delete_males)


Parameters <- c("log.lambda", "log.rho", "mgene", "newx", "log.kappa")
true_value2 <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
true_value <- c(-4.71, 0.804, 2.13, 1, 1.15)
True_df <- tibble::tibble(Parameters, true_value)

simulation_penmodel <- function(i, data_list, frailty.dist) {
  data_sim <- data_list[[i]]
  tryCatch({
    Complete_model <- quiet(penmodel(survival::Surv(time, status) ~ mgene + newx, cluster = "famID", 
                                     gvar = "mgene", 
                                     design = "pop+", base.dist = "Weibull", frailty.dist = frailty.dist, 
                                     agemin = 18, 
                                     data = data_sim,
                                     parms = true_value2))
    return(list(Estimates = Complete_model$estimates,
                SE = Complete_model$se))
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)
  })
}


results_LogNormal <- parallel::mclapply(1:500, simulation_penmodel,
                                        data_list = Complete_Data_LogNormal_kinship, 
                                        frailty.dist = "lognormal", mc.cores = n_cores)

LogNormal_est <- data.frame( do.call(rbind, lapply(results_LogNormal, function(x) x$Estimates)) )
LogNormal_est
LogNormal_SE <- data.frame( do.call(rbind, lapply(results_LogNormal, function(x) x$SE)) )
LogNormal_SE


colnames(LogNormal_SE) <- paste0(colnames(LogNormal_SE), "_SE")
combined_data <- cbind(LogNormal_est, LogNormal_SE)
clean_data <- na.omit(combined_data)
n_cols <- ncol(LogNormal_est)
LogNormal_est <- clean_data[,1:n_cols]
LogNormal_SE <- clean_data[,(n_cols + 1) : (2*n_cols)]
colnames(LogNormal_SE) <- colnames(LogNormal_est)

nsim_left_lognormal <- nrow(LogNormal_est)


print("Kinship using LogNormal - Males deleted")
nsim_left_lognormal

True_df |>
  rowwise() |>
  mutate(Parameters = Parameters, 
         Ave.Est = round( mean(LogNormal_est[[Parameters]], na.rm = TRUE), 3 ),
         #CI_Length = mean(LogNormal_est[[Parameters]] + 1.96 * LogNormal_est[[Parameters]], na.rm = TRUE) - mean(LogNormal_est[[Parameters]] - 1.96 * LogNormal_est[[Parameters]], na.rm = TRUE),
         Bias = round( mean(LogNormal_est[[Parameters]], na.rm = TRUE) - true_value, 3 ),
         Model_SE = round( mean(LogNormal_SE[[Parameters]], na.rm = TRUE), 3 ), 
         EmpSE = round( sd(LogNormal_est[[Parameters]], na.rm = TRUE), 3 ),
         RMSE = round( sqrt(mean((LogNormal_est[[Parameters]] - true_value)^2, na.rm = TRUE)), 3 ),
         Coverage = mean((LogNormal_est[[Parameters]] - 1.96 * LogNormal_SE[[Parameters]] <= true_value) & (true_value <= LogNormal_est[[Parameters]] + 1.96 * LogNormal_SE[[Parameters]]), na.rm = TRUE),
  ) |>
  print(width = Inf)
