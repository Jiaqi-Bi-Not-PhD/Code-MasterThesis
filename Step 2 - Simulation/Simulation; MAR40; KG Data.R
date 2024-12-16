Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

#############################################################################################################
############################################## Packages and Sources #########################################
#############################################################################################################

library(parallelly)
library(kinship2)
library(survival)
library(Matrix)
library(lme4)
library(parallel)
library(mice)
library(dplyr)
library(MASS)
library(lme4qtl)
#MI_FamEvent_K_H0T
#library(FamEvent)
library(truncnorm)
library(numDeriv)
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

source("find_closest.R")
source("Bias, MSE, Coverage Functions.R")
source("MI function; K log(H0(T)).R") # MI_FamEvent_K_H0T
source("MI function; No K log(H0(T)).R") # MI_FamEvent_noK_H0T
source("MI function; Random Intercept logT.R") # MI_FamEvent_noK_logT
source("MI function; with K log(T); df using Satterthwaite Approximation.R") # MI_FamEvent_K_logT
#source("MI function; with K log(T); variances adjusted.R") # MI_FamEvent_K_logT using lmer variances
source("Rubins Rule.R")
source("Suppress Cat.R")

#############################################################################################################
############################################ Simulated Data #################################################
#############################################################################################################

simulated_dataset_500_MAR_Data_Gamma20 <- readRDS("simulated_dataset_500_MAR40_Data_Kinship.RData")

Parameters <- c("log.lambda", "log.rho", "mgene", "newx_I", "log.kappa")
true_value <- c(-4.71, 0.804, 2.13, 1, log(1)) # log(1) - HighCorr; 1.15 - LowCorr
true_value2 <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(log(1))) # exp(log(1)) - HighCorr; exp(1.15) - LowCorr
True_df <- tibble::tibble(Parameters, true_value)
n_cores <- parallelly::availableCores()

frailty_dist <- "gamma"

#############################################################################################################
############################################ CCA ############################################################
#############################################################################################################
Parameters <- c("log.lambda", "log.rho", "mgene", "newx", "log.kappa")
True_df <- tibble::tibble(Parameters, true_value)
methods = "CCA"
simulation_penmodel <- function(i, data_list, frailty.dist, True_dfarg = True_df) {
  data_sim <- data_list[[i]]
  tryCatch({
    Complete_model <- quiet(penmodel(survival::Surv(time, status) ~ mgene + newx, cluster = "famID", 
                                     gvar = "mgene", 
                                     design = "pop+", base.dist = "Weibull", frailty.dist = frailty.dist, 
                                     agemin = 18, 
                                     data = data_sim,
                                     parms = true_value2))
    
    df <- tibble::tibble(
      Parameters = names(Complete_model$estimates),
      Estimates = Complete_model$estimates,
      SE = Complete_model$se,
      SimulationID = i,
      Method = "CCA"
    )
    
    merged_df <- merge(df, True_dfarg, by = "Parameters")
    return(merged_df)
  }, error = function(e) {
    message("Error occurred: ", e$message)
    return(NULL)
  })
}

results_Gamma <- parallel::mclapply(1:500, simulation_penmodel,
                                    data_list = simulated_dataset_500_MAR_Data_Gamma20, 
                                    frailty.dist = frailty_dist, mc.cores = n_cores)
results_Gamma <- do.call(rbind, results_Gamma)

results_CCA <- results_Gamma |>
  mutate(Bias = Estimates - true_value)
#saveRDS(results_Gamma, "Simulation Raw Results; CCA; Gamma using Log Normal 20; High Corr.RData")

n_sim <- length(unique(results_CCA$SimulationID))

print(methods)
results_CCA |>
  group_by(Parameters) |>
  summarise(MCSE_Bias = round(sd(Bias, na.rm = TRUE)/ sqrt(n_sim), 3),
            Bias = round(mean(Bias, na.rm = TRUE), 3),
            Ave.Est = round( mean(Estimates, na.rm = TRUE), 3),
            True_Est = mean(true_value),
            Model_SE = round( mean(SE, na.rm = TRUE), 3), 
            EmpSE = round( sd(Estimates, na.rm = TRUE), 3), 
            RMSE = round( sqrt(mean(Estimates - true_value)^2), 3), 
            Coverage = mean( (Estimates - 1.96 * SE <= true_value) & (true_value <= Estimates + 1.96 * SE )) )

#############################################################################################################
############################################ noK logT #######################################################
#############################################################################################################
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")
Parameters <- c("log.lambda", "log.rho", "mgene", "newx_I", "log.kappa")
True_df <- tibble::tibble(Parameters, true_value)

methods <- "MI-logT"
## Automated simulation function - lapply using parallel to the list of 1000 datasets
process_datasets_K_H0T <- function(datasets, true_value) {
  num_simulations <- length(datasets)
  results <- parallel::mclapply(seq_len(num_simulations), function(i) {
    data <- datasets[[i]]
    
    imputed_results <- MI_FamEvent_noK_logT(
      data, 
      option = "General", 
      M = 5, 
      frailty.dist = frailty_dist, 
      true_value = true_value2, 
      robust = FALSE
    )
    imputed_results <- imputed_results[[1]]
    
    if (any(is.na(imputed_results))) return(NULL)
    
    metrics <- calculate_metrics(imputed_results, true_value)
    metrics$SimulationID <- i
    return(metrics)
  }, mc.cores = n_cores, mc.silent = TRUE)
  results <- purrr::compact(results)
  results <- dplyr::bind_rows(results)
  return(results)
}

## Start running the simulation
sim_results <- process_datasets_K_H0T(
  datasets = simulated_dataset_500_MAR_Data_Gamma20,  ## <- Change name accordingly here! ##
  true_value = true_value
)

## On some cases, dataset is abandoned due to model divergence
#nsim_left <- nrow(sim_results) / length(true_value)
sim_results <- sim_results |>
  dplyr::select(SimulationID, everything())
sim_results_noKlogT <- sim_results |>
  dplyr::mutate(Method = methods)

n_sim <- length(unique(sim_results$SimulationID))

sim_results_toprint <- sim_results |> 
  dplyr::group_by(Parameters) |> 
  dplyr::summarise(Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 )
  )
sim_results_toprint <- merge(sim_results_toprint, True_df, by = "Parameters")
sim_results_toprint <- sim_results_toprint |>
  dplyr::select(c(Parameters, true_value, Ave.Est, Bias, MCSE_Bias, RMSE, MSE, EmpSE, Ave.ModSE, Coverage))

scenario <- "simulated_dataset_500_MAR_Data_Gamma20" ## Change name accordingly ##
final_results <- list(n_sim, scenario, sim_results_toprint)
print(methods)
print(final_results)

#############################################################################################################
############################################# noK H0T #######################################################
#############################################################################################################
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

methods <- "MI-H0T"
## Automated simulation function - lapply using parallel to the list of 1000 datasets
process_datasets_K_H0T <- function(datasets, true_value) {
  num_simulations <- length(datasets)
  results <- parallel::mclapply(seq_len(num_simulations), function(i) {
    data <- datasets[[i]]
    
    imputed_results <- MI_FamEvent_noK_H0T(
      data, 
      option = "General", 
      M = 5, 
      frailty.dist = frailty_dist, 
      true_value = true_value2, 
      robust = FALSE
    )
    imputed_results <- imputed_results[[1]]
    
    if (any(is.na(imputed_results))) return(NULL)
    
    metrics <- calculate_metrics(imputed_results, true_value)
    metrics$SimulationID <- i
    return(metrics)
  }, mc.cores = n_cores, mc.silent = TRUE)
  results <- purrr::compact(results)
  results <- dplyr::bind_rows(results)
  return(results)
}

## Start running the simulation
sim_results <- process_datasets_K_H0T(
  datasets = simulated_dataset_500_MAR_Data_Gamma20,  ## <- Change name accordingly here! ##
  true_value = true_value
)

## On some cases, dataset is abandoned due to model divergence
#nsim_left <- nrow(sim_results) / length(true_value)
sim_results <- sim_results |>
  dplyr::select(SimulationID, everything())
sim_results_noKH0T <- sim_results |>
  dplyr::mutate(Method = methods) ####### Change name here #######

n_sim <- length(unique(sim_results_noKH0T$SimulationID))

sim_results_toprint <- sim_results_noKH0T |> 
  dplyr::group_by(Parameters) |> 
  dplyr::summarise(Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 )
  )
sim_results_toprint <- merge(sim_results_toprint, True_df, by = "Parameters")
sim_results_toprint <- sim_results_toprint |>
  dplyr::select(c(Parameters, true_value, Ave.Est, Bias, MCSE_Bias, RMSE, MSE, EmpSE, Ave.ModSE, Coverage))

scenario <- "simulated_dataset_500_MAR_Data_Gamma20" ## Change name accordingly ##
final_results <- list(n_sim, scenario, sim_results_toprint)
print(methods)
print(final_results)

#############################################################################################################
############################################## K logT #######################################################
#############################################################################################################
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

methods <- "MI-K-logT"
## Automated simulation function - lapply using parallel to the list of 1000 datasets
process_datasets_K_H0T <- function(datasets, true_value) {
  num_simulations <- length(datasets)
  results <- parallel::mclapply(seq_len(num_simulations), function(i) {
    data <- datasets[[i]]
    
    imputed_results <- MI_FamEvent_K_logT(
      data, 
      option = "General", 
      M = 5, 
      frailty.dist = frailty_dist, 
      true_value = true_value2, 
      robust = FALSE
    )
    imputed_results <- imputed_results[[1]]
    
    if (any(is.na(imputed_results))) return(NULL)
    
    metrics <- calculate_metrics(imputed_results, true_value)
    metrics$SimulationID <- i
    return(metrics)
  }, mc.cores = n_cores, mc.silent = TRUE)
  results <- purrr::compact(results)
  results <- dplyr::bind_rows(results)
  return(results)
}

## Start running the simulation
sim_results <- process_datasets_K_H0T(
  datasets = simulated_dataset_500_MAR_Data_Gamma20,  ## <- Change name accordingly here! ##
  true_value = true_value
)

## On some cases, dataset is abandoned due to model divergence
#nsim_left <- nrow(sim_results) / length(true_value)
sim_results <- sim_results |>
  dplyr::select(SimulationID, everything())
sim_results_KlogT <- sim_results |>
  dplyr::mutate(Method = methods) ####### Change name here #######

n_sim <- length(unique(sim_results$SimulationID))

sim_results_toprint <- sim_results |> 
  dplyr::group_by(Parameters) |> 
  dplyr::summarise(Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 )
  )
sim_results_toprint <- merge(sim_results_toprint, True_df, by = "Parameters")
sim_results_toprint <- sim_results_toprint |>
  dplyr::select(c(Parameters, true_value, Ave.Est, Bias, MCSE_Bias, RMSE, MSE, EmpSE, Ave.ModSE, Coverage))

scenario <- "simulated_dataset_500_MAR_Data_Gamma20" ## Change name accordingly ##
final_results <- list(n_sim, scenario, sim_results_toprint)
print(methods)
print(final_results)

#############################################################################################################
############################################### K H0T #######################################################
#############################################################################################################
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1",
           LAPACK_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

methods <- "MI-K-H0T"
## Automated simulation function - lapply using parallel to the list of 1000 datasets
process_datasets_K_H0T <- function(datasets, true_value) {
  num_simulations <- length(datasets)
  results <- parallel::mclapply(seq_len(num_simulations), function(i) {
    data <- datasets[[i]]
    
    imputed_results <- MI_FamEvent_K_H0T(
      data, 
      option = "General", 
      M = 5, 
      frailty.dist = frailty_dist, 
      true_value = true_value2, 
      robust = FALSE
    )
    imputed_results <- imputed_results[[1]]
    
    if (any(is.na(imputed_results))) return(NULL)
    
    metrics <- calculate_metrics(imputed_results, true_value)
    metrics$SimulationID <- i
    return(metrics)
  }, mc.cores = n_cores, mc.silent = TRUE)
  results <- purrr::compact(results)
  results <- dplyr::bind_rows(results)
  return(results)
}

## Start running the simulation
sim_results <- process_datasets_K_H0T(
  datasets = simulated_dataset_500_MAR_Data_Gamma20,  ## <- Change name accordingly here! ##
  true_value = true_value
)

## On some cases, dataset is abandoned due to model divergence
#nsim_left <- nrow(sim_results) / length(true_value)
sim_results <- sim_results |>
  dplyr::select(SimulationID, everything())
sim_results_KH0T <- sim_results |>
  dplyr::mutate(Method = methods) ####### Change name here #######

n_sim <- length(unique(sim_results$SimulationID))

sim_results_toprint <- sim_results |> 
  dplyr::group_by(Parameters) |> 
  dplyr::summarise(Ave.Est = round( mean(Estimates, na.rm = TRUE), 3 ),
                   MCSE_Bias = round(sd(Bias, na.rm = TRUE) / sqrt(n_sim), 3),
                   Bias = round( mean(Bias, na.rm = TRUE), 3 ),
                   RMSE = round( sqrt(mean(MSE, na.rm = TRUE)), 3 ),
                   MSE = round( mean(MSE, na.rm = TRUE), 3 ),
                   Coverage = mean(Coverage, na.rm = TRUE),
                   EmpSE = round( sd(Estimates, na.rm = TRUE), 3 ) ,
                   Ave.ModSE = round( sqrt(mean(Total_Variance, na.rm = TRUE)), 3 )
  )
sim_results_toprint <- merge(sim_results_toprint, True_df, by = "Parameters")
sim_results_toprint <- sim_results_toprint |>
  dplyr::select(c(Parameters, true_value, Ave.Est, Bias, MCSE_Bias, RMSE, MSE, EmpSE, Ave.ModSE, Coverage))

scenario <- "simulated_dataset_500_MAR_Data_Gamma20" ## Change name accordingly ##
final_results <- list(n_sim, scenario, sim_results_toprint)
print(methods)
print(final_results)

#############################################################################################################
############################################### Combine results #############################################
#############################################################################################################

AllSimResults <- bind_rows(results_CCA, sim_results_noKlogT, 
                           sim_results_noKH0T, sim_results_KlogT,
                           sim_results_KH0T)
saveRDS(AllSimResults, "MAR40 KG Misspecified; Simulation Results.RData")