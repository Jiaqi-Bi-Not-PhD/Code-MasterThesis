library(parallel)
library(dplyr)
library(survival)
library(parallelly)
library(mice)
#library(FamEvent)
library(truncnorm)
library(MASS)
library(kinship2)
#source("cumhaz.R")
#source("dlaplace.R")
source("find_closest.R")
#source("gh_REVISED; May 27.R")
#source("hazards.R")
#source("hermite.R")
#source("laplace.R")
#source("loglik_frailty_REVISED; May 27.R")
#source("penmodel_REVISED; May 27.R")
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


## True thetas

true_value <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
true_value2 <- c(-4.71, 0.804, 2.13, 1, 1.15)

n_simulations <- 500
n_cores <- parallelly::availableCores()

#######################################################################################
################################ Generate MAR #########################################
#######################################################################################
Complete_Data_Gamma <- readRDS("SimData; 500Fams 500 Sims Gamma; HighCorr.RData")
Complete_Data_Gamma <- lapply(Complete_Data_Gamma,
                                  delete_males)
Complete_Data_LogNormal <- readRDS("SimData; 500Fams 500 Sims LogNormal; HighCorr.RData")
Complete_Data_LogNormal <- lapply(Complete_Data_LogNormal,
                              delete_males)
Complete_Data_Kinship <- readRDS("SimData; 500Fams 500 Sims Kinship; HighCorr.RData")
Complete_Data_Kinship <- lapply(Complete_Data_Kinship,
                                  delete_males)

proportions <- c(0.20, 0.40, 0.60, 0.80)

miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)

#parallel::clusterEvalQ(cl, library(mice))
## Missingness & proportion function
introduce_missingness <- function(index, dataset_list, prop, miss_pattern) {
  famx <- dataset_list[[index]]
  
  ampute_test <- mice::ampute(data = famx, prop = prop, patterns = miss_pattern, mech = "MAR")
  reasonable_weights <- ampute_test$weights
  reasonable_weights[1,] <- c(0, 0, 0, 0, 0, -7, 0, 0, 0, 0, 0, 0.05, -0.56, -3.27, 0, 0, 0, 0, 0)
  ampute_test <- mice::ampute(data = famx, prop = prop, patterns = miss_pattern, 
                              mech = "MAR", weights = reasonable_weights)
  miss_famx <- ampute_test$amp
  
  return(miss_famx)
}

## Function to run the missingness introduction for all proportions
run_for_all_proportions <- function(i, dataset_list, proportions, miss_pattern) {
  results <- lapply(proportions, function(prop) {
    introduce_missingness(i, dataset_list, prop, miss_pattern)
  })
  names(results) <- as.character(proportions)
  return(results)
}


#n_cores <- parallelly::availableCores()
#cl <- parallel::makeCluster(n_cores)


## List of complete datasets
complete_datasets <- list(
  Data_Gamma = Complete_Data_Gamma,
  Data_LogNormal = Complete_Data_LogNormal,
  Data_Kinship = Complete_Data_Kinship
)

#parallel::clusterExport(cl, c("introduce_missingness", "run_for_all_proportions", "miss_pattern", "proportions",
#                    "complete_datasets"))

mar_datasets <- list()

## Iterate over each complete dataset and apply the missing data mechanism
for (dataset_name in names(complete_datasets)) {
  dataset <- complete_datasets[[dataset_name]]
  # parallel::clusterExport(cl, "dataset")
  ## Apply the function in parallel
  results_parallel <- parallel::mclapply(1:n_simulations, function(i) {
    run_for_all_proportions(i, dataset, proportions, miss_pattern)
  }, mc.cores=n_cores)
  
  ## Store the results in the mar_datasets list
  for (prop in proportions) {
    mar_datasets[[paste0("simulated_dataset_500_MAR", prop * 100, "_", dataset_name)]] <- parallel::mclapply(1:n_simulations, function(i) results_parallel[[i]][[as.character(prop)]], mc.cores = n_cores)
  }
}

#stopCluster(cl)

### mar_datasets contains 24 lists of 1000 datasets -> 2400 datasets
### nomiss_datasets contains 6 lists of 1000 datasets -> 600 datasets
### mar_datasets
#names(mar_datasets)
#saveRDS(complete_datasets, file = "500 complete datasets lists.RData")
saveRDS(mar_datasets, file = "500_mar_lists.RData")

###### Split data ######
mar_list <- mar_datasets
list_names <- names(mar_list)
for(name in list_names) assign(name, mar_list[[name]])
for(name in list_names) saveRDS(mar_list[[name]], file = paste0(name, ".RData"))
