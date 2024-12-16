#true_value2 <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
source("Suppress Cat.R")
####################################################################################
####################################################################################
####################################################################################
############## Multiple Imputations WITHOUT kinship using log(H_0(T)) ##############
####################################################################################
####################################################################################
####################################################################################
MI_FamEvent_noK_H0T <- function(data, M = 5, option = "General", frailty.dist = "gamma", true_value = true_value2,
                                robust) {
  attempts <- 0
  max_attempts <- 1
  success <- FALSE
  while(attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
      ## Step X - Initial baseline cumulative hazard
      miss50_gamma_cca <- FamEvent::penmodel(survival::Surv(time, status) ~ mgene + newx, 
                                             cluster = "famID", gvar = "mgene",
                                             design = "pop+", base.dist = "Weibull", frailty.dist = "none", 
                                             agemin = 18, data = data,
                                             parms = c(exp(-4.71), exp(0.804), 2.13, 1))
      baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
      logalpha <- baseline_gammafr[1]
      loglambda <- baseline_gammafr[2]
      data <- data |>
        mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
  
  Y_obs <- data$newx
  data_imp <- list()
  for (i in 1:(M+2)) {
    ## Step 1 - Empirical estimates
   # imp_model <- lmer(newx ~ log(log(ageonset)) + status * H0 + mgene + (1|famID), data = data) 
    #X <- model.matrix( ~ log(log(ageonset)) + status * H0 + mgene, data = data)
    #imp_model <- lmer(newx ~ log(log(ageonset)) + status * H0 + mgene + (1|famID), data = data)
    #X <- model.matrix(~ log(log(ageonset)) + status * H0 + mgene, data = data)
    
    imp_model <- lme4::lmer(newx ~ log(ageonset) + status * H0 + mgene + proband + (1|famID) , data = data)
    # Changes made here
    X <- model.matrix(~ log(ageonset) + status * H0 + mgene + proband, data = data)
    
    betas <- lme4::fixef(imp_model)
    p_pred <- ncol(X)-1
    varcorr_list <- lme4::VarCorr(imp_model)
    sigma_e_hat <- attr(varcorr_list, "sc")
    sigma_b_hat <- attr(varcorr_list$famID, "stddev")
    V <- vcov(imp_model)
    data_obs <- data |> dplyr::filter(!is.na(newx))
    n_obs <- nrow(data_obs)
    n <- nrow(data)
    n_fam_obs <- length(unique(data_obs$famID))
    n_fam <- length(unique(data$famID))
    ## Step 2 - g ~ chi^2 nobs-p
    g <- rchisq(n = 1, df = n_obs - p_pred - n_fam_obs + 1)
    g_b <- rchisq(n = 1, df = n_fam_obs - 1)
    
    ## Step 3 - sigma star
    sigma_e_star <- sigma_e_hat * sqrt( (n_obs - p_pred - n_fam_obs + 1)/g )
    sigma_b_star <- sigma_b_hat * sqrt( (n_fam_obs - 1)/g_b )
    sigmastar <- sqrt( sigma_e_star^2 + sigma_b_star^2 )
    sigmahat <- sqrt( sigma_e_hat^2 + sigma_b_hat^2 )
    
    ## Step 4 - u1
    u1 <- rnorm(n = p_pred+1, mean = 0, sd = 1)
    
    ## Step 5
    betastar <- betas + (sigmastar/sigmahat) * u1 %*% chol(V)
    betastar <- as.vector(betastar)
    
    ## Step 6
    u_j <- rnorm(n = n_fam, mean = 0, sd = 1)
    u2i <- rnorm(n = nrow(data), mean = 0, sd = 1)
    newx_Imp <- X %*% betastar + u2i*sigma_e_star
    data$newx_Imp <- as.vector(newx_Imp)
    #if (option == "PMM") newx_Imp <- sapply(newx_Imp, find_closest, Y_obs)  # Comment out for regular MI
    
    famIDs <- unique(data$famID)
    u_j_df <- data.frame(famID = famIDs, u_j = u_j)
    
    ## Step 7
    data <- data |> 
      dplyr::mutate(newx_Imp_temp = newx_Imp) 
    
    if ("u_j" %in% names(data)) {
      data <- data |> dplyr::select(-u_j)
    }
    
    data <- left_join(data, u_j_df, by = "famID")
    
    data <- data |>
      dplyr::mutate(newx_Imp_temp = newx_Imp_temp + u_j * sigma_b_star)
  
    data <- data |>
      dplyr::mutate(newx_I = ifelse(is.na(newx), newx_Imp_temp, newx))

    if (i <= 2) {
      if (frailty.dist == "gamma") {
        fit <- survival::coxph(survival::Surv(time, status) ~ mgene + newx_I + frailty(famID, distribution = "gamma"), data = data)
      }
      else if (frailty.dist == "lognormal") {
        fit <- survival::coxph(survival::Surv(time, status) ~ mgene + newx_I + frailty(famID, distribution = "gaussian"), data = data)
      }
      baseline_haz <- survival::basehaz(fit, centered = FALSE)
      missing_times <- base::setdiff(data$time, baseline_haz$time)
      if (length(missing_times) > 0) {
        missing_haz <- data.frame(time = missing_times, hazard = 1)
        baseline_haz <- rbind(baseline_haz, missing_haz)
      }
      
      data <- data |>
        dplyr::left_join(baseline_haz, by = c("time" = "time")) |>
        dplyr::mutate(H0 = hazard) |>
        dplyr::select(-c("hazard"))
    }
    data_imp[[i]] <- data
  }
  
  ####################################################################################
  ################################# Analysis Step ####################################
  ####################################################################################
  data_imp <- data_imp[3:(M+2)]
  gamma_results <- list()
  for (i in 1:M) {
    gamma_results[[i]] <- penmodel(survival::Surv(time, status) ~ mgene + newx_I, cluster = "famID", 
                                   gvar = "mgene", design = "pop+", base.dist = "Weibull", 
                                   frailty.dist = frailty.dist, agemin = 18,
                                   data = data_imp[[i]], parms = true_value, robust = robust)
  }
  model_list <- gamma_results
  
  pooled_est <- Pooling(model = model_list, imputed_data = data_imp, robust = robust)
  return(pooled_est)
    }, error = function(e) {
      message("Error occurred: ", e$message, "on attempt ", attempts)
      if (attempts == max_attempts) {
        message("Max attempts reached, abandoning this dataset")
        return(NULL)
      }
    })
  }
}


