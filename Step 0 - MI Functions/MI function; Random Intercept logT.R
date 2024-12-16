#true_value2 <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
####################################################################################
####################################################################################
####################################################################################
################# Multiple Imputations WITHOUT kinship using log(T) ################
####################################################################################
####################################################################################
####################################################################################

####################################################################################
####################################################################################

MI_FamEvent_noK_logT <- function(data, M = 5, option = "General", frailty.dist = "gamma", true_value = true_value2,
                                 robust) {
  attempts <- 0
  max_attempts <- 1
  success <- FALSE
  while (attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
      ## Step 1 - Empirical estimates
      #imp_model <- lmer(newx ~ log(log(ageonset)) + status * log(time) + mgene + (1|famID) , data = data) # determined by user
      #X <- model.matrix( ~ log(log(ageonset)) + status * log(time) + mgene, data = data) # determined by user
      #imp_model <- lmer(newx ~ log(log(ageonset)) + status * log(time) + mgene + (1|famID), data = data)
      #X <- model.matrix(~ log(log(ageonset)) + status * log(time) + mgene , data = data)
      
      imp_model <- lme4::lmer(newx ~ log(ageonset) + status * log(time) + mgene + proband + (1|famID) , data = data)
      # Changes made here
      X <- model.matrix(~ log(ageonset) + status * log(time) + mgene + proband, data = data)
      
      betas <- lme4::fixef(imp_model)
      varcorr_list <- lme4::VarCorr(imp_model)
      sigma_e_hat <- attr(varcorr_list, "sc")
      sigma_b_hat <- attr(varcorr_list$famID, "stddev")
      V <- vcov(imp_model)

      Y_obs <- data$newx
      n_obs <- sum(!is.na(data$newx))
      n_fam_obs <- length(unique(data$famID[!is.na(data$newx)]))
      n_fam <- length(unique(data$famID))
      p_pred <- ncol(X)-1
      
      data_imp <- list()
      for (i in 1:M) {
        ## Step 2 - g ~ chi^2 n_obs - p_pred - n_fam_obs + 1
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
        
        famIDs <- unique(data$famID)
        u_j_df <- data.frame(famID = famIDs, u_j = u_j)
        
        ## Step 9
        data_update <- data |> 
          dplyr::mutate(newx_Imp_temp = newx_Imp) |>
          dplyr::left_join(u_j_df, by = "famID") |>
          dplyr::mutate(newx_Imp_temp = newx_Imp_temp + u_j * sigma_b_star)
        
        data_imp[[i]] <- data_update |>
          dplyr::mutate(newx_I = ifelse(is.na(newx), newx_Imp_temp, newx))
        
      }
      
      gamma_results <- list()
      for (i in 1:M) {
        gamma_results[[i]] <- penmodel(survival::Surv(time, status) ~ mgene + newx_I, cluster = "famID", 
                                       gvar = "mgene", design = "pop+", base.dist = "Weibull", 
                                       frailty.dist = frailty.dist, 
                                       agemin = 18, 
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
    }, warning = function(w) {
      message("Warning occurred: ", w$message, "on attempt ", attempts)
      if (attempts == max_attempts) {
        message("Max attempts reached, abandoning this dataset")
        return(NULL)
      }
    })
  }
}
