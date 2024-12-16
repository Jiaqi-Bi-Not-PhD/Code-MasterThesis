source("Calculating SE of Variance Components.R")

MI_FamEvent_K_logT <- function(data, M = 10, option = "General", frailty.dist = "gamma",
                               true_value, robust) {
  attempts <- 0 
  max_attempts <- 1
  success <- FALSE
  while(attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
      data_imp <- list()
      
      ## Step 1 - Kinship matrix
      kinship_mat <- with(data, kinship2::kinship(id = indID, dadid = fatherID, momid = motherID,
                                                  sex = gender))
      kinship_mat_sparse <- Matrix::Matrix(kinship_mat, sparse = TRUE)
      kinship_mat_sparse <- 2 * kinship_mat_sparse
      kinship_mat <- as.matrix(kinship_mat_sparse)
      
      Iden_mat <- diag(nrow = nrow(data)) # Identity matrix n*n
      Iden_mat_sparse <- Matrix::Matrix(Iden_mat, sparse = TRUE)
      
      Y_obs <- data$newx
      newx <- Y_obs
      
      ## Step 2 
      #formula <- newx ~ log(ageonset) + status * log(time) + mgene + (1 | indID) 
      #model_test <- lme4qtl::relmatLmer(newx ~ log(ageonset) + status * log(time) + mgene + (1|indID), data = data, relmat = list(indID = kinship_mat))
      #X <- model.matrix(~ log(ageonset) + status * log(time) + mgene  , data = data) # imputation model design matrix
      formula <- newx ~ log(ageonset) + status * log(time) + mgene  + proband + (1 | indID)  # changes made here
      model_test <- lme4qtl::relmatLmer(newx ~ log(ageonset) + status * log(time) + mgene  + proband + (1|indID) , data = data, relmat = list(indID = kinship_mat))
      # Changes made here
      X <- model.matrix(~ log(ageonset) + status * log(time) + mgene  + proband, data = data)
      
      V <- vcov(model_test)
      
      betas <- as.vector(lme4::fixef(model_test)) 
      sigma_g2 <- attr(lme4::VarCorr(model_test)$indID, "stddev")^2 # genetic variance
      sigma_e2 <- attr(lme4::VarCorr(model_test), "sc")^2 # residual variance

      p_pred <- ncol(X) - 1
      
      ## Total variance (This is not necessary...)
      #sigma2_hat <- sigma_e2 + sigma_g2
      #sigma_hat <- sqrt(sigma2_hat)
      
      n_obs <- sum(!is.na(Y_obs))  # Number of observed rows
      
      ## Bootstrap to Estimate SEs of Variance Components
      families <- unique(data$famID)
      n_families <- length(families)
      
      SEs <- SE_VarComp(model_test)
      
      se_sigma_g2 <- SEs["indID.(Intercept)"]
      #print(se_sigma_g2)
      se_sigma_e2 <- SEs["sigma"]
      #print(se_sigma_e2)
      ## df approximation
      df_sigma_g2 <- (sigma_g2) / (2 * se_sigma_g2^2)
      df_sigma_e2 <- (sigma_e2) / (2 * se_sigma_e2^2)
      
      #df_sigma_g2 <- (2*sigma_g2^2) / (se_sigma_g2^2)
      #df_sigma_e2 <- (2*sigma_e2^2) / (se_sigma_e2^2) #Either this or the above is the correct one
      
      df_sigma_g2 <- max(df_sigma_g2, 1)
      #print(df_sigma_g2)
      df_sigma_e2 <- max(df_sigma_e2, 1)
      #print(df_sigma_e2)
      
      newx <- data$newx
      
      ####################################################################################
      ################################# Imputation Step ##################################
      ####################################################################################
      for (m in 1:M) {
        ## Step 3 - Draw from chi-squared distributions
        s1 <- rchisq(1, df = df_sigma_g2)
        s2 <- rchisq(1, df = df_sigma_e2)
        
        ## Step 4 - Adjust variance estimates
        sigma_g_star <- sqrt( sigma_g2 * df_sigma_g2/s1 )
        sigma_e_star <- sqrt( sigma_e2 * df_sigma_e2/s2 ) 
        #sigma_star <- sqrt(sigma_g_star^2 + sigma_e_star^2)
        
        ## Step 5 - Draw random vector for coefficient adjustment
        w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
        
        ## Step 6 - Adjust fixed-effect coefficients
        #betastar <- betas + (sigma_star / sigma_hat) * (chol(V) %*% w_1)
        #betastar <- betas + (chol(V) %*% w_1)
        
        ## Step 7 - Compute predicted values
        #mu_star <- X %*% betastar
        
        ## Step 8 - Adjust covariance matrix
        
        
        Sigma <- sigma_g_star^2 * kinship_mat_sparse + sigma_e_star^2 * Iden_mat 
        V <- Matrix::solve(t(X) %*% Matrix::solve(Sigma) %*% X)
        ## %%%%%%%%%%% TODO Use Cholesky Decomposition here - Even faster %%%%%% ##
        betastar <- betas + chol(V) %*% w_1
        mu_star <- X %*% betastar
        
        cond_var_and_expect <- function(i) {
          y_minus_i <- newx[-i]
          mu_star_minus_i <- mu_star[-i]
          non_NA <- which(!is.na(y_minus_i))
          y_minus_i <- y_minus_i[non_NA]
          mu_star_minus_i <- mu_star_minus_i[non_NA]
          Sigma_nonNA <- Sigma[-i,-i][non_NA, non_NA]
          Sigmahat <- Sigma[i,-i][non_NA]
          
          ## Perform Cholesky decomposition for the non-NA subset 
          chol_Sigma_nonNA <- Matrix::Cholesky(Sigma_nonNA, LDL = FALSE)
          Sigmahat_Sigma_nonNA_inv <- Matrix::solve(chol_Sigma_nonNA, Sigmahat)
          conditional_var <- Sigma[i, i] - t(Sigmahat) %*% Sigmahat_Sigma_nonNA_inv
          conditional_Expect <- mu_star[i] + t(Sigmahat_Sigma_nonNA_inv) %*% (y_minus_i - mu_star_minus_i)
          
          return(list(var = conditional_var, expect = conditional_Expect))
        }
        
        results <- lapply(1:nrow(data), cond_var_and_expect)
        
        conditional_variances <- sapply(results, function(x) x$var)
        conditional_expectations <- sapply(results, function(x) x$expect)
        
        ## Step 8
        u <- rnorm(n = length(conditional_expectations), mean = 0, sd = 1)
        
        newx_star <- conditional_expectations + u * sqrt(conditional_variances)
        if (option == "PMM") newx_star <- sapply(newx_star, find_closest, Y_obs) # Comment out for regular MI
        data$newx_star <- newx_star
        data <- data |>
          mutate(newx_I = ifelse(!is.na(newx), newx, newx_star)) 
        
        data_imp[[m]] <- data
      }
      
      ####################################################################################
      ################################# Analysis Step ####################################
      ####################################################################################
      
      gamma_results <- list()
      for (i in 1:M) {
        gamma_results[[i]] <- quiet(penmodel(survival::Surv(time, status) ~ mgene + newx_I, cluster = "famID", 
                                             gvar = "mgene", 
                                             design = "pop+", base.dist = "Weibull", frailty.dist = frailty.dist, 
                                             agemin = 18, 
                                             data = data_imp[[i]],
                                             parms = true_value, robust = robust) )
      }
      model_list <- gamma_results
      pooled_est <- Pooling(model = model_list, imputed_data = data_imp, robust = robust)
      return(pooled_est)
    }, error = function(e) {
      message("Error occurred: ", e$message, " on attempt ", attempts)
      if (attempts == max_attempts) {
        message("Max attempts reached, abandoning this dataset")
        return(NULL)
      }
    })
  }
}

#time_start <- Sys.time()
#time_end <- Sys.time()
#time_operation <- time_end - time_start
#time_operation

