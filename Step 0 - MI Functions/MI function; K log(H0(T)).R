#true_value2 <- c(exp(-4.71), exp(0.804), 2.13, 1, exp(1.15))
source("Suppress Cat.R")
source("Calculating SE of Variance Components.R")
####################################################################################
####################################################################################
####################################################################################
################# Multiple Imputations WITH kinship using log(H_0(T)) ##############
####################################################################################
####################################################################################
####################################################################################
MI_FamEvent_K_H0T <- function(data, M = 10, option = "General", frailty.dist = "gamma", true_value,
                              robust = FALSE) {
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
      
      ## Initial baseline cumulative hazard
      miss50_gamma_cca <- quiet(FamEvent::penmodel(survival::Surv(time, status) ~ mgene + newx, 
                                                   cluster = "famID", gvar = "mgene",
                                                   design = "pop+", base.dist = "Weibull", frailty.dist = "none", 
                                                   agemin = 18, data = data,
                                                   parms = c(exp(-4.71), exp(0.804), 2.13, 1)))
      baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
      logalpha <- baseline_gammafr[1]
      loglambda <- baseline_gammafr[2]
      data <- data |>
        mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
      
      Y_obs <- data$newx 
      
      ####################################################################################
      ################################# Imputation Step ##################################
      ####################################################################################
      
      for (m in 1:M) {
        ## Step 2 - empirical estimates
        #formula <- newx ~ log(log(ageonset)) + status * H0 + mgene  + (1|indID)
        #model_test <- lme4qtl::relmatLmer(newx ~ log(log(ageonset)) + status * H0 + mgene  + (1|indID), data = data, relmat = list(indID = kinship_mat))
        #X <- model.matrix(~ log(log(ageonset)) + status * H0 + mgene  , data = data) # imputation model design matrix
        
        formula <- newx ~ log(ageonset) + status * H0 + mgene  + proband + (1 | indID)  # changes made here
        model_test <- lme4qtl::relmatLmer(newx ~ log(ageonset) + status * H0 + mgene  + proband + (1|indID) , data = data, relmat = list(indID = kinship_mat))
        # Changes made here
        X <- model.matrix(~ log(ageonset) + status * H0 + mgene + proband, data = data)
        
        V  <- vcov(model_test)
        
        betas <- as.vector(summary(model_test)$coefficients[,1]) # beta coefficients
        sigma_g2 <- attr(lme4::VarCorr(model_test)$indID, "stddev")^2 # genetic variance
        sigma_e2 <- attr(lme4::VarCorr(model_test), "sc")^2 # residual variance
        p_pred <- ncol(X) - 1
        
        ## Total variance (This is not necessary as long as the beta is not adjusted...)
        sigma2_hat <- sigma_e2 + sigma_g2
        sigma_hat <- sqrt(sigma2_hat)
        
        n_obs <- sum(!is.na(Y_obs))
        
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
        df_sigma_e2 <- max(df_sigma_e2, 1)
        
        ## Draw from chi-squared distributions
        s1 <- rchisq(1, df = df_sigma_g2)
        s2 <- rchisq(1, df = df_sigma_e2)
        
        ## Adjust variance estimates
        sigma_g_star <- sqrt( sigma_g2 * df_sigma_g2/s1 )
        sigma_e_star <- sqrt( sigma_e2 * df_sigma_e2/s2 ) 
        
        newx <- data$newx
        
        ## Draw random vector for coefficient adjustment
        w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
        
        #betastar <- betas + (sigma_star / sigma_hat) * (chol(V) %*% w_1)
        ## Adjust covariance matrix
        #Sigma_orig <- sigma_g2 * kinship_mat_sparse + sigma_e2 * Iden_mat
        Sigma <- sigma_g_star^2 * kinship_mat_sparse + sigma_e_star^2 * Iden_mat
        V <- Matrix::solve(t(X) %*% Matrix::solve(Sigma) %*% X)
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
        
        ## Step 9
        u <- rnorm(n = length(conditional_expectations), mean = 0, sd = 1)
        
        newx_star <- conditional_expectations + u * sqrt(conditional_variances)
        #if (option == "PMM") newx_star <- sapply(newx_star, find_closest, Y_obs) # Comment out for regular MI
        data$newx_star <- newx_star
        data <- data |>
          mutate(newx_I = ifelse(!is.na(newx), newx, newx_star))
        
        data_imp[[m]] <- data
        
        if (m <= 2) {
          ## update baseline cumulative hazard
          
          if (frailty.dist == "gamma") {
            fit <- survival::coxph(survival::Surv(time, status) ~ mgene + newx_I + frailty(famID, distribution = "gamma"), data = data)
          }
          else if (frailty.dist == "lognormal") {
            fit <- survival::coxph(survival::Surv(time, status) ~ mgene + newx_I + frailty(famID, distribution = "gaussian"), data = data)
          }
          baseline_haz <- survival::basehaz(fit, centered = FALSE)
          data <- data |>
            dplyr::left_join(baseline_haz, by = c("time" = "time")) |>
            dplyr::mutate(H0 = hazard) |>
            dplyr::select(-c("hazard"))
          data <- data |> dplyr::mutate(H0 = ifelse(is.na(H0), 0, H0))
        
        }
      }
      
      ####################################################################################
      ################################# Analysis Step ####################################
      ####################################################################################
      ## Analysis Gamma
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
#MI_FamEvent_K_H0T(data = data, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal")
