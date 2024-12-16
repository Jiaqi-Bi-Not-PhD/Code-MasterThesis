source("Archived/FamEvent/R/loglik_frailty.R")
source("FamEvent/R/cumhaz.R")
source("FamEvent/R/hazards.R")
source("FamEvent/R/gh.R")
source("FamEvent/R/penmodel.R")
source("FamEvent/R/dlaplace.R")
source("FamEvent/R/laplace.R")
source("Calculating SE of Variance Components.R")


brca1_prs_new <- brca1_prs |> mutate(missingindicator = ifelse(is.na(PRS), 1, 0))
model_BC <- glm(missingindicator ~ BC, family = binomial(), data = brca1_prs_new)
model_probandinfo <- glm(missingindicator ~ proband:currentage + proband, family = binomial(), data = brca1_prs_new)
model_mgene <- glm(missingindicator ~ mgeneI, family = binomial(), data = brca1_prs_new)
model_timeBC <- glm(missingindicator ~ timeBC, family = binomial(), data = brca1_prs_new)
summary(model_BC)
summary(model_probandinfo)
summary(model_mgene)
summary(model_timeBC)
###############################################################
###################### BRCA1 Data Application #################
####################### Complete Case Analysis ################
###############################################################
brca1_prs <- brca1_prs |> mutate(time = timeBC, 
                                 status = BC)
CCA_model_gamma <- penmodel(Surv(time, status) ~ mgeneI + PRS, cluster = "famID", 
                      gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
                      frailty.dist = "gamma", agemin = 18, 
                      data = brca1_prs, parms = c(1/41.41327,1,0,0, 1), robust = FALSE)
summary(CCA_model_gamma)
CCA_model_lognormal <- penmodel(Surv(time, status) ~ mgeneI + PRS, cluster = "famID", 
                                gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
                                frailty.dist = "lognormal", agemin = 18, 
                                data = brca1_prs, parms = c(1/41.41327,1,0,0, 1), robust = FALSE)
summary(CCA_model_lognormal)



true_value2 <- c(1/41.41327,1,0,0, 1)


###############################################################
###################### BRCA1 Data Application #################
####################### MI without K log(T) ###################
###############################################################
## Step 1 - Empirical estimates
MI_FamEvent_noK_logT <- function(data, M = 5, option = "General", frailty.dist = "gamma", true_value = true_value2,
                                 robust) {
  attempts <- 0
  max_attempts <- 1
  success <- FALSE
  while (attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
      ## Step 1 - Empirical estimates
      imp_model <- lmer(PRS ~ proband + log(time) * status + mgeneI + (1|famID)  , data = data) 
      X <- model.matrix( ~ proband + log(time) * status + mgeneI , data = data) 
      betas <- lme4::fixef(imp_model)
      varcorr_list <- lme4::VarCorr(imp_model)
      sigma_e_hat <- attr(varcorr_list, "sc")
      sigma_b_hat <- attr(varcorr_list$famID, "stddev")
      V <- vcov(imp_model)
      Y_obs <- data$PRS
      n_obs <- sum(!is.na(data$PRS))
      n_fam_obs <- length(unique(data$famID[!is.na(data$PRS)]))
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
        PRS_Imp <- X %*% betastar + u2i*sigma_e_star
        data$PRS_Imp <- as.vector(PRS_Imp)
        
        famIDs <- unique(data$famID)
        u_j_df <- data.frame(famID = famIDs, u_j = u_j)
        
        ## Step 9
        data_update <- data |> 
          dplyr::mutate(PRS_Imp_temp = PRS_Imp) |>
          dplyr::left_join(u_j_df, by = "famID") |>
          dplyr::mutate(PRS_Imp_temp = PRS_Imp_temp + u_j * sigma_b_star)
        
        data_imp[[i]] <- data_update |>
          dplyr::mutate(PRS_I = ifelse(is.na(PRS), PRS_Imp_temp, PRS))
        
      }
      
      gamma_results <- list()
      for (i in 1:M) {
        gamma_results[[i]] <- penmodel(survival::Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", 
                                       gvar = "mgeneI", design = "pop+", base.dist = "Weibull", 
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

BRCA1_noKlogT_gamma <- MI_FamEvent_noK_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "gamma", M=10)
# [[1]]
# A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
# 1 log.lambda    -4.52             0.00425          0.000164        0.00443  0          -4.65     -4.39 
# 2 log.rho        1.17             0.00176          0.00000861      0.00177  3.83e-150   1.09      1.26 
# 3 mgeneI          2.24             0.0602           0.00199         0.0624   6.96e- 19   1.75      2.73 
# 4 PRS_I          0.288            0.0189           0.00831         0.0281   8.90e-  2  -0.0449    0.622
# 5 log.kappa      1.11             0.232            0.00486         0.238    2.26e-  2   0.157     2.07 

# $`Pooled Var-Cov`
#             log.lambda       log.rho        mgeneI         PRS_I     log.kappa
# log.lambda  0.0044305287  0.0006626435 -0.011990757 -0.0001639542  0.0036745846
# log.rho     0.0006626435  0.0017720878  0.003442257  0.0004739397 -0.0087251476
# mgeneI      -0.0119907575  0.0034422567  0.062431407  0.0089224728 -0.0411338313
# PRS_I      -0.0001639542  0.0004739397  0.008922473  0.0280736935 -0.0009240991
# log.kappa   0.0036745846 -0.0087251476 -0.041133831 -0.0009240991  0.2377505054
BRCA1_noKlogT_lognormal <- MI_FamEvent_noK_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal", M=10)
# [[1]]
# A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.56             0.00511          0.000543        0.00571  1.99e-259  -4.71     -4.41 
# 2 log.rho        1.18             0.00180          0.00000286      0.00180  1.98e-148   1.09      1.26 
# 3 mgeneI          2.22             0.0599           0.00543         0.0659   2.73e- 17   1.72      2.73 
# 4 PRS_I          0.337            0.0187           0.0162          0.0364   8.57e-  2  -0.0497    0.724
# 5 log.kappa      1.18             0.234            0.00317         0.237    1.57e-  2   0.223     2.13 

# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI        PRS_I    log.kappa
# log.lambda  0.0057054232  0.0003024115 -0.01493816 0.0014274240  0.014069108
# log.rho     0.0003024115  0.0018047543  0.00344044 0.0004369887 -0.009276057
# mgeneI      -0.0149381552  0.0034404395  0.06587008 0.0068919102 -0.040217446
# PRS_I       0.0014274240  0.0004369887  0.00689191 0.0364333611  0.007589255
# log.kappa   0.0140691078 -0.0092760568 -0.04021745 0.0075892554  0.237084928


###############################################################
#################### BRCA1 Data Application ###################
#################### MI without K H0(T) #######################
###############################################################
MI_FamEvent_noK_H0T <- function(data, M = 5, option = "General", frailty.dist = "gamma", true_value = true_value2,
                                robust) {
  attempts <- 0
  max_attempts <- 1
  success <- FALSE
  while(attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
      ## Step X - Initial baseline cumulative hazard
      miss50_gamma_cca <- quiet(FamEvent::penmodel(survival::Surv(time, status) ~ mgeneI + PRS, 
                                                   cluster = "famID", gvar = "mgeneI",
                                                   design = "pop+", base.dist = "Weibull", 
                                                   frailty.dist = "none", 
                                                   agemin = 18, data = data,
                                                   parms = c(1/41.41327,1,0,0)))
      baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
      #print(baseline_gammafr) #####
      logalpha <- baseline_gammafr[1]
      loglambda <- baseline_gammafr[2]
      data <- data |>
        mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
      
      Y_obs <- data$PRS
      data_imp <- list()
      for (i in 1:(M+2)) {
        ## Step 1 - Empirical estimates
        # imp_model <- lmer(PRS ~ log(log(ageonset)) + status * H0 + mgeneI + (1|famID), data = data) 
        #X <- model.matrix( ~ log(log(ageonset)) + status * H0 + mgeneI, data = data)
        #imp_model <- lmer(PRS ~ log(log(ageonset)) + status * H0 + mgeneI + (1|famID), data = data)
        #X <- model.matrix(~ log(log(ageonset)) + status * H0 + mgeneI, data = data)
        
        imp_model <- lme4::lmer(PRS ~ proband + H0 * status + mgeneI + (1|famID) , data = data)
        # Changes made here
        X <- model.matrix(~ proband + H0 * status + mgeneI, data = data)
        
        betas <- lme4::fixef(imp_model)
        p_pred <- ncol(X)-1
        varcorr_list <- lme4::VarCorr(imp_model)
        sigma_e_hat <- attr(varcorr_list, "sc")
        sigma_b_hat <- attr(varcorr_list$famID, "stddev")
        V <- vcov(imp_model)
        data_obs <- data |> dplyr::filter(!is.na(PRS))
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
        PRS_Imp <- X %*% betastar + u2i*sigma_e_star
        data$PRS_Imp <- as.vector(PRS_Imp)
        #if (option == "PMM") PRS_Imp <- sapply(PRS_Imp, find_closest, Y_obs)  # Comment out for regular MI
        
        famIDs <- unique(data$famID)
        u_j_df <- data.frame(famID = famIDs, u_j = u_j)
        
        ## Step 7
        data <- data |> 
          dplyr::mutate(PRS_Imp_temp = PRS_Imp) 
        
        if ("u_j" %in% names(data)) {
          data <- data |> dplyr::select(-u_j)
        }
        
        data <- left_join(data, u_j_df, by = "famID")
        
        data <- data |>
          dplyr::mutate(PRS_Imp_temp = PRS_Imp_temp + u_j * sigma_b_star)
        
        data <- data |>
          dplyr::mutate(PRS_I = ifelse(is.na(PRS), PRS_Imp_temp, PRS))
        
        if (i <= 2) {
          if (frailty.dist == "gamma") {
            fit <- survival::coxph(survival::Surv(time, status) ~ mgeneI + PRS_I + frailty(famID, distribution = "gamma"), data = data)
          }
          else if (frailty.dist == "lognormal") {
            fit <- survival::coxph(survival::Surv(time, status) ~ mgeneI + PRS_I + frailty(famID, distribution = "gaussian"), data = data)
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
        gamma_results[[i]] <- penmodel(survival::Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", 
                                       gvar = "mgeneI", design = "pop+", base.dist = "Weibull", 
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


BRCA1_noKH0T_gamma <- MI_FamEvent_noK_H0T(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "gamma", M = 10)

# [[1]]
# A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.67            0.00130            0.000191      0.00151   8.58e-308   -4.75    -4.60 
# 2 log.rho        0.910           0.000814           0.0000165     0.000832  2.01e-181    0.853    0.966
# 3 mgeneI         2.22            0.0108             0.00647       0.0180    4.60e- 23    1.95     2.49 
# 4 PRS_I          0.144           0.00524            0.0256        0.0334    4.46e-  1   -0.253    0.541
# 5 log.kappa      1.11            0.0510             0.00181       0.0530    1.38e-  6    0.663    1.57 

# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI         PRS_I     log.kappa
# log.lambda  0.0015119685  0.0004660541 -0.002304388 -0.0005262026  0.0006230351
# log.rho     0.0004660541  0.0008322693  0.001342287  0.0006732604 -0.0017020301
# mgeneI     -0.0023043883  0.0013422866  0.017953431  0.0141412393 -0.0043958069
# PRS_I      -0.0005262026  0.0006732604  0.014141239  0.0333606367  0.0012628642
# log.kappa   0.0006230351 -0.0017020301 -0.004395807  0.0012628642  0.0529569948

BRCA1_noKH0T_lognormal <- MI_FamEvent_noK_H0T(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal", M = 10)

# [[1]]
# A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.72            0.00150            0.000571      0.00213   1.65e-101  -4.81     -4.62 
# 2 log.rho        0.920           0.000829           0.0000343     0.000867  9.79e-169   0.862     0.978
# 3 mgeneI         2.23            0.0106             0.00197       0.0128    1.29e- 54   2.01      2.46 
# 4 PRS_I          0.338           0.00540            0.0200        0.0274    6.12e-  2  -0.0183    0.694
# 5 log.kappa      1.06            0.0441             0.000712      0.0449    5.91e-  7   0.646     1.48 

# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI         PRS_I     log.kappa
# log.lambda  0.0021302707  0.0004758841 -0.002318791  0.0023178650  0.0030060558
# log.rho     0.0004758841  0.0008667426  0.001170275  0.0008484437 -0.0018170584
# mgeneI     -0.0023187909  0.0011702750  0.012770009  0.0032656665 -0.0045506741
# PRS_I       0.0023178650  0.0008484437  0.003265666  0.0274482496 -0.0007117185
# log.kappa   0.0030060558 -0.0018170584 -0.004550674 -0.0007117185  0.0448953596


###############################################################
###################### BRCA1 Data Application #################
####################### MI with K log(T) ######################
###############################################################
MI_FamEvent_K_logT <- function(data, M = 5, option = "General", frailty.dist = "gamma", true_value, 
                               robust) {
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
      
      Y_obs <- data$PRS
      
      #knots <- quantile(log(data$H0[data$status == 1]), prob = c(0.2, 0.4, 0.6, 0.8))
      
      #start_time <- Sys.time() # Starting time 
      
      ## Step 2 - empirical estimates
      model_test <- lme4qtl::relmatLmer(PRS ~ proband + status * log(time) + mgeneI  + (1|indID) + (1|famID), data = data, relmat = list(indID = kinship_mat), control = lmerControl(optimizer ="Nelder_Mead"))
      #summary(model_test)
      X <- model.matrix(~ proband + status * log(time) + mgeneI  , data = data) # imputation model design matrix
      V <- vcov(model_test)
      
      #sigmastar <- sigmahat*(sqrt((n_obs-p_pred)/g))
      #betas <- coef(model_test)
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
      
      PRS <- data$PRS
      
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
          y_minus_i <- PRS[-i]
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
        
        PRS_star <- conditional_expectations + u * sqrt(conditional_variances)
        if (option == "PMM") PRS_star <- sapply(PRS_star, find_closest, Y_obs) # Comment out for regular MI
        data$PRS_star <- PRS_star
        data <- data |>
          mutate(PRS_I = ifelse(!is.na(PRS), PRS, PRS_star)) 
        
        data_imp[[m]] <- data
      }
      
      ####################################################################################
      ################################# Analysis Step ####################################
      ####################################################################################
      
      gamma_results <- list()
      for (i in 1:M) {
        gamma_results[[i]] <- quiet(penmodel(survival::Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", 
                                             gvar = "mgeneI", 
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


BRCA2_K_LOGT_gamma <- MI_FamEvent_K_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "gamma", M=5)

# [[1]]
# A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.68            0.00130             0.000349      0.00172   2.35e-76   -4.76    -4.59 
# 2 log.rho        0.916           0.000817            0.000342      0.00123   1.42e-24    0.845    0.987
# 3 mgeneI         2.19            0.0106              0.0125        0.0256    1.73e- 8    1.84     2.54 
# 4 PRS_I          0.186           0.00568             0.305         0.372     7.77e- 1   -1.55     1.93 
# 5 log.kappa      1.08            0.0490              0.00132       0.0506    1.57e- 6    0.643    1.53 

# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI       PRS_I    log.kappa
# log.lambda  1.717539e-03  8.971586e-05 -0.003534895 -0.01114951  0.001109067
# log.rho     8.971586e-05  1.227137e-03  0.003429787  0.01220222 -0.002310984
# mgeneI     -3.534895e-03  3.429787e-03  0.025627868  0.07512227 -0.006606065
# PRS_I      -1.114951e-02  1.220222e-02  0.075122268  0.37178233 -0.018080979
# log.kappa   1.109067e-03 -2.310984e-03 -0.006606065 -0.01808098  0.050574562

BRCA2_K_LOGT_lognormal <- MI_FamEvent_K_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal", M=5)

# > BRCA2_K_LOGT_lognormal
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.73            0.00150             0.000971     0.00266    2.24e-28   -4.83    -4.62 
# 2 log.rho        0.924           0.000836            0.000118     0.000978   7.64e-70    0.862    0.986
# 3 mgeneI         2.25            0.0107              0.00664      0.0186     1.03e-13    1.97     2.54 
# 4 PRS_I          0.394           0.00498             0.148        0.183      4.09e- 1   -0.794    1.58 
# 5 log.kappa      1.04            0.0425              0.00314      0.0463     1.94e- 6    0.614    1.46 
# 
# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI       PRS_I    log.kappa
# log.lambda  0.0026617019  0.0007958518 -0.001883132 0.007711396  0.004306446
# log.rho     0.0007958518  0.0009781821  0.001307410 0.003982015 -0.001252676
# mgeneI     -0.0018831321  0.0013074096  0.018628135 0.030541357 -0.005646775
# PRS_I       0.0077113961  0.0039820152  0.030541357 0.182715974  0.005052835
# log.kappa   0.0043064465 -0.0012526757 -0.005646775 0.005052835  0.046287293




###############################################################
###################### BRCA1 Data Application #################
####################### MI with K H0(T) #######################
###############################################################
MI_FamEvent_K_H0T <- function(data, M = 5, option = "General", frailty.dist = "gamma", true_value = true_value2,
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
      
      ## Step X - Initial baseline cumulative hazard
      miss50_gamma_cca <- quiet(FamEvent::penmodel(survival::Surv(time, status) ~ mgeneI + PRS, 
                                                   cluster = "famID", gvar = "mgeneI",
                                                   design = "pop+", base.dist = "Weibull", 
                                                   frailty.dist = "none", 
                                                   agemin = 18, data = data,
                                                   parms = c(1/41.41327,1,0,0)))
      baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
      #print(baseline_gammafr) #####
      logalpha <- baseline_gammafr[1]
      loglambda <- baseline_gammafr[2]
      data <- data |>
        mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
      #fit <- survival::coxph(survival::Surv(time, status) ~ mgeneI + PRS, data = data)
      #baseline_haz <- survival::basehaz(fit, centered = FALSE)
      #data <- data |>
      #  dplyr::left_join(baseline_haz, by = c("time" = "time")) |>
      #  dplyr::mutate(H0 = hazard) |>
      #  dplyr::select(-c("hazard"))
      
      
      #H0 <- mice::nelsonaalen(data, time, status)
      #data$H0 <- H0
      
      
      
      
      Y_obs <- data$PRS # breslow
      
      #knots <- quantile(log(data$H0[data$status == 1]), prob = c(0.2, 0.4, 0.6, 0.8))
      
      #start_time <- Sys.time() # Starting time 
      ####################################################################################
      ################################# Imputation Step ##################################
      ####################################################################################
      
      for (m in 1:M) {
        ## Step 2 - empirical estimates
        #formula <- PRS ~ log(log(ageonset)) + status * H0 + mgeneI  + (1|indID)
        #model_test <- lme4qtl::relmatLmer(PRS ~ log(log(ageonset)) + status * H0 + mgeneI  + (1|indID), data = data, relmat = list(indID = kinship_mat))
        #X <- model.matrix(~ log(log(ageonset)) + status * H0 + mgeneI  , data = data) # imputation model design matrix
        
        formula <- PRS ~ proband + H0 * status + mgeneI  + (1 | indID)  # changes made here
        model_test <- lme4qtl::relmatLmer(PRS ~ proband + H0 * status + mgeneI + (1|indID) , data = data, relmat = list(indID = kinship_mat))
        # Changes made here
        X <- model.matrix(~ proband + H0 * status + mgeneI, data = data)
        
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
        
        PRS <- data$PRS
        
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
          y_minus_i <- PRS[-i]
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
        
        PRS_star <- conditional_expectations + u * sqrt(conditional_variances)
        #if (option == "PMM") PRS_star <- sapply(PRS_star, find_closest, Y_obs) # Comment out for regular MI
        data$PRS_star <- PRS_star
        data <- data |>
          mutate(PRS_I = ifelse(!is.na(PRS), PRS, PRS_star))
        
        data_imp[[m]] <- data
        
        if (m <= 2) {
          ## update baseline cumulative hazard
          
          if (frailty.dist == "gamma") {
            fit <- survival::coxph(survival::Surv(time, status) ~ mgeneI + PRS_I + frailty(famID, distribution = "gamma"), data = data)
          }
          else if (frailty.dist == "lognormal") {
            fit <- survival::coxph(survival::Surv(time, status) ~ mgeneI + PRS_I + frailty(famID, distribution = "gaussian"), data = data)
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
        gamma_results[[i]] <- quiet(penmodel(survival::Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", 
                                             gvar = "mgeneI", 
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

BRCA1_gamma_K_H0T <- MI_FamEvent_K_H0T(data = brca1_prs, frailty.dist = "gamma", robust = FALSE, M=10)

# [[1]]
# A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.66            0.00129             0.000930     0.00231   2.28e- 53   -4.76    -4.56 
# 2 log.rho        0.912           0.000815            0.000101     0.000925  1.95e-113    0.852    0.972
# 3 mgeneI         2.16            0.0103              0.00190      0.0124    2.32e- 53    1.94     2.37 
# 4 PRS_I          0.116           0.00469             0.0962       0.111     7.36e-  1   -0.636    0.867
# 5 log.kappa      1.10            0.0499              0.000169     0.0500    9.54e-  7    0.660    1.54 

# $`Pooled Var-Cov`
#             log.lambda       log.rho        mgeneI        PRS_I     log.kappa
# log.lambda  0.0023082874  0.0007931672 -0.0018804293 0.0094845536  0.0003482212
# log.rho     0.0007931672  0.0009251295  0.0009800761 0.0033033603 -0.0016451061
# mgeneI     -0.0018804293  0.0009800761  0.0123647992 0.0033401483 -0.0033017236
# PRS_I       0.0094845536  0.0033033603  0.0033401483 0.1105324592  0.0003615726
# log.kappa   0.0003482212 -0.0016451061 -0.0033017236 0.0003615726  0.0500374001

BRCA1_lognormal_K_H0T <- MI_FamEvent_K_H0T(data = brca1_prs, frailty.dist = "lognormal", robust = FALSE, M=5)

# > BRCA1_lognormal_K_H0T
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.73            0.00151            0.000428      0.00202   5.32e- 70   -4.82    -4.64 
# 2 log.rho        0.920           0.000831           0.0000448     0.000885  3.67e-136    0.861    0.978
# 3 mgeneI         2.16            0.0104             0.00444       0.0157    2.43e- 18    1.90     2.41 
# 4 PRS_I          0.101           0.00273            0.128         0.156     8.12e-  1   -1.02     1.22 
# 5 log.kappa      1.01            0.0421             0.00135       0.0438    1.40e-  6    0.603    1.42 
# 
# $`Pooled Var-Cov`
#             log.lambda       log.rho        mgeneI       PRS_I    log.kappa
# log.lambda  0.0020247021  0.0005627869 -0.0017932137 0.006733853  0.002435869
# log.rho     0.0005627869  0.0008850863  0.0008990831 0.001680047 -0.001927216
# mgeneI     -0.0017932137  0.0008990831  0.0156947923 0.017075916 -0.000944295
# PRS_I       0.0067338527  0.0016800469  0.0170759161 0.156272055  0.006115651
# log.kappa   0.0024358692 -0.0019272160 -0.0009442950 0.006115651  0.043770553




