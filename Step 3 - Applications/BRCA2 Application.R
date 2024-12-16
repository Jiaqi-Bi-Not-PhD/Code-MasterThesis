source("Archived/FamEvent/R/loglik_frailty.R")
source("FamEvent/R/cumhaz.R")
source("FamEvent/R/hazards.R")
source("FamEvent/R/gh.R")
source("FamEvent/R/penmodel.R")
source("FamEvent/R/dlaplace.R")
source("FamEvent/R/laplace.R")

###############################################################
###################### BRCA1 Data Application #################
####################### Complete Case Analysis ################
###############################################################
brca1_prs <- brca1_prs |> mutate(time = time, 
                                 status = status)
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
MI_FamEvent_noK_logT <- function(data, M = 10, option = "General", frailty.dist = "gamma", true_value = true_value2,
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

BRCA1_noKlogT_gamma <- MI_FamEvent_noK_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "gamma",
                                            M=5)

# > BRCA1_noKlogT_gamma
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda   -4.68              0.00165          0.000290        0.00200  2.29e-121   -4.77    -4.59 
# 2 log.rho       0.994             0.00106          0.00000593      0.00106  7.36e-166    0.930    1.06 
# 3 mgeneI        2.04              0.0154           0.00808         0.0251   7.84e- 13    1.71     2.36 
# 4 PRS_I         0.0399            0.00679          0.0212          0.0322   8.31e-  1   -0.394    0.474
# 5 log.kappa     0.967             0.0579           0.000313        0.0583   6.43e-  5    0.494    1.44 
# 
# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI         PRS_I    log.kappa
# log.lambda  0.0020032021  0.0004955123 -0.004514426 -0.0036315128  0.001652904
# log.rho     0.0004955123  0.0010643729  0.001398511  0.0004644649 -0.002215662
# mgeneI     -0.0045144258  0.0013985111  0.025127881  0.0197170315 -0.007875681
# PRS_I      -0.0036315128  0.0004644649  0.019717032  0.0322479647 -0.004143451
# log.kappa   0.0016529042 -0.0022156622 -0.007875681 -0.0041434509  0.058297886


BRCA1_noKlogT_lognormal <- MI_FamEvent_noK_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal", M=5)

# > BRCA1_noKlogT_lognormal
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
# 1 log.lambda   -4.73              0.00191           0.0000465      0.00197  0           -4.82    -4.64 
# 2 log.rho       0.996             0.00107           0.0000445      0.00113  1.19e-135    0.931    1.06 
# 3 mgeneI        1.98              0.0145            0.00300        0.0181   3.05e- 26    1.71     2.25 
# 4 PRS_I         0.0159            0.00682           0.0428         0.0582   9.50e-  1   -0.603    0.635
# 5 log.kappa     1.01              0.0538            0.00182        0.0560   1.98e-  5    0.550    1.48 
# 
# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI         PRS_I    log.kappa
# log.lambda  0.0019704627  0.0004854518 -0.002614141  0.0005778717  0.003751444
# log.rho     0.0004854518  0.0011266220  0.001547467  0.0017190696 -0.002682012
# mgeneI     -0.0026141412  0.0015474670  0.018092164  0.0160533059 -0.008352248
# PRS_I       0.0005778717  0.0017190696  0.016053306  0.0581607569 -0.010713816
# log.kappa   0.0037514440 -0.0026820115 -0.008352248 -0.0107138161  0.056036237



###############################################################
#################### BRCA1 Data Application ###################
#################### MI without K H0(T) #######################
###############################################################
MI_FamEvent_noK_H0T <- function(data, M = 10, option = "General", frailty.dist = "gamma", true_value = true_value2,
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

BRCA1_noKH0T_gamma <- MI_FamEvent_noK_H0T(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "gamma",
                                          M=5)

# > BRCA1_noKH0T_gamma
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda   -4.67              0.00155           0.000411       0.00204  2.26e- 74   -4.76    -4.58 
# 2 log.rho       0.998             0.00106           0.0000533      0.00112  2.34e-128    0.932    1.06 
# 3 mgeneI        2.03              0.0144            0.0207         0.0393   1.40e-  6    1.59     2.48 
# 4 PRS_I         0.0804            0.00695           0.0902         0.115    8.23e-  1   -0.832    0.993
# 5 log.kappa     0.933             0.0567            0.00512        0.0628   2.30e-  4    0.440    1.43 
# 
# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI        PRS_I    log.kappa
# log.lambda  0.0020441456  0.0003926553 -0.005758928 -0.007237546  0.002876749
# log.rho     0.0003926553  0.0011227305  0.002222972  0.002368874 -0.002799188
# mgeneI     -0.0057589282  0.0022229717  0.039262341  0.054678996 -0.017123419
# PRS_I      -0.0072375459  0.0023688742  0.054678996  0.115246368 -0.024538263
# log.kappa   0.0028767491 -0.0027991879 -0.017123419 -0.024538263  0.062805988

BRCA1_noKH0T_lognormal <- MI_FamEvent_noK_H0T(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal", M=5)

# > BRCA1_noKH0T_lognormal
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda   -4.74              0.00197          0.000178        0.00219  2.00e-254   -4.83    -4.65 
# 2 log.rho       0.996             0.00107          0.00000462      0.00108  1.02e-164    0.932    1.06 
# 3 mgeneI        2.02              0.0150           0.00436         0.0202   2.01e- 20    1.73     2.30 
# 4 PRS_I         0.0515            0.00773          0.0119          0.0221   7.36e-  1   -0.282    0.386
# 5 log.kappa     1.01              0.0535           0.000944        0.0547   1.62e-  5    0.553    1.47 
# 
# $`Pooled Var-Cov`
#             log.lambda       log.rho       mgeneI         PRS_I    log.kappa
# log.lambda  0.0021870901  0.0004398707 -0.003942784 -0.0022806259  0.004116570
# log.rho     0.0004398707  0.0010797861  0.001275160  0.0003502284 -0.002426047
# mgeneI     -0.0039427842  0.0012751596  0.020245680  0.0124097582 -0.006990336
# PRS_I      -0.0022806259  0.0003502284  0.012409758  0.0220541980 -0.002917155
# log.kappa   0.0041165704 -0.0024260474 -0.006990336 -0.0029171551  0.054655386



###############################################################
###################### BRCA1 Data Application #################
####################### MI with K log(T) ######################
###############################################################
MI_FamEvent_K_logT <- function(data, M = 10, option = "General", frailty.dist = "gamma", true_value, 
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

BRCA2_K_LOGT_gamma <- MI_FamEvent_K_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "gamma",
                                         M=5)

# > BRCA2_K_LOGT_gamma
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.68             0.00160           0.0000915      0.00171  0           -4.76    -4.60 
# 2 log.rho        0.998            0.00106           0.0000436      0.00111  8.79e-138    0.933    1.06 
# 3 mgeneI         2.11             0.0152            0.00755        0.0243   7.53e- 14    1.80     2.43 
# 4 PRS_I          0.217            0.00770           0.0299         0.0435   3.40e-  1   -0.298    0.731
# 5 log.kappa      0.963            0.0568            0.00184        0.0590   7.83e-  5    0.487    1.44 
# 
# $`Pooled Var-Cov`
#             log.lambda      log.rho       mgeneI        PRS_I    log.kappa
# log.lambda  0.001706938  0.000490538 -0.003231736 -0.002158819  0.001422738
# log.rho     0.000490538  0.001108387  0.001840013  0.001437964 -0.002460405
# mgeneI     -0.003231736  0.001840013  0.024272553  0.022611511 -0.008433104
# PRS_I      -0.002158819  0.001437964  0.022611511  0.043531300 -0.005761341
# log.kappa   0.001422738 -0.002460405 -0.008433104 -0.005761341  0.059045472

BRCA2_K_LOGT_lognormal <- MI_FamEvent_K_logT(brca1_prs, true_value = true_value2, robust = FALSE, frailty.dist = "lognormal", M=5)

# > BRCA2_K_LOGT_lognormal
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.` `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>      <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.73             0.00186           0.000167       0.00206  5.16e-261   -4.82    -4.64 
# 2 log.rho        0.999            0.00107           0.0000332      0.00111  2.09e-146    0.934    1.06 
# 3 mgeneI         2.07             0.0146            0.00815        0.0244   1.33e- 12    1.74     2.39 
# 4 PRS_I          0.157            0.00763           0.0344         0.0489   5.08e-  1   -0.396    0.710
# 5 log.kappa      1.01             0.0531            0.00105        0.0544   1.57e-  5    0.553    1.47 
# 
# $`Pooled Var-Cov`
#               log.lambda       log.rho       mgeneI        PRS_I    log.kappa
# log.lambda  0.0020624750  0.0004794677 -0.002262738  0.000522197  0.004049115
# log.rho     0.0004794677  0.0011135507  0.001765066  0.001290464 -0.002544249
# mgeneI     -0.0022627379  0.0017650657  0.024354964  0.023817826 -0.008730721
# PRS_I       0.0005221970  0.0012904644  0.023817826  0.048904901 -0.005522745
# log.kappa   0.0040491150 -0.0025442489 -0.008730721 -0.005522745  0.054405959




###############################################################
###################### BRCA1 Data Application #################
####################### MI with K H0(T) #######################
###############################################################
MI_FamEvent_K_H0T <- function(data, M = 10, option = "General", frailty.dist = "gamma", true_value = true_value2,
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
        model_test <- lme4qtl::relmatLmer(PRS ~ proband + H0 * status + mgeneI + (1|indID) , data = data, relmat = list(indID = kinship_mat), control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
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

BRCA1_gamma_K_H0T <- MI_FamEvent_K_H0T(data = brca1_prs, frailty.dist = "gamma", robust = FALSE, M=5)

# > BRCA1_gamma_K_H0T
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.`  `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>       <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.63             0.00149             0.0129       0.0170  0.000000726   -4.97     -4.29
# 2 log.rho        1.03             0.00107             0.00533      0.00746 0.0000477      0.810     1.25
# 3 mgeneI         2.05             0.0147              0.0430       0.0663  0.000138       1.43      2.67
# 4 PRS_I          0.244            0.00660             0.204        0.251   0.652         -1.16      1.65
# 5 log.kappa      0.938            0.0557              0.00420      0.0608  0.000162       0.453     1.42
# 
# $`Pooled Var-Cov`
#             log.lambda      log.rho       mgeneI       PRS_I    log.kappa
# log.lambda  0.017011211  0.009936632  0.005437379  0.03542657 -0.005603587
# log.rho     0.009936632  0.007458785  0.011397380  0.03243929 -0.007467211
# mgeneI      0.005437379  0.011397380  0.066317497  0.10859626 -0.018618952
# PRS_I       0.035426567  0.032439289  0.108596261  0.25109462 -0.035035091
# log.kappa  -0.005603587 -0.007467211 -0.018618952 -0.03503509  0.060763730

BRCA1_lognormal_K_H0T <- MI_FamEvent_K_H0T(data = brca1_prs, frailty.dist = "lognormal", robust = FALSE,
                                           M=5)

# > BRCA1_lognormal_K_H0T
# [[1]]
# # A tibble: 5 × 8
# Parameters Estimates `Within Imp. Var.` `Between Imp. Var.` `Total Var.`   `p values` Lower.CI Upper.CI
# <chr>          <dbl>              <dbl>               <dbl>        <dbl>        <dbl>    <dbl>    <dbl>
#   1 log.lambda    -4.70             0.00179             0.0140       0.0186  0.000000661    -5.05     -4.34
# 2 log.rho        1.03             0.00108             0.00493      0.00699 0.0000327       0.823     1.24
# 3 mgeneI         2.13             0.0146              0.0153       0.0329  0.0000000360    1.73      2.52
# 4 PRS_I          0.349            0.00685             0.150        0.187   0.465          -0.846     1.54
# 5 log.kappa      1.01             0.0532              0.00515      0.0594  0.0000412       0.535     1.49
# 
# $`Pooled Var-Cov`
#           log.lambda      log.rho       mgeneI       PRS_I     log.kappa
# log.lambda 0.018602958 0.0100330126  0.001495647 0.039327638  0.0109579742
# log.rho    0.010033013 0.0069949252  0.006305629 0.029192739  0.0007215737
# mgeneI     0.001495647 0.0063056287  0.032941125 0.051729170 -0.0108728450
# PRS_I      0.039327638 0.0291927391  0.051729170 0.187244100  0.0021998784
# log.kappa  0.010957974 0.0007215737 -0.010872845 0.002199878  0.0594037492

