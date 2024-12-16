update.lmerMod <- function(object,theta,...) {
  if (missing(theta)) return(update.default(object,...))
  object2 <- object
  ## deep-copy the (only) reference-class slots ...
  object2@pp <- object@pp$copy()
  object2@resp <- object@resp$copy()
  object2@pp$setTheta(theta)
  dd <- as.function(object2)
  dval <- dd(theta)  ## update internal structures
  ## object2@resp$updateMu(object2@resp$mu)  ## ?? not helping/not necessary
  
  mkMerMod(environment(dd),
           opt=list(par=theta,fval=dval,conv=0),
           fr=model.frame(object2),
           reTrms=getME(object2,
                        c("Zt","Lambdat","Lind","theta",
                          "lower","flist","cnms","Gp"))
  )
}

tn <- function(object) {
  c(names(getME(object,"theta")),"sigma")
} ## Extract all variance components including the residual variance
## %%%%%%%%% CONFIRMED: sigma, so need to perform a Delta method  %%%%%%%%%%%%

outfun <- function(object,t) {
  newmod <- update(object,theta=t)
  av <- as.data.frame(VarCorr(newmod),
                      order="lower.tri")
  r <- setNames(av[,"sdcor"],tn(object))
  return(r)
} ## 

#waldVar2 <- function(object) {
  ## test for/warn if ML fit?
#  dd <- lme4:::devfun2(object,useSc=TRUE,signames=FALSE)
#  nvp <- length(attr(dd,"thopt"))+1 ## variance parameters (+1 for sigma)
#  pars <- attr(dd,"optimum")[seq(nvp)] ## var params come first
#  hh <- hessian(dd,pars)
  ## factor of 2: deviance -> negative log-likelihood
#  vv <- 2*solve(hh)
#  nn <- tn(object)
#  dimnames(vv) <- list(nn,nn)
#  return(vv)
#}

waldVar1 <- function(object) {
  pp <- getME(object,"theta") ## GLMM: unlist(getME(dd,c("theta","beta")))    
  Jmat <- jacobian(outfun,pp,object=object)
  dimnames(Jmat) <- list(tn(object),paste0("theta",seq(ncol(Jmat))))
  dd <- as.function(object) ## deviance/REMLcrit function
  hh <- 1/2*hessian(dd,pp)  ## ... calculate information matrix ...
  ## 1/2 = deviance to log-lik scale)
  ## vv <- solve(hh)        ## invert to get Wald vars of theta parameters
  ## m2 <- Jmat %*% vv %*% t(Jmat)  ## delta method
  ## slightly better linear algebra:
  m1 <- Jmat %*% solve(hh,t(Jmat))
  return(m1)
}

SE_VarComp <- function(fm1) {
  #fm1.ML <- refitML(fm1)
  wsd1 <- sqrt(diag(waldVar1(fm1)))
  return(wsd1)
}

#SE_VarComp(model_test)
#res <- list()
#for (i in 1:10) {
#  data <- simulated_dataset_500_MAR_Data_Gamma20[[i]]
#  kinship_mat <- with(data, kinship2::kinship(id = indID, dadid = fatherID, momid = motherID,
#                                              sex = gender))
#  kinship_mat_sparse <- Matrix::Matrix(kinship_mat, sparse = TRUE)
#  kinship_mat_sparse <- 2 * kinship_mat_sparse
#  kinship_mat <- as.matrix(kinship_mat_sparse)
  
#  Iden_mat <- diag(nrow = nrow(data)) # Identity matrix n*n
#  Iden_mat_sparse <- Matrix::Matrix(Iden_mat, sparse = TRUE)
  
 # Y_obs <- data$newx
 # newx <- Y_obs
  
  ## Step 2 
#  formula <- newx ~ log(log(ageonset)) + status * log(time) + mgene + (1 | indID) 
#  model_test <- lme4qtl::relmatLmer(newx ~ log(log(ageonset)) + status * log(time) + mgene + (1|indID), data = data, relmat = list(indID = kinship_mat), REML = TRUE)
#  X <- model.matrix(~ log(log(ageonset)) + status * log(time) + mgene  , data = data) # imputation model design matrix
#  V <- vcov(model_test)
  
#  betas <- as.vector(lme4::fixef(model_test)) 
#  sigma_g2 <- attr(lme4::VarCorr(model_test)$indID, "stddev")^2 # genetic variance
#  sigma_e2 <- attr(lme4::VarCorr(model_test), "sc")^2 # residual variance
#  p_pred <- ncol(X) - 1
  
  ## Total variance (This is not necessary as long as the beta is not adjusted...)
#  sigma2_hat <- sigma_e2 + sigma_g2
#  sigma_hat <- sqrt(sigma2_hat)
  
 # n_obs <- sum(!is.na(Y_obs))  # Number of observed rows
  
  ## Bootstrap to Estimate SEs of Variance Components
#  families <- unique(data$famID)
#  n_families <- length(families)
  
#  SEs <- SE_VarComp(model_test)
  
#  res <- c(res, SEs)
#}

