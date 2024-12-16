penmodelEM <- function(formula, cluster="famID", gvar="mgene", parms, cuts=NULL, 
                       data, design="pop", base.dist="Weibull", agemin=NULL, robust=FALSE, method="data", 
                       mode="dominant", q=0.02){
  options(na.action='na.pass')
  
  agemin <- attr(data, "agemin")
  if(is.null(agemin)){
    agemin <- 0
    warning("agemin = 0 was used or assign agemin to attr(data, \"agemin\").")
  }
  
  if(sum(data$time <=  agemin, na.rm = TRUE) > 0) cat("Individuals with time < agemin (", agemin,") were removed from the analysis.\n")
  data <- data[data$time >=  agemin, ]
  
  data <- data[data$time >  agemin, ]
  data$famID.byuser <- data[, cluster]
  m <- model.frame(formula, data)
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  
  if (!inherits(Y, "Surv")) stop("Response must be a survival object.")
  type <- attr(Y, "type")
  if (type == "counting") stop("start-stop type Surv objects are not supported.")
  if (type == "mright" || type == "mcounting") stop("multi-state survival is not supported.")

  X <- model.matrix(Terms, m)
  
  n <- nrow(X)
  nvar <- ncol(X)-1
  var.names <- colnames(X)[-1]
  if(nvar==1) X <- matrix(X[,-1])
  else X <- X[,-1]
  
  #number of parameters for baseline
  nbase <- ifelse(base.dist=="logBurr", 3, ifelse(base.dist=="piecewise", length(cuts)+1, 2)) 

  colnames(X) <- var.names
  vbeta <- parms[-c(1:nbase)]

  if(length(parms) != (nvar+nbase) ) stop("The size of parms is incorrect.")
  if(is.null(data$weight)) data$weight <- 1

  newdata <- carrierprobgeno(data, method=method, mode=mode, q=q)
 
  X0 <- X1 <- X
  
  X0[, gvar] <- 0
  X1[, gvar] <- 1
  
  est0 <- est <- parms
  dd <- lval0 <- lval <- 1
  i <- 0
  lval <- loglikem(X=X, X0=X0, X1=X1, Y=Y, theta=est, theta0=est0, cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=FALSE)
  cat("Iterations = ")
  while(dd>0.00001){
    i <- i+1
    est0 <- est
    lval0 <- lval
    nlm.est <- optim(est0, loglikem, X=X, X0=X0, X1=X1, Y=Y, theta0=est0, cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=FALSE)
    lval <- nlm.est$value
    est <- nlm.est$par
    dd <- abs(lval0-lval)
    #dd <- abs(sum(est-est0))    
    #print(c(i, dd, lval, est))
    cat(i, " ")
  }

  logLik <- -lval
  EST <- nlm.est$par
  H <- hessian(loglikem, nlm.est$par, X=X, X0=X0, X1=X1, Y=Y, theta0=nlm.est$par, cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=FALSE)
  Var <- try(solve(H), TRUE)
  
  if(!is.null(attr(Var,"class"))) stop("Model didn't converge.\n  Try again with different initial values")
  else{  
    bparms.name <- c("log(lambda)","log(rho)", "log(eta)")
    if(base.dist=="lognormal") bparms.name[1] <- "lambda" 
    else if(base.dist=="piecewise") bparms.name <- paste0("log(q", 1:nbase,")")
    
    parms.cov <- Var
    parms.se <- sqrt(diag(parms.cov))
    parms.cov.robust <- NULL
    parms.se.robust <- NULL
    names(EST) <- names(parms.se)  <- rownames(parms.cov) <- colnames(parms.cov) <- c(bparms.name[1:nbase], colnames(X))
    
  	if(robust){
  	  grad <- jacobian(loglikem, nlm.est$par, X=X, X0=X0, X1=X1, Y=Y, theta0=nlm.est$par, cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=TRUE)
  	  Jscore <- t(grad)%*%grad
  	  parms.cov.robust <- Var%*%(Jscore)%*%Var
  	  parms.se.robust <- sqrt(diag(parms.cov.robust))
  	  rownames(parms.cov.robust) <- colnames(parms.cov.robust) <- c(bparms.name[1:nbase], colnames(X))
  	
  	}
  }
  
  aic = 2*length(EST) - 2*logLik
  
  out <- list(estimates = EST, varcov = parms.cov, varcov.robust = parms.cov.robust, se = parms.se, se.robust = parms.se.robust, logLik = logLik, AIC = aic)  
  class(out) <- "penmodel"
  attr(out, "design") <- design
  attr(out, "base.dist") <- base.dist
  attr(out, "agemin") <- agemin
  attr(out, "cuts") <- cuts
  attr(out, "nbase") <- nbase
  attr(out, "data") <- data
  attr(out, "iterations") <- i
  attr(out, "robust") <- robust
  attr(out, "gvar") <- gvar
  attr(out, "formula") <- formula
  attr(out, "X") <- X
  attr(out, "Y") <- Y
  
  invisible(out)
  
  }
