surv.dist_tvc <- function(t, base.dist, parms, xbeta, tvc.age, tvc.type, tvc.parms, alpha, res){

  ind1 <- which(t <= tvc.age)
  ind2 <- which(t > tvc.age)
  
  Haz <- rep(0, length(t))
  
  if(length(ind1)>0) Haz[ind1] <- cumhaz(base.dist, t[ind1], parms)*exp(xbeta+alpha)
  if(length(ind2)>0) {
  if(tvc.type == "PE") Htvc <- cumhaz_pe(tvc.age, t[ind2], tvc.parms, base.dist, parms)
  else if(tvc.type == "CO") Htvc <- cumhaz_co(tvc.age, t[ind2], tvc.parms, base.dist, parms)
  
  Haz[ind2] <- (cumhaz(base.dist, tvc.age, parms) + Htvc) * exp(xbeta+alpha) 
  }

  return(exp(-Haz) - res)
} 
  	