cumhaz_co <- function(a, b, tvc.parms, base.dist, base.parms) {
  beta <- tvc.parms[1]
  eta <- tvc.parms[2]
  eta0 <- tvc.parms[3]
  
  integrand <- function(t) cumhaz(dist = base.dist, t, parms = base.parms) * exp((beta * exp(-(t - a) * eta) + eta0) * (t > a))
  int <- try(integrate(integrand, lower = a, upper = b), silent = TRUE)
  if(inherits(int, "try-error")) {
    warning(as.vector(int))
    integrated <- NA_real_
  } else {
    integrated <- int$value
  }
  return(integrated)
}