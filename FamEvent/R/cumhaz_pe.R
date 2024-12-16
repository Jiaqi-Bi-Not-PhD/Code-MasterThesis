cumhaz_pe <- function(a, b, tvc.parms, base.dist, base.parms) {
  H1 <- cumhaz(dist = base.dist, a, parms = base.parms)
  H2 <- cumhaz(dist = base.dist, b, parms = base.parms)
  return((H2-H1)*exp(tvc.parms[1]))
}