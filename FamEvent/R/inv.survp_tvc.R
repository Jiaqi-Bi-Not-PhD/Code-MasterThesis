inv.survp_tvc <- function(val, base.dist, parms, tvc.type, tvc.parms, alpha){
	out<-try(uniroot(survp.dist_tvc, lower=0,upper=100000, base.dist=base.dist, parms=parms, 
	                 tvc.age = val[2], tvc.type = tvc.type, tvc.parms = tvc.parms,
	                 alpha=alpha, xbeta=val[1], currentage=val[3], res=val[4])$root)
  if(is.null(attr(out,"class"))) return(out)
	else print(c(parms, val))
	}
