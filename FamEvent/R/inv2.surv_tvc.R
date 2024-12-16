inv2.surv_tvc <- function(val, base.dist, parms, tvc.type, tvc.parms){
  
  out<-try(uniroot(surv.dist_tvc, lower=0,upper=100000, base.dist=base.dist, parms=parms, 
                   alpha=val[3], xbeta=val[1], tvc.age=val[2], tvc.type=tvc.type, 
                   tvc.parms=tvc.parms, res=val[4])$root)
  if(is.null(attr(out,"class"))) return(out)
	else print(c(parms, val))
	}
