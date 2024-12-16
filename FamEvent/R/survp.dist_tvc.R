survp.dist_tvc <- function(t, base.dist, currentage, parms, xbeta, 
                           tvc.age, tvc.type, tvc.parms, alpha, res){
  
  A = surv.dist_tvc(t=t, base.dist=base.dist, parms=parms, xbeta=xbeta, 
                    tvc.age=tvc.age, tvc.type=tvc.type, tvc.parms = tvc.parms, alpha=alpha, res=0)
  B = surv.dist_tvc(t=currentage, base.dist=base.dist, parms=parms, xbeta=xbeta, 
                    tvc.age=tvc.age, tvc.type=tvc.type, tvc.parms = tvc.parms, alpha=alpha, res=0)
  return((A-B)/(1-B)-res)  

} 
  	