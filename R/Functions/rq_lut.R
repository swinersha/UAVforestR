#-----------------------------------------------
# 
# A function that produces a height to radius lookup table from tau = 0.01 to tau=1
#
# Partner functions are included for using the lookup table to convert 
# from height to radius and vice versa
#
# Tom Swinfield
# 17-02-17
#
#-----------------------------------------------

rq_lut<-function(x, y, log){
  rq_tau<-function(x, y, tau, log){
    if(log){
     x<-log(x); y<-log(y) 
    }
      
    fit = quantreg::rq(y ~ x, tau = tau )
    a <- coef(fit)[1]
    b <- coef(fit)[2]
  
    params<-data.frame(tau=tau, a=a, b=b)
    rownames(params)<-NULL
    return(params)
  }
  lut<-lapply((1:100)/100, function(tau) rq_tau(x, y, tau, log))
  lut<-do.call(rbind, lut)
  return(lut)
}

allom_lookup<-function(x, lut, tau, antilog) 
{
  rad<-x
  a<-lut[tau,'a']
  b<-lut[tau,'b']
  if(antilog)
    rad<-(exp(a) * x^b) 
  else
    rad<-a + x*b 
  return(rad)
}