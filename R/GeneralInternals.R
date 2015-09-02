####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Functions takes the data:
# - est which is the estimates of the coefficients
# - eta is the eta for each logeps (within the suitable area)
# - logeps or log(ACER).
# - w is the weight factor
#
# GeneralParam(f/w) returns a the coefficients for a GEV model.
# 
# mseGeneral3(f/w) returns the sum of square error (using 3 a~,b,c since xi and q
# is linearly depentent of these)
# 
# mseGeneral5(f/w) returns the sum of sqare error using all 5 coefficients
#
# GradmseGeneral5(f/w) returns the gradients of GEV in the point est
####################################################################################

# Input est is a~ ,b ,c (from GEV) and output 
# a,b,c,q,xi (from GEV), where a~=a*xi (Frechet)
GeneralParamf<-function(est,eta,logeps,w){
  temp=1+est[1]*(eta-est[2])^est[3]
  if(sum(is.na(est))>1||min(temp)<0){return(rep(NA,5))}
  x<-log(temp)
  xbar<-mean(x)
  y<-logeps
  ybar<-mean(y)
  xi<--sum(w*(x-xbar)^2)/sum(w*(x-xbar)*(y-ybar))
  logq<-ybar+xbar/xi
  return(c(est[1]/xi,est[2:3],exp(logq),xi))
}
# (Weibull)
GeneralParamw<-function(est,eta,logeps,w){
  temp=1+est[1]*(eta-est[2])
  if(sum(is.na(est))>1||min(temp)<0){return(rep(NA,4))}
  x<-log(temp)
  xbar<-mean(x)
  y<-logeps
  ybar<-mean(y)
  xi<--sum(w*(x-xbar)^2)/sum(w*(x-xbar)*(y-ybar))
  logq<-ybar+xbar/xi
  return(c(est[1]/xi,est[2],exp(logq),xi))
}

# est=c(a~,b,c), a~=a*xi (Frechet)
mseGeneral3f<-function(est,eta,logeps,w){
  temp=1+est[1]*(eta-est[2])^est[3]
  if(sum(is.na(est))>0||min(temp)<0){return(Inf)}
  x<-log(temp) 
  xbar<-mean(x)
  y<-logeps
  ybar<-mean(y)
  xi<--sum(w*(x-xbar)^2)/sum(w*(x-xbar)*(y-ybar))
  logq<-ybar+xbar/xi
  return(sum(w*(y-logq+x/xi)^2))
}
#est=c(a~,b), a~=a*xi (Weibull)
mseGeneral3w<-function(est,eta,logeps,w){
  temp=1+est[1]*(eta-est[2])
  if(sum(is.na(est))>0||min(temp)<0){return(Inf)}
  x<-log(temp) 
  xbar<-mean(x)
  y<-logeps
  ybar<-mean(y)
  xi<--sum(w*(x-xbar)^2)/sum(w*(x-xbar)*(y-ybar))
  logq<-ybar+xbar/xi
  return(sum(w*(y-logq+x/xi)^2))
}

#est=c(a,b,c,q,xi), a~=a*xi (Frechet)
mseGeneral5f<-function(est,eta,logeps,w){
  temp=1+est[1]*est[5]*(eta-est[2])^est[3]
  if(sum(is.na(est))>0||min(temp)<0){return(Inf)}
  return(sum(w*(logeps-log(est[4])+log(temp)/est[5])^2))
}
#est=c(a,b,q,xi), a~=a*xi (Weibull)
mseGeneral5w<-function(est,eta,logeps,w){
  temp=1+est[1]*est[4]*(eta-est[2])
  if(sum(is.na(est))>0||min(temp)<0){return(Inf)}
  return(sum(w*(logeps-log(est[3])+log(temp)/est[4])^2))
}

#gradient of mse with 5 parameters (a,b,c,q,xi) where a~=a*xi (Frechet)
GradmseGeneral5f<-function(est,eta,logeps,w) {
  x<-1+est[1]*est[5]*(eta-est[2])^est[3]
  FX<-2*w*(logeps-log(est[4])+log(x)/est[5])
  Ga<-sum(FX*(eta-est[2])^est[3]/x)
  Gb<-sum(-FX*est[1]*est[3]*(eta-est[2])^(est[3]-1)/x)
  Gc<-sum(FX*est[1]*(eta-est[2])^est[3]*log(eta-est[2])/x)
  Gq<-sum(-FX/est[4])
  Gxi<-sum(FX*((est[1]/est[5]*(eta-est[2])^est[3])/x-log(x)/(est[5]^2)))
  return(c(Ga,Gb,Gc,Gq,Gxi))       
}
#Fow Weibull
GradmseGeneral5w<-function(est,eta,logeps,w) {
  x<-1+est[1]*est[4]*(eta-est[2])
  Fx<-2*w*(logeps-log(est[3])+log(x)/est[4])
  Ga<-sum(Fx*(eta-est[2])/x)
  Gb<-sum(-Fx*est[1]/x)
  Gq<-sum(-Fx/est[3])
  Gxi<-sum(Fx*((est[1]/est[4]*(eta-est[2]))/x-log(x)/(est[4]^2)))
  return(c(Ga,Gb,Gq,Gxi))       
}

# d=1e-4
# numbericGrad<-function(est,eta,logeps,w,d){
#  Ga<-(mseGeneral5(est+c(d,0,0,0,0),eta,logeps,w)-mseGeneral5(est,eta,logeps,w))/d
#  Gb<-(mseGeneral5(est+c(0,d,0,0,0),eta,logeps,w)-mseGeneral5(est,eta,logeps,w))/d
#  Gc<-(mseGeneral5(est+c(0,0,d,0,0),eta,logeps,w)-mseGeneral5(est,eta,logeps,w))/d
#  Gq<-(mseGeneral5(est+c(0,0,0,d,0),eta,logeps,w)-mseGeneral5(est,eta,logeps,w))/d
#  Gxi<-(mseGeneral5(est+c(0,0,0,0,d),eta,logeps,w)-mseGeneral5(est,eta,logeps,w))/d
#  return(c(Ga,Gb,Gc,Gq,Gxi)) 
# }