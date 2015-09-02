####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Functions takes the data:
# - est which is the estimates of the coefficients
# - eta is the eta for each logeps (within the suitable area)
# - logeps or log(ACER).
# - w is the weight factor
# - penalty opertunity if c is very close to 1
# - alpha is the penalty factor
#
# GumbelParam returns the four coefficients for the Gumbel model, when b and c is
# the input
# 
# mseGumbel2 returns the sum of square error (using 2 b,c since a and q
# is linearly depentent of these)
# 
# mseGumbel4 returns the sum of sqare error using all 4 coefficients
#
# GradmseGumbel4 returns the gradients of Gumbel in the point est
####################################################################################
GumbelParam<-function(est,eta,logeps,w){
  b=est[1]
  c=est[2]
  x=(eta-b)^c
  xbar=mean(x)
  y=logeps
  ybar=mean(y)
  a=-sum(w*(x-xbar)*(y-ybar))/sum(w*(x-xbar)^2)
  logq=ybar+a*xbar
  return(c(a,est,exp(logq)))
}

# est=c(b,c)
mseGumbel2<-function(est,eta,logeps,w,penalty=FALSE,alpha=0.05) {
  x=(eta-est[1])^est[2]
  xbar=mean(x)
  y=logeps
  ybar=mean(y)
  a=-sum(w*(x-xbar)*(y-ybar))/sum(w*(x-xbar)^2)
  logq=ybar+a*xbar
  if(penalty==FALSE){
    sol=sum(w*(y-logq+a*x)^2)
  }
  else if(penalty==TRUE){
    sol=exp(alpha*(est[2]^sign(est[2]-1)-1))*sum(w*(y-logq+a*x)^2)
  }
  return(sol)
}

# est=c(a,b,c,q)
mseGumbel4<-function(est,eta,logeps,w,penalty=FALSE,alpha=0.05) {
  if(penalty==FALSE){
    sol=sum(w*(logeps-log(est[4])+est[1]*(eta-est[2])^est[3])^2)
  }
  else if(penalty==TRUE){
    sol=exp(alpha*(est[3]^sign(est[3]-1)-1))*
      sum(w*(logeps-log(est[4])+est[1]*(eta-est[2])^est[3])^2)
  }
  return(sol)
}

GradmseGumbel4<-function(est,eta,logeps,w) {
  Ga=sum(2*w*(logeps-log(est[4])+est[1]*(eta-est[2])^est[3])*(eta-est[2])^est[3])
  Gb=sum(-2*w*(logeps-log(est[4])+est[1]*(eta-est[2])^est[3])*est[3]*est[1]*(eta-est[2])^(est[3]-1))
  Gc=sum(2*w*(logeps-log(est[4])+est[1]*(eta-est[2])^est[3])*est[1]*(eta-est[2])^est[3]*log((eta-est[2])))
  Gq=sum(-2*w*(logeps-log(est[4])+est[1]*(eta-est[2])^est[3])*1/est[4])
  return(c(Ga,Gb,Gc,Gq))
}