####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X an ACER class 
# - level is the confidence level
# - tdist.CI is logical. If TRUE t-distribudtion is used to 
#   calculate confidence intervall, if FALSE poisson is used.
#   For data with only one realization t-distribudtion is impossible
#   and poisson is done.
#
# Return a updated version of X (ACER class)
# level: Confidence Intervals of the ACER functions (default=0.95)
# tdist.CI: is TRUE if CI should be calculatied using t-statistics else poisson
####################################################################################
CI<- function(X,...) UseMethod("CI")
CI.ACER<-function(X,level=0.95,tdist.CI=TRUE,...){
  #Controling inputs
  if(level>=100||level<=0){
    print('level must be between 0 and 100')
    return(0)}
  else if(level>1&&level<100){level=level/100}
  if(!is.logical(tdist.CI)){print('tdist.CI must be logical')} 
  m<-X$dim[1]
  if(m!=1){
    if(tdist.CI){ #t-dist
      conf_coef<-qt((level+1)/2,m-1)
      CI<-conf_coef*X$acer.sd/sqrt(m)
      X$tdist.CI<-TRUE
    }else if(!tdist.CI){ #poisson
      conf_coef<-qnorm((level+1)/2,m-1)
      CI<-conf_coef*sqrt(X$acer)/sqrt(sum(X$RL)-X$k+1)
      X$tdist.CI<-FALSE
    }
  }else{
    if(tdist.CI){
      print('tdist.CI=TRUE, but Only one realization, ')
      print('automatically aplly poisson for CI => tdist.CI=FALSE')}
    conf_coef<-qnorm((level+1)/2,m-1)
    CI<-conf_coef*sqrt(X$acer)/sqrt(sum(X$RL)-X$k+1)
    X$tdist.CI<-FALSE
  }
  X$CI<-CI
  X$CIlevel<-level
  X$call<-match.call()
  return(X)
}