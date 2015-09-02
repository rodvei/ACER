####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X, k, stationary same as X in ?ACERm
# - CI is the same as CI.level in ?CI.
# - eta1 is the same as eta1 in ?updatecond
#
# Full process of filling out the ACER class
####################################################################################
ACER<- function(X, ...) UseMethod("ACER")
ACER.default<-function(X,k=1,CI=0.95, eta1=NULL, stationary=TRUE, method="general", 
                       check.weibull=TRUE, penalty=FALSE, alpha=0.05, neta=500,...){
  Y<-ACERm(X,k=k,stationary=stationary,neta=500)
  if(!(class(Y)=='ACER')){
    print('Error, << X >> input most be corrct. See ?ACERm')
    return(NULL)
  }
  if(!(method=='general'||method=='gumbel')){
    print('method must be either << general >> or << gumbel >>')
    return(NULL)
  }
  Xtemp=as.matrix(X) 
  if(min(dim(Xtemp))==1){tdist.CI<-FALSE
  }else{tdist.CI<-TRUE}
  
  if(is.null(eta1)){
    plot(Y,CI.level=CI,tdist.CI=tdist.CI)
  }else{
    Y<-CI(Y,level=CI,tdist.CI=tdist.CI)
    Y<-updatecond(Y,eta1=eta1)
  }
  Y<-MSEoptimization(Y,method=method, check.weibull=check.weibull, 
                     penalty=penalty,alpha=alpha)
  Y$call<-match.call()
  return(Y)
}