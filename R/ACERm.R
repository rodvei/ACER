####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is the data set or peaks 
# - k is an integer, which represent the k-dependence 
#   (see article description)
# - stationary is a logical value. Should be TRUE if timeseries 
#   is stationarry
# - neta is the number of barrier levels on the interval min to max
#   of the data set
#
# Return a ACER class, which is a list of values of interest.
# Includes the average conditional exceedance rates (at listname 
# acer),
####################################################################################
ACERm<- function(X, ...) UseMethod("ACERm")
ACERm.default<-function(X,k=1,stationary=TRUE,eta1=TRUE,neta=500,...){
  #Controling inputs
  if(!is.logical(stationary)){
    print('stationary must be logical')
    return(NULL)
  }
  # if vector, makes vector into matrix
  if(is.vector(X)){
    X<-X[!is.na(X)]
    n<-length(X)
    m<-1
    RL<-n
    X<-t(as.matrix(X))
  }else if (!is.matrix(X)){
    print('The inserted data set should be of class vector or matrix')
    return(NULL)
  }else{
    m<-dim(X)[1]
    RL<-RealizationLength(X)
    Xtemp=cbind(matrix(NA,m,k),X)
    Xtemp[1,1:k]=X[m,(RL[m]-k+1):RL[m]]
    for(i in 2:m){
      Xtemp[i,1:k]=X[i-1,(RL[i-1]-k+1):RL[i-1]]
    }
    X=Xtemp
    n<-dim(X)[2]
    RL<-RealizationLength(X)
  }
  if(!(RealizationByRows(X))){return(X)} #check that realization is by rows
  
  etaStart=min(X,na.rm=TRUE)
  etaEnd=max(X,na.rm=TRUE)
  eps<-.Machine$double.eps # smallest floating point.
  Deta<-(etaEnd-etaStart+2*eps)/(neta-1)
  barrierLevels<-seq(etaStart-eps,etaEnd+eps,Deta)
  Akj<-matrix(NA,m,neta)
  Bkj<-Akj
  eps_hat<-Akj
  eps_hat_<-Akj
  NDeta<-(X-barrierLevels[1])/Deta+1 # number of barier levels each peak contain
  NDeta[which(NDeta==max(NDeta,na.rm=TRUE),arr.ind=TRUE)]<-neta 
  # remove roundof problems
  NDeta[which(NDeta==min(NDeta,na.rm=TRUE),arr.ind=TRUE)]<-1        
  # remove roundof problems
  NDetaFloor<-floor(NDeta)           # number of barier levels each peak crosses
  NDetaCeiling<-ceiling(NDeta)       # number of barier levels not met by each peak 
  
  for(j in 1:m){
    num<-rep(0,length(barrierLevels))
    if(k==1){
      for(i in 1:RL[j]){
        counter<-1:NDetaFloor[j,i]
        num[counter]<-num[counter]+1
      }
      den<-rep(RL[j],length(barrierLevels))
    }
    else{ #k>1
      den<-num
      for(i in k:RL[j]){
        counterDen<-max(NDetaCeiling[j,(i-k+1):(i-1)])
        den[counterDen:neta]<-den[counterDen:neta]+1
        if(counterDen<=NDetaFloor[j,i]){
          counterDen<-counterDen:NDetaFloor[j,i]
          num[counterDen]<-num[counterDen]+1
        }
      }
      den[den==0]<-eps      
    }
    Akj[j,]<-num
    Bkj[j,]<-den
    eps_hat[j,]<-num/den
    eps_hat_[j,]<-num/(RL[j]-k+1)
  }
  rm(i,j) #remove i and j
  
  #Calculating Mean & SD of ACER functions
  if(m!=1){
    if(stationary){
      eps_hat_mean<-apply(eps_hat,2,mean)
      eps_hat_std<-apply(eps_hat,2,sd) 
    }else if(!stationary){
      eps_hat_mean<-apply(eps_hat_,2,mean)
      eps_hat_std<-apply(eps_hat_,2,sd)
    }
  }else{#only 1 realization => poisson
    if(stationary){
      eps_hat_mean<-as.vector(eps_hat)
      eps_hat_std<-NA
    }else if(!stationary){
      eps_hat_mean<-as.vector(eps_hat_)
      eps_hat_std<-NA    
    }
  }
  est<-list(eta=barrierLevels,acer=eps_hat_mean)
  est$call<-match.call()
  est$CI<-NA
  est$CIlevel<-NA
  est$tdist.CI<-NA
  est$acer.sd<-eps_hat_std
  est$coef<-NA
  est$upperCIcoef<-NA
  est$lowerCIcoef<-NA
  if(eta1){est$eta1est<-eta1(Akj,barrierLevels)
  }else{est$eta1est<-NA}
  est$k<-k
  est$dim<-c(m,n)
  est$RL<-RL
  est$eta1<-NA
  est$StartEnd<-NA
  est$equation<-NA
  
  class(est)<-"ACER"
  return(est)  
}