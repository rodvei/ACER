####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is the data set
# 
# Return length of series (until Nan)
# Calucalting the Realization Length (Length of the peaks 
# extracted from time series for each realization)
# Y is all data for j'th realisations.
####################################################################################

RealizationLength<-function(X){
  nrow<-length(X[,1]) #numbers of realization
  ncol<-length(X[1,]) #length of each realization
  RL<-rep(NA,nrow) #RL will hold the length of each realization
  for(j in 1:nrow){
    Y<-X[j,]
    if(!is.na(Y[ncol])){RL[j]<-ncol
    }else if(is.na(Y[1])){RL[j]<-0
    }else{
      i<-floor(ncol/2)
      ind<-TRUE;
      min<-1
      max<-ncol
      while(ind){
        if(is.na(Y[i])){
          max<-i
          i<-floor((max-min)/2)+min
        }else{
          if(is.na(Y[i+1])&&!is.na(Y[i])){
            ind<-FALSE
            RL[j]<-i}
          min<-i
          i<-ceiling((max-min)/2)+min
        }
      }
    }
  }
  return(RL)
}