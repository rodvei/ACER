####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is a array or matrix combined with the number of neighbors 
#   a point need to exceed to be considered a peak.
# - n is the distanse from the point where points is considered to 
#   be a neighbore
# Return an array or matrix which include all peacks from each 
# realizations
####################################################################################
peak<-function(X,n=1){
  # Extract peaks from matrix X where timeseries data is by colums and 
  # realizations buy raws 
  if(is.vector(X)){ #if vector, X get turned in to matrix
    X<-X[!is.na(X)] #remove NA's
    X<-t(as.matrix(X)) # make it a [1 x length(X)] matrix
  }
  if(!(RealizationByRows(X))){return(X)} #check that realization is by rows
  col<-dim(X)[2] #n
  row<-dim(X)[1] #m
  process<-matrix(NA,row,ceiling(col/(n+1)))
  ii<-t(rbind(rep(TRUE,row),diff(t(X))!=0))
  # diff each time series with itself and add TRUE in first colum 
  # and add TRUE in first colum (lost by diff)
  # ii purpos is to remove equal neighbors
  for(i in 1:row){
    ind<-0
    dummy<-X[i,ii[i,]]
    # Remove equal values wich are next to each other
    dummy<-dummy[!is.na(dummy)]
    target<-1+n
    while((target+n)<=length(dummy)){
      dummyMax<-max(dummy[(target-n):(target+n)])
      if(dummy[target]>=dummyMax){ #kept if larger then its n neighbors
        ind<-ind+1
        process[i,ind]<-dummy[target]
        target<-target+n+1
      }
      else{target<-target+1}
    }
  }
  return(process[,1:max(RealizationLength(process))])
}