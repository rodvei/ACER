####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - Akj number of data which exceed the barrier level (see ACERm)
# - eta barrier level to the corresponding Akj values
#
# Returns the eta value which is 1.5 standard divination away from
# the maximum amount of exceedance (the barrier level where we 
# have the highest number of independent exceedance)
####################################################################################
eta1<-function(Akj,eta){
  row<-dim(Akj)[1]
  col<-dim(Akj)[2]
  pos<-rep(NA,row)
  for(j in 1:row){
    positions<-which(Akj[j,]==max(Akj[j,]))
    pos[j]<-positions[length(positions)]
  }
  if(row==1){return(eta[ceiling(pos)])
  }else{return(eta[ceiling(mean(pos)+1.5*sd(pos))])}
}