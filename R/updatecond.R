####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is an ACER class
# - eta1 should be included if known
#
# Return a updated ACER class.
# Update the list by removing elemtes where CI>acer. 
# This because the resulting confidence interval would be
# smaller then zero ()
####################################################################################
updatecond<-function(X,...) UseMethod("updatecond")
updatecond.ACER<-function(X, eta1=NA,...){
  condition<-(X$CI<X$acer)
  StartEnd<-findcond(condition)
  if(!is.na(eta1)){
    StartEnd[1]<-which.min(abs(X$eta-eta1))[1]
    X$eta1<-eta1
  }
  X$StartEnd<-StartEnd
  X$call<-match.call()
  return(X)
}

