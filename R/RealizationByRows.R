####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is the data set
# 
# Check if realizations is by rows, while measurements is by columns
# (which is needed for further functions).
# Return TRUE if so, else the function print an problem massage 
# and returns FALSE 
####################################################################################
RealizationByRows<-function(X){
  if (dim(X)[1]>dim(X)[2]){print("Realizations should be by rows. Suggest transformation <X=t(X)>")}
  return(!(dim(X)[1]>dim(X)[2]))
}