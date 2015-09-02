####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - condition, an array where each element is logical. TRUE when
#   condition is met, else FALSE.
#
# Return an array of two elements, first start posistion, second
# last positon of which all elements in between is TRUE. 
####################################################################################
findcond<-function(condition){
  if(min(condition)==0){
    # condition=X # for testing purposees
    val=diff(condition) #Finding intermitten data
    poz=which(val!=0) #position of change
    if(length(poz)!=0){
      leng=diff(c(0,poz,length(condition)))
      N=length(leng)
      if(val[poz[1]]==-1){maxLeng=max(leng[seq(1,N,2)]) 
      # first block contains ones
      # seq step size 2 because ones go odd
      }else{maxLeng=max(leng[seq(2,N,2)])}
      # first block contains zeros => ones go even
      n=which(leng==maxLeng) # position or largest set containing ones
      start=sum(leng[1:n-1])+1
      end=sum(leng[1:n-1])+maxLeng  
    }else{print("ERROR: CI is always larger then ACER, try lower CI level")}
  }else if(sum(condition)==length(condition)){
    start=1
    end=length(condition)
  }
  return(c(start,end))
}