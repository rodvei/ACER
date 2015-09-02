####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - eta is the eta for each logeps (within the suitable area)
# - logeps or log(ACER).
# - w is the weight factor
# - eta1, is the beginning of regular tail behavior
# - mineta is the minimum eta used
#
# Returns a matrix which contains 4 different start values (for gumbel).
# The 4 guesses are C=1.1, 2, 3, 4.
####################################################################################
Guess<-function(eta,logeps, w, eta1, mineta, maxeta, method, check.weibull=TRUE){  
  # eps=q(1+a(n-b)^c)^-y
  # log(eps)=log(q)-y log(1+a(n-b)^c)
  if(method=="gumbel"){
    guess<-matrix(NA,4,4)
  }else if(method=="general"&&check.weibull==FALSE){guess<-matrix(NA,4,5)
  }else if(method=="general"&&check.weibull==TRUE){guess<-matrix(NA,6,5)}
  y<-logeps
  x<-eta
  x2<-x^2
  x3<-x^3
  x4<-x^4
  
  c<-1.1 #c=1.1 because if c=1 there are inf # of q that is corect
  mod1<-lm(y~x)
  a<--mod1$coef[2]
  b<-mod1$coef[1]/a
  coef1<-c(a,b,1.1,1)
  if(b>eta1){
    coef1<-GumbelParam(c(eta1,c),eta=x,logeps=y,w=w)
  }else if(b<mineta){coef1<-GumbelParam(c(mineta,c),eta=x,logeps=y,w=w)}
  guess[1,1:4]<-coef1
  rm(a,b,c)
  
  c<-2
  mod2<-lm(y~x+x2)
  a<--mod2$coef[3]
  b<-mod2$coef[2]/(2*a)
  q<-exp(mod2$coef[1]+a*b^2)
  coef2<-c(a,b,c,q)
  if(b>eta1){
    coef2<-GumbelParam(c(eta1,c),eta=x,logeps=y,w=w)
  }else if(b<mineta){coef2<-GumbelParam(c(mineta,c),eta=x,logeps=y,w=w)}
  guess[2,1:4]<-coef2
  rm(a,b,c,q)
  
  c<-3
  mod3<-lm(y~x+x2+x3)
  a<--mod3$coef[4]
  b<-mod3$coef[3]/(3*a)
  q<-exp(mod3$coef[1]-a*b^3)
  coef3<-c(a,b,c,q)
  if(b>eta1){
    coef3<-GumbelParam(c(eta1,c),eta=x,logeps=y,w=w)
  }else if(b<mineta){coef3<-GumbelParam(c(mineta,c),eta=x,logeps=y,w=w)}
  guess[3,1:4]<-coef3
  rm(a,b,c,q)
  
  c<-4
  mod4<-lm(y~x+x2+x3+x4)
  a<--mod4$coef[5]
  b<-mod4$coef[4]/(4*a)
  q<-exp(mod4$coef[1]+a*b^4)
  coef4<-c(a,b,c,q)
  if(b>eta1){
    coef4<-GumbelParam(c(eta1,c),eta=x,logeps=y,w=w)
  }else if(b<mineta){coef4<-GumbelParam(c(mineta,c),eta=x,logeps=y,w=w)}
  guess[4,1:4]<-coef4
  rm(a,b,c,q)
  
  if(method=="general"){
    # Convert guess from gumbel to general
    whichguess<-guess[1:4,1]>0
    if(check.weibull){whichguess<-c(whichguess,!is.na(guess[5:length(guess[,1])]))}
    guess[whichguess,]<-t(apply(guess[whichguess,1:3],1,GeneralParamf,eta=eta,logeps=logeps,w=w))
    guess[!whichguess,]<-NA
    # Weibull guess
    if(check.weibull==TRUE){
      c<-1
      b<-eta1+2*(maxeta-eta1)/3
      guess[5,]<-GeneralParamf(est=c(-1/(maxeta-b)*0.8,b,c),eta,logeps,w)
      rm(b)
      
      b<-eta1+(maxeta-eta1)/3
      guess[6,]<-GeneralParamf(est=c(-1/(maxeta-b)*0.9,b,c),eta,logeps,w)
    }
    colnames(guess)<-c('a','b','c','q','xi')
    rownames(guess)<-NULL
  }else if(method=="gumbel"){    
    colnames(guess)<-c('a','b','c','q')
    rownames(guess)<-NULL
  }
  return(guess)  
}