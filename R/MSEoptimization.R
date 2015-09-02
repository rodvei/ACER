####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is an ACER class 
# - method which the ACER should be estimated to (gumbel or general)
# - penalty factor for more robust estimation
# - alpha is the penalty factor (only used when penalty is TRUE)
#
# Return a updated ACER class with confidence intervall which also 
# include coeficients for general or gumbel model
####################################################################################
MSEoptimization<-function(X,...) UseMethod("MSEoptimization")
MSEoptimization.ACER<-function(X,method="general", check.weibull=TRUE, 
                               penalty=FALSE,alpha=0.05,...){
  #options(warn=-1) #Turn of warnings
  if(method!="general"&&method!="gumbel"){
    print('Error: method should be either "general" or "gumbel".')
    return(X)
  }
  start<-X$StartEnd[1]
  end<-X$StartEnd[2]
  if(is.na(sum(X$CI[start:end]))){
    print('Run CI and/or updatecond first, required for weights')
    return(X)
  }
  eta<-X$eta[start:end]
  logeps<-log(X$acer[start:end])
  lowerCI<-X$acer[start:end]-X$CI[start:end]
  upperCI<-X$acer[start:end]+X$CI[start:end]  
  W<-(log(upperCI)-log(lowerCI))^(-2)
  w<-W/sum(W)
  if(method=="gumbel"){
    guess<-Guess(eta=eta,logeps=logeps,w=w,eta1=X$eta1,mineta=min(X$eta),maxeta=max(X$eta),method=method)
    equation<-'eps(eta)=q*exp[-a(eta-b)^c]'
    sol<-GumbelOptim(guess,eta,logeps,w=w,eta1=X$eta1,mineta=min(X$eta),
                     penalty=TRUE,alpha=0.05)
    upperCI2<-Gumbel(eta,sol[1:4])+X$CI[start:end]
    lowerCI2<-Gumbel(eta,sol[1:4])-X$CI[start:end]
    if(min(lowerCI2)<=0){
      condition<-lowerCI2>0
      StartEndL<-findcond(condition)
      startL<-StartEndL[1];endL=StartEndL[2]
    }else{
      startL<-1;endL=(end-start+1)}
    guessUpper<-Guess(eta=eta,logeps=log(upperCI2),w=w,eta1=X$eta1,mineta=min(X$eta),maxeta=max(X$eta),method=method)
    guessLower<-Guess(eta=eta[startL:endL],logeps=log(lowerCI2[startL:endL]),
                      w=w[startL:endL], eta1=X$eta1,mineta=min(X$eta),maxeta=max(X$eta),method=method)
    solUpper<-GumbelOptim(guessUpper,eta,log(upperCI2),w=w,eta1=X$eta1,
                          mineta=min(X$eta),penalty=TRUE,alpha=0.05)
    solLower<-GumbelOptim(guessLower,eta[startL:endL],log(lowerCI2[startL:endL]),
                          w=w[startL:endL], eta1=X$eta1,mineta=min(X$eta),
                          penalty=TRUE,alpha=0.05)
  }else if(method=="general"){
    if(penalty==TRUE){
      print('Penalty is not necessary in the general case. penalty = FALSE')}
    equation<-'eps(eta)=q*[1+a*xi*(eta-b)^c]^(-1/xi)'
    guess<-Guess(eta=eta,logeps=logeps,w=w,eta1=X$eta1,mineta=min(X$eta),maxeta=max(X$eta),method=method, check.weibull=check.weibull)
    sol<-GeneralOptim(guess,eta,logeps,w,eta1=X$eta1,mineta=min(X$eta),maxeta=max(X$eta))
    upperCI2<-General(eta,sol[1:5])+X$CI[start:end]
    lowerCI2<-General(eta,sol[1:5])-X$CI[start:end]
    if(min(lowerCI2)<0){
      condition<-lowerCI2>0
      StartEndL<-findcond(condition)
      startL<-StartEndL[1];endL=StartEndL[2]
    }else{
      startL<-1;endL=(end-start+1)} 
    guessUpper<-Guess(eta=eta,logeps=log(upperCI2),w=w,eta1=X$eta1,mineta=min(X$eta),maxeta=max(X$eta),method=method)
    guessLower<-Guess(eta=eta[startL:endL],logeps=log(lowerCI2[startL:endL]),
                      w=w[startL:endL], eta1=X$eta1,mineta=min(X$eta),maxeta=max(eta),method=method)
    if(sol[5]<0){
      guessUpper<-guessUpper[guessUpper[,5]<0,]
      guessLower<-guessLower[guessLower[,5]<0,]
    }else{
      guessUpper<-guessUpper[guessUpper[,5]>0,]
      guessLower<-guessLower[guessLower[,5]>0,]
    }
    solUpper<-GeneralOptim(guessUpper,eta,log(upperCI2),w,eta1=X$eta1,
                           mineta=min(X$eta),maxeta=max(X$eta))
    solLower<-GeneralOptim(guessLower,eta[startL:endL],log(lowerCI2[startL:endL]),
                           w[startL:endL], eta1=X$eta1,mineta=min(X$eta),maxeta=max(eta))
  }  
  #options(warn=0)
  X$coef<-sol
  X$upperCIcoef<-solUpper
  X$lowerCIcoef<-solLower
  X$equation<-equation
  X$call<-match.call()
  return(X)
}