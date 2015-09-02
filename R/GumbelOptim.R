####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - guess is the startvalue of the iteration 
# - eta is the eta for each logeps (within the suitable area)
# - logeps or log(ACER).
# - w is the weight factor
# - eta1, is the beginning of regular tail behavior
# - mineta is the minimum eta used
# - penalty indicate if penalty should be used
# - alpha is the penalty factor
#
# Returns the optimized coeficients for the genral case together
# with the error term.
####################################################################################
GumbelOptim<-function(guess,eta,logeps,w,eta1,mineta,penalty=FALSE,alpha=0.05){
  options(warn=-1)
  guessLength<-dim(guess)[1]
  sol<-matrix(NA,guessLength,5)
  for(i in 1:guessLength){
    startpara<-guess[i,]
    tempsol<-matrix(NA,5,5)
    
    # minimizing mse using nlm with 2 parameters
    gumbel2nlm<-tryCatch(nlm(mseGumbel2,startpara[c(2,3)],eta=eta,logeps=logeps,w=w,
                        penalty=penalty,alpha=alpha)$estimate,
                        error=function(e) rep(NA,2))
    tempsol[1,]<-c(GumbelParam(gumbel2nlm,eta,logeps,w),
                  mseGumbel4(GumbelParam(gumbel2nlm,eta,logeps,w),eta,logeps,w))
    
    #minimizing mse using nlminb 4 param with box
    if(penalty==FALSE){
      gumbel4box<-tryCatch(nlminb(start=startpara,mseGumbel4,gradient=GradmseGumbel4,
                          eta=eta,logeps=logeps,w=w,lower=c(0,mineta,0,
                          .Machine$double.eps),upper=c(Inf,eta1,5,Inf),control=
                          list(iter.max=300,eval.max=400))$par,
                          error=function(e) rep(NA,4))
      tempsol[2,]<-c(gumbel4box,mseGumbel4(gumbel4box,eta,logeps,w))
    }else{
      gumbel4box<-tryCatch(nlminb(start=startpara,mseGumbel4,eta=eta,logeps=logeps,
                          penalty=penalty,alpha=alpha, w=w,lower=c(0,mineta,0,
                          .Machine$double.eps),upper=c(Inf,eta1,5,Inf),control=
                          list(iter.max=300,eval.max=400))$par,
                          error=function(e) rep(NA,4))
      tempsol[2,]<-c(gumbel4box,mseGumbel4(gumbel4box,eta,logeps,w))
    }
    
    #minimizing mse using LM 4 param with box
    if(penalty==FALSE){
      gumbel4LM<-tryCatch(coef(nlsLM(logeps~log(qq)-aa*(eta-bb)^cc, start = 
                         list(aa = startpara[1],bb =startpara[2], cc=startpara[3]
                         ,qq=startpara[4]),weights=w,lower=c(0,mineta,0,
                         .Machine$double.eps),upper=c(Inf,eta1,5,Inf),control=
                         list(maxiter=500))),error=function(e) rep(NA,4))
      tempsol[3,]<-c(gumbel4LM,mseGumbel4(gumbel4LM,eta,logeps,w))}
    
    #minimizing mse using LM 2 param with box
    if(penalty==FALSE){
      gumbel2LM<-tryCatch(coef(nlsLM(logeps~mean(logeps)+sum(w*((eta-bb)^cc-
                         mean((eta-bb)^cc))*(logeps-mean(logeps)))/sum(w*((eta-bb)^cc
                         -mean((eta-bb)^cc))^2)*((eta-bb)^cc-mean((eta-bb)^cc)),
                         start = list(bb = startpara[2], cc=startpara[3]),weights=w,
                         control=list(maxiter=500),lower=c(mineta,0),
                         upper=c(eta1,5))),error=function(e) rep(NA,2))
      tempsol[4,]<-c(GumbelParam(gumbel2LM,eta,logeps,w),
                    mseGumbel4(GumbelParam(gumbel2LM,eta,logeps,w),eta,logeps,w))}
    
    #minimizing mse using "L-BFGS-B" 4 parameters
    if(penalty==FALSE){
      gumbel4optim<-tryCatch(optim(par=startpara,fn=mseGumbel4,gr=GradmseGumbel4,
                            eta=eta,logeps=logeps,w=w,lower=c(0,mineta,0,
                            .Machine$double.eps),upper=c(Inf,eta1,5,Inf),control=
                            list(maxit=1000),method='L-BFGS-B')$par,
                            error=function(e) rep(NA,4))
      tempsol[5,]<-c(gumbel4optim,mseGumbel4(gumbel4optim,eta,logeps,w))
    }else{
      gumbel4optim<-tryCatch(optim(par=startpara,fn=mseGumbel4,eta=eta,logeps=logeps,
                            penalty=penalty,alpha=alpha,w=w,lower=c(0,mineta,0,
                            .Machine$double.eps),upper=c(Inf,eta1,5,Inf),control=
                            list(maxit=1000))$par,error=function(e) rep(NA,4))                                  
      tempsol[5,]<-c(gumbel4optim,mseGumbel4(gumbel4optim,eta,logeps,w))
    }
    
    conditionMatrix<-t(apply(tempsol[,1:4],1,function(x) 
      ((x<=c(Inf,eta1,5,Inf))*(x>=c(0,mineta,0,.Machine$double.eps)))))
    
    if(sum(is.na(tempsol[,5]))>=5){sol[i,]<-rep(NA,5)
    }else{
      condition<-floor(apply(conditionMatrix,1,sum)/4)
      condition[condition==0]<-NA
      sol[i,]<-tempsol[which.min(tempsol[,5]*condition),]
    }
  }
  options(warn=0)
  gumbelsol<-sol[which.min(sol[,5]),]
  names(gumbelsol)<-c('a','b','c','q','W.MSE')
  return(gumbelsol)
}