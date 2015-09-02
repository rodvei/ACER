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
#
# Returns the optimized coeficients for the genral case together
# with the error term.
####################################################################################

GeneralOptim<-function(guess, eta, logeps, w, eta1, mineta, maxeta){
  options(warn=-1)
  guessLength<-dim(guess)[1]
  guessindex<-(!is.na(guess)%*%rep(1,5))*(1:guessLength)
  guessindex<-guessindex[guessindex!=0]
  sol<-matrix(NA,guessLength,6)
  
  nlsLM3f<-logeps~mean(logeps)+sum(w*(log(1+aa*(eta-bb)^cc)-
              mean(log(1+aa*(eta-bb)^cc)))*(logeps-mean(logeps)))/
              sum(w*(log(1+aa*(eta-bb)^cc)-mean(log(1+aa*(eta-bb)^cc)))^2)*
              (log(1+aa*(eta-bb)^cc)-mean(log(1+aa*(eta-bb)^cc)))
  nlsLM5f<-logeps~log(qq)-log(1+aa*xxi*(eta-bb)^cc)/xxi
  if(guess[guessindex[length(guessindex)],5]<0){
    nlsLM3w<-logeps~mean(logeps)+sum(w*(log(1+aa*(eta-bb))-
                mean(log(1+aa*(eta-bb))))*(logeps-mean(logeps)))/
                sum(w*(log(1+aa*(eta-bb))-mean(log(1+aa*(eta-bb))))^2)*
                (log(1+aa*(eta-bb))-mean(log(1+aa*(eta-bb))))
    nlsLM5w<-logeps~log(qq)-log(1+aa*xxi*(eta-bb))/xxi    
  }
  
  for(i in guessindex){
    startpara<-guess[i,]
    if(startpara[5]>=0){
      lowerFull<-c(.Machine$double.eps,mineta,.Machine$double.eps,.Machine$double.eps,.Machine$double.eps)
      upperFull<-c(Inf,eta1,5,Inf,Inf)
      lowerReduced<-lowerFull[1:3]
      upperReduced<-upperFull[1:3]
      nlsLM3<-nlsLM3f
      nlsLM5<-nlsLM5f
      mseReduced<-mseGeneral3f
      mseFull<-mseGeneral5f
      GeneralParam<-GeneralParamf
      Gradient<-GradmseGeneral5f
      errorReturn<-list(rep(NA,5),rep(NA,3))
      
      startparaReduced<-c(startpara[1]*startpara[5],startpara[c(2,3)])
      startList<-list(aa = startpara[1],bb =startpara[2],
                      cc=startpara[3],qq=startpara[4],xxi=startpara[5])
      startListReduced<-list(aa=(startpara[1]*startpara[5]),bb = startpara[2], cc=startpara[3])
      tempsol<-matrix(NA,5,6)
    }else{
      lowerFull<-c(.Machine$double.eps,eta1,.Machine$double.eps,-Inf)
      upperFull<-c(Inf,maxeta,Inf,-.Machine$double.eps)
      lowerReduced<-c(-Inf,lowerFull[2])
      upperReduced<-c(-.Machine$double.eps,upperFull[2])
      nlsLM3<-nlsLM3w
      nlsLM5<-nlsLM5w
      mseReduced<-mseGeneral3w
      mseFull<-mseGeneral5w
      GeneralParam<-GeneralParamw
      Gradient<-GradmseGeneral5w
      errorReturn<-list(rep(NA,4),rep(NA,2))
      
      startpara<-c(startpara[1:2],startpara[4:5])
      startparaReduced<-c(startpara[1]*startpara[4],startpara[2])
      startList<-list(aa = startpara[1],bb =startpara[2],
                      qq=startpara[3],xxi=startpara[4])
      startListReduced<-list(aa=(startpara[1]*startpara[4]),bb = startpara[2])
      tempsol<-matrix(NA,5,5)
    }
    
    #values scaled: a,b,c,q,xi 
    # minimizing mse using nlm with 3 param a~,b,c where a=a~/xi
    
    general3nlm<-tryCatch(nlm(mseReduced,startparaReduced,eta=eta,logeps=logeps,
                              w=w)$estimate,error=function(e) errorReturn[[2]])
    if(sum(is.na(general3nlm))!=0||(sum(general3nlm>=lowerReduced)!=length(general3nlm)
                    ||sum(general3nlm<=upperReduced)!=length(general3nlm))){
        general3nlm=errorReturn[[2]]
    }
    tempsol[1,]<-c(GeneralParam(general3nlm,eta,logeps,w),
                   mseFull(GeneralParam(general3nlm,eta,logeps,w),eta,logeps,w))
    
    #minimizing mse using LM 3 param with box
    general3LM<-tryCatch(coef(nlsLM(nlsLM3,start=startListReduced,
                         weights=w, control=list(maxiter=500),lower=lowerReduced,
                         upper=upperReduced)),error=function(e) errorReturn[[2]])
    tempsol[2,]<-c(GeneralParam(general3LM,eta,logeps,w),
                  mseFull(GeneralParam(general3LM,eta,logeps,w),eta,logeps,w)) 
    
#    # The rest of the parameters will use coef with 5 paramteres (need to add xi
#    # to startpara). length(startpara)=4->length(startpara)=5
#    startpara<-GeneralParam(startpara[1:3],eta,logeps,w)
    
    # minimizing mse using nlminb with 5 param a,b,c,q,xi, a~=a*xi
    general5box<-tryCatch(nlminb(start=startpara,mseFull,gradient=Gradient,
                          eta=eta,logeps=logeps,w=w,control=list(iter.max=300,
                          eval.max=400),lower=lowerFull,
                          upper=upperFull)$par,error=function(e) errorReturn[[1]])
    tempsol[3,]<-c(general5box,mseFull(general5box,eta,logeps,w))
    
    #minimizing mse using LM 5 param with box
    general5LM<-tryCatch(coef(nlsLM(nlsLM5, 
                         start = startList,weights=w,
                         lower=lowerFull,upper=upperFull,control=list(maxiter=500))),
                         error=function(e) errorReturn[[1]])
    tempsol[4,]<-c(general5LM,mseFull(general5LM,eta,logeps,w))
    
    #minimizing mse using "L-BFGS-B" 5 parameters 
    general5optim<-tryCatch(optim(par=startpara,fn=mseFull,gr=Gradient,
                           eta=eta,logeps=logeps,w=w,lower=lowerFull,upper=upperFull,
                           control=list(maxit=1000),method='L-BFGS-B')$par,
                           error=function(e) errorReturn[[1]])
    tempsol[5,]<-c(general5optim,mseFull(general5optim,eta,logeps,w))
    
    if(length(tempsol[1,])==5){
      tempsol<-cbind(tempsol[,1:2],rep(1,5),tempsol[,3:5])
      upperFull<-c(upperFull[1:2],1,upperFull[3:4])
      lowerFull<-c(lowerFull[1:2],1,lowerFull[3:4])
    }

    conditionMatrix<-t(apply(tempsol[,1:5],1,function(x) 
      ((x<=upperFull)*(x>=lowerFull))))
    
    if(sum(!is.finite(tempsol[,6]))>=5){sol[i,]<-rep(NA,6)
    }else{
      condition<-floor(apply(conditionMatrix,1,sum)/5)
      condition[condition==0]<-NA
      sol[i,]<-tempsol[which.min(tempsol[,6]*condition),]
    }    
  }
  options(warn=0)
  generalsol=sol[which.min(sol[,6]),]
  names(generalsol)<-c('a','b','c','q','xi','W.MSE')
  return(generalsol)
}