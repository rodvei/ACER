####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is an ACER calss (after MSE optimization)
# - eta is the barrier level for where extrapolation should be done
# - prob is the epsilon in the general og Gumbel model
#
# Returns a matrix with known values, Confidence intervalls and extrapolated values
####################################################################################
predict.ACER<-function(object,eta=NULL,prob=NULL,weibull.limit=FALSE,...){
  if(is.na(sum(object$coef))){
    print('need to run MSEoptimization first')
    return(NULL)}
  general<-(length(object$coef)==6)
  if(!is.null(eta)&&!is.null(prob)||
       (is.null(eta)&&is.null(prob)&&(!weibull.limit||object$coef[5]>0))){
    cat('Error: either eta or prob should be NULL. ')
    cat('The variable which is NULL is predictied.')
    if(object$coef[5]<0){
      cat('If limits of the Weibull is desired eta and prob should be NULL') 
      cat('while weibull.limit=TRUE.')}
    return(NULL)
  }else if(!is.null(eta)){
    if(general){
      pred<-General(eta,object$coef)
      upperpred<-General(eta,object$upperCIcoef)
      lowerpred<-General(eta,object$lowerCIcoef)
      if(object$coef[5]<0&&eta>(object$coef[2]-1/(object$coef[1]*object$coef[5]))){
         cat('Weibull does not exist for eta>b-1/(xi*a) and such prediction should be avoided.')
         cat('Run predict(object,weibull.limit=FALSE) to investigate the Weibull limits.')
      }
    }else{
      pred<-Gumbel(eta,object$coef)
      lowerpred<-Gumbel(eta,object$lowerCIcoef)
      upperpred<-Gumbel(eta,object$upperCIcoef)
    }
    ret<-rbind(eta=eta,upperpred,PredictedProb=pred,lowerpred)
    ret[is.na(ret)]<-0
    rownames(ret)[c(2,4)]<-paste(c('Upper','Lower'),c('%CI'),sep=toString(object$CIlevel*100))
    colnames(ret)<-NULL
    return(ret)
  }else if(!is.null(prob)){
    if(general){
      pred<-reverseGeneral(prob,object$coef)
      lowerpred<-reverseGeneral(prob,object$lowerCIcoef)
      upperpred<-reverseGeneral(prob,object$upperCIcoef)
    }else{
      pred<-reverseGumbel(prob,object$coef)
      lowerpred<-reverseGumbel(prob,object$lowerCIcoef)
      upperpred<-reverseGumbel(prob,object$upperCIcoef)
    }
    ret<-rbind(prob=prob,upperpred,PredictedEta=pred,lowerpred)
    rownames(ret)[c(2,4)]<-paste(c('Upper','Lower'),c('CI%'),sep=toString(object$CIlevel*100))
    colnames(ret)<-NULL
    return(ret)
  }else if(weibull.limit&&object$coef[5]<0){
    pred<-object$coef[2]-1/(object$coef[1]*object$coef[5])
    upperpred<-object$upperCIcoef[2]-1/(object$upperCIcoef[1]*object$upperCIcoef[5])
    lowerpred<-object$lowerCIcoef[2]-1/(object$lowerCIcoef[1]*object$lowerCIcoef[5])
    cat('Weibull limits for prediction, upper and lower confidence estimates of eta.')
    ret<-rbind(upperpred,PredictedEtaLim=pred,lowerpred)
    rownames(ret)[c(1,3)]<-paste(c('Upper','Lower'),c('CI%Lim'),sep=toString(object$CIlevel*100))
    colnames(ret)<-NULL
    return(ret)
  }
}
#https://stat.ethz.ch/R-manual/R-patched/library/stats/html/predict.glm.html

Gumbel<-function(eta,est){return(est[4]*exp(-est[1]*(eta-est[2])^est[3]))}
reverseGumbel<-function(eps,est){return(est[2]+(-1/est[1]*log(eps/est[4]))^(1/est[3]))}

General<-function(eta,est){
  if(est[5]<0){eta[eta>(est[2]-1/(est[1]*est[5]))]=NA}
  return(est[4]*(1+est[1]*est[5]*(eta-est[2])^est[3])^(-1/est[5]))
}

reverseGeneral<-function(eps,est){
  return(est[2]+(1/(est[1]*est[5])*((eps/est[4])^(-est[5])-1))^(1/est[3]))
}