####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - x is an ACER object 
#
# Plor, returns a summary or plot a summary of an ACER class.
####################################################################################
print.ACER <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat('\nACER class with k=',x$k,'\n')
  if(x$call[1]=='ACERm.default()'){
    cat('\nnext run plot.ACER() or CI.ACER()\n')  
  }else if(x$call[1]=='CI.ACER()') {
    cat('\nnext run plot.ACER() or updatecond.ACER(X, eta1=<<suitableNumber>>,...)\n')
  }else if(x$call[1]=='updatecond.ACER()'){
    if(is.na(x$eta1)){
      cat('\nnext run plot.ACER() or updatecond.ACER(X, eta1=<<suitableNumber>>,...)')
    }else{cat('\nnext run plot.ACER() or MSEoptimization.ACER()')}
  }else if((x$call[1]=='ACER.default()')||(x$call[1]=='MSEoptimization.ACER()')){
    cat("\nCoefficients:\n")
    coef<-t(as.matrix(x$coef))
    rownames(coef)<-'Estimate'
    printCoefmat(coef)
  }
}

summary.ACER<-function(object,...){
  data<-c(length=sum(object$RL),realization=length(object$RL))
  df=data[2]-1
  TAB<-rbind(object$upperCIcoef,
             Estimate=object$coef,
             object$lowerCIcoef)
  rownames(TAB)[c(1,3)]=paste(c('Upper','Lower'),c('CI%'),sep=toString(object$CIlevel*100))
  range=object$StartEnd[1]:object$StartEnd[2]
  res<-list(call=object$call,eta=object$eta[range],acer=object$acer[range], k=object$k, coefficients=TAB, data=data,t.dist=object$tdist.CI, 
            eta1=object$eta1, numobs=diff(object$StartEnd),equation=object$equation, df=df)
  class(res)<-"summary.ACER"
  return(res)
}

print.summary.ACER<-function(x,...){
  cat("Call:\n")
  print(x$call)
  cat('\nACER class with k=',x$k,'.',sep = "")
  if(x$call[1]=='ACERm.default()'){
    cat('\nnext run plot.ACER() or CI.ACER()\n')  
  }else if(x$call[1]=='CI.ACER()') {
    cat('\nnext run plot.ACER() or updatecond.ACER(X, eta1=<<suitableNumber>>,...)\n')
  }else if(x$call[1]=='updatecond.ACER()'){
    if(is.na(x$eta1)){
      cat('\nnext run plot.ACER() or updatecond.ACER(X, eta1=<<suitableNumber>>,...)\n')
    }else{cat('\nnext run plot.ACER() or MSEoptimization.ACER()\n')}
  }else if((x$call[1]=='ACER.default()')||(x$call[1]=='MSEoptimization.ACER()')){
    if(length(x$coefficients[1,])==6){
      cat('\nEstimated general model: ', x$equation,'\n')
    } else if(length(x$coefficients[1,])==5){
      cat('\nEstimated Gumbel model: ', x$equation,'\n')
    }
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients)
  }
  if(!is.na(x$t.dist)){
    if(x$t.dist==TRUE){
      cat("\nConfidence intervals is calculated using t-distribution with",(x$df),"df.")
    }else if(x$t.dist==FALSE){
      cat("\nConfidence intervals is calculated using poisson-distribution")
    }
  }
  if((x$call[1]=='ACER.default()')||(x$call[1]=='MSEoptimization.ACER()')){
    cat("\nStart of regular tail behavior is set to eta1=", x$eta1,".\nThis gives ",sep = "")
    cat(x$numobs, "suitable ACER's used in the weighted MSE optimization.")
  }
  cat('\nThe ACER model was build using',x$data[1], 'data points with',x$data[2],'realizations.' ) 
}