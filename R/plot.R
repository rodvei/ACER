####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - x is an ACER class 
# - extrapolate is logic, whether the extrapolation is desired.
#   by defualt TRUE if coeficients have been optimized.
# - prob the disired level you want ACER to be extrapolated 
#   to defualt=1e-6
# - add.eta1 TRUE if model should add eta1 when its not presented 
# - alpha is how transparent the plot should be (range from 0 to 1 
#   where 0 is fully transparent)
#
# Does not return value, only shows a graphical plot
####################################################################################

plot.ACER<-function(x, conf.int=TRUE, CI.level, extrapolate, 
                    prob=1e-6, add.eta1=TRUE, xlab=expression(eta),
                    ylab=expression(hat(epsilon)[k](eta)), type='l', 
                    xlim, ylim, alpha=1, tdist.CI=TRUE, ...){
  if(missing(extrapolate)){
    extrapolate<-(!is.na(sum(x$coef)))
  }
  if(missing(CI.level)){
    if(!is.na(x$CIlevel)){
      CI.level=x$CIlevel
    }else{CI.level=0.95}
  }
    
  if(!extrapolate){
    xtemp<-x
    if(is.na(sum(xtemp$CI))||xtemp$CIlevel!=CI.level){ 
      #make CI.level=0.95 to just CI.level??
      xtemp<-CI(xtemp,level=CI.level,tdist.CI=tdist.CI)
      xtemp<-updatecond(xtemp)
    }else if(is.na(sum(xtemp$StartEnd))){xtemp<-updatecond(xtemp)}
    
    if(is.na(xtemp$eta1)&&missing(xlim)){
      start<-which(xtemp$acer==max(xtemp$acer))[1]
    }else if(!is.na(xtemp$eta1)&&missing(xlim)){
      start<-which.min(abs(xtemp$eta-xtemp$eta1))[1]
    }else if(is.na(xtemp$eta1)&&!missing(xlim)){
      start<-which(abs(xtemp$eta-xlim[1])==min(abs(xtemp$eta-xlim[1])))
    }else if(xlim[1]>xtemp$eta1){
      start=which(abs(xtemp$eta-xlim[1])==min(abs(xtemp$eta-xlim[1])))
    }else{
      start<-which.min(abs(xtemp$eta-xtemp$eta1))[1]            
    }
    end<-xtemp$StartEnd[2]
    absend<-length(xtemp$eta)
    eta<-xtemp$eta[start:end]
    acer<-xtemp$acer[start:end]
    par(mgp=c(2.5,1,0)) #Labals right position in plot
    if(missing(xlim)){
      xlim<-NULL
    }
    if(!conf.int){
      col<-adjustcolor(col=1, alpha.f=alpha)
      plot(xtemp$eta[start:absend],xtemp$acer[start:absend],xlab=xlab,ylab=ylab,
           type=type,xlim=xlim,log='y',col=col,...)
    }else{
      CI<-xtemp$CI[start:end]
      lowerCI<-acer-CI
      upperCI<-acer+CI
      if(missing(ylim)){
        ylim<-c(min(lowerCI),max(upperCI))
      }
      col<-adjustcolor(col=1, alpha.f=alpha)
      CIcol<-adjustcolor(col=4, alpha.f=alpha)
      plot(xtemp$eta[start:absend],xtemp$acer[start:absend],xlab=xlab,ylab=ylab,
           type=type,log='y',ylim=ylim,xlim=xlim,
           col=col,...)
      lines(eta,upperCI,col=CIcol,lty=2)
      lines(eta,lowerCI,col=CIcol,lty=2)
    }
    if(is.na(xtemp$eta1)){
      abline(v=x$eta1est,lty=2,col='indianred1')
      cat('Choose eta1, where plot indicate regular teilbehavior \n')
      cat('Program suggest regular tailbehavior for eta =',x$eta1est)
      if(add.eta1){
        eta1temp <- readline(prompt="Enter value of eta1: ")
        xtemp<-updatecond(xtemp,as.double(eta1temp))
        eval.parent(substitute(x<-xtemp))
      }
    }
    par(mgp=c(3,1,0)) # Labals back to defalt position
  }else if(extrapolate){
    if(is.na(x$eta1)){
      plot(x,conf.int=conf.int,CI.level=CI.level,tdist.CI=tdist.CI,type=type)   
    }
    start<-x$StartEnd[1]
    end<-x$StartEnd[2]
    eta<-x$eta[start:end]
    acer<-x$acer[start:end] 
    if(is.na(sum(x$coef))){x<-MSEoptimization(x)}
    if(!(length(x$coef)==5||length(x$coef)==6)){
      cat('Error:x$coef should have length = 5 or 6 dependent ')
      cat('if it is gumbel or general')
      return(0)
    }
    general<-(length(x$coef)==6)
    CI<-x$CI[start:end]
    upperCI<-acer+CI
    if(!is.null(prob)&&missing(ylim)){
      ylim<-c(prob*0.1,max(upperCI))
    }else if(is.null(prob)&&missing(ylim)){
      print('if prob=NULL, should spesify ylim')}
    if(missing(xlim)){
      xlim<-c(x$eta1,NA)
      if(general){
        xlim[2]<-reverseGeneral(prob*0.1,x$upperCIcoef)
#        cat("With General method and ACER (eps) =",prob,", the corresponding \n")
#        cat("estimated eta =",reverseGeneral(prob,x$coef),"\n")
      }else if(!general){
        xlim[2]<-reverseGumbel(prob*0.1,x$upperCIcoef)
#        cat("With Gumbel method and ACER (eps) =",prob,", the corresponding \n")
#        cat("estimated eta =",reverseGumbel(prob,x$coef)."\n")
      }
    }
    plot(x,extrapolate=FALSE,conf.int=conf.int,CI.level=CI.level,
         tdist.CI=tdist.CI, xlab=expression(eta),
         ylab=expression(hat(epsilon)[k](eta)),type='l',ylim=ylim, 
         xlim=xlim,alpha=0.4,...)
    xeta<-seq(max(xlim[1],x$eta1),xlim[2],length.out=100)
    if(x$coef[5]<0){
      etaLower<-reverseGeneral(prob*0.05,x$lowerCIcoef)
      etaPred<-reverseGeneral(prob*0.05,x$coef)
      etaUpper<-reverseGeneral(prob*0.05,x$upperCIcoef)
      xeta[which(abs(xeta-etaLower)==min(abs(xeta-etaLower)))]<-etaLower
      xeta[which(abs(xeta-etaPred)==min(abs(xeta-etaPred)))]<-etaPred
      xeta[which(abs(xeta-etaUpper)==min(abs(xeta-etaUpper)))]<-etaUpper
    }
    if(general){
      lines(xeta,General(xeta,x$coef),type='l')
      lines(xeta,General(xeta,x$upperCIcoef),lty=2,col=4)
      lines(xeta,General(xeta,x$lowerCIcoef),lty=2,col=4)
      if(x$coef[5]<0){
        abline(v=x$coef[2]-1/(x$coef[1]*x$coef[5]),lty=2)
      }
    }else if(!general){
      lines(xeta,Gumbel(xeta,x$coef),type='l')
      lines(xeta,Gumbel(xeta,x$upperCIcoef),lty=2,col=4)
      lines(xeta,Gumbel(xeta,x$lowerCIcoef),lty=2,col=4)
    }
    abline(h=prob,lty=2,col='indianred1')
  }
}