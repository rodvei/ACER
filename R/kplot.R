####################################################################################
# Kristoffer Kofoed Roedvei, NTNU, 2014.
#-----------------------------------------------------------------------------------
# Function wich take the data:
# - X is the data set or peaks 
# - k is a array of numbers, which represent the k-dependence 
#   (see article description)
# - stationary is a logical value. Should be TRUE if timeseries 
#   is stationarry
# - neta is the number of barrier levels on the interval min to max
#   of the data set
#
# Will make plots wich can be used to find the desired k for further
# analysis.
# Plots ACER function with different k, for k-selection
# Input dataset and vector containing k values.
####################################################################################

kplot<-function(X,k=1:6,stationary=TRUE,neta=500,xlab=expression(eta),
                ylab=expression(hat(epsilon)[k](eta)),...){
  k<-sort(k)
  dots=list(...)
  # Make plot collors (saved in cols)
  palette <- colorRampPalette(colors=c('#E2E2F2',rgb(t(col2rgb(colors()[30]))/255)))
  data<-ACERm(X,k[1],neta=neta,stationary=stationary,eta1=FALSE)
  par(mgp=c(2.5,1,0)) #change labels position so its inside plot
  start=1
  end=neta
  # how to handel xlim
  if(!is.null(dots$xlim)){
    start=tail(which(data$eta<=dots$xlim[1]),n=1)
    if(dots$xlim[2]<max(X,na.rm=TRUE)){end=which(data$eta>=dots$xlim[2])[1]
    }else{end=neta}
  }
  
  if(k[1]==1){
    cols<-c('#000000',palette(length(k)-1))
    plot(data$eta[start:end],data$acer[start:end],col=cols[1],type='l', 
         xlab=xlab, ylab=ylab,log='y',...)
  }
  else{
    cols<-palette(length(k))
    if(start<which(data$acer==max(data$acer))[1]){
      startTemp<-which(data$acer==max(data$acer))[1]
    }else{startTemp=start}
    plot(data$eta[startTemp:neta],data$acer[startTemp:neta],col=cols[1],
         type='l', xlab=xlab, ylab=ylab,log='y',...)
  }
  
  for(i in 2:length(k)){
    data<-ACERm(X,k[i],neta=neta,stationary=stationary,eta1=FALSE)
    if(start<which(data$acer==max(data$acer))[1]){
      startTemp<-which(data$acer==max(data$acer))[1]
    }else{startTemp=start}
    lines(data$eta[startTemp:neta],data$acer[startTemp:neta],col=cols[i])
  }
  if(length(k)<=7){legend_size=1}
  else{legend_size=length(k)/7}
  legend("bottomleft", inset=0,
         paste(rep('k=',length(k)), paste(k), sep = ""), fill=cols, horiz=FALSE,cex=1/legend_size)
  par(mgp=c(3,1,0)) # Labals back to defalt position
}