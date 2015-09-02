\name{ACER-package}
\alias{ACER-package}
\alias{ACERpackage}
\alias{ACER-packages}
\alias{ACERpackages}
\alias{ACERinfo}
\docType{package}
\title{ Average Conditional Exceedance Rate with Additional Useful Functions}
\description{
ACER is a package made for calculating average conditional exceedance rate (ACER). The package includes useful functions for data preparations, analysis and presentations of results.
}
\details{
\tabular{ll}{
Package: \tab ACER\cr
Type: \tab Package\cr
Title: \tab Average Conditional Exceedance Rate Method as a Tool for Extreme Value Statistics\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-20\cr
Author: \tab Kristoffer Kofoed Roedvei\cr
Maintainer: \tab Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>\cr
Description: \tab ACER is a package made for calculating average conditional exceedance rate (ACER). \cr
\tab The package includes useful functions for data preparations, analysis and presentations\cr
\tab of results.\cr
License: \tab GPL-3\cr
Depends: \tab minpack.lm\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\references{
For details about how to get starting useing the ACER package see examples below, \code{\link{ACER}} and \code{\link{ACERm}}.\cr \cr
More about the ACER method see \url{http://downloads.hindawi.com/journals/jps/2013/797014.pdf}. \cr
Some of the ACER R package functions is based on the MATLAB program \url{http://folk.ntnu.no/karpa/ACER/ACER_User_guide.pdf}.
}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}

\seealso{
\code{\link{ACER}} \cr 
\code{\link{ACERm}} \cr
\code{\link{CI.ACER}} \cr
\code{\link{Haltenbanken}} \cr
\code{\link{internalACER}} \cr
- \code{\link{eta1}} \cr
- \code{\link{RealizationByRows}} \cr
- \code{\link{RealizationLength}} \cr
- \code{\link{findcond}} \cr
- \code{\link{Gumbel}} \cr
- \code{\link{reverseGumbel}} \cr
- \code{\link{General}} \cr
- \code{\link{reverseGeneral}} \cr
\code{\link{internal-MSEoptimization}} \cr
- \code{\link{Guess}} \cr
- \code{\link{GeneralOptim}} \cr
- \code{\link{GeneralParamf}} \cr
- \code{\link{GeneralParamw}} \cr
- \code{\link{mseGeneral3f}} \cr
- \code{\link{mseGeneral3w}} \cr
- \code{\link{mseGeneral5f}} \cr
- \code{\link{mseGeneral5w}} \cr
- \code{\link{GradmseGeneral5f}} \cr
- \code{\link{GradmseGeneral5w}} \cr
- \code{\link{GumbelOptim}} \cr
- \code{\link{GumbelParam}} \cr
- \code{\link{mseGumbel2}} \cr
- \code{\link{mseGumbel4}} \cr
- \code{\link{GradmseGumbel4}} \cr 
\code{\link{kplot}} \cr
\code{\link{MSEoptimization}} \cr
\code{\link{peak}} \cr
\code{\link{plot.ACER}} \cr
\code{\link{predict.ACER}} \cr
\code{\link{print.ACER}} \cr
\code{\link{summary.ACER}} \cr
\code{\link{print.summary.ACER}} \cr
\code{\link{updatecond}} \cr
}

\examples{
\dontrun{
##################
# Gumbel Example #
##################
ar<-c(11/60,7/10)
set.seed(1989)
dataAR<-arima.sim(n = 131400, list(ar = ar))
#make it into 15 realizations
dataAR<-matrix(dataAR,nrow=15,byrow=TRUE)
library(ACER)
kplot(dataAR,k=1:5,xlim=c(0.5,max(dataAR)),main='k-plot')
#With Zoom plot
kplot(dataAR,k=1:5,xlim=c(3,5), main='Zoomed k-plot')
X<-ACER(dataAR,k=3,method='general',eta1=sort(as.vector(t(dataAR)))[131400-10000])
plot(X,main='Fitted plot with real data (transparent)')
# print the xi value. Gumbel is when xi get close to zero
summary(X)
(xi<-summary(X)$coef['Estimate','xi'])
# predict the 100 year senario
(predCI<-predict(X,prob=1/(100*365*24)))

# calculate 100 numbers of the 100 year return to compare it to the prediction
max100<-rep(NA,100)
for(i in 1:100){
  max100[i]<-max(arima.sim(n = 131400, list(ar = ar)))
}
mean(max100)
(realCI<-sort(max100)[c(5,95)])
#
#
###################
# Frechet Example #
###################
library(Quandl)
X<-Quandl('GOOG/NASDAQ_MSFT')

# Organize the data into the correct years (realization)
end<-FALSE
count<-1
leng<-length(X[,3])
M<-matrix(NA,28,260)
i<-1
j<-1
while(count<leng){
  year<-as.numeric(format(X[count,1],'\%Y'))
  if(year==2015||year==1986){
    count<-count+1
  }else{
    if(j<(2015-year)){
      i<-1
      j<-2015-year
    }
    if(is.na(X[count,2])){
      M[j,i]<-0
    }else{M[j,i]<-(X[count,5]-X[count-1,5])/X[count-1,5]}
    count<-count+1
    i<-i+1
  }
}
# M contain the daily return by column and realization by rows
Mvec<-as.vector(t(M))
Mvec<-Mvec[!is.na(Mvec)]
# Mvec is an vector of the return (easier to plot than a matrix)
ts.plot(Mvec,main='Timeseries',ylab='daily return')
# Mvec is non stationary in variance, as we can see from the plot
#ACER analysis
library(ACER)
kplot(M,main='k-plot')
# Looks like the acer's are indepentent of k in k-plot.
msft<-ACER(M,stationary=FALSE,method='general',eta1=-0.01)
plot(msft,prob=1e-05,main='Fitted plot with real data (transparent)')
summary(msft)
(xi<-summary(msft)$coef['Estimate','xi'])
#
#
#####################################################
# Weibull Example                                   #
# which does not work with the Subasymptotic method #
# shows how to identify such problems.              #
# should be corected in later versions              #
#####################################################
library(ACER)
set.seed(1989)
dataRN<-rnorm(131400)
dataRN[dataRN>3.55]<-NA
#make it into 15 realizations
dataRN<-matrix(dataRN,nrow=15,byrow=TRUE)
sum(is.na(dataRN))
# Since 23 << 131400, the 23 data is just removed without replacement
# without significant loss of information.
peakdataRN<-peak(dataRN)
kplot(peakdataRN,main='k-plot')
X<-ACER(peakdataRN,method='general',eta1=sort(as.vector(t(dataRN)))[131400-3000],CI=0.70)
plot(X,main='Fitted plot with real data (transparent)')
summary(X)
#
#
#####################
# Practical Example #
#####################
library(ACER)
data(Haltenbanken)
dim(Haltenbanken)
# Realization by row and data by column
hist(Haltenbanken,breaks=50)
ts.plot(apply(Haltenbanken,2,mean),main='Average at each time, over years',
        xlab='Time',ylab='Haltenbanken mean value')
# plots the mean at the same time over each year.
# The plot clearly shows the seasonality, 
# which means that the data is nonstationary
Hpeak<-peak(Haltenbanken)
# peak extracts peaks
kplot(Hpeak,k=1:8,stationary=FALSE,main='k-plot')
# Loks like the lines converge at arround 8.
plot(ACERm(Hpeak,k=8,stationary=FALSE),add.eta1=FALSE,
     main='ACER plot with confidence interval')
# eta1=7 seems like a good chose.
# prediction plot using 4 different CI interval
par(mfrow=c(2,2))
for(i in c(0.65,0.75,0.85,0.95)){
  Hacer<-ACERm(Hpeak,k=8,stationary=FALSE)
  Hacer<-CI(Hacer,level=i)
  Hacer<-updatecond(Hacer,eta1=7)
  Hacer<-MSEoptimization(Hacer,method='gumbel')
  plot(Hacer,main=paste(i*100,'\% CI interval', sep=''))
}
par(mfrow=c(1,1))
# Hacer is saved with the last value (CI.level=0.95)
summary(Hacer)
predict(Hacer,prob=1e-06)
predict(Hacer,eta=17)
}
}

