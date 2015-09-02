\name{ACERm}
\alias{ACERm}
\alias{ACERm.default}
\title{ Manually Calculate the Average Conditional Exceedance Rate (ACER)}
\description{
The \code{ACERm} function calculates the average conditional exceedance rate it's corresponding standard divination and an estimate of the barrier level where regular tail behavior starts. The function returns a \code{ACER} class which is constructed with some extra empty variables which will be added by other functions in the package. 
}
\usage{
ACERm(X, ...)
\method{ACERm}{default}(X, k = 1, stationary = TRUE, eta1 = TRUE, neta = 500, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ matrix or vecotor including data to be analysed. \code{X} is often the output from \code{\link{peak}} since \code{ACERm} becomes costly for large data sets.}
  \item{k}{ numeric, indicating the k-step memory approximations, which by default is equal to 1 corresponding to all data point being analyzed as if they where independent of each other. See details and \code{\link{kplot}}.}
  \item{stationary}{ logical. If \code{TRUE} the data is analyzed as a stationary data set. If \code{FALSE} the data is considered nonstationary.}
  \item{eta1}{ logical. If \code{TRUE} the function will estimate the beginning of regular tail behavior.}
  \item{neta}{ numeric. The numbers of barrier levels considered from \code{min} to \code{max} when calculating ACER.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}
\details{
\code{k} indicates the numbers of previous observations which should be considered when observing independent exceedence for a certain barrier level.\cr \cr
The \code{eta1} estimate should only be used as indication, further plots and analyzes should be considered before concluding the value of \code{eta1}. When \code{eta1=TRUE}, 
later plots will mark the estimate of \code{eta1}, which could be helpful for locating the beginning of regular tail behavior.\cr \cr
\code{neta} is adjustable for special usage when the function is to slow, or extra accuracy is needed, but for the average user 500 should be fine.\cr \cr
if any of the description wasn't satisfactory see references.
}
\value{
Returns an ACER class. Some of the returns are not important for the user. Extra data is added to the ACER class by later functions.
%%  If it is a LIST, use
\item{eta}{array of barrier levels.}
\item{acer}{array of average conditional exceedance rates of the corresponding \code{eta}'s.}
\item{acer.sd}{array of standard deviations of the corresponding \code{eta}'s.}
\item{eta1est}{a numeric estimate for the beginning of regular tail behavior if \code{eta1=TRUE} in call.}
\item{k}{integer from call.}
\item{dim}{the dimension of the data set \code{X}}
\item{RL}{the array of realization lengths}

%% ...
}
\references{
For more in depth description of the math behind the ACER see chapter \url{http://downloads.hindawi.com/journals/jps/2013/797014.pdf}. \cr
The \code{ACER} function is based on the functions from the \url{http://folk.ntnu.no/karpa/ACER/ACER_User_guide.pdf}.
}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
\code{\link{ACER}} \cr
For further package usage see \code{\link{plot.ACER}}, \code{\link{CI.ACER}}, \code{\link{updatecond.ACER}}, \code{\link{MSEoptimization}}\cr
(\code{\link{predict.ACER}}, \code{\link{print.ACER}}, \code{\link{summary.ACER}}, \code{\link{print.summary.ACER}}).
}
\examples{
# Load the data set
data(Haltenbanken)
#
# 1. example
X<-Haltenbanken[1:18,1:1000] 
# Haltenbanken is to large for reasonable fast  computation
# thats why X<-Haltenbanken[1:18,1:1000] is used instead of
# just X<-Haltenbanken
ACERm(X)
#
# 2. example
# A better way of making the data of Haltenbanken smaller 
# whithout loss of extremes, is by useing peak():
X<-peak(Haltenbanken)
Y<-ACERm(X)
#
# 3. example, non stationary and k=5
plot(colSums(Haltenbanken,na.rm=TRUE),type='l')
# We see that the time series is non stationary and we should instead do
X<-peak(Haltenbanken)
Y<-ACERm(X,k=5,stationary=FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
