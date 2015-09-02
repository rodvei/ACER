\name{CI}
\alias{CI}
\alias{CI.ACER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update Confidence Interval in ACER Class}
\description{
\code{CI} is often used after \code{\link{ACERm}}. The function is used for updating the confidence levels in the \code{ACER} class. The function is not necessary for most user since they directly use \code{\link{ACER}} or the confidence intervals are added by \code{plot.ACER}. But in cases when time consuming calculations have already been preformed and the user only want to change the confidence interval, the function can be useful, see \code{\link{ACER-package}} for practical example of this.\cr
\code{CI} should be used with a following \code{\link{updatecond}} to update the \code{StartEnd} condition.
}
\usage{
CI(X,...)
\method{CI}{ACER}(X, level = 0.95, tdist.CI = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{an \code{ACER} class, see \code{\link{ACER}} and \code{\link{ACERm}}.}
  \item{level}{the confidence level.}
  \item{tdist.CI}{ logical. If FALSE poisson approximation is used to calculate the confidence intervals.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}
\details{
\code{level} should be decimal (between 0 and 1) or percentage (between 0 and 100).\cr \cr
\code{tdist.CI} indicate if poisson or t-distribution is used to calculate the confidence interval. t-distribution uses the variation between realization to arrive at a confidence interval, while poisson dose not and is therefor automatically used for single realization cases.
}
\value{
returns an updated version of the \code{ACER} class. New elements added are:
\item{CI}{ an array of \code{1/2} confidence width. One confidence element for each \code{acer}.}
\item{CIlevel }{numeric, confidence interval level, equaling the \code{level} argument.}
}
\references{
For more details of how confidence interval is calculated using t-distribution or poisson see \url{http://downloads.hindawi.com/journals/jps/2013/797014.pdf} 
}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
\code{\link{ACERm}}\cr
\code{\link{ACER}} \cr
For further package usage see \code{\link{plot.ACER}}, \code{\link{updatecond.ACER}}, \code{\link{MSEoptimization.ACER}}\cr
(\code{\link{predict.ACER}}, \code{\link{print.ACER}}, \code{\link{summary.ACER}}, \code{\link{print.summary.ACER}}).
}
\examples{
# Load the data set
data(Haltenbanken)
#
# 1. example 
X<-peak(Haltenbanken)
Y<-ACERm(X)
# Add confidence level to the ACER object
Y<-CI(Y,level=0.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
