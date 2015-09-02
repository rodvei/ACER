\name{updatecond}
\alias{updatecond}
\alias{updatecond.ACER}

\title{ Update Area of Interest For the ACER Class}
\description{
update \code{StartEnd} condition for the ACER class. When using the \code{\link{ACER}} function, \code{updatecond} is not needed, but users who manually adds confidence interval (\code{CI}) and/or point of where regular tail behavior begins (\code{eta1}) should always run \code{updatecond} afterword. Users who manually use \code{ACERm}, but adds \code{CI} and \code{eta1} through \code{plot} does not need \code{updatecond}.
}
\usage{
updatecond(X,...)
\method{updatecond}{ACER}(X, eta1=NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{an ACER class}
  \item{eta1}{barrier level for where regular tail behavior starts.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}

\value{
Returns a updated version of the ACER class.
\item{StartEnd}{an array of two elements indicating the start and end position of interval of interest.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
\code{\link{ACERm}}\cr
\code{\link{CI.ACER}}\cr
\code{\link{ACER}} \cr
For further package usage see \code{\link{plot.ACER}}, \code{\link{MSEoptimization}}\cr
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
# After adding the Confidence level we need suitable area
Y<-updatecond(Y)

# 2. example beginning of regular tail behavior (eta1) is known
# from previous calculations or other reasons.
# Here we have eta1 approximately equal to 7
X<-peak(Haltenbanken)
Y<-ACERm(X,k=5)
# Add confidence level to the ACER object
Y<-CI(Y,level=0.95)
Y<-updatecond(Y,eta1=7)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
