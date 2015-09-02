\name{predict.ACER}
\alias{predict.ACER}
%- Also NEED an '\alias' for EACH other topic documented here.
  \title{Makes an Extrapolation of the ACER}
\description{
Makes an extrapolation of the eta or the probability for an ACER class. The prediction is made on eta if prob is given and vice versa.
}
\usage{
\method{predict}{ACER}(object, eta = NULL, prob = NULL, 
        weibull.limit = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of "\code{ACER}" calss.}
  \item{eta}{ numeric or vector of numbers. The number indcate the barrier level of where the extrapolation should be preformed.}
  \item{prob}{ numeric or vector of numbers each element should be between 0 and 1. The number indicate the probability or epsilon of where the extrapolation should be preformed.}
  \item{weibull.limit}{ logical. When \code{TRUE} and \code{X} is Weibull distributed, \code{eta} and \code{prob} are NULL, \code{predict} returns the Weibull eta limits.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}
\details{
predict uses the Gumbel or general model to extrapolate the barrier level or probability when the other is given.
}
\value{
Returns a matrix with the predicted values, confidence intervalls and the known input for prediction.\cr \cr
If \code{weibull.limit} is \code{TRUE}, \code{X} is Weibull distributed, \code{eta} and \code{prob} are NULL, the function returns upper confidence, estimated and lower confidence limits for each curves. This is possible since Weibull distributions are limited.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
(\code{\link{print.ACER}}, \code{\link{summary.ACER}}, \code{\link{print.summary.ACER}}).
}
\examples{
data(Haltenbanken)
X<-peak(Haltenbanken)
\donttest{Y<-ACER(X,k=7,stationary=FALSE)}\dontshow{Y<-ACER(X,k=7,stationary=FALSE,eta1=6.613026)} %check fails since user input is required for this line
# choose an eta1, example eta1 = 6.613026 like the program suggest
predict(Y,eta=20)
predict(Y,prob=1e-6)
# Multiple extrapolation
predict(Y,eta=c(17,18,19,20,21,22,23))
predict(Y,prob=c(1e-5,1e-6,1e-7))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
