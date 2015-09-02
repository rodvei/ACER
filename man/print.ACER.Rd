\name{print.ACER}
\alias{print.ACER}
\title{ printing ACER Fits}
\description{
\code{print} method for class "ACER".
}
\usage{
\method{print}{ACER}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of calss "\code{ACER}", usually, a result of a call to \code{\link{ACER}} or \code{\link{ACERm}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
prints a format of coefficients. If \code{\link{ACERm}} is used, the output will try to suggest the next logical step in the process.
}
\value{
No returns
}

\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
The model fitting functions \code{\link{ACER}}, \code{\link{ACERm}}.\cr
Similar functions \code{\link{summary.ACER}},\code{\link{print.summary.ACER}}.
}
\examples{
data(Haltenbanken)
X<-peak(Haltenbanken)
\donttest{Y<-ACER(X,k=7,stationary=FALSE)} \dontshow{Y<-ACER(X,k=7,stationary=FALSE,eta1=6.613026)} %check fails since user input is required for this line
# choose an eta1, for example eta1 = 6.613026 like the program suggest
print(Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
