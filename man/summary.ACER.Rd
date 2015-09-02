\name{summary.ACER}
\alias{summary.ACER}
\alias{print.summary.ACER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summarizing ACER Fits}
\description{ \code{summary} method for class "ACER" and \code{print} method for the summary.}
\usage{
\method{summary}{ACER}(object, ...)

\method{print}{summary.ACER}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of calss "\code{ACER}", usually, a result of a call to \code{\link{ACER}} or \code{\link{ACERm}}.}
  \item{x}{an object of calss "\code{summary.ACER}", usually, a result of call to \code{summary.ACER}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{print.summary.ACER prints a format of coefficients, confidence intervall and other informations about the object. If \code{\link{ACERm}} is used, the output will try to suggest the next step in the process.}
\value{
\code{summary.ACER} returns a list of calculated values useful for further analysis.
\item{call}{the last function call on the "\code{ACER}" class.}
\item{eta}{vector of the barrier levels (\eqn{\eta}) within the applicable area, meaning eta>eta1 and where lower confidence interval of acer is positive. For mor details see \code{\link{updatecond}}}
\item{acer}{vector of the average conditional exceedance rate (ACER) within the applicable area, see \code{eta} abow.}
\item{k}{the number used in the k-step memory approximations.}
\item{coefficients}{ matrix holding the coefficients for the estimated model, upper and lower confidence intervall. In the coefficient matrix there is also included the value of weighted mean sqare error (\code{W.MSE}) for each of the model estimated.}
\item{data}{vector of two elements where the first is the numbers of datapoint and the second is the numbers of realization used to calculate the ACER.}
\item{t.dist}{logical, if \code{TRUE} t-distribution is used to calculate the confidence interval. If \code{FALSE} poisson is used.}
\item{eta1}{the eta where the regular tail behavior start (as set by the user).}
\item{numobs}{the number of suitable ACER's used in the weighted MSE optimization. }
\item{equation}{string including the equation which the coefficients are optimized for.}
\item{df}{The degree of freedom for the t-distribution used to calculate the confidence interval.}


%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
The model fitting functions \code{\link{ACER}}, \code{\link{ACERm}}.\cr
Similar function \code{\link{print.ACER}}.
}
\examples{
data(Haltenbanken)
X<-peak(Haltenbanken)
\donttest{Y<-ACER(X,k=7,stationary=FALSE)} \dontshow{Y<-ACER(X,k=7,stationary=FALSE,eta1=6.613026)} %check fails since user input is required for this line
# choose an eta1, for example eta1 = 6.613026 like the program suggest
smry<-summary(Y)
smry
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
