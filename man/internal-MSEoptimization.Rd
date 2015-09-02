\name{internal-MSEoptimization}
\alias{internal-MSEoptimization}
\alias{Guess}
\alias{GeneralOptim}
\alias{GeneralParamf}
\alias{GeneralParamw}
\alias{mseGeneral3f}
\alias{mseGeneral3w}
\alias{mseGeneral5f}
\alias{mseGeneral5w}
\alias{GradmseGeneral5f}
\alias{GradmseGeneral5w}

\alias{GumbelOptim}
\alias{GumbelParam}
\alias{mseGumbel2}
\alias{mseGumbel4}
\alias{GradmseGumbel4}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal function used by MSEoptimization}
\description{
\code{Guess} calculate a matrix containing four or six starting values if \code{method} is \code{gumbel} or \code{general} respectively.\cr \cr
\code{GeneralOptim} returns the coefficients which minimizes the sum of squares error for the GEV distribution, see \code{\link{MSEoptimization}}.\cr
\code{GeneralParam(f/w)} uses the three first paramters as input and returns the five parameters used in the general case. The f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}<0 (Weibull). \cr 
\code{mseGeneral3(f/w)} returns the weighted sum of squares error using the three first parameters since \eqn{\xi} and q are linearly dependent of these. The f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}<0 (Weibull).\cr 
\code{mseGeneral5(f/w)} returns the weighted sum of squares error using all five parameters. The f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}<0 (Weibull).\cr 
\code{GradmseGeneral5(f/w)} returns the gradient of the weighted sum og squares error. The f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}<0 (Weibull).\cr 
\cr
\code{GumbelOptim} returns the coefficients which minimizes the sum of squares error for the Gumbel distribution, see \code{\link{MSEoptimization}}.\cr 
\code{GumbelParam} uses b and c paramters as input and returns the four parameters used in the Gumbel case.\cr
\code{mseGumbel2} returns the weighted sum of squares error using the b and c parameters since a and q are linearly dependent of these.\cr 
\code{mseGumbel4}  returns the weighted sum of squares error using all four parameters.\cr 
\code{GradmseGumbel4} returns the gradient of the weighted sum og squares error.
}


\usage{
Guess(eta,logeps, w, eta1, mineta, maxeta, method, check.weibull=TRUE)

GeneralOptim(guess, eta, logeps, w, eta1, mineta, maxeta)
GeneralParamf(est,eta,logeps,w)
GeneralParamw(est,eta,logeps,w)
mseGeneral3f(est,eta,logeps,w)
mseGeneral3w(est,eta,logeps,w)
mseGeneral5f(est,eta,logeps,w)
mseGeneral5w(est,eta,logeps,w)
GradmseGeneral5f(est,eta,logeps,w)
GradmseGeneral5w(est,eta,logeps,w)

GumbelOptim(guess, eta, logeps, w, eta1, mineta, penalty=FALSE,alpha=0.05)
GumbelParam(est,eta,logeps,w)
mseGumbel2(est,eta,logeps,w,penalty=FALSE,alpha=0.05)
mseGumbel4(est,eta,logeps,w,penalty=FALSE,alpha=0.05)
GradmseGumbel4(est,eta,logeps,w)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{vector, with the barrier levels of the elements position.}
  \item{logeps}{vector of log(\eqn{\epsilon}).}
  \item{est}{ vector, with the estimated coefficients (a, b, c, q, \eqn{\xi}). Dependent on which function it's used on \code{est} has length 2, 3, 4 or 5.}
  \item{w}{ vector, with the weight of the elements position.}
  \item{eta1}{ numeric, which indicate the beginning of regular tail behavior.}
  \item{mineta}{ numeric, equals the smallest measured barrier level.}
  \item{maxeta}{ numeric, equals the largest measured barrier level.}
  \item{method}{ string. Either \code{general} or \code{gumbel} dependent on method which ACER should look for start values.}
  \item{check.weibull}{logical. Will also look for Weibull solution when \code{check.weibull} is \code{TRUE} and \code{method} is \code{general}.}
  \item{guess}{ matrix, with starting values by row's. One row would be \code{c(a,b,c,q)} or \code{c(a,b,c,q,\eqn{\xi})} dependent on mehtod.}
  \item{penalty}{ logical. \code{TRUE} if panalty should be used.}
  \item{alpha}{ numeric, the penalty factor. Should be between \code{0} and \code{1}.}
 
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{Guess}, returns a matrix of four or six guesses.\cr
\code{GeneralOptim}, returns the best coefficients achieved for GEV.\cr
\code{GeneralParam(f/w)}, returns an estimated of the GEV paramteres. Where f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}<0 (Weibull).\cr
\code{mseGeneral3(f/w)}, returns the weighted sum of squares error. Where f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}<0 (Weibull).\cr
\code{mseGeneral5(f/w)}, returns the weighted sum of squares error. Where f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}).\cr
\code{GradmseGeneral5(f/w)}, returns the gradient of the weighted sum of squares error. Where f indicate \eqn{\xi}>0 (Frechet), while w indicate \eqn{\xi}).\cr
\cr
\code{GumbelOptim}, returns the best coefficients achieved for Gumbel.\cr
\code{GumbelParam}, returns an estimated of the Gumbel paramteres.\cr
\code{mseGumbel2},  returns the weighted sum of squares error.\cr 
\code{mseGumbel4},  returns the weighted sum of squares error.\cr 
\code{GradmseGumbel4}, returns the gradient of the weighted sum og squares error for the Gumbel.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MSEoptimization}}
}
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
