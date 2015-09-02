\name{MSEoptimization}
\alias{MSEoptimization}
\alias{MSEoptimization.ACER}
\title{MSE Optimization of ACER Through GEV or Gumbel Distribution}
\description{
Uses weighted mean sqare error (MSE) to fit a subasymptotic Gumbel or General Extreme Value (GEV) distribution. The return of \code{MSEoptimization} is an updated ACER class which include GEV or Gumbel coefficients for the ACER and confidence intervals.
}
\usage{
MSEoptimization(X,...)
\method{MSEoptimization}{ACER}(X, method = "general", check.weibull = TRUE, 
                penalty = FALSE, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{an \code{ACER} class, see \code{\link{ACER}} or \code{\link{ACERm}}.}
  \item{method}{string. Either \code{general} or \code{gumbel} dependent on method which the ACER should be optimized by. see details.}
  \item{check.weibull}{ logical. Will also look for Weibull solution when \code{check.weibull} is \code{TRUE} and \code{method} is \code{general}.}
  \item{penalty}{logical. Can only be used in the gumbel case (see details). When \code{TRUE} penalty is used, else normal weighted MSE is preformed. For the gumbel case when \code{c} get close to 1, \code{q} becomes unnecessary and we have an over fitted model (see details). Penalty option is added with the goal of reducing this behavior in the certain cases.}
  \item{alpha}{numeric. alpha is the penalty factor between \code{0} or \code{1}. see details.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}
\details{
Equation fitted for Gumbel distribution (\code{method="gumbel"}):\cr
\deqn{\epsilon(\eta)=q*exp[-a(\eta-b)^c]}
In the penalty case, the function get an extra penalty factor as you can see here:
\deqn{\epsilon_p(\eta)=exp(\alpha*c^{[sign(c-1)-1]}) * \epsilon(\eta)}
When \code{c} gets close to \code{1} in the gumbel case, this could be useful to provide over fitting.\cr \cr \cr
Equation fitted for GEV distribution (\code{method="general"}):\cr
\deqn{\epsilon(\eta)=q*[1+a \xi(\eta-b)^c]^{-1/\xi}}
Note that the gumbel distribution is a subset of the general distribution. This can be seen when \eqn{\xi} goes to zero, \eqn{\epsilon(\eta)} goes toward the gumbel case. See \url{http://downloads.hindawi.com/journals/jps/2013/797014.pdf}\cr \cr
\code{MSEoptimization} rely on internal functions, see \code{\link{internal-MSEoptimization}}.\cr
}
\value{
Return an ACER class:
\item{coef}{coefficients for the chosen model. 4 or 5 coefficients respectively if method is gumbel or general }
\item{upperCIcoef}{coefficients for the chosen model. 4 or 5 coefficients respectively if method is gumbel or general }
\item{lowerCIcoef}{coefficients for the chosen model. 4 or 5 coefficients respectively if method is gumbel or general }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{updatecond.ACER}}\cr
\code{\link{ACER}} \cr
For further package usage see \code{\link{plot.ACER}}\cr
(\code{\link{predict.ACER}}, \code{\link{print.ACER}}, \code{\link{summary.ACER}}, \code{\link{print.summary.ACER}}).
}
\examples{
# Load the data set
data(Haltenbanken)
# 1. example
X<-peak(Haltenbanken)
Y<-ACERm(X,k=5,stationary=FALSE)
\donttest{plot(Y)}\dontshow{Y<-CI(Y,level=0.95);Y<-updatecond(Y,eta1=6.087174)}%check fails since user input is required for this line
# choose an eta1, example eta1 = 6.087174 like the program suggest
Y<-MSEoptimization(Y,method="general")

# 2. example beginning of regular tail behavior (eta1) is known
# from previous calculations or other reasons.
# Here we have eta1 approximately equal to 7
X<-peak(Haltenbanken)
Y<-ACERm(X,k=5)
# Add confidence level to the ACER object
Y<-CI(Y,level=0.95)
Y<-updatecond(Y,eta1=7)
Y<-MSEoptimization(Y,method="gumbel",penalty=TRUE,alpha=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
