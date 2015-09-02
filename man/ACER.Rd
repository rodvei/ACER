\name{ACER}
\alias{ACER}
\alias{ACER.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Automatically Constructing a Complete ACER class}
\description{
The \code{ACER} function calculate the average conditional exceedance rate, together with all the variables for a complete ACER class. Probably the most useful function, since it is fully automatic and the easiest to use for ACER analysis. The class is ready for prediction, plotting and analysis. The function is build up by the manual functions  \code{\link{ACERm}}, \code{\link{plot.ACER}}, \code{\link{CI.ACER}},\code{\link{updatecond.ACER}}, \code{\link{MSEoptimization}, ect.}
}
\usage{
ACER(X, ...)
\method{ACER}{default}(X, k = 1, CI = 0.95, eta1 = NULL, stationary = TRUE, 
     method = "general", check.weibull = TRUE, penalty = FALSE, 
     alpha = 0.05, neta = 500,...)
}
\arguments{
  \item{X}{ a matrix or vecotor including data to be analysed. X is often the output from \code{\link{peak}} since ACER becomes costly for large data sets. Realizations should be by rows. See \code{\link{ACERm}} for more details.}
  \item{k}{numeric, indicating the k-step memory approximations, which by default is equal to 1 corresponding to all data point being analyzed as if they where independent of each other. See details, \code{\link{kplot}} and \code{\link{ACERm}} for more in-depth information.}
  \item{CI}{numeric, between 0 and 1. Indicates the confidence level. See \code{\link{CI}}.}
  \item{stationary}{ logical. If \code{TRUE} the data is analyzed as a stationary data set. If \code{FALSE} the data is considered nonstationary. See \code{\link{ACERm}} for more details.}
  \item{method}{string. Either \code{general} or \code{gumbel} dependent on which method ACER should be optimized by. See \code{\link{MSEoptimization}}.}
  \item{eta1}{ numeric, barrier level for where regular tail behavior starts. If not included, plot and suggestions is presented for easier choice of eta1. See \code{\link{updatecond}}}
  \item{check.weibull}{ logical. Will also look for Weibull solution when \code{check.weibull} is \code{TRUE} and \code{method} is \code{general}. See \code{\link{MSEoptimization}}.}
  \item{penalty}{logical. Can only be used in the gumbel case. When \code{TRUE} penalty is used, else normal weighted MSE is preformed. For the gumbel case when c get close to 1, q becomes unnecessary and we have an over fitted model. Penalty option is added with the goal of reducing this behavior in the certain cases. See \code{\link{MSEoptimization}}.}
  \item{alpha}{numeric. \code{alpha} is the penalty factor. see \code{\link{MSEoptimization}}.}
  \item{neta}{ numeric. The numbers of barrier levels considered from \code{min} to \code{max} when calculating ACER.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}
\details{
The function automatically preforms the average conditional exceedance rate (see \code{\link{ACERm}}, and \url{http://downloads.hindawi.com/journals/jps/2013/797014.pdf}) and MSE optimization of a subasymptotic Gumbel or GEV with confidence lines (see \code{\link{MSEoptimization}}).\cr
\code{k} indicates the numbers of previous observations which should be considered when observing independent exceedence for a certain barrier level.
}
\value{
Returns an updated and fully capable ACER class.
\item{eta}{array of barrier levels.}
\item{acer}{array of average conditional exceedance rates of the corresponding \code{eta}'s.}
\item{acer.sd}{array of standard deviations of the corresponding \code{eta}'s.}
\item{eta1est}{a numeric estimate for the beginning of regular tail behavior if \code{eta1=TRUE} in call.}
\item{k}{integer from call.}
\item{dim}{the dimension of the data set \code{X}}
\item{RL}{the array of realization lengths}
\item{CI}{ an array of \code{1/2} confidence width. One confidence element for each \code{acer}.}
\item{CIlevel }{confidence interval level}
\item{StartEnd}{an array of two elements indicating the start and end position of interval of interest.}
\item{coef}{coefficients for the choose model. 4 or 5 coefficients respectively if method is gumbel or general }
\item{upperCIcoef}{coefficients for the upper confidence line. 4 or 5 coefficients respectively if method is gumbel or general }
\item{lowerCIcoef}{coefficients for the lower confidence line. 4 or 5 coefficients respectively if method is gumbel or general }
}
\references{
For more in depth description of the math and statistics behind the ACER see chapter \url{http://downloads.hindawi.com/journals/jps/2013/797014.pdf}. \cr
The \code{ACER} function is based on the functions from the \url{http://folk.ntnu.no/karpa/ACER/ACER_User_guide.pdf}.
}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
For more information and examples see \code{\link{ACERm}}, \code{\link{MSEoptimization}}. \cr
For further package usage see \code{\link{plot.ACER}}\cr 
(\code{\link{predict.ACER}}, \code{\link{print.ACER}}, \code{\link{summary.ACER}},\code{\link{print.summary.ACER}}).
}
\examples{
data(Haltenbanken)
X=peak(Haltenbanken)
\donttest{Y=ACER(X,k=5,stationary=FALSE)} %check fails since user input is required for this line
# choose an eta1, for example eta1 = 6.087174 like the program suggest
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
