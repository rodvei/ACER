\name{plot.ACER}
\alias{plot.ACER}
\title{Plots an ACER Class}
\description{
Plots an object of ACER class. Makes a graphical plot of the average conditional exceedance rate vs barrier level. Good for estimating and adding \code{eta1} to the object. Can also plot extrapolation lines together with the ACER plot if \code{\link{ACER}} or \code{\link{MSEoptimization}} have been run first. If the user don't use \code{\link{ACER}} but still don't want to manually run \code{\link{CI}} and \code{\link{updatecond}}, \code{plot} will add confidence interval to the object and make an update. \cr
Black solid line represents the data, while the blue dashed lines are the confidence interval. When extrapolation is preformed the colors is the same for the data and estimated line, but the real data is more transparent than the fitted lines.
}
\usage{
\method{plot}{ACER}(x, conf.int = TRUE, CI.level, extrapolate, prob = 1e-6,
     add.eta1 = TRUE, xlab = expression(eta), 
     ylab = expression(hat(epsilon)[k](eta)), type = 'l', xlim, 
     ylim, alpha = 1, tdist.CI = TRUE, ...)
}
\arguments{
  \item{x}{an \code{ACER} class, see \code{\link{ACER}} or \code{\link{ACERm}}}
  \item{conf.int}{logical. Indicate if confidence interval should be plotted.}
  \item{CI.level}{numeric, the value of the confidence level, see \code{\link{CI}}.}%make CI.level=0.95 to just CI.level??
  \item{extrapolate}{logical. Whether the extrapolation is desired or not. By defualt \code{TRUE} if coeficients have been optimized by \code{\link{MSEoptimization}} or \code{\link{ACER}}.}
  \item{prob}{numeric, the disired level you want ACER to be graphically extrapolated to (defualt=1e-6)}
  \item{add.eta1}{logical. If the \code{ACER} class does not hold a value for eta1 and \code{add.eta1=TRUE}, after plotting the \code{ACER} the program will wait for the user to insert an estimate of eta1 (based on the plot) and add it to the \code{ACER} class.}
  \item{xlab}{string. Title for x-axis, see \code{\link{plot}}.}
  \item{ylab}{string. Title for y-axis, see \code{\link{plot}}.}
  \item{type}{string indicating the type of plot, see \code{\link{plot}}.}
  \item{xlim}{array of length 2 indicating min and max for x-labels, see \code{\link{plot}}.}
  \item{ylim}{array of length 2 indicating min and max for y-labels, see \code{\link{plot}}.}
  \item{alpha}{number between 0 and 1 indicating how transparent the plot should be. Where 0 is fully transparent and 1 is solid.}
  \item{tdist.CI}{logical. If object don't have confidence interval \code{plot} will run \code{\link{CI}} to add this.}
  \item{\dots}{additional arguments to be passed to the built in plot function, see \code{\link{plot}}}
}
\value{
The function don't return values. The function can still add eta1 to the class without return of the class, but this will be clear by the question asked in the console.
}
\examples{
data(Haltenbanken)
# 1. example
X<-peak(Haltenbanken)
\donttest{Y1<-ACER(X,5)} \dontshow{Y1<-ACER(X,5,eta1=6.087174)}%check fails since user input is required for this line
# choose an eta1, for example eta1 = 6.087174 like the program suggest
plot(Y1)

# Using the manual ACERm function
# 2. example 
X<-peak(Haltenbanken)
Y2<-ACERm(X,5)
\donttest{plot(Y2)} \dontshow{Y2<-CI(Y2,level=0.95);Y2<-updatecond(Y2,eta1=8)}%check fails since user input is required for this line
# Inset number for where regular tail behavior starts eksample eta1=8
plot(Y2)
#
# 3. example
Y2=MSEoptimization(Y2)
plot(Y2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
