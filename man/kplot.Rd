\name{kplot}
\alias{kplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting ACER for Different k-dependency}
\description{
\code{kplot} plots the ACER by barrier levels for different \code{k} dependency. The plot contain multiple \code{k} dependency in the same figure for comparison. The comparison should make it easier to choose a suitable \code{k} value for further analysis.
}
\usage{
kplot(X, k = 1:6, stationary = TRUE, 
      neta = 500, xlab = expression(eta), 
      ylab = expression(hat(epsilon)[k](eta)),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix or vecotor including data to be analysed. \code{X} is often the output from \code{\link{peak}} since \code{kplot} becomes costly for large data sets.}
  \item{k}{vector or numeric, containing every \code{k} for the k-step memory approximations which should be compared. See details.
}
  \item{stationary}{logical. If \code{TRUE} the data is analyzed as a stationary data set. If \code{FALSE} the data is considered nonstationary.}
  \item{neta}{integer. The numbers of barrier levels considered from \code{min} to \code{max} when calculating ACER.}
  \item{xlab}{string with a title for the x axis.}
  \item{ylab}{string with a title for the y axis.}
  \item{\dots}{additional arguments to be passed to lower level functions.}
}
\details{
\code{k} indicates the numbers of previous observations which should be considered when observing independent exceedence for a certain barrier level. \cr \cr
\code{neta} is adjustable for special usage when the function is to slow, or extra accuracy is needed, but for the average user 500 should be fine.\cr \cr
if any of the description wasn't satisfactory see references.
}
\value{
Shows a visual plot, but the function has no quantitative returns.
}
\references{
For more details about \code{kplot}, see \code{\link{ACERm}} or \code{\link{ACER}}.
}
\author{
Kristoffer Kofoed Roedvei <rodvei@stud.ntnu.no>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ACER-package}}\cr
\code{\link{ACER}} \cr
\code{\link{ACERm}}\cr
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
kplot(X)
#
# 2. example
# A better way of making the data of Haltenbanken smaller 
# whithout loss of extremes, is by useing peak():
X<-Haltenbanken
Y<-peak(X)
kplot(Y)
#
# 3. example
# We can also specify the k's which should be plotted:
kplot(Y,k=c(1,5,10,25,50))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
