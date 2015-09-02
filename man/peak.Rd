\name{peak}
\alias{peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting Peak Values from Data
}
\description{
Peak is used to extract values from time series which is higher than it's neighbor. It is also capable of extracting peaks from time series with multiple realizations.
}
\usage{
peak(X, n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ 
a data set in form of an array or matrix respectively if it has single or multiple realizations. The realizations of the set should be by rows.}
  \item{n}{
an integer for the distance which define a points neighbor (see below).
}
}
\details{
The neighborhood of a point A is defined as all the point which is less then \code{n} points away from A. If A is larger then its neighbors it will be extracted as a peak.
}
\value{
The peaks. Returns an array or matrix respectively if \code{X} is an array or matrix. The numbers of rows equal \code{X} but the numbers of column is smaller due to the removed elements which are not local maximas.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## With multiple realizatiosn
data(Haltenbanken)
X<-peak(Haltenbanken,3)
## X will only contain data which is larger than all it's 3 steps away neighbor
## End example}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line