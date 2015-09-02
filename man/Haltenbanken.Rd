\name{Haltenbanken}
\alias{Haltenbanken}
\docType{data}
\title{ Haltenbanken_Hs Data Set}
\description{A practical data set.}
\usage{Haltenbanken}
\format{
157784 observation over 18 realizations. One measurement each hour for 18 years.
}
\details{
One observations each hour for 18 years. There are some \code{NA} values from leap year.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(Haltenbanken)
X<-Haltenbanken
# X contain the data set, 
# and we can now for example plot the histogram
hist(X)
}
\keyword{datasets}
