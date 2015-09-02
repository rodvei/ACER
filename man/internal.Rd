\name{internal}
\alias{internalACER}
\alias{eta1}
\alias{RealizationByRows}
\alias{RealizationLength}
\alias{findcond}
\alias{Gumbel}
\alias{reverseGumbel}
\alias{General}
\alias{reverseGeneral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Internal ACER2 functions}
\description{
Internal functions in ACER2. \code{RealizationByRows} controls that the data set is correct, meaning that realization is by rows and data by columns. If each realization has different size (example after \code{peak}) \code{RealizationLength} returns the length of each realization. \code{eta1} gives an estimate of the barrier level for where regular tail behavior begins. For MSEoptimization internal function see \code{\link{internal-MSEoptimization}}.
}
\usage{
RealizationByRows(X)
RealizationLength(X)
eta1(Akj,eta)
findcond(condition)
Gumbel(eta,est)
reverseGumbel(eps,est)
General(eta,est)
reverseGeneral(eps,est)

}

\keyword{internal} 
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ the data set in matrix formate with realizations by rows. If a realization have any \code{NA} they should be placed in the last columns.}
  \item{Akj}{ the full vector of Akj extacted from the ACER class (see also \code{\link{ACERm}}).}
  \item{eta}{ the full vector of barrier levels extacted from the ACER class.}
  \item{condition}{condition, an array where each element is logical. TRUE when condition is met, else FALSE.}
  \item{eta}{ the barrier level}
  \item{eps}{ the probablility}
  \item{est}{ coefficients for hte model}
}
\details{
\code{X} Is often the data set from \code{peak} which returns a matrix where the first elements in each realizations is numbers while the last values are empty \code{NA}'s (peaks returns difference size for each realizations, but the return is formatted as matrix and to ensure equal length empty spots are filled with \code{NA}'s, see \code{\link{peak}}). For \code{RealizationLength}, \code{X} Should be formatted with numbers first, and \code{NA} last, else the returned value will become incorrect.
}
\value{
\code{RealizationByRows} returns a logical \code{TRUE} or \code{FALSE} respectively if realizations is by row as it should or if it's not.\cr
\code{RealizationLength} returns a vector. The i'th element in the vector respond to the length of the data set in the i'th realizations. \cr
\code{eta1} returns an estimate of the barrier level for where regular tail behavior begins.\cr
\code{StartEnd} Return an array of two elements, first start posistion, second last positon of which all elements in between is TRUE.\cr
\code{Gumbel} Returns the epsilon in the Gumbel case, see Gumbel details \code{\link{MSEoptimization}}.
\code{reverseGumbel} Returns the eta for a given epsilon in the general case, see Gumbel details \code{\link{MSEoptimization}}.
\code{General} Returns the epsilon in the general case, see general details \code{\link{MSEoptimization}}.
\code{reverseGeneral} Returns the eta for a given epsilon in the general case, see general details \code{\link{MSEoptimization}}.

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

\seealso{\code{\link{peak}}, \code{\link{ACER}}, \code{\link{ACERm}}, \code{\link{kplot}}}
\examples{
#no example for private functions
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
