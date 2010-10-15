\name{leeRPKM}
\alias{leeRPKM}
\docType{data}
\title{
supplemental data extract on RNA seq results in yeast
}
\description{
supplemental data extract on RNA seq results in yeast
}
\usage{data(leeRPKM)}
\format{
  A data frame with 6291 observations on the following 16 variables.
  \describe{
    \item{\code{chr}}{a numeric vector}
    \item{\code{strand}}{a numeric vector}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{name}}{a factor with levels \code{LSR1} \code{NME1} \code{YAL001C} \code{YAL002W} \code{YAL003W} ...}
    \item{\code{feature}}{a factor with levels \code{CDS} \code{CDS_unchar} \code{snRNA} \code{snoRNA}}
    \item{\code{orf_classification}}{a factor with levels \code{Uncharacterized} \code{Verified} \code{silenced_gene3AVerified}}
    \item{\code{gene}}{a factor with levels \code{AAC1} \code{AAC3} \code{AAD10} \code{AAD14} \code{AAD15} \code{AAD16} \code{AAD3} \code{AAD4} ...}
    \item{\code{wt.reads}}{a numeric vector}
    \item{\code{rrp.reads}}{a numeric vector}
    \item{\code{ski.reads}}{a numeric vector}
    \item{\code{xrn.reads}}{a numeric vector}
    \item{\code{wt.rpkm}}{a numeric vector}
    \item{\code{rrp.rpkm}}{a numeric vector}
    \item{\code{ski.rpkm}}{a numeric vector}
    \item{\code{xrn.rpkm}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
imported from supplemental data
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Lee et al PLOS genetics December 2008 ; Volume 4 ; Issue 12 ; e1000299
}
\examples{
data(leeRPKM)
leeRPKM[1:5,]
}
\keyword{datasets}
