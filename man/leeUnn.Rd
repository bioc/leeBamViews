\name{leeUnn}
\alias{leeUnn}
\docType{data}
\title{
supplemental data extracts on existing evidence of transcription in yeast
}
\description{
supplemental data extracts on existing evidence of transcription in yeast
}
\usage{data(leeUnn)}
\format{
  A data frame with 54822 observations on the following 11 variables.
  \describe{
    \item{\code{chr}}{a numeric vector}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
    \item{\code{strand}}{a numeric vector}
    \item{\code{lengthWithoutMask}}{a numeric vector}
    \item{\code{length}}{a numeric vector}
    \item{\code{lambda}}{a numeric vector}
    \item{\code{background5}}{a logical vector}
    \item{\code{background20}}{a logical vector}
    \item{\code{reads}}{a numeric vector}
    \item{\code{study}}{a factor with levels \code{David} \code{Davis} \code{Miura} \code{Nagalakshmi}}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
from Lee et al PLoS genetics December 2008  Volume 4  Issue 12  e1000299
supplemental data

information on unnannotated transcripts for which some evidence of
transcription was obtained in this experiment
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(leeUnn)
leeUnn[1:5,]
}
\keyword{datasets}
