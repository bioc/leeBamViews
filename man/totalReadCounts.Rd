\name{totalReadCounts}
\alias{totalReadCounts}
\alias{totalReadCounts,BamViews-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
scan BAM files for total read counts
}
\description{
scan BAM files for total read counts
}
\usage{
totalReadCounts(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{\link[Rsamtools]{BamViews-class}} instance
}
}
\details{
slow procedure -- does lightweight scan of entire file
}
\value{
named integer vector of read counts per sample
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
example(bs1)
totalReadCounts(bs1)
}
\keyword{ models }
