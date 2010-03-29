\name{tabulateReads}
\alias{tabulateReads}
\alias{tabulateReads,BamViews,character-method}
\alias{tabulateReads,BamViews,missing-method}
\alias{tabulateReads,BamViews,characterORNULL-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tabulate counts of alignments occurring in specified genomic regions
}
\description{
tabulate counts of alignments occurring in specified genomic regions
}
\usage{
tabulateReads(bv, strandmarker)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bv}{
\code{\link[Rsamtools]{BamViews-class}} instance
}
  \item{strandmarker}{
character atom: `+' or `-'; if missing, ignore strand
}
}
\details{
\code{\link[Rsamtools]{readBamGappedAlignments}} is the basic engine
for this task
}
\value{
annotated matrix with start, end, and samples as rows, regions
as columns, and read counts as cell entries
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
example(bs1)
#
# counts in a partition
#
myrn = RangedData(IRanges(start=seq(861250, 862750, 100), width=100),
    space="Scchr13")
myrn$name = paste("til", 1:nrow(myrn), sep=".")
bamRanges(bs1) = myrn
tabulateReads(bs1, "+")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
