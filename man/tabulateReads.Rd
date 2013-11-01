\name{tabulateReads}
\alias{tabulateReads}
\alias{tabulateReads,BamViews,characterORNULL,logical,function-method}
\alias{tabulateReads,BamViews,characterORNULL,missing,missing-method}
\alias{tabulateReads,BamViews,missing,missing,missing-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tabulate counts of alignments occurring in specified genomic regions
}
\description{
tabulate counts of alignments occurring in specified genomic regions
}
\usage{
tabulateReads(bv, strandmarker=NULL, as.GRanges=FALSE, applier=lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bv}{
\code{\link[Rsamtools]{BamViews-class}} instance
}
  \item{strandmarker}{
character atom: `+' or `-'; if missing, ignore strand
}
  \item{as.GRanges}{
logical directive to return a GRanges instance instead of a matrix
  }
 \item{applier}{lapply-like function; if unspecified and multicore is
attached will use mclapply}
}
\details{
\code{\link[Rsamtools]{readGAlignmentsFromBam}} is the basic engine
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
myrn = GRanges(IRanges(start=seq(861250, 862750, 100), width=100),
    seqnames="Scchr13", strand="+")

values(myrn)$name = paste("til", 1:length(myrn), sep=".")
bamRanges(bs1) = myrn
tabulateReads(bs1, "+")
#
# a related computation based on countBam
lapply(bamPaths(bs1)[1:2], function(x) 
    countBam(x,  param=ScanBamParam(which=bamRanges(bs1))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
