\name{bs1}
\alias{bs1}
\docType{data}
\title{
BamViews instance construction related to yeast RNA-seq
}
\description{
BamViews instance construction related to yeast RNA-seq
}
%\usage{data(bs1)}
\format{
  The format is:
Formal class 'BamViews' [package "Rsamtools"] with 5 slots \cr
  ..@ bamPaths     : chr [1:8] "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/isowt5_13e.bam" "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/isowt6_13e.bam" "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/rlp5_13e.bam" "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/rlp6_13e.bam" ... \cr
  ..@ bamIndicies  : chr [1:8] "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/isowt5_13e.bam" "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/isowt6_13e.bam" "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/rlp5_13e.bam" "/Users/stvjc/ExternalSoft/R-devel/library/leeBamViews/bam/rlp6_13e.bam" ... \cr
  ..@ bamSamples   :Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. ..@ rownames       : chr [1:8] "isowt.5" "isowt.6" "rlp.5" "rlp.6" ... \cr
  .. .. ..@ nrows          : int 8 \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ elementType    : chr "ANY \cr"
  .. .. ..@ metadata       : list() \cr
  .. .. ..@ listData       :List of 2 \cr
  .. .. .. ..$ geno: chr [1:8] "isowt" "isowt" "rlp" "rlp" ... \cr
  .. .. .. ..$ lane: chr [1:8] "5" "6" "5" "6" ... \cr
  ..@ bamRanges    :Formal class 'GRanges' [package "GenomicRanges"] with 7 slots \cr
  .. .. ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 5 slots \cr
  .. .. .. .. ..@ values         : Factor w/ 1 level "Scchr13": 1 \cr
  .. .. .. .. ..@ lengths        : int 27 \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "ANY" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. .. .. ..@ start          : int [1:27] 798517 801771 804455 808999 810465 811088 818826 820255 822762 832338 ... \cr
  .. .. .. .. ..@ width          : int [1:27] 2862 933 636 234 114 108 1122 2199 1869 915 ... \cr
  .. .. .. .. ..@ NAMES          : NULL \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ strand         :Formal class 'Rle' [package "IRanges"] with 5 slots \cr
  .. .. .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 \cr
  .. .. .. .. ..@ lengths        : int 27 \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "ANY" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ seqlengths     : Named int NA \cr
  .. .. .. ..- attr(*, "names")= chr "Scchr13" \cr
  .. .. ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. .. .. ..@ rownames       : NULL \cr
  .. .. .. .. ..@ nrows          : int 27 \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "ANY" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. ..@ listData       :List of 1 \cr
  .. .. .. .. .. ..$ name: chr [1:27] "YMR266W" "YMR267W" "YMR269W" "YMRWdelta20" ... \cr
  .. .. ..@ elementType    : chr "ANY" \cr
  .. .. ..@ metadata       : list() \cr
  ..@ bamExperiment:List of 1 \cr
  .. ..$ annotation: chr "org.Sc.sgd.db" \cr
}
\details{
Illumina short reads from a very small segment of yeast chr XIII have
been collected
}
\source{
FASTQ data available at \url{ftp://ftp.ncbi.nlm.nih.gov/sra/Studies/SRP000/SRP000632/}
}
\references{
Albert Lee and Kasper Daniel Hansen and James Bullard and Sandrine Dudoit and Gavin Sherlock,
Novel Low Abundance and Transient RNAs in Yeast Revealed by Tiling Microarrays and Ultra High--Throughput Sequencing Are Not Conserved Across Closely Related Yeast Species,
PLoS Genet,
v4,
e1000299,
Dec 2008
}
\examples{
library(leeBamViews)  # bam files stored in package
bpaths = dir(system.file("bam", package="leeBamViews"), full=TRUE, patt="bam$")
#
# extract genotype and lane information from filenames
#
gt = gsub(".*/", "", bpaths)
gt = gsub("_.*", "", gt)
lane = gsub(".*(.)$", "\\\\1", gt)
geno = gsub(".$", "", gt)
#
# format the sample-level information appropriately
#
pd = DataFrame(geno=geno, lane=lane, row.names=paste(geno,lane,sep="."))
prd = new("DataFrame")  # protocol data could go here
#
# create the views object, adding some arbitrary experiment-level information
#
bs1 = BamViews(bamPaths=bpaths, bamSamples=pd,
        bamExperiment=list(annotation="org.Sc.sgd.db"))
bs1
# add ranges and tabulate reads

START=c(861250, 863000)
END=c(862750, 864000)
exc = GRanges(IRanges(start=START, end=END), seqnames="Scchr13", strand="+")
values(exc)$name = c("intv1", "intv2")  # necessary
bamRanges(bs1) = exc
bs1
tabulateReads(bs1, "+")
}
\keyword{datasets}
