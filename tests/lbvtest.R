###################################################
### chunk number 1: lkd
###################################################
library(leeBamViews)  # bam files stored in package
bpaths = dir(system.file("bam", package="leeBamViews"), full=TRUE, patt="bam$")
#
# extract genotype and lane information from filenames
#
gt = gsub(".*/", "", bpaths)
gt = gsub("_.*", "", gt)
lane = gsub(".*(.)$", "\\1", gt)
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
#
# get some sample-level data
#
bamSamples(bs1)$geno


###################################################
### chunk number 2: lkc
###################################################
START=c(861250, 863000)
END=c(862750, 864000)
exc = GRanges(IRanges(start=START, end=END), seqnames="Scchr13", strand="+")
bamRanges(bs1) = exc

cov2baseTrack = function(rle, start, end,
   dp = DisplayPars(type="l", lwd=0.5, color="black"),
   countTx=function(x)log10(x+1)) {
 require(GenomeGraphs)
 if (!is(rle, "Rle")) stop("requires instance of Rle")
 dat = rle@values
 loc = cumsum(rle@lengths)
 ok = which(loc >= start & loc <= end)
 makeBaseTrack(base = loc[ok], value=countTx(dat[ok]),
    dp=dp)
}

library(biomaRt)
mart = useMart("ensembl", "scerevisiae_gene_ensembl")


###################################################
plotStrains = function(bs, query, start, end, snames, mart, martchr, seqname, strand="+") {
 mm = as.matrix(findOverlaps(bamRanges(bs), query))
 if (nrow(mm) < 1) stop("no overlap between query and input bamViews")
 filtbs = bs[mm[,"subjectHits"], ]
 cov = lapply(bamPaths(filtbs), function(x)coverage(readGappedAlignments(x))[[seqname]])
 covtrs = lapply(cov, function(x) cov2baseTrack(x, start, end,
   countTx = function(x) pmin(x,80)))
 names(covtrs) = snames
 gr = makeGeneRegion(start, end, chromosome=martchr,
       strand=strand, biomart=mart, dp=DisplayPars(plotId=TRUE,
       idRotation=0, idColor="black"))
 grm = makeGeneRegion(start, end, chromosome=martchr,
       strand="-", biomart=mart, dp=DisplayPars(plotId=TRUE,
       idRotation=0, idColor="black"))
 covtrs[[length(covtrs)+1]] = gr
 covtrs[[length(covtrs)+1]] = makeGenomeAxis()
 covtrs[[length(covtrs)+1]] = grm
 gdPlot( covtrs, minBase=start[1], maxBase=end[1] )
}


NN = GRanges(IRanges(start=START,end=END), seqnames="Scchr13")
plotStrains(bs1, NN, 800000, 900000, LETTERS[1:8], mart, "XIII", seqname="Scchr13")

