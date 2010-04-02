setMethod("strand", "logical", function(x) {
if (0L == length(x)) strand()
else strand(ifelse(x, "+", "-"))
})

.tabulateReads = function(bv, strandmarker=NULL, as.GRanges=FALSE) {
 if (!is.null(strandmarker) && !(strandmarker %in% c("+", "-")))
   stop("if non-missing, strandmarker must be either NULL, +, or -")
 br = bamRanges(bv)
 rnames = elementMetadata(br)$name 
 if (length(rnames)==0) stop("bamRanges(bv) must have a name element for table annotation")
 nregions = length(br)
# you need to use the file interface for now (March 25 2010)
 alignByFirstRange = function(bv) lapply(bamPaths(bv), function(x) readBamGappedAlignments(x, which=bamRanges(bv[1,])))
 als = lapply(1:nregions, function(i)try(alignByFirstRange(bv[i,]), silent=TRUE))
 ok = !sapply(als, inherits, "try-error")
 if (any(!ok)) warning(paste("readGappedAlignments failed for range(s)",
      paste(rnames[-which(ok)], collapse=" "), ", so dropping these"))
 als = als[which(ok)]
# strs = lapply(als, sapply, strand)
 strs = lapply(als, function(x)sapply(as(x, "list"), strand))
 if (is.null(strandmarker)) ans = sapply(strs, sapply, length)
 else ans = sapply(strs, sapply, function(x)sum(x==strandmarker))
 sta = start(br)[which(ok)]
 end = end(br)[which(ok)]
 ans = rbind(start=sta, end=end, ans)
 colnames(ans) = rnames[which(ok)]
 rownames(ans) = c("start", "end", rownames(bamSamples(bv)))
 dat = ans  # temporary hold for GRanges discipline introduction
 if (!as.GRanges) return(ans)
 ir = IRanges(sta, end)
 ans = GRanges(ir, seqnames=seqnames(br), strand=strandmarker)
 values(ans)$name = rnames[which(ok)]
 dat = DataFrame(t(dat)[,-c(1,2)])
 values(ans) = DataFrame(values(ans), dat)
 ans
}

setGeneric("tabulateReads", function(bv, strandmarker, as.GRanges)
 standardGeneric("tabulateReads"))

setMethod("tabulateReads", c("BamViews", "characterORNULL", "logical"), 
  function(bv, strandmarker=NULL, as.GRanges) {
  .tabulateReads(bv, strandmarker, as.GRanges)
})
setMethod("tabulateReads", c("BamViews", "characterORNULL", "missing"), 
  function(bv, strandmarker=NULL, as.GRanges) {
  .tabulateReads(bv, strandmarker, FALSE)
})

setMethod("tabulateReads", c("BamViews", "missing", "missing"), 
  function(bv, strandmarker=NULL) {
  .tabulateReads(bv, NULL, FALSE)
})
