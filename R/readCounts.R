
setGeneric("totalReadCounts", function(x)standardGeneric("totalReadCounts"))

setMethod("totalReadCounts", "BamViews", function(x) {
 tmp = lapply( bamPaths(x), scanBam, param=ScanBamParam(what="pos"))
 ans = sapply(tmp, function(x) length(x[[1]][[1]]))
 names(ans) = sapply( bamPaths(x), basename )
 ans
})
