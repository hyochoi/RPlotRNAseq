#' find_region(gpos=c(10614056,10596790),Ranges=geneRanges)
#' @import SCISSOR
#' @export
find_region = function(gpos,Ranges,ranges.outside=T) {
  if (missing(gpos)) {
    stop("gpos (genomic position) should be specified.")
  }
  if (missing(Ranges)) {
    stop("Ranges is missing. See get_Ranges.")
  }
  sapply(gpos,function(x) {find_region_single(gpos=x,Ranges=Ranges,ranges.outside=ranges.outside)})
}

find_region_single = function(gpos,Ranges,ranges.outside=T) {
  # find region in the locus provided by Ranges
  if (missing(gpos)) {
    stop("gpos (genomic position) should be specified.")
  }
  if (missing(Ranges)) {
    stop("Ranges is missing. See get_Ranges.")
  }
  x = as.numeric(gpos)
  ie = which(((x - Ranges$gRanges[,1])*(x - Ranges$gRanges[,4])) <= 0)  # which exon
  if (length(ie)==1) {
    ip = which(c(Ranges$gRanges[ie,1]:Ranges$gRanges[ie,4])==x)  # which position
    return(c(Ranges$lRanges[ie,1]:Ranges$lRanges[ie,4])[ip])
  } else {
    if (ranges.outside) {
      nexon = dim(Ranges$gRanges)[1]
      ib = which(((x - Ranges$gRanges[1:(nexon-1),1])*(x - Ranges$gRanges[2:nexon,4])) <= 0)  # between which parts
      if (length(ib)==1) {
        Ranges$lRanges[ib,4]
      } else {
        if (which.min(abs(x - Ranges$gRanges))==1) {
          Ranges$lRanges[1,1] - abs(Ranges$gRanges[1,1]-x)
        } else {
          Ranges$lRanges[nexon,4] + abs(Ranges$gRanges[nexon,4]-x)
        }
      }
    } else {
      warning("The provided position is outside the regions supported by Ranges. If needed, use ranges.outside=T. ")
      return("NA")
    }
  }
}
