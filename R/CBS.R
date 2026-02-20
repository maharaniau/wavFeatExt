#' @export
CBS <- function(obj, chr = rep(1L, length(obj))) {
  CNA.obj <- CNA(obj, chrom = chr, maploc = seq_along(obj), data.type = "binary")
  seg.obj <- segment(CNA.obj)
  num.seg <- length(seg.obj$segRows$startRow)
  cbs.seg <- numeric(length(obj))
  for (k in seq_len(num.seg)) {
    cbs.seg[seg.obj$segRows$startRow[k]:seg.obj$segRows$endRow[k]] <- seg.obj$output$seg.mean[k]
  }
  cbs.seg
}
