#' Normalize the single-cell count matrix
#'
#' @param count raw count matrix
#' @param meta meta information
#'
#' @return normalized count matrix
#' @import DESeq2 scran
#' @export
#'
#' @examples
normalize <- function(count, meta){
  data = DESeqDataSetFromMatrix(count, colData = meta, design = ~1)
  data = computeSumFactors(data)
  size_factor = as.vector(sizeFactors(data))
  count.norm = t(t(count)/size_factor)
  return(count.norm)
}
