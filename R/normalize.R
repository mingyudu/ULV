#' Normalize the single-cell count matrix
#'
#' @param count raw count matrix
#' @param meta meta information
#' @param option a character value representing the normalization method to apply to the count matrix
#'
#' @return normalized count matrix
#' @import DESeq2 scran
#' @export
#'
#' @examples
normalize <- function(count, meta, option = 'none'){
  if(option=='pooling'){
    data = DESeqDataSetFromMatrix(count, colData = meta, design = ~1)
    data = computeSumFactors(data)
    size_factor = as.vector(sizeFactors(data))
    count.norm = t(t(count)/size_factor)
  }else if(option=='clr'){
    count.norm = apply(count, 1, function(y){
      y = as.numeric(y)
      y[y==0] = 0.0001
      y = log(y) - mean(log(y))
      })
    count.norm = t(count.norm)
    colnames(count.norm) = colnames(count)
  }else if(option=='none'){
    count.norm = count
  }
  return(count.norm)
}
