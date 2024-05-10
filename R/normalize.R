#' Normalize the single-cell count matrix
#'
#' @param count raw count matrix.
#' @param meta meta information.
#' @param subject_name a character for subject name in \code{meta}.
#' @param option a character value representing the normalization method to apply to the count matrix.
#'
#' @return normalized count matrix
#' @import DESeq2 scran
#' @export
#'
#' @examples
normalize_data <- function(count, meta, subject_name, option = 'none'){
  if(option=='pooling'){
    df = data.frame(cell_name = colnames(count),
                    subject = meta[, subject_name])
    rownames(df) = df$cell_name

    df$size_factor = NA
    subject_list = unique(meta[, subject_name])
    for (i in 1:length(subject_list)) {
      subj = subject_list[i]
      # print(subj)
      dds = DESeqDataSetFromMatrix(count, colData = meta, design = ~1)
      dds = dds[,colData(dds)[,subject_name]==subj]
      dds = computeSumFactors(dds)
      size_factor = sizeFactors(dds)
      # cell_name = names(size_factor)
      # print(all(cell_name == df$cell_name[df$subject==subj]))
      size_factor = as.vector(size_factor)
      df$size_factor[df$subject==subj] = size_factor
    }
    size_factor = df$size_factor
    count.norm = t(t(count)/size_factor)
  }else if(option=='LogNormalize'){
    # LogNormalize from Seurat
    rds = colSums(count)
    med_rds = median(rds)
    count.norm = log1p(t(t(count)/rds)*med_rds)
  }else if(option=='RC'){
    # LogNormalize from Seurat
    rds = colSums(count)
    med_rds = median(rds)
    count.norm = t(t(count)/rds)*med_rds
  }else if(option=='clr'){
    # clr function from Seurat
    clr_function <- function(x) {
      return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
    }
    count.norm = apply(count, 2, clr_function)
  }else if(option=='none'){
    count.norm = count
  }
  return(count.norm)
}
