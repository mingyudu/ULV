# create pseudobulk
#' Create pseudobulk count matrix for DESeq2 method
#'
#' @param count count matrix.
#' @param meta meta data.
#' @param subject_name a character for the subject name in \code{meta}.
#' @param cond_name a character for the condition name in \code{meta}.
#' @param ctrl_cond a character for the control condition name.
#' @param case_cond a character for the case condition name.
#' @param aggregate a character for the aggregation method, either 'sum' or 'mean'.
#' @param covariate_name_list a vector of covariate name to adjust.
#'  Set to NULL if there is no covariate adjustment.
#'
#' @return A list containing the pseudobulk count matrix and metadata
#'  at pseudobulk level.
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix rowMeans
#' @export
#'
#' @examples
create_pseudobulk = function(count, meta, subject_name, cond_name,
                             ctrl_cond, case_cond,
                             aggregate, covariate_name_list=NULL){

  stopifnot(subject_name %in% colnames(meta),
            cond_name %in% colnames(meta))
  stopifnot(ctrl_cond %in% meta[, cond_name],
            case_cond %in% meta[, cond_name])
  if(!is.null(covariate_name_list)){
    stopifnot(all(covariate_name_list %in% colnames(meta)))
  }

  # only select two conditions
  ind = which(meta[, cond_name] %in% c(ctrl_cond, case_cond))
  count = count[, ind]
  meta = meta[ind, ]

  # extract subject IDs
  tab = table(meta[,subject_name], meta[,cond_name])
  ctrl_ids = rownames(tab)[tab[,ctrl_cond]!=0]
  case_ids = rownames(tab)[tab[,case_cond]!=0]
  message(paste('ctrl group:', colnames(tab)[1], length(ctrl_ids)))
  message(paste('case group:', colnames(tab)[2], length(case_ids)))
  all_ids = c(ctrl_ids, case_ids)

  # create pseudobulk matrix
  count.pseudo = matrix(NA, nrow=dim(count)[1], ncol=length(all_ids))
  cell_num = rep(NA, length(all_ids))
  for (j in 1:length(all_ids)) {
    id = all_ids[j]
    if(aggregate=='sum'){
      count.pseudo[,j] = Matrix::rowSums(count[, meta[, subject_name]==id])
    }else if (aggregate=='mean'){
      count.pseudo[,j] = Matrix::rowMeans(count[, meta[, subject_name]==id])
    }else{
      stop('aggregate can only be either sum or mean.')
    }
  }
  rownames(count.pseudo) = rownames(count)
  colnames(count.pseudo) = all_ids

  # create pseudobulk metadata
  meta.pseudo = unique(meta[,c(subject_name, cond_name, covariate_name_list)])
  rownames(meta.pseudo) = meta.pseudo[, subject_name]
  meta.pseudo = meta.pseudo[all_ids,]

  return(list(count.pseudo = count.pseudo, meta.pseudo = meta.pseudo))
}


#' DESeq2 at pseudobulk level
#'
#' @param count count matrix.
#' @param meta meta data.
#' @param subject_name a character for the subject name in \code{meta}.
#' @param cond_name a character for the codition name in \code{meta}.
#' @param ctrl_cond a character for the control condition name.
#' @param case_cond a character for the case condition name.
#' @param numerical_covar a vector of numerical covariate name to adjust.
#'  Set to NULL if there is no covariate adjustment.
#' @param categorical_covar a vector of categorical covariate name to adjust.
#'  Set to NULL if there is no covariate adjustment.
#'
#' @return A result table for DESeq2 result.
#' @import DESeq2
#' @export
#'
#' @examples
run_DESeq_pseudobulk = function(count, meta, subject_name, cond_name,
                     ctrl_cond, case_cond,
                     numerical_covar = NULL, categorical_covar = NULL){
  # use raw count!
  data.pseudo = create_pseudobulk(count, meta, subject_name, cond_name,
                                  ctrl_cond, case_cond, aggregate='sum',
                                  c(numerical_covar, categorical_covar))
  count.pseudo = data.pseudo$count.pseudo
  meta.pseudo = data.pseudo$meta.pseudo
  meta.pseudo[, subject_name] = as.factor(meta.pseudo[, subject_name])
  meta.pseudo[, cond_name] = factor(meta.pseudo[, cond_name], levels = c(ctrl_cond, case_cond))

  for (vr in categorical_covar){
    meta.pseudo[,vr] = factor(meta.pseudo[,vr])
  }
  for (vr in numerical_covar) {
    meta.pseudo[,vr] = scale(meta.pseudo[,vr])
  }

  # DESeq2
  if(is.null(c(numerical_covar, categorical_covar))){
    dds = DESeqDataSetFromMatrix(count.pseudo, meta.pseudo,
                               design = as.formula(paste0('~ ', cond_name)))
  }else{
    dds = DESeqDataSetFromMatrix(count.pseudo, meta.pseudo,
                                 design = as.formula(paste0('~ ',
                                                            paste0(c(numerical_covar, categorical_covar), sep = ' + ', collapse = ''),
                                                            cond_name)))
  }

  dds = DESeq(dds, minReplicatesForReplace = Inf)
  res = results(dds, contrast = c(cond_name, case_cond, ctrl_cond),
                cooksCutoff = FALSE, independentFiltering = FALSE)
  return(res)
}


#' DESeq2 at single cell level
#'
#' @param count count matrix.
#' @param meta meta data.
#' @param subject_name a character for subject name in \code{meta}.
#' @param cond_name a character for the codition name in \code{meta}.
#' @param ctrl_cond a character for the control condition name.
#' @param case_cond a character for the case condition name.
#' @param normalize_option either 'pooling' or 'RC' for the normalization method to apply on the single-cell data.
#'  Set to 'pooling' as default.
#' @param numerical_covar a vector of numerical covariate name to adjust.
#'  Set to NULL if there is no covariate adjustment.
#' @param categorical_covar a vector of categorical covariate name to adjust.
#'  Set to NULL if there is no covariate adjustment.
#'
#' @return A result table for DESeq2 result.
#' @import DESeq2 scran
#' @export
#'
#' @examples
run_DESeq_sc = function(count, meta, subject_name, cond_name,
                        ctrl_cond, case_cond, normalize_option = 'pooling',
                        numerical_covar = NULL, categorical_covar = NULL){

  meta[, cond_name] = factor(meta[, cond_name], levels = c(ctrl_cond, case_cond))

  for (vr in categorical_covar){
    meta[,vr] = factor(meta[,vr])
  }
  for (vr in numerical_covar) {
    meta[,vr] = scale(meta[,vr])
  }

  # obtain the size factor on a subject-by-subject basis
  if(normalize_option == 'pooling'){
    df = data.frame(cell_name = colnames(count),
                    subject = meta[, subject_name])
    rownames(df) = df$cell_name

    df$size_factor = NA
    subject_list = unique(meta[, subject_name])
    for (i in 1:length(subject_list)) {
      subj = subject_list[i]
      dds = DESeqDataSetFromMatrix(count, colData = meta, design = ~1)
      dds = dds[,colData(dds)[,subject_name]==subj]
      # use scran normalization suggested by DESeq2 tutorial
      # https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
      dds = computeSumFactors(dds)
      size_factor = sizeFactors(dds)
      size_factor = as.vector(size_factor)
      df$size_factor[df$subject==subj] = size_factor
    }
    size_factor = df$size_factor
  }else if(normalize_option == 'RC'){
    rds = colSums(count)
    med_rds = median(rds)
    size_factor = rds/med_rds
  }

  # DESeq2
  if(is.null(c(numerical_covar, categorical_covar))){
    dds = DESeqDataSetFromMatrix(count, meta,
                                 design = as.formula(paste0('~ ', cond_name)))
  }else{
    dds = DESeqDataSetFromMatrix(count, meta,
                                 design = as.formula(paste0('~ ',
                                                            paste0(c(numerical_covar, categorical_covar), sep = ' + ', collapse = ''),
                                                            cond_name)))
  }

  sizeFactors(dds) = size_factor

  # run DESeq in single-cell version
  dds = DESeq(dds, useT=TRUE, minmu=1e-6, minReplicatesForReplace = Inf)
  res = results(dds, contrast = c(cond_name, case_cond, ctrl_cond),
                cooksCutoff = FALSE, independentFiltering = FALSE)
  return(res)
}
