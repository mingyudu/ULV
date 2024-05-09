#' Wilcoxon test at subject level
#'
#' @param count count matrix.
#' @param meta metadata.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param subject_name a character for subject name in \code{meta}.
#' @param cond_name a character for condition name in \code{meta}.
#' @param case_cond a character for case name.
#' @param ctrl_cond a character for control name.
#' @param aggregate_method aggregate method.
#'
#' @return A result table containing log2-fold change, p-values, and FDR-adjusted p-values.
#'
#' @export
#'
#' @examples
wilcox_subject <- function(count, meta, normalize_option = 'pooling',
                           subject_name, cond_name,
                           case_cond, ctrl_cond,
                           aggregate_method = 1){
  #---------------------------------
  # pre-processing
  #---------------------------------
  if(normalize_option %in% c('pooling', 'LogNormalize', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize_data(count, meta, option = normalize_option)
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize_data(count, meta, option = normalize_option)
  }else {
    stop('The normalizatio option argument must be one of the following options: pooling, clr, or none.')
  }

  ngene = nrow(count)
  # only select the cells in the two groups
  count = count[, meta[,cond_name] %in% c(case_cond, ctrl_cond)]
  meta = meta[meta[,cond_name] %in% c(case_cond, ctrl_cond),]
  pval = rep(NA, ngene)
  fc = rep(NA, ngene)

  # find case and control subjects
  case_subjs = unique(meta[meta[,cond_name]==case_cond, subject_name])
  ctrl_subjs = unique(meta[meta[,cond_name]==ctrl_cond, subject_name])

  if(length(case_subjs)<2 | length(ctrl_subjs)<2){
    stop('At least one condition has less than 2 samples!!!')
  }

  if(!aggregate_method %in% c(1:3)){
    stop('aggregate_method argument must be from 1 to 3. For aggregate_method=1, it means average value;
           for aggregate_method=2, it means median value;
           for aggregate_method=3, it means proportion of values greater than 0.')
  }
  # Function to perform Wilcoxon test for each gene
  perform_wilcoxon_test <- function(y) {

    y <- as.numeric(y)
    subj_info <- meta[,subject_name]
    y.split <- split(y, subj_info)

    # calculate pseudo values based on aggregate method
    pseudo <- switch(
      aggregate_method,
      as.numeric(sapply(y.split, mean)),   # Mean
      as.numeric(sapply(y.split, median)), # Median
      as.numeric(sapply(y.split, function(x) mean(x > 0))) # Proportion of values > 0
    )

    # get ind for case and control groups
    case_ind <- which(names(y.split) %in% case_subjs)
    cases <- pseudo[case_ind]
    ctrl_ind <- which(names(y.split) %in% ctrl_subjs)
    ctrls <- pseudo[ctrl_ind]

    # wilcoxon
    test <- tryCatch(
      wilcox.test(cases, ctrls),
      error = function(e) list(p.value = NA, statistic = NA)  # Handle errors during test
    )

    if(is.na(test$p.value)) {
      message(g)
      message('cases:')
      print(cases)
      message('ctrls:')
      print(ctrls)
    }

    # Return test results including fold change FC and p-value
    c(fc = test$statistic / (length(case_subjs) * length(ctrl_subjs)),
      pval = test$p.value)
  }

  # apply the function across all genes
  res <- t(apply(count, 1, perform_wilcoxon_test))

  res <- as.data.frame(res)
  colnames(res)[which(colnames(res)=='fc.W')] = 'FC'

  res$padj = p.adjust(res$pval, method = 'BH')
  rownames(res) = rownames(count)

  return(res)
}


#' Wilcoxon test at cell level
#'
#' @param count count matrix.
#' @param meta metadata.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param subject_name a character for subject name in \code{meta}.
#' @param cond_name a character for condition name in \code{meta}.
#' @param case_cond a character for case name.
#' @param ctrl_cond a character for control name.
#'
#' @return A result table containing log2-fold change, p-values, and FDR-adjusted p-values.
#'
#' @export
#'
#' @examples
wilcox_cell <- function(count, meta, normalize_option = 'pooling',
                           subject_name, cond_name,
                           case_cond, ctrl_cond){
  #---------------------------------
  # pre-processing
  #---------------------------------
  if(normalize_option %in% c('pooling', 'LogNormalize', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize_data(count, meta, option = normalize_option)
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize_data(count, meta, option = normalize_option)
  }else {
    stop('The normalizatio option argument must be one of the following options: pooling, clr, or none.')
  }

  ngene = nrow(count)
  # only select the cells in the two groups
  count = count[, meta[,cond_name] %in% c(case_cond, ctrl_cond)]
  meta = meta[meta[,cond_name] %in% c(case_cond, ctrl_cond),]

  # wilcoxon at single cell level
  res = apply(count, 1, function(x){
    test = wilcox.test(x[meta[,cond_name]==case_cond], x[meta[,cond_name]==ctrl_cond])
    fc = test$statistic/(sum(meta[,cond_name]==case_cond)*sum(meta[,cond_name]==ctrl_cond))
    pval = test$p.value
    return(c(fc = fc, pval = pval))
  })
  res = as.data.frame(t(res))
  colnames(res) = c('FC', 'pval')

  res$padj = p.adjust(res$pval, method = 'BH')
  rownames(res) = rownames(count)
  return(res)
}
