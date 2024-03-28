#' Wilcoxon test at subject level
#'
#' @param count count matrix.
#' @param meta metadata.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param subj_name a character for subject name in \code{meta}.
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
                           subj_name, cond_name,
                           case_cond, ctrl_cond,
                           aggregate_method = 1){
  #---------------------------------
  # pre-processing
  #---------------------------------
  if(normalize_option %in% c('pooling', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize(count, meta, option = normalize_option)
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize(count, meta, option = normalize_option)
  }else {
    stop('The normalizatio option argument must be one of the following options: pooling, clr, or none.')
  }

  ngene = nrow(count)
  # only select the cells in the two groups
  count = count[, meta[,cond_name] %in% c(case_cond, ctrl_cond)]
  meta = meta[meta[,cond_name] %in% c(case_cond, ctrl_cond),]
  pval = rep(NA, ngene)
  logfc = rep(NA, ngene)

  # find case and control subjects
  case_subjs = unique(meta[meta[,cond_name]==case_cond, subj_name])
  ctrl_subjs = unique(meta[meta[,cond_name]==ctrl_cond, subj_name])

  if(length(case_subjs)<2 | length(ctrl_subjs)<2){
    stop('At least one condition has less than 2 samples!!!')
  }

  for (g in 1:ngene) {

    if(g %% 100 == 0){
      message('Wilcoxon testing: ', g)
    }
    y = as.numeric(count[g,])
    subj_info = meta[,subj_name]
    y.split = split(y, subj_info)

    if(aggregate_method == 1){
      pseudo = as.numeric(sapply(y.split, mean))
    }else if(aggregate_method == 2){
      pseudo = as.numeric(sapply(y.split, median))
    }else if(aggregate_method == 3){
      pseudo = as.numeric(sapply(y.split, function(x){mean(x>0)}))
    }else{
      stop('aggregate_method argument must be from 1 to 3. For aggregate_method=1, it means average value;
           for aggregate_method=2, it means median value;
           for aggregate_method=3, it means proportion of values greater than 0.')
    }

    case_ind = which(names(y.split) %in% case_subjs)
    cases = pseudo[case_ind]
    ctrl_ind = which(names(y.split) %in% ctrl_subjs)
    ctrls = pseudo[ctrl_ind]

    test = wilcox.test(cases, ctrls)
    pval[g] = test$p.value
    logfc[g] = log2(mean(cases+0.001)/mean(ctrls+0.001))
    if(is.na(test$p.value)){
      message(g)
      message('cases:')
      print(cases)
      message('ctrls:')
      print(ctrls)
    }
  }
  padj = p.adjust(pval, method = 'BH')
  res = data.frame(logfc=logfc, pval=pval, padj=padj)
  rownames(res) = rownames(count)
  return(res)
}
