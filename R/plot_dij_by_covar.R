#' Make a scatterplot to contrast the covariate difference and subject-level pairwise difference
#'
#' @param count count matrix.
#' @param meta a data frame of meta information.
#' @param res_table result table of ULV analysis.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param cond_name a character for condition name in \code{meta}.
#' @param subject_name a character for subject name in \code{meta}.
#' @param ctrl_cond a character for control name.
#' @param case_cond a character for case name.
#' @param feature a character for the feature to plot.
#' @param covariate_name a character for covariate name in \code{meta}.
#' @param xlab text label for x-axis.
#' @param ylab text label for y-axis.
#' @param title title of the plot.
#'
#' @return A \code{\link{ggplot2} object.
#' @import ggplot2
#' @export
#'
#' @examples
plot_dij_by_covar <- function(count, meta, res_table = NULL,
                              normalize_option = 'pooling',
                              cond_name, subject_name,
                              ctrl_cond, case_cond,
                              feature, covariate_name,
                              xlab = paste0(covariate_name, '_diff'),
                              ylab='d_ij', title = feature){
  #-------------------------------------------------
  # preprocessing
  #-------------------------------------------------

  if(normalize_option %in% c('pooling', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize(count, meta, option = normalize_option)
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize(count, meta, option = normalize_option)
  }else {
    stop('The normalizatio option argument must be one of the following options: pooling, clr, or none.')
  }
  gene_names = rownames(count)
  ngene = nrow(count)

  #-------------------------------------------------
  # extract condition and subject information
  #-------------------------------------------------

  ## only select the cells within the case or control condition group
  ind = which(meta[, cond_name] %in% c(case_cond, ctrl_cond))
  count = count[,ind]
  meta = meta[ind,]

  condition = meta[, cond_name]
  subject = meta[, subject_name]

  table = table(subject, condition)

  ## if there is no sample in either group, return 0
  if(!(ctrl_cond %in% colnames(table)) | !(case_cond %in% colnames(table))){
    message('At least in one condition we do not have any subject!')
    return(0)
  }

  ctrl_subjs = rownames(table)[table[,ctrl_cond]!=0]
  case_subjs = rownames(table)[table[,case_cond]!=0]
  n0 = length(ctrl_subjs); n1 = length(case_subjs)
  ## if there is only 1 sample in either group, return 0
  if(n0==1 | n1==1){
    message('At least in one condition we only have one subject!')
    return(0)
  }
  ind_fct = factor(subject, levels = c(ctrl_subjs, case_subjs))

  #--------------------------------------------------
  # check covariate information
  #--------------------------------------------------

  if (length(covariate_name)>1){
    stop('Only one covariate name can be provided.')
  }
  if (!covariate_name %in% colnames(meta)){
    stop('The covariate name is not in meta information.')
  }

  g = which(rownames(count)==feature)
  y = as.numeric(count[g,])
  y.split = split(y, ind_fct)

  #---------------------------------------------------------
  # comparison between case and control
  # use rank-base method to calculate probabilistic index
  #---------------------------------------------------------

  d_mat=matrix(0, n1, n0)
  for (i in (n0+1):(n0+n1)) { # case
    for (j in 1:n0) { # ctrl
      d_mat[i-n0,j]=wilcox.test(unlist(y.split[i]), unlist(y.split[j]))$statistic/
        length(unlist(y.split[i]))/length(unlist(y.split[j]))
    }
  }
  rownames(d_mat) = names(y.split)[(n0+1):(n0+n1)]
  colnames(d_mat) = names(y.split)[1:n0]
  d = matrix(t(d_mat), n1*n0, 1)

  id_case = rep((n0+1):(n0+n1), each=n0) # case
  id_ctrl = rep(1:n0, n1) # control
  d.latent = data.frame(d, id_case, id_ctrl)

  #-----------------------------------------------------
  # design matrix for covariate adjustment
  #-----------------------------------------------------

  if (class(meta[,covariate_name]) %in% c('numeric', 'integer')){
    diff = covarDiff(meta=meta, subject_name=subject_name,
                     covariate_name=covariate_name,
                     ctrl_subjs = ctrl_subjs, case_subjs = case_subjs,
                     id_ctrl = d.latent$id_ctrl, id_case = d.latent$id_case, scale=1)
    d.latent[,paste0(covariate_name,'_diff')] = diff
  } else if (class(meta[,covariate_name])=='character'){
    diff = covarDiffDummy(meta=meta, subject_name=subject_name,
                          covariate_name=covariate_name, ref_name='female',
                          ctrl_subjs=ctrl_subjs, case_subjs=case_subjs,
                          id_ctrl=d.latent$id_ctrl, id_case=d.latent$id_case)
    d.latent = cbind(d.latent, diff)
  } else {
    stop('class of covariates should be either numeric/integer(for numeric covariate)
           or character(for categorical covariate).')
  }
  #---------------------------------
  # create scatter plot
  #---------------------------------
  df = d.latent[,c(1,4)]
  colnames(df) = c('d', 'covar_diff')
  fig <- df %>%
    ggplot(aes(x = covar_diff, y = d)) +
    geom_point() +
    theme_bw() +
    geom_smooth(method = lm)

  if(!is.null(res_table)){
    covar_ind = which(str_detect(colnames(res_table), 'diff') & !str_detect(colnames(res_table), 'pval'))
    est_pi = res_table[feature, 'PI']
    est_coeff = res_table[feature, covar_ind]
    coeff_name = colnames(res_table)[covar_ind]
    fig <- fig +
      labs(title = title, x = xlab, y = ylab,
           caption = paste0('PI: ', signif(est_pi, 3), '\n',
                            coeff_name, ': ', signif(est_coeff, 3)))
  }else{
    fig <- fig +
      labs(title = title, x = xlab, y = ylab)
  }

  return(fig)
}
