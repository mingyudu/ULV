#' Fit ULV model
#'
#' \code{fit_ULV} fits a U-statistic based Latant Variable model that takes
#' advantage of the robustness of rank-based methods and the statistical
#' efficiency of parametric methods for small sample sizes.
#'
#' @param count count matrix.
#' @param meta a data frame of meta information.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param subject_name a character for subject name in \code{meta}.
#' @param cond_name a character for condition name in \code{meta}.
#' @param ctrl_cond a character for control name.
#' @param case_cond a character for case name.
#' @param weighted a binary variable indicating whether the analysis is
#'  weighted by varying cluster sizes.
#' @param covariate_name_list a vector of character of covariate names.
#'  If not, set it to NULL. The names must be included in the columns
#'  of \code{meta}.
#'
#' @return a result table summarizing the probabilistic index and p-value
#'  of each gene.
#' @import lme4 lmerTest plyr
#' @export
#'
#' @examples
#' library(ULV)
#' data('example_data')
#' count = example_data$count_matrix
#' count = count[1:10,]
#' meta = example_data$metadata
#'
#' res_table = fit_ULV(count, meta, normalize_option='pooling', subject_name = 'donor',
#'                 cond_name = 'group_per_sample', ctrl_cond = 'mild',
#'                 case_cond = 'severe', weighted = TRUE,
#'                 covariate_name_list=c('age_yr','sex'))

fit_ULV <- function(count, meta, normalize_option='none',
                    subject_name, cond_name,
                    ctrl_cond, case_cond,
                    weighted = FALSE,
                    covariate_name_list=NULL) {

  #-------------------------------------------------
  # preprocessing
  #-------------------------------------------------

  if(normalize_option %in% c('pooling', 'clr', 'none')){
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

  if(!is.null(covariate_name_list)){
    ncovar = length(covariate_name_list)
    for (c in 1:ncovar) {
      if (!covariate_name_list[c] %in% colnames(meta)){
        stop('At least one of covariate names are not in meta info.')
      }
    }
  }

  #---------------------------------------------------
  # ULV analysis
  #---------------------------------------------------

  res_table = Reduce(plyr::rbind.fill, lapply(1:ngene, function(g){

    # message('gene ', g)
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
    d = d - 0.5

    id_case = rep((n0+1):(n0+n1), each=n0) # case
    id_ctrl = rep(1:n0, n1) # control
    d.latent = data.frame(d, id_case, id_ctrl)

    #-----------------------------------------------------
    # design matrix for covariate adjustment
    #-----------------------------------------------------

    if(!is.null(covariate_name_list)){
      ncovar = length(covariate_name_list)
      for (c in 1:ncovar) {
        if (class(meta[,covariate_name_list[c]]) %in% c('numeric', 'integer')){
          diff = covarDiff(meta=meta, subject_name=subject_name,
                           covariate_name=covariate_name_list[c],
                           ctrl_subjs = ctrl_subjs, case_subjs = case_subjs,
                           id_ctrl = d.latent$id_ctrl, id_case = d.latent$id_case, scale=1)
          d.latent[,paste0(covariate_name_list[c],'_diff')] = diff
        } else if (class(meta[,covariate_name_list[c]])=='character'){
          diff = covarDiffDummy(meta=meta, subject_name=subject_name,
                                covariate_name=covariate_name_list[c], ref_name='female',
                                ctrl_subjs=ctrl_subjs, case_subjs=case_subjs,
                                id_ctrl=d.latent$id_ctrl, id_case=d.latent$id_case)
          d.latent = cbind(d.latent, diff)
        } else {
          stop('class of covariates should be either numeric/integer(for numeric covariate)
               or character(for categorical covariate).')
        }
      }
    }

    covariate_name_list_new = colnames(d.latent)[!colnames(d.latent) %in% c('d', 'id_case', 'id_ctrl')]

    #-----------------------------------------------
    # weighted by cell number of each subject
    #-----------------------------------------------
    cell_num = sapply(y.split, length)
    cell_case = rep(cell_num[(n0+1):(n0+n1)], each=n0)
    cell_ctrl = rep(cell_num[1:n0], n1)
    z_case = rep(1/sqrt(cell_num[(n0+1):(n0+n1)]), each=n0)
    z_ctrl = rep(1/sqrt(cell_num[1:n0]), n1)
    d.latent$cell_case = cell_case
    d.latent$cell_ctrl = cell_ctrl
    d.latent$z_case = z_case
    d.latent$z_ctrl = z_ctrl

    ###########################
    d.latent$id_case = as.factor(d.latent$id_case)
    d.latent$id_ctrl = as.factor(d.latent$id_ctrl)

    #------------------------------------------------------
    # fit mixed-effect model
    #------------------------------------------------------

    if (is.null(covariate_name_list)){
      if (weighted){
        model_fit = try(lmer(d ~ 1 + (z_case-1 | id_case) + (z_ctrl-1 | id_ctrl),
                             data = d.latent, REML = FALSE))
      }else{
        model_fit = try(lmer(d ~ 1 + (1 | id_case) + (1 | id_ctrl),
                             data = d.latent, REML = FALSE))
      }
    }else{

      if (weighted){
        model_fit = try(lmer(as.formula(paste0('d ~ 1 +',
                                               paste0(covariate_name_list_new,
                                                      sep = ' + ', collapse = ''),
                                               '(z_case-1 | id_case) + (z_ctrl-1 | id_ctrl)')),
                             control = lmerControl(optimizer = "Nelder_Mead"),
                             data = d.latent, REML = FALSE))
      }else{
        model_fit = try(lmer(as.formula(paste0('d ~ 1 +',
                                               paste0(covariate_name_list_new, sep = ' + ', collapse = ''),
                                               '(1 | id_case) + (1 | id_ctrl)')),
                             data = d.latent, REML = FALSE))
      }
    }

    ## check if there is any error or convergence issue in model fitting
    if('try-error' %in% class(model_fit)){

      message('feature ', g, ': an error occurred during model fitting!')
      res = data.frame(PI = NA,
                       PI.SE = NA,
                       vcov.case = NA,
                       vcov.ctrl = NA,
                       conv_info = 'fitting_error',
                       pval = NA)
    }else{

      w = try(model_fit@optinfo$conv$lme4$messages)
      conv_info = try(is.null(w))

      df = as.data.frame(VarCorr(model_fit))
      vcov1 = df[df$grp=='id_case','vcov']
      vcov2 = df[df$grp=='id_ctrl','vcov']

      if(conv_info){
        res = try(data.frame(PI = summary(model_fit)$coefficients[1,1]+0.5,
                             PI.SE = summary(model_fit)$coefficients[1,2],
                             vcov.case = vcov1,
                             vcov.ctrl = vcov2,
                             conv_info = 'converge',
                             pval = summary(model_fit)$coefficients[1,5]))
      }else{
        message('feature ', g, ': model fitting did not converged!')
        print(w)
        res = try(data.frame(PI = summary(model_fit)$coefficients[1,1]+0.5,
                             PI.SE = summary(model_fit)$coefficients[1,2],
                             vcov.case = vcov1,
                             vcov.ctrl = vcov2,
                             conv_info = 'not_converge',
                             pval = summary(model_fit)$coefficients[1,5]))
      }
    }
    res
  }))

  #---------------------------------------------------
  # create summary table
  #---------------------------------------------------

  res_table$padj = p.adjust(res_table$pval, 'BH')
  rownames(res_table) = gene_names
  # res_table = res_table[which(!is.na(res_table[,1])),]

  return(res_table)
}
