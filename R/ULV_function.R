#' ULV: U-statistic based Latant Variable model that takes advantage of the robustness of rank-based methods and the statistical efficiency of parametric methods for small sample sizes
#'
#' @param count count matrix
#' @param meta a data frame of meta information
#' @param normalize if TRUE, then conduct the normalization during ULV.
#' @param subject_name a character for subject name in meta
#' @param cond_name a character for condition name in meta
#' @param ctrl_cond a character for control name
#' @param case_cond a character for case name
#' @param weighted a binary variable indicating whether the analysis is weighted by varying cluster sizes
#' @param covariate_name_list a vector of character of covariate names. If not, set it to NULL
#'
#' @return a result table summarizing the probabilistic index and p-value of each gene
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
#' res_table = ULV(count, meta, normalize=TRUE, subject_name = 'donor', cond_name = 'group_per_sample', ctrl_cond = 'mild', case_cond = 'severe', weighted = TRUE, covariate_name_list=c('age_yr','sex'))

ULV <- function(count, meta, normalize=TRUE,
                    subject_name, cond_name,
                    ctrl_cond, case_cond,
                    weighted = FALSE,
                    covariate_name_list=NULL) {

  #-------------------------------------------------
  # preprocessing
  #-------------------------------------------------

  if(normalize){
    count = normalize(count, meta)
  }
  gene_names = rownames(count)
  ngene = nrow(count)

  #-------------------------------------------------
  # extract condition and subject information
  #-------------------------------------------------

  condition = meta[, cond_name]
  subject = meta[, subject_name]

  table = table(subject, condition)
  ctrl_subjs = rownames(table)[table[,ctrl_cond]!=0]
  case_subjs = rownames(table)[table[,case_cond]!=0]
  ind_fct = factor(subject, levels = c(ctrl_subjs, case_subjs))
  n0 = length(ctrl_subjs); n1 = length(case_subjs)

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

    message('gene ', g)
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
                           id_ctrl = d.latent$id_ctrl, id_case = d.latent$id_case, scale=20) # scale?
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

    res = try(data.frame(PI = summary(model_fit)$coefficients[1,1]+0.5,
                         pval = summary(model_fit)$coefficients[1,5]))
    if('try-error' %in% class(res)){
      message('An error occurred during model fitting for gene ', gene_names[g])
      res = data.frame(PI = NA)
    }
    res
  }))

  #---------------------------------------------------
  # create summary table
  #---------------------------------------------------

  res_table$padj = p.adjust(res_table$pval, 'BH')
  rownames(res_table) = gene_names
  res_table = res_table[which(!is.na(res_table[,1])),]

  return(res_table)
}
