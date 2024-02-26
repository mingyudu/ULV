#' Compute the subject-level categorical covariate difference for covariate adjustment
#'
#' @param meta a data frame of meta information
#' @param subject_name a character for subject name in meta
#' @param covariate_name a character for covariate name in meta
#' @param ref_name a character representing the reference level for the categorical covariate
#' @param ctrl_subjs a vector of character for all control subjects
#' @param case_subjs a vector of character for all case subjects
#' @param id_ctrl index of control subject to compare
#' @param id_case index of case subject to compare
#'
#' @return a data frame of covariate difference between a category and the reference category
#' @export
#'
#' @examples
covarDiffDummy <- function(meta, subject_name, covariate_name, ref_name=NULL,
                           ctrl_subjs, case_subjs, id_ctrl, id_case){

  subj_info = unique(meta[,c(subject_name, covariate_name)])
  rownames(subj_info) = subj_info[,subject_name]

  n0 = length(ctrl_subjs); n1 = length(case_subjs)
  covar_case = subj_info[case_subjs[id_case - n0], covariate_name]
  covar_ctrl = subj_info[ctrl_subjs[id_ctrl], covariate_name]

  cats = unique(c(covar_case, covar_ctrl))
  if(is.null(ref_name)){
    ref_name = cats[1]
  }
  if(!ref_name %in% cats){
    stop('User-defined reference level is not in the levels of categorical covariates.')
  }
  cats2 = cats[cats!=ref_name]

  # if there are M categories, then design matrix will have M-1 columns (except the reference level)
  mat_case = matrix(0, nrow = length(covar_case), ncol = length(cats)-1)
  for (i in 1:length(cats2)) {
    mat_case[,i] = ifelse(covar_case==cats2[i], 1, 0)
  }
  mat_ctrl = matrix(0, nrow = length(covar_ctrl), ncol = length(cats)-1)
  for (i in 1:length(cats2)) {
    mat_ctrl[,i] = ifelse(covar_ctrl==cats2[i], 1, 0)
  }
  mat_diff = as.data.frame(mat_case - mat_ctrl)
  colnames(mat_diff) = paste0(paste0(paste0(covariate_name, '_diff_'), cats2), '_vs_', ref_name)
  return(mat_diff)
}
