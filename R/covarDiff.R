#' Compute the subject-level covariate difference for covariate adjustment
#'
#' @param meta a data frame of meta information
#' @param subject_name a character for subject name in meta
#' @param covariate_name a character for covariate name in meta
#' @param ctrl_subjs a vector of character for all control subjects
#' @param case_subjs a vector of character for all case subjects
#' @param id_ctrl index of control subject to compare
#' @param id_case index of case subject to compare
#' @param scale scale for covariate difference. If it is NULL, we will scale the covariate difference by its standard deviation.
#'
#' @return a vector of numeric values representing the covariate differences
#' @export
#'
#' @examples
covarDiff <- function(meta, subject_name, covariate_name,
                      ctrl_subjs, case_subjs,
                      id_ctrl, id_case, scale=1){

  subj_info = unique(meta[,c(subject_name, covariate_name)])
  rownames(subj_info) = subj_info[,subject_name]

  n0 = length(ctrl_subjs); n1 = length(case_subjs)
  covar_case = subj_info[case_subjs[id_case - n0], covariate_name]
  covar_ctrl = subj_info[ctrl_subjs[id_ctrl], covariate_name]
  covar_diff = (covar_case - covar_ctrl)

  if(is.null(scale)){
    scale = sd(covar_diff)
    covar_diff = covar_diff/scale
  }else if(is.numeric(scale)){
    covar_diff = covar_diff/scale
  }else{
    stop('scale must be either NULL or a numeric value.')
  }
  return(covar_diff)
}
