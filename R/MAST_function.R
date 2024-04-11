#' Fit MAST mixed effect model
#'
#' @param count Raw count matrix.
#' @param meta a data frame of metadata.
#' @param subj_name a character for the subject name in \code{meta}.
#' @param cond_name a character for the condition name in \code{meta}.
#' @param case_cond a character for case name.
#' @param ctrl_cond a character for control name.
#' @param strictConvergence logical (default: TRUE) return results even when the model fitting
#'  complains about convergence
#' @param nAGQ nAGQ = 1 (default) will have more accurate results but longer runtime
#'  and a higher chance of convergence failures, otherwise use nAGQ = 0
#' @param categorical_covar a list of character values for the categorical covariate names
#' @param numerical_covar a list of character values for the numerical covariate names
#'
#' @return A data frame of MAST model fitting results containing log2-fold change, p-values, and BH-adjusted p-values.
#'
#' @import MAST dplyr
#' @importFrom SingleCellExperiment colData
#' @export
#'
#' @examples
#' data('example_data')
#' count = example_data$count_matrix
#' count = count[1:100,]
#' meta = example_data$metadata
#'
#' res_table = run_MAST_mixed(count, meta,
#'                            subj_name = 'donor',
#'                            cond_name = 'group_per_sample',
#'                            ctrl_cond = 'mild',
#'                            case_cond = 'severe',
#'                            strictConvergence = FALSE,
#'                            nAGQ = 0)
run_MAST_mixed = function(count, meta,
                          subj_name, cond_name,
                          case_cond, ctrl_cond,
                          strictConvergence = TRUE,
                          nAGQ = 1,
                          categorical_covar = NULL,
                          numerical_covar = NULL){

  # subset the cells in the two groups of comparison
  count = count[, meta[,cond_name] %in% c(case_cond, ctrl_cond)]
  meta = meta[meta[,cond_name] %in% c(case_cond, ctrl_cond),]

  # use raw count matrix as input
  if(!all(count == floor(count))){
    stop('Use raw count matrix for MAST.')
  }
  stopifnot(nAGQ %in% c(0,1))
  stopifnot(all(categorical_covar %in% colnames(meta)))
  stopifnot(all(numerical_covar %in% colnames(meta)))

  # count matrix normalization: log2(TPM+1)
  rds = colSums(count)
  med_rds = median(rds)
  count = t(t(count)/rds)*med_rds
  count = log1p(count)

  # create the SingleCellAssay (sca) object
  cell_id = colnames(count)
  gene_id = rownames(count)
  fData = data.frame(primerid = gene_id)
  cData = data.frame(wellKey = cell_id)
  sca = MAST::FromMatrix(as.matrix(count), cData, fData)

  # add variable to metadata
  cdr = scale(colSums(count>0)) # how to calculate cdr
  diagnosis = factor(meta[,cond_name], levels=c(ctrl_cond, case_cond))
  colData(sca)$cngeneson = as.numeric(cdr)
  colData(sca)$diagnosis = diagnosis
  colData(sca)$ind = as.factor(meta[,subj_name])

  # add covariate to metadata if you have any
  sca_colData = colData(sca)
  if(!is.null(categorical_covar)){
    for (i in 1:length(categorical_covar)) {
      var_name = categorical_covar[i]
      sca_colData[, var_name] = as.factor(meta[, var_name])
    }
  }
  if(!is.null(numerical_covar)){
    for (i in 1:length(numerical_covar)) {
      var_name = numerical_covar[i]
      sca_colData[, var_name] = scale(meta[, var_name])
    }
  }
  colData(sca) = sca_colData
  covariate_name_list = c(categorical_covar, numerical_covar)

  # run MAST-mixed effect model
  if(is.null(covariate_name_list)){
    b0 = zlm(formula = ~ diagnosis + cngeneson + ( 1 | ind ), sca = sca,
             strictConvergence = strictConvergence,
             fitArgsD = list(nAGQ = nAGQ),
             force = TRUE, silent = FALSE,
             method = 'glmer', ebayes = FALSE, parallel = TRUE)
  }else{
    b0 = zlm(formula = as.formula(paste0('~ diagnosis + cngeneson + (1 | ind) +',
                                         paste(covariate_name_list, collapse = '+'))),
             sca = sca,
             strictConvergence = strictConvergence,
             fitArgsD = list(nAGQ = nAGQ),
             force = TRUE, silent = FALSE,
             method = 'glmer', ebayes = FALSE, parallel = TRUE)
  }
  # lrt = MAST::lrTest(b0, 'diagnosis')

  lrt_term = paste0('diagnosis', case_cond)
  summary_cond = summary(b0, doLRT = lrt_term, fitArgsD = list(nAGQ = nAGQ))

  summary_Dt = summary_cond$datatable
  # fcHurdle = merge(summary_Dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
  #                  summary_Dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  # fcHurdle[,padj:=p.adjust(`Pr(>Chisq)`, 'BH')]
  hurdle_p_values <- summary_Dt %>%
    filter(contrast == lrt_term, component == 'H') %>%
    select(primerid, `Pr(>Chisq)`)

  logfc_coefficients <- summary_Dt %>%
    filter(contrast == lrt_term, component == 'logFC') %>%
    select(primerid, coef, ci.hi, ci.lo)

  fcHurdle <- merge(hurdle_p_values, logfc_coefficients, by = 'primerid')

  fcHurdle <- fcHurdle %>%
    mutate(padj = p.adjust(`Pr(>Chisq)`, method = 'BH'))

  fcHurdle_df = fcHurdle %>%
    as_tibble() %>%
    # arrange(padj) %>%
    dplyr::rename(feature = primerid,
                  pval = `Pr(>Chisq)`,
                  log2FC = coef)

  # pval = apply(lrt, 1, function(x){x[3,3]})
  return(fcHurdle_df)
}


#' Fit MAST glm effect model
#'
#' @param count Raw count matrix.
#' @param meta a data frame of metadata.
#' @param subj_name a character for the subject name in \code{meta}.
#' @param cond_name a character for the condition name in \code{meta}.
#' @param case_cond a character for case name.
#' @param ctrl_cond a character for control name.
#' @param categorical_covar a list of character values for the categorical covariate names
#' @param numerical_covar a list of character values for the numerical covariate names
#'
#' @return A data frame of MAST model fitting results containing log2-fold change, p-values, and BH-adjusted p-values.
#'
#' @import MAST dplyr
#' @importFrom SingleCellExperiment colData
#' @export
#'
#' @examples
run_MAST_glm = function(count, meta,
                          subj_name, cond_name,
                          case_cond, ctrl_cond,
                          categorical_covar = NULL,
                          numerical_covar = NULL){

  # subset the cells in the two groups of comparison
  count = count[, meta[,cond_name] %in% c(case_cond, ctrl_cond)]
  meta = meta[meta[,cond_name] %in% c(case_cond, ctrl_cond),]

  # use raw count matrix as input
  if(!all(count == floor(count))){
    stop('Use raw count matrix for MAST.')
  }

  stopifnot(all(categorical_covar %in% colnames(meta)))
  stopifnot(all(numerical_covar %in% colnames(meta)))

  # count matrix normalization: log2(TPM+1)
  rds = colSums(count)
  med_rds = median(rds)
  count = t(t(count)/rds)*med_rds
  count = log1p(count)

  # create the SingleCellAssay (sca) object
  cell_id = colnames(count)
  gene_id = rownames(count)
  fData = data.frame(primerid = gene_id)
  cData = data.frame(wellKey = cell_id)
  sca = MAST::FromMatrix(as.matrix(count), cData, fData)

  # add variable to metadata
  cdr = scale(colSums(count>0)) # how to calculate cdr
  diagnosis = factor(meta[,cond_name], levels=c(ctrl_cond, case_cond))
  colData(sca)$cngeneson = as.numeric(cdr)
  colData(sca)$diagnosis = diagnosis
  colData(sca)$ind = as.factor(meta[,subj_name])

  # add covariate to metadata if you have any
  sca_colData = colData(sca)
  if(!is.null(categorical_covar)){
    for (i in 1:length(categorical_covar)) {
      var_name = categorical_covar[i]
      sca_colData[, var_name] = as.factor(meta[, var_name])
    }
  }
  if(!is.null(numerical_covar)){
    for (i in 1:length(numerical_covar)) {
      var_name = numerical_covar[i]
      sca_colData[, var_name] = scale(meta[, var_name])
    }
  }
  colData(sca) = sca_colData
  covariate_name_list = c(categorical_covar, numerical_covar)

  # run MAST-glm model
  if(is.null(covariate_name_list)){
    b0 = zlm(formula = ~ diagnosis + cngeneson , sca = sca,
             force = TRUE, silent = FALSE,
             method = 'glm', ebayes = FALSE, parallel = TRUE)
  }else{
    b0 = zlm(formula = as.formula(paste0('~ diagnosis + cngeneson +',
                                         paste(covariate_name_list, collapse = '+'))),
             sca = sca,
             force = TRUE, silent = FALSE,
             method = 'glm', ebayes = FALSE, parallel = TRUE)
  }
  # lrt = MAST::lrTest(b0, 'diagnosis')

  lrt_term = paste0('diagnosis', case_cond)
  summary_cond = summary(b0, doLRT = lrt_term)

  summary_Dt = summary_cond$datatable
  # fcHurdle = merge(summary_Dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
  #                  summary_Dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  # fcHurdle[,padj:=p.adjust(`Pr(>Chisq)`, 'BH')]
  hurdle_p_values <- summary_Dt %>%
    filter(contrast == lrt_term, component == 'H') %>%
    select(primerid, `Pr(>Chisq)`)

  logfc_coefficients <- summary_Dt %>%
    filter(contrast == lrt_term, component == 'logFC') %>%
    select(primerid, coef, ci.hi, ci.lo)

  fcHurdle <- merge(hurdle_p_values, logfc_coefficients, by = 'primerid')

  fcHurdle <- fcHurdle %>%
    mutate(padj = p.adjust(`Pr(>Chisq)`, method = 'BH'))

  fcHurdle_df = fcHurdle %>%
    as_tibble() %>%
    # arrange(padj) %>%
    dplyr::rename(feature = primerid,
                  pval = `Pr(>Chisq)`,
                  log2FC = coef)

  # pval = apply(lrt, 1, function(x){x[3,3]})
  return(fcHurdle_df)
}


