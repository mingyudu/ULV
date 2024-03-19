#' Create boxplot of feature and subject-level summary statistic metric with covariate information
#'
#' @param count count matrix.
#' @param meta a data frame of meta information.
#' @param res_table result table of ULV analysis.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param cond_name a character for condition name in \code{meta}.
#' @param subj_name a character for subject name in \code{meta}.
#' @param feature a character for the feature to plot.
#' @param covar_name a character for covariate name in \code{meta}.
#' @param sub_cond a vector of conditions to plot.
#' @param xlab.p1 text label on x-axis of the first plot.
#' @param ylab.p1 text label on y-axis of the first plot.
#' @param title.p1 title of the first plot. Set to \code{feature} if not specified.
#' @param xlab.p2 text label on x-axis of the second plot.
#' @param ylab.p2 text label on y-axis of the second plot.
#' @param title.p2 title of the second plot.
#'
#' @return A \code{\link{ggplot2} object.
#' @import ggplot2 stringr patchwork
#' @export
#'
#' @examples
plot_dist_by_covar <- function(count, meta, res_table = NULL,
                               normalize_option = 'pooling',
                               cond_name, subj_name, feature, covar_name,
                               sub_cond = NULL, xlab.p1 = covar_name, ylab.p1 = NULL,
                               title.p1 = feature,
                               xlab.p2 = covar_name, ylab.p2 = 'Subject-level metric',
                               title.p2 = 'Subject-level metric'){

  if(is.null(rownames(count))){
    stop('The count input should have rownames.')
  }
  stopifnot(feature %in% rownames(count))
  stopifnot(cond_name %in% colnames(meta),
            subj_name %in% colnames(meta),
            covar_name %in% colnames(meta))
  #---------------------------
  # pre-processing: normalize
  #---------------------------

  if(normalize_option %in% c('pooling', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize(count, meta, option = normalize_option)
    if(is.null(ylab.p1)){
      ylab.p1 = 'Normalized expression'
    }
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize(count, meta, option = normalize_option)
    if(is.null(ylab.p1)){
      ylab.p1 = 'Expression'
    }
  }else {
    stop('The normalizatio option argument must be one of the following options: pooling, clr, or none.')
  }

  y = as.numeric(count[feature,])
  subj_ids = meta[subj_name][,1]
  cond_ids = meta[cond_name][,1]
  covar = meta[covar_name][,1]

  # create a data frame for the figure
  df = data.frame(y = y,
                  subj_ids = subj_ids,
                  cond_ids = cond_ids,
                  covar = covar)

  # only include sub_condition in the plot
  if(!is.null(sub_cond)){
    stopifnot(all(sub_cond %in% df$cond_ids))
    df = df[df$cond_ids %in% sub_cond,]
    df$cond_ids = factor(df$cond_ids, levels = sub_cond)
  }

  # whether or not to show subject ID on x-axis
  if(show_subj){
    axis = element_text(angle = 90, vjust = 0.5, hjust=1)
  }else{
    axis = element_blank()
  }

  #------------------
  # create boxplot
  #------------------
  if(class(covar) %in% c('numeric', 'integer')){
    fig1 <- df %>%
      ggplot(aes(x = covar, y = y, color = cond_ids, group = subj_ids)) +
      geom_boxplot() +
      theme_bw() +
      ggtitle(title.p1) +
      facet_wrap(.~cond_ids) +
      xlab(xlab.p1) +
      ylab(ylab.p1) +
      theme(legend.position = 'none')
  }else if(class(covar)=='character'){ ################## need to fix
    fig1 <- df %>%
      ggplot(aes(x = covar, y = y, color = cond_ids, group = subj_ids)) +
      geom_boxplot() +
      theme_bw() +
      ggtitle(title.p1) +
      facet_wrap(.~cond_ids) +
      xlab(xlab.p1) +
      ylab(ylab.p1) +
      theme(legend.position = 'none')
  }

  if(!is.null(res_table)){
    covar_ind = which(str_detect(colnames(res_table), 'diff') & !str_detect(colnames(res_table), 'pval'))
    est_pi = res_table[feature, 'PI']
    est_coeff = res_table[feature, covar_ind]
    coeff_name = colnames(res_table)[covar_ind]
    fig1 <- fig1 +
      labs(caption = paste0('PI: ', signif(est_pi, 3), '\n',
                            coeff_name, ': ', signif(est_coeff, 3)))
  }

  #--------------------------------------------
  # summary statistic plots: mean and median
  #--------------------------------------------
  df_summary1 = df %>%
    group_by(subj_ids) %>%
    summarize(metric = mean(y), .groups = 'drop') %>%
    left_join(df %>% distinct(subj_ids, cond_ids, covar), by = "subj_ids")
  df_summary1$summary = 'mean'

  df_summary2 = df %>%
    group_by(subj_ids) %>%
    summarize(metric = median(y), .groups = 'drop') %>%
    left_join(df %>% distinct(subj_ids, cond_ids, covar), by = "subj_ids")
  df_summary2$summary = 'median'
  df_summary = rbind(df_summary1, df_summary2)

  # Plotting
  fig2 = ggplot(df_summary, aes(x = covar, y = metric,
                                color = cond_ids, group = cond_ids)) +
    geom_point() +
    facet_wrap(.~summary) +
    geom_smooth(method = lm) +
    labs(title = title.p2,
         x = xlab.p2, y = ylab.p2) +
    theme_bw() +
    theme(legend.position = 'none')

  return(fig1/fig2)
}
