#' Create boxplot of feature from multiple subjects
#'
#' \code{plot_multi_subj_dist} takes the count matrix and cell-level metadata,
#'  and creates a boxplot to show the distribution of \code{feature}, grouped by
#'  condition and subject.
#'
#' @param count count matrix.
#' @param meta a data frame of meta information.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param cond_name a character for condition name in \code{meta}.
#' @param subj_name a character for subject name in \code{meta}.
#' @param feature a character for the feature to plot.
#' @param sub_cond a vector of conditions to plot.
#' @param y_range The range of y-axis in the plot. Set to NULL if not specified.
#' @param show_subj whether or not to show the subject ID on x-axis.
#' @param xlab text label on x-axis of the plot.
#' @param ylab text label on y-axis of the plot.
#' @param title title of the plot. Set to \code{feature} if not specified.
#'
#' @return A \code{\link{ggplot2} object.
#' @import ggplot2
#' @export
#'
#' @examples
#' data('example_data')
#' count = example_data$count_matrix
#' count = count[1:10,]
#' meta = example_data$metadata
#' feat = rownames(count)[1]
#'
#' plot_multi_subj_dist(count, meta, normalize_option = 'none',
#'                     cond_name = 'group_per_sample', subj_name = 'donor',
#'                     feature = feat)
plot_multi_subj_dist <- function(count, meta, normalize_option = 'none',
                                 cond_name, subj_name, feature, sub_cond = NULL,
                                 y_range = NULL, show_subj = FALSE,
                                 xlab = 'Subject', ylab = NULL, title = feature){

  if(is.null(rownames(count))){
    stop('The count input should have rownames.')
  }
  stopifnot(feature %in% rownames(count))
  stopifnot(cond_name %in% colnames(meta), subj_name %in% colnames(meta))
  #---------------------------
  # pre-processing: normalize
  #---------------------------

  if(normalize_option %in% c('pooling', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize_data(count, meta, option = normalize_option)
    if(is.null(ylab)){
      ylab = 'Normalized expression'
    }
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize_data(count, meta, option = normalize_option)
    if(is.null(ylab)){
      ylab = 'Expression'
    }
  }else {
    stop('The normalizatio option argument must be one of the following options: pooling, clr, or none.')
  }

  y = as.numeric(count[feature,])
  subj_ids = meta[subj_name][,1]
  cond_ids = meta[cond_name][,1]

  # create a data frame for the figure
  df = data.frame(y = y,
                  subj_ids = subj_ids,
                  cond_ids = cond_ids)

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

  #-------------------
  # create boxplots
  #-------------------
  if(is.null(y_range)){
    fig <-df %>%
            ggplot(aes(x = subj_ids, y = y, color = cond_ids)) +
            geom_boxplot() +
            scale_color_manual(values = c('#80c241', '#eb91bd')) +
            theme_bw() +
            ggtitle(title) +
            xlab(xlab) +
            ylab(ylab) +
            theme(axis.text.x = axis,
                  legend.position = 'none') +
            facet_wrap(.~cond_ids, scales = 'free')
  }else{
    stopifnot(length(y_range)==2)
    fig <- df %>%
              ggplot(aes(x = subj_ids, y = y, color = cond_ids)) +
              geom_boxplot() +
              scale_color_manual(values = c('#80c241', '#eb91bd')) +
              theme_bw() +
              ggtitle(title) +
              xlab(xlab) +
              ylab(ylab) +
              theme(axis.text.x = axis,
                    legend.position = 'none') +
              facet_wrap(.~cond_ids, scales = 'free') +
              ylim(y_range[1], y_range[2])
  }
  return(fig)
}

