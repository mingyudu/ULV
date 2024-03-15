#' Create volcano plots
#'
#' \code{plot_volcano} create volcano plots. This function takes the ULV
#' analysis result table and produces volcano plots with the horizontal
#' axis representing the estimated probabilistic index.
#'
#' @param res_table a data frame with the column 'PI' representing the
#'  probabilistic index, and the column 'padj' representing the adjusted
#'  p-values of ULV analysis. The rownames should be the names of the features
#'  (genes, proteins, etc.). The format is the same as the output from ULV
#'  function.
#' @param alpha the threshold of adjusted p-values. The default value is set
#'  as 0.1.
#' @param PI_thres the threshold of estimated probabilistic index.
#'  The default value is set as 0.05.
#' @param title Title of the volcano plot.
#' @param xlab text label on x-axis of the plot.
#' @param ylab text label on y-axis of the plot.
#' @param add_labels whether or not to add text labels for the significant
#'  features on the plot.
#' @param max.overlaps Equivalent of max.overlaps in ggrepel.
#'  Set to 'Inf' to always display all labels.
#' @import dplyr ggplot2 ggrepel
#'
#' @return A ggplot2 figure.
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
#'                     cond_name = 'group_per_sample', ctrl_cond = 'mild',
#'                     case_cond = 'severe', weighted = TRUE,
#'                     covariate_name_list=c('age_yr','sex'))
#' plot_volcano(res_table)
plot_volcano <- function(res_table, alpha = 0.1, PI_thres = 0.05,
                         title = 'Volcano Plot', xlab = 'Estimated probabilistic index',
                         ylab = '-log10(adjusted p-value)',
                         add_labels = TRUE, max.overlaps = 15){

  if(is.null(rownames(res_table))){
    stop('The input table for volcano plotting is not valid. It should have row names representing the features (genes, proteins, etc.).')
  }
  res_table$name = rownames(res_table)

  if (!all(c('padj', 'PI') %in% colnames(res_table))){
    stop('The column names of res_table should include padj and PI.')
  }
  # subset significant features
  significant_features = subset(res_table, padj < alpha & abs(PI - 0.5) > PI_thres)

  # create the volcano plot
  fig = res_table %>%
    as.data.frame() %>%
    mutate(significant = ifelse(padj < alpha & abs(PI - 0.5) > PI_thres,
                                'Sig', 'Not Sig')) %>%
    ggplot(aes(x = PI, y = -log10(padj))) +
    geom_point(aes(color = significant, size = significant, alpha = 0.5)) +
    scale_color_manual(values = c('#56B4E9', 'red')) +
    scale_size_manual(values = c(0.5, 0.6)) +
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw()

  # Add threshold lines
  fig = fig +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-PI_thres + 0.5, PI_thres + 0.5), linetype = "dashed", linewidth = 0.3) +
    theme(legend.position = 'none')

  # Add text labels for significant features
  if(add_labels){
    fig = fig +
      geom_label_repel(
        data = significant_features,
        aes(label = name),
        size = 2,
        label.size = 0.25,
        box.padding = 0.5,
        max.overlaps = max.overlaps,
        min.segment.length = 0,
        arrow = arrow(length = unit(0.01, 'npc')),
        segment.size = 0.3
      )
  }
  return(fig)
}
