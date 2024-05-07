#' Make a boxplot of feature at subject level and a heatmap of pairwise difference
#'
#' @param count count matrix.
#' @param meta a data frame of meta information.
#' @param normalize_option a character value to show the normalization method to
#'  apply to the count matrix.
#' @param cond_name a character for condition name in \code{meta}.
#' @param subject_name a character for subject name in \code{meta}.
#' @param ctrl_cond a character for control name.
#' @param case_cond a character for case name.
#' @param feature a character for the feature to plot.
#' @param xlab text label on x-axis for the boxplot.
#' @param ylab text label on y-axis for the boxplot.
#' @param title title for the boxplot.
#'
#' @return A \code{\link{ggplot2} object.
#' @import dplyr ggplot2 ComplexHeatmap circlize grid ggplotify patchwork
#' @export
#'
#' @examples
plot_dij_heatmap <- function(count, meta, normalize_option = 'pooling',
                             cond_name, subject_name, ctrl_cond, case_cond,
                             feature, xlab='', ylab='Normalized expression',
                             title = feature){
  # library(ComplexHeatmap)
  # library(circlize)
  # library(grid)
  # library(ggplotify)
  # library(patchwork)
  #-------------------------------------------------
  # preprocessing
  #-------------------------------------------------

  if(normalize_option %in% c('pooling', 'clr')){
    message('Normalizing the input count matrix using ', normalize_option)
    count = normalize_data(count, meta, option = normalize_option)
  }else if (normalize_option == 'none'){
    message('Use the input count matrix directly. No normalization was utilized.')
    count = normalize_data(count, meta, option = normalize_option)
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

  #-------------------------------
  # extract feature information
  #-------------------------------
  g = which(rownames(count)==feature)
  y = as.numeric(count[g,])
  y.split = split(y, ind_fct)

  #------------------------------------------
  # make boxplot for distribution of subject
  #------------------------------------------
  # df = data.frame(y = y,
  #                 subject = subject,
  #                 condition = condition)
  # max.y = max(df$y)
  # fig1 = df %>%
  #   ggplot(aes(x = subject, y = y, color = condition)) +
  #   geom_boxplot() +
  #   scale_color_manual(values = c('#80c241', '#eb91bd')) +
  #   facet_wrap(.~condition, scales = 'free') +
  #   ggtitle(title) +
  #   xlab(xlab) +
  #   ylab(ylab) +
  #   theme_bw(base_size = 22) +
  #   theme(panel.grid.minor = element_blank(),
  #         panel.grid.major.x = element_blank(),
  #         axis.line = element_line(colour = "black",linewidth = 1),
  #         axis.text = element_text(color = 'black'),
  #         strip.text = element_text(size = 22),
  #         strip.background = element_rect(linewidth = 2),
  #         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  #         legend.position = 'none') +
  #   ylim(0, max.y)

  #---------------------------------------------------------
  # compute pairwise difference d_ij
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

  #------------------------
  # make heatmap for d_ij
  #------------------------
  col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  fig2 <- Heatmap(d_mat,
            col = col_fun,
            name = 'Probabilistic Index',
            rect_gp = gpar(col = 'white', lwd = 1),
            row_names_gp = gpar(fontsize = 15, fontface = 'plain'),
            column_names_gp = gpar(fontsize = 15, fontface = 'plain'),
            column_names_rot = 60,
            heatmap_legend_param = list(grid_width = unit(0.6, 'cm'),
                                        legend_width = unit(8, "cm"),
                                        title_gp = gpar(fontface = "plain", fontsize = 18),
                                        labels_gp = gpar(fontsize = 15),
                                        direction = 'horizontal',
                                        title_position = 'topcenter'),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.2f", d_mat[i, j]), x, y, gp = gpar(fontsize = 18, col = 'black'))},
            show_row_dend=FALSE, show_column_dend = FALSE,
            cluster_rows = FALSE, cluster_columns = FALSE)

  #------------------------------
  # combine two plots together
  #------------------------------
  grob2 = grid.grabExpr(draw(fig2, heatmap_legend_side = "bottom"))
  gg2 = as.ggplot(grob2)
  # fig = fig1 / gg2 +
  #   plot_layout(heights = c(1,2))

  return(gg2)
}
