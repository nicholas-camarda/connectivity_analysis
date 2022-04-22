#' insert at around line 970 of heatmaps.R and add back in left annot in heatmap 
#' body call to get diff left hand annotation

base_cell_indices <- which(grp_color_df$grp_fac == base_cell_id)
base_grp_color <- unique(grp_color_df %>% 
                           filter(grp_fac == base_cell_id) %>%
                           .$color)
rest_cell_color <- "darkgray"
rest_cell_indices <- which(grp_color_df$grp_fac != base_cell_id)

# scale of the boxplot annotation
rg <- range(reordered_mat, na.rm = TRUE)
rg[1] = rg[1] - (rg[2] - rg[1])* 0.02
rg[2] = rg[2] + (rg[2] - rg[1])* 0.02

anno_multiple_boxplot = function(index) {
  # rows of matrix
  nr = length(index)
  # message(nr)
  pushViewport(viewport(xscale = rg, yscale = c(0.5, nr + 0.5)))
  for(i in seq_along(index)) {
    grid.rect(y = nr-i+1, height = 1, default.units = "native")
    # first group
    # grid.abline(intercept = 0, gp = gpar(lty = 3,lwd = 3, col = "red"))
    grid.boxplot(reordered_mat[ index[i], base_cell_indices] %>% na.omit(),
                 pos = nr-i+1 + 0.2,
                 box_width = 0.4, size = unit(1, "mm"),
                 gp = gpar(fill = base_grp_color),
                 direction = "horizontal")
    # second group
    grid.boxplot(reordered_mat[ index[i], rest_cell_indices] %>% na.omit(),
                 pos = nr-i+1 - 0.2,
                 box_width = 0.4,
                 size = unit(1, "mm"),
                 gp = gpar(fill = rest_cell_color),
                 direction = "horizontal")

  }
  grid.xaxis(gp = gpar(fontsize = 8))
  popViewport()
}

# save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-p100-heat.RData")
# load("debug/debug_dat/debug-p100-heat.RData")
# stop()

# get max of finite values
bar_plot_q_values <- signif_log_q_for_bar_plot
bar_plot_q_values[is.infinite(bar_plot_q_values)] = max(bar_plot_q_values[is.finite(bar_plot_q_values)], na.rm = TRUE)

left_ha <- rowAnnotation(
  boxplot = anno_multiple_boxplot,
  `bar` = anno_barplot(bar_plot_q_values,
                       gp = gpar(fontsize = 8, fill = "darkgray"),
                       axis_param = list(labels_rot = 0)),
  width = unit(4, "cm"),
  gap = unit(2, "mm"),
  # annotation_name_rot = 45,
  show_annotation_name = FALSE
  # show_annotation_names = FALSE
# signif_log_q_for_bar_plot)
  # width = unit(4, "cm"),
  # gap = unit(1.5, "mm")
  # annotation_legend_param = list(`bar` = list(title = "Q-value"))
)
draw(left_ha)
