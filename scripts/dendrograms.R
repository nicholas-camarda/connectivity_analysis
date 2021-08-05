# source("~/OneDrive - Tufts/phd/jaffe/workspace/ws/scripts/master-source.R")

##########################

#' @note generate a constant labeling scheme for cells
#' @note this function gets called in plot_pretty_dendrogram to generate the informative labels from hard code
create_informative_labels <- function(){
  # num_grps_cells = 3
  
  #' colors: 
  #' https://www.colorhexa.com/
  # colors <- brewer.pal(num_grps_cells, "Dark2") # cancer, not cancer / not vascular, vascular
  # colors <- c("#1eb6ff", "#ffd81e", "#000000", "#ff671e") # pal_igv("alternating")(2)
  colors <- c("#ff677f", "#ffe766", "#225ea8", "#41b6c4", "#9ba2ff", "#ffb781")
  # names(colors) <- c("blue", "yellow", "black", "orange")
  # "#a1dab4", green
  names(colors) <- c("light red", "yellow", "dark blue", "light blue", "light purple", "orange")
  # yellow = vascular
  # purple = huvecs
  # orange = haosmcs
  # light red = non-vascular
  # dark blue = other
  # light blue = cancer
  
  vascular_char_vec <- c("HUVEC","HAoSMC")
  other_char_vec <- c("MCF10A", "NPC")
  
  informative_labels_df <- tibble(lbs = c("HAoSMC", "HUVEC", "MCF10A", "PC3", "YAPC", "A375", "A549", "MCF7", "NPC"),
                                  new_lbs = c("HAoSMC", "HUVEC", "MCF10A  \n(fibrocystic)", "PC3    \n(prostate)",
                                              "YAPC    \n(pancreas)", "A375\n(skin)", "A549\n(lung)", "MCF7\n(breast)", "NPC  \n(neural)")) %>%
    mutate(grp = ifelse(lbs %in% vascular_char_vec, "Vascular cells", 
                        ifelse(lbs %in% other_char_vec, "Other cells", "Cancer cells")),
           grp = factor(grp, levels = c("Vascular cells", "Cancer cells", "Other cells"))) %>%
    mutate(grp_l = ifelse(lbs %in% vascular_char_vec, "Vascular cells", "Non-vascular cells"),
           grp_l = factor(grp_l, levels = c("Vascular cells", "Non-vascular cells"))) %>%
    mutate(grp_l2 = ifelse(lbs == "HUVEC", "HUVEC",
                           ifelse(lbs == "HAoSMC", "HAoSMC", "Non-vascular cells")),
           grp_l2 = factor(grp_l2, levels = c("HAoSMC", "HUVEC", "Non-vascular cells"))) %>%
    mutate(color = ifelse(lbs == vascular_char_vec[1], colors[5], # check huvec
                          ifelse(lbs == vascular_char_vec[2], colors[6], # check haosmc
                                 ifelse(lbs %in% other_char_vec, colors[3], colors[4]))), # check other vs cancer
           vasc_grp_color = ifelse(lbs %in% vascular_char_vec, colors[2], colors[1])) # check vasc vs non-vasc
  
  grp_labels_ordered_for_legend_df <- informative_labels_df %>% 
    distinct(grp_l, vasc_grp_color) %>% 
    arrange(desc(levels(grp_l)))
  
  res <- list(informative_labels_df, grp_labels_ordered_for_legend_df) %>% 
    set_names("informative_labels_df", "legend_dend_labels_df")
  return(res)
}


# data <- gcp_pert_data_dmso

#' @note 
#' @param args pert or moa specific lst-dataset, from *cluster3.R*
#' @param INSET_MARGIN percent of plotsize in x direction to adjust legend, usually negative
plot_pretty_dendrogram <- function(args, 
                                   INSET_MARGIN = -0.06, 
                                   rotate_dendrogram = FALSE) {
  # args <- my_clust_obj[2,]
  
  data <- args$cluster_lst
  
  base_output_dir <- force_natural(args$path)
  
  dataset <- force_natural(args$which_dat)
  fn_name_s <- force_natural(args[1])
  plot_title <- qq("Hierarchical clustering on connectivity matrix\n@{fn_name_s} || @{dataset}")
  dendro_output_fn <- file.path(base_output_dir, qq("@{fn_name_s}.eps"))
  
  cell_id_clusters_char_vec <- force_natural(data$cluster_assignments); # cell_id_clusters_char_vec
  pvclust_obj <- force_natural(data$clust_obj); # print(pvclust_obj)
  co_clust <- force_natural(data$co_clust_bool)
  max_n_clusts <- max(cell_id_clusters_char_vec)
  
  dend_height <- max(pvclust_obj$height) # height of the first split!, total height of dendro
  # print(dend_height)
  hclust_dend <- as.dendrogram(pvclust_obj)
  
  inf_lst <- create_informative_labels()
  informative_labels_df <- inf_lst$informative_labels_df
  grp_labels_ordered_for_legend_df <- inf_lst$legend_dend_labels_df
  
  # the last element of this vector is the rightmost on a vertically oriented dendrogram
  vasc_top_labels_df <- tibble(lbs = labels(hclust_dend)) %>% 
    mutate(dend_order = 1:nrow(.)) %>%
    left_join(informative_labels_df, by="lbs") 
  # right juxtapose huvec/haosmc
  dend_order_int_vec <- vasc_top_labels_df$dend_order
  last_2_positions <- c(length(dend_order_int_vec)-1, length(dend_order_int_vec))
  names_last_2_position <- vasc_top_labels_df$lbs[last_2_positions]
  
  vasc_dend_pos <- filter(vasc_top_labels_df, grp == "Vascular cells")$dend_order
  names_vasc_dend_pos <- vasc_top_labels_df$lbs[vasc_dend_pos]
  
  if ( !any(vasc_dend_pos %in% last_2_positions)  ) { #& co_clust
    message("\n Re-organizing labels to 'right-justify' vascular cells") 
    # ensure vasc is last
    
    # switch the names and positions of the last positions with 
    names_involved <-  c(names_vasc_dend_pos,names_last_2_position)
    to_join_original_df <- filter(vasc_top_labels_df, !(lbs %in% names_involved)) 
    # %>%
    #   mutate(new_dend_order = dend_order)
    vasc_top_labels_df <- tibble(lbs = c(names_vasc_dend_pos,names_last_2_position), 
                                 dend_order = c(last_2_positions,vasc_dend_pos)) %>% # this could be renamed to 'new_dend_labels' t
      left_join(vasc_top_labels_df %>% dplyr::select(-dend_order),by="lbs") %>% # give their lbs their properties back
      bind_rows(to_join_original_df) %>% # don't duplicate names when recombining
      arrange(dend_order)
  } 
  
  
  label_order_df <- tibble(lbs = labels(hclust_dend), 
                           dendro_order = 1:length(labels(hclust_dend)))
  
  # so this should be the same as vasc_top_df but with dendro_order column
  informative_labels_df_ordered <- informative_labels_df %>% 
    inner_join(label_order_df, by="lbs") %>%  # use inner, because if the cell isn't in this data, need to remove
    arrange(dendro_order); informative_labels_df_ordered
  # arrange(desc(top_vasc_order)); informative_labels_df_ordered
  
  # assign cluster vector with new labels
  # names(cell_id_clusters_char_vec) <- informative_labels_df_ordered %>% 
  #   arrange(lbs) %>% 
  #   .$new_lbs
  
  # normal labels and colors, according to dendro
  colors_char_vec <- as.character(informative_labels_df_ordered$color)
  labels_char_vec <- as.character(informative_labels_df_ordered$new_lbs)
  
  # rotated labels and colors to have vasc on top
  if (rotate_dendrogram){
    vasc_top_labels_df_rotated_order <- vasc_top_labels_df %>%
      arrange(dend_order)
    labels_char_vec <- as.character(vasc_top_labels_df_rotated_order$new_lbs)
    # print(vasc_top_labels_df_rotated_order)
  }
  
  
  ## color the leaves
  # https://stackoverflow.com/questions/10571266/colouring-branches-in-a-dendrogram-in-r
  # https://talgalili.github.io/dendextend/articles/dendextend.html
  
  # good for vertical dendrogram
  # pdf_height <- 10
  # pdf_width <- 18
  # INSET_MARGIN_VERT <- -0.06
  
  
  pdf_height <- 21
  pdf_width <- 35
  
  setEPS()
  postscript(file = dendro_output_fn, height = pdf_height, width = pdf_width, onefile = F)
  # pdf(file = dendro_output_fn, width = pdf_width, height = pdf_height, onefile = FALSE) # for some reason, without onefile = F, the pdf gets printed on the second page, with a blank page to start.. weird
  plot.new(); i <- 0
  # par(mar = c(2,2,2,2),
  par(mar = c(7,12,7,14), # c(bottom, left, top, right)
      cex = 2.2, # adjust plotting text size (not title or legend)
      xpd = TRUE)  # allow plotting outside window
  
  branch_lwd_set <- 8
  cell_dend_color_leave_edges <- hclust_dend %>% 
    assign_values_to_leaves_edgePar(hclust_dend, 
                                    value = informative_labels_df_ordered$color, 
                                    edgePar = "col") %>% # this function replaces colLab 
    set("branches_lwd", branch_lwd_set) %>%
    set("labels", informative_labels_df_ordered$new_lbs) 
  
  if (rotate_dendrogram){
    message("rotating dendro")
    cell_dend_color_leave_edges <- cell_dend_color_leave_edges %>% 
      rotate(labels_char_vec) # rotate the dendrogram according to labels_char_vec
  }
  
  # border colors; 1 = black, 2 = red, 3 = green, 4 = light blue, 5 = cyan, 6 =magenta, 7 = gold/yellow, 8 = gray, 9 = repeats back to black
  # lty : 1 = solid, 2 = short dash, 3 = dot 4 = dot and short dash alt, 5 = long dash, 6= tight short dash / dot, 7 = repeats back to solid
  rect_lty_set <- 2 
  border_color_set <- 8 
  rect_lwd_set <- branch_lwd_set / 1.5
  plot(cell_dend_color_leave_edges, 
       main = plot_title, 
       bty='L', xlim = c(2,0), # to help standardize x-axis for comparison
       cex.main = 2.5, # main title text size
       horiz = T, # horizontal
       axes = T)   # remove height with axes = F
  
  if (rotate_dendrogram){
    message("rotating rect boxes")
    rect.dendrogram(tree = cell_dend_color_leave_edges %>% rotate(labels_char_vec), # %>% rotate(labels_char_vec), 
                    k = max_n_clusts, # draw boxes around clusters # lower_rect = par("usr")[2L]*-0.25, # draw from right to left
                    border = border_color_set, horiz = T, # horizontal plot
                    lwd = rect_lwd_set, lty = rect_lty_set,
                    lower_rect = dend_height*-0.1718011 - 0.009845852,
                    # lower_rect = dend_height*-0.156641 + pdf_height*0.000286 - 0.020006,
                    xpd = NA)  # take the bottom border and increse slightly
    
  } else{
    rect.dendrogram(tree = cell_dend_color_leave_edges, # %>% rotate(labels_char_vec), 
                    k = max_n_clusts, # draw boxes around clusters # lower_rect = par("usr")[2L]*-0.25, # draw from right to left
                    border = border_color_set, horiz = T, # horizontal plot
                    lwd = rect_lwd_set, lty = rect_lty_set,
                    lower_rect = dend_height*-0.1718011 - 0.009845852,
                    # lower_rect = dend_height*-0.156641 + pdf_height*0.000286 - 0.020006,
                    xpd = NA)  # take the bottom border and increse slightly
  }
  
  # draw colored bars below dend
  if (rotate_dendrogram){
    message("rotating color bar")
    new_color_order <- tibble(new_lbs = labels(cell_dend_color_leave_edges), 
                              order = 1:length(labels(cell_dend_color_leave_edges)))
    color_bars <- suppressMessages(left_join(vasc_top_labels_df_rotated_order %>% 
                              dplyr::select(new_lbs, color, vasc_grp_color), new_color_order ) %>%
      arrange(order) %>% 
      dplyr::select(vasc_grp_color, color))
    colored_bars(dend = cell_dend_color_leave_edges %>% rotate(labels_char_vec), 
                 horiz = T,
                 colors = color_bars,
                 sort_by_labels_order = FALSE, add = T, rowLabels = "Cell group\nidentity")
  } else {
    color_bars <- informative_labels_df_ordered %>% 
      dplyr::select(vasc_grp_color, color)
    colored_bars(dend = cell_dend_color_leave_edges, horiz = T,
                 colors = color_bars,
                 sort_by_labels_order = FALSE, add = T, rowLabels = "Cell group\nidentity")
  }
  
  
  num_grps <- length(unique(informative_labels_df$grp)) + 1 # plus 1 for the 'non-vasc' group, essential unique sum of both columns grp, and grp_l
  
  # stop()
  non_vasc_color <- unique(informative_labels_df_ordered %>% filter(grp != "Vascular cells") %>% .$vasc_grp_color)
  vasc_color <- unique(informative_labels_df_ordered %>% filter(grp == "Vascular cells") %>% .$vasc_grp_color)
  
  # "topright"
  # x = max(dend_height) + 0.25, y = 9.55
  leg_temp <- bind_rows(informative_labels_df_ordered %>% distinct(grp_l, color, vasc_grp_color) %>% rename(grp = grp_l),
                        informative_labels_df_ordered %>% distinct(grp_l2, color, vasc_grp_color) %>% rename(grp = grp_l2),
                        informative_labels_df_ordered %>% distinct(grp, color, vasc_grp_color) ) %>%
    distinct(grp, color, vasc_grp_color) %>%
    mutate(grp = factor(grp, levels = c("HUVEC", "HAoSMC", "Cancer cells", "Other cells", "Vascular cells", "Non-vascular cells"))) %>%
    arrange(grp)
  leg_ <- bind_rows(leg_temp %>% slice(1:4) %>% dplyr::select(grp, color), 
                    leg_temp %>% slice(5:8) %>% distinct(grp, vasc_grp_color) %>% rename(color = vasc_grp_color))
  
  leg_final <- c(as.character(leg_$grp), "Cluster assignment")
  col_final <- c(leg_$color, border_color_set)
  
  # stop()
  
  xpos <- dend_height + 0.95
  ypos <- ifelse(rotate_dendrogram, xpos*1.5, xpos + 7.65)
  legend( x = xpos, y = ypos , title = "Legend", text.font = 3,
          legend = leg_final,
          col = col_final, 
          lty = c(rep(1,length(leg_final)-1), rect_lty_set), # line style 
          lwd = c(rep(8, length(leg_final)-1), 2.5), # line width
          cex = 2.1, inset = c(INSET_MARGIN)+1.2,
          yjust = 1.2) # legend size
  dev.off()
  
  # stop()
}

# i suspect that the dendrogram height scales linearly with the lower_rect_border
# fit a linear model to see what proportion to adjust
# pdf_height_vec <- seq(4, 20, length.out = 6)

# did this for a given pdf height (8)
# checked that these values look nice
# dend_height_vec <- c(1.385926, 1.439883, 1.560293, 1.35678, 1.241355, 0.8397211)
# lower_rect_border_vec <- c( -0.235, -0.25, -0.25, -0.23, -0.22, -0.14 )
# y <- lm(lower_rect_border_vec ~ dend_height_vec)
# b <- coef(y)[1]; b = -0.009845852 
# m <- coef(y)[2]; m # so lower rect border changes -0.1618011 per unit of dendrogram height
# nice fucking work man, this is sick
# p100_base_output_dir <- "~/Downloads"







# p100_data_lst <- p100_lst_obj; p100_analytes <- p100_data_lst$feature_set
# gcp_data_lst <- gcp_lst_obj; gcp_analytes <- gcp_data_lst$feature_set
# 
# 
# pert_interest_char <- "dmso"
# #p100
# i <- which(p100_res$pert$pert_iname %in% pert_interest_char)
# p100_pert_data_dmso <- p100_res$pert[i,] # choose row
# #gcp
# i <- which(gcp_res$pert$pert_iname %in% pert_interest_char)
# gcp_pert_data_dmso <- gcp_res$pert[i,] # choose row
# 
# 
# moa_interest_char <- "Kinase inhibitor"
# i <- which(p100_res$moa$pert_class %in% moa_interest_char)
# p100_moa_data_ki <- p100_res$moa[i,] # choose row
# 
# i <- which(gcp_res$moa$pert_class %in% moa_interest_char)
# gcp_moa_data_ki <- gcp_res$moa[i, ]
# 
# 
# moa_interest_char <- "Epigenetic" 
# i <- which(p100_res$moa$pert_class %in% moa_interest_char)
# p100_moa_data_epi <- p100_res$moa[i,] # choose row
# 
# i <- which(gcp_res$moa$pert_class %in% moa_interest_char)
# gcp_moa_data_epi <- gcp_res$moa[i, ]

# plot_pretty_dendrogram(data = p100_pert_data_dmso, base_output_dir = p100_base_output_dir,
#                        rotate_dendrogram = TRUE, INSET_MARGIN = -0.06)
# plot_pretty_dendrogram(data = gcp_pert_data_dmso, base_output_dir = gcp_base_output_dir,
#                        rotate_dendrogram = TRUE, INSET_MARGIN = -0.06)
# 
# plot_pretty_dendrogram(data = p100_moa_data_ki, base_output_dir = p100_base_output_dir, 
#                        rotate_dendrogram = TRUE, INSET_MARGIN = -0.06)
# plot_pretty_dendrogram(data = p100_moa_data_epi, base_output_dir = p100_base_output_dir, 
#                        rotate_dendrogram = TRUE, INSET_MARGIN = -0.1)
# 
# plot_pretty_dendrogram(data = gcp_moa_data_ki, base_output_dir = gcp_base_output_dir, 
#                        rotate_dendrogram = TRUE, INSET_MARGIN = -0.15)
# plot_pretty_dendrogram(data = gcp_moa_data_epi, base_output_dir = gcp_base_output_dir, 
#                        rotate_dendrogram = TRUE, INSET_MARGIN = -0.15)






























