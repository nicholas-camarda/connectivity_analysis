
# load("debug-p100-kinase-heat.RData")
# source(file.path("scripts", "init.R"))

analytes_ <- rownames(filtered_test_mat)
u_cell_id_ <- colnames(filtered_test_mat)
cluster_ids <- column_annots_df %>%
  dplyr::distinct(cluster, cell_id, pert_iname, cluster_name, grp_fac) ; 
# %>%
#   separate(col = u_cell_id, into = c("cell_id", "pert_iname", "pert_class_and_batch"), sep = "--") %>%
#   separate(col = pert_class_and_batch, into = c("pert_class", "replicate", "plate_id"), sep = "::") %>%
#   distinct(cluster, cell_id, pert_iname, cluster_name, grp_fac, plate_id); cluster_ids

reduced_long <- filtered_test_mat %>% 
  as_tibble() %>%
  rownames_to_column(var = "temp") %>%
  bind_cols(tibble(analytes = analytes_)) %>%
  dplyr::select(analytes, everything(), -temp) %>%
  pivot_longer(all_of(u_cell_id_), names_to = "u_cell_id") %>%
  separate(col = all_of(u_cell_id), into = c("cell_id", "pert_iname", "pert_class_and_batch"), sep = "--") %>%
  separate(col = pert_class_and_batch, into = c("pert_class", "replicate", "plate_id"), sep = "::") %>%
  left_join(cluster_ids, by = c("cell_id", "pert_iname"))

reduced_tbl <- reduced_long %>%
  group_by(analytes, cell_id, pert_iname) %>%
  summarize(median_value = median(value, na.rm = T), .groups = "keep") %>%
  ungroup() %>%
  mutate(u_cell_id = str_c(cell_id, pert_iname, sep = "::")) %>%
  pivot_wider(id_cols = analytes, names_from = u_cell_id, values_from = median_value)

reduced_mat <- reduced_tbl %>% dplyr::select(-analytes) %>% as.matrix()
rownames(reduced_mat) <- analytes_

n_max_clust <- length(unique(cluster_ids$cluster))
n_max_clust


NAME <- "Log2(value)"
# set color breaks at quantiles
breaks <- sapply(rev(c(.05, .10, .30, .50, .70, .90, .95)), function(q) quantile(reduced_mat, q, na.rm = T, names = F))
col_fun <- colorRamp2(breaks = breaks, colors = brewer.pal(n = length(breaks), name = "RdYlBu"))


ht_opt$COLUMN_ANNO_PADDING <- unit(2, "mm")
ht_opt$ROW_ANNO_PADDING <- unit(1, "mm")
ht_opt$annotation_border <- TRUE


# in the order of the grouping -> we need to adjust colors depending on group content

# stop()
informative_labels_lst <- create_informative_labels()
ildf <- informative_labels_lst$informative_labels_df %>%
  mutate(
    grp = map_chr(grp, function(x) str_split(x, pattern = " ", simplify = T)[, 1]),
    grp_l = map_chr(grp_l, function(x) str_split(x, pattern = " ", simplify = T)[, 1]),
    grp_l2 = map_chr(grp_l2, function(x) str_split(x, pattern = " ", simplify = T)[, 1])
  ) %>%
  left_join(cluster_ids %>% 
              dplyr::distinct(cell_id, cluster, cluster_name), 
            by = c("lbs" = "cell_id")) 

non_vasc_color <- ildf %>% 
  filter(grp_l == "Non-vascular") %>%
  distinct(vasc_grp_color) %>% 
  .$vasc_grp_color

dom_colors <- ildf %>% 
  filter(lbs %in% c("HUVEC", "HAoSMC")) %>%
  dplyr::select(dom = lbs, dom_color = cell_individual_color) %>%
  bind_rows(tibble(dom = "Non-vascular", dom_color = non_vasc_color))

dominant_colors <- ildf %>%
  distinct(cluster_name, lbs, cluster, cell_individual_color) %>%
  group_by(cluster) %>%
  summarize(
    which_in_clust = str_c(sort(unique(lbs)), collapse = ","),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  arrange(cluster) %>%
  mutate(clust_vec = map(which_in_clust, .f = function(x) str_c(str_split(string = x, pattern = ",", simplify = T)[1,]))) %>%
  mutate(dom = map(clust_vec, .f = function(x) {
    if ("HUVEC" %in% x){
      return("HUVEC")
    } else if ("HAoSMC" %in% x) {
      return("HAoSMC")
    } else {
      return("Non-vascular")
    }
  })) %>%
  unnest(cols = "dom") %>%
  distinct() %>%
  # mutate(dom = ifelse(is.na(dom), "none", dom)) %>%
  left_join(dom_colors, by= "dom")

vasc_cluster_id <- cluster_ids %>%
  filter(cell_id %in% c("HUVEC", "HAoSMC")) %>%
  distinct(cell_id, cluster) %>%
  arrange(cell_id)

# you know which cell line is the "base" by looking at the last row
# this must be a factor so we get order!
base_cell <- tail(cluster_ids, n = 1)$grp_fac
grp_factor_lvls <- levels(base_cell)

grp_colors <- cluster_ids %>%
  distinct(cell_id, cluster) %>%
  mutate(dom = ifelse(cluster %in% vasc_cluster_id$cluster[1], "HAoSMC",
                      ifelse(cluster %in% vasc_cluster_id$cluster[2], "HUVEC", "Non-vascular"))) %>%
  left_join(dominant_colors %>% dplyr::select(-cluster), by="dom") %>%
  na.omit() %>%
  distinct(cell_id, cluster, dom, dom_color) %>%
  arrange(cluster) %>%
  dplyr::select(everything(), color = dom_color) %>%
  left_join(ildf %>% distinct(lbs, cell_individual_color), by = c("cell_id" = "lbs")) %>%
  mutate(dom = factor(dom, levels = grp_factor_lvls))



#' @Note:
#' [1] grp_color_df defines the column ordering!!!

# generate the clustering color palette, and establish the column order in grp_color_df
grp_color_df <- left_join(cluster_ids,
                          grp_colors %>% distinct(cluster, dom, color), by = "cluster") %>%
  # in order of the column_df
  mutate(color = ifelse(is.na(color), "#b6b6b6", color)) %>%
  arrange(dom, cell_id, pert_iname) %>%
  mutate(unique_id = make.unique(str_c(cell_id, pert_iname, sep = "::")))
# grayish color for NA
col_fun_annot1 <- grp_color_df$color # 1st group being closest to the labels on the right

group_order <- levels(grp_color_df$dom)
grp_names_char_vec <- grp_color_df$dom
names(col_fun_annot1) <- grp_names_char_vec
block_col_fun_annot1 <- unique(col_fun_annot1)



# reorder the matrix with this -- this could also be done in the heatmap function;
# however, be careful about the variable you use to assign column order and column attributes!
# ESPECIALLY, column split
dim(reduced_mat)
length(analytes_)
length(grp_color_df$cell_id)

length(intersect(colnames(reduced_mat), grp_color_df$unique_id)) == length(grp_color_df$unique_id)

reordered_mat_temp <- reduced_mat[analytes_, grp_color_df$unique_id, drop = F]
reordered_mat <- apply(reordered_mat_temp, 2, FUN = as.numeric)

# perturbations, in the column order
perturbations_char_vec_temp <- str_split(string = grp_color_df$pert_iname, pattern = "_", simplify = TRUE)[, 1]

extra_pert_info <- read_excel(file.path(references_directory, "Drug Glossary_edited.xlsx")) %>%
  mutate(pert_iname = tolower(Drug), moa_simple = `MOA simplified`) %>%
  dplyr::select(pert_iname, moa_simple) %>%
  mutate(moa_simple = ifelse(is.na(moa_simple), "", moa_simple))

link_moa_to_pert <- tibble(pert_iname = perturbations_char_vec_temp) %>%
  left_join(extra_pert_info, by = "pert_iname") %>%
  rename(pert_iname_temp = pert_iname) %>%
  mutate(pert_iname = str_c(pert_iname_temp, " (", moa_simple, ")", sep = ""))

perturbations_char_vec <- link_moa_to_pert$pert_iname
n_u_perts <- length(unique(perturbations_char_vec))
n_u_perts

# generate pert color palette
col_palette_2_df <- generate_color_palette()
col_fun_annot2_df <- inner_join(col_palette_2_df,
                                link_moa_to_pert %>% distinct(pert_iname_temp, pert_iname),
                                by = c("pert_iname" = "pert_iname_temp")
) %>%
  rename(full_pert_name = pert_iname.y)

col_fun_annot2 <- col_fun_annot2_df$colors
names(col_fun_annot2) <- col_fun_annot2_df$full_pert_name

# generate the cell id color palette
cell_color_df <- grp_color_df %>%
  dplyr::select(cell_id) %>%
  left_join(grp_colors, by = "cell_id") %>%
  mutate(color = ifelse(is.na(color), "#b6b6b6", color)) # grayish color for NA
col_fun_annot3 <- cell_color_df$cell_individual_color # 1st group being closest to the labels on the right

cell_names_char_vec <- cell_color_df$cell_id
names(col_fun_annot3) <- cell_names_char_vec
# block_col_fun_annot3 <- unique(col_fun_annot3)

# top annotations
top_ha <- HeatmapAnnotation(
  `Cluster` = anno_block(gp = gpar(
    fill = block_col_fun_annot1,
    labels = group_order,
    labels_gp = gpar(col = "black", font = 2) #, fontsize = FONTSIZE + 2)
  )),
  `Cell line` = cell_names_char_vec,
  `Perturbations` = perturbations_char_vec,
  col = list(
    `Cluster` = col_fun_annot1,
    `Perturbations` = col_fun_annot2,
    `Cell line` = col_fun_annot3
  ),
  
  # gp = gpar(col = "black"), # this draws black lines in the perturbation annotation block
  annotation_name_gp = gpar(font = 2), #, fontsize = FONTSIZE + 2),
  gap = unit(1, "mm"),
  simple_anno_size = unit(0.5, "cm"),
  which = "column",
  # where do the names of the bars go?
  annotation_name_side = "left",
  
  # show_legend = c(TRUE),
  
  annotation_legend_param = list(`Perturbations` = list(ncol = 1))
)


# get top up and down
# needs to be relative to 1
top_up <- row_annots_df %>%
  filter(fc >= 1) %>%
  arrange(desc(`fc`))
top_up
top_down <- row_annots_df %>%
  filter(fc < 1) %>%
  # filter(logFC < LOGFC_CUTOFF) %>%
  arrange(desc(`fc`))
top_down
row_order_df <- bind_rows(top_up, top_down)
row_order_df

stopifnot(nrow(row_order_df) == nrow(filtered_test_mat))

# annotations for right side (by row)
signif_log_q_for_bar_plot <- row_order_df$`-log10(q)`
# signif_log_q <- row_order_df %>% mutate(`-log10(q)` = round(-log(p_val_bh, base = 10),5) ) %>% .$`-log10(q)`
signif_q_val <- round(row_order_df$p_val_bh, 5)
extra_analyte_info <- row_order_df$mark
d_statistics_vec <- round(row_order_df$d_stat_val, 3)
# fold change
fc_vec <- round(2^row_order_df$logFC, 3) # logFC is log_2
analytes_ordered <- row_order_df$analyte
color_fc_row_labels <- tibble(fc_vec) %>%
  mutate(color = ifelse(fc_vec < 1, "red", "black")) %>%
  .$color

right_ha <- HeatmapAnnotation(
  `d_statistic` = anno_text(d_statistics_vec,
                            just = "center", location = 0.5,
                            gp = gpar(
                              # fontsize = FONTSIZE,
                              border = "black"
                            ),
                            width = max_text_width(d_statistics_vec) * 1.075
  ),
  `FoldChange` = anno_text(fc_vec,
                           just = "center", location = 0.5, # log_fc_vec
                           gp = gpar(
                             # fontsize = FONTSIZE,
                             col = color_fc_row_labels,
                             border = "black"
                           ),
                           width = max_text_width(fc_vec) * 1.075
  ),
  `q-value` = anno_text(signif_q_val,
                        just = "center", location = 0.5,
                        gp = gpar(border = "black"), #, fontsize = FONTSIZE),
                        width = max_text_width(signif_q_val) * 1.075
  ), #  `-log10(q)`, anno_barplot(signif_log_q_for_bar_plot, extend = 0.25), # signif_log_q
  
  `site` = anno_text(extra_analyte_info,
                     just = "center", location = 0.5, # unit(1, "npc")
                     gp = gpar(border = "black"), #,fontsize = FONTSIZE),
                     width = max_text_width(extra_analyte_info) * 1.075
  ),
  `analyte` = anno_text(analytes_ordered,
                        just = "left", location = 0.1, # unit(1, "npc")
                        gp = gpar(font = 4), #, fontsize = FONTSIZE + 2),
                        width = max_text_width(extra_analyte_info) * 1.075
  ),
  which = "row",
  gap = unit(1.5, "mm") # gap between annotations
)

ht <- Heatmap(reordered_mat,
              name = NAME,
              col = col_fun,
              na_col = "gray",
              top_annotation = top_ha,
              right_annotation = right_ha,
              
              # row_title = "Analyte",
              row_title_side = "right",
              row_title_rot = 0,
              show_row_names = TRUE,
              # row_names_gp = gpar(fontsize = FONTSIZE + 2, font = 4),
              # column_names_gp = gpar(fontsize = FONTSIZE),
              
              cluster_rows = FALSE,
              # row_order = row_order_df$analyte,
              # cluster_rows = analyte_dend_colored, # logical or pass in hclust res from pvclust
              # row_split = max(analyte_cut_trees),
              
              cluster_columns = FALSE,
              column_names_rot = 90,
              show_column_names = FALSE,
              # column_labels = column_order_df$u_cell_id,
              # column_split = grp_color_df %>% arrange(grp_fac) %>% .$grp_fac,
              column_split = grp_color_df$dom, # split the columns by cluster! need an assignment per column
              
              column_title = " ", # no column title, otherwise defaults to groups
              # column_title_gp = gpar(fill = c("white"), font = 3), # if split column title, then fill can take on mulitple vals
              # column_title = clusters$group,
              # column_order = column_order_df$u_cell_id, # order the columns of this matrix with a new order, by u_cell_id (which matches matrix column names)
              # cluster_columns = col_dend_colored,
              # column_dend_reorder = TRUE,
              
              # column_title_gp = gpar(fontsize = 15, fontface = "bold"),
              # column_title = "Cell ID (unique)",
              # column_title_side = "bottom",
              
              border = TRUE, 
              heatmap_legend_param = list(direction = "horizontal"),
              width = ncol(reordered_mat)*unit(5, "mm"),
              height = nrow(reordered_mat)*unit(5, "mm"),
              raster_device = "png",
              rect_gp = gpar(col = "white", lwd = 0.7)) # outline heatmap with black line

draw(ht)
#   # set general plot params
#   height = unit(6, "in"),
#   width = unit(10, "in"),
#   # height = unit(10, "mm")*nrow(reordered_mat),
#   # width = unit(4, "mm")*ncol(reordered_mat),
#   # height = unit(5, "mm")*nrow(reordered_mat),
#   # width = unit(5, "mm")*ncol(reordered_mat), # scale the shape of the heatmap
#   raster_device = "png",
#   heatmap_legend_param = list(direction = "horizontal"),
#   rect_gp = gpar(col = "white", lwd = 0.7)
# )

# legend for cell lines

lgd_list <- list(Legend(labels = group_order, title = "Cluster", 
                        legend_gp = gpar(fill = block_col_fun_annot1)))

# calculate size
init_ht <- draw(ht,
                heatmap_legend_side = "top",
                annotation_legend_list = lgd_list,
                annotation_legend_side = "left")
# padding = unit(c(2, 2, 4, 2), "cm")
# )
size <- calc_ht_size(init_ht)

# draw the plot
# draw(init_ht)

## Save the plot
save_plot <- function() {
  # in pixels
  linear_height <- size[2] + 1
  linear_height # 9
  linear_width <- size[1] + 1
  linear_width # 16
  
  # linear_width <- 0.1969*ncol(reordered_mat); linear_width
  # linear_height <- unit(6 + 0.1969*nrow(reordered_mat), "mm"); linear_height
  
  setEPS()
  # par(mar=c(2, 2, 2, 2), oma=c(2,2,1,1))
  postscript(
    file = heatmap_output_fn, horizontal = F,
    height = linear_height, width = linear_width, onefile = F
  )
  
  draw(init_ht)
  
  right_annots_height <- unit(1.01, "npc") + unit(3, "mm")
  #####################
  # add title
  # decorate_heatmap_body(NAME, {
  #   grid.text(TITLE,
  #     x = unit(linear_width * 0.75, "cm"), y = unit(1.2, "npc"),
  #     just = "bottom",
  #     gp = gpar(font = 2) #, fontsize = FONTSIZE + 10)
  #   )
  # })
  
  # add the annotation titles
  decorate_annotation("site", {
    grid.text("Site",
              y = right_annots_height, just = "top",
              gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
    )
  })
  
  # add the annotation titles
  decorate_annotation(expression("q-value"), {
    grid.text("q-value",
              y = right_annots_height, just = "top", #
              gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
    )
  })
  
  
  # add the annotation titles
  decorate_annotation(expression("analyte"), {
    grid.text("Analyte",
              y = right_annots_height, just = "top", #
              gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
    )
  })
  
  
  decorate_annotation("FoldChange", {
    grid.text(expression(delta),
              y = right_annots_height, just = "top",
              gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
    )
  })
  
  decorate_annotation("d_statistic", {
    grid.text("D-stat",
              y = right_annots_height, just = "top",
              gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
    )
  })
  dev.off()
}
save_plot()