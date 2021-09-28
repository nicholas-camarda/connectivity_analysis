# source("~/OneDrive - Tufts/phd/jaffe/workspace/ws/scripts/master-source.R")

# library(heatmaply)
# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
# https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html#a-brief-description-of-following-chapters
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/

#'@note answering these questions
#' [#1] For the DMSO data, what P100 or GCP marks differentiate vascular from cancer cells at baseline?
#' [#1a] P100 marks
#' [#1b] GCP marks
#' [#2] For the P100 phosphoproteomic response, compare vascular cells to cancer cells and figure out which 
#' phosphorylation events differentiate the:
#' [#2a] p100 marks under epigenetic modifying drug moa
#' [#2b] p100 marks under the kinase inhibitor moa

# dataframes must be passed as lists
# args_lst <- list(  data = list(p100_dmso_res_obj), 
#                    analytes = list(p100_dmso_lst_obj$feature_set), 
#                    title_var = "pert_iname", 
#                    heat_map_base_output_dir = "~/Downloads/test_pwalk")
# pwalk(.l = args_lst, 
#       .f = plot_heatmap)


generate_color_palette <- function(){
  #' @note this is from master-source.R
  colors_ <- pal_igv("default")(51)[(1:nrow(all_drugs_mapping))+5]; length(colors_)
  all_drugs_mapping_final <- all_drugs_mapping %>% bind_cols(colors = colors_)
  return(all_drugs_mapping_final)
}



############## HEATMAP FUNCTION ####################
############## HEATMAP FUNCTION ####################

#' @param data dataframe of column-lists, output of *cluster3.R*
#' @param analytes character vector of analytes from original data
#' @param title_var column name of analysis you'd wish to run (e.g. pert_iname, drug_class, or all)
#' @param base_output_dir specific to dataset, e.g. ~/Downloads/p100 
plot_heatmap <- function(args){
  
  #' colors: 
  #' https://www.colorhexa.com/
  # stop()
  dataset <- force_natural(args$which_dat)
  grouping_var <- force_natural(args$grouping_var)
  output_dir_temp <- force_natural(args$path)
  spec_char <- args %>% pluck(1)
  # output_dir <- file.path(output_dir_temp, spec_char)
  # dir.create(output_dir, recursive = T, showWarnings = F)
  
  heatmap_output_fn <- force_natural(args$path)
  # stop()
  # dat_tbl <- force_natural(args$input_data); unique_clust_assignments <- force_natural(args$cluster_lst) %>% pluck(2); diff_ex_df <- force_natural(args$diffe_lst) %>% pluck(1)
  
  # load transposed matrix for condition
  dat_tbl <- args$input_data
  
  w_annot_tbl <- dat_tbl %>% 
    ungroup() %>% 
    dplyr::select(replicate_id, pert_iname, pert_class, master_id, pr_gene_symbol, value) %>%
    pivot_wider(names_from = pr_gene_symbol, 
                values_from = value, 
                values_fn = median) %>%
    ungroup()
  
  rnames <- w_annot_tbl$replicate_id
  analytes <- unique(dat_tbl$pr_gene_symbol)
  
  tbl <- w_annot_tbl %>% dplyr::select(replicate_id, all_of(analytes))
  t_mat <- tbl %>% dplyr::select(-replicate_id) %>% as.data.frame()
  rownames(t_mat) <- rnames
  mat <- t_mat %>% t() %>%as.data.frame()
  
  
  # get clust assignments and create dataframe of column annotations
  unique_clust_assignments <- args$clust_lst %>% pluck(2)
  clusters_by_cell_id_df <- tibble::enframe(x = unique_clust_assignments, name = "cell_id", value = "cluster")
  column_annots_df_temp <- sapply(colnames(mat),  
                                  function(x) str_split(x, "--",simplify = T)[,1]) %>%
    tibble::enframe(name = "u_cell_id", value = "cell_id") %>% 
    left_join(clusters_by_cell_id_df, by="cell_id") %>%
    mutate(pert_iname = w_annot_tbl$pert_iname,
           pert_class = w_annot_tbl$pert_class)
  
  
  # get diff_ex df
  if (class(args$diffe_lst) == "list") {
    diff_ex_df <- args$diffe_lst %>% pluck(1)
  } else { 
    diff_ex_df <- args$diffe_lst
  }
  
  # generate the row annotations
  row_annots_df_include_nonsig_all <- diff_ex_df %>% 
    mutate(`-log10(q)` = -log(p_val_boot_bh, base = 10)) %>%
    # only handles the boot statistic -- make it handle the regular stat
    mutate(d_stat_val = ifelse(directional_stat == "--" & !is.na(ks_boot_statistic), -1*ks_boot_statistic, ks_boot_statistic)) 
  
  # row_order_df <- 
  #'@NOTE: *TO-DO* need to be able to plot multiple heatmaps if 
  #' vasc cells don't cluster together
  # Get co-clust
  # co_clust <-  data$co_clust[[1]]; co_clust
  
  # this helps count the number of extra cells that clustered with the vasculars, or not!
  # clustered together with other cells
  # clustered together and alone
  # clustered separately with no others
  # clustered separately with others
  
  vasc_clust_id_df <- filter(clusters_by_cell_id_df, cell_id %in% vascular_char_vec); vasc_clust_id_df
  
  vasc_cluster_assignments_int <- unique(vasc_clust_id_df %>% .$cluster) # will be a doublet if they are in different clusters
  others_clustered_df <- filter(clusters_by_cell_id_df, 
                                (cluster %in% vasc_cluster_assignments_int)); others_clustered_df
  
  huvec_cluster <- filter(vasc_clust_id_df, cell_id == "HUVEC")$cluster
  haosmc_cluster <- filter(vasc_clust_id_df, cell_id == "HAoSMC")$cluster
  if (huvec_cluster == haosmc_cluster) { 
    co_clust <- T
  } else {
    co_clust <- F
  }
  
  if (nrow(others_clustered_df) > 2) {
    clustered_with_others <- T
  } else{
    clustered_with_others <- F
  }
  
  # stop()
  co_clust; clustered_with_others; vasc_clust_id_df; clusters_by_cell_id_df
  # stop()
  if (!co_clust & (clustered_with_others | !clustered_with_others)) {
    message("Vascular cells clustered separately, alone or with other cells. Separating plots...")
    
    # extract significant analytes per cluster
    row_annots_df1 <- filter(row_annots_df_include_nonsig_all, 
                             str_detect(string = base_clust_comp_name, pattern = "HUVEC"), signif); row_annots_df1
    row_annots_df2 <- filter(row_annots_df_include_nonsig_all, 
                             str_detect(string = base_clust_comp_name, pattern = "HAoSMC"), signif); row_annots_df2
    
    # get the character vector of the cluster name
    bccn1_char_vec <- str_extract(string = unique(row_annots_df1$base_clust_comp_name), pattern = "HUVEC")
    bccn2_char_vec <- str_extract(string = unique(row_annots_df2$base_clust_comp_name), pattern = "HAoSMC")
    
    cell_id_huvec_clust_char_vec <- as.character(unique(str_split(unique(row_annots_df1$base_clust_comp_name), ",", simplify = T))); cell_id_huvec_clust_char_vec
    cell_id_haosmc_clust_char_vec <- as.character(unique(str_split(unique(row_annots_df2$base_clust_comp_name), ",", simplify = T))); cell_id_haosmc_clust_char_vec
    
    column_annots_df1 <- column_annots_df_temp %>%
      mutate(group = ifelse(cell_id %in% vascular_char_vec, cell_id, "Non-vascular"), 
             cluster_name = ifelse(cell_id %in% cell_id_huvec_clust_char_vec, bccn1_char_vec, 
                                   ifelse(cell_id %in% bccn2_char_vec, bccn2_char_vec, "Non-vascular")),
             which_dat = dataset, 
             grp_fac = factor(cluster_name, 
                              levels = c("Non-vascular", bccn2_char_vec, bccn1_char_vec))); column_annots_df1
    
    filtered_test_mat1 <- mat[rownames(mat) %in% row_annots_df1$analyte, 
                              colnames(mat) %in% column_annots_df1$u_cell_id, 
                              drop=F ]
    
    
    
    column_annots_df2 <- column_annots_df_temp %>%
      mutate(group = ifelse(cell_id %in% vascular_char_vec, cell_id, "Non-vascular"), 
             cluster_name = ifelse(cell_id %in% cell_id_huvec_clust_char_vec, bccn1_char_vec, 
                                   ifelse(cell_id %in% bccn2_char_vec, bccn2_char_vec, "Non-vascular")),
             which_dat = dataset, 
             grp_fac = factor(cluster_name, levels = c("Non-vascular", bccn1_char_vec, bccn2_char_vec))); column_annots_df2
    
    filtered_test_mat2 <- mat[rownames(mat) %in% row_annots_df2$analyte, colnames(mat) %in% column_annots_df2$u_cell_id, drop=F ] 
    
    split_fn <- str_split(heatmap_output_fn, "\\.", simplify = T)
    heatmap_output_fn1 <- str_c(split_fn[,1], qq("-@{bccn1_char_vec}."), split_fn[,2], collapse = "")
    heatmap_output_fn2 <- str_c(split_fn[,1], qq("-@{bccn2_char_vec}."), split_fn[,2], collapse = "")
    
    # stop()
    ht1_lst <- organize_and_plot_heatmap_subfunction(filtered_test_mat = filtered_test_mat1, 
                                                     row_annots_df = row_annots_df1, column_annots_df = column_annots_df1,
                                                     heatmap_output_fn = heatmap_output_fn1, title_var = grouping_var)
    
    ht2_lst <- organize_and_plot_heatmap_subfunction(filtered_test_mat = filtered_test_mat2, 
                                                     row_annots_df = row_annots_df2, column_annots_df = column_annots_df2,
                                                     heatmap_output_fn = heatmap_output_fn2, title_var = grouping_var)
    
    return(list("HUVEC" = list(bare_ht = ht1_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all, 
                               filtered_row_annot = ht1_lst$row_order_df, col_annot = ht1_lst$column_order_df, 
                               # transformed_mat = row_minmax_filtered_test_mat,
                               plotted_matrix = ht2_lst$reordered_mat), 
                "HAoSMC" = list(bare_ht = ht2_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all, 
                                filtered_row_annot = ht2_lst$row_order_df, col_annot = ht2_lst$column_order_df, 
                                # transformed_mat = row_minmax_filtered_test_mat,
                                plotted_matrix = ht2_lst$reordered_mat)))
    
  } else if (co_clust & clustered_with_others) {
    message("Vascular cells clustered together, but with other cells. Plotting once...")
    
    #' @NOTE: 
    #' in here need to separate out HUVECs just for labeling purposes
    row_annots_df <- filter(row_annots_df_include_nonsig_all, 
                            str_detect(string = base_clust_comp_name, pattern = "HUVEC|HAoSMC"), signif); row_annots_df
    
    choose_analytes <- rownames(mat) %in% row_annots_df$analyte
    filtered_test_mat <- mat[choose_analytes,,drop=F ] 
    
    clust_name_char <- unique(row_annots_df$base_clust_comp_name); clust_name_char
    clust_name_char_split_vec <- unlist(map(clust_name_char, function(c) str_split(c, ",", simplify = T))); clust_name_char_split_vec
    
    column_annots_df <- column_annots_df_temp %>%
      mutate(group = ifelse(cell_id %in% clust_name_char_split_vec, 
                            clust_name_char,"Non-vascular"),
             cluster_name = group,
             which_dat = dataset) %>%
      mutate(grp_fac = factor(cluster_name, levels = c("Non-vascular", clust_name_char))); column_annots_df
    
    ht_lst <- organize_and_plot_heatmap_subfunction(filtered_test_mat = filtered_test_mat, 
                                                    row_annots_df = row_annots_df, column_annots_df = column_annots_df,
                                                    heatmap_output_fn = heatmap_output_fn, title_var = grouping_var)
    
    return("Combined_Vascular" = list(bare_ht = ht_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all, 
                                      filtered_row_annot = ht_lst$row_order_df, col_annot = ht_lst$column_order_df, 
                                      # transformed_mat = row_minmax_filtered_test_mat,
                                      plotted_matrix = ht_lst$reordered_mat))
    
  } else {
    message("Vascular cells clustered together with no other cells. Plotting once...")
    row_annots_df <- filter(row_annots_df_include_nonsig_all, 
                            str_detect(string = base_clust_comp_name, pattern = "HUVEC|HAoSMC"), signif); row_annots_df
    
    filtered_test_mat <- mat[rownames(mat) %in% row_annots_df$analyte,,drop=F ]
    
    column_annots_df <- column_annots_df_temp %>%
      mutate(group = ifelse(cell_id %in% c("HUVEC","HAoSMC"), 
                            "Vascular","Non-vascular"),
             cluster_name = group,
             which_dat = dataset) %>%
      mutate(grp_fac = factor(cluster_name, levels = c("Non-vascular", "Vascular"))); column_annots_df
    
    ht_lst <- organize_and_plot_heatmap_subfunction(filtered_test_mat = filtered_test_mat, 
                                                    row_annots_df = row_annots_df, column_annots_df = column_annots_df,
                                                    heatmap_output_fn = heatmap_output_fn, title_var = grouping_var)
    
    return("Vascular" = list(bare_ht = ht_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all, 
                             filtered_row_annot = ht_lst$row_order_df, col_annot = ht_lst$column_order_df, 
                             # transformed_mat = row_minmax_filtered_test_mat,
                             plotted_matrix = ht_lst$reordered_mat))
  }
}


calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}


#' @note helper function for plotting heatmaps
organize_and_plot_heatmap_subfunction <- function(filtered_test_mat, 
                                                  row_annots_df, column_annots_df,
                                                  heatmap_output_fn, title_var = "pert_iname"){
  
  column_order_df <- column_annots_df %>% 
    arrange(grp_fac, cluster_name, pert_iname, u_cell_id); column_order_df
  
  # clusters in the order desired for clusters
  group_order <- levels(column_order_df$grp_fac);group_order
  n_max_clust <-  length(group_order); n_max_clust
  
  # perturbations, in the column order
  perturbations_char_vec_temp <- column_order_df$pert_iname; unique(perturbations_char_vec_temp)
  
  extra_pert_info <- read_excel(file.path(references_directory, "Drug Glossary_edited.xlsx")) %>%
    mutate(pert_iname = tolower(Drug), moa_simple = `MOA simplified`) %>%
    dplyr::select(pert_iname, moa_simple) %>%
    mutate(moa_simple = ifelse(is.na(moa_simple), "", moa_simple))
  
  link_moa_to_pert <- tibble(pert_iname = perturbations_char_vec_temp) %>%
    left_join(extra_pert_info, by = "pert_iname") %>%
    rename(pert_iname_temp = pert_iname) %>%
    mutate(pert_iname = str_c(pert_iname_temp, " (", moa_simple, ")", sep = ""))
  
  perturbations_char_vec <- link_moa_to_pert$pert_iname
  n_u_perts <- length(unique(perturbations_char_vec)); n_u_perts
  
  # get top up and down
  # needs to be relative to 1
  top_up <- row_annots_df %>% filter(logFC > 1) %>% arrange(desc(`logFC`)); top_up
  top_down <- row_annots_df %>% filter(logFC < 1) %>% arrange(desc(`logFC`)); top_down
  row_order_df <- bind_rows(top_up,top_down) ; row_order_df
  
  stopifnot(nrow(row_order_df) == nrow(filtered_test_mat))
  
  # annotations for right side (by row)
  signif_log_q_for_bar_plot <-  row_order_df$`-log10(q)`;
  # signif_log_q <- row_order_df %>% mutate(`-log10(q)` = round(-log(p_val_bh, base = 10),5) ) %>% .$`-log10(q)`
  signif_q_val <- round(row_order_df$p_val_bh,5)
  extra_analyte_info <- row_order_df$mark
  d_statistics_vec <- round(row_order_df$d_stat_val, 3)
  fc_vec <- round(2^row_order_df$logFC, 3) # logFC is log_2
  analytes <- row_order_df$analyte
  color_fc_row_labels <- tibble(fc_vec) %>% 
    mutate(color = ifelse(fc_vec < 1, "red", "black"))  %>% 
    .$color
  
  # reorder the matrix first -- this could also be done in the heatmap function; 
  # however, be careful about the variable you use to assign column order and column attributes!
  # ESPECIALLY, column split
  reordered_mat_temp <- filtered_test_mat[analytes, column_order_df$u_cell_id, drop = F] 
  reordered_mat <- apply(reordered_mat_temp, 2, FUN = as.numeric)
  # tail(reordered_mat)
  ## PLOTTING
  FONTSIZE <- 7
  print_friendly_name <- str_replace(string = str_c(str_split(unique(column_order_df[title_var][[1]]), "_", simplify = TRUE), 
                                                    collapse = " "), pattern = "excl", replacement = "excluding")
  TITLE <-  str_c("Differential analytes using connectivity clusters\n", 
                  str_c(print_friendly_name, 
                        toupper(unique(column_order_df$which_dat[[1]])), sep = ' || ')); TITLE
  
  NAME <- "Log2(value)"; NAME
  # NAME <- "r/max(r)"; NAME # this is how i calculate each row's value, as fraction of max of row
  
  ht_opt$COLUMN_ANNO_PADDING <- unit(2, "mm")
  ht_opt$ROW_ANNO_PADDING <- unit(1, "mm")
  ht_opt$annotation_border <- TRUE
  
  # set color breaks at quantiles
  breaks <- sapply(rev(c(.05, .10,.30,.50,.70,.90, .95)), function(q) quantile(reordered_mat, q, na.rm = T, names = F))
  col_fun = colorRamp2(breaks = breaks, colors = brewer.pal(n = length(breaks), name = "RdYlBu"))
  
  # in the order of the grouping -> we need to adjust colors depending on group content
  
  # stop()
  informative_labels_lst <- create_informative_labels()
  ildf <- informative_labels_lst$informative_labels_df %>%
    mutate(grp = map_chr(grp, function(x) str_split(x, pattern = " ", simplify =T)[,1]),
           grp_l = map_chr(grp_l, function(x) str_split(x, pattern = " ", simplify =T)[,1]),
           grp_l2 = map_chr(grp_l2, function(x) str_split(x, pattern = " ", simplify =T)[,1]))
  grp_names <- unique(column_order_df$group)
  
  # you know which cell line is the "base" by looking at the last row 
  # this must be a factor so we get order!
  base_cell <- tail(column_order_df, n = 1)$grp_fac
  
  if (length(grp_names) == 3){
    grp_colors <- bind_rows(ildf %>% dplyr::select(grp_l2, color, lbs, cell_individual_color) %>% filter(grp_l2 %in% c("HAoSMC", "HUVEC")), 
                            ildf %>% dplyr::select(grp_l2, vasc_grp_color,lbs,  cell_individual_color) %>% filter(!(grp_l2 %in% c("HAoSMC", "HUVEC"))) %>% rename(color = vasc_grp_color)) %>%
      distinct() %>%
      mutate(grp_fac = factor(grp_l2, levels = levels(base_cell))) %>%
      arrange(grp_fac); grp_colors
  } else if (length(grp_names) == 2 & "Vascular" %in% grp_names){
    grp_colors <- ildf %>% 
      dplyr::select(grp_l, vasc_grp_color, lbs, cell_individual_color) %>%
      rename(grp_l2 = grp_l, color = vasc_grp_color) %>%
      distinct() %>%
      mutate(grp_fac = factor(grp_l2, levels = levels(base_cell))) %>%
      arrange(grp_fac); grp_colors
  } else {
    # clustered together but with other cells - factors won't match
    grp_colors <- ildf %>% 
      dplyr::select(grp_l, vasc_grp_color,lbs,  cell_individual_color) %>%
      rename(grp_l2 = grp_l, color = vasc_grp_color) %>%
      distinct() %>%
      mutate(grp_fac = factor(grp_l2)) %>%
      arrange(grp_fac); grp_colors
  }
  
  
  #' @TODO:
  #' [1] resolve huvec or haosmc in "all other" -> then the rest becomes non-vascular
  
  grp_color_df <- column_order_df %>% 
    dplyr::select(group) %>%
    rename(grp_l2 = group) %>%
    left_join(grp_colors, by = "grp_l2") %>%
    mutate(color = ifelse(is.na(color), "#b6b6b6", color)) # grayish color for NA
  col_fun_annot1 <- grp_color_df$color  # 1st group being closest to the labels on the right
  
  grp_names_char_vec <- grp_color_df$grp_l2 
  names(col_fun_annot1) <- grp_names_char_vec
  block_col_fun_annot1 <- unique(col_fun_annot1)
  
  col_palette_2_df <- generate_color_palette() 
  col_fun_annot2 <- col_palette_2_df %>% 
    filter(pert_iname %in% unique(link_moa_to_pert$pert_iname_temp)) %>% 
    .$colors
  
  names(col_fun_annot2) <- unique(perturbations_char_vec)
  
  cell_color_df <- column_order_df %>% 
    dplyr::select(cell_id) %>%
    rename(lbs = cell_id) %>%
    left_join(grp_colors, by = "lbs") %>%
    mutate(color = ifelse(is.na(color), "#b6b6b6", color)) # grayish color for NA
  col_fun_annot3 <- cell_color_df$cell_individual_color  # 1st group being closest to the labels on the right
  
  cell_names_char_vec <- cell_color_df$lbs 
  names(col_fun_annot3) <- cell_names_char_vec
  # block_col_fun_annot3 <- unique(col_fun_annot3)
  
  # top annotations
  top_ha <- HeatmapAnnotation(
    `Cluster` = anno_block(gp = gpar(fill = block_col_fun_annot1,
                                     labels = group_order,
                                     labels_gp = gpar(col = "black", font = 2, fontsize = FONTSIZE + 2))),
    `Cell line` = cell_names_char_vec,
    
    `Perturbations` = perturbations_char_vec,
    
    col = list(`Cluster` = col_fun_annot1,
               `Perturbations` = col_fun_annot2,
               `Cell line` = col_fun_annot3),
    
    # gp = gpar(col = "black"), # this draws black lines in the perturbation annotation block
    annotation_name_gp = gpar(fontsize = FONTSIZE + 2, font = 2), 
    gap = unit(1, "mm"),
    
    simple_anno_size = unit(0.5, "cm"),
    
    annotation_name_side = "left",
    which = "column",
    # show_legend = c(TRUE),
    
    annotation_legend_param = list(`Perturbations` = list(ncol = 2))
  )
  
  # stop()
  # right annotations
  
  right_ha <- HeatmapAnnotation(
    `d_statistic` = anno_text(d_statistics_vec, just = "center",location = 0.5,
                              gp = gpar(fontsize = FONTSIZE,
                                        border = "black"),
                              width = max_text_width(d_statistics_vec)*1.075),
    
    `FoldChange` = anno_text(fc_vec, just = "center",location = 0.5, # log_fc_vec
                             gp = gpar(fontsize = FONTSIZE, 
                                       col = color_fc_row_labels,
                                       border = "black"),
                             width = max_text_width(fc_vec)*1.075),
    
    `q-value` = anno_text(signif_q_val, just = "center",location = 0.5,
                          gp = gpar(fontsize = FONTSIZE, border = "black"),
                          width = max_text_width(signif_q_val)*1.075), #  `-log10(q)`, anno_barplot(signif_log_q_for_bar_plot, extend = 0.25), # signif_log_q
    
    `site` = anno_text(extra_analyte_info, just = "center",location = 0.5, # unit(1, "npc")
                       gp = gpar(fontsize = FONTSIZE, border = "black"),
                       width = max_text_width(extra_analyte_info)*1.075),
    
    `analyte` = anno_text(analytes, just = "left", location = 0.1, # unit(1, "npc")
                          gp = gpar(fontsize = FONTSIZE + 2, font = 4),
                          width = max_text_width(extra_analyte_info)*1.075),
    
    
    which = "row",
    gap = unit(1.5, "mm") # gap between annotations
  ) 
  
  #' @To-Do: consider:
  # order the columns by their cluster identity, in addition to non-vascular vs vascular?
  # stop()
  
  ht <- Heatmap(reordered_mat,
                name = NAME, 
                col = col_fun, 
                na_col = "gray",
                
                top_annotation = top_ha, 
                right_annotation = right_ha,
                
                # row_title = "Analyte",
                row_title_side = "right",
                row_title_rot = 0, 
                show_row_names = FALSE,
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
                column_split = column_order_df$grp_fac, # split the columns by cluster! need an assignment per column
                
                column_title = " ", # no column title, otherwise defaults to groups
                # column_title_gp = gpar(fill = c("white"), font = 3), # if split column title, then fill can take on mulitple vals
                # column_title = clusters$group,
                # column_order = column_order_df$u_cell_id, # order the columns of this matrix with a new order, by u_cell_id (which matches matrix column names)
                # cluster_columns = col_dend_colored,
                # column_dend_reorder = TRUE,
                
                # column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                # column_title = "Cell ID (unique)",
                # column_title_side = "bottom",
                
                border = TRUE, # outline heatmap with black line
                
                # set general plot params
                height = unit(6, "in"),
                width = unit(10, "in"),
                # height = unit(10, "mm")*nrow(reordered_mat), 
                # width = unit(4, "mm")*ncol(reordered_mat),
                # height = unit(5, "mm")*nrow(reordered_mat),
                # width = unit(5, "mm")*ncol(reordered_mat), # scale the shape of the heatmap
                raster_device = "png",
                heatmap_legend_param = list(direction = "horizontal"),
                rect_gp = gpar(col = "white", lwd = 0.7));draw(ht)
  # legend for cell lines
  
  lgd_list = list(Legend(labels = group_order, title = "Group", legend_gp = gpar(fill = block_col_fun_annot1)))
  
  size = calc_ht_size(ht)
  
  ## Save the plot
  save_plot <- function(){
    
    linear_height <- size[2] + 1.5; linear_height # 9 
    linear_width <- size[1] + 1.5; linear_width # 16
    
    # linear_width <- 0.1969*ncol(reordered_mat); linear_width
    # linear_height <- unit(6 + 0.1969*nrow(reordered_mat), "mm"); linear_height
    
    setEPS()
    # par(mar=c(2, 2, 2, 2), oma=c(2,2,1,1))
    postscript(file = heatmap_output_fn, horizontal = F,
               height = linear_height, width = linear_width, onefile = F)
    draw(ht,
         heatmap_legend_side = "left", 
         annotation_legend_list = lgd_list,
         annotation_legend_side = "left",
         padding = unit(c(2, 2, 4, 2), "cm"))
    #####################
    # add title
    decorate_heatmap_body(NAME, {
      grid.text(TITLE, x = unit(linear_width*0.75, "cm"), y = unit(1.2, "npc"), 
                just = "bottom",
                gp = gpar(font = 2, fontsize = FONTSIZE + 10))
    })
    
    # add the annotation titles
    decorate_annotation("site", {
      grid.text("Site", y = unit(1, "npc") + unit(3, "mm"), just = "top",
                gp = gpar(font = 2, fontsize = FONTSIZE + 2))
    })
    
    # add the annotation titles
    decorate_annotation(expression("q-value"), {
      grid.text("q-value",  y = unit(1, "npc") + unit(3, "mm"), just = "top",  #
                gp = gpar(font = 2, fontsize = FONTSIZE + 2 ))
    })
    
    
    # add the annotation titles
    decorate_annotation(expression("analyte"), {
      grid.text("Analyte",  y = unit(1, "npc") + unit(3, "mm"), just = "top",  #
                gp = gpar(font = 2, fontsize = FONTSIZE + 2 ))
    })
    
    
    decorate_annotation("FoldChange", {
      grid.text(expression(delta),  y = unit(1, "npc") + unit(3, "mm"), just = "top", 
                gp = gpar(font = 2, fontsize = FONTSIZE + 2 ))
    })
    
    decorate_annotation("d_statistic", {
      grid.text("D-stat",  y = unit(1, "npc") + unit(3, "mm"), just = "top", 
                gp = gpar(font = 2, fontsize = FONTSIZE + 2 ))
    })
    dev.off()
    
  }
  save_plot()
  
  # stop()
  message("Done!")
  return(list(ht, row_order_df, column_order_df, reordered_mat) %>% 
           set_names("ht","row_order_df", "column_order_df", "reordered_mat"))
}


############## HEATMAP FUNCTION ####################
############## HEATMAP FUNCTION ####################





