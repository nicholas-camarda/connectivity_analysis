# source("~/OneDrive - Tufts/phd/jaffe/workspace/ws/scripts/master-source.R")

# library(heatmaply)
# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
# https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html#a-brief-description-of-following-chapters
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/

#' @note answering these questions
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
#' @note generate color palette for heatmaps
#' @note this needs to be here because all_drugs_mapping hasn't been made yet
generate_color_palette <- function() {
  #' @note this is from master-source.R
  # need divergent colors
  # colors_ <- pal_igv("default")(51)[(1:nrow(all_drugs_mapping)) + 5]
  # set.seed(1234)
  # colors_ <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  # sampled_colors <- colors_[sample(x = 1:nrow(all_drugs_mapping), 
  #                                  size = nrow(all_drugs_mapping), replace = F)]
  n <- nrow(all_drugs_mapping)
  palette <- distinctColorPalette(n)
  all_drugs_mapping_final <- all_drugs_mapping %>% 
    bind_cols(colors = palette)
  return(all_drugs_mapping_final)
}

############## HEATMAP FUNCTION ####################
############## HEATMAP FUNCTION ####################

#' @param data dataframe of column-lists, output of *cluster3.R*
#' @param analytes character vector of analytes from original data
#' @param title_var column name of analysis you'd wish to run (e.g. pert_iname, drug_class, or all)
#' @param base_output_dir specific to dataset, e.g. ~/Downloads/p100
plot_heatmap <- function(args) {
  
  #' colors:
  #' https://www.colorhexa.com/
  # stop()
  
  dataset <- force_natural(unique(args$which_dat))
  grouping_var <- force_natural(unique(args$grouping_var))
  spec_char <- force_natural(args$dirname_)
  heatmap_output_fn <- force_natural(args$path)
  
  # load transposed matrix for condition
  dat_tbl <- unique(args$input_data)
  
  w_annot_tbl <- dat_tbl %>%
    ungroup() %>%
    dplyr::select(replicate_id, pert_iname, pert_class, master_id, pr_gene_symbol, value) %>%
    pivot_wider(
      names_from = pr_gene_symbol,
      values_from = value,
      values_fn = median
    ) %>%
    ungroup()
  
  rnames <- w_annot_tbl$replicate_id
  analytes <- unique(dat_tbl$pr_gene_symbol)
  
  tbl <- w_annot_tbl %>% dplyr::select(replicate_id, all_of(analytes))
  t_mat <- tbl %>%
    dplyr::select(-replicate_id) %>%
    as.data.frame()
  rownames(t_mat) <- rnames
  mat <- t_mat %>%
    t() %>%
    as.data.frame()
  
  
  
  # get clust assignments and create dataframe of column annotations
  unique_clust_assignments <- args$clust_lst %>%
    pluck(2)
  
  clusters_by_cell_id_df <- tibble::enframe(x = unique_clust_assignments, 
                                            name = "cell_id", 
                                            value = "cluster")
  column_annots_df_temp <- sapply(
    colnames(mat),
    function(x) str_split(x, "--", simplify = T)[, 1]
  ) %>%
    tibble::enframe(name = "u_cell_id", value = "cell_id") %>%
    left_join(clusters_by_cell_id_df, by = "cell_id") %>%
    mutate(
      pert_iname = w_annot_tbl$pert_iname,
      pert_class = w_annot_tbl$pert_class
    )
  
  
  # get diff_ex df
  diff_ex_df <- args$diffe_lst %>%
    pluck(1)
  
  # generate the row annotations
  row_annots_df_include_nonsig_all <- diff_ex_df %>%
    mutate(`-log10(q)` = -log(p_val_boot_bh, base = 10)) %>%
    # only handles the boot statistic -- make it handle the regular stat
    mutate(d_stat_val = ifelse(directional_stat == "--" & !is.na(ks_boot_statistic), 
                               -1 * ks_boot_statistic, 
                               ks_boot_statistic))
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-p100-heat-pre.RData")
  # load("debug/debug_dat/debug-p100-heat-pre.RData")
  # stop()
  
  
  # row_order_df <-
  #' @NOTE: *TO-DO* need to be able to plot multiple heatmaps if
  #' vasc cells don't cluster together
  # Get co-clust
  # co_clust <-  data$co_clust[[1]]; co_clust
  
  # this helps count the number of extra cells that clustered with the vasculars, or not!
  # clustered together with other cells
  # clustered together and alone
  # clustered separately with no others
  # clustered separately with others
  
  vasc_clust_id_df <- filter(clusters_by_cell_id_df, cell_id %in% vascular_char_vec)
  vasc_clust_id_df
  
  vasc_cluster_assignments_int <- unique(vasc_clust_id_df %>% .$cluster) # will be a doublet if they are in different clusters
  others_clustered_df <- filter(
    clusters_by_cell_id_df,
    (cluster %in% vasc_cluster_assignments_int)
  )
  others_clustered_df
  
  huvec_cluster <- filter(vasc_clust_id_df, cell_id == "HUVEC")$cluster
  haosmc_cluster <- filter(vasc_clust_id_df, cell_id == "HAoSMC")$cluster
  if (huvec_cluster == haosmc_cluster) {
    co_clust <- T
  } else {
    co_clust <- F
  }
  
  if (nrow(others_clustered_df) > 2) {
    clustered_with_others <- T
  } else {
    clustered_with_others <- F
  }
  
  # stop()
  # co_clust
  # clustered_with_others
  # vasc_clust_id_df
  # clusters_by_cell_id_df
  
  
  if (!co_clust & (clustered_with_others | !clustered_with_others)) {
    message("Vascular cells clustered separately, alone or with other cells. Separating plots...")
    
    # extract significant analytes per cluster
    row_annots_df1 <- filter(
      row_annots_df_include_nonsig_all,
      str_detect(string = base_clust_comp_name, pattern = "HUVEC"), signif
    ); row_annots_df1
    row_annots_df2 <- filter(
      row_annots_df_include_nonsig_all,
      str_detect(string = base_clust_comp_name, pattern = "HAoSMC"), signif
    ); row_annots_df2
    
    # get the character vector of the cluster name
    bccn1_char_vec <- str_extract(string = unique(row_annots_df1$base_clust_comp_name),
                                  pattern = "HUVEC")
    bccn2_char_vec <- str_extract(string = unique(row_annots_df2$base_clust_comp_name), 
                                  pattern = "HAoSMC")
    
    cell_id_huvec_clust_char_vec <- as.character(unique(str_split(unique(row_annots_df1$base_clust_comp_name), ",", simplify = T)))
    cell_id_huvec_clust_char_vec
    cell_id_haosmc_clust_char_vec <- as.character(unique(str_split(unique(row_annots_df2$base_clust_comp_name), ",", simplify = T)))
    cell_id_haosmc_clust_char_vec
    
    column_annots_df1 <- column_annots_df_temp %>%
      mutate(
        group = ifelse(cell_id %in% vascular_char_vec, cell_id, "Non-vascular"),
        cluster_name = ifelse(cell_id %in% cell_id_huvec_clust_char_vec, bccn1_char_vec,
                              ifelse(cell_id %in% bccn2_char_vec, bccn2_char_vec, "Non-vascular")
        ),
        which_dat = dataset,
        grp_fac = factor(cluster_name,
                         levels = c("Non-vascular", bccn2_char_vec, bccn1_char_vec)
        )
      ) %>%
      arrange(grp_fac)
    column_annots_df1
    
    filtered_test_mat1 <- mat[rownames(mat) %in% row_annots_df1$analyte,
                              colnames(mat) %in% column_annots_df1$u_cell_id,
                              drop = F
    ]
    
    column_annots_df2 <- column_annots_df_temp %>%
      mutate(
        group = ifelse(cell_id %in% vascular_char_vec, cell_id, "Non-vascular"),
        cluster_name = ifelse(cell_id %in% cell_id_huvec_clust_char_vec, bccn1_char_vec,
                              ifelse(cell_id %in% bccn2_char_vec, bccn2_char_vec, "Non-vascular")
        ),
        which_dat = dataset,
        grp_fac = factor(cluster_name, levels = c("Non-vascular", bccn1_char_vec, bccn2_char_vec))
      ) %>%
      arrange(grp_fac)
    column_annots_df2
    
    filtered_test_mat2 <- mat[rownames(mat) %in% row_annots_df2$analyte, 
                              colnames(mat) %in% column_annots_df2$u_cell_id, 
                              drop = F]
    
    split_fn <- str_split(heatmap_output_fn, "\\.", simplify = T)
    heatmap_output_fn1 <- str_c(split_fn[, 1], qq("-@{bccn1_char_vec}."), split_fn[, 2], collapse = "")
    heatmap_output_fn2 <- str_c(split_fn[, 1], qq("-@{bccn2_char_vec}."), split_fn[, 2], collapse = "")
    
    # stop()
    ht1_lst <- organize_and_plot_heatmap_subfunction(base_cell_id = "HUVEC",
                                                     filtered_test_mat = filtered_test_mat1,
                                                     row_annots_df = row_annots_df1, column_annots_df = column_annots_df1,
                                                     heatmap_output_fn = heatmap_output_fn1, title_var = grouping_var
    )
    
    ht2_lst <- organize_and_plot_heatmap_subfunction(base_cell_id = "HAoSMC",
                                                     filtered_test_mat = filtered_test_mat2,
                                                     row_annots_df = row_annots_df2, column_annots_df = column_annots_df2,
                                                     heatmap_output_fn = heatmap_output_fn2, title_var = grouping_var
    )
    
    return(list(
      "HUVEC" = list(
        bare_ht = ht1_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all,
        filtered_row_annot = ht1_lst$row_order_df, col_annot = ht1_lst$column_order_df,
        # transformed_mat = row_minmax_filtered_test_mat,
        plotted_matrix = ht2_lst$reordered_mat
      ),
      "HAoSMC" = list(
        bare_ht = ht2_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all,
        filtered_row_annot = ht2_lst$row_order_df, col_annot = ht2_lst$column_order_df,
        # transformed_mat = row_minmax_filtered_test_mat,
        plotted_matrix = ht2_lst$reordered_mat
      )
    ))
  } else {
    message("Vascular cells clustered together with no other cells. Plotting once...")
    row_annots_df <- filter(
      row_annots_df_include_nonsig_all,
      str_detect(string = base_clust_comp_name, pattern = "HUVEC|HAoSMC"), signif
    )
    row_annots_df
    
    column_annots_df <- column_annots_df_temp %>%
      mutate(
        group = ifelse(cell_id %in% c("HUVEC", "HAoSMC"),
                       "Vascular", "Non-vascular"
        ),
        cluster_name = group,
        which_dat = dataset
      ) %>%
      mutate(grp_fac = factor(cluster_name, levels = c("Non-vascular", "Vascular"))) %>%
      arrange(grp_fac)
    column_annots_df
    
    filtered_test_mat <- mat[rownames(mat) %in% row_annots_df$analyte, 
                             colnames(mat) %in% column_annots_df$u_cell_id,
                             drop = F]
    
    ht_lst <- organize_and_plot_heatmap_subfunction(base_cell_id = "Vascular",
                                                    filtered_test_mat = filtered_test_mat,
                                                    row_annots_df = row_annots_df, column_annots_df = column_annots_df,
                                                    heatmap_output_fn = heatmap_output_fn, title_var = grouping_var
    )
    
    return("Vascular" = list(
      bare_ht = ht_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all,
      filtered_row_annot = ht_lst$row_order_df, col_annot = ht_lst$column_order_df,
      # transformed_mat = row_minmax_filtered_test_mat,
      plotted_matrix = ht_lst$reordered_mat
    ))
  }
  
  #' else if (co_clust & clustered_with_others) {
  #'   message("Vascular cells clustered together, but with other cells. Plotting once...")
  #'   
  #'   #' @NOTE:
  #'   #' in here need to separate out HUVECs just for labeling purposes
  #'   row_annots_df <- filter(
  #'     row_annots_df_include_nonsig_all,
  #'     str_detect(string = base_clust_comp_name, pattern = "HUVEC|HAoSMC"), signif
  #'   )
  #'   row_annots_df
  #'   
  #'   choose_analytes <- rownames(mat) %in% row_annots_df$analyte
  #'   
  #'   clust_name_char <- unique(row_annots_df$base_clust_comp_name)
  #'   clust_name_char
  #'   clust_name_char_split_vec <- unlist(map(clust_name_char, function(c) str_split(c, ",", simplify = T)))
  #'   clust_name_char_split_vec
  #'   
  #'   column_annots_df <- column_annots_df_temp %>%
  #'     mutate(
  #'       group = ifelse(cell_id %in% clust_name_char_split_vec,
  #'                      clust_name_char, "Non-vascular"
  #'       ),
  #'       cluster_name = group,
  #'       which_dat = dataset
  #'     ) %>%
  #'     mutate(grp_fac = factor(cluster_name, levels = c("Non-vascular", clust_name_char))) %>%
  #'     arrange(grp_fac)
  #'   column_annots_df
  #'   
  #'   filtered_test_mat <- mat[choose_analytes, 
  #'                            colnames(mat) %in% column_annots_df$u_cell_id,
  #'                            drop = F]
  #'   
  #'   ht_lst <- organize_and_plot_heatmap_subfunction(base_cell_id = clust_name_char,
  #'                                                   filtered_test_mat = filtered_test_mat,
  #'                                                   row_annots_df = row_annots_df, column_annots_df = column_annots_df,
  #'                                                   heatmap_output_fn = heatmap_output_fn, title_var = grouping_var
  #'   )
  #'   
  #'   return("Combined_Vascular" = list(
  #'     bare_ht = ht_lst$ht, complete_row_annot = row_annots_df_include_nonsig_all,
  #'     filtered_row_annot = ht_lst$row_order_df, col_annot = ht_lst$column_order_df,
  #'     # transformed_mat = row_minmax_filtered_test_mat,
  #'     plotted_matrix = ht_lst$reordered_mat
  #'   ))
  #' } 
}

#' @note returns size of plot in pixels, w x h
calc_ht_size <- function(ht, unit = "inches") {
  pdf(NULL)
  ht <- draw(ht)
  w <- ComplexHeatmap:::width(ht)
  w <- convertX(w, unit, valueOnly = TRUE)
  h <- ComplexHeatmap:::height(ht)
  h <- convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}


#' @note helper function for plotting heatmaps
organize_and_plot_heatmap_subfunction <- function(base_cell_id = NA, 
                                                  filtered_test_mat,
                                                  row_annots_df, column_annots_df,
                                                  heatmap_output_fn, title_var = "pert_iname") {
  
  analytes_ <- rownames(filtered_test_mat)
  u_cell_id_ <- colnames(filtered_test_mat)
  cluster_ids <- column_annots_df %>%
    dplyr::distinct(cluster, cell_id, pert_iname, cluster_name, grp_fac) ; 
  
  # message(tail(cluster_ids))
  
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
    summarize(median_value = mean(value, na.rm = T), .groups = "keep") %>%
    ungroup() %>%
    mutate(u_cell_id = str_c(cell_id, pert_iname, sep = "::")) %>%
    pivot_wider(id_cols = analytes, names_from = u_cell_id, values_from = median_value)
  
  reduced_mat <- reduced_tbl %>% dplyr::select(-analytes) %>% as.matrix()
  rownames(reduced_mat) <- analytes_
  
  n_max_clust <- length(unique(cluster_ids$cluster))
  n_max_clust
  
  
  NAME <- "Mean Log2(value)"
  which_dat <- unique(column_annots_df$which_dat[[1]])
  print_friendly_name <- str_replace(string = str_c(str_split(unique(column_annots_df[title_var][[1]]), "_", simplify = TRUE),
                                                    collapse = " "
  ), pattern = "excl", replacement = "excluding")
  TITLE <- str_c(
    "Differential analytes using connectivity clusters\n",
    str_c(print_friendly_name,
          toupper(which_dat),
          sep = " || "
    )
  )
  # set color breaks at quantiles
  breaks <- sapply(rev(c(.05, .10, .30, .50, .70, .90, .95)), function(q) quantile(reduced_mat, q, na.rm = T, names = F))
  col_fun <- colorRamp2(breaks = breaks, colors = brewer.pal(n = length(breaks), name = "RdYlBu"))
  
  
  ht_opt$COLUMN_ANNO_PADDING <- unit(2, "mm")
  ht_opt$ROW_ANNO_PADDING <- unit(1, "mm")
  ht_opt$annotation_border <- TRUE
  
  
  # in the order of the grouping -> we need to adjust colors depending on group content
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-gcp-heat.RData")
  # load("debug/debug_dat/debug-gcp-heat.RData")
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
  
  vasc_color <- ildf %>% 
    filter(grp_l == "Vascular") %>%
    distinct(vasc_grp_color) %>% 
    .$vasc_grp_color
  
  dom_colors <- ildf %>% 
    filter(lbs %in% c("HUVEC", "HAoSMC")) %>%
    dplyr::select(dom = lbs, dom_color = cell_individual_color) %>%
    bind_rows(tibble(dom = "Non-vascular", dom_color = non_vasc_color)) %>%
    bind_rows(tibble(dom = "Vascular", dom_color = vasc_color))
  
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
    mutate(dom = map_chr(clust_vec, .f = function(x) {
      if ("HUVEC" %in% x & !("HAoSMC" %in% x)) {
        return("HUVEC")
      } else if ("HAoSMC" %in% x & !("HUVEC" %in% x)) {
        return("HAoSMC")
      } else if ("HAoSMC" %in% x & "HUVEC" %in% x){
        return("Vascular")
      } else {
        return("Non-vascular")
      }
    })) %>%
    distinct() %>%
    # mutate(dom = ifelse(is.na(dom), "none", dom)) %>%
    left_join(dom_colors, by= "dom") %>%
    na.omit() %>%
    mutate(dom = str_c(dom, cluster, sep = "::")) %>%
    mutate(which_in_clust_label = map(clust_vec, function(x) {
      init_chr <- str_c(x, collapse = " ")
      res_temp <- str_wrap(init_chr, width = 14); res_temp
      res_temp1 <- str_split(res_temp, pattern = " ", simplify = TRUE); res_temp1
      res_temp2 <- str_c(res_temp1, collapse = ","); res_temp2
      res <- gsub(x = res_temp2, pattern = "\n", replacement = "\n   "); res
      return(res)
    })) %>%
    mutate(dom_info_label = str_c(dom, " ", which_in_clust_label)); dominant_colors
  
  vasc_cluster_id <- cluster_ids %>%
    filter(cell_id %in% c("HUVEC", "HAoSMC")) %>%
    distinct(cell_id, cluster) %>%
    arrange(cell_id); vasc_cluster_id
  
  base_cell <- dominant_colors %>%
    mutate(clust_vec = map(which_in_clust, .f = function(x) str_c(str_split(string = x, pattern = ",", simplify = T)[1,]))) %>%
    mutate(match_base = map_chr(clust_vec, .f = function(x) {
      if ("HUVEC" %in% x & !("HAoSMC" %in% x)) {
        return("HUVEC")
      } else if ("HAoSMC" %in% x & !("HUVEC" %in% x)) {
        return("HAoSMC")
      } else if ("HAoSMC" %in% x & "HUVEC" %in% x){
        return("Vascular")
      } else {
        return("Non-vascular")
      }
    })) %>%
    filter(match_base == base_cell_id) %>%
    .$dom %>%
    unique(); base_cell
  
  grp_factor_lvls_temp <- sort(relevel(factor(unique(dominant_colors$dom)), ref = base_cell)); grp_factor_lvls_temp
  grp_factor_lvls <- sort(factor(grp_factor_lvls_temp, levels = rev(grp_factor_lvls_temp))); grp_factor_lvls
  
  unnested_dom_colors <- dominant_colors %>% 
    unnest(cols = c(clust_vec)); unnested_dom_colors
  
  
  grp_colors <- cluster_ids %>%
    distinct(cell_id, cluster) %>%
    left_join(unnested_dom_colors, by=c("cell_id" = "clust_vec", "cluster")) %>%
    na.omit() %>%
    arrange(cluster) %>%
    dplyr::select(everything(), color = dom_color) %>%
    # indiv cell colors
    left_join(ildf %>% distinct(lbs, cell_individual_color), by = c("cell_id" = "lbs"))
  
  #' @Note:
  #' [1] grp_color_df defines the column ordering!!!
  # grp_factor_lvls <- levels(grp_colors$dom); grp_factor_lvls
  # generate the clustering color palette, and establish the column order in grp_color_df
  grp_color_df <- left_join(cluster_ids,
                            grp_colors %>% distinct(cluster, dom, dom_info_label, color), by = "cluster") %>%
    # in order of the column_df
    mutate(color = ifelse(is.na(color), "#b6b6b6", color)) %>%
    # arrange(dom, pert_iname, cell_id) %>%
    mutate(unique_id = make.unique(str_c(cell_id, pert_iname, sep = "::"))) %>%
    mutate(unique_id = str_split(unique_id, pattern = "_excl_", simplify = T)[,1]) %>%
    mutate(grp_fac_char = as.character(grp_fac)) %>%
    # mutate(dom = ifelse(dom != grp_fac_char, grp_fac_char, dom)) %>%
    mutate(dom = factor(dom, levels = grp_factor_lvls),
           dom = droplevels(dom)) %>%
    arrange(dom, cell_id); grp_color_df
  
  # grayish color for NA
  col_fun_annot1 <- grp_color_df$color # 1st group being closest to the labels on the right
  
  group_order <- unique(grp_color_df$dom_info_label); group_order
  grp_names_char_vec <- grp_color_df$dom
  names(col_fun_annot1) <- grp_names_char_vec
  block_col_fun_annot1 <- grp_color_df %>% distinct(dom, color) %>% .$color
  
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-gcp-dmso-heat.RData")
  # load("debug/debug_dat/debug-gcp-dmso-heat.RData")
  # stop()
  
  {
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
      # set up the colored bars above heatmap
      `Cluster` = anno_block(gp = gpar(
        fill = block_col_fun_annot1,
        labels = group_order,
        labels_gp = gpar(col = "black", font = 2) #, fontsize = FONTSIZE + 2)
      )),
      `Cell line` = cell_names_char_vec,
      `Perturbations` = perturbations_char_vec,
      # color parameters
      col = list(
        `Cluster` = col_fun_annot1,
        `Perturbations` = col_fun_annot2,
        `Cell line` = col_fun_annot3
      ),
      gp = gpar(col = "black"), # this draws black lines in the perturbation annotation block
      annotation_name_gp = gpar(font = 2), #, fontsize = FONTSIZE + 2),
      gap = unit(1, "mm"),
      simple_anno_size = unit(0.5, "cm"),
      show_annotation_name = c(TRUE, TRUE, TRUE),
      show_legend = c(FALSE, FALSE, FALSE),
      annotation_name_side = "left",
      annotation_name_align = TRUE,
      which = "column",
      
      # only need this for simple annotations, like annot_block
      annotation_legend_param = list(`Perturbations` = list(title = "Perturbations"))
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
    
    setdiff(rownames(filtered_test_mat), row_order_df$analyte)
    setdiff(row_order_df$analyte,rownames(filtered_test_mat))
    stopifnot(nrow(row_order_df) == nrow(filtered_test_mat))
    
    # annotations for right side (by row)
    # don't use the boot bh pvalue, use the regular - the boot one is rounded
    signif_log_q_for_bar_plot <- row_order_df$neg_log10_p_val_bh
    # signif_log_q <- row_order_df %>% mutate(`-log10(q)` = round(-log(p_val_bh, base = 10),5) ) %>% .$`-log10(q)`
    signif_q_val <- format(row_order_df$p_val_bh, scientific = TRUE)
    extra_analyte_info_init <- row_order_df$mark
    
    if (which_dat == "gcp") {
      extra_analyte_info_temp <- str_split(extra_analyte_info_init, 
                                           pattern = "me[0-9]*|ac[0-9]*|ph[0-9]*", 
                                           simplify = T) %>%
        as_tibble(.name_repair = "unique") %>%
        suppressMessages() 
      
      # this is going to go in the place of analytes for gcp
      mod_table <- tibble(full_sites = extra_analyte_info_init) %>%
        mutate(s_sites = str_split(extra_analyte_info_init, 
                                   pattern = "(me[0-9]*)|(ac[0-9]*)|(ph[0-9]*)")) %>%
        mutate(mods = str_extract_all(string = full_sites, 
                                      pattern = "(me[0-9]*)|(ac[0-9]*)|(ph[0-9]*)")) %>%
        mutate(mods_idx = gregexpr(text = full_sites, 
                                   pattern = "(me[0-9]*)|(ac[0-9]*)|(ph[0-9]*)")) %>%
        mutate(mod_counts = pmap_chr(list(s_sites, mods, mods_idx), .f = function(s, m, i){
          # idx = 1
          # s = mod_table$s_sites[[idx]]; m = mod_table$mods[[idx]]; i = mod_table$mods_idx[[idx]]
          
          s <- s[s != ""]
          
          diff_i <- tibble(i = i) %>%
            mutate(spacing = i - min(i),
                   unit_space = spacing/3,
                   unit_diff = diff(c(0, unit_space))) %>%
            mutate(together = ifelse(lead(unit_diff) <= 1 | unit_diff == 1, T, F)) %>%
            rownames_to_column(var = "id") %>%
            # group_by(id, together) %>%
            distinct(id, together) %>%
            mutate(together = ifelse(is.na(together), F, together)); diff_i
          
          if (any(diff_i$together)) {
            m_ids_df <- tibble(m) %>%
              rownames_to_column("id") %>%
              left_join(diff_i, by= "id") %>%
              group_by(together) %>%
              # row ids denote the id of the s character vector
              summarize(marks = str_c(m, collapse = ","),
                        .groups = "keep") %>%
              ungroup() %>%
              dplyr::select(-together); m_ids_df
            m_ids <- m_ids_df$marks
            #names(m_ids) <- s
            
            res <- tibble(s, m_ids) %>%
              mutate(cmbd = str_c(s, ": ",  m_ids)) %>%
              summarize(final_cmbd = str_c(cmbd, collapse = " | ")) %>%
              .$final_cmbd; res
          } else {
            m_ids <- m
            res <- tibble(s, m_ids) %>%
              mutate(cmbd = str_c(s, ": ",  m_ids)) %>%
              summarize(final_cmbd = str_c(cmbd, collapse = " | ")) %>%
              .$final_cmbd %>%
              unique(); res
          }
          
          return(res)
        }))
      
      # this goes in the 'site' column
      colnames_ <- colnames(extra_analyte_info_temp)
      extra_analyte_info_df <- extra_analyte_info_temp %>%
        unite(col = colnames_, sep = "")
      extra_analyte_info <- extra_analyte_info_df$colnames_
      
    } else {
      # p100 is easy
      extra_analyte_info <- extra_analyte_info_init
    }
    
    d_statistics_vec <- round(row_order_df$d_stat_val, 3)
    # fold change
    fc_vec <- round(2^row_order_df$logFC, 3) # logFC is log_2
    analytes_ordered <- row_order_df$analyte
    color_fc_row_labels <- tibble(fc_vec) %>%
      mutate(color = ifelse(fc_vec < 1, "red", "black")) %>%
      .$color
    
    # order of rows
    analytes_reordered_df <- row_order_df %>%
      mutate(new_analyte_label = ifelse(signif_and_fold, str_c(analyte, " **"), analyte)) 
    
    analytes_reordered_labels <- analytes_reordered_df$new_analyte_label
    analytes_reordered <- analytes_reordered_df$analyte
    
    if (which_dat == "gcp") {
      analytes_reordered_labels <- row_order_df %>% 
        left_join(mod_table, by = c("analyte" = "full_sites")) %>%
        mutate(new_analyte_label = ifelse(signif_and_fold, str_c(mod_counts, " **"), mod_counts)) %>%
        .$new_analyte_label
    }
    
    right_ha <- HeatmapAnnotation(
      `d_statistic` = anno_text(d_statistics_vec,
                                just = "center", location = 0.5,
                                gp = gpar(# font = "Arial",
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
      ), #  
      # `-log10(q)` = anno_barplot(signif_log_q_for_bar_plot), # signif_log_q
      
      `site` = anno_text(extra_analyte_info,
                         just = "center", location = 0.5, # unit(1, "npc")
                         gp = gpar(border = "black"), #,fontsize = FONTSIZE),
                         width = max_text_width(extra_analyte_info) * 1.075
      ),
      `analyte` = anno_text(analytes_reordered_labels,
                            just = "left", location = 0.1, # unit(1, "npc")
                            gp = gpar(font = 4), #, fontsize = FONTSIZE + 2),
                            width = max_text_width(extra_analyte_info) * 1.075
      ),
      which = "row",
      gap = unit(1.5, "mm") # gap between annotations
    )
    
    
    # reorder the matrix with this -- this could also be done in the heatmap function;
    # however, be careful about the variable you use to assign column order and column attributes!
    # ESPECIALLY, column split
    dim(reduced_mat)
    length(analytes_reordered)
    length(grp_color_df$unique_id)
    
    reordered_mat_temp <- reduced_mat[analytes_reordered, grp_color_df$unique_id, drop = F]
    reordered_mat <- apply(reordered_mat_temp, 2, FUN = as.numeric)
    


    # base_cell_indices <- which(grp_color_df$grp_fac == base_cell_id)
    # base_grp_color <- unique(grp_color_df %>% filter(grp_fac == base_cell_id) %>% .$color)
    # rest_cell_color <- "darkgray"
    # rest_cell_indices <- which(grp_color_df$grp_fac != base_cell_id)
    
    # # scale of the boxplot annotation
    # rg <- range(reordered_mat, na.rm = TRUE)
    # rg[1] = rg[1] - (rg[2] - rg[1])* 0.02
    # rg[2] = rg[2] + (rg[2] - rg[1])* 0.02
    # 
    # anno_multiple_boxplot = function(index) {
    #   # rows of matrix
    #   nr = length(index)
    #   # message(nr)
    #   pushViewport(viewport(xscale = rg, yscale = c(0.5, nr + 0.5)))
    #   for(i in seq_along(index)) {
    #     grid.rect(y = nr-i+1, height = 1, default.units = "native")
    #     # first group
    #     # grid.abline(intercept = 0, gp = gpar(lty = 3,lwd = 3, col = "red"))
    #     grid.boxplot(reordered_mat[ index[i], base_cell_indices] %>% na.omit(),
    #                  pos = nr-i+1 + 0.2,
    #                  box_width = 0.4, size = unit(1, "mm"),
    #                  gp = gpar(fill = base_grp_color),
    #                  direction = "horizontal")
    #     # second group
    #     grid.boxplot(reordered_mat[ index[i], rest_cell_indices] %>% na.omit(),
    #                  pos = nr-i+1 - 0.2,
    #                  box_width = 0.4,
    #                  size = unit(1, "mm"),
    #                  gp = gpar(fill = rest_cell_color),
    #                  direction = "horizontal")
    #     
    #   }
    #   grid.xaxis(gp = gpar(fontsize = 8))
    #   popViewport()
    # }

    # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-p100-heat.RData")
    # load("debug/debug_dat/debug-p100-heat.RData")
    # stop()
    
    # get max of finite values
    # bar_plot_q_values <- signif_log_q_for_bar_plot
    # bar_plot_q_values[is.infinite(bar_plot_q_values)] = max(bar_plot_q_values[is.finite(bar_plot_q_values)], na.rm = TRUE)
    # 
    # left_ha <- rowAnnotation(
    #   boxplot = anno_multiple_boxplot,
    #   `bar` = anno_barplot(bar_plot_q_values, 
    #                        gp = gpar(fontsize = 8, fill = "darkgray"),
    #                        axis_param = list(labels_rot = 0)),
    #   width = unit(4, "cm"),
    #   gap = unit(2, "mm"),
    #   # annotation_name_rot = 45,
    #   show_annotation_name = FALSE
    #   # show_annotation_names = FALSE
    # # signif_log_q_for_bar_plot)
    #   # width = unit(4, "cm"),
    #   # gap = unit(1.5, "mm")
    #   # annotation_legend_param = list(`bar` = list(title = "Q-value"))
    # )
    # draw(left_ha)
    
    
    ht_opt$message = FALSE
    
    ht <- Heatmap(reordered_mat,
                  name = NAME,
                  col = col_fun,
                  na_col = "gray",
                  top_annotation = top_ha,
                  right_annotation = right_ha,
                  # left_annotation = left_ha,
                  
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
    
    # draw(ht)
    # legend for cell lines
    lgd_cluster <- Legend(labels = group_order, 
                          title = "Cluster", 
                          legend_gp = gpar(fill = block_col_fun_annot1))
    lgd_cell <- Legend(labels = unique(cell_names_char_vec), 
                       title = "Cell Line", 
                       legend_gp = gpar(fill = unique(col_fun_annot3)),  
                       direction = "horizontal", 
                       ncol = 2)
    lgd_pert <- Legend(labels = unique(perturbations_char_vec), 
                       title = "Perturbations", 
                       legend_gp = gpar(fill = unique(col_fun_annot2)),  
                       direction = "horizontal", 
                       ncol = 2)
    
    # lgd_box <- Legend(labels = c(base_cell_id, "All other clusters"), title = "Boxplots",
    #                   legend_gp = gpar(fill = c(base_grp_color, rest_cell_color)))
    
    lgd_lst <- packLegend(lgd_cluster, lgd_cell, lgd_pert, # lgd_box,
                          max_width = unit(20, "cm"), 
                          direction = "horizontal")

    
    init_ht <- draw(ht,
                    heatmap_legend_side = "bottom",
                    annotation_legend_list = lgd_lst,
                    annotation_legend_side = "left",
                    column_title = TITLE, 
                    column_title_gp = gpar(font = 2),
                    merge_legend = TRUE)
    
    size <- calc_ht_size(init_ht)
  }
  
  ## Save the plot
  # in pixels
  linear_height <- size[2] + 2
  linear_height # 9
  linear_width <- size[1] + 5
  linear_width # 16
  
  # get a new filename for pdf outoput
  heatmap_output_fn_pdf <- file.path(str_c(str_split(heatmap_output_fn, pattern = "\\.", simplify = TRUE)[,1], ".pdf", sep = ""))
  
  to_plot <- function() {
    
    draw(init_ht)
    right_annots_height <- unit(1.01, "npc") + unit(3, "mm")
    height_top_dec <- convertX(ComplexHeatmap:::height(top_ha), "mm", valueOnly = TRUE)
    
    # add the annotation titles
    decorate_annotation("site", {
      grid.text("Site",
                y = right_annots_height, just = "top",
                gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
      )
    })
    
    # decorate_annotation("bar", {
    #   # grid.text("barplot", unit(-10, "mm"), just = "bottom", rot = 90)
    #   grid.lines(y=c(0, 1), unit(c(1, 1), "native"), 
    #              gp = gpar(lty = 2, lwd = 2, col = "red"))
    # })
    # decorate_annotation("boxplot", {
    #   # grid.text("barplot", unit(-10, "mm"), just = "bottom", rot = 90)
    #   grid.lines(y=c(0, 1), unit(abs(min(reordered_mat, na.rm = T)), "native"), 
    #              gp = gpar(lty = 2, lwd = 2, col = "red"))
    # })
    
    
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
      grid.text("fc",
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
    
  }
  
  plot_eps <- function(){
    setEPS()
    # par(mar=c(2, 2, 2, 2), oma=c(2,2,1,1))
    postscript(bg = "white",
               file = heatmap_output_fn, horizontal = F,
               height = linear_height, width = linear_width, onefile = F
    )
    # plot
    to_plot()
    # turn off device
    dev.off()
  }
  plot_pdf <- function(){
    # par(mar=c(2, 2, 2, 2), oma=c(2,2,1,1))
    pdf(file = heatmap_output_fn_pdf, height = linear_height, width = linear_width)
    # plot
    to_plot()
    # turn off device
    dev.off()
  }
  # stop()
  
  plot_eps()
  plot_pdf()
  message("Done!")
  return(list(ht, row_order_df, column_annots_df, reordered_mat) %>%
           set_names("ht", "row_order_df", "column_order_df", "reordered_mat"))
}


############## HEATMAP FUNCTION ####################
############## HEATMAP FUNCTION ####################

#' n_max_clust <- length(unique(column_annots_df$cluster))
#' n_max_clust
#' 
#' # perturbations, in the column order
#' perturbations_char_vec_temp <- str_split(string = column_annots_df$pert_iname, pattern = "_", simplify = TRUE)[, 1]
#' 
#' extra_pert_info <- read_excel(file.path(references_directory, "Drug Glossary_edited.xlsx")) %>%
#'   mutate(pert_iname = tolower(Drug), moa_simple = `MOA simplified`) %>%
#'   dplyr::select(pert_iname, moa_simple) %>%
#'   mutate(moa_simple = ifelse(is.na(moa_simple), "", moa_simple))
#' 
#' link_moa_to_pert <- tibble(pert_iname = perturbations_char_vec_temp) %>%
#'   left_join(extra_pert_info, by = "pert_iname") %>%
#'   rename(pert_iname_temp = pert_iname) %>%
#'   mutate(pert_iname = str_c(pert_iname_temp, " (", moa_simple, ")", sep = ""))
#' 
#' perturbations_char_vec <- link_moa_to_pert$pert_iname
#' n_u_perts <- length(unique(perturbations_char_vec))
#' n_u_perts
#' 
#' # get top up and down
#' # needs to be relative to 1
#' top_up <- row_annots_df %>%
#'   filter(fc >= 1) %>%
#'   arrange(desc(`fc`))
#' top_up
#' top_down <- row_annots_df %>%
#'   filter(fc < 1) %>%
#'   # filter(logFC < LOGFC_CUTOFF) %>%
#'   arrange(desc(`fc`))
#' top_down
#' row_order_df <- bind_rows(top_up, top_down)
#' row_order_df
#' 
#' stopifnot(nrow(row_order_df) == nrow(filtered_test_mat))
#' 
#' # annotations for right side (by row)
#' signif_log_q_for_bar_plot <- row_order_df$`-log10(q)`
#' # signif_log_q <- row_order_df %>% mutate(`-log10(q)` = round(-log(p_val_bh, base = 10),5) ) %>% .$`-log10(q)`
#' signif_q_val <- round(row_order_df$p_val_bh, 5)
#' extra_analyte_info <- row_order_df$mark
#' d_statistics_vec <- round(row_order_df$d_stat_val, 3)
#' # fold change
#' fc_vec <- round(2^row_order_df$logFC, 3) # logFC is log_2
#' analytes <- row_order_df$analyte
#' color_fc_row_labels <- tibble(fc_vec) %>%
#'   mutate(color = ifelse(fc_vec < 1, "red", "black")) %>%
#'   .$color
#' 
#' 
#' # tail(reordered_mat)
#' ## PLOTTING
#' FONTSIZE <- 7
# print_friendly_name <- str_replace(string = str_c(str_split(unique(column_annots_df[title_var][[1]]), "_", simplify = TRUE),
#                                                   collapse = " "
# ), pattern = "excl", replacement = "excluding")
# TITLE <- str_c(
#   "Differential analytes using connectivity clusters\n",
#   str_c(print_friendly_name,
#         toupper(unique(column_annots_df$which_dat[[1]])),
#         sep = " || "
#   )
# )
#' TITLE
#' 
#' NAME <- "Log2(value)"
#' NAME
#' 
#' ht_opt$COLUMN_ANNO_PADDING <- unit(2, "mm")
#' ht_opt$ROW_ANNO_PADDING <- unit(1, "mm")
#' ht_opt$annotation_border <- TRUE
#' 
#' # set color breaks at quantiles
#' breaks <- sapply(rev(c(.05, .10, .30, .50, .70, .90, .95)), function(q) quantile(reordered_mat, q, na.rm = T, names = F))
#' col_fun <- colorRamp2(breaks = breaks, colors = brewer.pal(n = length(breaks), name = "RdYlBu"))
#' 
#' # in the order of the grouping -> we need to adjust colors depending on group content
#' 
#' # stop()
#' informative_labels_lst <- create_informative_labels()
#' ildf <- informative_labels_lst$informative_labels_df %>%
#'   mutate(
#'     grp = map_chr(grp, function(x) str_split(x, pattern = " ", simplify = T)[, 1]),
#'     grp_l = map_chr(grp_l, function(x) str_split(x, pattern = " ", simplify = T)[, 1]),
#'     grp_l2 = map_chr(grp_l2, function(x) str_split(x, pattern = " ", simplify = T)[, 1])
#'   ) %>%
#'   left_join(column_annots_df %>% 
#'               dplyr::distinct(cell_id, cluster, cluster_name), 
#'             by = c("lbs" = "cell_id")) 
#' 
#' non_vasc_color <- ildf %>% 
#'   filter(grp_l == "Non-vascular") %>%
#'   distinct(vasc_grp_color) %>% 
#'   .$vasc_grp_color
#' 
#' dom_colors <- ildf %>% 
#'   filter(lbs %in% c("HUVEC", "HAoSMC")) %>%
#'   dplyr::select(dom = lbs, dom_color = cell_individual_color) %>%
#'   bind_rows(tibble(dom = "Non-vascular", dom_color = non_vasc_color))
#' 
#' dominant_colors <- ildf %>%
#'   distinct(cluster_name, lbs, cluster, cell_individual_color) %>%
#'   group_by(cluster) %>%
#'   summarize(
#'     which_in_clust = str_c(sort(unique(lbs)), collapse = ","),
#'     .groups = "keep"
#'   ) %>%
#'   ungroup() %>%
#'   arrange(cluster) %>%
#'   mutate(clust_vec = map(which_in_clust, .f = function(x) str_c(str_split(string = x, pattern = ",", simplify = T)[1,]))) %>%
#'   mutate(dom = map(clust_vec, .f = function(x) {
#'     if ("HUVEC" %in% x){
#'       return("HUVEC")
#'     } else if ("HAoSMC" %in% x) {
#'       return("HAoSMC")
#'     } else {
#'       return("Non-vascular")
#'     }
#'   })) %>%
#'   unnest(cols = "dom") %>%
#'   distinct() %>%
#'   # mutate(dom = ifelse(is.na(dom), "none", dom)) %>%
#'   left_join(dom_colors, by= "dom")
#' 
#' vasc_cluster_id <- column_annots_df %>%
#'   filter(cell_id %in% c("HUVEC", "HAoSMC")) %>%
#'   distinct(cell_id, cluster) %>%
#'   arrange(cell_id)
#' 
#' # you know which cell line is the "base" by looking at the last row
#' # this must be a factor so we get order!
#' base_cell <- tail(column_annots_df, n = 1)$grp_fac
#' grp_factor_lvls <- levels(base_cell)
#' 
#' grp_colors <- column_annots_df %>%
#'   distinct(cell_id, cluster) %>%
#'   mutate(dom = ifelse(cluster %in% vasc_cluster_id$cluster[1], "HAoSMC",
#'                       ifelse(cluster %in% vasc_cluster_id$cluster[2], "HUVEC", "Non-vascular"))) %>%
#'   left_join(dominant_colors %>% dplyr::select(-cluster), by="dom") %>%
#'   na.omit() %>%
#'   distinct(cell_id, cluster, dom, dom_color) %>%
#'   arrange(cluster) %>%
#'   dplyr::select(everything(), color = dom_color) %>%
#'   left_join(ildf %>% distinct(lbs, cell_individual_color), by = c("cell_id" = "lbs")) %>%
#'   mutate(dom = factor(dom, levels = grp_factor_lvls))
#' 
#' 
#' 
#' #' @Note:
#' #' [1] grp_color_df defines the column ordering!!!
#' 
#' # generate the clustering color palette, and establish the column order in grp_color_df
#' grp_color_df <- left_join(column_annots_df,
#'                           grp_colors %>% distinct(cluster, dom, color), by = "cluster") %>%
#'   # in order of the column_df
#'   mutate(color = ifelse(is.na(color), "#b6b6b6", color)) %>%
#'   arrange(dom, cell_id, pert_iname)
#' # grayish color for NA
#' col_fun_annot1 <- grp_color_df$color # 1st group being closest to the labels on the right
#' 
#' group_order <- levels(grp_color_df$dom)
#' grp_names_char_vec <- grp_color_df$dom
#' names(col_fun_annot1) <- grp_names_char_vec
#' block_col_fun_annot1 <- unique(col_fun_annot1)
#' 
#' 
#' # reorder the matrix with this -- this could also be done in the heatmap function;
#' # however, be careful about the variable you use to assign column order and column attributes!
#' # ESPECIALLY, column split
#' reordered_mat_temp <- filtered_test_mat[analytes, grp_color_df$u_cell_id, drop = F]
#' reordered_mat <- apply(reordered_mat_temp, 2, FUN = as.numeric)
#' 
#' 
#' 
#' # generate pert color palette
#' col_palette_2_df <- generate_color_palette()
#' col_fun_annot2_df <- inner_join(col_palette_2_df,
#'                                 link_moa_to_pert %>% distinct(pert_iname_temp, pert_iname),
#'                                 by = c("pert_iname" = "pert_iname_temp")
#' ) %>%
#'   rename(full_pert_name = pert_iname.y)
#' 
#' col_fun_annot2 <- col_fun_annot2_df$colors
#' names(col_fun_annot2) <- col_fun_annot2_df$full_pert_name
#' 
#' # generate the cell id color palette
#' cell_color_df <- grp_color_df %>%
#'   dplyr::select(cell_id) %>%
#'   left_join(grp_colors, by = "cell_id") %>%
#'   mutate(color = ifelse(is.na(color), "#b6b6b6", color)) # grayish color for NA
#' col_fun_annot3 <- cell_color_df$cell_individual_color # 1st group being closest to the labels on the right
#' 
#' cell_names_char_vec <- cell_color_df$cell_id
#' names(col_fun_annot3) <- cell_names_char_vec
#' # block_col_fun_annot3 <- unique(col_fun_annot3)
#' 
#' # top annotations
#' top_ha <- HeatmapAnnotation(
#'   `Cluster` = anno_block(gp = gpar(
#'     fill = block_col_fun_annot1,
#'     labels = group_order,
#'     labels_gp = gpar(col = "black", font = 2) #, fontsize = FONTSIZE + 2)
#'   )),
#'   `Cell line` = cell_names_char_vec,
#'   `Perturbations` = perturbations_char_vec,
#'   col = list(
#'     `Cluster` = col_fun_annot1,
#'     `Perturbations` = col_fun_annot2,
#'     `Cell line` = col_fun_annot3
#'   ),
#'   
#'   # gp = gpar(col = "black"), # this draws black lines in the perturbation annotation block
#'   annotation_name_gp = gpar(font = 2), #, fontsize = FONTSIZE + 2),
#'   gap = unit(1, "mm"),
#'   simple_anno_size = unit(0.5, "cm"),
#'   which = "column",
#'   # where do the names of the bars go?
#'   annotation_name_side = "left",
#'   
#'   # show_legend = c(TRUE),
#'   
#'   annotation_legend_param = list(`Perturbations` = list(ncol = 1))
#' )
#' 
#' # stop()
#' # right annotations
#' 
#' right_ha <- HeatmapAnnotation(
#'   `d_statistic` = anno_text(d_statistics_vec,
#'                             just = "center", location = 0.5,
#'                             gp = gpar(
#'                               # fontsize = FONTSIZE,
#'                               border = "black"
#'                             ),
#'                             width = max_text_width(d_statistics_vec) * 1.075
#'   ),
#'   `FoldChange` = anno_text(fc_vec,
#'                            just = "center", location = 0.5, # log_fc_vec
#'                            gp = gpar(
#'                              # fontsize = FONTSIZE,
#'                              col = color_fc_row_labels,
#'                              border = "black"
#'                            ),
#'                            width = max_text_width(fc_vec) * 1.075
#'   ),
#'   `q-value` = anno_text(signif_q_val,
#'                         just = "center", location = 0.5,
#'                         gp = gpar(border = "black"), #, fontsize = FONTSIZE),
#'                         width = max_text_width(signif_q_val) * 1.075
#'   ), #  `-log10(q)`, anno_barplot(signif_log_q_for_bar_plot, extend = 0.25), # signif_log_q
#'   
#'   `site` = anno_text(extra_analyte_info,
#'                      just = "center", location = 0.5, # unit(1, "npc")
#'                      gp = gpar(border = "black"), #,fontsize = FONTSIZE),
#'                      width = max_text_width(extra_analyte_info) * 1.075
#'   ),
#'   `analyte` = anno_text(analytes,
#'                         just = "left", location = 0.1, # unit(1, "npc")
#'                         gp = gpar(font = 4), #, fontsize = FONTSIZE + 2),
#'                         width = max_text_width(extra_analyte_info) * 1.075
#'   ),
#'   which = "row",
#'   gap = unit(1.5, "mm") # gap between annotations
#' )