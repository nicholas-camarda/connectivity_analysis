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
  
  
  color_fn <- file.path(data_directory, "color_mapping.rds")
  if (!file.exists(color_fn)) {
    n <- nrow(all_drugs_mapping)
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    palette <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                             rownames(qual_col_pals)))[1:n]
    
    # palette <- distinctColorPalette(k = n)
    # palette <- pals::alphabet(n = n)
    # palette <- colorspace::qualitative_hcl(n, palette = "Dynamic")
    all_drugs_mapping_final <- all_drugs_mapping %>% 
      bind_cols(colors = palette)
    write_rds(all_drugs_mapping_final, file = color_fn)
  } else {
    all_drugs_mapping_final <- read_rds(color_fn)
  }
  
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
  
  # stop()
  # generate the row annotations
  row_annots_df_include_nonsig_all <- diff_ex_df %>%
    mutate(`-log10(q)` = -log(p_val_boot_bh, base = 10)) %>%
    # only handles the boot statistic -- make it handle the regular stat
    mutate(d_stat_val = ifelse(directional_stat == "--" & !is.na(ks_boot_statistic), 
                               -1 * ks_boot_statistic, 
                               ks_boot_statistic))
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-p100-pre-heat.RData")
  # load("debug/debug_dat/debug-p100-pre-heat.RData")
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
    message("Vascular cells clustered separately, alone or with other cells. Separating heatmap plots...")
    
    # extract significant analytes per cluster
    row_annots_df1 <- filter(
      row_annots_df_include_nonsig_all,
      str_detect(string = base_clust_comp_name, pattern = "HUVEC"), signif
    ); row_annots_df1
    
    row_annots_df2 <- filter(
      row_annots_df_include_nonsig_all,
      str_detect(string = base_clust_comp_name, pattern = "HAoSMC"), signif
    ); row_annots_df2
    
    
    # get the 'vasc-centric' character vector of the cluster name, i.e. just HUVEC, or HAoSMC
    bccn1_char_vec <- str_extract(string = unique(row_annots_df1$base_clust_comp_name),
                                  pattern = "HUVEC"); bccn1_char_vec
    bccn2_char_vec <- str_extract(string = unique(row_annots_df2$base_clust_comp_name), 
                                  pattern = "HAoSMC"); bccn2_char_vec
    
    # this will be the full cluster name, i.e. A375, HUVEC, YAPC
    cell_id_huvec_clust_char_vec <- as.character(unique(str_split(unique(row_annots_df1$base_clust_comp_name), ",", simplify = T)))
    cell_id_huvec_clust_char_vec
    
    cell_id_haosmc_clust_char_vec <- as.character(unique(str_split(unique(row_annots_df2$base_clust_comp_name), ",", simplify = T)))
    cell_id_haosmc_clust_char_vec
    
    column_annots_df1 <- column_annots_df_temp %>%
      mutate(
        group = ifelse(cell_id %in% vascular_char_vec, cell_id, "Non-vascular"),
        cluster_name = ifelse(cell_id %in% cell_id_huvec_clust_char_vec, bccn1_char_vec,
                              ifelse(cell_id %in% cell_id_haosmc_clust_char_vec, bccn2_char_vec, "Non-vascular")
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
    
    #' #' @CHECK check og table
    #' dat_tbl %>%
    #' filter(cell_id == "HUVEC", pr_gene_symbol == "pS12 EIF4A3") %>%
    #'   mutate(median_value = median(value, na.rm = TRUE)) %>%
    #'   dplyr::select(pr_gene_symbol, cell_id, value, median_value)
    #' 
    #' #' @CHECK check mat
    #' mat %>%
    #'   rownames_to_column() %>%
    #'   rename(analyte = rowname) %>%
    #'   as_tibble() %>%
    #'   pivot_longer(cols = colnames(.)[-1], names_to = 'joint') %>%
    #'   separate(col = joint, into = c("joint2","replicate_id", "plate_id"), sep = "::") %>%
    #'   separate(col = joint2, into = c("cell_id", "pert_iname","pert_class")) %>%
    #'   filter(cell_id == "HUVEC", analyte == "pS12 EIF4A3") %>%
    #'   mutate(median_value = median(value, na.rm = TRUE))
    #' 
    #' #' @CHECK check filtered_mat
    #' filtered_test_mat1 %>%
    #'  rownames_to_column() %>%
    #'   rename(analyte = rowname) %>%
    #'   as_tibble() %>%
    #'   pivot_longer(cols = colnames(.)[-1], names_to = 'joint') %>%
    #'   separate(col = joint, into = c("joint2","replicate_id", "plate_id"), sep = "::") %>%
    #'   separate(col = joint2, into = c("cell_id", "pert_iname","pert_class")) %>%
    #'   filter(cell_id == "HUVEC", analyte == "pS12 EIF4A3") %>%
    #'   mutate(median_value = median(value, na.rm = TRUE))
    
    
    column_annots_df2 <- column_annots_df_temp %>%
      mutate(
        group = ifelse(cell_id %in% vascular_char_vec, cell_id, "Non-vascular"),
        cluster_name = ifelse(cell_id %in% cell_id_huvec_clust_char_vec, bccn1_char_vec,
                              ifelse(cell_id %in% cell_id_haosmc_clust_char_vec, bccn2_char_vec, "Non-vascular")
        ),
        which_dat = dataset,
        grp_fac = factor(cluster_name, levels = c("Non-vascular", bccn1_char_vec, bccn2_char_vec))
      ) %>%
      arrange(grp_fac)
    column_annots_df2
    
    filtered_test_mat2 <- mat[rownames(mat) %in% row_annots_df2$analyte, 
                              colnames(mat) %in% column_annots_df2$u_cell_id, 
                              drop = F]
    
    # name the output FNs
    split_fn <- str_split(heatmap_output_fn, "\\.", simplify = T)
    heatmap_output_fn1 <- str_c(split_fn[, 1], qq("-@{bccn1_char_vec}."), split_fn[, 2], collapse = "")
    heatmap_output_fn2 <- str_c(split_fn[, 1], qq("-@{bccn2_char_vec}."), split_fn[, 2], collapse = "")
    
    # stop()
    ht1_lst <- organize_and_plot_heatmap_subfunction(base_cell_id = "HUVEC",
                                                     filtered_test_mat = filtered_test_mat1,
                                                     row_annots_df = row_annots_df1, 
                                                     column_annots_df = column_annots_df1,
                                                     heatmap_output_fn = heatmap_output_fn1, 
                                                     title_var = grouping_var
    )
    
    ht2_lst <- organize_and_plot_heatmap_subfunction(base_cell_id = "HAoSMC",
                                                     filtered_test_mat = filtered_test_mat2,
                                                     row_annots_df = row_annots_df2, 
                                                     column_annots_df = column_annots_df2,
                                                     heatmap_output_fn = heatmap_output_fn2, 
                                                     title_var = grouping_var
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
    message("Vascular cells clustered together with no other cells. Plotting heatmap once...")
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
                                                  filtered_test_mat = NA,
                                                  row_annots_df = NA, column_annots_df = NA,
                                                  heatmap_output_fn = NA, title_var = NA,
                                                  plot_cluster_bar = FALSE,
                                                  add_D_stat = FALSE) {
  message("Plotting for: ", base_cell_id)
  cluster_ids <- column_annots_df %>%
    dplyr::distinct(cluster, cell_id, pert_iname, cluster_name, grp_fac) ; 
  
  # message(tail(cluster_ids))
  reduced_long <- filtered_test_mat %>% 
    rownames_to_column() %>%
    rename(analytes = rowname) %>%
    as_tibble() %>%
    pivot_longer(cols = colnames(.)[-1], names_to = 'joint') %>%
    separate(col = joint, into = c("joint2","replicate_id", "plate_id"), sep = "::") %>%
    separate(col = joint2, into = c("cell_id", "pert_iname","pert_class"), sep = "--")  %>%
    left_join(cluster_ids, by = c("cell_id", "pert_iname"))

  reduced_tbl <- reduced_long %>%
    group_by(analytes, cell_id, pert_iname, cluster) %>%
    summarize(median_value = median(value, na.rm = T), .groups = "keep") %>%
    ungroup() %>%
    # summarize(mean_value = mean(value, na.rm = T), .groups = "keep") %>%
    mutate(u_cell_id = str_c(cell_id, pert_iname, sep = "::")) %>%
    dplyr::select(-cell_id, -pert_iname) %>%
    pivot_wider(id_cols = analytes, 
                names_from = u_cell_id, 
                values_from = median_value); reduced_tbl
  
  reduced_mat <- reduced_tbl %>% dplyr::select(-analytes) %>% as.matrix()
  rownames(reduced_mat) <- reduced_tbl$analytes
  # reduced_mat["pS12 EIF4A3",]
  
  n_max_clust <- length(unique(cluster_ids$cluster))
  n_max_clust
  
  
  NAME <- "Median Log2(value)"
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
  ); TITLE
  # set color breaks at quantiles
  # qnts <- seq(0, 1, 0.125)
  qnts <- c(0.05, 0.10, 0.30, 0.50, 0.70, 0.90, 0.95)
  # qnts <- c(0.025, 0.05, 0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.975)
  # qnts <- c(0.001, 0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975, 0.999)
  # qnts <- c(0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.95)
  #  c(0.001, 0.025, 0.05, 0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.975, 0.999)
  # qnts <- c(0.001, 0.025, 0.05, 0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.975, 0.999)
  # palette_colors <- diverge_hcl(n = length(qnts), h = c(260, 0), c = 81, l = c(30, 90), power = 1.5)
  # palette_colors <- wes_palette("Zissou1", length(qnts), type = "continuous")
  
  breaks <- sapply(rev(qnts), function(q) quantile(reduced_mat, q, na.rm = T, names = F))
  palette_colors <- brewer.pal(n = length(breaks), name = "RdYlBu")
  col_fun <- colorRamp2(breaks = breaks, colors = palette_colors)
  
  # stop()
  # in the order of the grouping -> we need to adjust colors depending on group content
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-p100-heat.RData")
  # load("debug/debug_dat/debug-p100-heat.RData")
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
    arrange(dom, pert_iname) %>%
    left_join(ildf %>% distinct(lbs, new_lbs),
              by = c("cell_id" = "lbs")); grp_color_df
  
  # grayish color for NA
  col_fun_annot1 <- grp_color_df$color # 1st group being closest to the labels on the right
  
  group_order <- unique(grp_color_df$dom_info_label); group_order
  
  grp_names_char_vec <- grp_color_df$dom
  names(col_fun_annot1) <- grp_names_char_vec
  block_col_fun_annot1 <- grp_color_df %>% distinct(dom, color) %>% .$color
  
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-heat-2.RData")
  # load("debug/debug_dat/debug-heat-1.RData")
  # load("debug/debug_dat/debug-heat-2.RData")
  # stop()
  
  # perturbations, in the column order
  # must be non-unique
  perturbations_char_vec_temp <- str_split(string = grp_color_df$pert_iname, 
                                           pattern = "_", 
                                           simplify = TRUE)[, 1]
  
  extra_pert_info <- read_excel(file.path(references_directory, 
                                          "Drug Glossary_edited.xlsx")) %>%
    mutate(pert_iname = tolower(Drug), moa_simple = `MOA simplified`) %>%
    mutate(pert_iname = ifelse(pert_iname == "dmso", toupper(pert_iname), pert_iname)) %>%
    dplyr::select(pert_iname, moa_simple) %>% # Class
    mutate(moa_simple = ifelse(is.na(moa_simple), "", moa_simple)) %>%
    distinct()
  
  link_moa_to_pert <- tibble(pert_iname = perturbations_char_vec_temp) %>%
    left_join(extra_pert_info, by = "pert_iname") %>%
    rename(pert_iname_temp = pert_iname) %>%
    mutate(pert_iname = str_c(pert_iname_temp, " (", moa_simple, ")", 
                              sep = "")) 
  
  perturbations_char_vec <- link_moa_to_pert$pert_iname
  n_u_perts <- length(unique(perturbations_char_vec))
  n_u_perts
  
  # generate pert color palette
  col_palette_2_df <- generate_color_palette()
  col_fun_annot2_df <- left_join(link_moa_to_pert %>% 
                                   dplyr::select(pert_iname_temp, pert_iname),
                                 col_palette_2_df,
                                 by = c("pert_iname_temp" = "pert_iname")
  ) %>%
    rename(full_pert_name = pert_iname); col_fun_annot2_df
  
  col_fun_annot2 <- col_fun_annot2_df$colors
  names(col_fun_annot2) <- col_fun_annot2_df$full_pert_name
  
  # generate the cell id color palette
  cell_color_df <- grp_color_df %>%
    dplyr::select(cell_id, new_lbs) %>%
    # fix this to render properly
    mutate(new_lbs = map(new_lbs, .f = function(s) {
      str_split(string = s, pattern = "\n", simplify = TRUE)})) %>%
    mutate(new_lbs = map_chr(new_lbs, function(s){
      str_c(str_trim(s, side = "both"), collapse = " ")})) %>%
    mutate(new_lbs = map2_chr(cell_id, new_lbs, .f = function(s, y) {
      res <- ifelse(str_detect(string = s, pattern = "HAoSMC"), 
                    str_c(s, " (vascular SMC)", collapse = ""), 
                    ifelse(str_detect(string = s, pattern = "HUVEC"), 
                           str_c(s, " (vascular EC)", collapse = ""), 
                           y))
      return(res)
    })) %>%
    left_join(grp_colors, by = "cell_id") %>%
    mutate(color = ifelse(is.na(color), "#b6b6b6", color)) # grayish color for NA
  col_fun_annot3 <- cell_color_df$cell_individual_color # 1st group being closest to the labels on the right
  
  cell_names_char_vec <- cell_color_df$new_lbs
  names(col_fun_annot3) <- cell_names_char_vec
  # block_col_fun_annot3 <- unique(col_fun_annot3)
  
  if (plot_cluster_bar) {
    # top annotations
    top_ha <- HeatmapAnnotation(
      # set up the colored bars above heatmap
      `Cluster` = anno_block(gp = gpar(
        fill = block_col_fun_annot1,
        labels = group_order,
        labels_gp = gpar(col = "black", font = 2) #, fontsize = FONTSIZE + 2)
      )),
      `Cell line` = cell_names_char_vec,
      `Perturbation` = perturbations_char_vec,
      # color parameters
      col = list(
        `Cluster` = col_fun_annot1,
        `Perturbation` = col_fun_annot2,
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
      annotation_legend_param = list(`Perturbation` = list(title = "Perturbation"))
    ) 
  } else {
    # top annotations
    top_ha <- HeatmapAnnotation(
      # set up the colored bars above heatmap
      `Cell line` = cell_names_char_vec,
      `Perturbation` = perturbations_char_vec,
      # color parameters
      col = list(
        `Perturbation` = col_fun_annot2,
        `Cell line` = col_fun_annot3
      ),
      gp = gpar(col = "black"), # this draws black lines in the perturbation annotation block
      annotation_name_gp = gpar(font = 2), #, fontsize = FONTSIZE + 2),
      gap = unit(1, "mm"),
      simple_anno_size = unit(0.5, "cm"),
      show_annotation_name = c(TRUE, TRUE),
      show_legend = c(FALSE, FALSE),
      annotation_name_side = "left",
      annotation_name_align = TRUE,
      which = "column",
      
      # only need this for simple annotations, like annot_block
      annotation_legend_param = list(`Perturbation` = list(title = "Perturbation"))
    ) 
  }
  
  
  # save(list = ls(all.names = TRUE),
  #            file = "debug/debug_dat/debug-p100.RData")
  # stop()
  
  # get top up and down
  # needs to be relative to 1
  top_up <- row_annots_df %>%
    filter(fc > 1 & signif) %>%
    arrange(desc(`fc`))
  top_up
  top_down <- row_annots_df %>%
    filter(fc < 1 & signif) %>%
    # filter(logFC < LOGFC_CUTOFF) %>%
    arrange(desc(`fc`))
  top_down
  row_order_df <- bind_rows(top_up, top_down)
  row_order_df
  
  matrix_analyte_check <- sort(rownames(filtered_test_mat))
  row_order_analyte_check <- sort(row_order_df$analyte)
  
  setdiff(matrix_analyte_check, row_order_analyte_check)
  setdiff(row_order_analyte_check, matrix_analyte_check)
  
  stopifnot(nrow(row_order_df) == nrow(filtered_test_mat))
  
  # annotations for right side (by row)
  # don't use the boot bh pvalue, use the regular - the boot one is rounded
  signif_log_q_for_bar_plot <- row_order_df$neg_log10_p_val_bh
  # signif_log_q <- row_order_df %>% mutate(`-log10(q)` = round(-log(p_val_bh, base = 10),5) ) %>% .$`-log10(q)`
  signif_q_val <- format(row_order_df$p_val_bh, digits = 3, scientific = TRUE)
  extra_analyte_info_init <- row_order_df$mark
  
  # edit analyte information when talking about GCP
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
  fc_vec <- round(2^row_order_df$logFC, 3) # logFC is median log_2
  
  color_fc_row_labels_df <- row_order_df %>%
    mutate(color = ifelse(fc < 1, "blue", "red")) %>%
    mutate(border_color = ifelse(signif_and_fold, "white", "white"), # can switch to light orange #FFD580
           border_thickness = 1.9) # ifelse(signif_and_fold, 1.9, 1))
  color_fc_row_labels <- color_fc_row_labels_df  %>% .$color
  color_fc_border_row_labels <- color_fc_row_labels_df %>% .$border_color
  color_fc_border_row_thickness <- color_fc_row_labels_df %>% .$border_thickness
  
  # change the way the analyte is represented
  analytes_reordered <- row_order_df$analyte
  
  if (which_dat == "gcp") {
    analytes_reordered_labels <- row_order_df %>% 
      left_join(mod_table, by = c("analyte" = "full_sites")) %>%
      mutate(new_analyte_label = mod_counts) %>%
      .$new_analyte_label
  } else {
    analytes_reordered_labels <- str_split(string = gsub(x = row_order_df$analyte, 
                                                         pattern = "_[0-9]", replacement = "*"), 
                                           pattern = " ",simplify = TRUE)[,2] 
  }
  
  
  if (add_D_stat & which_dat == "gcp") {
    # add d-stat, no site
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
                            just = "center", location = 0.5, # location = 0.1;  unit(1, "npc")
                            gp = gpar(font = 2,  border = "black",
                                      fill = color_fc_border_row_labels, 
                                      lwd = color_fc_border_row_thickness), #, fontsize = FONTSIZE + 2),
                            width = max_text_width(analytes_reordered_labels) * 1.075
      ),
      which = "row",
      gap = unit(1.5, "mm") # gap between annotations
    )
  } else if (add_D_stat & which_dat != "gcp") {
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
                            just = "center", location = 0.5, # location = 0.1;  unit(1, "npc")
                            gp = gpar(font = 2,  border = "black",
                                      fill = color_fc_border_row_labels, 
                                      lwd = color_fc_border_row_thickness), #, fontsize = FONTSIZE + 2),
                            width = max_text_width(analytes_reordered_labels) * 1.075
      ),
      which = "row",
      gap = unit(1.5, "mm") # gap between annotations
    )
  } else if (!add_D_stat & which_dat == "gcp") {
    # no d_stat, no site
    right_ha <- HeatmapAnnotation(
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
      
      `analyte` = anno_text(analytes_reordered_labels,
                            just = "center", location = 0.5, # location = 0.1;  unit(1, "npc")
                            gp = gpar(font = 2,  border = "black",
                                      fill = color_fc_border_row_labels, 
                                      lwd = color_fc_border_row_thickness), #, fontsize = FONTSIZE + 2),
                            width = max_text_width(analytes_reordered_labels) * 1.075
      ),
      which = "row",
      gap = unit(1.5, "mm") # gap between annotations
    )
  } else {
    # no d_stat, add site
    right_ha <- HeatmapAnnotation(
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
                            just = "center", location = 0.5, # location = 0.1;  unit(1, "npc")
                            gp = gpar(font = 2, border = "black",
                                      fill = color_fc_border_row_labels, 
                                      lwd = color_fc_border_row_thickness), #, fontsize = FONTSIZE + 2),
                            width = max_text_width(analytes_reordered_labels) * 1.075
      ),
      which = "row",
      gap = unit(1.5, "mm") # gap between annotations
    )
  }
  
  
  
  # reorder the matrix with this -- this could also be done in the heatmap function;
  # however, be careful about the variable you use to assign column order and column attributes!
  # ESPECIALLY, column split
  dim(reduced_mat)
  length(analytes_reordered)
  length(grp_color_df$unique_id)
  
  ### ORDER HERE!!! ###
  reordered_mat <- reduced_mat[analytes_reordered, grp_color_df$unique_id, drop = F]

  ### insert code for left-annotation here if you want to use it ###
  
  ht_opt$message = FALSE
  ht_opt$COLUMN_ANNO_PADDING <- unit(2, "mm")
  ht_opt$ROW_ANNO_PADDING <- unit(1, "mm")
  ht_opt$annotation_border <- TRUE
  
  # this IS the same as matrix_for_diffe, which is what differential analytes is computed on
  # t2_reduced_long_ <- reduced_long %>% 
  #   filter(analytes == "pS12 EIF4A3") %>%
  #   arrange(cell_id) %>%
  #   pluck("value")
  
  # check incoming matrix
  
  ht <- Heatmap(reordered_mat,
                name = NAME,
                col = col_fun,
                na_col = "gray",
                top_annotation = top_ha,
                right_annotation = right_ha,
                
                # left_annotation = left_ha,
                
                #' @note turn this on to draw in the values of each cell!
                # cell_fun = function(j, i, x, y, width, height, fill) {
                #   grid.text(sprintf("%.1f", reordered_mat[i, j]), x, y, gp = gpar(fontsize = 5))
                # },
                
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
  
  draw(ht)
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
                     title = "Perturbation", 
                     legend_gp = gpar(fill = unique(col_fun_annot2)),  
                     direction = "horizontal", 
                     ncol = 2)
  
  # lgd_box <- Legend(labels = c(base_cell_id, "All other clusters"), title = "Boxplots",
  #                   legend_gp = gpar(fill = c(base_grp_color, rest_cell_color)))
  if (plot_cluster_bar){
    lgd_lst <- packLegend(lgd_cluster, lgd_cell, lgd_pert, # lgd_box,
                          max_width = unit(20, "cm"), 
                          direction = "horizontal")
  }  else {
    lgd_lst <- packLegend(lgd_cell, lgd_pert, # lgd_box,
                          max_width = unit(20, "cm"), 
                          direction = "horizontal")
  }
  
  
  
  init_ht <- draw(ht,
                  heatmap_legend_side = "bottom",
                  annotation_legend_list = lgd_lst,
                  annotation_legend_side = "left",
                  column_title = TITLE, 
                  column_title_gp = gpar(font = 2),
                  merge_legend = TRUE)
  
  size <- calc_ht_size(init_ht)
  
  
  ## Save the plot
  # in pixels
  linear_height <- size[2] + 2
  linear_height # 9
  linear_width <- size[1] + 5
  linear_width # 16
  
  # get a new filename for pdf outoput
  heatmap_output_fn_pdf <- file.path(str_c(str_split(heatmap_output_fn, 
                                                     pattern = "\\.", 
                                                     simplify = TRUE)[,1], ".pdf", sep = ""))
  
  to_plot <- function() {
    
    draw(init_ht)
    right_annots_height <- unit(1.01, "npc") + unit(3, "mm")
    height_top_dec <- convertX(ComplexHeatmap:::height(top_ha), "mm", valueOnly = TRUE)
    
    # add the annotation titles
    if (which_dat != "gcp") {
      # remove this for GCP, not helpful
      decorate_annotation("site", {
        grid.text("Site",
                  y = right_annots_height, just = "top",
                  gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
        )
      })
    }
    
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
      grid.text("mFC",
                y = right_annots_height, just = "top",
                gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
      )
    })
    
    if (add_D_stat) {
      decorate_annotation("d_statistic", {
        grid.text("D-stat",
                  y = right_annots_height, just = "top",
                  gp = gpar(font = 2) #, fontsize = FONTSIZE + 2)
        )
      })
    }
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
