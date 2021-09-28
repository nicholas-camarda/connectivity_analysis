

#' @note helper function for compute_boot_pvclust and morpheus plotter
#' @param long_df to be converted to matrix -- [assumes long df is a connectivity df!]
#' @return numeric matrix, symmetric
make_numeric_mat <- function(long_df) {
  if ("p_val" %in% colnames(long_df)) {
    long_df <- long_df %>% dplyr::select(-p_val)
  }
  x <- long_df %>%
    pivot_wider(names_from = group_b, 
                values_from = conn, 
                values_fn = median) %>%
    as.data.frame() %>%
    dplyr::select(-group_a) %>%
    as.matrix()
  rownames(x) <- colnames(x)
  return(x)
}

run_corr_lst <- function(lst, tie_method = "average") {
  message("Computing correlation...")
  res <- map(lst, run_corr, tie_method = tie_method)
  message("Done.")
  return(res)
}

#' @note run correlation
#' @param x long_form df from melt_gct(parse_gctx(.x))
#' @param tie_method for computing the ranks, whether to do average or dense
#' @return list obj containing matrix / long tbl forms of correlation result
run_corr <- function(x, tie_method = "average") {
  # x <- full_splt_lst[[1]]
  x_ungrouped <- x %>%
    ungroup()
  # 
  # stop('debug')
  
  wide_x <- x_ungrouped %>%
    dplyr::select(replicate_id, pr_gene_symbol, value) %>%
    pivot_wider(names_from = pr_gene_symbol, 
                values_from = value, 
                values_fn = median); wide_x
  
  transposed_x <- wide_x %>%
    dplyr::select(-replicate_id) %>%
    t() %>%
    # unnaming the large matrix results
    # in massive performance upgrade
    unname()
  
  # data has been log transformed, so must be spearman bc rho doesn't 
  # change with monotonic transfromation; pearson DOES!
  res <- cor(transposed_x, method = "spearman", use = "pairwise.complete.obs")
  
  rownames(res) <- wide_x$replicate_id
  colnames(res) <- wide_x$replicate_id
  
  res_final1 <- res %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble()
  
  colnames(res_final1) <- c("unique_id_a", wide_x$replicate_id)
  
  look_up_df <- x_ungrouped %>%
    distinct(master_id, replicate_id)
  
  res_final2 <- res_final1 %>%
    pivot_longer(
      cols = colnames(res_final1)[-1],
      names_to = "unique_id_b", values_to = "corr"
    ) %>%
    right_join(look_up_df %>%
                 dplyr::rename(group_a = master_id, unique_id_a = replicate_id),
               by = "unique_id_a"
    ) %>%
    right_join(look_up_df %>%
                 dplyr::rename(group_b = master_id, unique_id_b = replicate_id),
               by = "unique_id_b"
    )
  return(list("matrix" = res, "tibble" = res_final2))
}

#' @note get unique names from matrix
#' @param M matrix with row and column names
get_grp_names_from_matrix <- function(M, unique_cs, sep_ = "--") {
  grp_names <- str_split(unique_cs, pattern = sep_, simplify = T)[, 1]
  return(grp_names)
}

#' @note collapse connectivity between pairs by taking the median
#' @NOTE you need to get a unique ID for the 3+ replicates!!
#' @param cell_conn [matrix]!! of connectivities among cells
#' @return symmetric matrix of connectivities
calc_conn_by_median <- function(m) {
  if (!is.matrix(m)) {
    m <- make_numeric_mat(m)
  }
  
  comps <- length(colnames(m))
  # here conn is really med_conn
  df <- tibble() # p_val/conn for compatability issues
  for (i in 1:comps) {
    for (j in 1:comps) {
      n1 <- colnames(m)[i]
      n2 <- colnames(m)[j]
      v1 <- m[i, j]
      v2 <- m[j, i]
      median_val <- median(c(v1, v2))
      new_df <- tibble(group_a = n1, group_b = n2, conn = median_val)
      df <- bind_rows(df, new_df)
    }
  }
  symm_med_m <- make_numeric_mat(df)
  return(symm_med_m)
}

#' @note compute for cell A vs cell B and cell B vs all others (not B)
#' @param a character identifier for 'group a'
#' @param b character identifier for 'group b'
#' @param use_bootstrap logical to calculate boostrapped p-value for the ks.test
calc_conn <- function(a, b, my_mat, my_tib, use_bootstrap = FALSE) {
  # a <- looper$group_a[5];
  # b <- looper$group_b[5]
  # my_mat = mat;  my_tib = cs_tib
  # message(a,b)
  # A vs B
  grp_col_idx_a <- my_tib %>% filter(grp_names == a) %>% .$cs_idx
  grp_col_idx_b <- my_tib %>% filter(grp_names == b) %>% .$cs_idx
  
  # B vs all (except itself)
  grp_col_idx_b <- my_tib %>% filter(grp_names == b) %>% .$cs_idx
  grp_col_idx_all <- my_tib %>% filter(!(grp_names %in% c(a, b))) %>% .$cs_idx
  
  
  test <- my_mat[grp_col_idx_a, grp_col_idx_b]
  bkg <- my_mat[grp_col_idx_b, grp_col_idx_all]
  
  if (use_bootstrap) {
    message("Using bootstrap...")
    res <- ks.boot(test, bkg,
                   nboots = 1000, alternative = "two.sided"
    )
    connectivity <- res$ks$statistic
    p_val <- res$ks.boot.pvalue
  } else {
    res <- ks.test(test, bkg, alternative = "two.sided")
    connectivity <- res$statistic
    p_val <- res$p.value
  }
  
  if (median(test, na.rm = T) < median(bkg, na.rm = T)) {
    connectivity <- -1 * connectivity
  }
  
  conn_sub_obj <- tibble(conn = connectivity, p_val = p_val)
  return(conn_sub_obj)
}

#' @note run connectivity analysis using the helper function [calc_conn]
#' @param m matrix of correlation values, must be numeric and symmetric
#' @param use_bootstrap logical to calculate boostrapped p-value for the ks.test
#' @return connectivity tbl
run_conn <- function(mat, use_bootstrap = FALSE) {
  # mat <- sample_corr_lst$Epigenetic$matrix
  if (!is.matrix(mat)) mat <- make_numeric_mat(mat)
  
  unique_cs <- unique(rownames(mat))
  grp_names <- get_grp_names_from_matrix(mat, unique_cs, sep_ = "--")
  
  cs_idx <- seq_len(length(rownames(mat)))
  cs_tib <- tibble(unique_cs, cs_idx, grp_names)
  looper <- expand.grid(
    group_a = unique(grp_names),
    group_b = unique(grp_names)
  )
  
  if (use_bootstrap) message("Using ks.boot")
  conn_df <- looper %>%
    mutate(map2_df(
      .x = group_a, .y = group_b,
      .f = calc_conn,
      my_mat = mat, my_tib = cs_tib, 
      use_bootstrap = use_bootstrap
    )) %>%
    as_tibble()
  
  median_mat <- calc_conn_by_median(conn_df)
  rnames <- rownames(median_mat)
  # stop("debug")
  conn_tbl_median <- median_mat %>% 
    as.data.frame() %>%
    rownames_to_column(var = "group_a") %>% 
    as_tibble() %>%
    pivot_longer(cols = all_of(rnames), 
                 names_to = "group_b", values_to = "median_conn")
  
  
  return(list("conn_tbl" = conn_df, 
              "conn_median_mat" = median_mat,
              "conn_tbl_median" = conn_tbl_median))
}

#' @note runn connectivity analysis over a list of matrices
#' @param lst list of correlation matrices
#' @return connectivity long_df
run_conn_lst <- function(lst, use_bootstrap = FALSE) {
  # global assign
  message("Computing connectivity...")
  # need the matrix form, not tibble --> transpose
  mats <- purrr::transpose(lst)$matrix
  connectivity_lst_obj <- map(mats, run_conn, use_bootstrap = use_bootstrap)
  message("Done.")
  return(connectivity_lst_obj)
}


#' @note run clustering analysis over list of matrices
#' @param lst list of conn_lst obj
#' @return clustering object results
run_clust_lst <- function(lst, use_parallel = FALSE) {
  message("Computing clusters...")
  mat <- purrr::transpose(lst)$conn_median_mat 
  clust_lst_obj <- map(mat, run_clust, use_parallel = use_parallel)
  message("Done.")
  return(clust_lst_obj)
}

#' @note compute connectivity clustering analysis on conn matrix
#' @param mat connectivity matrix
run_clust <- function(mat, use_parallel = FALSE) {
  # mat <- sample_conn_obj$conn_lst[[1]]$conn_median_mat
  mat <- force_natural(mat)
  
  # stop("debug")
  # print(Heatmap(mat))
  
  pv_obj <- compute_boot_pvclust(
    x = mat,
    parallel_flag = use_parallel, 
    n_boot = 1000
  )
  
  # dendro_cut_thresh is global
  named_clusters <- generate_clusters(
    x = pv_obj,
    thresh = dendro_cut_thresh
  )
  
  co_clust_bool <- co_cluster(named_clusters)
  
  res <- list(
    "clust_obj" = pv_obj,
    "cluster_assignments" = named_clusters,
    "co_clust_bool" = co_clust_bool
  )
  return(res)
}

#' @note compute clusters via Lev's method
#' @param x numeric matrix for clustering
#' @param dname names of drugs (groups) -- only used for naming the result plot
#' @param base_output_dir name of base output dir for plots
compute_boot_pvclust <- function(x, parallel_flag = FALSE, n_boot = 1000) {
  if (!(nrow(x) > 2)) {
    message("Not enough data to perform clustering analysis.")
    return(NULL)
  }
  
  # stop("debug")
  res <- pvclust(x,
                 method.dist = "correlation",
                 method.hclust = "average", 
                 # weight = TRUE,
                 quiet = TRUE,
                 nboot = n_boot,
                 iseed = 2334,
                 parallel = parallel_flag
  )
  
  # db <- fpc::dbscan(x, eps = 0.35, MinPts = 5)
  # plot(db, x, main = "DBSCAN", frame = FALSE)
  # dbscan::kNNdistplot(x, k =  5)
  # abline(h = 0.3, lty = 2)
  # plot(res)
  # # msplot(res)
  # pvpick(res)
  
  return(res)
}

#' @note extract cluster assignments but cutting dendrograms at desired fraction of max tree height
#' @param x pvclust object
#' @param thresh percent of max height of dendrogram to cut
#' @return cluster assignments
generate_clusters <- function(x, thresh = 0.6) {
  if (is.null(x)) {
    return(NULL)
  }
  tree <- x$hclust
  cut_height <- thresh * max(tree$height)
  clusters <- cutree(tree, h = cut_height)
  return(clusters)
}

#' @note return whether or not, e.g. HUVECs and SMCs, clustered together
#' @param cut_tree_obj named vector of cluster assignments
#' @param name1 first var you want to compare
#' @param name2 second var you want to compare
co_cluster <- function(cut_tree_obj, name1 = "HAoSMC", name2 = "HUVEC") {
  if (is.null(cut_tree_obj)) {
    return(NA)
  }
  names_split <- map_chr(names(cut_tree_obj), function(ct) str_split(string = ct, pattern = "--", simplify = T)[, 1])
  name1_obj_idx <- which(name1 == names_split)
  name2_obj_idx <- which(name2 == names_split)
  c1 <- cut_tree_obj[name1_obj_idx]
  c2 <- cut_tree_obj[name2_obj_idx]
  ni <- intersect(c1, c2)
  if (length(ni) != 1) {
    return(FALSE)
  }
  res <- unname(c1 == c2 & table(cut_tree_obj)[ni] == 2)
  return(res)
}


run_diffe_lst <- function(lst, clust_lst) {
  message("\n")
  diffe_lst_obj <- purrr::pmap(list(lst, clust_lst, names(lst)), 
                        .f = run_diffe)
  message("Done.")
  return(diffe_lst_obj)
}

#' @note run differential expression on dat
#' @param dat tbl containing raw data
#' @param cob list obj containing results of clustering analysis
#' @param dname grouping_var, e.g. Kinase inhibitor or dmso, for naming the plot
run_diffe <- function(dat, cob, dname) {
  # dat <- lst$`Kinase inhibitor`; cob <- clust_lst$`Kinase inhibitor`; dname <- "Kinase inhibitor"; 
  
  # base_output_dir <- lst$data
  # which_dat <- "P100"
  
  
  # stop()
  which_dat <- unique(force_natural(dat$which_dat))
  message(qq("Computting differential expression on @{dname}, @{which_dat}"))
  
  mat_tbl <- dat %>% 
    dplyr::select(replicate_id, cell_id, pert_iname, pert_class, pr_gene_symbol, value) %>%
    ungroup() 
  
  clust_assignments <- cob$cluster_assignments
  # join clusters and annotations
  ca_df_temp <- create_ca_df(ca = clust_assignments) %>%
    rename(base_clust_comp = cluster); ca_df_temp
  
  # # which cluster(s) are vascular cells in
  # which_cluster_vasc_df <- ca_df_temp %>% 
  #   filter(cell_id %in% c("HUVEC", "HAoSMC"))
  
  # get cluster labels
  clust_label_df <- ca_df_temp %>% 
    group_by(base_clust_comp) %>%
    summarize(base_clust_comp_name = str_c(sort(unique(cell_id)), 
                                           collapse = ","), 
              .groups = "keep") %>%
    ungroup() %>%
    arrange(base_clust_comp) 
  
  ca_df <- left_join(ca_df_temp, clust_label_df, by="base_clust_comp")
  
  dat_ <- dat %>% 
    left_join(ca_df, by=c("cell_id")) 
  
  clusts_int_vec <- clust_label_df$base_clust_comp
  
  
  matrix_for_diffe <- mat_tbl %>%
    left_join(ca_df, by = "cell_id") %>%
    dplyr::select(replicate_id, cell_id, pert_iname, pert_class, 
                  base_clust_comp, base_clust_comp_name, value, pr_gene_symbol) %>%
    pivot_wider(names_from = pr_gene_symbol, 
                values_from = value, 
                values_fn = median) # need this in case dups with different values
  
  feature_names <- dat_ %>% ungroup() %>% .$pr_gene_symbol %>% unique()
  
  # print(dat_)
  # print(cob)
  # print(dname)
  # stop("help")
  p2 <- progressr::progressor(steps = length(clusts_int_vec)*length(feature_names) )
  message("Computing differential analytes..")
  diffe_by_clust_df <- map_df(clusts_int_vec, function(ith_cluster){
    # ith_cluster <- clusts_int_vec[1]
    
    base_clust_comp_name <- clust_label_df %>% 
      dplyr::select(base_clust_comp, base_clust_comp_name) %>% # take the correct cluster naming convention column
      filter(base_clust_comp == ith_cluster) %>% # take the correct name according to the cluster integer
      distinct(base_clust_comp_name) %>% # non-vascular clusters need to get lumped together, we just want the unique cluster assign
      pluck(1) # make it a character vector 
    
    message(str_c("\nCluster ID:",ith_cluster, "\n#", base_clust_comp_name, sep =" "))
    p2()
    
    cluster_res <- map2_df(feature_names, ith_cluster, function(k, i){
      # k <- feature_names[1]; i <- ith_cluster[1]
      # 
      # DEBUG: 
      # message(k)
      
      c1 <- matrix_for_diffe %>% filter(base_clust_comp == i) %>% pluck(k)
      c2 <- matrix_for_diffe %>% filter(base_clust_comp != i) %>% pluck(k)
      
      # dens_dat <- tibble(c = c(c1, c2), g = c(rep("c1", length(c1)), rep("c2", length(c2))))
      # ggecdf(data = dens_dat, x = "c", color = "g", ggtheme = theme_bw(), palette = "jco")
      # 
      n_non_na_vec1 <- length(c1[!is.na(c1)]); n_non_na_vec1
      n_non_na_vec2 <- length(c2[!is.na(c2)]); n_non_na_vec2
      
      if (n_non_na_vec1 < 1 | n_non_na_vec2 < 1) {
        # message("Not enough data to compute diffe for ", k)
        res <- tibble(analyte = k, 
                      logFC = NA, base_clust_comp = as.integer(i), 
                      base_clust_comp_name = base_clust_comp_name, 
                      ks_statistic = NA, 
                      ks_boot_statistic = NA, 
                      signal_to_noise_statistic = NA,
                      directional_stat = NA, 
                      p_val = NA, p_val_boot = NA,
                      p_val_bh = NA,
                      p_val_boot_bh = NA, 
                      signif = NA) %>%
          mutate(k_clust_dat = list(c1), 
                 all_others_dat = list(c2))
        
        # cat("*")
        return(res)
      }
      
      ks <- ks.test(c1, c2, alternative = "two.sided", tol = 1e-8); ks
      ks_boot <- ks.boot(c1, c2, nboots = 1000, alternative = "two.sided"); ks_boot
      
      
      # LINCS espouses the concept of making different data levels available for public use.  Different data levels correspond different steps along our processing workflow.  The LINCS PCCSE levels are defined as follows:
      # Level 0 - Raw Mass Spectrometry Data (LCMS) - will be available through a chorusproject.org repository in the future
      # Level 1 - Probe Reads (SKY) - Curated Skyline documents; available on this website, including metadata
      # Level 2 - Raw Numerical Data (RPT) - Matrix data of extracted signal ratios of endogenous probes vs. internal standards (log2 transformed); available on this website, including metadata
      # Level 3 - Normalized and QC'ed Numerical Data (QCNORM) - Matrix data derived from Level 2 after automated processing and normalization
      # Level 4 - Differential Quantification (DIFF) - Matrix data of Level 3 with plate-wide median ratio of each analyte subtracted from each sample; available on this website, including metadata 
      
      # compute ks stat
      #quantifies the distance between the empirical distribution of given two samples
      stat <- ks$statistic
      stat_boot <- ks_boot$ks$statistic
      
      c1_median <- median(c1, na.rm = T); c1_median
      c2_median <- median(c2, na.rm = T); c2_median
      
      #' @note assumes data are already log transformed!!! 
      #' So a fold change is simply a difference, not a division!!
      #' is this the correct way to calculate diffex??
      #' 
      logfc <- median(c1, na.rm = T) -  median(c2, na.rm = T) ; logfc
      d_stat <- ifelse(logfc > 0, "++", "--")
      
      min_x <- min(c(min(c1, na.rm = T), min(c2, na.rm = T)))
      max_x <- max(c(max(c1, na.rm = T), max(c2, na.rm = T)))
      
      c1_names_df <- matrix_for_diffe %>% 
        filter(base_clust_comp == i) %>% 
        mutate(id = str_c(base_clust_comp_name, pert_iname, sep = "-"),
               plot_clust_id = 1) %>%
        dplyr::select(plot_clust_id, id, all_of(k)) %>%
        rename(value := !!sym(k)); c1_names_df
      
      c2_names_df <- matrix_for_diffe %>% 
        filter(base_clust_comp != i) %>% 
        mutate(id = str_c(base_clust_comp_name, pert_iname, sep = "-"),
               plot_clust_id = 2) %>%
        dplyr::select(plot_clust_id, id, all_of(k)) %>%
        rename(value := !!sym(k)); c2_names_df
      
      # compute signal to noise
      mean_grp_diff <- mean(c1, na.rm = T) - mean(c2, na.rm = T)
      sum_grp_sd <- sd(c2, na.rm = T) + sd(c1, na.rm = T)
      signal_to_nose <- mean_grp_diff / sum_grp_sd
      
      #'@note I'm pretty sure we want to adjust the p-values per cluster comparison... 
      res <- tibble(analyte = k,  logFC = logfc, base_clust_comp = i, 
                    base_clust_comp_name = base_clust_comp_name, 
                    ks_statistic = stat, 
                    ks_boot_statistic = stat_boot, 
                    signal_to_noise_statistic = signal_to_nose,
                    directional_stat = d_stat, 
                    p_val = ks$p.value, p_val_boot = ks_boot$ks.boot.pvalue,
                    p_val_bh = p.adjust(ks$p.value, method = "BH"),
                    p_val_boot_bh = p.adjust(ks_boot$ks.boot.pvalue, method = "BH")) %>%
        mutate(signif = p_val_boot_bh < bh_thresh_val) %>%
        mutate(k_clust_dat = list(c1), 
               all_others_dat = list(c2))
      
      # if (k == "CLTA" & base_clust_comp_name == "A375,HUVEC,PC3") stop("debug")
      # na_omit_res <- res %>% na.omit()
      p2()
      # cat(".")
      return(res)
    }) 
    # cat(".Done\n")
    
    return(cluster_res)
  }) 
  
  stopifnot(nrow(diffe_by_clust_df) > 0)
  
  phosphosite_meta <- dat_ %>% 
    ungroup() %>%
    dplyr::distinct( pr_gene_symbol, pr_gene_id, mark) %>%
    rename(analyte = pr_gene_symbol)
  
  diffe_final_res <- diffe_by_clust_df  %>%
    left_join(phosphosite_meta, by = "analyte") %>%
    dplyr::select(analyte, pr_gene_id, mark, everything()) %>%
    mutate(label_ = str_c("p",str_c(mark, analyte,sep =" ")), sep="") %>%
    mutate(signif_and_fold = ifelse(signif & logFC > LOGFC_CUTOFF, T, ifelse(signif & logFC < -LOGFC_CUTOFF, T, F))) %>% 
    dplyr::select(signif_and_fold, signif, logFC, everything()) %>%
    mutate(fc = 2^logFC); diffe_final_res

  signif_df <- diffe_final_res %>% filter(signif_and_fold); nrow(signif_df)
  # save this plot
  
  # stop()
  diffe_g <- ggplot(diffe_final_res %>% na.omit()) + 
    geom_point(aes(x = fc, y = -log10(p_val_boot_bh), color = directional_stat), size = 2) +
    geom_vline(xintercept=c(2^-LOGFC_CUTOFF, 2^LOGFC_CUTOFF), col="purple", alpha = 0.5) +
    geom_hline(yintercept=-log10(0.1), col="red", alpha = 0.5) +
    geom_label_repel(signif_df,
                     mapping = aes(x = fc, y = -log10(p_val_boot_bh), label = label_), force = 5,
                     min.segment.length = unit(0, 'lines')) +
    scale_color_manual(values = c("blue","red")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 10, angle = 45)) +
    guides(size = FALSE) + 
    labs(color = "Directional change", 
         caption = "BH.q.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1 \nDifference between groups calculated using ks.test()") +
    ylab("-Log10(BH.q.val)") +
    xlab("Fold change") + 
    ggtitle(label = qq("@{toupper(which_dat)}, @{dname}"))
  
  if (nrow(clust_label_df) > 1) {
    diffe_g <- diffe_g +
      facet_grid(~base_clust_comp_name, scales = "free_x") 
  }
  diffe_g
  # message("Done")

  return(list(diffe_final_res, diffe_g) %>% 
           setNames(c("diffe_final_res", "diffe_ggplot")))
}

#' @note run connectivity and clustering analysis, with progress updates
#' @param lst a list of tbl's split by the grouping structure for the given analysis
#' @return result containing the input data, corr_lst (correlation results),
#' conn_lst (connectivity results), and clust_lst (clustering results)
run_analysis <- function(lst, tie_method = "average",
                         use_bootstrap = FALSE, use_parallel = FALSE) {
  p <- progressr::progressor(steps = 3)
  p(message = "Computing correlation")
  corr_lst <- run_corr_lst(lst, tie_method = tie_method)
  
  p(message = "Computing connectivity")
  conn_lst <- run_conn_lst(corr_lst, use_bootstrap = use_bootstrap)
  
  p(message = "Clustering")
  clust_lst <- run_clust_lst(conn_lst, use_parallel = use_parallel)
  
  # message("Computing differential expression")
  diffe_lst <- run_diffe_lst(lst, clust_lst)
  message("Done")
  
  res <- list(
    "input_data" = lst, 
    "output_results" = list("corr_lst" = corr_lst,
                            "conn_lst" = conn_lst, 
                            "clust_lst" = clust_lst,
                            "diffe_lst" = diffe_lst)
  )
  return(res)
}



update_dataset_with_clustering <- function(data, clust_lst) {
  cell_ids <- names(clust_lst$cluster_assignments)
  clust_id_df <- tibble(cell_id = cell_ids, cut_trees = clust_lst$cluster_assignments )
  return(data %>% left_join(clust_id_df, by = "cell_id"))
}







#' #' @NOTE: this is SO SO slow
#' #' @note compute for cell A vs cell B and cell B vs all others (not A)
#' #' @param cell_corr matrix of correlations among cells
#' #' @param use_bootstrap logical to calculate boostrapped p-value for the ks.test
#' compute_connectivity <- function(M, use_bootstrap = FALSE) {
#'   
#'   # M <- corr_mat_reps %>% slice(1) %>% .$cell_corr_r_mat ; M <- M[[1]]
#'   
#'   unique_cs <- unique(rownames(M))
#'   grp_names <- get_grp_names_from_matrix(M, unique_cs)
#'   unique_cs_idx <- 1:length(rownames(M))
#'   
#'   cs_tib <- tibble(unique_cs, unique_cs_idx, grp_names)
#'   
#'   grp_idx <- bind_cols(tibble(
#'     grp_names = unique(grp_names),
#'     grp_idx = 1:length(unique(grp_names))
#'   ))
#'   
#'   cs_tib <- left_join(cs_tib, grp_idx, by = "grp_names")
#'   looper <- expand.grid(grp_namesA = unique(cs_tib$grp_names), grp_namesB = unique(cs_tib$grp_names))
#'   
#'   if (use_bootstrap) message("Using ks.boot")
#'   conn_df <- looper %>%
#'     mutate(map2_df(grp_namesA, grp_namesB, function(x, y) {
#'       # x <- looper$grp_namesA[5]; y <- looper$grp_namesB[5]
#'       
#'       grpA_idx1 <- cs_tib %>%
#'         filter(grp_names == x) %>%
#'         .$unique_cs_idx
#'       
#'       grpB_idx2 <- cs_tib %>%
#'         filter(grp_names == y) %>%
#'         .$unique_cs_idx
#'       
#'       test <- M[grpA_idx1, grpB_idx2]
#'       
#'       grpA_idx2 <- cs_tib %>%
#'         filter(grp_names == y) %>%
#'         .$unique_cs_idx
#'       
#'       grpB_idx2 <- cs_tib %>%
#'         filter(!(grp_names %in% y)) %>%
#'         .$unique_cs_idx
#'       
#'       bkg <- M[grpA_idx2, grpB_idx2]
#'       
#'       if (use_bootstrap) {
#'         res <- ks.boot(test, bkg, nboots = 1000, alternative = "two.sided")
#'         connectivity <- res$ks$statistic
#'         p.val <- res$ks.boot.pvalue
#'         
#'       } else {
#'         res <- ks.test(test, bkg, alternative = "two.sided")
#'         connectivity <- res$statistic
#'         p.val <- res$p.val
#'       }
#'       
#'       if (median(test, na.rm = T) < median(bkg, na.rm = T)) connectivity <- -1 * connectivity
#'       
#'       cat(".")
#'       return(tibble(conn = connectivity, p_val = p.val))
#'       
#'       
#'     })) %>%
#'     bind_rows()
#' }
#' 
# stop("debug here")
