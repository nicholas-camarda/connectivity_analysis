
#' @note load cached objects helper function
#' @param output_paths takes this dataframe of output paths and results mapping for reading
#' @param grouping_var takes in the grouping var and ensures correct directory placement
#' @return loaded cached objects
load_cached_objs <- function(output_paths, grouping_var) {
  message("Loading cached obj's...")
  res_paths_tbl_temp <- output_paths %>%
    mutate(match_temp = str_extract(
      string = path,
      pattern = "clust_clust_obj--cluster_assignments--co_clust_bool|conn_conn_tbl--conn_median_mat--conn_tbl_median|corr_matrix--tibble|diffe_diffe_final_res--signif_df"
    )) %>%
    na.omit() %>%
    mutate(match = str_c(str_extract(match_temp, "clust|conn|corr|diffe"), "_lst", sep = ""))
  
  res_paths_tbl <- res_paths_tbl_temp %>%
    mutate(lst = map(path, read_rds)) %>%
    dplyr::select(
      extract_path, path, !!grouping_var,
      match_temp, match, lst, everything()
    ) %>%
    mutate(match_temp = map_chr(
      extract_path,
      function(x) str_extract(string = x, pattern = "corr|conn|clust|diffe")
    )) %>%
    mutate(match = str_c(match_temp, "_lst")) %>%
    dplyr::select(match, lst, path, !!grouping_var)
  message("Done.")
  return(res_paths_tbl)
}



#' @note helper function for compute_boot_pvclust and morpheus plotter
#' @param long_df to be converted to matrix -- [assumes long df is a connectivity df!]
#' @return numeric matrix, symmetric
make_numeric_mat <- function(long_df) {
  if ("p_val" %in% colnames(long_df)) {
    long_df <- long_df %>% dplyr::select(-p_val)
  }
  x <- long_df %>%
    pivot_wider(
      names_from = group_b,
      values_from = conn,
      values_fn = median
    ) %>%
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
    pivot_wider(
      names_from = pr_gene_symbol,
      values_from = value,
      values_fn = median
    )
  wide_x
  
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
  grp_col_idx_a <- my_tib %>%
    filter(grp_names == a) %>%
    .$cs_idx
  grp_col_idx_b <- my_tib %>%
    filter(grp_names == b) %>%
    .$cs_idx
  
  # B vs all (except itself)
  grp_col_idx_b <- my_tib %>%
    filter(grp_names == b) %>%
    .$cs_idx
  grp_col_idx_all <- my_tib %>%
    filter(!(grp_names %in% c(a, b))) %>%
    .$cs_idx
  
  
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
    pivot_longer(
      cols = all_of(rnames),
      names_to = "group_b", values_to = "median_conn"
    )
  
  
  return(list(
    "conn_tbl" = conn_df,
    "conn_median_mat" = median_mat,
    "conn_tbl_median" = conn_tbl_median
  ))
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
  mat <- as.matrix(force_natural(mat))
  
  # stop("debug")
  print(Heatmap(mat))
  
  pv_obj <- compute_boot_pvclust(
    x = mat,
    parallel_flag = use_parallel,
    n_boot = 1000
  )
  
  print(plot(pv_obj))
  
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
  } else if (nrow(x) == 3) {
    message("Distance matrix doesn't have enough comparisons. Using Hclust, euclidean, average")
    x <- dist(x)
    res <- list(hclust = hclust(x, method = "average"))
  } else {
    message("Using pvclust..")
    # hclust(x)
    # stop("debug")
    res <- pvclust(x,
                   method.dist = "correlation",
                   method.hclust = "average",
                   # weight = TRUE,
                   quiet = FALSE,
                   nboot = n_boot,
                   iseed = 2334,
                   parallel = parallel_flag
    )
  }
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
                               .f = run_diffe
  )
  message("Done.")
  return(diffe_lst_obj)
}

#' @note run differential expression on dat
#' @param dat tbl containing raw data
#' @param cob list obj containing results of clustering analysis
#' @param dname grouping_var, e.g. Kinase inhibitor or dmso, for naming the plot
run_diffe <- function(dat, cob, dname) {
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-diffe.RData")
  # load("debug/debug_dat/debug-diffe.RData")
  # stop()
  
  which_dat <- unique(force_natural(dat$which_dat))
  message(qq("Computting differential expression on @{dname}, @{which_dat}"))
  
  clust_assignments <- cob$cluster_assignments; clust_assignments
  # join clusters and annotations
  ca_df_temp <- create_ca_df(ca = clust_assignments) %>%
    rename(base_clust_comp = cluster); ca_df_temp
  
  # get cluster labels
  clust_label_df <- ca_df_temp %>%
    group_by(base_clust_comp) %>%
    summarize(
      base_clust_comp_name = str_c(sort(unique(cell_id)),
                                   collapse = ","
      ),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    arrange(base_clust_comp);clust_label_df
  
  ca_df <- left_join(ca_df_temp, clust_label_df, by = "base_clust_comp"); ca_df
  
  # assign the clustering details to the original dataset, with a new name now
  dat_ <- dat %>%
    left_join(ca_df, by = c("cell_id"))
  
  mat_tbl <- dat_ %>%
    dplyr::select(
      replicate_id, cell_id, pert_iname, pert_class,
      base_clust_comp, base_clust_comp_name, value, pr_gene_symbol
    ) %>%
    ungroup(); mat_tbl
  
  # calculate the correct logFC's on mat_tbl
  #' @param df_ is the mat_tbl df
  #' @param i_clust is the ith cluster to be base_clust_comp, or the cluster
  #' we are comparing everything else to
  #' @return df containing logFC, fc, d_stat, and signif info
  calculate_logFC_grouped <- function(df_ = NA, i_clust = NA){
    # message(i_clust)
    median_df <- df_ %>%
      mutate(helper_cluster = ifelse(base_clust_comp == i_clust, 1, 0)) %>%
      arrange(cell_id); median_df
  
    # get some extra info
    extra_info_temp <- median_df %>%
      group_by(helper_cluster, cell_id, pr_gene_symbol) %>%
      summarize(gene_clust_median_val = median(value, na.rm = TRUE), .groups = "keep")
    
    # wide summarized info
    # this is the median value of each cell_ID!
    summarized_median_df_extr_info <- extra_info_temp %>%
      pivot_wider(id_cols = pr_gene_symbol, 
                  names_from = c(cell_id, helper_cluster), names_sep = "__", 
                  values_from = gene_clust_median_val) %>%
      ungroup(); summarized_median_df_extr_info
    
    # then make final df, combining cluster-specific medians and cell_id specific medians
    summarized_median_df <- median_df %>%
      group_by(helper_cluster, cell_id, pr_gene_symbol) %>%
      summarize(gene_clust_median_val = median(value, na.rm = TRUE), .groups = "keep") %>%
      pivot_wider(id_cols = pr_gene_symbol, 
                  names_from = helper_cluster, names_sep = "__", 
                  values_from = gene_clust_median_val,
                  values_fn = function(x) median(x, na.rm= TRUE),
                  names_prefix = "cluster_") %>%
      ungroup() %>%
      # join the extra info to main
      left_join(summarized_median_df_extr_info, by = "pr_gene_symbol") %>%
      rename(cluster_1_median = cluster_1, cluster_0_median = cluster_0) %>%
      # base cluster minus everything ELSE
      mutate(logFC = cluster_1_median - cluster_0_median,
             fc = 2^logFC) %>%
      rename(analyte = pr_gene_symbol) %>%
      mutate(directional_stat = ifelse(logFC > 0, "++", "--")) 
    # summarized_median_df
    
    # res <- list(logfc_df = summarized_median_df)
    return(summarized_median_df)
  }
  
  logFC_df <- clust_label_df %>%
    mutate(x = list(mat_tbl)) %>% # just rep a tbl as a row in another tbl
    mutate(res =  map2(.x = x, 
                       .y = base_clust_comp, 
                       .f = calculate_logFC_grouped)) %>%
    dplyr::select(-x) %>%
    unnest(c(res)) %>%
    arrange(base_clust_comp, analyte); 
  
  
  # then pivot for diffe ks.test calculations
  matrix_for_diffe <- mat_tbl %>%
    pivot_wider(
      names_from = pr_gene_symbol,
      values_from = value,
      values_fn = median # need this in case dups with different values, which there shouldn't be but...
    ) ; matrix_for_diffe
  
  feature_names <- dat_ %>%
    ungroup() %>%
    .$pr_gene_symbol %>%
    unique(); feature_names
  
  clusts_int_vec <- clust_label_df$base_clust_comp; clusts_int_vec
  p2 <- progressr::progressor(steps = length(clusts_int_vec) * length(feature_names))

  diffe_by_clust_df <- map_df(clusts_int_vec, function(ith_cluster) {
    # ith_cluster <- clusts_int_vec[2]
    
    base_clust_comp_name <- clust_label_df %>%
      filter(base_clust_comp == ith_cluster) %>% # take the correct name according to the cluster integer
      distinct(base_clust_comp_name) %>% # non-vascular clusters need to get lumped together, we just want the unique cluster assign
      pluck(1) ; base_clust_comp_name# make it a character vector
    
    message(str_c("\nCluster ID:", ith_cluster, "\n#", base_clust_comp_name, sep = " "))
    p2()
    
    # i want this result to be bound by rows! since each feature-cluster pair will be an 'observation'
    cluster_res <- future_map2_dfr(feature_names, ith_cluster, function(k, i) {
      # k <- feature_names[35]; i <- ith_cluster
      # DEBUG:
      # message(k)
      
      c1 <- matrix_for_diffe %>%
        filter(base_clust_comp == i) %>%
        pluck(k); c1
      
      c2 <- matrix_for_diffe %>%
        filter(base_clust_comp != i) %>%
        pluck(k); c2
      
      c1_length <- length(c1)
      c2_length <- length(c2)
      n_non_na_vec1 <- length(c1[!is.na(c1)]); n_non_na_vec1
      n_non_na_vec2 <- length(c2[!is.na(c2)]); n_non_na_vec2
      
      # if at least 1/3 of the results are missing, then quit
      # if (n_non_na_vec1/c1_length < (1/3) | n_non_na_vec2/c2_length < (1/3)) {
      if (n_non_na_vec1 < 1 | n_non_na_vec2 < 1) {
        # message("Not enough data to compute diffe for ", k)
        res <- tibble(
          analyte = k, 
          base_clust_comp = as.integer(i),
          base_clust_comp_name = base_clust_comp_name,
          ks_statistic = NA,
          ks_boot_statistic = NA,
          p_val = NA, p_val_boot = NA,
          p_val_bh = NA,
          p_val_boot_bh = NA,
          signif = NA
        ) %>%
          mutate(
            k_clust_dat = list(c1),
            all_others_dat = list(c2)
          )
        return(res)
      }
      
      ks <- ks.test(c1, c2, alternative = "two.sided")
      ks_boot <- ks.boot(c1, c2, nboots = 1000, alternative = "two.sided")
      
      # LINCS espouses the concept of making different data levels available for public use.  Different data levels correspond different steps along our processing workflow.  The LINCS PCCSE levels are defined as follows:
      # Level 0 - Raw Mass Spectrometry Data (LCMS) - will be available through a chorusproject.org repository in the future
      # Level 1 - Probe Reads (SKY) - Curated Skyline documents; available on this website, including metadata
      # Level 2 - Raw Numerical Data (RPT) - Matrix data of extracted signal ratios of endogenous probes vs. internal standards (log2 transformed); available on this website, including metadata
      # Level 3 - Normalized and QC'ed Numerical Data (QCNORM) - Matrix data derived from Level 2 after automated processing and normalization
      # Level 4 - Differential Quantification (DIFF) - Matrix data of Level 3 with plate-wide median ratio of each analyte subtracted from each sample; available on this website, including metadata
      
      # compute ks stat
      # quantifies the distance between the empirical distribution of given two samples
      stat <- ks$statistic
      stat_boot <- ks_boot$ks$statistic
  
      c1_names_df <- matrix_for_diffe %>%
        filter(base_clust_comp == i) %>%
        mutate(
          id = str_c(base_clust_comp_name, pert_iname, sep = "-"),
          plot_clust_id = 1
        ) %>%
        dplyr::select(plot_clust_id, id, all_of(k)) %>%
        rename(value := !!sym(k))
      c1_names_df
      
      c2_names_df <- matrix_for_diffe %>%
        filter(base_clust_comp != i) %>%
        mutate(
          id = str_c(base_clust_comp_name, pert_iname, sep = "-"),
          plot_clust_id = 2
        ) %>%
        dplyr::select(plot_clust_id, id, all_of(k)) %>%
        rename(value := !!sym(k))
      c2_names_df
      
      # compute signal to noise
      sum_grp_sd <- sd(c2, na.rm = T) + sd(c1, na.rm = T) # this is sd of log-transformed data
      # signal_to_nose <- mean_logfc / sum_grp_sd
      
      #' @note I'm pretty sure we want to adjust the p-values per cluster comparison...
      res <- tibble(
        analyte = k, 
        base_clust_comp = i,
        base_clust_comp_name = base_clust_comp_name,
        ks_statistic = stat,
        ks_boot_statistic = stat_boot,
        p_val = ks$p.value, p_val_boot = ks_boot$ks.boot.pvalue,
        p_val_bh = p.adjust(ks$p.value, method = "BH"),
        p_val_boot_bh = p.adjust(ks_boot$ks.boot.pvalue, method = "BH")
      ) %>%
        # determine significance of comparison based on bh_thresh
        mutate(signif = p_val_bh < bh_thresh_val) %>%
        mutate(
          k_clust_dat = list(c1),
          all_others_dat = list(c2)
        ); res
      
      p2()
      # cat(".")
      return(res)
    })
    # cat(".Done\n")
    
    return(cluster_res)
  })
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-gcp-epi-diffe.RData")
  # load("debug/debug_dat/debug-gcp-epi-diffe.RData")
  # stop()
  
  stopifnot(nrow(diffe_by_clust_df) > 0)
  
  phosphosite_meta <- dat_ %>%
    ungroup() %>%
    dplyr::distinct(pr_gene_symbol, pr_gene_id, mark) %>%
    rename(analyte = pr_gene_symbol)
  
  diffe_final_res_temp <- diffe_by_clust_df %>%
    left_join(phosphosite_meta, by = "analyte") %>%
    dplyr::select(analyte, pr_gene_id, mark, everything()) %>%
    # so that log doesn't cause inf
    # mutate(p_val_boot_bh = p_val_boot_bh + .Machine$double.xmin) %>%
    mutate(neg_log10_p_val_bh = -log10(p_val_bh),
           neg_log10_p_val_boot_bh = -log10(p_val_boot_bh)) %>%
    # add in the correctly calculated logFC 
    left_join(logFC_df, by = c("analyte", "base_clust_comp", "base_clust_comp_name")) %>%
    # use the p-value significance and fold change to compute signif diff and fold change
    mutate(signif_and_fold = ifelse(signif & (fc >= FC_UPPER_BOUND), TRUE,
                                    ifelse(signif & (fc <= FC_LOWER_BOUND), TRUE, FALSE)))
  
  diffe_final_res <- diffe_final_res_temp %>% 
    mutate(label_ = analyte) %>%
    dplyr::select(signif_and_fold, signif, logFC, fc, analyte, everything()); diffe_final_res
  
  signif_df <- diffe_final_res %>%
    filter(signif_and_fold)
  nrow(signif_df)
  
  # save this data to plot later
  final_lst <- list(diffe_final_res, signif_df) %>%
    setNames(c("diffe_final_res", "signif_df"))
  
  return(final_lst)
}

#' @note plot differential expression as volcano plot
#' @param dname name of the perturbation in question, or class of perturbations
#' @param diffe_final_res from the output of run_diffe()
#' @param signif_df from the output of run_diffe(), significant analytes
#' @param which_dat which dataset is this coming from? P100 or GCP
#' @return list of ggplots, combined as a facetted plot and singles
plot_diffe_results <- function(args){
  # DEBUG: diffe_final_res <- args$diffe_lst[[1]][[1]]; signif_df <- args$diffe_lst[[1]][[2]]
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-diffe.RData")
  # load("debug/debug_dat/debug-diffe.RData")
  # stop()
  
  dname <- force_natural(args$dirname_)
  diffe_final_res <- args$diffe_lst[[1]]
  signif_df <- args$diffe_lst[[2]] 
  which_dat <- force_natural(args$which_dat)
  
  dname_splt_temp <- str_split(dname, "_excl_", simplify = T)
  dname_splt <- dname_splt_temp[dname_splt_temp != ""]; dname_splt
  if (!(length(dname_splt) == 1)) {
    dname_title <- str_c(dname_splt, collapse = " excluding ")
  } else {
    dname_title <- dname
  }
  
  #' *DECIDE WHETHER TO USE MEAN OR MEDIAN* - here we're using median if comments are active
  diffe_final_res <- diffe_final_res # %>% mutate(fc = mean_fc) # for mean
  signif_df_final <- signif_df # %>% mutate(fc = mean_fc) # for mean
  
  sig_up_df <- signif_df_final %>% 
    filter(neg_log10_p_val_bh >= 1 & 
             (fc >= FC_UPPER_BOUND))
  sig_down_df <- signif_df_final %>% 
    filter(neg_log10_p_val_bh >= 1 & 
             (fc <= FC_LOWER_BOUND))
  
  to_plot_signif_df_final <- bind_rows(sig_up_df,sig_down_df) %>%
    mutate(label_ = str_replace(string = label_, pattern = "_[0-9]", replacement = "*"))
  
  # pretty breaks!! pretty_breaks
  # https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
  plot_diffe_all <- function(diffe_final_res, to_plot_signif_df_final, which_dat, dname_title) {
    
    min_x_lim <- min(to_plot_signif_df_final$fc, na.rm = T) - 0.1
    max_x_lim <- max(to_plot_signif_df_final$fc, na.rm = T) + 0.1
    
    filt_signif <- to_plot_signif_df_final %>% 
      filter(neg_log10_p_val_bh < Inf & neg_log10_p_val_bh > -Inf) %>% 
      .$neg_log10_p_val_bh
    max_y_lim <- max(filt_signif, na.rm = T)
    
    rel_size_label <- rel(4.5)
    rel_segment_size <- rel(1)
    rel_size_point <- rel(3)
    
    diffe_g <- ggplot(diffe_final_res) +
      geom_point(aes(x = fc, y = neg_log10_p_val_bh), # size = neg_log10_p_val_bh
                 size = rel_size_point, shape = 21, color = "black", # filled circle with outline
      ) +
      geom_point(data = to_plot_signif_df_final,
                 mapping = aes(x = fc, y = neg_log10_p_val_bh, fill = fc), 
                 size = rel_size_point, shape = 21, color = "black") +
      geom_vline(
        xintercept = as.numeric(c(FC_LOWER_BOUND, FC_UPPER_BOUND)),
        col = "purple"
      ) + # , alpha = 0.5
      geom_hline(
        yintercept = as.numeric(-log10(0.1)),
        col = "orange", linetype = 6,
      ) +
      geom_label_repel(to_plot_signif_df_final %>% 
                         filter(neg_log10_p_val_bh >= 1 & 
                                  (fc >= FC_UPPER_BOUND)),
                       mapping = aes(x = fc, y = neg_log10_p_val_bh, 
                                     label = label_),
                       seed = 42,
                       alpha = 0.7,
                       direction = "both", # direction to adjust labels, x, y, both
                       verbose = FALSE, # for debugging
                       size = rel_size_label,
                       max.iter = 1e9,
                       segment.size = rel_segment_size,
                       # max.time = 10,
                       # segment.shape = -1,
                       segment.curvature = -0.6,
                       segment.square = TRUE,
                       segment.color = 'red',
                       box.padding = 0.2,
                       min.segment.length = 0,
                       point.padding = 0.8,
                       # force_pull = 0.5,
                       force = 100,
                       nudge_x = 0.01,
                       nudge_y = 0.01,
                       xlim = c(1, max_x_lim),
                       ylim = c(0, max_y_lim+1),
                       arrow = arrow(length = unit(0.0075, "npc"))) +
      
      geom_label_repel(to_plot_signif_df_final %>%
                         filter(neg_log10_p_val_bh >= 1 & 
                                  (fc <= FC_LOWER_BOUND)),
                       mapping = aes(x = fc, y = neg_log10_p_val_bh,
                                     label = label_),
                       seed = 42,
                       alpha = 0.7,
                       direction = "both", # direction to adjust labels, x, y, both
                       verbose = FALSE, # for debugging
                       size = rel_size_label,
                       max.iter = 1e9,
                       segment.size = rel_segment_size,
                       segment.curvature = -0.6,
                       segment.square = TRUE,
                       segment.color = 'blue',
                       box.padding = 0.2,
                       min.segment.length = 0,
                       point.padding = 0.8,
                       # force_pull = 0.5,
                       force = 100,
                       nudge_x = -0.01,
                       nudge_y = 0.01,
                       xlim = c(min_x_lim-0.25, 1),
                       ylim = c(0, max_y_lim+1),
                       arrow = arrow(length = unit(0.0075, "npc"))) +
      
      scale_x_continuous(limits = c(min_x_lim-0.25, max_x_lim),
                         breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      # increase the space on the y axis so that the labels don't get clipped
      expand_limits(y = c(0, max_y_lim+1)) +
      scale_fill_gradient2(midpoint = 1.0, low = "blue", mid = "white", high = "red",
                           na.value = NA, space = "Lab",
                           name = "Fold Change") +
      theme_bw(base_size = 7) +
      theme(aspect.ratio = 1, # this is hard to pick
            legend.title = element_text(size = rel(1.75)),
            legend.text = element_text(size = rel(1.5)), 
            strip.text = element_text(size = rel(2)),
            plot.title = element_text(hjust = 0.5, size = rel(3), face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = rel(2.5), face = "bold"),
            axis.title = element_text(size = rel(2.75)),
            axis.text.y = element_text(size = rel(2.5), angle = 45),
            axis.text.x = element_text(size = rel(2.5), angle = 45,  vjust = 0.5)
      ) +
      labs(
        caption = "BH.q.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1 \nDifference between groups calculated using ks.test()"
      ) +
      ylab("-Log10(BH.q.val)") +
      xlab("Ratio of Median Fold Change\nbetween Clusters") +
      ggtitle(label = qq("@{toupper(which_dat)}, @{dname_title}")) 
    return(diffe_g)
  }
  
  diffe_g <- plot_diffe_all(diffe_final_res, to_plot_signif_df_final, which_dat, dname_title) + 
    facet_grid(~base_clust_comp_name, scales = "free"); diffe_g
  
  
  # for plotting just vascular cells later on
  # discovered a weird silent bug that was only apparent on windows...
  new_main <- diffe_final_res %>%
    mutate(keep_ = map_lgl(base_clust_comp_name, .f = function(x) {
      r <- ifelse(any(!is.na(str_locate(string = x, pattern = "HUVEC|HAoSMC"))), T, F)
      return(r)
    })) %>%
    filter(keep_); new_main
  
  new_sig <- to_plot_signif_df_final %>%
    mutate(keep_ = map_lgl(base_clust_comp_name, .f = function(x) {
      r <- ifelse(any(!is.na(str_locate(string = x, pattern = "HUVEC|HAoSMC"))), T, F)
      return(r)
    })) %>% 
    filter(keep_); new_sig

  groups <- unique(new_main$base_clust_comp_name)
  
  message("Plotting volcano, aligning labels to avoid overlap...")
  diffe_g_vasc <- future_map(groups, .f = function(g) {
    ndf <- filter(new_main, base_clust_comp_name == g)
    nns <- filter(new_sig, base_clust_comp_name == g)
    res <- plot_diffe_all(diffe_final_res = ndf, to_plot_signif_df_final = nns, 
                          which_dat = which_dat, dname_title = dname_title) 
    return(res)
  }); diffe_g_vasc
  
  group_name <- map_chr(groups, .f = function(g) {
    r <- str_extract_all(g, pattern = "HAoSMC|HUVEC", simplify = TRUE)
    return(str_c(r, collapse = ","))
  })
  # group_name <- ifelse(group_name == "HAoSMC,HUVEC", "All", group_name)
  names(diffe_g_vasc) <- str_c(dname, "-", group_name)
  
  init_lst <- list(diffe_g) %>%
    setNames("diffe_full_g")
  next_lst <- list(diffe_g_vasc) %>% setNames("diffe_singles")
  final_lst_gg <- c(init_lst, next_lst)
  return(final_lst_gg)
}

#' @note run connectivity and clustering analysis, with progress updates
#' @param lst a list of tbl's split by the grouping structure for the given analysis
#' @return result containing the input data, corr_lst (correlation results),
#' conn_lst (connectivity results), and clust_lst (clustering results)
run_analysis <- function(lst, tie_method = "average",
                         use_bootstrap = FALSE, use_parallel = FALSE) {
  # lst <- full_splt_lst; tie_method = "average"; use_bootstrap = FALSE; use_parallel = TRUE
  p <- progressr::progressor(steps = 3)
  p(message = "Computing correlation")
  corr_lst <- run_corr_lst(lst, tie_method = tie_method)
  
  p(message = "Computing connectivity")
  conn_lst <- run_conn_lst(corr_lst, use_bootstrap = use_bootstrap)
  
  p(message = "Clustering")
  clust_lst <- run_clust_lst(conn_lst, use_parallel = use_parallel)
  
  diffe_lst <- run_diffe_lst(lst, clust_lst)
  message("Done")
  
  res <- list(
    "input_data" = lst,
    "output_results" = list(
      "corr_lst" = corr_lst,
      "conn_lst" = conn_lst,
      "clust_lst" = clust_lst,
      "diffe_lst" = diffe_lst
    )
  )
  return(res)
}



update_dataset_with_clustering <- function(data, clust_lst) {
  cell_ids <- names(clust_lst$cluster_assignments)
  clust_id_df <- tibble(cell_id = cell_ids, cut_trees = clust_lst$cluster_assignments)
  return(data %>% left_join(clust_id_df, by = "cell_id"))
}





#' @note moved this out of main function
#' @note needs to be reworked...
get_debug_stats <- function(output_directory, dataset_type){
  message("Plotting pert-cell distribution data...")
  debug_plot_directory <- file.path(
    output_directory,
    dataset_type, "plots", "summaries"
  )
  dir.create(debug_plot_directory, recursive = T, showWarnings = F)
  
  prop_pert_df <- tibble(
    filter_vars = dir_name_df$new_filter_var,
    data = full_splt_lst
  ) %>%
    mutate(
      dirname_ = .[[1]],
      !!grouping_var := str_split(
        string = dirname_, pattern = "_",
        simplify = TRUE
      )[, 1]
    ) %>%
    mutate(prop_df = map(data, function(d) {
      res <- d %>%
        group_by(cell_id, pert_iname) %>%
        dplyr::summarize(n = n(), .groups = "keep") %>%
        group_by(cell_id) %>%
        mutate(prop = n / sum(n))
      return(res)
    })) %>%
    mutate(gplot = map(prop_df, function(d) {
      gplot <- ggbarplot(
        data = d, x = "cell_id", y = "prop",
        fill = "pert_iname",
        palette = "igv", ggtheme = theme_bw()
      ) +
        ggtitle("Proportion of therapies used in each cell type")
      return(gplot)
    })) %>%
    mutate(
      gplot_base_dir = debug_plot_directory,
      gplot_path = file.path(
        gplot_base_dir,
        str_c(dirname_, "--prop-of-drugs-per-cell.pdf")
      )
    )
  
  walk2(as.list(prop_pert_df$gplot),
        as.list(prop_pert_df$gplot_path),
        .f = function(x, y) {
          ggsave(x, filename = y, device = "pdf", width = 8, height = 10)
        }
  )
}


#' @note still working on this...
cor_strct <- function(args) {
  # DEBUG: args <- analysis_dat[1,]; dat_ <- args$data %>% pluck(1)
  dat_ <- args$data
  
  analytes_ <- unique(dat_$pr_gene_symbol)
  header_cols <- c("master_id", "replicate_id", "cell_id") # "master_id", "replicate_id", "pert_iname"
  master_short_dat <- dat_ %>% 
    dplyr::select(replicate_id, master_id) %>% 
    distinct()
  
  df_temp <- dat_ %>%
    group_by(replicate_id) %>%
    pivot_wider(replicate_id, pr_gene_symbol, values_fn = median) %>%
    left_join(master_short_dat) %>%
    dplyr::select(master_id, replicate_id, pert_iname, cell_id, everything()) %>%
    group_by(pert_iname) %>% # cell_id # master_id
    nest(data = all_of(c(header_cols, analytes_))) %>%
    ungroup()
  
  all_dat_df <- right_join(master_short_dat,
                           bind_rows(df_temp$data)) %>%
    nest(all_dat = colnames(.)) %>%
    suppressMessages()
  
  
  
  df <- df_temp %>%
    bind_cols(all_dat_df) %>%
    slice(1:10) %>%
    mutate(res = pmap(.l = list(pert_iname, data, all_dat),  # cell_id
                      .f = function(id, d, all_d) {
                        
                        # id = df$pert_iname[1]; d = df$data[1] %>% pluck(1); all_d = df$all_dat[1] %>% pluck(1)
                        message(id)
                        # title_ <-                   
                        # removes pert from "all"
                        bkgrd <- anti_join(all_d, d) %>%
                          ungroup() %>%
                          suppressMessages()
                        test <- d 
                        
                        mat_test <- as.matrix(test %>% dplyr::select(all_of(analytes_) ))
                        rownames(mat_test) <- d$replicate_id
                        cor_mat_test <- round(cor(t(mat_test), use = "pairwise.complete.obs"), 3)
                        
                        df_corr_test <- cor_mat_test %>%
                          as.data.frame() %>%
                          rownames_to_column(var = "cor1") %>%
                          as_tibble() %>%
                          mutate(grp = "test", .before = 1) %>%
                          pivot_longer(cols = all_of(d$replicate_id), names_to = "cor2")
                        
                        mat_bkgrd <- as.matrix(bkgrd %>% dplyr::select(all_of(analytes_)))
                        rownames(mat_bkgrd) <- bkgrd$replicate_id
                        cor_mat_bkgrd <- round(cor(t(mat_bkgrd), use = "pairwise.complete.obs"), 3)
                        
                        df_corr_bkgrd <- cor_mat_bkgrd %>%
                          as.data.frame() %>%
                          rownames_to_column(var = "cor1") %>%
                          as_tibble() %>%
                          mutate(grp = "bkgrd", .before = 1) %>%
                          pivot_longer(cols = all_of(bkgrd$replicate_id), names_to = "cor2")
                        
                        cmbd_df_temp <- bind_rows(df_corr_test, df_corr_bkgrd) %>%
                          dplyr::select(-cor1, -cor2)
                        
                        n_cmbd_df <- cmbd_df_temp %>%
                          group_by(grp) %>%
                          tally() %>%
                          mutate(n = format(n, big.mark=",")) %>%
                          mutate(n = str_trim(n))
                        
                        cmbd_df <- cmbd_df_temp %>%
                          left_join(n_cmbd_df, by = "grp") %>%
                          mutate(grp = ifelse(grp == "bkgrd", "Non-replicate", "Replicate")) %>%
                          mutate(grp = str_c(grp, " (n=",n,")"))
                        
                        # g_dens <- ggdensity(cmbd_df, x = "value", fill = "grp", 
                        #           xlab = "Spearman Correlation") +
                        #   ggtitle(id)
                        # 
                        
                        # median_cor_mat <- median(cor_mat, na.rm = TRUE)
                        # ggcorrplot(cor_mat, hc.order = TRUE, outline.col = "white")
                        # res <- list(cmbd_df, g_dens) %>%
                        #   setNames(c("res", "g_dens"))
                        
                        return(cmbd_df)
                      }))
  df
  ggdensity(df$res[[3]], x = "value", fill = "grp", add = "mean")
  
  return(df)
}
