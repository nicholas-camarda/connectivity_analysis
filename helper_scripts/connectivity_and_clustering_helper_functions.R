#' @note helper function for compute_boot_pvclust and morpheus plotter
#' @param long_df to be converted to matrix -- [assumes long df is a connectivity df!]
#' @return numeric matrix, symmetric
make_numeric_mat <- function(long_df) {
  if ("p_val" %in% colnames(long_df)) {
    long_df <- long_df %>% dplyr::select(-p_val)
  }
  x <- long_df %>%
    pivot_wider(names_from = group_b, values_from = conn) %>%
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
  wide_x <- x_ungrouped %>%
    dplyr::select(master_id, replicate_id, pr_gene_symbol, value) %>%
    pivot_wider(names_from = pr_gene_symbol, values_from = value)

  transposed_x <- wide_x %>%
    dplyr::select(-master_id, -replicate_id) %>%
    t() %>%
    # unnaming the large matrix results
    # in massive performance upgrade
    unname()

  # first convert to ranks, 
  # then compute pearson correlation with pairwise.complete.obs
  # try dense ties.method, might switch to average
  # rank <- matrixStats::rowRanks(transposed_x, ties.method = tie_method)
  # res <- cor(rank, method = "pearson", use = "pairwise.complete.obs")

  # https://stackoverflow.com/questions/10711395/spearman-correlation-and-ties
  # for lots of ties in data, use the kendall coefficient
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
  # message(a,b)

  # A vs B
  grp_col_idx_a <- my_tib %>% filter(grp_names == a) %>% .$cs_idx
  grp_col_idx_b <- my_tib %>% filter(grp_names == b) %>% .$cs_idx

  # B vs all (except itself)
  grp_col_idx_b <- my_tib %>% filter(grp_names == b) %>% .$cs_idx
  grp_col_idx_all <- my_tib %>% filter(grp_names != b) %>% .$cs_idx

  m <- matrix(NA,
    nrow = nrow(my_tib), ncol = nrow(my_tib)
  )
  m[grp_col_idx_a, grp_col_idx_b] <- 5
  m[grp_col_idx_b, grp_col_idx_all] <- 1
  m_outer <- outer(seq_len(nrow(m)), seq_len(ncol(m)),
    FUN = "paste", sep = ","
  )
  # temp <- matrix(1:20, nrow = 4, ncol = 5)
  # temp_outer <- outer(seq_len(nrow(temp)), seq_len(ncol(temp)),
  #   FUN = "paste", sep = ","
  # )

  # visualize which parts of the matrix are being used for connectivity

  # fill in (color 1) for A vs B
  # fill in (color 2) for B vs all
  # ht <- ComplexHeatmap::Heatmap(m, name = "mat vis",
  #   na_col = "white", border = TRUE,
  #   cell_fun = function(j, i, x, y, width, height, fill) {
  #     if (!is.na(m[i, j])) {
  #       grid.text(sprintf("%s", m_outer[i, j]), x, y, gp = gpar(fontsize = 10))
  #     }
  #   },
  #   cluster_rows = FALSE, cluster_columns = FALSE, use_raster = TRUE,
  # )
  # draw(ht)
  # pdf(qq("~/Downloads/heatmap_@{a}-@{b}.pdf"), width = 8, height = 8)
  # dev.off()

  
  test <- my_mat[grp_col_idx_a, grp_col_idx_b]
  bkg <- my_mat[grp_col_idx_b, grp_col_idx_all]

  if (use_bootstrap) {
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

  
  return(list("conn_tbl" = conn_df, "conn_median_mat" = median_mat))
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
  pv_obj <- compute_boot_pvclust(
    x = mat,
    parallel_flag = use_parallel, n_boot = 1000
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
  res <- pvclust(x,
    method.dist = "cor",
    method.hclust = "average",
    quiet = TRUE,
    nboot = n_boot,
    iseed = 2334,
    parallel = parallel_flag
  )
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

#' @note run connectivity and clustering analysis, with progress updates
#' @param lst a list of tbl's split by the grouping structure for the given analysis
#' @return result containing the input data, corr_lst (correlation results),
#' conn_lst (connectivity results), and clust_lst (clustering results)
run_analysis <- function(lst, tie_method = "average",
                         use_bootstrap = FALSE, use_parallel = FALSE) {
  p <- progressr::progressor(steps = 4)
  p(message = "Computing correlation")
  corr_lst <- run_corr_lst(lst, tie_method = tie_method)
  
  p(message = "Computing connectivity")

  # DEBUG REMOVED
  stop("Debug here")

  conn_lst <- run_conn_lst(corr_lst, use_bootstrap = use_bootstrap)
  p(message = "Clustering")
  clust_lst <- run_clust_lst(conn_lst, use_parallel = use_parallel)
  message("Done")

  res <- list(
    "input_data" = lst, "output_results" = list("corr_lst" = corr_lst,
    "conn_lst" = conn_lst, "clust_lst" = clust_lst)
  )
  return(res)
}

