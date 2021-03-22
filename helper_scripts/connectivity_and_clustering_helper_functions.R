
#' @note compute correlation off of data
#' @param x is a numeric matrix (transpose before import if you want rows -> columns and vice versa)
compute_correlation <- function(x, remove = TRUE) {
  c <- cor(x, method = "spearman", use = "pairwise.complete.obs")
  if (remove) {
    cx <- remove_x_perc_NA(dat = c, thresh_row = 0.5, thresh_col = 0.5, d = "both")$mat
    return(cx)
  }
  return(c)
}


#' @note get unique names from matrix
#' @param M matrix with row and column names
get_grp_names_from_matrix <- function(M, unique_cs, sep_ = "--") {
  grp_names <- str_split(unique_cs, pattern = sep_, simplify = T)[, 1]
  return(grp_names)
}

#' @NOTE: this is SO SO slow
#' @note compute for cell A vs cell B and cell B vs all others (not A)
#' @param cell_corr matrix of correlations among cells
#' @param use_bootstrap logical to calculate boostrapped p-value for the ks.test
compute_connectivity <- function(M, use_bootstrap = FALSE) {
  # M <- corr_mat_reps %>% slice(1) %>% .$cell_corr_r_mat ; M <- M[[1]]

  # here, we don't care above replicates because they get collapsed down!
  # so no need for the indexing
  unique_cs <- unique(rownames(M))
  grp_names <- get_grp_names_from_matrix(M, unique_cs, sep_ = "--")

  unique_cs_idx <- 1:length(rownames(M))
  cs_tib <- tibble(unique_cs, unique_cs_idx, grp_names)

  grp_idx <- bind_cols(tibble(
    grp_names = unique(grp_names),
    grp_idx = 1:length(unique(grp_names))
  ))

  cs_tib <- left_join(cs_tib, grp_idx, by = "grp_names")
  looper <- expand.grid(grp_namesA = unique(cs_tib$grp_names), grp_namesB = unique(cs_tib$grp_names))

  if (use_bootstrap) message("Using ks.boot")
  conn_df <- looper %>%
    mutate(map2_df(grp_namesA, grp_namesB, function(x, y) {
      # x <- looper$grp_namesA[5]; y <- looper$grp_namesB[5]
      # message(x,y)
      grpA_idx1 <- cs_tib %>%
        filter(grp_names == x) %>%
        .$unique_cs_idx
      grpB_idx2 <- cs_tib %>%
        filter(grp_names == y) %>%
        .$unique_cs_idx

      test <- M[grpA_idx1, grpB_idx2]

      grpA_idx2 <- cs_tib %>%
        filter(grp_names == y) %>%
        .$unique_cs_idx
      grpB_idx2 <- cs_tib %>%
        filter(!(grp_names %in% y)) %>%
        .$unique_cs_idx

      bkg <- M[grpA_idx2, grpB_idx2]

      if (use_bootstrap) {
        #' estimate ks stat on a population by sampling data with replacement
        #' provides correct coverage even when the data are not continuous
        #' (similar values in the two groups) ties are allowed in this test,
        #' unlike the original
        #'
        #' https://stats.stackexchange.com/questions/222294/understanding-kolmogorov-smirnov-test-in-r
        #'  The test statistic is the maximum distance between the ECDF's of the
        #'  two samples. The p-value is the probability of seeing a test
        #'  statistic as high or higher than the one observed if the two samples
        #'  were drawn from the same distribution. (It is not the "probability
        #'  that var1 = var2". And furthermore, 1-p_value
        #' is NOT the that probability either.)
        res <- ks.boot(test, bkg, nboots = 1000, alternative = "two.sided")
        connectivity <- res$ks$statistic
        p.val <- res$ks.boot.pvalue
      } else {
        res <- ks.test(test, bkg, alternative = "two.sided")
        connectivity <- res$statistic
        p.val <- res$p.val
      }

      if (median(test, na.rm = T) < median(bkg, na.rm = T)) connectivity <- -1 * connectivity
      cat(".")
      return(tibble(conn = connectivity, p_val = p.val))
    })) %>%
    as_tibble()
  cat("\n")

  return(conn_df)
}


#' @note collapse connectivity between pairs by taking the median
#' @NOTE you need to get a unique ID for the 3+ replicates!!
#' @param cell_conn [matrix]!! of connectivities among cells
#' @return symmetric matrix of connectivities
collapse_connectivity_by_median <- function(M) {

  # conn_mat_reps <- conn_reps %>%  mutate(conn_mat = map(conn, make_numeric_mat)) %>% mutate(rc_conn_mat = map(conn_mat, replace_rc_names, from = "well", to = "cell_id", rc = c("r", "c")))

  if (!is.matrix(M)) {
    M <- make_numeric_mat(M)
  }

  comps <- length(colnames(M))
  # here conn is really med_conn
  df <- tibble() # p_val/conn for compatability issues
  for (i in 1:comps) {
    for (j in 1:comps) {
      n1 <- colnames(M)[i]
      n2 <- colnames(M)[j]
      v1 <- M[i, j]
      v2 <- M[j, i]
      # v <- sort(c(v1,v2),decreasing = F) # important for median calculation
      median_val <- median(c(v1, v2))
      median_val
      # message(n1, ", ", i, ", ",n2, ", ",j, ", ", mat) # debug:
      new_df <- tibble(group_a = n1, group_b = n2, conn = median_val)
      df <- bind_rows(df, new_df)
    }
  }
  symm_med_M <- make_numeric_mat(df)
  return(symm_med_M)
}


#' @note PLOT AND *SAVE* pvclust results
#' @param res the pvclust result
#' @param dname the drug
#' @param co_clust true or false boolean indicating whether vars of interest co-clustered
#' @param output_dir output directory
#' @param thresh percent of max height to cut tree (visualize clusters)
plot_pvclust <- function(res, dname, co_clust, base_output_dir, thresh = 0.6) {
  if (is.null(res)) {
    return(NULL)
  }
  output_dir <- file.path(base_output_dir, "connectivity_pvclust", qq("h_@{thresh}"))
  add <- ifelse(co_clust, "CO_CLUST/", "OTHER/")
  dir.create(file.path(output_dir, "CO_CLUST"), recursive = T, showWarnings = F)
  dir.create(file.path(output_dir, "OTHER"), recursive = T, showWarnings = F)
  # Clusters with AU larger than 95% are highlighted by rectangles!! Strongly supported by the data
  par(mfrow = c(1, 1))
  pdf(file.path(output_dir, qq("@{add}@{dname}__h_@{thresh}.pdf")), width = 10, height = 12)
  plot(res, main = qq("@{dname} - cut at @{thresh} max height"))
  pvrect(res, alpha = 0.95)
  abline(h = thresh * max(res$hclust$height), col = "dark blue")
  dev.off()
}

#' @note PLOT AND *DISPLAY* pvclust results
#' @param res the pvclust result
#' @param dname the drug
#' @param thresh percent of max height to cut tree (visualize clusters)
display_pvclust <- function(res, dname, dataset, co_clust, thresh = 0.6) {
  if (is.null(res)) {
    return(NULL)
  }
  sub_col <- ifelse(co_clust, "red", "black")
  plot(res, main = qq("@{dname}"), sub = qq("@{dataset}\nCo-clust = @{co_clust}"), cex.sub = 2, col.sub = sub_col)
  pvrect(res, alpha = 0.95)
  abline(h = thresh * max(res$hclust$height), col = "dark blue")
}


#' @note PLOT AND DISPLAY pvclust results
#' @param res the pvclust result
#' @param ct cutree result
#' @param dname the drug
#' @param thresh percent of max height to cut tree (visualize clusters)
display_pvclust_poster <- function(res, ct, dname, dataset) {
  if (is.null(res)) {
    return(NULL)
  }
  k <- length(unique(unname(ct)))
  # print(qq("@{dname}, @{k}"))
  plot(res, main = qq("@{dname}"), sub = qq("@{dataset}"), cex.sub = 2)
  # pvrect(res, alpha=0.95, border = 16)
  rect.hclust(res$hclust, k = k, border = 2:(k + 1))
}


#' @note compute clusters via Lev's method
#' @param x numeric matrix for clustering
#' @param dname names of drugs (groups) -- only used for naming the result plot
#' @param base_output_dir name of base output dir for plots
compute_boot_pvclust_legacy <- function(x, dname, n_boot = 1000) {
  message(dname)
  if (!(nrow(x) > 2)) {
    message("Not enough data to perform clustering analysis. Degenerate clusters.")
    return(NULL)
  }
  res <- pvclust(x, method.dist = "cor", method.hclust = "average", nboot = n_boot, iseed = 2334, parallel = T)
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
    iseed = 2334, parallel = parallel_flag
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


#' @note unsupervised analysis
#' split cell x analyte x [group] -- group could be MOA, drug name, or All the drugs at once!
#' HCLUST by similarity -- call morpheus
#' Note that, conclusions about the proximity of two observations can be drawn only based on the height where
#' branches containing those two observations first are fused.
#' We cannot use the proximity of two observations along the horizontal axis as a criteria of their similarity.
#'
#' https://uc-r.github.io/hc_clustering
#' https://rlbarter.github.io/superheat/missing-data.html
#' https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html#methods-for-correlation-analyses

#' @note run connectivity analysis
#' @param obj from [load_data] OR an AVG connectivity object
#' @param base_output_dir for output files
#' @param avg_toggle whether or not combining datasets
#' @param dataset for naming, e.g. differentiation between gcp and p100
#' @param group_var grouping variable passed from batch_exec
#' @param suffix ending of rds file ("", "_moa", "_all)
#' @param rerun default FALSE, whether or not to rerun pvclust
run_connectivity_analysis <- function(obj, avg_toggle = FALSE, base_output_dir_conn_final,
                                      dataset, group_var = "pert_iname", suffix = "", rerun = FALSE,
                                      dendro_cut_thresh = 0.6, plt_morph_tog = TRUE, master_key = NULL) {

  # ensure string matching with lowercase
  group_var <- tolower(group_var)
  stopifnot(group_var %in% c("pert_iname", "drug_class", "all"))

  # obj = p100_obj_pert; base_output_dir; dataset
  clustering_output_fn <- file.path(base_output_dir_conn_final, qq("@{dataset}_conn_clust@{suffix}.rds"))
  if (!file.exists(clustering_output_fn)) {
    message("Clustering .rds doesn't exist. Running clustering...")
    rerun <- TRUE
  }

  if (rerun) {
    if (avg_toggle) {
      message("AVG connectivity already computed. Plotting..")
      # here, obj is made from joining two or more previous results

      message("Running AVG clustering...")
      if (group_var == "all") group_var <- "org"

      conn_clust <- obj %>%
        mutate(pvclust = map2(med_conn, !!sym(group_var), .f = compute_boot_pvclust))

      # plot connectivity
      walk2(
        .x = obj$med_conn, .y = obj[group_var] %>% pull(),
        function(x, y) {
          plot_morpheus(
            dat = x, drug_name = y, dataset = dataset, base_output_dir = base_output_dir_conn_final,
            subdir = "conn_reps_med", minr = -0.75, maxr = 0.75, sym = TRUE
          )
        }
      )

      conn_clust_trees <- conn_clust %>%
        mutate(which_dat = dataset) %>%
        mutate(cut_trees = map(pvclust, .f = generate_clusters, thresh = dendro_cut_thresh)) %>%
        mutate(co_clust = map(cut_trees, co_cluster))
      # conn_clust_trees$which_dat[[1]]
      # conn_clust_trees$cut_trees[[1]]
      # conn_clust_trees$co_clust[[1]]

      # plot dendrogram
      message("Drawing clusters and plotting basic dendrograms...")
      pwalk(list(conn_clust_trees$pvclust, conn_clust_trees[group_var] %>% pull(), conn_clust_trees$co_clust),
        .f = plot_pvclust,
        base_output_dir = base_output_dir_conn_final, thresh = dendro_cut_thresh
      )

      message("Writing .rds...")
      write_rds(conn_clust_trees, clustering_output_fn, compress = "gz")
    } else {
      message("Gathering replicates...")
      mat_reps_lst <- get_replicates(dat_obj = obj, grouping_var = group_var)
      mat_reps <- mat_reps_lst$post

      message("Computing correlation...")
      corr_mat_reps_temp <- mat_reps %>%
        mutate(
          corr_r = map(t_dataframe, compute_correlation),
          corr_nr = map(t_dataframe, compute_correlation, remove = F)
        )

      # stop()
      # define cell_id coding by grouping variable up front
      if (group_var != "all") {
        key_idx <- master_key %>%
          distinct(well, cell_id, pert_iname, drug_class) %>%
          group_by(!!sym(group_var)) %>%
          mutate(cell_id = make.unique(cell_id))
      } else {
        key_idx <- master_key %>%
          distinct(well, cell_id, pert_iname, drug_class) %>%
          mutate(cell_id = make.unique(cell_id))
      }

      corr_mat_reps <- corr_mat_reps_temp %>%
        mutate(match_corr_lst = map(corr_r,
          .f = replace_rc_names,
          match_df = key_idx
        )) %>%
        mutate(cell_corr_r_mat = map(match_corr_lst, function(l) l$new_mat))
      # M <- corr_mat_reps$cell_corr_r[[1]]

      message("Computing connectivity among cell lines...")
      conn_reps <- corr_mat_reps %>%
        mutate(conn = map(cell_corr_r_mat, .f = compute_connectivity))

      conn_mat_reps <- conn_reps %>%
        mutate(conn_mat = map(conn, make_numeric_mat)) %>%
        mutate(med_conn_mat = map(conn_mat, collapse_connectivity_by_median))

      # hack to make the 'all' condition work out
      if (group_var == "all") group_var <- "org"

      if (plt_morph_tog) {
        message("\nPlotting correlation and connectivity matrices...")

        walk2(
          .x = corr_mat_reps$cell_corr_r_mat, .y = corr_mat_reps[group_var] %>% pull(),
          function(x, y) {
            plot_morpheus(
              dat = x, drug_name = y, dataset = dataset, base_output_dir = base_output_dir_conn_final,
              subdir = "corr_reps", minr = -0.75, maxr = 0.75, sym = TRUE
            )
          }
        )

        walk2(
          .x = conn_mat_reps$conn_mat, .y = conn_mat_reps[group_var] %>% pull(),
          function(x, y) {
            plot_morpheus(
              dat = x, drug_name = y, dataset = dataset, base_output_dir = base_output_dir_conn_final,
              subdir = "conn_reps", minr = -0.75, maxr = 0.75, sym = TRUE
            )
          }
        )

        walk2(
          .x = conn_mat_reps$med_conn_mat, .y = conn_mat_reps[group_var] %>% pull(),
          function(x, y) {
            plot_morpheus(
              dat = x, drug_name = y, dataset = dataset, base_output_dir = base_output_dir_conn_final,
              subdir = "conn_reps_med", minr = -0.75, maxr = 0.75, sym = TRUE
            )
          }
        )
      }

      message("\n\nPerforming clustering...")
      conn_clust <- conn_mat_reps %>%
        mutate(pvclust = map2(med_conn_mat, !!sym(group_var), .f = compute_boot_pvclust))

      message(qq("\nCutting tree at @{dendro_cut_thresh}.."))
      conn_clust_trees <- conn_clust %>%
        mutate(which_dat = dataset) %>%
        mutate(cut_trees = map(pvclust, .f = generate_clusters, thresh = dendro_cut_thresh)) %>%
        mutate(co_clust = map(cut_trees, co_cluster))
      # conn_clust_trees$which_dat[[1]]
      # conn_clust_trees$cut_trees[[1]]
      # conn_clust_trees$co_clust[[1]]

      # plot dendrogram
      message("Drawing clusters and plotting basic dendrograms...")
      pwalk(list(conn_clust_trees$pvclust, conn_clust_trees[group_var] %>% pull(), conn_clust_trees$co_clust),
        .f = plot_pvclust,
        base_output_dir = base_output_dir_conn_final, thresh = dendro_cut_thresh
      )

      message("Writing .rds...")
      write_rds(conn_clust_trees, clustering_output_fn, compress = "gz")
    }
  } else {
    message("Reading clustering .rds...")
    conn_clust_trees <- read_rds(clustering_output_fn)
  }

  message("Done!")
  return(conn_clust_trees)
}




#### EXTRA



#' @note return grouping for poster
#' @param cut_tree_obj named vector of cluster assignments
#' @param name1 first var you want to compare
#' @param name2 second var you want to compare
poster_grouping <- function(cut_tree_obj, name1 = "HAoSMC", name2 = "HUVEC") {
  if (is.null(cut_tree_obj)) {
    return(NA)
  }
  name1_obj_idx <- which(name1 == names(cut_tree_obj))
  name2_obj_idx <- which(name2 == names(cut_tree_obj))
  c1 <- cut_tree_obj[name1_obj_idx]
  c2 <- cut_tree_obj[name2_obj_idx]

  if (c1 != c2) {
    return("3")
  }
  ni <- intersect(c1, c2)
  g <- which(cut_tree_obj == ni)
  if (length(g) > 2) {
    return("2")
  } else {
    return("1")
  }
}


dend_flowchart <- function(pvclust_obj1, pvclust_obj2) {
  # pv1 <- pvclust_obj1$pvclust[[4]]; pv2 <- pvclust_obj2$pvclust[[4]]
  pv1_dend <- as.dendrogram(pv1$hclust)
  pv2_dend <- as.dendrogram(pv2$hclust)
  dend_list <- dendlist(pv1_dend, pv2_dend)
  tanglegram(pv1_dend, pv2_dend,
    common_subtrees_color_lines = TRUE,
    common_subtrees_color_branches = TRUE,
    main = paste("entanglement =", round(entanglement(dend_list), 2))
  )
}

#' @Note cophonetic distance to compute similarities across hierarchical clustering outputs
#' The cophenetic distance between two observations that have been clustered is defined to be the intergroup
#' dissimilarity at which the two observations are first combined into a single cluster. Note that this distance
#' has many ties and restrictions.
#' It can be argued that a dendrogram is an appropriate summary of some data if the correlation between the original distances
#' and the cophenetic distances is high.
#' @param obj_list e.g. perts$p100
#' @return cophenetic distance matrix. The value can range between -1 to 1.
#' With near 0 values meaning that the two trees are not statistically similar.

# BROKEN NEED TO FIXE DIMENSION ISSUES -- ugh
dend_corr <- function(obj_list, members_thresh = 8) {
  # obj_lst <- map(perts, function(x) filter(x, pert_iname %in% p100_group1))$p100
  which.pvclust.keep <- which(!sapply(obj_list$pvclust, is.null))

  complete_labs <- obj_list$pert_iname[which.pvclust.keep]
  left_out_num_cells <- c()


  dends <- dendlist()
  mat_labels <- c()
  for (i in which.pvclust.keep) {
    t <- obj_list$pvclust[[i]]$hclust

    if (length(t$labels) != members_thresh) {
      left_out_num_cells <- c(left_out_num_cells, length(t$labels))
    }
    y <- as.dendrogram(t)

    dends <- dendlist(dends, y)
    lab <- obj_list$pert_iname[[i]]
    mat_labels <- c(mat_labels, lab)
  }
  res <- cor.dendlist(dends, method = "cophenetic")
  colnames(res) <- mat_labels
  rownames(res) <- mat_labels

  cat(qq("\n\n Correlation Structure across Hierarchical Clustering Results \n\n"))
  cor_sig <- cor.mtest(res, conf.level = .95)
  corrplot(res, "pie", "lower")
  cat("\n\n\\pagebreak\n")
  corrplot(res,
    p.mat = cor_sig$p, method = "color", type = "lower",
    sig.level = c(.001, .01, .05), pch.cex = .9,
    insig = "label_sig", pch.col = "red"
  )
  cat("\n\n\\pagebreak\n")

  return(res)
}

#' @note which drugs were left out of the clustering similarity analysis
#' @param mat_labs names of drugs that are being compared
#' @param complete_labs all the drugs
left_out <- function(mat_labs, complete_labs, num_cells) {
  tbl <- tibble(`Left Out` = setdiff(complete_labs, mat_labs), `Number of cells` = num_cells)
  print(kable(tbl, booktabs = TRUE, caption = qq("Summary of perturbation distribution for @{dataset}")) %>%
    kable_styling(latex_options = "hold_position"))
}


#' @note creates a cytoscape network from a connectivity matrix
#' @param mat a connectivity matrix (n x n)
#' @param top_x_perc top x percent of connections to retain in the network, default is 5%. [note: should be between 0 and 1]
#' @param dname group naming for file
#' @param base_output_dir name output directory for graph
create_cytoscape_network <- function(mat, top_x_perc = 0.20, dname, dataset, base_output_dir) {

  # mat <- conn_mat_reps$conn_mat[[23]] ; top_x_perc=0.1
  # use quantile to extract the top x percent of correlations
  message(qq("@{dname}"))
  edges <- data.frame(mat) %>%
    flattenMatrix() %>%
    select(source = i, target = j, temp_weight = cor) %>%
    mutate(temp_weight = map_dbl(temp_weight, function(x) x + exp(runif(1, log(1e-9), log(1))))) %>% # add small correction
    mutate(interaction = ifelse(temp_weight > 0, "connected", "disconnected")) %>%
    mutate(pos_weight = abs(temp_weight)) %>%
    filter(pos_weight > quantile(pos_weight, prob = 1 - top_x_perc)) %>%
    mutate(weight = exp(temp_weight)) %>%
    mutate(color = ifelse(temp_weight > 0, "green", "red"))

  nodes <- data.frame(id = dplyr::union(edges$source, edges$target))

  # gD <- igraph::simplify(igraph::graph.data.frame(edges, directed=FALSE))
  # gD <- igraph::set.edge.attribute(gD, "weight", index = igraph::E(gD), value = igraph::E(gD)$weight)

  createNetworkFromDataFrames(nodes, edges, title = qq("@{dname}"), collection = qq("@{dataset}"))
  # createNetworkFromIgraph(gD, title=qq("@{dname}"), collection=qq("@{dataset}"))

  style.name <- "myStyle"
  defaults <- list(
    NODE_SHAPE = "diamond",
    NODE_SIZE = 30,
    EDGE_TRANSPARENCY = 120,
    NODE_LABEL_POSITION = "W,E,c,0.00,0.00"
  )
  nodeLabels <- mapVisualProperty(visual.prop = "node label", table.column = "id", mapping.type = "p")
  arrowShapes <- mapVisualProperty("Edge Target Arrow Shape", "interaction", "d", c("connected", "disconnected"), c("Arrow", "T"))
  edgeWidth <- mapVisualProperty(visual.prop = "edge width", table.column = "weight", mapping.type = "p")
  # edgeColors <- mapVisualProperty(visual.prop = 'edge color', table.column = "color", table.column.values = c("green", "red"),
  #                     visual.prop.values = col2hex(c("green","red")), mapping.type = "d")

  createVisualStyle(style.name, defaults, list(nodeLabels, arrowShapes, edgeWidth)) # edgeColors,
  setVisualStyle(style.name)

  output_dir <- file.path(base_output_dir, "cytoscape")
  dir.create(output_dir, showWarnings = F, recursive = T)
  full.path <- file.path(output_dir, qq("@{dname}.pdf"))
  exportImage(full.path, "PDF") # .pdf
}


#' [ARCHIVED]
#' @note kmeans on connectivity matrix
do_kmeans <- function(conn_mat, dname, base_output_dir, subdir = "kmeans_cell_conn") {
  output_dir <- file.path(base_output_dir, "kmeans", subdir)
  dir.create(output_dir, recursive = T)

  sil <- fviz_nbclust(conn_mat, kmeans, method = "silhouette", k.max = nrow(conn_mat) - 1)
  numK <- which.max(sil$data$y)
  k3 <- kmeans(conn_mat, centers = numK, iter.max = 1000, nstart = 25)
  viz <- fviz_cluster(object = k3, data = conn_mat) + ggtitle(dname)
  ggexport(viz, filename = file.path(output_dir, str_c(dname, ".pdf")))
  return(k3)
}