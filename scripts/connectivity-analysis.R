
winos <- ifelse(grepl("windows", Sys.info()["sysname"], ignore.case = T), 1, 0)
if (winos == 1) {
  source("C:\\Users\\ncama\\OneDrive - Tufts\\phd\\ws\\scripts\\master-source.R")
} else {
  source("/Users/Nicholas/OneDrive - Tufts/phd/ws/proteomics/scripts/master-source.R")
}

cat("\nReading and merging data...")
#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
#' the values in this data are Z-scores!!!

collect_args <- function(str_char) {
  res <- as.vector(str_split(string = str_char, pattern = "--", simplify = T))
  return(res)
}
# TODO change this to [analysis_args.csv] when ready
analysis_fn <- file.path(DATASETS_DIRECTORY, "test_args.csv")
analysis_dat_temp <- read_csv(analysis_fn, comment = "#") %>%
  mutate_all(str_trim) %>%
  mutate(filter_vars = map(filter_vars, collect_args))

#' @note read in both GCP and P100 since it's cheap, then analyze accordingly
p100_fn <- file.path(
  DATASETS_DIRECTORY, "Inherited Data",
  "1st gen data", "P100", "P100-All-Cell-Lines.gct"
)
gcp_fn <- file.path(
  DATASETS_DIRECTORY,
  "Inherited Data", "1st gen data",
  "GCP", "GCP All Cell Lines.gct"
)
merged_obj_choice <- tibble(gct = map(c(p100_fn, gcp_fn), parse_gctx)) %>%
  bind_cols(tibble(dataset_type = c("P100", "GCP")))

analysis_dat <- left_join(
  analysis_dat_temp,
  merged_obj_choice,
  by = "dataset_type"
)

#' @note convenience function to 
#' help with annoying list obj with purrr
force_natural <- function(obj) {
  if (is.list(obj)) {
    return(obj[[1]])
  } else {
    return(obj)
  }
}

# OPTIMIZE all of these could probably be made faster
# OPTIMIZE with matrix mult
#' @note run correlation
run_corr <- function(x) {
  # x <- full_splt_lst[[1]]
  wide_x <- x %>%
    ungroup() %>%
    select(master_id, replicate_id, pr_gene_symbol, value) %>%
    pivot_wider(names_from = pr_gene_symbol, values_from = value)

  transposed_x <- wide_x %>%
    select(-master_id, -replicate_id) %>%
    t()

  res <- compute_correlation(transposed_x)

  rownames(res) <- wide_x$replicate_id
  colnames(res) <- wide_x$replicate_id

  res_final1 <- res %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble()
  colnames(res_final1) <- c("unique_id_a", wide_x$replicate_id)

  look_up_df <- x %>% distinct(master_id, replicate_id)
  res_final2 <- res_final1 %>%
    pivot_longer(
      cols = colnames(res_final1)[-1],
      names_to = "unique_id_b", values_to = "corr"
    ) %>%
    right_join(look_up_df %>%
      rename(group_a = master_id, unique_id_a = replicate_id),
    by = "unique_id_a"
    ) %>%
    right_join(look_up_df %>%
      rename(group_b = master_id, unique_id_b = replicate_id),
    by = "unique_id_b"
    )

  res_final2
  cat(".")
  return(list("matrix" = res, "tibble" = res_final2))
}

# connectivity is computed across replicates!!
#' @note connectiivty and clustering super function
#' @param corr_lst sample-level correlation values
run_conn <- function(corr_lst) {
  # NOTE:
  # suppressMessages(library(tidyverse))
  # suppressMessages(library(Matching))
  # suppressMessages(library(pvclust))
  # suppressWarnings(suppressMessages(source("./scripts/master-source.R")))
  # corr_lst <- sample_corr_lst$`1271738-62-5`

  corr_tbl <- corr_lst$tibble
  groups <- corr_tbl %>%
    distinct(group_a) %>%
    .$group_a

  # is this a multidyplr moment?
  # this is SO slow
  # i think i can apply some fancy matrix stuff to get this to go faster
  res <- map_df(groups, function(g1) {
    hold_g1 <- map_df(groups, function(g2) {
      # g1 <- groups[1]; g2 <- groups[2]

      groups_not_in_g1 <- groups[(groups != g1)]
      test <- corr_tbl %>%
        filter(group_a == g1 & group_b == g2) %>%
        pluck("corr")
      background <- corr_tbl %>%
        filter(group_a == g2, group_b %in% groups_not_in_g1) %>%
        pluck("corr")

      if (is.null(background)) {
        sub_res <- tibble(
          conn = NA, p_val = NA,
          group_a = g1, group_b = g2,
          test = list(test), background = list(background)
        )
        return(sub_res)
      }

      ks_res <- ks.boot(test, background,
        nboot = 1000, alternative = "two.sided"
      )

      stat <- ks_res$ks$statistic
      if (median(test, na.rm = T) < median(background, na.rm = T)) {
        stat <- -1 * stat
      }

      p_val <- ks_res$ks.boot.pvalue
      sub_res <- tibble(
        conn = stat, p_val = p_val,
        group_a = g1, group_b = g2,
        test = list(test), background = list(background)
      )
      return(sub_res)
    })

    return(hold_g1)
  })

  res_coded_for_median <- res %>%
    mutate(
      group_a = as.character(group_a),
      group_b = as.character(group_b)
    ) %>%
    mutate(grouping_var_code = map2_chr(group_a, group_b, function(a, b) {
      vec_a <- str_split(string = a, pattern = "--", simplify = T)
      a_code <- vec_a[1] # c(1, length(vec_a))
      vec_b <- str_split(string = b, pattern = "--", simplify = T)
      b_code <- vec_b[1] # c(1, length(vec_b))
      # grouping by this code should create pairs
      coded <- str_c(sort(c(a_code, b_code)), collapse = "_")
      return(coded)
    })) %>%
    group_by(grouping_var_code) %>%
    arrange(grouping_var_code)

  res_median <- res_coded_for_median %>%
    # should guarantee that no NA's in data
    summarize(median_conn = median(conn, na.rm = T), .groups = "drop") %>%
    right_join(res_coded_for_median, by = "grouping_var_code") %>%
    dplyr::select(-test, -background, -p_val) %>%
    ungroup() %>%
    dplyr::select(-grouping_var_code, -conn) %>%
    spread(group_b, median_conn) %>%
    # replace all i == j entries with 1
    mutate_at(vars(-group_cols()), ~ replace(., is.na(.), 1)) %>%
    as.data.frame() %>%
    dplyr::select(-group_a) %>%
    as.matrix()
  rownames(res_median) <- colnames(res_median)
  cat(".")

  return(list(
    "res" = res,
    "res_median_tibble" = res_coded_for_median,
    "res_median_matrix" = res_median
  ))
}

# clustering
run_clust <- function(conn_lst, pvclust_parallel = FALSE) {
  # extract median matrix from conn_lst
  # DEBUG:
  res_median <- conn_lst$res_median_matrix

  # converts cell-drug-class names into just cell names
  cell_names_from_long <- str_split(
    string = rownames(res_median),
    pattern = "--", simplify = T
  )[, 1]
  pretty_name_res_median <- res_median
  rownames(pretty_name_res_median) <- cell_names_from_long
  colnames(pretty_name_res_median) <- cell_names_from_long

  clust_obj <- compute_boot_pvclust(
    x = pretty_name_res_median,
    parallel_flag = pvclust_parallel,
    n_boot = 1000
  )
  clust_assignments <- generate_clusters(
    x = clust_obj,
    thresh = DENDRO_CUT_THRESH
  )
  co_clust_bool <- co_cluster(cut_tree_obj = clust_assignments)

  return(list(
    "clust_obj" = clust_obj,
    "clust_assignments" = clust_assignments,
    "co_clust_bool" = co_clust_bool
  ))
}


# analysis apply function
analysis_res <- apply(analysis_dat, 1, function(args) {
  # DEBUG: args = analysis_dat[1,]

  dataset_type <- args$dataset_type
  grouping_var <- args$grouping_var
  filter_vars <- force_natural(args$filter_vars)
  merged_obj <- force_natural(args$gct)

  # convenience fnc

  message("Creating top output directory: ")
  SPECIFIC_OUTPUT_DIR <- file.path(OUTPUT_DIRECTORY, dataset_type, grouping_var)
  dir.create(SPECIFIC_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  message(qq("@{SPECIFIC_OUTPUT_DIR}"))

  drugs_moa_df <- create_my_drugs_df(ref_dir = REFERENCES_DIRECTORY)

  # rank_obj <- rank_gct(merged_obj)

  
  my_obj <- melt_gct(merged_obj) %>%
    as_tibble() %>%
    rename(
      row_id = id.x,
      column_id = id.y,
      temp_pr_gene_symbol = pr_gene_symbol,
    ) %>%
    mutate(pert_iname = tolower(pert_iname)) %>%
    ## inner join to only do drugs that I pick!!
    left_join(drugs_moa_df, by = "pert_iname") %>%
    mutate(master_id = str_c(cell_id, pert_iname, pert_class, sep = "--")) %>%
    mutate(replicate_id = str_c(master_id, column_id, sep = "::")) %>%
    select(master_id, replicate_id, everything()) %>%
    # non-standard evaluation to find column for group 
    # and filter
    filter(!!sym(grouping_var) %in% filter_vars) %>%
    group_by(master_id, replicate_id) %>%
    # some gene names are the same
    # so you need to make them unique and reference
    # the phosphosite later
    mutate(pr_gene_symbol = make.unique(temp_pr_gene_symbol, sep = "_"))

  summary_read_in <- my_obj %>%
    ungroup() %>%
    distinct(!!sym(grouping_var))

  message("Data filtered down to contain perts/classes: ")
  print(summary_read_in)
  message("Cell IDs in this analysis: ")
  print(my_obj$cell_id %>% unique())

  full_splt_lst <- split(my_obj, f = my_obj[[sym(grouping_var)]])
  stopifnot(length(full_splt_lst) == length(filter_vars))

  message("Computing correlation...")
  # OPTIMIZE could this be faster with matrix mult??
  sample_corr_lst <- map(full_splt_lst, run_corr)
  cat(".Done.\n")


  connctivity_output_dir <- file.path(SPECIFIC_OUTPUT_DIR, "conn")
  dir.create(connctivity_output_dir, showWarnings = FALSE, recursive = TRUE)
  clust_output_dir <- file.path(SPECIFIC_OUTPUT_DIR, "clust")
  dir.create(connctivity_output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Computing connectivity, clustering...")
  sample_conn_clust_obj <- tibble(conn_lst = map(sample_corr_lst, run_conn)) %>%
    mutate(
      clust_lst = map(conn_lst, run_clust, pvclust_parallel = TRUE),
      name = names(clust_lst)
    ) %>%
    mutate(
      s_connctivity_output_dir = file.path(connctivity_output_dir, name),
      s_clust_output_dir = file.path(clust_output_dir, name)
    )
  message("Done.")

  # TODO: PLOT PVCLUST
  # TODO: HEATMAPS
  # TODO: CLUSTER HIERARCHY

  anon_write_to_file = function(obj, dir_) {
    dir.create(dir_, recursive = T, showWarnings = F)
    fn <- file.path(dir_, qq("@{basename(dir_)}_clust.rds"))
    write_rds(x = obj, file = fn, compress = "gz")
  }
  message("Writing to files...")
  # write to dirs
  walk2(
    sample_conn_clust_obj$conn_lst,
    sample_conn_clust_obj$s_connctivity_output_dir, 
    anon_write_to_file
  )

  walk2(
    sample_conn_clust_obj$clust_lst,
    sample_conn_clust_obj$s_clust_output_dir,
    anon_write_to_file
  )
  message("Done.")

  message("Caching full obj...")
  complete_obj <- list(
    "data" = my_obj,
    "corr_lst" = sample_corr_lst,
    "conn_clust_lst" = sample_conn_clust_obj
  )
  cache_output_dir <- file.path(SPECIFIC_OUTPUT_DIR, "cache")
  dir.create(cache_output_dir, showWarnings = F, recursive = T)
  cache_fn <- file.path(cache_output_dir, qq("cache-@{Sys.Date()}.rds"))
  write_rds(x = complete_obj, file = cache_fn)
  message("Done.\n")
  return(complete_obj)
})

# TODO: diff ex, heatmaps dendros



# library(Matrix)

# mat1 <- matrix(data = 1:8, nrow = 4, ncol = 4)
# mat2 <- matrix(data = 9:16, nrow = 4, ncol = 4)
# # diag returns vector, whereas constructor Matrix::Diagonal() returns matrix
# Diagonal(4, diag(mat1))
# mat1 %o% mat2



  # response <- readline(prompt = "Overwrite [y/n]: ")
  # if (response == "y") {
  #   message("\nComputing connectivity and performing clustering...")
  #   multi_obj_save <- run_conn_clust_par(my_obj, sample_corr_lst)
  # } else if (response == "n") {
  #   if (!file.exists(res_cache_obj_fn)) {
  #     message("Cache obj doesn't exist. Computing conn/clust...")
  #     message("\nComputing connectivity and performing clustering...")
  #     multi_obj_save <- run_conn_clust_par(my_obj, sample_corr_lst)
  #   } else {
  #     message("Loading from cache...")
  #     multi_obj_save <- read_rds(file.path(connctivity_output_dir, "data.rds"))
  #   }
  # } else {
  #   message("Incorrect input. Will attempt to load from cache.")
  #   multi_obj_save <- tryCatch(
  #     {
  #       read_rds(file.path(connctivity_output_dir, "data.rds"))
  #       message("Successfully loaded from cache.")
  #     },
  #     error = function(cond) {
  #       message("Loading failed. Running conn/clust analysis...")
  #       multi_obj_save <- run_conn_clust_par(my_obj, sample_corr_lst)
  #     }
  #   )
  # }

 #' @note run parallel connectivity and clustering analysis
  #' @param my_obj
  #' @param sample_corr_lst sample correlation list obj
  #' @return [multi_obj_save] combined output of analysis so far
  # run_conn_clust_par <- function(my_obj, sample_corr_lst) {
    # clock1 <- proc.time()
  #   n_cores <- registered()$SnowParam$.clusterargs$spec
  #   n_workers <- ifelse(length(sample_corr_lst) < n_cores,
  #     length(sample_corr_lst),
  #     n_cores
  #   )
    
  #   param <- SnowParam(
  #     # number of cores
  #     workers = n_workers,
  #     # socket cluster can be used with windows backend
  #     type = "SOCK",
  #   )

  #   bplog(param) <- TRUE
  #   bpthreshold(param) <- "INFO"
  #   log_dir <- file.path(OUTPUT_DIRECTORY, "log")
  #   dir.create(log_dir)
  #   bplogdir(param) <- log_dir

  #   # dmso <- compute_connectivity(M = sample_corr_lst$dmso$matrix)
  #   # dmso_median <- collapse_connectivity_by_median(dmso %>%
  #   #   rename(group_a = grp_namesA, group_b = grp_namesB))
  #   # dmso_clust_obj <- compute_boot_pvclust(dmso_median)
  #   # dmso_clust_assign <- generate_clusters(x = dmso_clust_obj,
  #   # thresh = DENDRO_CUT_THRESH)
  #   #  co_clust_bool <- co_cluster(cut_tree_obj = dmso_clust_assign)

  #   
  #  
  #     X <- list(1, 2, "3", 4, 5, 6)
  
  #   # parallel computation over each pert_iname or pert_class
  #   sample_conn_lst <- bplapply(
  #     X = sample_corr_lst,
  #     FUN = run_conn_clust,
  #     BPPARAM = param
  #   )

  #   clock2 <- proc.time() - clock1
  #   message("Done")
  #   message(qq("Completed processing @{length(sample_conn_lst)} items"))
  #   message(qq("Total time: in @{unname(clock2[3])} seconds"))


  #   return(sample_conn_lst)
  # }

  