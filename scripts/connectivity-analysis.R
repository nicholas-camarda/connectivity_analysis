# source the init.R file first
source(file.path("scripts/init.R"))

cat("\nReading and merging data...")
#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
#' the values in this data are Z-scores!!!

# TODO change this to [analysis_args.csv] when ready
analysis_fn <- file.path(data_directory, "test_args.csv")
analysis_dat_temp <- read_csv(analysis_fn, comment = "#") %>%
  mutate_all(str_trim) %>%
  mutate(filter_vars = map(filter_vars, collect_args))

#' @note read in both GCP and P100 since it's cheap, then analyze accordingly
p100_fn <- "P100-All-Cell-Lines.gct"
gcp_fn <- "GCP All Cell Lines.gct"
merged_obj_choice <- tibble(gct = map(
  .x = c(
    file.path(datasets_directory, p100_fn),
    file.path(datasets_directory, gcp_fn)
  ),
  .f = parse_gctx
)) %>%
  bind_cols(tibble(dataset_type = c("P100", "GCP")))

analysis_dat <- left_join(
  analysis_dat_temp,
  merged_obj_choice,
  by = "dataset_type"
)

# analysis apply function
analysis_res <- apply(analysis_dat, 1, function(args) {
  # DEBUG: args = analysis_dat[1,]

  dataset_type <- args$dataset_type
  grouping_var <- args$grouping_var
  filter_vars <- force_natural(args$filter_vars)
  merged_obj <- force_natural(args$gct)

  drugs_moa_df <- create_my_drugs_df(ref_dir = references_directory)
  my_obj <- melt_gct(merged_obj) %>%
    as_tibble() %>%
    dplyr::rename(
      row_id = id.x,
      column_id = id.y,
      temp_pr_gene_symbol = pr_gene_symbol,
    ) %>%
    mutate(pert_iname = tolower(pert_iname)) %>%
    ## inner join to only do drugs that I pick!!
    inner_join(drugs_moa_df, by = "pert_iname") %>%
    mutate(master_id = str_c(cell_id, pert_iname, pert_class, sep = "--")) %>%
    mutate(replicate_id = str_c(master_id, column_id, sep = "::")) %>%
    dplyr::select(master_id, replicate_id, everything())

  sub_obj <- my_obj %>%
  # non-standard evaluation to find column for group 
  # and filter
  filter(!!sym(grouping_var) %in% filter_vars) %>%
  group_by(master_id, replicate_id) %>%
  # some gene names are the same
  # so you need to make them unique and reference
  # the phosphosite later
  mutate(pr_gene_symbol = make.unique(temp_pr_gene_symbol, sep = "_"))

  #' print some summary information about filtered obj
  print_helper_info(sub_obj, grouping_var)

  full_splt_lst <- split(sub_obj, f = sub_obj[[sym(grouping_var)]])
  stopifnot(length(full_splt_lst) == length(filter_vars))

  output_dirs_lst <- create_od_str(filter_vars,
    output_directory,
    dataset_type, grouping_var
  ) 

  analysis_result_lst <- run_analysis(full_splt_lst,
    tie_method = "average", use_bootstrap = FALSE, use_parallel = TRUE
  )

  conn_clust_obj <- analysis_result_lst$output_results
  tbl_conn_clust_obj <- as_tibble(conn_clust_obj) %>%
    mutate(name = filter_vars) %>% # can do this bc always in order...
    rename(!!grouping_var := name) %>%
    pivot_longer(corr_lst:clust_lst, names_to = "match", values_to = "lst")

  res_paths_tbl <- suppressMessages(left_join(
    output_dirs_lst,
    tbl_conn_clust_obj
  ))
  stopifnot(nrow(res_paths_tbl) == nrow(tbl_conn_clust_obj))

  
  cache_objects(
    dir_tbl = res_paths_tbl,
    target_col = "lst", path_col = "path"
  )


  # TODO: PLOT PVCLUST
  # TODO: HEATMAPS
  # TODO: CLUSTER HIERARCHY

 
  message("Done.")

  message("Caching full obj...")
  complete_obj <- list(
    "data" = sub_obj,
    "corr_lst" = sample_corr_lst,
    "conn_clust_lst" = sample_conn_clust_obj
  )
  cache_output_dir <- file.path(specific_output_dir, "cache")
  dir.create(cache_output_dir, showWarnings = F, recursive = T)
  cache_fn <- file.path(cache_output_dir, qq("cache-@{Sys.Date()}.rds"))
  write_rds(x = complete_obj, file = cache_fn)
  message("Done.\n")
  return(complete_obj)
})

# TODO: diff ex, heatmaps dendros
