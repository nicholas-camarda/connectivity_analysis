# source the init.R file first
source(file.path("scripts", "init.R"))

cat("\nReading and merging data...")
#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
#' the values in this data are Z-scores!!!

#' TODO change this to [analysis_args.csv] when ready
analysis_fn <- file.path(data_directory, "test_args.csv")
analysis_dat_temp <- read_csv(analysis_fn, comment = "#", 
                              show_col_types = FALSE) %>%
  mutate_all(str_trim) %>%
  mutate(filter_vars = map(filter_vars, collect_args),
         exclude = map(exclude, collect_args)) %>%
  left_join(dir_tbl, by = "dataset_type")

my_data <- tibble(fns = dir(file.path(datasets_directory, "LINCS-data"), 
                            full.names = T, recursive = T)) %>%
  distinct() %>% 
  mutate(gct = map(
    .x = fns,
    .f = parse_gctx
  )) %>% 
  mutate(dataset_type = str_extract(string = fns, pattern = "GCP|P100")) %>%
  arrange(desc(dataset_type)) %>%
  mutate(ranked_gct = map(.x = gct, .f = cmapR::rank_gct, dim = "col"),
         mat = map(.x = gct, .f = cmapR::mat),
         ranked_mat = map(.x = ranked_gct, .f = cmapR::mat)) %>%
  mutate(data = map(.x = gct, .f = melt_gct)) %>%
  inner_join(analysis_dat_temp, by= "dataset_type")

drugs_moa_df <- create_my_drugs_df(ref_dir = references_directory)

# read in shared plated IDs from Srila's data
# run debug-dataset.R first!
shared_plate_id_df <- read_tsv(file.path(working_directory, "debug/shared_plate_ids.tsv"), 
                               show_col_types = FALSE) %>% 
  arrange(shared_plate_ids)
shared_plate_ids <- shared_plate_id_df$shared_plate_ids


# organize what was read in
my_data_lst <- my_data %>% 
  split(.$dataset_type) 
dataset_type_col <- names(my_data_lst)

message("Reading and summarizing data...")
my_data_obj_final <- tibble(dataset_type = dataset_type_col, 
                            # combine the data rows into a single dataset, and nest it
                            data = map2(my_data_lst, dataset_type, .f = read_and_summarize_data))

# grab drugs from BOTH P100 and GCP
my_temp_obj <- bind_rows(my_data_obj_final$data) %>% 
  ungroup() 

HUVEC_HAoSMC_perts <- my_temp_obj %>% 
  filter(cell_id %in% vascular_char_vec) %>% 
  ungroup() %>% 
  dplyr::select(pert_iname) %>% 
  .$pert_iname %>% 
  unique()
other_perts <- my_temp_obj %>% 
  ungroup() %>%
  filter(!(cell_id %in% vascular_char_vec)) %>% 
  distinct(pert_iname) %>% 
  .$pert_iname
my_perts <- intersect(HUVEC_HAoSMC_perts, other_perts)
all_drugs_mapping <- my_temp_obj %>% 
  distinct(pert_iname) %>%
  filter(pert_iname %in% my_perts)

# bind final object
analysis_dat <- inner_join(
  analysis_dat_temp,
  my_data_obj_final, # my_data_obj_final
  by = "dataset_type"
)


# analysis apply function
analysis_res <- apply(analysis_dat, 1, function(args) {
  # DEBUG: args = analysis_dat[1,]; 
  
  dataset_type <- tolower(args$dataset_type)
  grouping_var <- args$grouping_var
  filter_vars <- force_natural(args$filter_vars)
  if(any(is.na(force_natural(args$exclude)))) {
    exclude <- ""
  } else { 
    exclude <- force_natural(args$exclude)
  }
  my_obj <- args$data
  message("Excluding: ", exclude)
  
  # DEBUG:  my_obj <- force_natural(args$data)
  
  # this is convenient to filter the drugs out here... 
  # but probably not the most intelligent design
  # need to filter out by group of analysis, e.g. per kinase inhibitor, per epigenetics, etc
  HUVEC_HAoSMC_perts <- my_obj %>% 
    filter(cell_id %in% vascular_char_vec) %>% 
    ungroup() %>% 
    dplyr::select(pert_iname) %>% 
    .$pert_iname %>% 
    unique()
  
  other_perts <- my_obj %>% 
    ungroup() %>%
    filter(!(cell_id %in% vascular_char_vec)) %>% 
    distinct(pert_iname) %>% 
    .$pert_iname
  
  my_perts <- intersect(HUVEC_HAoSMC_perts, other_perts)
  
  
  sub_obj_temp <- my_obj %>%
    # non-standard evaluation to find column for group 
    # and filter
    # filter for only perts that show in both cancer and vascular
    filter(pert_iname %in% my_perts) %>%
    filter(!!sym(grouping_var) %in% filter_vars) %>%
    # if we are excluding perts, make sure this runs
    filter(!(pert_iname %in% exclude))
  
  #' print some summary information about filtered obj
  print_helper_info(sub_obj_temp, grouping_var)
  
  dir_name_df <- tibble(!!grouping_var := filter_vars) %>%
    mutate(new_filter_var = map_chr(filter_vars, str_c, qq("_excl_@{exclude}"))) 
  
  sub_obj <- sub_obj_temp %>% 
    left_join(dir_name_df) %>%
    ungroup() %>% 
    dplyr::select(-!!grouping_var) %>%
    rename(!!grouping_var := new_filter_var) %>% 
    suppressMessages()
  
  full_splt_lst <- split(sub_obj, f = sub_obj[[sym(grouping_var)]])
  stopifnot(length(full_splt_lst) == length(filter_vars))
  
  output_dirs_lst <- create_od_str(filter_vars = dir_name_df$new_filter_var,
                                   output_directory = output_directory, 
                                   dataset_type = dataset_type, 
                                   grouping_var = grouping_var
  ) %>%
    mutate(dirname_ = .[[1]],
           !!grouping_var := str_split(string = dirname_, pattern = "_", simplify = TRUE )[,1]); output_dirs_lst
  
  message("Plotting pert-cell distribution data...")
  debug_plot_directory <- file.path(output_directory, dataset_type, "plots", "summaries")
  dir.create(debug_plot_directory, recursive = T, showWarnings = F)
  
  prop_pert_df <- tibble(filter_vars = dir_name_df$new_filter_var,
                         data = full_splt_lst) %>%
    mutate(dirname_ = .[[1]],
           !!grouping_var := str_split(string = dirname_, pattern = "_", simplify = TRUE )[,1]) %>%
    mutate(prop_df = map(data, function(d){
      res <- d %>% group_by(cell_id, pert_iname) %>% 
        dplyr::summarize( n = n(), .groups = "keep") %>%
        group_by(cell_id) %>% 
        mutate(prop = n / sum(n))
      return(res)
    })) %>% 
    mutate(gplot = map(prop_df, function(d) {
      gplot <- ggbarplot(data = d, x = "cell_id", y="prop", fill = "pert_iname",
                         palette = "igv", ggtheme = theme_bw()) + 
        ggtitle("Proportion of therapies used in each cell type")
      return(gplot)
    })) %>%
    mutate(gplot_base_dir = debug_plot_directory,
           gplot_path = file.path(gplot_base_dir, str_c( dirname_, "--prop-of-drugs-per-cell.pdf")))
  
  walk2(as.list(prop_pert_df$gplot), as.list(prop_pert_df$gplot_path), .f = function(x,y) {
    ggsave(x, filename = y, device = 'pdf', width = 8, height = 10)
  })
  
  # check output dir for results
  # stop()
  cache_output_directory <- file.path(output_directory, dataset_type, "cache", grouping_var)
  directories_to_search <- file.path(cache_output_directory, unique(dir_name_df$new_filter_var))
  search_string <- "clust_clust_obj|conn_conn_tbl|corr_matrix|diffe_diffe_final_res"
  split_search_string <- str_split(string = search_string, pattern = "\\|", simplify = T)[1,]
  output_paths <- tibble(path = dir(directories_to_search, 
                                    full.names = T, recursive = T)) %>%
    mutate(temp_col = str_extract(string = path, pattern = str_c(filter_vars, collapse = "|"))) %>%
    rename(!!grouping_var := temp_col) %>%
    mutate(extract_path = str_extract(string = path,  
                                      pattern = search_string)) 
  
  output_path_check <- nrow(output_paths) == length(split_search_string)*length(directories_to_search); output_path_check
  
  # looking for 4 files per condition
  if (!output_path_check) {
    message("Did not detect cached output... Starting computation.")
    analysis_result_lst <- run_analysis(full_splt_lst,
                                        tie_method = "average", 
                                        use_bootstrap = FALSE, 
                                        use_parallel = TRUE
    )
    #' how to read pvclust results
    #' The plot should be read from bottom to top. There are three numbers around each node. 
    #' The number below each node specifies the rank of the cluster (here, from 1 to 13, i.e. 
    #' from the 1st generated cluster at the bottom to the 13th at the top). The two numbers 
    #' above each node indicate two types of p-values, which are calculated via two different 
    #' bootstrapping algorithms: AU and BP.1
    
    res_obj <- analysis_result_lst$output_results
    tbl_res_obj <- as_tibble(res_obj) %>%
      mutate(name = names(res_obj$clust_lst)) %>% 
      pivot_longer(corr_lst:diffe_lst, 
                   names_to = "match", values_to = "lst") %>%
      rename(!!grouping_var := name) %>%
      mutate(dirname_ = .[[1]],
             !!grouping_var := str_split(string = dirname_, pattern = "_", simplify = TRUE )[,1])
    
    res_paths_tbl_temp <- left_join(
      output_dirs_lst,
      tbl_res_obj
    ) %>% suppressMessages()
    
    stopifnot(nrow(res_paths_tbl_temp) == nrow(tbl_res_obj))
    
    # TODO: write ifelse to skip analysis if already cached
    cache_objects(
      dir_tbl = res_paths_tbl_temp,
      target_col = "lst", 
      path_col = "path"
    )
    
    res_paths_tbl <- suppressMessages(left_join(tbl_res_obj, 
                                                res_paths_tbl_temp %>% 
                                                  dplyr::select(-lst))) %>%
      dplyr::select(match, lst, path, everything())
    
  } else {
    message("Loading cached obj's...")
    res_paths_tbl_temp <- output_paths %>%
      mutate(match_temp = str_extract(string = path, 
                                      pattern = "clust_clust_obj--cluster_assignments--co_clust_bool|conn_conn_tbl--conn_median_mat--conn_tbl_median|corr_matrix--tibble|diffe_diffe_final_res--diffe_ggplot")) %>%
      na.omit() %>%
      mutate(match = str_c(str_extract(match_temp, "clust|conn|corr|diffe"), "_lst", sep = "")) 
    
    res_paths_tbl <- res_paths_tbl_temp %>% 
      mutate(lst = map(path, read_rds)) %>%
      dplyr::select(extract_path, path, !!grouping_var,  match_temp, match, lst, everything()) %>%
      mutate(match_temp = map_chr(extract_path, function(x) str_extract(string = x, pattern = "corr|conn|clust|diffe"))) %>%
      mutate(match = str_c(match_temp,  "_lst")) %>%
      dplyr::select(match, lst, path, !!grouping_var); res_paths_tbl
  }
  message("Done.")
  
  my_heatmap_and_dendro_obj_temp <- res_paths_tbl %>% 
    pivot_wider(id_cols = !!grouping_var, names_from = match, values_from = lst) %>%
    arrange(grouping_var) %>%
    mutate(base_path = file.path(output_directory, dataset_type),
           which_dat = dataset_type) %>%
    left_join(output_dirs_lst) %>%
    suppressMessages()
  
  # figure out how to get paths easily here
  # TODO FIX FIX FIX
  my_dendro_obj <- my_heatmap_and_dendro_obj_temp %>%
    mutate(dir_path = file.path(base_path, "plots", "prettydendros"),
           path = file.path(dir_path, str_c(dirname_, ".eps"))) %>%
    distinct(path, .keep_all = TRUE) %>%
    dplyr::select(-match)
  walk(unique(my_dendro_obj$dir_path), dir.create, showWarnings = F, recursive = T)
  
  apply(my_dendro_obj, 1, FUN = function(args_) {
    plot_pretty_dendrogram(args = args_, rotate_dendrogram = TRUE)
  })
  
  input_dat_wide <- tibble(dirname_ = names(full_splt_lst), 
                           input_data = full_splt_lst) %>%
    mutate(dirname_ = .[[1]],
           !!grouping_var := str_split(string = dirname_, pattern = "_", simplify = TRUE )[,1])
  
  my_heatmap_obj <- suppressMessages(left_join(my_heatmap_and_dendro_obj_temp,
                                               input_dat_wide)) %>%
    mutate(dir_path = file.path(base_path,"plots",  "heatmaps"),
           path = file.path(dir_path, str_c(dirname_, ".eps"))) %>%
    distinct(path, .keep_all = TRUE)%>%
    dplyr::select(-match) %>%
    mutate(grouping_var = grouping_var)
  
  walk(unique(my_heatmap_obj$dir_path), dir.create, showWarnings = F, recursive = T)
  
  apply(X = my_heatmap_obj, 1, FUN = function(args_) {
    plot_heatmap(args = args_)
  })
  
  # TODO move diffe plot out of other function
  diffe_ggplot_outputpath <- my_heatmap_and_dendro_obj_temp %>%
    mutate(dir_path = file.path(base_path, "plots", "volcano_plots"),
           path = file.path(dir_path, str_c(dirname_, ".eps"))) %>%
    distinct(path, .keep_all = TRUE) %>%
    dplyr::select(-match) %>%
    mutate(grouping_var = grouping_var)
  
  walk(diffe_ggplot_outputpath$dir_path, ~ dir.create(.x, recursive = T, showWarnings = F))
  
  ggplots_diffe <- purrr::transpose(diffe_ggplot_outputpath$diffe_lst)$diffe_ggplot
  walk2(as.list(ggplots_diffe), as.list(diffe_ggplot_outputpath$path), ~ ggsave(filename = .y, plot = .x, width = 14, height = 8))
  
  message("Done with that batch!!")
})


message("\nDone with everything!")
gc()

