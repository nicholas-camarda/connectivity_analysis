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
  mutate(filter_vars = map(filter_vars, collect_args)) %>%
  left_join(dir_tbl, by = "dataset_type")

# read in all 
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


# # plot z-score vs rank for the first 25 genes (rows)
# ranked_m <- mat(my_ds_rank_by_column)
# print(plot(ranked_m[1:25, ],
#            m[1:25, ],
#            xlab="rank",
#            ylab="differential expression score",
#            main="score vs. rank"))

# create drug df from excel file of drugs
drugs_moa_df <- create_my_drugs_df(ref_dir = references_directory)

# read in shared plated IDs from Srila's data
# run debug-dataset.R first!
shared_plate_id_df <- read_tsv("/Users/ncamarda/OneDrive - Tufts/phd/ws/proteomics/debug/shared_plate_ids.tsv", 
                               show_col_types = FALSE) %>% 
  arrange(shared_plate_ids)
shared_plate_ids <- shared_plate_id_df$shared_plate_ids


# organize what was read in
my_data_lst <- my_data %>% 
  split(.$dataset_type) 
dataset_type_col <- names(my_data_lst)
my_data_obj_final <- tibble(dataset_type = dataset_type_col, 
                            # combine the data rows into a single dataset, and nest it
                            data = map2(my_data_lst, dataset_type, .f = read_and_summarize_data))

analysis_dat <- inner_join(
  analysis_dat_temp,
  my_data_obj_final, # my_data_obj_final
  by = "dataset_type"
)

# analysis apply function
analysis_res <- apply(analysis_dat, 1, function(args) {
  # DEBUG: args = analysis_dat[1,]; 
  
  dataset_type <- args$dataset_type
  grouping_var <- args$grouping_var
  filter_vars <- force_natural(args$filter_vars)
  my_obj <- args$data
  
  # DEBUG:  my_obj <- force_natural(args$data)
  
  # this is convenient to filter the drugs out here... 
  # but probably not the most intelligent design
  # need to filter out by group of analysis, e.g. per kinase inhibitor, per epigenetics, etc
  vascular_cells <- c("HUVEC", "HAoSMC")
  HUVEC_HAoSMC_perts <- my_obj %>% 
    filter(cell_id %in% vascular_cells) %>% 
    ungroup() %>% 
    dplyr::select(pert_iname) %>% 
    .$pert_iname %>% 
    unique()
  other_perts <- my_obj %>% 
    ungroup() %>%
    filter(!(cell_id %in% vascular_cells)) %>% 
    distinct(pert_iname) %>% 
    .$pert_iname
  my_perts <- intersect(HUVEC_HAoSMC_perts, other_perts)
  
  sub_obj <- my_obj %>%
    # non-standard evaluation to find column for group 
    # and filter
    # filter for only perts that show in both cancer and vascular
    filter(pert_iname %in% my_perts) %>%
    filter(!!sym(grouping_var) %in% filter_vars)
  
  
  #' print some summary information about filtered obj
  print_helper_info(sub_obj, grouping_var)
  
  full_splt_lst <- split(sub_obj, f = sub_obj[[sym(grouping_var)]])
  stopifnot(length(full_splt_lst) == length(filter_vars))
  
  output_dirs_lst <- create_od_str(filter_vars,
                                   output_directory,
                                   dataset_type, grouping_var
  )
  
  
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
    pivot_longer(corr_lst:clust_lst, 
                 names_to = "match", values_to = "lst") %>%
    rename(!!grouping_var := name)
  
  res_paths_tbl <- suppressMessages(left_join(
    output_dirs_lst,
    tbl_res_obj
  ))
  stopifnot(nrow(res_paths_tbl) == nrow(tbl_res_obj))
  
  # TODO: write ifelse to skip analysis if already cached
  cache_objects(
    dir_tbl = res_paths_tbl,
    target_col = "lst", 
    path_col = "path"
  )
  
  # plot dendrograms and heatmaps
  
  my_heatmap_and_dendro_obj_temp <- res_paths_tbl %>% 
    filter(match == "clust_lst") %>%
    arrange(!!grouping_var) %>%
    rename(cluster_lst = lst)
  
  # needs to take into account split by grouping_var
  dendro_dirs <- tibble(path = file.path(output_directory,
                                         dataset_type, 
                                         "prettydendros"), 
                        which_dat = dataset_type)
  
  dir.create(dendro_dirs$path %>% unique(), 
             recursive = T, showWarnings = F)
  my_dendro_obj <- suppressMessages(bind_cols(my_heatmap_and_dendro_obj_temp %>% 
                                               dplyr::select(-path), 
                                             dendro_dirs))
  
  apply(my_dendro_obj, 1, FUN = function(args_) {
    plot_pretty_dendrogram(args = args_, rotate_dendrogram = TRUE)
  })
  
  # TODO: HEATMAPS
  heatmap_dirs <- tibble(path = file.path(output_directory,
                                          dataset_type, 
                                          "heatmaps"), 
                         which_dat = dataset_type)
  
  my_heatmap_obj <- suppressMessages(bind_cols(my_heatmap_and_dendro_obj_temp %>% 
                                dplyr::select(-path), 
                              heatmap_dirs) %>% 
    left_join(tibble(input_data = analysis_result_lst$input_data, 
                      temp_col = names(analysis_result_lst$input_data)) %>%
                 rename(!!grouping_var := temp_col) %>%
                 dplyr::select(!!grouping_var, input_data)))
                              
  plot_heatmap()
  apply(X = , 1, FUN = function(args_) {
    
  })
  message("Done.")
  
  # return(complete_obj)
})

