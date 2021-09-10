# source the init.R file first
options(error = recover)
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
shared_plate_id_df <- read_tsv(file.path(working_directory, "debug/shared_plate_ids.tsv"), 
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

# grab drugs
my_temp_obj <- bind_rows(my_data_obj_final$data) %>% 
  ungroup() 
vascular_cells <- c("HUVEC", "HAoSMC")
HUVEC_HAoSMC_perts <- my_temp_obj %>% 
  filter(cell_id %in% vascular_cells) %>% 
  ungroup() %>% 
  dplyr::select(pert_iname) %>% 
  .$pert_iname %>% 
  unique()
other_perts <- my_temp_obj %>% 
  ungroup() %>%
  filter(!(cell_id %in% vascular_cells)) %>% 
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
  # DEBUG: args = analysis_dat[3,]; 
  
  dataset_type <- tolower(args$dataset_type)
  grouping_var <- args$grouping_var
  filter_vars <- force_natural(args$filter_vars)
  if(any(is.na(force_natural(args$exclude)))){
    exclude <- ""
  } else { 
    exclude <- force_natural(args$exclude)
  }
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
  dir_name <- tibble(!!grouping_var := filter_vars) %>%
    mutate(new_filter_var = map_chr(filter_vars, str_c, qq("_excl_@{exclude}"))) 
  if (grouping_var == "pert_class"){
    sub_obj <- suppressMessages(sub_obj_temp %>% 
                                  left_join(dir_name) %>%
                                  dplyr::select(-!!grouping_var) %>%
                                  rename(!!grouping_var := new_filter_var))
  } else {
    sub_obj <- sub_obj_temp
  }
 
  
  full_splt_lst <- split(sub_obj, f = sub_obj[[sym(grouping_var)]])
  stopifnot(length(full_splt_lst) == length(filter_vars))
  
  
  output_dirs_lst <- create_od_str(dir_name$new_filter_var,
                                   output_directory,
                                   dataset_type, grouping_var
  )
  


  
  # check output dir for results
  # stop()
  directories_to_search <- file.path(output_directory, dataset_type, grouping_var, dir_name$new_filter_var)
  search_string <- "clust_clust_obj|conn_conn_tbl|corr_matrix|diffe_diffe_final_res"
  split_search_string <- str_split(string = search_string, pattern = "\\|", simplify = T)[1,]
  output_paths <- tibble(path = dir(directories_to_search, 
                                    full.names = T,recursive = T)) %>%
    mutate(temp_col = str_extract(string = path, pattern = str_c(filter_vars, collapse = "|"))) %>%
    rename(!!grouping_var := temp_col) %>%
    mutate(extract_path = str_extract(string = path,  
                                      pattern = search_string)) %>%
    na.omit() 
  
  output_path_check <- nrow(output_paths) == length(split_search_string); output_path_check
  
  # looking for 4 files
  if (!output_path_check) {
    message("Did not detected cached output... Starting computation.")
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
      rename(!!grouping_var := name)
    
    res_paths_tbl_temp <- suppressMessages(left_join(
      output_dirs_lst,
      tbl_res_obj
    ))
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
  # %>%

  
  my_heatmap_and_dendro_obj_temp <- res_paths_tbl %>% 
    pivot_wider(id_cols = !!grouping_var, names_from = match, values_from = lst)
  
  # needs to take into account split by grouping_var
  dendro_dirs <- tibble(path = file.path(output_directory,
                                         dataset_type, 
                                         "prettydendros"), 
                        which_dat = dataset_type)
  

  dir.create(dendro_dirs$path %>% unique(), 
             recursive = T, showWarnings = F)
  my_dendro_obj <- suppressMessages(bind_cols(my_heatmap_and_dendro_obj_temp,
                                              dendro_dirs)) 
  
  apply(my_dendro_obj, 1, FUN = function(args_) {
    plot_pretty_dendrogram(args = args_, rotate_dendrogram = TRUE)
  })
  
  # TODO: HEATMAPS
  heatmap_dirs <- tibble(path = file.path(output_directory,
                                          dataset_type, 
                                          "heatmaps"), 
                         which_dat = dataset_type)
  
  input_dat_wide <- tibble(input_data = full_splt_lst, 
                           temp_col = names(full_splt_lst)) %>%
    rename(!!grouping_var := temp_col)
  my_heatmap_obj <- suppressMessages(bind_cols(my_heatmap_and_dendro_obj_temp, 
                                               heatmap_dirs) %>% 
                                       left_join(input_dat_wide) %>%
    mutate(grouping_var = grouping_var))
  
  
  apply(X = my_heatmap_obj, 1, FUN = function(args_) {
    # plot_heatmap_and_clustering(args = args_)
    plot_heatmap(args = args_)
  })
  
  
  
  # stop("Debug")
  message("Done.")
  
  # return(complete_obj)
})


plot_heatmap_and_clustering <- function(args){
  
  # Debug
  stop()
  # args$input_data
  gene_names <- unique(args$input_data$pr_gene_symbol)
  mat_temp <- args$input_data %>%
    dplyr::select(replicate_id, pert_iname, cell_id, pr_gene_symbol, value) %>%
    pivot_wider(replicate_id:cell_id, pr_gene_symbol)  %>%
    ungroup()
  
  mat_no_replicates <- mat_temp %>% 
    group_by(pert_iname, cell_id) %>% 
    summarize(across(.cols = all_of(gene_names), ~ median(.x, na.rm = T)), .groups = "keep") %>%
    ungroup()
  
  # rnames <- mat_temp$replicate_id
  rnames <- mat_no_replicates$cell_id
  
  mat <- t(apply(mat_no_replicates %>%  # mat_temp
                   dplyr::select(all_of(gene_names)) %>% 
                   as.matrix(), 2, as.numeric))
  colnames(mat) <- rnames
  clustering_obj <- args$clust_lst$clust_obj$hclust
  
  cluster_tb <- left_join(mat_no_replicates, # mat_temp
                          tibble::enframe(args$clust_lst$cluster_assignments,
                                          name =  "cell_id", value = "base_clust_comp"), 
                          by="cell_id") ; cluster_tb
  cluster_tb_name <- left_join(cluster_tb, 
                               args$diffe_lst %>% 
                                 force_natural() %>% 
                                 distinct(base_clust_comp_name, base_clust_comp),
                               by = "base_clust_comp") %>%
    dplyr::select(cell_id, base_clust_comp, base_clust_comp_name)  %>% # replicate_id, 
    arrange(base_clust_comp) ; cluster_tb_name
    
  
  cluster_tb_name_distinct <- cluster_tb_name %>%
    distinct(cell_id, base_clust_comp, base_clust_comp_name) 
  
  cluster_tb_name_vec <- cluster_tb_name %>% 
    distinct(base_clust_comp, base_clust_comp_name) %>% 
    .$base_clust_comp_name
  
  clusters <- cluster_tb_name %>% 
    dplyr::select(base_clust_comp, base_clust_comp_name) %>%
    .$base_clust_comp
  
  # analytes_order <- args$diffe_lst[[1]] %>% arrange(logFC)
  # set color breaks at quantiles
  breaks <- sapply(rev(c(.05, .10,.30,.50,.70,.90, .95)), function(q) quantile(mat, q, na.rm = T, names = F))
  col_fun = colorRamp2(breaks = breaks, colors = brewer.pal(n = length(breaks), name = "RdYlBu"))
  # col_fun <- colorRamp2(c(-1,0,1), colors =c("red", "white", "blue"))
  # col_fun(seq(-3,3))
  library(yarrr)
  cols_annot1 <- piratepal(length.out = length(unique(cluster_tb$base_clust_comp_name)), trans = 0.1,  palette = "basel")
  cols_annot2 <- left_join(cluster_tb_name, tibble(color = piratepal(length.out = length(unique(cluster_tb$cell_id)), trans = 0.1,  palette = "basel")) %>%
    bind_cols(cluster_tb_name_distinct) %>%
      mutate(n = 1:n())) 
  
  #' sort by cell name and then by cluster
  names(cols_annot2$color) <- cols_annot2$n
  
  # cols_annot2_final <- 
  # fh = function(x) fastcluster::hclust(dist(x))
  ha <- HeatmapAnnotation(`cell id` = anno_simple(x = cols_annot2$n, col = cols_annot2$color))
  ht <- Heatmap(mat,
                na_col = "gray",
                col = col_fun,
                top_annotation = ha,
                name = 'foo', 
                cluster_rows = FALSE,
                cluster_columns = clustering_obj,
                column_labels = cluster_tb_name$cell_id,
                column_split = max(clusters), # clusters
                # cluster_columns = cluster_within_group(mat, factor = cluster_tb_name$base_clust_comp),
                border = TRUE, # outline heatmap with black line
                
                # set general plot params
                height = unit(10, "in"),
                width = unit(6, "in"),

                raster_device = "png",
                row_names_gp =  gpar(fontsize = 6, fontface = "bold"),
                heatmap_legend_param = list(direction = "horizontal"),
                rect_gp = gpar(col = "white", lwd = 0.7)); draw(ht)
                # # row_title = "Analyte",
                # row_title_side = "right",
                # row_title_rot = 0, 
                # show_row_names = TRUE,
                # cluster_rows = FALSE,
                # # cluster_columns = fh, 
                # # cluster_columns = TRUE,
                # column_names_rot = 90,
                # show_column_names = FALSE, 
                # cluster_columns = FALSE,
                # # column_labels = unique(cluster_tb_name_distinct$base_clust_comp_name),
                # # column_title = cluster_tb_name_vec,
                # # column_split = cluster_tb$base_clust_comp,
                # # cluster_columns = clustering_obj,
                # # column_sp
                # # cluster_columns = clustering_obj,
                # # column_split = cluster_tb$base_clust_comp, # split the columns by cluster! need an assignment per column
                # 
                # # column_title = " ", 
                # # no column title, otherwise defaults to groups
                # # column_title_gp = gpar(fill = c("white"), font = 3), # if split column title, then fill can take on mulitple vals
                # # column_title = clusters$group,
                # column_order = cluster_tb_name$replicate_id, # order the columns of this matrix with a new order, by u_cell_id (which matches matrix column names)
                # # cluster_columns = as.dendrogram(clustering_obj),
                # # column_dend_reorder = TRUE,
                # 
                # # column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                # # column_title = "Cell ID (unique)",
                # # column_title_side = "bottom",
                # 
                # border = TRUE, # outline heatmap with black line
                # 
                # # set general plot params
                # height = unit(6, "in"),
                # width = unit(10, "in"),
                # # height = unit(10, "mm")*nrow(reordered_mat),
                # # width = unit(4, "mm")*ncol(reordered_mat),
                # # height = unit(5, "mm")*nrow(reordered_mat),
                # # width = unit(5, "mm")*ncol(reordered_mat), # scale the shape of the heatmap
                # raster_device = "png",
                # row_names_gp =  gpar(fontsize = 6, fontface = "bold"),
                # heatmap_legend_param = list(direction = "horizontal"),
                # rect_gp = gpar(col = "white", lwd = 0.7)); draw(ht)
  
  
}
