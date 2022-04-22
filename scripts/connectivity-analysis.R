# source the init.R file firstso
source(file.path("scripts", "init.R"))

# analysis apply function
tic()
analysis_res <- apply(analysis_dat, 1, function(args) {
  # DEBUG: args = analysis_dat[1,];
  # DEBUG:  my_obj <- force_natural(args$data)
  
  dataset_type <- tolower(args$dataset_type)
  grouping_var <- args$grouping_var
  filter_vars <- force_natural(args$filter_vars)
  if (any(is.na(force_natural(args$exclude)))) {
    exclude <- ""
  } else {
    exclude <- force_natural(args$exclude)
  }
  my_obj <- args$data
  exclude_msg <- ifelse(exclude == "", "None", exclude)
 
  
  sub_obj_temp <- my_obj %>%
    # non-standard evaluation to find column for group
    # filter for only my_perts that show in both cancer and vascular, 
    # which is defined in init.R
    filter(pert_iname %in% my_perts) %>%
    filter(!!sym(grouping_var) %in% filter_vars) %>%
    # if we are excluding perts, make sure this runs
    filter(!(pert_iname %in% exclude))
  
  #' print some summary information about filtered obj
  print_helper_info(sub_obj_temp, grouping_var)
  message("Excluding perturbations?: ", exclude_msg)
  message()
  
  dir_name_df <- tibble(!!grouping_var := filter_vars) %>%
    mutate(new_filter_var = map_chr(
      filter_vars, .f = function(x) {
        if (exclude != ""){
          return(str_c(x, qq("_excl_@{exclude}"), sep = ""))
        } else {
          return(x)
        }
      }
    )); dir_name_df
  
  sub_obj <- sub_obj_temp %>%
    left_join(dir_name_df) %>%
    ungroup() %>%
    dplyr::select(-!!grouping_var) %>%
    rename(!!grouping_var := new_filter_var) %>%
    suppressMessages()
  
  full_splt_lst <- split(sub_obj, f = sub_obj[[sym(grouping_var)]])
  stopifnot(length(full_splt_lst) == length(filter_vars))
  
  output_dirs_lst <- create_od_str(
    filter_vars = dir_name_df$new_filter_var,
    output_directory = output_directory,
    dataset_type = dataset_type,
    grouping_var = grouping_var
  ) %>%
    mutate(
      dirname_ = .[[1]], # dirname is the first column of this df always
      # grouping var is correct only if '_excl_' is removed
      !!grouping_var := str_split(
        string = dirname_,
        pattern = "_", simplify = TRUE
      )[, 1]
    )
  output_dirs_lst
  
  # check output dir for results
  # stop()
  cache_output_directory <- file.path(output_directory, dataset_type, "cache", grouping_var)
  directories_to_search <- file.path(cache_output_directory, unique(dir_name_df$new_filter_var))
  search_string <- "clust_clust_obj|conn_conn_tbl|corr_matrix|diffe_diffe_final_res"
  split_search_string <- str_split(string = search_string, pattern = "\\|", simplify = T)[1, ]
  
  output_paths <- tibble(path = dir(directories_to_search,
                                    full.names = T, recursive = T)) %>%
    mutate(temp_col = str_extract(
      string = path,
      pattern = str_c(filter_vars, collapse = "|")
    )) %>%
    rename(!!grouping_var := temp_col) %>%
    mutate(extract_path = str_extract(
      string = path,
      pattern = search_string
    ))
  expected_files_length <- length(split_search_string) * length(directories_to_search)
  output_path_check <- nrow(output_paths) == expected_files_length
  output_path_check
  
  # looking for 4 files per condition
  if (!output_path_check) {
    message("Did not detect cached output... Starting computation.")
    analysis_result_lst <- run_analysis(full_splt_lst,
                                        tie_method = "average",
                                        use_bootstrap = FALSE,
                                        use_parallel = TRUE
    )
    
    # DEBUG:
    # lst = full_splt_lst; tie_method = "average"; use_bootstrap = FALSE; use_parallel = TRUE
    
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
                   names_to = "match", values_to = "lst"
      ) %>%
      rename(!!grouping_var := name) %>%
      # if we excluded anything from a class analysis, 
      # then the directory needs to be different
      mutate(
        dirname_ = .[[1]],
        !!grouping_var := str_split(string = dirname_, pattern = "_", simplify = TRUE)[, 1]
      )
    
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
    
    res_paths_tbl <- suppressMessages(left_join(
      tbl_res_obj,
      res_paths_tbl_temp %>%
        dplyr::select(-lst)
    )) %>%
      dplyr::select(match, lst, path, everything())
    
  } else {
    
    res_paths_tbl <- load_cached_objs(output_paths, 
                                      grouping_var = grouping_var)
    
  }
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-main.RData")
  # load("debug/debug_dat/debug-main.RData")
  # stop()
  
  my_heatmap_and_dendro_obj_temp <- res_paths_tbl %>%
    arrange(grouping_var) %>%
    mutate(
      base_path = file.path(output_directory, dataset_type),
      which_dat = dataset_type
    ) %>%
    left_join(output_dirs_lst %>% dplyr::select(-path)) %>%
    suppressMessages()
  
  my_dendro_obj <- my_heatmap_and_dendro_obj_temp %>%
    filter(match == "clust_lst") %>%
    mutate(
      dir_path = file.path(base_path, "plots", "prettydendros"),
      path = file.path(dir_path, str_c(dirname_, ".eps"))
    ) %>%
    distinct(path, .keep_all = TRUE)
  
  # create output directory
  walk(unique(my_dendro_obj$dir_path), dir.create, showWarnings = F, recursive = T)
  
  apply(my_dendro_obj, 1, FUN = function(args_) {
    plot_pretty_dendrogram(args = args_, rotate_dendrogram = TRUE)
  })
  
  input_dat_wide <- tibble(dirname_ = names(full_splt_lst), 
                           input_data = full_splt_lst) %>%
    mutate(dirname_ = .[[1]], 
           !!grouping_var := str_split(string = dirname_, 
                                       pattern = "_", 
                                       simplify = TRUE)[, 1])
  
  my_heatmap_obj <- left_join(
    my_heatmap_and_dendro_obj_temp,
    input_dat_wide
  ) %>%
    mutate(
      dir_path = file.path(base_path, "plots", "heatmaps"),
      # pdf saving is done in ploting function... too hard to rework now
      path = file.path(dir_path, str_c(dirname_, ".eps"))
    ) %>%
    # distinct(path, .keep_all = TRUE) %>%
    mutate(grouping_var = grouping_var) %>%
    suppressMessages() %>%
    filter(match %in% c("diffe_lst", "clust_lst")) %>%
    pivot_wider(id_cols = c(path, !!grouping_var, dirname_, base_path, which_dat, 
                            input_data, dir_path, grouping_var),
                names_from = match, 
                values_from = lst)
  
  # create output directory
  walk(unique(my_heatmap_obj$dir_path), dir.create, showWarnings = F, recursive = T)
  
  heatmap_res <- apply(X = my_heatmap_obj, 1, FUN = function(args_) {
    plot_heatmap(args = args_)
  })
  
  my_diffe_obj <- my_heatmap_and_dendro_obj_temp %>%
    filter(match == "diffe_lst") %>%
    mutate(
      dir_path_group = file.path(base_path, "plots", "volcano_plots", "faceted"),
      dir_path_single = file.path(base_path, "plots", "volcano_plots", "single"),
      path_group_eps = file.path(dir_path_group, str_c(dirname_, ".eps")),
      path_single_eps = file.path(dir_path_single, str_c(dirname_, ".eps")),
      path_group_pdf = file.path(dir_path_group, str_c(dirname_, ".pdf")),
      path_single_pdf = file.path(dir_path_single, str_c(dirname_, ".pdf"))
    ) %>%
    distinct(path_group_eps, path_single_eps, 
             path_group_pdf, path_single_pdf, .keep_all = TRUE) %>%
    mutate(grouping_var = grouping_var) %>%
    filter(match %in% c("diffe_lst", "clust_lst")) %>%
    pivot_wider(id_cols = c(path, !!grouping_var, dirname_, base_path, which_dat, 
                            dir_path_group, dir_path_single, grouping_var,
                            path_group_eps, path_single_eps,
                            path_group_pdf, path_single_pdf),
                names_from = match, 
                values_from = lst)
  
  # make the directories
  walk(
    my_diffe_obj$dir_path_group,
    ~ dir.create(.x, recursive = T, showWarnings = F)
  )
  walk(
    my_diffe_obj$dir_path_single,
    ~ dir.create(.x, recursive = T, showWarnings = F)
  )
  
  
  # plot things
  my_diffe_plots <- apply(my_diffe_obj, 1, FUN = function(args_) {
    plot_diffe_results(args = args_)
  }) %>%
    setNames(force_natural(my_diffe_obj[1])) %>%
    purrr::transpose()
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-diffe.RData")
  # load("debug/debug_dat/debug-diffe.RData")
  # stop()
  
  
  ggplots_diffe <- my_diffe_plots$diffe_full_g
  ggplots_diffe_singles <- my_diffe_plots$diffe_singles %>%
    purrr::flatten()
  
  # ugly part to correct paths
  base_path_chr <- unique(my_diffe_obj$base_path)
  dirname_chr <- unique(my_diffe_obj$dirname_)
  names_tbl <- tibble(names_ = names(ggplots_diffe_singles)) %>%
    mutate(new_name = ifelse(is.na(names_), lag(names_), names_)) %>%
    mutate(unique_name = make.unique(new_name, sep = "_")) %>%
    mutate(dir_path_group = file.path(base_path_chr, "plots", "volcano_plots", "faceted"),
           dir_path_single = file.path(base_path_chr, "plots", "volcano_plots", "single"), 
           path_group_eps = file.path(dir_path_group, str_c(unique_name, ".eps")),
           path_single_eps = file.path(dir_path_single, str_c(unique_name, ".eps")),
           path_group_pdf = file.path(dir_path_group, str_c(unique_name, ".pdf")),
           path_single_pdf = file.path(dir_path_single, str_c(unique_name, ".pdf")))
  names(ggplots_diffe_singles) <- names_tbl$unique_name
  
  # stop()
  # plot the complete faceted diffe plots
  write_diffe_objs_to_file(ggplots_diffe, ggplots_diffe_singles, names_tbl, lvl4_bool_data)
  
  # extract the diffe result and retrieve drugs that were used
  
  # str(heatmap_res, max.level = 2)
  temp <- purrr::transpose(heatmap_res) %>% flatten()
  if (length(temp) != 5) {
    # another transposition if it's nested
    temp <- purrr::transpose(temp)
  }
  df_drugs_used <- bind_rows(lapply(temp$col_annot, 
                                    FUN = function(l) l)) %>% 
    distinct(pert_iname, pert_class) %>%
    mutate(pert_class = str_split(pert_class, "_excl_", simplify = TRUE)[,1])
  
  drug_dose_info_df <- my_obj %>% 
    ungroup() %>%
    distinct(cell_id, pert_iname, pert_class, pert_dose, pert_dose_unit, det_plate) %>% 
    arrange(pert_iname, cell_id) %>%
    filter(pert_iname %in% df_drugs_used$pert_iname)
  drug_dose_info_df_distinct <- drug_dose_info_df %>% 
    distinct(pert_iname)
  
  dirname_char <- unique(dir_name_df$new_filter_var)
  drug_info_dir <- file.path(output_directory, dataset_type, "dataset_info")
  drug_dose_info_fn_name <- file.path(drug_info_dir, str_c(dirname_char, ".tsv"))
  drug_dose_info_fn_name_distinct <- file.path(drug_info_dir, str_c(dirname_char, "_distinct.tsv"))
  dir.create(drug_info_dir, showWarnings = F, recursive = T)
  
  write_tsv(drug_dose_info_df, file = drug_dose_info_fn_name)
  write_tsv(drug_dose_info_df_distinct, file = drug_dose_info_fn_name_distinct)
  
  
  message("Done with that batch!\n")
  
  return(list(my_dendro_obj, my_heatmap_obj, my_diffe_obj, res_paths_tbl) %>% setNames(c("dendro_res","heatmap_res", "diffe_res", "cache")))
})

message("Writing analysis res..")
names(analysis_res) <- str_c(analysis_dat$dataset_type, analysis_dat$filter_vars, analysis_dat$grouping_var, sep = "-")
write_rds(analysis_res, file = file.path(output_directory, "final_result.rds"), compress = "gz")

message("\nDone with everything!")
gc()
toc()

