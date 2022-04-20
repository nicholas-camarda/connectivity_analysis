
#' @note extract genes and sites based on list of genes, string
extract_sites <- function(analysis_dat, lst_of_genes) {
  
  res <- analysis_dat %>% 
    filter(dataset_type == "P100") %>%
    unnest(cols = c(data)) %>%
    dplyr::select(non_unique_pr_gene_symbol, mark) %>%
    distinct() %>%
    filter(non_unique_pr_gene_symbol %in% lst_of_genes) %>%
    rename(pr_gene_symbol = non_unique_pr_gene_symbol); res
  write_tsv(res, "~/Downloads/sites.tsv")
  return(res)
}


#' @note write diffe objects to file
#' @param ggplots_diffe list of faceted ggplots to plot
#' @param ggplots_diffe_singles list of single ggplots to plot
#' @param names_tbl file name and ggplot mapping table
write_diffe_objs_to_file <- function(ggplots_diffe, ggplots_diffe_singles, names_tbl, lvl4_bool_data){
  
  if (lvl4_bool_data) {
    width_ <- 8
    height_ <- 8
  } else {
    # could change this if you wanted
    width_ <- 14
    height_ <- 8
  }
  walk2(
    as.list(ggplots_diffe),
    as.list(names_tbl$path_group_eps),
    ~ ggsave(filename = .y, plot = .x, width = width_, height = height_, device = cairo_ps)
  )
  
  walk2(
    as.list(ggplots_diffe),
    as.list(names_tbl$path_group_pdf),
    ~ ggsave(filename = .y, plot = .x, width = width_, height = height_)
  )
  
  # plot the singles
  # for (i in )
  walk2(
    as.list(ggplots_diffe_singles), 
    as.list(names_tbl$path_single_eps),
    ~ ggsave(filename = .y, plot = .x, width = width_, height = height_, device = cairo_ps)
  )
  
  walk2(
    as.list(ggplots_diffe_singles), 
    as.list(names_tbl$path_single_pdf),
    ~ ggsave(filename = .y, plot = .x, width = width_, height = height_)
  )
}

#' @note plot missing values
#' @param sub_obj sub_obj filtered 
plot_missing <- function(sub_obj){
  wide <- sub_obj %>%
    dplyr::select(replicate_id, pert_iname, cell_id, pr_gene_symbol, value) %>%
    pivot_wider(names_from =  pr_gene_symbol, id_cols = replicate_id:cell_id )
  summary_ <- wide %>% 
    group_by(replicate_id, pert_iname, cell_id) %>% 
    miss_var_summary() %>%
    summarize(mean_pct_pres = 100-mean(pct_miss)) %>%
    # filter(mean_pct_pres >= 70) %>%
    ungroup() %>%
    arrange(mean_pct_pres) %>%
    group_by(pert_iname, cell_id) %>%
    mutate(n = 1:n()) %>% # number teh replicates
    mutate(plot_id = str_c(cell_id, pert_iname,n, sep = "-")) %>%
    ungroup(); summary_
  
  ggbarplot(summary_ %>% slice(1:50), x = "plot_id", y = "mean_pct_pres", ggtheme = theme_bw()) +
    theme(axis.text.x = element_text(angle = 90, vjust = -0.0001))
}


find_duplicates <- function(tbl) {
  #' only works on pivot_wider objects
  apply(tbl[,-1], 2, function(x) {
    if(any(which(x > 1))){
      return(x)
    }
  })
}



#' @note use this in an apply call to merge all data
#' @param l list of already-read-in GCT files
#' @param dtype_ the dataset type of the GCT, e.g. P100 or GCP
#' @return a merged object that contains all P100 or GCP data
read_and_summarize_data <- function(l, dtype_) {
  # l <- my_data_lst[[1]]; dtype_ = "GCP"
  # l <- my_data_lst[[2]]; dtype_ = "P100"
  # message(dtype_)
  res_temp <- data.table::rbindlist(l$data, fill = TRUE, use.names = TRUE) %>%
    as_tibble() %>%
    mutate(which_dat = dtype_) %>%
    mutate(
      cell_id = ifelse(cell_id == "Pericytes", "Pericyte", cell_id),
      pert_iname = tolower(pert_iname)
    ) %>%
    mutate(pert_iname = ifelse(pert_iname == "dmso", toupper(pert_iname), pert_iname)) %>%
    mutate(pert_iname = str_trim(str_replace(pert_iname, pattern = "\\(chembl1797936\\)", ""), "both")) %>%
    # do not include MCF10A
    dplyr::filter(cell_id != "MCF10A") %>%
    mutate(det_normalization_group_vector = as.character(det_normalization_group_vector)) %>%
    ## inner join to only do drugs that I pick!!
    left_join(drugs_moa_df, by = "pert_iname") %>%
    dplyr::rename( 
      row_id = id.x,
      column_id = id.y,
    ) %>%
    group_by(cell_id, pert_iname) %>%
    mutate(master_id = str_c(cell_id, pert_iname, pert_class, sep = "--")) %>%
    # this is critical. must make replicate id by plate!
    mutate(replicate_id = str_c(master_id, column_id, det_plate,  sep = "::")) %>%
    dplyr::select(master_id, replicate_id, everything()) %>%
    # group_by(replicate_id, pr_gene_symbol) %>%
    distinct(replicate_id, pr_gene_symbol, value, .keep_all = TRUE) 
  message("Removed MCF10A cells...")
  
  peptide_mod_dir <- file.path(data_directory, "peptide_modification_data")
  dir.create(peptide_mod_dir, showWarnings = F, recursive = T)
  if ("pr_gcp_histone_mark" %in% colnames(res_temp)) {
    # no need to do unique names for the histones...
    
    res <- res_temp %>%
      mutate(pr_gcp_histone_mark = str_trim(pr_gcp_histone_mark,"both")) %>%
      mutate(mark = pr_gcp_histone_mark, 
             pr_gene_symbol_temp = pr_gene_symbol,
             pr_gene_symbol = mark) %>% 
      dplyr::select(master_id, replicate_id, column_id, 
                    row_id, value, pr_gene_symbol, mark, everything()) %>%
      ungroup()
    
    write_tsv(res %>% distinct(pr_gene_symbol, pr_gcp_histone_mark), 
              file = file.path(peptide_mod_dir, "GCP-genes_with_marks_and_peptides.tsv"))
    
  } else {
    
    unique_gene_names_df <- res_temp %>% 
      ungroup() %>% 
      distinct(pr_gene_symbol, pr_p100_phosphosite, pr_p100_modified_peptide_code) %>% 
      rename(mark = pr_p100_phosphosite) %>% 
      # group_by(pr_gene_symbol, pr_p100_phosphosite) %>%
      mutate(pr_gene_symbol_u = make.unique(str_c("p", str_c(mark, pr_gene_symbol, sep = " "), 
                                                  sep = ""), 
                                            sep = "_")) %>%
      ungroup(); unique_gene_names_df
    
    write_tsv(unique_gene_names_df, 
              file = file.path(peptide_mod_dir, "P100-genes_with_phospho_marks_and_peptides.tsv"))
    
    res_temp2 <- res_temp %>%
      rename(mark = pr_p100_phosphosite) %>%
      left_join(unique_gene_names_df, 
                by = c("pr_gene_symbol", "mark", "pr_p100_modified_peptide_code")) %>%
      rename(non_unique_pr_gene_symbol = pr_gene_symbol) %>%
      rename(pr_gene_symbol = pr_gene_symbol_u); res_temp2
    
    res <- res_temp2 %>%
      dplyr::select(master_id, replicate_id, row_id, column_id, 
                    value, pr_gene_symbol, mark, everything()) %>%
      ungroup()
    
  }
  
  # save(list = ls(all.names = TRUE), file = "debug/debug_dat/debug-summary.RData")
  # load("debug/debug_dat/debug-summary.RData")
  # stop()
  
  message("Plotting summary data across pr_gene_symbol...")
  g <- ggplot(res) +
    geom_boxplot(aes(x=pr_gene_symbol, y=value)) + # 
    theme_bw(base_size = 7) +
    theme(axis.text.x = element_text(size = rel(1.5), angle=90, hjust=1, vjust=1)) +
    ggtitle(dtype_) 
  
  plot_dir <- file.path(output_directory, dtype_, "summary")
  dir.create(plot_dir, recursive = T, showWarnings = F)
  ggsave(g, filename = file.path(plot_dir, qq("@{dtype_}-gene_values-all.png")),
         width = 20, height = 10)
  
  return(res)
}


#' @note anonymous function for writing rds objects to file
#' @param obj to write to rds
#' @param dir_ full path to output directory, on which filename is made
anon_write_to_file <- function(obj, dir_) {
  dir.create(dir_, recursive = T, showWarnings = F)
  type_ <- str_c(names(obj), collapse = "--")
  fn <- file.path(dir_, qq("@{basename(dir_)}_@{type_}.rds"))
  write_rds(x = obj, file = fn, compress = "gz")
}

#' @note function to cache objects using walk2
#' @param dir_tbl df containing a column of list-obj to write,
#' and a column of paths associated with their destinations
#' @param target_col character string identifying obj-list column
#' @param path_col character string identifying column with paths for 
#' associated obj-list column
cache_objects <- function(dir_tbl, target_col = "lst",
                          path_col = "path") {
  message("Caching...")
  # write to dirs
  walk2(
    .x = dir_tbl %>% pluck(target_col),
    .y = dir_tbl %>% pluck(path_col),
    .f = anon_write_to_file
  )
  message("Done.")
}

#' @note create cluster assignments df for plotting
#' @param ca cluster assignments NAMED array!! this is the "cut_trees" field in diff_ex
#' @param cluster_column_id name of output cluster id
create_ca_df <- function(ca, cluster_column_id = "cluster") {
  t <- tibble::enframe(ca, name = "cell_id", value = cluster_column_id) %>%
    distinct()
  return(t)
}

#' @note convenience function for creating output directory structures
#' @param filter_vars from [filter_vars] 
#' @param output_directory from global env
#' @param dataset_type P100 or GCP
#' @param grouping_var pert_iname, pert_class,
#' @param dirs_to_make vector of dirs under this nest
#' etc for establishing sub directories
#' @return lst containing specific_output_dir (top level),
#' conncetivity_output_dir, and clust_output_dir
create_od_str <- function(filter_vars = c("Epigenetic", "Kinase Inhibitor"),
                          output_directory = "~/output",
                          dataset_type  = "P100", grouping_var = "pert_class", 
                          dirs_to_make = c("corr", "conn", "clust", "diffe")) {
  message("Creating top output directory: ")
  
  specific_output_dir <- file.path(output_directory, 
                                   tolower(dataset_type), 
                                   "cache", 
                                   grouping_var)
  
  dir.create(specific_output_dir, showWarnings = FALSE, recursive = TRUE)
  message(qq("@{specific_output_dir}"))
  
  message("Creating subdirs: ")
  msg <- str_c(dirs_to_make, collapse = ' | ')
  message(qq("@{msg}"))
  
  output_dirs_lst <- map(
    filter_vars,
    function(c) {
      map_chr(dirs_to_make, function(d) {
        file.path(specific_output_dir, c, d)
      })
    }
  ) %>%
    setNames(filter_vars) %>%
    bind_rows(.) %>%
    pivot_longer(everything(),
                 names_to = grouping_var, values_to = "path"
    ) %>%
    mutate(match = str_c(basename(path), "_lst"))
  
  walk(.x = output_dirs_lst$path, .f = dir.create,
       showWarnings = FALSE, recursive = TRUE
  )
  return(output_dirs_lst)
}

#' @note collects arguments from analysis_args.csv
#' @param str_chr presumably a line from .csv that needs to be split
#' @return split
collect_args <- function(str_char) {
  res <- as.vector(str_split(string = str_char, pattern = "--", simplify = T))
  return(res)
}

# print helper information
print_helper_info <- function(sub_obj, grouping_var) {
  summary_read_in <- sub_obj %>%
    ungroup() %>%
    distinct(!!sym(grouping_var))
  
  message("Data filtered down to contain perts/classes: ")
  print(summary_read_in)
  message("Cell IDs in this analysis: ")
  print(sub_obj$cell_id %>% unique())
}


#' @note convenience function to 
#' help with annoying list obj with purrr
#' @param obj some list-column obj
force_natural <- function(obj) {
  if (is.list(obj)) {
    return(obj[[1]])
  } else {
    return(obj)
  }
}

#' @note expects ref_dir to contain a file named [Drug Glossary_edited.xlsx]
#' @param ref_dir REFERENCES_DIRECTORY
create_my_drugs_df <- function(ref_dir = REFERENCES_DIRECTORY) {
  drug_classes_fn <- file.path(ref_dir, "Drug Glossary_edited.xlsx")
  cancer_drug_moa_df <- read_excel(path = drug_classes_fn, sheet = 1) %>%
    filter(!is.na(Drug)) %>%
    mutate(pert_iname = Drug, pert_class = Class, pert_category = "cancer")
  
  cv_drug_moa_df <- read_excel(path = drug_classes_fn, sheet = 2) %>%
    filter(!is.na(`Drug (Generic)`)) %>%
    mutate(
      pert_iname = `Drug (Generic)`,
      pert_class = `Class`, pert_category = "cv"
    )
  
  drugs_moa_df <- bind_rows(cancer_drug_moa_df, cv_drug_moa_df) %>%
    mutate(pert_iname = tolower(pert_iname)) %>%
    distinct() %>%
    mutate(pert_iname = ifelse(pert_iname == "dmso", toupper(pert_iname), pert_iname))
  return(drugs_moa_df)
}


#' @note anonymous function GCT merge for merging two gct's
#' @param dtf1 gctx1
#' @param dtf2 gctx2
anonymous_gct_merge <- function(dtf1, dtf2) {
  res <- merge_gct(dtf1, dtf2, dim = "column", matrix_only = FALSE)
  return(res)
}

#' @note extraction function for read_and_merge_gcts
#' @param chr is a string to split
#' @param idx is the position of the split string to return, defaulting to 1
#' @param sep_ is the pattern on which to split, defaulting to '-'
my_extract <- function(chr, idx = 1, sep_ = "-") {
  res <- str_split(string = chr, pattern = sep_, simplify = T)[, idx]
  return(res)
}

#' @note reads in data from RAW GCT folder
#' @param parent_dir_fn folder path for all the GCTs
#' @param dataset_grp P100 or GCT, depending on datatype
read_and_merge_gcts <- function(parent_dir_fn,
                                dataset_grp = "P100") {
  # NOTE: read all the data from RAW GCT directory and
  #' merge everything using library functions, not excel
  
  # TEST: parent_dir_fn <- file.path(DATASETS_DIRECTORY, "RAW GCT")
  all_data_fns <- tibble(parent_dir = parent_dir_fn) %>%
    mutate(fn = map(parent_dir, function(pr) list.files(pr))) %>%
    unnest(cols = c(fn)) %>%
    mutate(full_path_fn = map2_chr(fn, parent_dir, function(fn, pd) {
      return(file.path(pd, fn))
    })) %>%
    mutate(ext = map_chr(fn, function(fn) {
      res <- str_split(string = fn, pattern = "\\.", simplify = T)
      len <- length(res)
      return(res[len])
    })) %>%
    filter(ext == "gct") %>%
    arrange(full_path_fn)
  
  data <- all_data_fns %>%
    dplyr::select(fn, full_path_fn) %>%
    mutate(obj = map(full_path_fn, function(f) {
      res <- suppressMessages(parse_gctx(f))
      cat(".")
      return(res)
    })) %>%
    dplyr::select(obj, everything()) %>%
    mutate(
      cell_type = map_chr(fn, my_extract, idx = 1),
      dataset_type = map_chr(fn, my_extract, idx = 2)
    )
  cat("Done.\n")
  
  p100_data <- data %>% filter(dataset_type == dataset_grp)
  data_lst <- as.list(p100_data$obj) %>% setNames(p100_data$fn)
  # https://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list
  #' recursively merge list
  merged_obj <- suppressMessages(Reduce(
    f = anonymous_gct_merge,
    x = data_lst
  ))
  return(merged_obj)
}


