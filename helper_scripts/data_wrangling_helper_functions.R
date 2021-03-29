
#' @note convenience function for creating output directory structures
#' @param output_directory from global env
#' @param dataset_type P100 or GCP
#' @param grouping_var pert_iname, pert_class,
#' etc for establishing sub directories
#' @return lst containing specific_output_dir (top level),
#' conncetivity_output_dir, and clust_output_dir
create_output_directory_str <- function(output_directory,
                                        dataset_type, grouping_var) {
  message("Creating top output directory: ")
  specific_output_dir <- file.path(output_directory, dataset_type, grouping_var)
  dir.create(specific_output_dir, showWarnings = FALSE, recursive = TRUE)
  message(qq("@{specific_output_dir}"))

  connectivity_output_dir <- file.path(specific_output_dir, "conn")
  dir.create(connectivity_output_dir, showWarnings = FALSE, recursive = TRUE)
  clust_output_dir <- file.path(specific_output_dir, "clust")
  dir.create(connectivity_output_dir, showWarnings = FALSE, recursive = TRUE)
  return(list(
    "specific_output_dir" = specific_output_dir,
    "connectivity_output_dir" = connectivity_output_dir,
    "clust_output_dir" = clust_output_dir
  ))
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
    transmute(pert_iname = Drug, pert_class = Class, pert_category = "cancer")

  cv_drug_moa_df <- read_excel(path = drug_classes_fn, sheet = 2) %>%
    filter(!is.na(`Drug (Generic)`)) %>%
    transmute(
      pert_iname = `Drug (Generic)`,
      pert_class = `Class`, pert_category = "cv"
    )

  drugs_moa_df <- bind_rows(cancer_drug_moa_df, cv_drug_moa_df) %>%
    mutate(pert_iname = tolower(pert_iname)) %>%
    distinct()
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
    select(fn, full_path_fn) %>%
    mutate(obj = map(full_path_fn, function(f) {
      res <- suppressMessages(parse_gctx(f))
      cat(".")
      return(res)
    })) %>%
    select(obj, everything()) %>%
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