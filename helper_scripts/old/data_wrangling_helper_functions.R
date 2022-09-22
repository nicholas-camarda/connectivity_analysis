
#' @note print path name to easily copy into file.path(), universal paths
#' @param path any path character string
print_path <- function(path) {
  options(useFancyQuotes = FALSE)
  res <- sapply(strsplit(path, "/|\\\\"), function(x) toString(sQuote(x, q = FALSE)))
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
  p100_data_lst <- as.list(p100_data$obj) %>% setNames(p100_data$fn)
  # https://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list
  #' recursively merge list
  merged_p100_obj <- suppressMessages(Reduce(
    f = anonymous_gct_merge,
    x = p100_data_lst
  ))
  return(merged_p100_obj)
}


#' @note get a combined dataframe ...
#' @param mdf is a dataframe of column-dataframes
get_combined_tbl <- function(mdf) {
  mdf2 <- mdf %>%
    mutate(dataframe_for_analytes = map(t_dataframe, function(d) {
      as.data.frame(d) %>% rownames_to_column(var = "analyte")
    })) %>%
    mutate(dataframe_for_samples = map(dataframe, function(d) {
      as.data.frame(d) %>% rownames_to_column(var = "samples")
    }))


  # plot how many na's per analyte, left join by analyte
  # weird that bind_cols doesn't work
  analytes_across_samples_df <- plyr::join_all(mdf2$dataframe_for_analytes,
    type = "left",
    by = "analyte"
  ) %>%
    as_tibble()

  key_idx_cell_id <- key_idx_special %>%
    select(well, cell_id) %>%
    mutate(nu_cell_id = map_chr(cell_id, function(chr) str_split(string = chr, pattern = "-", simplify = T)[, 1])) %>%
    ungroup()

  samples_across_analytes_unfiltered_df <- bind_rows(mdf2$dataframe_for_samples) %>%
    as_tibble() %>%
    left_join(key_idx_cell_id %>%
      select(well, cell_id, nu_cell_id), by = c("samples" = "well")) %>%
    select(cell_id, nu_cell_id, everything(), -samples)

  #' @note then filter the analytes out that we discovered had too many NAs
  samples_across_analytes_df <- samples_across_analytes_unfiltered_df %>%
    select(cell_id, nu_cell_id, everything())

  return(list(analytes_across_samples_df, samples_across_analytes_df) %>% set_names("row_analytes", "row_samples"))
}


#' @note log pvalue to pval
#' @param log_p_val
logp_to_p <- Vectorize(function(logp_val) {
  10^-logp_val
})

#' @note convert a windwos path to a mac/unix acceptable R path
#' @param path any path
my.file.path <- function(path) {
  return(chartr("\\", "/", path))
}

#' @note load gcp dataset
#' @param first_gen_dataset path to first gen dataset in [GCT format]
load_gcp <- function(first_gen_dataset) {
  gcp_r <- read_tsv(first_gen_dataset, skip = 2)
  gcp_a <- gcp_r[c(1, 4, 10, 14, 21:nrow(gcp_r)), c(1, 4, 9:ncol(gcp_r))] %>%
    mutate(pr_gene_symbol = make.unique(pr_gcp_histone_mark, sep = "_")) %>%
    select(-pr_gcp_histone_mark) %>%
    select(id, pr_gene_symbol, everything())
  return(gcp_a)
}

#' @note load p100 dataset
#' @param first_gen_dataset path to first gen dataset in [GCT format]
load_p100 <- function(first_gen_dataset) {
  p100_r <- read_tsv(first_gen_dataset, skip = 2)
  p100a <- p100_r[c(1, 3, 11, 15, 23:nrow(p100_r)), c(1, 3, 13:ncol(p100_r))] %>%
    mutate(pr_gene_symbol = make.unique(pr_gene_symbol, sep = "_")) # remove metadata columns/rows
  return(p100a)
}

load_p100_alt <- function(first_gen_dataset) {
  p100_r <- read_tsv(first_gen_dataset, skip = 2)
  p100a <- p100_r[c(1, 3, 14, 18, 29:nrow(p100_r)), c(1, 3, 13:ncol(p100_r))] %>%
    mutate(pr_gene_symbol = make.unique(pr_gene_symbol, sep = "_")) # remove metadata columns/rows
  return(p100a)
}


#' @note load in data
#' @note more info about the .GCT file format here: https://clue.io/connectopedia/gct_format
#' @param data data frame that has been processed by [load_p100] or [load_gcp] functions (ABOVE)
#' @param drug_classes_fn excel file for drug classifications
load_data <- function(data, drug_classes_fn = "~/Drug Glossary_edited.xlsx", dataset = "p100",
                      filter_by_vasc_perts = FALSE) {

  # DEBUG:
  # data <- p100b
  # drug_classes_fn = file.path(references_directory, "Drug Glossary_edited.xlsx")

  cancer_drug_moa_df <- read_excel(path = drug_classes_fn, sheet = 1) %>%
    filter(!is.na(Drug)) %>%
    transmute(pert_iname = Drug, drug_class = Class, drug_cat = "cancer")
  dim(cancer_drug_moa_df)
  cv_drug_moa_df <- read_excel(path = drug_classes_fn, sheet = 2) %>%
    filter(!is.na(`Drug (Generic)`)) %>%
    transmute(pert_iname = `Drug (Generic)`, drug_class = `Class`, drug_cat = "cv")
  dim(cv_drug_moa_df)
  drugs_moa_df <- bind_rows(cancer_drug_moa_df, cv_drug_moa_df) %>%
    mutate(pert_iname = tolower(pert_iname)) %>%
    distinct()

  message(qq("Loading dataset: @{dataset}"))
  # p100a <- load_p100(first_gen_dataset = file.path(DATASETS_DIRECTORY, "P100","P100 All Cell Lines.gct"))
  # data <- p100a
  # in GCTs, Typically, each column represents a specific experiment (e.g. treatment of cell line MCF7 with a small-molecule drug)
  # and each row represents features (e.g. genes) that are measured in the assay.

  wells <- colnames(data)[3:ncol(data)]
  unique_wells_n <- length(unique(wells))
  message(qq("Unique experimental conditions (wells): @{unique_wells_n}"))

  cell_ids <- data[1, ] %>%
    gather(well, cell_id, wells) %>%
    select(well, cell_id)

  cell_cat <- cell_ids %>%
    distinct(cell_id) %>%
    mutate(cell_cat = ifelse(cell_id %in% c("HUVEC", "HAoSMC"), "vascular",
      ifelse(cell_id %in% c("NPC", "MCF10A"), "other", "cancer")
    ))

  message(qq("\nUnique cells: "))
  print(cell_cat)

  # map pert_id to well number
  det_norm_group_vect <- data[2, ] %>%
    gather(well, det_normalization_group_vector, wells) %>%
    select(well, det_normalization_group_vector) %>%
    rename(grp = det_normalization_group_vector) %>%
    distinct()
  dim(det_norm_group_vect)

  replicates <- data[3, ] %>%
    gather(well, pert_batch_internal_replicate, wells) %>%
    select(well, pert_batch_internal_replicate) %>%
    rename(rep = pert_batch_internal_replicate) %>%
    distinct()
  dim(replicates)

  message("\nTable of replicates by normalization group: ")
  dn <- det_norm_group_vect$grp
  rp <- replicates$rep
  print(table(rp, dn))

  message(qq("Observations (experimental conditions, e.g. cell + drug) by normalization group:")) # total 'observations' in each replicate
  print(colSums(table(rp, dn)))
  message(qq("Total observations (rowSum): @{sum(colSums(table(rp, dn)))}")) # sample size - 1904 for gcp
  # datb %>% filter(rep == 6, grp == "1,2")
  # >   table(dn,rp)
  # rp
  # dn        1     2     3     4     5     6
  # 1   27440 27048 27195     0     0     0
  # 1,1   931   931   931     0     0     0
  # 1,2  1029  1029   686    49    49    49
  # 1,3   588   588   588     0     0     0
  # 2     637   882   686     0     0     0
  # 2,3   539   539   882     0     0     0

  # map pert_id to well number
  pert_names <- data[4, ] %>%
    gather(well, pert_iname, wells) %>%
    select(well, pert_iname) %>%
    mutate(
      pert_iname = tolower(pert_iname),
      pert_iname = map_chr(pert_iname, function(x) str_split(string = x, pattern = "_", simplify = T)[, 1])
    )

  # # just for safe keeping, in case we want later
  phosphosite_long_df <- data %>%
    slice(-c(1:4)) %>%
    select(pr_gene_symbol, id, 3:ncol(data)) %>%
    gather(3:ncol(data), key = "well", value = pex) %>%
    select(well, pr_gene_symbol, pid = id) %>%
    distinct(pr_gene_symbol, pid)
  dim(phosphosite_long_df)

  # basically check if gcp
  if (!any(str_detect(pattern = "BI100", string = unique(phosphosite_long_df$pid)))) {
    phosphosite_long_df <- phosphosite_long_df %>% mutate(phosphosite = str_split(pid, "_", simplify = T)[, 3])
  } else {
    phosphosite_long_df <- phosphosite_long_df %>% rename(phosphosite = pid)
  }
  dim(phosphosite_long_df)

  # make "long" dataset - gene, well, cell_id, expression
  datb <- data %>%
    slice(-c(1, 2, 3, 4)) %>% # take out cell_id, det_normalization_group_vector, pert_batch_internal_replicate, and pert_iname
    select(pr_gene_symbol, wells) %>%
    gather(well, pex, wells) %>%
    left_join(cell_ids, by = "well") %>% # add cell_Id back in correctly
    left_join(pert_names, by = "well") %>%
    left_join(det_norm_group_vect, by = "well") %>%
    left_join(replicates, by = "well") %>%
    left_join(phosphosite_long_df, by = c("pr_gene_symbol")) %>% # , "well"
    left_join(cell_cat, by = "cell_id") %>%
    arrange(cell_id, pr_gene_symbol, pert_iname) %>%
    mutate(pex = as.numeric(pex)) %>%
    left_join(drugs_moa_df, by = "pert_iname")
  dim(datb)

  # take median of replicates for simplified dataset
  # filter(pert_iname == "1271738-62-5", cell_id == "A375", grp == "1,3")

  #' @note maybe to-do: should check to see whether replicates are any good...?
  datbb <- datb %>%
    # testbb <- datb %>%
    group_by(pr_gene_symbol, cell_id, pert_iname) %>%
    dplyr::summarize(pex_m = median(pex, na.rm = T), .groups = "keep") %>% # median to collapse replicates
    ungroup() %>%
    left_join(drugs_moa_df, by = "pert_iname") %>%
    left_join(cell_cat, by = "cell_id") %>%
    left_join(phosphosite_long_df, by = c("pr_gene_symbol"))

  # to get a dataset from a single well (e.g) cell and drug condition, you need to filter down to the well, pert, cell
  # heat <- datb %>%
  #   select(pr_gene_symbol, pert_iname, pex_m) %>%
  #   arrange(pert_iname) %>%
  #   mutate(unique_id = str_c("i",1:nrow(.))) %>%
  #   group_by(pert_iname)

  cell_lines_identities <- datbb %>%
    distinct(cell_id, cell_cat) %>%
    mutate(is_cancer = ifelse(cell_cat == "cancer", T, F)) %>%
    select(-cell_cat)
  cell_lines_identities

  key_fn <- qq("~/Downloads/master_tsv-@{dataset}-@{filter_by_vasc_perts}.tsv")
  message(qq("Writing master key file to @{key_fn} for @{dataset}"))
  write_tsv(datb, key_fn)


  # give cell line ids to final dat. datc because we have collapsed replicates by median
  if (filter_by_vasc_perts) {

    # get only drugs that were tested in CV cell lines
    tested_drugs_in_cv_cells <- datbb %>%
      select(-drug_class, -drug_cat, -cell_cat, -phosphosite) %>%
      filter(cell_id %in% c("HUVEC", "HAoSMC"))
    tested_drugs_in_smcs <- tested_drugs_in_cv_cells %>%
      filter(cell_id == "HAoSMC") %>%
      .$pert_iname %>%
      unique()
    tested_drugs_in_ecs <- tested_drugs_in_cv_cells %>%
      filter(cell_id == "HUVEC") %>%
      .$pert_iname %>%
      unique()

    perturbation_feature_set <- intersect(tested_drugs_in_smcs, tested_drugs_in_ecs)
    feature_set <- datbb %>%
      distinct(pr_gene_symbol) %>%
      .$pr_gene_symbol

    unique_perts_n <- length(unique(pert_names$pert_iname))
    message(qq("Unique total perturbations: @{unique_perts_n}"))
    unique_perts_shared_n <- length(unique(perturbation_feature_set))
    message(qq("Unique perturbations common to both HUVECs and HAoSMCs: @{unique_perts_shared_n}"))
    unique_features_n <- length(feature_set)
    message(qq("Total features: @{unique_features_n}"))


    message("Filtering datasets (with/witout replicates) to contain only those drugs also tested in vascular cell lines")
    message("Note: data still need to be filtered down to shared perts between datasets (e.g. P100 vs GCP)")
    no_replicates_dat <- datbb %>%
      filter(pert_iname %in% perturbation_feature_set) %>% # filter to drugs used in the
      left_join(cell_lines_identities, by = "cell_id") %>%
      select(cell_id, pert_iname, drug_class, pr_gene_symbol, pex_m, drug_cat, cell_cat, phosphosite, is_cancer)
    no_replicates_dat


    # # A tibble: 11,025 x 9
    # pr_gene_symbol   cell_id pert_iname    pex_m drug_class       drug_cat is_cancer pid     phosphosite
    # <chr>            <chr>   <chr>         <dbl> <chr>            <chr>    <lgl>     <chr>   <chr>
    #   1 H3.3K27me0K36me0 A375    ar a014418     0.11 Kinase inhibitor cancer   TRUE      BI10051 BI10051
    # 2 H3.3K27me0K36me0 A375    curcumin       0.05 Other            cancer   TRUE      BI10051 BI10051
    # 3 H3.3K27me0K36me0 A375    decitabine    -0.47 Epigenetic       cancer   TRUE      BI10051 BI10051
    # 4 H3.3K27me0K36me0 A375    dexamethasone  0.06 Other            cancer   TRUE      BI10051 BI10051
    # 5 H3.3K27me0K36me0 A375    dmso           0.11 control          cancer   TRUE      BI10051 BI10051
    # 6 H3.3K27me0K36me0 A375    geldanamycin   0.03 Epigenetic       cancer   TRUE      BI10051 BI10051
    # 7 H3.3K27me0K36me0 A375    gsk-j4        -0.43 Epigenetic       cancer   TRUE      BI10051 BI10051
    # 8 H3.3K27me0K36me0 A375    gsk126         0.99 Epigenetic       cancer   TRUE      BI10051 BI10051
    # 9 H3.3K27me0K36me0 A375    jq1-s         -0.02 Epigenetic       cancer   TRUE      BI10051 BI10051

    # give cell line ids to datb
    dat_with_reps <- datb %>%
      filter(pert_iname %in% perturbation_feature_set) %>%
      left_join(cell_lines_identities, by = "cell_id")

    # # A tibble: 37,828 x 13
    # pr_gene_symbol   well           pex cell_id pert_iname   rep   grp   pid    phosphosite cell_cat is_cancer drug_class    drug_cat
    # <chr>            <chr>        <dbl> <chr>   <chr>        <chr> <chr> <chr>  <chr>       <chr>    <lgl>     <chr>         <chr>
    #   1 H3.3K27me0K36me0 GA5-421DF-0…  0.11 A375    ar a014418   1     1     BI100… BI10051     cancer   TRUE      Kinase inhib… cancer
    # 2 H3.3K27me0K36me0 GA5-421DF-0…  0.13 A375    ar a014418   2     1     BI100… BI10051     cancer   TRUE      Kinase inhib… cancer
    # 3 H3.3K27me0K36me0 GA5-421DF-0… -0.09 A375    ar a014418   3     1     BI100… BI10051     cancer   TRUE      Kinase inhib… cancer
    # 4 H3.3K27me0K36me0 GA5-11373-0…  0.09 A375    curcumin     1     1     BI100… BI10051     cancer   TRUE      Other         cancer
    # 5 H3.3K27me0K36me0 GA5-11373-0…  0.01 A375    curcumin     2     1     BI100… BI10051     cancer   TRUE      Other         cancer
    # 6 H3.3K27me0K36me0 GA5-11373-0…  0.05 A375    curcumin     3     1     BI100… BI10051     cancer   TRUE      Other         cancer
    # 7 H3.3K27me0K36me0 GA5-30E17-0…  0.14 A375    decitabine   1     1,3   BI100… BI10051     cancer   TRUE      Epigenetic    cancer
    # 8 H3.3K27me0K36me0 GA5-30E17-0… -0.59 A375    decitabine   2     1,3   BI100… BI10051     cancer   TRUE      Epigenetic    cancer
    # 9 H3.3K27me0K36me0 GA5-30E17-0… -0.47 A375    decitabine   3     1,3   BI100… BI10051     cancer   TRUE      Epigenetic    cancer
    # 10 H3.3K27me0K36me0 GA5-11373-0…  0.23 A375    dexamethaso… 1     1     BI100… BI10051     cancer   TRUE      Other         cancer
    #
  } else {
    message(" !! Not !! filtering perturbations!")
    message("Note: data still need to be filtered down to shared perts between datasets (e.g. P100 vs GCP)")

    no_replicates_dat <- datbb %>%
      # filter(pert_iname %in% perturbation_feature_set) %>%
      left_join(cell_lines_identities, by = "cell_id") %>%
      select(cell_id, pert_iname, drug_class, pr_gene_symbol, pex_m, drug_cat, cell_cat, phosphosite, is_cancer)
    no_replicates_dat

    dat_with_reps <- datb %>%
      # filter(pert_iname %in% perturbation_feature_set) %>%
      left_join(cell_lines_identities, by = "cell_id")

    perturbation_feature_set <- datbb %>%
      distinct(pert_iname) %>%
      .$pert_iname
    feature_set <- datbb %>%
      distinct(pr_gene_symbol) %>%
      .$pr_gene_symbol
  }

  return(list(data = no_replicates_dat, data_with_reps = dat_with_reps, feature_set = feature_set, perturbation_feature_set = perturbation_feature_set))
}


#' @note merge data obj, created by *load_data*
merge_data <- function(data1, data2) {
  # DEBUG:
  # data1 = p100_obj_a
  # data2 = p100_obj_b

  merged_obj_lst <- map2(data1, data2, function(x, y) {
    # x <- data1[[1]]; y <- data2[[1]]
    if (typeof(x) == "list" & typeof(y) == "list") {
      r <- bind_rows(x %>% mutate(tag = "1st gen"), y %>% mutate(tag = "2nd gen")) %>% distinct()
      # duplicates between 1st and 2nd gen handled in get_replicates function
    } else {
      r <- unique(c(x, y))
    }
  })

  return(merged_obj_lst)
}


#' @Note convenience wrapper for loading the 'ready to analyze' data
#' @param datasets_directory base location of GCP and P100 1st gen data
#' @param references_directory location of *Drug Glossary_edited.xlsx*
load_wrapper <- function(datasets_directory = NULL, references_directory = NULL,
                         load_from_file = FALSE, filter_common_by_data_set = FALSE) {

  #' @note first, go to https://clue.io/cmapPy/build.html#install and install cmapPy to manipulate gct's
  #' @note then check out this tutorial from Oana: https://github.com/cmap/cmapPy/blob/master/tutorials/cmapPy_pandasGEXpress_tutorial.ipynb

  # DEBUG:
  # datasets_directory = DATASETS_DIRECTORY;
  # references_directory = REFERENCES_DIRECTORY;
  # load_from_file = FALSE;
  # filter_common_by_data_set = FALSE
  #
  output_fn <- file.path(datasets_directory, "ready-data.rds")
  if (load_from_file) {
    if (file.exists(output_fn)) {
      return(read_rds(output_fn))
    } else {
      message("Rds doesn't exist. Creating analysis-ready data object...")
    }
  }

  gcp_a <- load_gcp(first_gen_dataset = file.path(datasets_directory, "1st gen data", "GCP", "GCP All Cell Lines.gct"))
  gcp_obj <- load_data(
    data = gcp_a, drug_classes_fn = file.path(references_directory, "Drug Glossary_edited.xlsx"), dataset = "gcp",
    filter_by_vasc_perts = TRUE
  )

  p100a <- load_p100(first_gen_dataset = file.path(datasets_directory, "1st gen data", "P100", "P100 All Cell Lines.gct"))

  p100_obj <- load_data(
    data = p100a, drug_classes_fn = file.path(references_directory, "Drug Glossary_edited.xlsx"), dataset = "p100",
    filter_by_vasc_perts = TRUE
  )

  #' @TO-DO: you need to merge in the new HUVEC and HAoSMC data
  # p100b <- load_p100_alt(first_gen_dataset = file.path(datasets_directory, "2nd gen data", "P100", "P100-2nd-gen_reps-removed.gct"))

  # p100_obj_a <- load_data(data = p100a, drug_classes_fn = file.path(references_directory, "Drug Glossary_edited.xlsx"), dataset="p100",
  #                         filter_by_vasc_perts = TRUE)
  # p100_obj_b <- load_data(data = p100b, drug_classes_fn = file.path(references_directory, "Drug Glossary_edited.xlsx"), dataset="p100",
  #                         filter_by_vasc_perts = FALSE)
  #
  # perts_1st_gen <- p100_obj_a$data_with_reps$pert_iname %>% unique(); perts_1st_gen
  # analytes_1st_gen <- p100_obj_a$data_with_reps$pr_gene_symbol %>% unique(); analytes_1st_gen
  # redone_perts_for_huvecs_and_smcs_2nd_gen <- p100_obj_b$data_with_reps$pert_iname %>% unique(); redone_perts_for_huvecs_and_smcs_2nd_gen
  # redone_analytes_for_huvecs_and_smcs_2nd_gen <- p100_obj_b$data_with_reps$pr_gene_symbol %>% unique(); redone_analytes_for_huvecs_and_smcs_2nd_gen
  #
  # new_reps_a <- p100_obj_a$data_with_reps %>%
  #   filter(!(pert_iname %in% redone_perts_for_huvecs_and_smcs_2nd_gen) &
  #          !(cell_id %in% c("HUVEC", "HAoSMC")))
  #
  # p100_obj_a$data_with_reps <- new_reps_a
  #
  # p100_obj <- merge_data(data1 = p100_obj_a, data2 = p100_obj_b)


  # common perts and moas ACROSS datasets!!
  if (filter_common_by_data_set) {
    common_perts <- intersect(p100_obj$perturbation_feature_set, gcp_obj$perturbation_feature_set)
    common_moas <- intersect(p100_obj$data_with_reps$drug_class, gcp_obj$data_with_reps$drug_class)
    n_common_pert <- length(common_perts)
    n_common_moas <- length(common_moas)
  } else {
    common_perts <- NULL
    common_moas <- NULL
    n_common_pert <- "NA"
    n_common_moas <- "NA"
  }

  gcp_data_ready <- create_obj_lst(gcp_obj, dataset_common_perts = common_perts)
  p100_data_ready <- create_obj_lst(p100_obj, dataset_common_perts = common_perts)

  p100_key <- read_tsv("~/Downloads/master_tsv-p100-TRUE.tsv") %>%
    # bind_rows(read_tsv("~/Downloads/master_tsv-p100-FALSE.tsv")) %>%
    mutate(dataset = "p100")
  gcp_key <- read_tsv("~/Downloads/master_tsv-gcp-TRUE.tsv") %>%
    mutate(dataset = "gcp")

  master_key <- bind_rows(p100_key, gcp_key) %>%
    distinct()
  write_tsv(master_key, file.path(datasets_directory, "master-key.tsv"))

  message(qq("Total number of common drugs used: @{n_common_pert} "))
  message(qq("Total common MOAs among these drugs: @{n_common_moas} "))
  message(qq("Number of P100 perturbations: @{length(p100_obj$perturbation_feature_set)} "))
  message(qq("Number of GCP perturbations: @{length(gcp_obj$perturbation_feature_set)} "))

  lst_dat <- list(gcp = gcp_data_ready, p100 = p100_data_ready, master_key = master_key)
  message(qq("\nWriting to:\n @{output_fn} \nfor easy loading.\n"))
  write_rds(lst_dat, output_fn, compress = "gz")

  return(lst_dat)
}


#' @note wrapper to create the lst object for processing [perts], [moa], [all]
#' @param obj e.g. p100_obj, created from function *load_data*
#' @param dataset_common_perts the common perturbations between gcp and p100, calculated in *load_data*
create_obj_lst <- function(obj, dataset_common_perts = c("staurosporine")) {
  if (is.null(dataset_common_perts)) {
    message("!!! Not filtering perturbations between datasets !!!")
    dat1 <- obj$data
    dat2 <- obj$data_with_reps
  } else {
    message("Filtering perturbations by common drugs between datasets")
    dat1 <- obj$data %>%
      filter(pert_iname %in% dataset_common_perts)
    dat2 <- obj$data_with_reps %>%
      filter(pert_iname %in% dataset_common_perts)
  }

  dmap <- dat1 %>%
    distinct(pert_iname, drug_class)

  obj_ <- list(
    data = dat1,
    data_with_reps = dat2,
    feature_set = obj$feature_set,
    perturbation_feature_set = obj$perturbation_feature_set,
    drug_mapping = dmap
  )


  return(obj_)
}


#' @note create matrices ready to be analyzed from an object produced by [load_data]
#' @param dat_obj is a produced by load_data
#' @param grouping_var determines how the reps are collapsed - here we are interested in either each pert,
#' by drug_class, or by all of the perts together
#'    @NOTE if [grouping_var] is set to *all*, we REMOVE the grouping variable!!! so that we can combine everything
get_replicates <- function(dat_obj, grouping_var = "pert_iname") {

  # perturbations <- dat_obj$perturbation_feature_set; perturbations
  features <- dat_obj$feature_set
  features
  drug_classes <- dat_obj$drug_mapping
  drug_classes

  group_var_sym <- sym(grouping_var) # use symbol to programmatically evaluate

  mat_reps_temp <- dat_obj$data_with_reps %>% # using the replicates data
    select(pr_gene_symbol, well, cell_id, pert_iname, drug_class, pex)


  # data_by_group

  # grouping_var = "pert_iname"; grouping_var = "drug_class"; grouping_var = "all"
  if (tolower(grouping_var) == "pert_iname") {
    print(grouping_var)
    mat_reps_temp_a <- mat_reps_temp %>%
      group_by(!!group_var_sym) %>% # using !! here to evaluate group_var_sym to get the grouping variable
      nest(data_by_tx_group = c(pr_gene_symbol, well, cell_id, pex))

    mat_reps_temp_b <- mat_reps_temp_a %>%
      # if not 'all', then data2 needs to include pert_iname and drug class so that this info can be propogated to df at the end
      mutate(data_by_cell = pmap(.l = list(data_by_tx_group, pert_iname, drug_class), function(x, y, z) {
        # x <- mat_reps_temp$data[[1]]; y= mat_reps_temp$pert_iname[[1]]; z = mat_reps_temp$drug_class[[1]]
        pc <- x %>%
          mutate(
            pert_iname = y, # give drug class and pert_iname a column outside of the nested datasets
            drug_class = z
          ) %>%
          group_by(cell_id, drug_class, pert_iname) %>%
          nest(data_per_analyte = c(pr_gene_symbol, well, pex)) # each replicates data [longform of gene, pex], for each cell and pert
        return(pc)
      })) %>%
      mutate(org = grouping_var) %>%
      select(org, pert_iname, drug_class, data_by_tx_group, data_by_cell)

    # this inspects dmso for
    # mat_reps_temp_b$data_by_cell[[4]] # notice that there are 3x wells except for huvecs and haosmcs
    # table(mat_reps_temp_b$data_by_cell[[4]]$data_per_analyte[[1]]$well) # inspect one
  } else if (tolower(grouping_var) == "drug_class") {
    print(grouping_var)
    mat_reps_temp_a <- mat_reps_temp %>% # using !! here to evaluate group_var_sym to get the grouping variable
      group_by(!!group_var_sym) %>%
      nest(data_by_tx_group = c(pr_gene_symbol, well, pert_iname, cell_id, pex))

    mat_reps_temp_b <- mat_reps_temp_a %>%
      # x <- mat_reps_temp_a$data[[2]]; y <- mat_reps_temp_a$drug_class[[2]]
      mutate(data_by_cell = map2(data_by_tx_group, drug_class, function(x, y) {
        pc <- x %>%
          mutate(drug_class = y) %>% # give drug class a column outside of the nested datasets
          group_by(cell_id, drug_class) %>%
          nest(data_per_analyte = c(pr_gene_symbol, well, pex)) # each replicates data [longform of gene, pex], for each cell and pert
        return(pc)
      })) %>%
      mutate(org = grouping_var) %>%
      select(org, drug_class, data_by_tx_group, data_by_cell)
  } else {
    print(grouping_var)
    mat_reps_temp_a <- mat_reps_temp %>% # using !! here to evaluate group_var_sym to get the grouping variable
      nest(data_by_tx_group = c(pr_gene_symbol, well, pert_iname, drug_class, cell_id, pex))

    mat_reps_temp_b <- mat_reps_temp_a %>%
      # x <- mat_reps$data[[1]]
      mutate(data_by_cell = map(data_by_tx_group, function(x) {
        pc <- x %>%
          group_by(cell_id) %>% # no columns needed outside nested dfs, other than org
          nest(data_per_analyte = c(pr_gene_symbol, well, pex)) # each replicates data [longform of gene, pex], for each cell and pert
        return(pc)
      })) %>%
      mutate(org = grouping_var) %>%
      select(org, data_by_tx_group, data_by_cell)
  }

  mat_reps <- mat_reps_temp_b %>%
    # cell_dat_df <- mat_reps_temp_b$data_by_cell[[1]]
    mutate(df = map(data_by_cell, function(cell_dat_df) {

      # print(unique(cell_dat_df$drug_class))
      cell_dat <- cell_dat_df$data_per_analyte
      cell_dat
      nrows_df <- cell_dat_df %>%
        mutate(n_rows = map_int(data_per_analyte, function(x) nrow(x)))
      nrows_df
      # num_perturbations <- length(perturbations); num_perturbations

      master_cell_dat <- nrows_df %>%
        mutate(cell_ids = map2(data_per_analyte, cell_id, function(d, c) {
          cell_ids <- rep(c, nrow(d %>% distinct(well)))
        })) %>%
        mutate(pert_names = map2(data_per_analyte, pert_iname, function(d, p) {
          pert_names <- rep(p, nrow(d %>% distinct(well)))
        })) %>%
        mutate(d_classes = map2(data_per_analyte, drug_class, function(d, dc) {
          d_classes <- rep(dc, nrow(d %>% distinct(well)))
        }))
      master_cell_dat
      # basically shows that some cells were treated with a drug more than 3 times (some 9 times)
      # check this to debug below, should depend on the number of reps!

      # n_series_vec <- master_cell_dat$num_series
      # cell_vec <- master_cell_dat$cell_id
      # pert_vec <- master_cell_dat$pert_iname
      # moa_vec <- master_cell_dat$drug_class

      # repeat these names number of times defined by num_series
      # cell_ids <- rep(cell_vec, n_series_vec)
      # pert_names <- rep(pert_vec, n_series_vec)
      # d_classes <- rep(moa_vec,n_series_vec)
      cell_ids <- master_cell_dat$cell_ids %>% unlist()
      pert_names <- master_cell_dat$pert_names %>% unlist()
      d_classes <- master_cell_dat$d_classes %>% unlist()

      # go through each drug for this cell line and spread to get it into cell_id/well x gene/analyte with value pex
      well_by_analyte_df_plus_annots <- map_df(cell_dat, function(analyte_df) {
        # analyte_df <- cell_dat[[3]]
        # gives a well x analyte matrix
        wide_analyte_df <- analyte_df %>%
          group_by(well, pr_gene_symbol) %>% # if there are dup wells, compress by median
          summarize(pex = median(pex, na.rm = T), .groups = "keep") %>%
          spread(pr_gene_symbol, pex)
        wide_analyte_df
        return(wide_analyte_df)
      }) %>%
        ungroup() %>%
        mutate(
          cell_id = cell_ids,
          pert_iname = pert_names,
          drug_class = d_classes
        ) %>%
        select(well, cell_id, pert_iname, drug_class, everything())
      # now you have a nice matrix of well x analyte

      return(well_by_analyte_df_plus_annots)
    }))

  final_mat_reps <- mat_reps %>%
    select(org, everything(), -data_by_tx_group, -data_by_cell) %>% # new organization variable, 'org'
    mutate(
      dataframe = map(df, .f = convert_to_dataframe, keep_cols_vec = features),
      t_dataframe = map(dataframe, function(x) t(x))
    )

  # pre for debugging, post for downstream analysis
  return(list(mat_reps, final_mat_reps) %>% set_names("pre", "post"))
}

#' @note match dfs
#' @param cur_mat the matrix you would like to replace the rows and column names of
#' @param master_df made by reading and combining files written function *load_data*
#' @param from name of column in master_df you'd like to convert from
#' @param to name of column in master_df you'd liek to convert to
#' @param rc takes "r" or "c" or both as a character vector, whether you'd like to replace the columns or rows
match_function <- function(cur_mat, master_df = NULL,
                           from = "well", to = "cell_id", rc = c("r", "c")) {
  which_axis <- which(rc %in% c("r", "c"))
  ordered_names_df <- tibble(
    rnames = rownames(cur_mat),
    cnames = colnames(cur_mat)
  ) %>%
    select(all_of(which_axis))

  # rename(!!new_col_name_[1] := rnames)
  match_df_temp <- master_df %>%
    select(from, to) %>%
    distinct()

  from_idx <- which(colnames(match_df_temp) %in% from)
  to_idx <- which(colnames(match_df_temp) %in% to)

  stopifnot(length(unique(match_df_temp[from_idx])) == length(unique(match_df_temp[to_idx])))

  match_col_name_from <- eval(from)
  match_col_name_from_ <- str_c(match_col_name_from, 1:length(which_axis))

  match_col_name_to <- eval(to)
  match_col_name_to_ <- str_c(match_col_name_to, 1:length(which_axis))

  match_df <- map(1:length(which_axis), function(i) {
    i_names_df_match_from <- ordered_names_df %>% rename(!!match_col_name_from_[i] := !!colnames(ordered_names_df)[i])
    mdf_match_from <- match_df_temp %>% rename(!!match_col_name_from_[i] := !!colnames(match_df_temp)[from_idx])

    f_df <- suppressMessages(left_join(i_names_df_match_from, mdf_match_from) %>%
      rename(!!match_col_name_to_[i] := !!colnames(.)[length(colnames(.))]))
  }) %>%
    bind_cols() %>%
    select(all_of(match_col_name_from_), all_of(match_col_name_to_))

  match_df[match_col_name_to_[1]] <- make.unique(match_df[, match_col_name_to_[1]] %>% pull())
  match_df[match_col_name_to_[2]] <- make.unique(match_df[, match_col_name_to_[1]] %>% pull())

  return(match_df)
}


#' @note this function replaces row and/or column names with a target vector, from and to must be exchangeable!!
#' @param cur_mat the matrix you would like to replace the rows and column names of
#' @param match_df unique well, cell_id, pert_iname, drug_class tibble
#' @param new_match_id col name of match_df that you'd like to replace the rows and cols of cur_mat with
replace_rc_names <- function(cur_mat, match_df, rc = c("r", "c"), new_match_id = "cell_id") {
  # cur_mat <- corr_mat_reps$corr_r[[1]];
  # match_df = master_key; rc <- c("r", "c")

  new_mat <- cur_mat
  new_mat

  if ("r" %in% rc) {
    # match up by rows
    convert_rows <- tibble(old = rownames(cur_mat)) %>%
      left_join(match_df, by = c("old" = "well")) %>%
      rename("well" = "old")
    rownames(new_mat) <- convert_rows[new_match_id] %>% pull()
  } else {
    convert_rows <- NULL
    rownames(new_mat) <- rownames(cur_mat)
  }

  if ("c" %in% rc) {
    # match up by colnames
    convert_cols <- tibble(old = colnames(cur_mat)) %>%
      left_join(match_df, by = c("old" = "well")) %>%
      rename("well" = "old")
    colnames(new_mat) <- convert_cols[new_match_id] %>% pull()
  } else {
    convert_cols <- NULL
    colnames(new_mat) <- colnames(cur_mat)
  }

  return(list("new_mat" = new_mat, "row_match" = convert_rows, "col_match" = convert_cols))
}


#' @note convert a tibble to a dataframe
#' @param x assumes that x is a tibble where the first column will become the row names of the new data.frame
convert_to_dataframe <- function(x, keep_cols_vec) {
  x_df <- as.data.frame(x)
  rownames(x_df) <- x[, 1, drop = T]
  x_df_c <- x_df %>% select(all_of(keep_cols_vec))
  return(x_df_c)
}

#' @note helper function for compute_boot_pvclust and morpheus plotter
#' @param long_df to be converted to matrix -- [assumes long df is a connectivity df!]
make_numeric_mat <- function(long_df) {
  if ("p_val" %in% colnames(long_df)) {
    long_df <- long_df %>% select(-p_val)
  }
  x <- long_df %>%
    spread(group_b, conn) %>%
    as.data.frame() %>%
    select(-group_a) %>%
    as.matrix()
  rownames(x) <- colnames(x)
  return(x)
}

#' @note create cluster assignments df for plotting
#' @param ca cluster assignments NAMED array!! this is the "cut_trees" field in diff_ex
create_ca_df <- function(ca) {
  t <- tibble::enframe(ca, name = "cell_id", value = "cluster") %>%
    distinct(cell_id, cluster)
  return(t)
}

#' @note flattens a correlation matrix
#' @param m a correlation matrix
flattenMatrix <- function(m) {
  ut <- upper.tri(m)
  lt <- lower.tri(m)
  bind_rows(
    data.frame(
      i = rownames(m)[row(m)[ut]],
      j = rownames(m)[col(m)[ut]],
      cor = t(m)[ut],
      stringsAsFactors = F
    ),
    data.frame(
      i = rownames(m)[row(m)[lt]],
      j = rownames(m)[col(m)[lt]],
      cor = t(m)[lt],
      stringsAsFactors = F
    )
  )
}


#' @note Remove columns/rows/both that have an NA percentage greater than or equal to the specified threshold
#' @param dat matrix
#' @param thresh_row percent of NA values along rows that you want to remove
#' @param thresh_col percent of NA values along cols that you want to remove
#' @param d rows and columns? possible values include c("both", "cols", "rows")
remove_x_perc_NA <- function(dat, thresh_row = 1.0, thresh_col = 1.0, d) {
  if (as.logical(match(d, "both"))) {
    # message("both")
    rsd <- which(rowMeans(is.na(dat)) >= thresh_row) # message("Removing: ", length(rsd), " rows\n")
    rs <- setdiff(1:nrow(dat), rsd)
    csd <- which(colMeans(is.na(dat)) >= thresh_col)
    cs_n <- paste(names(csd), sep = "", collapse = " ") # message(qq("Removing cols: @{cs_n}"))
    cs <- setdiff(1:ncol(dat), csd)
    new <- as.matrix(dat[rs, cs, drop = F])
    return(list(mat = new, removed = list(cols = names(csd), rows = rsd)))
  } else if (as.logical(match(d, "cols"))) {
    # message("cols")
    csd <- which(colMeans(is.na(dat)) >= thresh_col)
    cs_n <- paste(names(csd), sep = " ", collapse = "") # message(qq("Removing cols: @{cs_n}"))
    cs <- setdiff(1:ncol(dat), csd)
    new <- as.matrix(dat[, cs, drop = F])
    return(list(mat = new, removed = list(cols = names(csd), rows = NA)))
  } else {
    # message("rows")
    rsd <- which(rowMeans(is.na(dat)) >= thresh_row) # message("Removing: ", length(rsd), " rows")
    rs <- setdiff(1:nrow(dat), rsd)
    new <- as.matrix(dat[rs, , drop = F])
    return(list(mat = new, removed = list(cols = NA, rows = rsd)))
  }
}


#' @note compute the average of two connectivity matrices and output a single matrix
#' @param conn_mat1 p100 or gcp connectivity matrix, could also be an MOA connectivity matrix
#' @param conn_mat2 p100 or gcp connectivity matrix
compute_avg_conn_mat <- function(conn_mat1, conn_mat2) {
  # test <- left_join(gcp_conn_clust %>% filter(pert_iname %in% common_perts), p100_conn_clust %>% filter(pert_iname %in% common_perts), by="pert_iname", suffix=c("_gcp","_p100"))
  # conn_mat1 <- test$med_conn_mat_p100[[4]]; conn_mat2 <- test$med_conn_mat_gcp[[4]]
  all_cols <- union(colnames(conn_mat1), colnames(conn_mat2))
  cols <- intersect(colnames(conn_mat1), colnames(conn_mat2))
  msg <- str_c(setdiff(all_cols, cols), collapse = ",")
  qq("Dropping @{msg}...")
  conn_mat1_c <- conn_mat1[cols, cols, drop = F]
  conn_mat2_c <- conn_mat2[cols, cols, drop = F]
  X <- list(conn_mat1_c, conn_mat2_c)
  v <- reduce(X, `+`) / length(X)
  return(v)
}


#' @note combine data across analytes for a replicate; assumes df is in [replicate x analyte] form
#' @param df1 first df to combine
#' @param df2 second df to combine
make_avg_df <- function(df1, df2, vars_to_join = c("well", "cell_id", "pert_iname", "drug_class")) {

  # message(qq("Joining by vars: \n@{str_c(vars_to_join,collapse=', ')}"))
  # we must create ID system by which to compare cell-drug replicates with each other! then inner join on this
  # otherwise, the two systems are not comparable

  # does it matter how we decide to combine them?? this is probably wrong
  temp1 <- df1 %>%
    mutate(
      twell = "avg",
      idx = make.unique(str_c(pert_iname, drug_class, sep = "-"))
    )

  temp2 <- df2 %>%
    mutate(
      twell = "avg",
      idx = make.unique(str_c(pert_iname, drug_class, sep = "-"))
    )

  # mutate(twell = "avg",
  #        idx = make.unique(str_c(pert_iname, drug_class, order_names[2], sep = "-")))
  #
  join1 <- temp1 %>%
    mutate(well = temp1$idx) %>%
    select(-twell, -idx)
  join2 <- temp2 %>%
    mutate(well = temp2$idx) %>%
    select(-twell, -idx)

  # we must bind by drug-cell combinations... this is probably wrong
  df <- inner_join(join1, join2, by = c(vars_to_join))

  return(df)
}



#' @note function to plot pvclust and diffex by pert
#' @param obj returned by batch_exec, must contain column ca_df
plot_pvclust_and_diffex <- function(obj, dend_thresh = 0.6, bh_thresh = 0.1) {
  pvclust_res <- obj$pvclust[[1]]
  diff_ex_res <- obj$diff_ex[[1]]
  dname <- obj$pert_iname
  ca_df <- obj$ca_df[[1]]

  p1 <- suppressWarnings(suppressMessages(display_pvclust(pvclust_res, dname = dname, thresh = dend_thresh)))
  p2 <- suppressWarnings(suppressMessages(plot_diff_analyte_results(diffe_by_clust = diff_ex_res, cluster_assignments_df = ca_df, dname = dname, BH_THRESH = bh_thresh)))
  print(p2)
}
