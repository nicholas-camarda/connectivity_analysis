
winos <- ifelse(grepl("windows", Sys.info()["sysname"], ignore.case = T), 1, 0)
if (winos == 1) {
  source("C:\\Users\\ncama\\OneDrive - Tufts\\phd\\ws\\scripts\\master-source.R")
} else {
  source("/Users/Nicholas/OneDrive - Tufts/phd/ws/proteomics/scripts/master-source.R")
}

my_extract <- function(chr, idx) {
  res <- str_split(string = chr, pattern = "-", simplify = T)[, idx]
  return(res)
}
# NOTE: read all the data from RAW GCT directory and
#' merge everything using library functions, not excel

cat("\nReading and merging data...")
parent_dir_fn <- file.path(DATASETS_DIRECTORY, "RAW GCT")
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

anonymous_gct_merge <- function(dtf1, dtf2) {
  res <- merge_gct(dtf1, dtf2, dim = "column", matrix_only = FALSE)
  return(res)
}

#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
dataset <- "P100"
grouping_var <- "pert_iname"
SPECIFIC_OUTPUT_DIR <- file.path(OUTPUT_DIRECTORY, dataset, grouping_var)

p100_data <- data %>% filter(dataset_type == dataset)
p100_data_lst <- as.list(p100_data$obj) %>% setNames(p100_data$fn)
# https://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list
#' recursively merge list
merged_p100_obj <- suppressMessages(Reduce(
  f = anonymous_gct_merge,
  x = p100_data_lst
))

# (p100_se <- as(merged_p100_obj, "SummarizedExperiment"))

# my_mat <- mat(merged_p100_obj) %>% as.data.frame() %>%
#   rownames_to_column("id") my_row_meta <- meta(merged_p100_obj, dimension =
#   "row") joined <- left_join(my_mat, my_row_meta) %>% as_tibble() %>%
#   select(pr_gene_symbol, `PA5-1B07E-001A01`:`PYC-35F57-196H12`) rbm17 <-
#   joined %>% filter(pr_gene_symbol == "RBM17")


# my meta data
drug_classes_fn <- file.path(REFERENCES_DIRECTORY, "Drug Glossary_edited.xlsx")
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

mat_ <- mat(merged_p100_obj)


# melt the GCT into a long a table (thank god for this function)
melted_merged_p100_obj <- melt_gct(merged_p100_obj) %>%
  as_tibble()


p100_obj <- melted_merged_p100_obj %>%
  rename(
    row_id = id.x,
    column_id = id.y
  ) %>%
  mutate(pert_iname = tolower(pert_iname)) %>%
  ## inner join to only do drugs that I pick!!
  inner_join(drugs_moa_df, by = "pert_iname") %>%
  mutate(master_id = str_c(cell_id, pert_iname, pert_class, sep = "--")) %>%
  mutate(replicate_id = str_c(master_id, column_id, sep = "::")) %>%
  select(master_id, replicate_id, everything()) %>%
  group_by(master_id, replicate_id) %>%
  # some gene names are the same, so you need to make them unique and reference
  # the phosphosite later
  mutate(u_pr_gene_symbol = make.unique(pr_gene_symbol, sep = "_"))


full_splt_lst <- split(p100_obj, f = p100_obj[[sym(grouping_var)]])
length(full_splt_lst)

cat("\nComputing correlation...")
sample_corr_lst <- map(full_splt_lst, function(x) {
  # x <- full_splt_lst[[1]]

  wide_x <- x %>%
    ungroup() %>%
    select(master_id, replicate_id, u_pr_gene_symbol, value) %>%
    pivot_wider(names_from = u_pr_gene_symbol, values_from = value)

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
    pivot_longer(cols = colnames(res_final1)[-1], 
    names_to = "unique_id_b", values_to = "corr") %>%
    right_join(look_up_df %>%
      rename(group_a = master_id, unique_id_a = replicate_id), 
      by = "unique_id_a") %>%
    right_join(look_up_df %>%
      rename(group_b = master_id, unique_id_b = replicate_id), 
      by = "unique_id_b")

  res_final2
  cat(".")
  return(list("matrix" = res, "tibble" = res_final2))
})
cat(".Done.\n")


# connectivity is computed across replicates!!
run_conn_clust <- function(corr_lst) {
  # NOTE: 
  suppressMessages(library(tidyverse))
  suppressMessages(library(Matching))
  suppressMessages(library(pvclust))
  suppressWarnings(suppressMessages(source("./scripts/master-source.R")))
  # corr_lst <- sample_corr_lst$`1271738-62-5` 

  corr_tbl <- corr_lst$tibble
  groups <- corr_tbl %>%
    distinct(group_a) %>%
    .$group_a

  # is this a multidyplr moment?
  res <- map_df(groups, function(g1) {
    df2 <- map_df(groups, function(g2) {
      # g1 <- groups[1]; g2 <- groups[2]
      #' this is by definition all of the comparisons between A and B
      not_ <- groups[!(groups %in% c(g1, g2))]
      test <- corr_tbl %>%
        filter(group_a == g1 & group_b == g2) %>%
        pluck("corr")
      background <- corr_tbl %>%
        filter(group_a == g2, group_b %in% not_) %>%
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
    return(df2)
  })
  cat(".")

  res_coded_for_median <- res %>%
    mutate(group_a = as.character(group_a),
           group_b = as.character(group_b)) %>%
    mutate(grouping_var_code = map2_chr(group_a, group_b, function(a,b) {
      vec_a <- str_split(string = a, pattern = "--", simplify = T)
      a_code <- vec_a[1] # c(1, length(vec_a))
      vec_b <- str_split(string = b, pattern = "--", simplify = T)
      b_code <- vec_b[1] # c(1, length(vec_b))
      # grouping by this code should create pairs
      coded <- str_c(sort(c(a_code, b_code)), collapse =  "_")
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
    mutate_at(vars(-group_cols()), ~replace(., is.na(.), 1)) %>%
    as.data.frame() %>%
    dplyr::select(-group_a) %>%
    as.matrix()
  rownames(res_median) <- colnames(res_median)

  cell_names_from_long <- str_split(
    string = rownames(res_median),
    pattern = "--", simplify = T

  )[, 1]
  pretty_name_res_median <- res_median
  rownames(pretty_name_res_median) <- cell_names_from_long
  colnames(pretty_name_res_median) <- cell_names_from_long

  clust_obj <- compute_boot_pvclust(
    x = pretty_name_res_median,
    parallel_flag = FALSE,
    n_boot = 1000
  )
  clust_assignments <- generate_clusters(
    x = clust_obj,
    thresh = DENDRO_CUT_THRESH
  )
  co_clust_bool <- co_cluster(cut_tree_obj = clust_assignments)
  clustering_output_lst <- list(clust_obj, clust_assignments, co_clust_bool) %>%
    setNames(c("clust_obj", "clust_assignments", "co_clust_bool"))

 # complete, and collapsed by median
 return(list(
   "res_tibble" = res_coded_for_median,
   "res_median_matrix" = res_median,
   "clust_obj" = clust_obj,
   "clust_assignments" = clust_assignments,
   "co_clust_bool" = co_clust_bool
 ))
}


message("\nComputing connectivity...")
clock1 <- proc.time()
param <- SnowParam(
  workers = registered()$SnowParam$.clusterargs$spec,
  type = "SOCK"
)
sample_conn_lst <- bplapply(
  X = sample_corr_lst,
  FUN = run_conn_clust,
  BPPARAM = param
)

clock2 <- proc.time() - clock1
message("Done")
message(qq("Completed processing @{length(sample_conn_lst)} items"))
message(qq("Total time: in @{unname(clock2[3])} seconds"))

connctivity_output_dir <- file.path(SPECIFIC_OUTPUT_DIR, "conn")
dir.create(connctivity_output_dir, showWarnings = FALSE, recursive = T)
message(qq("Writing output to path:\n@{connctivity_output_dir}"))
write_rds(
  x = sample_conn_lst,
  file = file.path(connctivity_output_dir, "conn.rds"),
  compress = "gz"
)

# TODO: diff ex, heatmaps dendros


