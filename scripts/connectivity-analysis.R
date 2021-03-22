
winos <- ifelse(grepl("windows", Sys.info()["sysname"], ignore.case = T), 1, 0)
if (winos == 1) {
  source("C:\\Users\\ncama\\OneDrive - Tufts\\phd\\ws\\scripts\\master-source.R")
} else {
  source("/Users/Nicholas/OneDrive - Tufts/phd/ws/proteomics/scripts/master-source.R")
}

cat("\nReading and merging data...")
#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
dataset <- "P100"
grouping_var <- "pert_iname"
SPECIFIC_OUTPUT_DIR <- file.path(OUTPUT_DIRECTORY, dataset, grouping_var)

merged_p100_obj <- parse_gctx(file.path(
  DATASETS_DIRECTORY, "Inherited Data", "1st gen data",
  "P100", "P100-All-Cell-Lines.gct"
))

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
#' @note connectiivty and clustering super function
#' @param corr_lst sample-level correlation values
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

#' @note run parallel connectivity and clustering analysis
#' @param sample_corr_lst sample correlation list obj
run_conn_clust_par <- function(sample_corr_lst) {
  clock1 <- proc.time()
  param <- SnowParam(
    # number of cores
    workers = registered()$SnowParam$.clusterargs$spec,
    # socket cluster can be used with windows backend
    type = "SOCK"
  )
  # dmso <- compute_connectivity(M = sample_corr_lst$dmso$matrix)
  # dmso_median <- collapse_connectivity_by_median(dmso %>%
  #   rename(group_a = grp_namesA, group_b = grp_namesB))
  # dmso_clust_obj <- compute_boot_pvclust(dmso_median)
  # dmso_clust_assign <- generate_clusters(x = dmso_clust_obj,
  # thresh = DENDRO_CUT_THRESH)
  #  co_clust_bool <- co_cluster(cut_tree_obj = dmso_clust_assign)

  sample_conn_lst <- bplapply(
    X = sample_corr_lst,
    FUN = run_conn_clust,
    BPPARAM = param
  )

  clock2 <- proc.time() - clock1
  message("Done")
  return(sample_conn_lst)
}

connctivity_output_dir <- file.path(SPECIFIC_OUTPUT_DIR, "cache")
dir.create(connctivity_output_dir, showWarnings = FALSE, recursive = TRUE)
res_cache_obj_fn <- file.path(connctivity_output_dir, "data.rds")


response <- readline(prompt = "Overwrite [y/n] (yes or load from cache): ")
if (response == "y") {
  message("\nComputing connectivity and performing clustering...")
  sample_conn_lst <- run_conn_clust_par(sample_corr_lst)
} else if (response == "n") {
  if (!file.exists(res_cache_obj_fn)) {
    message("Cache obj doesn't exist. Computing conn/clust...")
    message("\nComputing connectivity and performing clustering...")
    sample_conn_lst <- run_conn_clust_par(sample_corr_lst)
  }
  else {
    message("Loading from cache...")
    sample_conn_lst <- read_rds(file.path(connctivity_output_dir, "data.rds"))
  }
} else {
  message("Incorrect input. Will attempt to load from cache.")
  sample_conn_lst <- tryCatch({
    read_rds(file.path(connctivity_output_dir, "data.rds"))
    message("Successfully loaded from cache.")
  },
    error = function(cond) {
      message("Loading failed. Running conn/clust analysis...")
      sample_conn_lst <- run_conn_clust_par(sample_corr_lst)
    }
  )
}
message("Done.")

message(qq("Completed processing @{length(sample_conn_lst)} items"))
message(qq("Total time: in @{unname(clock2[3])} seconds"))


message(qq("Writing output to path:\n@{res_cache_obj_fn}"))

multi_obj_save <- list(
  "p100_obj" = p100_obj,
  "corr_obj" = sample_corr_lst,
  "conn_obj" = sample_conn_lst
)

write_rds(
  x = multi_obj_save,
  file = res_cache_obj_fn,
  compress = "gz"
)

# DEBUG:
# sample_conn_lst <- read_rds(file.path(connctivity_output_dir, "data.rds"))

co_clust_lst <- map(sample_conn_lst, function(lst) {
  if (is.na(lst$co_clust_bool)) {
    return(NULL)
  }
  if (lst$co_clust_bool == TRUE) {
    return(lst)
  } else { 
    return(NULL)
  }
})


# TODO: diff ex, heatmaps dendros


