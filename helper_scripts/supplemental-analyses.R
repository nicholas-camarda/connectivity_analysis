library(patchwork)
library(umap)
library(uwot)
library(ggbiplot)
library(factoextra)

###' These are *directory* and *input data* related variables
#' name of *output directory*, not the full path (this will be appended automatically)
my_output_directory <- "final" # "test-LVL4" # All-LINCS-data-LVL4

#' name of *data directory*, should contain "-LVL3|4" info so that script knows what kind of data we are using
specific_data_directory <- "All-LINCS-data-LVL4"

# the name of the input file containing analysis run structure; i.e.
# which grouping structure (drug name or class), which drugs / classes to filter in, 
# which dataset, and whether to exclude any perturbations specifically from the analysis
args_fn_name <- "all_args.csv" 


###' These are *result filtering* related variables
# these are the vascular cell types
vascular_char_vec <- c("HUVEC", "HAoSMC", "Pericyte")
# whether to filter out analytes that aren't well represented across cell types and perturbations
filter_analytes <- FALSE

# this variable is only important if filter_analytes = TRUE
# this percentage denotes the amount of usable data tolerated; 
# i.e. if this is 0.75, that means an analyte may not have more than 25% missing data
percent_analyte_na_to_throwout <- 0.75

# cut trees at x% the maximum height of the tree
dendro_cut_thresh <- 0.6

data_directory <- file.path("data")
datasets_directory <- file.path(data_directory, "datasets")
my_output_directory <- "final" # "test-LVL4" # All-LINCS-data-LVL4

# filter out analytes that are strictly above this threshold
bh_thresh_val <- 0.1

source(file.path("scripts", "load.R"))

output_directory_temp <- file.path(str_c("output", my_output_directory, sep = "_"))
no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)
# analysis dat loaded into memory with this source call and a bunch of functions

supplemental_output_dir <- file.path(output_directory, "supplemental_output")
dir.create(supplemental_output_dir, showWarnings = F, recursive = T)

informative_lbls <- create_informative_labels()$informative_labels_df 
pert_color_map_init <- generate_color_palette()

# SWITCH THIS OFF FOR NO DEBUG
# analysis_dat_temp <- analysis_dat %>% slice(1)


conn_results <- apply(analysis_dat, MARGIN = 1, function(args){
  # args <- analysis_dat[1,]
  d <- args$data
  # d <- args$data[[1]]
  grouping_id <- args$grouping_var
  filter_id <- args$filter_vars
  # filter_id <- args$filter_vars[[1]]
  
  print(grouping_id)
  print(filter_id)
  # stop()
  dat <- d %>%
    filter(!!sym(grouping_id) == filter_id) %>%
    filter(pert_iname %in% my_perts)
  
  dataset_type <- unique(dat$which_dat)
  
  # dat <- dat %>% filter(cell_id %in% c("MCF7", "HAoSMC", "HUVEC"))
  
  # now do this with connectivity values
  cor_dat <- run_corr(dat)$matrix
  unique_cs <- unique(rownames(cor_dat))
  grp_names <- get_grp_names_from_matrix(cor_dat, unique_cs, sep_ = "--")
  
  cs_idx <- seq_len(length(rownames(cor_dat)))
  # collapse replicates by drug - cell pair, not just cell
  cs_tib <- tibble(unique_cs, cs_idx, grp_names)
  # %>%
  #   mutate(cell_ = str_split(unique_cs, "--", simplify = T)[,1],
  #          pert_ = str_split(unique_cs, "--", simplify = T)[,2],
  #          g_cell_pert_ = str_c(cell_, pert_, sep = "--"))
  # 
  # to do cell-pert replicates
  # looper <- expand.grid(
  #   rep_a = unique(cs_tib$g_cell_pert_), 
  #   rep_b = unique(cs_tib$g_cell_pert_)
  # ) %>%
  #   as_tibble()
  
  # to do it like i did in the original analysis
  looper <- expand.grid(
    group_a = unique(grp_names), 
    group_b = unique(grp_names)
  ) %>%
    as_tibble()
  
  # connectivity of each cell drug (all 3 replicates) vs all others
  
  with_progress(expr = {
    conn_progress <- progressr::progressor(steps = nrow(looper))
    conn_df <- map_dfr(1:nrow(looper), .f = function(z){
      # message(looper$rep_a[z], "  ", looper$rep_b[z])
      res <- calc_conn(a = looper$group_a[z], b = looper$group_b[z],
                                  my_mat = cor_dat, my_tib = cs_tib,
                                  use_bootstrap = F)
      conn_progress()
      return(res)
    }) %>%
      bind_rows() 
    gc(verbose = F)
  })

  replicate_matrix <- conn_df %>%
    dplyr::select(group_a, group_b, conn) %>%
    pivot_wider(id_cols = group_a,
                names_from = group_b,
                values_from = conn,
                values_fn = function(x) median(x, na.rm = TRUE))
  
  numeric_replicate_matrix <- replicate_matrix %>%
    column_to_rownames("group_a") %>%
    dplyr::select(where(is.numeric))
  
  # scale and 0-center the data
  scaled_numeric_replicate_matrix <- scale(numeric_replicate_matrix)
  
  final_result <- tibble(grouping_id = grouping_id, filter_id = filter_id, dataset_type = dataset_type,
                         dat = list(dat),
                         conn_df = list(conn_df),
                         replicate_matrix = list(replicate_matrix),
                         numeric_replicate_matrix = list(numeric_replicate_matrix),
                         scaled_numeric_replicate_matrix = list(scaled_numeric_replicate_matrix))
  return(final_result)
}) %>%
  bind_rows()

write_rds(conn_results, file = file.path(supplemental_output_dir, "conn_results-cell_cell.rds"))
conn_res <- read_rds(file.path(supplemental_output_dir, "conn_results-cell_cell.rds"))

n_neighbors_tib <- tribble(~dataset_type, ~filter_id, ~n_neighbors,
                           "P100", "Kinase inhibitor", 4,
                           "P100", "DMSO", 4,
                           "GCP", "Epigenetic", 4,
                           "GCP", "DMSO", 4)
custom.config <- umap.defaults
updated_conn_res <- left_join(conn_res, n_neighbors_tib, by = c("filter_id", "dataset_type"))

umap_uwot_pca_results <- apply(updated_conn_res, MARGIN = 1, FUN = function(args){
  n_neighbors_uwot <- args$n_neighbors
  custom.config$n_neighbors <- n_neighbors_uwot
  
  dat <- args$dat # dat <- args$dat[[1]]
  grouping_id <- unique(args$grouping_id)
  filter_id <- unique(args$filter_id)
  dataset_type <- unique(args$dataset_type)
  
  replicate_matrix_dat <- dat %>%
    dplyr::select(replicate_id, value, pr_gene_symbol) %>%
    pivot_wider(id_cols = c(replicate_id),
                names_from = pr_gene_symbol,
                values_from = value,
                values_fn = function(x) median(x, na.rm = TRUE)) %>%
    column_to_rownames("replicate_id")
  # UMAP on raw data
  
  dat_distance <- as.matrix(dist(replicate_matrix_dat))
  umap_fit <- umap::umap(dat_distance, input = "dist", config=custom.config)
  
  umap_tbl <- tibble(UMAP_1 = umap_fit$layout[,1],
                     UMAP_2 = umap_fit$layout[,2],
                     replicate_id = rownames(dat_distance),
                     cell_id = str_split(replicate_id, "--", simplify = T)[,1]) %>%
                     # pert_iname = str_split(replicate_id, "--", simplify = T)[,2]) %>%
    left_join(informative_lbls, by= c("cell_id" = "lbs")) %>%
    mutate(cell_id = factor(cell_id)) %>%
    mutate(filter_id, grouping_id, dataset_type)
  
  numeric_replicate_matrix <- args$numeric_replicate_matrix # numeric_replicate_matrix <- args$numeric_replicate_matrix[[1]]
  # perform UMAP and PCA
  
  umap_conn_replicates_fit <- umap::umap(numeric_replicate_matrix, config=custom.config)
  res_pca <- prcomp(numeric_replicate_matrix, scale. = F, center = F)
  
  scaled_uwot_conn_replicates_fit <- uwot::umap(X = numeric_replicate_matrix,
                                                n_neighbors = n_neighbors_uwot, min_dist = 0.001,
                                                # metric = "correlation",
                                                verbose = TRUE, n_threads = no_cores) # no_cores defined in load.R
  
  
  scaled_numeric_replicate_matrix <- args$scaled_numeric_replicate_matrix # scaled_numeric_replicate_matrix <- args$scaled_numeric_replicate_matrix[[1]]
  # perform UMAP and PCA on SCALED data
  scaled_umap_conn_replicates_fit <- umap::umap(scaled_numeric_replicate_matrix, config=custom.config)
  scaled_res_pca <- prcomp(scaled_numeric_replicate_matrix, scale. = F, center = F)
  
  # stop()
  conn_replicates_tbl <- tibble(UMAP_1 = umap_conn_replicates_fit$layout[,1],
                                UMAP_2 = umap_conn_replicates_fit$layout[,2],
                                
                                PCA_1 = res_pca$x[,1],
                                PCA_2 = res_pca$x[,2],
                                
                                UWOT_1 = scaled_uwot_conn_replicates_fit[,1],
                                UWOT_2 = scaled_uwot_conn_replicates_fit[,2],
                                
                                UMAP_1_scaled = scaled_umap_conn_replicates_fit$layout[,1],
                                UMAP_2_scaled = scaled_umap_conn_replicates_fit$layout[,2],
                                
                                PCA_1_scaled = scaled_res_pca$x[,1],
                                PCA_2_scaled = scaled_res_pca$x[,2],
                                
                                replicate_id = rownames(numeric_replicate_matrix),
                                cell_id = str_split(replicate_id, "--", simplify = T)[,1]) %>%
                                # pert_iname = str_split(replicate_id, "--", simplify = T)[,2]) %>%
    left_join(informative_lbls, by= c("cell_id" = "lbs")) %>%
    mutate(cell_id = factor(cell_id)) %>%
    mutate(filter_id, grouping_id, dataset_type)
  
  final_res <- tibble(grouping_id = grouping_id, filter_id = filter_id,
                      n_neighbors = n_neighbors_uwot,
                      dataset_type = dataset_type,
                      matrix_raw_dat = list(replicate_matrix_dat),
                      umap_raw = list(umap_tbl),
                      umap_uwot_pca_connectivity = list(conn_replicates_tbl),
                      umap_obj = list(umap_conn_replicates_fit),
                      prcomp_obj = list(res_pca),
                      scaled_uwot_obj = list(scaled_uwot_conn_replicates_fit), 
                      umap_scaled_obj = list(scaled_umap_conn_replicates_fit),
                      prcomp_obj_scaled = list(scaled_res_pca))
  return(final_res)
}) %>%
  bind_rows()

cmbd_dat_and_res <- left_join(conn_res, umap_uwot_pca_results, by = c("grouping_id", "filter_id", "dataset_type"))
write_rds(cmbd_dat_and_res, file = file.path(supplemental_output_dir, "umap_uwot_pca_results-cell_cell.rds"))

umap_uwot_pca_data_and_results <- read_rds(file.path(supplemental_output_dir, "umap_uwot_pca_results-cell_cell.rds"))

pca_plots <- apply(umap_uwot_pca_data_and_results, MARGIN = 1, FUN = function(args){
  
  # args <- umap_uwot_pca_data_and_results[1,]
  filter_id <- unique(args$filter_id)
  grouping_id <- unique(args$grouping_id)
  dataset_type <- unique(args$dataset_type)
  
  numeric_replicate_matrix <- args$numeric_replicate_matrix # numeric_replicate_matrix <- args$numeric_replicate_matrix[[1]]
  ut_cells <- args$conn_df # ut_cells <- ut_cells[[1]]
  ut <- args$umap_uwot_pca_connectivity # ut <- args$umap_uwot_pca_connectivity[[1]]
  prcomp_obj <- args$prcomp_obj # prcomp_obj <- args$prcomp_obj[[1]]
  
  cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
  ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
  color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
  # pert_color_map <- deframe(pert_color_map_init %>% filter(pert_iname %in% ut$pert_iname))
  
  explained_var <- round(summary(prcomp_obj)$importance[2,c(1,2)]*100, 2)
  
  by_cell <- ggplot(ut, aes(PCA_1, PCA_2)) +
    geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 5,alpha = 0.75) + 
    theme_bw() +
    scale_fill_manual(values = color_map) +
    labs(fill = "Cell ID") +
    ggtitle(qq("PCA on @{filter_id} | @{dataset_type}")) +
    labs(caption = qq("@{dataset_type} connectivity profiles\nDefault PCA parameters")) +
    xlab(qq("standarized PC1 (@{explained_var[1]}% explained var.)")) +
    ylab(qq("standarized PC2 (@{explained_var[2]}% explained var.)")); by_cell
  
  pca_proper_plot <- fviz_pca_ind(prcomp_obj,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  ) + 
    ggtitle(qq("PCA on @{filter_id} | @{dataset_type}")) +
    labs(caption = qq("@{dataset_type} connectivity profiles\nDefault PCA parameters"))
  
  # by_pert <- ggplot(ut, aes(PCA_1, PCA_2)) +
  #   geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3,alpha = 0.75) + 
  #   theme_bw() +
  #   scale_fill_manual(values = pert_color_map) +
  #   labs(fill = "Perturbation") +
  #   ggtitle(qq("PCA on @{filter_id} | @{dataset_type}")) +
  #   labs(caption = qq("@{dataset_type} connectivity profiles\nDefault PCA parameters")) +
  #   xlab(qq("standarized PC1 (@{explained_var[1]}% explained var.)")) +
  #   ylab(qq("standarized PC2 (@{explained_var[2]}% explained var.)")); by_pert
  # 
  
  return(tibble(filter_id, grouping_id, dataset_type) %>% mutate(cell_g = list(by_cell), pca_proper = list(pca_proper_plot))) # , pert_g = list(by_pert)
}) %>%
  bind_rows() %>%
  pivot_longer(cols = cell_g:pca_proper, names_to = "plot_type") %>% # :pert_g
  mutate(dir_type = "pca",
         plot_dir = file.path(supplemental_output_dir, dir_type),
         fn_name = file.path(plot_dir, str_c(str_c(filter_id, dataset_type, plot_type, sep = "-"),".png", sep = "")))
# (pca_plots$cell_g[[1]] | pca_plots$cell_g[[2]]) /
#   (pca_plots$cell_g[[3]] | pca_plots$cell_g[[4]])

# (pca_plots$cell_g[[1]] | pca_plots$pert_g[[1]])
# (pca_plots$cell_g[[3]] | pca_plots$pert_g[[3]])

umap_plots <- apply(umap_uwot_pca_data_and_results, MARGIN = 1, FUN = function(args){
  # args <- umap_uwot_pca_data_and_results[1,]
  filter_id <- unique(args$filter_id)
  grouping_id <- unique(args$grouping_id)
  dataset_type <- unique(args$dataset_type)
  
  ut_cells <- args$conn_df # ut_cells <- args$conn_df[[1]]
  ut <- args$umap_uwot_pca_connectivity #  ut <- args$umap_uwot_pca_connectivity[[1]]
  n_neighbors <- args$n_neighbors
  
  cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
  ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
  color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
  # pert_color_map <- deframe(pert_color_map_init %>% filter(pert_iname %in% ut$pert_iname))
  
  by_cell <- ggplot(ut, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 3,alpha = 0.75) + 
    theme_bw() +
    scale_fill_manual(values = color_map) +
    labs(fill = "Cell ID") +
    ggtitle(qq("UMAP on @{filter_id} | @{dataset_type}")) +
    labs(caption = qq("Euclidean distance on @{dataset_type} connectivity profiles\nDefault UMAP parameters\nn_neighbors = @{n_neighbors}"))
  
  # by_pert <- ggplot(ut, aes(UMAP_1, UMAP_2)) +
  #   geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3,alpha = 0.75) +
  #   theme_bw() +
  #   scale_fill_manual(values = pert_color_map) +
  #   labs(fill = "Perturbation") +
  #   ggtitle(qq("UMAP on @{filter_id} | @{dataset_type}")) +
  #   labs(caption = qq("Euclidean distance on @{dataset_type} connectivity profiles\nDefault UMAP parameters\nn_neighbors = @{n_neighbors}"))
  # 
  return(tibble(filter_id, grouping_id, dataset_type) %>% mutate(cell_g = list(by_cell))) # , pert_g = list(by_pert)
}) %>%
  bind_rows()%>%
  pivot_longer(cols = cell_g, names_to = "plot_type") %>% #:pert_g
  mutate(dir_type = "umap",
         plot_dir = file.path(supplemental_output_dir, dir_type),
         fn_name = file.path(plot_dir, str_c(str_c(filter_id, dataset_type, plot_type, sep = "-"),".png", sep = "")))

uwot_plots <- apply(umap_uwot_pca_data_and_results, MARGIN = 1, FUN = function(args){
  filter_id <- unique(args$filter_id)
  grouping_id <- unique(args$grouping_id)
  dataset_type <- unique(args$dataset_type)
  
  ut_cells <- args$conn_df
  ut <- args$umap_uwot_pca_connectivity
  n_neighbors <- args$n_neighbors
  
  cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
  ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
  color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
  # pert_color_map <- deframe(pert_color_map_init %>% filter(pert_iname %in% ut$pert_iname))
  
  by_cell <- ggplot(ut, aes(UWOT_1, UWOT_2)) +
    geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 3,alpha = 0.75) + 
    theme_bw() +
    scale_fill_manual(values = color_map) +
    labs(fill = "Cell ID") +
    ggtitle(qq("UWOT on @{filter_id} | @{dataset_type}")) +
    labs(caption = qq("Correlation distance on @{dataset_type} connectivity profiles\nDefault UWOT parameters\nn_neighbors = @{n_neighbors}"))
  
  # by_pert <- ggplot(ut, aes(UWOT_1, UWOT_2)) +
  #   geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3,alpha = 0.75) +
  #   theme_bw() +
  #   scale_fill_manual(values = pert_color_map) +
  #   labs(fill = "Perturbation") +
  #   ggtitle(qq("UWOT on @{filter_id} | @{dataset_type}")) +
  #   labs(caption = qq("Correlation distance on @{dataset_type} connectivity profiles\nDefault UWOT parameters\nn_neighbors = @{n_neighbors}"))
  return(tibble(filter_id, grouping_id, dataset_type) %>% mutate(cell_g = list(by_cell))) # , pert_g = list(by_pert)
}) %>%
  bind_rows() %>%
  pivot_longer(cols = cell_g, names_to = "plot_type") %>% # pert_g
  mutate(dir_type = "uwot",
         plot_dir = file.path(supplemental_output_dir, dir_type),
         fn_name = file.path(plot_dir, str_c(str_c(filter_id, dataset_type, plot_type, sep = "-"),".png", sep = "")))

dir.create(unique(pca_plots$plot_dir), recursive = T, showWarnings = F)
dir.create(unique(umap_plots$plot_dir), recursive = T, showWarnings = F)
dir.create(unique(uwot_plots$plot_dir), recursive = T, showWarnings = F)


walk2(.x = as.list(pca_plots$fn_name), .y = pca_plots$value, .f = function(fn, plot_) ggsave(filename = fn, plot = plot_, device = "png"))
walk2(.x = as.list(umap_plots$fn_name), .y = umap_plots$value, .f = function(fn, plot_) ggsave(filename = fn, plot = plot_, device = "png"))
walk2(.x = as.list(uwot_plots$fn_name), .y = uwot_plots$value, .f = function(fn, plot_) ggsave(filename = fn, plot = plot_, device = "png"))


# pca_plots$value[[1]] | umap_plots$value[[1]] | uwot_plots$value[[1]]


# # determine the number of neighbors...?
# uwot_plots <- apply(conn_res , MARGIN = 1, FUN = function(args){
#   # args <- conn_res[1,]
#   dat <- args$dat # dat <- args$dat[[1]]
#   grouping_id <- unique(args$grouping_id)
#   filter_id <- unique(args$filter_id)
#   dataset_type <- unique(args$dataset_type)
#   
#   numeric_replicate_matrix <- args$numeric_replicate_matrix # numeric_replicate_matrix <- args$numeric_replicate_matrix[[1]]
#   
#   ut_cells <- args$conn_df # ut_cells <- args$conn_df[[1]] 
#   
#   cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
#   ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
#   color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
#   
#   numeric_raw_replicate_matrix <- dat %>%
#     dplyr::select(replicate_id, pr_gene_symbol, value) %>%
#     pivot_wider(id_cols = replicate_id,
#                 names_from = pr_gene_symbol,
#                 values_from = value,
#                 values_fn = function(x) median(x, na.rm = TRUE)) %>%
#     column_to_rownames("replicate_id")
#   numeric_cor_replicate_matrix <- cor(t(numeric_raw_replicate_matrix), use = "pairwise.complete.obs")
#   
#   uwot_cor_replicates_fit <- map_dfr(c(15, 25, 50, 75), .f = function(i) {
#     metric_ <- "euclidean" #"correlation"
#     uwot_i_res_temp <- uwot::umap(X = numeric_cor_replicate_matrix, 
#                                   n_neighbors = i, min_dist = 0.001,
#                                   metric = metric_,
#                                   verbose = TRUE, n_threads = no_cores)  # no_cores defined in load.R
#     replicate_id <- rownames(uwot_i_res_temp)
#     cell_id <- str_split(replicate_id, "--", simplify = T)[,1]
#     pert_iname <- str_split(replicate_id, "--", simplify = T)[,2]
#     
#     uwot_i_res <- tibble(UWOT_1 = uwot_i_res_temp[,1], 
#                          UWOT_2 = uwot_i_res_temp[,2],
#                          replicate_id = replicate_id,
#                          cell_id = cell_id,
#                          pert_iname = pert_iname)
#     res <- tibble(n_neighbor = i, metric = metric_, uwot_res = list(uwot_i_res))
#     return(res)
#   })
#   
#   # numeric_replicate_matrix <- numeric_replicate_matrix[[1]]
#   uwot_conn_replicates_fit <- map_dfr(c(15, 25, 50, 75), .f = function(i) {
#     metric_ <- "euclidean" #"correlation"
#     uwot_i_res_temp <- uwot::umap(X = numeric_replicate_matrix, 
#                                   n_neighbors = i, min_dist = 0.001,
#                                   metric = metric_,
#                                   verbose = TRUE, n_threads = no_cores)  # no_cores defined in load.R
#     replicate_id <- rownames(uwot_i_res_temp)
#     cell_id <- str_split(replicate_id, "--", simplify = T)[,1]
#     pert_iname <- str_split(replicate_id, "--", simplify = T)[,2]
#     
#     uwot_i_res <- tibble(UWOT_1 = uwot_i_res_temp[,1], 
#                          UWOT_2 = uwot_i_res_temp[,2],
#                          replicate_id = replicate_id,
#                          cell_id = cell_id,
#                          pert_iname = pert_iname)
#     res <- tibble(n_neighbor = i, metric = metric_, uwot_res = list(uwot_i_res))
#     return(res)
#   })
#   
#   by_cell_lst <- apply(uwot_conn_replicates_fit, MARGIN = 1, FUN = function(r){
#     ut <- r$uwot_res
#     n_neighbors <- unique(r$n_neighbor)
#     metric <- unique(r$metric)
#     by_cell <- ggplot(ut, aes(UWOT_1, UWOT_2)) +
#       geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 3) + 
#       theme_bw() +
#       scale_fill_manual(values = color_map) +
#       labs(fill = "Cell ID") +
#       ggtitle(qq("UWOT on @{filter_id} | @{dataset_type}")) +
#       labs(caption = qq("Connectivity profiles\nDefault UWOT parameters\nN_neighbors = @{n_neighbors}\nMetric = @{metric}"))
#     return(by_cell)
#   })
#   (by_cell_lst[[1]] | by_cell_lst[[2]]) / (by_cell_lst[[3]] | by_cell_lst[[4]]) 
#   
#   by_pert_lst <- apply(uwot_conn_replicates_fit, MARGIN = 1, FUN = function(r){
#     ut <- r$uwot_res
#     n_neighbors <- unique(r$n_neighbor)
#     metric <- unique(r$metric)
#     by_cell <- ggplot(ut, aes(UWOT_1, UWOT_2)) +
#       geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3) + 
#       theme_bw() +
#       # scale_fill_manual(values = color_map) +
#       labs(fill = "Cell ID") +
#       ggtitle(qq("UWOT on @{filter_id} | @{dataset_type}")) +
#       labs(caption = qq("Connectivity profiles\nDefault UWOT parameters\nN_neighbors = @{n_neighbors}\nMetric = @{metric}"))
#     return(by_cell)
#   })
#   (by_pert_lst[[1]] | by_pert_lst[[2]]) / (by_pert_lst[[3]] | by_pert_lst[[4]]) 
#   
#   
#   
#   corr_pca <- prcomp(numeric_cor_replicate_matrix, scale. = FALSE, center = FALSE)
#   summary(corr_pca)
# ggbiplot(corr_pca, choices = c(1,2),
#          groups = str_split(rownames(numeric_cor_replicate_matrix), "--", simplify = T)[,1],
#          var.axes=FALSE, ellipse = TRUE) +
#   scale_color_manual(values = color_map) +
#   theme_bw()
#   
#   conn_pca <- prcomp(numeric_replicate_matrix, scale. = FALSE, center = FALSE)
#   summary(conn_pca)
#   ggbiplot(conn_pca, choices = c(1,2),
#            groups = str_split(rownames(numeric_replicate_matrix), "--", simplify = T)[,1], 
#            var.axes=FALSE, ellipse = TRUE) + 
#     scale_color_manual(values = color_map) +
#     theme_bw()
#   
#   
#   # conn_replicates_tbl <- tibble(UWOT_1 = scaled_uwot_conn_replicates_fit[,1],
#   #                               UWOT_2 = scaled_uwot_conn_replicates_fit[,2],
#   #                               replicate_id = rownames(numeric_replicate_matrix),
#   #                               cell_id = str_split(replicate_id, "--", simplify = T)[,1],
#   #                               pert_iname = str_split(replicate_id, "--", simplify = T)[,2]) %>%
#   #   left_join(informative_lbls, by= c("cell_id" = "lbs")) %>%
#   #   mutate(cell_id = factor(cell_id)) %>%
#   #   mutate(filter_id, grouping_id, dataset_type)
#   # 
#   # final_res <- tibble(grouping_id = grouping_id, filter_id = filter_id,
#   #                     dataset_type = dataset_type,
#   #                     matrix_raw_dat = list(replicate_matrix_dat),
#   #                     umap_raw = list(umap_tbl),
#   #                     umap_uwot_pca_connectivity = list(conn_replicates_tbl))
#   # return(final_res)
#   
# }) %>% 
#   bind_rows()

# conn_ALL_results <- apply(analysis_dat[c(1,3),], MARGIN = 1, function(args){
#   # args <- analysis_dat[1,]
#   d <- args$data
#   # d <- args$data[[1]]
#   grouping_id <- "ALL"
#   filter_id <- "ALL"
#   # filter_id <- args$filter_vars[[1]]
#   
#   print(grouping_id)
#   print(filter_id)
#   # stop()
#   dat <- d 
#   
#   dataset_type <- unique(dat$which_dat)
#   
#   # dat <- dat %>% filter(cell_id %in% c("MCF7", "HAoSMC", "HUVEC"))
#   
#   # now do this with connectivity values
#   message("Computing correlation...")
#   cor_dat <- run_corr(dat)$matrix
#   unique_cs <- unique(rownames(cor_dat))
#   grp_names <- get_grp_names_from_matrix(cor_dat, unique_cs, sep_ = "--")
#   
#   cs_idx <- seq_len(length(rownames(cor_dat)))
#   # collapse replicates by drug - cell pair, not just cell
#   cs_tib <- tibble(unique_cs, cs_idx, grp_names) %>%
#     mutate(cell_ = str_split(unique_cs, "--", simplify = T)[,1],
#            pert_ = str_split(unique_cs, "--", simplify = T)[,2],
#            g_cell_pert_ = str_c(cell_, pert_, sep = "--"))
#   
#   looper <- expand.grid(
#     rep_a = unique(cs_tib$g_cell_pert_), 
#     rep_b = unique(cs_tib$g_cell_pert_)
#   ) %>%
#     as_tibble()
#   
#   # connectivity of each cell drug (all 3 replicates) vs all others
#   ggcorr(cor_dat)
#   
#   with_progress(expr = {
#     conn_progress <- progressr::progressor(steps = nrow(looper))
#     conn_df <- future_map_dfr(1:nrow(looper), .f = function(z){
#       # message(looper$rep_a[z], "  ", looper$rep_b[z])
#       res <- calc_conn_replicates(a = looper$rep_a[z], b = looper$rep_b[z],
#                                   my_mat = cor_dat, my_tib = cs_tib,
#                                   use_bootstrap = F)
#       conn_progress()
#       return(res)
#     }) %>%
#       bind_rows() %>%
#       mutate(cell_a = str_split(column_a, pattern = "--", simplify = T)[,1],
#              cell_b = str_split(column_b, pattern = "--", simplify = T)[,1], 
#              pert_a = str_split(column_a, pattern = "--", simplify = T)[,2],
#              pert_b = str_split(column_a, pattern = "--", simplify = T)[,2])
#     gc(verbose = F)
#   })
#   
#   replicate_matrix <- conn_df %>%
#     dplyr::select(column_a, column_b, connectivity) %>%
#     pivot_wider(id_cols = column_a,
#                 names_from = column_b,
#                 values_from = connectivity,
#                 values_fn = function(x) median(x, na.rm = TRUE))
#   
#   numeric_replicate_matrix <- replicate_matrix %>%
#     column_to_rownames("column_a") %>%
#     dplyr::select(where(is.numeric))
#   
#   # scale and 0-center the data
#   scaled_numeric_replicate_matrix <- scale(numeric_replicate_matrix)
#   
#   final_result <- tibble(grouping_id = grouping_id, filter_id = filter_id, dataset_type = dataset_type,
#                          dat = list(dat),
#                          conn_df = list(conn_df),
#                          replicate_matrix = list(replicate_matrix),
#                          numeric_replicate_matrix = list(numeric_replicate_matrix),
#                          scaled_numeric_replicate_matrix = list(scaled_numeric_replicate_matrix))
#   return(final_result)
# }) %>%
#   bind_rows()


# make_gct <- apply(analysis_dat, MARGIN = 1, function(args){
#   # args <- analysis_dat[1,]
#   d <- args$data
#   # d <- args$data[[1]]
#   grouping_id <- args$grouping_var
#   filter_id <- args$filter_vars
#   # filter_id <- args$filter_vars[[1]]
#   
#   # stop()
#   dat <- d %>%
#     filter(!!sym(grouping_id) == filter_id) %>%
#     filter(pert_iname %in% my_perts)
#   
#   dataset_type <- unique(dat$which_dat)
#   
#   mat_gct <- dat %>%
#     dplyr::select(replicate_id, pr_gene_symbol, value) %>%
#     pivot_wider(id_cols = pr_gene_symbol, names_from = replicate_id, values_from = value,
#                 values_fn = function(x) median(x, na.rm = TRUE)) %>%
#     column_to_rownames("pr_gene_symbol") %>%
#     as.matrix()
#   
#   rdesc <- data.frame(pr_gene_symbol = rownames(mat_gct))
#   cdesc <- data.frame(cnames = colnames(mat_gct)) %>%
#     mutate(cell_id = str_split(cnames, pattern = "--", simplify = T)[,1],
#            pert_iname = str_split(cnames, pattern = "--", simplify = T)[,2],
#            pert_class = str_split(str_split(cnames, pattern = "--", simplify = T)[,3], pattern = "::", simplify = T)[,1],
#            batch = str_split(cnames, pattern = "::", simplify = T)[,2],
#            det_plate = str_split(cnames, pattern = "::", simplify = T)[,3])
#   gct_file_path <- file.path(supplemental_output_dir, str_c(str_c(grouping_id, filter_id, dataset_type, sep = "--"), ".gct"))
#   
#   new_gct <- new("GCT", mat = mat_gct, rdesc = rdesc , cdesc = cdesc)
#   write_gct(new_gct, ofile = gct_file_path)
# })