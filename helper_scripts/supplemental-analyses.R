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
# dendro_cut_thresh <- 0.75

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

version_ <- 2
dir_name_supp <- qq("supplemental_output_dmso_subtracted_v@{version_}-thresh-@{dendro_cut_thresh}")
# dir_name_supp <- "supplemental_output_dmso_subtracted_v2-thresh-060"
supplemental_output_dir <- file.path(output_directory, dir_name_supp)
dir.create(supplemental_output_dir, showWarnings = F, recursive = T)

informative_lbls <- create_informative_labels()$informative_labels_df 
pert_color_map_init <- generate_color_palette()

# SWITCH THIS OFF FOR NO DEBUG
# analysis_dat_temp <- analysis_dat %>% slice(1)

#' @note make a heatmap out of correlation or connectivity matrix
#' @param m matrix like object, e.g. correlation or connectivity matrix
#' @param pert_arg character string, like DMSO
#' @param filter_arg character string, like "DMSO" or "Kinase inhibitor"
#' @param dataset_arg character string, like "P100" or "GCP"
#' @param height_unit numeric, height of heatmap cell
#' @param width_unit numeric, width of heatmap cell
#' @param use_cell_fun boolean, put values in cell heatmap or not
#' @param clustering_method character, one of the clustering methods valid for 'hclust'
#' @param cut_tree_thresh numeric, percent of the maximum height of the tree to cut, to create clusters
#' @return heatmap object
make_easy_heatmap <- function(m, pert_arg = "DMSO", filter_arg = "DMSO", dataset_arg = "P100",
                              height_unit = 2, width_unit = 2, use_cell_fun = FALSE, 
                              clustering_method  = "euclidean", cut_tree_thresh = 0.6){
  # DEBUG:
  # m = cor_dat; pert_arg = pert_arg_; filter_arg = filter_id;height_unit = 7.5; width_unit = 7.5; use_cell_fun = TRUE; dataset_arg = dataset_type
  unique_rnames_dat <- rownames(m)
  unique_colnames_dat <- colnames(m)
  
  sorted_unique_rnames_dat <- sort(unique_rnames_dat)
  sorted_unique_colnames_dat <- sort(unique_colnames_dat)
  
  rnames_split <- tibble(sorted_unique_rnames_dat) %>%
    mutate(rsplit_ = str_split(sorted_unique_rnames_dat, pattern = "\\.|-", simplify = TRUE)[,1])
  cnames_split <- tibble(sorted_unique_colnames_dat) %>%
    mutate(csplit_ = str_split(sorted_unique_colnames_dat, pattern = "\\.|-", simplify = TRUE)[,1])
  rsplit <- factor(rnames_split$rsplit_)
  rsplit_int <- as.integer(rsplit); rsplit_int
  csplit <- factor(cnames_split$csplit_)
  csplit_int <- as.integer(csplit); csplit_int
  
  my_fontsize <- 6
  
  #' @note convenience function to get all dendro info
  #' @param mat numeric matrix for which clustering will be performed on the rows
  #' @param clustering_method This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall"
  get_dendro_lst <- function(mat, clst_mtd = "euclidean") {
    
    clst <- hclust(d = factoextra::get_dist(mat, method = clst_mtd), method = "average")
    dend_temp <- as.dendrogram(clst)
    
    ct <- cutree(dend_temp, h = cut_tree_thresh*max(clst$height))
    
    clust_palette <- hcl.colors(n = max(ct), palette = "ag_Sunset")
    
    cluster_identities <- enframe(ct)  %>%
      mutate(cluster_name = str_split(name, "--", simplify = TRUE)[,1]) %>%
      group_by(value) %>%
      summarize(
        cluster_name = str_c(sort(unique(cluster_name)),
                             collapse = ","
        ),
        .groups = "keep"
      ) %>%
      ungroup() %>%
      mutate(clust_colors = clust_palette) 
    
    full_clst_id <- enframe(ct)  %>%
      mutate(cluster_name = str_split(name, "--", simplify = TRUE)[,1]) %>%
      dplyr::select(-cluster_name) %>%
      left_join(cluster_identities, by = c("value")) %>%
      mutate(cell_id = str_split(name, "--", simplify = TRUE)[,1]) %>%
      left_join(informative_lbls, by = c("cell_id" = "lbs"))
    
    block_clust_vals <- cluster_identities$clust_colors
    names(block_clust_vals) <- cluster_identities$cluster_name
    
    block_cell_vals <- full_clst_id$cell_individual_color
    names(block_cell_vals) <- full_clst_id$cell_id
    
    block_rep_vals_df <- cluster_identities %>%
      distinct(cluster_name, clust_colors)
    block_rep_vals <- block_rep_vals_df$clust_colors
    names(block_rep_vals) <- block_rep_vals_df$cluster_name
    
    dend <- color_branches(dend_temp, h = cut_tree_thresh*max(clst$height), col = block_rep_vals)
    
    # from https://stackoverflow.com/questions/68837845/map-colors-to-labels-in-ggplot-based-on-background-color
    # super helpful to choose label color
    hcl <- farver::decode_colour(block_clust_vals, "rgb", "hcl")
    block_clust_labs <- ifelse(hcl[, "l"] > 50, "black", "white")
    names(block_clust_labs) <- cluster_identities$cluster_name
    
    final_split <- max(ct)
    return(list(clst, dend, ct, block_clust_vals, block_clust_labs, block_cell_vals, final_split, cluster_identities, full_clst_id) %>%
             set_names("hclust_obj", "dendrogram", "cut_values", "block_clust_vals", "block_clust_labs", "block_cell_vals", "final_split", "cluster_identities", "full_clst_ids"))
  }
  
  row_dendro_lst <- get_dendro_lst(mat = m, clst_mtd = clustering_method)
  col_dendro_lst <- get_dendro_lst(mat = t(m), clst_mtd = clustering_method)
  
  # get label orders from dendrogram -- this isn't necessary right now but was hard to come up with, so can't delete
  # row_cell_labels <- left_join(tibble(new_order = row_dendro_lst$hclust_obj$labels[row_dendro_lst$hclust_obj$order]) %>%
  #                                mutate(name = new_order), row_dendro_lst$full_clst_ids, by= "name") %>%
  #   dplyr::select(value, name, cell_id, cell_individual_color) 
  # 
  # col_cell_labels <- left_join(tibble(new_order = col_dendro_lst$hclust_obj$labels[col_dendro_lst$hclust_obj$order]) %>%
  #                                mutate(name = new_order), col_dendro_lst$full_clst_ids, by= "name") %>%
  #   dplyr::select(value, name, cell_id, cell_individual_color) 
  
  # unit_height_width <- min(c(width_unit, height_unit))
  
  name_ht <- str_c(dataset_arg, pert_arg,  toupper(substr(clustering_method, 1, 1)), sep ="-") # "matrix" 
  
  my_cell_fun <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", m[i, j]), x, y, gp = gpar(fontsize = my_fontsize))
  }
  
  block_fill <- col_dendro_lst$full_clst_ids$cell_individual_color
  names(block_fill) <- col_dendro_lst$full_clst_ids$cell_id
  
  top_annot <- HeatmapAnnotation(`Cell Type` = names(block_fill),
                                 col = list(`Cell Type` = block_fill),
                                 gp = gpar(col = "black"),
                                 which = "col",
                                 `Cell Type` = list(title = "Cell Type"),
                                 # annotation_legend_param = list(`Cell Type` = list(direction = "horizontal")), 
                                 show_legend = TRUE, show_annotation_name = FALSE)
  # 
  left_annot <- rowAnnotation(`Cell Type` = names(block_fill),
                              col = list(`Cell Type` = block_fill),
                              gp = gpar(col = "black"),
                              `Cell Type` = list(title = "Cell Type"),
                              show_legend = FALSE, show_annotation_name = FALSE)
  
  
  
  if (use_cell_fun){
    ht <- Heatmap(m, 
                  name = name_ht, 
                  top_annotation = top_annot, 
                  left_annotation = left_annot,
                  cell_fun = my_cell_fun, 
                  
                  row_split = row_dendro_lst$final_split,
                  column_split = col_dendro_lst$final_split,
                  
                  row_dend_reorder = FALSE,
                  column_dend_reorder = FALSE,
                  row_title_rot = 0, column_title_rot = 90,
                  
                  cluster_rows = row_dendro_lst$dend,
                  cluster_columns = col_dendro_lst$dend,
                  
                  border = TRUE,
                  show_heatmap_legend = TRUE,
                  column_names_gp = grid::gpar(fontsize = my_fontsize),
                  row_names_gp = grid::gpar(fontsize = my_fontsize),
                  rect_gp = gpar(col = "white", lwd = 0.25),
                  width = ncol(sorted_cor_dat)*unit(width_unit, "mm"), 
                  height = nrow(sorted_cor_dat)*unit(height_unit, "mm"),
                  raster_device = "png",  heatmap_legend_param = list(legend_direction = "horizontal")); ht
  } else {
    ht <- Heatmap(m, 
                  name = name_ht, 
                  top_annotation = top_annot, 
                  left_annotation = left_annot,
                  # cell_fun = my_cell_fun, 
                  
                  row_split = row_dendro_lst$final_split,
                  column_split = col_dendro_lst$final_split,
                  
                  row_dend_reorder = FALSE,
                  column_dend_reorder = FALSE,
                  row_title_rot = 0, column_title_rot = 90,
                  
                  cluster_rows = row_dendro_lst$dend,
                  cluster_columns = col_dendro_lst$dend,
                  
                  border = TRUE,
                  show_heatmap_legend = TRUE,
                  column_names_gp = grid::gpar(fontsize = my_fontsize),
                  row_names_gp = grid::gpar(fontsize = my_fontsize),
                  rect_gp = gpar(col = "white", lwd = 0.25),
                  width = ncol(sorted_cor_dat)*unit(width_unit, "mm"), 
                  height = nrow(sorted_cor_dat)*unit(height_unit, "mm"),
                  raster_device = "png",  heatmap_legend_param = list(legend_direction = "horizontal")); ht
  }
  
  return(ht)
}


#' @note make a heatmap out of correlation or connectivity matrix, no clustering
#' @param m matrix like object
#' @param pert_arg character string, like DMSO
#' @param filter_arg character string, like "DMSO" or "Kinase inhibitor"
#' @param dataset_arg character string, like "P100" or "GCP"
#' @param height_unit height of heatmap cell
#' @param width_unit width of heatmap cell
#' @param use_cell_fun put values in cell heatmap
#' @return heatmap object
make_easy_heatmap_no_clust <- function(m, pert_arg = "DMSO", filter_arg = "DMSO", dataset_arg = "P100",
                                       height_unit = 2, width_unit = 2, use_cell_fun = FALSE){
  
  unique_rnames_dat <- rownames(m)
  unique_colnames_dat <- colnames(m)
  
  sorted_unique_rnames_dat <- sort(unique_rnames_dat)
  sorted_unique_colnames_dat <- sort(unique_colnames_dat)
  
  rnames_split <- tibble(sorted_unique_rnames_dat) %>%
    mutate(rsplit_ = str_split(sorted_unique_rnames_dat, pattern = "\\.|-", simplify = TRUE)[,1])
  cnames_split <- tibble(sorted_unique_colnames_dat) %>%
    mutate(csplit_ = str_split(sorted_unique_colnames_dat, pattern = "\\.|-", simplify = TRUE)[,1])
  rsplit <- factor(rnames_split$rsplit_)
  rsplit_int <- as.integer(rsplit); rsplit_int
  csplit <- factor(cnames_split$csplit_)
  csplit_int <- as.integer(csplit); csplit_int
  
  my_fontsize <- 6
  # if (filter_arg != "DMSO") my_fontsize <- my_fontsize /1.5
  
  # honestly, just always sort the matrix before...
  sorted_cor_dat <- m[sorted_unique_rnames_dat, sorted_unique_colnames_dat, drop = FALSE]
  
  name_ht <- str_c(dataset_arg, pert_arg, sep ="-") # "matrix"
  
  my_cell_fun <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", sorted_cor_dat[i, j]), x, y, gp = gpar(fontsize = my_fontsize))
  }
  
  
  if (use_cell_fun){
    ht <- Heatmap(sorted_cor_dat, 
                  name = name_ht, 
                  cell_fun = my_cell_fun,
                  row_split = rsplit,
                  column_split = csplit,
                  row_title_rot = 0, column_title_rot = 90,
                  cluster_rows = F,
                  cluster_columns = F,
                  border = TRUE,
                  column_names_gp = grid::gpar(fontsize = my_fontsize),
                  row_names_gp = grid::gpar(fontsize = my_fontsize),
                  rect_gp = gpar(col = "white", lwd = 0.25),
                  width = ncol(sorted_cor_dat)*unit(width_unit, "mm"), 
                  height = nrow(sorted_cor_dat)*unit(height_unit, "mm"),
                  raster_device = "png", heatmap_legend_param = list(legend_direction = "horizontal"));
  } else {
    ht <- Heatmap(sorted_cor_dat, 
                  name = name_ht, 
                  row_split = rsplit,
                  column_split = csplit,
                  cell_fun = my_cell_fun,
                  row_title_rot = 0, column_title_rot = 90,
                  cluster_rows = F,
                  cluster_columns = F,
                  border = TRUE,
                  column_names_gp = grid::gpar(fontsize = my_fontsize),
                  row_names_gp = grid::gpar(fontsize = my_fontsize),
                  rect_gp = gpar(col = "white", lwd = 0.25),
                  width = ncol(sorted_cor_dat)*unit(width_unit, "mm"), 
                  height = nrow(sorted_cor_dat)*unit(height_unit, "mm"),
                  raster_device = "png",  heatmap_legend_param = list(legend_direction = "horizontal"));
  }
  
  return(ht)
}

#' @note write the heatmap to a file, using supplemental_output_dir as base dir
#' @param plot_df corr_results or conn_results
#' @param top_directory corr_ht or conn_ht, to separate out the plots into directories by type (correlation vs connectivity)
write_plot_to_file <- function(plot_df, top_directory = "corr_ht"){
  x_df <- plot_df %>%
    mutate(output_dir = file.path(supplemental_output_dir, top_directory, dataset_type, ht_type, clust_methods),
           output_fn = file.path(output_dir, str_c(str_c(dataset_type, pert_iname, filter_id, ht_type, top_directory, sep = "-"), ".png")))
  
  walk(unique(x_df$output_dir), .f = function(x) dir.create(x, showWarnings = F, recursive = T))
  
  tic()
  pwalk(.l = list(as.list(x_df$ht), as.list(x_df$output_fn), 
                  as.list(x_df$pert_iname), as.list(x_df$dataset_type), 
                  as.list(x_df$clust_methods)), 
        .f = function(ht_, path_, pert_name_, d_type_, clust_method_){
          # debug: lst <- list(as.list(x_df$ht), as.list(x_df$output_fn), as.list(x_df$pert_iname), as.list(x_df$dataset_type), as.list(x_df$clust_methods))
          # i_ = 19; ht_ <- lst[[1]][[i_]]; path_ = lst[[2]][[i_]]; pert_name_ = lst[[3]][[i_]]; d_type_ = lst[[4]][[i_]]; clust_method_ = lst[[5]][[i_]]; 
          
          print(path_)
          d_ht_ <- draw(ht_) # , padding = unit(c(2, 2, 15, 2), "mm")
          width_ <- ComplexHeatmap:::width(d_ht_) + unit(35, "mm"); width_
          height_ <- ComplexHeatmap:::height(d_ht_) + unit(50, "mm"); height_
          
          # heatmap_title_ <- str_c(pert_name_, " | ", d_type_, " | ", clust_method_, sep = "")
          
          png(file = path_, units = "mm", width = width_, height = height_, res = 300)
          draw(ht_, heatmap_legend_side = "bottom") #  padding = unit(c(2, 2, 15, 2), "mm")
          # decorate_heatmap_body(heatmap = names(d_ht_), {
          #   grid.text(heatmap_title_, y = height_/3.85, x = width_/3, just = "bottom", gp = gpar(fontsize = 20, font = 2))
          # })
          dev.off()
          
        })
  toc()
}

#' @note takes a df of corr or conn matrices and combines them into a single matrix
#' @param single_res_df like conn_results or corr_results
#' @param filter_mark the filter_id to select, e.g. Epigenentic or Kinase inhibitor
#' @param column_to_select the column name containing the cor_dat or conn_dat info
get_combined_matrix <- function(single_res_df, filter_mark = "Kinase Inhibitor", column_to_select = "cor_dat"){
  matrices <- single_res_df %>%
    filter(filter_id == filter_mark) %>% 
    pluck(column_to_select)
  
  pert_names <- single_res_df %>%
    filter(filter_id == filter_mark) %>%
    .$pert_iname
  
  # need to fill out the 
  combined_matrix <- lapply(1:length(matrices), FUN = function(i){
    mat <- matrices[[i]]
    pert_name <- pert_names[i]
    
    abbrev_colnames <- make.unique(str_split(colnames(mat), "--", simplify = TRUE)[,1])
    abbrev_rownames <- make.unique(str_split(rownames(mat), "--", simplify = TRUE)[,1])
    
    colnames(mat) <- str_c(abbrev_colnames, pert_name, sep = "-")
    rownames(mat) <- abbrev_rownames
    df <- as.data.frame(mat) %>%
      rownames_to_column() %>%
      as_tibble() %>%
      pivot_longer(colnames(mat), names_to = "colname")
    return(df)
  }) %>%
    bind_rows() %>%
    pivot_wider(id_cols = rowname, names_from = colname, values_from = value) %>%
    column_to_rownames("rowname") %>%
    as.matrix()
  
  return(combined_matrix)
}

# subtract DMSO from pert
new_analysis_dat <- analysis_dat %>% 
  mutate(dmso_subtracted_data = map(data, .f = function(d){
    # d <- analysis_dat$data[[1]]
    #  d <- analysis_dat$data[[3]]
    temp_veh_d <- d %>%
      filter(pert_iname == "DMSO") %>%
      # dplyr::select(abbrev_replicate_id, cell_id, pert_iname, det_plate, pr_gene_symbol, value) %>%
      rename(dmso_value = value, dmso_pert_name = pert_iname, 
             dmso_abbrev_replicate_id = abbrev_replicate_id) ; temp_veh_d
    
    summarized_temp_veh_d <- temp_veh_d %>%
      group_by(cell_id, det_plate, pr_gene_symbol) %>%
      summarize(median_cell_id_plate_dmso_value = median(dmso_value, na.rm= TRUE), .groups = "keep") 
    
    # wierdly, this plate doesn't have a DMSO group - take most closely related plate and add in
    if ("G-0016R" %in% unique(summarized_temp_veh_d$det_plate)){
      summarized_temp_veh_d <- summarized_temp_veh_d %>% 
        filter(det_plate == "G-0016R") %>%
        mutate(det_plate = "G-0016") %>%
        bind_rows(summarized_temp_veh_d) %>%
        distinct()
    }
    
    temp_pert_d_1 <- d %>%
      filter(pert_iname != "DMSO") %>%
      # dplyr::select(abbrev_replicate_id, cell_id, det_plate, pert_iname, pr_gene_symbol, value) %>%
      rename(pert_value = value, pert_abbrev_replicate_id = abbrev_replicate_id); temp_pert_d_1
    
    # setdiff(summarized_temp_veh_d  %>% ungroup() %>% distinct(det_plate) %>% .$det_plate %>% sort(), 
    #         temp_pert_d_1 %>% ungroup() %>% distinct(det_plate) %>% .$det_plate %>% sort())
    # setdiff(temp_pert_d_1 %>% ungroup() %>% distinct(det_plate) %>% .$det_plate %>% sort(), 
    #         summarized_temp_veh_d  %>% ungroup() %>% distinct(det_plate) %>% .$det_plate %>% sort())
    
    # nrow(temp_pert_d_1) + nrow(temp_veh_d) == nrow(analysis_dat$data[[1]])
    # TRUE
    
    cmbd <- left_join(temp_pert_d_1, summarized_temp_veh_d, 
                      by = c("cell_id", "det_plate", "pr_gene_symbol")) %>%
      mutate(value = pert_value - median_cell_id_plate_dmso_value, .before = 4) %>%
      rename(abbrev_replicate_id = pert_abbrev_replicate_id) %>%
      dplyr::select(master_id, replicate_id, abbrev_replicate_id, 
                    pert_value, median_cell_id_plate_dmso_value, value, everything()); cmbd
    
    # all.equal(nrow(cmbd), nrow(temp_pert_d_1))
    # TRUE
    return(cmbd)
  })) %>%
  mutate(unsubtracted_data = data, 
         data = dmso_subtracted_data) %>%
  unnest(filter_vars) %>%
  filter(filter_vars != "DMSO")


hclust_methods_ <- c("euclidean", "pearson", "spearman")

# correlation
corr_results <- apply(new_analysis_dat, MARGIN = 1, function(args){
  #  args <- new_analysis_dat[1,]
  # args <- analysis_dat[2,]
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
  
  # now do this with connectivity values
  lst_dat <- dat %>% split(.$pert_iname)
  cor_dat_lst <- lapply(X = lst_dat, FUN = function(l){
    cor_dat <- run_corr(l)$matrix
    return(cor_dat)
  })
  
  
  lst_ht <- lapply(1:length(cor_dat_lst), FUN = function(i){
    
    cor_dat <- cor_dat_lst[[i]]
    pert_arg_ <- names(cor_dat_lst)[[i]]
    print(pert_arg_)
    
    unit_height_width <- 7.5
    ht <- make_easy_heatmap_no_clust(cor_dat, 
                                     pert_arg = pert_arg_, 
                                     filter_arg = filter_id, 
                                     height_unit = unit_height_width, 
                                     width_unit = unit_height_width,
                                     # cluster_arg = FALSE,
                                     dataset_arg = dataset_type)
    
    ht_1_clust_lst <- lapply(1:length(hclust_methods_), FUN = function(j){
      ht_1_clust <- make_easy_heatmap(m = cor_dat, 
                                      pert_arg = pert_arg_, 
                                      filter_arg = filter_id, 
                                      height_unit = unit_height_width, 
                                      width_unit = unit_height_width,
                                      use_cell_fun = TRUE, 
                                      clustering_method = hclust_methods_[j],
                                      dataset_arg = dataset_type, 
                                      cut_tree_thresh = dendro_cut_thresh)
      return(ht_1_clust)
    }) %>%
      set_names(hclust_methods_)
    
    ht_1_clust_no_cell_lst <- lapply(1:length(hclust_methods_), FUN = function(j){
      ht_1_clust_no_cell <- make_easy_heatmap(m = cor_dat, 
                                              pert_arg = pert_arg_, 
                                              filter_arg = filter_id, 
                                              height_unit = unit_height_width, 
                                              width_unit = unit_height_width,
                                              use_cell_fun = FALSE, 
                                              clustering_method = hclust_methods_[j],
                                              dataset_arg = dataset_type, 
                                              cut_tree_thresh = dendro_cut_thresh)
      return(ht_1_clust_no_cell)
    }) %>%
      set_names(hclust_methods_)
    
    names_df <- tibble(rnames = rownames(cor_dat), cnames = colnames(cor_dat)) %>%
      mutate(r_cell_types = str_split(rnames, pattern = "--", simplify = TRUE)[,1], 
             unique_r_cell_types = make.unique(r_cell_types),
             c_cell_types = str_split(cnames, pattern = "--", simplify = TRUE)[,1],
             unique_c_cell_types = make.unique(r_cell_types)) 
    
    fx_cor_dat <- cor_dat 
    rownames(fx_cor_dat) <- names_df$unique_r_cell_types
    colnames(fx_cor_dat) <- names_df$unique_c_cell_types
    
    long_fx_cor_dat <- fx_cor_dat %>%
      as.data.frame() %>%
      rownames_to_column("rownames") %>%
      as_tibble() %>%
      pivot_longer(cols = names_df$unique_c_cell_types, names_to = "colnames") %>%
      mutate(r_cell_types = str_split(rownames, pattern = "\\.", simplify = TRUE)[,1],
             c_cell_types = str_split(colnames, pattern = "\\.", simplify = TRUE)[,1]) %>%
      group_by(r_cell_types, c_cell_types); long_fx_cor_dat
    
    summarized_long_fx_cor_dat <- long_fx_cor_dat %>%
      summarize(mean_corr = mean(value), .groups = "keep"); summarized_long_fx_cor_dat
    
    summarized_cor_dat <- summarized_long_fx_cor_dat %>%
      pivot_wider(id_cols = r_cell_types, names_from = c_cell_types, 
                  values_from = mean_corr) %>%
      column_to_rownames("r_cell_types") %>%
      as.matrix()
    
    ht_2_clust_lst <- lapply(1:length(hclust_methods_), FUN = function(k){
      ht_2_clust <- make_easy_heatmap(summarized_cor_dat, 
                                      pert_arg = pert_arg_, 
                                      filter_arg = filter_id, 
                                      dataset_arg = dataset_type, 
                                      clustering_method = hclust_methods_[k],
                                      height_unit = 10, 
                                      width_unit = 10, 
                                      use_cell_fun = TRUE, 
                                      cut_tree_thresh = dendro_cut_thresh)
      return(ht_2_clust)
    }) %>%
      set_names(hclust_methods_)
    
    res <- tibble(dataset_type, filter_id, cor_dat = list(cor_dat), 
                  summarized_cor_dat = list(summarized_cor_dat),
                  ht = list(ht), clust_methods = list(names(ht_1_clust_lst)), 
                  ht_1_clust = list(ht_1_clust_lst), 
                  ht_1_clust_no_cell = list(ht_1_clust_no_cell_lst),
                  ht_2_clust = list(ht_2_clust_lst))
    return(res)
  }) %>% 
    bind_rows() %>%
    mutate(pert_iname = names(cor_dat_lst), .before = 1) %>%
    unnest(c(clust_methods, ht_1_clust, ht_1_clust_no_cell, ht_2_clust))
  
  return(lst_ht)
}) %>%
  bind_rows()

corr_results_long <- corr_results %>% 
  pivot_longer(cols = c(ht, ht_1_clust, ht_1_clust_no_cell, ht_2_clust), names_to = "ht_type", values_to = "ht") %>%
  mutate(ht_type = ifelse(ht_type == "ht_2_clust", "clustered-summarized", 
                          ifelse(ht_type == "ht_1_clust", "clustered-full", 
                                 ifelse(ht_type == "ht_1_clust_no_cell", "clustered-no_cell-full", "full"))))


curr_plot_long_df <- corr_results_long %>% 
  filter(ht_type == "clustered-full" | ht_type == "clustered-no_cell-full")

write_plot_to_file(plot_df = curr_plot_long_df, top_directory = "corr_ht")

# write_plot_to_file(corr_results_long %>% filter(ht_type == "full"))
# write_plot_to_file(corr_results_long %>% filter(ht_type == "clustered-summarized"))



# regular corr
combined_ki_corr <- get_combined_matrix(single_res_df = corr_results, 
                                        filter_mark = "Kinase inhibitor")
combined_ki_corr_ht <- make_easy_heatmap(combined_ki_corr, 
                                         pert_arg = "Kinase inhibitor", 
                                         filter_arg = "Kinase inhibitor", 
                                         dataset_arg = "P100", 
                                         height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_ki_corr_ht
combined_ki_corr_no_clust_ht <- make_easy_heatmap(combined_ki_corr, 
                                                  pert_arg = "Kinase inhibitor", 
                                                  filter_arg = "Kinase inhibitor", 
                                                  dataset_arg = "P100", 
                                                  height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_ki_corr_no_clust_ht

# summarized KI corr
combined_ki_summarized_corr <- get_combined_matrix(single_res_df = corr_results, 
                                                   filter_mark = "Kinase inhibitor", 
                                                   column_to_select = "summarized_cor_dat")
combined_ki_summarized_corr_ht <- make_easy_heatmap(combined_ki_summarized_corr, 
                                                    pert_arg = "Kinase inhibitor", 
                                                    filter_arg = "Kinase inhibitor", 
                                                    dataset_arg = "P100", 
                                                    height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_ki_summarized_corr_ht
combined_ki_summarized_corr_no_clust_ht <- make_easy_heatmap(combined_ki_summarized_corr, 
                                                             pert_arg = "Kinase inhibitor", 
                                                             filter_arg = "Kinase inhibitor", 
                                                             dataset_arg = "P100", 
                                                             height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_ki_summarized_corr_no_clust_ht

# regular epi corr
combined_epi_corr <- get_combined_matrix(single_res_df = corr_results, filter_mark = "Epigenetic")
combined_epi_corr_ht <- make_easy_heatmap(combined_epi_corr, 
                                          pert_arg = "Epigenetic", 
                                          filter_arg = "Epigenetic", 
                                          dataset_arg = "GCP", 
                                          height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_epi_corr_ht
combined_epi_corr_no_clust_ht <- make_easy_heatmap(combined_epi_corr, 
                                                   pert_arg = "Epigenetic", 
                                                   filter_arg = "Epigenetic", 
                                                   dataset_arg = "GCP", 
                                                   height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_epi_corr_no_clust_ht

# summarized epi corr
combined_epi_summarized_corr <- get_combined_matrix(single_res_df = corr_results, 
                                                    filter_mark = "Epigenetic", 
                                                    column_to_select = "summarized_cor_dat")
combined_epi_summarized_corr_ht <- make_easy_heatmap(combined_epi_summarized_corr, 
                                                     pert_arg = "Epigenetic", 
                                                     filter_arg = "Epigenetic", 
                                                     dataset_arg = "GCP", 
                                                     height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_epi_summarized_corr_ht
combined_epi_summarized_corr_no_clust_ht <- make_easy_heatmap(combined_epi_summarized_corr, 
                                                              pert_arg = "Epigenetic", 
                                                              filter_arg = "Epigenetic", 
                                                              dataset_arg = "GCP", 
                                                              height_unit = 3.5, width_unit = 2, use_cell_fun = FALSE); # combined_epi_summarized_corr_no_clust_ht

# put corr dat into a tribble and plot it
global_corr_res <- tribble(
  ~dataset_type, ~pert_iname, ~filter_id, ~conn_dat, ~ht_type, ~ht,
  "P100", "Kinase inhibitor", "Kinase inhibitor", combined_ki_corr, "agg", combined_ki_corr_ht,
  "P100", "Kinase inhibitor", "Kinase inhibitor", combined_ki_corr, "agg-no_clust", combined_ki_corr_no_clust_ht,
  "P100", "Kinase inhibitor", "Kinase inhibitor", combined_ki_summarized_corr, "summary-agg", combined_ki_summarized_corr_ht,
  "P100", "Kinase inhibitor", "Kinase inhibitor", combined_ki_summarized_corr, "summary-agg-no_clust", combined_ki_summarized_corr_no_clust_ht,
  "GCP", "Epigenetic", "Epigenetic", combined_epi_corr, "agg", combined_epi_corr_ht,
  "GCP", "Epigenetic", "Epigenetic", combined_epi_corr, "agg-no_clust", combined_epi_corr_no_clust_ht,
  "GCP", "Epigenetic", "Epigenetic", combined_epi_summarized_corr, "summary-agg", combined_epi_summarized_corr_ht,
  "GCP", "Epigenetic", "Epigenetic", combined_epi_summarized_corr, "summary-agg-no_clust", combined_epi_summarized_corr_no_clust_ht,
)
write_plot_to_file(global_corr_res, top_directory = "global_corr_ht")


# connectivity !
conn_results <- apply(corr_results, MARGIN = 1, FUN = function(args){
  # args <- corr_results[9,]
  dataset_type <- args$dataset_type
  filter_id <- args$filter_id
  cor_dat <- args$cor_dat # cor_dat <- args$cor_dat[[1]]
  pert_arg_ <- args$pert_iname
  
  print(pert_arg_)
  
  unique_cs <- unique(rownames(cor_dat))
  grp_names <- get_grp_names_from_matrix(cor_dat, unique_cs, sep_ = "--")
  
  cs_idx <- seq_len(length(rownames(cor_dat)))
  # collapse replicates by drug - cell pair, not just cell
  cs_tib <- tibble(unique_cs, cs_idx, grp_names)
  
  # to do it like i did in the original analysis
  looper <- expand.grid(
    group_a = unique(grp_names), 
    group_b = unique(grp_names)
  ) %>%
    as_tibble()
  
  # conn_progress <- progressr::progressor(steps = nrow(looper))
  tic()
  conn_df <- map_dfr(1:nrow(looper), .f = function(z){
    # message(looper$rep_a[z], "  ", looper$rep_b[z])
    res <- calc_conn(a = looper$group_a[z], b = looper$group_b[z],
                     my_mat = cor_dat, my_tib = cs_tib,
                     use_bootstrap = F)
    # conn_progress()
    return(res)
  }) %>%
    bind_rows() 
  toc()
  
  replicate_matrix <- conn_df %>%
    dplyr::select(group_a, group_b, conn) %>%
    pivot_wider(id_cols = group_a,
                names_from = group_b,
                values_from = conn,
                values_fn = function(x) median(x, na.rm = TRUE))
  
  numeric_replicate_matrix <- replicate_matrix %>%
    column_to_rownames("group_a") %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()
  
  # should I do this?
  symmetric_numeric_replicate_matrix <- calc_conn_by_median(m = numeric_replicate_matrix)
  
  ht <- make_easy_heatmap(symmetric_numeric_replicate_matrix, 
                          pert_arg = pert_arg_, 
                          filter_arg = filter_id, 
                          dataset_arg = dataset_type, cluster_arg = TRUE, 
                          height_unit = 10, width_unit = 10, use_cell_fun = TRUE); ht
  
  res <- tibble(dataset_type, pert_iname = pert_arg_, filter_id, conn_dat = list(symmetric_numeric_replicate_matrix), ht = list(ht))
  return(res)
}) %>%
  bind_rows()

conn_results_long <- conn_results %>%
  pivot_longer(cols = ht, names_to = "ht_type", values_to = "ht") %>%
  mutate(ht_type = ifelse(ht_type == "ht_2", "summarized", "full"))

write_plot_to_file(conn_results_long, top_directory = "conn_ht") 

# Combined ki Conn
combined_ki_conn <- get_combined_matrix(single_res_df = conn_results, filter_mark = "Kinase inhibitor", column_to_select = "conn_dat")
combined_ki_conn_ht <- make_easy_heatmap(combined_ki_conn, 
                                         pert_arg = "Kinase inhibitor", 
                                         filter_arg = "Kinase inhibitor", 
                                         dataset_arg = "P100", cluster_arg = TRUE, 
                                         height_unit = 8, width_unit = 6.5, use_cell_fun = TRUE); combined_ki_conn_ht
combined_ki_conn_no_clust_ht <- make_easy_heatmap(combined_ki_conn, 
                                                  pert_arg = "Kinase inhibitor", 
                                                  filter_arg = "Kinase inhibitor", 
                                                  dataset_arg = "P100", cluster_arg = FALSE, 
                                                  height_unit = 8, width_unit = 6.5, use_cell_fun = TRUE); combined_ki_conn_no_clust_ht
# Combined epi Conn
combined_epi_conn <- get_combined_matrix(single_res_df = conn_results, filter_mark = "Epigenetic", column_to_select = "conn_dat")
combined_epi_conn_ht <- make_easy_heatmap(combined_epi_conn, 
                                          pert_arg = "Epigenetic", 
                                          filter_arg = "Epigenetic", 
                                          dataset_arg = "GCP", cluster_arg = TRUE, 
                                          height_unit = 6.5, width_unit = 5, use_cell_fun = TRUE); combined_epi_conn_ht
combined_epi_conn_no_clust_ht <- make_easy_heatmap(combined_epi_conn, 
                                                   pert_arg = "Epigenetic", 
                                                   filter_arg = "Epigenetic", 
                                                   dataset_arg = "GCP", cluster_arg = FALSE, 
                                                   height_unit = 6.5, width_unit = 5, use_cell_fun = TRUE); combined_epi_conn_no_clust_ht

# put combined conn into tribble and plot
global_conn_res <- tribble(
  ~dataset_type, ~pert_iname, ~filter_id, ~conn_dat, ~ht_type, ~ht,
  "P100", "Kinase inhibitor", "Kinase inhibitor", combined_ki_conn, "agg", combined_ki_conn_ht,
  "P100", "Kinase inhibitor", "Kinase inhibitor", combined_ki_conn, "agg-no_clust", combined_ki_conn_no_clust_ht,
  "GCP", "Epigenetic", "Epigenetic", combined_epi_conn, "agg", combined_epi_conn_ht,
  "GCP", "Epigenetic", "Epigenetic", combined_epi_conn, "agg-no_clust", combined_epi_conn_no_clust_ht
)
write_plot_to_file(global_conn_res, top_directory = "global_conn_ht")


# library(missMDA)
# library(FactoMineR)
# 
# # estimate the number of components from incomplete data
# combined_conn_pca <- apply(global_conn_res, MARGIN = 1, FUN = function(args){
#   tic()
#   title_ <- str_c(args$pert_iname, args$datset_type, sep = "-")
#   message(title_)
#   mat <- args %>% pluck("conn_dat")
#   nb <- missMDA::estim_ncpPCA(mat, method.cv = "Kfold", verbose = FALSE) 
#   # iterativePCA algorithm
#   res.comp <- missMDA::imputePCA(mat, ncp = nb$ncp) 
#   imp <- res.comp$completeObs
#   res.pca <- FactoMineR::PCA(imp, ncp = nb$ncp, graph=FALSE) # quanti.sup = 1, quali.sup = 12, 
#   plot(res.pca, title = title_); # hab=12, lab="quali"
#   toc()
#   return(res.pca)
# })


# write_rds(conn_results, file = file.path(supplemental_output_dir, "conn_results-cell_cell.rds"))
# conn_res <- read_rds(file.path(supplemental_output_dir, "conn_results-cell_cell.rds"))
# 
# n_neighbors_tib <- tribble(~dataset_type, ~filter_id, ~n_neighbors,
#                            "P100", "Kinase inhibitor", 4,
#                            "P100", "DMSO", 4,
#                            "GCP", "Epigenetic", 4,
#                            "GCP", "DMSO", 4)
# custom.config <- umap.defaults
# updated_conn_res <- left_join(conn_res, n_neighbors_tib, by = c("filter_id", "dataset_type"))
# 
# umap_uwot_pca_results <- apply(updated_conn_res, MARGIN = 1, FUN = function(args){
#   n_neighbors_uwot <- args$n_neighbors
#   custom.config$n_neighbors <- n_neighbors_uwot
#   
#   dat <- args$dat # dat <- args$dat[[1]]
#   grouping_id <- unique(args$grouping_id)
#   filter_id <- unique(args$filter_id)
#   dataset_type <- unique(args$dataset_type)
#   
#   replicate_matrix_dat <- dat %>%
#     dplyr::select(replicate_id, value, pr_gene_symbol) %>%
#     pivot_wider(id_cols = c(replicate_id),
#                 names_from = pr_gene_symbol,
#                 values_from = value,
#                 values_fn = function(x) median(x, na.rm = TRUE)) %>%
#     column_to_rownames("replicate_id")
#   # UMAP on raw data
#   
#   dat_distance <- as.matrix(dist(replicate_matrix_dat))
#   umap_fit <- umap::umap(dat_distance, input = "dist", config=custom.config)
#   
#   umap_tbl <- tibble(UMAP_1 = umap_fit$layout[,1],
#                      UMAP_2 = umap_fit$layout[,2],
#                      replicate_id = rownames(dat_distance),
#                      cell_id = str_split(replicate_id, "--", simplify = T)[,1]) %>%
#     # pert_iname = str_split(replicate_id, "--", simplify = T)[,2]) %>%
#     left_join(informative_lbls, by= c("cell_id" = "lbs")) %>%
#     mutate(cell_id = factor(cell_id)) %>%
#     mutate(filter_id, grouping_id, dataset_type)
#   
#   numeric_replicate_matrix <- args$numeric_replicate_matrix # numeric_replicate_matrix <- args$numeric_replicate_matrix[[1]]
#   # perform UMAP and PCA
#   
#   umap_conn_replicates_fit <- umap::umap(numeric_replicate_matrix, config=custom.config)
#   res_pca <- prcomp(numeric_replicate_matrix, scale. = F, center = F)
#   
#   scaled_uwot_conn_replicates_fit <- uwot::umap(X = numeric_replicate_matrix,
#                                                 n_neighbors = n_neighbors_uwot, min_dist = 0.001,
#                                                 # metric = "correlation",
#                                                 verbose = TRUE, n_threads = no_cores) # no_cores defined in load.R
#   
#   
#   scaled_numeric_replicate_matrix <- args$scaled_numeric_replicate_matrix # scaled_numeric_replicate_matrix <- args$scaled_numeric_replicate_matrix[[1]]
#   # perform UMAP and PCA on SCALED data
#   scaled_umap_conn_replicates_fit <- umap::umap(scaled_numeric_replicate_matrix, config=custom.config)
#   scaled_res_pca <- prcomp(scaled_numeric_replicate_matrix, scale. = F, center = F)
#   
#   # stop()
#   conn_replicates_tbl <- tibble(UMAP_1 = umap_conn_replicates_fit$layout[,1],
#                                 UMAP_2 = umap_conn_replicates_fit$layout[,2],
#                                 
#                                 PCA_1 = res_pca$x[,1],
#                                 PCA_2 = res_pca$x[,2],
#                                 
#                                 UWOT_1 = scaled_uwot_conn_replicates_fit[,1],
#                                 UWOT_2 = scaled_uwot_conn_replicates_fit[,2],
#                                 
#                                 UMAP_1_scaled = scaled_umap_conn_replicates_fit$layout[,1],
#                                 UMAP_2_scaled = scaled_umap_conn_replicates_fit$layout[,2],
#                                 
#                                 PCA_1_scaled = scaled_res_pca$x[,1],
#                                 PCA_2_scaled = scaled_res_pca$x[,2],
#                                 
#                                 replicate_id = rownames(numeric_replicate_matrix),
#                                 cell_id = str_split(replicate_id, "--", simplify = T)[,1]) %>%
#     # pert_iname = str_split(replicate_id, "--", simplify = T)[,2]) %>%
#     left_join(informative_lbls, by= c("cell_id" = "lbs")) %>%
#     mutate(cell_id = factor(cell_id)) %>%
#     mutate(filter_id, grouping_id, dataset_type)
#   
#   final_res <- tibble(grouping_id = grouping_id, filter_id = filter_id,
#                       n_neighbors = n_neighbors_uwot,
#                       dataset_type = dataset_type,
#                       matrix_raw_dat = list(replicate_matrix_dat),
#                       umap_raw = list(umap_tbl),
#                       umap_uwot_pca_connectivity = list(conn_replicates_tbl),
#                       umap_obj = list(umap_conn_replicates_fit),
#                       prcomp_obj = list(res_pca),
#                       scaled_uwot_obj = list(scaled_uwot_conn_replicates_fit), 
#                       umap_scaled_obj = list(scaled_umap_conn_replicates_fit),
#                       prcomp_obj_scaled = list(scaled_res_pca))
#   return(final_res)
# }) %>%
#   bind_rows()
# 
# cmbd_dat_and_res <- left_join(conn_res, umap_uwot_pca_results, by = c("grouping_id", "filter_id", "dataset_type"))
# write_rds(cmbd_dat_and_res, file = file.path(supplemental_output_dir, "umap_uwot_pca_results-cell_cell.rds"))
# 
# umap_uwot_pca_data_and_results <- read_rds(file.path(supplemental_output_dir, "umap_uwot_pca_results-cell_cell.rds"))
# 
# pca_plots <- apply(umap_uwot_pca_data_and_results, MARGIN = 1, FUN = function(args){
#   
#   # args <- umap_uwot_pca_data_and_results[1,]
#   filter_id <- unique(args$filter_id)
#   grouping_id <- unique(args$grouping_id)
#   dataset_type <- unique(args$dataset_type)
#   
#   numeric_replicate_matrix <- args$numeric_replicate_matrix # numeric_replicate_matrix <- args$numeric_replicate_matrix[[1]]
#   ut_cells <- args$conn_df # ut_cells <- ut_cells[[1]]
#   ut <- args$umap_uwot_pca_connectivity # ut <- args$umap_uwot_pca_connectivity[[1]]
#   prcomp_obj <- args$prcomp_obj # prcomp_obj <- args$prcomp_obj[[1]]
#   
#   cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
#   ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
#   color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
#   # pert_color_map <- deframe(pert_color_map_init %>% filter(pert_iname %in% ut$pert_iname))
#   
#   explained_var <- round(summary(prcomp_obj)$importance[2,c(1,2)]*100, 2)
#   
#   by_cell <- ggplot(ut, aes(PCA_1, PCA_2)) +
#     geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 5,alpha = 0.75) + 
#     theme_bw() +
#     scale_fill_manual(values = color_map) +
#     labs(fill = "Cell ID") +
#     ggtitle(qq("PCA on @{filter_id} | @{dataset_type}")) +
#     labs(caption = qq("@{dataset_type} connectivity profiles\nDefault PCA parameters")) +
#     xlab(qq("standarized PC1 (@{explained_var[1]}% explained var.)")) +
#     ylab(qq("standarized PC2 (@{explained_var[2]}% explained var.)")); by_cell
#   
#   pca_proper_plot <- fviz_pca_ind(prcomp_obj,
#                                   col.ind = "cos2", # Color by the quality of representation
#                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                                   repel = TRUE     # Avoid text overlapping
#   ) + 
#     ggtitle(qq("PCA on @{filter_id} | @{dataset_type}")) +
#     labs(caption = qq("@{dataset_type} connectivity profiles\nDefault PCA parameters"))
#   
#   # by_pert <- ggplot(ut, aes(PCA_1, PCA_2)) +
#   #   geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3,alpha = 0.75) + 
#   #   theme_bw() +
#   #   scale_fill_manual(values = pert_color_map) +
#   #   labs(fill = "Perturbation") +
#   #   ggtitle(qq("PCA on @{filter_id} | @{dataset_type}")) +
#   #   labs(caption = qq("@{dataset_type} connectivity profiles\nDefault PCA parameters")) +
#   #   xlab(qq("standarized PC1 (@{explained_var[1]}% explained var.)")) +
#   #   ylab(qq("standarized PC2 (@{explained_var[2]}% explained var.)")); by_pert
#   # 
#   
#   return(tibble(filter_id, grouping_id, dataset_type) %>% mutate(cell_g = list(by_cell), pca_proper = list(pca_proper_plot))) # , pert_g = list(by_pert)
# }) %>%
#   bind_rows() %>%
#   pivot_longer(cols = cell_g:pca_proper, names_to = "plot_type") %>% # :pert_g
#   mutate(dir_type = "pca",
#          plot_dir = file.path(supplemental_output_dir, dir_type),
#          fn_name = file.path(plot_dir, str_c(str_c(filter_id, dataset_type, plot_type, sep = "-"),".png", sep = "")))
# # (pca_plots$cell_g[[1]] | pca_plots$cell_g[[2]]) /
# #   (pca_plots$cell_g[[3]] | pca_plots$cell_g[[4]])
# 
# # (pca_plots$cell_g[[1]] | pca_plots$pert_g[[1]])
# # (pca_plots$cell_g[[3]] | pca_plots$pert_g[[3]])
# 
# umap_plots <- apply(umap_uwot_pca_data_and_results, MARGIN = 1, FUN = function(args){
#   # args <- umap_uwot_pca_data_and_results[1,]
#   filter_id <- unique(args$filter_id)
#   grouping_id <- unique(args$grouping_id)
#   dataset_type <- unique(args$dataset_type)
#   
#   ut_cells <- args$conn_df # ut_cells <- args$conn_df[[1]]
#   ut <- args$umap_uwot_pca_connectivity #  ut <- args$umap_uwot_pca_connectivity[[1]]
#   n_neighbors <- args$n_neighbors
#   
#   cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
#   ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
#   color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
#   # pert_color_map <- deframe(pert_color_map_init %>% filter(pert_iname %in% ut$pert_iname))
#   
#   by_cell <- ggplot(ut, aes(UMAP_1, UMAP_2)) +
#     geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 3,alpha = 0.75) + 
#     theme_bw() +
#     scale_fill_manual(values = color_map) +
#     labs(fill = "Cell ID") +
#     ggtitle(qq("UMAP on @{filter_id} | @{dataset_type}")) +
#     labs(caption = qq("Euclidean distance on @{dataset_type} connectivity profiles\nDefault UMAP parameters\nn_neighbors = @{n_neighbors}"))
#   
#   # by_pert <- ggplot(ut, aes(UMAP_1, UMAP_2)) +
#   #   geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3,alpha = 0.75) +
#   #   theme_bw() +
#   #   scale_fill_manual(values = pert_color_map) +
#   #   labs(fill = "Perturbation") +
#   #   ggtitle(qq("UMAP on @{filter_id} | @{dataset_type}")) +
#   #   labs(caption = qq("Euclidean distance on @{dataset_type} connectivity profiles\nDefault UMAP parameters\nn_neighbors = @{n_neighbors}"))
#   # 
#   return(tibble(filter_id, grouping_id, dataset_type) %>% mutate(cell_g = list(by_cell))) # , pert_g = list(by_pert)
# }) %>%
#   bind_rows()%>%
#   pivot_longer(cols = cell_g, names_to = "plot_type") %>% #:pert_g
#   mutate(dir_type = "umap",
#          plot_dir = file.path(supplemental_output_dir, dir_type),
#          fn_name = file.path(plot_dir, str_c(str_c(filter_id, dataset_type, plot_type, sep = "-"),".png", sep = "")))
# 
# uwot_plots <- apply(umap_uwot_pca_data_and_results, MARGIN = 1, FUN = function(args){
#   filter_id <- unique(args$filter_id)
#   grouping_id <- unique(args$grouping_id)
#   dataset_type <- unique(args$dataset_type)
#   
#   ut_cells <- args$conn_df
#   ut <- args$umap_uwot_pca_connectivity
#   n_neighbors <- args$n_neighbors
#   
#   cell_ids_present <- str_split(ut_cells[,1, drop = TRUE] %>% unique(), "--", simplify = TRUE)[,1]
#   ilbs <- informative_lbls %>% filter(lbs %in% cell_ids_present)
#   color_map <- deframe(ilbs %>% dplyr::select(lbs, cell_individual_color)) 
#   # pert_color_map <- deframe(pert_color_map_init %>% filter(pert_iname %in% ut$pert_iname))
#   
#   by_cell <- ggplot(ut, aes(UWOT_1, UWOT_2)) +
#     geom_point(aes(fill = cell_id), pch = 21, color = "black", size = 3,alpha = 0.75) + 
#     theme_bw() +
#     scale_fill_manual(values = color_map) +
#     labs(fill = "Cell ID") +
#     ggtitle(qq("UWOT on @{filter_id} | @{dataset_type}")) +
#     labs(caption = qq("Correlation distance on @{dataset_type} connectivity profiles\nDefault UWOT parameters\nn_neighbors = @{n_neighbors}"))
#   
#   # by_pert <- ggplot(ut, aes(UWOT_1, UWOT_2)) +
#   #   geom_point(aes(fill = pert_iname), pch = 21, color = "black", size = 3,alpha = 0.75) +
#   #   theme_bw() +
#   #   scale_fill_manual(values = pert_color_map) +
#   #   labs(fill = "Perturbation") +
#   #   ggtitle(qq("UWOT on @{filter_id} | @{dataset_type}")) +
#   #   labs(caption = qq("Correlation distance on @{dataset_type} connectivity profiles\nDefault UWOT parameters\nn_neighbors = @{n_neighbors}"))
#   return(tibble(filter_id, grouping_id, dataset_type) %>% mutate(cell_g = list(by_cell))) # , pert_g = list(by_pert)
# }) %>%
#   bind_rows() %>%
#   pivot_longer(cols = cell_g, names_to = "plot_type") %>% # pert_g
#   mutate(dir_type = "uwot",
#          plot_dir = file.path(supplemental_output_dir, dir_type),
#          fn_name = file.path(plot_dir, str_c(str_c(filter_id, dataset_type, plot_type, sep = "-"),".png", sep = "")))
# 
# dir.create(unique(pca_plots$plot_dir), recursive = T, showWarnings = F)
# dir.create(unique(umap_plots$plot_dir), recursive = T, showWarnings = F)
# dir.create(unique(uwot_plots$plot_dir), recursive = T, showWarnings = F)
# 
# 
# walk2(.x = as.list(pca_plots$fn_name), .y = pca_plots$value, .f = function(fn, plot_) ggsave(filename = fn, plot = plot_, device = "png"))
# walk2(.x = as.list(umap_plots$fn_name), .y = umap_plots$value, .f = function(fn, plot_) ggsave(filename = fn, plot = plot_, device = "png"))
# walk2(.x = as.list(uwot_plots$fn_name), .y = uwot_plots$value, .f = function(fn, plot_) ggsave(filename = fn, plot = plot_, device = "png"))


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