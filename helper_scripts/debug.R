rm(list = ls())
library(morpheus) # heat maps
library(tidyverse)
library(readxl)
library(purrr)
library(amap) # clustering of correlations according to Lev's paper
library(pvclust) # cluster stability, Lev's paper
library(dendextend)
library(viridis)
library(GetoptLong)
library(gridExtra)
library(ggplotify)
library(gridGraphics)
library(ggrepel)
library(ggpubr)
library(doMC)
library(conflicted)
library(showtext)
# devtools::install_github('cmap/morpheus.R')

set.seed(42) 
ncores <- 6
registerDoMC(cores=ncores)

conflict_prefer("rbind", "BiocGenerics")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("summarize","dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("union", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("intersect","dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("qq","GetoptLong")
conflict_prefer("which", "Matrix")
conflict_prefer("map","purrr")
conflict_prefer("layout", "plotly")
conflict_prefer("flatten","purrr")

font_add_google("Noto Sans")

## IMPORTANT DIRECTORIES ##
WORKING_DIRECTORY <- "/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws"
OUTPUT_DIRECTORY <- file.path(WORKING_DIRECTORY, "output")

CARDONC_DIRECTORY <- file.path("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/Cardio-oncology")
REFERENCES_DIRECTORY <- file.path(CARDONC_DIRECTORY, "References")
DATASETS_DIRECTORY <- file.path(CARDONC_DIRECTORY, "Datasets", "1st gen data")

## load helper functions ##
source(file.path(WORKING_DIRECTORY, "helper_scripts", "differential_analyte_helper_functions.R"))
source(file.path(WORKING_DIRECTORY, "helper_scripts", "connectivity_and_clustering_helper_functions.R"))
source(file.path(WORKING_DIRECTORY, "helper_scripts", "data_wrangling_helper_functions.R"))
source(file.path(WORKING_DIRECTORY, "helper_scripts", "morpheus_helper_functions.R"))



###############################################################################################################################################################
###################################################### LOAD and RESHAPE DATA ###############################################################################
###############################################################################################################################################################

#' @note first, go to https://clue.io/cmapPy/build.html#install and install cmapPy to manipulate gct's
#' @note then check out this tutorial from Oana: https://github.com/cmap/cmapPy/blob/master/tutorials/cmapPy_pandasGEXpress_tutorial.ipynb
gcp_a <- load_gcp(first_gen_dataset = file.path(DATASETS_DIRECTORY, "GCP","GCP All Cell Lines.gct"))
gcp_obj <- load_data(data = gcp_a, drug_classes_fn = file.path(REFERENCES_DIRECTORY, "Drug Glossary_edited.xlsx"), dataset="gcp")
gcp_c <- gcp_obj$data
gcp_feature_set <- gcp_obj$feature_set
gcp_perturbation_feature_set <- gcp_obj$perturbation_feature_set


p100a <- load_p100(first_gen_dataset = file.path(DATASETS_DIRECTORY, "P100","P100 All Cell Lines.gct"))
p100_obj <- load_data(data = p100a, drug_classes_fn = file.path(REFERENCES_DIRECTORY, "Drug Glossary_edited.xlsx"), dataset="p100")
p100c <- p100_obj$data
p100_feature_set <- p100_obj$feature_set
p100_perturbation_feature_set <- p100_obj$perturbation_feature_set

common_perts <- intersect(p100_obj$perturbation_feature_set, gcp_obj$perturbation_feature_set)
common_moas <- intersect(p100_obj$data_with_reps$drug_class, gcp_obj$data_with_reps$drug_class)

# p100_key <- read_tsv("~/Downloads/master_tsv-p100.tsv") %>% mutate(dataset = "p100")
# gcp_key <- read_tsv("~/Downloads/master_tsv-gcp.tsv") %>% mutate(dataset = "gcp")
# MASTER_KEY <- bind_rows(p100_key, gcp_key) %>%
#   distinct(well, pr_gene_symbol, phosphosite, cell_id, pert_iname, cell_cat, dataset)

message(qq("Total number of drugs used: @{length(common_perts)} "))
message(qq("Number of P100 analytes: @{length(p100_feature_set)} "))
message(qq("Number of GCP analytes: @{length(gcp_feature_set)} "))


old_gcp_moa <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/old_output/gcp/moa/gcp_diff_ex-h_0.6-q_0.1.rds")
new_gcp_moa <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/gcp/moa/diffex/cancer_vs_noncancer/gcp_diff_ex-h_0.6-q_0.1.rds")

  new_gcp_moa; old_gcp_moa
  
new_other <- new_gcp_moa$df[[2]]$well
old_other <- old_gcp_moa$df[[2]]$well

heatmap(new_gcp_moa$cell_corr_r_mat[[1]])
heatmap(old_gcp_moa$cell_corr_nr[[1]])
wells_in_old_not_new <- setdiff(old_other, new_other); wells_in_old_not_new

gcp_obj$data_with_reps %>% filter(well %in% wells_in_old_not_new) 


# DEBUG: 
gcp_lst_obj <- create_obj_lst(gcp_obj, dataset_common_perts = common_perts)

base_output_dir <- file.path(OUTPUT_DIRECTORY, "gcp")
dir.create(base_output_dir, recursive = T, showWarnings = F)

obj_lst = gcp_lst_obj; dataset = "gcp"; 
AVG_TOGGLE = FALSE; RERUN = FALSE; DENDRO_CUT_THRESH = 0.6; bh_thresh_val = 0.1; 
cancer_vs_non_cancer = TRUE

i = 2# i = 1 == pert
sub_dirs <- c("pert", "moa", "all")
suffixes <- c("", "_moa", "_all")

base_output_dir_batch <- file.path("~/Downloads/test")
dir.create(base_output_dir_batch)
base_output_dir_conn_final <- file.path(base_output_dir_batch, "clustering")


base_output_dir_diffex <- file.path(base_output_dir_batch, "diffex")
if (cancer_vs_non_cancer) {
  base_output_dir_diffex_final <- file.path(base_output_dir_diffex, "cancer_vs_noncancer") # separate out
} else {
  base_output_dir_diffex_final <- file.path(base_output_dir_diffex,"complete")
}

message(base_output_dir_conn_final)
message(base_output_dir_diffex_final)


# in run connectivity  function
obj = gcp_lst_obj
# 
mat_reps <- get_replicates(obj) # dat_obj = obj
corr_mat_reps <- mat_reps %>%
  mutate(corr_r = map(t_dataframe, compute_correlation),
         corr_nr = map(t_dataframe, compute_correlation, remove = F))  %>%
  mutate(match_df = map(corr_r, match_function, from = "well", to = "cell_id", rc = c("r", "c"))) %>%
  mutate(cell_corr_r_mat = map2(corr_r, match_df, replace_rc_names))

# walk2(.x = corr_mat_reps$cell_corr_r_mat, .y = corr_mat_reps$pert_iname,
#       function(x,y) plot_morpheus(dat = x, drug_name = y, dataset=dataset, base_output_dir = base_output_dir_conn_final,
#                                   subdir = "corr_reps", minr = -0.75, maxr = 0.75, sym = TRUE))

# compute connectivity for each
# conn_reps <- corr_mat_reps %>%
#   mutate(conn = map(cell_corr_r_mat, compute_connectivity))
# 
# conn_mat_reps <- conn_reps %>%
#   mutate(conn_mat = map(conn, make_numeric_mat)) %>%
#   mutate(med_conn_mat = map(conn_mat, collapse_connectivity_by_median))
# write_rds(conn_mat_reps, qq("~/Downloads/@{dataset}-test-@{suffixes[i]}.rds"),compress = "gz")
conn_mat_reps <- read_rds( qq("~/Downloads/@{dataset}-test-@{suffixes[i]}.rds"))

# res <- compute_boot_pvclust(x = conn_mat_reps$med_conn_mat[[1]], conn_mat_reps$pert_iname[[1]])
# cut_trees <- generate_clusters(res, thresh = dend_thresh)
# co_clust <- co_cluster(cut_trees)

# conn_clust <- conn_mat_reps %>%
#   mutate(pvclust = map2(med_conn_mat, pert_iname, .f = compute_boot_pvclust))
# write_rds(conn_clust, qq("~/Downloads/@{dataset}-test-conn-@{suffixes[i]}.rds"),compress = "gz")
conn_clust <- read_rds(qq("~/Downloads/@{dataset}-test-conn-@{suffixes[i]}.rds"))


# in run differential expression function
conn_clust_obj = conn_clust

dendro_cut_thresh <- DENDRO_CUT_THRESH; bh_thresh <- bh_thresh_val

conn_clust_trees <- conn_clust_obj %>% 
  mutate(which_dat = dataset) %>%
  mutate(cut_trees = map(pvclust, .f = generate_clusters, thresh = dendro_cut_thresh)) %>%
  mutate(co_clust = map(cut_trees, co_cluster)) 

# conn_clust_trees_test <- conn_clust_trees %>% filter(pert_iname == "Epigenetic")
# conn_clust_trees_test <- conn_clust_trees %>% filter(pert_iname == "ly-294002")
conn_clust_trees_test <- conn_clust_trees %>% filter(pert_iname == "vorinostat")

j = 1;
master_key <- MASTER_KEY; clust_assignments = conn_clust_trees_test$cut_trees[[j]]; dat = conn_clust_trees_test$df[[j]]; 
dname = conn_clust_trees_test$pert_iname[[j]]; co_clust = conn_clust_trees_test$co_clust[[j]]; which.dat = dataset
cancer_vs_non_cancer = TRUE


diff_ex_obj <- conn_clust_trees_test %>%
  mutate(diff_ex = pmap(list(clust_assignments = cut_trees, dat = df, dname = pert_iname, co_clust = co_clust), 
                        .f = calc_differential_analytes, 
                        which.dat = dataset, gid = c("HUVEC","HAoSMC"),
                        cancer_vs_non_cancer = cancer_vs_non_cancer,
                        base_output_dir = base_output_dir,
                        dendro_cut_thresh = dendro_cut_thresh,
                        bh_thresh = bh_thresh, master_key = master_key))

pwalk(.l = list(diff_ex_obj$diff_ex, diff_ex_obj$pert_iname, diff_ex_obj$cut_trees, diff_ex_obj$co_clust),
      .f = plot_diff_analyte_results_CVs_html, 
      which.dat = dataset, 
      gid = c("HUVEC", "HAoSMC"), 
      cancer_vs_non_cancer = cancer_vs_non_cancer, 
      base_output_dir = base_output_dir, 
      dendro_cut_thresh = dendro_cut_thresh, bh_thresh = bh_thresh)








(x = matrix(c(-2,-1,0,1,2,1.5,2,0,1,2,NA,NA,0,1,2),5))
(x[1:3, 1:3])
(cor(x[1:3, 1:3]))

# default
(cor(x, use = "everything")) 

# expect perfect correlation since the only rows with no observations are exactly identical
(cor(x, use = "complete"))

# different parts of columns are used to compute correlation if the column has a missing value
(cor(x, use="pairwise.complete.obs"))


ht1_lst <- organize_and_plot_heatmap_subfunction(filtered_test_mat = filtered_test_mat1, 
                                                 row_annots_df = row_annots_df1, column_annots_df = column_annots_df1,
                                                 heatmap_output_fn = heatmap_output_fn1, title_var = title_var)




