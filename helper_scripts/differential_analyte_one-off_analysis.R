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


# assumes you have run cluster3.R
p100_clust <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/p100/pert/p100_conn_clust.rds")
gcp_clust <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/gcp/pert/gcp_conn_clust.rds")
avg_clust <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/avg/pert/avg_conn_clust.rds")


# extract and return the top 10 most differential analytes
get_most_diff <- function(obj, direction="++"){
  # direction needs to be ++ or --
  x <- obj$diff_ex[[1]] %>% 
    arrange(p_val) %>% 
    filter(base_clust_comp == 2, d_stat == direction) %>% 
    slice_head(n = 10)
  return(x)
}

# For the DMSO data, 
# what P100 or GCP marks differentiate vascular from cancer cells at baseline?

## p100 separate
p100_base_output_dir <- file.path(OUTPUT_DIRECTORY, "p100")
p100_cancer_vs_noncancer_output_dir <- file.path(p100_base_output_dir, "cancer_vs_noncancer", "pert")

p100_clust_dmso <- p100_clust %>% filter(pert_iname == "dmso")
p100_diffe_dmso <- run_differential_analyte_analysis(conn_clust_obj = p100_clust_dmso, dataset = "p100",base_output_dir_specific = p100_cancer_vs_noncancer_output_dir,
                                                     cancer_vs_non_cancer = TRUE)
p100_dmso_most_pos <- get_most_diff(p100_diffe_dmso, direction = "++")
p100_dmso_most_neg <- get_most_diff(p100_diffe_dmso, direction = "--")
p_p100 <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/p100/cancer_vs_noncancer/diff_exp/h_0.6/q_0.25/CO_CLUST/dmso_diff_exp__q_0.25.rds")


## gcp separate
gcp_base_output_dir <- file.path(OUTPUT_DIRECTORY, "gcp")
gcp_cancer_vs_noncancer_output_dir <- file.path(gcp_base_output_dir, "cancer_vs_noncancer", "pert")

gcp_clust_dmso <- gcp_clust %>% filter(pert_iname == "dmso")
gcp_diffe_dmso <- run_differential_analyte_analysis(conn_clust_obj = gcp_clust_dmso, dataset = "gcp",base_output_dir_specific = gcp_cancer_vs_noncancer_output_dir,
                                                     cancer_vs_non_cancer = TRUE)
gcp_dmso_most_pos <- get_most_diff(gcp_diffe_dmso, direction = "++")
gcp_dmso_most_neg <- get_most_diff(gcp_diffe_dmso, direction = "--")
p_gcp <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/gcp/cancer_vs_noncancer/diff_exp/h_0.6/q_0.25/CO_CLUST/dmso_diff_exp__q_0.25.rds")


## avg 
avg_base_output_dir <- file.path(OUTPUT_DIRECTORY, "avg")
avg_cancer_vs_noncancer_output_dir <- file.path(avg_base_output_dir, "cancer_vs_noncancer", "pert")

avg_clust_dmso <- avg_clust %>% filter(pert_iname == "dmso")
avg_diffe_dmso <- run_differential_analyte_analysis(conn_clust_obj = avg_clust_dmso, dataset = "avg",base_output_dir_specific = avg_cancer_vs_noncancer_output_dir,
                                                    cancer_vs_non_cancer = TRUE)
avg_dmso_most_pos <- get_most_diff(avg_diffe_dmso, direction = "++")
avg_dmso_most_neg <- get_most_diff(avg_diffe_dmso, direction = "--")
p_avg <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/avg/cancer_vs_noncancer/diff_exp/h_0.6/q_0.25/CO_CLUST/dmso_diff_exp__q_0.25.rds")


bind_cols("p100" = p100_dmso_most_pos$analyte, "gcp" = gcp_dmso_most_pos$analyte, "avg" = avg_dmso_most_pos$analyte)
bind_cols("p100" = p100_dmso_most_neg$analyte, "gcp" = gcp_dmso_most_neg$analyte, "avg" = avg_dmso_most_neg$analyte)

# And for the P100 phosphoproteomic response, 
# can you compare vascular cells to cancer cells and figure out which phosphorylation events 
# differentiate the 1) epigenetic modifying drug response and then 2) the kinase inhibitor response? 

# p100 epi
p100_moa_clust <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/p100/moa/p100_conn_clust_moa.rds")
p100_moa_base_output_dir <- file.path(OUTPUT_DIRECTORY, "p100")
p100_moa_epi_cancer_vs_noncancer_output_dir <- file.path(p100_moa_base_output_dir, 
                                                     "cancer_vs_noncancer", "epi")

p100_moa_clust_epi <- p100_moa_clust %>% filter(pert_iname == "Epigenetic")
p100_moa_diffe_epi <- run_differential_analyte_analysis(conn_clust_obj = p100_moa_clust_epi, dataset = "p100",
                                                    base_output_dir_specific = p100_moa_epi_cancer_vs_noncancer_output_dir,
                                                    cancer_vs_non_cancer = TRUE)
p100_moa_diffe_epi_most_pos <- get_most_diff(p100_moa_diffe_epi, direction = "++")
p100_moa_diffe_epi_most_neg <- get_most_diff(p100_moa_diffe_epi, direction = "--")
p_p100_moa_epi <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/p100/cancer_vs_noncancer/epi/diff_exp/h_0.6/q_0.25/CO_CLUST/Epigenetic_diff_exp__q_0.25.rds")

bind_cols("pos" = p100_moa_diffe_epi_most_pos$analyte, "neg" = p100_moa_diffe_epi_most_neg$analyte)

# p100 ki

p100_moa_ki_cancer_vs_noncancer_output_dir <- file.path(p100_moa_base_output_dir, 
                                                     "cancer_vs_noncancer", "ki")

p100_moa_clust_ki <- p100_moa_clust %>% filter(pert_iname == "Kinase inhibitor")
p100_moa_diffe_ki <- run_differential_analyte_analysis(conn_clust_obj = p100_moa_clust_ki, dataset = "p100",
                                                        base_output_dir_specific = p100_moa_ki_cancer_vs_noncancer_output_dir,
                                                        cancer_vs_non_cancer = TRUE)
p100_moa_diffe_ki_most_pos <- get_most_diff(p100_moa_diffe_ki, direction = "++")
p100_moa_diffe_ki_most_neg <- get_most_diff(p100_moa_diffe_ki, direction = "--")
p_p100_moa_ki <- read_rds("/Users/Nicholas/OneDrive - Tufts/phd/jaffe/workspace/ws/output/p100/cancer_vs_noncancer/ki/diff_exp/h_0.6/q_0.25/CO_CLUST/Kinase inhibitor_diff_exp__q_0.25.rds")

bind_cols("pos" = p100_moa_diffe_ki_most_pos$analyte, "neg" = p100_moa_diffe_ki_most_neg$analyte)
# I think if you did that, you might have 6-7 figures that could be a paper for you to write in the fall.


