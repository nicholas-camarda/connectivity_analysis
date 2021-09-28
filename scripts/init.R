library(tidyverse)
library(ggrepel)
library(readxl)
library(ggsci)
library(ggpubr)

library(xml2)
library(XML)

library(RColorBrewer)

# bioconductor pacakges that need to be installed manually
library(cmapR)
library(circlize)
library(ComplexHeatmap)
library(BiocParallel)

# clustering
library(cluster)
library(factoextra)
library(dbscan)
library(fpc)

library(progressr)

library(minerva) # non-linear correlation analysis with mutual information theory
library(matrixStats)
library(pvclust) # cluster stability, Lev's paper
library(Matching)

library(dendextend)

library(GetoptLong)
library(gridExtra)


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("ComplexHeatmap", "BiocParallel", "cmapR", "circlize"))

## progress bar ##
handlers(global = TRUE) # no need to wrap every call with_progress
handlers("progress")
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
  width    = 60,
  complete = "+"
))
## progress bar ##


winos <- ifelse(grepl("windows", Sys.info()["sysname"], ignore.case = T), 1, 
                ifelse(grepl("linux", Sys.info()["sysname"], ignore.case = T), 2, 0))
if (winos == 1) {
  working_directory <- file.path(
    "C:", "Users", "ncama",
    "OneDrive - Tufts", "phd", "ws"
  )
} else if (winos == 0) {
  working_directory <- file.path(
    "", "Users", "ncamarda",
    "OneDrive - Tufts", "phd", "ws", "proteomics"
  )
} else {
  working_directory <- file.path(
    "", "home", "ncamarda93",
    "OneDrive - Tufts", "phd", "ws", "proteomics"
  )
}
output_directory <- file.path(working_directory, "output")
data_directory <- file.path(working_directory, "data")

setwd(working_directory)

references_directory <- file.path(data_directory, "references")
datasets_directory <- file.path(data_directory, "datasets")

print_important_directories <- function() {
  message(qq("\n\nWorking directory [working_directory] = @{working_directory}"))
  message(qq("Output directory [output_directory] = @{output_directory}"))
  message(qq("Cardio-oncology directory [data_directory] = @{data_directory}"))
  message(qq("References directory [references_directory] = @{references_directory}"))
  message(qq("Datasets directory [datasets_directory] = @{datasets_directory}"))
}
print_important_directories()

## load helper functions ##
source(file.path(working_directory, "scripts", "dendrograms.R"), local = T)
source(file.path(working_directory, "scripts", "heatmaps.R"), local = T)
source(file.path(working_directory, "helper_scripts", "differential_analyte_helper_functions.R"), local = T)
source(file.path(working_directory, "helper_scripts", "connectivity_and_clustering_helper_functions.R"), local = T)
source(file.path(working_directory, "helper_scripts", "data_wrangling_helper_functions.R"), local = T)
# source(file.path(working_directory, "helper_scripts", "morpheus_helper_functions.R"), local = T)


#' @note constants
drug_bank_database_fn <- file.path(references_directory, "full database.xml")
rerun_clustering <- TRUE; rerun_diffe <- TRUE;
set_run_organization <- c("pert_iname", "drug_class", "all")
cancer_vs_non_cancer <- TRUE; plot_morpheus_toggle <- FALSE
dendro_cut_thresh <- 0.6; bh_thresh_val <- 0.1
message(qq("\nLoaded env variables: \nrerun_clustering = @{rerun_clustering}\nrerun_diffe = @{rerun_diffe}\nset_run_organization = @{str_c(set_run_organization, collapse=',')}\ncancer_vs_non_cancer = @{cancer_vs_non_cancer}\nplot_morpheus_toggle = @{plot_morpheus_toggle}\ndendro_cut_thresh = @{dendro_cut_thresh}\nbh_thresh_val = @{bh_thresh_val}\n\n"))

#' @note load output directories
p100_base_output_dir <- file.path(output_directory, "p100")
dir.create(p100_base_output_dir, recursive = T, showWarnings = F)
gcp_base_output_dir <- file.path(output_directory, "gcp")
dir.create(gcp_base_output_dir, recursive = T, showWarnings = F)
avg_base_output_dir <- file.path(output_directory, "avg")
dir.create(avg_base_output_dir, recursive = T, showWarnings = F)

dir_tbl <- tribble(~dataset_type, ~output_dir,
                   "P100", p100_base_output_dir,
                   "GCP", gcp_base_output_dir,
                   "AVG", avg_base_output_dir)

vascular_char_vec <- c("HUVEC", "HAoSMC")

# p100_fn <- file.path("combined-datasets", "P100-All-Cell-Lines.gct")
# gcp_fn <- file.path("combined-datasets", "GCP All Cell Lines.gct")

