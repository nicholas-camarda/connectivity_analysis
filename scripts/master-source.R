library(tidyverse)
library(data.table)
library(readxl)
library(cmapR)
library(amap) # clustering of correlations according to Lev's paper
library(pvclust) # cluster stability, Lev's paper
library(Matching)
library(dendextend)

library(GetoptLong)
library(gridExtra)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)

library(morpheus) # heat maps
library(viridis)
library(plotly)
library(gplots)
library(htmlwidgets)
library(leaflet)
library(ggplotify)
library(gridGraphics)
library(ggrepel)
library(ggpubr)

# pathway analysis
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(limma)


library(BiocParallel)
library(conflicted)
library(showtext)

set.seed(42) 
# ncores <- 6
# registerDoMC(cores=ncores)
conflict_prefer("paste", "base")
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
conflict_prefer("transpose", "purrr")
conflict_prefer("layout", "plotly")
conflict_prefer("set", "dendextend")
conflict_prefer("rotate", "dendextend")
conflict_prefer("rename", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("reduce", "purrr")
conflict_prefer("rename", "dplyr")
conflict_prefer("clusterExport", "parallel")
conflict_prefer("parLapply", "parallel")
conflict_prefer("clusterEvalQ", "parallel")
# font_add_google("Noto Sans")

## IMPORTANT DIRECTORIES ##

OUTPUT_DIR_NAME <- "output"

if (winos == 1){
  WORKING_DIRECTORY <- "C:\\Users\\ncama\\OneDrive - Tufts\\phd\\ws"
  OUTPUT_DIRECTORY <- file.path("C:\\Users\\ncama\\phd\\proteomics", OUTPUT_DIR_NAME)
  CARDONC_DIRECTORY <- file.path("C:\\Users\\ncama\\OneDrive - Tufts\\phd\\Cardio-oncology")
} else {
  WORKING_DIRECTORY <- "/Users/Nicholas/OneDrive - Tufts/phd/ws/proteomics"
  OUTPUT_DIRECTORY <- file.path("/Users/Nicholas/Desktop/coding/proteomics", OUTPUT_DIR_NAME)
  CARDONC_DIRECTORY <- file.path("/Users/Nicholas/OneDrive - Tufts/phd/Cardio-oncology")
}

setwd(WORKING_DIRECTORY)



message(qq("\n\nWORKING DIRECTORY [WORKING_DIRECTORY] = @{WORKING_DIRECTORY}"))
message(qq("OUTPUT DIRECTORY [OUTPUT_DIRECTORY] = @{OUTPUT_DIRECTORY}"))

REFERENCES_DIRECTORY <- file.path(CARDONC_DIRECTORY, "References")
DATASETS_DIRECTORY <- file.path(CARDONC_DIRECTORY, "Datasets")

message(qq("Cardio-oncology directory [CARDONC_DIRECTORY] = @{CARDONC_DIRECTORY}"))
message(qq("References directory [REFERENCES_DIRECTORY] = @{REFERENCES_DIRECTORY}"))
message(qq("Dasets directory [DATASETS_DIRECTORY] = @{DATASETS_DIRECTORY}"))

## load helper functions ##
source(file.path(WORKING_DIRECTORY, "scripts", "dendrograms.R"), local = T)
source(file.path(WORKING_DIRECTORY, "scripts", "heatmaps.R"), local = T)
source(file.path(WORKING_DIRECTORY, "helper_scripts", "differential_analyte_helper_functions.R"), local = T)
source(file.path(WORKING_DIRECTORY, "helper_scripts", "connectivity_and_clustering_helper_functions.R"), local = T)
source(file.path(WORKING_DIRECTORY, "helper_scripts", "data_wrangling_helper_functions.R"), local = T)
source(file.path(WORKING_DIRECTORY, "helper_scripts", "morpheus_helper_functions.R"), local = T)


#' @note constants
RERUN_CLUSTERING <- TRUE; RERUN_DIFFE <- TRUE;
SET_RUN_ORGANIZATION <- c("pert_iname", "drug_class","all")
CANCER_VS_NON_CANCER <- TRUE; PLOT_MORPHEUS_TOGGLE <- FALSE
DENDRO_CUT_THRESH <- 0.6; BH_THRESH_VAL <- 0.1
message(qq("\nLoaded env variables: \nRERUN_CLUSTERING = @{RERUN_CLUSTERING}\nRERUN_DIFFE = @{RERUN_DIFFE}\nSET_RUN_ORGANIZATION = @{str_c(SET_RUN_ORGANIZATION, collapse=',')}\nCANCER_VS_NON_CANCER = @{CANCER_VS_NON_CANCER}\nPLOT_MORPHEUS_TOGGLE = @{PLOT_MORPHEUS_TOGGLE}\nDENDRO_CUT_THRESH = @{DENDRO_CUT_THRESH}\nBH_THRESH_VAL = @{BH_THRESH_VAL}\n\n"))

#' @note load output directories
p100_base_output_dir <- file.path(OUTPUT_DIRECTORY, "p100")
dir.create(p100_base_output_dir, recursive = T, showWarnings = F)
gcp_base_output_dir <- file.path(OUTPUT_DIRECTORY, "gcp")
dir.create(gcp_base_output_dir, recursive = T, showWarnings = F)
avg_base_output_dir <- file.path(OUTPUT_DIRECTORY,"avg")
dir.create(avg_base_output_dir, recursive = T, showWarnings = F)

#' #' @note load data
# lst_dat <- load_wrapper(datasets_directory = DATASETS_DIRECTORY,
#'                         references_directory = REFERENCES_DIRECTORY, 
#'                         load_from_file = TRUE,
#'                         filter_common_by_data_set = FALSE) 
#' 
#' p100_lst_obj <- lst_dat$p100; gcp_lst_obj <- lst_dat$gcp; MASTER_KEY <- lst_dat$master_key
#' message("Loaded data: \np100_lst_obj, gcp_lst_obj, MASTER_KEY")
#' 
#' p100_final_res_fn <- file.path(p100_base_output_dir, "p100-pert_iname_drug_class_all-final_res-clusters-diffex.rds")
#' if (file.exists(p100_final_res_fn)) {
#'   message("\nFound a cached p100 final res object. \nLoading p100_res ...")
#'   p100_res <- read_rds(p100_final_res_fn)
#' } else {
#'   p100_res <- NULL
#' }
#' 
#' gcp_final_res_fn <- file.path(gcp_base_output_dir, "gcp-pert_iname_drug_class_all-final_res-clusters-diffex.rds")
#' if (file.exists(gcp_final_res_fn)) {
#'   message("\nFound a cached gcp final res object. \nLoading gcp_res ...")
#'   gcp_res <- read_rds(gcp_final_res_fn)
#' } else {
#'   gcp_res <- NULL
#' }

# avg_final_res_fn <- file.path(avg_base_output_dir, "avg-final_res-clusters-diffex.rds")
# if (file.exists(gcp_final_res_fn)) {
#   message("\nFound a cached avg final res object. \nLoading avg_res ...")
#   avg_res <- read_rds(avg_final_res_fn)
# }



