# main libraries
library(data.table)
library(tidyverse)
library(readxl)

library(ggrepel)
library(ggsci)
library(ggforce)
library(ggpubr)
library(scales)

# colors
library(RColorBrewer)
library(randomcoloR)
library(pals)
library(colorspace)

# bioconductor
library(cmapR)
library(circlize)
library(ComplexHeatmap)

library(extrafont)
# do this once
# font_import()
extrafont::fonts()

# progress bar
library(progressr)

# clustering and trees
library(pvclust) # cluster stability, Lev's paper
library(Matching)
library(dendextend)

# formatting
library(GetoptLong)
library(gridExtra)


# multicore configuration
library(furrr)
no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)

# random seed set for global session
set.seed(25)

# name of output directory
my_output_directory <- "test-LVL4" # "All-LINCS-data-LVL4"
# name of data directory, should contain LVL3|4 info!
specific_data_directory <- "All-LINCS-data-LVL4"
# args.csv 
args_fn_name <- "test_args.csv" # "all_args.csv" 


## progress bar ##
options(ggrepel.max.overlaps = Inf, error = recover)
handlers(global = TRUE) # no need to wrap every call with_progress
handlers("progress")
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
  width    = 60,
  complete = "+"
))
## progress bar ##

## set up for multiple OS
winos <- ifelse(grepl("windows", Sys.info()["sysname"], ignore.case = T), 1, 
                ifelse(grepl("linux", Sys.info()["sysname"], ignore.case = T), 2, 0))
if (winos == 1) {
  # windows
  working_directory <- file.path(
    "C:", "Users", "ncama",
    "OneDrive - Tufts", "phd", "ws", "proteomics"
  )
} else if (winos == 0) {
  # mac
  working_directory <- file.path(
    "", "Users", "ncamarda",
    "OneDrive - Tufts", "phd", "ws", "proteomics"
  )
} else {
  # linux
  working_directory <- file.path(
    "", "home", "ncamarda93",
    "OneDrive - Tufts", "phd", "ws", "proteomics"
  )
}

# The level of data we are processing
lvl4_bool_data <- stringr::str_detect(string = specific_data_directory, 
                             pattern = "LVL4") 

if (lvl4_bool_data) {
  FC_UPPER_BOUND <- 1.1
  FC_LOWER_BOUND <- 0.9
} else {
  FC_UPPER_BOUND <- 3
  FC_LOWER_BOUND <- 0.33
}
message(qq("FC upper bound: @{FC_UPPER_BOUND}"))
message(qq("FC lower bound: @{FC_LOWER_BOUND}"))
message(qq("\nRunning analysis on @{specific_data_directory} data\n"))
if (lvl4_bool_data) {
  message("LVL4 data has been row normalized to the PLATE MEDIAN")
} else {
  message("LVL3 data has ")
}

output_directory <- file.path(working_directory, str_c("output", my_output_directory, sep = "_"))
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
source(file.path(working_directory, "helper_scripts", "connectivity_and_clustering_helper_functions.R"), local = T)
source(file.path(working_directory, "helper_scripts", "data_wrangling_helper_functions.R"), local = T)

#' @note constants
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

dir_tbl <- tribble(~dataset_type, ~output_dir,
                   "P100", p100_base_output_dir,
                   "GCP", gcp_base_output_dir)

vascular_char_vec <- c("HUVEC", "HAoSMC", "Pericyte")

##########################################


p100_fn <- file.path("combined-datasets", "P100-All-Cell-Lines.gct")
gcp_fn <- file.path("combined-datasets", "GCP All Cell Lines.gct")

cat("\nReading and merging data...")
#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
#' the values in this data are Z-scores!!!

#' # test_args.csv
analysis_fn <- file.path(data_directory, args_fn_name)
analysis_dat_temp <- read_csv(analysis_fn,
                              comment = "#",
                              show_col_types = FALSE
) %>%
  mutate_all(str_trim) %>%
  mutate(
    filter_vars = map(filter_vars, collect_args),
    exclude = map(exclude, collect_args)
  ) %>%
  left_join(dir_tbl, by = "dataset_type")

stopifnot(nrow(analysis_dat_temp) >= 1)

my_data <- tibble(fns = dir(file.path(datasets_directory, specific_data_directory),
                            full.names = T, recursive = T
)) %>%
  distinct() %>%
  mutate(gct = map(
    .x = fns,
    .f = parse_gctx
  )) %>%
  mutate(dataset_type = str_extract(string = fns, pattern = "GCP|P100")) %>%
  arrange(desc(dataset_type)) %>%
  mutate(data = map(.x = gct, .f = melt_gct)) %>%
  inner_join(analysis_dat_temp, by = "dataset_type")

drugs_moa_df <- create_my_drugs_df(ref_dir = references_directory)

# write summary file of data paths
file_summary_dat <- my_data %>%
  distinct(dataset_type, fns) %>%
  dplyr::select(dataset_type, fns) %>%
  mutate(fns = basename(fns)) %>%
  mutate(plate = str_extract(string = fns, pattern = "[P|p]late[0-9]*[a-z]*")) %>%
  dplyr::select(dataset_type, plate, fns)
write_tsv(file_summary_dat, file = file.path(data_directory, "fn_list.tsv"))

# organize what was read in
my_data_lst <- my_data %>%
  split(.$dataset_type)
dataset_type_col <- names(my_data_lst)


message("Reading and summarizing data...")
#' If you want to edit the raw data, edit in *read_and_summarize_data*
my_data_obj_final <- tibble(
  dataset_type = dataset_type_col,
  # combine the data rows into a single dataset, and nest it
  data = map2(my_data_lst, dataset_type,
              .f = read_and_summarize_data
  )
)

# grab drugs from BOTH P100 and GCP
my_temp_obj <- bind_rows(my_data_obj_final$data) %>%
  ungroup()

HUVEC_HAoSMC_perts <- my_temp_obj %>%
  filter(cell_id %in% vascular_char_vec) %>%
  ungroup() %>%
  dplyr::select(pert_iname) %>%
  .$pert_iname %>%
  unique()

other_perts <- my_temp_obj %>%
  ungroup() %>%
  filter(!(cell_id %in% vascular_char_vec)) %>%
  distinct(pert_iname) %>%
  .$pert_iname

my_perts <- intersect(HUVEC_HAoSMC_perts, other_perts); my_perts


all_drugs_mapping <- my_temp_obj %>%
  distinct(pert_iname) %>%
  filter(pert_iname %in% my_perts)


# bind final object
analysis_dat <- inner_join(
  analysis_dat_temp,
  my_data_obj_final, # my_data_obj_final
  by = "dataset_type"
)

grouping_var_for_summary <- unlist(analysis_dat_temp$grouping_var)
filter_vars_for_summary <- unlist(analysis_dat_temp$filter_vars)
get_drug_and_cell_summary_data(analysis_dat = analysis_dat, 
                               output_dir = output_dir,
                               filter_vars_for_summary, 
                               grouping_var_for_summary)
