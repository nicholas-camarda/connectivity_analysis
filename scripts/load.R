# main libraries
library(data.table)
library(tidyverse)
library(readxl)

library(ggrepel)
library(ggsci)
library(ggforce)
library(ggpubr)
library(scales)
library(ggradar)

# colors
library(RColorBrewer)
library(randomcoloR)
library(pals)
library(colorspace)

# bioconductor
library(cmapR)
library(circlize)
library(ComplexHeatmap)

# progress bar
library(progressr)

# clustering and trees
library(pvclust) # cluster stability, Lev's paper
library(Matching)
library(dendextend)

# formatting
library(GetoptLong)
library(gridExtra)
library(tictoc)

# multicore configuration
library(furrr)
no_cores <- availableCores() - 1
plan(multisession, workers = no_cores)

# random seed set for global session
set.seed(25)



################### SESSION INFO ##################
message("#############################################################")
message("########## Running cell-cell connectivity analysis ##########")
message("#############################################################\n")
message(Sys.time())
message()

#' @note constants
# source(file.path("scripts", "environment_constants.R"))

# session options
options(ggrepel.max.overlaps = Inf)

## progress bar ##
handlers(global = TRUE) # no need to wrap every call with_progress
handlers("progress")
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
  width    = 60,
  complete = "+"
))
## progress bar ##

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
message(qq("\nRunning analysis on data contained in @{specific_data_directory} folder"))
if (lvl4_bool_data) {
  message("**LVL4** data has been row normalized to the PLATE MEDIAN\n")
} else {
  message("**LVL3** data\n")
}

output_directory_temp <- file.path(str_c("output", my_output_directory, sep = "_"))
data_directory <- file.path("data")

references_directory <- file.path(data_directory, "references")
datasets_directory <- file.path(data_directory, "datasets")

print_important_directories <- function() {
  # message(qq("\n\nWorking directory [working_directory] = @{getwd()}"))
  message(qq("Output directory [output_directory] = @{output_directory_temp}"))
  message(qq("Cardio-oncology directory [data_directory] = @{data_directory}"))
  message(qq("References directory [references_directory] = @{references_directory}"))
  message(qq("Datasets directory [datasets_directory] = @{datasets_directory}"))
}
print_important_directories()

## load helper functions ##
source(file.path("scripts", "dendrograms.R"), local = T)
source(file.path("scripts", "heatmaps.R"), local = T)
source(file.path("helper_scripts", "connectivity_and_clustering_helper_functions.R"), local = T)
source(file.path("helper_scripts", "data_wrangling_helper_functions.R"), local = T)
source(file.path("helper_scripts", "get_analytes.R"), local = T)
source(file.path("helper_scripts", "get_drugs.R"), local = T)


#' @note load output directories
p100_base_output_dir <- file.path(output_directory_temp, "p100")
dir.create(p100_base_output_dir, recursive = T, showWarnings = F)
gcp_base_output_dir <- file.path(output_directory_temp, "gcp")
dir.create(gcp_base_output_dir, recursive = T, showWarnings = F)

dir_tbl <- tribble(~dataset_type, ~output_dir,
                   "P100", p100_base_output_dir,
                   "GCP", gcp_base_output_dir)

##########################################

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
# make sure we read in the analysis files correctly
stopifnot(nrow(analysis_dat_temp) >= 1)

# grab the inital listing of drugs and mechanism of action
drugs_moa_df <- create_my_drugs_df(ref_dir = references_directory)

# find the combined dataset, if not made - make it; if made, read it
obj_final_fn <- file.path(data_directory, "datasets", "combined-datasets", 
                          "analysis_ready-formatted_datasets.rds")
# note what the filenames are that are going into the combined dataset
fn_data_lst_dir <- file.path(data_directory, "datasets", specific_data_directory)
if (!file.exists(obj_final_fn)) {
  message("Reading gct* files and melting into long form...")
  my_data <- tibble(fns = dir(file.path(datasets_directory, specific_data_directory),
                              pattern = "gct*",
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
  
  # write summary file of data paths
  file_summary_dat <- my_data %>%
    distinct(dataset_type, fns) %>%
    dplyr::select(dataset_type, fns) %>%
    mutate(fns = basename(fns)) %>%
    mutate(plate = str_extract(string = fns, pattern = "[P|p]late[0-9]*[a-z]*")) %>%
    dplyr::select(dataset_type, plate, fns)
  
  write_tsv(file_summary_dat, file = file.path(fn_data_lst_dir, "fn_list.tsv"))
  message(qq("Wrote data file names to fn_lst.tsv in:\n@{fn_data_lst_dir}\n"))
  
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
  
  message(qq("Writing combined dataset file to: \n@{obj_final_fn}\n"))
  write_rds(my_data_obj_final, obj_final_fn)
} else {
  message(qq("Data files used in analysis stored in: @{fn_data_lst_dir}/fn_list.tsv"))
  message(qq("Reading combined dataset file stored in: @{obj_final_fn}\n"))
  my_data_obj_final <- read_rds(obj_final_fn)
}

# grab drugs from BOTH P100 and GCP
all_data <- bind_rows(my_data_obj_final$data) %>%
  ungroup()

# run get_drugs.R first to generate the appropriate files in the appropriate place!
pert_fn <- file.path(data_directory, "perturbation_data", "my_perts.rds")
if (file.exists(pert_fn)) {
  message("Reading my_perts.rds ...")
  my_perts_df <- read_rds(pert_fn)
  message("Done.\n")
} else {
  message("Couldn't find my_perts_df... preparing for conditions Epigenetic")
  
  ### apply functions and get output ---
  p100_perts <- get_analysis_perturbations(data_ = all_data, dataset_type = "P100")
  p100_perts_edit <- p100_perts %>%
    dplyr::select(pert_iname, pert_class, dataset_type)
  gcp_perts <- get_analysis_perturbations(data_ = all_data, dataset_type = "GCP")
  gcp_perts_edit <- gcp_perts %>%
    dplyr::select(pert_iname, pert_class, dataset_type)
  
  my_perts_df_temp <- bind_rows(p100_perts_edit, gcp_perts_edit) %>%
    pivot_wider(id_cols = pert_iname:pert_class, 
                names_from = dataset_type, 
                values_from = dataset_type) %>%
    # edit this area to pick which classes of perturbations we retain
    mutate(keep = ifelse(pert_class == "Epigenetic" & !is.na(GCP), T,
                         ifelse(pert_class == "Kinase inhibitor" & !is.na(P100), T, 
                                ifelse(pert_class == "control", T, F))))
  # my_perts! ---
  my_perts_df <- my_perts_df_temp %>%
    filter(keep) %>%
    arrange(pert_class); my_perts_df
  
  # write all this info to files
  pert_data_dir <- file.path(data_directory, "perturbation_data")
  dir.create(pert_data_dir, showWarnings = F, recursive = T)
  
  # this is the one we'll use for the analysis
  write_rds(my_perts_df, file.path(pert_data_dir, "my_perts.rds"))
  
  # these are helper files in case we want to look back
  write_tsv(p100_perts %>% dplyr::select(-data), file.path(pert_data_dir, "p100_pert.tsv"))
  write_rds(p100_perts, file.path(pert_data_dir, "p100_pert.rds"))
  write_tsv(gcp_perts %>% dplyr::select(-data), file.path(pert_data_dir, "gcp_pert.tsv"))
  write_rds(gcp_perts, file.path(pert_data_dir, "gcp_pert.rds"))
  message("Done.\n")
}

my_perts <- my_perts_df$pert_iname

# make perturbation mapping df now, then filter out drugs in analysis loop in connectivity-analysis.R
all_drugs_mapping <- all_data %>%
  distinct(pert_iname) %>%
  filter(pert_iname %in% my_perts)


# bind final object
analysis_dat_pre <- inner_join(
  analysis_dat_temp,
  my_data_obj_final, # my_data_obj_final
  by = "dataset_type"
)

if (filter_analytes){
  message("Running analysis on P100 and GCP ***with analyte filtering***!")

  analysis_dat <- get_analysis_dat_with_filtered_analytes(data_split_ = analysis_dat_pre, 
                                                          my_perts_df_ = my_perts_df, 
                                                          perc_trash = percent_analyte_na_to_throwout) %>%
    rename(unfiltered_data = data, 
           data = pert_and_analyte_filtered_data)
  output_directory <- str_c(output_directory_temp, "_filtered_analytes", sep = "")
  message(qq("New output directory [output_directory] = @{output_directory}"))
} else {
  message("Running analysis on P100 and GCP, as is - ***no analyte filtering***!")
  analysis_dat <- analysis_dat_pre
  output_directory <- output_directory_temp
  message(qq("Unchanged output directory [output_directory] = @{output_directory}"))
}

message("\nFiles accumulated, ready for analysis\n")
