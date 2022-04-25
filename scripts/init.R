#' @NOTE *this is the scripts that runs EVERYTHING*
#' Change (or keep) these values before sourcing connectivity-analysis.R
# Flow of scripts:: init.R -> load.R -> connectivity-analysis.R

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

# filter out analytes that are strictly above this threshold
bh_thresh_val <- 0.1


message("##############################")
message(qq("\nLoaded env variables:\nfilter_analytes = @{filter_analytes}\ndendro_cut_thresh = @{dendro_cut_thresh}\nbh_thresh_val = @{bh_thresh_val}\n"))
message("##############################")

source(file.path("scripts", "load.R"))
source(file.path("scripts", "connectivity-analysis.R"))
