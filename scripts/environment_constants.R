# make sure to set these before analysis
# these are the defaults!
filter_analytes <- FALSE
dendro_cut_thresh <- 0.6 # cut trees at 60% the maximum height
bh_thresh_val <- 0.1 # filter out analytes that are strictly above this threshold

message("##############################")
message(qq("\nLoaded env variables:\nfilter_analytes = @{filter_analytes}\ndendro_cut_thresh = @{dendro_cut_thresh}\nbh_thresh_val = @{bh_thresh_val}\n\n"))
message("##############################")

### deprecated ###
# set_run_organization <- c("pert_iname", "drug_class", "all")
# rerun_clustering <- TRUE
# rerun_diffe <- TRUE
# cancer_vs_non_cancer <- TRUE 
# plot_morpheus_toggle <- FALSE 
# message(qq("\nLoaded env variables:\nfilter_analytes = @{filter_analytes}\nrerun_clustering = @{rerun_clustering}\nrerun_diffe = @{rerun_diffe}\nset_run_organization = @{str_c(set_run_organization, collapse=',')}\ncancer_vs_non_cancer = @{cancer_vs_non_cancer}\nplot_morpheus_toggle = @{plot_morpheus_toggle}\ndendro_cut_thresh = @{dendro_cut_thresh}\nbh_thresh_val = @{bh_thresh_val}\n\n"))
