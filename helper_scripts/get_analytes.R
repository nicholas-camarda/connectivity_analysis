## THIS REQUIRES get_drugs.R output

#' @note this function returns the 'valid' analytes to be used in the analysis
#' @param data_split_ all data, split by dataset type, attached with exclusion, grouping variable info
#' @param my_pert_df from get_drugs.R, df containing the perturbations for the analysis
#' @param perc_trash threshold; if the percent of NA's is greater than this number for an analyte, remove it
#' @return analysis_dat with an extra column *filtered_data* which should be renamed to *data* in init.R when filtered_analytes flag = TRUE
get_analysis_dat_with_filtered_analytes <- function(data_split_ = NA, my_perts_df_ = NA, perc_trash = 0.75){
  # data_split_ = analysis_dat_pre; my_perts_df_ = my_perts_df; perc_trash = 0.75
  
  filt_lst_res <- apply(data_split_, 1, FUN = function(args){
    # args <- data_split_[1,]; my_obj <- args$data[[1]]
    grouping_var <- args$grouping_var
    filter_vars <- force_natural(args$filter_vars)
    if (any(is.na(force_natural(args$exclude)))) {
      exclude <- ""
    } else {
      exclude <- force_natural(args$exclude)
    }
    my_obj <- args$data
    
    new_pert_df <- my_perts_df_ %>% 
      filter(!(pert_iname %in% exclude)) %>%
      filter(!!sym(grouping_var) == filter_vars); new_pert_df
    
    res <- my_obj %>% filter(pert_iname %in% new_pert_df$pert_iname)
    return(res)
  }) 
  
  analysis_dat_filtered_all <- data_split_ %>% 
    bind_cols(tibble(perturbation_filtered_data = filt_lst_res))
  
  # i_ = 4
  analyte_and_pert_filtered_dat <- analysis_dat_filtered_all %>%
    mutate(pert_and_analyte_filtered_data = map(perturbation_filtered_data, .f = function(d){
      # d <- analyte_and_pert_filtered_dat$perturbation_filtered_data[[3]]
      
      wide_d <- d %>%
        dplyr::select(replicate_id, cell_id, pert_iname, 
                      pert_class, pr_gene_symbol, value) %>%
        pivot_wider(id_cols = replicate_id:pert_class,
                    names_from = pr_gene_symbol, values_from = value, 
                    values_fn = function(x) median(x, na.rm = TRUE))
      
      gene_cols <- colnames(wide_d)[grepl(x = colnames(wide_d), pattern = "^p[A-Z]|^H3")]
      
      analytes_in_huvec <- wide_d %>% filter(cell_id == "HUVEC") 
      remove_analytes_huvec <- tibble(frac_na = map_dbl(analytes_in_huvec, function(x) sum(is.na(x))/length(x)),
                             names = colnames(analytes_in_huvec)) %>%
        filter(str_detect(names, "^p[A-Z]|^H[3-4]")) %>%
        arrange(desc(frac_na)) %>%
        filter(frac_na > perc_trash); remove_analytes_huvec
      
      # smooth
      analytes_in_haosmc <- wide_d %>%
        mutate(temp = replicate_id, .before = 1) %>%
        separate(temp, into = c("cell_id", "pert_iname", "pert_class"), sep = "--") %>%
        filter(cell_id %in% "HAoSMC")
      
      remove_analytes_haosmc <- tibble(frac_na = map_dbl(analytes_in_haosmc, function(x) sum(is.na(x))/length(x)),
                                       names = colnames(analytes_in_haosmc)) %>%
        filter(str_detect(names, "^p[A-Z]|^H[3-4]")) %>%
        arrange(desc(frac_na)) %>%
        filter(frac_na > perc_trash); remove_analytes_haosmc
      
      
      # non-vascular
      analytes_in_nonvasc <- wide_d %>%
        mutate(temp = replicate_id, .before = 1) %>%
        separate(temp, into = c("cell_id", "pert_iname", "pert_class"), sep = "--") %>%
        filter(!(cell_id %in% vascular_char_vec))
      
      remove_analytes_nonvasc <- tibble(frac_na = map_dbl(analytes_in_nonvasc, function(x) sum(is.na(x))/length(x)),
                                        names = colnames(analytes_in_nonvasc)) %>%
        filter(str_detect(names, "^p[A-Z]|^H[3-4]")) %>%
        arrange(desc(frac_na)) %>%
        filter(frac_na > perc_trash); remove_analytes_nonvasc
      
      remove_analytes <- bind_rows(remove_analytes_haosmc, remove_analytes_huvec,
                                   remove_analytes_nonvasc) %>%
        distinct(names)
      
      res_d <- d %>%
        filter(!(pr_gene_symbol %in% remove_analytes$names))
      return(res_d)
    }))
  
  return(analyte_and_pert_filtered_dat)
}

