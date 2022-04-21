
#' @note this function returns the 'valid' perturbations to be used in the analysis
#' @param data_ all data, combined into one LONG dataframe 
#' @param dataset_type to filter between P100 and GCP, since these are combined in data_
get_analysis_perturbations <- function(data_ = NA, dataset_type = NA){
  res <- data_ %>%
    filter(which_dat == dataset_type) %>% 
    dplyr::select(cell_id, pert_iname, pert_class, pr_gene_symbol, value) %>%
    group_by(pert_iname) %>%
    mutate(cell_id = factor(cell_id),
           cell_type = factor(ifelse(cell_id %in% vascular_char_vec, "vascular", "non_vascular")),
           pr_gene_symbol = factor(pr_gene_symbol)) %>%
    nest(data = c(cell_id, cell_type, pr_gene_symbol, value)) %>%
    mutate(map_df(data, .f = function(d){
      
      d_ <- d %>%
        complete(cell_id = cell_id, 
                 cell_type = cell_type, 
                 pr_gene_symbol = pr_gene_symbol) ; d_
      
      r_long <- d_ %>%
        group_by(cell_type, pr_gene_symbol) %>%
        summarize(is_not_na = sum(!is.na(value)), .groups = "keep"); r_long
      
      r <- r_long %>%
        pivot_wider(pr_gene_symbol, 
                    names_from = cell_type, 
                    values_from = is_not_na, 
                    names_prefix = "n_NOT_na_") ; r
      r2 <- r %>%
        ungroup() %>%
        summarize(n_NOT_na_non_vascular = sum(n_NOT_na_non_vascular, na.rm= TRUE),
                  n_NOT_na_vascular = sum(n_NOT_na_vascular, na.rm=TRUE)) %>%
        arrange(n_NOT_na_vascular, n_NOT_na_non_vascular); r2
      return(r2)
    })) %>%
    mutate(dim_og = map_int(data, nrow)) %>%
    filter(n_NOT_na_vascular > 0 & n_NOT_na_non_vascular > 0) %>%
    mutate(dataset_type) %>%
    ungroup()
  return(res)
}


