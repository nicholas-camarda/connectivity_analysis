library(tidyverse)


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


### Body ---

all_data_split <- read_rds(file.path("data", "datasets", "combined-datasets", "analysis_ready-formatted_datasets.rds"))
all_data <- bind_rows(all_data_split$data) %>% ungroup()

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
  mutate(keep = ifelse(pert_class == "Epigenetic" & !is.na(GCP), T,
                       ifelse(pert_class == "Kinase inhibitor" & !is.na(P100), T, 
                              ifelse(pert_class == "control", T, F))))
# my_perts! ---
my_perts_df <- my_perts_df_temp %>%
  filter(keep) %>%
  arrange(pert_class); my_perts_df

# write to file
pert_data_dir <- file.path(data_directory, "perturbation_data")
dir.create(pert_data_dir, showWarnings = F, recursive = T)

# this is the one we'll use for the analysis
write_rds(my_perts_df, file.path(pert_data_dir, "my_perts.rds"))

# these are helper files in case we want to look back
write_tsv(p100_perts %>% dplyr::select(-data), file.path(pert_data_dir, "p100_pert.tsv"))
write_rds(p100_perts, file.path(pert_data_dir, "p100_pert.rds"))
write_tsv(gcp_perts %>% dplyr::select(-data), file.path(pert_data_dir, "gcp_pert.tsv"))
write_rds(gcp_perts, file.path(pert_data_dir, "gcp_pert.rds"))
