# source the init.R file firstso
source(file.path("scripts", "init.R"))


# lst_of_genes <- c("DDX54", "AHNAK", "NUP214", "RBM17", "BRD4", 
#                   "HN1", "RBBP6", "BAT2", "ZC3HC1", "PPP1R10", 
#                   "C17orf85", "LIMA1")
# extract_sites(analysis_dat, lst_of_genes)

my_obj <- analysis_dat$data[[1]] %>%
  filter(pert_class == "VEGFR inhibitor")

df_temp <- my_obj %>% 
  ungroup() %>%
  dplyr::select(master_id, pr_gene_symbol,cell_id, mark, pert_iname, value) %>%
  pivot_wider(id_cols = c(master_id, cell_id, pert_iname), 
              names_from = pr_gene_symbol, 
              values_from = value, 
              values_fn = median)

ids <- df_temp$master_id
mat <- as.matrix(df_temp[-c(1:3)]) %>% t()
colnames(mat) <- ids

corr <- cor(mat, use = "pairwise.complete.obs")
Heatmap(corr)