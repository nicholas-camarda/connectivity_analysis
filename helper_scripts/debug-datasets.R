# source the init.R file first
source(file.path("scripts", "init.R"))

cat("\nReading and merging data...")
#' @note cmapR help functions https://rdrr.io/github/cmap/cmapR/man/
#' the values in this data are Z-scores!!!

#' TODO change this to [analysis_args.csv] when ready
analysis_fn <- file.path(data_directory, "all_args.csv")
analysis_dat_temp <- read_csv(analysis_fn, comment = "#") %>%
  mutate_all(str_trim) %>%
  mutate(filter_vars = map(filter_vars, collect_args)) %>%
  left_join(dir_tbl, by = "dataset_type")

dir_tbl
#' @note read in both GCP and P100 since it's cheap, then analyze accordingly

srilas_data <- tibble(fns = c( file.path(datasets_directory, p100_fn),
                               file.path(datasets_directory, gcp_fn)), 
                      gct = map(
                        .x = fns,
                        .f = parse_gctx
                      )) %>% 
  mutate(dataset_type = str_extract(string = fns, pattern = "GCP|P100")) %>%
  mutate(data = map(gct, .f = function(x) as_tibble(melt_gct(x)))) 

srilas_data_final <- srilas_data$data[[1]]%>%
  mutate(det_plate = str_extract(string = det_plate, pattern = "P-[0-9]*"))

# read in all 
my_data <- tibble(fns = dir(file.path(datasets_directory, "All-LINCS-data-LVL4"), full.names = T, recursive = T)) %>%
  distinct() %>% 
  mutate(gct = map(
    .x = fns,
    .f = parse_gctx
  )) %>% 
  mutate(dataset_type = str_extract(string = fns, pattern = "GCP|P100")) %>%
  arrange(desc(dataset_type)) %>%
  mutate(data = map(.x = gct, .f = melt_gct)) %>%
  filter(dataset_type == "P100")

message("Merging... Matching to Srila's data")
my_data_temp_final <- data.table::rbindlist(my_data$data, fill = TRUE, use.names = TRUE) %>% 
  as_tibble() %>%
  mutate(cell_id = ifelse(cell_id == "Pericytes", "Pericyte", cell_id)) %>%
  mutate(det_plate = str_extract(string = det_plate, pattern = "P-[0-9]*"))
  # mutate(det_plate = str_replace(string = det_plate, pattern = "G", replacement = "P")) 

setdiff(srilas_data_final$det_plate, my_data_temp_final$det_plate)
shared_plate_ids <- intersect(srilas_data_final$det_plate, my_data_temp_final$det_plate)
all_plate_ids <- union(srilas_data_final$det_plate, my_data_temp_final$det_plate)

my_data_final <- my_data_temp_final %>% filter(det_plate %in% shared_plate_ids)
dir.create("/Users/ncamarda/OneDrive - Tufts/phd/ws/proteomics/debug")
write_tsv(tibble(shared_plate_ids = shared_plate_ids), 
          file = "/Users/ncamarda/OneDrive - Tufts/phd/ws/proteomics/debug/shared_plate_ids.tsv")
# less data probably means more filtered?

# create drug df from excel file of drugs
drugs_moa_df <- create_my_drugs_df(ref_dir = references_directory)

all_together <- bind_rows(srilas_data_final %>% mutate(data_id = "Srila") %>% 
                            dplyr::mutate(across(.cols = c(det_normalization_group_vector, det_well_enrichment_score, 
                                                           pert_batch_internal_compound_enumerator, pert_batch_internal_replicate, 
                                                           pert_time, pubchem_cid), as.factor)), 
                          my_data_final %>% mutate(data_id = "Panorama") %>% 
                            dplyr::mutate(across(.cols = c(det_normalization_group_vector, det_well_enrichment_score, 
                                                           pert_batch_internal_compound_enumerator, pert_batch_internal_replicate, 
                                                           pert_time, pubchem_cid), as.factor))) %>%
  arrange(cell_id) %>%
  # filter(!(det_plate %in% c("P-0057", "P-0069"))) %>%
  filter(det_plate %in% all_plate_ids) %>%
  mutate(pert_iname = tolower(pert_iname)) %>%
  inner_join(drugs_moa_df, by = "pert_iname")

vascular_cells <- c("HUVEC", "HAoSMC")
HUVEC_HAoSMC_perts <- all_together %>% filter(cell_id %in% vascular_cells) %>% distinct(pert_iname) %>% .$pert_iname
other_perts <- all_together %>% filter(!(cell_id %in% vascular_cells)) %>% distinct(pert_iname) %>% .$pert_iname
my_perts <- intersect(HUVEC_HAoSMC_perts, other_perts)

my_theme <- theme(axis.text.x = element_text(angle = 60, vjust = 0.5))


g1 <- ggboxplot(all_together, x = "cell_id", y = "value", color = "data_id", 
                ggtheme = theme_bw(), add = "jitter", 
                palette = "jco") +
  ggtitle("All data") +
  ylab("Normalized phosphosite expression") + 
  xlab("Cell line") + 
  my_theme

g1a <- ggboxplot(all_together %>% filter(pert_iname %in% my_perts),
                 x = "cell_id", y = "value", color = "data_id", 
                 ggtheme = theme_bw(), add = "jitter", 
                 palette = "jco", facet.by = "pert_iname") +
  ggtitle("Data distribution by perturbations") +
  labs(caption = "(only those in both vascular cells and cancer cells)") +
  ylab("Normalized phosphosite expression") + 
  xlab("Cell line") + 
  my_theme; g1a

# kolmogorov smirnov to test whether the two datasets are drawn from different distributions by cell_id and pert_iname
different_dist <- all_together %>% 
  filter(pert_iname %in% my_perts) %>%
  group_by(cell_id, data_id) %>%
  dplyr::select(cell_id,pert_iname, data_id, value) %>%
  pivot_wider(id_cols = pert_iname:cell_id, data_id, 
              values_fn = function(x) median(x, na.rm = T)) %>%
  group_by(pert_iname) %>%
  mutate(p_val_ks = map2_dbl(.x = Srila, .y = Panorama, .f = function(x,y) {
    if (length(x) > 2 & length(y) > 2) {
      r <- ks.test(x,y)$p.value 
    } else {
      r <- NA
    }
    return(r)
  })) %>%
  # if significant, these two were draw from separate distribtions
  mutate(different_distributions = ifelse(p_val_ks < 0.05, T, F))

ecdf_dat <- all_together %>% 
  filter(pert_iname %in% my_perts) %>%
  # separate by cell and perturbation, because these are different files
  group_by(cell_id, pert_iname) %>%
  dplyr::select(value, cell_id, pert_iname, data_id) %>%
  arrange(cell_id, pert_iname) %>%
  nest(d = c(value, data_id)) %>%
  left_join(different_dist %>% dplyr::select(-Srila, -Panorama)) %>%
  left_join(drugs_moa_df); ecdf_dat 

plot_ecdfs <- function(data_to_plot, parent_dir = "all"){
  apply(data_to_plot, 1, FUN = function(dat){
    test_stat <- str_c(str_c("ks.test p-val: ", round(dat$p_val_ks, 5)), 
                       str_c("\ndifferent dist?: ", dat$different_distributions))
    p <- unique(dat$pert_iname)
    c <- unique(dat$cell_id)
    n <- str_c(p, c, sep = " || ")
    message(n)
    g <- ggecdf(dat$d, x = "value",
                color = "data_id", linetype = "data_id",
                palette = c("#00AFBB", "#E7B800")) +
      ggtitle(n) + 
      labs(caption = test_stat) 
    same_dist <- ifelse(unique(dat$different_distributions) == TRUE, "diff", "same")
    ecdf_dir <- file.path(analysis_dat_temp$output_dir, "ecdf", parent_dir, same_dist, c)
    dir.create(ecdf_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(g, file= file.path(ecdf_dir, qq("@{n}.pdf")), width = 5, height = 5)
  })
}
plot_ecdfs(data_to_plot = ecdf_dat %>% filter(pert_class == "Kinase inhibitor"), parent_dir = "kinase-inhibitor")
plot_ecdfs(data_to_plot = ecdf_dat %>% filter(pert_class == "Epigenetic"), parent_dir = "epigenetics")
plot_ecdfs(data_to_plot = ecdf_dat, parent_dir = "all")




g2 <- ggboxplot(all_together %>% filter(cell_id %in% c("HUVEC", "HAoSMC") ), 
                x = "cell_id", y = "value", color = "data_id",
                ggtheme = theme_bw(), add = "jitter", palette = "jco")+
  ggtitle("Only HUVECs and HAoSMCS") +
  ylab("Normalized phosphosite expression") + 
  xlab("Cell line")+ 
  my_theme

g3b <- ggboxplot(all_together %>% filter(pert_class == "Kinase inhibitor"), 
                x = "det_plate", y = "value", color = "data_id",
                ggtheme = theme_bw(), add = "jitter")+
  ggtitle("All cell lines, only Kinase inhibitors") +
  ylab("Normalized phosphosite expression") + 
  xlab("Cell line")+ 
  facet_grid(cell_id ~ pert_iname, scales = "free_x") +
  my_theme; g3b

g3 <- ggboxplot(all_together %>% filter(pert_class == "Kinase inhibitor"), 
                x = "cell_id", y = "value", color = "data_id",
                ggtheme = theme_bw(), add = "jitter")+
  ggtitle("All cell lines, only Kinase inhibitors") +
  ylab("Normalized phosphosite expression") + 
  xlab("Cell line") +
  my_theme; g3

g4 <- ggboxplot(all_together %>% filter(pert_iname == "dmso"), 
                x = "cell_id", y = "value", color = "data_id",
                ggtheme = theme_bw(), add = "jitter", palette = "jco")+
  ggtitle("All cell lines, only dmso") +
  ylab("Normalized phosphosite expression") + 
  xlab("Cell line")+ 
  my_theme

debug_datasets_dir <- file.path(analysis_dat_temp$output_dir, "debug-datasets")
dir.create(debug_datasets_dir, showWarnings = FALSE, recursive = TRUE)
glst <- ggarrange(g1, g2, g3, g4, nrow = 2, ncol = 2)
ggsave(glst, file = file.path(debug_datasets_dir, "glst.pdf"), width = 12, height = 12)
ggsave(g1a, file = file.path(debug_datasets_dir, "g1a.pdf"), width = 12, height = 12)

