source(file.path("scripts", "init.R"))

library(magick)
library(pdftools)

p100_pdfs <- dir(file.path(output_directory, "p100", "plots"), 
                 pattern = ".pdf", recursive = T, full.names = TRUE)
gcp_pdfs <- dir(file.path(output_directory, "gcp", "plots"), 
                pattern = ".pdf", recursive = T, full.names = TRUE)

create_pdf_mapping <- function(pdfs) {
  # DEBUG: pdfs <- p100_pdfs
  res_tbl <- tibble(pdfs = pdfs) %>%
    # don't include the faceted plots
    mutate(keep = str_detect(pdfs, pattern = "faceted",negate = TRUE)) %>%
    filter(keep) %>%
    dplyr::select(-keep) %>%
    mutate(which_dat = str_extract(pdfs, pattern = "p100|gcp|P100|GCP"),
           type_data = str_extract(pdfs, pattern = "heatmaps|prettydendros|volcano_plots"),
           name_temp = str_split(basename(pdfs), "\\.", simplify = TRUE)[,1],
           cell_type = sort(str_split(name_temp, "-", simplify = TRUE)[,2]),
           cell_type = ifelse(cell_type == "", "HAoSMC,HUVEC", cell_type),
           name_ = str_split(name_temp, pattern = "-", simplify = TRUE)[,1],
           exclude_ = str_split(name_temp, pattern = "_excl_", simplify = TRUE)[,2],
           exclude_ = ifelse(exclude_ == "", "None", exclude_)) %>%
    arrange(name_) %>%
    mutate(name_copy = name_,
           type_data_copy = type_data) %>%
    group_by(name_, type_data) %>%
    nest(data = c(name_copy, type_data_copy, name_temp, exclude_, cell_type, 
                 which_dat, pdfs)); res_tbl
  
  a <- res_tbl$data[[2]]$pdfs %>% magick::image_read_pdf()
  b1 <- res_tbl$data[[1]]$pdfs[1] %>% magick::image_read_pdf()
  b2 <- res_tbl$data[[1]]$pdfs[2] %>% magick::image_read_pdf()
  c1 <- res_tbl$data[[3]]$pdfs[1] %>% magick::image_read_pdf()
  c2 <- res_tbl$data[[3]]$pdfs[2] %>% magick::image_read_pdf()
  
  top <- a %>% image_append() %>% plot()
  # savePlot("~/Downloads/top.pdf")
  dev.off()
  
}