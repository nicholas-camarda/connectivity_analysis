library(tidyverse)
library(readxl)
library(cmapR)

# none of these contain the VEGFR inhibitors
fn3 <- "/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/Cardio-oncology/Datasets/Inherited Data/3rd gen data/HUVEC 2+3 P100 - TKI and CV.txt"
fn2 <- "/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/Cardio-oncology/Datasets/Inherited Data/2nd gen data/P100/HUVEC2 P100.gct"
fn1 <- "/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/Cardio-oncology/Datasets/Inherited Data/1st gen data/P100/P100-All-Cell-Lines.gct"

complete_table1 <- read_delim(fn1, delim = "\t", skip = 2, show_col_types = F) 
complete_table2 <- read_delim(fn2, delim = "\t", skip = 2, show_col_types = F)
complete_table3 <- read_delim(fn3, delim = "\t", skip = 2, show_col_types = F)

meta <- complete_table3 %>%
  dplyr::select(id, `PHS-1FE9D-040D04`:`PHC-139BE-069F09`) %>%
  dplyr::slice(1:32)
meta_long <- meta %>%
  pivot_longer()

data <- complete_table1 %>%
  dplyr::select(c(id:pr_uniprot_id)) %>%
  dplyr::slice(33:123) 



