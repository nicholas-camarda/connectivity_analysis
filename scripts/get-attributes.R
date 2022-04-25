library(tidyverse)
library(readxl)
library(stringdist)

source(file.path("helper_scripts", "data_wrangling_helper_functions.R"))
references_directory <- file.path("data", "references")

#' @note find drug names available based on our master excel spreadsheet *NOT ALL INCLUSIVE*
#' @param fuzzy_drug_name a string that can partially match a drug name in the dataset
#' @param print_all_unique toggle to TRUE to print all the unique drug names in this excel file
#' @return either a list of all the unique drug names OR a dataframe sorted by the similarity of your query 
#' to drug names in the excel file
find_drug_name <- function(fuzzy_drug_name = "ly", print_all_unique = FALSE){
  drugs_df <- create_my_drugs_df(ref_dir = references_directory)
  if (print_all_unique){
    return(drugs_df$Drug %>% unique())
  }
  
  drugs_sim <- create_my_drugs_df(ref_dir = references_directory) %>%
    dplyr::select(Drug, Class, MOA, `MOA simplified`) %>%
    mutate(query = fuzzy_drug_name, .before = 1,
           drug = tolower(Drug)) %>%
    mutate(similarity_ = stringsim(query, drug), .before = 1) %>%
    arrange(desc(similarity_))
  return(drugs_sim)
}


#' @note find drug classes available based on our master excel spreadsheet *NOT ALL INCLUSIVE*
#' @param fuzzy_drug_class a string that can partially match a drug class in the dataset
#' @param print_all_unique toggle to TRUE to print all the unique drug classes in this excel file
#' @return either a list of all the unique drug classes OR a dataframe sorted by the similarity of your query 
#' to drug classes in the excel file
find_drug_class <- function(fuzzy_drug_class = "Epigen", print_all_unique = FALSE){
  drugs_df <- create_my_drugs_df(ref_dir = references_directory)
  if (print_all_unique){
    return(drugs_df$Class %>% unique())
  }
  drugs_class_sim <- drugs_df %>%
    dplyr::select(Drug, Class, MOA, `MOA simplified`) %>%
    mutate(query = fuzzy_drug_class, .before = 1,
           class_ = tolower(Class)) %>%
    mutate(similarity_ = stringsim(query, class_), .before = 1) %>%
    arrange(desc(similarity_))
  return(drugs_class_sim)
}

# find_drug_class(fuzzy_drug_class = "Epi")
find_drug_class(print_all_unique = TRUE)

find_drug_name(fuzzy_drug_name = "bos")
find_drug_name(print_all_unique = TRUE)
