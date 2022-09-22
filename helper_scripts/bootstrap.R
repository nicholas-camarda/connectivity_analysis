library(boot)
library(tidyverse)
library(readxl)
library(reprex)
library(progressr)
# library(styler)

## progress bar ##
handlers(global = TRUE) # no need to wrap every call with_progress
handlers("progress")
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
  width    = 60,
  complete = "+"
))
## progress

set.seed(1)

# P100
library(reprex)
weights <- c(1/6, 1/6, 1/6, 1/8, 1/5, 1/6, 1/8, 1/6)
my_function <- function(){
  res <- rbinom(n = 8, size = 1, prob = weights)
  return(sum(res)/length(res) == 5/8)
}
reps <- replicate(n = 10000, expr = my_function(), simplify = TRUE)
sum(reps)/length(reps)

reprex()

# GCP
library(reprex)
weights2 <- c(1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8)
my_function <- function(){
  res <- rbinom(n = 8, size = 1, prob = weights2)
  return(sum(res)/length(res) == 5/8)
}
reps <- replicate(n = 10000, expr = my_function(), simplify = TRUE)
sum(reps)/length(reps)

reprex()
# set.seed(1)
# bootstrap_p100_tabs <- boot(p100_tabs, my_function, R=10000)
# summary(bootstrap_p100_tabs)


library(tidyverse)
library(readxl)
library(progressr)
fn <- file.path('/Users/ncamarda/Library/CloudStorage/OneDrive-Tufts/phd/projects/proteomics-paper/supplemental/Clustering Tabs.xlsx')
p100_tabs <- read_excel(fn, sheet = 1); p100_tabs
gcp_tabs <- read_excel(fn, sheet = 2); gcp_tabs

# try label swapping
#' @note this function swaps the labels randomly of the 1st column of a dataframe
label_swap <- function(data, i){
  new_data_temp <- data %>%
    # sample without replacement, name the column after the bootstrap replicate number
    mutate(!!sym(as.character(i)) := sample(data[,1, drop = TRUE], size = 8, replace = FALSE)) %>% #, weights_
    dplyr::select(!!sym(as.character(i)), 2) %>%
    rename(Perturbation = as.character(i), 
           !!sym(as.character(i)) := 2)
  
  return(new_data_temp)
}

#' @note this function executes the swapping procedure over n_boot
exc_lbl_swap_fn <- function(data_init, n_boot = 1000){ # , my_weights_ = NA
  swap_progress <- progressr::progressor(steps = n_boot)
  for (i_ in 1:n_boot){
    res <- label_swap(data_init, i = i_) # weights_ = my_weights_
    # save the whether each perturbation segregated for each i_
    data_init <- left_join(data_init, res, by= "Perturbation")
    # message(i_)
    swap_progress()
  }
  names_to_pivot <- colnames(data_init)[-1]
  # summarize the number of times each pert segregated, and calculate the probability
  final_res <- data_init %>% 
    pivot_longer(cols = names_to_pivot, names_to = "rep", values_to = "value") %>%
    group_by(Perturbation) %>%
    summarize(num_sgrgt = sum(value)) %>%
    mutate(n_boot = n_boot,
           ratio = num_sgrgt/n_boot)
  gc()
  return(final_res)
}

#' @results
#' @note p100 results, euclidean
euc_p100 <- exc_lbl_swap_fn(data_init = p100_tabs,
                            n_boot = 1000); euc_p100

#' @note gcp results, euclidean
euc_gcp <- exc_lbl_swap_fn(data_init = gcp_tabs,
                            n_boot = 1000); euc_gcp

reprex(comment = "#>")
