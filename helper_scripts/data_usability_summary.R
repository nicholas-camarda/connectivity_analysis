# summarize data usability
source(file.path("scripts","init.R"))

#' @note this function plots radar of usable data across cell lines for a given perturbation class
#' @param data_ all_data from init.R
#' @param dataset_type character string of either P100 or GCP
#' @param conditions_to_show vector of the perturbation classes in question, e.g. c("Kinase inhibitor")
get_usable_data_plot <- function(data_ = NA, dataset_type = "P100", 
                                 conditions_to_show = c("Kinase inhibitor", "Epigenetic", "control")){
  
  # get number of NON NA's per drug, by cell id
  # essentially a summary of the workable data for each cell type and drug
  usable_data_df <- data_ %>%
    filter(which_dat == dataset_type) %>%
    dplyr::select(cell_id, pert_iname, pert_class, pr_gene_symbol, value) %>%
    group_by(pert_iname) %>%
    mutate(cell_id = factor(cell_id),
           pr_gene_symbol = factor(pr_gene_symbol)) %>%
    nest(data = c(cell_id, pr_gene_symbol, value)) %>%
    mutate(map_df(data, .f = function(d){
      
      d_ <- d %>%
        complete(cell_id = cell_id, 
                 pr_gene_symbol = pr_gene_symbol) ; d_
      
      r_long <- d_ %>%
        group_by(cell_id, pr_gene_symbol) %>%
        summarize(is_not_na = sum(!is.na(value)), .groups = "keep"); r_long
      
      r <- r_long %>%
        pivot_wider(pr_gene_symbol, 
                    names_from = cell_id, 
                    values_from = is_not_na) ; r
      cell_ids <- levels(d$cell_id)
      r2 <- r %>%
        ungroup() %>%
        summarize(across(.cols = all_of(cell_ids), .fns = sum));  r2
      return(r2)
    })) %>%
    mutate(dataset_type) %>%
    ungroup()
  
  # scale the data for use in radar plot, proportion of usable data
  meta_dat <- usable_data_df %>% 
    dplyr::select(pert_iname, pert_class, data, dataset_type)
  scaled_dat <- usable_data_df %>% 
    dplyr::select(-c(pert_iname, pert_class, data, dataset_type)) %>%
    mutate(across(everything(), scales::rescale))
  
  radar_plot_df_master <- bind_cols(meta_dat, scaled_dat) %>%
    filter(pert_class %in% conditions_to_show) %>%
    split(.$pert_class)
  
  grad_lst <- tibble(pert_class = names(radar_plot_df_master)) %>% 
    mutate(gplot = future_map(radar_plot_df_master, .f = function(df_){
      
      radar_plot_df <- df_ %>% 
        dplyr::select(-pert_class, -data, -dataset_type) %>%
        filter(HAoSMC > 0 | HUVEC > 0)
      
      dataset_type_ <- unique(df_$dataset_type)
      class_type <- unique(df_$pert_class)
      # print(class_type)
      
      # whole graphs together
      grad <- ggradar(radar_plot_df, 
                      legend.title = class_type,
                      grid.label.size = 4,
                      axis.label.size = 4,
                      group.point.size = 2,
                      group.line.width = 0.8,
                      legend.text.size = 10) +
        labs(title = qq("Potentially usable data across\ncell lines in @{dataset_type_}")); grad
      
      
      # split graphs
      split_radar_plot_df <- radar_plot_df %>%
        split(.$pert_iname) 
      grad_split <- tibble(pert_iname = names(split_radar_plot_df)) %>%
        mutate(gplot = map(split_radar_plot_df, .f = function(r){
          pert_name <- unique(r$pert_iname)
          # print(pert_iname)
          ggradar(r, 
                  grid.label.size = 4,
                  axis.label.size = 4,
                  group.point.size = 2,
                  group.line.width = 0.8,
                  legend.text.size = 10) +
            labs(title = qq("Potentially usable data across\ncell lines in @{dataset_type_}"),
                 subtitle = pert_name)
        }))
      
      final_result <- list(grad, grad_split) %>% set_names("main", "by_pert")
      return(final_result)
    })) %>%
    unnest(c(gplot)) %>%
    mutate(name = names(gplot)) %>%
    pivot_wider(id_cols = pert_class, 
                names_from = name, values_from = gplot)
  
  main_output_directory <- file.path(output_directory, dataset_type, "summary")
  dir.create(main_output_directory, showWarnings = FALSE, recursive = TRUE)
  by_pert_output_directory <- file.path(output_directory, dataset_type, "summary", "by_perturbation")
  dir.create(by_pert_output_directory, showWarnings = FALSE, recursive = TRUE)
  
  message("Plotting main...")
  grad_lst_main <- grad_lst %>%
    dplyr::select(pert_class, main) %>%
    mutate(output_fn = file.path(main_output_directory,
                                 str_c(str_c(dataset_type,pert_class,sep = "-"),".pdf")))
  future_walk2(grad_lst_main$main, as.list(grad_lst_main$output_fn), .f = function(x, y) ggsave(plot = x, filename = y, width = 10))
  
  message("Plotting by perturbation...")
  grad_lst_pert <- grad_lst %>%
    dplyr::select(pert_class, by_pert) %>%
    unnest(c(by_pert)) %>%
    mutate(output_fn = file.path(by_pert_output_directory,
                                 str_c(str_c(dataset_type, pert_iname,sep = "-"),".pdf")))
  
  future_walk2(grad_lst_pert$gplot, as.list(grad_lst_pert$output_fn), .f = function(x, y) ggsave(plot = x, filename = y, width = 10))
  
  return(grad_lst)
}

all_data_temp <- read_rds(obj_final_fn)
all_data <- all_data_temp$data %>% bind_rows() %>% ungroup()

message(qq("Plotting usable data radar graph in output directory: @{my_output_directory}/{gcp|p100}/summary"))
gcp_radar <- get_usable_data_plot(data_ = all_data, dataset_type = "GCP")
p100_radar <- get_usable_data_plot(data_ = all_data, dataset_type = "P100")

