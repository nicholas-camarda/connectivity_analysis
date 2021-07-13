#' @note run differential analyte analysis
#' @param conn_clust_obj from *run_connectivity_analysis* 
#' @param base_output_dir for output files
#' @param dataset for naming, e.g. differentiation between gcp and p100
#' @param rerun logical to rerun or load existing analysis file
#' @param group_var grouping variable passed from batch_exec
#' @param dendro_cut_thresh threshold to cut the dendrograms produced in clustering step previously
#' @param bh_thresh threshold under which analyte is called significant using benjamini-hochberg
run_differential_analyte_analysis <- function(conn_clust_obj, bh_thresh = 0.1, dendro_cut_thresh = dendro_cut_thresh,
                                              dataset = "p100", group_var = "pert_iname", rerun = F, 
                                              base_output_dir_specific = "~/Downloads", master_key = NULL){
  # DEBUG: 
  # read in final rds and debug
  # conn_clust_obj <- p100_only_staurosporine_pert %>% select(-which_dat, cut_trees, co_clust, diff_ex); 
  # group_var <- "pert_iname"
  # dataset <- "p100"
  # bh_thresh <- 0.1
  # dendro_cut_thresh <- 0.6
  # base_output_dir_specific <- "~/Downloads/test"
  # master_key <- MASTER_KEY
  
  # ensure string matching with lowercase
  group_var <- tolower(group_var)
  stopifnot(group_var %in% c("pert_iname","drug_class","all")); group_var
  if (group_var == "all") group_var <- "org"
  
  # DEBUG:
  
  # args <- list(clust_assignments = conn_clust_trees$cut_trees, dat = conn_clust_trees$df,
  #      dname = conn_clust_trees[group_var][[1]], co_clust = conn_clust_trees$co_clust)
  # # gid <- c("HUVEC","HAoSMC")
  # # cancer_vs_non_cancer <- T
  # res <- pmap(.l = args, 
  #             .f = calc_differential_analytes, 
  #             which.dat = dataset, gid = c("HUVEC","HAoSMC"),
  #             cancer_vs_non_cancer = cancer_vs_non_cancer,
  #             base_output_dir = base_output_dir_specific,
  #             dendro_cut_thresh = dendro_cut_thresh,
  #             bh_thresh = bh_thresh, master_key = master_key)
  
  # differential expression
  
  out_fn <- file.path(base_output_dir_specific, qq("@{dataset}_diff_ex-h_@{dendro_cut_thresh}-q_@{bh_thresh}.rds"))
  if (!file.exists(out_fn)) {
    message("Differential analyte .rds doesn't exist..")
    rerun <- TRUE
  }
  
  if (rerun){
    message("\nStarting differential expression analysis..")
    diff_ex_obj <- conn_clust_obj %>%
      mutate(diff_ex = pmap(list(clust_assignments = cut_trees, dat = df, 
                                 dname = !!sym(group_var), co_clust = co_clust), 
                            .f = calc_differential_analytes, 
                            dendro_cut_thresh = dendro_cut_thresh,
                            which.dat = dataset, gid = c("HUVEC","HAoSMC"),
                            # cancer_vs_non_cancer = cancer_vs_non_cancer,
                            base_output_dir = base_output_dir_specific,
                            bh_thresh = bh_thresh, master_key = master_key))
    
    out_fn <- file.path(base_output_dir_specific, qq("@{dataset}_diff_ex-h_@{dendro_cut_thresh}-q_@{bh_thresh}.rds"))
    message(qq("\nWriting combined output to: \n@{out_fn}"))
    write_rds(diff_ex_obj, out_fn, compress = "gz")
    
    message("Plotting volcanos....")
    pwalk(.l = list(diff_ex_obj$diff_ex, diff_ex_obj[group_var] %>% pull, diff_ex_obj$cut_trees, diff_ex_obj$co_clust),
          .f = plot_diff_analyte_results_CVs_html, 
          which.dat = dataset, 
          gid = c("HUVEC", "HAoSMC"), dendro_cut_thresh = dendro_cut_thresh,
          # cancer_vs_non_cancer = cancer_vs_non_cancer, 
          base_output_dir = base_output_dir_specific, bh_thresh = bh_thresh)
  } else {
    message("Reading diffe .rds ... ")
    diff_ex_obj <- read_rds(out_fn)
  }

  message("Done!")
  
  return(diff_ex_obj)
}


#' @note Using cluster assignments, calculate which analytes are differentially expressed between clusters
#' So, compare 1 vs others, 2 vs others, 3 vs others
#' @param clust_assignments named vector of cell_id cluster assignments
#' @param dat dataframe of analyte values (replicates x analyte)
#' @param dname grouping variable
#' @param co_clust true or false boolean indicating whether vars of interest co-clustered
#' @param base_output_dir output directory for summary plot
calc_differential_analytes <- function(clust_assignments, dat, dname, co_clust, which.dat = "p100",
                                       base_output_dir, 
                                       gid = c("HUVEC","HAoSMC"),
                                       dendro_cut_thresh = 0.6, bh_thresh = 0.1, 
                                       master_key = NULL){
  
  # DEBUG:
  
  # a = 1
  # clust_assignments = conn_clust_trees$cut_trees[[a]];
  # dat = conn_clust_trees$df[[a]];
  # dname = conn_clust_trees[group_var][[1]][a];
  # co_clust = conn_clust_trees$co_clust[[a]];
  # which.dat <- "p100"
  # gid <- c("HUVEC","HAoSMC")
  # cancer_vs_non_cancer <- T
  
  # attach cluster assignments to data
  message(qq("@{dname}, p.val cutoff = @{bh_thresh}"))
  if (is.null(clust_assignments)) {
    message("No clusters to compute differential analytes. Skipping.")
    return(NULL)
  }
  
   # join clusters and annotations
  ca_df_temp <- create_ca_df(ca = clust_assignments) %>%
    rename(base_clust_comp = cluster); ca_df_temp
  
  # which cluster(s) are vascular cells in
  which_cluster_vasc_df <- ca_df_temp %>% 
    filter(cell_id %in% gid)
  
  # get cluster labels
  clust_label_df <- ca_df_temp %>% 
    group_by(base_clust_comp) %>%
    summarize(base_clust_comp_name = str_c(unique(cell_id), collapse = ","), .groups = "keep") %>%
    ungroup() %>%
    arrange(base_clust_comp); clust_label_df

  ca_df <- left_join(ca_df_temp, clust_label_df, by="base_clust_comp")
  
  dat_ <- dat %>% 
    left_join(ca_df, by=c("cell_id")) 
  
  # start by getting column classes and finding where the data is
  # order the columns with the numeric, 'matrix' columns at the end
  ordered_column_names_df <- tibble(col_name = colnames(dat_), 
                           col_class = sapply(dat_, class)) %>% # we are using class here, so expect numeric, not double
    arrange(col_class) ; ordered_column_names_df
  
  # note - cluster assingment is NOT a numeric! it is an integer!!
  start_col_idx <- which(ordered_column_names_df$col_class == "numeric")[1]; start_col_idx
  # grabs the first column that's a numeric
  # separate out our data from annots
  
  features_idx <- start_col_idx:nrow(ordered_column_names_df); tail(features_idx) # i.e. cols to analyze
  features_names <- ordered_column_names_df$col_name[features_idx]; tail(features_names)

  clusts_int_vec <- dat_ %>% 
    select(base_clust_comp) %>% # grab the correct cluster column
    pluck(1) %>% # convert this 1-d dataframe into a vector
    unique() %>% # unique it and sort it
    sort(); clusts_int_vec
  
  # order the df with annots in front and data in back
  ordered_dat_ <- select(dat_, all_of(ordered_column_names_df$col_name)); ordered_dat_
  # get a 'matrix' for diffe 
  matrix_for_diffe <- select(ordered_dat_, base_clust_comp, all_of(features_names)); matrix_for_diffe
  # 1:num clusters
  
  # get a unique df to look up cluster names and integer assignments 
  # that incoudes both base_clust_comp and cvnc_cluster names and integer cluster assignments
  unique_clust_labels_df <- ca_df %>% 
    select(-cell_id) %>% 
    distinct()
  
  # stop()
  
  message("Computing differential analytes..")
  message("Note: '*' indicates that at least one group didn't have enough data for comparison")
  diffe_by_clust_df <- map_df(clusts_int_vec, function(ith_cluster){
    # ith_cluster <- clusts_int_vec[3]
    
    base_clust_comp_name <- unique_clust_labels_df %>% 
      select(base_clust_comp, base_clust_comp_name) %>% # take the correct cluster naming convention column
      filter(base_clust_comp == ith_cluster) %>% # take the correct name according to the cluster integer
      distinct(base_clust_comp_name) %>% # non-vascular clusters need to get lumped together, we just want the unique cluster assign
      pluck(1) # make it a character vector 
    
    
    
    cat("\n",ith_cluster, base_clust_comp_name, "vs all others \n")
    
    fail <- numeric(length(features_names))
    cluster_res <- map2_df(features_names, ith_cluster, function(k, i){
      # k <- features_names[1]
      # DEBUG: 
      # message(k)

      c1 <- matrix_for_diffe %>% filter(base_clust_comp == i) %>% pluck(k)
      c2 <- matrix_for_diffe %>% filter(base_clust_comp != i) %>% pluck(k)
      
      n_non_na_vec1 <- length(c1[!is.na(c1)]); n_non_na_vec1
      n_non_na_vec2 <- length(c2[!is.na(c2)]); n_non_na_vec2
      
      if (n_non_na_vec1 < 1 | n_non_na_vec2 < 1) {
        # message("Not enough data to compute diffe for ", k)
        res <- tibble(analyte = k, 
                      logFC = NA, base_clust_comp = as.integer(i), 
                      base_clust_comp_name = base_clust_comp_name, 
                      ks_statistic = NA, 
                      ks_boot_statistic = NA, 
                      signal_to_noise_statistic = NA,
                      directional_stat = NA, 
                      p_val = NA, p_val_boot = NA,
                      p_val_bh = NA,
                      p_val_boot_bh = NA, 
                      signif = NA) %>%
          mutate(k_clust_dat = list(c1), 
                 all_others_dat = list(c2))
        
        cat("*")
        return(res)
      }
      
      ks <- ks.test(c1, c2, alternative = "two.sided", tol = 1e-8); ks
      ks_boot <- ks.boot(c1, c2, nboots = 1000, alternative = "two.sided"); ks_boot
   
      
      # compute ks stat
      #quantifies the distance between the empirical distribution of given two samples
      stat <- ks$statistic
      stat_boot <- ks_boot$ks$statistic
      
      c1_median <- median(c1, na.rm = T); c1_median
      c2_median <- median(c2, na.rm = T); c2_median
       
      #' @note assumes data are already log transformed!!! So a fold change is simply a difference, not a division!!
      logfc <- median(c1, na.rm = T) - median(c2, na.rm = T)
      d_stat <- ifelse(logfc > 0, "++", "--");
      
      
      # plot the ecdf to look at cluster separation
      min_x <- min(c(min(c1, na.rm = T), min(c2, na.rm = T)))
      max_x <- max(c(max(c1, na.rm = T), max(c2, na.rm = T)))

      sub_dir_name <- str_replace_all(string = base_clust_comp_name, 
                                      pattern = ',', replacement = '-');sub_dir_name
      ecdf_dir <- file.path(base_output_dir, "ecdf", dname, sub_dir_name); # ecdf/dmso/cluster1 .. cluster2... 
      dir.create(ecdf_dir, recursive = T, showWarnings = F)
      
      c1_names_df <- ordered_dat_ %>% 
        filter(base_clust_comp == i) %>% 
        mutate(id = str_c(base_clust_comp_name, pert_iname, sep = "-"),
               plot_clust_id = 1) %>%
        select(plot_clust_id, id, all_of(k)) %>%
        rename(value := !!sym(k)); c1_names_df
      
      # c1_ecdf <- ecdf(c1_names_df$value)
      
      c2_names_df <- ordered_dat_ %>% 
        filter(base_clust_comp != i) %>% 
        mutate(id = str_c(base_clust_comp_name, pert_iname, sep = "-"),
               plot_clust_id = 2) %>%
        select(plot_clust_id, id, all_of(k)) %>%
        rename(value := !!sym(k)); c2_names_df
      # c2_ecdf <- ecdf(c2_names_df$value)
      
      ecdf_names_df <- bind_rows(c1_names_df, c2_names_df) %>%
        mutate(plot_clust_id = factor(plot_clust_id, levels=c(1,2))); ecdf_names_df
      
      ecdf_plot <- ggplot(ecdf_names_df) + 
        stat_ecdf(aes(x = value, color = plot_clust_id)) +
        scale_color_manual(values = c("#00AFBB", "#E7B800"),  labels = c(qq("Base cluster @{i}"), "All others")) + 
        labs(title = qq("@{dname} eCDFs of @{base_clust_comp_name} vs all\nAnalyte @{k}"), 
             y = "f(value)", 
             color = "cluster") +
        theme_minimal()
      ggsave(filename = file.path(ecdf_dir, qq("@{k}_ecdfs.eps")), plot = ecdf_plot, width = 5, height = 5, device = "eps")
  

      
      # compute signal to noise
      mean_grp_diff <- mean(c2, na.rm = T) - mean(c1, na.rm = T)
      sum_grp_sd <- sd(c2, na.rm = T) + sd(c1, na.rm = T)
      signal_to_nose <- mean_grp_diff / sum_grp_sd
      
      #'@note I'm pretty sure we want to adjust the p-values per cluster comparison... 
      res <- tibble(analyte = k,  logFC = logfc, base_clust_comp = i, 
                    base_clust_comp_name = base_clust_comp_name, 
                    ks_statistic = stat, 
                    ks_boot_statistic = stat_boot, 
                    signal_to_noise_statistic = signal_to_nose,
                    directional_stat = d_stat, 
                    p_val = ks$p.value, p_val_boot = ks_boot$ks.boot.pvalue,
                    p_val_bh = p.adjust(ks$p.value, method = "BH"),
                    p_val_boot_bh = p.adjust(ks_boot$ks.boot.pvalue, method = "BH")) %>%
        mutate(signif = p_val_bh < bh_thresh) %>%
        mutate(k_clust_dat = list(c1), 
               all_others_dat = list(c2))
      
      # na_omit_res <- res %>% na.omit()
      cat(".")
      return(res)
    }) 
    
    cat(".Done\n")
    
    return(cluster_res)
  })
  

 
  # i th clust
  # k th analyte
  # i <- 1
  # k <- "ATXN2L" # 1
  
  phospho_master_df <- master_key %>% 
    distinct(pr_gene_symbol, phosphosite) 
  diffe_by_clust_w_phos_df <- diffe_by_clust_df %>% 
    left_join(phospho_master_df, by=c("analyte"="pr_gene_symbol")) %>%
    select(analyte, phosphosite, everything()); diffe_by_clust_w_phos_df
  
  # diffe_by_clust = diffe_by_clust_w_phos; which.dat = "p100"; gid = c("HUVEC", "HAoSMC"); cluster_assignments_df = ca_df; dname = dname; dendo_thresh = DENDRO_CUT_THRESH; bh_thresh = BH_THRESH_VAL
  
  morpheus_folder <- file.path(base_output_dir, "data_for_morpheus_plotting", dname)
  dir.create(morpheus_folder, recursive = T, showWarnings = F)
  
  morpheus_fn <- file.path(morpheus_folder, qq("@{dname}-for_morpheus.csv"))
  annot_fn <- file.path(morpheus_folder, qq("@{dname}-annots.csv"))
  raw_fn <- file.path(morpheus_folder, qq("@{dname}-raw.csv"))
  
  morpheus_lst <- write_morpheus_diffex_file(data = ordered_dat_, dname = dname, gid = gid,
                                             morpheus_fn = morpheus_fn, 
                                             start_col_idx = start_col_idx)
  
  write_csv(morpheus_lst$raw_dat, path = raw_fn)
  write_csv(morpheus_lst$annots, path = annot_fn)
  
  
  return(diffe_by_clust_w_phos_df)
}



#' @note this functions wraps the verbose method of converting a matrix into a 'gct'-like matrix for 
#' plotting in morpheus 
#' @param data differential expression, essentially cells x analyte
#' @param dname drug name (or group)
#' @param morpheus_fn output fn
#' @param start_col_idx essentially the index of data where analytes start, or the number of attributes (cells, wells, groups)
#' @return transposed matrix with analytes in rows and cells (and all attributes) in columns
write_morpheus_diffex_file <- function(data, dname = dname, gid = c("HUVEC","HAoSMC"),
                                       morpheus_fn = morpheus_fn, 
                                       start_col_idx = start_col_idx){
  
  # message("\nWriting rds diff_ex_obj, and writing csv for morpheus...")
  # transpose this matrix!!
  
  
  informative_labels_lst <- create_informative_labels() # this is from the dendrograms.R script
  il_df <- informative_labels_lst$informative_labels_df %>% 
    select(lbs, grp) %>% 
    mutate(grp = as.character(grp)) %>%
    rename(cell_id = lbs, group = grp)
  
  for_morpheus <- data %>% 
    group_by(cell_id) %>%
    left_join(il_df, by = "cell_id") %>%
    mutate(u_cell_id = make.unique(cell_id)) %>% # added one new column here
    ungroup() %>%
    select(well, cell_id, u_cell_id, group, pert_iname, drug_class, 
           base_clust_comp, base_clust_comp_name,
           # cvnc_cluster, cvnc_cluster_name, 
           everything())
  
  
  new_start_col_idx <- (start_col_idx+2); new_start_col_idx # because we added "u_cell_id, group"
  g_for_morpheus <- for_morpheus %>%
    gather(analyte, pex, all_of(new_start_col_idx:ncol(.))); g_for_morpheus
  
  # create 'gct' kind of document
  well_cell_info <- for_morpheus %>% distinct(well, u_cell_id); well_cell_info
  clust_info <- for_morpheus %>% distinct(u_cell_id, base_clust_comp, base_clust_comp_name); clust_info
  group_info <- for_morpheus %>% distinct(u_cell_id, group); group_info
  pert_info <- for_morpheus %>% select(drug_class,pert_iname); pert_info
  
  cname_total <- nrow(for_morpheus) + 1
  wells_row <- tibble(well = c("well",  well_cell_info$well), cnames = 1:cname_total) %>%
    pivot_wider(names_from = cnames, values_from = well);wells_row
  cells_row <- tibble(cell_id =  c("cell_id",  well_cell_info$u_cell_id), cnames = 1:cname_total) %>%
    pivot_wider(names_from = cnames, values_from = cell_id)
  clust_row <- tibble(base_clust_comp = c("cluster", clust_info$base_clust_comp), cnames = 1:cname_total) %>%
    pivot_wider(names_from = cnames, values_from = base_clust_comp)
  clust_name_row <-tibble(cluster = c("cluster_name", clust_info$base_clust_comp_name), cnames = 1:cname_total) %>%
    pivot_wider(names_from = cnames, values_from = cluster)
  group_row <- tibble(group = c("group", group_info$group), cnames = 1:cname_total ) %>%
    pivot_wider(names_from = cnames, values_from = group)

  gpert_row <- tibble(drug_class = c("drug_class", pert_info$drug_class), cnames = 1:cname_total) %>%
    pivot_wider(names_from = cnames, values_from = drug_class)
  ipert_row <- tibble(pert_iname = c("pert_iname", pert_info$pert_iname), cnames = 1:cname_total) %>%
    pivot_wider(names_from = cnames, values_from = pert_iname)
  
  annots <- bind_rows(wells_row, cells_row, clust_row, clust_name_row, group_row, gpert_row, ipert_row)
  
  to_bind_with_annots <- g_for_morpheus %>%
    select(-cell_id, -well, -base_clust_comp, -base_clust_comp_name, -group, 
           # -cvnc_cluster, -cvnc_cluster_name, 
           -pert_iname, -drug_class) %>%
    spread(u_cell_id,pex)
  
  colnames(annots) <- colnames(to_bind_with_annots)
  
  # bind annots on top
  t_for_morpheus <- rbindlist(list(annots, to_bind_with_annots))
  write_csv(t_for_morpheus, morpheus_fn)

  message("..Done.\n")
  
  return(list("for_morpheus" = t_for_morpheus, "annots" = annots, "raw_dat" = to_bind_with_annots))
}
         


#' @note plot differential analyte results using a faceted volcano plot 
#' AND displays most differential analytes -- for html only
#' @param diffe_by_clust differential analyte results with cluster assignments
#' @param which.dat which dataset is primary
#' @param cluster_assignments_df tbl of cluster assignments
#' @param dname grouping variable, e.g. name of a perturbation ID (passed in here for labeling the graph)
#' @param bh_thresh bh threshold to plot analyte names (passed in here for labeling the graph)
plot_diff_analyte_results_CVs_html <- function(diffe_by_clust, dname, cut_trees, co_clust,
                                               which.dat, gid = c("HUVEC", "HAoSMC"), 
                                               base_output_dir = base_output_dir, 
                                               dendro_cut_thresh = dendro_cut_thresh, bh_thresh = bh_thresh){
  
  if (is.null(diffe_by_clust)) return("No clustering was performed for this condition.")
  
  # diffe_by_clust <- diff_ex_obj$diff_ex[[1]]; cut_trees <- diff_ex_obj$cut_trees[[1]]; dname <- diff_ex_obj$pert_iname[[1]]; gid = c("HUVEC", "HAoSMC")
  cluster_assignments_df <- create_ca_df(ca = cut_trees) 
  
  ca_df_gid <- cluster_assignments_df %>% 
    filter(cell_id %in% gid)
  
  which_cluster_vasc_df <- ca_df_gid %>%
    distinct(cluster)
  
  clust_label <- cluster_assignments_df %>% 
    group_by(cluster) %>%
    summarize(name = str_c(unique(cell_id), collapse = ",")) %>%
    ungroup() %>%
    rename(base_clust_comp = cluster)
  
  # do this twice - all and cancer vs non cancer
  
  # cancer vs non-cancer
  clust_label1 <- clust_label %>%
    semi_join(ca_df_gid, by=c("base_clust_comp" = "cluster")) %>% # if this is left_join, then we get all combinations of clusters
    distinct(name, base_clust_comp)
  
  # this doesn't really work for plotly....
  # all 
  clust_label2 <- clust_label %>%
    left_join(ca_df_gid, by=c("base_clust_comp" = "cluster")) %>% # if this is left_join, then we get all combinations of clusters
    distinct(name, base_clust_comp)
  
    list_clust_labels <- list(clust_label1, clust_label2) # these are dataframes
    clust_types <- list("grps", "all") # these are strings, obviously
    
    # get plots for each type 
    walk2(.x = list_clust_labels, .y = clust_types, function(clust_label, clust_type){
      
      to_plot <- diffe_by_clust %>% 
        left_join(which_cluster_vasc_df, by=c("base_clust_comp" = "cluster")) %>%
        left_join(clust_label, by="base_clust_comp") %>%
        mutate(`-log(p_val_bh)` = -log(p_val_bh,base = 10)) %>%
        na.omit() %>%
        mutate(plot_analyte_name = str_c(analyte, phosphosite, sep=" "))
      
      # message("Checking if plotting df has no NAs...")
      stopifnot(nrow(to_plot) > 0 )
      # message("Good!")
      
      min_max <- to_plot %>% 
        arrange(p_val_bh,logFC) %>%
        filter(signif)
      
      plot_subtitle <- str_c(clust_label$name, collapse=' and '); plot_subtitle
      
      # plotly interactive html
      plotly_p <- to_plot %>%
        plot_ly(x = ~ logFC, y = ~ `-log(p_val_bh)`, color = ~ `-log(p_val_bh)`) %>%
        add_markers(hoverinfo = "text", #  text = ~ plot_analyte_name
                    text = ~paste('</br> analyte: ', plot_analyte_name,
                                  '</br> -log10(q-val): ', round(`-log(p_val_bh)`,3),
                                  '</br> q-val: ', round(p_val_bh,3),
                                  '</br> log2 fold change: ', round(logFC,3))) %>%
        
        layout(shapes = list(plotly_hline(y = -log(bh_thresh, base=10), 
                                          name = qq("BH threshold: q < @{bh_thresh}"), 
                                          color = "black", 
                                          dash = "dot")),
               title = list(text = qq("@{toupper(which.dat)} | @{dname} \n@{plot_subtitle} \nvs all other clusters"), y =2),
               xaxis = list(title = "Log2(Fold Change)"),
               yaxis = list(title = "-Log10(BH.q.val)"),
               annotations = 
                 list(x = 1.57, y = -0.06, #position of text adjust as needed 
                      text = qq("BH.q.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1 \nDifference between groups calculated using ks.test()"), 
                      showarrow = F, xref='paper', yref='paper', 
                      xanchor='right', yanchor='auto', xshift=0, yshift=0,
                      font=list(size=5, color="black"))) %>%
        rangeslider(); plotly_p
      
      # static ggplot pdf
      p <- ggplot(to_plot) + 
        geom_point(mapping = aes(x = logFC, y = `-log(p_val_bh)`, size = `-log(p_val_bh)`, color = `-log(p_val_bh)`))  +
        geom_label_repel(aes(x = logFC, y = `-log(p_val_bh)`, family="Courier", label = ifelse(signif, as.character(plot_analyte_name), "")),
                         box.padding = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') +
        scale_color_viridis(direction = 1) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 10)) +
        guides(size = FALSE) + 
        labs(color = "-log10(BH.p.val)", 
             caption = "BH.q.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1 \nDifference between groups calculated using ks.test()") +
        ylab("-Log10(BH.q.val)") +
        xlab("Log2(Fold Change)") + 
        ggtitle(label = qq("@{toupper(which.dat)}, @{dname}"),
                subtitle = qq("@{plot_subtitle} vs all other clusters"))
      
      
      if (nrow(which_cluster_vasc_df) > 1 | clust_type == "all") {
        # message("Building facet...")
        # do a facet plot for the different clusters if huvecs and smcs clustered into different groups
        p <- p + facet_grid(cols = vars(name)); p
        
        # have to build new plot for plotly conversion, without geom_repel
        # include aesthetics in main call to ggplot
        p_for_plotly <- ggplot(to_plot, aes(x = logFC, y = `-log(p_val_bh)`, size = `-log(p_val_bh)`, color = `-log(p_val_bh)`, label = plot_analyte_name)) + 
          geom_point() +
          scale_color_viridis(direction = 1) + 
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 10)) +
          guides(size = FALSE) + 
          labs(color = "-log10(BH.p.val)", 
               caption = "BH.q.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1 \nDifference between groups calculated using ks.test()") +
          ylab("-Log10(BH.q.val)") +
          xlab("Log2(Fold Change)") + 
          ggtitle(label = qq("@{toupper(which.dat)}, @{dname}"),
                  subtitle = qq("@{plot_subtitle} vs all other clusters")) + 
          facet_grid(cols = vars(name))
        plotly_p <- suppressWarnings(ggplotly(p_for_plotly, tooltip = c("label", "x", "y"))); plotly_p
      } 
      
      output_dir <- file.path(base_output_dir, "diff", qq("@{clust_type}-h_@{dendro_cut_thresh}-q_@{bh_thresh}"))
      dir.create(file.path(output_dir, "CO_CLUST"), recursive = T, showWarnings = F)
      dir.create(file.path(output_dir, "OTHER"), recursive = T, showWarnings = F)
      add <- ifelse(co_clust, "CO_CLUST/","OTHER/")
      
      plt_name <- qq("@{add}@{dname}_diff_exp__q_@{bh_thresh}")
      saveWidget(plotly_p, file = file.path(output_dir, qq("@{plt_name}.html")), selfcontained = T)
      # message(qq("Writing plot html and pdf to @{file.path(output_dir,plt_name)}"))
      ggexport(p, filename = file.path(output_dir, qq("@{plt_name}.pdf")), width = 10, height = 12)
      
      # return(p)
    })
  
}                       

# helper 
plotly_vline <- function(x = 0, name = "test line", color = "red", dash = "dash") {
  list(name = name,
       type = "line", 
       y0 = 0, 
       y1 = 1, 
       yref = "paper",
       x0 = x, 
       x1 = x, 
       line = list(color = color, dash = dash)
  )
}

plotly_hline <- function(y = 0, name = "test line", color = "blue", dash = "dash") {
  list(name = name,
       type = "line", 
       x0 = 0, 
       x1 = 1, 
       xref = "paper",
       y0 = y, 
       y1 = y, 
       line = list(color = color, dash = dash)
  )
}






# p <- plot_diff_analyte_results_CVs_html(diffe_by_clust = diffe_by_clust_w_phos, 
#                                         which.dat = which.dat, 
#                                         gid = c("HUVEC", "HAoSMC"), 
#                                         dname = dname, 
#                                         cancer_vs_non_cancer = cancer_vs_non_cancer, 
#                                         base_output_dir = base_output_dir, 
#                                         dendro_cut_thresh = dendro_cut_thresh, bh_thresh = bh_thresh)

#############################################################################################

#' #' @note plot differential analyte results using a faceted volcano plot 
#' #' AND displays most differential analytes
#' #' @param diffe_by_clust differential analyte results with cluster assignments
#' #' @param which.dat which dataset is primary
#' #' @param cluster_assignments_df tbl of cluster assignments
#' #' @param dname grouping variable, e.g. name of a perturbation ID (passed in here for labeling the graph)
#' #' @param bh_thresh bh threshold to plot analyte names (passed in here for labeling the graph)
#' plot_diff_analyte_results_CVs <- function(diffe_by_clust, which.dat, cluster_assignments_df, 
#'                                           dname, BH_THRESH, gid = c("HUVEC", "HAoSMC")){
#'   if (is.null(diffe_by_clust)) return("No clustering was performed for this condition.")
#'   # diffe_by_clust = test$diff_ex[[1]]; which.dat = "p100"; cluster_assignments_df = test$ca_df[[1]]
#'   
#'   ca_df <- cluster_assignments_df %>% 
#'     filter(cid %in% gid)
#'   
#'   clust_label <- cluster_assignments_df %>% 
#'     group_by(cluster) %>%
#'     summarize(name = str_c(cid, collapse = ",")) %>%
#'     ungroup() %>%
#'     rename(base_clust_comp = cluster) %>%
#'     semi_join(ca_df, by=c("base_clust_comp" = "cluster")) %>%
#'     distinct(name, base_clust_comp)
#'   
#'   which_cluster_vasc_df <- ca_df %>%
#'     distinct(cluster)
#'   
#'   to_plot <- diffe_by_clust %>% 
#'     semi_join(which_cluster_vasc_df, by=c("base_clust_comp" = "cluster")) %>%
#'     left_join(clust_label, by="base_clust_comp") %>%
#'     mutate(`-log(p_val_bh)` = -log(p_val_bh,base = 10)) %>%
#'     na.omit()
#'   
#'   min_max <- to_plot %>% 
#'     arrange(p_val_bh,logFC) %>%
#'     filter(signif)
#'   
#'   if (nrow(min_max) > 1){
#'     top_analytes <- tibble()
#'     for (i in unique(min_max$base_clust_comp)){
#'       df <- filter(min_max, base_clust_comp == i)
#'       max <- df[which.max(df$logFC),] %>% mutate(status = "Max logFC")
#'       min <- df[which.min(df$logFC),] %>% mutate(status = "Min logFC")
#'       top_analytes <- bind_rows(top_analytes, max, min)
#'     }
#'     tbl <- top_analytes %>% 
#'       select(`cluster` = name, analyte, logFC, `BH.p.val` = p_val_bh, status)
#'     
#'     # cat("\n\n\\pagebreak\n")
#'     print(kable(tbl,
#'                 booktabs=TRUE, caption = "Maximally Differential Analytes",align = 'r') %>%
#'             kable_styling(latex_options = "hold_position",font_size = 7, full_width = F))
#'   }
#'   
#'   LOG_PLOT = TRUE
#'   
#'   if (LOG_PLOT){
#'     plot_subtitle <- str_c(clust_label$name, collapse=' and ')
#'     p <- ggplot(to_plot) + 
#'       geom_point(mapping = aes(x = logFC, y = `-log(p_val_bh)`, size = `-log(p_val_bh)`, color = `-log(p_val_bh)`))  +
#'       geom_label_repel(aes(x = logFC, y = `-log(p_val_bh)`, family="Courier", label = ifelse(signif, as.character(analyte), "")),
#'                        box.padding = 0.35, 
#'                        point.padding = 0.5,
#'                        segment.color = 'grey50') +
#'       scale_color_viridis(direction = 1) + 
#'       theme_bw() +
#'       theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
#'             plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
#'             axis.title = element_text(size = 15),
#'             axis.text = element_text(size = 10)) +
#'       guides(size = FALSE) + 
#'       labs(color = "-log10(BH.p.val)", 
#'            caption = "BH.p.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1") +
#'       ylab("-Log10(BH.p.val)") +
#'       xlab("Log2(Fold Change)") + 
#'       ggtitle(label = qq("@{toupper(which.dat)}"),
#'               subtitle = qq("@{plot_subtitle}"))
#'   } else {
#'     plot_subtitle <- str_c(clust_label$name, collapse=' and ')
#'     p <- ggplot(to_plot) + 
#'       geom_point(mapping = aes(x = logFC, y = `p_val_bh`, size = `p_val_bh`, color = `p_val_bh`))  +
#'       geom_label_repel(aes(x = logFC, y = `p_val_bh`, family="Courier", label = ifelse(signif, as.character(analyte), "")),
#'                        box.padding = 0.35, 
#'                        point.padding = 0.5,
#'                        segment.color = 'grey50') +
#'       scale_color_viridis(direction = 1) + 
#'       theme_bw() +
#'       theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
#'             plot.subtitle = element_text(hjust = 0.5, size = 15, face = "bold"),
#'             axis.title = element_text(size = 15),
#'             axis.text = element_text(size = 10)) +
#'       guides(size = FALSE) + 
#'       labs(color = "BH.p.val", 
#'            caption = "BH.p.val = Benjamini-Hochberg-Corrected P-value (q) \nCut-off for displaying label: q < 0.1") +
#'       ylab("BH.p.val") +
#'       xlab("Log2(Fold Change)") + 
#'       ggtitle(label = qq("@{toupper(which.dat)}"),
#'               subtitle = qq("@{plot_subtitle}"))
#'   }
#'   
#'   if (nrow(which_cluster_vasc_df) > 1) {
#'     # do a facet plot for the different clusters if huvecs and smcs clustered into different groups
#'     return(p + facet_grid(cols = vars(name)))
#'   } else{
#'     return(p)
#'   }
#' }
#' 
#' #' @note plot differential analyte results using a faceted volcano plot
#' #' @param diffe_by_clust differential analyte results with cluster assignments
#' #' @param cluster_assignments_df tbl of cluster assignments
#' #' @param dname grouping variable, e.g. name of a perturbation ID (passed in here for labeling the graph)
#' #' @param BH_THRESH bh threshold to plot analyte names (passed in here for labeling the graph)
#' plot_diff_analyte_results <- function(diffe_by_clust, cluster_assignments_df, plot_only_vascular = cancer_vs_non_cancer, dname, BH_THRESH){
#'   # do a facet plot !!
#'   message("Plotting analytes that are differentially enriched compared to all other clusters")
#'   g_names <- cluster_assignments_df %>%
#'     group_by(cluster) %>%
#'     summarize(name = str_c(cid, collapse = ",")) %>%
#'     ungroup() %>%
#'     rename(base_clust_comp = cluster)
#'   
#'   
#'   if (plot_only_vascular){
#'     # plotting vascular only!! if cancer vs vascular, we only need one plot
#'     g_names <- g_names %>% filter(name == "HAoSMC,HUVEC")
#'     diffe_by_clust <- diffe_by_clust %>% semi_join(g_names)
#'   }
#'   
#'   to_plot <- diffe_by_clust %>% left_join(g_names,by="base_clust_comp") 
#'   p <- ggplot(to_plot) + 
#'     geom_point(mapping = aes(x = logFC, y = -log(p_val), color = -log(p_val_bh))) + 
#'     facet_grid(name ~ .) +
#'     geom_label_repel(data = to_plot %>% filter(signif),
#'                      aes(x = logFC, y = -log(p_val), family="Courier", 
#'                          label = ifelse(signif, unique(as.character(analyte)), "")),
#'                      box.padding = 0.35, 
#'                      point.padding = 0.5,
#'                      segment.color = 'grey50') +
#'     scale_color_viridis(direction = 1) + 
#'     theme_bw() + # base_size = 4, base_family = "Noto Sans"
#'     ggtitle(qq("@{dname}, p.val cutoff = @{BH_THRESH}"))
#'   return(p)
#' }
#' 