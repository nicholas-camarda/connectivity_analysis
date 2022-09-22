#' @note creates side by side plots of pvclust results, then plots differential analyte results. For use in the RMarkdown document
#' @note I realize this is very messy code. i'm working on it :(
#' @param lst_obj list of diff_ex obj produced in rmarkdown document. [MUST BE NAMED!!!]
#' @param which.dat c(p100, gcp, avg) for organizing which dataset results to plot. If CV cell lines clustered together for perturbations in P100, then the clustering for those perturbations are compared to the other datasets  
#' @param filter_co_clust if TRUE, plot just the conditions in which CV cell lines co-clustered. FALSE for "ALL" condition.
#' @param plot_coclust_diff_ex if FALSE, just plot the main diff_ex volcano plot (if which.dat == p100, then just diff_ex for p100 is plotted). If TRUE, plot for all datasets in which CVs co-clustered for that perturbation group
plot_everything <- function(lst_obj, which.dat, filter_co_clust = TRUE, plot_coclust_diff_ex = FALSE,
                            bh_thresh = 0.1, dendro_thresh = 0.6){
  
  main_obj_idx <- which(names(lst_obj) == which.dat)
  mapping <- seq(1:3)
  names(mapping) <- names(lst_obj)
  
  pick <- setdiff(seq(1:3), main_obj_idx)
  if (filter_co_clust){
    obj_1 <- lst_obj[[main_obj_idx]] %>% filter(co_clust) %>% arrange(pert_iname)
  } else {
    obj_1 <- lst_obj[[main_obj_idx]] %>% arrange(pert_iname)
  }
  
  obj_2 <- lst_obj[[pick[1]]] %>% filter(pert_iname %in% obj_1$pert_iname) %>% arrange(pert_iname)
  obj_3 <- lst_obj[[pick[2]]] %>% filter(pert_iname %in% obj_1$pert_iname) %>% arrange(pert_iname)
  
  f_obj_lst <- list(obj_1, obj_2, obj_3)
  names(f_obj_lst) <- c(names(mapping[main_obj_idx]), names(mapping[pick[1]]), names(mapping[pick[2]]))
  
  # perts_1 always main
  for (i in 1:nrow(obj_1)){
    # i = 3
    if (is.null(obj_1$pvclust[[i]])) next
    
    par(mfrow = c(1,3))
    dname <- obj_1$pert_iname[[i]]; dname
    cat(qq("### @{toupper(dname)}\n\n"))
    cat(qq("\n\n PVCLUST Results \n\n"))
    
    # thresh = dendrogram cut-off threshold
    cls_1 <- obj_1 %>% filter(pert_iname == dname) 
    cls_2 <- obj_2 %>% filter(pert_iname == dname); if (nrow(cls_2) < 1) next 
    cls_3 <- obj_3 %>% filter(pert_iname == dname); if (nrow(cls_3) < 1) next 
    
    display_pvclust(cls_1$pvclust[[1]], dname = dname, dataset = toupper(names(mapping[main_obj_idx])), 
                    co_clust = cls_1$co_clust[[1]], thresh = dendro_thresh)
    display_pvclust(cls_2$pvclust[[1]], dname = dname, dataset = toupper(names(mapping[pick[1]])), 
                    co_clust = cls_2$co_clust[[1]], thresh = dendro_thresh)
    display_pvclust(cls_3$pvclust[[1]], dname = dname, dataset = toupper(names(mapping[pick[2]])), 
                    co_clust = cls_3$co_clust[[1]], thresh = dendro_thresh)
    

    par(mfrow = c(1,1))
    if (!filter_co_clust){
      # this only happens for "ALL" - we want to show the diff_ex analysis for whichever dataset is main, regardless of co-cluster
      g <- plot_diff_analyte_results_CVs(diffe_by_clust = obj_1$diff_ex[[i]], which.dat = which.dat,
                                     cluster_assignments_df = obj_1$ca_df[[i]], 
                                     dname = dname, BH_THRESH = bh_thresh) 
      if (!is.null(g)) {
        cat(qq("\n\n Differential Analytes Results \n\n"))
        print(g)
        cat("\n\n\\pagebreak\n")
      }
      next # this is important -- we don't want to plot anything else! Just for our main dataset. The others will be plotted later
    }

    # get boolean vector for plotting co-clust pairs
    vec <- c(obj_1[i,]$co_clust, obj_2[i,]$co_clust, obj_3[i,]$co_clust)
    for (j in 1:length(vec)){
      oj <- f_obj_lst[[j]] # for the jth dataset object
      if (vec[j]){
        # if co-clust, plot the ith perturbation condition
        g <- plot_diff_analyte_results_CVs(diffe_by_clust = oj$diff_ex[[i]], which.dat = names(f_obj_lst[j]),
                                       cluster_assignments_df = oj$ca_df[[i]], 
                                       dname = dname, BH_THRESH = bh_thresh)
        if (!is.null(g)) {
          cat(qq("\n\n Differential Analytes Results \n\n"))
          print(g)
          cat("\n\n ")
        }
      } else { next }
    }
    
    cat("\n\n\\pagebreak\n")
  }
}

#' @note creates side by side plots of pvclust results, then plots differential analyte results. For use in the RMarkdown document
#' @note I realize this is very messy code. i'm working on it :(
#' @param lst_obj list of diff_ex obj produced in rmarkdown document. [MUST BE NAMED!!!]
#' @param which.dat c(p100, gcp, avg) for organizing which dataset results to plot. If CV cell lines clustered together for perturbations in P100, then the clustering for those perturbations are compared to the other datasets  
#' @param plot_coclust_diff_ex if FALSE, just plot the main diff_ex volcano plot (if which.dat == p100, then just diff_ex for p100 is plotted). If TRUE, plot for all datasets in which CVs co-clustered for that perturbation group
plot_poster_figures <- function(lst_obj, which.dat="p100", plot_coclust_diff_ex = FALSE,
                            bh_thresh = 0.1, dendro_thresh = 0.6){
  
  obj_1 <- lst_obj[[which.dat]] %>% arrange(pert_iname)
 
  # for every drug, moa, etc..
  for (i in 1:nrow(obj_1)){
    if (is.null(obj_1$pvclust[[i]])) return(NULL)
    par(mfrow = c(1,1))
    dname <- obj_1$pert_iname[[i]]; dname
    
    cat(qq("### @{toupper(dname)}\n\n"))
    cat(qq("\n\n PVCLUST Results \n\n"))
    
    cls_1 <- obj_1 %>% filter(pert_iname == dname) 
    
    display_pvclust_poster(cls_1$pvclust[[1]], ct = cls_1$cut_trees[[1]], dname = dname, dataset = toupper(which.dat))

    # this only happens for "ALL" - we want to show the diff_ex analysis for whichever dataset is main, regardless of co-cluster
    g <- plot_diff_analyte_results_CVs(diffe_by_clust = cls_1$diff_ex[[1]], which.dat = which.dat,
                                       cluster_assignments_df = cls_1$ca_df[[1]], 
                                       dname = dname, BH_THRESH = bh_thresh) 
    if (!is.null(g)) {
      cat(qq("\n\n Differential Analytes Results \n\n"))
      print(g)
      cat("\n\n\\pagebreak\n")
    }
  }
}

#' @note plot summary of drugs used
#' @param fn_name name of the file containing drug/class info
summarize_drugs <- function(fn_name = "common_drugs_classes.tsv", 
                            dataset = "P100",
                            REFERENCES_DIRECTORY = file.path("/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/Cardio-oncology/References")){
  drugs <- read_tsv(file.path(REFERENCES_DIRECTORY, fn_name)) %>%
    arrange(drug_class) %>%
    rename(Perturbation = pert_iname, `Drug Class` = drug_class) 
  sum_drugs <- drugs %>% 
    group_by(`Drug Class`) %>% 
    summarize(Count = n()) %>% 
    mutate(Percent = percent(Count/sum(Count)))
  
  print(kable(drugs, booktabs=TRUE, caption = qq("Summary of perturbation distribution for @{dataset}")) %>%
    kable_styling(latex_options = "hold_position"))
  return(list("drugs" = drugs, "sum_drugs" = sum_drugs))
}

#' @note for quickly expanding a single string into a vector
check_list <- function(str1){
  return(str_trim(str_split(str1, ",", simplify = T)))
}

#' @note 
#' @param obj_dat takes, for example, perts$p100
extract_cell_lines_per_pert <- function(obj_dat, dataset = "p100"){
  res <- obj_dat %>% 
    select(pert_iname, df) %>% 
    mutate(cell_lines = map(obj_dat$df, function(x) sort(unique(sapply(x$cell_id, function(y) str_split(y, "\\.",simplify = T)[,1]))))) %>%
    group_by(pert_iname) %>% 
    summarize(`Cell Lines` = map(cell_lines, function(x) str_c(x, collapse = ", "))) %>% 
    select(`Perturbation`=pert_iname, `Cell Lines`) %>%
    unnest()
  
  res_tbl <- kable(res, booktabs=TRUE, caption = qq("Summary of cell line distribution for @{dataset}")) %>%
          kable_styling(latex_options = "hold_position")
  
  return(res_tbl)
}
