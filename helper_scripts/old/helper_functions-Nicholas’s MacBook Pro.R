# install.packages(c("tidyverse","readxl","pvclust","gplots","viridis",
#                    "GetoptLong","gridExtra","ggrepel","ggrepel","doMC",
#                    "conflicted"))

##############################################
### Functions ################################
##############################################

# load gcp
load_gcp <- function(first_gen_datasets = "/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/Cardio-oncology/Datasets/1st gen data",
                     gcp_dir = "GCP"){
  gcp_f <- file.path(first_gen_datasets, gcp_dir, "GCP All Cell Lines.gct")
  
  gcp_r <- read_tsv(gcp_f, skip = 2)
  gcp_a <- gcp_r[c(1,4,10,14,21:nrow(gcp_r)),c(1,4,9:ncol(gcp_r))] %>%
    mutate(pr_gene_symbol = make.unique(pr_gcp_histone_mark,sep="_")) %>%
    select(-pr_gcp_histone_mark) %>%
    select(id, pr_gene_symbol,everything())
  return(gcp_a)
}

# load p100
load_p100 <- function(first_gen_datasets, p100_dir = "P100"){
  
  p100_f <- file.path(first_gen_datasets, p100_dir,"P100 All Cell Lines.gct")
  p100_r <- read_tsv(p100_f, skip = 2)
  p100a <- p100_r[c(1,3,11,15,23:nrow(p100_r)),c(1,3,13:ncol(p100_r))] %>%
    mutate(pr_gene_symbol = make.unique(pr_gene_symbol,sep="_")) # remove metadata columns/rows
  return(p100a)
}


# load data
load_data <- function(data,
                      drug_classes_fn = "/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/Cardio-oncology/References/Drug Glossary_edited.xlsx",
                      directory = "P100"){

  # missing_data_plot_dir <- file.path(first_gen_datasets, directory, "plots", "missing")
  # dir.create(missing_data_plot_dir,recursive = T)
  # 
  
  #' @note reading in the data
  #' @note more info about the .GCT file format here: https://clue.io/connectopedia/gct_format
  
  
  cancer_drug_moa_df <- read_excel(path = drug_classes_fn,sheet = 1) %>% filter(!is.na(Drug)) %>% transmute(pert_iname = Drug, drug_class = Class, drug_cat = "cancer"); dim(cancer_drug_moa_df)
  cv_drug_moa_df <- read_excel(path = drug_classes_fn,sheet = 2) %>% filter(!is.na(`Drug (Generic)`)) %>% transmute(pert_iname = `Drug (Generic)`, drug_class = `Class`, drug_cat = "cv"); dim(cv_drug_moa_df)
  drugs_moa_df <- bind_rows(cancer_drug_moa_df,cv_drug_moa_df) %>%
    mutate(pert_iname = tolower(pert_iname))
  
  # in GCTs, Typically, each column represents a specific experiment (e.g. treatment of cell line MCF7 with a small-molecule drug) 
  # and each row represents features (e.g. genes) that are measured in the assay.
  
  wells <- colnames(data)[3:ncol(data)]
  cell_ids <- data[1,] %>%
    gather(well, cell_id, wells) %>%
    select(well, cell_id) 
  
  cell_cat <- cell_ids %>% 
    distinct(cell_id) %>%
    mutate(cell_cat = ifelse(cell_id %in% c("HUVEC","HAoSMC"), "cv","cancer"))
  
  # map pert_id to well number
  det_norm_group_vect <- data[2,] %>%
    gather(well, det_normalization_group_vector, wells) %>%
    select(well, det_normalization_group_vector) 
  
  replicates <- data[3,] %>%
    gather(well, pert_batch_internal_replicate, wells) %>%
    select(well, pert_batch_internal_replicate)
  
  # map pert_id to well number
  pert_names <- data[4,] %>%
    gather(well, pert_iname, wells) %>%
    select(well, pert_iname) %>% 
    mutate(pert_iname = tolower(pert_iname))
  
  # # just for safe keeping, in case we want later
  # phosphosite_long_df <- data %>% 
  #   slice(-c(1,2,3)) %>% 
  #   select(pr_gene_symbol,pr_p100_phosphosite, `PM7-46D43-001A01`:`PHC-44828-080H10`)
  
  # make "long" dataset - gene, well, cell_id, expression
  datb <- data %>% 
    slice(-c(1,2,3,4)) %>% # take out cell_id, replicate, group_vect, and pert_iname
    select(pr_gene_symbol,wells) %>%
    gather(well, pex, wells) %>%
    left_join(cell_ids, by="well") %>% # add cell_Id back in correctly
    left_join(pert_names, by="well") %>%
    left_join(replicates, by="well") %>%
    left_join(det_norm_group_vect, by="well") %>%
    left_join(cell_cat, by="cell_id") %>%
    arrange(cell_id,pr_gene_symbol,pert_iname) %>%
    mutate(pex = as.numeric(pex)) %>%
    rename(grp = det_normalization_group_vector,
              rep = pert_batch_internal_replicate)
  
  
  # take median of replicates for simplified dataset
  datbb <- datb %>%
    group_by(pr_gene_symbol,cell_id,pert_iname) %>% 
    dplyr::summarize(pex_m = median(pex, na.rm = T)) %>% # median to collapse replicates
    ungroup() %>%
    left_join(drugs_moa_df, by="pert_iname")
  # p100c %>% distinct(pr_gene_symbol) %>% write_tsv("/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/ws/clustering/proteins.tsv")
  
  tested_drugs_in_cv_cells <- datbb %>% select(-drug_class, -drug_cat) %>% 
    filter(cell_id %in% c("HUVEC","HAoSMC")) %>% 
    spread(pert_iname, pex_m) %>%
    group_by(cell_id) %>%
    nest() %>%
    mutate(tested_drugs = map(data, function(x) colnames(x[,!apply(is.na(x), 2, all)])[-1]))
  
  perturbation_feature_set <- intersect(tested_drugs_in_cv_cells$tested_drugs[[1]],tested_drugs_in_cv_cells$tested_drugs[[2]])
  dir.create("/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/ws/feature_sets",recursive = T)
  write_tsv(tibble(perturbation_feature_set),'/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/ws/feature_sets/pert_set.tsv')
  feature_set <- datbb %>% distinct(pr_gene_symbol) %>% .$pr_gene_symbol
  write_tsv(tibble(feature_set),'/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/ws/feature_sets/protein_set.tsv')
  
  cell_lines_identities <- datbb %>% 
    distinct(cell_id) %>% 
    mutate(is_cancer = !(cell_id %in% c("HUVEC","HAoSMC")))
  
  # give cell line ids to final dat
  datc <- datbb %>% 
    filter(pert_iname %in% perturbation_feature_set) %>%
    left_join(cell_lines_identities, by="cell_id")
  
  # give cell line ids, drug MOAs, to datb
  dat_with_reps <- datb %>% 
    filter(pert_iname %in% perturbation_feature_set) %>%
    left_join(cell_lines_identities, by="cell_id") %>%
    left_join(drugs_moa_df, by="pert_iname")
  
  return(list(data = datc, data_with_reps = dat_with_reps, feature_set = feature_set, perturbation_feature_set = perturbation_feature_set))
}


#' @note Remove columns/rows/both that have an NA percentage greater than or equal to the specified threshold
#' @param dat matrix
#' @param thresh_row percent of NA values along rows that you want to remove
#' @param thresh_col percent of NA values along cols that you want to remove
#' @param d rows and columns? possible values include c("both", "cols", "rows")
remove_x_perc_NA <- function(dat, thresh_row = 1.0, thresh_col = 1.0, d){
  if (as.logical(match(d, "both"))){
    # message("both")
    rsd <- which(rowMeans(is.na(dat)) >= thresh_row);# message("Removing: ", length(rsd), " rows\n")
    rs <- setdiff(1:nrow(dat),rsd)
    csd <- which(colMeans(is.na(dat)) >= thresh_col); cs_n <- paste(names(csd), sep="", collapse=" "); #message(qq("Removing cols: @{cs_n}"))
    cs <- setdiff(1:ncol(dat),csd)
    new <- as.matrix(dat[rs,cs,drop=F])
    return(list(mat = new, removed = list(cols = names(csd), rows = rsd)))
  } else if (as.logical(match(d, "cols"))){
    # message("cols")
    csd <- which(colMeans(is.na(dat)) >= thresh_col); cs_n <- paste(names(csd), sep=" ", collapse=""); #message(qq("Removing cols: @{cs_n}"))
    cs <- setdiff(1:ncol(dat),csd)
    new <- as.matrix(dat[,cs,drop=F])
    return(list(mat = new, removed = list(cols = names(csd), rows = NA)))
  } else {
    # message("rows")
    rsd <- which(rowMeans(is.na(dat)) >= thresh_row); #message("Removing: ", length(rsd), " rows")
    rs <- setdiff(1:nrow(dat),rsd)
    new <- as.matrix(dat[rs,,drop=F])
    return(list(mat = new, removed = list(cols = NA, rows = rsd)))
  }
}



#' @param mat_labs long data frame of median phosphorylation values
#' @param p_base_output_dir base directory for output
#' @param dataset p100 or gcp, just for naming output specific to dataset
#' @param feature_set set of proteins or set of histone markers
cluster_analysis_pert <- function(mat_labs, 
                                  p_base_output_dir = "/Users/Nicholas/OneDrive - Tufts/phd/labs/jaffe/workspace/ws/clustering", 
                                  dataset = "p100",
                                  feature_set = p100_feature_set){
  
  base_output_dir <- file.path(p_base_output_dir, dataset)
  message(qq("Writing to @{base_output_dir}"))
  
  methods <- c( "average", "single", "complete", "ward")
  names(methods) <- c( "average", "single", "complete", "ward")
  
  # convert data to a "data frame), impute median of column for missing values and set up for clustering methods that best suit the data
  # do I want to just do the same clustering method for all of them?? probably...
  # also want to remove NA clumns
  message("\nPreparing data for analysis....")
  clust <- mat_labs %>% 
    group_by(pert_iname) %>%
    nest(cell_id, feature_set) %>%
    mutate(dataframe = map(data, function(x) {
      x_df <- as.data.frame(x)
      rownames(x_df) <- x[,1,drop=T]
      x_df_c <- x_df[,-1]
      x_df_c_na <- remove_x_perc_NA(dat = x_df_c, thresh_row = 1.0, thresh_col = 0.5, d = "both")$mat
      return(x_df_c_na)
    })) 
  
  # correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables
  # Kendall -- non-parametric test that measures the strength of dependence between two variables. 
  # https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html#methods-for-correlation-analyses
  
  message("\nComputing correlation matrices off of raw data....")
  clust2 <- clust %>% 
    mutate(correlation = map(dataframe, compute_correlation),
           correlation_cell = map(dataframe, function(x) {
             t_d <- t(dataframe)
             t_c <- compute_correlation(t_d)
           }))
  
  
  # generate heatmaps!
  
  message("\nComputing connectivities between cells.. and clustering using pvclust...")
  clust3 <- clust2 %>% 
    mutate(cell_conn = map(correlation_cell, compute_connectivity_archive)) %>%
    mutate(cell_conn_mat = map(cell_conn, make_numeric_mat)) %>%
    mutate(pvclust_cell = map2(cell_conn_mat, pert_iname, .f = compute_boot_pvclust, base_output_dir = base_output_dir))
  
  walk2(clust3$cell_conn_mat, clust3$pert_iname, plot_morpheus, base_output_dir = base_output_dir, subdir = "cell_conn", minr = -1, maxr = 1, sym = TRUE)
  
  message("\nPerforming kmeans")
  clust3 <- clust3 %>% 
    mutate(kmeans_cell_conn = map2(cell_conn_mat, clust3$pert_iname, .f = do_kmeans, subdir = "kmeans_cell_conn"))
  
  # write out clust for safe keeping
  write_rds(clust3, path = file.path(base_output_dir,qq("@{dataset}_clust3.rds")))
  return(clust3)
}


## creating and saving MORPHEUS objects
# from https://github.com/cmap/morpheus.R/blob/master/R/morpheus.R

dendToTree <- function(dend) {
  tree <- c(
    as.list(attributes(dend)[c('height')])
  )
  
  # Recursively add children
  if (! is.leaf(dend)) {
    tree$children <- lapply(dend, dendToTree)
  }
  tree
}

is.dendrogram <- function (x) { inherits(x, "dendrogram")}


create.payload <- function(x,
                           labRow = rownames(x),
                           labCol = colnames(x),
                           Rowv = TRUE,
                           Colv=if (symm)"Rowv" else TRUE,
                           distfun = dist,
                           hclustfun = hclust,
                           dendrogram = c("both", "row", "column", "none"),
                           reorderfun = function(d, w) reorder(d, w),
                           symm = FALSE,
                           na.rm = TRUE,
                           rowAnnotations=NULL,
                           columnAnnotations=NULL,
                           colorScheme = NULL,
                           rowSize = 13,
                           columnSize = 13,
                           drawGrid = TRUE,
                           gridColor = "#808080",
                           gridThickness = 0.1,
                           drawValues = FALSE,
                           width = NULL,
                           height = NULL) {
  name <- deparse(substitute(x))
  ## x is a matrix!
  if (! is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (! is.matrix(x)) stop("x must be a matrix")
  options(expressions= 500000)
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  ddc <- NULL
  ddr <- NULL
  if (! inherits(Rowv, "dendrogram")) {
    if (((is.logical(Rowv) && ! isTRUE(Rowv)) || (is.null(Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
      if (dendrogram == "both")
        dendrogram <- "column"
      else dendrogram <- "none"
    }
  }
  if (! inherits(Colv, "dendrogram")) {
    if (((is.logical(Colv) && ! isTRUE(Colv)) || (is.null(Colv))) &&
        (dendrogram %in% c("both", "column"))) {
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
      if (dendrogram == "both")
        dendrogram <- "row"
      else dendrogram <- "none"
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
    if (length(rowInd) > nr || any(rowInd < 1 | rowInd > nr))
      stop("Rowv dendrogram doesn't match size of x")
    if (length(rowInd) < nr)
      nr <- length(rowInd)
  }
  else if (is.integer(Rowv)) {
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (! isTRUE(Rowv)) {
    rowInd <- nr : 1
    ddr <- as.dendrogram(hclust(dist(diag(nr))))
  }
  else {
    rowInd <- nr : 1
    ddr <- as.dendrogram(Rowv)
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
    if (length(colInd) > nc || any(colInd < 1 | colInd > nc))
      stop("Colv dendrogram doesn't match size of x")
    if (length(colInd) < nc)
      nc <- length(colInd)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    distc <- distfun(if (symm)
      x
      else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    distc <- distfun(if (symm)
      x
      else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (! isTRUE(Colv)) {
    colInd <- 1 : nc
    ddc <- as.dendrogram(hclust(dist(diag(nc))))
  }
  else {
    colInd <- 1 : nc
    ddc <- as.dendrogram(Colv)
  }
  
  ddr <- rev(ddr)
  rowInd <- rev(rowInd) # reverse to match order of R heat maps
  if(!is.null(rowAnnotations)){
    rowAnnotations <- rowAnnotations[rowInd,]
  }
  if(!is.null(columnAnnotations)){
    columnAnnotations <- columnAnnotations[colInd,]
  }
  
  ## Labels for Row/Column
  rownames(x) <- labRow %||% paste(1 : nrow(x))
  colnames(x) <- labCol %||% paste(1 : ncol(x))
  options(htmlwidgets.TOJSON_ARGS = list(dataframe = "column"))
  morpheusOptions <- list()
  morpheusOptions$colorScheme <- colorScheme
  morpheusOptions$rowSize <- rowSize
  morpheusOptions$columnSize <- columnSize
  morpheusOptions$drawGrid <- drawGrid
  morpheusOptions$gridColor <- gridColor
  morpheusOptions$gridThickness <- gridThickness
  morpheusOptions$drawValues <- drawValues
  
  x <- x[rowInd, colInd]
  if (!is.null(morpheusOptions$colorScheme$colors)) {
    morpheusOptions$colorScheme$colors <- lapply(morpheusOptions$colorScheme$colors, function(color){
      rgb <- col2rgb(color)
      paste("rgb(", rgb[1], ",", rgb[2], ",", rgb[3], ")", sep='')
    })
  }
  
  if (!is.null(morpheusOptions$colorScheme$colors) && is.null(morpheusOptions$colorScheme$values)) {
    rng = range(x)
    nvals <- length(morpheusOptions$colorScheme$colors)
    fractionStep <- 1/(nvals-1)
    values <- vector("list", nvals)
    dataRange <- rng[2] - rng[1]
    values[1] <- rng[1]
    values[nvals] <- rng[2]
    
    for(i in 2:nvals-1) {
      fraction <-fractionStep*(i-1)
      values[i] <- rng[1] + fraction*dataRange
    }
    morpheusOptions$colorScheme$values <- values
  }
  
  columnDendrogram <- if (!is.null(ddc) &&
                          is.dendrogram(ddc) &&
                          dendrogram %in% c("both", "column")) dendToTree(ddc) else NULL
  rowDendrogram <- if (!is.null(ddr) &&
                       is.dendrogram(ddr) &&
                       dendrogram %in% c("both", "row")) dendToTree(ddr) else NULL
  
  
  
  rowVectors <- list()
  rowVectors[[1]] = list(name='id', array=rownames(x))
  
  if(!is.null(rowAnnotations)) {
    for (i in 1:ncol(rowAnnotations)) {
      rowVectors[[i+1]] = list(name= names(rowAnnotations)[i], array=rowAnnotations[,i])
    }
  }
  rowMetadataModel <- list(vectors=rowVectors)
  columnVectors <- list()
  columnVectors[[1]] = list(name='id', array=colnames(x))
  if(!is.null(columnAnnotations)) {
    for (i in 1:ncol(columnAnnotations)) {
      columnVectors[[i+1]] = list(name= names(columnAnnotations)[i], array=columnAnnotations[,i])
    }
  }
  columnMetadataModel <- list(vectors=columnVectors)
  dataset <- list(seriesNames=list(name), rows = nrow(x), columns = ncol(x),  seriesDataTypes=list('number'),
                  seriesArrays=list(x), rowMetadataModel=rowMetadataModel, columnMetadataModel=columnMetadataModel)
  
  morpheusOptions$dataset = dataset
  morpheusOptions$name = name
  payload <- list(rowDendrogram = rowDendrogram, columnDendrogram = columnDendrogram, options=morpheusOptions)
  return(payload)
}



## my code to extract

input_null_for_NA <- function(x){
  x %>% mutate_all(~replace(., is.na(.), "null"))
}

get_morpheus_json_object_from_drug_name <- function(lst_object, drug_name, Rowv = TRUE,
                                     Colv=TRUE,
                                     distfun = dist,
                                     hclustfun = hclust,
                                     dendrogram = c("both", "row", "column", "none"),
                                     reorderfun = function(d, w) reorder(d, w),
                                     symm = FALSE,
                                     na.rm = TRUE, # NAs are removed by default!!!???
                                     rowAnnotations=NULL,
                                     columnAnnotations=NULL,
                                     colorScheme = NULL,
                                     rowSize = 13,
                                     columnSize = 13,
                                     drawGrid = TRUE,
                                     gridColor = "#808080",
                                     gridThickness = 0.1,
                                     drawValues = FALSE,
                                     width = NULL,
                                     height = NULL,...) {
  stopifnot(class(x) == "matrix")
  
  xl <- lst_object %>% filter(pert_iname == drug_name) %>% .$dataframe
  x <- xl[[1]] %>% as.matrix() # unpack list object
  
  labRow = rownames(x); labCol = colnames(x)
  
  payload <- create.payload(x,
                            labRow = labRow,
                            labCol = labCol,
                            Rowv = Rowv,
                            Colv=Colv,
                            distfun = distfun,
                            hclustfun = hclustfun,
                            dendrogram = dendrogram,
                            reorderfun = reorderfun,
                            symm = symm,
                            na.rm = na.rm,
                            rowAnnotations=rowAnnotations,
                            columnAnnotations=columnAnnotations,
                            colorScheme = colorScheme,
                            rowSize = rowSize,
                            columnSize = columnSize,
                            drawGrid = drawGrid,
                            gridColor = gridColor,
                            gridThickness = gridThickness,
                            drawValues = drawValues,
                            width = width,
                            height = height)
  
  json_object <- jsonlite::toJSON(payload$options,dataframe = "columns", null = "null", na = "null", auto_unbox = TRUE,
                                  digits = getOption("shiny.json.digits", 16), use_signif = TRUE, force = TRUE,
                                  POSIXt = "ISO8601", UTC = TRUE, rownames = FALSE, keep_vec_names = TRUE,
                                  strict_atomic = TRUE)
  return(json_object)
}


get_morpheus_json_object <- function(x, Rowv = TRUE,
                                        Colv=TRUE,
                                        distfun = dist,
                                        hclustfun = hclust,
                                        dendrogram = c("both", "row", "column", "none"),
                                        reorderfun = function(d, w) reorder(d, w),
                                        symm = FALSE,
                                        na.rm = TRUE, # NAs are removed by default!!!???
                                        rowAnnotations=NULL,
                                        columnAnnotations=NULL,
                                        colorScheme = NULL,
                                        rowSize = 13,
                                        columnSize = 13,
                                        drawGrid = TRUE,
                                        gridColor = "#808080",
                                        gridThickness = 0.1,
                                        drawValues = FALSE,
                                        width = NULL,
                                        height = NULL,...) {
 
  
  new_x <- remove_x_perc_NA(dat = x,thresh_row = 1.0, thresh_col = 0.5, d = "both")$mat
  labRow = rownames(new_x); labCol = colnames(new_x)
  # new_x <- x[!apply(is.na(x),1,all), !apply(is.na(x),2,all)]
  # new_x <- input_null_for_NA(x)
  
  
  payload <- create.payload(new_x,
                            labRow = labRow,
                            labCol = labCol,
                            Rowv = Rowv,
                            Colv=Colv,
                            distfun = distfun,
                            hclustfun = hclustfun,
                            dendrogram = dendrogram,
                            reorderfun = reorderfun,
                            symm = symm,
                            na.rm = na.rm,
                            rowAnnotations=rowAnnotations,
                            columnAnnotations=columnAnnotations,
                            colorScheme = colorScheme,
                            rowSize = rowSize,
                            columnSize = columnSize,
                            drawGrid = drawGrid,
                            gridColor = gridColor,
                            gridThickness = gridThickness,
                            drawValues = drawValues,
                            width = width,
                            height = height)
  
  json_object <- jsonlite::toJSON(payload$options,dataframe = "columns", null = "null", na = "null", auto_unbox = TRUE,
                                  digits = getOption("shiny.json.digits", 16), use_signif = TRUE, force = TRUE,
                                  POSIXt = "ISO8601", UTC = TRUE, rownames = FALSE, keep_vec_names = TRUE,
                                  strict_atomic = TRUE)
  return(json_object)
}

