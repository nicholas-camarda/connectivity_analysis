#' @note [DOES NOT WORK]
#' @param dat data to plot
#' @param drug_name for labeling
#' @param base_output_dir base output directory for plots
#' @param subdir is the subdir location, relative to [base_output_dir]
#' @param minr min plotting scale
#' @param maxr max plotting scale
#' @param sym impose symmetry on heatmap (for plotting correlations/connectivities)
#' @param master_key for linking annotations, 
#' @param plt_col ...
plot_morpheus_w_annotations <- function(dat, drug_name, base_output_dir, dataset, subdir = "raw", 
                                        minr = -0.75, maxr = 0.75, sym = TRUE, match_key = match_df){
  # dat <- corr_mat_reps$cell_corr_r_mat[[1]]; drug_name = corr_mat_reps$pert_iname[[1]];; minr = -0.75; maxr = 0.75; sym = TRUE;  subdir = "raw"
  # match_key <- corr_mat_reps$match_df[[1]]
 
  unique_rows1 <- match_key %>% pull(3) %>% make.unique()
  unique_cols2 <- match_key %>% pull(4) %>% make.unique()
  match_key_u <- match_key %>% mutate(unique_rows1, unique_cols2)
  
  # could rewrite this to allow for row / column annotations..
  #' @Note if you're getting small files... make sure that you have the heatmap package morpheus loaded, 
  #' not the CRAN stats package...
  
  output_dir <- file.path(base_output_dir, "morpheus_output_matrices", subdir)
  write_tsv(match_key_u, file.path(output_dir, qq("@{drug_name}-match.tsv")))
  
  rowAnnotations <- match_key_u %>% select(contains("1"))
  colAnnotations <- match_key_u %>% select(contains("2"))
  
  dir.create(output_dir, recursive = T)
  message(drug_name)
  kolors <- 3
  vals <- seq(from=minr, to=maxr,length.out = kolors)
  payload <- create.payload(x = dat, 
                            columnAnnotations = colAnnotations,
                            rowAnnotations = rowAnnotations,
                            distfun = function(x) Dist(x, method = "spearman"),
                            hclustfun = hclust,
                            symm = sym,
                            colorScheme = list(scalingMode="fixed", 
                                             stepped = FALSE, 
                                             values = vals, 
                                             colors=hcl.colors(kolors))) # using euclidian distance
  payload$options$name <- str_c(drug_name, subdir, dataset, sep= "-")
  
  widget <- htmlwidgets::createWidget(
    name = 'morpheus',
    payload,
    package = 'morpheus', 
    sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE)
  ); widget
  
  output_fn_n <- file.path(normalizePath(output_dir), str_c(drug_name, ".html"))
  htmlwidgets::saveWidget(widget,file = output_fn_n, selfcontained = TRUE)
  # message(output_fn_n)
}



#' @note save morpheus heatmaps as interactive html files
#' @param dat is a numeric matrix that will be plotted
#' @param drug_name for directory naming
#' @param base_output_dir base output directory for plots
#' @param subdir is the subdir location, relative to [base_output_dir]
#' @param minr min plotting scale
#' @param maxr max plotting scale
#' @param sym impose symmetry on heatmap (for plotting correlations/connectivities)
plot_morpheus <- function(dat, drug_name, base_output_dir, dataset, subdir = "raw", minr = -0.75, maxr = 0.75, sym = TRUE){
  # dat <- corr_mat_reps$cell_corr_r[[1]]; drug_name = corr_mat_reps$pert_iname[[1]]; minr = -0.75; maxr = 0.75; sym = TRUE;  subdir = "raw"
  
  # could rewrite this to allow for row / column annotations..
  #' @Note if you're getting small files... make sure that you have the heatmap package morpheus loaded, 
  #' not the CRAN stats package...
  
  output_dir <- file.path(base_output_dir, "morpheus_output_matrices", subdir)
  dir.create(output_dir, recursive = T)
  message(drug_name)
  kolors <- 3
  vals <- seq(from=minr, to=maxr,length.out = kolors)
  payload <- create.payload(x = dat, 
                            # distfun = function(x) Dist(x, method = "spearman"),
                            # hclustfun = hclust,
                            # symm = sym,
                            colorScheme=list(scalingMode="fixed", stepped = FALSE, values = vals, colors=hcl.colors(kolors))) # using euclidian distance
  payload$options$name <- str_c(drug_name, subdir, dataset, sep= "-")
  
  widget <- htmlwidgets::createWidget(
    name = 'morpheus',
    payload,
    package = 'morpheus',
    sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE)
  ); widget
  
  output_fn_n <- file.path(normalizePath(output_dir), str_c(drug_name, ".html"))
  htmlwidgets::saveWidget(widget,file = output_fn_n, selfcontained = TRUE)
  message(output_fn_n)
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



###### DEPRECATED
#' @note display morpheus heatmaps as interactive html files
#' @param dat is a numeric matrix that will be plotted
#' @param drug_name for naming plot
#' @param minr min plotting scale
#' @param maxr max plotting scale
#' @param sym impose symmetry on heatmap (for plotting correlations/connectivities)
display_morpheus <- function(dat, drug_name, dataset, minr = -0.75, maxr = 0.75, sym = TRUE){
  
  # could rewrite this to allow for row / column annotations..
  #' @Note if you're getting small files... make sure that you have the heatmap package morpheus loaded, 
  #' not the CRAN stats package...
  
  kolors <- 3
  vals <- seq(from=minr, to=maxr,length.out = kolors)
  payload <- create.payload(x = dat, 
                            distfun = function(x) Dist(x, method = "spearman"),
                            hclustfun = hclust,
                            symm = sym,
                            colorScheme=list(scalingMode="fixed", stepped = FALSE, values = vals, colors=hcl.colors(kolors))) # using euclidian distance
  payload$options$name <- str_c(drug_name, dataset, sep= "-")
  widget <- htmlwidgets::createWidget(
    name = 'morpheus',
    payload,
    package = 'morpheus',
    sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE)
  ); widget
  frameWidget(widget, width='90%')
}
