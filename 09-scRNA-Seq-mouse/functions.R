library(glue)

plot_cutoffs <- function(obj, 
                         metrics_list = c("nFeature_RNA","nCount_RNA","percent.mt"),
                         filterset){
  
  # extract feature data from seurat to plot
  df <- FetchData(obj, vars = metrics_list) %>% 
    mutate(sample = obj@project.name)
  
  # course filter to highlight lower region
  if("percent.mt" %in% metrics_list) df <- filter(df, percent.mt < 50)
  
  # grow a list of cutoffs
  cutlines <- tibble(metric = names(unlist(filterset)),
                     value = unlist(filterset)) %>% 
    mutate(metric = gsub("\\.[A-za-z]+$","",metric))
  
  df  %>%
    pivot_longer(cols = -sample,
                 names_to = "metric",
                 values_to = "value") %>%
    ggplot(aes(x=sample,y=value)) +
    geom_jitter(size=0.5) +
    geom_violin(color="red", fill="pink",alpha=0.5) +
    geom_hline(data=cutlines, aes(yintercept = value)) +
    facet_wrap(~metric, scales = "free") +
    theme_minimal() +
    theme(axis.title = element_blank())
}


# filterset is a nested list of meta.data column names with either min or max values
# filterset$nFeature_RNA$min=1200
# filterset$percent.mt$max=10
apply_filterset <- function(obj, filterset){
  filtered_obj <- obj
  
  # sequentially apply each filter using min/max limits
  for(f in names(filterset)) {
    if(f %in% names(filtered_obj@meta.data)){
      # set defaults for missing limits
      min_val <- filterset[[f]]$min %||% 0
      max_val <- filterset[[f]]$max %||% Inf
      
      expr <- FetchData(filtered_obj, vars = f)
      filtered_obj <- filtered_obj[, which(x = expr > min_val & expr < max_val)]
    }else{
      warning(glue::glue('Attempting to filter using feature name not in this object.
                       \tfilter name: {f}; seurat object: {filtered_obj@project.name}'))
      return()
    }
  }
  # report number of cells removed from this object
  message(glue::glue("{dim(obj)[[2]] - dim(filtered_obj)[[2]]} cells filtered from {filtered_obj@project.name}"))
  return(filtered_obj)
}

downsample_seurat_within_clusters <- function(obj, 
                                              n_cells=NULL, 
                                              fraction_cells=NULL, 
                                              cluster_col = 'seurat_clusters'){
  
  if( cluster_col %in% names(obj@meta.data) ){
    
    cluster_names <- as.character(unique(obj@meta.data$seurat_clusters))
    
    # convert n_cells to percentage so we can apply it within each cluster
    # if both are specified this will override fraction_cells
    if(!is.null(n_cells)){
      fraction_cells <- n_cells / dim(obj)[[2]]
    }
    
    # grab subsample of cells from within each cluster if fraction_cells is 0-1
    if( fraction_cells > 0 & fraction_cells <= 1 ){
      cell_subset <- sapply(cluster_names, function(x){
        cluster_cell_names <- rownames(obj@meta.data)[obj@meta.data[cluster_col]==x]
        sample(cluster_cell_names, size = (length(cluster_cell_names) * fraction_cells))
      })
      cell_subset <- unlist(cell_subset)
    }else{
      warning('n_cells or fraction_cells not appropriate for this seurat object')
      return()
    }
    message(glue('Subset Seurat object "{seur@project.name}" to {length(cell_subset)} total cells'))
    subset(seur, cells = cell_subset)
    
  }else{
    warning('cluster_col not found in this seurat object')
    return()
  }
}


### Code for extracting slingshot curves for ggplot2 ----------------------
### Point on curve function ----
points_on_curve <- function(curve, lambda, ...) {
  UseMethod("points_on_curve", curve)
}

points_on_curve.principal_curve <- function(curve, lambda, ...) {
  if (nrow(curve$s) == length(curve$lambda)) { # didn't use approx_points
    S <- apply(curve$s, 2, function(sjj) {
      return(approx(
        x = curve$lambda[curve$ord],
        y = sjj[curve$ord],
        xout = lambda, ties = "ordered"
      )$y)
    })
  } else {
    if (all(curve$ord == seq_along(curve$lambda))) { # used approx_points
      curvelambda <- seq(min(curve$lambda), max(curve$lambda), length.out = nrow(curve$s))
      S <- apply(curve$s, 2, function(sjj) {
        return(approx(
          x = curvelambda,
          y = sjj,
          xout = lambda, ties = "ordered"
        )$y)
      })
    }
  }
  return(S)
}

points_on_curve.SlingshotDataSet <- function(curve, lambda, ...) {
  locs <- lapply(slingCurves(curve), function(crv) {
    points_on_curve(crv, lambda, ...)
  })
  locs <- do.call('rbind', locs)
  colnames(locs) <- paste0("dim", seq_len(ncol(locs)))
  return(as.data.frame(locs))
}

### Extend ggplot function

#' Plot the gene in reduced dimension space
#'
#' @param sds The output from a lineage computation
#' @param col The assignation of each cell to a label. If none is provided, 
#' default to the cluster labels provided when creating the \code{\link{SlingshotDataSet}}
#' @param title Title for the plot.
#' @param lineSize Size of the curve lineages. Default to 1.
#' @param ... Other options passed to \code{\link{geom_point}}
#' @return A \code{\link{ggplot}} object
#' @examples
#' data('slingshotExample')
#' sds <- slingshot(rd, cl)
#' gg_plot(sds)
#' 
#' ## Change point size and transparency
#' gg_plot(sds, size = 2, alpha = .5)
#' 
#' ## Use grey background
#' 
#' gg_plot(sds) + theme_grey()
#' 
#' ## Color by gene expression
#' gene_count <- sample(0:10, nrow(reducedDims(sds)), replace = TRUE) 
#' gg_plot(sds, col = gene_count)
#' 
#' ## Add a marker of pseudotime
#' gg_plot(sds) + geom_point(data = points_on_curve(sds, 10), size = 3)
#' @importFrom slingshot slingPseudotime slingCurves reducedDim slingClusterLabels
#' @import ggplot2
#' @export

gg_plot <- function(sds, col = NULL, title = NULL, lineSize = 1, ...) {
  rd <- reducedDim(sds)
  
  if (is.null(col)) {
    cl <- slingClusterLabels(sds)
    if ("matrix" %in% is(cl)) {
      cl <- apply(cl, 1, which.max)
      cl <- as.character(cl)
    }
  } else {
    cl <- col
  }
  
  # Getting the main plot
  df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cl)
  p <- ggplot(df, aes(x = dim1, y = dim2, col = col)) +
    geom_point(...) +
    theme_classic() +
    labs(title = title, col = "")
  
  # Adding the curves
  for (i in seq_along(slingCurves(sds))) {
    curve_i <- slingCurves(sds)[[i]]
    curve_i <- curve_i$s[curve_i$ord, ]
    colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1)
  }
  return(p)
}

