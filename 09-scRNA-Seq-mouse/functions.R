

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