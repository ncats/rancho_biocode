library(tidyverse)
library(dyno)
library(dynwrap)
seur <- readRDS("objects/filtered-merged-umap.RDS")

# dyno tools don't use seurat objects, so we extract just what is needed
# only retain variable genes
variable_counts_matrix <- as.matrix(seur@assays$RNA@counts[seur@assays$RNA@var.features,])
variable_raw_matrix <- as.matrix(seur@assays$RNA@data[seur@assays$RNA@var.features,])

dataset <- wrap_expression(
  counts = t(variable_counts_matrix),
  expression= t(variable_raw_matrix),
  cell_info = data.frame(cell_id = rownames(seur@meta.data),
                         cell_info = seur@meta.data$seurat_clusters) )

# config <- create_docker_config()
# set_default_config(config)
# model.slingshot <- infer_trajectory(dataset = task, 
#                                     method = 'slingshot', 
#                                     backend = "docker")

# saveRDS(model.slingshot, "objects/model.slingshot.RDS")
model <- readRDS('ncats_09_scRNA_seq_mouse/data/objects/model.slingshot.RDS')
model <- model %>% 
  add_dimred(dyndimred::dimred_mds, 
             expression_source = dataset$expression)
saveRDS(model, "objects/model-dyno-slingshot-mds.RDS")

# library(hexbin)
# plot_dimred(model, dimred = "mds")
# 
# plot_dimred( model, 
#              expression_source = dataset$expression,
#              grouping = dataset$cell_info)



# plot_dimred(
#   model, 
#   expression_source = task$expression, 
#   grouping = grouping,
#   groups=grouping,
#   color_cells="auto"
# )