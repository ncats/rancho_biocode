#' Run iPSC profiler to calc % of cells that fall w/in
#' each category of modules to broadly assign cell type.
#' 
#' Make sure you id what version of Seurat
#' you are using as v4 just came out

library(Seurat)
library(tidyverse)

s <- readRDS("objects/filtered-merged-umap.RDS")

# 1248 > 1187 > 1185
markers <- read_tsv("reference/iPSC-mouse-markers.tsv") %>% 
  # filter for valid symbols
  filter(!is.na(mouse)) %>% 
  filter(mouse %in% rownames(s@assays$RNA) ) %>% 
  
  # clean up module names
  mutate(name = make.names(name)) %>% 
  
  # group into named modules with gene lists
  group_by(name) %>% 
  summarise(genes = list(mouse))

new <- s
new <- AddModuleScore(new, features = markers['genes'][[1]], name = markers['name'][[1]])  

apply(markers[1:2,],1,function(x){
  print(x['genes'])
  new <- AddModuleScore(new, features = x['genes'], name = x['name'])  
})

names(new@meta.data)
FeaturePlot(new, features = names(new@meta.data)[7:15])
ggsave("data/plots/ipsc-feature-plots-1.png", width = 10, height = 8)

FeaturePlot(new, features = names(new@meta.data)[16:24])
ggsave("data/plots/ipsc-feature-plots-2.png", width = 10, height = 8)

FeaturePlot(new, features = names(new@meta.data)[25:33])
ggsave("data/plots/ipsc-feature-plots-3.png", width = 10, height = 8)

FeaturePlot(new, features = names(new@meta.data)[34:42])
ggsave("data/plots/ipsc-feature-plots-4.png", width = 10, height = 8)

FeaturePlot(new, features = names(new@meta.data)[43:44])
ggsave("data/plots/ipsc-feature-plots-5.png", width = 8, height = 4)







# pull marker list and associated genes
markers <- as.list(openxlsx::readWorkbook("reference/iPSC_Profiler_curated_markers.xlsx", sheet = 1))
markers <- lapply(markers, function(x) {
  # filter out na's from each profiler list
  y <- x[!is.na(x)]
  # keep only genes found in the seurat obj
  y <- y[y %in% rownames(seur)]
})
# remove marker list elements that are empty
markers <- markers[sapply(markers, length) > 0]
names(markers) <- gsub("_", replacement = "-", x = names(markers))

# run module score with marker list
seuratData_withModuleScore <- AddModuleScore(seur, features = markers)
seur_plots <- seuratData_withModuleScore
# append "Cluster" to the modules (confusing, it's acutally module)
moduleNamesRefVec <- paste0("Cluster", 1:length(markers))
# set names as module names
names(moduleNamesRefVec) <- names(markers)
names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
seur_nam <- moduleNamesRefVec
# transpose the metadata form seurat obj w/module
# score
moduleScoreMat <- t(seuratData_withModuleScore@meta.data[, paste0("Cluster", 1:length(markers))])
#add rownames based on module names
rownames(moduleScoreMat) <- names(markers)
# create a new seurat object with the matrix
seuratData_ModuleScore <- CreateSeuratObject(moduleScoreMat)
# label column names based on meta.data in original 
# seurat object
for(i in colnames(seur@meta.data)) {
  seuratData_ModuleScore@meta.data[, i] <- seur@meta.data[, i]
}
# store in a diff obj
mod_score_seurat <- seuratData_ModuleScore

######################################################
# plot heatmap for all modules - if you want to look #
# at specific modules, you'll have to adjust code    #
######################################################
groupbyVal <- "seurat_clusters"
labAngle <- 45
hjust <- 0
rmLegend <- NULL
plotMargins <- c(1,1,1,1,1)
colLabs <- F

g <- DoHeatmap(mod_score_seurat, slot = "counts",
               features = names(markers),
               group.by = groupbyVal,
               label = colLabs,
               angle = labAngle,
               draw.lines = F,
               hjust = hjust) +
  theme(plot.margin = grid::unit(plotMargins, "cm")) +
  labs(fill = "Module score", col = groupbyVal) +
  scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))+
  rmLegend
g$data$Feature <- gsub("-", replacement = " ", x = g$data$Feature)
g

######################################################
# plot umap for all clusters by selected module      #
######################################################
xLab <- "UMAP 1"
yLab <- "UMAP 2"

FeaturePlot(seur_plots, reduction = "umap",
            features = seur_nam[1],
            cols = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) +
  xlab(xLab) + ylab(yLab) +
  labs(title = names(seur_nam)[1],
       color = "Module score")

######################################################
# plot umap for all clusters just as side by side w/ #
# above                                              #
######################################################
groupbyVal <- "seurat_clusters"
DimPlot(seur, reduction = "umap",
        group.by = groupbyVal) +
  xlab(xLab) + ylab(yLab) +
  labs(title = "Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

######################################################
# plot violin plot for all clusters by selected mod- #
# ule                                                #
######################################################
hjust <- NULL
VlnPlot(seur_plots, features = seur_nam[1],
        group.by = groupbyVal) +
  theme(axis.text.x = element_text(angle = labAngle, hjust = hjust)) +
  NoLegend() +
  xlab("Cluster") +
  ylab("Module score") +
  labs(title = names(seur_nam[1]))

rm(list = ls())
gc()
