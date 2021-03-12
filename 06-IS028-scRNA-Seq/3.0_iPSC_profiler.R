#' Run iPSC profiler to calc % of cells that fall w/in
#' each category of modules to broadly assign cell type.
#' 
#' Make sure you id what version of Seurat
#' you are using as v4 just came out

library(Seurat)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

#get isolation method, can be: Accumax, Ctube, or nuclei_introns
isoMethod <- "Ctube"

#initialize directories
dir.create(paste("./figures/", isoMethod, "/umap/", sep = ""), showWarnings = F)
dir.create(paste("./figures/", isoMethod, "/heatmap/", sep = ""), showWarnings = F)
dir.create(paste("./figures/", isoMethod, "/iPSC/", sep = ""), showWarnings = F)


###########################################################


# load seurat obj after you select your resolution
seur <- readRDS(paste("./data/seurat_files/",isoMethod, ".rds", sep = ""))
  
# pull marker list and associated genes
#markers <- as.list(readWorkbook("./data/iPSC_Profiler_curated_markers.xlsx", sheet = 1))
markers <- as.list(readWorkbook("./data/Nociceptor_subtype plan_2020_clean2.xlsx", sheet = 1))

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

# transpose the metadata form seurat obj w/module score
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

p <- DoHeatmap(mod_score_seurat, slot = "counts",
               features = names(markers),
               group.by = groupbyVal,
               label = colLabs,
               angle = labAngle,
               draw.lines = F,
               hjust = hjust) +
  #theme(plot.margin = grid::unit(plotMargins, "cm")) +
  labs(fill = "Module score", col = groupbyVal, title = isoMethod) +
  scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))+
  rmLegend + 
  theme(plot.title = element_text(hjust = 0.5)) +
  mytheme

p$data$Feature <- gsub("-", replacement = " ", x = p$data$Feature) 
p

#save the figure
plotName <- paste("./figures/", isoMethod, "/heatmap/", isoMethod, "_heatmap_clusters.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot = p, 
       height = 8, 
       width = 9, 
       units = "in")


######################################################
# plot umap for all clusters by selected module      #
######################################################
xLab <- "UMAP 1"
yLab <- "UMAP 2"

for(i in 1:length(names(seur_nam))){
p <- FeaturePlot(seur_plots, reduction = "umap",
                 features = seur_nam[i],
                 cols = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) +
                 xlab(xLab) + ylab(yLab) +
                 labs(title = paste(names(seur_nam)[i], ": ", isoMethod, sep = ""),
                      color = "Module score") +
                 mytheme
p

#save the figure
plotName <- paste("./figures/", isoMethod, "/umap/", isoMethod, "_", names(seur_nam[i]), "_umap_continuous.png", sep = "")
plotName <- gsub(" ", "_", plotName)
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")
}

######################################################
# plot umap for all clusters just as side by side w/ #
# above                                              #
######################################################
groupbyVal <- "seurat_clusters"
p <- DimPlot(seur, reduction = "umap",
             group.by = groupbyVal) +
             xlab(xLab) + ylab(yLab) +
             labs(title = paste("Clusters: ", isoMethod, sep = "")) +
             theme(plot.title = element_text(hjust = 0.5)) +
             mytheme

p

#save the figure
plotName <- paste("./figures/", isoMethod, "/umap/", isoMethod, "_umap_by_cluster.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")

######################################################
# plot violin plot for all clusters by selected mod- #
# ule                                                #
######################################################
hjust <- NULL
i <- 1
for(i in 1:length(names(seur_nam))){
p <- VlnPlot(seur_plots, features = seur_nam[i],
             group.by = groupbyVal) +
             theme(axis.text.x = element_text(angle = labAngle, hjust = hjust)) +
             NoLegend() +
             xlab("Cluster") +
             ylab("Module score") +
             labs(title = paste(names(seur_nam[i]), ": ", isoMethod, sep = "")) +
             mytheme
p


#save the figure
plotName <- paste("./figures/", isoMethod, "/iPSC/", isoMethod, "_", names(seur_nam[i]), "_module_score.png", sep = "")
plotName <- gsub(" ", "_", plotName)
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")
}

###########################################################

#clean up namespace
rm(list = ls())
gc()
