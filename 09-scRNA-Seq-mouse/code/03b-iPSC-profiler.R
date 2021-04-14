#' Run iPSC profiler to calc % of cells that fall w/in
#' each category of modules to broadly assign cell type.
#' 
#' Make sure you id what version of Seurat
#' you are using as v4 just came out

library(Seurat)
library(RColorBrewer)
library(tidyverse)

#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))


###########################################################


# load seurat obj after you select your resolution
seur <- readRDS("objects/filtered-merged-umap.RDS")

# pull marker list and associated genes
ipsc_db <- markers <- read_tsv("iPSC-modules/iPSC-mouse-markers.tsv") %>% 
  # filter for valid symbols
  filter(!is.na(mouse)) %>% 
  filter(mouse %in% rownames(s@assays$RNA) ) %>% 
  
  # clean up module names
  mutate(name = make.names(name)) %>% 
  # prefix ipsc module names
  mutate(name = str_to_upper( glue("ipsc_{name}")) ) %>% 
  transmute(gs_name = name, gene_symbol = mouse)

markers <- read_tsv("iPSC-modules/iPSC-mouse-markers.tsv", 
                               col_types = cols( name = col_character(),
                                                 human = col_character(),
                                                 mouse = col_character()) ) %>% 
  
  # filter for valid symbols
  filter(!is.na(mouse)) %>% 
  filter(mouse %in% rownames(s@assays$RNA) ) %>% 
  
  # clean up module names
  mutate(name = make.names(name)) %>% 
  mutate(name = gsub("_","-",name)) %>% 
  
  # group into named modules with gene lists
  group_by(name) %>% 
  summarise(genes = list(mouse))

s_processed <- AddModuleScore(seur, features = markers$genes)

# rename modules in meta.data table
moduleNamesRefVec <- markers$name
names(moduleNamesRefVec) <- paste0(rep("Cluster",length(moduleNamesRefVec)),seq(1,length(moduleNamesRefVec)))
names(s_processed@meta.data)[(names(s_processed@meta.data) %in% names(moduleNamesRefVec))] <-  moduleNamesRefVec[names(s_processed@meta.data)[(names(s_processed@meta.data) %in% names(moduleNamesRefVec))]]
moduleNamesRefVec


# run module score with marker list
# seuratData_withModuleScore <- AddModuleScore(seur, features = markers)
# seur_plots <- s_processed
# module_names <- names(s_processed@meta.data)[7:44]


# append "Cluster" to the modules (confusing, it's acutally module)
# moduleNamesRefVec <- paste0("Cluster", 1:length(markers))
# 
# # set names as module names
# names(moduleNamesRefVec) <- names(markers)
# names(moduleNamesRefVec) <- gsub("-", replacement = " ", x = names(moduleNamesRefVec))
# seur_nam <- moduleNamesRefVec

# transpose the metadata from seurat obj w/module score
moduleScoreMat <- t(s_processed@meta.data[, moduleNamesRefVec])

# #add rownames based on module names
# rownames(moduleScoreMat) <- names(moduleNamesRefVec)

# # create a new seurat object with the matrix
seuratData_ModuleScore <- CreateSeuratObject(moduleScoreMat)
# 
# # label column names based on meta.data in original 
# # seurat object
for(i in colnames(seur@meta.data)) {
  seuratData_ModuleScore@meta.data[, i] <- seur@meta.data[, i]
}
# 
# # store in a diff obj
# mod_score_seurat <- seuratData_ModuleScore

######################################################
# plot heatmap for all modules - if you want to look #
# at specific modules, you'll have to adjust code    #
######################################################


mod_score_seurat <- seuratData_ModuleScore[, sample(colnames(seuratData_ModuleScore), size = 10000, replace=F)]


mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

groupbyVal <- "seurat_clusters"
labAngle <- 45
hjust <- 0
rmLegend <- NULL
plotMargins <- c(1,1,1,1,1)
colLabs <- F

library(RColorBrewer)
p <- DoHeatmap(mod_score_seurat, 
               slot = "counts",
               features = moduleNamesRefVec,
               group.by = groupbyVal,
               label = colLabs,
               angle = labAngle,
               draw.lines = F,
               hjust = hjust) +
  # scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))+
  rmLegend +
  theme(plot.title = element_text(hjust = 0.5)) 

p$data$Feature <- gsub("-", replacement = " ", x = p$data$Feature) 
p




#save the figure
ggsave(filename = "plots/ipsc-module-heatmap.png", 
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
