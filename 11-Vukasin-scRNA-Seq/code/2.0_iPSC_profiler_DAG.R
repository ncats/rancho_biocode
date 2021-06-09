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

#set isolation method
celltype <- "all_celltypes"

#initialize directories
dir.create(paste("./figures/", celltype, "/heatmap/", sep = ""), showWarnings = F)
dir.create(paste("./figures/", celltype, "/iPSC/", sep = ""), showWarnings = F)
dir.create(paste("./figures/", celltype, "/umap/", sep = ""), showWarnings = F)
dir.create(paste("./figures/", celltype, "/signatures/", sep = ""), showWarnings = F)


###########################################################

# load seurat obj after you select your resolution
#seur <- readRDS(paste("./data/seurat_files/",celltype, ".rds", sep = ""))
seur <- readRDS("./data/seurat_files/allcells.rds") #loads all cells seurat object, takes a long time 
  
# pull marker list and associated genes
#markers <- as.list(readWorkbook("./data/iPSC_Profiler_curated_markers.xlsx", sheet = 1))
markers <- as.list(readWorkbook("./data/astrocyte_profiler.xlsx", sheet = 1))

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
  labs(fill = "Module score", col = groupbyVal) +
  scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))+
  rmLegend + 
  theme(plot.title = element_text(hjust = 0.5)) +
  mytheme

p$data$Feature <- gsub("-", replacement = " ", x = p$data$Feature) 
p

#save the figure
plotName <- paste("./figures/", celltype, "/heatmap/", celltype, "_heatmap_clusters.png", sep = "")
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
                 features = seur_nam[i]) +
                 #cols = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) +
                 xlab(xLab) + ylab(yLab) +
                 labs(title = paste(names(seur_nam)[i], sep = ""),
                      color = "Module score") +
                 mytheme
p

#save the figure
plotName <- paste("./figures/", celltype, "/umap/", celltype, "_", names(seur_nam[i]), "_umap_continuous.png", sep = "")
plotName <- gsub(" ", "_", plotName)
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")
}


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
             labs(title = paste(names(seur_nam[i]), ": ", celltype, sep = "")) +
             mytheme
p


#save the figure
plotName <- paste("./figures/", celltype, "/iPSC/", celltype, "_", names(seur_nam[i]), "_module_score.png", sep = "")
plotName <- gsub(" ", "_", plotName)
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")
}


###########################################################

#to create a feature plot based on a specific gene
p <- FeaturePlot(seur, reduction = "umap", features = "AGT") +
  xlab(xLab) + ylab(yLab) +
  #labs(title = paste("Clusters: ", celltype, sep = "")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  mytheme

p

###########################################################
# examining heatmap with a few gene signatures
###########################################################

#high-fidelity astrocyte genes (Kelley et al 2018)
markers_Astrocyte <- c("AGT", "EDNRB", "GJA1", "PON2", "BBOX1", "SLC1A2", "TP53BP2", "ALDH6A1", "SPON1", "NTRK2", "MERTK", "RAB31", "GPAM", "NOTCH2", "ARHGEF26", "GRAMD1C", "RAB34", "PPAP2B", "MLC1", "GPR125", "BMPR1B", "SLC1A3", "PDLIM5", "PLSCR4", "HEPH", "AQP4", "ATP1A2", "FGFR3", "ALDH2", "SLC4A4", "SLC25A18", "SOX2", "TRIL", "ETNPPL", "AHCYL1", "METTL7A")
#Glutamate neurons should show high-expression of high-fidelity neuron genes: 
markers_Glutamate_Neurons <- c("NSF", "CNR1", "TRIM37", "SNAP25", "MAPK9", "PPP3CB", "ATP2A2")
#Additional useful markers as positive controls might be considered from Polioudakis et al 2019 
#oRG
markers_Outer_Radial_Glia <- c("PTN", "VIM", "PEA15", "SFRP1", "PTPRZ1", "FABP7", "SOX2", "SLC1A3", "HES1", "FOS", "ID4", "JUN", "DBI", "FABP5", "NRG1", "C1orf61")
#OPC
markers_Oligodendrocyte_Precursors <- c("PTPRZ1", "BCAN", "SCRG1", "OLIG1", "PMP2", "PDGFRA", "DBI", "EGR1", "LIMA1", "PTN", "APOD", "EPN2", "S100B")
#vRG 
markers_Ventricular_Radial_Glia <- c("EGR1", "HES1", "CYR61", "SFRP1", "ID4", "SOX9", "ZFP36L1", "IER2")
#intermediate progenitor 
markers_Intermediate_Progenitor <- c("PPP1R17", "SSTR2", "HES6", "SOX4", "CCND2", "EZR", "EOMES", "ENC1", "PRDX1", "CORO1C", "LINC01158", "SYNE2", "HNRNPA1", "TMSB4X", "NHLH1", "MARCKS", "PENK")
#maturing excitatory neurons 
markers_Maturing_Excitatory_Neurons <- c("SATB2", "STMN2", "LIMCH1", "SYT4", "GAP43", "SLA", "NEUROD2", "ZBTB18", "PLXNA4", "GPM6A", "ARPP21", "ANK2", "NEUROD6", "VCAN", "MEF2C")
#maturing upper-layer excitatory neurons 
markers_Maturing_Upper_Layer_Excitatory_Neurons <- c("NEFM", "NEFL", "RUNX1T1", "CELF2", "MYT1L", "LMO4", "NELL2", "FABP7", "R3HDM1", "MAP1B", "GAP43", "DAAM1")
#migrating neurons 
markers_Migrating_Neurons <- c("ENC1", "SOX4", "MEIS2", "EZR", "ROBO2", "SOX11", "GRIA2", "POU3F2", "MLLT3", "NFIB", "DDAH2", "SDCBP", "EPHB6", "ZNF704")

#make the heatmap for each signature above
#you can also do this with the cell lines instead of clusters
#replace groupbyval with orig.ident
for(i in ls(pattern = "markers_")){
  nam <- gsub("markers_", "", i)
  groupbyVal <- "seurat_clusters"
  labAngle <- 45
  hjust <- 0
  rmLegend <- NULL
  plotMargins <- c(1,1,1,1,1)
  colLabs <- F
  
  DoHeatmap(object = seur,
            features = mget(i)[[1]],
            size = 3, disp.max = 1.4, disp.min = -1.5) + 
    scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                          mid = "white", 
                          high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill") 
  
  
  p <- DoHeatmap(seur, 
                 #slot = "counts",
                 features = mget(i)[[1]],
                 group.by = groupbyVal,
                 label = colLabs,
                 angle = labAngle,
                 draw.lines = F,
                 hjust = hjust) +
    #theme(plot.margin = grid::unit(plotMargins, "cm")) +
    labs(fill = "Expression", 
         col = "Cluster",
         title = gsub("_", " ", nam)) +
    scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                          mid = "white", 
                          high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")  +
    rmLegend + 
    theme(plot.title = element_text(hjust = 0.5)) +
    mytheme
  
  p$data$Feature <- gsub("-", replacement = " ", x = p$data$Feature) 
  p
  
  #save the figure
  plotName <- paste("./figures/", celltype, "/signatures/", nam, "_heatmap_cluster.png", sep = "")
  ggsave(filename = plotName, 
         device = "png", 
         plot = p, 
         height = 8, 
         width = 9, 
         units = "in")

}

#clean up namespace
rm(list = ls())
gc()
