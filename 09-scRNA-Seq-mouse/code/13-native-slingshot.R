# Since the tradeSeq method downstream requires export from a native slingshot 
# object we are running again here after an initial exploration with dyne-wrapped version

library(Seurat)
library(tidyverse)
library(slingshot)
source("ncats_09_scRNA_seq_mouse/functions.R")
DATA_DIR <- "ncats_09_scRNA_seq_mouse/data"

seur <- readRDS(file.path(DATA_DIR, "objects/filtered-merged-umap-modules.RDS"))
sling <- slingshot(Embeddings(seur, "umap"), clusterLabels = seur$seurat_clusters)
saveRDS(sling, file.path(DATA_DIR, "objects/filtered-merged-umap-slingshot.RDS"))
