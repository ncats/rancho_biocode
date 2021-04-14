
## Method for pK Identification without ground-truth calls
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

# devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(tidyverse)
library(glue)
library(future)
options(future.globals.maxSize = +Inf)
plan(strategy = "multicore")

# Performs pN-pK parameter sweeps on a 10,000-cell subset of a
# pre-processed Seurat object. Will use all cells if Seurat object
# contains less than 10,000 cells. Results are fed into
# 'summarizeSweep' and 'find.pK' functions during optimal pK parameter
# selection workflow. Parameters tested: pN = 0.05-0.3, pK =
# 0.0005-0.3
s <- readRDS("110220DRG-processed.RDS")
sweep_res_list <- paramSweep_v3(s, PCs = 1:10, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
bcmvn <- find.pK(sweep_stats)

# view plot/table and find pK = 0.01 maximum
best_pk <- 0.005

library(clustree)
s <- FindNeighbors(s, dims = 1:10)
s <- FindClusters(s, resolution = 0.4)
DimPlot(s, reduction = "umap")

head(s@meta.data)
homotypic.prop <- modelHomotypic(s@meta.data$RNA_snn_res.0.4) ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(s@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# with homotypic population estimates
xx <- doubletFinder_v3(s, PCs = 1:10, pN = 0.25, pK = best_pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(xx@meta.data)
table(xx@meta.data$DF.classifications_0.25_0.005_827)
dbl <- FetchData(xx,vars = "DF.classifications_0.25_0.005_827")
DimPlot(xx, reduction = "umap", group.by  = "DF.classifications_0.25_0.005_827")

# without homotypic calls
hh <- doubletFinder_v3(s, PCs = 1:10, pN = 0.25, pK = best_pk, nExp = 0, reuse.pANN = FALSE, sct = FALSE)
head(hh@meta.data)
table(hh@meta.data$DF.classifications_0.25_0.005_0)
dbl <- FetchData(hh,vars = "DF.classifications_0.25_0.005_0")
DimPlot(hh, reduction = "umap", group.by  = "DF.classifications_0.25_0.005_0")
