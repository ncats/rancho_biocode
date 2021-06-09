library(Seurat)
library(tidyverse)
library(glue)

# read in individual, pre-filtered samples and merge for downstream analysis
samples <- c('082120DRG','101920DRG','103020DRG','110220DRG','111220DRG')
objs <- lapply(samples, 
               function(sample)readRDS(glue("objects/{sample}-filtered.RDS")) )
seur <- merge(objs[[1]], y = objs[2:length(objs)], add.cell.ids = samples, project = "DRG")

# plot post filter merged metrics
print('violin plot of aggregated and filtered dataset')
VlnPlot(seur, features = c("nFeature_RNA","nCount_RNA","percent.mt"))
ggsave("plots/post-filter-scatterplots.png", width = 12, height = 8)


# score cells using mouse cell cycle markers
mouse_map <- read_tsv("reference/HMD_HumanPhenotype.tsv", 
                      col_names = c("human", "ensg", "entrez", 
                                    "id", "mouse", "accession", "note"))
symbol_map <- mouse_map$mouse
names(symbol_map) <- mouse_map$human

mouse_s_genes <- symbol_map[cc.genes.updated.2019$s.genes]
mouse_g2m_genes <- symbol_map[cc.genes.updated.2019$g2m.genes]

print('cell cycle scoring')
seur <- CellCycleScoring(seur, 
                       s.features = mouse_s_genes,
                       g2m.features = mouse_g2m_genes, 
                       set.ident = F)

print('normalize')
seur <- NormalizeData(seur, 
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

print('variable features')
seur <- FindVariableFeatures(seur, 
                             selection.method = "vst", 
                             nfeatures = 2000)
print('scale')
seur <- ScaleData(seur, 
                  features = rownames(seur),
                  vars.to.regress = c("S.Score", "G2M.Score"), 
                  verbose = T)
print('pca')
seur <- RunPCA(seur, 
               features = VariableFeatures(object = seur), 
               nfeatures.print = 10, ndims.print = 1:5)
print('plot elbow')
ElbowPlot(seur, ndims = 50)
ggsave("plots/elbow.png", width = 10, height = 8)

# from elbow plot we want to cluster using 28 PCs
DIM <- 28

print('umap')
seur <- RunUMAP(seur, dims = 1:DIM)
DimPlot(seur, reduction = "umap")
ggsave("plots/umap-by-sample.png", width = 10, height = 8)

# lower end resolution testing
# identify clusters by clustree
print('clustree')
library(clustree)
# copy object for optimizing resolution
res_test <- seur
res_test <- FindNeighbors(res_test, dims = 1:DIM)
res_test <- FindClusters(res_test, resolution = seq(0.01,0.1,0.01))
clustree(res_test, prefix = "RNA_snn_res.")
ggsave("plots/clustree-plot_low-end.png", width = 14, height = 7)

RES <- 0.03
seur <- FindNeighbors(seur, dims = 1:DIM)
seur <- FindClusters(seur, resolution = RES)
DimPlot(seur, reduction = 'umap', label = T)
ggsave(glue("plots/umap-{RES}.png"), width = 8, height = 8)


table(seur@meta.data$seurat_clusters)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19 
# 5878 5090 4916 2818 2717 2540 2420 2106 1893 1783 1589 1342 1336  778  594  568  535  416  101   84 

# s_test <- subset(seur, seurat_clusters %in% c(18, 19), invert=T)
# DimPlot(s_test, reduction = "umap", label = T) + theme(legend.position = "none")

saveRDS(seur, "objects/filtered-merged-umap.RDS")


