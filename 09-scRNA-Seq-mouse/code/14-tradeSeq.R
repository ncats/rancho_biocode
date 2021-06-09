# since we need an r slingshot object to run tradeSeq, I'll rerun the TI method here
library(Seurat)
seur <- readRDS("data/objects/filtered-merged-umap.RDS")
sling <- slingshot(Embeddings(seur, "umap"), clusterLabels = seur$seurat_clusters)
saveRDS(sling, "data/objects/filtered-merged-umap-slingshot.RDS")

# install tradeSeq
# BiocManager::install("tradeSeq")
library(tradeSeq)
set.seed(5)
icMat <- evaluateK(counts = seur@assays$RNA@counts, sds = sling, k = 3:10, nGenes = 200, verbose = T)

# fit a general additive model (GAM) using a negative binomial noise distribution to model relationships between gene expression and pseudotime
set.seed(7)
pseudotime <- slingPseudotime(sling, na = FALSE)
cellWeights <- slingCurveWeights(sling)
sce <- fitGAM(counts = seur@assays$RNA@counts, pseudotime = pseudotime, cellWeights = cellWeights, nknots = 6, verbose = FALSE)

# test for significant associations between expression and pseudotime using the associationTest

# build a heatmap to visualize these differences over time