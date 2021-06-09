# when attempting to run tradeSeq's fitGAM method to the entire 50k gene x 40k cell dataset
# failed repeatedly after 5-7 days of runtime. Here we test various methods of down sampling
# to determine the best methods for achieving results from a representative set of features

# "tradeseq-micro" experiments used the following method:
# 	
#		slingshot method (performed on previous script <13-native-slingshot.R>)
#			generated with all cells
#			native slingshot method (not dyne)
#			no root start, end or milestone genes provided
#			results saved  in <filtered-merged-umap-slingshot.RDS>
#			
#		seurat
#			start with filtered-merged-umap-modules.RDS or presampled filtered-merged-umap-modules-10k.RDS
#			downsample_seurat_within_clusters() method to 10, 100, 1000, 10000, 20000 or all (40000) cells
#			all methods performed with genes subset to use only 2000 variable genes
#			
#		tradeSeq
#		  we alternatively used either lineage-filtered to only pseudotime/cellWeight for lineage 1, cluster 0->1 curves
#			this also requires re-filtering cells to only those with weights on lineage 1
#	
# "tradeseq-variable" experiments:
#	    all steps same as above, but uses all 5 lineages
#	      tradeSeq no longer requires re-filtering based on weight since all lineages included
#
# Conclusions
#   turns out that the only reduction required to get this to complete in <1day was 
#   filtering of genes to 2000 variable genes. All methods above completed 


# the code saved here represents the final "variable40000" method which 
# produces <tradeseq-variable40000-cells.RDS>

library(Seurat)
library(tidyverse)
library(slingshot)
library(tradeSeq)
source("ncats_09_scRNA_seq_mouse/functions.R")
DATA_DIR <- "ncats_09_scRNA_seq_mouse/data"

set.seed(5)
seur <- readRDS(file.path(DATA_DIR, "objects/filtered-merged-umap-modules.RDS"))
# seur <- downsample_seurat_within_clusters(seur, n_cells = 1000)
sling <- readRDS(file.path(DATA_DIR, "objects/filtered-merged-umap-slingshot.RDS"))


# subset slingshot objects to just the cells in out subset seurat
subset_cell_names <- rownames(seur@meta.data)

# subset cells weights for our lineage
# since we are only using one lineage, need to remove weight=0 cells
all_weights <- slingCurveWeights(sling)
cellWeights <- all_weights[subset_cell_names,]
#cellWeights <- cellWeights[cellWeights>0]

# grab pseudotime on this lineage for remaining cells
pseudotime <- slingPseudotime(sling, na = FALSE)[rownames(cellWeights),]

# subset to only variable genes
variable_counts <- seur@assays$RNA@counts[seur@assays$RNA@var.features,]
# we also need to re-subset for only remaining cells
variable_counts <- variable_counts[,rownames(cellWeights)]

# icMat <- evaluateK(counts = variable_counts,
#                    pseudotime = pseudotime,
#                    cellWeights = cellWeights,
#                    k = 3:10, 
#                    nGenes = 50, verbose = T)

dim(variable_counts)

# run the fitGAM method and save object for downstream analysis
# 6 knots was used from previous evaluateK determination
sce <- fitGAM(counts = variable_counts,
              pseudotime = pseudotime,
              cellWeights = cellWeights,
              nknots = 6,
              verbose = TRUE)
saveRDS(sce, file.path(DATA_DIR, "objects/tradeseq-variable40000-cells.RDS"))
