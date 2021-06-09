#load libraries
library(ggplot2)
library(clustifyr)
library(Seurat)
library(cowplot)
library(ggcorrplot)
library(tidyverse)

#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

#initialize directories
dir.create("./figures/all_celltypes/h_tree/", showWarnings = F)

###########################################################

#load the significant hits
df <- read.delim("./results/DEGs/all_celltypes_DEGs.tsv", sep = "\t")

#load seurat object
seur <- readRDS("./data/seurat_files/allcells.rds") #loads all cells seurat object, takes a long time 

#build cluster tree
seur <- BuildClusterTree(seur)

#set PDF
pdf(file = "./figures/all_celltypes/h_tree/allCelltypes_Hcluster.pdf",
    width = 8,
    height = 5) 

#plot the cluster tree
PlotClusterTree(seur)

#save file
dev.off()


celltype.averages <- AverageExpression(seur, group.by = "orig.ident")
celltype.averages <- as.matrix(celltype.averages$RNA)
celltype.averages[1:5, 1:5]

seur@meta.data$classified <- seur@meta.data$orig.ident
# Calculate correlation coefficients for each cluster (spearman by default)
vargenes <- df$gene

res <- clustify(
  input = seur@assays$RNA[,], # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  metadata = seur@meta.data, # meta.data table containing cell clusters
  cluster_col = "seurat_clusters", # name of column in meta.data containing cell clusters
  ref_mat = celltype.averages, # matrix of RNA-seq expression data for each cell type
  query_genes = vargenes # list of highly varible genes identified with Seurat
)


# Call cell types
res2 <- cor_to_call(
  cor_mat = res,                  # matrix correlation coefficients
  cluster_col = "seurat_clusters" # name of column in meta.data containing cell clusters
)
res2[1:5, ]


seur_meta2 <- call_to_metadata(
  res = res2,                     # data.frame of called cell type for each cluster
  metadata = seur@meta.data,           # original meta.data table containing cell clusters
  cluster_col = "seurat_clusters" # name of column in meta.data containing cell clusters
)


#set PDF
pdf(file = "./figures/all_celltypes/all_celltypes_cluster_celltype_correlations.pdf",
    width = 7,
    height = 7) 

#plot the correlation heatmap 
plot_cor_heatmap(cor_mat = res)

#save file
dev.off()


#######################################################################
#Correlation between samples

#get the average expression for each sample and convert to a matrix
cluster.exp <- AverageExpression(seur, group.by = "orig.ident")
head(cluster.exp[["RNA"]][, 1:5])
cluster.exp <- cluster.exp[["RNA"]]

#select only significant genes
cluster.exp.sig <- cluster.exp[rownames(cluster.exp) %in% df$gene,]

#generate correlation matrix from the genes listed above
corr <- round(cor(cluster.exp.sig, use = "complete.obs", method = "pearson"), 2)

#plot correlation matrix
p <- ggcorrplot(corr,
                hc.order = TRUE,
                type = "lower",
                lab = TRUE) + 
  labs(title = "Correlation Between Samples") +
  mytheme

p

#save file
filename <- paste("./figures/all_celltypes/all_celltypes_sample_correlations.png", sep = "")
ggsave(filename = filename,
       device   = "png", 
       plot     =  p, 
       height   = 6, 
       width    = 7, 
       units    = "in")
