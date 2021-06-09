#load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(clustree)
library(tidyverse)

#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

#initialize directories
dir.create(paste("./figures/all_celltypes/", sep = ""), showWarnings = F)
dir.create(paste("./results/", sep = ""), showWarnings = F)
dir.create("./results/DEGs", showWarnings = F)
dir.create("./results/DEGs/all_celltypes_pairwise/", showWarnings = F)
dir.create("./results/DEGs/all_celltypes_pairwise/cluster/", showWarnings = F)
dir.create("./results/DEGs/all_celltypes_pairwise/sample/", showWarnings = F)


#load seurat objects
#df.combined <- readRDS("./data/seurat_files/allcells.rds")

###########################################################

#get a list of the samples in the data directory
data_files <- list.files("./data/")

#remove excel and seurat files and from the list
data_files <- data_files[!str_detect(data_files, pattern=".xlsx|seurat")]

for(i in data_files){
  print(i)
  #set data file location
  data.file <- paste("data/", i, "/raw_feature_bc_matrix/", sep = "")
  # Load the dataset
  df.data <- Read10X(data.dir = data.file)
  # Initialize the Seurat object with the raw (non-normalized data).
  assign(paste(i,"seur", sep = "."), CreateSeuratObject(counts = df.data, project = i, min.cells = 3, min.features = 200))
}

df.combined <- merge(x = D21.seur, 
                     y = c(D50.seur, `hPSC-FCDI-Astro.seur`, `hPSC-WA09-Astro-2.seur`, 
                           `PHA-Astro.seur`, `Glutamate-neurons.seur`), 
                     add.cell.ids = data_files)


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# here we add the percent mitochondrial for each cell
df.combined[["percent.mt"]] <- PercentageFeatureSet(df.combined, pattern = "^MT-")

#remove large unused seurat objects free up space
rm(list = ls(pattern = ".seur"), df.data)

#############################################
# QC plots
#############################################

# Visualize QC metrics as a violin plot
p <- VlnPlot(df.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#output the plot
p

#save the figure
plotName <- paste("./figures/all_celltypes/all_celltypes_QC_violin.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 12, 
       units = "in")


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# this will plot the data prior to any sort of QC filtering with the correlation at the top
plot1 <- FeatureScatter(df.combined, feature1 = "nCount_RNA", feature2 = "percent.mt") + mytheme
plot2 <- FeatureScatter(df.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + mytheme
p <- plot1 + plot2

#output the figure
p

#save the figure
plotName <- paste("./figures/all_celltypes/all_celltypes_QC_FeatureScatter.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 12, 
       units = "in")


#############################################
# Subset Data
#############################################

# next is sort of a manual step to determine potential cutoffs for each of the QC metrics 
# based on the plots. Examine them to see if there are any particular cutoffs that make sense
# apply filtering based on previous violin plots based on IQR

#get IQR
iqr <- summary(df.combined$nFeature_RNA)
minCutoff <- unname(iqr)[1]
maxCutoff <- unname(iqr)[5]

#doing 5% and 90% instead of first and third quartiles
Qdf <- quantile(df.combined$nFeature_RNA, c(0.05, 0.9))
minCutoff <- unname(Qdf)[1]
maxCutoff <- unname(Qdf)[2]

#remove cells that do not pass the feature count thresholds and mitochondrial percentage cutoff
df.combined <- subset(df.combined, subset = nFeature_RNA > minCutoff & nFeature_RNA < maxCutoff & percent.mt < 10) 

#############################################
# Normalize the data
#############################################

# Normalize the data using a global-scaling normalization method "LogNormalize" that normalizes the 
# feature expression measurements for each cell by the total expression, multiplies this by a scale factor 
# (10,000 by default), and log-transforms the result. Normalized values are stored in df[["RNA"]]@data
df.combined <- NormalizeData(df.combined)



#############################################
# Feature Selection
#############################################

# Calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in
# some cells, and lowly expressed in others). Focusing on these genes  helps to highlight biological signal in single-cell datasets.
# directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures function. 
#By default, it returns 2,000 features per dataset. These will be used in downstream analysis, like PCA.
df.combined <- FindVariableFeatures(df.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(df.combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(df.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#create and save plot
p <- plot2 + 
  labs(title = "Expression Vs. Variance in All Cell Types")
p

#save the figure
plotName <- paste("./figures/all_celltypes/all_celltypes_VariableFeaturePlot.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 9, 
       units = "in")



#############################################
# Scaling the Data
#############################################

# this applies a linear transformation ('scaling') that is a standard pre-processing 
# step prior to dimensional reduction techniques like PCA. This will shift the expression 
# of each gene, so that the mean expression across cells is 0 and scales the expression of 
# each gene, so that the variance across cells is 1. This step gives equal weight in 
# downstream analyses, so that highly-expressed genes do not dominate. The results of 
# this are stored in df.combined[["RNA"]]@scale.data

#get gene list from data
all.genes <- rownames(df.combined)

#scale the data for all genes
df.combined <- ScaleData(df.combined, features = all.genes)



#############################################
#Perform linear dimensional reduction
#############################################

# Next we perform PCA on the scaled data. By default, only the previously determined 
#variable features are used as input, but can be defined using features argument 
#if you wish to choose a different subset.

#run the PCA and add it to the Seurat Object
df.combined <- RunPCA(df.combined, features = VariableFeatures(object = df.combined))


# Examine and visualize PCA results a few different ways
# this will print the top 5 positive and negative genes
# associated with each PC
print(df.combined[["pca"]], dims = 1:5, nfeatures = 5)

# plot the top loadings for genes associated with PC1 and PC2
p <- VizDimLoadings(df.combined, dims = 1:4, reduction = "pca")
p

#plot the PCA using PC1 and PC2
p <- DimPlot(df.combined, reduction = "pca") + 
  labs(x = "PC1",
       y = "PC2",
       title = "PCA for all cell types") +
  mytheme

p

#save the figure
plotName <- paste("./figures/all_celltypes/all_celltypes_PCA.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 9, 
       units = "in")



#############################################
# Determine the 'dimensionality' of the dataset
#############################################

#examine the dimensionality of the data using an elbow plot
#You are looking to find the dim where the plot starts to level off
#10 dims with nuclei_introns
ElbowPlot(df.combined, ndims = 30)


#set dims based on the elbow plot above. For this dataset, it seems
#like 10 is a good cutoff for each isolation method
dims <- 10 

# will need to adjust dims based on above
df.combined <- FindNeighbors(df.combined, dims = 1:dims) 

# this iterates through different resolutions & runs
# clustree to help with identifying optimal resolution
p <- clustree(FindClusters(df.combined, resolution = seq(0.01,0.1,0.01))) + 
  labs(title = "All Cell Types") +
  theme(plot.title = element_text(hjust = 0.5))

#output clustree diagram
p 

#save the figure
plotName <- paste("./figures/all_celltypes/all_celltypes_clustree.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 10, 
       width = 12, 
       units = "in")


#will need to select the appropriate resolution based on clustree
#diagram. You are looking for the point where the columns do not 
#change much
#Combined resolution = 0.05
res = 0.04
df.combined <- FindClusters(df.combined, resolution = res)

#create the data for UMAP
df.combined <- RunUMAP(df.combined, dims = 1:dims) 

#visualize clusters using UMAP
p <- DimPlot(df.combined, reduction = "umap", group.by = "seurat_clusters") + 
        labs(title = "UMAP by Cluster") +
        mytheme

p

plotName <- paste("./figures/all_celltypes/all_combined_umap_cluster", res, ".png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")


#visualize cell types using UMAP
p <- DimPlot(df.combined, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP by Cell Type") +
  mytheme

p

plotName <- paste("./figures/all_celltypes/all_combined_umap_celltype", res, ".png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")

# Look at cluster IDs of the first 5 cells
head(Idents(df.combined), 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
# this take a long time to run
df.combined.markers <- FindAllMarkers(df.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
df.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#save all significant hits
write.table(df.combined.markers, file = paste("./results/DEGs/all_celltypes_DEGs.tsv", sep = ""), sep = "\t", row.names = T)

### create a heatmap ###
#Explore the top 50 genes for each cluster and create a heatmap
top100 <- df.combined.markers %>% top_n(n = 50, wt = avg_log2FC)

#create heatmap
p <- DoHeatmap(df.combined, features = top100$gene, group.by = "orig.ident") + 
  #NoLegend() +
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                       mid = "white", 
                       high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                       midpoint = 0, guide = "colourbar", aesthetics = "fill", na.value = "white") +
  mytheme

p

plotName <- paste("./figures/all_celltypes/all_combined_top50_heatmap.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 10, 
       width = 10, 
       units = "in")


####################################################
### dge with all pairwise comparisons by cluster ###
####################################################
# identify the unique clusters
clusters <- unique(Idents(df.combined))
# create combinations of those clusters
pairwise <- combn(clusters, 2)

# run dge of cluster x vs y & save each comp
# to a csv file named with the comparison
object_clust_x_vs_y <- lapply(1:ncol(pairwise), function(x) {
  y <- FindMarkers(df.combined, 
                   ident.1 = pairwise[1, x], 
                   ident.2 = pairwise[2, x],
                   logfc.threshold = 1,
                   test.use = "wilcox", 
                   min.pct = 0.1)
  nam <- paste("cluster", pairwise[, x], collapse = "_vs_", sep = "")
  write.csv(y, file = paste0("./results/DEGs/all_celltypes_pairwise/cluster/all_celltypes_", nam, "_dge.csv"), row.names = T)
  return(y)
})


###################################################
### dge with all pairwise comparisons by cell line ###
###################################################

### dge with all pairwise comparisons ###
# identify the unique clusters
clusters <- unique(df.combined@meta.data$orig.ident)
# create combinations of those clusters
pairwise <- combn(clusters, 2)

#set the cell line name as the identity
Idents(df.combined.samples) <- "orig.ident"

# run dge of cluster x vs y & save each comp
# to a csv file named with the comparison
object_clust_x_vs_y <- lapply(1:ncol(pairwise), function(x) {
  y <- FindMarkers(df.combined.samples, 
                   ident.1 = pairwise[1, x], 
                   ident.2 = pairwise[2, x],
                   logfc.threshold = 1,
                   test.use = "wilcox", 
                   min.pct = 0.1)
  nam <- paste(pairwise[, x], collapse = "_vs_", sep = "")
  write.csv(y, file = paste0("./results/DEGs/all_celltypes_pairwise/sample/all_celltypes_", nam, "_dge.csv"), row.names = T)
  return(y)
})

#creating a new seurat object to have the identities set
#as the cell line instead of the clusters
#df.combined.samples <- df.combined
#Idents(object = df.combined.samples) <- df.combined.samples@meta.data$sample


#save seurat object
saveRDS(df.combined, file = paste("./data/seurat_files/allcells.rds", sep = ""))
###########################################################

#clean up workspace
rm(list = ls())
gc()
