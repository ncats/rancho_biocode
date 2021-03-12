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

#get isolation method, can be: Accumax, Ctube, or nuclei_introns    
isoMethod <- "Ctube"
data.file <- paste("data/Lonza-Noci-D32-", isoMethod, "/outs/filtered_feature_bc_matrix/", sep = "")

#initialize directories
dir.create(paste("./figures/", isoMethod, sep = ""), showWarnings = F)
dir.create(paste("./figures/", isoMethod, "/QC/", sep = ""), showWarnings = F)
dir.create(paste("./results/", isoMethod, sep = ""), showWarnings = F)
dir.create("./data/DEGs", showWarnings = F)

#load seurat objects
#load("./data/seurat_files/Lonza-Noci-D32-Accumax.seurat.Rdata")

#df <- `Lonza-Noci-D32-Accumax`
###########################################################

# Load the dataset
df.data <- Read10X(data.dir = data.file)

# Initialize the Seurat object with the raw (non-normalized data).
df <- CreateSeuratObject(counts = df.data, project = isoMethod, min.cells = 3, min.features = 200)

# get a summary of the Seurat object
df

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# here we add the percent mitochondrial for each cell
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")



#############################################
# QC plots
#############################################

# Visualize QC metrics as a violin plot
p <- VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#output the plot
p

#save the figure
plotName <- paste("./figures/", isoMethod, "/QC/", isoMethod, "_QC_violin.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# this will plot the data prior to any sort of QC filtering with the correlation at the top
plot1 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt") + mytheme
plot2 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + mytheme
p <- plot1 + plot2

#output the figure
p

#save the figure
plotName <- paste("./figures/", isoMethod, "/QC/", isoMethod, "_QC_FeatureScatter.png", sep = "")
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
iqr <- summary(df$nFeature_RNA)
minCutoff <- unname(iqr)[1]
maxCutoff <- unname(iqr)[5]

#remove cells that do not pass the feature count thresholds and mitochondrial percentage cutoff
df <- subset(df, subset = nFeature_RNA > minCutoff & nFeature_RNA < maxCutoff & percent.mt < 10) 

#############################################
# Normalize the data
#############################################

# Normalize the data using a global-scaling normalization method "LogNormalize" that normalizes the 
# feature expression measurements for each cell by the total expression, multiplies this by a scale factor 
# (10,000 by default), and log-transforms the result. Normalized values are stored in df[["RNA"]]@data
df <- NormalizeData(df)



#############################################
# Feature Selection
#############################################

# Calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in
# some cells, and lowly expressed in others). Focusing on these genes  helps to highlight biological signal in single-cell datasets.
# directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures function. 
#By default, it returns 2,000 features per dataset. These will be used in downstream analysis, like PCA.
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(df), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#create and save plot
p <- plot2 + 
  labs(title = paste("Expression Vs. Variance in ", isoMethod, sep = ""))
p

#save the figure
plotName <- paste("./figures/", isoMethod, "/", isoMethod, "_VariableFeaturePlot.png", sep = "")
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
# this are stored in df[["RNA"]]@scale.data

#get gene list from data
all.genes <- rownames(df)

#scale the data for all genes
df <- ScaleData(df, features = all.genes)



#############################################
#Perform linear dimensional reduction
#############################################

# Next we perform PCA on the scaled data. By default, only the previously determined 
#variable features are used as input, but can be defined using features argument 
#if you wish to choose a different subset.

#run the PCA and add it to the Seurat Object
df <- RunPCA(df, features = VariableFeatures(object = df))


# Examine and visualize PCA results a few different ways
# this will print the top 5 positive and negative genes
# associated with each PC
print(df[["pca"]], dims = 1:5, nfeatures = 5)

# plot the top loadings for genes associated with PC1 and PC2
# only plotting dims 1:4 because it looks better 
p <- VizDimLoadings(df, dims = 1:4, reduction = "pca")
p

#plot the PCA using PC1 and PC2
p <- DimPlot(df, reduction = "pca") + 
  labs(x = "PC1",
       y = "PC2",
       title = paste("PCA for ", isoMethod, sep = "")) +
  mytheme

p

#save the figure
plotName <- paste("./figures/", isoMethod, "/", isoMethod, "_PCA.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 9, 
       units = "in")



#############################################
#Create Heat map
#############################################

# DimHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, 
# and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores. Setting cells to a number 
# plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for 
# large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for 
# exploring correlated feature sets.

#plot heatmap for PC1
DimHeatmap(df, dims = 1, cells = 500, balanced = TRUE)


# if you want to plot heatmaps for the top 15 PCs
DimHeatmap(df, dims = 1:10, cells = 500, balanced = TRUE)



#############################################
# Determine the 'dimensionality' of the dataset
#############################################

#examine the dimensionality of the data using an elbow plot
#You are looking to find the dim where the plot starts to level off
#10 dims with nuclei_introns
ElbowPlot(df, ndims = 30)


#set dims based on the elbow plot above. For this dataset, it seems
#like 10 is a good cutoff for each isolation method
dims <- 10 

# will need to adjust dims based on above
df <- FindNeighbors(df, dims = 1:dims) 

# this iterates through different resolutions 
# & appends it to the meta.data in the seurat
# object
seur_obj <- list()
for (i in seq(0, 1, by = .1)) {
  seurat_resolution <- 0 + i
  # ^ iteratively incrementing resolution parameter
  obj <- tryCatch(Seurat::FindClusters(df, resolution = seurat_resolution),
                  error = function(e) NULL)
  if(is.null(obj)) next
  if(length(unique(obj$seurat_clusters)) == 1) next
  seur_obj[[paste0("res.", seurat_resolution)]] <- obj
}


#prepare the data for clustree by combining the resolution
#data generated above into one dataframe to compare
#initialize the dataframe with the base information
ct <- df

#add in resolution cluster data 
for (i in 1:length(seur_obj)) {
  tmp_obj <- seur_obj[[i]]@meta.data %>%
    select(5)
  ct@meta.data <- cbind(ct@meta.data, tmp_obj)
}

#run clustree to help with identifying optimal resolution
p <- clustree(ct) + 
  labs(title = isoMethod) +
  theme(plot.title = element_text(hjust = 0.5))

#output clustree diagram
p 

#save the figure
plotName <- paste("./figures/", isoMethod, "/", isoMethod, "_clustree", dims,".png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 8.7, 
       width = 9, 
       units = "in")

#will need to select the appropriate resolution based on clustree
#diagram. You are looking for the point where the columns do not 
#change much, for 
#Ctube: 0.5
#Accumax: 0.6
#Nuclei: 0.5
so <- seur_obj[[4]]

#create the data for UMAP
so <- RunUMAP(so, dims = 1:dims) #for Ctube

#visualize clusters using UMAP
DimPlot(so, reduction = "umap")

#save Seurat object
saveRDS(so, file = paste("./data/seurat_files/", isoMethod,".rds", sep = ""))

# Look at cluster IDs of the first 5 cells
head(Idents(so), 5)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(so, ident.1 = 1, min.pct = 0.25, logfc.threshold = 1)
head(cluster1.markers, n = 5)
dim(cluster1.markers)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster2.markers <- FindMarkers(so, ident.1 = 1, ident.2 = c(0, 3), min.pct = 0.25, logfc.threshold = 1)
head(cluster2.markers, n = 5)
dim(cluster2.markers)


# find markers for every cluster compared to all remaining cells, report only the positive ones
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
so.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#save all significant hits
write.table(so.markers, file = paste("./data/DEGs/", isoMethod,"_DEGs.tsv", sep = ""), sep = "\t")

#see antisense genes
cluster1.markers %>%
  rownames_to_column("gene") %>%
  dplyr::filter(grepl("-AS", gene)) %>%
  distinct(gene) %>% 
  dim ()

VlnPlot(so, features = c("SERPINE2", "PRSS23"))


# you can plot raw counts as well
#VlnPlot(so, features = c("S100B", "S100B"), slot = "counts", log = TRUE)

#plots expression of given genes on the umap plot
FeaturePlot(so, features = c("SERPINE2", "PRSS23"))

#Explore the top 10 genes for each cluster and create a heatmap
top10 <- so.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(so, features = top10$gene) + NoLegend()



###########################################################

#clean up workspace
rm(list = ls())
gc()
