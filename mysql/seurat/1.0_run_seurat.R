#' Process scRNA-seq data to build a Seurat object
#' to push to MySQL db for complex app

# load libraries
x <- c("Seurat", "tidyverse", "patchwork", "clustree")
sapply(x, library, character.only = TRUE)

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

###########################################################

# UPDATE PATH AND FILENAME
# load the dataset
df.data <- Read10X(data.dir = "path/to/10x-data/filename")

# UPDATE PROEJCT NAME
# initialize the seurat object with
# the raw (non-normalized data)
df <- CreateSeuratObject(counts = df.data, project = "name-your-seurat",
                         min.cells = 3, min.features = 200)

# get a summary of the seurat object
df

# calculate the % of mitochondrial genes
# and store it in a column in the meta.data
# of the seurat object
df$percent.mt <- PercentageFeatureSet(df, pattern = "^MT-")

#############################################
# qc data
#############################################

# visualize qc metrics as a violin plot
p <- VlnPlot(df,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3)

#output the plot
p

#save the figure
plotName <- "./qc_violin.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")

# FeatureScatter will plot the data
# prior to qc filtering with the
# correlation value included at the top
plot1 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt") + mytheme
plot2 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + mytheme
p <- plot1 + plot2

# output the figure
p

# save the figure
plotName <- "./qc_featurescatter.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 12, 
       units = "in")

#############################################
# subset data
#############################################

# UPDATE CUT-OFFS BASED ON YOUR DATA
# the data needs to be subset based on basic
# metrics: number of features (genes) and
# typically % mitochondrial gene content
# this uses the interquartile range of the
# number of features and a hard cut-off of
# 10% to filter the data

# get IQR
iqr <- summary(df$nFeature_RNA)
minCutoff <- unname(iqr)[1]
maxCutoff <- unname(iqr)[5]

#remove cells that do not pass the feature count thresholds and mitochondrial percentage cutoff
df <- subset(df, subset = nFeature_RNA > minCutoff & nFeature_RNA < maxCutoff & percent.mt < 10) 

#############################################
# normalize the data
#############################################

# scRNA-seq data needs to be normalized:
# we use global-scaling log normalization;
# it normalizes feature expression measurements
# for each cell by total expresion, multiplies
# by a scaling factor (10,000 default) and
# log transforms - data stored in the seurat object
#  under df[["RNA"]]@data
df <- NormalizeData(df)

#############################################
# select features (highly variable genes)
#############################################

# we calculate a subset of features that
# exhibit high cell-to-cell variation in
# the dataset (i.e, highly expressed in
# some cells and lowly expressed in others).
# focusing on these genes  helps to highlight
# biological signal in single-cell datasets.
# default = 2,000 features per dataset
df <- FindVariableFeatures(df,
                           selection.method = "vst",
                           nfeatures = 2000)

# identify the 10 most highly variable genes
top10 <- head(VariableFeatures(df), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#create and save plot
p <- plot2 + 
  labs(title = "Expression Vs. Variance")
p

#save the figure
plotName <- "./variablefeatureplot.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 9, 
       units = "in")

#############################################
# scaling the data
#############################################

# we apply a linear transformation ('scaling')
# that is a standard pre-processing
# step prior to dimensional reduction techniques
# like PCA scale (mean = 0, variance = 1)
# data is stored in df[["RNA"]]@scale.data

#get gene list from data
all.genes <- rownames(df)

#scale the data for all genes
df <- ScaleData(df, features = all.genes)

#############################################
# perform linear dimensional reduction
#############################################

# we perform pca on scaled data using only
# highly variable features - can update
# features by selecting targeted genes
df <- RunPCA(df, features = VariableFeatures(object = df))

# visualize pca results:
# this will print the top 5 positive
# and negative genes associated with
# each pc
print(df[["pca"]], dims = 1:5, nfeatures = 5)

# plot the top loadings for genes
# associated with PC1 and PC2; adjust
# dims as you would like
p <- VizDimLoadings(df, dims = 1:4, reduction = "pca")
p

# plot pca using PC1 and PC2
p <- DimPlot(df, reduction = "pca") + 
  labs(x = "PC1",
       y = "PC2") +
  mytheme

# plot the pca plot
p

#save the figure
plotName <- "./pca.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 9, 
       units = "in")

#############################################
# determine the 'dimensionality' of the dataset
#############################################

# use an elbow plot to determine optimal
# dimensions to reduce data
ElbowPlot(df, ndims = 30)


# UPDATE dims BASED ON YOUR DATA!
# set dims based on the elbow plot above
# choose dims based on the "bend" in the
# elbow plot.
dims <- 20

# will need to adjust dims based on above
df <- FindNeighbors(df, dims = 1:dims)

# UPDATE dims BASED ON YOUR DATA!
# determine what the "best" granularity of
# clustering for your dataset
seurat_res <- seq(0.4, 2.8, by = 0.4)

# this clusters based on each resolution
# and appends it to the meta.data df in
# the seurat object
df <- FindClusters(df, resolution = seurat_res)

# UPDATE dims BASED ON YOUR DATA!
# run umap dimensional reduction based
# on optimal dims
df <- RunUMAP(df, dims = 1:dims)

# target all the columns for different
# resolutions run in seurat
col_select <- names(df@meta.data)[grepl("RNA_snn_", names(df@meta.data))]
# remove the number for resolution: req'd to
# run clustree and target all columns
cols <- unique(gsub("[0-9].*", "", col_select))

# run clustree to identify optimal resolution
# ideally, you are looking for nodes where you
# don't have multiple incoming edges (arrows),
# which indicates overclustering
p <- clustree(df, prefix = cols) +
  theme(plot.title = element_text(hjust = 0.5))

# output clustree diagram
p 

# save the plot
plotName <- "./clustree.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 8.7, 
       width = 9, 
       units = "in")


# in conjunction with clustree, you need to
# explore dimensional reduction plots to select
# the optimal resolution
res <- names(df@meta.data)[grepl("res.", names(df@meta.data))]

# visualize clusters using UMAP
# use this in conjunction with
# the clustree plot to determine
# optimal resolution for clustering
p <- DimPlot(df, reduction = "umap", group.by = res)

# print all umap plots by resolution
p

# save the plot
plotName <- "./umap.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 8.7, 
       width = 9, 
       units = "in")

# UPDATE SEURAT OBJECT ONCE YOU HAVE
# REVIEWED CLUSTREE AND UMAP PLOTS
# SPECIFY YOUR OPTIMAL RESOLUTION UNDER
# opt_res
opt_res <- "0.4"
others <- names(df@meta.data)[grepl("RNA_snn_res.", names(df@meta.data))]
others <- others[!grepl(opt_res, others)]
keep_cols<- names(df@meta.data)[!names(df@meta.data) %in% others]

df@meta.data <- df@meta.data[, keep_cols]
df@meta.data$seurat_clusters <- df@meta.data[, grepl("RNA_snn_", names(df@meta.data))]

#save seurat object to push to mysql db
saveRDS(df, file = paste0("./name-your-dataset", "_seurat.RDS"))

#clean up workspace
rm(list = ls())
gc()
