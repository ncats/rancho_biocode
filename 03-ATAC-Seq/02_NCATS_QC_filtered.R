#load libraries
library(magrittr)
library(ggplot2)
library(Rtsne)
library(variancePartition)
library(tidyverse)


#CONSTANTS
#set a constant to attach to file names
outfileBase <- "03_NCATS_ATAC"

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))


#read in the filtered data containing the normalized counts from Active Motif
df <- read.delim("./data/deseq/014JNIH_ATAC_mergedregs_NormCount_Filtered.txt", header = T, row.names = 1, sep = "\t")

#create a design to tell DESeq2 which samples 
#belong to which treatments
doe <- data.frame("Sample" = colnames(df))
doe$Treatment <- sub("\\.[^.]*$", "", doe$Sample)
rownames(doe) <- doe$Sample
doe %<>% select(Treatment)

#create a factor for treatment to color PCA and tSNE
group_cols <- apply(doe, 2, function(x){length(unique(x))})


#run a PCA on the normalized data
pca <- prcomp(t(df), scale = TRUE, center = TRUE, retx = TRUE)

#get information for the top components of the PCA
if(ncol(pca$x) > 10){
  pcaNum <- 10
  topComponents <- as.data.frame(pca$x[,1:pcaNum])
} else {
  pcaNum <- ncol(pca$x)
  topComponents <- as.data.frame(pca$x[,1:pcaNum])
}

#verify that the information in the doe is the same as in PCA 
all(rownames(doe) == rownames(topComponents))

#merge treatment information into the PCA data
topComponents <- cbind(topComponents, doe)
topComponents$Treatment <- factor(topComponents$Treatment, levels = c("D0", "A1", "LSB", "D30", "D50"))

#get percentages from PCA
summary(pca)

#for each factor (only 1 here) plot the two two components and save the file
for(i in names(group_cols)){
  p <- ggplot(topComponents, aes_string("PC1", "PC2", color = i)) +
    geom_point(size=2) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_hline(yintercept = 0, color = "gray") +
    labs(x = "PC1 (56.78%)",
         y = "PC2 (14.79%)")+
    theme_classic() +
    mytheme
  print(p)
  
  plotName <- paste(outfileBase, "PCA", i, "png", sep = ".")
  ggsave(filename = paste("./figures/pca/", plotName, sep = ""), 
         plot =  p, 
         height = 6, 
         width = 8, 
         units = "in")
}


#take an iterative approach to select the highest possible perplexity for the data
#by starting at 50 and going down to 2, using a suggested range from the published data
#then run tSNE with that perplexity.
for(p in 50:2){
  successOrFailure <- try({
    tsneOut <- Rtsne(t(df), perplexity = p, check_duplicates = FALSE, pca = FALSE, max_iter = 5000, verbose = TRUE)
  }, silent = TRUE)
  
  if(class(successOrFailure) != "try-error"){
    if(interactive()) print(paste("tSNE perplexity = ", p, sep = ""))
    break
  }
}

#get tSNE data out of the tSNE object and combine it with the treatment information
tsne <- as.data.frame(tsneOut$Y)
rownames(tsne) <- row.names(doe)
tsne <- cbind(tsne, doe)
all(rownames(doe) == rownames(tsne))

for(i in names(group_cols)){
  p <- ggplot(tsne, aes_string("V1", "V2", color = i)) +
    geom_point(size=2) +
    labs(x = "tSNE-1", y = "tSNE-2") +
    theme_classic() +
    mytheme
  print(p)
  
  plotName <- paste(outfileBase, "tSNE", i, "png", sep = ".")
  ggsave(filename = paste("./figures/tsne/", plotName, sep = ""),
         plot =  p, 
         height = 6, 
         width = 8, 
         units = "in")
}


#convert the data to a matrix to feed into rowVars to identify
#the variance for each merged region
df <- as.matrix(df)

#get variance estimates for reach merged region
v <- matrixStats::rowVars(df)

#need to convert back to a data frame to add the variance estimates
#into the data frame, and filter for low variance genes because
#variance partition cannot handle 0 variance genes. Here the cutoff
#is set at 0.000001, this is not set at > 0 because there could still be
#rounding errors that may set it to 0. 
df <- as.data.frame(df)
df$v <- v
df %<>%
  filter(v > 0.000001) %>%
  select(-v)

#variance partition needs to have a formula to identify the variance 
#attributed to treatment. If the variable is categorical, it needs to 
#be have a 1| before it. 
form <- ~(1|Treatment)

#variance partition to identify the percentage of variance that can be contributed to
#treatment. This step takes a lot of time, especially if you include more variables
out <- capture.output(
  varPart <- fitExtractVarPartModel(df, form, doe, useWeights = F)
)

#plot the output from variancePartition
plotVarPart(sortCols(varPart))
