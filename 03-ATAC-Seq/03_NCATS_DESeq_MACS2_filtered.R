#load libraries
library(magrittr)
library(DESeq2)
library(openxlsx)
library(limma)
library(edgeR)
library(tidyverse)

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

#load merged region data to annotate reads
#this can take a bit to read in
#NOTE: you'll need to point this to the correct directory as I did not 
#want to download the whole dataset, just the relevant files. It can be found 
#in sctl/Vukasin/Vuk_ATAC/AM NIH ATAC-Seq 36311/
mRegs <- read.xlsx("./data/mergedRegions/014JNIH_ATAC_mergedregs.xlsx")

#select only the merged region, chromosome, and start and end positions
#select only those merged regions are in the 
mRegs %<>% 
  select(Merged.Region, Chromosome, Start, End) %>%
  filter(Merged.Region %in% rownames(df)) %>%
  mutate(Merged.Region = as.character(Merged.Region))

#create a design to tell DESeq2 which samples 
#belong to which treatments

#create a column with the sample names
doe <- data.frame("Sample" = colnames(df))

#remove the sample information after the first period to leave
#only the treatment information
doe$Treatment <- sub("\\.[^.]*$", "", doe$Sample)

#convert to a factor
doe$Treatment <- as.factor(doe$Treatment)

#set rownames as the sample and remove the sample column
rownames(doe) <- doe$Sample
doe %<>% select(Treatment)

#set NCRM5 as the reference group
doe$Treatment <- relevel(doe$Treatment, ref = "D0")

#manually specify contrasts because we are not interested in ALL
#pairwise comparisons
contrasts <- list()
contrasts[[1]] <- c("Treatment", "A1", "D0")
contrasts[[2]] <- c("Treatment", "LSB", "A1")
contrasts[[3]] <- c("Treatment", "LSB", "D0")
contrasts[[4]] <- c("Treatment", "D30", "D0")
contrasts[[5]] <- c("Treatment", "D50", "D0")
contrasts[[6]] <- c("Treatment", "D50", "D30")
contrasts[[7]] <- c("Treatment", "D50", "A1")
contrasts[[8]] <- c("Treatment", "D30", "A1")

#set up the DESeq object
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = doe, 
                              design= ~ Treatment)

#set a filter for the minimum number of reads per sample
#if we filter for minimum of 10 reads/sample, total is 73758
#if we filter for minimum of 20 reads/sample, total is 33364
#We should probably use 10
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]

#setting sizeFactors to 1 because this has already been done by Active Motif
#in their normalized count data file
sizeFactors(dds) <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

#run DESeq
dds <- DESeq(object = dds)

#initialize a list to add data
mylist <- list()

#loop through all of the contrasts
for (i in 1:length(contrasts)) {
  #extract results table from the DESeq analysis
  #res <- results(dds, name = i, alpha = 0.05)
  res <- results(dds, contrast = contrasts[[i]], alpha = 0.05)
  #add shrunken log2 fold change and SE to the results table
  #using an adaptive shrinkage estimator (ashr) instead of apeglm, due to using contrasts
  #this should not affect the data much, if at all, and was used because apeglm
  #could not handle contrasts
  resShrink <- lfcShrink(dds, contrast = contrasts[[i]], type = "ashr", res = res)
  #add the data to the list object created earlier
  mylist[[i]] <- as.data.frame(resShrink)
  
  #add merged region annotation to the results files and create
  #a consistent format to combine results
  resShrink %<>%
    as.data.frame() %>%
    rownames_to_column("Merged.Region") %>%
    mutate(Merged.Region = as.character(Merged.Region)) %>%
    left_join(mRegs, ., by = "Merged.Region") %>%
    filter_at(.vars = 9, all_vars(. < 0.05)) %>%
    select(Chromosome, 
           Start, 
           End, 
           starts_with("baseMean"), 
           starts_with("log2"), 
           starts_with("lfc"), 
           starts_with("pvalue"), 
           starts_with("padj")) %>%
    mutate(Chromosome = paste("chr", Chromosome, sep = ""))
  
  #set contrast name
  cn <- paste(contrasts[[i]], collapse = "_")
  
  #set file name
  fn <- paste("./results/deseq/contrasts/", outfileBase, "_", cn, "_DESeqMACS.txt", sep = "")
  
  #write the results from each contrast to their own individual files
  write.table(resShrink, file = fn, row.names = F, sep = "\t", quote = F)
}

#set the names in the list with the contrast names
for (i in seq(1:length(mylist))){
  names(mylist)[i] <- paste(contrasts[[i]], collapse = "_")
}

#combine information from all contrasts into one large data frame
ResAllContrasts <- do.call("cbind",mylist)

#remove "Treatment_" from the column names that is pasted in when combining
#the elements of the list above
colnames(ResAllContrasts) <- gsub("Treatment_", "", colnames(ResAllContrasts))

#the baseMean column from each treatment should be identical because it is
#calculated from all samples. This converts the name of the first one to
#"BaseMean" and removes the rest, then annotate with the merged region info
ResAllContrasts %<>%
  dplyr::rename(BaseMean = A1_D0.baseMean) %>%
  select(BaseMean, contains("log2FoldChange"), contains("lfcSE"), contains("pvalue"), contains("padj")) %>%
  rownames_to_column("Merged.Region") %>%
  right_join(mRegs, ., by = "Merged.Region")

#write the data to a file
write.table(ResAllContrasts, "./results/deseq/DESeqMACS2filtered.txt", row.names = T, sep = "\t")
