library(openxlsx)
library(magrittr)
library(ggplot2)
library(tidyverse)


#CONSTANTS
outfileBase <- "03_NCATS_ATAC"

mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

#get merged region data
#NOTE: you'll need to point this to the correct directory as I did not 
#want to download the whole dataset, just the relevant files. It can be found 
#in sctl/Vukasin/Vuk_ATAC/AM NIH ATAC-Seq 36311/
df <- read.xlsx("./data/mergedRegions/014JNIH_ATAC_mergedregs.xlsx")

#adjust treatment names to differentiate from NCRM5 later
colnames(df) <- gsub("NCRM5-D30", "D30", colnames(df))
colnames(df) <- gsub("NCRM5-D50", "D50", colnames(df))
colnames(df) <- gsub("NCRM5", "D0", colnames(df))


#get presence absence for each bio rep, then count number of present
#peaks for each treatment group. Then select only the peaks that are
#present in at least two bio reps
MRegs <- df %>%
  select(Merged.Region, Chromosome, Start, End, ends_with("Present")) %>%
  gather(key = Sample, value = Value, -Merged.Region, -Chromosome, -Start, -End) %>%
  mutate(Sample = gsub("_ATAC_hg38::1 Present", "", Sample),
         Treatment = gsub("-.*", "", Sample),
         Treatment = gsub("^.*_", "", Treatment)) %>%
  group_by(Treatment, Merged.Region) %>%
  mutate(PresSum = sum(Value),
         Remove = any(PresSum < 2)) %>% 
  ungroup() %>%
  filter(Remove == TRUE) %>%
  pull(Merged.Region) %>%
  unique()

#256,058 merged regions total
#227,119 merged regions after filtering for at least 2 samples in all treatments

#remove df to clean up workspace, it is not longer needed
rm(df)

#create a data frame to make a QC plot to show the number
#of merged regions that were removed after the filtering
df <- data.frame(PreFilter = 256058, PostFilter = 227119)
rownames(df) <- "NumberOfGenes"
df <- as.data.frame(t(df))
df$Filter <- rownames(df)

#plot the number of regions pre and post filtering to 
#see how many were removed
require(scales)
p <- ggplot(df, aes(x = Filter, y = NumberOfGenes)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0), labels = comma) +
  labs(x = "",
       y = "Number of Merged Regions") +
  theme_classic() + 
  mytheme

#putput plot
p

#save plot
plotName <- paste(outfileBase, "NumberAfterFiltering.png", sep = ".")
ggsave(filename = paste("./figures/qc/", plotName, sep = ""),
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")

#read in normalized counts from active motif 
#NOTE: you'll need to point this to the correct directory as I did not 
#want to download the whole dataset, just the relevant files. It can be found 
#in sctl/Vukasin/Vuk_ATAC/AM NIH ATAC-Seq diff analysis 36311/
df <- read.xlsx("./data/deseq/014JNIH_ATAC_mergedregs_DESeq2.xlsx")


#select only columns containing normalized counts and 
#merged region information, then remove the peaks 
#that did not pass the filter above
df %<>%
  select(Merged.Region, ends_with("Norm.Counts")) %>%
  filter(Merged.Region %in% MRegs) %>%
  as.data.frame()

#clean up column names by removing unnecessary information
#and replacing _ and - with . to ensure common formatting
colnames(df) <- gsub("_ATAC_hg38::1.Norm.Counts","", colnames(df))
colnames(df) <- gsub(".*_","", colnames(df))
colnames(df) <- gsub("_|-", ".", colnames(df))
colnames(df) <- gsub("NCRM5.D", "D", colnames(df))
colnames(df) <- gsub("NCRM5", "D0", colnames(df))

#write table
write.table(df, "./data/deseq/014JNIH_ATAC_mergedregs_NormCount_Filtered.txt", row.names = F, sep = "\t")
