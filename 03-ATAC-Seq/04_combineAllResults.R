#load libraries
library(magrittr)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)
library(VennDetail)
library(tidyverse)
library(data.table)

GenesToRemove <- fread("./data/refs/Genes_To_Remove_filter.txt")

#CONSTANTS
outfileBase <- "03_NCATS_ATAC_"

#get contrasts based on files that are available
#this is not needed if you are just using the list below
contrasts <- list.files("./data/combined/")
contrasts <- contrasts[grepl("_atacseq_genes.txt", contrasts)]
contrasts <- gsub("_atacseq_genes.txt", "", contrasts)

#so that all files are not needed for contrast labels
contrasts <- c("A1_vs_NCRM5", "D30_vs_A1", "D30_vs_NCRM5", "D50_vs_D30", "D50_vs_NCRM5", "LSB_vs_A1", "LSB_vs_NCRM5")

mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

#the files above are large, just use this section
#write.table(annot, "./data/NCATS_annot.txt", row.names = F, sep = "\t")
annot <- read.delim("./data/refs/NCATS_annot.txt", header = T, sep = "\t")

#remove chr, start, and end, because they can be confused with the 
#merged region information, which will be different than gene annotation
annot %<>%
  dplyr::select(-Chr, -Start, -End)

#initialize data frame 
allDf <- data.frame()

#combine all results from each of the contrasts into one long file
#this includes data from ChipDuo's combine output
for(contrast in contrasts) {
  
  ###################################
  #Read in data for each contrast   
  ###################################
  
  #load data from specific contrast
  df <- read.delim(file = paste("./data/combined/", contrast, "_atacseq_genes.txt", sep = ""), header = T, sep = "\t")
  
  
  #remove parantheses around genes and distances over 10kb away
  #we can filter these later
  df %<>%
    mutate_if(is.factor, as.character) %>%
    separate_rows(Distance, Gene, sep = ",", convert = TRUE) %>%
    mutate(Gene = gsub("[()]", "", Gene),
           Distance = gsub("[()]", "", Distance),
           Distance = as.numeric(Distance),
           Source = case_when(Source == "ChipDuo, MACS2/DESeq2" ~ "Both",
                              Source == "MACS2/DESeq2, ChipDuo" ~ "Both",
                              TRUE ~ Source)) %>%
    filter(abs(Distance) <= 10000, 
           #Source != "ChipDuo", 
           !Gene %in% GenesToRemove$Gene) 
  
  #Take the most significant hit for each gene
  sigDf <- df %>%
    arrange(AdjP.ScoreInv) %>%
    mutate(Significant = case_when(AdjP.ScoreInv < 0.05 & AbsLog2FC > 1 ~ 1,
                                   TRUE ~ 0),
           Significant = as.factor(Significant)) %>%
    left_join(., annot, by = c("Gene")) %>%
    mutate(Source = case_when(Significant == 0 ~ "NS",
                              TRUE ~ Source),
           Contrast = contrast)
  
  allDf <- rbind(allDf, sigDf)
  allDf$Contrast <- gsub("NCRM5.D", "D", allDf$Contrast)
  allDf$Contrast <- gsub("NCRM5", "D0", allDf$Contrast)  
}

#write table with all results
write.table(allDf, paste("./results/", outfileBase, "AllCombinedResults.txt", sep = ""), row.names = F, sep = "\t", quote = F)
