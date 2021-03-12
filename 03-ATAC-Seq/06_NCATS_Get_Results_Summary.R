#load libraries
library(magrittr)
library(tidyverse)

#CONSTANTS
#set a constant to attach to file names
outfileBase <- "03_NCATS_ATAC_"

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

#read in all results file
df <- read.delim("./results/03_NCATS_ATAC_AllCombinedResults.txt", header = T, sep = "\t")


#################################################
#Region level
#################################################

#calculate summary by region
resSum <- df %>%
  filter(abs(Distance) <= 10000 & AdjP.ScoreInv <= 0.05) %>%
  group_by(Contrast) %>%
  mutate(Up = sum(Log2FC > 1),
         Down = sum(Log2FC < -1)) %>%
  dplyr::select(Contrast, Up, Down) %>%
  distinct(.keep_all = T) %>%
  as.data.frame()

resSum

#write table
write.table(resSum, paste("./results/", outfileBase, "_summary.txt", sep = ""), row.names = F, sep = "\t", quote = F)



#################################################
#Gene level
#################################################

#calculate summary by gene
resSum <- df %>%
  filter(abs(Distance) <= 10000 & AdjP.ScoreInv <= 0.05) %>%
  group_by(Contrast) %>%
  arrange(Gene, AdjP.ScoreInv) %>% 
  distinct(Gene, .keep_all = T) %>%
  mutate(Up = sum(Log2FC > 1),
         Down = sum(Log2FC < -1)) %>%
  dplyr::select(Contrast, Up, Down) %>%
  distinct(.keep_all = T) %>%
  as.data.frame()

resSum


#write table
write.table(resSum, paste("./results/", outfileBase, "_gene_level_summary.txt", sep = ""), row.names = F, sep = "\t", quote = F)
