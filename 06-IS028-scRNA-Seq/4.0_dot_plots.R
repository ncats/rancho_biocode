#load libraries
library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)


#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

#get isolation method, can be: Accumax, Ctube, or nuclei_introns. 
isoMethod <- "Ctube"

#initialize directories
dir.create(paste("./figures/", isoMethod, "/dotplots/", sep = ""), showWarnings = F)

###########################################################


# load seurat obj after you select your resolution
seur <- readRDS(paste("./data/seurat_files/",isoMethod, ".rds", sep = ""))

#chose appropriate genes for dot plot
features <- c("APOE", "FABP7", "MPZ", "MBP", "NGFR", "TUBB3", "POU4F1", "S100B",
              "SOX10", "ISL1", "TAC1", "CALCA", "P2RX3", "NTRK1", "NEFH", "TRPV1", 
              "MRGPRD", "SST", "RET", "TRPA1", "SCN10A", "SCN9A")

#make the dot plot
p <- DotPlot(seur, features = features) + 
        coord_flip() + 
        labs(title = paste("Dotplot: ", isoMethod, sep = "")) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        mytheme

p

#save the figure
plotName <- paste("./figures/", isoMethod, "/dotplots/", isoMethod, "_dotplot.png", sep = "")
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 8, 
       width = 7, 
       units = "in")



###########################################################

#clean up workspace
rm(list = ls())
gc()
