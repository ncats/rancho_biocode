#load libraries
library(magrittr)
library(pheatmap)
library(tidyverse)


#CONSTANTS
outfileBase <- "03_NCATS_ATAC_"

#create a list of the contrasts to loop through
contrasts <- c("A1_vs_NCRM5", "D30_vs_A1", "D30_vs_NCRM5", "D50_vs_D30", "D50_vs_NCRM5", "LSB_vs_A1", "LSB_vs_NCRM5")

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))


#function to save heatmap from pheatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  cairo_pdf(filename, width = width, height = height, family = "NotoSans_Condensed")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#make a heatmap for top 100 genes of all contrast
for(contrast in contrasts) {
  
  #############################
  #Get Data
  #############################
  
  df <- read.delim(file = paste("./data/chipduo/heatmap/counts/", contrast, "_ChipDuoCounts.txt", sep = ""), header = T, row.names = 1, sep = "\t")
  
  #change headers
  #colnames(df) <- gsub("_", ".", colnames(df))
  colnames(df) <- gsub(".*JNIH_", "", colnames(df))
  colnames(df) <- gsub("_ATAC.*", "", colnames(df))
  colnames(df) <- gsub("NCRM5.D", "D", colnames(df))
  colnames(df) <- gsub("NCRM5", "D0", colnames(df))
  
  #############################
  #Create Design
  #############################
  
  doe <- data.frame("Sample" = colnames(df))
  doe$Treatment <- gsub("\\..*", "", doe$Sample)
  rownames(doe) <- doe$Sample
  doe %<>% dplyr::select(Treatment)
  doe$Treatment <- as.factor(doe$Treatment)
  doe$Treatment <- relevel(doe$Treatment, ref = "D0")
  
  #get treatment name from the contrast name
  contrast <- gsub("NCRM5", "D0", contrast)
  treat <- strsplit(contrast, "_vs_")
  
  #get sample names for the contrast of interest
  doeContrast <- doe %>%
    filter(Treatment %in% treat[[1]])
  
  #pull samples from df
  df %<>%
    dplyr::select(rownames(doeContrast))
  
  names(doe) <- " "
  
  
  #############################
  #Heatmap
  #############################
  
  p <- pheatmap(log2(df + 1), 
                color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
                annotation_col = doe, 
                border_color = NA,
                show_rownames = F, 
                show_colnames = F,
                #cellheight = 9,
                cellwidth = 45,
                main = paste("Top 100 Significant Genes for ", gsub("_", " ", contrast), sep = ""))
  
  filename <- paste("./figures/heatmap/", outfileBase, contrast, "_heatmap.pdf", sep = "")
  save_pheatmap_pdf(p, filename)
  
}

