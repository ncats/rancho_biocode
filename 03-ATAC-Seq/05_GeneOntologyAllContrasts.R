#load libraries
library(magrittr)
library(clusterProfiler)
library(openxlsx)
library(org.Hs.eg.db)
library(VennDetail)
library(tidyverse)

#CONSTANTS
#set a constant to attach to file names
#outfileBase <- "NCATS_ATAC_DESeqOnly_" #used if filtering out ChipDuo
outfileBase <- "03_NCATS_ATAC_"

#get contrasts based on files that are available
#this is not needed if you are just using the list below
#contrasts <- list.files("./data/bam/combinedResults")
#contrasts <- contrasts[grepl("_atacseq_genes.txt", contrasts)]
#contrasts <- gsub("_atacseq_genes.txt", "", contrasts)

#so that all files are not needed for contrast labels
contrasts <- c("A1_vs_D0", "D30_vs_A1", "D30_vs_D0", "D50_vs_D30", "D50_vs_D0", "LSB_vs_A1", "LSB_vs_D0")

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

#get annotation to use as a background for gene ontology
annot <- read.delim("./data/refs/NCATS_annot.txt", header = T, sep = "\t")

#initialize an excel workbook
xlwb <- createWorkbook()

#read in all results file
df <- read.delim("./results/03_NCATS_ATAC_AllCombinedResults.txt", header = T, sep = "\t")

#run enrichment analyses for each contrast
for(contrast in contrasts) {
  
  ###################################
  #select data for each contrast   
  ###################################
  
  #Take the most significant hit for each gene
  sigDf <- df %>%
    filter(Contrast == contrast & Significant == 1) %>%
    arrange(AdjP.ScoreInv)
  
  
  ###################################
  #Volcano Plot
  ###################################
  
  #volcano plot
  p <- ggplot(sigDf, aes(x = Log2FC, y = -log10(AdjP.ScoreInv), color = Source)) +
    geom_point(size = 2, alpha = 0.5) + 
    geom_hline(yintercept = -log10(0.01), color="gray7", linetype="dashed") +
    geom_vline(xintercept = c(-1,1), color="gray7", linetype="dashed") + 
    labs(x = "Log2 Fold Change",
         y = "-Log10 FDR/ScoreInversion",
         title = paste(contrast, " Volcano Plot", sep = "")) +
    #scale_color_manual(values = c("lightGray", "red")) + 
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "gray80")) +
    mytheme
  
  p
  
  #save plot
  filename <- paste("./figures/volcano/", outfileBase, contrast, "_volcano.png", sep = "")
  ggsave(filename = filename,
         plot     =  p, 
         height   = 6, 
         width    = 8, 
         units    = "in")
  
  
  #filter to remove low fc genes for plotting
  sigDf1 <- sigDf %>%
    filter(AbsLog2FC > 1)
  
  ###################################
  #Upset Plot
  ###################################
  
  #get chipduo hits
  chipDuo <- sigDf1 %>%
    filter(Source %in% c("ChipDuo", "Both")) %>%
    mutate(id = paste(Gene, Start, sep = ".")) %>%
    pull(id)
  
  
  #get macs2 hits
  macs <- sigDf1 %>%
    filter(Source %in% c("MACS2/DESeq2", "Both")) %>%
    mutate(id = paste(Gene, Start, sep = ".")) %>%
    pull(id)
  
  #create upset plot
  ven <- venndetail(list("ChipDuo" = chipDuo, "MACS2" = macs))
  
  
  #initialize the PDf
  #dev.off()
  cairo_pdf(filename = paste("./figures/venn/", outfileBase, contrast, "_upset.pdf", sep = ""), 
            family = "NotoSans-Condensed",
            width = 9, height = 5)
  print(plot(ven, type = "upset"))
  graphics.off()
  
  
  
  ###################################
  #Enrichment Analysis GO
  ###################################
  
  #run gene set enrichment using GO
  ego <- enrichGO(gene          = sigDf1$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
                  universe      = annot$ENSEMBL, 
                  keyType       = 'ENSEMBL',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01)
  
  dim(as.data.frame(ego))
  
  #make dot plot
  p <- dotplot(ego, showCategory = 10) + 
    labs(title = paste(contrast, " Gene Ontology", sep = "")) +
    mytheme
  p
  
  #save plot
  filename <- paste("./figures/GO/", outfileBase, contrast, "_dotplot.png", sep = "")
  ggsave(filename = filename, 
         plot     =  p, 
         height   = 6, 
         width    = 9, 
         units    = "in")
  
  #emapplot throws an error if there are no significant hits  
  if(nrow(as.data.frame(ego)) > 1){  
    #plot enrichment map
    p <- emapplot(ego, showCategory = 10) +
      theme(axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()) +
      mytheme +
      ylab("") + xlab("")
    p
    
    #save plot
    filename <- paste("./figures/GO/", outfileBase, contrast, "_enrichMap.png", sep = "")
    ggsave(filename = filename,
           plot     =  p, 
           height   = 6, 
           width    = 8, 
           units    = "in")
  }
  
  addWorksheet(xlwb, paste(contrast, "_Enrichment", sep = ""))
  writeData(xlwb, 
            sheet = paste(contrast, "_Enrichment", sep = ""), 
            x = as.data.frame(ego))
  
  
  
  ###################################
  #Enrichment Analysis KEGG
  ###################################
  
  #run gene set enrichment using KEGG
  kk <- enrichKEGG(gene          = sigDf1$EntrezID,
                   organism      = 'hsa',
                   universe      = as.character(annot$EntrezID), 
                   pvalueCutoff  = 0.05) 
  
  dim(as.data.frame(kk))
  
  
  #make dot plot
  p <- dotplot(kk, showCategory = 10) + 
    labs(title = paste(contrast, " KEGG", sep = "")) +
    mytheme
  p
  
  #save plot
  filename <- paste("./figures/GO/", outfileBase, contrast, "_KEGG_dotplot.png", sep = "")
  ggsave(filename = filename,
         plot     =  p, 
         height   = 6, 
         width    = 8, 
         units    = "in")
  
  
  #emapplot throws an error if there are no significant hits  
  if(nrow(as.data.frame(kk)) > 1){  
    #plot enrichment map
    p <- emapplot(kk, showCategory = 10) +
      theme(axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()) +
      mytheme +
      ylab("") + xlab("")
    p
    
    #save plot
    filename <- paste("./figures/GO/", outfileBase, contrast, "_KEGG_enrichMap.png", sep = "")
    ggsave(filename = filename,
           plot     =  p, 
           height   = 6, 
           width    = 8, 
           units    = "in")
  }
  
  addWorksheet(xlwb, paste(contrast, "_KEGG_Enrichment", sep = ""))
  writeData(xlwb, 
            sheet = paste(contrast, "_KEGG_Enrichment", sep = ""), 
            x = as.data.frame(kk))
  
}

#save excel workbook with Enrichment results
saveWorkbook(xlwb, paste("./results/", outfileBase, "Enrichment.xls", sep = ""), overwrite = T)
