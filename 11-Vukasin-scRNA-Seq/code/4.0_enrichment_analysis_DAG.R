#load libraries
library(magrittr)
library(clusterProfiler)
library(Seurat)
library(org.Hs.eg.db)
library(openxlsx)
library(msigdbr)
library(tidyverse)

#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))

#initialize directories
dir.create("./figures/all_celltypes/enrich/", showWarnings = F)
dir.create("./results/enrich/", showWarnings = F)

###########################################################

#load the significant hits
df <- read.delim("./results/DEGs/all_celltypes_DEGs.tsv", sep = "\t")

#load seurat object
seur <- readRDS("./data/seurat_files/allcells.rds") #loads all cells seurat object

#get a list of genes found in the seurat object for a background list
genes <- rownames(seur@assays$RNA[,])

#convert gene symbols to ensembl
annot <- bitr(genes, fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID", "UNIPROT"),
                OrgDb = org.Hs.eg.db)

#build db for KEGG and GO
kegg_db <- msigdbr(species = "Homo sapiens", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, human_entrez_gene) 

go_db <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  filter(gs_subcat != "HPO") %>%
  dplyr::select(gs_name, human_entrez_gene)

#combine into one db
#gs_query <- bind_rows(kegg_db,go_db)

#initialize excel workbook
xlwb <- createWorkbook()

#run enrichment analyses for each cluster
for(cluster in unique(df$cluster)) {
  
  ###################################
  #select data for each cluster   
  ###################################
  
  #Take the most significant hit for each gene and add
  #ensembl gene IDs 
  sigDf <- df %>%
    filter(cluster == !!cluster & avg_log2FC >= 1) %>%
    arrange(p_val_adj) %>%
    left_join(., annot, by = c("gene" = "SYMBOL"))
  

  ###################################
  #Enrichment Analysis GO
  ###################################
  
  #run gene set enrichment using GO
  ego <- enricher(gene = sigDf$ENTREZID, 
                  TERM2GENE = go_db,
                  universe = annot$ENTREZID,
                  pvalueCutoff = 0.05)
  
  dim(as.data.frame(ego))
  
  #if using enrichGO, simplify the output by removing overlapping categories
  #ego.bp <- simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min)
  
  #make dot plot
  p <- dotplot(ego, showCategory = 10) + 
    labs(title = paste("Cluster ", cluster, " Gene Ontology", sep = "")) +
    mytheme
  p
  
  #save plot
  filename <- paste("./figures/all_celltypes/enrich/all_celltypes_cluster", cluster, "_GO_dotplot.png", sep = "")
  ggsave(filename = filename, 
         plot     =  p, 
         height   = 6, 
         width    = 11, 
         units    = "in")
  
  #if using enrichplot newer than 1.8.0, uncomment this
  #ego <- enrichplot::pairwise_termsim(ego)
  
  #emapplot throws an error if there are no significant hits  
  if(nrow(as.data.frame(ego)) > 1){  
    #plot enrichment map
    p <- emapplot(ego, showCategory = 15) +
      theme(axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()) +
      mytheme +
      ylab("") + xlab("")
    p
    
    #save plot
    filename <- paste("./figures/all_celltypes/enrich/all_celltypes_cluster", cluster, "_GO_emap.png", sep = "")
    ggsave(filename = filename,
           plot     =  p, 
           height   = 6, 
           width    = 8, 
           units    = "in")
  }
  
  addWorksheet(xlwb, paste("cluster", cluster, "_GO_Enrichment", sep = ""))
  writeData(xlwb, 
            sheet = paste("cluster", cluster, "_GO_Enrichment", sep = ""), 
            x = as.data.frame(ego))
  

  
  ###################################
  #Enrichment Analysis KEGG
  ###################################
  
  kk <- enricher(gene = sigDf$ENTREZID, 
                 TERM2GENE = kegg_db,
                 universe = annot$ENTREZID,
                 pvalueCutoff = 0.05)

  dim(as.data.frame(kk))
  
  #if using enrichplot newer than 1.8.0, uncomment this
  #kk <- enrichplot::pairwise_termsim(kk)
  
  #make dot plot
  p <- dotplot(kk, showCategory = 15) + 
    labs(title = paste(cluster, " KEGG", sep = "")) +
    mytheme
  p
  
  #save plot
  filename <- paste("./figures/all_celltypes/enrich/all_celltypes_cluster", cluster, "_KEGG_dotplot.png", sep = "")
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
    filename <- paste("./figures/all_celltypes/enrich/all_celltypes_cluster", cluster, "_KEGG_emap.png", sep = "")
    ggsave(filename = filename,
           plot     =  p, 
           height   = 6, 
           width    = 8, 
           units    = "in")
  }
  
  addWorksheet(xlwb, paste("cluster", cluster, "_KEGG_Enrichment", sep = ""))
  writeData(xlwb, 
            sheet = paste("cluster", cluster, "_KEGG_Enrichment", sep = ""), 
            x = as.data.frame(kk))
  
}

#save excel workbook with Enrichment results
saveWorkbook(xlwb, "./results/enrich/all_celltypes_Enrichment.xls", overwrite = T)
