#' Perform GSEA for GO/KEGG across all contrasts

# req'd libs
x <- c("magrittr", "clusterProfiler", "openxlsx",
       "msigdbr", "tidyverse", "Cairo", "data.table")
sapply(x, library, character.only = T)

#initialize an excel workbook
xlwb <- createWorkbook()

# pull GO/KEGG genesets from msigdbr
anot <- msigdbr("Homo sapiens") %>%
  filter(gs_subcat %in% c("GO:MF", "GO:BP", "GO:CC", "CP:KEGG"))
write.csv(anot, file = "./adj_data/go_kegg_anot.csv", row.names = F)

# pull descriptions for gene sets for merging
# post gsea
descrip <- anot %>%
  dplyr::select(c(gs_exact_source, gs_name)) %>%
  distinct(.)
names(descrip) <- c("ID", "description")

# term2gene: df of term & gene
anot <- anot %>%
  dplyr::select(c(gs_exact_source, human_gene_symbol))
names(anot) <- c("term", "gene")

# load deseq files
all <- list.files("./adj_data/deseq/filtered/", full.names = T)

# remove unncessary items from naming contrasts
all_nam <- gsub(".*\\/|_deseq.*", "", all)
all_nam <- gsub("-", ".", all_nam)

# read in deseq files, filter padj < 0.05;
# pull symbol/padj; pull unique genes from
# deseq output; recombine
# THIS PULLS ANNOTATIONS FROM ACTIVE MOTIF ONLY!
all <- lapply(all, function(x) {
  y <- fread(x) %>%
    filter(padj < 0.05) %>%
    dplyr::select(c(MergedRegion, gene_active, padj)) %>%
    arrange(padj)
  
  genes <- y %>%
    group_by(MergedRegion) %>%
    separate_rows(gene_active, sep = "; ", convert = TRUE) %>%
    filter(gene_active != "NA") %>%
    distinct(.) %>%
    mutate(symbol = gene_active) %>%
    dplyr::pull(symbol)
})

# run gsea using enricher + GO/KEGG;
# put results in df; merge w/descrip
gsea_all <- lapply(1:length(all), function(x) {
  en <- enricher(gene = all[[x]], TERM2GENE = anot)
  res <- data.frame(en) %>%
    left_join(., descrip) %>%
    mutate(GeneRatio = as.numeric(gsub("\\/.*", "", GeneRatio))/ as.numeric(gsub(".*\\/", "", GeneRatio)))
})

# dotplot of sign gene sets by contrast
lapply(1:length(gsea_all), function(x) {
  y <- gsea_all[[x]]
  if (nrow(y) == 0) {
    return(NULL)
  } else {
    y <- y[c(1:10),]
    CairoPDF(file = paste("./adj_data/plots/gsea/gsea_", gsub("_", "_vs_", gsub("\\.", "-", all_nam[[x]])), "_dotplot.pdf", sep = ""),
             width = 15, height = 12,
             family = "Noto Sans Cond")
    print(ggplot(y, aes(x = reorder(description, Count), y = Count, fill = p.adjust)) +
            geom_point(aes(size = GeneRatio), shape = 21) + 
            scale_size_continuous(range = c(3, 6)) + 
            coord_flip() +
            scale_fill_gradient("Adjusted p-value", low = "Red", high = "Royal blue") +
            xlab("GO/KEGG pathways") + ylab("Count of genes in GO/KEGG gene sets") +
            theme_minimal() + 
            ggtitle(paste("Significant GSEA for ", gsub("_", "_vs_", gsub("\\.", "-", all_nam[[x]])), " from GO/KEGG gene sets", sep = "")) +
            theme(text = element_text(size = 12)))
    graphics.off()
  }
})

# save gsea results for ea contrast into a
# workbook
lapply(1:length(gsea_all), function(x) {
  nam <- gsub("J98i\\.", "", all_nam[[x]])
  addWorksheet(xlwb, paste(gsub("_", "_vs_", nam), "_Enrichment", sep = ""))
  writeData(xlwb, 
            sheet = paste(gsub("_", "_vs_", nam), "_Enrichment", sep = ""), 
            x = as.data.frame(gsea_all[[x]]))
})

#save excel workbook with all gsea results
saveWorkbook(xlwb, paste("./adj_data/final/", "gsea_all_contrasts.xlsx", sep = ""), overwrite = T)

rm(list = ls())
gc()
