# overexpression analysis

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(glue)

s <- readRDS("objects/filtered-merged-umap.RDS")

# msigdbr_show_species()
kegg_db <- msigdbr(species = "Mus musculus", subcategory = "CP:KEGG") %>%
  select(gs_name, gene_symbol) 

go_db <- msigdbr(species = "Mus musculus", category = "C5") %>%
  select(gs_name, gene_symbol) 

ipsc_db <- markers <- read_tsv("iPSC-modules/iPSC-mouse-markers.tsv") %>% 
  # filter for valid symbols
  filter(!is.na(mouse)) %>% 
  filter(mouse %in% rownames(s@assays$RNA) ) %>% 
  
  # clean up module names
  mutate(name = make.names(name)) %>% 
  # prefix ipsc module names
  mutate(name = str_to_upper( glue("ipsc_{name}")) ) %>% 
  transmute(gs_name = name, gene_symbol = mouse)

gs_query <- bind_rows(kegg_db,go_db,ipsc_db)

run_enrichr <- function(x, cats){
  e <- enricher(x, TERM2GENE = cats)
  e %>% as.data.frame()
}

enrich_results <- read_tsv(glue("dex-tables/merged-dex-results.tsv"), 
                           col_types = cols(.default = col_double(),
                                            gene = col_character())) %>% 
  group_by(cluster) %>%
  summarise(hits = n(),
            genes = list(c(gene))) %>%
  rowwise() %>% 
  mutate(results = list(run_enrichr(genes, cats = gs_query)) ) %>%
  unnest(results)

enrich_results %>% 
  select(-genes) %>% 
  write_tsv("dex-tables/enriched-sets.tsv")



ann <- enrich_results %>% 
  mutate(name = glue("Cluster{cluster}")) %>% 
  # mutate(name = cluster) %>% 
  select(cluster, name) %>% 
  distinct() %>% 
  column_to_rownames("name")

# KEGG 
data <- enrich_results %>% 
  mutate(name = glue("Cluster{cluster}"))  %>% 
  select(name,ID,p.adjust) %>% 
  filter(grepl("KEGG",ID)) %>% 
  pivot_wider(names_from = name, values_from = p.adjust, values_fill = 1)  %>% 
  column_to_rownames("ID") %>%
  as.matrix() 

pheatmap::pheatmap(data,
                   annotation_col = ann, 
                   cluster_rows = T,
                   color = c("blue","white"),
                   angle_col = 45, 
                   # border_color = NA,
                   fontsize_col = 6,
                   fontsize_row = 5,
                   height = 10, width = 10, 
                   filename = "plots/kegg-heatmap.png")

# GO 
# way too bug!
data <- enrich_results %>% 
  mutate(name = glue("Cluster{cluster}"))  %>% 
  select(name,ID,p.adjust) %>% 
  filter(grepl("GO",ID)) %>% 
  pivot_wider(names_from = name, values_from = p.adjust, values_fill = 1)  %>% 
  column_to_rownames("ID") %>%
  as.matrix() 

pheatmap::pheatmap(data,
                   annotation_col = ann, 
                   cluster_rows = T,
                   color = c("blue","white"),
                   angle_col = 45, 
                   # border_color = NA,
                   fontsize_col = 6,
                   fontsize_row = 5,
                   height = 10, width = 10, 
                   filename = "plots/go-heatmap.png")

