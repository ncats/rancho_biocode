# Once an optimal Seurat resolution has been selected, we will run differential
# gene expression analysis in Seurat using the FindMarkers function in Seurat.
# We will also utilize a list of targeted genes to produce a dot plot of the
# average normalized expression for each gene across clusters. We will follow-up
# the DGE analysis with overrepresentation analysis using msigdb and
# clusterprofiler in R exploring the GO gene sets and KEGG.
library(Seurat)
library(tidyverse)

s <- readRDS("objects/filtered-merged-umap.RDS")

results <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write_tsv(results, "dex-tables/merged-dex-results.tsv")

results %>% 
  filter(p_val_adj<=0.05) %>% 
  mutate(sig = if_else(p_val_adj<=0.0001,"***","*" )) %>% 
  ggplot(aes(x=cluster, fill=sig)) +
  geom_bar() +
  labs(x="Cluster#",
       y="Count of significant DEX genes",
       fill="Significance level") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("plots/dex-gene-counts.png")