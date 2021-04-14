# dot plots for dex genes
library(Seurat)
library(tidyverse)
library(glue)

s <- readRDS(glue("objects/filtered-merged-umap.RDS"))

toptable <- read_tsv(glue("dex-tables/merged-dex-results.tsv")) %>% 
  arrange(p_val_adj, desc(avg_log2FC)) 

summary(toptable$p_val_adj)
top20 <- toptable[1:20,]$gene

DimPlot(s, reduction = "umap", label = T) + theme(legend.position = "none")

DotPlot(s, features = top20) + 
  coord_flip() + 
  labs(title=glue("Top 20 differentially expressed features"),
       y="Cluster Identity",
       x="Gene name" )+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(glue("plots/top20-dotplot_merged.png"), height = 5, width = 8)
