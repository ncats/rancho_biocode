# dot plots for dex genes
library(Seurat)
library(tidyverse)
library(glue)
library(patchwork)

seur <- readRDS(glue("ncats_09_scRNA_seq_mouse/data/objects/filtered-merged-umap.RDS"))

top_dotplot <- function(obj = seur, gene_list, name){
  p <- DotPlot(seur, features = gene_list) + 
    coord_flip() + 
    labs(title=glue("Top differentially expressed features"),
         y="Cluster Identity",
         x="Gene name") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = glue("ncats_09_scRNA_seq_mouse/data/plots/dex-top-dotplot_cluster-{unique(name)}.png"), 
         plot = p,
         height = 8, width = 8)
  list(p)
}

top_plots <- read_tsv(glue("ncats_09_scRNA_seq_mouse/data/tables/merged-dex-results.tsv")) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>% 
  select(cluster, gene) %>% 
  group_by(cluster) %>% 
  slice(1:20) %>% 
  group_by(cluster) %>% 
  summarize(dot_plot = top_dotplot(obj = seur, gene_list=gene, name=cluster))

wrap_plots(top_plots$dot_plot) + patchwork::plot_layout(ncol = 3)
ggsave("ncats_09_scRNA_seq_mouse/data/plots/dex-top-dotplot_COMBINED.png", height = 12, width = 20)

