# dot plots for dex genes
library(Seurat)
library(slingshot)
library(glue)

s <- readRDS("objects/filtered-merged-umap.RDS")

s <- slingshot(Embeddings(s, "umap"), clusterLabels = s$seurat_clusters)

saveRDS(s, "objects/filtered-merged-umap-slingshot.RDS")

# 
# plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
# lines(sds, lwd = 2, type = 'lineages', col = 'black')
# 
# DimPlot(s, reduction = "umap")

# toptable <- read_tsv(glue("dex-tables/{sample}.tsv")) %>% 
#   arrange(p_val_adj, desc(avg_log2FC)) %>% 
#   slice(1:20)
# 
# DimPlot(s, reduction = "umap", label = T) + theme(legend.position = "none")
# 
# p <- DotPlot(s, features = toptable$gene) + 
#   coord_flip() + 
#   labs(title=glue("Top 20 differentially expressed\nfeatures from {sample}"),
#        y="Cluster Identity",
#        x="Gene name" )+
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave(plot = p, glue("plots/top20-dotplot_{sample}.png"), height = 6)
# p
