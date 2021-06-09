# tradeSeq analysis

library(Seurat)
library(tidyverse)
library(slingshot)
library(tradeSeq)
source("ncats_09_scRNA_seq_mouse/functions.R")
DATA_DIR <- "ncats_09_scRNA_seq_mouse/data"

# read in previously prepared objects to obtain counts and curves
seur <- readRDS(file.path(DATA_DIR, "objects/filtered-merged-umap-modules.RDS"))
sling <- readRDS(file.path(DATA_DIR, "objects/filtered-merged-umap-slingshot.RDS"))

# comparison of downsampling --------------------------------------------------
# 
# We ended up using all cells in the final dataset, but I 
# still wanted to check the cell subsampled datasets to 
# determine how conserved the results were across experiments
# 
import_tradeseq <- function(n){
  print(n)
  sce <- readRDS(file.path(DATA_DIR, glue("objects/tradeseq-{n}-cells.RDS")))
  associationTest(sce) %>% 
    mutate(name=n) %>% 
    rownames_to_column('gene')
}
expts <- c('micro10','micro10000','micro20000','micro40000',
           'variable1000','variable10000','variable20000')

data <- lapply(expts, import_tradeseq) 

data %>% 
  bind_rows() %>% 
  filter(grepl("micro",name)) %>% 
  group_by(name) %>% 
  arrange(pvalue,-meanLogFC) %>% 
  mutate(rank = 1:n()) %>% 
  pivot_wider(id_cols = gene, names_from = name, values_from = rank) %>% 
  column_to_rownames("gene") %>% 
  cor(method = 'spearman') %>% 
  corrplot::corrplot.mixed()

data %>% 
  bind_rows() %>% 
  filter(grepl("var",name)) %>% 
  group_by(name) %>% 
  arrange(pvalue,-meanLogFC) %>% 
  mutate(rank = 1:n()) %>% 
  pivot_wider(id_cols = gene, names_from = name, values_from = rank) %>% 
  column_to_rownames("gene") %>% 
  cor(method = 'spearman') %>% 
  corrplot::corrplot.mixed()


# slingshot analysis 
# plot pseudotime for each lineage curve ---------------------------------------
for(this_curve in seq(1,5)){
  
  # compile a "cells" table with pseudotime, weight and umap coordiantes
  pseudotimes <- slingPseudotime(sling, na = FALSE) %>% 
    as.data.frame() %>% 
    rename_with(~ paste0(.x, "_pseudotime")) %>% 
    rownames_to_column('id')
  
  weights <- slingCurveWeights(sling) %>% 
    as.data.frame() %>% 
    rename_with(~ paste0(.x, "_weight")) %>% 
    rownames_to_column('id')
  
  cells <- seur@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% 
    rownames_to_column('id') %>% 
    inner_join(pseudotimes, by = "id") %>% 
    inner_join(weights, by = "id") 
  
  # save cells that are unweighted for this lineage
  background_cells <- filter(cells, 
                             across(all_of(glue("curve{this_curve}_weight")), 
                                    function(x) x<=0.01 ) )
  
  # main umap point plot with pseudotime
  p <- ggplot(cells, aes(UMAP_1,UMAP_2)) +
    geom_point(aes_string(color=glue("curve{this_curve}_pseudotime"), 
                          alpha=glue("curve{this_curve}_weight") ), 
               shape=19, size=2) +
    scale_alpha_identity() +
    scale_color_distiller(palette = "Spectral", 
                          direction = 1) +
    
    # add background points to show unweighted cluster context
    geom_point(data = background_cells, 
               aes(UMAP_1,UMAP_2), 
               color='grey', 
               alpha=0.6, 
               shape=19, 
               size=1.5) + 
    
    
    theme_classic()
  
  # add lineage curves to the plot
  for (i in seq_along(slingCurves(sling))) {
    curve_i <- slingCurves(sling)[[i]]
    curve_i <- curve_i$s[curve_i$ord, ]
    if(i==this_curve){
      p <- p + geom_path(data = as.data.frame(curve_i), 
                         col = "black", 
                         size = 1, 
                         linetype='solid',
                         arrow = arrow() )
    }else{
      p <- p + geom_path(data = as.data.frame(curve_i), 
                         col = "darkgray", 
                         size = 1, 
                         linetype='dotted')
    }
  }
  p
  ggsave(filename = file.path(DATA_DIR, glue("plots/ti-pseudotime-curve{this_curve}.png")), width=9, height = 8)
}


# top lineage-associated genes -------------------------------------------------
ts_obj <- readRDS(file.path(DATA_DIR, glue("objects/tradeseq-variable40000-cells.RDS")))

# write the sorted top table output
top_table <- associationTest(ts_obj, lineages=F, l2fc=log2(2)) %>% 
  rownames_to_column('gene') %>% 
  arrange(-waldStat)
write_tsv(top_table, 
          file.path(DATA_DIR, glue("tables/tradeseq-associationTest-wald-sorted.tsv")))

# plot points for wald stats of significant genes
top_table %>% 
  filter(!is.na(waldStat)) %>% 
  arrange(-waldStat) %>% 
  mutate(n=1:n()) %>% 
  mutate(l = if_else(n<10,gene,"")) %>% 
  mutate(l = if_else(meanLogFC>1000,gene,l)) %>% 
  ggplot( aes(n,waldStat)) + 
  geom_point(aes(size=meanLogFC), shape=21) +
  ggrepel::geom_label_repel(aes(label=l), nudge_y = 1, nudge_x = 20) +
  scale_y_log10() +
  theme_classic() +
  labs(title="Genes showing significant expression differences across pseudotime across all lineages",
       y = "Wald statistic",
       x="Ordered Top Results")
ggsave(file.path(DATA_DIR, glue("plots/tradeseq-associationTest-wald-results.png")), width=8, height = 5)

# heatmap of top100 significant genes
Seurat::DoHeatmap(seur, cells = sample(rownames(seur@meta.data),1000), features = top_table$gene[1:100]) 
ggsave(file.path(DATA_DIR, glue("plots/tradeseq-heatmap-associationTest-top100-wald-sorted.png")), width=12, height = 8)

# count significant genes
top_table %>% 
  group_by(pvalue<0.05) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(pct = (n/sum(n))*100)

# plot top genes smoothers by lineage ------------------------------------------
FeaturePlot(seur, reduction = 'umap', features = top_table$gene[1:12])
ggsave(file.path(DATA_DIR, glue("plots/tradeseq-featureplots.png")), 
       width=16, height = 10)

library(patchwork)
out <- lapply(seq(1,12), function(i){ 
  plotSmoothers(ts_obj, 
                counts = seur@assays$RNA@counts, 
                gene = top_table$gene[[i]]) + 
    labs(title = top_table$gene[[i]])})

wrap_plots(out) + patchwork::plot_layout(ncol = 4, guides = 'collect')
ggsave(file.path(DATA_DIR, glue("plots/tradeseq-plotSmoothers.png")), 
       width=16, height = 10)


# median cluster weights by lineage --------------------------------------------
## lets take a look at how the curves are weighted 
## based on current cluster definitions this isn't 
## included in any analysis results, just FYI
clusters <- seur@meta.data %>% 
  rownames_to_column("id") %>% 
  select(id,seurat_clusters)

pseudotimes <- slingPseudotime(sling, na = FALSE) %>% 
  as.data.frame() %>% 
  rename_with(~ paste0(.x, "_pseudotime")) %>% 
  rownames_to_column('id')

weights <- slingCurveWeights(sling) %>% 
  as.data.frame() %>% 
  rename_with(~ paste0(.x, "_weight")) %>% 
  rownames_to_column('id')

cells <- seur@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% 
  rownames_to_column('id') %>% 
  inner_join(pseudotimes, by = "id") %>% 
  inner_join(weights, by = "id") %>% 
  inner_join(clusters, by = "id") 

cells %>% 
  pivot_longer(cols = ends_with("weight")) %>% 
  group_by(seurat_clusters,name) %>% 
  summarise(median = median(value)) %>%
  ggplot(aes(name, median)) +
  geom_point(aes(color=seurat_clusters), size=3) +
  ggrepel::geom_label_repel(aes(label=seurat_clusters)) +
  theme_classic()

ggsave(file.path(DATA_DIR, glue("plots/slingshot-median-lineage-weights-by-cluster.png")), width=10, height = 6)


