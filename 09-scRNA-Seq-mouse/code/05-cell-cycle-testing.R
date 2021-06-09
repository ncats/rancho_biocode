# explore cell cycle effect
# https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html

mouse_map <- read_tsv("reference/HMD_HumanPhenotype.tsv", 
                      col_names = c("human", "ensg", "entrez", 
                                    "id", "mouse", "accession", "note"))
symbol_map <- mouse_map$mouse
names(symbol_map) <- mouse_map$human

mouse_s_genes <- symbol_map[cc.genes.updated.2019$s.genes]
mouse_g2m_genes <- symbol_map[cc.genes.updated.2019$g2m.genes]



seur <- readRDS("objects/filtered-merged-umap.RDS")

# check out if any PCs have loadings contributed by cell cycle genes
# Looks like PC1-29 are significant
ElbowPlot(seur, ndims = 50)
pca_loads <- Loadings(seur, reduction = 'pca')[,1:29]
# convert to absolute loadings 
pca_loads <- abs(pca_loads)
# pivot and plot cell cycle markers to see if any are abnormally contributing to variance
pca_loads[rownames(pca_loads) %in% unlist(c(mouse_s_genes, mouse_g2m_genes)),] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene) %>%
  mutate(name = factor(name, levels = paste("PC",seq(1,50), sep = "_"))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_boxplot() +
  labs(title = "PC-10 shows enrichment for cell cycle gene after remediation",
       y= "Absolute PCA loading of cell cycle genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.x = element_blank())


# what are the top cell cycle genes in loading DIM
DIM<-10
pcs <- Loadings(seur, reduction = 'pca')[,DIM]

# grab top 4 cell cycle genes in PC10 for further inspection
pcs <- pcs[names(pcs)%in%c(mouse_s_genes)]
top_pcs <- pcs[ order(abs(pcs), decreasing = T) ][1:4]

# check out if this expression is driving PC10, yep
FeaturePlot(seur, reduction = "pca", dims = c(1,DIM), features = names(top_pcs)) 

str(seur@reductions$pca)




## check out if these variables were actually regressed out?
new_seur <- CellCycleScoring(seur,
                             s.features = mouse_s_genes,
                             g2m.features = mouse_g2m_genes,
                             set.ident = F)

tibble(cell = colnames(seur),
       scaled_s = new_seur@meta.data$G2M.Score,
       scaled_g2m = new_seur@meta.data$S.Score,
       original_s = seur@meta.data$G2M.Score,
       original_g2m = seur@meta.data$S.Score) %>% 
  pivot_longer(-cell) %>% 
  separate(name, into=c("version","phase")) %>% 
  mutate(version = factor(version, levels = c("original","scaled"))) %>% 
  ggplot(aes(x=version, y=value))+
  geom_point(alpha=0.2) +
  geom_line(aes(group=cell), alpha=0.2)+
  facet_wrap(~phase, scales = "free") +
  labs(title="Effect of normalization with cell cycle regression",
       y="Cell cycle score value",
       x="") +
  theme_bw()
     

# we should also see if this effect is driving any clustering
FeaturePlot(seur, reduction = "umap", features = "Birc5") 

new_seur <- AddModuleScore(seur, features = list(mouse_g2m_genes, mouse_s_genes), name = c("mouse_g2m_genes","mouse_s_genes"))
FeaturePlot(new_seur, reduction = "umap", features = "mouse_g2m_genes1") 
FeaturePlot(new_seur, reduction = "umap", features = "mouse_s_genes2") 
