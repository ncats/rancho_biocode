# Ready for QC
#' Create updated plots for ATAC- and MeDIP-Seq analysis

# req'd pkgs
x <- c("openxlsx", "tidyverse", "gtools", 
       "pheatmap", "scales", "data.table", 
       "broom", "Cairo", "readxl", "gdata",
       "extrafont", "extrafontdb")
sapply(x, library, character.only = TRUE)
loadfonts()

# load custom fxns
source("./functions/save_pheatmap_pdf.R")

# set cairo text backend for cairo_pdf
CairoFonts(
  regular="NotoSans-Condensed:style=Regular",
  bold="NotoSans-Condensed:style=Bold",
  italic="NotoSans-Condensed:style=Italic",
  bolditalic="NotoSans-Condensed:style=Bold Italic, BoldItalic",
  symbol="Symbol"
)

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

#### DOTPLOTS FOR ATAC-SEQ ####
# genes to plot
genes <- c("FABP7", "FAT1", "SOX9", "CD44", 
           "SPARC", "SLC1A2", "GFAP", "AQP4")

# load merged regions and raw count file
mr <- readWorkbook("./data/014JNIH_ATAC_mergedregs_DESeq2.xlsx") %>%
  dplyr::select(c(1, 6:20, 36:38))
mr_sep <- mr %>%
  dplyr::select(c(1, 17:19)) %>%
  separate_rows(c(Gene.List, Dist.to.Start, Position), sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start))
mr_gene <- mr_sep %>%
  filter(Position == "in gene")
mr_up <- mr_sep %>%
  filter((Position == "upstream" & as.numeric(Dist.to.Start) >= -1000))
mr_down <- mr_sep %>%
  filter((Position == "downstream" & as.numeric(Dist.to.Start) <= 1000))
mr_sep <- rbind(mr_gene, mr_up, mr_down) %>%
  mutate(MR_gene = paste(Merged.Region, Gene.List, sep = "_"))
rm(mr, mr_gene, mr_up, mr_down)

# load filtered normalized count atac-seq
# peaks
counts <- fread("./data/014JNIH_ATAC_mergedregs_NormCount_Filtered.txt")  %>%
  right_join(., mr_sep, by = "Merged.Region") %>%
  filter(!is.na(NCRM5.1)) %>%
  distinct(.) 

# adjust names for samples
names(counts) <- gsub("NCRM5\\.", "", names(counts))
names(counts)[c(2:4)] <- paste("D0.", seq(1,3), sep = "")

# keep MergedRegion, sample count data &
# gene symbol; combine MergedRegion & gene
# symbol
sub_counts <- counts %>% 
  mutate(mr_gene = paste(Merged.Region, Gene.List, sep = "_")) %>%
  dplyr::select(c(2:16, 20)) %>% 
  distinct(.)

# generate a pca matrix with MergedRegion +
# gene symbol as rownames, transpose for pca
pca_matrix <- sub_counts %>% 
  # make the "gene" column become the rownames of the table
  column_to_rownames("MR_gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t() 

# perform pca
pca_matrix <- prcomp(pca_matrix)
# pull variances
pc_eigenvalues <- pca_matrix$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# pull pca scores for all dims from prcomp obj
pc_scores <- pca_matrix$x
# convert samples into column
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# plot sampels to look at distribution
# for pc1 vs pc2
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = gsub("\\..*", "", sample)))

# extract pca dims by MergedRegion_gene symbol
pc_loadings <- pca_matrix$rotation

mr_sub <- mr_sep[, -c(1:2)]
# filter to retain genes for exploring
# specific variances
pc_dot <- pc_loadings %>% 
  as_tibble(rownames = "MR_gene") %>%
  data.frame(.) %>%
  mutate(mr = gsub("_.*", "", MR_gene),
         gene = gsub(".*_", "", MR_gene)) %>%
  separate_rows(MR_gene, sep = ", ") %>%
  data.frame(.) %>%
  left_join(., mr_sub) %>%
  filter(gene %in% genes)

# visually explore pc1 vs pc2 for peaks associated
# with these targeted genes
ggplot(data = pc_dot) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "purple") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            #nudge_y = 0.005, 
            check_overlap = T,
            size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) + 
  xlab("PC1") + ylab("PC2") + theme(text = element_text(family = "NotoSans-Condensed"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.line = element_line(colour = "black"))

# extract the peak associated with the largest
# variance and gene
sub_pca_dot <- pc_dot[, c(2:3, 1)] %>%
  mutate(pc1_per = abs(PC1 * 100),
         pc2_per = abs(PC2 * 100),
         gene = gsub(".*_", "", MR_gene)) %>%
  arrange(desc(pc1_per)) %>%
  group_by(gene) %>%
  slice(1) %>%
  mutate(Merged.Region = gsub("_.*", "", MR_gene)) %>%
  dplyr::select(c(Merged.Region, gene))

# filter norm counts to select only those genes
# and peak w/highest variance
counts_dot <- counts %>%
  filter(Merged.Region %in% sub_pca_dot$Merged.Region) %>%
  separate_rows(Gene.List, sep = ", ") %>%
  filter(Gene.List %in% sub_pca_dot$gene)

# pull sample and norm counts; gather into long
# format for plotting; merge with annotation info
col_counts <- counts_dot %>%
  dplyr::select(-MR_gene) %>%
  gather(Replicate, Counts, -Merged.Region, -Gene.List, -Dist.to.Start, -Position) %>%
  mutate(Sample = gsub("\\..*", "", Replicate),
         Counts = as.numeric(Counts)) 

# extract IQR for norm count data by gene
iqr <- col_counts %>% group_by(Gene.List) %>% 
  summarise_at(vars(Counts),
               list(min=min, Q1=~quantile(., probs = 0.25),
                    median=median, Q3=~quantile(., probs = 0.75),
                    max=max))

# plot distribution facetted by gene
# not normal; go w/non-parametric stats
ggplot(col_counts, aes(Counts, fill = Sample)) + geom_histogram() + facet_wrap(~ Gene.List)

# remove LSB, left join with IQR data;
# specify as white, plum2, or purple 
# depending on accessibility & IQR cut-off
counts_access <- col_counts %>%
  filter(Sample != "LSB") %>%
  left_join(., iqr) %>%
  group_by(Sample, Gene.List, Merged.Region, Dist.to.Start, Position) %>%
  summarise(access = ifelse(median(Counts) < Q1, "white",
                            ifelse(median(Counts) > Q3, "purple",
                                   ifelse(median(Counts) >= Q1 & median(Counts) <= Q3, "plum2", "NA")))) %>%
  arrange(Sample, access) %>%
  mutate(annot = paste("In gene ", Gene.List, sep = ""))

# for specifying order in plots
counts_access$Sample <- factor(counts_access$Sample, 
                               levels = c("D0", "A1", "D30", "D50"))
# manually specified this by looking at the data
ord <- paste("In gene ", unique(counts_access$Gene.List), sep = "")
ord <- ord[c(1, 4, 6, 3, 5, 7, 2)]
counts_access$annot <- factor(counts_access$annot, 
                                  levels = ord)
counts_access$access <- factor(counts_access$access, 
                               levels = c("white", "plum2", "purple"))

###### plot atac-seq data dotplot #######
p <- ggplot(counts_access, aes(y = annot, x = Sample)) + 
  geom_point(aes(fill = access), color = "black", shape = 21, size = 6) +
  scale_fill_identity() + ylab("Gene") + xlab("Sample") +
  theme_classic() + mytheme +
  scale_y_discrete(limits = rev(levels(counts$Gene.List))) +
  scale_fill_manual(name = "Accessibility", 
                    values =c('white', 'plum2', 'purple'),
                    labels = c( 'white' = 'Inaccessible', 
                                'plum2' = 'Mid-accessible', 
                                'purple' = 'Accessible')
                    ) 
ggsave(filename = "./results/atac_seq_dotplot.png", plot = p, 
       height = 7, width = 7, units = "in", dpi = 300)
rm(mr_sep, p, pc_dot, sub_counts, sub_pca_dot, 
   col_counts, counts_access, counts_dot,
   iqr, mr_sub, genes, ord)
      
#### HEATMAP NOWAKOWSKI ET AL 2017 FOR ATAC-SEQ ####
# genes to plot
genes <- c("ZIC4", "CLIP3", "NLGN2", "CHD2", "PSAP", "FZD5", "PMEL",
           "EPHA7", "NUSAP1", "BNIP3", "BIRC5", "ASPM", "CDH6", "ZNF124",
           "KIF20A", "HIST1H4J", "IER5L", "MRPL30", "SCARNA7", "HIST1H1B",
           "HIST1H2BK", "HIST1H2BC", "ZNF788P", "HIST1H1E", "HIST1H1C",
           "HIST1H2BD", "HIST1H3J", "HIST2H2AC", "HIST1H2AG", "HIST2H2AC",
           "ACBD6", "HIST1H2AI", "HIST1H2BJ", "HIST1H1D", "FAT1", "HIST2H2BE",
           "SERINC5", "SH3BP4", "BCAT1", "CACHD1", "KCNQ2", "LDHA", "HMMR",
           "TRIM24")

# load merged regions and raw count file
mr <- readWorkbook("./data/014JNIH_ATAC_mergedregs_DESeq2.xlsx") %>%
  dplyr::select(c(1:4, 6:20, 36:38))
mr <- mr %>%
  separate_rows(c(Gene.List, Dist.to.Start, Position), sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start))
mr_gene <- mr %>%
  filter(Position == "in gene")
mr_up <- mr %>%
  filter((Position == "upstream" & Dist.to.Start<= -1000))
mr_down <- mr %>%
  filter((Position == "downstream" & Dist.to.Start <= 1000))
mr <- rbind(mr_gene, mr_up, mr_down)
names(mr) <- gsub("\\_ATAC\\_hg38\\:\\:1\\.Norm.Counts|^[0-9]{2}_", "", names(mr))
names(mr) <- gsub(".*-D", "D", names(mr))
names(mr) <- gsub("NCRM5", "D0", names(mr))
names(mr) <- gsub("-", ".", names(mr))
mr <- data.frame(mr) %>%
  mutate(mr_gene = paste(Merged.Region, Gene.List, sep = "_")) %>%
  dplyr::select(-c(Merged.Region, Gene.List))
mr_genes <- mr %>%
  dplyr::select(c(21, 19:20))
rm(mr, mr_gene, mr_up, mr_down)

# filter to retain genes for exploring
# specific variances
pc_now <- pc_loadings %>% 
  as_tibble(rownames = "mr_gene") %>%
  data.frame(.) %>%
  mutate(mr = gsub("_.*", "", mr_gene),
         gene = gsub(".*_", "", mr_gene)) %>%
  separate_rows(gene, sep = ", ") %>%
  data.frame(.) %>%
  left_join(., mr_genes) %>%
  #filter(mr_gene %in% mr_gene) %>%
  filter(gene %in% genes) %>%
  mutate(mr_gene = paste(mr, gene, sep = "_")) %>%
  filter(!is.na(Position))

# visually explore pc1 vs pc2 for peaks associated
# with these targeted genes
ggplot(data = pc_now) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "purple") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            #nudge_y = 0.005, 
            check_overlap = T,
            size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) + 
  xlab("PC1") + ylab("PC2") + theme(text = element_text(family = "NotoSans-Condensed"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.line = element_line(colour = "black"),
                                    panel.background = element_rect(fill = "white"))

# extract the peak associated with the largest
# variance and gene
sub_pca_now <- pc_now[, c(1:3)] %>%
  mutate(pc1_per = abs(PC1 * 100),
         pc2_per = abs(PC2 * 100),
         gene = gsub(".*_", "", mr_gene)) %>%
  arrange(desc(pc1_per)) %>%
  group_by(gene) %>%
  slice(1) %>%
  mutate(Merged.Region = gsub("_.*", "", mr_gene)) %>%
  dplyr::select(c(Merged.Region, gene))

# filter norm counts to select only those genes
# and peak w/highest variance
counts_now <- counts %>%
  filter(Merged.Region %in% sub_pca_now$Merged.Region) %>%
  separate_rows(Gene.List, sep = ", ") %>%
  filter(Gene.List %in% sub_pca_now$gene) %>%
  mutate(mr_gene = paste(Merged.Region, Gene.List, sep = "_")) %>%
  left_join(., mr_genes)

y <- counts_now %>%
  dplyr::select(c(2:10, 17)) %>%
  column_to_rownames("Gene.List")
y <- y[, c(1:3, 7:9, 4:6)]

# define the annotation by sample, treatment
annot = data.frame(
  Sample = gsub("\\.[0-9]{1}", "", names(y))
)

# set rownames of annotation df so 
# colors will plot correctly
rownames(annot) <- names(y)

gene_type <- data.frame(gsub("i", "I", counts_now$Position))
names(gene_type) <- "Location"
rownames(gene_type) <- rownames(y)

# add annotation for RG markers overlap
# label <- data.frame(
#   RG.markers.overlap = c(rep("Pan-RG", 4), rep("Pan-RG", 2), rep("oRG", 4), 
#                          rep("Pan-RG", 38), rep("Pan-RG", 18), rep("Pan-RG", 1), 
#                          rep("Pan-RG", 1), rep("RGdiv1", 1), rep("Pan-RG", 1), 
#                          rep("RGdiv1", 2), rep("RGdiv1", 2), rep("Pan-RG", 2), 
#                          rep("Pan-RG", 1), rep("Pan-RG", 4), rep("RGdiv1", 7))
#)
#nam <- rownames(y)[grepl("ASPM|BNIP3|BIRC5|CACHD1|FAT1|HIST1H1B|HIST1H1D|HIST1H1E|HIST2H2BE|HMMR|KIF20A|LDHA|NUSAP1|PSAP|TRIM24", rownames(y))]
#rownames(label) <- nam

# specify color of annotation
# specify annotation colors
treat_samp <- list(
  Location = c(`In gene` = "#00B0F6", `downstream` = "#E76BF3"),
  Sample = c('D0' = '#F8766D', 'A1' = '#A3A500', 'LSB' = '#00BF7D')
)
#set names of elements in the list
#names(treat_samp[[1]]) <- unique(gene_type$Location)
#names(treat_samp[[2]]) <- unique(annot$Sample)

# Replace title of legend w/nothing
# names(treat_samp) <- " "
# names(annot) <- " "
# names(gene_type) <- "  "
# names(treat_samp)[c(1:2)] <- "   "

# plot heatmap w/the parameters:
# blue -> red heatmap; treatment
# annotation; treatment colors;
# no border, row, colnames
p <- pheatmap(log2(y + 1), 
              scale = "column",
              legend = T,
              legend_breaks = c(2, 1, 0, -1, -2, max(2)),
              legend_labels = c(2, 1, 0, -1, -2, "Row Z-score\n\n"),
              fontsize = 12,
              color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
              cluster_cols = F,
              cluster_rows = T,
              annotation_row = gene_type,
              annotation_col = annot, 
              annotation_colors = treat_samp,
              cellwidth = 35,
              treeheight_row = 8, 
              treeheight_col = 10,
              border_color = NA,
              show_rownames = T,
              show_colnames = F)
filename <- paste("./results/atac_nowakowsi_heatmap.pdf", sep = "")
save_pheatmap_pdf(p, filename, width = 12, height = 8)
rm(annot, counts_now, gene_type, p, pc_now, 
   sub_pca_now, treat_samp, y, genes)

#### HEATMAP FOR ZHANG et al 2016: ATAC-SEQ ####
genes <- readWorkbook("./data/Copy of Fetal vs adult astrocytes 1-s2.0-S0896627315010193-mmc5.xlsx") %>%
  dplyr::select(1) %>%
  distinct(.)
names(genes) <- genes$X1[1]
genes <- genes[-1,1]

# filter to retain genes for exploring
# specific variances
pc_zhang <- pc_loadings %>% 
  as_tibble(rownames = "mr_gene") %>%
  data.frame(.) %>%
  mutate(mr = gsub("_.*", "", mr_gene),
         gene = gsub(".*_", "", mr_gene)) %>%
  separate_rows(gene, sep = ", ") %>%
  data.frame(.) %>%
  left_join(., mr_genes) %>%
  #filter(mr_gene %in% mr_gene) %>%
  filter(gene %in% genes) %>%
  mutate(mr_gene = paste(mr, gene, sep = "_")) %>%
  filter(!is.na(Position))

# visually explore pc1 vs pc2 for peaks associated
# with these targeted genes
ggplot(data = pc_zhang) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "purple") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            #nudge_y = 0.005, 
            check_overlap = T,
            size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) + 
  xlab("PC1") + ylab("PC2") + theme(text = element_text(family = "NotoSans-Condensed"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.line = element_line(colour = "black"),
                                    panel.background = element_rect(fill = "white"))

# extract the peak associated with the largest
# variance and gene
sub_pca_zhang <- pc_zhang[, c(1:3)] %>%
  mutate(pc1_per = abs(PC1 * 100),
         pc2_per = abs(PC2 * 100),
         gene = gsub(".*_", "", mr_gene)) %>%
  arrange(desc(pc1_per)) %>%
  group_by(gene) %>%
  slice(1) %>%
  mutate(Merged.Region = gsub("_.*", "", mr_gene)) %>%
  dplyr::select(c(Merged.Region, gene))

# filter norm counts to select only those genes
# and peak w/highest variance
counts_zhang <- counts %>%
  mutate(mr_gene = paste(Merged.Region, Gene.List, sep = "_")) %>%
  filter(mr_gene %in% paste(sub_pca_zhang$Merged.Region, sub_pca_zhang$gene, sep = "_"))
counts_zhang <- counts_zhang[, c(1:4, 11:20)]

y <- counts_zhang %>%
  dplyr::select(c(2:11)) %>%
  column_to_rownames("Gene.List")

# define the annotation by sample, treatment
annot = data.frame(
  Sample = gsub("\\.[0-9]{1}", "", names(y))
)

# set rownames of annotation df so 
# colors will plot correctly
rownames(annot) <- names(y)

gene_type <- data.frame(gsub("i", "I", counts_zhang$Position))
names(gene_type) <- "Location"
rownames(gene_type) <- rownames(y)

# add annotation for RG markers overlap
# label <- data.frame(
#   RG.markers.overlap = c(rep("Pan-RG", 4), rep("Pan-RG", 2), rep("oRG", 4), 
#                          rep("Pan-RG", 38), rep("Pan-RG", 18), rep("Pan-RG", 1), 
#                          rep("Pan-RG", 1), rep("RGdiv1", 1), rep("Pan-RG", 1), 
#                          rep("RGdiv1", 2), rep("RGdiv1", 2), rep("Pan-RG", 2), 
#                          rep("Pan-RG", 1), rep("Pan-RG", 4), rep("RGdiv1", 7))
#)
#nam <- rownames(y)[grepl("ASPM|BNIP3|BIRC5|CACHD1|FAT1|HIST1H1B|HIST1H1D|HIST1H1E|HIST2H2BE|HMMR|KIF20A|LDHA|NUSAP1|PSAP|TRIM24", rownames(y))]
#rownames(label) <- nam

# specify color of annotation
#samp <- hue_pal()(3)
#gen <- hue_pal()(4)
#gen <- gen[4]
# specify annotation colors
treat_samp <- list(
  Location = c(`In gene` = "#00B0F6", `downstream` = "#E76BF3"),
  Sample = c('D0' = '#F8766D', 'D30' = '#A3A500', 'D50' = '#00BF7D')
)
#set names of elements in the list
#names(treat_samp[[1]]) <- unique(gene_type$Location)
#names(treat_samp[[2]]) <- unique(annot$Sample)

# Replace title of legend w/nothing
# names(treat_samp) <- " "
# names(annot) <- " "
# names(gene_type) <- "  "
# names(treat_samp)[c(1:2)] <- "   "

# plot heatmap w/the parameters:
# blue -> red heatmap; treatment
# annotation; treatment colors;
# no border, row, colnames
p <- pheatmap(log2(y + 1), 
              scale = "column",
              main = "",
              legend = T,
              legend_breaks = c(2, 1, 0, -1, -2, max(2)),
              legend_labels = c(2, 1, 0, -1, -2, "Row Z-score\n\n"),
              fontsize = 12,
              color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
              cluster_cols = F,
              cluster_rows = T,
              annotation_row = gene_type,
              annotation_col = annot, 
              annotation_colors = treat_samp,
              cellwidth = 35,
              treeheight_row = 8, 
              treeheight_col = 8,
              border_color = NA,
              show_rownames = T,
              show_colnames = F)
filename <- paste("./results/atac_zhang_heatmap.pdf", sep = "")
save_pheatmap_pdf(p, filename, width = 12, height = 58)
keep(save_pheatmap_pdf, mytheme, sure = TRUE)

#### DOTPLOTS FOR MEDIP-SEQ ####

# genes to plot
genes <- c("EGFR", "PDGFRB", "NCOR2", "NOTCH1", "NOTCH2", "NOTCH3", "HES1",
           "HES2", "HES3", "HES4", "HES5", "HES6", "STAT1", "STAT2", "STAT3",
           "STAT4", "STAT5", "STAT6", "MTOR", "GFAP")

# load merged regions - all
mr <- readWorkbook("./data/4075NIH_Lonza_MeDIP_mergedregs.xlsx") %>%
  dplyr::select(Merged.Region, Gene.List, Dist.to.Start, Position)

# load filtered medip-seq peaks
filt_mr <- fread("./data/final_peaks_annotatr.csv")
names(filt_mr)[1] <- "Merged.Region"
filt_mr <- filt_mr %>% pull(Merged.Region)

# load filtered normalized count atac-seq
# peaks
counts <- readRDS("./data/medipseq_counts.RDS")
names(counts)[1] <- "Merged.Region"
counts <- counts %>%
  dplyr::select(c(1:15)) %>%
  left_join(., mr, by = "Merged.Region") %>%
  distinct(.) 

# adjust names for samples
names(counts) <- gsub("NCRM5\\.", "", names(counts))
names(counts)[c(2:4)] <- paste(seq(1,3), ".", "D0", sep = "")
names(counts) <- gsub("^0", "", names(counts))

# filter mr
mr <- mr %>%
  filter(Merged.Region %in% filt_mr)
mr <- mr %>%
  separate_rows(c(Gene.List, Dist.to.Start, Position), sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start))
mr_gene <- mr %>%
  filter(Position == "in gene")
mr_up <- mr %>%
  filter((Position == "upstream" & Dist.to.Start<= -1000))
mr_down <- mr %>%
  filter((Position == "downstream" & Dist.to.Start <= 1000))
mr <- rbind(mr_gene, mr_up, mr_down) 
rm(mr_gene, mr_up, mr_down, filt_mr)

counts <- counts %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start))
counts <- counts %>%
  filter(!is.na(Gene.List)) %>%
  right_join(., mr)
counts_comb <- counts %>%
  group_by(Merged.Region) %>%
  summarise(Gene.List = paste(Gene.List, collapse = ", "),
            Dist.to.Start = paste(Dist.to.Start, collapse = ", "),
            Position = paste(Position, collapse = ", "))
counts <- counts[, c(1:15)] %>%
  distinct(.) %>%
  left_join(., counts_comb)

# keep MergedRegion, sample count data &
# gene symbol; combine MergedRegion & gene
# symbol
sub_counts <- counts %>% 
  dplyr::select(c(1:16)) %>% 
  distinct(.) %>%
  mutate(MR_gene = paste(Merged.Region, Gene.List, sep = "_")) %>%
  dplyr::select(-c(Merged.Region, Gene.List))

# generate a pca matrix with MergedRegion +
# gene symbol as rownames, transpose for pca
pca_matrix <- sub_counts %>% 
  # make the "gene" column become the rownames of the table
  column_to_rownames("MR_gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t() 

# perform pca
pca_matrix <- prcomp(pca_matrix)
# pull std dev ^ 2 (ie eigenvalues)
pc_eigenvalues <- pca_matrix$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# pull pca scores for all dims from prcomp obj
pc_scores <- pca_matrix$x
# convert samples into column
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# plot sampels to look at distribution
# for pc1 vs pc2
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = gsub(".*\\.", "", sample)))

# extract pca dims by MergedRegion_gene symbol
pc_loadings <- pca_matrix$rotation
# filter to retain genes for exploring
# specific variances
pc_medip_dot <- pc_loadings %>% 
  as_tibble(rownames = "gene") %>%
  data.frame(.) %>%
  mutate(mr = gsub("_.*", "", gene),
         gene = gsub(".*_", "", gene)) %>%
  separate_rows(gene, sep = ", ") %>%
  data.frame(.) %>%
  filter(gene %in% genes) %>%
  mutate(mr_gene = paste(mr, gene, sep = "_"))

# visually explore pc1 vs pc2 for peaks associated
# with these targeted genes
ggplot(data = pc_medip_dot) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "purple") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            #nudge_y = 0.005, 
            check_overlap = T,
            size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) + 
  xlab("PC1") + ylab("PC2") + theme(text = element_text(family = "NotoSans-Condensed"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.line = element_line(colour = "black"))

# extract the peak associated with the largest
# variance and gene
sub_pca_medip_dot <- pc_medip_dot[, c(2:3, 17)] %>%
  mutate(pc1_per = abs(PC1 * 100),
         pc2_per = abs(PC2 * 100),
         gene = gsub(".*_", "", mr_gene)) %>%
  arrange(desc(pc1_per)) %>%
  group_by(gene) %>%
  slice(1) %>%
  mutate(Merged.Region = gsub("_.*", "", mr_gene)) %>%
  dplyr::select(c(Merged.Region, gene))

# filter norm counts to select only those genes
# and peak w/highest variance
counts_medip_dot <- counts %>%
  filter(Merged.Region %in% sub_pca_medip_dot$Merged.Region) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  filter(Gene.List %in% sub_pca_medip_dot$gene)

# pull sample and norm counts; gather into long
# format for plotting; merge with annotation info
col_counts <- counts_medip_dot %>%
  gather(Replicate, Counts, -Merged.Region, -Gene.List, -Dist.to.Start, -Position) %>%
  mutate(Sample = gsub(".*\\.", "", Replicate))

# extract IQR for norm count data by gene
iqr <- col_counts %>% group_by(Gene.List) %>% 
  summarise_at(vars(Counts),
               list(min=min, Q1=~quantile(., probs = 0.25),
                    median=median, Q3=~quantile(., probs = 0.75),
                    max=max))

# plot distribution facetted by gene
# not normal; go w/non-parametric stats
ggplot(col_counts, aes(Counts, fill = Sample)) + geom_histogram() + facet_wrap(~ Gene.List)

# remove LSB, left join with IQR data;
# specify as white, plum2, or purple 
# depending on accessibility & IQR cut-off
counts_access <- col_counts %>%
  filter(Sample != "LSB") %>%
  left_join(., iqr) %>%
  group_by(Sample, Gene.List, Merged.Region, Dist.to.Start, Position) %>%
  summarise(access = ifelse(median(Counts) < Q1, "white",
                            ifelse(median(Counts) > Q3, "purple",
                                   ifelse(median(Counts) >= Q1 & median(Counts) <= Q3, "plum2", "NA")))) %>%
  arrange(Sample, access) %>%
  mutate(annot = paste("In gene ", Gene.List, sep = "")) %>%
  data.frame(.)

# for specifying order in plots
counts_access$Sample <- factor(counts_access$Sample, 
                               levels = c("D0", "A1", "D30", "D50"))
# manually specified this by looking at the data
ord <- unique(paste("In gene ", counts_access$Gene.List, sep = ""))
ord <- ord[c(8:9, 12, 13, 14, 15, 16:18, 3, 6, 11, 2, 4, 1, 5, 7, 10)]
counts_access$annot <- factor(counts_access$annot, 
                                  levels = ord)
counts_access$access <- factor(counts_access$access, 
                               levels = c("white", "plum2", "purple"))

###### plot medip-seq data dotplot #######
p <- ggplot(counts_access, aes(y = annot, x = Sample)) + 
  geom_point(aes(fill = access), color = "black", shape = 21, size = 6) +
  scale_fill_identity() + ylab("Gene") + xlab("Sample") +
  theme_classic() + mytheme +
  scale_y_discrete(limits = rev(levels(counts$Gene.List))) +
  scale_fill_manual(name = "Accessibility", 
                    values =c('white', 'plum2', 'purple'),
                    labels = c( 'white' = 'Hypo-methylated', 
                                'plum2' = 'Mid-methylated', 
                                'purple' = 'Hyper-Methylated')
  ) 
ggsave(filename = "./results/medip_seq_dotplot.png", plot = p, 
       height = 7, width = 7, units = "in", dpi = 300)
rm(genes, col_counts, counts_access, counts_comb,
   iqr, p, sub_counts, sub_pca_medip_dot, ord,
   counts_medip_dot)

#### HEATMAP ZHANG ET AL 2017 FOR MEDIP-SEQ ####

genes <- readWorkbook("./data/Copy of Fetal vs adult astrocytes 1-s2.0-S0896627315010193-mmc5.xlsx") %>%
  dplyr::select(1) %>% 
  distinct(.)
names(genes) <- genes$X1[1]
genes <- genes[-1,1]

mr$MR_gene <- paste(mr$Merged.Region, mr$Gene.List, sep = "_")
# filter to retain genes for exploring
# specific variances
pc_medip_zhang <- pc_loadings %>% 
  as_tibble(rownames = "MR_gene") %>%
  data.frame(.) %>%
  mutate(mr = gsub("_.*", "", MR_gene),
         gene = gsub(".*_", "", MR_gene)) %>%
  separate_rows(gene, sep = ", ") %>%
  left_join(., mr) %>%
  data.frame(.) %>%
  filter(gene %in% genes) %>%
  mutate(mr_gene = paste(mr, gene, sep = "_")) %>%
  filter(!is.na(Gene.List))

# visually explore pc1 vs pc2 for peaks associated
# with these targeted genes
ggplot(data = pc_medip_zhang) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "purple") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            #nudge_y = 0.005, 
            check_overlap = T,
            size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) + 
  xlab("PC1") + ylab("PC2") + theme(text = element_text(family = "NotoSans-Condensed"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.line = element_line(colour = "black"),
                                    panel.background = element_rect(fill = "white"))

# extract the peak associated with the largest
# variance and gene
sub_pca_medip_zhang <- pc_medip_zhang[, c(2:3, 22)] %>%
  mutate(pc1_per = abs(PC1 * 100),
         pc2_per = abs(PC2 * 100),
         gene = gsub(".*_", "", mr_gene)) %>%
  arrange(desc(pc1_per)) %>%
  group_by(gene) %>%
  slice(1) %>%
  mutate(Merged.Region = gsub("_.*", "", mr_gene)) %>%
  dplyr::select(c(Merged.Region, gene))

# filter norm counts to select only those genes
# and peak w/highest variance
counts_medip_zhang <- counts %>%
  filter(Merged.Region %in% sub_pca_medip_zhang$Merged.Region) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  filter(Gene.List %in% sub_pca_medip_zhang$gene)

y <- counts_medip_zhang %>%
  dplyr::select(c(2:4, 11:16)) %>%
  column_to_rownames("Gene.List")

# define the annotation by sample, treatment
annot = data.frame(
  Sample = gsub("[0-9]{1,2}\\.", "", names(y))
)

# set rownames of annotation df so 
# colors will plot correctly
rownames(annot) <- names(y)

gene_type <- data.frame(gsub("i", "I", counts_medip_zhang$Position))
names(gene_type) <- "Location"
rownames(gene_type) <- rownames(y)

treat_samp <- list(
  Location = c(`In gene` = "#00B0F6", `upstream` = "#E76BF3"),
  Sample = c('D0' = '#F8766D', 'D30' = '#A3A500', 'D50' = '#00BF7D')
)

# plot heatmap w/the parameters:
# blue -> red heatmap; treatment
# annotation; treatment colors;
# no border, row, colnames
p <- pheatmap(log2(y + 1), 
               scale = "column",
               main = "",
               legend = T,
               legend_breaks = c(4, 2, 0, -2, -4, max(4)),
               legend_labels = c(4, 2, 0, -2, -4, "Row Z-score\n"),
               fontsize = 12,
               color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
               cluster_cols = F,
               cluster_rows = T,
               annotation_row = gene_type,
               annotation_col = annot, 
               annotation_colors = treat_samp,
               cellwidth = 35,
               treeheight_row = 8, 
               treeheight_col = 10,
               border_color = NA,
               show_rownames = T,
               show_colnames = F)
filename <- paste("./results/medip_zhang_heatmap.pdf", sep = "")
save_pheatmap_pdf(p, filename, width = 12, height = 58)

rm(list = ls())
gc()
