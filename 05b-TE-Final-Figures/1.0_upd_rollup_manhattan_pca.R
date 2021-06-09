#' Generate final tables (deseq), updated Manhattan plot
#' and updated PCA plots

# req'd pkgs
x <- c("tidyverse", "data.table", "openxlsx", 
       "Gviz", "biomaRt", "showtext", "ggplot2",
       "ggrepel", "gtools", "scales", "DESeq2")
sapply(x, library, character.only = TRUE)
font_add("noto_cond", "./data/Noto-Sans-Condensed/NotoSans-Condensed.ttf")
font_add("noto_bold", "./data/Noto-Sans-Condensed/NotoSans-CondensedBold.ttf")
showtext_auto()

# load custom fxns
source("./functions/upd_functions.R")
source("./functions/adj_pca.R")

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "noto_bold"), 
                 axis.text = element_text(size = 14, family = "noto_cond"),
                 axis.title = element_text(colour = "Black", size = 16, family = "noto_bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "noto_cond"),
                 legend.title = element_text(colour = "Black", size = 14, family = "noto_cond"))

###### ROLLUP DESEQ/GSEA DATA ######
# j98i
deseq <- readWorkbook("./data/deseq_filtered_all_contrasts.xlsx", sheet = 2) %>%
  dplyr::select(c(1:4, 15:20)) %>% 
  mutate(contrast = "iPSC vs TEexp",
         cell_line = "J98i")
names(deseq)[1:4] <- c("peak", "chr", "st", "end")

# load initial peaks to merge in annotation from 
# active motif
orig <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx") %>%
  dplyr::select(c(Merged.Region, Length, IntervalCount,
                  CGIslandCount, PromoterCount, GeneCount, 
                  Gene.List, Dist.to.Start, Position))
names(orig)[1] <- "peak"

deseq <- deseq %>%
  left_join(., orig) %>%
  dplyr::select(c(11:12, 2:4, 1, 13, 18:20, 14:17, 5:10))

deseq_final <- deseq %>%
  dplyr::select(c(1:10, 15:20)) %>%
  arrange(padj)
names(deseq_final) <- c("Contrast", "Cell line", "Chromosome",
                        "Start", "End", "MeDIP-Seq peak", "Length",
                        "Gene annotation", "Distance to start", 
                        "Location", "Base Mean", "Log2FC", "Log2FC SE",
                        "Statistic", "P-value", "Adjusted p-value")

# j91i
deseq_j91 <- fread("./data/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  dplyr::select(c(1:11, 16:18)) %>% 
  mutate(contrast = "iPSC vs TEexp",
         cell_line = "J91i")
deseq_j91 <- deseq_j91[, c(15:16, 8:10, 1, 11,
                           12:14, 2:7)]
names(deseq_j91) <- c("Contrast", "Cell line", "Chromosome",
                                            "Start", "End", "MeDIP-Seq peak", "Length",
                                            "Gene annotation", "Distance to start", 
                                            "Location", "Base Mean", "Log2FC", "Log2FC SE",
                                            "Statistic", "P-value", "Adjusted p-value")
deseq_final <- rbind(deseq_final, deseq_j91)

gsea <- readWorkbook("./data/gsea_all_contrasts.xlsx", sheet = 2) %>%
  mutate(contrast = "iPSC vs TEexp",
         cell_line = "J98i")
gsea <- gsea %>%
  dplyr::select(c(11:12, 1, 10, 3:9))
names(gsea) <- c("Contrast", "Cell line", "GO ID",
                 "Description", "Gene ratio", "Bg Ratio",
                 "P-value", "Adjusted p-value", "Q-value",
                 "Genes", "Count")

# roll deseq into a single workbook
xlwb <- createWorkbook()
addWorksheet(xlwb, "iPSC_vs_TEexp")
writeData(xlwb,
          sheet = "iPSC_vs_TEexp",
          x = deseq_final)
saveWorkbook(xlwb, "./adj_data/deseq_final_j98i_j91i.xlsx", overwrite = T)

# roll gsea into a single workbook (j98i only)
xlwb <- createWorkbook()
addWorksheet(xlwb, "iPSC_vs_TEexp_J98i")
writeData(xlwb,
          sheet = "iPSC_vs_TEexp_J98i",
          x = gsea)
saveWorkbook(xlwb, "./adj_data/gsea_final_j98i.xlsx", overwrite = T)
rm(deseq_final, deseq, deseq_j91, gsea, xlwb, x)
gc()

###### UPDATED MANHATTAN PLOTS #####

# load sign deseq (shrink data)
dat <- fread("./data/J98i-iPSC_TEexp_deseq_shrink_filtered.csv") %>%
  dplyr::select(c(1:10)) %>%
  dplyr::rename("peak" = "MergedRegion") %>%
  arrange(padj, abs(log2FoldChange)) %>%
  left_join(., orig) %>%
  mutate(type = ifelse(log2FoldChange < 0, "TEexp",
                       ifelse(log2FoldChange > 0, "iPSC", "NA")))
names(dat)[c(8:10)] <- c("chr", "st", "end")
# build updated manhattan plot
upd_build_manhattan(dat, "Upd_manhattan_TEexp_vs_iPSC")

rm(dat, orig)
gc()

##### UPDATED PCA PLOT: J98i and J91i #####

##### COMBINED PCA PLOT (J98i and J91i) #####
# j98 sign deseq contrasts
deseq <- readWorkbook("./data/deseq_filtered_all_contrasts.xlsx", sheet = 2) %>%
  dplyr::select(c(1:4)) %>% 
  mutate(contrast = "iPSC vs TEexp",
         cell_line = "J98i")
names(deseq)[1:4] <- c("peak", "chr", "st", "end")

# j91 deseq contrasts - all
deseq_j91l <- fread("./data/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  filter(padj < 0.05) %>%
  dplyr::select(c(1, 8:10)) %>% 
  mutate(contrast = "iPSC vs TEexp",
         cell_line = "J91l")
names(deseq_j91l)[1:4] <- c("peak", "chr", "st", "end")
names(deseq_j91l) <- gsub("^[0-9]{1,2}_", "", names(deseq_j91l))

# pull chr, st, end from j98 data
j98 <- deseq %>%
  dplyr::select(c(chr, st, end)) %>%
  arrange(chr, st)
j98_nam <- paste(j98$chr, j98$st, j98$end, sep = "_")
# put it into a genomicranges obj
j98 <- GenomicRanges::GRanges(
  seqnames = j98$chr,
  ranges = IRanges(start = j98$st,
                   end = j98$end)
)
names(j98) <- j98_nam

# pull chr, st, end from j91 data
j91 <- deseq_j91l %>%
  dplyr::select(c(chr, st, end)) %>%
  arrange(chr, st)
j91_nam <- paste(j91$chr, j91$st, j91$end, sep = "_")
# put it into a genomicranges obj
j91 <- GenomicRanges::GRanges(
  seqnames = j91$chr,
  ranges = IRanges(start = j91$st,
                   end = j91$end)
)
names(j91) <- j91_nam

# find overlap in significant peak regions b/t j98 and j91l
# for chromosome, start and end position
overlap <- GenomicRanges::findOverlaps(j98, j91)

# pull chr, st, end from j98, j91 overlap
# rename peak regions for merging
overlap2 <- data.frame("chr_j98" = names(j98)[queryHits(overlap)],
                       "chr_j91" = names(j91)[subjectHits(overlap)]) %>%
  separate(chr_j98, c("chr_j98", "st_j98", "end_j98")) %>%
  separate(chr_j91, c("chr_j91", "st_j91", "end_j91")) %>%
  distinct(.) %>%
  mutate(new_peak = c(1:length(chr_j98)))

# deseq for j98 with new peaks only
j98_peaks <- overlap2[, c(1:3, 7)] %>%
  mutate(comb = paste(chr_j98, st_j98, end_j98, sep = "_"))
# deseq for j91 with new peaks only
j91_peaks <- overlap2[, c(4:7)] %>%
  mutate(comb = paste(chr_j91, st_j91, end_j91, sep = "_"))

# pull deseq peaks from j98
deseq <- readWorkbook("./data/deseq_filtered_all_contrasts.xlsx", sheet = 2) %>%
  dplyr::select(c(1:4))
names(deseq)[1:4] <- c("peak", "chr", "st", "end")

# pull deseq peaks from j91
deseq_j91l <- fread("./data/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  filter(padj < 0.05) %>%
  dplyr::select(c(1, 8:10, 20:25))
names(deseq_j91l)[1:4] <- c("peak", "chr", "st", "end")
names(deseq_j91l) <- gsub("^[0-9]{1,2}_", "", names(deseq_j91l))

# load initial peaks to merge in annotation & cts
# from active motif
orig <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx")
names(orig)[1] <- "peak"

# merge sign deseq peaks from j98 with count/annotation data
deseq <- deseq %>%
  left_join(., orig)

# filter for norm count data and remove TEp0
cols <- names(deseq)[grepl("Norm.Counts", names(deseq))]
cols <- cols[!grepl("TEp0", cols)]

deseq <- deseq %>%
  dplyr::select(c("peak", "chr", "st", "end", cols))
names(deseq) <- gsub("^[0-9]{1,2}_|_MeDIP_.*", "", names(deseq))

# merge deseq data w/new peak naming j98
deseq <- deseq %>%
  mutate(comb = paste(chr, st, end, sep = "_")) %>%
  right_join(., j98_peaks, by = "comb") %>%
  dplyr::select(c(15, 5:10))
# merge deseq data w/new peak naming j91
deseq_j91l <- deseq_j91l %>%
  mutate(comb = paste(chr, st, end, sep = "_")) %>%
  right_join(., j91_peaks, by = "comb") %>%
  dplyr::select(c(15, 5:10))

# join both deseq from j98 and j
deseq <- deseq %>%
  full_join(., deseq_j91l, by = "new_peak") %>%
  distinct(.)

# subset to keep only count data
cts <- data.frame(deseq[, c(2:13)])
nam <- deseq$new_peak
rownames(cts) <- nam

cols <- names(cts)
# rename contrasts for ease downstream
coldata <- data.frame(treatment = gsub("\\.", " ", cols))
coldata$treatment <- gsub("\\s[0-9]{1}$", "", coldata$treatment)
rownames(coldata) <- colnames(cts)
cts <- as.matrix(cts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(cts, coldata, ~treatment)
# variance stabilizing transformation
rld <- vst(dds, blind = T)

# plot pca plot with grouping of treatment
p <- adj_pca(rld, intgroup = "treatment", level = c("J91i iPSC", "J98i iPSC", "J91i TEexp", "J98i TEexp")) + theme_minimal() + guides(col=guide_legend("Treatment")) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") + 
  mytheme
ggsave(filename = "./adj_data/plots/pca_j98_j91.png", plot = p,
       width = 6, height = 4, units = "in")

##### individual pca plots for (j98 & j91) #####
# load medip/annotation data
data <- fread("./data/j98_final_peaks_annotatr.csv")
counts <- readRDS("./data/j98_medipseq_counts.RDS") %>%
  dplyr::select(-c(Chromosome, Start, End))
# merge w/count data
data <- data %>% left_join(., counts, by = "MergedRegion")

# subset to keep only count data
counts <- data.frame(data[, c(24:26, 30:32)])
nam <- data$MergedRegion
names(counts) <- gsub("^X[0-9]{1,2}_", "", names(counts))
rownames(counts) <- nam

cols <- names(counts)
# rename contrasts for ease downstream
coldata <- data.frame(treatment = ifelse(grepl("iPSC", cols), "J98i-iPSC",
                                                ifelse(grepl("TEexp", cols), "J98i-TEexp", "NA")))
rownames(coldata) <- colnames(counts)
counts <- as.matrix(counts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(counts, coldata, ~treatment)
rld <- vst(dds, blind = T)
p <- adj_pca(rld, intgroup = "treatment", level = c("J98i-iPSC", "J98i-TEexp")) + 
  theme_minimal() + guides(col=guide_legend("Treatment")) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") + 
  mytheme
ggsave(filename = "./adj_data/plots/pca_j98_only.png", plot = p,
       width = 6, height = 4, units = "in")
pdf(file = "./adj_data/plots/pca_j98_only.pdf", width = 6, height = 4)
p
graphics.off()

# load medip/annotation data
data <- readRDS("./data/j91_medipseq_counts.RDS")

# subset to keep only count data
cols <- names(data)[grepl("J91i", names(data))]
counts <- data.frame(data[, cols])
counts <- counts[, c(1:3, 7:9)]
nam <- data$MergedRegion
names(counts) <- gsub("^X", "", names(counts))
rownames(counts) <- nam

cols <- names(counts)
# rename contrasts for ease downstream
coldata <- data.frame(treatment = ifelse(grepl("iPSC", cols), "J91i-iPSC",
                                                ifelse(grepl("TEexp", cols), "J91i-TEexp", "NA")))
rownames(coldata) <- colnames(counts)
counts <- as.matrix(counts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(counts, coldata, ~treatment)
rld <- vst(dds, blind = T)
p <- adj_pca(rld, intgroup = "treatment", level = c("J91i-iPSC", "J91i-TEexp")) + 
  theme_minimal() + guides(col=guide_legend("Treatment")) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") + 
  mytheme
ggsave(filename = "./adj_data/plots/pca_j91_only.png", plot = p,
       width = 6, height = 4, units = "in")
pdf(file = "./adj_data/plots/pca_j91_only.pdf", width = 6, height = 4)
p
graphics.off()

rm(list = ls())
gc()