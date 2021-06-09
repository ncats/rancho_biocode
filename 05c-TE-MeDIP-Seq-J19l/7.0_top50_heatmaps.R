#' Build top 50 heatmaps from normalized counts
#' across all contrasts

# req'd pkgs
x <- c("magrittr", "pheatmap", "tidyverse", "data.table",
       "RColorBrewer", "scales", "Cairo",  "extrafontdb",
       "extrafont", "openxlsx")
sapply(x, library, character.only = TRUE)
loadfonts()

# set cairo text backend for cairo_pdf
CairoFonts(
  regular="NotoSans-Condensed:style=Regular",
  bold="NotoSans-Condensed:style=Bold",
  italic="NotoSans-Condensed:style=Italic",
  bolditalic="NotoSans-Condensed:style=Bold Italic, BoldItalic",
  symbol="Symbol"
)

# source custom fxns
source("./functions/save_pheatmap_pdf.R")

# load initial peaks to merge in annotation from 
# active motif
orig <- readWorkbook("./data/4441NIH_J91i-MeDIP_mergedregs.xlsx") %>%
  dplyr::select(c(Merged.Region, Length, IntervalCount,
                  CGIslandCount, PromoterCount, GeneCount, 
                  Gene.List, Dist.to.Start, Position))
names(orig)[1] <- "peak"

# load deseq raw data
dat <- fread("./adj_data/deseq/filtered/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  arrange(padj, abs(log2FoldChange))
dat <- dat[c(1:50), ]
cols <- names(dat)[grepl("^[0-9]{2}\\_", names(dat))]
dat <- dat %>% 
  dplyr::select(c(MergedRegion, cols)) %>%
  dplyr::rename("peak" = "MergedRegion") %>%
  left_join(., orig) %>%
  dplyr::select(Gene.List, cols) %>%
  mutate(Gene.List = ifelse(Gene.List == "NA", NA, Gene.List)) %>%
  arrange(Gene.List)

# make 'KIF1A' peaks unique
idx <- grep('KIF1A', dat$Gene.List)
dat$Gene.List[idx[1]] <- 'KIF1A p1'
dat$Gene.List[idx[2]] <- 'KIF1A p2'

uni_na <- seq(1, sum(is.na(dat$Gene.List)))
upd_nam <- paste0("UNANNOTATED LOCUS ", uni_na)
upd_nam <- c(dat$Gene.List[c(1:46)], upd_nam)
# keep the MergedRegion & counts of
# the contrasts
dat <- dat %>%
  mutate(Gene.List = upd_nam) %>%
  column_to_rownames("Gene.List")
#dat <- dat[, c(4:6, 1:3)]

# address one treatment w/three samples
nam <- data.frame("Treatment" = gsub("^[0-9]{2}\\_|\\-[0-9]{1}$", "", names(dat)))
# define the annotation by sample, treatment
annot = data.frame(
  Treatment = gsub("\\.", "-", nam$Treatment)
)

# set rownames of annotation df so 
# colors will plot correctly
rownames(annot) <- names(dat)

# specify color of annotation
treat <- hue_pal()(2)
# specify annotation colors
treat_samp <- list(
  Treatment = treat
)
# set names of elements in the list
names(treat_samp[[1]]) <- sort(unique(annot$Treatment))
#names(treat_samp) <- " "
#names(annot) <- "  "

# adjust title for contrast
new_title <- "iPSC vs TEexp in J91i"

# plot heatmap w/the parameters:
# blue -> red heatmap; treatment
# annotation; treatment colors;
# no border, row, colnames
p <- pheatmap(log2(dat + 1), 
              color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
              annotation_col = annot, 
              annotation_colors = treat_samp,
              border_color = NA,
              show_rownames = T, 
              show_colnames = F,
              cluster_cols = F
              #main = paste("Top 50 Significant MeDIP-Seq Peak regions for\n", new_title, sep = "")
)

filename <- paste("./adj_data/plots/", new_title, "_heatmap.pdf", sep = "")
save_pheatmap_pdf(p, filename)

rm(list = ls())
gc()
