#' Generate updated heatmap (j98 only), top 10 locus zoom
#' plots (j98 only) and targeted genes locus zoom plots
#' both (j98 and j91)

# req'd pkgs
x <- c("tidyverse", "data.table", "openxlsx", "ggplot2",
       "gtools", "scales", "RColorBrewer", "rtracklayer",
       "magrittr", "pheatmap", "Cairo",  "extrafontdb",
       "extrafont", "Gviz", "biomaRt")
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

# load custom fxns
source("./functions/upd_functions.R")

###### UPDATED HEATMAP PLOTS #####

# load initial peaks to merge in annotation from 
# active motif
orig <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx") %>%
  dplyr::select(c(Merged.Region, Length, IntervalCount,
                  CGIslandCount, PromoterCount, GeneCount, 
                  Gene.List, Dist.to.Start, Position))
names(orig)[1] <- "peak"

# load deseq raw data
dat <- fread("./data/J98i-iPSC_TEexp_deseq_filtered.csv") %>%
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

# label unannotated loci
uni_na <- seq(1, sum(is.na(dat$Gene.List)))
upd_nam <- paste0("UNANNOTATED LOCUS ", uni_na)
upd_nam <- c(dat$Gene.List[c(1:47)], upd_nam)
# keep the MergedRegion & counts of
# the contrasts
dat <- dat %>%
  mutate(Gene.List = upd_nam) %>%
  column_to_rownames("Gene.List")

# create annotation for col-side colors
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

# adjust title for contrast
new_title <- "iPSC vs TEexp in J98i"

# plot heatmap w/the parameters:
# blue -> red heatmap; treatment
# annotation; treatment colors;
# no border, colnames
p <- pheatmap(log2(dat + 1), 
              color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
              annotation_col = annot, 
              annotation_colors = treat_samp,
              border_color = NA,
              show_rownames = T, 
              show_colnames = F,
              cluster_cols = F)

filename <- paste("./adj_data/plots/", new_title, "_heatmap.pdf", sep = "")
save_pheatmap_pdf(p, filename)
rm(annot, dat, nam, p, treat_samp, cols, new_title, treat,
   uni_na, upd_nam, x)
gc()

###### UPDATED LOCUS-ZOOM PLOTS - TOP 10 for J98i ONLY #####
# J91i is under 05c

# for j98i only
top10 <- fread("./data/J98i-iPSC_TEexp_deseq_filtered.csv") %>%
  arrange(padj, abs(log2FoldChange))
# pull top 10 differentially methylated peaks ("in gene")
top10 <- top10[c(1:10), ]
top10 <- top10 %>% 
  dplyr::rename("peak" = "MergedRegion") %>%
  left_join(., orig) %>%
  dplyr::select(c(1,8:10, 33:35)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  filter(Position == "in gene")
names(top10)[2:7] <- c("chr", "st", "end", "gene", "dist_st", "position")

# load annotation data from active motif to 
# merge to deseq data (cleaner)
loc <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx") %>%
  data.frame(.) %>%
  dplyr::select(c(Merged.Region, Chromosome, Start, End,
                  Gene.List, Dist.to.Start, Position)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  distinct(.) %>%
  filter(Gene.List %in% top10$gene) %>%
  distinct(.) 
names(loc) <- c("mr", "chr", "start", "end", "gene",
                "dist_start", "position")
loc$dist_start <- as.numeric(loc$dist_start)

peak_top10 <- top10

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:nrow(top10), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste0("chr", top10$chr[x]))
)

########## Biomart Track #############
########## This plots the biomart annotation for
########## the selected feature

# create a mart object
ensembl <- useMart("ensembl")
# pull the homo sapiens dataset from the mart object
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
# extract data from mart object
gene_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 
                                   'start_position', 'end_position', "strand"),
                      filters = 'hgnc_symbol',
                      values = top10$gene,
                      mart = ensembl) %>%
  tibble(.) %>%
  filter(!grepl("_", chromosome_name)) 
# subtract 1000 bp up/downstream of the targeted peak
df <- gene_ensembl %>%
  mutate(start_position = start_position - 1000,
         end_position = end_position + 1000)
names(df) <- c("gene", "ensembl", "chr", "st", "end", "strand")
df<- df %>%
  mutate(chr = paste0("chr", chr))

# merge biomart data to top10 peak data
top10 <- top10 %>%
  dplyr::select(-c(chr, st, end)) %>%
  full_join(., df)

# create biomart track
bmTrack <- lapply(1:nrow(top10), function(x)
  BiomartGeneRegionTrack(genome = "hg38",
                         chromosome = top10$chr[x],
                         name = NULL,
                         stacking = "squish",
                         collapseTranscripts = "meta",
                         filters = list("ensembl_gene_id" = unlist(strsplit(top10$ensembl[x], "; "))),
                         #using fill color from UCSC genome browser
                         fill = "#0C0C78",
                         col = "#0C0C78",
                         col.line = "black",
                         background.title = "white",
                         col.title = "black",
                         fontsize = 14,
                         lwd = 1.5,
                         min.height = 10,
                         min.width = 5
  )
)

########## Annotation Track #############
########## This plots the peak features for each
########## peak region

peaks <- lapply(1:nrow(top10), function(x)
  loc %>%
    filter(position == "in gene") %>%
    filter(paste0("chr", chr) %in% top10$chr[x]) %>%
    filter(gene %in% top10$gene[x]) %>%
    filter(end <= top10$end[x])
)

aTrack <- lapply(1:length(peaks), function(x) 
  aTrack <- anno_track(st = peaks[[x]]$start, 
                       chr = paste0("chr", peaks[[x]]$chr[1]),
                       plot_name = " ",
                       wid = abs(peaks[[x]]$start - peaks[[x]]$end),
                       gen = "hg38",
                       title = "Peaks",
                       feat = "Peak ",
                       len = length(peaks[[x]]$start))
)

########## Gene Axis Track #############
########## This plots the start/end base pair locations
########## for targeted peak region

# extract st/end locations for each peak region
range_all <- lapply(1:nrow(top10), function(x) {
  range <- top10[x, ] %>%
    dplyr::slice(c(1, n())) %>%
    dplyr::select(st, end)
  range <- c(range$st[1], range$end[2])
})

# build genome axis track
gtrack <- lapply(1:length(range_all), function(x)
  GenomeAxisTrack(range=IRanges(start = range_all[[x]][[1]],
                                end = range_all[[x]][2]
  ),
  fontsize = 14,
  fontcolor = "black",
  fill.range = "black",
  col = "black",
  col.line = "black"
  ))

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull bam files
bam_files <- list.files("./data/BAM_MeDIP3/", pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("TEp0", bam_files)]
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

# calculate coverage within each peak region
# and bam track file
all_cov <- lapply(1:length(range_all), function(x) {
  unlist(lapply(c(1:7), function(y) {
    build_cov(chr = top10$chr[x],
              range1 = range_all[[x]][1],
              range2 = range_all[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})

nam <- bam_nam
nam <- ifelse(grepl("iPSC", nam), "J98i iPSC",
          ifelse(grepl("TEexp", nam), "J98i TEexp", 
            ifelse(grepl("input2", nam), "J98i input", "NA")))
nam[c(2, 4:5, 7)] <- ""

# generate data tracks
all_dat <- multi_comp_dat_t(all_cov, bam_nam, nam)

########## Highlight Track ############
########## This overlays highlighted tracks over promoter
########## and significant peak regions in the data tracks.

# extract entire st/end site of peak region (gene)
sub_medip <- top10 %>%
  dplyr::select(c(gene, chr, st, end, strand)) %>%
  dplyr::rename("full_st" = "st") %>%
  dplyr::rename("full_end" = "end") %>%
  mutate(chr = as.numeric(gsub("chr", "", chr)))
# calculate width of peak region and merge to
# full peak region above
peak_top10 <- peak_top10 %>%
  mutate(width = end - st,
         chr = as.numeric(chr)) %>%
  left_join(., sub_medip)

# overlay highlight track on data track
ht1 <- lapply(1:nrow(peak_top10), function(x) {
  if (peak_top10$strand[x] > 0) {
    HighlightTrack(trackList = all_dat[[x]], 
                   genome = "hg38",
                   start = c(peak_top10$st[x], peak_top10$full_st[x]), 
                   width = c(peak_top10$width[x], 1000), 
                   chromosome = peak_top10$chr[x],
                   lty = c("solid", "dashed"),
                   col = c("black", "black"),
                   fill = c(NA, NA),
                   inBackground = T)
  } else {
    HighlightTrack(trackList = all_dat[[x]], 
                   genome = "hg38",
                   start = c(peak_top10$st[x], peak_top10$full_end[x] - 1000), 
                   width = c(peak_top10$width[x], 1000), 
                   chromosome = peak_top10$chr[x],
                   lty = c("solid", "dashed"),
                   col = c("black", "black"),
                   fill = c(NA, NA),
                   inBackground = T)
  }
})

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- top10$gene

# plot top 10 peaks for ea contrast
plot_data2(itrack, aTrack, ht1, gtrack, bmTrack, title, "top10_j98")
rm(all_cov, all_dat, aTrack, bmTrack,
   ensembl, gene_ensembl, gtrack, ht1,
   loc, peak_top10, peaks, range_all,
   sub_medip, top10, df, loc, bam_files,
   bam_nam, nam, title, itrack)
gc()

##### TARGETED GENES FOR BOTH J98I AND J91I#####

# based on request ilyas/jaro
goi <- c("POU5F1", "SOX2", "TGF1", "UTF1",
         "ELF5", "GATA3", "TFAP2A", "CDX2",
         "TFAP2C", "YAP1", "TP63", "XAGE3",
         "HAVCR1", "MAGEA10", "HLA-A", "HLA-B")

# filter to retain peaks that are in gene, 
#  up/downstream within 10 kb
targ <- fread("./data/J98i-iPSC_TEexp_deseq_filtered.csv") %>%
  arrange(padj, abs(log2FoldChange)) %>% 
  dplyr::rename("peak" = "MergedRegion") %>%
  left_join(., orig) %>%
  dplyr::select(c(1,8:10, 33:35)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start)) %>%
  filter(Gene.List %in% goi) %>%
  filter(Position == "in gene" | (Position == "downstream" & Dist.to.Start <= 10000) | (Position == "upstream" & Dist.to.Start >= -10000)) %>%
  # these two were annotated in gene, but were not in the gene
  filter(peak != 113813) %>%
  filter(peak != 102166) %>%
  arrange(Gene.List)
names(targ)[2:7] <- c("chr", "st", "end", "gene", "dist_st", "position")

# filter to retain sign peaks that are in gene, 
#  up/downstream within 10 kb
sign_targ <- fread("./data/J98i-iPSC_TEexp_deseq_filtered.csv") %>%
  filter(padj < 0.05) %>%
  arrange(padj, abs(log2FoldChange)) %>% 
  dplyr::rename("peak" = "MergedRegion") %>%
  left_join(., orig) %>%
  dplyr::select(c(1,8:10, 33:35)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start)) %>%
  filter(Gene.List %in% goi)%>%
  filter(Position == "in gene" | (Position == "downstream" & Dist.to.Start <= 10000) | (Position == "upstream" & Dist.to.Start >= -10000)) %>%
  filter(peak != 102166) %>%
  arrange(Gene.List)
names(sign_targ)[2:7] <- c("chr", "st", "end", "gene", "dist_st", "position")

# use this to pull clean active motif annotation
loc <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx") %>%
  data.frame(.) %>%
  dplyr::select(c(Merged.Region, Chromosome, Start, End,
                  Gene.List, Dist.to.Start, Position)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  distinct(.) %>%
  filter(Gene.List %in% goi) %>%
  distinct(.) 
names(loc) <- c("mr", "chr", "start", "end", "gene",
                "dist_start", "position")
loc$dist_start <- as.numeric(loc$dist_start)

uni_genes <- targ %>%
  dplyr::select(c("chr", "gene")) %>%
  distinct(.)
  
########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:nrow(uni_genes), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste0("chr", uni_genes$chr[x]))
)

########## Biomart Track #############
########## This plots the biomart annotation for
########## the selected feature

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
gene_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 
                                   'start_position', 'end_position', "strand"),
                      filters = 'hgnc_symbol',
                      values = uni_genes$gene,
                      mart = ensembl) %>%
  tibble(.) %>%
  filter(!grepl("_", chromosome_name)) 
#df <- data.frame(gene_ensembl) 
df <- gene_ensembl %>%
  mutate(start_position = start_position - 10000,
         end_position = end_position + 10000)
names(df) <- c("gene", "ensembl", "chr", "st", "end", "strand")
df<- df %>%
  mutate(chr = paste0("chr", chr))

uni_genes <- uni_genes %>%
  mutate(chr = paste0("chr", chr)) %>%
  full_join(., df)

bmTrack <- lapply(1:nrow(uni_genes), function(x)
  BiomartGeneRegionTrack(genome = "hg38",
                         chromosome = uni_genes$chr[x],
                         name = NULL,
                         #name = uni_genes$gene[x], #to shrink the track
                         stacking = "squish",
                         collapseTranscripts = "meta",
                         filters = list("ensembl_gene_id" = unlist(strsplit(uni_genes$ensembl[x], "; "))),
                         #using fill color from UCSC genome browser
                         fill = "#0C0C78",
                         col = "#0C0C78",
                         col.line = "black",
                         background.title = "white",
                         col.title = "black",
                         fontsize = 14,
                         lwd = 1.5,
                         min.height = 10#,
                         #min.width = 5
  )
)

########## Annotation Track #############
########## This plots the peak features for each
########## peak region

#peaks <- split(targ, targ$gene)
peaks <- lapply(1:nrow(uni_genes), function(x)
  targ %>%
    filter(paste0("chr", chr) %in% uni_genes$chr[x]) %>%
    filter(gene %in% uni_genes$gene[x]) %>%
    filter(st >= uni_genes$st[x]) %>%
    filter(end <= uni_genes$end[x])
)

aTrack <- lapply(1:length(peaks), function(x) {
  if (nrow(peaks[[x]]) == 0) {
    return(NULL)
  } else {
    aTrack <- anno_track(st = peaks[[x]]$st, 
                         chr = paste0("chr", peaks[[x]]$chr[1]),
                         plot_name = " ",
                         wid = abs(peaks[[x]]$st - peaks[[x]]$end),
                         gen = "hg38",
                         title = "Peaks",
                         feat = "Peak ",
                         len = length(peaks[[x]]$st))
  }
})

########## Gene Axis Track #############
########## This plots the start/end base pair locations
########## for targeted peak region

range_all <- lapply(1:nrow(uni_genes), function(x) {
  range <- uni_genes[x, ] %>%
    dplyr::slice(c(1, n())) %>%
    dplyr::select(st, end)
  range <- c(range$st[1], range$end[2])
})

# build genome axis track for h9
gtrack <- lapply(1:length(range_all), function(x)
  GenomeAxisTrack(range=IRanges(start = range_all[[x]][[1]],
                                end = range_all[[x]][2]
  ),
  fontsize = 14,
  fontcolor = "black",
  fill.range = "black",
  col = "black",
  col.line = "black"
  ))

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull coverage from bam files
bam_files <- list.files("./data/BAM_MeDIP3/", pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("TEp0", bam_files)]
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

all_cov <- lapply(1:length(range_all), function(x) {
  unlist(lapply(c(1:7), function(y) {
    build_cov(chr = uni_genes$chr[x],
              range1 = range_all[[x]][1],
              range2 = range_all[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})

nam <- bam_nam
nam <- ifelse(grepl("iPSC", nam), "J98i iPSC",
              ifelse(grepl("TEexp", nam), "J98i TEexp",
                     ifelse(grepl("input", nam), "J98i input", "NA")))
nam[c(2, 4:5, 7)] <- ""

all_dat <- multi_comp_dat_t(all_cov, bam_nam, nam)

########## Highlight Track ############
########## This overlays highlighted tracks over promoter
########## and significant peak regions in the data tracks.
sub <- gene_ensembl[, c(2, 4:5)]
names(sub) <- c("ensembl", "gene_st", "gene_end")
sub_medip <- uni_genes %>%
  dplyr::select(c(gene, chr, st, end, strand, ensembl)) %>%
  dplyr::rename("full_st" = "st") %>%
  dplyr::rename("full_end" = "end") %>%
  mutate(chr = as.numeric(gsub("chr", "", chr))) %>%
  left_join(., sub)
sign_targ <- sign_targ %>%
  mutate(width = end - st,
         chr = as.numeric(chr)) %>%
  full_join(., sub_medip)
sign_targ <- split(sign_targ, sign_targ$gene)

ht1 <- lapply(1:length(sign_targ), function(x) {
  if (nrow(sign_targ[[x]]) >= 1 & !is.na(sign_targ[[x]]$peak[1])) {
    if (sign_targ[[x]]$strand[1] > 0) {
      st <- c(sign_targ[[x]]$st, sign_targ[[x]]$gene_st[1] -1000)
      width <- c(sign_targ[[x]]$width, 1000)
      col <- c(rep("black", length(st) - 1), "black")
      fill <- c(rep(NA, length(st) - 1), NA)
      lty <- c(rep("solid", length(st) - 1), "dashed")
      HighlightTrack(trackList = all_dat[[x]], 
                     genome = "hg38",
                     start = st, 
                     width = width, 
                     chromosome = sign_targ[[x]]$chr[1],
                     lty = lty,
                     col = col,
                     fill = fill,
                     inBackground = T)
    } else {
      st <- c(sign_targ[[x]]$st, sign_targ[[x]]$gene_end[1])
      width <- c(sign_targ[[x]]$width, 1000)
      col <- c(rep("black", length(st) - 1), "black")
      fill <- c(rep(NA, length(st) - 1), NA)
      lty <- c(rep("solid", length(st) - 1), "dashed")
      HighlightTrack(trackList = all_dat[[x]], 
                     genome = "hg38",
                     start = st, 
                     width = width, 
                     chromosome = sign_targ[[x]]$chr[1],
                     lty = lty,
                     col = col,
                     fill = fill,
                     inBackground = T)
    }
  } else {
    st <- sign_targ[[x]]$gene_end[1]
    width <- 1000
    col <- "black"
    fill <- NA
    lty <- "dashed"
    HighlightTrack(trackList = all_dat[[x]], 
                   genome = "hg38",
                   start = st, 
                   width = width, 
                   chromosome = sign_targ[[x]]$chr[1],
                   lty = lty,
                   col = col,
                   fill = fill,
                   inBackground = T)
  }
})

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- uni_genes$gene

# plot top 10 peaks for ea contrast
plot_data2(itrack, aTrack, ht1, gtrack, bmTrack, title, "targeted_j98")

rm(all_cov, all_dat, aTrack, bmTrack, df, ensembl,
   gene_ensembl, gtrack, ht1, itrack, loc, loc,
   peaks, range_all, sign_targ, sub_medip, targ,
   uni_genes, bam_files, bam_nam, nam, title)
gc()

# j91i
# filter to retain peaks that are in gene, 
#  up/downstream within 10 kb
targ <- fread("./data/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  arrange(padj, abs(log2FoldChange)) %>% 
  dplyr::rename("peak" = "MergedRegion") %>%
  dplyr::select(c(1,8:10, 16:18)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start)) %>%
  filter(Gene.List %in% goi) %>%
  filter(Position == "in gene" | (Position == "downstream" & Dist.to.Start <= 10000) | (Position == "upstream" & Dist.to.Start >= -10000)) %>%
  arrange(Gene.List)
names(targ)[2:7] <- c("chr", "st", "end", "gene", "dist_st", "position")

# filter to retain sign peaks that are in gene, 
#  up/downstream within 10 kb
sign_targ <- fread("./data/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  filter(padj < 0.05) %>%
  arrange(padj, abs(log2FoldChange)) %>% 
  dplyr::rename("peak" = "MergedRegion") %>%
  dplyr::select(c(1,8:10, 16:18)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  mutate(Dist.to.Start = as.numeric(Dist.to.Start)) %>%
  filter(Gene.List %in% goi)%>%
  filter(Position == "in gene" | (Position == "downstream" & Dist.to.Start <= 10000) | (Position == "upstream" & Dist.to.Start >= -10000)) %>%
  #filter(peak != 102166) %>%
  arrange(Gene.List)
names(sign_targ)[2:7] <- c("chr", "st", "end", "gene", "dist_st", "position")

# use this to pull clean active motif annotation
loc <- readWorkbook("./data/4441NIH_J91i-MeDIP_mergedregs.xlsx") %>%
  data.frame(.) %>%
  dplyr::select(c(Merged.Region, Chromosome, Start, End,
                  Gene.List, Dist.to.Start, Position)) %>%
  separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  distinct(.) %>%
  filter(Gene.List %in% goi) %>%
  distinct(.) 
names(loc) <- c("mr", "chr", "start", "end", "gene",
                "dist_start", "position")
loc$dist_start <- as.numeric(loc$dist_start)

uni_genes <- targ %>%
  dplyr::select(c("chr", "gene")) %>%
  distinct(.)

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:nrow(uni_genes), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste0("chr", uni_genes$chr[x]))
)

########## Biomart Track #############
########## This plots the biomart annotation for
########## the selected feature

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
gene_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 
                                   'start_position', 'end_position', "strand"),
                      filters = 'hgnc_symbol',
                      values = uni_genes$gene,
                      mart = ensembl) %>%
  tibble(.) %>%
  dplyr::filter(!grepl("_", chromosome_name)) 
df <- gene_ensembl %>%
  mutate(start_position = start_position - 10000,
         end_position = end_position + 10000)
names(df) <- c("gene", "ensembl", "chr", "st", "end", "strand")
df<- df %>%
  mutate(chr = paste0("chr", chr))

uni_genes <- uni_genes %>%
  mutate(chr = paste0("chr", chr)) %>%
  full_join(., df)

bmTrack <- lapply(1:nrow(uni_genes), function(x)
  BiomartGeneRegionTrack(genome = "hg38",
                         chromosome = uni_genes$chr[x],
                         name = NULL,
                         #name = uni_genes$gene[x], #to shrink the track
                         stacking = "squish",
                         collapseTranscripts = "meta",
                         filters = list("ensembl_gene_id" = unlist(strsplit(uni_genes$ensembl[x], "; "))),
                         #using fill color from UCSC genome browser
                         fill = "#0C0C78",
                         col = "#0C0C78",
                         col.line = "black",
                         background.title = "white",
                         col.title = "black",
                         fontsize = 14,
                         lwd = 1.5,
                         min.height = 10#,
                         #min.width = 5
  )
)

########## Annotation Track #############
########## This plots the peak features for each
########## peak region

peaks <- lapply(1:nrow(uni_genes), function(x)
  targ %>%
    filter(paste0("chr", chr) %in% uni_genes$chr[x]) %>%
    filter(gene %in% uni_genes$gene[x]) %>%
    filter(st >= uni_genes$st[x]) %>%
    filter(end <= uni_genes$end[x])
)

aTrack <- lapply(1:length(peaks), function(x) {
  if (nrow(peaks[[x]]) == 0) {
    return(NULL)
  } else {
    aTrack <- anno_track(st = peaks[[x]]$st, 
                         chr = paste0("chr", peaks[[x]]$chr[1]),
                         plot_name = " ",
                         #plot_name = "", # if you want to remove the names 
                         wid = abs(peaks[[x]]$st - peaks[[x]]$end),
                         gen = "hg38",
                         title = "Peaks",
                         feat = "Peak ",
                         len = length(peaks[[x]]$st))
  }
})

########## Gene Axis Track #############
########## This plots the start/end base pair locations
########## for targeted peak region

range_all <- lapply(1:nrow(uni_genes), function(x) {
  range <- uni_genes[x, ] %>%
    dplyr::slice(c(1, n())) %>%
    dplyr::select(st, end)
  range <- c(range$st[1], range$end[2])
})

# build genome axis track for h9
gtrack <- lapply(1:length(range_all), function(x)
  GenomeAxisTrack(range=IRanges(start = range_all[[x]][[1]],
                                end = range_all[[x]][2]
  ),
  fontsize = 14,
  fontcolor = "black",
  fill.range = "black",
  col = "black",
  col.line = "black"
  ))

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull coverage from bam files
bam_files <- list.files("./data/BAM_MeDIP3_j91l/", pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("TEp0", bam_files)]
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

all_cov <- lapply(1:length(range_all), function(x) {
  unlist(lapply(c(1:7), function(y) {
    build_cov(chr = uni_genes$chr[x],
              range1 = range_all[[x]][1],
              range2 = range_all[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})

nam <- bam_nam
nam <- ifelse(grepl("iPSC", nam), "J91i iPSC",
              ifelse(grepl("TEexp", nam), "J91i TEexp",
                     ifelse(grepl("Input", nam), "J91i input", "NA")))
nam[c(2, 4:5, 7)] <- ""

all_dat <- multi_comp_dat_t(all_cov, bam_nam, nam)

########## Highlight Track ############
########## This overlays highlighted tracks over promoter
########## and significant peak regions in the data tracks.
sub <- gene_ensembl[, c(2, 4:5)]
names(sub) <- c("ensembl", "gene_st", "gene_end")
sub_medip <- uni_genes %>%
  dplyr::select(c(gene, chr, st, end, strand, ensembl)) %>%
  dplyr::rename("full_st" = "st") %>%
  dplyr::rename("full_end" = "end") %>%
  mutate(chr = as.numeric(gsub("chr", "", chr))) %>%
  left_join(., sub)
sign_targ <- sign_targ %>%
  mutate(width = end - st,
         chr = as.numeric(chr)) %>%
  full_join(., sub_medip) 
sign_targ <- split(sign_targ, sign_targ$gene)
sign_targ <- lapply(sign_targ, function(x){
  x %>%
    mutate(peak = ifelse(st < full_st, NA, peak),
           st = ifelse(st < full_st, NA, st),
           end = ifelse(st < full_st, NA, end),
           dist_st = ifelse(st < full_st, NA, dist_st),
           position = ifelse(st < full_st, NA, position),
           width = ifelse(st < full_st, NA, width))
})

ht1 <- lapply(1:length(sign_targ), function(x) {
  if (nrow(sign_targ[[x]]) >= 1 & !is.na(sign_targ[[x]]$peak[1])) {
    if (sign_targ[[x]]$strand[1] > 0) {
      st <- c(sign_targ[[x]]$st, sign_targ[[x]]$gene_st[1] -1000)
      width <- c(sign_targ[[x]]$width, 1000)
      col <- c(rep("black", length(st) - 1), "black")
      fill <- c(rep(NA, length(st) - 1), NA)
      lty <- c(rep("solid", length(st) - 1), "dashed")
      HighlightTrack(trackList = all_dat[[x]], 
                     genome = "hg38",
                     start = st, 
                     width = width, 
                     chromosome = sign_targ[[x]]$chr[1],
                     lty = lty,
                     col = col,
                     fill = fill,
                     inBackground = T)
    } else {
      st <- c(sign_targ[[x]]$st, sign_targ[[x]]$gene_end[1])
      width <- c(sign_targ[[x]]$width, 1000)
      col <- c(rep("black", length(st) - 1), "black")
      fill <- c(rep(NA, length(st) - 1), NA)
      lty <- c(rep("solid", length(st) - 1), "dashed")
      HighlightTrack(trackList = all_dat[[x]], 
                     genome = "hg38",
                     start = st, 
                     width = width, 
                     chromosome = sign_targ[[x]]$chr[1],
                     lty = lty,
                     col = col,
                     fill = fill,
                     inBackground = T)
    }
  } else {
    st <- sign_targ[[x]]$gene_end[1]
    width <- 1000
    col <- "black"
    fill <- NA
    lty <- "dashed"
    HighlightTrack(trackList = all_dat[[x]], 
                   genome = "hg38",
                   start = st, 
                   width = width, 
                   chromosome = sign_targ[[x]]$chr[1],
                   lty = lty,
                   col = col,
                   fill = fill,
                   inBackground = T)
  }
})

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- uni_genes$gene

# plot top 10 peaks for ea contrast
plot_data2(itrack, aTrack, ht1, gtrack, bmTrack, title, "targeted_j91")

rm(list = ls())
gc()
