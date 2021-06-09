#' Build locus zoom plots using Gviz for top 10
#' significant peaks

# req'd pkgs
x <- c("Gviz", "GenomicRanges", "rtracklayer", 
       "Rsamtools", "GenomicAlignments",
       "dplyr", "data.table", "Cairo", "extrafontdb",
       "extrafont", "splitstackshape", "openxlsx", "biomaRt",
       "scales")
sapply(x, library, character.only = TRUE)
loadfonts()

source("./functions/upd_build_locus_zoom.R")

# set cairo text backend for cairo_pdf
CairoFonts(
  regular="NotoSans-Condensed:style=Regular",
  bold="NotoSans-Condensed:style=Bold",
  italic="NotoSans-Condensed:style=Italic",
  bolditalic="NotoSans-Condensed:style=Bold Italic, BoldItalic",
  symbol="Symbol"
)

###### UPDATED LOCUS-ZOOM PLOTS #####
top10 <- fread("./adj_data/deseq/filtered/J91i-iPSC_TEexp_deseq_filtered.csv") %>%
  arrange(padj, abs(log2FoldChange))
top10 <- top10 %>% 
  dplyr::rename("peak" = "MergedRegion") %>%
  tidyr::separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
  filter(Position == "in gene") %>%
  dplyr::select(c(1,8:10,16:18))
names(top10)[2:7] <- c("chr", "st", "end", "gene", "dist_st", "position")
top10 <- top10[c(1:10), ]

loc <- readWorkbook("./data/4441NIH_J91i-MeDIP_mergedregs.xlsx") %>%
  data.frame(.) %>%
  dplyr::select(c(Merged.Region, Chromosome, Start, End,
                  Gene.List, Dist.to.Start, Position)) %>%
  tidyr::separate_rows(Gene.List, Dist.to.Start, Position, sep = ", ") %>%
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

####### Gene_Region_Track ######

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
gene_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 
                                   'start_position', 'end_position', "strand"),
                      filters = 'hgnc_symbol',
                      values = top10$gene,
                      mart = ensembl) %>%
  tibble(.) %>%
  filter(!grepl("_", chromosome_name)) 
df <- gene_ensembl %>%
  mutate(start_position = start_position - 1000,
         end_position = end_position + 1000)
names(df) <- c("gene", "ensembl", "chr", "st", "end", "strand")
df<- df %>%
  mutate(chr = paste0("chr", chr))

top10 <- top10 %>%
  dplyr::select(-c(chr, st, end)) %>%
  full_join(., df)

bmTrack <- lapply(1:nrow(top10), function(x)
  BiomartGeneRegionTrack(genome = "hg38",
                         chromosome = top10$chr[x],
                         #name = top10$gene[x], #to shrink the track
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
                         min.height = 10
  )
)

########## Annotation Track ############

peaks <- lapply(1:nrow(top10), function(x)
  loc %>%
    filter(position == "in gene") %>%
    filter(paste0("chr", chr) %in% top10$chr[x]) %>%
    filter(gene %in% top10$gene[x]) %>%
    filter(start >= top10$st[x]) %>%
    filter(end <= top10$end[x])
)

aTrack <- lapply(1:length(peaks), function(x) 
  aTrack <- anno_track(st = peaks[[x]]$start, 
                       chr = paste0("chr", peaks[[x]]$chr[1]),
                       plot_name = " ",
                       #plot_name = "", # if you want to remove the names 
                       wid = abs(peaks[[x]]$start - peaks[[x]]$end),
                       gen = "hg38",
                       title = "Peaks",
                       feat = "Peak ",
                       len = length(peaks[[x]]$start))
)

########## Genome Axis Track ############

range_all <- lapply(1:nrow(top10), function(x) {
  range <- top10[x, ] %>%
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
bam_files <- bam_files[!grepl('J98i', bam_files)]
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

all_cov <- lapply(1:length(range_all), function(x) {
  unlist(lapply(c(1:7), function(y) {
    build_cov(chr = top10$chr[x],
              range1 = range_all[[x]][1],
              range2 = range_all[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})

nam <- bam_nam
nam <- ifelse(grepl("Input", nam), "J91i input",
              ifelse(grepl("iPSC", nam), "J91i iPSC",
              ifelse(grepl("TEexp", nam), "J91i TEexp", "NA")))
nam[c(2, 4:5, 7 )] <- ""

all_dat <- multi_comp_dat_t(all_cov, bam_nam, nam)

########## highlight peak var & ######
########## promoter region      ######
sub_medip <- top10 %>%
  dplyr::select(c(gene, chr, st, end, strand)) %>%
  dplyr::rename("full_st" = "st") %>%
  dplyr::rename("full_end" = "end") %>%
  mutate(chr = as.numeric(gsub("chr", "", chr)))
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
plot_data2(itrack, aTrack, ht1, gtrack, bmTrack, title, bam_nam)
rm(list = ls())
gc()
