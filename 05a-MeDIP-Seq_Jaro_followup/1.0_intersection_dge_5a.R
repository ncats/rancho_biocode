#' Find in common significant peaks across
#' all samples from deseq binomial comparisons.
#' Plot all samples in locus-zoom plots.

# req'd pkgs
x <- c("Gviz", "GenomicRanges", "rtracklayer", 
       "Rsamtools", "GenomicAlignments",
       "tidyverse", "data.table", "Cairo", 
       "splitstackshape", "biomaRt", "extrafont",
       "extrafontdb")
sapply(x, library, character.only = TRUE)

# source locus zoom fxns
source("./functions/build_locus_zoom_imp.R")

CairoFonts(
  regular="NotoSans-Condensed:style=Regular",
  bold="NotoSans-Condensed:style=Bold",
  italic="NotoSans-Condensed:style=Italic",
  bolditalic="NotoSans-Condensed:style=Bold Italic,BoldItalic",
  symbol="Symbol"
)

###### UNSUPERVISED PLOTS #######
# load deseq objects for all contrasts
all <- list.files("./adj_data/deseq/filtered/", pattern = "deseq_filtered", full.names = T)
all <- lapply(all, function(x) fread(x))
all <- lapply(all, function(x) {
  y <- x %>%
    dplyr::select(MergedRegion, log2FoldChange, padj,
                  chromosome_name_medip, start_position_medip,
                  end_position_medip, gene_active) %>%
    filter(padj < 0.05) %>%
    distinct(.) %>%
    mutate(comb = paste(chromosome_name_medip, start_position_medip, end_position_medip, sep = "_")) %>%
    arrange(padj)
})

# common unique genes across all 3 datasets
region <- unique(all[[3]]$comb)[all[[3]]$comb %in% all[[2]]$comb]
region <- unique(region)[region %in% all[[1]]$comb]

# filter each deseq2 dataset to retain overlapping genes, 
# sort by padj
top <- lapply(1:length(all), function(x) {
  y <- all[[x]]
  y <- y %>%
    filter(comb %in% region) %>%
    group_by(comb) %>%
    arrange(padj, abs(log2FoldChange), comb)
})

# selected top10 overlap by the first time pt
ord10 <- top[[1]]$comb[1:10]
top10 <- lapply(1:length(all), function(x) {
  all[[x]] %>% filter(comb %in% ord10) %>%
    arrange(chromosome_name_medip)
})

# extract MergedRegion, chr, start, end
# cleaned up gene symbol for top 10 loci
top10 <- top10[[1]] %>% 
  dplyr::select(c(1, 4:6, 7)) %>%
  separate_rows(gene_active, sep = "; ") %>%
  distinct(.) %>%
  group_by(MergedRegion, chromosome_name_medip, start_position_medip,
           end_position_medip) %>%
  summarise(gene = paste(gene_active, collapse = "; ")) %>%
  ungroup() %>%
  mutate(locus = c(1:10)) %>%
  separate_rows(gene, sep = "; ")
names(top10)[2:4] <- c("chr", "start", "end")

filt <- unlist(strsplit(top10$gene, "; "))

# pull actual peaks from file
loc <- fread("./adj_data/final_peaks_annotatr.csv") %>%
  data.frame(.) %>%
  separate_rows(gene_active, sep = "; ") %>%
  distinct(.) %>%
  filter(!grepl("NA", gene_active)) %>%
  filter(gene_active %in% filt) %>%
  dplyr::select(c(MergedRegion, chromosome_name_medip, 
                  start_position_medip, end_position_medip,
                  gene_active)) %>%
  distinct(.) %>%
  mutate(locus = top10$locus[match(gene_active, top10$gene)])
names(loc)[2:5] <- c("chr", "start", "end", "gene")
loc <- loc %>%
  mutate(width = (end - start) + 1) %>%
  arrange(gene)

annot <- loc %>%
  group_by(locus, chr) %>%
  summarise(start = min(start),
            end = max(end)) %>%
  data.frame(.)

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:nrow(annot), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste("chr", annot$chr[x], sep = ""))
)

########## Annotation Track ############
# extract MergedReg (medip-seq peaks) from 
# filtered bed file

peaks <- lapply(1:nrow(annot), function(x)
  loc %>%
    filter(chr %in% annot$chr[x]) %>%
    filter(locus %in% annot$locus[x])
)

aTrack <- lapply(1:length(peaks), function(x) 
  aTrack <- anno_track(st = peaks[[x]]$start, 
                       chr = peaks[[x]]$chr[1],
                       plot_name = " ",
                       #plot_name = "", # if you want to remove the names 
                       wid = abs(peaks[[x]]$start - peaks[[x]]$end),
                       gen = "hg38",
                       title = "MeDIP-Seq Peaks",
                       feat = "Peak ",
                       len = length(peaks[[x]]$start))
)

########## Genome Axis Track ############

range_all <- lapply(1:nrow(annot), function(x) {
  range <- annot[x, ] %>%
    dplyr::slice(c(1, n())) %>%
    dplyr::select(start, end)
  range <- c(range$start[1], range$end[2])
})

# build genome axis track for h9
gtrack <- lapply(1:length(range_all), function(x)
  GenomeAxisTrack(range=IRanges(start = range_all[[x]][[1]],
                                end = range_all[[x]][2])
  ))

####### Gene_Region_Track ######

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
genes <- top10$gene
gene_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name'),
                     filters = 'hgnc_symbol',
                     values = genes,
                     mart = ensembl) %>%
 data_frame(.) %>%
 filter(!grepl("_", chromosome_name)) 
df <- data.frame(gene_ensembl)
names(df) <- c("gene", "ensembl", "chr")
add <- top10 %>% filter(grepl("LOC", gene)) %>%
  dplyr::select(c(gene, chr)) %>%
  mutate(ensembl = "NA")
df <- rbind(df, add)
df$locus <- top10$locus[match(df$gene, top10$gene)]

annot <- annot %>%
  full_join(., df)
annot <- annot %>%
  group_by(locus, chr, start, end) %>%
  summarise(gene = paste(gene, collapse ="; "),
         ensembl = paste(ensembl, collapse = "; "))

track <- top10 %>% dplyr::select(c(MergedRegion, chr, start, end, locus)) %>%
  distinct(.) %>%
  mutate(gene = annot$gene)

bmTrack <- lapply(1:nrow(track), function(x)
  BiomartGeneRegionTrack(genome = "hg38",
                         chromosome = track$chr[x], 
                         #start = track$start[x], 
                         #end = track$end[x],
                         name = annot$gene[x], #to shrink the track
                         stacking = "squish",
                         collapseTranscripts = "meta",
                         filters = list("ensembl_gene_id" = unlist(strsplit(annot$ensembl[x], "; "))),
                         #using fill color from UCSC genome browser
                         fill = "#0C0C78",
                         col = "#0C0C78"
  )
)

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull coverage from bam files
bam_files <- list.files("./BAM_MeDIP3/sort", pattern = ".bam$", full.names = T)
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

all_cov <- lapply(1:length(range_all), function(x) {
  unlist(lapply(c(1:9), function(y) {
    build_cov(chr = paste("chr", peaks[[x]]$chr[1], sep = ""),
              range1 = range_all[[x]][1],
              range2 = range_all[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})
all_dat <- multi_comp_dat(all_cov, bam_nam)

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- lapply(1:nrow(annot), function(x) {
  y <- unlist(strsplit(annot$gene[x], "; "))
  y <- unique(y)
  y <- y[!grepl("NA", y)]
  
  if (length(y) == 0) {
    b <- "Unannotated_locus"
  } else if (length(y) >= 2) {
    z <- paste(y, collapse = "_")
    b <- paste("Annotated_as_", z, "_loci", sep = "")
  } else {
    b <- y
  }
})

title <- lapply(1:length(title), function(x) {
  z <- unlist(title[[x]])
  if (sum(z == "Unannotated_locus") > 1) {
    count <- sum(z == "Unannotated_locus")
    count <- seq(1, count)
    count <- paste("Unannotated_locus_", count, sep = "")
    idx <- which(grepl("Unannot", z))
    z[idx] <- count
  }
  return(z)
})

# plot top 10 peaks for ea contrast
plot_data2(itrack, aTrack, all_dat, gtrack, bmTrack, top10, "TE_all_comps", title)

####### SUPERVISED LOCUS ZOOM PLOTS ########
targ_genes <- c("POU5F1", "NANOG", "SOX2", "NKX1-2",
                "ZSCAN10", "WNT6", "WNT10A", "TFAP2A",
                "TFAP2B", "TFAP2C", "TP63", "NR2F2",
                "ELF5")

loc <- fread("./adj_data/final_peaks_annotatr.csv") %>%
  separate_rows(gene_active, sep = "; ") %>%
  distinct(.) %>%
  #filter(!grepl("NA", gene_active)) %>%
  filter(gene_active %in% targ_genes) %>%
  dplyr::select(c(MergedRegion, chromosome_name_medip, 
                  start_position_medip, end_position_medip,
                  gene_active)) %>%
  distinct(.) 
names(loc)[2:5] <- c("chr", "start", "end", "gene")
loc <- loc %>%
  mutate(width = (end - start) + 1) %>%
  arrange(gene)

annot <- loc %>%
  group_by(gene, chr) %>%
  summarise(start = min(start),
            end = max(end)) %>%
  data.frame(.)

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:nrow(annot), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste("chr", annot$chr[x], sep = ""))
)

########## Annotation Track ############
# extract MergedReg (medip-seq peaks) from 
# filtered bed file

peaks <- lapply(1:nrow(annot), function(x)
  loc %>%
    filter(chr %in% annot$chr[x]) %>%
    filter(gene %in% annot$gene[x])
)

nam <- lapply(1:length(peaks), function(x) {
  y <- peaks[[x]] %>%
    pull(gene)
  y <- unlist(sapply(1:length(y), function(z)
    paste(y[[z]], collapse = "\n")
  ))
  return(y)
})

aTrack <- lapply(1:length(peaks), function(x) 
  aTrack <- anno_track(st = peaks[[x]]$start, 
                       chr = peaks[[x]]$chr[1],
                       plot_name = nam[[x]],
                       #plot_name = "", # if you want to remove the names 
                       wid = abs(peaks[[x]]$start - peaks[[x]]$end),
                       gen = "hg38",
                       title = "MeDIP-Seq Peaks",
                       feat = "Peak ",
                       len = length(peaks[[x]]$start))
)

########## Genome Axis Track ############

range_all <- lapply(1:nrow(annot), function(x) {
  range <- annot[x, ] %>%
    dplyr::slice(c(1, n())) %>%
    dplyr::select(start, end)
  range <- c(range$start[1], range$end[2])
})

# build genome axis track for h9
gtrack <- lapply(1:length(range_all), function(x)
  GenomeAxisTrack(range=IRanges(start = range_all[[x]][[1]],
                                end = range_all[[x]][2])
  ))

####### Gene_Region_Track ######

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
genes <- unlist(strsplit(annot$gene, "; "))
gene_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name'), 
                      filters = 'hgnc_symbol', 
                      values = genes, 
                      mart = ensembl) %>%
  data_frame(.) %>%
  filter(!grepl("_", chromosome_name)) %>%
  dplyr::select(ensembl_gene_id)
df <- data.frame("gene" = genes,
                 "ensembl" = gene_ensembl)
annot <- annot %>%
  left_join(., df) %>%
  group_by(chr, start, end)

bmTrack <- lapply(1:nrow(annot), function(x)
  BiomartGeneRegionTrack(genome = "hg38",
                         chromosome = annot$chr[x], 
                         start = annot$start[x], 
                         end = annot$end[x],
                         name = "", #to shrink the track
                         collapseTranscripts = "meta",
                         filters = list("ensembl_gene_id" = annot$ensembl_gene_id[x]),
                         #using fill color from UCSC genome browser
                         fill = "#0C0C78",
                         col = "#0C0C78"
  )
)

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

all_cov <- lapply(1:length(range_all), function(x) {
    unlist(lapply(c(1:9), function(y) {
      build_cov(chr = paste("chr", peaks[[x]]$chr[1], sep = ""),
                range1 = range_all[[x]][1],
                range2 = range_all[[x]][2],
                file = bam_files[y])
    }), recursive = F)
  })

all_dat <- multi_comp_dat(all_cov, bam_nam)


########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- lapply(1:nrow(annot), function(x) {
  z <- annot$gene[x]
})

# plot top 10 peaks for ea contrast
plot_data2(itrack, aTrack, all_dat, gtrack, bmTrack, annot, "Supervised", title)



rm(list = ls())
gc()
