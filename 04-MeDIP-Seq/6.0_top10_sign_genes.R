#' Build locus zoom plots using Gviz for top 10
#' significant peaks

# req'd pkgs
x <- c("Gviz", "GenomicRanges", "rtracklayer", 
       "Rsamtools", "GenomicAlignments",
       "dplyr", "data.table", "Cairo", 
       "splitstackshape", "tidyr", "biomaRt",
       "ggplot2")
sapply(x, library, character.only = TRUE)

CairoFonts(
  regular="NotoSans-Condensed:style=Regular",
  bold="NotoSans-Condensed:style=Bold",
  italic="NotoSans-Condensed:style=Italic",
  bolditalic="NotoSans-Condensed:style=Bold Italic,BoldItalic",
  symbol="Symbol"
)

# source locus zoom fxns
source("./functions/build_locus_zoom.R")

# load deseq objects for all contrasts
all <- list.files("./adj_data/deseq/filtered/", pattern = "deseq_filtered", full.names = T)
all <- lapply(all, function(x) fread(x))
all <- lapply(all, function(x) {
  y <- x %>%
    arrange(padj) %>%
    dplyr::select(MergedRegion, log2FoldChange, padj,
                  chromosome_name_medip, start_position_medip,
                  end_position_medip, strand, width, id,
                  gene_anotatr
                  #symbol
                  )
  y <- y[c(1:10), ]
  })

# extract MergedRegion, chr, start, end
# cleaned up gene symbol for top 10 loci
all <- lapply(1:length(all), function(x) {
  y <- all[[x]]
  y <- y %>% 
    dplyr::select(c(1, 4:6, 10)) %>%
    separate_rows(gene_anotatr, sep = "; ") %>%
    distinct(.) %>%
    group_by(MergedRegion, chromosome_name_medip, start_position_medip,
             end_position_medip) %>%
    summarise(gene = paste(gene_anotatr, collapse = "; ")) %>%
    ungroup() %>%
    mutate(locus = c(1:10)) %>%
    separate_rows(gene, sep = "; ")
  names(y)[2:4] <- c("chr", "start", "end")
  return(y)
})

filt <- lapply(all, function(x) {
  y <- unique(unlist(strsplit(x$gene, "; ")))
  y <- y[y != "NA"]
})
peak <- fread("./adj_data/final_peaks_annotatr.csv")

# pull actual peaks from file
loc <- lapply(1:length(all), function(x) {
  y <- peak %>%
    data.frame(.) %>%
    separate_rows(gene_anotatr, sep = "; ") %>%
    separate_rows(gene_anotatr, sep = ", ") %>%
    distinct(.) %>%
    filter(!grepl("NA", gene_anotatr)) %>%
    filter(gene_anotatr %in% filt[[x]]) %>%
    dplyr::select(c(MergedRegion, chromosome_name_medip, 
                    start_position_medip, end_position_medip,
                    gene_anotatr)) %>%
    distinct(.) %>%
    mutate(locus = all[[x]]$locus[match(gene_anotatr, all[[x]]$gene)])
  names(y)[2:5] <- c("chr", "start", "end", "gene")
  y <- y %>%
    mutate(width = (end - start) + 1) %>%
    arrange(gene)
  
  na <- all[[x]] %>%
    filter(gene == "NA") %>%
    mutate(width = abs(end - start))
  
  y <- rbind(y, na) %>% arrange(locus)
  
  return(y)
})

annot <- lapply(1:length(loc), function(x) 
  loc[[x]] %>%
  group_by(locus, chr) %>%
  summarise(start = min(start),
            end = max(end)) %>%
  data.frame(.))

# name that contrast
nam <- gsub("_de.*", "", list.files("./adj_data/deseq/", pattern = "deseq2.csv"))
nam <- gsub("hPSCs_day0", "D0", nam)
nam <- gsub("NCRM5-", "", nam)

# now read in filtered merged peaks
#medip_peaks <- import.bed("./adj_data/MergedPeaks.bed")

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:length(annot), function(x) 
  lapply(1:nrow(annot[[x]]), function(y)
    ideo_track(gen = 'hg38',
               chr = paste("chr", annot[[x]]$chr[y], sep = ""))
  )
)

########## Annotation Track ############
########## Utilizes .bed file from Active Motif
########## plot specified medipseq peaks using an annotation track ##########

#all_peaks <- data.frame(medip_peaks)
#all_peaks$name <- gsub("MergedReg_", "Peak", all_peaks$name)

# peaks <- lapply(1:length(annot), function(x) {
#   y <- lapply(1:nrow(annot[[x]]), function(y)
#     loc[[x]] %>%
#       filter(chr %in% annot[[x]]$chr[y]) %>%
#       filter(locus %in% annot[[x]]$locus[y])
#   )
#   #y <- unlist(y, recursive = F)
# })
# 
# aTrack <- lapply(1:length(peaks), function(x) 
#   lapply(1:length(peaks[[x]]), function(y) {
#     aTrack <- anno_track(st = peaks[[x]][[y]]$start, 
#                          chr = peaks[[x]][[y]]$chr[1],
#                          plot_name = " ", # if you want to remove the names 
#                          wid = abs(peaks[[x]][[y]]$start - peaks[[x]][[y]]$end),
#                          gen = "hg38",
#                          title = "MeDIP-Seq Peaks",
#                          feat = "Peak ",
#                          len = length(peaks[[x]][[y]]$start)
#                          )
#   })
# )

aTrack <- lapply(1:length(loc), function(x) {
  y <- loc[[x]]
  y <- split(y, y$locus)
  lapply(1:length(y), function(t) {
    aTrack <- anno_track(st = y[[t]]$start, 
                         chr = y[[t]]$chr[1],
                         plot_name = " ", # if you want to remove the names 
                         wid = abs(y[[t]]$start - y[[t]]$end),
                         gen = "hg38",
                         title = "MeDIP-Seq Peaks",
                         feat = "Peak ",
                         len = length(y[[t]]$start)
    )
  })
})

########## Genome Axis Track ############

range_all <- lapply(1:length(annot), function(x) {
  lapply(1:nrow(annot[[x]]), function(y) {
    range <- annot[[x]][y, ] %>%
      dplyr::slice(c(1, n())) %>%
      dplyr::select(start, end)
    range <- c(range$start[1], range$end[2])
  })
})

# build genome axis track for h9
gtrack <- lapply(1:length(range_all), function(x) {
  lapply(1:length(range_all[[x]]), function(y) {
    GenomeAxisTrack(range=IRanges(start = range_all[[x]][[y]][[1]],
                                  end = range_all[[x]][[y]][2]))
  })
})

####### Gene_Region_Track ######

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
genes <- lapply(all, function(x) x$gene)
gene_ensembl <- lapply(1:length(genes), function(x) 
  getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name'),
                      filters = 'hgnc_symbol',
                      values = genes[[x]],
                      mart = ensembl) %>%
  data_frame(.) %>%
  filter(!grepl("_", chromosome_name)))
df <- lapply(gene_ensembl, function(x) {
  y <- data.frame(x)
  names(y) <- c("gene", "ensembl", "chr")
  return(y)
  })

add <- lapply(all, function(x) 
  x %>% 
    filter(grepl("LOC|NA", gene)) %>%
    dplyr::select(c(gene, chr)) %>%
    mutate(ensembl = "NA"))
df <- lapply(1:length(df), function(x) rbind(df[[x]], add[[x]]))

locus <- lapply(1:length(df), function(x) {
  loc <- all[[x]]$locus[match(df[[x]]$gene, all[[x]]$gene)]
  #df[[x]]$locus <- loc
  #return(df)
  })
df <- lapply(1:length(df), function(x) df[[x]] %>% mutate(locus = locus[[x]]))

annot <- lapply(1:length(annot), function(x)
  annot[[x]] %>%
  full_join(., df[[x]]) %>%
  group_by(locus, chr, start, end) %>%
  summarise(gene = paste(gene, collapse ="; "),
            ensembl = paste(ensembl, collapse = "; "))
  )
annot[[5]] <- annot[[5]] %>% filter(!is.na(start))

bmTrack <- lapply(1:length(annot), function(x) {
  lapply(1:nrow(annot[[x]]), function(y) {
    BiomartGeneRegionTrack(genome = "hg38",
                           chromosome = annot[[x]]$chr[y], 
                           #start = track$start[x], 
                           #end = track$end[x],
                           name = annot[[x]]$gene[y], #to shrink the track
                           stacking = "squish",
                           collapseTranscripts = "meta",
                           filters = list("ensembl_gene_id" = unlist(strsplit(annot[[x]]$ensembl[y], "; "))),
                           #using fill color from UCSC genome browser
                           fill = "#0C0C78",
                           col = "#0C0C78"
    )
  })
})

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull coverage from bam files
bam_files <- list.files("./BAM", pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("Test", bam_files)]
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)
bam_nam <- bam_nam[!grepl("Test", bam_nam)]

# split out into contrasts:
d0_d30 <- bam_nam[grepl("D0|D30", bam_nam)]
d0_d30_files <- bam_files[grepl("D0|D30", bam_nam)]
# d0 vs d30 coverage/data track
d0_d30_cov <- multi_comp_cov(annot[[6]], d0_d30, range_all[[6]], d0_d30_files)
d0_d30_dat <- multi_comp_dat(d0_d30_cov, d0_d30)

d0_d50 <- bam_nam[grepl("D0|D50", bam_nam)]
d0_d50_files <- bam_files[grepl("D0|D50", bam_nam)]
# d0 vs d50 coverage/data track
d0_d50_cov <- multi_comp_cov(annot[[7]], d0_d50, range_all[[7]], d0_d50_files)
d0_d50_dat <- multi_comp_dat(d0_d50_cov, d0_d50)

d0_a1 <- bam_nam[grepl("D0|A1", bam_nam)]
d0_a1_files <- bam_files[grepl("D0|A1", bam_nam)]
# a1 vs d0 coverage/data track
d0_a1_cov <- multi_comp_cov(annot[[4]], d0_a1, range_all[[4]], d0_a1_files)
d0_a1_dat <- multi_comp_dat(d0_a1_cov, d0_a1)

a1_lsb <- bam_nam[grepl("A1|LSB", bam_nam)]
a1_lsb_files <- bam_files[grepl("A1|LSB", bam_nam)]
# a1 vs lsb coverage/data track
a1_lsb_cov <- multi_comp_cov(annot[[1]], a1_lsb, range_all[[1]], a1_lsb_files)
a1_lsb_dat <- multi_comp_dat(a1_lsb_cov, a1_lsb)

a1_d30 <- bam_nam[grepl("A1|D30", bam_nam)]
a1_d30_files <- bam_files[grepl("A1|D30", bam_nam)]
# a1 vs d30 coverage/data track
a1_d30_cov <- multi_comp_cov(annot[[2]], a1_d30, range_all[[2]], a1_d30_files)
a1_d30_dat <- multi_comp_dat(a1_d30_cov, a1_d30)

a1_d50 <- bam_nam[grepl("A1|D50", bam_nam)]
a1_d50_files <- bam_files[grepl("A1|D50", bam_nam)]
# a1 vs d50 coverage/data track
a1_d50_cov <- multi_comp_cov(annot[[3]], a1_d50, range_all[[3]], a1_d50_files)
a1_d50_dat <- multi_comp_dat(a1_d50_cov, a1_d50)

d0_lsb <- bam_nam[grepl("D0|LSB", bam_nam)]
d0_lsb_files <- bam_files[grepl("D0|LSB", bam_nam)]
# lsb vs d0 coverage/data track
d0_lsb_cov <- multi_comp_cov(annot[[5]], d0_lsb, range_all[[5]], d0_lsb_files)
d0_lsb_dat <- multi_comp_dat(d0_lsb_cov, d0_lsb)

d30_d50 <- bam_nam[grepl("D30|D50", bam_nam)]
d30_d50_files <- bam_files[grepl("D30|D50", bam_nam)]
# d30 vs d50 coverage/data track
d30_d50_cov <- multi_comp_cov(annot[[8]], d30_d50, range_all[[8]], d30_d50_files)
d30_d50_dat <- multi_comp_dat(d30_d50_cov, d30_d50)

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- lapply(1:length(all), function(x) {
  y <- all[[x]] %>%
    group_by(locus) %>%
    summarise(gene = paste(gene, collapse = "; ")) %>%
    pull(gene)
  
  y <- sapply(1:length(y), function(t) {
    if (y[[t]] == "NA") {
      b <- "Unannotated locus"
    } else {
      b <- y[[t]]
    }
  })
  #y <- unique(unlist(y))
})

title <- lapply(1:length(title), function(x) {
  z <- unlist(title[[x]])
  if (sum(z == "Unannotated locus") > 1) {
    count <- sum(z == "Unannotated locus")
    count <- seq(1, count)
    count <- paste("Unannotated locus ", count, sep = "")
    idx <- which(grepl("Unannot", z))
    z[idx] <- count
  }
  z <- unique(z)
  return(z)
})

# plot top 10 peaks for ea contrast
plot_data(itrack[[1]], aTrack[[1]], a1_lsb_dat, gtrack[[1]], annot[[1]], bmTrack[[1]], "A1_vs_LSB", title[[1]])
plot_data(itrack[[2]], aTrack[[2]], a1_d30_dat, gtrack[[2]], annot[[2]], bmTrack[[2]], "A1_vs_D30", title[[2]])
plot_data(itrack[[3]], aTrack[[3]], a1_d50_dat, gtrack[[3]], annot[[3]], bmTrack[[3]], "A1_vs_D50", title[[3]])
plot_data(itrack[[4]], aTrack[[4]], d0_a1_dat, gtrack[[4]], annot[[4]], bmTrack[[4]], "D0_vs_A1", title[[4]])
plot_data(itrack[[5]], aTrack[[5]], d0_lsb_dat, gtrack[[5]], annot[[5]], bmTrack[[5]], "D0_vs_LSB", title[[5]])
plot_data(itrack[[6]], aTrack[[6]], d0_d30_dat, gtrack[[6]], annot[[6]], bmTrack[[6]], "D0_vs_D30", title[[6]])
plot_data(itrack[[7]], aTrack[[7]], d0_d50_dat, gtrack[[7]], annot[[7]], bmTrack[[7]], "D0_vs_D50", title[[7]])
plot_data(itrack[[8]], aTrack[[8]], d30_d50_dat, gtrack[[8]], annot[[8]], bmTrack[[8]], "D30_vs_D50", title[[8]])


rm(list = ls())
gc()
