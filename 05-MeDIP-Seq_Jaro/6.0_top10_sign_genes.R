#' Build locus zoom plots using Gviz for top 10
#' significant peaks

# req'd pkgs
x <- c("Gviz", "GenomicRanges", "rtracklayer", 
       "Rsamtools", "GenomicAlignments",
       "dplyr", "data.table", "Cairo", 
       "splitstackshape")
sapply(x, library, character.only = TRUE)

# source locus zoom fxns
source("./functions/build_locus_zoom_imp.R")

# load deseq objects for all contrasts
all <- list.files("./adj_data/deseq/filtered/", pattern = "deseq_filtered", full.names = T)
all <- lapply(all, function(x) fread(x))
all <- lapply(all, function(x) {
  y <- x %>%
    arrange(padj) %>%
    dplyr::select(MergedRegion, log2FoldChange, padj,
                  chromosome_name_medip, start_position_medip,
                  end_position_medip, strand, width, id,
                  gene_anotatr, gene_active)
  y <- y[c(1:10),]
  })

# name that contrast
nam <- gsub("_de.*", "", list.files("./adj_data/deseq/", pattern = "deseq2.csv"))

# now read in filtered merged peaks
medip_peaks <- import.bed("./adj_data/MergedPeaks.bed")

########## Annotation Track ############
########## Utilizes .bed file from Active Motif
########## plot specified medipseq peaks using an annotation track ##########

all_peaks <- data.frame(medip_peaks)
all_peaks$name <- gsub("MergedReg_", "Peak", all_peaks$name)

# build annotation track for all comps
aTrack <- lapply(1:length(all), function(x) {
  t <- all[[x]]
  nam <- paste("Peak", t$MergedRegion, sep = "")
  # nam <- all_peaks %>%
  #   filter(name %in% paste("Peak", t$MergedRegion, sep = "")) %>%
  #   pull(name)
  width <- as.numeric(unlist(sapply(1:length(t$width), function(x) {
    y <- unlist(strsplit(t$width[x], "; "))[1]
  })))
  aTrack <- lapply(1:nrow(t), function(y) 
    anno_track(st = t[y]$start_position_medip,
               wid = width[y], 
               chr = t[y]$chromosome_name_medip[1],
               plot_name = nam[y],
               gen = "hg38",
               title = "MeDIP-Seq Peaks",
               feat = "Peak ",
               len = length(t[y]$start_position_medip))
  )
})

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

itrack <- lapply(1:length(all), function(x) 
  unlist(lapply(1:nrow(all[[x]]), function(y)
    ideo_track(gen = 'hg38',
               chr = paste("chr", all[[x]]$chromosome_name_medip[y], sep = ""))
    ), recursive = F)
  
)

########## Genome Axis Track ############

range_all <- lapply(1:length(all), function(x) {
  t <- all[[x]]
  lapply(1:nrow(t), function(y) {
    range <- t[y, ] %>%
      dplyr::slice(c(1, n())) %>%
      dplyr::select(start_position_medip, end_position_medip)
    range <- c(range$start[1], range$end[2])
  })
})

# build genome axis track for h9
gtrack <- lapply(1:length(range_all), function(x)
  lapply(1:length(range_all[[x]]), function(y)
    GenomeAxisTrack(range=IRanges(start = range_all[[x]][[y]][[1]],
                                  end = range_all[[x]][[y]][2])
  ))
)

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull coverage from bam files
bam_files <- list.files("./BAM_MeDIP3//", pattern = ".bam$", full.names = T)
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

# split out into contrasts:
ipsc_tep0 <- bam_nam[grepl("iPSC|TEp0", bam_nam)]
ipsc_tep0_files <- bam_files[grepl("iPSC|TEp0", bam_nam)]
# d0 vs d30 coverage/data track
ipsc_tep0_cov <- multi_comp_cov(all[[1]], ipsc_tep0, range_all[[1]], ipsc_tep0_files)
ipsc_tep0_dat <- multi_comp_dat(ipsc_tep0_cov, ipsc_tep0)

ipsc_teexp <- bam_nam[grepl("iPSC|TEexp", bam_nam)]
ipsc_teexp_files <- bam_files[grepl("iPSC|TEexp", bam_nam)]
# d0 vs d50 coverage/data track
ipsc_teexp_cov <- multi_comp_cov(all[[2]], ipsc_teexp, range_all[[2]], ipsc_teexp_files)
ipsc_teexp_dat <- multi_comp_dat(ipsc_teexp_cov, ipsc_teexp)

tep0_teexp <- bam_nam[grepl("TEp0|TEexp", bam_nam)]
tep0_teexp_files <- bam_files[grepl("TEp0|TEexp", bam_nam)]
# a1 vs d0 coverage/data track
tep0_teexp_cov <- multi_comp_cov(all[[3]], tep0_teexp, range_all[[3]], tep0_teexp_files)
tep0_teexp_dat <- multi_comp_dat(tep0_teexp_cov, tep0_teexp)

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- lapply(1:length(all), function(x) {
  y <- strsplit(all[[x]]$gene_anotatr, "; ")
  y <- lapply(y, function(t) {
    z <- unique(t)
    z <- z[!grepl("NA", z)]
    })
  
  z <- strsplit(all[[x]]$gene_active, "; ")
  z <- lapply(z, function(t) {
    u <- unique(t)
    u <- u[!grepl("NA", u)]
  })
  
  y <- lapply(1:length(y), function(x) 
    unique(c(y[[x]], z[[x]]))
    )
  
  y <- lapply(y, function(t) {
    if (length(t) == 0) {
      b <- "Unannotated_locus"
    } else if (length(t) >= 2) {
      z <- paste(t, collapse = "_")
      b <- paste("Annotated_as_", z, "_loci", sep = "")
    } else {
      b <- t
    }
  })
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
plot_data(itrack[[1]], aTrack[[1]], ipsc_tep0_dat, gtrack[[1]], all[[1]], "J98i-iPSC_vs_J98i-TEp0", title[[1]], "top10")
plot_data(itrack[[2]], aTrack[[2]], ipsc_teexp_dat, gtrack[[2]], all[[2]], "J98i-iPSC_vs_TEexp", title[[2]], "top10")
plot_data(itrack[[3]], aTrack[[3]], tep0_teexp_dat, gtrack[[3]], all[[3]], "J98i-TEp0_vs_TEexp", title[[3]], "top10")


rm(list = ls())
gc()
