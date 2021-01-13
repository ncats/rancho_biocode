#' Build locus zoom plots using Gviz for top 10
#' significant peaks

# req'd pkgs
x <- c("Gviz", "GenomicRanges", "rtracklayer", 
       "Rsamtools", "GenomicAlignments",
       "dplyr", "data.table", "Cairo", 
       "splitstackshape")
sapply(x, library, character.only = TRUE)

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
                  symbol)
  y <- y[c(1:10), ]
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
  nam <- all_peaks %>%
    filter(name %in% paste("Peak", t$MergedRegion, sep = "")) %>%
    pull(name)
  width <- as.numeric(unlist(sapply(1:length(t$width), function(x) {
    y <- unlist(strsplit(t$width[x], "; "))[1]
  })))
  aTrack <- lapply(1:length(t), function(y) 
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
bam_files <- list.files("./BAM/", pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("Test", bam_files)]
# re-name bam files
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)
bam_nam <- bam_nam[!grepl("Test", bam_nam)]

# split out into contrasts:
d0_d30 <- bam_nam[grepl("NCRM5-[0-9]{1}|D30", bam_nam)]
d0_d30_files <- bam_files[grepl("NCRM5-[0-9]{1}|D30", bam_nam)]
d# d0 vs d30 coverage/data track
d0_d30_cov <- multi_comp_cov(all[[6]], d0_d30, range_all[[6]], d0_d30_files)
d0_d30_dat <- multi_comp_dat(d0_d30_cov, d0_d30)

d0_d50 <- bam_nam[grepl("NCRM5-[0-9]{1}|D50", bam_nam)]
d0_d50_files <- bam_files[grepl("NCRM5-[0-9]{1}|D50", bam_nam)]
d# d0 vs d50 coverage/data track
d0_d50_cov <- multi_comp_cov(all[[7]], d0_d50, range_all[[7]], d0_d50_files)
d0_d50_dat <- multi_comp_dat(d0_d50_cov, d0_d50)

d0_a1 <- bam_nam[grepl("NCRM5-[0-9]{1}|A1", bam_nam)]
d0_a1_files <- bam_files[grepl("NCRM5-[0-9]{1}|A1", bam_nam)]
# a1 vs d0 coverage/data track
d0_a1_cov <- multi_comp_cov(all[[4]], d0_a1, range_all[[4]], d0_a1_files)
d0_a1_dat <- multi_comp_dat(d0_a1_cov, d0_a1)

a1_lsb <- bam_nam[grepl("A1|LSB", bam_nam)]
a1_lsb_files <- bam_files[grepl("A1|LSB", bam_nam)]
# a1 vs lsb coverage/data track
a1_lsb_cov <- multi_comp_cov(all[[1]], a1_lsb, range_all[[1]], a1_lsb_files)
a1_lsb_dat <- multi_comp_dat(a1_lsb_cov, a1_lsb)

a1_d30 <- bam_nam[grepl("A1|D30", bam_nam)]
a1_d30_files <- bam_files[grepl("A1|D30", bam_nam)]
# a1 vs d30 coverage/data track
a1_d30_cov <- multi_comp_cov(all[[2]], a1_d30, range_all[[2]], a1_d30_files)
a1_d30_dat <- multi_comp_dat(a1_d30_cov, a1_d30)

a1_d50 <- bam_nam[grepl("A1|D50", bam_nam)]
a1_d50_files <- bam_files[grepl("A1|D50", bam_nam)]
# a1 vs d50 coverage/data track
a1_d50_cov <- multi_comp_cov(all[[3]], a1_d50, range_all[[3]], a1_d50_files)
a1_d50_dat <- multi_comp_dat(a1_d50_cov, a1_d50)

d0_lsb <- bam_nam[grepl("NCRM5-[0-9]{1}|LSB", bam_nam)]
d0_lsb_files <- bam_files[grepl("NCRM5-[0-9]{1}|LSB", bam_nam)]
# lsb vs d0 coverage/data track
d0_lsb_cov <- multi_comp_cov(all[[5]], d0_lsb, range_all[[5]], d0_lsb_files)
d0_lsb_dat <- multi_comp_dat(d0_lsb_cov, d0_lsb)

d30_d50 <- bam_nam[grepl("D30|D50", bam_nam)]
d30_d50_files <- bam_files[grepl("D30|D50", bam_nam)]
# d30 vs d50 coverage/data track
d30_d50_cov <- multi_comp_cov(all[[8]], d30_d50, range_all[[8]], d30_d50_files)
d30_d50_dat <- multi_comp_dat(d30_d50_cov, d30_d50)

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# select one unique title name per locus
title <- lapply(1:length(all), function(x) {
  y <- strsplit(all[[x]]$symbol, "; ")
  y <- lapply(y, function(t) {
    z <- unique(t)
    z <- z[!grepl("NA", z)]
    })
  
  y <- lapply(y, function(t) {
    if (length(t) == 0) {
      b <- "Unannotated locus"
    } else if (length(t) >= 2) {
      z <- paste(t, collapse = ", ")
      b <- paste("Annotated as: ", z, " loci", sep = "")
    } else {
      b <- t
    }
  })
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
  return(z)
})

# plot top 10 peaks for ea contrast
plot_data(itrack[[1]], aTrack[[1]], a1_lsb_dat, gtrack[[1]], all[[1]], "A1 vs LSB at ", title[[1]])
plot_data(itrack[[2]], aTrack[[2]], a1_d30_dat, gtrack[[2]], all[[2]], "A1 vs Day 30 at ", title[[2]])
plot_data(itrack[[3]], aTrack[[3]], a1_d50_dat, gtrack[[3]], all[[3]], "A1 vs Day 50 at ", title[[3]])
plot_data(itrack[[4]], aTrack[[4]], d0_a1_dat, gtrack[[4]], all[[4]], "Day 0 vs A1 at ", title[[4]])
plot_data(itrack[[5]], aTrack[[5]], d0_lsb_dat, gtrack[[5]], all[[5]], "Day 0 vs LSB at ", title[[5]])
plot_data(itrack[[6]], aTrack[[6]], d0_d30_dat, gtrack[[6]], all[[6]], "Day 0 vs Day 30 at ", title[[6]])
plot_data(itrack[[7]], aTrack[[7]], d0_d50_dat, gtrack[[7]], all[[7]], "Day 0 vs Day 50 at ", title[[7]])
plot_data(itrack[[8]], aTrack[[8]], d30_d50_dat, gtrack[[8]], all[[8]], "Day 30 vs Day 50 at ", title[[8]])


rm(list = ls())
gc()
