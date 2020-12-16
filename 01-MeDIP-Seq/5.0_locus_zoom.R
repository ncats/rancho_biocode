# Build locus zoom plots using Gviz for 7/10
# targeted genes/imprinted genes
# goi: elf5, dlk1, ube3a, igf2,
# zfat, cdkn1c, proser2-as1

# req'd pkgs
x <- c("Gviz", "GenomicRanges", "rtracklayer", 
       "Rsamtools", "GenomicAlignments",
       "dplyr", "data.table", "Cairo")
sapply(x, library, character.only = TRUE)

# source locus zoom functions
source("./functions/build_locus_zoom.R")

# load anotatr file to pull goi loci
annot <- data.frame(readRDS("./adj_data/anotatr_annotations.RDS"))
# create vector of goi left
goi <- toupper(c("elf5", "dlk1", "ube3a", "igf2",
                 "zfat", "cdkn1c", "proser2-as1"))
# filter to retain only goi and the first/last
# sites and retain only chr 1 - 22, X and Y
# from anotatr file
annot <- annot %>%
  filter(symbol %in% goi) %>%
  arrange(symbol, start) %>%
  group_by(seqnames, symbol) %>%
  dplyr::slice(c(1, n())) %>%
  filter(!grepl("_", seqnames)) 
# split annotation file by symbol
annot <- split(annot, annot$symbol)
# pull the start/end sites from that goi
# and add upstream and downstream 1000 bp
# this targets the region of interest for
# a goi
annot <- rbindlist(lapply(1:length(annot), function(x) {
  y <- annot[[x]]
  y <- data.frame(seqnames = y$seqnames[1],
                  start = y$start[1] - 1000,
                  end = y$end[2] + 1000,
                  symbol = y$symbol[1],
                  strand = y$strand[1]
                  )
}))

# now read in filtered merged peaks
medip_peaks <- import.bed("./adj_data/MergedPeaks.bed")

########## Annotation Track ############
########## Utilizes .bed file from Active Motif
########## plot specified medipseq peaks using an annotation track ##########

all_peaks <- data.frame(medip_peaks)
# extract targeted roi from goi from
# MergedReg (medip-seq peaks) filtered
# bed file
peaks <- lapply(1:nrow(annot), function(x)
  all_peaks %>%
    filter(as.character(seqnames) == as.character(annot$seqnames[x])) %>%
    filter(start >= annot$start[x]) %>%
    filter(end <= annot$end[x]) %>%
    mutate(strand = annot$strand[x])
)

# if there are > 3 peaks, split out
# peaks into elements of a list contain-
# ing 3 peaks in ea
peaks <- lapply(peaks, function(x) {
  y <- x
  if (nrow(y) > 3) {
    len <- ceiling(nrow(y)/3)
    id <- rep(c(1:len), each = 3)
    y$id <- id[c(1:nrow(y))]
    y <- split(y, y$id)
    y <- lapply(y, function(z)
    z %>% dplyr::select(-id))
  } else {
    z <- y
  }
})

# re-name MergedReg_ to Peak for all 
# elements of the list at the same time
nam <- lapply(1:length(peaks), function(x) {
  t <- peaks[[x]]
  if (!is.data.frame(t)) {
    y <- lapply(1:length(t), function(y) {
      y <- t[[y]] %>%
        pull(name)
      y <- gsub("MergedReg_", "Peak ", y)
    })
    return(y)
  } else {
    y <- t %>%
      pull(name)
    y <- gsub("MergedReg_", "Peak ", y)
    return(y)
  }
})

# build the annotation track for all
aTrack <- unlist(lapply(1:length(peaks), function(x) {
  t <- peaks[[x]]
  nam <- nam[[x]]
  if (!is.data.frame(t)) {
    aTrack <- unlist(lapply(1:length(t), function(y) 
      anno_track(st = t[[y]]$start,
                 wid = t[[y]]$width, 
                 chr = t[[y]]$seqnames[1],
                 str = t[[y]]$strand,
                 plot_name = nam[[y]],
                 gen = "hg38",
                 title = "MeDIP-Seq Peaks",
                 feat = "Peak ",
                 len = length(t[[y]]$start))
    ))
  } else {
    aTrack <- anno_track(st = t$start,
                         wid = t$width, 
                         chr = t$seqnames[1],
                         str = t$strand,
                         plot_name = nam,
                         gen = "hg38",
                         title = "MeDIP-Seq Peaks",
                         feat = "Peak ",
                         len = length(t$start))
  }
}))

########## Ideogram Track ############
########## This plots the chromosome feature at
########## the top of the plot - the red line is
########## the loc of the zoomed plots below (i.e.
########## peaks and coverage)

# for genes w/ many medip-seq peaks,
# generate multiple rows for these goi
len <- rbindlist(lapply(1:length(nam), function(x) {
  t <- nam[[x]]
  if (is.list(t)) {
    y <- annot[x,]
    len <- length(t)
    y <- rbind(y, y[rep(1, len - 1), ])
  } else {
    return(annot[x,])
  }
}))

# create ideogram track for each row in len
itrack <- lapply(1:nrow(len), function(x) 
  ideo_track(gen = 'hg38',
             chr = len$seqnames[x])
)

########## Genome Axis Track ############

# collapse lists w/in lists into a vector
f <- function(x) {
  if (is.atomic(x)) {
    list(x)
  } else {
    x
  }
}

# pull st/end site for a given roi in a given
# goi for all medip-seq peaks
range <- lapply(1:length(peaks), function(x) {
  t <- peaks[[x]]
  if (!is.data.frame(t)) {
    range <- lapply(1:length(t), function(y) {
      range <- t[[y]] %>%
        dplyr::slice(c(1, n())) %>%
        dplyr::select(start, end)
      range <- c(range$start[1], range$end[2])
    })
  } else {
    range <- t %>%
      dplyr::slice(c(1, n())) %>%
      dplyr::select(start, end)
    range <- c(range$start[1], range$end[2])
  }
})
# collapse them into a single long list = 14
# ties out to nrow(len)
range <- unlist(lapply(range, f), recursive = F)

# build genome axis track for each
gtrack <- lapply(1:length(range), function(x)
  GenomeAxisTrack(#seqnames = peaks[[x]]$seqnames[1], 
    range=IRanges(start = range[[x]][1],
                  end = range[[x]][2])
  )
)

########## Data Track ############
########## Utilizes .bam files from Active Motif or Rsubread
########## as input to calc coverage & plot ##########

# pull coverage from bam files
bam_files <- list.files("./BAM/", pattern = ".BAM$", full.names = T)
# re-name bam files
bam_nam <- gsub(".*NIH_|_medip.*", "", bam_files)
# split out into h9 contrasts and wa17 contrasts
bam_h9 <- bam_nam[1:6]
bam_wa <- bam_nam[7:12]

# for h9 vs h9-te
# calc coverage for a given goi, roi
all_h9 <- lapply(1:length(range), function(x) {
  unlist(lapply(1:6, function(y) {
    build_cov(chr = gsub("chr", "", len$seqnames[x]),
              range1 = range[[x]][1],
              range2 = range[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})

# build the data track for all goi in roi
# for h9-te vs h9
all_dtrack_h9 <- lapply(1:length(all_h9), function(x) {
  y <- all_h9[[x]]
  y <- rbindlist(lapply(1:length(y), function(t) {
    z <- y[[t]]
    z <- data.frame(z)
  }))
  rang <- c(min(y$coverage), max(y$coverage))
  unlist(lapply(1:6, function(y) 
    dat_track(all_h9[[x]][[y]],
              gen = "hg38",
              style = "histogram",
              title = bam_h9[y],
              chr = gsub("chr", "", len$seqnames[x]),
              lim = c(rang[1], rang[2])
    )
  ), recursive = F)
})

# for wa17 vs wa17-te
# calc coverage for a given goi, roi

all_wa <- lapply(1:length(range), function(x) {
  unlist(lapply(7:12, function(y) {
    build_cov(chr = gsub("chr", "", len$seqnames[x]),
              range1 = range[[x]][1],
              range2 = range[[x]][2],
              file = bam_files[y])
  }), recursive = F)
})

# build the data track for all goi in roi
# for wa17-te vs wa17
all_dtrack_wa <- lapply(1:length(all_wa), function(x) {
  y <- all_wa[[x]]
  y <- rbindlist(lapply(1:length(y), function(t) {
    z <- y[[t]]
    z <- data.frame(z)
  }))
  rang <- c(min(y$coverage), max(y$coverage))
  unlist(lapply(1:6, function(y) 
    dat_track(all_wa[[x]][[y]],
              gen = "hg38",
              style = "histogram",
              title = bam_wa[y],
              chr = gsub("chr", "", len$seqnames[x]),
              lim = c(rang[1], rang[2])
    )
  ), recursive = F)
})

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf

# pool itrack and atrack
first <- lapply(1:length(itrack), function(x)
  list(itrack[[x]], aTrack[[x]]))
# then add data track and genome axis track
comb_h9 <- lapply(1:length(all_dtrack_h9), function(x)
  c(first[[x]], all_dtrack_h9[[x]], gtrack[[x]]))

# for goi that have more than 3 medip-seq peaks
# create a title with _2, or _3 so the plots
# are not written over
# id dup goi symbol in len
ind_dup <- which(duplicated(len$symbol)) 
# put into df
dup <- data.frame("symbol" = len$symbol[duplicated(len$symbol)])
# determine the count of duplication of 
# a given goi
ct <- dup %>% group_by(symbol) %>% summarize(ct = n())
# iterate +1 to each count for each goi duplicated
ct <- unlist(sapply(1:nrow(ct), function(x) seq(ct$ct[x]) + 1))
# merge counts into duplicate df, add _2, _3, etc
dup <- dup %>%
  cbind(., ct) %>%
  mutate(symbol = paste(symbol, "_", ct, sep = "")) %>%
  dplyr::select(-ct)
len$symbol <- as.character(as.factor(len$symbol))
# add duplicated name _2, _3 etc to len
len[ind_dup, 4] <- dup$symbol

# create title for h9 plots
title_h9 <- lapply(1:nrow(len), function(x)
  paste("H9 vs H9-TE at ", len$symbol[x], " locus", sep = ""))

# create plot name for h9 plots
plot_nam_h9 <- lapply(1:nrow(len), function(x)
  paste("./adj_data/imprinted_genes/h9_", as.character(len$symbol[x]), ".pdf", sep = "")
  )

# create final combined pdf for h9 plots
lapply(1:length(comb_h9), function(x) {
  CairoPDF(file = plot_nam_h9[[x]], width = 20, height = 20, family = "NotoSans-Condensed")
  plotTracks(comb_h9[[x]],
             main = title_h9[[x]],
             # specify color of feature for peaks
             `Peak ` = "red",
             # set fill for y-axis (title)
             background.title = "white",
             # set font color for y-axis (title)
             col.title = "black",
             # set font family for y-axis (title)
             fontfamily.title = "NotoSans-Condensed",
             # set font size for y-axis (title)
             fontsize = 16,
             sizes=c(0.5, 0.5, 1, 1, 1, 1, 1, 1, 0.5)
  )
  graphics.off()
})

# create title for wa17 plots
title_wa <- lapply(1:nrow(len), function(x)
  paste("WA17 vs WA17-TE at ", len$symbol[x], " locus", sep = ""))

# create plot name for wa17 plots
plot_nam_wa <- lapply(1:nrow(len), function(x)
  paste("./adj_data/imprinted_genes/wa17_", as.character(len$symbol[x]), ".pdf", sep = "")
)

# add data and genome axis tracks to ideogram
# and annotation track
comb_wa <- lapply(1:length(all_dtrack_wa), function(x)
  c(first[[x]], all_dtrack_wa[[x]], gtrack[[x]]))

# create final combined pdf for wa17 plots
lapply(1:length(comb_wa), function(x) {
  CairoPDF(file = plot_nam_wa[[x]], width = 20, height = 20, family = "NotoSans-Condensed")
  plotTracks(comb_wa[[x]],  
             main = title_wa[[x]],
             # specify color of feature for peaks
             # linked to group feature in annotation
             # track
             `Peak ` = "red",
             # set fill for y-axis (title)
             background.title = "white",
             # set font color for y-axis (title)
             col.title = "black",
             # set font family for y-axis (title)
             fontfamily.title = "NotoSans-Condensed",
             # set font size for y-axis (title)
             fontsize = 16,
             sizes = c(0.5, 0.5, 1, 1, 1, 1, 1, 1, 0.5)
  )
  graphics.off()
})

rm(list = ls())
gc()
