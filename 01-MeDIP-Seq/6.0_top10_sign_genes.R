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

# load deseq objects for h9 or wa17
h9 <- fread("./adj_data/h9_deseq2.csv") %>%
  arrange(padj) %>%
  dplyr::select(c(GeneID, log2FoldChange, padj,
                  chromosome_name_medip, start_position_medip,
                  end_position_medip, strand, width, id,
                  symbol))
# pull top 10 by padj
h9 <- h9[c(1:10),]
# split h9 deseq objects by gene (row)
h9 <- split(h9, seq(nrow(h9)))

# extract start/end for each goi, gene
# symbol and strand for h9
h9 <- rbindlist(lapply(1:length(h9), function(x) {
  y <- h9[[x]]
  st <- unlist(strsplit(y$start_position_medip, "; "))
  end <- unlist(strsplit(y$end_position_medip, "; "))
  sym <- unlist(strsplit(y$symbol, "; "))
  y <- data.frame("seqnames" = rep(y$chromosome_name_medip, length(st)),
                  "start" = as.numeric(st),
                  "end" = as.numeric(end),
                  "symbol" = sym,
                  "strand" = rep(y$strand, length(st))
  ) %>%
    mutate(symbol = gsub("-NA", "", symbol))
  y <- y %>%
    dplyr::slice(c(1, n()))
  y <- data.frame(seqnames = y$seqnames[1],
                  start = y$start[1],
                  end = y$end[2],
                  symbol = y$symbol[1],
                  strand = y$strand[1])
}))

# extract start/end for each goi, gene
# symbol and strand for wa17
wa <- fread("./adj_data/wa17_deseq2.csv") %>%
  arrange(padj) %>%
  dplyr::select(c(GeneID, log2FoldChange, padj,
                  chromosome_name_medip, start_position_medip,
                  end_position_medip, strand, width, id,
                  symbol))
wa <- wa[c(1:10),]
wa <- split(wa, seq(nrow(wa)))
wa <- rbindlist(lapply(1:length(wa), function(x) {
  y <- wa[[x]]
  st <- unlist(strsplit(y$start_position_medip, "; "))
  end <- unlist(strsplit(y$end_position_medip, "; "))
  sym <- unlist(strsplit(y$symbol, "; "))
  y <- data.frame("seqnames" = rep(y$chromosome_name_medip, length(st)),
                  "start" = as.numeric(st),
                  "end" = as.numeric(end),
                  "symbol" = sym,
                  "strand" = rep(y$strand, length(st))
  ) %>%
    mutate(symbol = gsub("-NA", "", symbol))
  y <- y %>%
    dplyr::slice(c(1, n()))
  y <- data.frame(seqnames = y$seqnames[1],
                  start = y$start[1],
                  end = y$end[2],
                  symbol = y$symbol[1],
                  strand = y$strand[1])
}))

# now read in filtered merged peaks
medip_peaks <- import.bed("./adj_data/MergedPeaks.bed")

########## Annotation Track ############
########## Utilizes .bed file from Active Motif
########## plot specified medipseq peaks using an annotation track ##########

all_peaks <- data.frame(medip_peaks)

# extract targeted roi from goi from
# MergedReg (medip-seq peaks) filtered
# bed file for h9
peaks_h9 <- lapply(1:nrow(h9), function(x) {
  y <- all_peaks %>%
    filter(as.character(gsub("chr", "", seqnames)) == as.character(h9$seqnames[x])) %>%
    filter(start >= as.numeric(h9$start[x])) %>%
    filter(end <= as.numeric(h9$end[x])) %>%
    mutate(strand = as.character(h9$strand[x]))
 })

# extract targeted roi from goi from
# MergedReg (medip-seq peaks) filtered
# bed file for wa17
peaks_wa <- lapply(1:nrow(wa), function(x)
  all_peaks %>%
    filter(as.character(gsub("chr", "", seqnames)) == as.character(wa$seqnames[x])) %>%
    filter(start >= as.numeric(wa$start[x])) %>%
    filter(end <= as.numeric(wa$end[x])) %>%
    mutate(strand = as.character(wa$strand[x]))
)

# if there are > 3 peaks, split out
# peaks into elements of a list contain-
# ing 3 peaks in ea for h9
peaks_h9 <- lapply(peaks_h9, function(x) {
  y <- x
  if (nrow(y) > 3) {
    len <- ceiling(nrow(y)/3)
    id <- rep(c(1:len), each = 3)
    y$id <- id[c(1:nrow(y))]
    y <- split(y, y$id)
    y <- lapply(y, function(z)
    z %>% dplyr::select(-id))
  } else {
    return(y)
  }
})

# if there are > 3 peaks, split out
# peaks into elements of a list contain-
# ing 3 peaks in ea for wa17
peaks_wa <- lapply(peaks_wa, function(x) {
  y <- x
  if (nrow(y) > 3) {
    len <- ceiling(nrow(y)/3)
    id <- rep(c(1:len), each = 3)
    y$id <- id[c(1:nrow(y))]
    y <- split(y, y$id)
    y <- lapply(y, function(z)
      z %>% dplyr::select(-id))
  } else {
    return(y)
  }
})

# re-name MergedReg_ to Peak for all 
# elements of the list at the same time
# for h9
nam_h9 <- lapply(1:length(peaks_h9), function(x) {
  t <- peaks_h9[[x]]
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

# re-name MergedReg_ to Peak for all 
# elements of the list at the same time
# for wa17
nam_wa <- lapply(1:length(peaks_wa), function(x) {
  t <- peaks_wa[[x]]
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

# build the annotation track for h9
aTrack_h9 <- unlist(lapply(1:length(peaks_h9), function(x) {
  t <- peaks_h9[[x]]
  nam <- nam_h9[[x]]
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

# build the annotation track for wa17
aTrack_wa <- unlist(lapply(1:length(peaks_wa), function(x) {
  t <- peaks_wa[[x]]
  nam <- nam_wa[[x]]
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
# for h9
len_h9 <- rbindlist(lapply(1:length(nam_h9), function(x) {
  t <- nam_h9[[x]]
  if (is.list(t)) {
    y <- h9[x,]
    len <- length(t)
    y <- rbind(y, y[rep(1, len - 1), ])
  } else {
    return(h9[x,])
  }
}))

# for genes w/ many medip-seq peaks,
# generate multiple rows for these goi
# for wa17
len_wa <- rbindlist(lapply(1:length(nam_wa), function(x) {
  t <- nam_wa[[x]]
  if (is.list(t)) {
    y <- wa[x,]
    len <- length(t)
    y <- rbind(y, y[rep(1, len - 1), ])
  } else {
    return(wa[x,])
  }
}))

# create ideogram track for each row in len_h9
itrack_h9 <- lapply(1:nrow(len_h9), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste("chr", len_h9$seqnames[x], sep = ""))
)

# create ideogram track for each row in len_wa
itrack_wa <- lapply(1:nrow(len_wa), function(x) 
  ideo_track(gen = 'hg38',
             chr = paste("chr", len_wa$seqnames[x], sep = ""))
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
# goi for all medip-seq peaks for h9
range_h9 <- lapply(1:length(peaks_h9), function(x) {
  t <- peaks_h9[[x]]
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
# collapse them into a single long list
# ties out to nrow(len_h9)
range_h9 <- unlist(lapply(range_h9, f), recursive = F)

# pull st/end site for a given roi in a given
# goi for all medip-seq peaks for wa17
range_wa <- lapply(1:length(peaks_wa), function(x) {
  t <- peaks_wa[[x]]
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
# collapse them into a single long list
# ties out to nrow(len_wa)
range_wa <- unlist(lapply(range_wa, f), recursive = F)

# build genome axis track for h9
gtrack_h9 <- lapply(1:length(range_h9), function(x)
  GenomeAxisTrack(#seqnames = peaks_h9[[x]]$seqnames[1], 
    range=IRanges(start = range_h9[[x]][1],
                  end = range_h9[[x]][2])
  )
)

# build genome axis track for wa17
gtrack_wa <- lapply(1:length(range_wa), function(x)
  GenomeAxisTrack(#seqnames = peaks[[x]]$seqnames[1],
    range=IRanges(start = range_wa[[x]][1],
                  end = range_wa[[x]][2])
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
all_h9 <- lapply(1:length(range_h9), function(x) {
  unlist(lapply(1:6, function(y) {
    build_cov(chr = gsub("chr", "", len_h9$seqnames[x]),
              range1 = range_h9[[x]][1],
              range2 = range_h9[[x]][2],
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
              chr = gsub("chr", "", len_wa$seqnames[x]),
              lim = c(rang[1], rang[2])
    )
  ), recursive = F)
})

# for wa17 vs wa17-te
# calc coverage for a given goi, roi
all_wa <- lapply(1:length(range_wa), function(x) {
  unlist(lapply(7:12, function(y) {
    build_cov(chr = gsub("chr", "", len_wa$seqnames[x]),
              range1 = range_wa[[x]][1],
              range2 = range_wa[[x]][2],
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
              chr = gsub("chr", "", len_wa$seqnames[x]),
              lim = c(rang[1], rang[2])
    )
  ), recursive = F)
})

########## plot all data! ############
########## combine itrack, atrack, gtrack and data tracks
########## into a single pdf
# pool itrack and atrack
first_h9 <- lapply(1:length(itrack_h9), function(x)
  list(itrack_h9[[x]], aTrack_h9[[x]]))
# then add data track and genome axis track for h9
comb_h9 <- lapply(1:length(all_dtrack_h9), function(x)
  c(first_h9[[x]], all_dtrack_h9[[x]], gtrack_h9[[x]]))

# for goi that have more than 3 medip-seq peaks
# create a title with _2, or _3 so the plots
# are not written over
# id dup goi symbol in len
ind_dup <- which(duplicated(len_h9$symbol)) 
# put into df
dup <- data.frame("symbol" = len_h9$symbol[duplicated(len_h9$symbol)])
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
len_h9$symbol <- as.character(as.factor(len_h9$symbol))
# add duplicated name _2, _3 etc to len
len_h9[ind_dup, 4] <- dup$symbol

# create title for h9 plots
title_h9 <- lapply(1:nrow(len_h9), function(x) {
  paste("H9 vs H9-TE at ", len_h9$symbol[x], " locus", sep = "")
})

# create plot name for h9 plots
plot_nam_h9 <- lapply(1:nrow(len_h9), function(x)
  paste("./adj_data/top10/h9_", as.character(len_h9$symbol[x]), ".pdf", sep = "")
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

# see annotation above - just applied to wa17
ind_dup <- which(duplicated(len_wa$symbol)) 
dup <- data.frame("symbol" = len_wa$symbol[duplicated(len_wa$symbol)])
ct <- dup %>% group_by(symbol) %>% summarize(ct = n())
ct <- unlist(sapply(1:nrow(ct), function(x) seq(ct$ct[x]) + 1))
dup <- dup %>%
  cbind(., ct) %>%
  mutate(symbol = paste(symbol, "_", ct, sep = "")) %>%
  dplyr::select(-ct)
len_wa$symbol <- as.character(as.factor(len_wa$symbol))
len_wa[ind_dup, 4] <- dup$symbol

plot_nam_wa <- lapply(1:nrow(len_wa), function(x)
  paste("./adj_data/top10/wa17_", as.character(len_wa$symbol[x]), ".pdf", sep = "")
)

title_wa <- lapply(1:nrow(len_wa), function(x)
  paste("WA17 vs WA17-TE at ", len_wa$symbol[x], " locus", sep = ""))

first_wa <- lapply(1:length(itrack_wa), function(x)
  list(itrack_wa[[x]], aTrack_wa[[x]]))
comb_wa <- lapply(1:length(all_dtrack_wa), function(x)
  c(first_wa[[x]], all_dtrack_wa[[x]], gtrack_wa[[x]]))

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
