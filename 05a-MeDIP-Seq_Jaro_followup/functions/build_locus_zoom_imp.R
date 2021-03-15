#' functions for building a locus zoom plot
#' using gviz: ideogram, annotation, genomic
#' axis and data track (coverage)

# build an ideogram track
ideo_track <- function(gen, chr) {
  itrack <- IdeogramTrack(
                          # specify which genome/version you want
                          genome = gen, 
                          # specify chromosome you are targeting
                          chromosome = chr,
                          # font
                          fontfamily = "NotoSans-Condensed",
                          # black font color
                          fontcolor = "#000000")
}

# build an annotation track: in this case,
# it's a track of peaks (ATAC, MeDIP, etc)
anno_track <- function(st = c(0, 1, 2),
                       wid = c(4, 5, 6),
                       chr,
                       #str = NULL,
                       plot_name,
                       gen,
                       title,
                       feat,
                       len
                       ) {
  track <- AnnotationTrack(
    # specify start locs for peaks
    start = st,
    # specify the width of the peak
    width = wid, 
    # specify chromosome
    chromosome = chr, 
    # specify strand info if it's avail
    #strand = str,
    # this groups peaks into 1 row
    group = plot_name,
    # specify annotation track
    genome = gen, 
    # specify y-axis label (title)
    name = title,
    # specify grouping annotation
    groupAnnotation = "group", 
    # specify loc of text annotation for peak
    just.group = "above",
    # sets font color for annotation of peaks to black
    fontcolor.group = "#000000",
    # sets font size for annotation of peak 
    # (diff than title size)
    fontsize.group = 24,
    # specify font family
    fontfamily.group = "NotoSans-Condensed",
    # specify border color for peak
    col = "black",
    # specify type of shape for annotation
    shape = "box")
  # specify a feature in elementMetadata
  # this allows you to color grouped annotations
  # such as medip-seq peaks
  feature(track) <- rep(feat, len)
  return(track)
}

# load bam file and calc coverage for plotting
build_cov <- function(chr, range1, range2, file) {
  require(Rsamtools)
  # specify what chromosome and range of values 
  # to pull from bam (instead of pulling the entire
  # bam file, which is large & takes a long time)
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(start = range1, end = range2))
  # specify a bam flag
  flag <- scanBamFlag(isUnmappedQuery = FALSE)
  # specify columns to pull from bam file (all of them)
  what <- c("qname", "flag", "rname", "strand", 
            "pos", "qwidth", "mapq", "cigar",
            "mrnm", "mpos", "isize", "seq", "qual")
  # set parameters for pulling  bam file: which chromo-
  # some and range, what columns from the bam file and
  # the flag
  param <- ScanBamParam(which = gr, what = what, flag = flag)
  # pull data from the bam file using parameters specified;
  # requires [[1]]
  bam <- scanBam(file,
                 param = param)[[1]]
  
  if (length(bam$qname) == 0) {
    comb <- GRanges(seqnames = seqnames(gr)[1],
                    ranges = IRanges(start = 0, end = 0),
                    coverage = 0)
    return(comb)
  } else {
    # create a IRanges object storing bam start and width data
    ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
    # create a GRanges object storing bam genomic locations and annotations
    ranges <- GRanges(seqnames=Rle(bam$rname), 
                      ranges=ranges, 
                      strand=Rle(bam$strand),
                      flag=bam$flag, 
                      readid=bam$rname )
    
    # calc coverage for a given chr, range
    cov <- coverage(ranges(ranges), width = end(gr))
    
    # sort start/end coverages
    st <- sort(unique(c(start(cov))))
    end <- sort(unique(c(end(cov))))
    
    # create an updated GRanges object with selected
    # chromosome, range1 & range2 limits w/two metadata
    # columns specifying coverage for plus/minus strands
    comb <- GRanges(seqnames = seqnames(gr)[1],
                    ranges = IRanges(start = st[-1], end = end[-1]),
                    coverage = as.numeric(cov[st[-1]]))
    return(comb)
  }
}

# build a data track based on coverage from
# bam files
dat_track <- function(file, gen, style, title, chr, lim) {
  DataTrack(
            # this is the GRanges object containing coverage
            # for plus & minus strands
            range = file,
            # specify genome and version
            genome = gen,
            # specify type of data track (i.e. line, histogram, etc)
            type = style,
            # set title of track (y-axis)
            name = title,
            # sets a fixed window for the track
            window = "auto",
            # specify chromosome
            chromosome = chr,
            #col = "#000000"
            # set font text for all items in track to black
            col.axis = "black",
            ylim = lim)
}

# these 2 fxns allow you to iterate thru many 
# comparisons of the build_cov and dat_track 
# functions (i.e. so your code isn't super long
# and redundant)
multi_comp_cov <- function(file, bam_comps, ranges, bam_file) {
  y <- lapply(1:nrow(file), function(x) {
    unlist(lapply(1:length(bam_comps), function(y) {
      build_cov(chr = paste("chr", file$chromosome_name_medip[x], sep = ""),
                range1 = ranges[[x]][[1]],
                range2 = ranges[[x]][[2]],
                file = bam_file[y])
    }), recursive = F)
  })
}

# these 2 fxns allow you to iterate thru many 
# comparisons of the build_cov and dat_track 
# functions (i.e. so your code isn't super long
# and redundant)
multi_comp_cov2  <- function(file, bam_comps, ranges, bam_file) {
  y <- lapply(1:nrow(file), function(x) {
    unlist(lapply(1:length(bam_comps), function(y) {
      build_cov(chr = paste("chr", file$chr[x], sep = ""),
                range1 = ranges[[x]][[1]],
                range2 = ranges[[x]][[2]],
                file = bam_file[y])
    }), recursive = F)
  })
}

multi_comp_dat <- function(file, bam_comps) {
  a1_lsb_dat <- lapply(1:length(file), function(x) {
    y <- file[[x]]
    y <- rbindlist(lapply(1:length(y), function(t) {
      z <- y[[t]]
      z <- data.frame(z)
    }))
    rang <- c(min(y$coverage), max(y$coverage))
    unlist(lapply(1:length(bam_comps), function(y) 
      dat_track(file[[x]][[y]],
                gen = "hg38",
                style = "histogram",
                title = bam_comps[y],
                chr = gsub("chr", "", data.frame(file[[x]][[y]])$seqnames),
                lim = c(rang[1], rang[2])
      )
    ), recursive = F)
  })
}

# combine data input for plot track across
# multi-comparisons
plot_data <- function(itrack_files, atrack_files, dat_files, gtrack_files,
                      raw_data, contrast, titles, loc) {
  # pool itrack and atrack
  first <- lapply(1:length(itrack_files), function(x)
    list(itrack_files[[x]], atrack_files[[x]]))
  # then add data track and genome axis track
  comb <- lapply(1:length(dat_files), function(x)
    c(first[[x]], dat_files[[x]], gtrack_files[[x]]))
  
  #con <- gsub(" at ", "", contrast)
  dir.exists(file.path("./adj_data/plots/", paste(loc, "/", contrast, sep = "")))
  
  ifelse(!dir.exists(file.path("./adj_data/plots/", paste(loc, "/", contrast, sep = ""))),
         dir.create(file.path("./adj_data/plots/", paste(loc, "/", contrast, sep = ""))), FALSE)
  
  # create plot name
  plot_nam <- lapply(1:length(titles), function(x)
    paste(paste("./adj_data/plots/", loc, "/", contrast, "/", sep = ""), contrast, "_", titles[[x]], ".pdf", sep = "")
  )
  
    # create final combined pdf
    lapply(1:length(comb), function(x) {
      CairoPDF(file = plot_nam[[x]], width = 20, height = 20, family = "NotoSans-Condensed")
      plotTracks(comb[[x]],
                 main = titles[[x]],
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
}

# combine data input for plot track across
# multi-comparisons
plot_data2 <- function(itrack_files, atrack_files, dat_files, gtrack_files,
                      btrack_files, raw_data, contrast, titles) {
  # pool itrack and atrack
  first <- lapply(1:length(itrack_files), function(x)
    list(itrack_files[[x]], atrack_files[[x]], gtrack_files[[x]], btrack_files[[x]]))
  # then add data track and genome axis track
  comb <- lapply(1:length(dat_files), function(x)
    c(first[[x]], dat_files[[x]]))
  
  # create plot name
  plot_nam <- lapply(1:length(titles), function(x)
    paste(paste("./adj_data/plots/all_top10/", contrast, "/", sep = ""), contrast, "_", titles[[x]], ".pdf", sep = "")
  )
  
  # create final combined pdf
  lapply(1:length(comb), function(x) {
    CairoPDF(file = plot_nam[[x]], width = 28, height = 28)
    plotTracks(comb[[x]],
               main = gsub("_", " ", titles[[x]]),
               # specify color of feature for peaks
               #`Peak ` = "red",
               # set fill for y-axis (title)
               background.title = "white",
               # set font color for y-axis (title)
               col.title = "black",
               # set font family for y-axis (title)
               fontfamily.title = "NotoSans-Condensed",
               fontfamily = "NotoSans-Condensed",
               # set font size for y-axis (title)
               fontsize = 20,
               sizes=c(0.5, 1.2, 1.4, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
    )
    graphics.off()
  })
}
