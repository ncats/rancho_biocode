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
    cex = 1,
    # font
    #fontfamily = "NotoSans-Condensed",
    fontface = "bold",
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
    cex.title = 1,
    # sets font size for annotation of peak 
    # (diff than title size)
    #fontsize.group = 28,
    # specify font family
    fontfamily.group = "NotoSans-Condensed",
    # specify border color for peak
    col = "black",
    col.line = "black",
    # specify type of shape for annotation
    shape = "box",
    background.title = "white",
    col.title = "black")
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
dat_track_t <- function(file, gen, style, title, chr, lim, col1, bam_comps) {
  if (title == '') {
    DataTrack(
      # this is the GRanges object containing coverage
      # for plus & minus strands
      range = file,
      # specify genome and version
      genome = gen,
      # specify type of data track (i.e. line, histogram, etc)
      type = style,
      # set title of track (y-axis)
      #name = title,
      # sets a fixed window for the track
      window = "auto",
      # specify chromosome
      chromosome = chr,
      #col = "#000000"
      background.title = "white",
      cex.title = 1,
      # set font text for all items in track to black
      col.axis = "black",
      #col.title = col1,
      ylim = lim,
      col.histogram = col1,
      fill.histogram = col1,
      groups = bam_comps)
  } else {
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
      background.title = "white",
      cex.title = 1,
      # set font text for all items in track to black
      col.axis = "black",
      col.title = col1,
      # cex.title = 20,
      ylim = lim,
      col.histogram = col1,
      fill.histogram = col1,
      groups = bam_comps)
  }
}


multi_comp_dat_t <- function(file, bam_comps, adj_nam) {
  a1_lsb_dat <- lapply(1:length(file), function(x) {
    y <- file[[x]]
    y <- rbindlist(lapply(1:length(y), function(t) {
      z <- y[[t]]
      z <- data.frame(z)
    }))
    rang <- c(min(y$coverage), max(y$coverage))
    
    bar_col <- c("#7CAE00", "#F8766D", "#00BFC4")
    
    unlist(lapply(1:length(bam_comps), function(y) {
      if (grepl("input|Input", bam_comps[y])) {
        col1 <- bar_col[1]
      } else if (grepl("iPSC", bam_comps[y])) {
        col1 <- bar_col[2]
      } else if(grepl("TEexp", bam_comps[y])) {
        col1 <- bar_col[3]
      }
      dat_track_t(file[[x]][[y]],
                  gen = "hg38",
                  style = "histogram",
                  title = adj_nam[y],
                  chr = gsub("chr", "", file[[x]][[y]]@seqnames@values),
                  lim = c(rang[1], rang[2]),
                  col1 = col1,
                  bam_comps = bam_comps[y]
      )
      
    }), recursive = F)
  })
}

# combine data input for plot track across
# multi-comparisons
plot_data2 <- function(itrack_files, atrack_files, dat_files, gtrack_files,
                       btrack_files, titles, bam_comps) {
  # pool itrack and atrack
  first <- lapply(1:length(itrack_files), function(x)
    list(itrack_files[[x]], atrack_files[[x]], gtrack_files[[x]], btrack_files[[x]]))
  # then add data track and genome axis track
  comb <- lapply(1:length(dat_files), function(x)
    c(first[[x]], dat_files[[x]]))
  
  # create plot name
  plot_nam <- lapply(1:length(titles), function(x) {
    dir.create('./adj_data/plots/top10/', showWarnings = FALSE)
    paste("./adj_data/plots/top10/", titles[[x]], ".pdf", sep = "") 
  })
  
  # create final combined pdf
  lapply(1:length(comb), function(x) {
    CairoPDF(file = plot_nam[[x]], width = 8, height = 11)
    plotTracks(comb[[x]],
               main = gsub("_", " ", titles[[x]]),
               background.title = "white",
               fontfamily = "NotoSans-Condensed",
               # set font size for y-axis (title)
               #fontsize = 20,
               sizes=c(0.8, 0.8, 0.8, 0.8, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7),
               stacking = "squish"
    )
    graphics.off()
  })
}

