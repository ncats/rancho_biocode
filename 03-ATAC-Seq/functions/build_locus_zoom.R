#' functions for building a locus zoom plot
#' using gviz

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
                       str,
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
    strand = str,
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
    shape = "box",
    #remove connecting line
    col.line = "#FFFFFF")
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
  # create a IRanges object storing bam start and width data
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  # create a GRanges object storing bam genomic locations and annotations
  ranges <- GRanges(seqnames=Rle(bam$rname), 
                    ranges=ranges, 
                    strand=Rle(bam$strand),
                    flag=bam$flag, 
                    readid=bam$rname )
  
  # split GRanges object by strand type (+ or -)
  # grs <- split(ranges, strand(ranges))
  # calc coverage for each ranges and associated
  # range 1 and 2
  # cov <- lapply(grs[c("+", "-")], function(y)
  #   coverage(ranges(y),
  #            width=end(gr)
  #   ))
  cov <- coverage(ranges(ranges), width = end(gr))
  
  # sort start/end coverages for positive strand
  # pos_st <- sort(unique(unlist(lapply(cov, function(y) c(start(y))))))
  # pos_end <- sort(unique(unlist(lapply(cov, function(y) c(end(y))))))
  
  pos_st <- sort(unique(c(start(cov))))
  pos_end <- sort(unique(c(end(cov))))
  
  # # if there are no coverage lengths on the positive strand
  # if(length(pos) == 0){
  #   
  #   mcols(gr) <- DataFrame(plus = 0, minus = 0)
  #   return(gr)
  # } else{
  
  # create an updated GRanges object with selected
  # chromosome, range1 & range2 limits w/two metadata
  # columns specifying coverage for plus/minus strands
  comb <- GRanges(seqnames = seqnames(gr)[1],
                  ranges = IRanges(start = pos_st[-1], end = pos_end[-1]),
                  coverage = as.numeric(cov[pos_st[-1]]))
  # }
  return(comb)
}

# build a data track based on coverage from
# bam files
dat_track <- function(file, gen, style, title, chr, ylim) {
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
    #window = width(range(all_cont[[x]][[y]])),
    #window = "auto",
    window = -1,
    windowSize = 5, 
    # specify chromosome
    chromosome = chr,
    #set color for histogram as all gray
    fill.histogram = "#929292",
    col.histogram = "#929292",
    # set font text for all items in track to black
    col.axis = "black",
    #add a ylim to make all plots of consistent size
    ylim = ylim)
}