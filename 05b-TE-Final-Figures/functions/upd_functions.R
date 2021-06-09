#' Build manhattan plot from DESeq2 output from 
#' MeDIP-Seq data

upd_build_manhattan <- function(data, filename) {
  # all chromo
  df <- data %>%
    # select gene symbol, chromosome,
    # start of gene and adj p-value
    dplyr::select(c(Gene.List, chr, st, padj, type)) %>%
    # mutate X and Y chr to numeric
    mutate(chr = ifelse(grepl("X", chr), gsub("X", "23", chr),
                        ifelse(grepl("Y", chr), gsub("Y", "24", chr), gsub("chr", "", chr)))) %>%
    # convert factor columns to character
    mutate_if(is.factor, as.character) %>%
    mutate_if(is.numeric, as.character)
  
  # re-name cols (SNP is NOT a SNP, pulled
  # from someone else's code)
  names(df)[1:4] <- c("SNP", "CHR", "BP", "P")
  
  goi <- df$SNP[1:20]
  
  df <- df %>%
    mutate_at(c(1:2), as.character) %>%
    mutate_at(c(3:4), as.numeric) %>%
    #arrange(CHR) %>%
    mutate(CHR = as.numeric(CHR))
  
  cols <- hue_pal()(2)
  
  # this code is used to plot the chrom
  # correctly in manhattan plot (taken from 
  # other ppl's code)
  df <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = as.numeric(max(BP))) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(df, .) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot) %>%
    distinct(.) %>%
    arrange(P) %>%
    mutate(is_top = c(rep("yes", 20), rep("no", length(df$SNP) - 20)),
           is_annotate = ifelse(P < 0.05, "yes", "no")) %>%
    mutate(is_top = ifelse(is_top == "yes" & type == "iPSC", cols[1],
                      ifelse(is_top == "yes" & type == "TEexp", cols[2], NA))) %>%
    mutate(upd = ifelse(is_top %in% cols, SNP, NA))
  df$SNP <- ifelse(is.na(df$SNP), "NA", df$SNP)

  sub <- df %>%
    group_by(CHR, is_annotate) %>%
    summarise(n = n())
  sub$col <- rep(c("gray32", "gray72", "dark blue", "blue"), times = round(nrow(sub)/4))[1:nrow(sub)]
  sub <- unlist(sapply(1:nrow(sub), function(x) rep(sub$col[x], sub$n[x])))
  df <- df %>% 
    arrange(CHR, is_annotate) %>%
    mutate(all_col = sub) %>%
    arrange(CHR)
  
  # determine axis plotting - taken from
  # someone else's code
  axisdf <- df %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  df$CHR <- factor(df$CHR, levels = c("1", "2", "3", "4", "5", 
                                      "6", "7", "8", "9", "10",
                                      "11", "12", "13", "14", "15",
                                      "16", "17", "18", "19", "20",
                                      "21", "22", "23", "24"))
  
  # this is for all chromosomes
  p <- ggplot(df, aes(x = BPcum, y = -log10(P))) + 
    geom_point( aes(color = df$all_col), alpha=0.8, size=1.3) +
    scale_color_identity("TE differentiation stage", 
                         guide = "legend", 
                         labels = c("Methylated in iPSC", "Methylated in TEexp", rep("", 4)),
                         breaks = c(cols[1], cols[2], rep(NA, 4))) +
    # custom X axis: 1 - 22, X and Y
    scale_x_continuous( label = c(seq(1,22), "X", "Y"),
                        breaks= axisdf$center ) +
    # drop axis at 0,0 pos on y-axis
    scale_y_continuous(expand = c(0, 0)) +
    
    # color goi with red points
    geom_point(data = df, aes(color = is_top), size = 3) +
    
    # repel labels so text does not overlap
    geom_label_repel(data = df, aes(label = upd), size = 4, max.overlaps = 20) +
    
    # customize the theme:
    theme_bw() +
    theme( 
      legend.position = "bottom",
      panel.border = element_blank(),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 14, color = "black")
    ) + xlab("Chromosome") + ylim(c(0, max(-log10(df$P)) + 1)) +
    mytheme
  
  ggsave(paste("./adj_data/plots/", filename, ".png", sep = ""), plot = p, width = 14, height = 10, units = "in")
  pdf(file = paste("./adj_data/plots/", filename, ".pdf", sep = ""), width = 12, height = 8)
  print(p)
  graphics.off()
}

#' Save pheatmap fxn
save_pheatmap_pdf <- function(x, filename, width=12, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  cairo_pdf(filename, width = width, height = height, family = "NotoSans-Condensed")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

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
      # sets a fixed window for the track
      window = "auto",
      # specify chromosome
      chromosome = chr,
      background.title = "white",
      cex.title = 1,
      # set font text for all items in track to black
      col.axis = "black",
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
      background.title = "white",
      cex.title = 1,
      # set font text for all items in track to black
      col.axis = "black",
      col.title = col1,
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
                  chr = file[[x]][[y]]$seqnames,
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
                       btrack_files, titles, loc) {
  # pool itrack and atrack
  first <- lapply(1:length(itrack_files), function(x)
    list(itrack_files[[x]], atrack_files[[x]], gtrack_files[[x]], btrack_files[[x]]))
  # then add data track and genome axis track
  comb <- lapply(1:length(dat_files), function(x)
    c(first[[x]], dat_files[[x]]))
  
  # create plot name
  plot_nam <- lapply(1:length(titles), function(x)
    paste0("./adj_data/plots/", loc, "/", titles[[x]], ".pdf", sep = "")
  )
  
  # create final combined pdf
  lapply(1:length(comb), function(x) {
    CairoPDF(file = plot_nam[[x]], width = 8, height = 11)
    plotTracks(comb[[x]],
               main = gsub("_", " ", titles[[x]]),
               background.title = "white",
               fontfamily = "NotoSans-Condensed",
               # set font size for y-axis (title)
               #fontsize = 12,
               sizes=c(0.8, 0.8, 0.8, 0.8, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7),
               stacking = "squish"
    )
    graphics.off()
  })
}

