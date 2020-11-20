#' Build manhattan plot from DESeq2 output from 
#' MeDIP-Seq data

build_manhattan <- function(data, filename) {
  # all chromo
  df <- data %>%
    # select gene symbol, chromosome,
    # start of gene and adj p-value
    dplyr::select(c(symbol, seqnames, start, padj)) %>%
    # mutate X and Y chr to numeric
    mutate(seqnames = ifelse(grepl("chrX", seqnames), gsub("chrX", "23", seqnames),
                             ifelse(grepl("chrY", seqnames), gsub("chrY", "24", seqnames), gsub("chr", "", seqnames)))) %>%
    # convert factor columns to character
    mutate_if(is.factor, as.character)
  
  # separate genes onto separate lines
  df <- rbindlist(lapply(1:nrow(df), function(x) {
    y <- df[x,]
    if (!grepl(";", y$symbol)) {
      return(y)
    } else {
      sym <- unlist(strsplit(y$symbol, "; "))
      seq <- unlist(strsplit(y$seqnames, "; "))
      st <- unlist(strsplit(y$start, "; "))
      p <- rep(y$padj, length(sym))
      z <- data.frame("symbol" = sym,
                      "seqnames" = seq,
                      "start" = st,
                      "padj" = p) %>%
        distinct(.)
      return(z)
    }
  }))
  
  # re-name cols (SNP is NOT a SNP, pulled
  # from someone else's code)
  names(df) <- c("SNP", "CHR", "BP", "P")
  
  # ggplot alternative: easier to adapt
  # update goi to remove: phlda2, h19, peg3, peg10
  # as they are not in the output
  goi <- toupper(c("elf5", "dlk1", "ube3a", "igf2",
                   "zfat", "cdkn1c", "proser2-as1"))
  
  df <- df %>%
    mutate_at(c(1:2), as.character) %>%
    mutate(BP = as.numeric(BP))
  
  df <- df %>%
    mutate(CHR = as.numeric(CHR)) %>%
    arrange(CHR)
  
  
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
    mutate(is_highlight = ifelse(SNP %in% goi, "yes", "no")) %>%
    mutate(is_annotate=ifelse(P < 0.05, "yes", "no"))
  # UBE3 is connected to another gene in annotation
  df$is_highlight <- ifelse(grepl("UBE3A", df$SNP), "yes", df$is_highlight)
  df$SNP <- ifelse(grepl("UBE3A", df$SNP), "UBE3A", df$SNP)
  
  # determine axis plotting - taken from
  # someone else's code
  axisdf <- df %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # this labels goi only once, rather than 
  # multiple times in the plot
  labs <- subset(df, is_highlight == "yes") %>%
    arrange(SNP)
  sum <- labs %>% group_by(SNP) %>% summarise(n = n())
  vec <- unlist(lapply(1:length(sum$SNP), function(x) {
    y <- c(sum$SNP[x], rep("", sum$n[x] - 1))
  }))
  labs$upd <- vec
  
  # this is for all chromosomes
  p <- ggplot(df, aes(x = BPcum, y = -log10(P))) + 
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    # color points alternating grey and blue 
    # by chr
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis: 1 - 22, X and Y
    scale_x_continuous( label = c(seq(1,22), "X", "Y"),
                        #label = axisdf$CHR, 
                        breaks= axisdf$center ) +
    # drop axis at 0,0 pos on y-axis
    scale_y_continuous(expand = c(0, 0)) +
    
    # color goi with orange points
    geom_point(data=subset(df, is_highlight=="yes"), color="orange", size=2) +
    
    # color sign genes with red (padj < 0.05)
    geom_point(data = subset(df, is_annotate == "yes"), color = "red", size = 2) +
    # repel labels so text does not overlap
    geom_label_repel(data = labs, aes(label = upd), size=4) +
    
    # customize the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 14, color = "black")
    ) + xlab("Chromosome") + ylim(c(0, max(-log10(df$P)) + 1))
  ggsave(paste("./adj_data/", filename, sep = ""), plot = p, width = 11, height = 8, units = "in")
}
