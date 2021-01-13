#' Build manhattan plot from DESeq2 output from 
#' MeDIP-Seq data

build_manhattan <- function(data, filename) {
  # all chromo
  df <- data %>%
    # select gene symbol, chromosome,
    # start of gene and adj p-value
    dplyr::select(c(symbol, chromosome_name_medip, start_position_medip, padj)) %>%
    # mutate X and Y chr to numeric
    mutate(chromosome_name_medip = ifelse(grepl("X", chromosome_name_medip), gsub("X", "23", chromosome_name_medip),
                             ifelse(grepl("Y", chromosome_name_medip), gsub("Y", "24", chromosome_name_medip), gsub("chr", "", chromosome_name_medip)))) %>%
    # convert factor columns to character
    mutate_if(is.factor, as.character) %>%
    mutate_if(is.numeric, as.character)
  
  # separate genes onto separate lines
  df <- rbindlist(lapply(1:nrow(df), function(x) {
    y <- df[x,]
    if (!grepl(";", y$symbol)) {
      names(y) <- c("symbol", "seqnames", "start", "padj")
      return(y)
    } else {
      sym <- unlist(strsplit(y$symbol, "; "))
      seq <- unlist(strsplit(y$chromosome_name_medip, "; "))
      st <- unlist(strsplit(y$start_position_medip, "; "))
      p <- rep(y$padj, length(sym))
      z <- data.frame("symbol" = as.character(sym),
                      "seqnames" = as.character(seq),
                      "start" = as.numeric(st),
                      "padj" = as.numeric(p)) %>%
        distinct(.)
      return(z)
    }
  }))
  
  # re-name cols (SNP is NOT a SNP, pulled
  # from someone else's code)
  names(df) <- c("SNP", "CHR", "BP", "P")
  
  df <- df %>%
    mutate_at(c(1:2), as.character) %>%
    mutate_at(c(3:4), as.numeric) %>%
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
    mutate(is_annotate=ifelse(P < 0.05, "yes", "no"))
  
  # determine axis plotting - taken from
  # someone else's code
  axisdf <- df %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
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
    
    # color sign genes with red (padj < 0.05)
    geom_point(data = subset(df, is_annotate == "yes"), color = "red", size = 2) +
    
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
  ggsave(paste("./adj_data/plots/manhattan/", filename, sep = ""), plot = p, width = 14, height = 10, units = "in")
}
