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
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.key.size = unit(2.5, 'cm'),
      panel.border = element_blank(),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 14, color = "black")
    ) + 
    guides(color = guide_legend(override.aes = list(size = 4))) +
    xlab("Chromosome") + ylim(c(0, max(-log10(df$P)) + 1))# +
    #mytheme
  
  ggsave(paste("./adj_data/plots/", filename, ".png", sep = ""), plot = p, width = 14, height = 9, units = "in")
  pdf(file = paste("./adj_data/plots/", filename, ".pdf", sep = ""), width = 14, height = 9)
  print(p)
  graphics.off()
}
