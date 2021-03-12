#' Build manhattan plot from DESeq2 output from 
#' MeDIP-Seq data

build_manhattan <- function(data, filename) {
  # all chromo
  df <- data %>%
    # select gene symbol, chromosome,
    # start of gene and adj p-value
    dplyr::select(c(gene_anotatr, gene_active, chromosome_name_medip, start_position_medip, padj)) %>%
    # mutate X and Y chr to numeric
    mutate(chromosome_name_medip = ifelse(grepl("X", chromosome_name_medip), gsub("X", "23", chromosome_name_medip),
                             ifelse(grepl("Y", chromosome_name_medip), gsub("Y", "24", chromosome_name_medip), gsub("chr", "", chromosome_name_medip)))) %>%
    # convert factor columns to character
    mutate_if(is.factor, as.character) %>%
    mutate_if(is.numeric, as.character)
  
  # separate genes onto separate lines
  df <- rbindlist(lapply(1:nrow(df), function(x) {
    y <- df[x,]
    anot <- unique(unlist(strsplit(y$gene_anotatr, "; ")))
    act <- unique(unlist(strsplit(y$gene_active, "; ")))
    all <- c(anot, act)
    all <- all[all != "NA"]
    
    if (length(all) == 1) {
      y <- y %>%
        dplyr::select(-gene_active) %>%
        mutate(gene_anotatr = all)
      names(y) <- c("symbol", "seqnames", "start", "padj")
      return(y)
    } else {
      seq <- rep(y$chromosome_name_medip, length(all))
      st <- rep(y$start_position_medip, length(all))
      p <- rep(y$padj, length(all))
      z <- data.frame("symbol" = as.character(all),
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
  
  # ggplot alternative: easier to adapt
  # update goi to remove: phlda2, h19, peg3, peg10
  # as they are not in the output
  goi <- toupper(c("elf5", "dlk1", "ube3a", "igf2",
                   "zfat", "cdkn1c", "proser2-as1",
                   "phlda2", "h19", "peg3", "peg10"))
  
  df <- df %>%
    mutate_at(c(1:2), as.character) %>%
    mutate_at(c(3:4), as.numeric) %>%
    arrange(CHR) %>%
    mutate(CHR = as.numeric(CHR))
  
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
    mutate(is_highlight = ifelse(SNP %in% goi & P < 0.05, "yes", "no"),
           is_annotate=ifelse(P < 0.05, "yes", "no"))
  # df$CHR <- factor(df$CHR, levels = c("1", "2", "3", "4", "5", 
  #                                     "6", "7", "8", "9", "10",
  #                                     "11", "12", "13", "14", "15",
  #                                     "16", "17", "18", "19", "20",
  #                                     "21", "22", "23", "24"))
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
    geom_point( aes(color = df$all_col), alpha=0.8, size=1.3) +
    scale_color_identity(guide = "none") +
    # color points alternating grey and blue 
    # by chr
    #scale_color_manual(values = rep(c("gray72", "blue"), 24)) +
    
    # custom X axis: 1 - 22, X and Y
    scale_x_continuous( label = c(seq(1,22), "X", "Y"),
                        breaks= axisdf$center ) +
    # drop axis at 0,0 pos on y-axis
    scale_y_continuous(expand = c(0, 0)) +
    
    # color goi with orange points
    geom_point(data=subset(df, is_highlight=="yes"), color = "red", size = 2, alpha = 0.5) +
    
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
  ggsave(paste("./adj_data/plots/manhattan/", filename, sep = ""), plot = p, width = 14, height = 10, units = "in")
}
