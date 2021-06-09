#' Create manhattan plots across all chromosomes 
#' for each contrast

# req'd pkgs
x <- c("ggplot2", "tidyverse", "qqman", "data.table",
       "ggrepel", "gtools", "scales")
sapply(x, library, character.only = TRUE)

source("./functions/upd_build_manhattan.R")

# load shrinkage data
ipsc_teexp <- fread("./adj_data/deseq_shrink/filtered/J91i-iPSC_TEexp_deseq_shrink_filtered.csv") %>%
  # not sure this is needed
  #dplyr::select(c(1:10)) %>%
  dplyr::rename("peak" = "MergedRegion") %>%
  arrange(padj, abs(log2FoldChange)) %>%
  mutate(type = ifelse(log2FoldChange < 0, "TEexp",
                       ifelse(log2FoldChange > 0, "iPSC", "NA")))
names(ipsc_teexp)[c(8:10)] <- c("chr", "st", "end")


# build manhattan plots
upd_build_manhattan(ipsc_teexp, "Manhattan_TEexp_vs_iPSC_J91l")

rm(list = ls())
gc()
