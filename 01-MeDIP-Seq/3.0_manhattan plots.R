#' Create manhattan plots for targeted genes:
#' elf5, dlk1, ube3a, igf2,
#' zfat, cdkn1c, proser2-as1
#' These are the only ones left of the 11

# req'd pkgs
x <- c("ggplot2", "tidyverse", "qqman", "data.table",
       "ggrepel")
sapply(x, library, character.only = TRUE)

source("./functions/build_manhattan.R")

# load shrinkage data
wa17 <- fread("./adj_data/wa17_deseq2_shrink.csv")
h9 <- fread("./adj_data/h9_deseq2_shrink.csv")

# build manhattan plots
build_manhattan(wa17, "MeDIP-Seq_manhattan_wa17_te_vs_wa17.png")
build_manhattan(h9, "MeDIP-Seq_manhattan_h9_te_vs_h9.png")

rm(list = ls())
gc()
