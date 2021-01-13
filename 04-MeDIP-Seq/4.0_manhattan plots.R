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
day0_a1 <- fread("./adj_data/deseq_shrink/filtered/hPSCs_day0_A1_shrink_filtered.csv")
day0_lsb <- fread("./adj_data/deseq_shrink/filtered/hPSCs_day0_LSB_shrink_filtered.csv")
day0_d30 <- fread("./adj_data/deseq_shrink/filtered/hPSCs_day0_NCRM5-D30_shrink_filtered.csv")
day0_d50 <- fread("./adj_data/deseq_shrink/filtered/hPSCs_day0_NCRM5-D50_shrink_filtered.csv")
a1_lsb <- fread("./adj_data/deseq_shrink/filtered/A1_LSB_shrink_filtered.csv")
a1_d30 <- fread("./adj_data/deseq_shrink/filtered/A1_NCRM5-D30_shrink_filtered.csv")
a1_d50 <- fread("./adj_data/deseq_shrink/filtered/A1_NCRM5-D50_shrink_filtered.csv")
d30_d50 <- fread("./adj_data/deseq_shrink/filtered/NCRM5-D30_NCRM5-D50_shrink_filtered.csv")

# build manhattan plots
build_manhattan(day0_a1, "MeDIP-Seq_manhattan_A1_vs_hPSCs.png")
build_manhattan(day0_lsb, "MeDIP-Seq_manhattan_LSB_vs_hPSCs.png")
build_manhattan(day0_d30, "MeDIP-Seq_manhattan_day30_vs_hPSCs.png")
build_manhattan(day0_d50, "MeDIP-Seq_manhattan_day50_vs_hPSCs.png")
build_manhattan(a1_lsb, "MeDIP-Seq_manhattan_LSB_vs_A1.png")
build_manhattan(a1_d30, "MeDIP-Seq_manhattan_day30_vs_A1.png")
build_manhattan(a1_d50, "MeDIP-Seq_manhattan_day50_vs_A1.png")
build_manhattan(d30_d50, "MeDIP-Seq_manhattan_day50_vs_day30.png")

rm(list = ls())
gc()
