#' Create manhattan plots across all chromosomes 
#' for each contrast

# req'd pkgs
x <- c("ggplot2", "tidyverse", "qqman", "data.table",
       "ggrepel", "gtools")
sapply(x, library, character.only = TRUE)

source("./functions/build_manhattan.R")

# load shrinkage data
ipsc_tep0 <- fread("./adj_data/deseq_shrink/filtered/J98i-iPSC_J98i-TEp0_shrink_filtered.csv")
ipsc_teexp <- fread("./adj_data/deseq_shrink/filtered/J98i-iPSC_TEexp_shrink_filtered.csv")
tep0_teexp <- fread("./adj_data/deseq_shrink/filtered/J98i-TEp0_TEexp_shrink_filtered.csv")


# build manhattan plots
build_manhattan(ipsc_tep0, "Manhattan_TEp0_vs_iPSC.png")
build_manhattan(ipsc_teexp, "Manhattan_TEexp_vs_iPSC.png")
build_manhattan(tep0_teexp, "Manhattan_TEexp_vs_TEp0.png")

rm(list = ls())
gc()
