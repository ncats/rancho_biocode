#' Filter out genes that should not be included
#' in analysis downstream (CM's list)

# req'd pkgs
x <- c("tidyverse", "data.table", "parallel")
sapply(x, library, character.only = TRUE)

source("./functions/filter_genes.R")

# load genes to filter out by CM
fg <- fread("./data/Genes_To_Remove_filter.txt", header = T)
fg <- fg[-1,]

# deseq files
deseq_files <- list.files("./adj_data/deseq/", pattern = "deseq2.csv$", full.names = T)
deseq_nam <- gsub(".*\\/|_deseq.*", "", deseq_files)

# filter out genes in deseq
dir.create('./adj_data/deseq/filtered/', showWarnings = FALSE)
lapply(1:length(deseq_files), function(x) {
 filter_genes(deseq_files[x], "./adj_data/deseq/filtered/", deseq_nam[x], "_deseq") 
})

# shrinkage files
shrink_files <- list.files("./adj_data/deseq_shrink/", pattern = "deseq2_shrink.csv$", full.names = T)
shrink_nam <- gsub(".*\\/|_deseq.*", "", shrink_files)

# filter out genes in deseq_shrink
dir.create('./adj_data/deseq_shrink/filtered/', showWarnings = FALSE)
lapply(1:length(shrink_files), function(x) {
  filter_genes(shrink_files[[x]], "./adj_data/deseq_shrink/filtered/", deseq_nam[x], "_deseq_shrink") 
})

rm(list = ls())
gc()
