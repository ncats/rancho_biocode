#' Combine all DESeq2 results into a single
#' excel file.

# req'd libs
x <- c("magrittr", "tidyverse", "openxlsx", "data.table")
sapply(x, library, character.only = T)

source("./functions/remove_redundancy.R")

# load deseq files
all <- list.files("./adj_data/deseq/filtered/", full.names = T)

# remove unncessary items from naming contrasts
all_nam <- gsub(".*\\/|_deseq.*", "", all)

# roll up medip-seq peak filtered, 10 kb up/down
# filtered data into a single excel file
# create an empty excel workbook
xlwb <- createWorkbook()

# remove redundancies in each row (i.e. location, distance, etc)
lapply(1:length(all), function(x) remove_redundancy(all[[x]], all_nam[x]))

#save excel workbook with rolled-up, non-redundant deseq results
saveWorkbook(xlwb, paste("./adj_data/final/", "deseq_filtered_all_contrasts.xlsx", sep = ""), overwrite = T)

# create an empty excel workbook
xlwb <- createWorkbook()

# load filtered, 10 kb removed deseq shrink contrasts
all <- list.files("./adj_data/deseq_shrink/filtered/", full.names = T)
lapply(1:length(all), function(x) remove_redundancy(all[[x]], all_nam[x]))

#save excel workbook with rolled-up, non-redundant deseq shrink results
saveWorkbook(xlwb, paste("./adj_data/final/", "deseq_shrink_filtered_all_contrasts.xlsx", sep = ""), overwrite = T)

rm(list = ls())
gc()
