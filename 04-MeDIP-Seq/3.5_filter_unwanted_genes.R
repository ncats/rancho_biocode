#' Filter out genes that should not be included
#' in analysis downstream (CM's list)

# THIS CODE TAKES A LONG TIME TO RUN, SO YOU 
# SHOULD RUN IT ON THE EC2 AND SPECIFY 7 CORES;
# IT STILL TAKES 10 MIN OR SO TO RUN FOR EACH
# LIST

# req'd pkgs
x <- c("tidyverse", "data.table", "parallel")
sapply(x, library, character.only = TRUE)

# load genes we want to filter out - list
# from CM
# load genes to filter out by CM
fg <- fread("./data/Genes_To_Remove_filter.txt", header = T)
fg <- fg[-1,]

# function to split up gene symbols from ActiveMotif
# and anotatr and filter out genes from list Claire 
# provided
filter_genes <- function(deseq) {
  y <- split(deseq, deseq$MergedRegion)
  
  # filter out genes from ActiveMotif annotation
  y <- rbindlist(lapply(1:length(y), function(x) {
    y <- y[[x]]
    nam <- names(y)
    sym <- unlist(strsplit(y$symbol, "; "))
    y <- y %>% dplyr::select(-symbol)
    y <- y[rep(seq_len(nrow(y)), each = length(sym)), ]
    y$symbol <- sym
    y <- y %>%
      filter(!symbol %in% fg$Gene) %>%
      data.frame(.)
    names(y) <- gsub("^X", "", names(y))
    y <- y[, nam]

    if (nrow(y) > 1) {
      sym <- paste(y$symbol, collapse = "; ")
      y <- y[1,]
      y$symbol <- sym
      return(y)
    } else {
      return(y)
    }
  }))
  
  z <- split(y, y$MergedRegion)
  
  # filter out genes from anotatr annotation
  z <- rbindlist(lapply(1:length(z), function(x) {
    z <- z[[x]]
    nam <- names(z)
    sym <- unlist(strsplit(z$gene_anotatr, ", "))
    z <- z %>% dplyr::select(-gene_anotatr)
    z <- z[rep(seq_len(nrow(z)), each = length(sym)), ]
    z$gene_anotatr <- sym
    z <- z %>%
      filter(!gene_anotatr %in% fg$Gene) %>%
      data.frame(.)
    names(z) <- gsub("^X", "", names(z))
    z <- z[, nam]
    
    if (nrow(z) > 1) {
      sym <- paste(z$gene_anotatr, collapse = "; ")
      z <- z[1,]
      z$gene_anotatr <- sym
      return(z)
    } else {
      return(z)
    }
  }))
  return(z)
}

# load deseq & shrinkage data
deseq_files <- list.files("./adj_data/deseq/", pattern = "deseq2.csv$", full.names = T)
deseq_nam <- gsub(".*\\/|_deseq.*", "", deseq_files)
deseq <- lapply(deseq_files, function(x) fread(x))
shrink_files <- list.files("./adj_data/deseq_shrink/", pattern = "deseq2_shrink.csv$", full.names = T)
shrink_nam <- gsub(".*\\/|_deseq.*", "", shrink_files)
shrink <- lapply(shrink_files, function(x) fread(x))

deseq <- mclapply(1:length(deseq), function(x) {
  y <- filter_genes(deseq[[x]])
  y <- y %>% mutate(contrast = deseq_nam[x])
}, mc.cores = 7)
lapply(1:length(deseq), function(x) 
  write.csv(deseq[[x]], file = paste("./adj_data/deseq/filtered/", deseq_nam[x], "_deseq_filtered.csv", sep = ""), row.names = F))

shrink <- mclapply(1:length(shrink), function(x) {
  y <- filter_genes(shrink[[x]])
  y <- y %>% mutate(contrast = shrink_nam[x])
}, mc.cores = 7)

lapply(1:length(shrink), function(x) 
  write.csv(shrink[[x]], file = paste("./adj_data/deseq_shrink/filtered/", shrink_nam[x], "_shrink_filtered.csv", sep = ""), row.names = F))

rm(list = ls())
gc()