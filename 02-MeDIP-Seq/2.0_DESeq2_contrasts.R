#' Re-run DESeq-2 on filtered MeDIP-Seq peaks
#' found in at least 2 replicates

# req'd pkgs
x <- c("ggplot2", "tidyverse", "DESeq2", 
       "data.table", "tidyverse", 
       "GenomicRanges", "IRanges")
sapply(x, library, character.only = TRUE)

# load count files per exon from rsubread
# after overlap w/medipseq & anotatr
data <- fread("./adj_data/final_peaks_annotatr.csv")
data$GeneID <- as.character(data$GeneID)

# multiple rows for annotation/counts, so 
# roll me-dip peaks into one row per entrez id
# as that's what is used in DESeq2
data <- split(data, data$GeneID)
data <- rbindlist(lapply(1:length(data), function(x) {
  y <- data[[x]]
  if (nrow(y) == 1) {
    return(y)
  } else {
    chr <- paste(unique(y$chromosome_name_medip), collapse = "; ")
    st <- paste(unique(y$start_position_medip), collapse = "; ")
    ed <- paste(unique(y$end_position_medip), collapse = "; ")
    seq <- paste(y$seqnames, collapse = "; ")
    st2 <- paste(y$start, collapse = "; ")
    ed2 <- paste(y$end, collapse = "; ")
    symb <- if (grepl("; ", paste(y$symbol, collapse = ""))) {
              a <- ifelse(grepl("; ", y$symbol), gsub("; ", "-", y$symbol), y$symbol)
              a <- paste(a, collapse = "; ")
            } else {
              a <- paste(y$symbol, collapse = "; ")
            }
    z <- data.frame(y[1, c(7:23, 25)])
    z <- z %>%
      mutate(chromosome_name_medip = chr,
             start_position_medip = st,
             end_position_medip = ed,
             seqnames = seq,
             start = st2,
             end = ed2,
             symbol = symb)
    z <- z[, c(19:24, 1:17, 25, 18)]
  }
}))

# subset to keep only count data
counts <- data.frame(data[, c(7:19)])
nam <- counts$GeneID
names(counts) <- gsub("^X", "", names(counts))
counts$GeneID <- NULL
rownames(counts) <- nam

cols <- names(counts)
# rename contrasts for ease downstream
coldata <- data.frame(treatment = ifelse(grepl("TE", cols) & grepl("H9", cols), "H9_TE",
                                         ifelse(!grepl("TE", cols) & grepl("H9", cols), "H9_hPSC",
                                                ifelse(grepl("TE", cols) & grepl("WA17", cols), "WA17_TE", "WA17"))))
rownames(coldata) <- colnames(counts)
counts <- as.matrix(counts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(counts, coldata, ~treatment)
# run DESeq2
dds <- DESeq(dds)

# pull the correct contrasts
wa17_res <- data.frame(results(dds, contrast = c("treatment", "WA17_TE", "WA17"))) %>%
  rownames_to_column("GeneID") %>%
  left_join(., data) %>%
  distinct(.)
h9_res <- data.frame(results(dds, contrast = c("treatment", "H9_TE", "H9_hPSC"))) %>%
  rownames_to_column("GeneID") %>%
  left_join(., data) %>%
  distinct(.)

# calc shrinkage log2fc values for plotting
resLFC_h9 <- data.frame(lfcShrink(dds, contrast = c("treatment", "H9_TE", "H9_hPSC"), type = "normal"))
resLFC_h9 <- resLFC_h9 %>%
  rownames_to_column("GeneID") %>%
  left_join(., data) %>%
  distinct(.)
resLFC_wa17 <- data.frame(lfcShrink(dds, contrast = c("treatment", "WA17_TE", "WA17"), type = "normal"))
resLFC_wa17 <- resLFC_wa17 %>%
  rownames_to_column("GeneID") %>%
  left_join(., data) %>%
  distinct(.)

# save output
write.csv(wa17_res, file = "./adj_data/wa17_deseq2.csv", row.names = F)
write.csv(h9_res, file = "./adj_data/h9_deseq2.csv", row.names = F)
write.csv(resLFC_wa17, file = "./adj_data/wa17_deseq2_shrink.csv", row.names = F)
write.csv(resLFC_h9, file = "./adj_data/h9_deseq2_shrink.csv", row.names = F)

rm(list = ls())
gc()
