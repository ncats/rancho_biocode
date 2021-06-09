#' Re-run DESeq-2 on filtered MeDIP-Seq peaks
#' found in at least 2 replicates

# req'd pkgs
x <- c("ggplot2", "tidyverse", "DESeq2", 
       "data.table", "tidyverse", 
       "GenomicRanges", "IRanges", "Cairo")
sapply(x, library, character.only = TRUE)

source("./functions/adj_pca.R")

# load medip/annotation data
data <- readRDS("./adj_data/medipseq_counts.RDS")

# subset to keep only count data
cols <- names(data)[grepl("J91i", names(data))]
counts <- data.frame(data[, cols])
nam <- data$MergedRegion
names(counts) <- gsub("^X", "", names(counts))
rownames(counts) <- nam

cols <- names(counts)
# rename contrasts for ease downstream
coldata <- data.frame(treatment = ifelse(grepl("iPSC", cols), "J91i-iPSC",
                                         ifelse(grepl("TEp0", cols), "J91i-TEp0",
                                                ifelse(grepl("TEexp", cols), "TEexp", "NA"))))
rownames(coldata) <- colnames(counts)
counts <- as.matrix(counts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(counts, coldata, ~treatment)
varstabilised <- vst(dds, blind = T)
dir.create('./adj_data/plots/', showWarnings = FALSE)
dir.create('./adj_data/plots/pca/', showWarnings = FALSE)
pdf(file = "./adj_data/plots/pca/pca_all.pdf", width = 5, height = 5)
print(adj_pca(varstabilised, intgroup = "treatment") + theme_minimal())
graphics.off()
# run DESeq2
dds <- DESeq(dds)

# pull the correct contrasts for deseq
# and deseq shrink data
dir.create('./adj_data/deseq/', showWarnings = FALSE)
dir.create('./adj_data/deseq_shrink/', showWarnings = FALSE)
contrasts_deseq <- function(con1, con2) {
  res <- data.frame(results(dds, contrast = c("treatment", con1, con2))) %>%
    rownames_to_column("MergedRegion") %>%
    mutate(MergedRegion = as.numeric(MergedRegion)) %>%
    left_join(., data) %>%
    distinct(.)
  
  nam <- paste(con1, con2, sep = "|")
  nam <- gsub("\\-", ".", nam)
  cols <- names(res)[grepl(nam, names(res))]
  
  res <- res %>%
    dplyr::select(c(1:21, cols))
  write.csv(res, file = paste("./adj_data/deseq/", con1, "_", con2, "_deseq2.csv", sep = ""), row.names = F)
  
  
  resLFC <- data.frame(lfcShrink(dds, contrast = c("treatment", con1, con2), type = "normal"))
  resLFC <- resLFC %>%
    rownames_to_column("MergedRegion") %>%
    mutate(MergedRegion = as.numeric(MergedRegion)) %>%
    left_join(., data) %>%
    distinct(.)
  resLFC <- resLFC %>%
    dplyr::select(c(1:21, cols))
  write.csv(res, file = paste("./adj_data/deseq_shrink/", con1, "_", con2, "_deseq2_shrink.csv", sep = ""), row.names = F)
  
}

contrasts_deseq("J91i-iPSC", "TEexp")

rm(list = ls())
gc()
