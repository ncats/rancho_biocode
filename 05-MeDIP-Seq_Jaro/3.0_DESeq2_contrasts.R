#' Re-run DESeq-2 on filtered MeDIP-Seq peaks
#' found in at least 2 replicates

# req'd pkgs
x <- c("ggplot2", "tidyverse", "DESeq2", 
       "data.table", "tidyverse", 
       "GenomicRanges", "IRanges", "Cairo")
sapply(x, library, character.only = TRUE)

source("./functions/adj_pca.R")

# load medip/annotation data
data <- fread("./adj_data/final_peaks_annotatr.csv")
counts <- readRDS("./adj_data/medipseq_counts.RDS") %>%
  dplyr::select(-c(Chromosome, Start, End))
# merge w/count data
data <- data %>% left_join(., counts, by = "MergedRegion")

# subset to keep only count data
counts <- data.frame(data[, c(24:32)])
nam <- data$MergedRegion
names(counts) <- gsub("^X", "", names(counts))
rownames(counts) <- nam

cols <- names(counts)
# rename contrasts for ease downstream
coldata <- data.frame(treatment = ifelse(grepl("iPSC", cols), "J98i-iPSC",
                                         ifelse(grepl("TEp0", cols), "J98i-TEp0",
                                                ifelse(grepl("TEexp", cols), "TEexp", "NA"))))
rownames(coldata) <- colnames(counts)
counts <- as.matrix(counts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(counts, coldata, ~treatment)
rld <- vst(dds, blind = T)
CairoPDF(file = "./adj_data/plots/pca/pca_all.png", width = 5, height = 5, family = "NotoSans-Condensed")
print(adj_pca(rld, intgroup = "treatment") + theme_minimal())
graphics.off()
# run DESeq2
dds <- DESeq(dds)

# pull the correct contrasts for deseq
# and deseq shrink data
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

contrasts_deseq("J98i-iPSC", "J98i-TEp0")
contrasts_deseq("J98i-iPSC", "TEexp")
contrasts_deseq("J98i-TEp0", "TEexp")

rm(list = ls())
gc()
