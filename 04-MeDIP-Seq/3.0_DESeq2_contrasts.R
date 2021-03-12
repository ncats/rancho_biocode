#' Re-run DESeq-2 on filtered MeDIP-Seq peaks
#' found in at least 2 replicates

# req'd pkgs
x <- c("ggplot2", "tidyverse", "DESeq2", 
       "data.table", "tidyverse", 
       "GenomicRanges", "IRanges", "extrafont",
       "extrafontdb")
sapply(x, library, character.only = TRUE)

loadfonts()

source("./functions/adj_pca.R")

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20),
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))


# load medip/annotation data
data <- fread("./adj_data/final_peaks_annotatr.csv")
counts <- readRDS("./adj_data/medipseq_counts.RDS") %>%
  dplyr::select(-c(Chromosome, Start, End))
# merge w/count data
data <- data %>% left_join(., counts, by = "MergedRegion")

# subset to keep only count data
counts <- data.frame(data[, c(15:28)])
nam <- data$MergedRegion
names(counts) <- gsub("^X", "", names(counts))
rownames(counts) <- nam

cols <- names(counts)
# rename contrasts for ease downstream
coldata <- data.frame(Treatment = ifelse(grepl("NCRM5.D30", cols), "D30",
                                         ifelse(grepl("NCRM5.D50", cols), "D50",
                                                ifelse(grepl("A1", cols), "A1",
                                                       ifelse(grepl("LSB", cols), "LSB",
                                                              ifelse(grepl("NCRM5$", cols), "D0", "NA"))))))
rownames(coldata) <- colnames(counts)
counts <- as.matrix(counts)

# set up deseq obj
dds <- DESeqDataSetFromMatrix(counts, coldata, ~Treatment)
rld <- vst(dds, blind = T)
p <- adj_pca(rld, intgroup = "Treatment")
p <- p +
  geom_point(size = 2 ) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  labs(x = "PC1 (95%)",
       y = "PC2 (2%)") +
  theme_classic() +
  mytheme
p$labels$colour <- "Treatment"

ggsave(filename = "./adj_data/plots/pca/pca_all.png",
       plot =  p, 
       height = 5, 
       width = 5, 
       units = "in")

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
  if (grepl("day0", nam)) {
    nam <- gsub("hPSCs_day0", "NCRM5$", nam)
    cols <- names(res)[grepl(nam, names(res))]
  } else {
    cols <- names(res)[grepl(nam, names(res))]
  }
  
  res <- res %>%
    dplyr::select(c(1:20, cols))
  write.csv(res, file = paste("./adj_data/deseq/", con1, "_", con2, "_deseq2.csv", sep = ""), row.names = F)
  
  
  resLFC <- data.frame(lfcShrink(dds, contrast = c("treatment", con1, con2), type = "normal"))
  resLFC <- resLFC %>%
    rownames_to_column("MergedRegion") %>%
    mutate(MergedRegion = as.numeric(MergedRegion)) %>%
    left_join(., data) %>%
    distinct(.)
  resLFC <- resLFC %>%
    dplyr::select(c(1:20, cols))
  write.csv(res, file = paste("./adj_data/deseq_shrink/", con1, "_", con2, "_deseq2_shrink.csv", sep = ""), row.names = F)
  
}

contrasts_deseq("hPSCs_day0", "LSB")
contrasts_deseq("hPSCs_day0", "A1")
contrasts_deseq("hPSCs_day0", "NCRM5-D30")
contrasts_deseq("hPSCs_day0", "NCRM5-D50")
contrasts_deseq("A1", "LSB")
contrasts_deseq("A1", "NCRM5-D30")
contrasts_deseq("A1", "NCRM5-D50")
contrasts_deseq("NCRM5-D30", "NCRM5-D50")

rm(list = ls())
gc()
