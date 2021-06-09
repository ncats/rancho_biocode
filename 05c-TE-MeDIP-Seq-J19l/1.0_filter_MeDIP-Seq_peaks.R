#' Filter to retain peaks in at least two biological replicates
#' for MeDIP-Seq data from Active Motif.

# req'd pkgs
x <- c("openxlsx", "magrittr", "tidyverse",
       "data.table", "GenomicRanges", "IRanges")
sapply(x, library, character.only = TRUE)

# load mergedregs file containing deseq2/peak info
# from active motif 
df <- readWorkbook("./data/4441NIH_J91i-MeDIP_mergedregs.xlsx")
names(df)[1] <- "MergedRegion"

# retain only those MergedRegions that have at least
# 2 biological replicates supporting it
MRegs <- df %>%
  dplyr::select(MergedRegion, Chromosome, Start, End, ends_with("Present")) %>%
  gather(key = Sample, value = Value, -MergedRegion, -Chromosome, -Start, -End) %>%
  mutate(Treatment =  ifelse(grepl("iPSC", Sample), "J91i-iPSC",
                             ifelse(grepl("TEp0", Sample), "J91i-TEp0",
                                    ifelse(grepl("TEexp", Sample), "TEexp", "NA")))) %>%
  group_by(Treatment, MergedRegion) %>%
  mutate(PresSum = sum(Value),
         Remove = any(PresSum < 2)) %>% 
  ungroup() %>%
  filter(Remove == TRUE) %>%
  pull(MergedRegion) %>%
  unique()
dir.create('./adj_data/', showWarnings = FALSE)
write.csv(MRegs, file = "./adj_data/filtered_medipseq_peaks.csv", row.names = F)

# retain only merged regions w/2>= bio reps
# supporting it
df %<>%
  dplyr::select(MergedRegion, Chromosome, Start, End, Length, 
                IntervalCount, CGIslandCount, PromoterCount, 
                GeneCount, Gene.List, Dist.to.Start, Position,
                UCSC.Link, ends_with("Counts")) %>%
  filter(MergedRegion %in% MRegs) %>%
  as.data.frame()
names(df) <- gsub("_MeDIP_.*", "", names(df))
saveRDS(df, file = "./adj_data/medipseq_counts.RDS")

rm(list = ls())
gc()
