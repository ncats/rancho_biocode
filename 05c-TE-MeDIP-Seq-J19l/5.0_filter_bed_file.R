#' Filter bed file to retain only peaks that 
#' were found in at least 2 bio reps

# req'd pkgs
x <- c("data.table", "dplyr")
sapply(x, library, character.only = TRUE)

# pull merged peaks from .bed file
t <- fread("./data/4441NIH_J91i-MeDIP_ucsctracks_bw.bed", skip = 3, fill = T)
t <- t %>% filter(!grepl("track", t$V1))
t <- t[-c(146254:nrow(t))]

# deseq contrasts to pull merged regions
mreg <- readRDS("./adj_data/medipseq_counts.RDS")
mreg$x <- paste0("MergedReg_", mreg$MergedRegion)

# filter it to remove those not represented by
# at least two bio reps
t <- t %>%
  filter(V4 %in% mreg$x)
write.table(t, file = "./adj_data/MergedPeaks.bed", sep = "\t", row.names = F, col.names = F, quote = F)

rm(list = ls())
gc()
