#' Filter bed file to retain only peaks that 
#' were found in at least 2 bio reps

# req'd pkgs
x <- c("data.table", "dplyr")
sapply(x, library, character.only = TRUE)

# pull merged peaks from .bed file
t <- fread("./data/4075NIH_5meC_ucsctracks_bw_noheader.bed", fill = T)
t <- t %>% filter(!grepl("track", t$V1))
t <- t[-c(177877:nrow(t))]

# load wa17 or h9 deseq contrasts to pull merged regions
mreg <- fread("./adj_data/filtered_medipseq_peak.csv") %>%
  mutate(x = paste("MergedReg_", x, sep = ""))

# filter it to remove those not represented by
# at least two bio reps
t <- t %>%
  filter(V4 %in% mreg$x)
write.table(t, file = "./adj_data/MergedPeaks.bed", sep = "\t", row.names = F, col.names = F, quote = F)

rm(list = ls())
gc()
