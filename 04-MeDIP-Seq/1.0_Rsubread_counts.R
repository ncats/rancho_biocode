#' Run Rsubread to obtain count data for re-running DESeq2 for
#' filtered MeDIP-Seq peaks - run this on linux ec2 for RAM > 16 cores.

# req'd pkgs
# NOTE: To run Rsubread for indexing of genome,
# alignment and counting, you need to run on a
# machine w/ > 16 GB RAM
x <- c("openxlsx", "magrittr", "tidyverse",
       "data.table", "annotatr", "GenomicRanges", 
       "IRanges", "Rsubread", "org.Hs.eg.db")
sapply(x, library, character.only = TRUE)

# file paths are related to ec2: ec2-user@ec2-3-84-14-241.compute-1.amazonaws.com
# load fastq files
fastq.files <- list.files(path = "./FASTQ_MeDIP2", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files <- fastq.files[-1]

nam <- gsub(".*\\/", "", fastq.files)
out <- paste("./BAM/", nam, "_subread_results.bam", sep = "")

# load dir of bam files
bam <- list.files("./BAM/", pattern = ".bam$", full.names = T)

# create a saf file based on mergedregions
# identified by active motif
mr <- readWorkbook("./data/4075NIH_Lonza_MeDIP_mergedregs.xlsx")
saf_mr <- mr[, c(1:4)] %>%
    mutate(Strand = "*")
names(saf_mr) <- c("GeneID", "Chr", "Start", "End", "Strand")
saveRDS(saf_mr, file = "./adj_data/mr_saf.RDS")

# this step takes time, so run on a screen
# calculate counts using featureCounts, SAF annotation
# file from above
mr_counts <- featureCounts(bam, annot.ext = saf_mr, isGTFAnnotationFile = F, useMetaFeatures = T)
saveRDS(mr_counts, file = "./adj_data/mr_counts.RDS")

rm(list = ls())
gc()
