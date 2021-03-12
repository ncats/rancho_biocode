#' Filter to retain peaks in at least two biological replicates
#' for MeDIP-Seq data from Active Motif.

# req'd pkgs
x <- c("openxlsx", "magrittr", "tidyverse",
       "data.table", "annotatr", "GenomicRanges", 
       "IRanges")
sapply(x, library, character.only = TRUE)

# load count files based on medipseq peaks
#mr_comb <- readRDS("./adj_data/mr_counts.RDS") 

# create anotator annotation
# pull annotations from anotatr pkg
hg38 <- builtin_annotations()[grepl("hg38", builtin_annotations())]
# keep promoters, exons, introns, intergenic, cpg islands
# cpg shores, shelves, inter, and cpgs
sub_hg <- hg38[c(2, 4:5, 7, 10:15, 19)]
# build annotations for hg38 - include all
annot_hg38 <- build_annotations(genome = 'hg38', annotations = sub_hg)
saveRDS(annot_hg38, file = "./adj_data/anotator_hg38.RDS")

# # fxn to extract count, annotation data & merge
# # them into a df
# count_annot_merge <- function(file) {
#   # extract count data only
#   comb <- data.frame(file$counts) %>%
#     rownames_to_column("GeneID") %>%
#     dplyr::select(-c(2)) %>%
#     mutate(GeneID = gsub("^X", "", GeneID)) %>% 
#     filter(GeneID != "NA.")
#   # extract annotation data
#   an <- data.frame(file$annotation) %>%
#     mutate(GeneID = as.character(GeneID)) %>%
#     filter(GeneID != "NA")
#   # combine count and annotation data
#   comb <- comb %>% 
#     left_join(., an, by = "GeneID") %>%
#     distinct(.)
# }
# 
# # merge annotation w/count data
# mr_comb <- count_annot_merge(mr_comb)
# mr_comb <- mr_comb %>%
#   dplyr::select(-c(14:15))
# names(mr_comb) <- gsub("^X|4441NIH\\.|\\.[0-9]{1}\\.MeDIP.*", "", names(mr_comb))
# names(mr_comb)[c(1, 11)] <- c("MergedRegion", "Chromosome")
# mr_comb$MergedRegion <- as.numeric(mr_comb$MergedRegion)
# 
# # filter out rows where counts are 0 across all samples
# mr_comb <- mr_comb[rowSums(!as.matrix(mr_comb[, c(2:10)])) < ncol(mr_comb[, c(2:10)]), ]
# saveRDS(mr_comb, file = "./adj_data/medipseq_counts.RDS")

# load mergedregs file containing deseq2/peak info
# from active motif
df <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx")
names(df)[1] <- "MergedRegion"

# retain only those MergedRegions that have at least
# 2 biological replicates supporting it
MRegs <- df %>%
  dplyr::select(MergedRegion, Chromosome, Start, End, ends_with("Present")) %>%
  gather(key = Sample, value = Value, -MergedRegion, -Chromosome, -Start, -End) %>%
  mutate(Treatment =  ifelse(grepl("iPSC", Sample), "J98i-iPSC",
                             ifelse(grepl("TEp0", Sample), "J98i-TEp0",
                                    ifelse(grepl("TEexp", Sample), "TEexp", "NA")))) %>%
  group_by(Treatment, MergedRegion) %>%
  mutate(PresSum = sum(Value),
         Remove = any(PresSum < 2)) %>% 
  ungroup() %>%
  filter(Remove == TRUE) %>%
  pull(MergedRegion) %>%
  unique()
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
# # merge peak data w/count data
# df <- df %>%
#   left_join(., mr_comb) %>%
#   filter(!is.na(`01.NCRM5`))

# create a df w/chr, st, end site from filtered
# medip-seq peaks
query <- df %>%
  dplyr::select(c(Chromosome, Start, End)) %>%
  arrange(Chromosome, Start)
query_nam <- paste(query$Chromosome, query$Start, query$End, sep = "_")
# put it into a genomicranges obj
query <- GenomicRanges::GRanges(
  seqnames = query$Chromosome,
  ranges = IRanges(start = query$Start,
                   end = query$End)
)
names(query) <- query_nam

# extract chr, st, end site from count data
anot <- data.frame(annot_hg38) %>%
  dplyr::select(c(seqnames, start, end)) %>%
  mutate_at(c(2:3), as.numeric) %>%
  mutate(seqnames = gsub("chr", "", seqnames)) %>%
  arrange(seqnames, start) %>%
  filter(!(grepl("_", seqnames)))
anot_nam <- paste(anot$seqnames, as.numeric(anot$start), as.numeric(anot$end), sep = "_")
# put it in a genomicranges obj
anot <- GenomicRanges::GRanges(
  seqnames = anot$seqnames,
  ranges = IRanges(start = anot$start,
                   end = anot$end)
)
names(anot) <- anot_nam

# need to find overlap b/t count data and medip peaks
# annot == count; query == medipseq peaks
overlap <- GenomicRanges::findOverlaps(query, anot)

# find overlap b/t counts from medipseq peaks/counts &
# anotatr annotation
overlap2 <- data.frame("chromosome_name_medip" = names(query)[queryHits(overlap)],
                       "chromosome_name_anot" = names(anot)[subjectHits(overlap)]) %>%
  # fix e+ issue in R in a few rows
  mutate(chromosome_name_anot = gsub("5e\\+05", "500000", chromosome_name_anot)) %>%
  mutate(chromosome_name_anot = gsub("4\\.2e\\+07$", "42000000", chromosome_name_anot)) %>%
  separate(chromosome_name_medip, c("chromosome_name_medip", "start_position_medip", "end_position_medip")) %>%
  separate(chromosome_name_anot, c("chromosome_name_anot", "start_position_anot", "end_position_anot")) %>%
  distinct(.) %>%
  mutate(start_position_medip = as.numeric(start_position_medip),
         end_position_medip = as.numeric(end_position_medip),
         start_position_anot = as.numeric(start_position_anot),
         end_position_anot = as.numeric(end_position_anot)
  )
# rename cols in count file for
# easy merging
annot_hg38 <- data.frame(annot_hg38)
names(annot_hg38)[1:3] <- c("chromosome_name_anot", "start_position_anot", "end_position_anot")
annot_hg38 <- annot_hg38 %>%
  mutate(chromosome_name_anot = gsub("chr", "", chromosome_name_anot))
# merge overlap b/t medipseq peaks and counts
# and count data
overlap2 <- overlap2 %>%
  left_join(., annot_hg38) %>%
  distinct(.)

df_sub <- df[, c(1:5, 10:12)]
names(df_sub)[c(2:4)] <- c("chromosome_name_medip", "start_position_medip", "end_position_medip")
overlap2 <- overlap2 %>%
  left_join(., df_sub) %>%
  distinct(.)
overlap2 <- overlap2[, c(14, 1:3, 8:13, 15:18)]

overlap2 <- overlap2 %>%
  mutate(comb = paste(MergedRegion, chromosome_name_medip, start_position_medip, end_position_medip, sep = "_"))
# split df into elements of a list for each medip-seq peak
overlap2 <- split(overlap2, overlap2$comb)
# combine annotation for a given medip-seq peak
# into a single row by catenating multiple annotations
# with the same chr, start, end site
overlap2 <- rbindlist(lapply(1:length(overlap2), function(x) {
  y <- overlap2[[x]]
  
  # these are from anotatr
  st <- paste(y$strand, collapse = "; ")
  id <- paste(y$id, collapse = "; ")
  tx <- paste(y$tx_id, collapse = "; ")
  gene <- paste(y$gene_id, collapse = "; ")
  sym <- paste(y$symbol, collapse = "; ")
  type <- paste(y$type, collapse = "; ")
  len <- paste(y$Length, collapse = "; ")
  
  # these are from active motif
  active_gene <- paste(unlist(strsplit(y$Gene.List, ", ")), collapse = "; ")
  dis <- paste(unlist(strsplit(y$Dist.to.Start, ", ")), collapse = "; ")
  pos <- paste(unlist(strsplit(y$Position, ", ")), collapse = "; ")

  z <- data.frame(y[1, c(1:4)])
  z <- z %>%
    mutate(strand = st,
           id = id,
           tx_id = tx,
           gene = gene,
           gene_anotatr = sym,
           type = type,
           width = len,
           gene_active = active_gene,
           location = pos,
           distance_to_start = dis
           )
}))
write.csv(overlap2, file = "./adj_data/final_peaks_annotatr.csv", row.names = F)

rm(list = ls())
gc()
