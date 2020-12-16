#' Run Rsubread to obtain count data for re-running DESeq2 for
#' filtered MeDIP-Seq peaks - run this on linux ec2 for RAM > 16 cores.
#' Filter to retain peaks in at least two biological replicates
#' for MeDIP-Seq data from Active Motif.

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
fastq.files <- list.files(path = "./FASTQ/", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files <- fastq.files[-1]

# build ref hg38
buildindex(basename = "hg38", reference = "./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa")

# pull annotations from anotatr pkg
hg38 <- builtin_annotations()[grepl("hg38", builtin_annotations())]
# keep promoters, exons, introns, intergenic, cpg islands
# cpg shores, shelves, inter, and cpgs
sub_hg <- hg38[c(2, 5, 7, 11:15, 19)]
# build annotations for hg38 - include all
annot_hg38 <- build_annotations(genome = 'hg38', annotations = sub_hg)
saveRDS(annot_hg38, file = "./adj_data/anotatr_annotations.RDS")
annot_hg38 <- data.frame(annot_hg38) %>% 
  dplyr::select(c(gene_id, seqnames, start, end, strand))
names(annot_hg38) <- c("GeneID", "Chr", "Start", "End", "Strand")
annot_hg38 <- unique(annot_hg38)
saveRDS(annot_hg38, file = "./adj_data/annotatr_saf.RDS")

# align to the index created above - should specify threads!!!
# this step takes time, so run on a screen
out <- paste(fastq.files, "subread_results.bam", sep = "")
align(index = "./FASTQ/ref_hg38/hg38", readfile1 = fastq.files, output_file = out, nthreads = 6)

# load dir of bam files
bam <- list.files("./FASTQ/", pattern = ".BAM$", full.names = T)

# need to run this code in samtools on the ec2
# this takes the BAM files generated in the align-
# ment on line 40 and sorts and indexes them
# this is used for building plots in Gviz
# replace what's in < > with the name of bam file 
# samtools sort <.BAM$> -o <.sort.BAM$> -@ 6 # threads
# samtools align <.sort.BAM$>

# load saf file from anotatr
anot <- readRDS("./adj_data/annotatr_saf.RDS")
# swap order of genes on negative strand from
# end to start & start to end; then take abs
# value to feed as input into Rsubread
anot <- anot %>% 
  mutate(add = ifelse(Start < 0, End, Start),
         End = ifelse(End < 0, Start, End))
  mutate(Start = add) %>% 
  mutate(Start = abs(Start), 
         End = abs(End)) %>% 
  dplyr::select(-add)

# this step takes time, so run on a screen
# calculate counts using featureCounts, SAF annotation
# file from above
counts <- featureCounts(bam, annot.ext = anot, isGTFAnnotationFile = F,
                        useMetaFeatures = T)

# FROM HERE DOWN, NO LARGE RAM REQ'TS
# extract count data only
comb <- data.frame(counts$counts) %>%
  rownames_to_column("GeneID") %>%
  mutate(GeneID = gsub("^X", "", GeneID)) %>% 
  filter(GeneID != "NA.")
# extract annotation data
an <- data.frame(counts$annotation) %>%
  mutate(GeneID = as.character(GeneID)) %>%
  filter(GeneID != "NA")
# combine count and annotation data
comb <- comb %>% 
  left_join(., an, by = "GeneID") %>%
  distinct(.)

# filter out rows where counts are 0 across all samples
comb <- comb[rowSums(!as.matrix(comb[, c(2:13)])) < ncol(comb[, c(2:13)]), ]

# split out chr, start, end & strand into separate rows
comb <- rbindlist(lapply(1:nrow(comb), function(x) {
  y <- comb[x,]
  chr <- unlist(strsplit(y$Chr, ";"))
  st <- unlist(strsplit(y$Start, ";"))  
  end <- unlist(strsplit(y$End, ";"))
  str <- unlist(strsplit(y$Strand, ";"))
  z <- y[, c(1:13)][rep(seq_len(nrow(y)), each = length(chr)), ]
  z <- z %>%
    mutate(Chr = chr,
           Start = st,
           End = end,
           Strand = str) %>%
    mutate(Length = as.numeric(End) - as.numeric(Start))
}))
names(comb) <- gsub("^X|\\.medip.*", "", names(comb))

# load mergedregs file containing deseq2/peak info
# from active motif
df <- readWorkbook("./data/4075NIH_5meC_mergedregs.xlsx")
names(df)[1] <- "MergedRegion"

# retain only those MergedRegions that have at least
# 2 biological replicates supporting it
MRegs <- df %>%
  dplyr::select(MergedRegion, Chromosome, Start, End, ends_with("Present")) %>%
  gather(key = Sample, value = Value, -MergedRegion, -Chromosome, -Start, -End) %>%
  mutate(Sample = gsub("-p.*|-i.*|-r.*", "", Sample),
         Treatment =  ifelse(grepl("TE", Sample) & grepl("H9", Sample), "H9-TE",
                             ifelse(!grepl("TE", Sample) & grepl("H9", Sample), "H9-hPSC",
                                    ifelse(grepl("TE", Sample) & grepl("WA17", Sample), "WA17-TE", "WA17")))) %>%
  group_by(Treatment, MergedRegion) %>%
  mutate(PresSum = sum(Value),
         Remove = any(PresSum < 2)) %>% 
  ungroup() %>%
  filter(Remove == TRUE) %>%
  pull(MergedRegion) %>%
  unique()
write.csv(MRegs, file = "./adj_data/filtered_medipseq_peak.csv", row.names = F)

# retain only merged regions w/2>= bio reps
# supporting it
df %<>%
  dplyr::select(MergedRegion, Chromosome, Start, End, Length, 
                IntervalCount, CGIslandCount, PromoterCount, 
                GeneCount, Gene.List, Dist.to.Start, Position,
                UCSC.Link) %>%
  filter(MergedRegion %in% MRegs) %>%
  as.data.frame()

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
cts <- comb %>%
  dplyr::select(c(Chr, Start, End)) %>%
  mutate_at(c(2:3), as.numeric) %>%
  mutate(Chr = gsub("chr", "", Chr)) %>%
  arrange(Chr, Start)
cts_nam <- paste(cts$Chr, as.numeric(cts$Start), as.numeric(cts$End), sep = "_")
# this was producing a weird error downstream,
# this fixes it
cts_nam <- ifelse(grepl("15_41999917_", cts_nam), "15_41999917_42000000", cts_nam)
cts_nam <- ifelse(grepl("11_5e", cts_nam), "11_500000_500483", cts_nam)
cts_nam <- ifelse(grepl("2.3e", cts_nam), "8_22999679_23000000", cts_nam)
# put it in a genomicranges obj
cts <- GenomicRanges::GRanges(
  seqnames = cts$Chr,
  ranges = IRanges(start = cts$Start,
                   end = cts$End)
)
names(cts) <- cts_nam

# need to find overlap b/t count data and medip peaks
# annot == count; query == medipseq peaks
overlap <- GenomicRanges::findOverlaps(query, cts)

# find overlap b/t counts from featureCounts &
# medip seq peaks
overlap2 <- data.frame("chromosome_name_medip" = names(query)[queryHits(overlap)],
                       "chromosome_name_cts" = names(cts)[subjectHits(overlap)]) %>%
  separate(chromosome_name_medip, c("chromosome_name_medip", "start_position_medip", "end_position_medip")) %>%
  separate(chromosome_name_cts, c("chromosome_name_cts", "start_position_cts", "end_position_cts")) %>%
  distinct(.) %>%
  mutate(start_position_medip = as.numeric(start_position_medip),
         end_position_medip = as.numeric(end_position_medip),
         start_position_cts = as.numeric(start_position_cts),
         end_position_cts = as.numeric(end_position_cts)
  )
# rename cols in count file for
# easy merging
names(comb)[14:16] <- c("chromosome_name_cts", "start_position_cts", "end_position_cts")
comb <- comb %>%
  mutate_at(c(15:16), as.numeric) %>%
  mutate(chromosome_name_cts = gsub("chr", "", chromosome_name_cts))
# merge overlap b/t medipseq peaks and counts
# and count data
overlap2 <- overlap2 %>%
  left_join(., comb) %>%
  distinct(.)
overlap2 <- overlap2 %>% dplyr::select(-Length) %>% distinct(.)

# load anotatr annotation
anotatr_anot <- readRDS("./adj_data/anotatr_annotations.RDS")
anotatr_anot <- data.frame(anotatr_anot)
# rename cols for easy merging
names(overlap2)[c(4:6, 20)] <- c("seqnames", "start", "end", "strand")
overlap2$seqnames <- paste("chr", overlap2$seqnames, sep = "")

# combine annotation from anotatr with
# overlap b/t medip-seq + counts and
# count data
overlap2 <- overlap2 %>% left_join(., anotatr_anot)
overlap2 <- overlap2 %>%
  mutate(comb = paste(chromosome_name_medip, start_position_medip, end_position_medip, sep = "_"))
# split df into elements of a list for each medip-seq peak
overlap2 <- split(overlap2, overlap2$comb)
# combine annotation for a given medip-seq peak
# into a single row by catenating multiple annotations
# with the same chr, start, end site
overlap2 <- rbindlist(lapply(1:length(overlap2), function(x) {
  y <- overlap2[[x]]
  st <- paste(unique(y$strand), collapse = "; ")
  wid <- paste(unique(y$width), collapse = "; ")
  gene <- paste(unique(y$GeneID), collapse = "; ")
  id <- paste(unique(y$id), collapse = "; ")
  tx <- paste(unique(y$tx_id), collapse = "; ")
  sym <- paste(unique(y$symbol), collapse = "; ")
  type <- paste(unique(y$type), collapse = "; ")
  z <- data.frame(y[1, c(1:19)])
  z <- z %>%
    mutate(strand = st,
           width = wid,
           id = id,
           tx_id = tx,
           symbol = sym,
           type = type)
}))
write.csv(overlap2, file = "./adj_data/final_peaks_annotatr.csv", row.names = F)

# NOTES: PHLDA2, H19, PEG10 were all found in regions
# upstream or downstream of these genes not in exonic
# promoter, etc. regions; PEG3 only had 1 replicate w
# a peak, so it's excluded in the analysis

rm(list = ls())
gc()
