#' Global DNA methylation plots

# req'd pkgs
x <- c("tidyverse", "data.table", "openxlsx", 
       "showtext", "ggplot2", "scales", "GenomicRanges",
       "Rsubread", "bamsignals", "RColorBrewer", "Cairo",
       "extrafont", "extrafontdb")
sapply(x, library, character.only = TRUE)

# so there's no scientific notation
options(scipen=999)

# add correct font
font_add("noto_cond", "./data/Noto-Sans-Condensed/NotoSans-Condensed.ttf")
font_add("noto_bold", "./data/Noto-Sans-Condensed/NotoSans-CondensedBold.ttf")
showtext_auto()
loadfonts()

# set cairo text backend for cairo_pdf
CairoFonts(
  regular="NotoSans-Condensed:style=Regular",
  bold="NotoSans-Condensed:style=Bold",
  italic="NotoSans-Condensed:style=Italic",
  bolditalic="NotoSans-Condensed:style=Bold Italic, BoldItalic",
  symbol="Symbol"
)

# load custom fxns
source("./functions/global_dna.R")

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "noto_bold"), 
                 axis.text = element_text(size = 14, family = "noto_cond"),
                 axis.title = element_text(colour = "Black", size = 16, family = "noto_bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "noto_cond"),
                 legend.title = element_text(colour = "Black", size = 14, family = "noto_cond"))

##### CREATE GLOBAL METHYLATION PLOTS #####

# j98i
orig <- readWorkbook("./data/4441NIH_J98i-MeDIP_mergedregs.xlsx") 
cols <- names(orig)[grepl("TEp0-[0-9]{1}", names(orig))]
orig <- orig %>%
  dplyr::select(-c(cols))
orig <- orig %>%
  dplyr::select(c(Merged.Region, Chromosome, Start, End,
                  CGIslandCount, PromoterCount, GeneCount))
names(orig)[1] <- "peak"

#j91i
orig_j91 <- readWorkbook("./data/4441NIH_J91i-MeDIP_mergedregs.xlsx") 
cols <- names(orig_j91)[grepl("TEp0-[0-9]{1}", names(orig_j91))]
orig_j91 <- orig_j91 %>%
  dplyr::select(-c(cols))
orig_j91 <- orig_j91 %>%
  dplyr::select(c(Merged.Region, Chromosome, Start, End,
                  CGIslandCount, PromoterCount, GeneCount))
names(orig_j91)[1] <- "peak"

# j98
deseq <- readWorkbook("./data/deseq_filtered_all_contrasts.xlsx", sheet = 2)
sub_deseq <- deseq %>%
  dplyr::select(`MeDIP-Seq.peak`)
names(sub_deseq)[1] <- "peak"

# j91
deseq_j91 <- fread("./data/J91i-iPSC_TEexp_deseq_filtered.csv")
sub_deseq_j91 <- deseq_j91 %>%
  filter(padj < 0.05) %>%
  dplyr::select(MergedRegion)
names(sub_deseq_j91)[1] <- "peak"

orig <- orig %>%
  filter(peak %in% sub_deseq$peak)

orig_j91 <- orig_j91 %>%
  filter(peak %in% sub_deseq_j91$peak)

# j98
bam_files <- list.files("./data/BAM_MeDIP3", pattern = ".bam$", full.names = T)
bam_files <- bam_files[!grepl("TEp0", bam_files)]
bam_nam <- gsub(".*NIH_|_MeDIP.*", "", bam_files)

# j91
bam_files_j91 <- list.files("./data/BAM_MeDIP3_j91l/", pattern = ".bam$", full.names = T)
bam_files_j91 <- bam_files_j91[!grepl("TEp0", bam_files_j91)]
bam_nam_j91 <- gsub(".*NIH_|_MeDIP.*", "", bam_files_j91)

##### CODE FOR CPG ISLANDS #####

# filter to retain only CPG island counts &
# upstream/downstream 1000 bp
cpg <- orig %>%
  filter(CGIslandCount > 0) %>%
  dplyr::select(-c(PromoterCount, GeneCount)) %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Start - 1000,
         End = End + 1000)
names(cpg)[2:4] <- c("chr", "start", "end")

cpg_j91 <- orig_j91 %>%
  filter(CGIslandCount > 0) %>%
  dplyr::select(-c(PromoterCount, GeneCount)) %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Start - 1000,
         End = End + 1000)
names(cpg_j91)[2:4] <- c("chr", "start", "end")

# split out cpg by peak
cpg_adj <- split(cpg, cpg$peak)
cpg_adj_j91 <- split(cpg_j91, cpg_j91$peak)

# bin the st/end sites into 10 bins
cpg_adj <- lapply(1:length(cpg_adj), function(x) 
  bin_methyl(cpg_adj[[x]], bin = 10))
cpg_adj_j91 <- lapply(1:length(cpg_adj_j91), function(x) 
  bin_methyl(cpg_adj_j91[[x]], bin = 10))

# put peak ranges into granges object
cpg_adj <- lapply(cpg_adj, function(x) makeGRangesFromDataFrame(x))
cpg_adj_j91 <- lapply(cpg_adj_j91, function(x) makeGRangesFromDataFrame(x))

# extract counts for peaks
cts <- lapply(1:length(cpg_adj), function(x) {
  lapply(1:length(bam_files), function(y) {
    bamCount(bam_files[y], cpg_adj[[x]], verbose = T)
  })
})
cts_j91 <- lapply(1:length(cpg_adj_j91), function(x) {
  lapply(1:length(bam_files_j91), function(y) {
    bamCount(bam_files_j91[y], cpg_adj_j91[[x]], verbose = T)
  })
})

# re-organize data for plotting
cts_cpg <- lapply(1:length(cts), function(x) 
  counts_reformat(cts[[x]], cpg_adj[[x]], split = 10))
cts_cpg <- rbindlist(cts_cpg)
names(cts_cpg)[2] <- "Sample-replicate"

cts_cpg_j91 <- lapply(1:length(cts_j91), function(x) 
  counts_reformat(cts_j91[[x]], cpg_adj_j91[[x]], split = 10))
cts_cpg_j91 <- rbindlist(cts_cpg_j91)
names(cts_cpg_j91)[2] <- "Sample-replicate"
saveRDS(cts_cpg_j91, file = "./adj_data/cts_cpg_j91.RDS")

# plot dna methylation for CpG Islands
plot_methyl(cts_cpg, 
            brks = c(1, 4, 7, 10),
            labs = c("-1kbp", "0bp", "0bp", "1kbp"),
            int = c(4, 7),
            name = "CpG Island Location",
            filename = "sign_cpg_j98",
            k = 10,
            h = 5,
            w = 9)
plot_methyl(cts_cpg_j91, 
            brks = c(1, 4, 7, 10),
            labs = c("-1kbp", "0bp", "0bp", "1kbp"),
            int = c(4, 7),
            name = "CpG Island Location",
            filename = "sign_cpg_j91",
            k = 10,
            h = 5,
            w = 9)
rm(cts, cts_cpg, cpg_adj, cpg_adj_j91, cpg_j91, cts_j91,
   cpg, cts_cpg_j91)
gc()

##### CODE FOR PROMOTERS #####
prom <- orig %>%
  filter(PromoterCount > 0) %>%
  dplyr::select(-c(CGIslandCount, GeneCount)) %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Start - 5000,
         End = End + 5000)
names(prom)[2:4] <- c("chr", "start", "end")

prom_j91 <- orig_j91 %>%
  filter(PromoterCount > 0) %>%
  dplyr::select(-c(CGIslandCount, GeneCount)) %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Start - 5000,
         End = End + 5000)
names(prom_j91)[2:4] <- c("chr", "start", "end")

# split out promoters by peak
prom_adj <- split(prom, prom$peak)
prom_adj_j91 <- split(prom_j91, prom_j91$peak)

# bin the st/end sites into 10 bins
prom_adj <- lapply(1:length(prom_adj), function(x) 
  bin_methyl(prom_adj[[x]], bin = 40))
prom_adj_j91 <- lapply(1:length(prom_adj_j91), function(x) 
  bin_methyl(prom_adj_j91[[x]], bin = 40))

# put peak ranges into granges object
prom_adj <- lapply(prom_adj, function(x) makeGRangesFromDataFrame(x))
prom_adj_j91 <- lapply(prom_adj_j91, function(x) makeGRangesFromDataFrame(x))

# extract counts for peaks
cts <- lapply(1:length(prom_adj), function(x) {
  lapply(1:length(bam_files), function(y) {
    bamCount(bam_files[y], prom_adj[[x]], verbose = T)
  })
})

cts_j91 <- lapply(1:length(prom_adj_j91), function(x) {
  lapply(1:length(bam_files), function(y) {
    bamCount(bam_files[y], prom_adj_j91[[x]], verbose = T)
  })
})

# re-organize data for plotting
cts_prom <- lapply(1:length(cts), function(x) 
  counts_reformat(cts[[x]], prom_adj[[x]], split= 40))
cts_prom <- rbindlist(cts_prom)
names(cts_prom)[2] <- "Sample-replicate"

cts_prom_j91 <- lapply(1:length(cts_j91), function(x) 
  counts_reformat(cts_j91[[x]], prom_adj_j91[[x]], split= 40))
cts_prom_j91 <- rbindlist(cts_prom_j91)
names(cts_prom_j91)[2] <- "Sample-replicate"
saveRDS(cts_prom_j91, file = "./adj_data/cts_prom_j91.RDS")

# plot dna methylation for promoters
plot_methyl(cts_prom, 
            brks = c(1, 20, 40),
            labs = c("-5kbp", "0bp", "5kbp"),
            int = c(20),
            name = "Promoter Location",
            filename = "sign_promoters_j98",
            k = 20,
            h = 5, 
            w = 9)

plot_methyl(cts_prom_j91, 
            brks = c(1, 20, 40),
            labs = c("-5kbp", "0bp", "5kbp"),
            int = c(20),
            name = "Promoter Location",
            filename = "sign_promoters_j91",
            k = 20,
            h = 5, 
            w = 9)
rm(cts_prom, cts_prom_j91, prom, prom_adj,
   prom_adj_j91, prom_j91)
gc()

##### CODE FOR GENE BODIES #####
gb <- orig %>%
  filter(GeneCount > 0) %>%
  dplyr::select(-c(CGIslandCount, PromoterCount)) %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Start - 2000,
         End = End + 2000)
names(gb)[2:4] <- c("chr", "start", "end")

gb_j91 <- orig_j91 %>%
  filter(GeneCount > 0) %>%
  dplyr::select(-c(CGIslandCount, PromoterCount)) %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Start = Start - 2000,
         End = End + 2000)
names(gb_j91)[2:4] <- c("chr", "start", "end")

# split out cpg by peak
gb_adj <- split(gb, gb$peak)
gb_adj_j91 <- split(gb_j91, gb_j91$peak)

# bin the st/end sites into 10 bins
gb_adj <- lapply(1:length(gb_adj), function(x) 
  bin_methyl(gb_adj[[x]], bin = 40))
gb_adj_j91 <- lapply(1:length(gb_adj_j91), function(x) 
  bin_methyl(gb_adj_j91[[x]], bin = 40))

# put peak ranges into granges object
gb_adj <- lapply(gb_adj, function(x) makeGRangesFromDataFrame(x))
gb_adj_j91 <- lapply(gb_adj_j91, function(x) makeGRangesFromDataFrame(x))

# extract counts for peaks
cts <- lapply(1:length(gb_adj), function(x) {
  lapply(1:length(bam_files), function(y) {
    bamCount(bam_files[y], gb_adj[[x]], verbose = T)
  })
})

cts_j91 <- lapply(1:length(gb_adj_j91), function(x) {
  lapply(1:length(bam_files), function(y) {
    bamCount(bam_files[y], gb_adj_j91[[x]], verbose = T)
  })
})

# re-organize data for plotting
cts_gb <- lapply(1:length(cts), function(x) 
  counts_reformat(cts[[x]], gb_adj[[x]], split = 40))
cts_gb <- rbindlist(cts_gb)
names(cts_gb)[2] <- "Sample-replicate"

cts_gb_j91 <- lapply(1:length(cts_j91), function(x) 
  counts_reformat(cts_j91[[x]], gb_adj_j91[[x]], split = 40))
cts_gb_j91 <- rbindlist(cts_gb_j91)
names(cts_gb_j91)[2] <- "Sample-replicate"
saveRDS(cts_gb_j91, file = "./adj_data/cts_gb_j91.RDS")

# plot dna methylation for promoters
plot_methyl(cts_gb, 
            brks = c(1, 5, 35, 40),
            labs = c("-2kbp", "0bp", "0bp", "2kbp"),
            int = c(5, 35),
            name = "Gene Body Location",
            filename = "sign_genebodies_j98",
            k = 20,
            h = 5,
            w = 9)

plot_methyl(cts_gb_j91, 
            sp = 0.2,
            brks = c(1, 5, 35, 40),
            labs = c("-2kbp", "0bp", "0bp", "2kbp"),
            int = c(5, 35),
            name = "Gene Body Location",
            filename = "sign_genebodies_j91",
            k = 20,
            h = 5,
            w = 9)
rm(cts, cts_gb, cts_gb_j91, cts_j91, gb, 
   gb_adj, gb_adj_j91, gb_j91)
gc()

##### HISTOGRAM FOR THE LONG ARM OF CHR 21 (CHR 21) #####
chr21 <- orig %>%
  filter(Chromosome == 21) %>%
  dplyr::select(c(1:4))
names(chr21)[2:4] <- c("chr", "start", "end")

chr21_j91 <- orig_j91 %>%
  filter(Chromosome == 21) %>%
  dplyr::select(c(1:4))
names(chr21_j91)[2:4] <- c("chr", "start", "end")

# split out cpg by peak
chr21 <- split(chr21, chr21$peak)
chr21_j91 <- split(chr21_j91, chr21_j91$peak)


# bin the st/end sites into 10 bins
chr21 <- lapply(1:length(chr21), function(x) 
  bin_methyl(chr21[[x]], bin = 40))
chr21_j91 <- lapply(1:length(chr21_j91), function(x) 
  bin_methyl(chr21_j91[[x]], bin = 40))

chr21 <- lapply(chr21, function(x) x %>% mutate(chr = paste0("chr", chr)))
chr21_j91 <- lapply(chr21_j91, function(x) x %>% mutate(chr = paste0("chr", chr)))

# put peak ranges into granges object
chr21 <- lapply(chr21, function(x) makeGRangesFromDataFrame(x))
chr21_j91 <- lapply(chr21_j91, function(x) makeGRangesFromDataFrame(x))

# extract counts for peaks
cts <- lapply(1:length(chr21), function(x) {
  lapply(1:length(bam_files), function(y) {
    bamCount(bam_files[y], chr21[[x]], verbose = T)
  })
})

cts_j91 <- lapply(1:length(chr21_j91), function(x) {
  lapply(1:length(bam_files_j91), function(y) {
    bamCount(bam_files_j91[y], chr21_j91[[x]], verbose = T)
  })
})

# re-organize data for plotting
cts_chr21 <- lapply(1:length(cts), function(x) 
  counts_reformat(cts[[x]], chr21[[x]], bam_nam, split = 40))
cts_chr21 <- rbindlist(cts_chr21)
names(cts_chr21)[2] <- "Sample-replicate"

cts_chr21_j91 <- lapply(1:length(cts_j91), function(x) 
  counts_reformat(cts_j91[[x]], chr21_j91[[x]], bam_nam_j91, split = 40))
cts_chr21_j91 <- rbindlist(cts_chr21_j91)
names(cts_chr21_j91)[2] <- "Sample-replicate"

# pull length of long arm of chr21
chr <- fread("./data/cytoBand.txt")
names(chr) <-  c("chrom","chromStart","chromEnd","name","gieStain")
chr <- chr[ , .(length = sum(chromEnd - chromStart)), 
     by = .(chrom, arm = substring(name, 1, 1)) ] %>%
  filter(chrom == "chr21" & arm == "q")

# plot dna methylation for chr21
p <- ggplot(cts_chr21) +
  geom_bar(aes(x = cut, y = n, color = `Sample-replicate`, fill = `Sample-replicate`), stat = "identity") +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40),
                     labels = c(0, round(chr$length[1] * (1/5)),  round(chr$length[1] * (3/5)), round(chr$length[1] * (4/5)), chr$length[1]),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("firebrick1", "tomato", "lightsalmon",
                               "blue1", "deepskyblue", "lightblue1")) +
  scale_color_manual(values = c("firebrick1", "tomato", "lightsalmon",
                                "blue1", "deepskyblue", "lightblue1")) +
  xlab("Long arm of chromosome 21 (bp)") + ylab("Normalized Counts") +
  theme_classic() +
  mytheme 
ggsave(filename = "./adj_data/plots/chr21_q_j98_global_methylation_plot.png", 
       units = "in",
       height = 5, width = 8, plot = p)
pdf(file = "./adj_data/plots/chr21_q_j98_global_methylation_plot.pdf",
    height = 5, width = 8)
print(p)
graphics.off()

p <- ggplot(cts_chr21_j91) +
  geom_bar(aes(x = cut, y = n, color = `Sample-replicate`, fill = `Sample-replicate`), stat = "identity") +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40),
                     labels = c(0, round(chr$length[1] * (1/5)),  round(chr$length[1] * (3/5)), round(chr$length[1] * (4/5)), chr$length[1]),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("firebrick1", "tomato", "lightsalmon",
                               "blue1", "deepskyblue", "lightblue1")) +
  scale_color_manual(values = c("firebrick1", "tomato", "lightsalmon",
                                "blue1", "deepskyblue", "lightblue1")) +
  xlab("Long arm of chromosome 21 (bp)") + ylab("Normalized Counts") +
  theme_classic() +
  mytheme 
ggsave(filename = "./adj_data/plots/chr21_q_j91_global_methylation_plot.png", 
       units = "in",
       height = 5, width = 8, plot = p)
pdf(file = "./adj_data/plots/chr21_q_j91_global_methylation_plot.pdf",
    height = 5, width = 8)
print(p)
graphics.off()


rm(list = ls())
gc()
