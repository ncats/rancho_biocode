# Ready for Review
#' Code provided by CM to update Vukasin's heatmaps 
#' for manuscript on astrocytes.
#' Working in Rv4.0.4

# req'd pkgs
x <- c("data.table", "RColorBrewer", "gplots", "ggplot2",
       "stringr", "dendextend", "ComplexHeatmap", "gtools",
       "DESeq2", "openxlsx", "dplyr", "scales", "Cairo",
       "ggtext", "extrafont", "tibble" , "tidyr")
sapply(x, library, character.only = T)

source("./functions/adj_pca.R")
source("./functions/heatmap_vuk.R")

# pull fetal genes (from CM)
fetal.genes <- as.data.table(readxl::read_xlsx('./data/Copy of Fetal vs adult astrocytes 1-s2.0-S0896627315010193-mmc5.xlsx'))
fetal.genes <- unique(fetal.genes$`...1`)[1:100]
fetal.genes <- fetal.genes[-1]

# pull adult genes (from CM)
adult.genes <- as.data.table(readxl::read_xlsx('./data/Copy of Fetal vs adult astrocytes 1-s2.0-S0896627315010193-mmc5.xlsx', sheet=2))
adult.genes <- unique(adult.genes$`...1`)[1:100]
adult.genes <- adult.genes[-1]

# load v5 dds obj
load("./data/Vukasin.DDS.RData")
old <- dds
old <- old[ , (old$condition %in% c('NCRM5_astro_D0', 'NCRM5_astro_D30', 'NCRM5_astro_D50'))]
load("./data/Vukasin_SRAv5.DDS.RData")
cts_old <- data.frame(assay(old))%>%
  rownames_to_column("gene")
cts_new <- data.frame(assay(dds)) %>%
  rownames_to_column("gene")
cts <- cts_old %>% 
  full_join(., cts_new, by = "gene")
rownames(cts) <- cts$gene
cts$gene <- NULL
cts <- cts %>% mutate_all(funs(replace_na(.,0)))
old_meta <- data.frame(colData(old)) %>%
  filter(grepl("NCRM5_astro_", condition))
new_meta <- data.frame(colData(dds))
old_meta <- rbind(old_meta, new_meta) %>%
  dplyr::select(-c(W_1, sizeFactor))
dds <- DESeqDataSetFromMatrix(countData = cts, colData = old_meta, design = ~ condition)

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
dds.genes <- row.names(dds.stabilized)

####### PLOT1 ########
# subset deseq2 obj to remove adult astrocytes
# plot adult/fetal genes
genes <- c(fetal.genes, adult.genes)
dds_no_adult <- dds.stabilized[ , !(dds.stabilized$condition %in% c("NCRM5_astro_D0", "adult_astrocytes_TL_Zhang", 'adult_astrocyte_HC_Zhang'))]
all_plot <- heatmap_vuk(genes)
png(file = "./results/zhang_all_genes.png", width = 3000, height = 4000, units = "px", res = 300)
print(all_plot)
graphics.off()

# found <- dds.genes[dds.genes %in% genes]
# length(found)
# length(genes)
# length(genes[!genes %in% found] ) # 100
# genes[!genes %in% found] # [1] 


####### PLOT2 ########
# subset deseq2 obj to remove adult astrocytes
# plot fetal genes only
fetal <- fetal.genes
fetal_plot <- heatmap_vuk(fetal)
png(file = "./results/zhang_fetal_genes_only.png", width = 3000, height = 4000, units = "px", res = 300)
print(fetal_plot)
graphics.off()

####### PLOT3 #######
# create pca colored appropriately
#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20),
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(size = 16, face = "italic", family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))
dds_no_adult <- dds[ , !(dds$condition %in% c("adult_astrocytes_TL_Zhang", 'adult_astrocyte_HC_Zhang'))]
dds_no_adult$condition <- ifelse(grepl("^NCRM5_astro_D30$|^NCRM5_astro_D50$", dds_no_adult$condition), "Jovanovic et al.\n sf",
                                 ifelse(grepl("^Astro", dds_no_adult$condition), "Jovanovic et al.\n s",
                                        ifelse(grepl("fetal", dds_no_adult$condition), "Primary astrocytes",
                                               ifelse(grepl("iPSC_derived", dds_no_adult$condition), "Santos et al.",
                                                      ifelse(grepl("hiPSC", dds_no_adult$condition), "TCW et al.",
                                                             ifelse(grepl("Tchieu", dds_no_adult$condition), "Tchieu et al.",
                                                                    ifelse(grepl("^NCRM5_astro_D0", dds_no_adult$condition), "NCRM5 D0", "NA")))))))
rld <- vst(dds_no_adult, blind = T)
p <- adj_pca(rld, intgroup = "condition")
p <- p +
  geom_vline(xintercept = 0, color = "gray") +
  geom_hline(yintercept = 0, color = "gray") +
  labs(x = "PC1 (35%)",
       y = "PC2 (25%)") 
        

p <- p +
  scale_color_manual(labels = c("Jovanovic et al.\n s", "Jovanovic et al.\n sf",
                                "NCRM5 D0", "Primary astrocytes", "Tchieu et al.", 
                                "TCW et al.", "Santos et al."),
                     values = c("plum2", "purple", "coral", "gray", "red", "deepskyblue", "forestgreen")) +
  coord_fixed(ratio = 1, xlim = c(-60, 60), ylim = c(-60, 60))
p$labels$colour <- ""

ggsave(filename = "./results/pca_all.png",
       plot =  p)

rm(list = ls())
gc()
