#' Build top 50 heatmaps from normalized counts
#' across all contrasts

# req'd pkgs
x <- c("magrittr", "pheatmap", "tidyverse", "data.table",
       "RColorBrewer", "scales")
sapply(x, library, character.only = TRUE)

# source custom fxns
source("./functions/save_pheatmap_pdf.R")

# load deseq files
all <- list.files("./adj_data/deseq/filtered/", full.names = T)

# remove unncessary items from naming contrasts
all_nam <- gsub(".*\\/|_deseq.*", "", all)
all_nam <- gsub("-", ".", all_nam)

# read in deseq files, sort by padj, keep top
# 50 medip-seq peaks & keep MergedRegion
# and selected contrast
all <- lapply(all, function(x) {
  y <- fread(x) %>%
    arrange(padj)
  y <- y[c(1:50), ]
  cols <- names(y)[grepl("^[0-9]{2}\\_", names(y))]
  y <- y %>% dplyr::select(c(MergedRegion, cols))
  return(y)
})

# massive lapply loop to build heatmaps - 
# detailed below:
all <- lapply(1:length(all), function(x) {
  # adjust contrast names so the correct
  # contrasts are pulled for each plot
  nam <- all_nam[x]
  nam <- gsub("_", "|", nam)
  cols <- names(all[[x]])[grepl(nam, names(all[[x]]))]
  
  # keep the MergedRegion & counts of
  # the contrasts
  y <- all[[x]] %>%
    dplyr::select(c(1, cols)) %>%
    column_to_rownames("MergedRegion")
  
  # address one treatment w/three samples
  nam <- data.frame("Treatment" = gsub("^[0-9]{2}\\_|\\.[0-9]{1}$", "", names(y)))
  # define the annotation by sample, treatment
  annot = data.frame(
    Treatment = gsub("\\.", "-", nam$Treatment)
  )
  
  # set rownames of annotation df so 
  # colors will plot correctly
  rownames(annot) <- names(y)

  # specify color of annotation
  treat <- hue_pal()(2)
  # specify annotation colors
  treat_samp <- list(
    Treatment = treat
  )
  # set names of elements in the list
  names(treat_samp[[1]]) <- unique(annot$Treatment)
  
  # adjust title for contrast
  new_title <- gsub("_", "_vs_", all_nam[[x]])
  new_title <- gsub("\\.", "-", new_title)
  
  # plot heatmap w/the parameters:
  # blue -> red heatmap; treatment
  # annotation; treatment colors;
  # no border, row, colnames
  p <- pheatmap(log2(y + 1), 
                color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
                annotation_col = annot, 
                annotation_colors = treat_samp,
                border_color = NA,
                show_rownames = F, 
                show_colnames = F,
                main = paste("Top 50 Significant MeDIP-Seq Peak regions for ", new_title, sep = ""))
  
  filename <- paste("./adj_data/plots/heatmap/", new_title, "_heatmap.pdf", sep = "")
  save_pheatmap_pdf(p, filename)
})

#### CREATE A SINGLE HEATMAP PLOT ####
# load deseq files
all <- list.files("./adj_data/deseq/filtered/", full.names = T)

# read in deseq files; keep cols for given
# contrast; filter padj < 0.05; keep only
# peak # and cols for given contrast
all <- lapply(1:length(all), function(x) {
  y <- fread(all[[x]]) 
  cols <- names(y)[grepl("[0-9]{2}_", names(y))]
  y <- y %>%
    dplyr::select(c(c(MergedRegion, log2FoldChange, padj, 
                    chromosome_name_medip, start_position_medip,
                    end_position_medip), all_of(cols))) %>%
    filter(padj < 0.05) %>%
    dplyr::select(c(MergedRegion, all_of(cols)))
  return(y)
})

# filter to retain only overlap in peak #
# across all contrasts
all[[1]] <- all[[1]] %>% filter(MergedRegion %in% all[[2]]$MergedRegion & MergedRegion %in% all[[3]]$MergedRegion)
all[[2]] <- all[[2]] %>% filter(MergedRegion %in% all[[1]]$MergedRegion & MergedRegion %in% all[[3]]$MergedRegion)
all[[3]] <- all[[3]] %>% filter(MergedRegion %in% all[[1]]$MergedRegion & MergedRegion %in% all[[2]]$MergedRegion)

# keep only TEexp column names
all[[2]] <- all[[2]] %>% dplyr::select(c(1, 5:7))
# keep peak #, iPSC, TEp0, TEexp
all <- all[[1]] %>% left_join(., all[[2]])
# put peak # in row name
all <- all %>%
  column_to_rownames("MergedRegion")

# address one treatment w/three samples
nam <- data.frame("Sample" = gsub("^[0-9]{2}\\_|\\.[0-9]{1}$", "", names(all)))
# define the annotation by sample, treatment
annot = data.frame(
  Sample = gsub("\\.", "-", nam$Sample)
)

# set rownames of annotation df so 
# colors will plot correctly
rownames(annot) <- names(all)

# specify color of annotation
treat <- hue_pal()(3)
# specify annotation colors
treat_samp <- list(
  Sample = treat
)
# set names of elements in the list
names(treat_samp[[1]]) <- unique(annot$Sample)

# adjust title for contrast
new_title <- "All_treat_sign_peaks_heatmap"

# Replace title of legend w/nothing
names(treat_samp) <- " "
names(annot) <- " "

# plot heatmap w/the parameters:
# blue -> red heatmap; treatment
# annotation; treatment colors;
# no border, row, colnames
p <- pheatmap(log2(all + 1), 
              cutree_rows = 4,
              color = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100),
              annotation_col = annot, 
              annotation_names_col = F,
              annotation_colors = treat_samp,
              border_color = NA,
              show_rownames = F, 
              show_colnames = F,
              main = "Differentially methylated peaks for all samples")

filename <- paste("./adj_data/plots/heatmap/", new_title, "_heatmap.pdf", sep = "")
save_pheatmap_pdf(p, filename)

rm(list = ls())
gc()
