#' Build top 50 heatmaps from normalized counts
#' across all contrasts

# req'd pkgs
x <- c("magrittr", "pheatmap", "tidyverse", "data.table",
       "RColorBrewer", "scales", "grid")
sapply(x, library, character.only = TRUE)

# source custom fxns
source("./functions/save_pheatmap_pdf.R")

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20),
                 axis.text = element_text(size = 14, family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

# load deseq shrink files
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
  cols <- names(y)[grepl("^[0-9]{2}\\.", names(y))]
  y <- y %>% dplyr::select(c(MergedRegion, cols))
  names(y) <- gsub("NCRM5$", "D0", names(y))
  names(y) <- gsub("NCRM5.D", "D", names(y))
  return(y)
})

# massive lapply loop to build heatmaps - 
# detailed below:
all <- lapply(1:length(all), function(x) {
  # adjust contrast names so the correct
  # contrasts are pulled for each plot
  nam <- all_nam[x]
  nam <- gsub("hPSCs_|NCRM5.", "", nam)
  nam <- gsub("_", "|", nam)
  if (grepl("day0", nam)) {
    nam <- gsub("day0", "D0", nam)
    cols <- names(all[[x]])[grepl(nam, names(all[[x]]))]
  } else {
    cols <- names(all[[x]])[grepl(nam, names(all[[x]]))]
  }
  
  # keep the MergedRegion & counts of
  # the contrasts
  y <- all[[x]] %>%
    dplyr::select(c(1, cols)) %>%
    column_to_rownames("MergedRegion")
  # re-name names of cols for ease in heatmap
  #names(y) <- ifelse(grepl("D30|D50", names(y)), gsub("NCRM5\\.", "", names(y)), names(y))
  #names(y) <- ifelse(grepl("^NCRM5$", names(y)), "hPSC day0", names(y))
  
  # address one treatment w/only two
  # samples
  if (ncol(y) == 5) {
    nam <- data.frame("name" = gsub("^[0-9]{2}\\.", "", names(y)))
    nam <- nam %>% group_by(name) %>% summarise(len = n())
    
    # define the annotation by sample, treatment
    annot = data.frame(
      Treatment = unlist(sapply(1:length(nam$name), function(x) rep(nam$name[x], nam$len[x])))
    )
    # address one treatment w/three samples
  } else if (ncol(y) == 6) {
    nam <- data.frame("Treatment" = gsub("^[0-9]{2}\\.", "", names(y)))
    # define the annotation by sample, treatment
    annot = data.frame(
      Treatment = nam
    )
  }
  
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
  new_title <- gsub("hPSCs_|NCRM5\\.", "", all_nam[[x]])
  new_title <- gsub("day0", "D0", new_title)
  new_title <- gsub("_", "_vs_", new_title)
  
  # Replace title of legend w/nothing
  names(treat_samp) <- " "
  names(annot) <- " "
  
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
                cellheight = 9,
                cellwidth = 50,
                main = paste("Top 50 Significant MeDIP-Seq Peak regions for ", gsub("_", " ", new_title), sep = ""))
  filename <- paste("./adj_data/plots/heatmap/", new_title, "_heatmap.pdf", sep = "")
  save_pheatmap_pdf(p, filename)
})

rm(list = ls())
gc()
