# Ready for Review
#' custom heatmap fxn

heatmap_vuk <- function(genes) {
  idx <- which(row.names(dds_no_adult) %in% genes)
  alt <- dds_no_adult[idx,]
  
  coldata <- as.data.frame(alt@colData)
  mat <- as.data.frame(assay(alt))
  names(mat) <- coldata$condition
  mat <- mat[, c(1:6, 48:50, 28:47, 12:27, 7:11)]
  names(mat) <- ifelse(grepl("fetal", names(mat)), "Primary astrocytes",
                       ifelse(grepl("iPSC_derived", names(mat)), "Santos et al.",
                              ifelse(grepl("hiPSC", names(mat)), "TCW et al.",
                                     ifelse(grepl("Tchieu", names(mat)), "Tchieu et al.",
                                            ifelse(grepl("^Astro", names(mat)), "Jovanovic et al.\n s",
                                                   ifelse(grepl("^NCRM", names(mat)), "Jovanovic et al.\n sf", "NA"))))))
  
  
  cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red
  scaled_mat <- t(scale(t(mat)))
  rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
  
  ha = HeatmapAnnotation(Sample = names(mat),
                         col = list(Sample = c(
                           "Jovanovic et al.\n s" = "plum2",
                           "Jovanovic et al.\n sf" = "purple",
                           "Tchieu et al." = "deepskyblue",
                           "Santos et al." = "forestgreen", 
                           "TCW et al." = "red",
                           "Primary astrocytes" = "gray"), gp = gpar(fontsize = 16)),
                         annotation_legend_param = list(title = "",
                                                        Sample = list(at = c("Jovanovic et al.\n s", "Jovanovic et al.\n sf",
                                                                             "Tchieu et al.", "Santos et al.", "TCW et al.",
                                                                             "Primary astrocytes")),
                                                        labels_gp = gpar(fontface = "italic", fontsize = 16)),
                         show_annotation_name = F)
  
  suppressPackageStartupMessages(library(ComplexHeatmap))
  ht_opt(heatmap_row_names_gp = gpar(fontfamily = "NotoSans-Condensed", fontsize = 6),
         legend_title_gp = gpar(fontfamily = "NotoSans-Condensed", fontsize = 16),
         legend_labels_gp = gpar(fontfamily = "NotoSans-Condensed", fontsize = 16)
  )
  
  ht <- Heatmap(as.matrix(na.omit(rescaled_mat)),
                top_annotation = ha,
                row_names_side = "right",
                column_names_side = "top",
                col = cols.use, 
                show_column_names = F, 
                show_row_names = T,
                cluster_rows = T, 
                cluster_columns = T,
                heatmap_legend_param = list(legend_height = unit(4, "cm"), title='Row Z-Score'))
  return(ht)
}
