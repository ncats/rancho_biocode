#' Remove redundancy in row data

remove_redundancy <- function(file, sheet) {
  # for ea contrast: re-order cols; collapse
  # unique values for each row in deseq output;
  # rename cols; store ea contrast in a sep
  # excel sheet
    y <- fread(file) %>%
      data.frame(.) %>%
      filter(padj < 0.05)
    # reorder cols
    y <- y %>%
      dplyr::select(c(MergedRegion, chromosome_name_medip,
                      start_position_medip, end_position_medip, gene_anotatr,
                      location, distance_to_start, id, tx_id, gene,
                      gene_active, type, width, strand, baseMean, log2FoldChange,
                      lfcSE, stat, pvalue, padj))

    # for ea medip-seq peak: remove redunant info
    # such as gene, location, etc from ea row
    y <- rbindlist(lapply(1:nrow(y), function(t) {
      # pull row
      z <- y[t, ]
      
      # collapse data into a single string fxn
      collapse <- function(col) {
        val <- unique(unlist(strsplit(col, "; ")))
        
        if (length(val) > 1) {
          val <- paste(val, collapse = "; ")
        } else {
          val <- val
        }
      }
      
      # pull unique cell data, collapse into one string
      gene_sym <- collapse(z$gene_anotatr)
      loc <- collapse(z$location)
      dis <- collapse(z$distance_to_start)
      id <- collapse(z$id)
      tx <- collapse(z$tx_id)
      gene2 <- collapse(z$gene)
      sym <- collapse(z$gene_active)
      type <- collapse(z$type)
      wid <- collapse(z$width)
      st <- collapse(z$strand)
      
      # re-create df for ea row removing redundant info
      z <- z %>%
        mutate(gene_active = sym,
               location = loc,
               distance_to_start = dis,
               id = id,
               tx_id = tx,
               gene = gene2,
               gene_anotatr = gene_sym,
               type = type,
               width = wid,
               strand = st
        )
    }))
    
    y <- data.frame(y)
    y <- y[, c(1:4, 11, 6:7, 8:10, 5, 12:14, 15:ncol(y))]
    # rename cols of deseq df
    names(y) <- c("MeDIP-Seq peak", "Chromosome",
                  "Start", "End", "ActiveMotif annotation", 
                  "Location", "Distance to start", "Anotatr id",
                  "Ensembl transcript id", "Entrez id", 
                  "Anotatr symbol", "Anotatr type", "Anotatr width",
                  "Anotatr strand", "baseMean", "log2FC", "lfcSE", 
                  "stat", "pvalue", "padj")
    
    # store in excel sheet
    addWorksheet(xlwb, paste(gsub("_", "_vs_", sheet), sep = ""))
    writeData(xlwb,
             sheet = gsub("_", "_vs_", sheet),
             x = as.data.frame(y))
}
