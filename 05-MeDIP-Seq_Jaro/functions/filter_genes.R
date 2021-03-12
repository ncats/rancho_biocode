#' Filter out unwanted genes

filter_genes <- function(file, loc, name, output) {
  # read in deseq output for given contrast
  df <- fread(file)
  
  y <- df %>%
    # separate data from anotatr's annotation
    # into separate rows
    separate_rows(strand, id, tx_id, gene, gene_anotatr, type, width, sep = "; ", convert = TRUE) %>% 
    # remove genes from anotatr annotation
    filter(!gene_anotatr %in% fg$Gene) %>%
    group_by(MergedRegion) %>%
    # recombine rows
    mutate(strand = paste(strand, collapse = "; "),
           id = paste(id, collapse = "; "),
           tx_id = paste(tx_id, collapse = "; "),
           gene = paste(gene, collapse = "; "),
           gene_anotatr = paste(gene_anotatr, collapse = "; "),
           type = paste(type, collapse = "; "),
           width = paste(width, collapse = "; ")) %>% 
    # take first slice of data grouped by MergedRegion
    slice(1) %>%
    ungroup() %>%
    # separate data from ActiveMotif's annotation
    # into separate rows
    separate_rows(gene_active, location, distance_to_start, sep = "; ", convert = TRUE) %>% 
    # remove genes from ActiveMotif's annotation
    filter(!gene_active %in% fg$Gene) %>%
    group_by(MergedRegion) %>%
    # recombine rows
    mutate(gene_active = paste(gene_active, collapse = "; "),
           location = paste(location, collapse = "; "),
           distance_to_start = paste(distance_to_start, collapse = "; ")) %>%
    # take first slice of data grouped by MergedRegion
    slice(1) %>%
    ungroup() %>%
    # add specified contrast
    mutate(contrast = name)
  
  # rename columns
  names(y) <- gsub("^X", "", names(y))
  
  # write out to file
  write.csv(y, file = paste(loc, name, output, "_filtered.csv", sep = ""), 
            row.names = F)
}
