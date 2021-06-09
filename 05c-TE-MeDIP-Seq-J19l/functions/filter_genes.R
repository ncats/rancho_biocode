#' Filter out unwanted genes

filter_genes <- function(file, loc, name, output) {
  # read in deseq output for given contrast
  df <- fread(file)
  
  y <- df %>%
    group_by(MergedRegion) %>%
    # take first slice of data grouped by MergedRegion
    dplyr::slice(1) %>%
    ungroup() %>%
    group_by(MergedRegion) %>%
    # take first slice of data grouped by MergedRegion
    dplyr::slice(1) %>%
    ungroup() %>%
    # add specified contrast
    mutate(contrast = name)
  
  # rename columns
  names(y) <- gsub("^X", "", names(y))
  
  # write out to file
  write.csv(y, file = paste(loc, name, output, "_filtered.csv", sep = ""), 
            row.names = F)
}
