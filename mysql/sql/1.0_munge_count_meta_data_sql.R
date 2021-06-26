#' Extract count matrices and metadata from DESeq2
#' dds object or seurat object, clean it up and write
#' it to a .csv (bulk RNA-Seq) or .RDS file (scRNA-Seq).
#' For scRNA-Seq, if the count matrix is > 500e6, then
#' split up the data into files of max 500e6 rows ea.

x <- c("RMySQL", "data.table", "dplyr", "tidyr",
       "Seurat", "tibble", "DESeq2")
sapply(x, library, character.only = TRUE)

options(scipen=999)

######## MUNGE ALL RNA-SEQ DATA FOR COMPLEX APP ########

###BULK###
new_bulk <- list.files("/home/ubuntu/bulk", pattern = ".RData|.rdata|.Rdata", full.names = T)
nam <- gsub(".*\\/|_.*|\\..*", "", new_bulk)
if (length(new_bulk) > 0) {
  for (i in 1:length(new_bulk)) {
    load(new_bulk[i])
    cts <- counts(dds) %>%
      data.frame(.) %>%
      rownames_to_column("gene") %>%
      tidyr::gather(sample_id, counts, -gene) %>%
      mutate_all(funs(gsub("Lonza|LONZA|iPSC|IPSC", "LiPSC.GR1.1", .)))
    col <- colData(dds) %>%
      data.frame(.) %>%
      rownames_to_column("sample_id")
    rownames(col) <- NULL
    names(col) <- tolower(names(col))
    
    cols_exist <- names(col)[grepl("sizefactor|W_1|replicate", names(col))]
    
    if (length(cols_exist) > 0) {
        col <- col %>% select(-cols_exist)
    }
    
    if (sum(grepl("treatment", names(col))) == 1) {
        col$treatment <- gsub("^Y$|y27|y27.*", "Y27", col$treatment)
    }
    
    if (sum(grepl("media", names(col))) == 1) {
        col$media <- gsub("", "None", col$media)
    }
    
    if (sum(grepl("day", names(col))) == 1) {
        col$day <- gsub("d", "D", col$day)
    }
    
    if (sum(grepl("cell_line", names(col))) == 1) {
        col$cell_line <- gsub("iPSC|Lonza|lonza|LONZA|IPSC", "LiPSC.GR1.1", col$cell_line)
    }
    
    if (sum(grepl("mutation", names(col))) == 1) {
        col$mutation <- gsub("", "None", col$mutation)
    }
    
    if (sum(grepl("dosage", names(col))) == 1) {
        col$dosage <- gsub("", "None", col$dosage)
    }
    
    if (sum(grepl("differentiation", names(col))) == 1) {
        col$differentiation <- gsub("", "None", col$differentiation)
    }
    
    if (sum(grepl("condition", names(col))) == 1) {
        col$condition <- gsub("\\-|\\.", "_", col$condition)
    } else if (sum(grepl("condition", names(col))) == 0) {
        cols <- names(col)
        cols <- cols[-sample_id]
        col <- col %>% mutate(condition = paste(cols, sep = "_"))
    }
    
    write.csv(cts, file = paste0("/home/ubuntu/adj_data/bulk/", nam[i], "_counts.csv"), row.names = F)
    write.csv(col, file = paste0("/home/ubuntu/adj_data/bulk/", nam[i], "_meta.csv"), row.names = F)
    rm(dds, cts, col)
  }
} else {
  print("There are no new bulk RNA-Seq datasets to munge!")
}
rm(list = ls())

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

### SC-RNA-SEQ###
new_sc <- list.files("/home/ubuntu/sc", pattern = ".RData|.rdata|.Rdata", full.names = T)
nam_short <- gsub(".*\\/|_.*|\\..*", "", new_sc)
if (length(new_sc) > 0) {
  for (i in 1:length(new_sc)) {
    seur <- loadRData(new_sc[i])
    seur2 <- seur
    
    # update for seurat obj
    counts <- GetAssayData(seur, slot = "counts")
    counts <- data.frame(counts)
    nam <- rownames(counts)
    counts <- counts %>% 
      mutate(gene = nam) %>%
      gather(cell_id, counts, -gene) %>%
      mutate_all(funs(gsub("Lonza|LONZA|iPSC|IPSC", "LiPSC.GR1.1", .)))
    
    # update for seurat obj
    sub <- seur@meta.data
    nam <- gsub("-", ".", rownames(sub))
    sub <- data.frame("cell_id" = nam, 
                      "sample" = sub$sample)
    
    counts <- counts %>% left_join(., sub, by = "cell_id")
    counts <- counts[, c(2, 4, 3, 1)]
    
    # extract total rows
    nr <- nrow(counts)
    
    if (nr > 500e+6) {
      # split data into 500e6 rows 
      len <- seq(1, nr, by = 500000000)
      len <- sapply(1:length(len), function(x) {
        if (len[x] == min(len)) {
          a <- c(len[x], len[x+1] - 1)
        } else if (len[x] == max(len)) {
          a <- c(len[x], nr)
        } else {
          a <- c(len[x], len[x+1] - 1)
        }
      })
      
      # save a csv for df by 500e6 rows
      lapply(1:ncol(len), function(x) {
        st <- len[1,x]
        end <- len[2,x]
        int <- paste0(st, ":", end)
        n <- x
        saveRDS(counts[int,], file = paste0("/home/ubuntu/adj_data/sc/", nam_short[i], "_part", n,".RDS"))
      })
      
      rm(sub, nam, seur, len, nr, counts)
      
    } else {
      saveRDS(counts, file = paste0("/home/ubuntu/adj_data/sc/", nam_short[i], "_counts.RDS"))
      rm(sub, nam, seur, nr, counts)
      gc()
    }
  
  # need to add required data & make sure all columns
  # are the same as the db
    keep_cols <- c("cell_id", "orig.ident", "nFeature_RNA",
                   "nCount_RNA", "seurat_clusters", "sample",
                   "treatment", "day", "differentiation",
                   "cell_line", "method", "media",
                   "phenotype", "location")
    meta <- seur2@meta.data
    res <- names(meta)[grepl("RNA_snn_res", names(meta))]
    keep <- keep_cols[keep_cols %in% names(meta)]
    meta <- meta %>% 
      rownames_to_column("cell_id") %>%
      mutate(cell_id = gsub("-", ".", cell_id))
      
    if (sum(names(meta) %in% "sample") == 1) {
        meta <- meta %>%
                select(c(keep, res))
    } else {
        excl <- c("cell_id", "orig.ident", "nFeature_RNA",
                  "nCount_RNA", "seurat_clusters")
        cond <- names(meta)[!names(meta) %in% c(excl, res)]
        
        my_func <- function(dat, vars){
             .vars <- rlang::syms(vars)

             result <- dat %>%
                  mutate(sample = paste(!!!.vars, sep = "_" ))
             return(result)
        }
        
        meta <- meta %>%
                dplyr::select(c(keep, res))
        meta <- my_func(meta, cond)
    }
    
    if (sum(grepl("treatment", names(meta))) == 1) {
        meta$treatment <- gsub("^Y$|y27|y27.*|^y", "Y27", meta$treatment)
    }
    if (sum(grepl("day", names(meta))) == 1) {
        meta$day <- gsub("d", "D", meta$day)
    }
    if (sum(grepl("cell_line", names(meta))) == 1) {
        meta$cell_line <- gsub("iPSC|Lonza|lonza|LONZA|IPSC", "LiPSC.GR1.1", meta$cell_line)
    }
    if (sum(grepl("differentiation", names(meta))) == 1) {
        meta$differentiation <- gsub("", "None", meta$differentiation)
    }
    if (sum(grepl("method", names(meta))) == 1) {
        meta$method <- gsub("", "None", meta$method)
    }
    if (sum(grepl("media", names(meta))) == 1) {
        meta$media <- gsub("", "None", meta$media)
    }
    if (sum(grepl("phenotype", names(meta))) == 1) {
        meta$phenotype <- gsub("", "None", meta$phenotype)
    }
    if (sum(grepl("location", names(meta))) == 1) {
        meta$location <- gsub("", "None", meta$location)
    }
      
    saveRDS(meta, file = paste0("/home/ubuntu/adj_data/sc/", nam_short[i], "_meta.RDS"))
    rm(meta)
    gc()
  }
} else {
  print("There are no new scRNA-Seq datasets to munge!")
}

rm(list = ls())
gc()
