#' Push cleaned and munged bulk and scRNA-Seq data to
#' mysql rds.
#' For scRNA-Seq, if the count matrix is > 500e6, then
#' data is pushed one file at a time (first overwriting,
#' then appending) until the entire count matrix is pushed
#' to the mysql rds.

x <- c("RMySQL", "data.table", "dplyr", "tidyr",
       "RMySQL", "tibble")
sapply(x, library, character.only = TRUE)

options(scipen=999)

######## UPDATE SC RNA-SEQ DATASETS IN SCDB ########

###BULK###
new_bulk <- list.files("/home/ubuntu/adj_data/bulk", pattern = ".csv", full.names = T)
name <- gsub(".*\\/|\\.csv", "", new_bulk)
uni_name <- unique(gsub("_.*", "", name))

if (length(new_bulk) > 0) {
  if (length(uni_name) > 1) {
    for (i in 1:length(uni_name)) {
      # create conn to SQL db BULKDB
      conn <- dbConnect(MySQL(), user = 'blk-app-admin', 
                        password = 'Bisi9adm1n', 
                        dbname = 'BULKDB', 
                        host = 'sctl-dev-rds-mysql-01.ceyknq0yekb3.us-east-1.rds.amazonaws.com', 
                        port = 3306)
      files <- new_bulk[grepl(uni_name[i], new_bulk)]
      cts <- files[grepl("counts", files)]
      meta <- files[grepl("meta", files)]
      for (t in 1:length(cts)) {
        tab <- fread(cts[t])
        if (t == 1) {
          dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_counts"), overwrite = TRUE, row.names = F)
          rm(tab)
        } else {
          dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_counts"), append = TRUE, row.names = F)
          rm(tab)
        }
      }
      tab <- fread(meta)
      dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_meta"), overwrite = TRUE, row.names = F)
      rm(tab)
    }
  } else {
    files <- new_bulk[grepl(uni_name[i], new_bulk)]
    cts <- files[grepl("counts", files)]
    meta <- files[grepl("meta", files)]
    tab <- fread(cts)
    dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_counts"), overwrite = TRUE, row.names = F)
    tab <- readRDS(meta)
    dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_meta"), overwrite = TRUE, row.names = F)
    rm(tab)
  }
} else {
  print("There are no count or meta data matrices to push to SQL!")
}
rm(list = ls())
gc()

# create conn to SQL db SCDB
new_sc <- list.files("/home/ubuntu/adj_data/sc", pattern = ".RDS", full.names = T)
name <- gsub(".*\\/|\\.RDS", "", new_sc)
uni_name <- unique(gsub("_.*", "", name))

if (length(uni_name) > 1) {
  for (i in 1:length(uni_name)) {
    # create conn to SQL db BULKDB
    conn <- dbConnect(MySQL(), user = 'sc-app-admin', 
                      password = 'Bisi9adm1n', 
                      dbname = 'SCDB', 
                      host = 'sctl-dev-rds-mysql-01.ceyknq0yekb3.us-east-1.rds.amazonaws.com', 
                      port = 3306)
    files <- new_sc[grepl(uni_name[i], new_sc)]
    cts <- files[grepl("counts", files)]
    meta <- files[grepl("meta", files)]
    for (t in 1:length(cts)) {
      tab <- readRDS(cts[t])
      if (t == 1) {
        dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_counts"), overwrite = TRUE, row.names = F)
        rm(tab)
      } else {
        dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_counts"), append = TRUE, row.names = F)
        rm(tab)
      }
    }
    tab <- readRDS(meta)
    dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_meta"), overwrite = TRUE, row.names = F)
    rm(tab)
  }
} else {
  files <- new_sc[grepl(uni_name[i], new_sc)]
  cts <- files[grepl("counts", files)]
  meta <- files[grepl("meta", files)]
  tab <- readRDS(cts)
  dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_counts"), overwrite = TRUE, row.names = F)
  tab <- readRDS(meta)
  dbWriteTable(conn, value = tab, name = paste0(uni_name[i], "_meta"), overwrite = TRUE, row.names = F)
  rm(tab)
}

rm(list = ls())
gc()
