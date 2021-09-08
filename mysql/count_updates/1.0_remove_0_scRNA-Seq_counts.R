#' Remove rows w/0's from scRNA-Seq MySQL count
#' matrices

# req'd pkgs
x <- c("dplyr", "data.table", "RMySQL")
sapply(x, library, character.only = T)

# create connection to SQL db
conn <- dbConnect(MySQL(), user = 'sc-app-admin',
                  password = 'Bisi9adm1n',
                  dbname = 'SCDB',
                  host = 'sctl-dev-rds-mysql-01.ceyknq0yekb3.us-east-1.rds.amazonaws.com',
                  port = 3306)

# read in sc table
all <- dbReadTable(conn, "sc")
uni_name <- all$unique_table
len <- nrow(all)

# do this one at a time b/c the tables are large
# and there's not enough space on the ec2 to store
# the data

# this pulls number of rows & lets us know if we
# need to chunk the data before pushing it
tbl <- paste0(uni_name[6], "_counts")
query <- paste0("SELECT COUNT(*) FROM ", tbl)
rows <- dbSendQuery(conn, query)
rows <- rows@Id[[1]]
max_rows <- max(rows)
rows <- seq(1, rows, by = 500000000)
rows <- sapply(1:length(rows), function(x) {
  if (rows[x] == min(rows)) {
    a <- c(rows[x], rows[x+1] - 1)
  } else if (x == max(length(rows))) {
    a <- c(rows[x], max_rows)
  } else {
    a <- c(rows[x], rows[x+1] - 1)
  }
})

if (max_rows > 500000000) {
    lapply(1:ncol(rows), function(x) {
        query <- paste0("SELECT * FROM ", tbl, " limit ", rows[1, x], ",", rows[2, x])
        t <- dbGetQuery(conn, query) %>% filter(counts != 0)
        print(nrow(t))
        fwrite(t, file = paste0("~/zero/", tbl, "_counts_zero_", x, ".csv"), row.names = F)
        
        if (x == 1) {
            dbWriteTable(conn, value = t, name = tbl, overwrite = TRUE, row.names = F)
        } else {
            dbWriteTable(conn, value = t, name = tbl, append = TRUE, row.names = F)
        }
    })
} else {
    t <- dbReadTable(conn, paste0(uni_name[1], "_counts")) %>% filter(counts != 0)
    fwrite(t, file = paste0("~/zero/", uni_name[1], "_counts_zero.csv"), row.names = F)

    dbWriteTable(conn, value = t, name = paste0(uni_name[1], "_counts"), overwrite = TRUE, row.names = F)
}


