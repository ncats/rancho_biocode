library(RMySQL)
library(data.table)
library(tidyverse)

options(scipen=999)

# create conn to SQL db SCDB
conn <- dbConnect(MySQL(), user = 'sc-app-admin', 
                  password = 'Bisi9adm1n', 
                  dbname = 'SCDB', 
                  host = 'sctl-dev-rds-mysql-01.ceyknq0yekb3.us-east-1.rds.amazonaws.com', 
                  port = 3306)

# split out huge scRNA-seq datasets into columns
query <- 'select * from SCDB.IS010_counts where sample = "H9_CE"'
h9_ce <- dbGetQuery(conn, query)
write.csv(h9_ce, file = "~/is010_h9_ce.csv", row.names = F)

query <- 'select * from SCDB.IS010_counts where sample = "H9_CEPT"'
h9_cept <- dbGetQuery(conn, query)
write.csv(h9_cept, file = "~/is010_h9_cept.csv", row.names = F)

query <- 'select * from SCDB.IS010_counts where sample = "H9_Y27632"'
h9_y27 <- dbGetQuery(conn, query)
write.csv(h9_y27, file = "~/is010_h9_y27.csv", row.names = F)

query <- 'select * from SCDB.IS010_counts where sample = "Lonza_CE"'
ipsc_ce <- dbGetQuery(conn, query)
ipsc_ce <- ipsc_ce %>%
mutate(cell_id = gsub("Lonza", "LiPSC.GR1.1", cell_id),
sample = gsub("Lonza", "LiPSC.GR1.1", sample))
write.csv(ipsc_ce, file = "~/is010_ipsc_ce.csv", row.names = F)

query <- 'select * from SCDB.IS010_counts where sample = "Lonza_CEPT"'
ipsc_cept <- dbGetQuery(conn, query)
ipsc_cept <- ipsc_cept %>%
mutate(cell_id = gsub("Lonza", "LiPSC.GR1.1", cell_id),
sample = gsub("Lonza", "LiPSC.GR1.1", sample))
write.csv(ipsc_cept, file = "~/is010_ipsc_cept.csv", row.names = F)

query <- 'select * from SCDB.IS010_counts where sample = "Lonza_Y27632"'
ipsc_y27 <- dbGetQuery(conn, query)
ipsc_y27 <- ipsc_y27 %>%
mutate(cell_id = gsub("Lonza", "LiPSC.GR1.1", cell_id),
sample = gsub("Lonza", "LiPSC.GR1.1", sample))
write.csv(ipsc_y27, file = "~/is010_ipsc_y27.csv", row.names = F)

tab <- fread("~/is010_h9_ce.csv")
dbWriteTable(conn, value = tab, name = 'IS010_counts', overwrite = TRUE, row.names = F)
tab <- fread("~/is010_h9_cept.csv")
dbWriteTable(conn, value = tab, name = 'IS010_counts', append = TRUE, row.names = F)
tab <- fread("~/is010_h9_y27.csv")
dbWriteTable(conn, value = tab, name = 'IS010_counts', append = TRUE, row.names = F)
tab <- fread("~/is010_ipsc_ce.csv")
dbWriteTable(conn, value = tab, name = 'IS010_counts', append = TRUE, row.names = F)
tab <- fread("~/is010_ipsc_cept.csv")
dbWriteTable(conn, value = tab, name = 'IS010_counts', append = TRUE, row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "NCCIT"'
nccit <- dbGetQuery(conn, query)
write.csv(nccit, file = "~/is017_nccit.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "NTERA2"'
ntera2 <- dbGetQuery(conn, query)
write.csv(ntera2, file = "~/is017_ntera2.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "LONZA"'
lonza <- dbGetQuery(conn, query)
lonza <- lonza %>%
mutate(cell_id = gsub("Lonza", "LiPSC.GR1.1", cell_id),
sample = gsub("Lonza", "LiPSC.GR1.1", sample))
write.csv(lonza, file = "~/is017_lonza.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "PA-1"'
pa1 <- dbGetQuery(conn, query)
write.csv(pa1, file = "~/is017_pa1.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "GM23720"'
gm23720 <- dbGetQuery(conn, query)
write.csv(gm23720, file = "~/is017_gm23720.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "GM23279"'
gm23279 <- dbGetQuery(conn, query)
write.csv(gm23279, file = "~/is017_gm23279.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "WA26"'
wa26 <- dbGetQuery(conn, query)
write.csv(wa26, file = "~/is017_wa26.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "WA14"'
wa14 <- dbGetQuery(conn, query)
write.csv(wa14, file = "~/is017_wa14.csv", row.names = F)

query <- 'select * from SCDB.IS017_counts where sample = "NCRM5"'
ncrm5 <- dbGetQuery(conn, query)
write.csv(ncrm5, file = "~/is017_ncrm5.csv", row.names = F)

tab <- fread("~/is017_nccit.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', overwrite = TRUE, row.names = F)
tab <- fread("~/is017_ntera2.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_lonza.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_pa1.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_gm23720.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_gm23279.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_wa26.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_wa14.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)
tab <- fread("~/is017_ncrm5.csv")
dbWriteTable(conn, value = tab, name = 'IS017_counts', append = TRUE, row.names = F)

query <- 'select * from SCDB.IS027_Tao_counts LIMIT 0,200000000'
a <- dbGetQuery(conn, query)
write.csv(a, file = "~/is027_a.csv", row.names = F)

query <- 'select * from SCDB.IS027_Tao_counts LIMIT 200000000,400000000'
b <- dbGetQuery(conn, query)
write.csv(b, file = "~/is027_b.csv", row.names = F)

query <- 'select * from SCDB.IS027_Tao_counts LIMIT 400000000,600000000'
c <- dbGetQuery(conn, query)
c <- c %>%
mutate(cell_id = gsub("Lonza", "LiPSC.GR1.1", cell_id),
sample = gsub("Lonza", "LiPSC.GR1.1", sample))
write.csv(c, file = "~/is027_c.csv", row.names = F)


tab <- fread("~/is027_a.csv")
dbWriteTable(conn, value = tab, name = 'IS027_Tao_counts', overwrite = TRUE, row.names = F)
tab <- fread("~/is027_b.csv")
dbWriteTable(conn, value = tab, name = 'IS027_Tao_counts', append = TRUE, row.names = F)
tab <- fread("~/is027_c.csv")
dbWriteTable(conn, value = tab, name = 'IS027_Tao_counts', append = TRUE, row.names = F)


a <- fread("/opt/sctl-data-local/adj_data/sc/IS009_counts.csv")
a <- a %>%
mutate(cell_id = gsub("lonza", "LiPSC.GR1.1", cell_id),
sample = gsub("lonza", "LiPSC.GR1.1", sample))
write_csv(a, file = "~/is009.RDS")

tab <- fread("~/is009.csv")
dbWriteTable(conn, value = tab, name = 'IS009_counts', append = TRUE, row.names = F)
