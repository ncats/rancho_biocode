#' Custom fxn to remove medip-seq peaks > 10kb 
#' upstream/downstream

# read in deseq files, filter padj < 0.05;
# keep all cols except counts; separate out
# medip-seq peak row annotation into individ-
# ual rows; split into up/downstream annotated
# genes; filter out medip-seq peaks up/down-
# stream of 10 kb; also, summarize count of
# hypo/hyper methylated peaks in ea contrast
remove_10kb <- function(file, contrast, file_type) {
  # filter padj < 0.05
  y <- fread(file) %>%
    filter(padj < 0.05)
  
  # filter out count data
  sub <- y %>%
    dplyr::select(c(1:20, ncol(y))) %>%
    data.frame(.)
  
  # pull location, distance from each medip-seq
  # peak and put on sep row in df
  expand_df <- function(idx) {
    loc <- gsub("\\;.*", "", sub$location[idx])
    loc <- unlist(strsplit(loc, ", "))
    dis <- gsub("\\;.*", "", sub$distance_to_start[idx])
    dis <- unlist(strsplit(dis, ", "))
    
    z <- data.frame(`MergedRegion` = rep(sub$MergedRegion[idx], length(dis)),
                    `chromosome_name_medip` = rep(sub$chromosome_name_medip[idx], length(dis)),
                    `start_position_medip` = rep(sub$start_position_medip[idx], length(dis)),
                    `end_position_medip` = rep(sub$end_position_medip[idx], length(dis)),
                    `location` = loc,
                    `distance` = dis
    )
  }
  
  # expand df & collapse into df
  m <- rbindlist(lapply(1:nrow(sub), function(x) expand_df(x)))
  
  # sep df into upstream elements; 
  # filter out medip-seq peaks 10 kb
  # upstream
  up <- m %>% 
    filter(location == "upstream") %>%
    mutate(distance = as.numeric(as.character(distance))) %>%
    filter(distance < -10000)
  
  # sep df into downstream elements; 
  # filter out medip-seq peaks 10 kb
  # downstream
  down <- m %>% 
    filter(location == "downstream") %>%
    mutate(distance = as.numeric(as.character(distance))) %>%
    filter(distance > 10000)
  
  remove <- rbind(up, down)
  remove <- remove %>%
    mutate(comb = paste(MergedRegion, chromosome_name_medip, start_position_medip,
                        end_position_medip, location, distance, sep = "_"))
  
  # filter out 10 kb up/downstream medip-seq peaks
  sub <- sub %>%
    mutate(comb = paste(MergedRegion, chromosome_name_medip, start_position_medip,
                        end_position_medip, location, distance_to_start, sep = "_")) %>%
    filter(!comb %in% remove$comb)
  
  write.csv(sub, file = paste("./adj_data/deseq_10kb/", contrast, "_", file_type, "_10kb_up_down_remove.csv", sep = ""), row.names = F)
  return(sub)
}