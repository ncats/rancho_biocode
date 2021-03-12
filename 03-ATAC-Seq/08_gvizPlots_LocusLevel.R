#load libraries
library(Gviz)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)
library(data.table)
library(Cairo)
library(magrittr)
library(biomaRt)
library(tidyverse)
library(extrafont)

#get functions
source("./code/functions/build_locus_zoom_DAG.R")


##CONSTANTS

#create a list of the contrasts to loop through
contrasts <- c("A1_vs_D0", "D30_vs_A1", "D30_vs_D0", "D50_vs_D30", "D50_vs_D0", "LSB_vs_A1", "LSB_vs_D0")

#read the results file
dfall <- read.delim("./results/03_NCATS_ATAC_AllCombinedResults.txt", header = T, sep = "\t")

#get annotation
ann <- read.delim("./data/refs/NCATS_annot.txt", header = T, sep = "\t")

all_dtrack_cont <- readRDS("~/Desktop/all_dtrack_cont.RDS")

for(i in 1:length(contrasts)) {
  
  #filter data for only significant genes and clean up the contrast column
  df <- dfall %>%
    filter(Contrast == contrasts[i],
           Significant == 1)
  
  #split the contrast name into the two treatments    
  contrast <- as.list(strsplit(contrasts[i], "_vs_"))
  
  #identify genes of interest
  #top 5 genes
  goi <- df %>%
    arrange(AdjP.ScoreInv, Distance) %>%
    distinct(Chr, Start, End, .keep_all = T) %>%
    head(n=5) %>%
    arrange(Gene, Start)
  
  #filter annotation
  annot <- ann %>%
    filter(Gene %in% goi$Gene) %>%
    arrange(Gene, Start)
  
  #increase the range of annotation using df to include peaks just outside of gene region
  annot <- df %>%
    dplyr::select(Gene, Distance) %>%
    group_by(Gene) %>%
    mutate(maxUp = min(Distance),
           maxDown = max(Distance)) %>%
    ungroup() %>%
    left_join(annot, ., by = "Gene") %>%
    mutate(maxUp = ifelse(maxUp > 0, 0, maxUp),
           maxDown = ifelse(maxDown < 0, 0, maxDown),
           #expand the window by 500 in both directions to make the plotting look nicer
           Start = Start + maxUp - 500,
           End = End + maxDown + 500) %>%
    dplyr::select(-Distance, -maxUp, -maxDown) %>%
    distinct(.keep_all = T)
  
  
  ####### Ideogram track #######
  itrack <- lapply(1:nrow(annot), function(x) 
    ideo_track(gen = 'hg38',
               chr = annot[x,1])
  )
  
  ########## Annotation Track ############
  # extract MergedReg (medip-seq peaks) from 
  # filtered bed file
  all_peaks <- data.frame(df)
  
  peaks <- lapply(1:nrow(annot), function(x)
    all_peaks %>%
      filter(as.character(Chr) == as.character(annot$Chr[x])) %>%
      #filter(Start >= annot$Start[x]) %>%
      #filter(End <= annot$End[x]) %>%
      filter(Gene %in% annot$Gene[x]) %>%
      mutate(Strand = annot$Strand[x])
  )
  
  nam <- lapply(1:length(peaks), function(x) {
    y <- peaks[[x]] %>%
      pull(Gene)
    y <- unlist(sapply(1:length(y), function(z)
      paste(y[[z]], collapse = "\n")
    ))
    return(y)
  })
  
  aTrack <- lapply(1:length(peaks), function(x) 
    aTrack <- anno_track(st = peaks[[x]]$Start,
                         chr = peaks[[x]]$Chr[1],
                         plot_name = nam[[x]],
                         #plot_name = "", # if you want to remove the names 
                         wid = abs(peaks[[x]]$Start - peaks[[x]]$End),
                         gen = "hg38",
                         title = "ATAC Peaks",
                         feat = "ATAC",
                         len = length(peaks[[x]]$Start))
  )
  
  
  ########## Genome Axis Track ############
  range <- lapply(1:nrow(annot), function(x) {
    range <- annot[x,] %>%
      dplyr::slice(c(1, n())) %>%
      dplyr::select(Start, End)
    range <- c(range$Start[1], range$End[2])
  })
  
  gtrack <- lapply(1:length(range), function(x)
    GenomeAxisTrack(#\chr = peaks[[x]]$chr[1],
      range=IRanges(start = range[[x]][1],
                    end = range[[x]][2])
    )
  )
  
  
  
  
  
  ####### Gene_Region_Track ######
  bmTrack <- lapply(1:nrow(annot), function(x)
    BiomartGeneRegionTrack(genome = "hg38",
                           chromosome = annot[x,"Chr"], 
                           start = annot[x,"Start"], 
                           end = annot[x,"End"],
                           #name = "ENSEMBL",
                           name = "", #to shrink the track
                           collapseTranscripts = "meta",
                           filters = list(ensembl_gene_id = annot[x,"ENSEMBL"]),
                           #using fill color from UCSC genome browser
                           fill = "#0C0C78",
                           col = "#0C0C78")
  )
  
  
  ########## Data Track ############
  # # pull coverage from bam files
  # bam_files <- list.files("./bam/")
  # bam_files <- bam_files[grepl("\\.bam$", bam_files)]
  # bam_files <- paste("./bam/", bam_files, sep = "")
  # 
  # 
  # #change the - to an _ to seperate D30 and D50 samples from NCRM5
  # #this will be changed back later
  # bam_files <- gsub("-D30", "_D30", bam_files)
  # bam_files <- gsub("-D50", "_D50", bam_files)
  # 
  # #create a name for the bam based on the sample to be used 
  # #as a title for the track in the visualization
  # bam_nam <- gsub(".*JNIH_", "", bam_files)
  # bam_nam <- gsub("_ATAC_.*", "", bam_nam)
  # 
  # 
  # #reduce files
  # contrast_files <- bam_files[grep(paste(contrast[[1]],"-", sep = "", collapse = "|"), bam_files)]
  # contrast_nam <- bam_nam[grep(paste(contrast[[1]],"-", sep = "", collapse = "|"), bam_nam)]
  # 
  # #changing this back to access correct files
  # contrast_files <- gsub("_D30", "-D30", contrast_files)
  # contrast_files <- gsub("_D50", "-D50", contrast_files)
  # 
  # 
  # # this has been reduced to only the target contrast, can do 1:6
  # all_cont <- lapply(1:length(range), function(x) {
  #   unlist(lapply(c(1:6), function(y) {
  #     build_cov(chr = peaks[[x]]$Chr[1],
  #               range1 = range[[x]][1],
  #               range2 = range[[x]][2],
  #               file = contrast_files[y])
  #   }), recursive = F)
  # })
  
  ####### Combine Tracks #######
  # all_dtrack_cont <- lapply(1:length(all_cont), function(x) {
  #   #get max coverage for consistent scales
  #   max_cov <- c()
  #   for(i in 1:length(all_cont[[x]])){
  #     tmp_max <- max(all_cont[[x]][[i]]$coverage)
  #     max_cov <- append(max_cov, tmp_max)
  #   }
  #   max_cov <- max(max_cov)
  #   ylims <- range(0, max_cov)
  #   #create data track
  #   unlist(lapply(c(1:6), function(y) 
  #     dat_track(all_cont[[x]][[y]],
  #               gen = "hg38",
  #               style = "histogram",
  #               title = contrast_nam[y],
  #               chr = gsub("chr", "", peaks[[x]]$Chr[1]),
  #               ylim = ylims
  #     )
  #   ), recursive = F)
  # })
  
  dtrack <- all_dtrack_cont[[i]]
  
  ########## plot all data! ############
  ########## combine itrack, atrack and data tracks
  ########## into a single png
  first <- lapply(1:length(itrack), function(x)
    list(itrack[[x]], aTrack[[x]], bmTrack[[x]]))
  
  comb_cont <- lapply(1:length(dtrack), function(x)
    c(first[[x]], dtrack[[x]], gtrack[[x]]))
  
  #set a contrast name
  cont <- paste(contrast[[1]], collapse = " vs ")
  contfn <- paste(contrast[[1]], collapse = "_vs_")
  
  title_cont <- lapply(1:length(itrack), function(x)
    paste(cont, " at ", annot$Gene[x], " Locus", sep = ""))
  
  plot_nam_cont <- lapply(1:nrow(annot), function(x)
    paste("./figures/gviz/Locus/Locus_", contfn , "_", as.character(annot$Gene[x]), ".pdf", sep = "")
  )
  
  lapply(1:length(comb_cont), function(x) {
    cairo_pdf(filename = plot_nam_cont[[x]], 
              width = 10, height = 12, 
              family = "NotoSans-Condensed")
    #CairoPDF(file = plot_nam_cont[[x]], width = 10, height = 15, family = "NotoSans-Condensed")
    #reordering the tracks so the genome access track is near the top
    plotTracks(comb_cont[[x]][c(1, 10, 2:9)],
               main = title_cont[[x]],
               # set fill for y-axis (title)
               background.title = "white",
               # set font color for y-axis (title)
               col.title = "black",
               # set font family for y-axis (title)
               fontfamily.title = "NotoSans-Condensed",
               fontfamily = "NotoSans-Condensed",
               # set font size for y-axis (title)
               fontsize = 12,
               #sizes=c(0.5, 0.5, 0.5, 0.7, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75),
               sizes=c(0.5, 0.5, 1, 0.3, 1, 1, 1, 1, 1, 1),
               from = annot$Start[x], 
               to = annot$End[x] 
    )
    graphics.off()
  })
  
}
