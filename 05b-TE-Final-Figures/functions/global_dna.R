#' global dna methylation functions for building
#' plots

# bin methylation data (CpG islands, TSS, and GB's)
bin_methyl <- function(file, bin) {
  df <- file
  cut <- as.character(unique(cut_interval(df$start:df$end, bin, dig.lab=10)))
  cut <- strsplit(gsub("\\(|\\]|\\[", "", cut), ",")
  
  if (ncol(df) == 5) {
    cut <- rbindlist(lapply(1:length(cut), function(x)
      data.frame("peak" = df$peak,
                 "chr" = df$chr,
                 "start" = round(as.numeric(cut[[x]][1])),
                 "end" = round(as.numeric(cut[[x]][2])),
                 "Count" = file[1, 5])
    ))
  } else {
    cut <- rbindlist(lapply(1:length(cut), function(x)
      data.frame("peak" = df$peak,
                 "chr" = df$chr,
                 "start" = round(as.numeric(cut[[x]][1])),
                 "end" = round(as.numeric(cut[[x]][2])))
    ))
  }
}

# re-format methylation data for plotting
counts_reformat <- function(counts, bin, bam_nam_adj, split) {
  y <- counts
  z <- lapply(1:length(y), function(t) {
    b <- data.frame(y[[t]])
  })
  z <- do.call("cbind", z)
  names(z) <- bam_nam_adj
  
  st <- data.frame(bin)
  z <- cbind(st, z)
  z <- z %>%
    mutate(cut = as.numeric(cut_number(z$start, split))) %>%
    gather(sample, n, -seqnames, -start, -end, -width, -strand, -cut) %>%
    group_by(cut, sample) %>%
    summarise(n = mean(n, na.rm = T))
}

# plot global methylation plots
plot_methyl <- function(file, sp, brks, int, labs, name, filename, h, w, k) {
  p <- ggplot(file, aes(x = cut, y = n, color = `Sample-replicate`)) +
    stat_smooth(method = "gam",  
                formula = y ~ s(x, k = k), 
                n = 365,
                se = TRUE, 
                aes(fill = `Sample-replicate`), 
                #alpha = 0.3
                ) + 
    scale_x_continuous(breaks = brks,
                       labels = labs) +
    scale_fill_manual(values = c("firebrick1", "tomato", "lightsalmon",
                                 "blue1", "deepskyblue", "lightblue1")) +
    scale_color_manual(values = c("firebrick1", "tomato", "lightsalmon",
                                  "blue1", "deepskyblue", "lightblue1")) +
    xlab(name) + ylab("Normalized Counts") +
    theme_classic() +
    mytheme + 
    geom_vline(xintercept = int, linetype="dotted", 
               color = "black", size = 0.8)
  ggsave(filename = paste0("./adj_data/plots/", filename, "_global_methylation_plot.png"), 
         units = "in",
         height = h, width = w, plot = p)
  pdf(file = paste0("./adj_data/plots/", filename, "_global_methylation_plot.pdf"),
                    height = h, width = w)
  print(p)
  graphics.off()
}
