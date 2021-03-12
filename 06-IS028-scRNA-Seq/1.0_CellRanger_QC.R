#load Libraries
library(openxlsx)
library(tidyverse)

#CONSTANTS

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, face = "bold", size = 20),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(face = "bold", colour = "Black", size = 16),
                 legend.text = element_text(colour = "Black", size = 12),
                 legend.title = element_text(colour = "Black", size = 14))


##################################


#load data
CRsummary <- read.xlsx("./results/cellRangerSummary.xlsx")

#separate filenames into columns to get isolation methods and introns
CRsummary <- CRsummary %>%
  mutate(Sample = gsub("Lonza-Noci-D32-", "", Sample)) %>%
  separate(Sample, sep = "_", into = c("Iso_Method", "Introns"), fill = "right") %>%
  replace_na(list(Introns = "NoIntrons"))


#plot estimated number of cells########################################################

p <- ggplot(CRsummary, aes(x = Iso_Method, y = Estimated.Number.of.Cells, fill = Introns)) + 
        geom_col(color = "black", position = "dodge") + 
        scale_y_continuous(expand = c(0,0)) + 
        labs(x = "Isolation Method",
             y = "Number of Cells",
             title = "Number of Cells For Each Isolation Method") +
        theme_classic() + 
        mytheme

p


plotName <- "./figures/NumberCells.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")


#plot reads per cell########################################################

p <- ggplot(CRsummary, aes(x = Iso_Method, y = Mean.Reads.per.Cell, fill = Introns)) + 
  geom_col(color = "black", position = "dodge") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Isolation Method",
       y = "Number of Reads",
       title = "Number of Reads Per Cell") +
  theme_classic() + 
  mytheme

p


plotName <- "./figures/ReadsPerCell.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")



#plot Genes per cell########################################################

p <- ggplot(CRsummary, aes(x = Iso_Method, y = Median.Genes.per.Cell, fill = Introns)) + 
  geom_col(color = "black", position = "dodge") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Isolation Method",
       y = "Median Number of Genes",
       title = "Median Number of Genes Per Cell") +
  theme_classic() + 
  mytheme

p


plotName <- "./figures/GenesPerCell.png"
ggsave(filename = plotName, 
       device = "png", 
       plot =  p, 
       height = 6, 
       width = 8, 
       units = "in")


#clean up workspace
rm(list = ls())
gc()
