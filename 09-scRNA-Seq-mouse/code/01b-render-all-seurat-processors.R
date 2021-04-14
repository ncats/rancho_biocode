
rend <- function(sample){
  rmarkdown::render( 'seurat-processor.Rmd', 
                     params = list(sample=sample),
                     knit_root_dir = "/home/danr/ncats-09_scRNA-Seq-mouse",
                     output_file = glue::glue("filtering-{sample}.html")) }

lapply(c('082120DRG','101920DRG','103020DRG','110220DRG','111220DRG'), rend)

