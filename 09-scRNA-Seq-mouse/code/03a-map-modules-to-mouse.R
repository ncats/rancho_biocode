
#  [4. Mouse/Human Orthology with Phenotype Annotations (tab-delimited)](http://www.informatics.jax.org/downloads/reports/index.html#homology)
mouse_map <- read_tsv("iPSC-modules/HMD_HumanPhenotype.tsv", 
                      col_names = c("human", "ensg", "entrez", 
                                    "id", "mouse", "accession", "note")) 

df <- readxl::read_xlsx("iPSC-modules/iPSC_Profiler_curated_markers.xlsx") %>% 
  pivot_longer(cols = everything()) %>% 
  filter(!is.na(value)) %>% 
  left_join(mouse_map, by = c("value"="human"))

write_tsv(df, "iPSC-modules/iPSC-mouse-markers.tsv")

df %>% 
  group_by(name) %>% 
  summarise(ortholog_found = sum(!is.na(mouse)),
            no_ortholog = sum(is.na(mouse))) %>% 
  write_tsv("iPSC-modules/iPSC-count-of-mapped-homology.tsv")
