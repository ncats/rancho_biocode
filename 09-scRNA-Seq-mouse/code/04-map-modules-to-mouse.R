
#  [4. Mouse/Human Orthology with Phenotype Annotations (tab-delimited)](http://www.informatics.jax.org/downloads/reports/index.html#homology)
mouse_map <- read_tsv("ncats_09_scRNA_seq_mouse/reference/HMD_HumanPhenotype.tsv", 
                      col_names = c("human", "ensg", "entrez", 
                                    "id", "mouse", "accession", "note")) 

df <- readxl::read_xlsx("reference/iPSC_Profiler_curated_markers.xlsx") %>% 
  pivot_longer(cols = everything()) %>% 
  filter(!is.na(value)) %>% 
  left_join(mouse_map, by = c("value"="human"))

write_tsv(df, "ncats_09_scRNA_seq_mouse/reference/iPSC-mouse-markers.tsv")

df %>% 
  group_by(name) %>% 
  summarise(ortholog_found = sum(!is.na(mouse)),
            no_ortholog = sum(is.na(mouse))) %>% 
  write_tsv("ncats_09_scRNA_seq_mouse/reference/iPSC-count-of-mapped-homology.tsv")

# convert new Satellite Glia markers from Pei-hsuan
read_tsv("ncats_09_scRNA_seq_mouse/reference/SatGlia.tsv") %>% 
  left_join(mouse_map, by = c("human_satglia"="human")) %>% 
  filter(!is.na(mouse)) %>% 
  write_tsv("ncats_09_scRNA_seq_mouse/reference/satglia-mapped-markers.tsv")
