library(vroom)
library(tidyverse)
library(coloc)


# Apply helper function 1 and 2 to actual COLOC colocalisation
trait_1_gwassumstat <- vroom::vromm("path_to_trait1_gwassumstats")
trait_2_gwassumstat <- vroom::vroom("path_to_trait2_gwassumstats")

# Information of testing regions (chr, start, stop) are extracted from LAVA significant regions
# For example, here I tested chr2:191051955-193033982 region
coloc_outcome <- run_coloc(
  gwas_trait_1 = trait_1_gwassumstat,
  gwas_trait_2 = trait_2_gwassumstat,
  chr = 2,
  start = 191051955,
  stop = 193033982)


# Export COLOC outcome
write.table(coloc_outcome, 
            file= paste0("trai1_trait2_chr2_191051955_193033982_COLOC.txt"), 
            sep= "\t", 
            row.names = TRUE, quote= FALSE)
