library(vroom)
library(tidyverse)
library(coloc)


# HELPER FUNCTION 1: Define scanning region
# Based on LAVA bivariate output to define scanning regions between pair of phenotypes
extract_shared_snps_in_region <- 
  function(gwas_trait_1, gwas_trait_2, chr, start, stop) {
    
  shared_snps <- intersect(gwas_trait_1$MarkerName, gwas_trait_2$MarkerName) # Get shared SNPs
  
  # Select SNPs from trait 1
  trait_1_snps <- gwas_trait_1 %>% 
    filter(CHR == as.numeric(chr), 
           BP >= as.numeric(start), BP < as.numeric(stop))   
  
  # Select SNPs from trait 2
  trait_2_snps <- gwas_trait_2 %>% 
    filter(CHR == as.numeric(chr), 
           BP >= as.numeric(start), BP < as.numeric(stop))
  
  # Get unique shared SNPs
  shared_snps <- intersect(trait_1_snps$MarkerName, trait_2_snps$MarkerName) %>% unique() 
  
  return(shared_snps)
}

# HELPER FUNCTION 2: Run Coloc
run_coloc <- 
  function(gwas_trait_1, gwas_trait_2, chr, start, stop) {
  
  # Extract shared snps between two traits
  shared_snps <- extract_shared_snps_in_region(gwas_trait_1, gwas_trait_2, chr, start, stop)
  
  print(paste0("There are ", shared_snps %>% length(), " SNPs in region ", chr, "_", start, "_", stop))
  
  trait_1_snp <- gwas_trait_1 %>% 
    filter(MarkerName %in% shared_snps) %>% 
    distinct(MarkerName, .keep_all = TRUE)
  
  trait_2_snp <- gwas_trait_2 %>% 
    filter(MarkerName %in% shared_snps) %>% 
    distinct(MarkerName, .keep_all = TRUE)
  
  # Prepare dataset 1 (from Trait 1) for COLOC
  dataset1 <- list(
    beta = trait_1_snp$Beta,
    varbeta= trait_1_snp$Var,
    snp = trait_1_snp$MarkerName,
    position= trait_1_snp$BP,
    type = "cc" # case-control type
    ) 
  
  # Prepare dataset 2 (from Trait 2) for COLOC
  dataset2 <- list(
    beta = trait_2_snp$Beta,
    varbeta= trait_2_snp$Var,
    snp = trait_2_snp$MarkerName,
    position= trait_2_snp$BP,
    type = "cc" 
    ) 
  
  # Sense check (optional-- but recommended)
  check_dataset(dataset1)
  check_dataset(dataset2)
  
  par(mfrow= c(2,1))
  par(mar= c(4,4,1,2))
  plot_dataset(dataset1)
  plot_dataset(dataset2)
  
  # Run colocalisation
  my.res <- coloc.abf(dataset1, dataset2)
  
  return(my.res$summary)
}

###########

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
