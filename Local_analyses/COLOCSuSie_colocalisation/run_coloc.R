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
