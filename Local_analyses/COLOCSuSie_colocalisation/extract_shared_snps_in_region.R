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
