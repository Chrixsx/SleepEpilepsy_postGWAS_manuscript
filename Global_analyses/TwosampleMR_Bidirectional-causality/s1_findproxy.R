# Script to run and extract data from TwosampleMR R package
library(vroom)
library(tidyverse)
library(proxysnps)
library(LDlinkR)

# Step 1: Find instrumental SNPs


exposure_trait_name <- " " # Name of exposure trait
outcome_trait_name <- " " # Name of outcome trait

# Load GWAS sumstat of exposure and outcome traits
exposure_trait_gwas <- vroom::vromm("path_to_exposure_gwassumstats")
outcome_trait_gwas <- vroom::vroom("path_to_outcome_gwassumstats")

# Load lead SNPs of exposure trait (extracted by FUMA)
exposure_leadsnp <- vroom::vroom("path_to_exposure_leadSNPs")

# Calculate F-statistics for SNP instruments
gwas_size = 446118 # GWAS sample size

exposure_instrument <- exposure_trait_gwas %>% 
  filter(RSID %in% exposure_leadsnp$rsID) %>% 
  mutate(A1 = toupper(A1), A2 = toupper(A2),
         R2 = 2*BETA^2*FREQ1*(1-FREQ1)/(2*BETA^2*FREQ1*(1-FREQ1) + SE^2*2*gwas_size*FREQ1*(1-FREQ1)),
         F_stat = R2*(gwas_size - 1- 1)/(1-R2)*1) %>% 
  filter(F_stat >= 10) # only chose SNPs with F-stat > 10

# Find if chosen SNP instruments are in outcome_trait GWAS sumstat
available_snp <- outcome_trait_gwas %>% filter(RSID %in% exposure_instrument$RSID)

missing_snp <- setdiff(exposure_instrument$RSID, available_snp$RSID)

# For exposure SNP instrument that are not in outcome GWAS, find proxy SNPs
result <- data.frame(ori_rsID = character(), proxy_rsID = character(), stringsAsFactors = FALSE)


# Iterate over each missing SNPs
for (query_snp in missing_snp) {
  # Use tryCatch to handle errors during LDproxy call
  my_proxies <- tryCatch({
    # Call LDproxy to find proxies
    LDproxy(
      snp = query_snp,
      pop = "EUR",
      r2d = "r2",
      genome_build = "grch37",
      token = " " # Get access Token from https://ldlink.nih.gov/?tab=apiaccess
    ) %>% filter(R2 >= 0.8)
  }, error = function(e) {
    # If an error occurs, return NULL
    return(NULL)
  })

  # Check if proxies were retrieved successfully
  if (is.null(my_proxies) || nrow(my_proxies) == 0) {
    # If no proxies are found or an error occurred, record NA for proxy_rsID
    result <- rbind(result, data.frame(ori_rsID = query_snp, proxy_rsID = NA))
  } else {
    # Check if any proxy SNP is found in outcome gwas
    proxy_found <- FALSE
    for (i in seq_len(nrow(my_proxies))) {
      proxy_snp_rsid <- my_proxies$RS_Number[i]
      if (proxy_snp_rsid %in% outcome_trait_gwas$RSID) {
        result <- rbind(result, data.frame(ori_rsID = query_snp, proxy_rsID = proxy_snp_rsid))
        proxy_found <- TRUE
        break
      }
    }

    # If no matching proxy is found in outcome.gwas, record NA
    if (!proxy_found) {
      result <- rbind(result, data.frame(ori_rsID = query_snp, proxy_rsID = NA))
    }
  }
}


# Check if any of the found proxy is available in outcome_trait gwas sumstat
available_snp <- outcome_trait_gwas %>% filter(RSID %in% exposure_instrument$RSID)
missing_snp2 <- result$ori_rsID # Lost SNP, would include both No found and Proxy

proxy.rsid <-na.omit(result)$proxy_rsID

avai.data <-  exposure_leadsnp %>% 
  filter(RSID %in% available_snp) %>%
  mutate(avai.type = "TRUE")

proxy.data <- exposure_trait_gwas %>% 
  filter(RSID %in% proxy.rsid) %>% 
  mutate(A1 = toupper(A1), A2 = toupper(A2),
         R2 = 2*BETA^2*FREQ1*(1-FREQ1)/(2*BETA^2*FREQ1*(1-FREQ1) + SE^2*2*size*FREQ1*(1-FREQ1)),
         F_stat = R2*(size - 1- 1)/(1-R2)*1,
         avai.type = "TRUE,PROXY") %>% 
  filter(F_stat > 10)

# Add outcome (to trace back) because proxy SNP was found based on Exp vs Outcome gwas
twosmr_snp <- rbind(avai.data, proxy.data) %>%
  mutate(EXPOSURE = exposure_trait_name,
         OUTCOME = outcome_trait_name)

# Save outcome
vroom_write(twosmr_snp, file.path("path_storage_SNPinstruments", paste0(exposure_trait_name,"_vs_", outcome_trait_name,"_with_proxy.txt")), delim="\t")
