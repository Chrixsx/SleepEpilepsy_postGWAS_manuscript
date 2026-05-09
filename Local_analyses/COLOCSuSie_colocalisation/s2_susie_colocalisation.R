library(vroom)
library(tidyverse)
library(coloc)
library(Rfast) # Make Susie run faster
library(ieugwasr)

# Import GWAS summary statistics for two traits
trait_1_gwassumstat <- vroom::vromm("path_to_trait1_gwassumstats")
trait_2_gwassumstat <- vroom::vroom("path_to_trait2_gwassumstats")

# Import information of the testing region
# Information of testing regions (chr, start, stop) are extracted from LAVA significant regions
# For example, here I tested chr2:191051955-193033982 region
chr <- 2
start <- 191051955
stop <- 193033982


# To run SuSiE, we need a Linkage disequilibrium matrix of shared SNPs.
# To do that, we need to find:
  # (1) shared SNPs between two traits, (2)
  # (2) extract the shared SNP from LD reference panel. However, not all shared SNPs will available in LD reference panel
  # (3) therefore, we need to do 3-way intersect (trait1-trait2-LD reference) 
  # (4) Once the above steps are done, prepare data and run SuSiE 


# (1) Find shared SNPs from two GWAS sumstats of two traits
shared_snps <- extract_shared_snps_in_region(trait_1_gwassumstat, trait_2_gwassumstat, chr, start, stop)
print(paste0("There are ", shared_snps %>% length(), " SNPs in region ", chr, "_", start, "_", stop))


# (2) Extract available shared SNPs from LD reference panel
# Reformat shared SNPs from Trait 1  
trait_1_snp <- trait_1_gwassumstat %>% 
  filter(MarkerName %in% shared_snps) %>% 
  distinct(MarkerName, .keep_all = TRUE) %>%
  
  rename(trait1_Allele1 = Allele1,
          trait1_Allele2 = Allele2,
          trait1_Beta = Beta,
          trait1_Var = Var,
          trait1_Chr = CHR,
          trait1_BP = BP) %>% 
  
  mutate(trait1_Allele1 = toupper(trait1_Allele1),
         trait1_Allele2 = toupper(trait1_Allele2)) %>% 
  
  arrange(trait1_Chr, trait1_BP)


# Reformat shared SNPs from Trait 2  
trait_2_snp <- trait_2_gwassumstat %>% 
  filter(MarkerName %in% shared_snps) %>% 
  distinct(MarkerName, .keep_all = TRUE) %>%
  
  rename(trait2_Allele1 = Allele1,
          trait2_Allele2 = Allele2,
          trait2_Beta = Beta,
          trait2_Var = Var,
          trait2_Chr = CHR,
          trait2_BP = BP) %>% 
  
  mutate(trait2_Allele1 = toupper(trait2_Allele1),
         trait2_Allele2 = toupper(trait2_Allele2)) %>% 
  
  arrange(trait2_Chr, trait2_BP)



# Align the order of those common SNPs in both traits, else, it will create error for LD matrix
# Check if the rsID in the same order between 2 traits
identical(trait_1_snp$MarkerName, trait_2_snp$MarkerName) # return TRUE/FALSE
re_shared_snps <- trait_1_snp$MarkerName # Extracted rsID in specific order


# Create LD matrix according to the ordered SNP
plink_path <- "path_to_plink_file/plink.exe"
bfile_path <- "path_to_1000G_LDpanel/g1000_eur"

# Create LD matrix
ld.matrix <- ieugwasr::ld_matrix(
    variants = re_shared_snps,
    with_alleles = FALSE,
    pop = "EUR",
    plink_bin = plink_path,
    bfile = bfile_path)

# (3) Checking LD matrix
# Not all SNPs are captured in LD matrix
snp_in_LD.matrix <- colnames(ld.matrix)
share_snps=intersect(snp_in_LD.matrix, re_shared_snps)
non_shared_snps <- setdiff(re_shared_snps, snp_in_LD.matrix)

# Filter SNP again
trait_1_snp <- trait_1_snp %>% filter(!MarkerName %in% non_shared_snps)
trait_2_snp <- trait_2_snp %>% filter(!MarkerName %in% non_shared_snps)


# (4) Prepare data and run SuSiE to find shared causal SNPs
dataset1 <- list(
    beta = trait_1_snp$trait1_Beta,
    varbeta= trait_1_snp$trait1_Var,
    snp = trait_1_snp$MarkerName,
    position= trait_1_snp$trait1_BP,
    type = "quant", # if trait is case-control (use cc), if continuous (use quant)
    sdY=1,
    LD= ld.matrix,
    N=as.integer("403195") # sample size of trait 1 GWAS
    )

dataset2 <- list(
    beta = trait_2_snp$trait2_Beta,
    varbeta= trait_2_snp$trait2_Var,
    snp = trait_2_snp$MarkerName,
    position= trait_2_snp$trait2_BP,
    type = "cc", # if trait is case-control (use cc), if continuous (use quant)
    LD= ld.matrix,
    N=as.integer("49388") # sample size of trait 2 GWAS
    )

# Sense check the dataset (optional -- but recommended)
check_dataset(dataset1, req="LD")
check_dataset(dataset2, req="LD") 


# Extract credible SNP from each trait
trait1_credible= runsusie(dataset1,estimate_prior_variance = FALSE)
summary(trait1_credible)

trait2_credible= runsusie(dataset2,estimate_prior_variance = FALSE)
summary(trait2_credible)

# Check if any of the credit SNPs are shared between two traits
pairwise_susie= coloc.susie(trait1_credible,trait2_credible)

# Export result
vroom_write(pairwise_susie$results, 
            file= file= paste0("trai1_trait2_chr2_191051955_193033982_SUSIE.txt"),
            delim="\t")
