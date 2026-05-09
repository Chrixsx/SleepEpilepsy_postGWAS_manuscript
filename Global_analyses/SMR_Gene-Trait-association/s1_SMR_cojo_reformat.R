# Step 1: Reformatting GWAS to COJO format
library(vroom)
library(tidyverse)

gwas_file <- " " # full path to gwas sumstats
output_dir <- " " # folder to store reformatted gwas


original_gwas <- vroom(gwas_file, show_col_types = F, delim= "\t")

cojo_gwas <- original_gwas %>%
  select(c("MarkerName", "Allele1", "Allele2", "Freq1", "Beta", "SE", "P-value", "Effective_N", "BP",  "CHR")) %>%
  rename(
    SNP = "MarkerName",
    A1 = "Allele1",
    A2 = "Allele2",
    freq = "Freq1",
    b = "Beta",
    se = "SE",
    p = "P-value",
    N = "Effective_N",
    pos = "BP",
    chr = "CHR"
  )

vroom_write(x = cojo_gwas, file = file.path(output_dir, "gwas_COJO_reformatted.tsv.gz"), delim= "\t")

# Step 2: Upload GWAS summary statistic in COJO format to https://yanglab.westlake.edu.cn/smr-portal/ 
