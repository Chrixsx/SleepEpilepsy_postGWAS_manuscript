# Script to run and collate data/image from Twosample Mendelian randomisation

library(vroom)
library(tidyverse)
library(TwoSampleMR)
library(gridExtra)
library(ggplot2)

result_dir <- " " # path to results storage folder
plot_dir <- " " # path to images/plots storage folder
exposure_trait_name <- " " # Name of exposure trait
outcome_trait_name <- " " # Name of outcome trait


twosmr_snp_file <- paste0(expos,"_vs_", outcome.trait,"_with_proxy.txt")


#Step 1: Prepare Data
exposure.dat= read_exposure_data(filename = file.path("path_storage_SNPinstruments", twosmr_snp_file),
                                 sep ="\t",
                                 phenotype_col = "EXPOSURE",
                                 snp_col = "RSID",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "FREQ1",
                                 pval_col = "PVAL",
                                 chr_col = "CHR",
                                 pos_col = "POS")


outcome.dat = read_outcome_data(
  snps= exposure.dat$SNP,
  filename = file.path("path_to_outcome_gwassumstats"),
  sep= "\t",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FREQ1",
  pval_col = "PVAL",
  chr_col = "CHR",
  pos_col = "POS")
outcome.dat$outcome = outcome_trait_name


#Step 2: Harmonisation
dat= harmonise_data(exposure_dat = exposure.dat, 
                    outcome_dat = outcome.dat)


# Step 3: Performing MR
mr_results= mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(mr_results)


p1= mr_scatter_plot(mr_results, dat)
res_single= mr_singlesnp(dat)

p2= mr_forest_plot(res_single)

res_1out= mr_leaveoneout(dat)
p3= mr_leaveoneout_plot(res_1out)
p4= mr_funnel_plot(res_single)

# Step 4: Plot all 3 plots together, and export
plot <- grid.arrange(p1[[1]], p2[[1]], p3[[1]], p4[[1]],  nrow = 2, ncol=2)
ggsave(filename = file.path(plot_dir,
                            paste0("2MR_",exposure_trait_name, "_vs_", outcome_trait_name,".png")), 
       plot = plot,
       width = 50,
       height = 30,
       units = "cm",
       dpi=300)


# Step 5: Sensitivity analysis
plei_res= mr_pleiotropy_test(dat) #Pleiotroppy
heteo_res= mr_heterogeneity(dat) #Heterogeneity

# Step 6: Merge results from Twosample MR and Sensitivity check, and export
combine_res=plyr::rbind.fill(mr_results, plei_res, heteo_res)
vroom_write(x = combine_res,
            file = file.path(result_dir, paste0("2MR_",exposure_trait_name, "_vs_", outcome_trait_name, ".txt")),
            delim = "\t")
