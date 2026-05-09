# Step 3 in LDSC is to estimate pairwise correlation between any two traits
# Note: LDSC estimated correlation coefficient (rg) is NOT bounded within [-1:1]. 
# Note: If estimated rg is out of bound, please check output log file, and refer to Issue in https://github.com/bulik/LDSC

conda activate ldsc

APP_DIR=" " # full path to ldsc software folder
REF_DIR=" " # full path to LD reference panel
SUMSTATS_DIR=" " # full path to gwas sumstat file
MUNG_DIR=" " # path to storage of munge output
HERI_DIR=" " # path to storage of snp-based heritability output
COR_DIR=" " # path to storage of correlation output

${APP_DIR}/ldsc.py \
--rg "${MUNG_DIR}/trait1_munge_filename",\
"${MUNG_DIR}/trait2_munge_filename" \
--ref-ld-chr ${REF_DIR}/ \
--w-ld-chr ${REF_DIR}/ \
--out "${COR_DIR}/correlation_trait1_trait2"
