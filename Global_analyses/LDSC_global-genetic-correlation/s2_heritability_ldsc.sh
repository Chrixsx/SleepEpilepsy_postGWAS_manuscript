# Step 2 in LDSC is to estimate SNP-based heritability from GWAS file

conda activate ldsc

APP_DIR=" " # full path to ldsc software folder
REF_DIR=" " # full path to LD reference panel
SUMSTATS_DIR=" " # full path to gwas sumstat file
MUNG_DIR=" " # path to storage of munge output
HERI_DIR=" " # path to storage of snp-based heritability output
COR_DIR=" " # path to storage of correlation output

${APP_DIR}/ldsc.py \
--h2 "${MUNG_DIR}/munge_filename" \
--ref-ld-chr ${REF_DIR}/ \
--w-ld-chr ${REF_DIR}/ \
--out "${HERI_DIR}/snph2_filename"
