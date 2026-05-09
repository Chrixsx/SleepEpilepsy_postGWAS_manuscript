# Step 1 in LDSC is to convert original sumstat to LDSC suitable format (munge)

conda activate ldsc

APP_DIR=" " # full path to ldsc software folder
REF_DIR=" " # full path to LD reference panel
SUMSTATS_DIR=" " # full path to gwas sumstat file
MUNG_DIR=" " # path to storage of munge output
HERI_DIR=" " # path to storage of snp-based heritability output
COR_DIR=" " # path to storage of correlation output

${APP_DIR}/munge_sumstats.py \
--sumstats "{$SUMSTATS_DIR}" \
--snp SNP \
--a1 A1 \
--a2 A2 \
--p P \
--chunksize 500000 \
#--N n_size \ #if use fixed number
--N-col N \
--signed-sumstats BETA,0 \
--merge-alleles "${REF_DIR}/w_hm3.snplist" \
--out "${MUNG_DIR}/munge_filename"
