#!/bin/bash

#SBATCH --job-name=fit1
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5000M
#SBATCH --nodes=1
#SBATCH --array=1-20

# approx 4 hrs to run
module load plink
module load python/3.7.0
set -o errexit # exit on errors


APP_DIR=" " # full path to mixer software folder
REF_DIR=" " # full path to 1000G plink files
SUMSTATS_DIR=" " # full path to gwas sumstat file
OUT_DIR=" " # output storage
REP=1 # repeat for REP in 1..20

# Univariate MiXeR
python ${APP_DIR}/precimed/mixer.py fit1 \
--trait1-file "${SUMSTATS_DIR}/trait1_gwas_file.gz" \
--out ${OUT_DIR}/trai1.fit.rep${SLURM_ARRAY_TASK_ID} \
--extract ${REF_DIR}/1000G.EUR.QC.prune.rep${SLURM_ARRAY_TASK_ID}.snps \
--bim-file "${REF_DIR}/1000G.EUR.bim" \
--ld-file "${REF_DIR}/1000G.EUR.QC.ld" \
--lib  ${APP_DIR}/src/build/lib/libbgmg.so
