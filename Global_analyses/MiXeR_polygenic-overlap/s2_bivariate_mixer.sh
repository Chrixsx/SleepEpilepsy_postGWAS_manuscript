#!/bin/bash

#SBATCH --job-name=test2.4
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5000M
#SBATCH --nodes=1
#SBATCH --array=1-20

# approx 3.5 hrs
module load python/3.7.0
set -o errexit # exit on errors

APP_DIR=" " # full path to mixer software folder
REF_DIR=" " # full path to 1000G plink files
SUMSTATS_DIR=" " # full path to gwas sumstat file
OUT_DIR=" " # output storage

# Bivariate MiXeR
python3 ${APP_DIR}/precimed/mixer.py test2 \
--trait1-file "${SUMSTATS_DIR}/trait1_gwas_file.gz" \
--trait2-file "${SUMSTATS_DIR}/trait2_gwas_file.gz" \
--load-params-file ${OUT_DIR}/combined/trait1_vs_trait2.fit.rep${SLURM_ARRAY_TASK_ID}.json \
--out ${OUT_DIR}/combined/trait1_vs_trait2.test.rep${SLURM_ARRAY_TASK_ID} \
--bim-file "${REF_DIR}/1000G.EUR.bim" \
--ld-file "${REF_DIR}/1000G.EUR.QC.ld" \
--lib  ${APP_DIR}/src/build/lib/libbgmg.so
