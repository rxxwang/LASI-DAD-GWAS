#!/bin/bash 
#SBATCH --job-name=gemmasub
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

chr=$1
chunk=$2

BASE=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1

GENO=${BASE}/genotype/chr${chr}_filtered
SNPS=${BASE}/SNPs/chr${chr}_${chunk}.txt
GRM=${BASE}/grm.txt
COVAR=${BASE}/covariate.txt
OUTDIR=${BASE}/gwas_results

[[ ! -d ${OUTDIR} ]] && mkdir -p ${OUTDIR}

# phenotype index and name
declare -a PHENOS=("gfap" "nfl" "ptau" "totaltau" "abeta_ratio")

for i in {1..5}; do
    pheno=${PHENOS[$((i-1))]}

    gemma -bfile ${GENO} \
          -k ${GRM} \
          -c ${COVAR} \
          -snps ${SNPS} \
          -n ${i} \
          -lmm 1 \
          -o ${pheno}_chr${chr}_${chunk} \
          -outdir ${OUTDIR}
done