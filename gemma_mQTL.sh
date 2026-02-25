#!/bin/bash 
#SBATCH --job-name=gemma
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

# Turn vcf file to bed file
# for chr in {1..22}; do
#     plink2 --bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr} \
#            --keep /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/keep.ids.txt \
#            --make-bed \
#            --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype_mqtl/chr${chr}_filtered
# done

# BASE="/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/"
# if [ ! -d "${BASE}/mQTL_results" ]; then
#   mkdir "${BASE}/mQTL_results"
# fi
# if [ ! -d "${BASE}/SNPs_mQTL" ]; then
#   mkdir "${BASE}/SNPs_mQTL"
# fi

# Rscript --vanilla /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/get_SNP.R

# echo "Finished revision"

BASE=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1
# mkdir -p ${BASE}/manifest_chunks
# # Split the manifest into 100 files, prefixed with 'chunk_'
# split -d -n l/100 ${BASE}/manifest_final.txt ${BASE}/manifest_chunks/chunk_

CHUNK_ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
CHUNK_FILE="${BASE}/manifest_chunks/chunk_${CHUNK_ID}"

GRM="${BASE}/grm_mqtl.txt"
COVAR="${BASE}/covariate_mqtl.txt"
OUTDIR="${BASE}/mQTL_results"

mkdir -p logs
while read -r CpG chr pos pheno_idx; do

    GENO="${BASE}/genotype_mqtl/chr${chr}_filtered"
    SNPS="${BASE}/SNPs_mQTL/${CpG}_chr${chr}.txt"

    # Only run if SNP list exists and isn't empty
    if [ -s "$SNPS" ]; then
        gemma -bfile ${GENO} \
            -k ${GRM} \
            -c ${COVAR} \
            -snps ${SNPS} \
            -n ${pheno_idx} \
            -lmm 1 \
            -o ${CpG} \
            -outdir ${OUTDIR}
    else
        echo "No SNPs for ${CpG} on chr${chr}. Skipping."
    fi

done < "$CHUNK_FILE"