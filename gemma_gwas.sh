#!/bin/bash 
#SBATCH --job-name=gemma
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

# # Turn vcf file to bed file
# for chr in {1..22}; do
# #     plink2 --vcf /net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/QC_UPenn/vcf_new_ids/chr${chr}.vcf.gz \
# #            --make-bed \
# #            --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}
# #     echo "Finished converting chr${chr}"
#     plink2 --bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr} \
#            --keep /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno.ids \
#            --make-bed \
#            --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}_filtered
#     awk 'NR==FNR{
#         if(NR>1) {  # skip header
#             pheno[$2] = $0   # store the whole line keyed by IID
#         }
#         next
#      }
#      {
#         split(pheno[$2], pvals)
#         $6 = pvals[3]   # gfap
#         $7 = pvals[4]   # nfl
#         $8 = pvals[5]   # ptau
#         $9 = pvals[6]   # totaltau
#         $10 = pvals[7]  # abeta_ratio
#         print
#      }' /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno_gwas.txt /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}_filtered.fam > /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}_filtered_new.fam
#     mv /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}_filtered_new.fam /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}_filtered.fam
#     echo "Finished revising chr${chr}"
# done

BASE=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1
GENO=${BASE}/genotype
SNPDIR=${BASE}/SNPs

# Create SNP directory if it does not exist
[[ ! -d ${SNPDIR} ]] && mkdir ${SNPDIR}

for chr in {1..22}; do
    BIM=${GENO}/chr${chr}_filtered.bim

    if [[ ! -f ${BIM} ]]; then
        echo "Missing ${BIM}, skipping chr${chr}"
        continue
    fi

    nsnp=$(wc -l < ${BIM})
    nchunk=$(( (nsnp + 999999) / 1000000 ))

    echo "Chr${chr}: ${nsnp} SNPs → ${nchunk} chunk(s)"

    for ((chunk=1; chunk<=nchunk; chunk++)); do
        start=$(( (chunk - 1) * 1000000 + 1 ))
        end=$(( chunk * 1000000 ))
        (( end > nsnp )) && end=${nsnp}

        echo "  chr${chr} chunk ${chunk}: SNP ${start}-${end}"

        sed -n "${start},${end}p" ${BIM} | awk '{print $2}' \
            > ${SNPDIR}/chr${chr}_${chunk}.txt

        sbatch gemma_job.sh ${chr} ${chunk}
    done
done

# # Combine all the bed file together
# plink2 --bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr1 \
#        --pmerge-list /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merge_list.txt \
#        --make-bed \
#        --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno

# # Get the samples that have biomarker data
# awk 'NR>1 {print $1, $2}' /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno_gwas.txt | sort > /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno.ids

# # Filter the samples
# plink2 --bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno \
#        --keep /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/phenos.ids \
#        --make-bed \
#        --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered

# # Check the sampleid
# head /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered.fam
# head /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno_gwas.txt

# # Add the phenotype into .fam file
# awk 'NR==FNR{
#         if(NR>1) {  # skip header
#             pheno[$2] = $0   # store the whole line keyed by IID
#         }
#         next
#      }
#      {
#         split(pheno[$2], pvals)
#         $6 = pvals[3]   # gfap
#         $7 = pvals[4]   # nfl
#         $8 = pvals[5]   # ptau
#         $9 = pvals[6]   # totaltau
#         $10 = pvals[7]  # abeta_ratio
#         print
#      }' /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno_gwas.txt /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered.fam > /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered_new.fam
# mv /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered_new.fam /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered.fam

# # GWAS via GEMMA
# gemma -bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr1_filtered \
#       -k /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/grm.txt \
#       -c /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/covariate.txt \
#       -snps /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/first100k_snps_chr1.txt \
#       -n 1 \
#       -lmm 1 \
#       -o gfap \
#       -outdir /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/gwas_results

# # GWAS via GEMMA
# gemma -bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered \
#       -k /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/grm.txt \
#       -c /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/covariate.txt \
#       -snps /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/first100k_snps.txt \
#       -n 1 \
#       -lmm 4 \
#       -o /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/gwas_results
