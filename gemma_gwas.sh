#!/bin/bash 
#SBATCH --job-name=gemma
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=7g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

# Turn vcf file to bed file
for chr in {1..22}; do
    plink2 --vcf /net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/QC_UPenn/vcf_new_ids/chr${chr}.vcf.gz \
           --make-bed \
           --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr${chr}
    echo "Finished converting chr${chr}"
done

# Combine all the bed file together
plink2 --bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/genotype/chr1 \
       --pmerge-list /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merge_list.txt \
       --make-bed \
       --out /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno

# Get the samples that have biomarker data
awk 'NR>1 {print $1, $2}' pheno_gwas.txt | sort > pheno.ids

# Filter the samples
plink2 --bfile merged_geno \
       --keep phenos.ids \
       --make-bed \
       --out merged_geno_filtered

# Check the sampleid
head merged_geno_filtered.fam
head pheno_gwas.txt

# Add the phenotype into .fam file
awk 'NR==FNR{
        if(NR>1) {  # skip header
            pheno[$2] = $0   # store the whole line keyed by IID
        }
        next
     }
     {
        split(pheno[$2], pvals)
        $6 = pvals[3]   # gfap
        $7 = pvals[4]   # nfl
        $8 = pvals[5]   # ptau
        $9 = pvals[6]   # totaltau
        $10 = pvals[7]  # abeta_ratio
        print
     }' pheno_gwas.txt merged_geno_filtered.fam > merged_geno_filtered_new.fam
mv merged_geno_filtered_new.fam merged_geno_filtered.fam

# GWAS via GEMMA
gemma -bfile /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/merged_geno_filtered \
      -k /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/grm.txt \
      -c /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/covariate.txt \
      -n 1 \
      -lmm 4 \
      -o /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/gwas_results
