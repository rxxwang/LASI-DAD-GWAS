#!/bin/bash 
#SBATCH --job-name=gwassumm
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

BASE=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1
GWAS=${BASE}/gwas_results
GENO=${BASE}/genotype

declare -a PHENOS=("gfap" "nfl" "ptau" "totaltau" "abeta_ratio")

for pheno in "${PHENOS[@]}"; do
    echo "Processing phenotype: ${pheno}"

    final_pheno_file=${GWAS}/${pheno}.assoc.txt
    rm -f "${final_pheno_file}"

    first_chr=true

    for chr in {1..22}; do
        BIM=${GENO}/chr${chr}_filtered.bim
        [[ ! -f ${BIM} ]] && continue

        nsnp=$(wc -l < "${BIM}")
        nchunk=$(( (nsnp + 999999) / 1000000 ))

        chr_out=${GWAS}/${pheno}_chr${chr}.assoc.txt
        rm -f "${chr_out}"

        first_chunk=true

        for ((chunk=1; chunk<=nchunk; chunk++)); do
            chunk_file=${GWAS}/${pheno}_chr${chr}_${chunk}.assoc.txt
            [[ ! -f ${chunk_file} ]] && continue

            if ${first_chunk}; then
                cat "${chunk_file}" > "${chr_out}"
                first_chunk=false
            else
                tail -n +2 "${chunk_file}" >> "${chr_out}"
            fi
        done

        echo "${pheno}: Chr${chr} done"

        # Now append chr file to phenotype-wide file
        if [[ -f ${chr_out} ]]; then
            if ${first_chr}; then
                cat "${chr_out}" > "${final_pheno_file}"
                first_chr=false
            else
                tail -n +2 "${chr_out}" >> "${final_pheno_file}"
            fi
        fi
    done
done
