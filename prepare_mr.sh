#!/bin/bash 
#SBATCH --job-name=pre_mr
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --mem=7g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

BASE=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1
GENO_BASE="${BASE}/genotype_mqtl"
RESULTS_DIR="${BASE}/mQTL_results"
CLEANED_RESULTS_DIR="${BASE}/mQTL_cleaned_results"
COR_MATRIX_DIR="${BASE}/cor_matrix"
PLINK_EXEC="./plink1.9_dir/plink"

CHUNK_ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
CHUNK_FILE="${BASE}/manifest_chunks/chunk_${CHUNK_ID}"

# if [ -d "$COR_MATRIX_DIR" ]; then
#     echo "Directory $COR_MATRIX_DIR exists. Cleaning..."
#     rm -rf "$COR_MATRIX_DIR"
# fi
# mkdir -p "$COR_MATRIX_DIR"

# if [ -d "$CLEANED_RESULTS_DIR" ]; then
#     echo "Directory $CLEANED_RESULTS_DIR exists. Cleaning..."
#     rm -rf "$CLEANED_RESULTS_DIR"
# fi
# mkdir -p "$CLEANED_RESULTS_DIR"

awk '{print $1, $2}' ${CHUNK_FILE} | while read CpG chr; do
    BFILE="${GENO_BASE}/chr${chr}_filtered"
    CLEANED_ASSOC="${CLEANED_RESULTS_DIR}/${CpG}.cleaned.assoc.txt"
    ASSOC_FILE="${RESULTS_DIR}/${CpG}.assoc.txt"

    if [[ -f "$ASSOC_FILE" ]]; then
        echo "Processing CpG: ${CpG} (Chr: ${chr})"

        # Filter variants that are SNPs and AF > 0.01
        awk '
        BEGIN {FS=" "; OFS="\t"} 
        NR == 1 {print; next} # Keep the header
        {
            # Check if Col 5 and Col 6 are only A, C, T, or G
            if ($5 ~ /^[ATCG]$/ && $6 ~ /^[ATCG]$/) {
                # Check if Col 7 (af) is > 0.01
                if ($7 > 0.01) {
                    print $0
                }
            }
        }' "$ASSOC_FILE" > "$CLEANED_ASSOC"

        awk 'NR > 1 && $3 != "" {print $2}' "$CLEANED_ASSOC" > "temp_snplist_${CpG}.txt"
    
        # Generate the LD matrix
        $PLINK_EXEC --bfile "$BFILE" \
                    --extract "temp_snplist_${CpG}.txt" \
                    --r square \
                    --out "${COR_MATRIX_DIR}/${CpG}_matrix" \
                    --write-snplist \
                    --allow-extra-chr \
                    --allow-no-sex \
                    --silent
                    
        rm "temp_snplist_${CpG}.txt"
    else
        echo "Error: ${ASSOC_FILE} not found. Skipping ${CpG}."
    fi
done