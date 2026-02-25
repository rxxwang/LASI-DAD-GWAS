#!/bin/bash 
#SBATCH --job-name=ld_clump
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
CLEANED_RESULTS_DIR="${BASE}/mQTL_clumped_results"
COR_MATRIX_DIR="${BASE}/cor_matrix"
OUTPUT_DIR="${BASE}/clumped_results"
PLINK_EXEC="./plink1.9_dir/plink"

CHUNK_ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
CHUNK_FILE="${BASE}/manifest_chunks/chunk_${CHUNK_ID}"

if [ -d "$OUTPUT_DIR" ]; then
    echo "Directory $OUTPUT_DIR exists. Cleaning..."
    rm -rf "$OUTPUT_DIR"
fi
mkdir -p "$OUTPUT_DIR"

if [ -d "$CLEANED_RESULTS_DIR" ]; then
    echo "Directory $CLEANED_RESULTS_DIR exists. Cleaning..."
    rm -rf "$CLEANED_RESULTS_DIR"
fi
mkdir -p "$CLEANED_RESULTS_DIR"

# --- Main Loop ---
# Col 1 = CpG, Col 2 = Chromosome (No header in manifest)
awk '{print $1, $2}' ${CHUNK_FILE} | while read CpG chr; do

    BFILE="${GENO_BASE}/chr${chr}_filtered"
    ASSOC_FILE="${RESULTS_DIR}/${CpG}.assoc.txt"

    if [[ -f "$ASSOC_FILE" ]]; then
        echo "Processing CpG: ${CpG} (Chr: ${chr})"

        # 1. Identify duplicates for this specific run's genotype file
        # We save this in the output dir to keep the working dir clean
        awk '{print $2}' "${BFILE}.bim" | sort | uniq -d > "temp_dups_${CpG}.txt"

        # 2. Create the clumping input
        echo -e "SNP\tP" > "temp_clump_${CpG}.txt"
        awk 'NR > 1 {print $2, $12}' "$ASSOC_FILE" >> "temp_clump_${CpG}.txt"

        # 3. Run PLINK Clumping
        # --exclude removes the duplicates found in step 1
        $PLINK_EXEC --bfile "$BFILE" \
              --exclude "temp_dups_${CpG}.txt" \
              --clump "temp_clump_${CpG}.txt" \
              --clump-field P \
              --clump-snp-field SNP \
              --clump-p1 1e-3 \
              --clump-r2 0.01 \
              --clump-kb 500 \
              --out "${OUTPUT_DIR}/${CpG}_clumped" \
              --silent

        # 4. Cleanup temp files for this CpG
        rm "temp_clump_${CpG}.txt"
        rm "temp_dups_${CpG}.txt"

        CLUMP_FILE="${OUTPUT_DIR}/${CpG}_clumped.clumped"

        CLEAN_OUT="${CLEANED_RESULTS_DIR}/${CpG}_leads_clean.txt"

        # Get the LD clumped mQTL results for each CpG
        if [[ -f "$CLUMP_FILE" ]]; then
            echo "Cleaning association file for ${CpG}..."
            awk 'NR > 1 {print $3}' "$CLUMP_FILE" > "temp_leads_${CpG}.txt"
            awk 'FNR==NR {leads[$1]; next} (FNR==1 || $2 in leads)' \
                "temp_leads_${CpG}.txt" "$ASSOC_FILE" > "$CLEAN_OUT"
            rm "temp_leads_${CpG}.txt"
        else
            echo "No clumped leads found for ${CpG}."
        fi
    else
        echo "Error: ${ASSOC_FILE} not found. Skipping ${CpG}."
    fi
done

echo "Workflow complete. Clumped results are in $OUTPUT_DIR."