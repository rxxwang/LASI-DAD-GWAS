#!/bin/bash 
#SBATCH --job-name=gemmasum
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=20g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

BASE=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1
MANIFEST="${BASE}/manifest_final.txt"
RESULTS_DIR="${BASE}/mQTL_results"
SUMMARY_FILE="${BASE}/mqtl_summary_report.txt"

# 1. Write the header
echo -e "CpG\tChr\tPos\tTotal_SNPs\tCount_P05\tCount_P1e3\tTop_SNP\tTop_P\tTop_AF" > $SUMMARY_FILE

# 2. Loop through the manifest
while read -r CpG chr pos pheno_idx; do
    FILE="${RESULTS_DIR}/${CpG}.assoc.txt"

    if [ -f "$FILE" ]; then
        # AWK Logic:
        # $2 = SNP name, $7 = Allele Frequency, $12 = p_wald
        awk -v cpg="$CpG" -v ch="$chr" -v ps="$pos" '
            BEGIN { 
                min_p = 1.1; 
                top_snp = "NA"; 
                top_af = "NA";
                p05 = 0; 
                p1e3 = 0; 
                count = 0; 
                OFS="\t" 
            }
            NR > 1 { # Skip header
                count++;
                p = $12; 
                snp = $2; 
                af = $7;
                
                if (p < 0.05) p05++;
                if (p < 0.001) p1e3++;
                
                if (p < min_p) {
                    min_p = p;
                    top_snp = snp;
                    top_af = af;
                }
            }
            END {
                if (count > 0)
                    print cpg, ch, ps, count, p05, p1e3, top_snp, min_p, top_af;
                else
                    print cpg, ch, ps, 0, 0, 0, "NA", "NA", "NA";
            }
        ' "$FILE" >> $SUMMARY_FILE
    else
        echo -e "${CpG}\t${chr}\t${pos}\t0\t0\t0\tNA\tNA\tNA" >> $SUMMARY_FILE
    fi
done < "$MANIFEST"

echo "Summary report generated with AF: $SUMMARY_FILE"