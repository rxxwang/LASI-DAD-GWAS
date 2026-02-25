#!/bin/bash 
#SBATCH --job-name=finemap
#SBATCH --partition=main
#SBATCH --time=120:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --mem=25g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

CHUNK_ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
Rscript --vanilla /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/FinemappingMR.R $CHUNK_ID