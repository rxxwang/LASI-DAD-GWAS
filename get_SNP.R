# R Script: process_manifest.R
library(data.table)

# Load data
manifest <- fread("manifest.txt", header = FALSE, col.names = c("CpG", "chr", "pos"))
pheno <- fread("pheno_mqtl.txt")

# 1. Handle chr0: Filter out unknowns (unless you have a specific lookup table)
manifest <- manifest[chr != 0 & !is.na(chr) & !is.na(pos)]

# 2. Generate SNP lists by chromosome
# We will loop through each chromosome to read the .bim file once (efficiency)
for (c in unique(manifest$chr)) {
  bim_file <- paste0("genotype_mqtl/chr", c, "_filtered.bim")
  if(!file.exists(bim_file)) next
  
  bim <- fread(bim_file, header = FALSE, col.names = c("chr", "snp", "morgans", "pos", "a1", "a2"))
  chr_manifest <- manifest[chr == c]
  
  for (i in 1:nrow(chr_manifest)) {
    cpg_name <- chr_manifest$CpG[i]
    cpg_pos <- chr_manifest$pos[i]
    
    # Define 1Mbp window (1,000,000 bp)
    window_snps <- bim[pos >= (cpg_pos - 1000000) & pos <= (cpg_pos + 1000000), snp]
    
    # Save SNP list
    fwrite(list(window_snps), 
           file = paste0("SNPs_mQTL/", cpg_name, "_chr", c, ".txt"), 
           col.names = FALSE)
  }
}

fwrite(manifest, "manifest_cleaned.txt", sep = "\t", col.names = FALSE)

# R Script: update_fam_and_manifest.R
# Load manifest processed from Step 1
manifest <- fread("manifest_cleaned.txt") 
pheno <- fread("pheno_mqtl.txt") # Columns: FID, IID, CpG1, CpG2...

# Reorder pheno to match a template .fam if necessary, but here we assume FID/IID match
# We will update the manifest with the 'column index' for GEMMA (-n flag)
cpg_names <- colnames(pheno)[3:ncol(pheno)]
manifest$pheno_idx <- match(manifest$CpG, cpg_names)

# Update .fam files for each chromosome
for (c in 1:22) {
  fam_path <- paste0("genotype_mqtl/chr", c, "_filtered.fam")
  if(!file.exists(fam_path)) next
  
  fam <- fread(fam_path, header = FALSE)
  # Standard .fam has 6 cols. We merge pheno data starting from column 7+
  # GEMMA -n 1 will look at the 6th column, -n 2 at the 7th, etc.
  
  # Ensure IDs match
  updated_fam <- merge(fam[, 1:5], pheno, by.x = c("V1", "V2"), by.y = c("FID", "IID"), all.x = TRUE)
  
  # Save back to .fam (Note: this overwrites the 6th column and appends others)
  fwrite(updated_fam, file = fam_path, sep = " ", col.names = FALSE)
}

fwrite(manifest, "manifest_final.txt", sep = "\t", col.names = FALSE)