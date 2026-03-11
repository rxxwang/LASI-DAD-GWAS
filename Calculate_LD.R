library(bigsnpr)
library(data.table)
library(magrittr)

CpG = snakemake@params[["CpG"]]
chr = snakemake@params[["chr"]]

# --- 2. Function to Calculate LD Matrix ---
calculate_ld_in_r <- function(cpg_id, chr) {
  
  # Load your cleaned association results to get the SNP list
  assoc_path <- file.path(snakemake@input[["exposure"]])
  df_stats <- fread(assoc_path)
  snp_list <- df_stats$rs  # Ensure this matches your ID column name
  
  bed_orig <- snakemake@input[["genotype_bed"]]
  bim_orig <- snakemake@input[["genotype_bim"]]
  fam_orig <- snakemake@input[["genotype_fam"]]
  
  # 1. Manually check how many individuals are actually in your FAM file
  # We use 'wc -l' via system command to get the true line count
  n_individuals <- as.numeric(system(paste("wc -l <", fam_orig), intern = TRUE))
  message("True number of individuals: ", n_individuals)

  # 2. Create a truly clean dummy FAM file
  # We only need 6 columns. We'll just generate IDs 1 to N.
  clean_fam <- data.frame(
    FID = 1:n_individuals,
    IID = 1:n_individuals,
    PID = 0, MID = 0, SEX = 0, PHENO = -9
  )

  # 3. Setup the temp directory again
  tmp_dir <- tempdir()
  tmp_base <- file.path(tmp_dir, "repaired_geno")
  fwrite(clean_fam, paste0(tmp_base, ".fam"), sep = " ", col.names = FALSE)

  # 4. Link the original BED and BIM
  file.symlink(normalizePath(bed_orig), paste0(tmp_base, ".bed"))
  file.symlink(normalizePath(bim_orig), paste0(tmp_base, ".bim"))

  # 5. NOW read the Bed. 
  # It will see 'n_individuals' lines in the FAM and create the correct N x P matrix
  rds_path <- snp_readBed(paste0(tmp_base, ".bed"))
  obj.bigSNP <- snp_attach(rds_path)

  # Check the dimensions again. It should be [n_individuals, 2082079]
  print(dim(obj.bigSNP$genotypes))
  
  # 6. Extract only the genotypes for the SNPs we need
  # Match SNPs in the genotype file with your summary stats
  idx_in_geno <- match(df_stats$rs, obj.bigSNP$map$marker.ID)

  # 7. Filter for SNPs that actually exist in the genotypes
  valid_mask <- !is.na(idx_in_geno)
  valid_indices <- idx_in_geno[valid_mask]
  
  # Get the genotype matrix (n_samples x n_snps)
  G <- obj.bigSNP$genotypes[, valid_indices]
  
  # 8. Calculate the Correlation (LD) Matrix
  # 'use = "pairwise.complete.obs"' is critical if there is missing genotype data
  R <- cor(G, use = "pairwise.complete.obs")
  
  return(R)
}

R <- calculate_ld_in_r(CpG, chr)
saveRDS(R, file = snakemake@output[["results"]])