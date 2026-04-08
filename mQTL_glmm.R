library(glmmTMB)
library(gaston)
library(coxme)
library(Matrix)

# --- Path Setup ---
BASE <- "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1"
cpg = "cg00000622"
OUTDIR <- file.path(BASE, "mQTL_results_glmm")
CHUNK_FILE <- file.path(BASE, "manifest_final.txt")

# --- 1. Load the GRM ---
# Reading as square matrix
grm_mat <- as.matrix(read.table(file.path(BASE, "grm_mqtl.txt"), header = FALSE))
n_samples <- nrow(grm_mat)
# Ensure it is a formal symmetric matrix for the solver
grm_sparse <- as(grm_mat, "dgCMatrix") 
sample_ids <- as.character(1:nrow(grm_mat))
rownames(grm_mat) <- sample_ids
colnames(grm_mat) <- sample_ids

# --- 2. Load Phenotypes & Covariates ---
meth_df <- read.table(file.path(BASE, "meth_pheno.txt"), header = TRUE, check.names = FALSE, row.names = 1)
unmeth_df <- read.table(file.path(BASE, "unmeth_pheno.txt"), header = TRUE, check.names = FALSE, row.names = 1)
covar <- read.table(file.path(BASE, "covariate_mqtl.txt"), header = FALSE)
colnames(covar) <- paste0("C", 1:ncol(covar))

# --- 3. Process the Manifest ---
manifest <- read.table(CHUNK_FILE, header = FALSE)
colnames(manifest) <- c("CpG", "chr", "pos")

for (i in 1:nrow(manifest)) {
  cpg <- manifest$CpG[i]
  chr <- manifest$chr[i]
  
  M <- meth_df[[cpg]]
  U <- unmeth_df[[cpg]]
  log_total <- log(M + U)
  M[M <= 0] <- 0.001 
  
  geno_path <- file.path(BASE, "genotype_mqtl", paste0("chr", chr, "_filtered"))
  snps_path <- file.path(BASE, "SNPs_mQTL", paste0(cpg, "_chr", chr, ".txt"))
  
  if (file.exists(snps_path) && file.size(snps_path) > 0) {
    bed <- read.bed.matrix(geno_path)
    target_snps <- read.table(snps_path)$V1
    
    for (s_id in target_snps) {
      # Single data frame for the 'data' argument
      snp_idx <- which(bed@snps$id == s_id)
      df_analysis <- data.frame(
        ID = as.factor(1:n_samples), # Grouping factor for random effect
        meth = M,
        log_total = log_total,
        SNP = as.numeric(as.matrix(bed[, snp_idx])),
        log_beta= log(M / (M + U))
      )
      df_analysis <- cbind(df_analysis, covar)
      
      # --- The Correct glmmTMB Call ---
      # We use the 'equalto' logic by defining the structure in the formula
      # and providing the matrix via the sparse covariance link
      covar_names <- colnames(covar)[-1]
      formula1 <- paste0("meth ~ SNP + ", paste(covar_names, collapse = " + "), " + offset(log_total) + (1 | ID)")
      formula2 <- as.formula(paste("log_beta ~ SNP +", paste(covar_names, collapse = " + "), "+ (1 | ID)"))

      tight_control <- glmmTMBControl(
        optimizer = nlminb,
        optArgs = list(iter.max = 2000, eval.max = 2000),
        profile = TRUE # This helps with Hessian stability in Gamma models
      )
      try({
        # Since 'equalto' is experimental, the most robust 'basic' way 
        # to ensure the GRM is used is often via the 'mapped' covariance 
        # but here we use the formula you identified:
        model_glmm <- glmmTMB(as.formula(formula1),
                         data = df_analysis,
                         family = Gamma(link = "log"),
                         # Here is where we "plug in" the GRM matrix structure
                         # to replace the default identity matrix for (1|ID)
                         start = list(theta = log(chol(grm_sparse)@x[1:1])), # Rough initialization
                         map = list(theta = factor(rep(NA, 1)))) # Fixes the variance structure

        model_lmm <- lmekin(as.formula(formula2), data = df_analysis, varlist = grm_mat)
        
        # Stats
        res <- summary(model)$coefficients$cond["SNP", ]
        
        # Write Out
        write.table(data.frame(CpG=cpg, SNP=s_id, Beta=res[1], P=res[4]),
                    file.path(OUTDIR, paste0(cpg, "_res.txt")), 
                    append=T, row.names=F, sep="\t", quote=F)
      }, silent = TRUE)
    }
  }
}