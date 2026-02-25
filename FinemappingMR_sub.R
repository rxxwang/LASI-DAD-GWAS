library(finemappingMR)
library(susieR)
library(dplyr)
setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/")

CpG = snakemake@params[["CpG"]]
bio = snakemake@params[["biomarker"]]
chr = snakemake@params[["chr"]]
exposure = read.table(snakemake@input[["exposure"]], header = T)
cor_matrix = read.table(snakemake@input[["cor_matrix"]])
outcome_file = read.table(snakemake@input[["outcome"]], header = T) 
outcome = outcome_file %>% filter(rs %in% exposure$rs) %>% arrange(ps)
overlap_rs = outcome$rs
cor_matrix = cor_matrix[exposure$rs %in% overlap_rs, exposure$rs %in% overlap_rs]
exposure = exposure %>% filter(rs %in% overlap_rs) %>% arrange(ps)
if(sum(is.na(as.matrix(cor_matrix))) > 0){
    na_cols <- which(colSums(is.na(cor_matrix)) > 0)
    removed_count <- length(na_cols)
    cor_matrix_clean <- cor_matrix[-na_cols, -na_cols]
    message(paste("Removed", removed_count, "rows and columns due to NA values."))
    R = as.matrix(Matrix::nearPD(as.matrix(cor_matrix_clean))$mat)
    outcome = outcome[-na_cols,]
    exposure = exposure[-na_cols,]
} else{
    R = as.matrix(Matrix::nearPD(as.matrix(cor_matrix))$mat)
}

covariate_gwas = read.table(snakemake@input[["covariate_gwas"]], header = F, sep = " ")
covariate_mqtl = read.table(snakemake@input[["covariate_mqtl"]], header = F, sep = " ")
pheno_gwas = read.table(snakemake@input[["pheno_gwas"]], header = T, sep = " ")
pheno_mqtl = read.table(snakemake@input[["pheno_mqtl"]], header = T, sep = " ")
covariate_gwas = covariate_gwas[pheno_gwas$IID %in% pheno_mqtl$IID,]
pheno_gwas = pheno_gwas[pheno_gwas$IID %in% pheno_mqtl$IID,]
colnames(pheno_mqtl)[3:ncol(pheno_mqtl)] = sub("_.*", "", colnames(pheno_mqtl)[3:ncol(pheno_mqtl)])

stopifnot(all.equal(exposure$rs, outcome$rs))
dim(cor_matrix)
dim(outcome)
dim(exposure)

Zx = exposure$beta / exposure$se
Zy = outcome$beta / outcome$se
n_x = 1129
n_y = 2224
Z_x = Zx * sqrt(n_x-1) / sqrt(Zx^2 + n_x - 2)
Z_y = Zy * sqrt(n_y-1) / sqrt(Zy^2 + n_y - 2)

print(Sys.time())
res_susie <- susie_suff_stat(XtX = (n_x - 1)*R, Xty = sqrt(n_x - 1)*Z_x, yty = n_x - 1, n = n_x, estimate_residual_variance = FALSE, max_iter = 1000)
if (res_susie$niter >= 1000){
    num_cs = -1
} else {
    num_cs <- length(res_susie$sets$cs)
}
print(Sys.time())


if(num_cs > 0){
    sig_SNP = c()
    for(k in 1:num_cs){
        sig_SNP = c(sig_SNP, res_susie$sets$cs[[k]])
    }
    ps_SNP = outcome[sig_SNP,3]
    if(num_cs > 1){
        sorted_SNP <- sort(ps_SNP)
        gaps <- diff(sorted_SNP)
        large_gaps_exist <- any(gaps > 1000000)
    } else {
        large_gaps_exist = FALSE
    }
    if(!large_gaps_exist){
        min_ps = max(min(ps_SNP)-500000, min(outcome$ps))
        max_ps = min(max(ps_SNP)+500000, max(outcome$ps))
        new_SNP = outcome$ps <= max_ps & outcome$ps >= min_ps
    } else{
        gap_num = which(gaps > 1000000)
        ps1 = max(min(ps_SNP)-500000, min(outcome$ps))
        ps2 = ps_SNP[gap_num]+500000
        ps3 = ps_SNP[gap_num+1]-500000
        ps4 = min(max(ps_SNP)+500000, max(outcome$ps))
        new_SNP = (outcome$ps <= ps4 & outcome$ps >= ps3) | (outcome$ps <= ps2 & outcome$ps >= ps1)
    }

    exposure = exposure[new_SNP,]
    outcome = outcome[new_SNP,]
    R = R[new_SNP, new_SNP]
    Zx = exposure$beta / exposure$se
    Zy = outcome$beta / outcome$se
    Z_x = Zx * sqrt(n_x-1) / sqrt(Zx^2 + n_x - 2)
    Z_y = Zy * sqrt(n_y-1) / sqrt(Zy^2 + n_y - 2)

    y_bio = pheno_gwas[,bio]
    x_bio = as.matrix(covariate_gwas[,2:ncol(covariate_gwas)])
    nabio = !is.na(y_bio)
    y_cpg = pheno_mqtl[,CpG]
    x_cpg = as.matrix(covariate_mqtl[,2:ncol(covariate_mqtl)])
    nacpg = !is.na(y_cpg)
    na = nabio*nacpg == 1
    y_bio = y_bio[na]
    x_bio = x_bio[na,]
    y_cpg = y_cpg[na]
    x_cpg = x_cpg[na,]
    model1 = lm(y_bio ~ x_bio)
    model2 = lm(y_cpg ~ x_cpg)
    rho_obs = cor(model1$residuals, model2$residuals)
    rho = rho_obs * sqrt(sum(na)^2/sum(nabio)/sum(nacpg))

    n_x = sum(nacpg)
    n_y = sum(nabio)

    res_func_attempt <- tryCatch({
        res_func <- finemappingMR_sampleOverlap(Z_x=Z_x, Z_y=Z_y, R=R, rho=0, n_x=n_x, n_y=n_y)
        if (res_func$res$iter >= 1000) {
            message(paste("CpG", CpG, bio,": No converging!"))
            "No converging" # Return NULL so we don't save bad data
        } else{
            res_func$res
        }
    }, error = function(e) {
        message(paste("!! Error in CpG", CpG, bio, ":", e$message))
        NULL # Return NULL on error to keep the loop moving
    })
    if (!is.null(res_func_attempt)) {
        cpg_bio = paste0(CpG, "_", bio)
        message(cpg_bio, ": Succeed!")
        results <- res_func_attempt
    } else {
        results = "NULL"
    }
} else {
    print(paste0(CpG, ": No variables in credible sets of SuSiE."))
    results = "No variables"
}
print(Sys.time())
print(lobstr::mem_used())

saveRDS(list(result = results, SuSiE = num_cs), file = snakemake@output[["results"]])




