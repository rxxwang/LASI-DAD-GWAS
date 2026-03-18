library(finemappingMR)
library(susieR)
library(dplyr)
setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/")

args <- commandArgs(trailingOnly = TRUE)
chunk <- args[1]
manifest = read.table(paste0("manifest_chunks/chunk_", chunk))
pval = readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/pval.model2.RDS")
pmat_filtered <- pval[rowSums(pval < 1e-5, na.rm = TRUE) > 0, ]
sig_cpg = sub("_.*", "", rownames(pmat_filtered))
results = list()
sig_cpg_bio = as.data.frame(matrix(NA, nrow = sum(pmat_filtered < 1e-5), ncol = 3))
biomarker_mapping = data.frame(
    old = c("abeta_ratio", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final"),
    new = c("abeta_ratio", "gfap", "nfl", "ptau", "totaltau")
)
k = 1
manifest = readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/manifest4.RDS") %>% mutate(chr = substr(Chr, start = 4, stop = nchar(Chr)))
for(i in 1:length(sig_cpg)){
    CpG = sig_cpg[i]
    chr = unique(as.numeric(manifest[manifest$CPG == CpG, 6]))
    if(chr == 0){
        next
    }
    biomarker = colnames(pmat_filtered)[which(pmat_filtered[which(sig_cpg == CpG),] <1e-5)]
    for(j in biomarker){
        bio = biomarker_mapping$new[biomarker_mapping$old == j]
        sig_cpg_bio[k,] = c(CpG, bio, chr)
        k = k+1
    }
}
sig_cpg_bio = sig_cpg_bio[1:(k-1),]
colnames(sig_cpg_bio) = c("CpG", "biomarker", "chr")
write.table(sig_cpg_bio, file = "sig_cpgs.txt", row.names = F, col.names = T, sep = "\t", quote = FALSE)

covariate_gwas = read.table("covariate.txt", header = F, sep = " ")
covariate_mqtl = read.table("covariate_mqtl.txt", header = F, sep = " ")
pheno_gwas = read.table("pheno_gwas.txt", header = T, sep = " ")
pheno_mqtl = read.table("pheno_mqtl.txt", header = T, sep = " ")
covariate_gwas = covariate_gwas[pheno_gwas$IID %in% pheno_mqtl$IID,]
pheno_gwas = pheno_gwas[pheno_gwas$IID %in% pheno_mqtl$IID,]
colnames(pheno_mqtl)[3:ncol(pheno_mqtl)] = sub("_.*", "", colnames(pheno_mqtl)[3:ncol(pheno_mqtl)])

for(i in 1:nrow(manifest)){
    # CpG = "cg00001687"

    CpG = manifest$V1[i]
    chr = manifest$V2[i]

    biomarker = colnames(pmat_filtered)[which(pmat_filtered[which(sig_cpg == CpG),] <1e-5)]
    for(j in biomarker){
        bio = biomarker_mapping$new[biomarker_mapping$old == j]
        exposure = read.table(paste0("mQTL_cleaned_results/",CpG,".cleaned.assoc.txt"), header = T)
        cor_matrix = read.table(paste0("cor_matrix/",CpG,"_matrix.ld"))
        chr = exposure$chr[1]
        outcome_file = read.table(paste0("gwas_results/", bio,"_chr", chr, ".assoc.txt"), header = T) 
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

        
        res_susie <- susie_suff_stat(XtX = (n_x - 1)*R, Xty = sqrt(n_x - 1)*Z_x, yty = n_x - 1, n = n_x, estimate_residual_variance = FALSE, max_iter = 1000)
        if (res_susie$niter >= 1000){
            num_cs = -1
        } else {
            num_cs <- length(res_susie$sets$cs)
        }

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
                large_gaps_exist = F
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
                res_func <- finemappingMR_sampleOverlap(Z_x=Z_x, Z_y=Z_y, R=R, rho=rho, n_x=n_x, n_y=n_y, verbose = TRUE)
                if (res_func$res$iter >= 1000) {
                    message(paste("CpG", CpG, bio,": No converging!"))
                    NULL # Return NULL so we don't save bad data
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
                results[[cpg_bio]] <- res_func_attempt
            }
        } else {
            print(paste0(CpG, ": No variables in credible sets of SuSiE."))
        }
        print(Sys.time())
        print(lobstr::mem_used())
    }
}

saveRDS(results, file = paste0("MR_result/chunk", chunk, ".RDS"))




