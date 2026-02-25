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


res_susie <- susie_suff_stat(XtX = (n_x - 1)*R, Xty = sqrt(n_x - 1)*Z_x, yty = n_x - 1, n = n_x, estimate_residual_variance = FALSE)
num_cs <- length(res_susie$sets$cs)
print(num_cs)
print(Sys.time())
print(lobstr::mem_used())

saveRDS(num_cs, file = snakemake@output[["SuSiE_results"]])




