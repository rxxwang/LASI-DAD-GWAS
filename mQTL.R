grm = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/Relatedness_check/Relatedness_check_zheng/grm_2680_genesis.csv")
pcs = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/PCA/PCA_analysis/LASI_only/LASI_PCs_2680_kinThresh0.044.csv")
#unrel = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/PCA/PCA_analysis/LASI_only/unrel_smps_cutoff0.044_N2569_pcair.txt")
mapping = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/ID_mapping/WGS_methy_map.csv")
library(haven)
library(dplyr)
library(stringr)
# Find the significant cpgs
pval = readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/pval.model2.RDS")
pmat_filtered <- pval[rowSums(pval < 1e-5, na.rm = TRUE) > 0, ]
sig_cpg = rownames(pmat_filtered)
# Get beta value of cpgs
cpgs = readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/cpg_list4.RDS")
meth_cpg = paste0(sig_cpg, "_M")
unmeth_cpg = paste0(sig_cpg, "_U")
pca_result = readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/pca_result4.RDS")
manifest = readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/manifest4.RDS")

for(i in 1:500){
  data = readRDS(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/data4/meth_pheno_data.",i,".RDS"))
  if(i == 1){
    meth = data[,colnames(data) %in% meth_cpg]
    unmeth = data[,colnames(data) %in% unmeth_cpg]
    col_cpg = colnames(data)[colnames(data) %in% unmeth_cpg]
  } else{
    meth = cbind(meth, data[,colnames(data) %in% meth_cpg])
    unmeth = cbind(unmeth, data[,colnames(data) %in% unmeth_cpg])
    col_cpg = c(col_cpg, colnames(data)[colnames(data) %in% unmeth_cpg])
  }

  print(paste(i, ":", sum(colnames(data) %in% meth_cpg)))
}
col_cpg = substr(col_cpg, 1, nchar(col_cpg) - 2)
col_cpg = sub("_.*", "", col_cpg)
colnames(meth) = col_cpg
rownames(meth) = data$Sample_name
meth = as.matrix(meth)
colnames(unmeth) = col_cpg
rownames(unmeth) = data$Sample_name
unmeth = as.matrix(unmeth)
beta = log(meth/(meth+unmeth))
beta = as.data.frame(beta)
beta$barcode_w1 = rownames(meth)
beta$sex = data$ragender
beta$age = data$r1hagey
beta$CD8T = data$CD8T
beta$CD4T = data$CD4T
beta$NK = data$NK
beta$Bcell = data$Bcell
beta$Mono = data$Mono

beta <- beta %>% 
  mutate(PC1 = scale(pca_result$x[,1]), PC2 = scale(pca_result$x[,2]), PC3 = scale(pca_result$x[,3]), PC4 = scale(pca_result$x[,4])) %>%
  left_join(mapping[,c(2,3)], by = "barcode_w1") %>%
  filter(!is.na(sampleid))
sample_overlap = beta$sampleid

# Take the overlap sample of biomarker and wgs
pcs = pcs[rownames(pcs) %in% sample_overlap,]
pcs = pcs[order(rownames(pcs)),]
grm = grm[rownames(grm) %in% sample_overlap, rownames(grm) %in% sample_overlap]
grm = grm[order(rownames(grm)),order(colnames(grm))]
beta = beta %>% filter(sampleid %in% rownames(grm))
beta = beta[order(beta$sampleid),]

# pheno data (4642 cpgs)
pheno = cbind(
  data.frame(FID = rep(0, nrow(beta)), IID = beta$sampleid),
  beta[,1:length(sig_cpg)])
pheno <- pheno[match(rownames(pcs), pheno$IID), ]

# covariate: age, gender, wbc, ewas PC1-PC4, and genetic PC1-PC10
pcs.df = data.frame(intercept = rep(1, nrow(pcs)),
                    age = beta$age, sex = beta$sex, CD8T = beta$CD8T, CD4T = beta$CD4T, NK = beta$NK,
                    Bcell = beta$Bcell, Mono = beta$Mono,
                    ewasPC1 = beta$PC1, ewasPC2 = beta$PC2, ewasPC3 = beta$PC3, ewasPC4 = beta$PC4, 
                    PC1 = pcs$PC1, PC2 = pcs$PC2, PC3 = pcs$PC3, PC4 = pcs$PC4, PC5 = pcs$PC5,
                    PC6 = pcs$PC6, PC7 = pcs$PC7, PC8 = pcs$PC8, PC9 = pcs$PC9, PC10 = pcs$PC10)

# keep ID
keep.ids = data.frame(FID = rep(0, nrow(beta)), IID = beta$sampleid)

# manifest
manifest = manifest %>% filter(CPG %in% col_cpg) %>% 
  mutate(chr = as.numeric(substr(Chr, 4, nchar(Chr)))) %>%
  dplyr::select(CPG, chr, Position)

write.table(
  pheno,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno_mqtl.txt",
  quote = FALSE,   
  row.names = FALSE,
  col.names = TRUE, 
  sep = " "      
)
write.table(
  grm,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/grm_mqtl.txt",
  quote = FALSE,   
  row.names = FALSE,
  col.names = FALSE, # GEMMA do not need colnames
  sep = " "   
)
write.table(
  pcs.df,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/covariate_mqtl.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = " " 
)
write.table(
  keep.ids,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/keep.ids.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = " " 
)
write.table(
  manifest,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/manifest.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = " " 
)
