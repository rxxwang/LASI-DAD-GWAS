grm = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/Relatedness_check/Relatedness_check_zheng/grm_2680_genesis.csv")
pcs = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/PCA/PCA_analysis/LASI_only/LASI_PCs_2680_kinThresh0.044.csv")
#unrel = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/PCA/PCA_analysis/LASI_only/unrel_smps_cutoff0.044_N2569_pcair.txt")
mapping = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/ID_mapping/WGS_methy_map.csv")
library(haven)
library(dplyr)
biomarker = read_dta("/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_methylID.dta")
adbio_col = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final")
# biomarker data is using barcode. Add sampleid from mapping file
biomarker <- biomarker %>% filter(inw1_adbio == 1 & barcode_w1 != "") %>% 
  dplyr::select(barcode_w1, w1abeta40_final, w1abeta42_final, w1gfap_final, w1nfl_final, w1ptau_final, w1totaltau_final) %>%
  left_join(mapping[,c(2,3)], by = "barcode_w1") %>%
  filter(!is.na(sampleid)) %>%
  filter(!if_all(all_of(adbio_col), is.na))
sample_overlap = biomarker$sampleid

# Take the overlap sample of biomarker and wgs
pcs = pcs[rownames(pcs) %in% sample_overlap,]
grm = grm[rownames(grm) %in% sample_overlap, rownames(grm) %in% sample_overlap]
biomarker = biomarker %>% filter(sampleid %in% rownames(pcs))

# read gender and age data
gender = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/Phenotype/Pheno_Jan_19_2022_WGS/Wei_gws_01_2022_gender_correct.csv") %>%
    mutate(barcode_w1 = barcode,
           age = r1agey,
           gender = case_when(ragender == "2.Female" ~ 1,
                              ragender == "1.Male" ~ 0)) %>%
    left_join(biomarker, by = "barcode_w1") %>%
    mutate(IID = sampleid) %>%
    select(IID, age, gender) %>%
    filter(!is.na(IID))

# pheno data (5 biomarkers)
pheno = data.frame(FID = rep(0, nrow(biomarker)), IID = biomarker$sampleid, 
                   gfap = biomarker$w1gfap_final, nfl = biomarker$w1nfl_final, 
                   ptau = biomarker$w1ptau_final, totaltau = biomarker$w1totaltau_final,
                   abeta_ratio = biomarker$w1abeta42_final/biomarker$w1abeta40_final) 
pheno <- pheno[match(rownames(pcs), pheno$IID), ]

# covariate: gender and PC1-PC10
pcs.df = data.frame(FID = rep(0, nrow(pcs)), IID = rownames(pcs), intercept = rep(1, nrow(pcs)),
                 PC1 = pcs$PC1, PC2 = pcs$PC2, PC3 = pcs$PC3, PC4 = pcs$PC4, PC5 = pcs$PC5,
                 PC6 = pcs$PC6, PC7 = pcs$PC7, PC8 = pcs$PC8, PC9 = pcs$PC9, PC10 = pcs$PC10) %>%
    left_join(gender, by = "IID")
pcs.df = pcs.df[,-(1:2)]

write.table(
  pheno,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/pheno_gwas.txt",
  quote = FALSE,   
  row.names = FALSE,
  col.names = TRUE, 
  sep = " "      
)
write.table(
  grm,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/grm.txt",
  quote = FALSE,   
  row.names = FALSE,
  col.names = FALSE, # GEMMA do not need colnames
  sep = " "   
)
write.table(
  pcs.df,
  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/covariate.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE, # GEMMA do not need colnames
  sep = " " 
)
