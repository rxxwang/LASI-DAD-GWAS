grm = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/Relatedness_check/Relatedness_check_zheng/grm_2680_genesis.csv")
pcs = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/PCA/PCA_analysis/LASI_only/LASI_PCs_2680_kinThresh0.044.csv")
#unrel = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/PCA/PCA_analysis/LASI_only/unrel_smps_cutoff0.044_N2569_pcair.txt")
mapping = read.csv("/net/orion/skardia_lab/clubhouse/research/projects/LASI/WGS_data_Dec_2021/QC_UM/ID_mapping/WGS_methy_map.csv")
library(haven)
library(dplyr)
library(lme4)
library(tidyr)
biomarker = read_dta("/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_methylID.dta")
adbio_col = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final")
# biomarker data is using barcode. Add sampleid from mapping file
biomarker <- biomarker %>% filter(inw1_adbio == 1 & barcode_w1 != "") %>% 
  dplyr::select(barcode_w1, w1abeta40_final, w1abeta42_final, w1gfap_final, w1nfl_final, w1ptau_final, w1totaltau_final,
                w1abeta40_plate, w1abeta42_plate, w1gfap_plate, w1nfl_plate, w1ptau_plate, w1totaltau_plate) %>%
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

col_biomarker <-data.frame(plate = c("w1abeta40_plate", "w1abeta42_plate", "w1gfap_plate", "w1nfl_plate", "w1ptau_plate", "w1totaltau_plate"),
                           biomarker = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final"))
for(i in 1:6){
  data = data.frame(biomarker = biomarker[,col_biomarker[i,2]],
                    plate = biomarker[,col_biomarker[i,1]],
                    rownum = 1:nrow(biomarker))
  data = data[complete.cases(data),]
  colnames(data) = c("biomarker", "plate", "rownum")
  data$plate = ifelse(data$plate < 2, 1, data$plate)
  data$plate = factor(data$plate)
  model <- lmer(log(biomarker) ~ (1 | plate), data)
  res <- data.frame(res = scale(residuals(model)), rownum = data$rownum) %>%
    complete(rownum = 1:nrow(biomarker))
  biomarker[[paste0(col_biomarker[i,2], "_res")]] = res$res
}

# pheno data (5 biomarkers)
pheno = data.frame(FID = rep(0, nrow(biomarker)), IID = biomarker$sampleid, 
                   gfap = biomarker$w1gfap_final_res, nfl = biomarker$w1nfl_final_res, 
                   ptau = biomarker$w1ptau_final_res, totaltau = biomarker$w1totaltau_final_res,
                   abeta_ratio = biomarker$w1abeta42_final_res-biomarker$w1abeta40_final_res) 
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
