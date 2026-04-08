dir = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/"
import pandas as pd
config_df = pd.read_csv("sig_cpgs.txt", sep="\t")

rule all:
    input:
        expand(dir + "MR_rho_results2/{CpG}.{biomarker}.chr{chr}.RDS", 
               zip, # zip ensures we take values row-by-row, not a cross-product
               CpG = config_df['CpG'], 
               biomarker = config_df['biomarker'], 
               chr = config_df['chr'])
        # expand(dir + "MR_rho_results3/{CpG}.{biomarker}.chr{chr}.RDS", 
        #        zip, # zip ensures we take values row-by-row, not a cross-product
        #        CpG = config_df['CpG'], 
        #        biomarker = config_df['biomarker'], 
        #        chr = config_df['chr'])

rule CalculateLD:
  input:
    exposure = dir + "mQTL_cleaned_results/{CpG}.cleaned.assoc.txt",
    genotype_bed = dir + "genotype_mqtl/chr{chr}_filtered.bed",
    genotype_bim = dir + "genotype_mqtl/chr{chr}_filtered.bim",
    genotype_fam = dir + "genotype_mqtl/chr{chr}_filtered.fam"
  output: 
    results = dir + "cor_matrix_R/{CpG}.chr{chr}.RDS"
  params:
    CpG = lambda wildcards: wildcards.CpG,
    chr = lambda wildcards: wildcards.chr
  resources:
    mem_mb = 50000,
    runtime = 5000
  script:
    "Calculate_LD.R"

rule FinemappingMR:
  input:
    exposure = dir + "mQTL_cleaned_results/{CpG}.cleaned.assoc.txt",
    cor_matrix = dir + "cor_matrix_R/{CpG}.chr{chr}.RDS",
    outcome = dir + "gwas_results/{biomarker}_chr{chr}.assoc.txt",
    covariate_gwas = dir + "covariate.txt",
    covariate_mqtl = dir + "covariate_mqtl.txt",
    pheno_gwas = dir + "pheno_gwas.txt",
    pheno_mqtl = dir + "pheno_mqtl.txt"
  output: 
    results = dir + "MR_rho_results2/{CpG}.{biomarker}.chr{chr}.RDS"
  params:
    CpG = lambda wildcards: wildcards.CpG,
    biomarker = lambda wildcards: wildcards.biomarker,
    chr = lambda wildcards: wildcards.chr
  resources:
    mem_mb = 20000,
    runtime = 5000
  script:
    "FinemappingMR_sub.R"

# rule CalculateLD2:
#   input:
#     exposure = dir + "mQTL_cleaned_results/{CpG}.cleaned.assoc.txt",
#     genotype_bed = dir + "genotype/chr{chr}_filtered.bed",
#     genotype_bim = dir + "genotype/chr{chr}_filtered.bim",
#     genotype_fam = dir + "genotype/chr{chr}_filtered.fam"
#   output: 
#     results = dir + "cor_matrix_R2/{CpG}.chr{chr}.RDS"
#   params:
#     CpG = lambda wildcards: wildcards.CpG,
#     chr = lambda wildcards: wildcards.chr
#   resources:
#     mem_mb = 50000,
#     runtime = 2000
#   script:
#     "Calculate_LD.R"

# rule FinemappingMR2:
#   input:
#     exposure = dir + "mQTL_cleaned_results/{CpG}.cleaned.assoc.txt",
#     cor_matrix = dir + "cor_matrix_R2/{CpG}.chr{chr}.RDS",
#     outcome = dir + "gwas_results/{biomarker}_chr{chr}.assoc.txt",
#     covariate_gwas = dir + "covariate.txt",
#     covariate_mqtl = dir + "covariate_mqtl.txt",
#     pheno_gwas = dir + "pheno_gwas.txt",
#     pheno_mqtl = dir + "pheno_mqtl.txt"
#   output: 
#     results = dir + "MR_rho_results3/{CpG}.{biomarker}.chr{chr}.RDS"
#   params:
#     CpG = lambda wildcards: wildcards.CpG,
#     biomarker = lambda wildcards: wildcards.biomarker,
#     chr = lambda wildcards: wildcards.chr
#   resources:
#     mem_mb = 30000,
#     runtime = 5000
#   script:
#     "FinemappingMR_sub.R"