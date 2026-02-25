library(dplyr)
library(purrr)
library(ggplot2)
setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/")

# 1. Load your index file
sig_cpgs <- read.table("sig_cpgs.txt", header = TRUE, stringsAsFactors = FALSE)
results_dir <- "MR_rho_results/"

# 2. Main processing function
process_results <- function(df) {
  results_list <- pmap(list(df$CpG, df$biomarker, df$chr), function(cpg, bio, ch) {
    # Construct the filename using the period format
    file_path <- file.path(results_dir, paste0(cpg, ".", bio, ".chr", ch, ".RDS"))

    # Handle missing files (in case the pipeline skipped some)
    if (!file.exists(file_path)) {
      return(list(status = "Time Limit", data = NULL, info = c(cpg, bio, ch)))
    }
    
    res_list <- readRDS(file_path)
    res = res_list[["result"]]
    SuSiE = res_list[["SuSiE"]]
    
    # Categorize
    if(SuSiE == -1) {
      return(list(status = "SuSiE No converging", data = NULL, info = c(cpg, bio, ch, SuSiE)))
    } else {
      if (identical(res, "NULL")) {
        return(list(status = "NULL (Error)", data = NULL, info = c(cpg, bio, ch, SuSiE)))
      } else if (identical(res, "No converging")) {
        return(list(status = "No converging", data = NULL, info = c(cpg, bio, ch, SuSiE)))
      } else if (identical(res, "No variables")) {
        return(list(status = "No variables", data = NULL, info = c(cpg, bio, ch, SuSiE)))
      } else if (is.data.frame(res)) {
        # Add the metadata to the result row
        res_full <- cbind(CpG = cpg, biomarker = bio, chr = ch, SuSiE = SuSiE, res)
        return(list(status = "Success", data = res_full, info = c(cpg, bio, ch, SuSiE)))
      } else{
        return(list(status = "Unknown", data = NULL, info = c(cpg, bio, ch, SuSiE)))
      }
    }

  })
  
  return(results_list)
}

# Run the processing
all_outcomes <- process_results(sig_cpgs)

# --- Task 1: Summary Table ---
status_vec <- map_chr(all_outcomes, ~ .x[["status"]] %||% "Unknown2")
status_counts <- table(status_vec)
print("Summary of Outcomes:")
print(status_counts)

# --- Task 2: Error Data Frame ---
error_df <- map_dfr(all_outcomes, function(x) {
  if (x$status == "NULL (Error)") {
    return(data.frame(CpG = x$info[1], biomarker = x$info[2], chr = x$info[3]))
  }
})

# --- Task 3: Success Data Frame & Plotting ---
success_df <- map_dfr(all_outcomes, ~ .x$data)
success_df <- success_df %>%
  mutate(
    se = sqrt(gamma_total_var),
    z_score = gamma / se,
    p_value = 2 * pnorm(abs(z_score), lower.tail = FALSE)
  )

# Histogram of gamma_total_var
p <- ggplot(success_df, aes(x = -log10(p_value))) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Gamma Total Variance",
       x = "-log10 p-value",
       y = "Count")
       

print(p)
