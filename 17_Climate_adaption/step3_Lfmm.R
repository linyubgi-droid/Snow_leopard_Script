library(lfmm)
library(dplyr)

# 1. Paths and parameters
BASE_DIR   <- "/EClab/Project/Threatened_Species/landscape/snowleopard"
INPUT_GEN  <- file.path(BASE_DIR, "02.lfmm/filtered_genotype_matrix.txt")
INPUT_META <- file.path(BASE_DIR, "02.lfmm/indi_df_matched_sorted.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "02.lfmm/result")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

env_vars  <- c("elev", "bio9", "bio11", "bio16", "bio2")
K_value   <- 2    
FDR_THRES <- 0.05 

# 2. Load and align data
gen     <- read.table(INPUT_GEN, sep = " ", row.names = 1, header = TRUE)
indi_df <- read.csv(INPUT_META, row.names = 1)

indi_df <- indi_df[rownames(gen), , drop = FALSE]
Y_mat   <- t(gen)

# 3. LFMM Analysis loop
all_candidates <- list()

for (env in env_vars) {
  X_vec <- as.numeric(indi_df[[env]])
  
  # Run LFMM Ridge and calibrate P-values
  mod_lfmm <- lfmm_ridge(Y = Y_mat, X = X_vec, K = K_value)
  pv_calib <- lfmm_test(Y = Y_mat, X = X_vec, lfmm = mod_lfmm, calibrate = "gif")
  
  pvalues        <- pv_calib$calibrated.pvalue
  names(pvalues) <- colnames(gen)
  qvalues        <- p.adjust(pvalues, method = "BH")
  
  # Format results
  res_df <- data.frame(
    SNP       = names(pvalues),
    Env_Var   = env,
    p.value   = pvalues,
    q.value   = qvalues,
    row.names = NULL
  )
  
  write.csv(res_df, file.path(OUTPUT_DIR, paste0("lfmm_all_", env, ".csv")), row.names = FALSE)
  
  # Store significant SNPs directly (even if 0 rows, R handles it gracefully)
  all_candidates[[env]] <- res_df[res_df$q.value < FDR_THRES, ]
}

# 4. Aggregate and save final candidate list
final_candidates <- do.call(rbind, all_candidates)
write.csv(final_candidates, file.path(OUTPUT_DIR, "lfmm_candidate_snps_fdr0.05.csv"), row.names = FALSE)
