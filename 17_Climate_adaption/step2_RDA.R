library(vegan)     
library(robust)    

# 1. Paths and parameters
BASE_DIR   <- "/EClab/Project/Threatened_Species/landscape/snowleopard"
INPUT_GEN  <- file.path(BASE_DIR, "02.lfmm/filtered_genotype_matrix.txt")
INPUT_META <- file.path(BASE_DIR, "03.RDA/result/indi_df_with_popid.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "03.RDA/result")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

env_vars <- c("elev", "bio9", "bio11", "bio16", "bio2")
K_axes   <- length(env_vars) 

# 2. Load and align data
gen     <- read.table(INPUT_GEN, sep = " ", header = TRUE, row.names = 1)
indi_df <- read.csv(INPUT_META, row.names = 1)

indi_df <- indi_df[colnames(gen), , drop = FALSE]
gen.imp <- t(gen)
pred    <- indi_df[, env_vars, drop = FALSE]

# 3. Run RDA model
rda_model <- rda(gen.imp ~ ., data = pred, scale = TRUE)
load_rda  <- scores(rda_model, choices = 1:K_axes, display = "species")
write.table(load_rda, file.path(OUTPUT_DIR, "rda_snp_loadings.txt"), quote = FALSE)

# 4. Significance testing (Mahalanobis distance)
resscale <- apply(load_rda, 2, scale) 
resmaha  <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist

lambda   <- median(resmaha) / qchisq(0.5, df = K_axes) 
p_values <- pchisq(resmaha / lambda, K_axes, lower.tail = FALSE)

res_rda  <- data.frame(SNP = rownames(load_rda), p.value = p_values)
write.csv(res_rda, file.path(OUTPUT_DIR, "rda_all_snps_pvalues.csv"), row.names = FALSE)

# 5. Extract adaptive loci (P < 0.05)
candidates <- res_rda[res_rda$p.value < 0.05, ]
write.csv(candidates, file.path(OUTPUT_DIR, "rda_candidate_snps_p0.05.csv"), row.names = FALSE)
