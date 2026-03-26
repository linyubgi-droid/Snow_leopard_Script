library(gradientForest)
library(raster)

# 1. Define Paths (Modify BASE_DIR for local environment)
BASE_DIR   <- "/EClab/Project/Threatened_Species/landscape"
BIO_DIR    <- file.path(BASE_DIR, "bio/worldclim_current/2.5m")
GF_DIR     <- file.path(BASE_DIR, "snowleopard/04.gf")

env_select <- c("elev", "bio9", "bio11", "bio16", "bio2") 

# 2. Load and Format Environmental Raster Stack
bio_files <- list.files(BIO_DIR, pattern = "^bio[0-9]+\\.tif$", full.names = TRUE)
elev_file <- file.path(BIO_DIR, "elevation.tif")
BIO       <- stack(c(bio_files, elev_file))

names(BIO) <- tools::file_path_sans_ext(basename(c(bio_files, elev_file)))

# 3. Process Coordinates and Extract Environments
coords <- read.csv(file.path(GF_DIR, "sample_coords.csv"), row.names = 1)[, c("lon", "lat")]
gf_env <- na.omit(data.frame(raster::extract(BIO, coords)))
rownames(gf_env) <- rownames(coords)

# 4. Process PLINK Genomic Dosages
raw_data <- read.table(file.path(GF_DIR, "snowleopard_dosages.raw"), header = TRUE, row.names = NULL)
gf_alle  <- raw_data[, 7:ncol(raw_data)]
rownames(gf_alle) <- raw_data$IID
colnames(gf_alle) <- gsub("_[A-Z]$", "", colnames(gf_alle))

# 5. Align Datasets
common_samples <- intersect(rownames(gf_env), rownames(gf_alle))
gf_env  <- gf_env[common_samples, env_select]
gf_alle <- gf_alle[common_samples, ]

# Save aligned input data for downstream offset calculations
saveRDS(list(gf_env = gf_env, gf_alle = gf_alle, BIO = BIO), 
        file.path(GF_DIR, "gf_input_data.rds"))

# 6. Train Gradient Forest Model
set.seed(123)
gf_model <- gradientForest(
  cbind(gf_env, gf_alle),
  predictor.vars = colnames(gf_env),
  response.vars  = colnames(gf_alle), 
  ntree          = 400,
  maxLevel       = 0.05, 
  trace          = FALSE, 
  corr.threshold = 0.40
)

# Save the trained model object
saveRDS(gf_model, file.path(GF_DIR, "gf_result_object.rds"))
