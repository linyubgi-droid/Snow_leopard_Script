library(gradientForest)
library(raster)

# 1. Paths and Data Loading
gf_work_dir <- "/EClab/Project/Threatened_Species/landscape/snowleopard/04.gf/"
base_fut_dir <- "/EClab/Project/Threatened_Species/landscape/bio/worldclim_future/2.5m/"

# Load prepared data and GF model
input_data  <- readRDS(paste0(gf_work_dir, "gf_input_data.rds"))
gf_result   <- readRDS(paste0(gf_work_dir, "gf_result_object.rds"))
BIO_current <- input_data$BIO

# Set up habitat mask (SDM suitability > 0.2)
SDM       <- raster(paste0(gf_work_dir, "snowleopard.asc"))
pts_mask  <- rasterToPoints(SDM)
point     <- as.data.frame(pts_mask[pts_mask[,3] > 0.2, ])
env_select <- c("elev", "bio9", "bio11", "bio16", "bio2")

# Predict current genomic composition once
SDM_env_curr <- data.frame(raster::extract(BIO_current[[env_select]], point[, c("x", "y")]))
df_current   <- predict(gf_result, SDM_env_curr)

# 2. Batch Processing for 4 Scenarios
scenarios <- c("ssp126_2050", "ssp126_2090", "ssp585_2050", "ssp585_2090")

for (scen in scenarios) {
  # Setup paths for each scenario
  fut_dir <- file.path(base_fut_dir, scen)
  out_dir <- file.path(gf_work_dir, scen)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load future climate stack
  fut_files <- list.files(fut_dir, pattern = "\\.tif$", full.names = TRUE)
  BIO_future <- stack(fut_files)
  names(BIO_future) <- paste0("BIO", 1:nlayers(BIO_future))
  
  # Extract and predict future genomic composition
  SDM_env_fut <- data.frame(raster::extract(BIO_future[[env_select]], point[, c("x", "y")]))
  
  # Ensure NAs are handled (only points valid in both current and future)
  valid_idx  <- complete.cases(SDM_env_fut)
  curr_pred  <- df_current[valid_idx, ]
  fut_pred   <- predict(gf_result, SDM_env_fut[valid_idx, ])
  pts_valid  <- point[valid_idx, c("x", "y")]
  
  # Calculate Genomic Offset (Euclidean distance)
  offset_val <- sqrt(rowSums((fut_pred - curr_pred)^2))
  df_offset  <- cbind(pts_valid, offset = offset_val)
  
  # Output 1: CSV for raw data
  write.csv(df_offset, file.path(out_dir, paste0("genetic_offset_", scen, ".csv")), row.names = FALSE)
  
  # Output 2: GeoTIFF for ArcGIS visualization
  go_raster <- raster(SDM)
  values(go_raster) <- NA
  go_raster[cellFromXY(go_raster, df_offset[, 1:2])] <- df_offset$offset
  
  writeRaster(go_raster, file.path(out_dir, paste0("genomic_offset_", scen, ".tif")), 
              overwrite = TRUE, datatype = 'FLT4S', format = "GTiff", options = "COMPRESS=LZW")
}
