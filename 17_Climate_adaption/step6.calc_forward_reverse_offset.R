library(gradientForest)
library(raster)
library(geosphere)  
library(doParallel)
library(fields)     

# 1. Define Paths and Parameters
BASE_DIR   <- "/EClab/Project/Threatened_Species/landscape"
GF_DIR     <- file.path(BASE_DIR, "snowleopard/04.gf")
SCENARIO   <- "ssp585_2090" # This can be parameterized via commandArgs()
FUT_DIR    <- file.path(BASE_DIR, "bio/worldclim_future/2.5m", SCENARIO)

env_select <- c("elev", "bio9", "bio11", "bio16", "bio2")
sdm_thresh <- 0.2

# 2. Data Loading & Preparation
input_data  <- readRDS(file.path(GF_DIR, "gf_input_data.rds"))
gf_result   <- readRDS(file.path(GF_DIR, "gf_result_object.rds"))
BIO_current <- input_data$BIO

SDM       <- raster(file.path(GF_DIR, "snowleopard.asc")) 
pts_mask  <- rasterToPoints(SDM)
point     <- pts_mask[pts_mask[, 3] > sdm_thresh, 1:2] 

# 3. Extract Climate Data and Project into Genomic Space
env_curr <- data.frame(extract(BIO_current[[env_select]], point))

fut_files <- list.files(FUT_DIR, pattern = "\\.tif$", full.names = TRUE)
BIO_future <- stack(fut_files)
names(BIO_future) <- paste0("BIO", 1:nlayers(BIO_future))
env_fut <- data.frame(extract(BIO_future[[env_select]], point))

# Clean NAs to ensure perfect matching
valid_idx <- complete.cases(env_curr) & complete.cases(env_fut)
pts_valid <- point[valid_idx, ]

# Project environments to Gradient Forest genomic space
df_curr_trans <- predict(gf_result, env_curr[valid_idx, ])
df_fut_trans  <- predict(gf_result, env_fut[valid_idx, ])

# 4. Local Offset Calculation
# Local: Euclidean distance between current and future at the exact same location
local_offset <- sqrt(rowSums((df_fut_trans - df_curr_trans)^2))

# 5. Forward & Reverse Offset Calculation (Parallel Processing)
numCores <- min(detectCores() - 1, 8) 
cl <- makeCluster(max(1, numCores))
registerDoParallel(cl)

# A. Reverse Offset: From each Future Cell back to the best Current Cell
reverse_offset_list <- foreach(i = 1:nrow(df_fut_trans), 
                               .packages=c("fields","geosphere")) %dopar% {
    
    # Calculate genomic distance from this specific future cell to ALL current cells
    gfOffset <- c(fields::rdist(df_fut_trans[i, , drop=FALSE], df_curr_trans)) 
    
    # Find the current cell(s) with the minimum genomic distance
    min_idx <- which(gfOffset == min(gfOffset))
    
    # If multiple minimums exist, break tie by minimum geographic distance
    if (length(min_idx) > 1) {
        geo_dists <- geosphere::distGeo(pts_valid[i, , drop=FALSE], pts_valid[min_idx, , drop=FALSE])
        min_idx <- min_idx[which.min(geo_dists)]
    }
    
    return(gfOffset[min_idx])
}

# B. Forward Offset: From each Current Cell forward to the best Future Cell
forward_offset_list <- foreach(i = 1:nrow(df_curr_trans), 
                               .packages=c("fields","geosphere")) %dopar% {
                               
    gfOffset <- c(fields::rdist(df_curr_trans[i, , drop=FALSE], df_fut_trans)) 
    min_idx <- which(gfOffset == min(gfOffset))
    
    if (length(min_idx) > 1) {
        geo_dists <- geosphere::distGeo(pts_valid[i, , drop=FALSE], pts_valid[min_idx, , drop=FALSE])
        min_idx <- min_idx[which.min(geo_dists)]
    }
    
    return(gfOffset[min_idx])
}

stopCluster(cl)

# Compile results
LFR_df <- data.frame(
    lon     = pts_valid[, "x"],
    lat     = pts_valid[, "y"],
    Local   = local_offset,
    Forward = unlist(forward_offset_list),
    Reverse = unlist(reverse_offset_list)
)

write.csv(LFR_df, file.path(GF_DIR, paste0("LFR_Offset_", SCENARIO, ".csv")), row.names = FALSE)

# 6. Normalize and Generate RGB GeoTIFF
# Normalizes a vector to 0-255 range for RGB mapping
norm_to_rgb <- function(v) {
    if(max(v, na.rm=T) == min(v, na.rm=T)) return(rep(255, length(v)))
    round((v - min(v, na.rm=T)) / (max(v, na.rm=T) - min(v, na.rm=T)) * 255)
}

rgb_df <- data.frame(
    x = pts_valid[, "x"], y = pts_valid[, "y"],
    R = norm_to_rgb(LFR_df$Local),
    G = norm_to_rgb(LFR_df$Forward),
    B = norm_to_rgb(LFR_df$Reverse)
)

# Convert to Raster Stack
r_layer <- raster(SDM); values(r_layer) <- NA; r_layer[cellFromXY(r_layer, rgb_df[, 1:2])] <- rgb_df$R
g_layer <- raster(SDM); values(g_layer) <- NA; g_layer[cellFromXY(g_layer, rgb_df[, 1:2])] <- rgb_df$G
b_layer <- raster(SDM); values(b_layer) <- NA; b_layer[cellFromXY(b_layer, rgb_df[, 1:2])] <- rgb_df$B

rgb_stack <- stack(r_layer, g_layer, b_layer)
names(rgb_stack) <- c("Local_R", "Forward_G", "Reverse_B")

writeRaster(rgb_stack, file.path(GF_DIR, paste0("LFR_RGB_", SCENARIO, ".tif")), 
            overwrite = TRUE, datatype='INT2U', format = "GTiff", options="COMPRESS=LZW")
