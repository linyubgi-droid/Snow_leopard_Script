# -----------------------------------------------------------------------------
# Ensemble Species Distribution Modeling (SDM)
# -----------------------------------------------------------------------------

rm(list = ls())

library(sdm)
library(usdm)      
library(terra)     
library(rJava)     

# 1. Define working directories
BASE_DIR <- "/EClab/Project/Threatened_Species/landscape/snowleopard"
DIR_CURRENT <- file.path(BASE_DIR, "../bio/worldclim_current/2.5m")
DIR_FUTURE <- file.path(BASE_DIR, "../bio/worldclim_future/2.5m")

# 2. Configure MaxEnt environment
jar_dir <- file.path(system.file(package = "sdm"), "java")
dir.create(jar_dir, showWarnings = FALSE)
file.copy(file.path(BASE_DIR, "03.SDM/maxent/maxent.jar"), file.path(jar_dir, "maxent.jar"), overwrite = TRUE)

# 3. Load and format current environmental predictors
elev_cur  <- rast(file.path(DIR_CURRENT, "elevation.tif"))
bio9_cur  <- rast(file.path(DIR_CURRENT, "bio9.tif"))
bio11_cur <- rast(file.path(DIR_CURRENT, "bio11.tif"))
bio16_cur <- rast(file.path(DIR_CURRENT, "bio16.tif"))
bio2_cur  <- rast(file.path(DIR_CURRENT, "bio2.tif"))

env_current <- c(elev_cur, bio9_cur, bio11_cur, bio16_cur, bio2_cur)
names(env_current) <- c("elev", "bio9", "bio11", "bio16", "bio2")

# 4. Load occurrence data and prepare SDM data object
occ_clean <- read.csv(file.path(BASE_DIR, "snow_leopard_thinned.csv")) 
d <- sdmData(formula = Species ~ ., train = occ_clean, predictors = env_current, bg = list(n = 10000))

# 5. Train ensemble model (5 algorithms, 5-fold cross-validation)
m <- sdm(Species ~ ., data = d, 
         methods = c('svm', 'brt', 'gam', 'rf', 'maxent'), 
         replication = 'cv', cv.folds = 5, n = 1)

# 6. Project to current climate
p_current <- ensemble(m, newdata = env_current, setting = list(method = 'weighted', stat = 'TSS'))
writeRaster(p_current, file.path(BASE_DIR, "03.SDM/Ensemble_Current.tif"), overwrite = TRUE)

# 7. Loop through future scenarios and project
folders <- c("ssp126_2050", "ssp126_2090", "ssp585_2050", "ssp585_2090")

for(folder in folders) {
  
  file_prefix <- toupper(folder)
  path_prefix <- file.path(DIR_FUTURE, folder, paste0(file_prefix, "_"))
  
  bio9_fut  <- rast(paste0(path_prefix, "9.tif"))
  bio11_fut <- rast(paste0(path_prefix, "11.tif"))
  bio16_fut <- rast(paste0(path_prefix, "16.tif"))
  bio2_fut  <- rast(paste0(path_prefix, "2.tif"))
  
  env_future <- c(elev_cur, bio9_fut, bio11_fut, bio16_fut, bio2_fut)
  names(env_future) <- c("elev", "bio9", "bio11", "bio16", "bio2")
  
  p_future <- ensemble(m, newdata = env_future, setting = list(method = 'weighted', stat = 'TSS'))
  
  output_filename <- file.path(BASE_DIR, "03.SDM", paste0("Ensemble_", file_prefix, ".tif"))
  writeRaster(p_future, output_filename, overwrite = TRUE)
}
