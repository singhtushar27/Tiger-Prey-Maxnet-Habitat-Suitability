# ============================================================
# TIGER PREY RESOURCE HABITAT SUITABILITY USING MAXNET
# ============================================================
# Purpose:
#   Reproducible R implementation of the MaxEnt-style workflow
#   used to model combined tiger-prey resource suitability.
#
# Important methodological note:
#   The original analysis used standalone MaxEnt. This script
#   implements a reproducible R equivalent using maxnet and
#   should not be described as the original standalone-MaxEnt run.
#
# Model settings:
#   Feature classes : lqph
#   Regularization  : 1
#   Output          : cloglog suitability (0-1)
#   Background      : 10,000 complete environmental cells
#   Validation      : 5-fold spatial cross-validation
# ============================================================

# -----------------------------
# 1. PACKAGES
# -----------------------------
library(terra)
library(maxnet)
library(pROC)

# -----------------------------
# 2. PROJECT PATHS
# -----------------------------
project_dir <- "D:/Tushar/Research_Papers/Animal_Corridor_Pandu_da/R_MaxEnt"

env_file <- file.path(
  project_dir, "environmental_layers", "selected_environmental_stack.tif"
)
occ_file <- file.path(
  "D:/Tushar/Research_Papers/Animal_Corridor_Pandu_da",
  "Species", "Book1.csv"
)
study_area_file <- file.path(
  "D:/Tushar/Research_Papers/Animal_Corridor_Pandu_da",
  "Projected_study_area", "Study_Area_UTM.shp"
)

model_dir <- file.path(project_dir, "models")
prediction_dir <- file.path(project_dir, "predictions")
results_dir <- file.path(project_dir, "results")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. ENVIRONMENTAL PREDICTORS
# -----------------------------
env <- rast(env_file)

expected_layers <- c(
  "bio2", "bio3", "bio7", "bio9", "bio10", "bio12", "bio17", "bio18",
  "NDVI", "slope", "dist_river", "dist_road", "dist_railway", "LULC"
)

stopifnot(identical(names(env), expected_layers))

# -----------------------------
# 4. OCCURRENCE DATA
# -----------------------------

occ <- read.csv(
  occ_file,
  stringsAsFactors = FALSE
)

stopifnot(all(c("Animal", "Lat", "Long") %in% names(occ)))

occ <- occ[complete.cases(occ[, c("Long", "Lat")]), ]

# Book1.csv contains projected coordinates:
# Long = X and Lat = Y, despite the column names.

# Create point vector from X/Y coordinates
occ_vect <- terra::vect(
  occ,
  geom = c("Long", "Lat"),
  crs = terra::crs(env),
  keepgeom = FALSE
)

# Preserve animal identity
occ_vect$Animal <- occ$Animal

# -----------------------------
# 5. ONE RECORD PER ENVIRONMENTAL CELL
# -----------------------------
occ_cells <- cellFromXY(env[[1]], crds(occ_vect))
valid <- !is.na(occ_cells)

occ_vect <- occ_vect[valid, ]
occ_cells <- occ_cells[valid]

keep <- !duplicated(occ_cells)
prey_vect_cell <- occ_vect[keep, ]
presence_cells <- occ_cells[keep]

cat("Original records:", nrow(occ), "\n")
cat("Unique modelling cells:", length(presence_cells), "\n")

# -----------------------------
# 6. PRESENCE ENVIRONMENTAL VALUES
# -----------------------------
presence_env <- terra::extract(env,prey_vect_cell)
presence_df <- presence_env[, -1]

if (anyNA(presence_df)) {
  stop("Presence cells contain missing environmental values.")
}

# -----------------------------
# 7. BACKGROUND SAMPLING
# -----------------------------
study_area <- vect(study_area_file)
study_area <- project(study_area, crs(env))

study_mask <- rasterize(
  study_area, env[[1]], field = 1, background = NA
)

available_cells <- which(!is.na(values(study_mask)))
complete_cells <- which(complete.cases(values(env)))

background_candidates <- intersect(available_cells, complete_cells)
background_candidates <- setdiff(background_candidates, presence_cells)

set.seed(123)
n_background <- min(10000, length(background_candidates))
background_cells <- sample(background_candidates, n_background, replace = FALSE)

background_xy <- xyFromCell(env[[1]], background_cells)
background_vect <- vect(background_xy, type = "points", crs = crs(env))

background_env <- terra::extract(env, background_vect)
background_df <- background_env[, -1]

if (anyNA(background_df)) {
  stop("Background cells contain missing environmental values.")
}

# -----------------------------
# 8. MAXNET DATASET
# -----------------------------
# LULC is categorical; all other predictors remain continuous.
lulc_levels <- sort(unique(c(
  as.character(presence_df$LULC),
  as.character(background_df$LULC)
)))

presence_df$LULC <- factor(presence_df$LULC, levels = lulc_levels)
background_df$LULC <- factor(background_df$LULC, levels = lulc_levels)

presence_df$presence <- 1
background_df$presence <- 0

maxent_data <- rbind(presence_df, background_df)

cat("\nPresence/background counts:\n")
print(table(maxent_data$presence))

# -----------------------------
# 9. RANDOM 70/30 VALIDATION
# -----------------------------

set.seed(123)

presence_id <- sample(
  seq_len(nrow(presence_df)),
  size = round(0.70 * nrow(presence_df))
)

background_id <- sample(
  seq_len(nrow(background_df)),
  size = round(0.70 * nrow(background_df))
)

train_random <- rbind(
  presence_df[presence_id, ],
  background_df[background_id, ]
)

test_random <- rbind(
  presence_df[-presence_id, ],
  background_df[-background_id, ]
)

x_train <- train_random[, setdiff(names(train_random), "presence")]
y_train <- train_random$presence

x_test <- test_random[, setdiff(names(test_random), "presence")]
y_test <- test_random$presence

x_train$LULC <- factor(
  x_train$LULC,
  levels = lulc_levels
)

x_test$LULC <- factor(
  x_test$LULC,
  levels = lulc_levels
)

set.seed(123)

mx_model_eval <- maxnet(
  p = y_train,
  data = x_train,
  f = maxnet.formula(
    p = y_train,
    data = x_train,
    classes = "lqph"
  ),
  regmult = 1
)

random_pred <- as.numeric(
  predict(
    mx_model_eval,
    x_test,
    type = "cloglog"
  )
)

random_auc <- as.numeric(
  auc(
    roc(
      y_test,
      random_pred,
      quiet = TRUE
    )
  )
)

# Find threshold maximizing TSS

thresholds <- seq(0, 1, by = 0.001)

tss_values <- sapply(
  thresholds,
  function(threshold) {
    
    predicted_class <- ifelse(
      random_pred >= threshold,
      1,
      0
    )
    
    TP <- sum(predicted_class == 1 & y_test == 1)
    TN <- sum(predicted_class == 0 & y_test == 0)
    FP <- sum(predicted_class == 1 & y_test == 0)
    FN <- sum(predicted_class == 0 & y_test == 1)
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    
    sensitivity + specificity - 1
  }
)

best_random <- which.max(tss_values)

random_threshold <- thresholds[best_random]
random_tss <- tss_values[best_random]

cat("\nRandom validation AUC:",
    round(random_auc, 4), "\n")

cat("Random validation TSS:",
    round(random_tss, 4), "\n")

cat("Random validation threshold:",
    random_threshold, "\n")

# -----------------------------
# 10. FINAL MAXNET MODEL
# -----------------------------
x_final <- maxent_data[, setdiff(names(maxent_data), "presence")]
y_final <- maxent_data$presence

set.seed(123)
mx_model_final <- maxnet(
  p = y_final,
  data = x_final,
  f = maxnet.formula(
    p = y_final,
    data = x_final,
    classes = "lqph"
  ),
  regmult = 1
)

saveRDS(
  mx_model_final,
  file.path(model_dir, "Tiger_Prey_Maxnet_final.rds")
)

# -----------------------------
# 11. FIVE-FOLD SPATIAL VALIDATION
# -----------------------------
set.seed(123)

presence_xy <- crds(prey_vect_cell)
background_xy <- crds(background_vect)

# Cluster presence locations into five spatial groups.
presence_clusters <- kmeans(
  presence_xy, centers = 5, nstart = 50
)

presence_df$spatial_fold <- presence_clusters$cluster

nearest_cluster <- function(x, centers) {
  d <- sapply(seq_len(nrow(centers)), function(i) {
    (x[, 1] - centers[i, 1])^2 +
      (x[, 2] - centers[i, 2])^2
  })
  max.col(-d, ties.method = "first")
}

background_df$spatial_fold <- nearest_cluster(
  background_xy, presence_clusters$centers
)

spatial_results <- data.frame()
spatial_predictions <- list()

for (fold_id in 1:5) {

  train_sp <- rbind(
    presence_df[presence_df$spatial_fold != fold_id, ],
    background_df[background_df$spatial_fold != fold_id, ]
  )

  test_sp <- rbind(
    presence_df[presence_df$spatial_fold == fold_id, ],
    background_df[background_df$spatial_fold == fold_id, ]
  )

  x_train <- train_sp[, setdiff(names(train_sp), c("presence", "spatial_fold"))]
  y_train <- train_sp$presence
  x_test <- test_sp[, setdiff(names(test_sp), c("presence", "spatial_fold"))]
  y_test <- test_sp$presence

  x_train$LULC <- factor(x_train$LULC, levels = lulc_levels)
  x_test$LULC <- factor(x_test$LULC, levels = lulc_levels)

  set.seed(123 + fold_id)
  mx_fold <- maxnet(
    p = y_train,
    data = x_train,
    f = maxnet.formula(p = y_train, data = x_train, classes = "lqph"),
    regmult = 1
  )

  pred <- as.numeric(predict(mx_fold, x_test, type = "cloglog"))

  auc_value <- as.numeric(auc(roc(y_test, pred, quiet = TRUE)))
  
  # Calculate TSS across thresholds
  
  thresholds <- seq(0, 1, by = 0.001)
  
  tss_values <- sapply(
    thresholds,
    function(threshold) {
      
      predicted_class <- ifelse(
        pred >= threshold,
        1,
        0
      )
      
      TP <- sum(predicted_class == 1 & y_test == 1)
      TN <- sum(predicted_class == 0 & y_test == 0)
      FP <- sum(predicted_class == 1 & y_test == 0)
      FN <- sum(predicted_class == 0 & y_test == 1)
      
      sensitivity <- TP / (TP + FN)
      specificity <- TN / (TN + FP)
      
      sensitivity + specificity - 1
    }
  )
  
  best_tss <- which.max(tss_values)
  
  fold_threshold <- thresholds[best_tss]
  fold_tss <- tss_values[best_tss]

  spatial_results <- rbind(
    spatial_results,
    data.frame(
      fold = fold_id,
      test_presence = sum(y_test == 1),
      test_background = sum(y_test == 0),
      AUC = auc_value,
      threshold = fold_threshold,
      TSS = fold_tss
    )
  )

  spatial_predictions[[fold_id]] <- data.frame(
    observed = y_test,
    prediction = pred
  )
}

print(spatial_results)

mean_spatial_auc <- mean(spatial_results$AUC)
sd_spatial_auc <- sd(spatial_results$AUC)
mean_spatial_tss <- mean(spatial_results$TSS)
sd_spatial_tss <- sd(spatial_results$TSS)

cat("\nMean spatial AUC:", round(mean_spatial_auc, 4), "\n")
cat("SD spatial AUC:", round(sd_spatial_auc, 4), "\n")
cat("Mean spatial TSS:", round(mean_spatial_tss, 4), "\n")
cat("SD spatial TSS:", round(sd_spatial_tss, 4), "\n")

write.csv(
  spatial_results[, c(
    "fold",
    "test_presence",
    "test_background",
    "AUC"
  )],
  file.path(
    results_dir,
    "spatial_validation_AUC.csv"
  ),
  row.names = FALSE
)

write.csv(
  spatial_results[, c(
    "fold",
    "test_presence",
    "test_background",
    "threshold",
    "TSS"
  )],
  file.path(
    results_dir,
    "spatial_validation_TSS.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# 12. SAVE MODEL SUMMARY
# -----------------------------

model_summary <- data.frame(
  Metric = c(
    "Presence cells",
    "Background cells",
    "Number of predictors",
    "Random test AUC",
    "Random test TSS",
    "Spatial mean AUC",
    "Spatial SD AUC",
    "Spatial mean TSS",
    "Spatial SD TSS"
  ),
  Value = c(
    length(presence_cells),
    n_background,
    nlyr(env),
    random_auc,
    random_tss,
    mean_spatial_auc,
    sd_spatial_auc,
    mean_spatial_tss,
    sd_spatial_tss
  )
)

write.csv(
  model_summary,
  file.path(
    results_dir,
    "model_summary.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# 13. CONTINUOUS SUITABILITY MAP
# -----------------------------
prey_suitability <- terra::predict(
  env,
  mx_model_final,
  type = "cloglog",
  na.rm = TRUE
)

prey_suitability_file <- file.path(
  prediction_dir,
  "Tiger_Prey_Resource_Suitability_Maxnet.tif"
)

writeRaster(
  prey_suitability,
  prey_suitability_file,
  overwrite = TRUE,
  wopt = list(
    datatype = "FLT4S",
    gdal = c("COMPRESS=LZW")
  )
)

# -----------------------------
# 14. OUTPUT SUMMARY
# -----------------------------

cat("\n============================================\n")
cat("TIGER PREY MAXNET WORKFLOW COMPLETE\n")
cat("============================================\n")

cat("Presence cells:", length(presence_cells), "\n")
cat("Background cells:", n_background, "\n")
cat("Predictors:", nlyr(env), "\n")

cat("\nRandom validation:\n")
cat("AUC:", round(random_auc, 4), "\n")
cat("TSS:", round(random_tss, 4), "\n")

cat("\nSpatial validation:\n")
cat("Mean AUC:", round(mean_spatial_auc, 4), "\n")
cat("SD AUC:", round(sd_spatial_auc, 4), "\n")
cat("Mean TSS:", round(mean_spatial_tss, 4), "\n")
cat("SD TSS:", round(sd_spatial_tss, 4), "\n")

cat("\nSuitability raster:\n")
cat(prey_suitability_file, "\n")

cat("\nFinal model:\n")
cat(
  file.path(
    model_dir,
    "Tiger_Prey_Maxnet_final.rds"
  ),
  "\n"
)

cat("============================================\n")
