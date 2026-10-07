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
  project_dir, "Species", "Book1.csv"
)
study_area_file <- file.path(
  "D:/Tushar/Research_Papers/Animal_Corridor_Pandu_da",
  "Projected_study_area", "Study_Area_UTM.shp"
)

model_dir <- file.path(project_dir, "models")
prediction_dir <- file.path(project_dir, "predictions")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prediction_dir, recursive = TRUE, showWarnings = FALSE)

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

occ_file <- "D:/Tushar/Research_Papers/Animal_Corridor_Pandu_da/Species/Book1.csv"

occ <- read.csv(
  occ_file,
  stringsAsFactors = FALSE
)

stopifnot(all(c("Animal", "Lat", "Long") %in% names(occ)))

occ <- occ[
  complete.cases(occ[, c("Long", "Lat")]),
]
# Book1.csv contains projected coordinates:
# Long = X and Lat = Y, despite the column names.
occ <- occ[complete.cases(occ[, c("Long", "Lat")]), ]

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
cattle_vect_cell <- occ_vect[keep, ]
presence_cells <- occ_cells[keep]

cat("Original records:", nrow(occ), "\n")
cat("Unique modelling cells:", length(presence_cells), "\n")

# -----------------------------
# 6. PRESENCE ENVIRONMENTAL VALUES
# -----------------------------
presence_env <- terra::extract(env, cattle_vect_cell)
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
# 9. FINAL MAXNET MODEL
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
# 10. FIVE-FOLD SPATIAL VALIDATION
# -----------------------------
set.seed(123)

presence_xy <- crds(cattle_vect_cell)
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

  spatial_results <- rbind(
    spatial_results,
    data.frame(
      fold = fold_id,
      test_presence = sum(y_test == 1),
      test_background = sum(y_test == 0),
      AUC = auc_value
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

cat("\nMean spatial AUC:", round(mean_spatial_auc, 4), "\n")
cat("SD spatial AUC:", round(sd_spatial_auc, 4), "\n")

# -----------------------------
# 11. CONTINUOUS SUITABILITY MAP
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
# 12. OUTPUT SUMMARY
# -----------------------------
cat("\n============================================\n")
cat("TIGER PREY MAXNET WORKFLOW COMPLETE\n")
cat("============================================\n")
cat("Presence cells:", length(presence_cells), "\n")
cat("Background cells:", n_background, "\n")
cat("Predictors:", nlyr(env), "\n")
cat("Spatial AUC:", round(mean_spatial_auc, 4),
    "+/-", round(sd_spatial_auc, 4), "\n")
cat("Suitability raster:\n", prey_suitability_file, "\n")
cat("Final model:\n",
    file.path(model_dir, "Tiger_Prey_Maxnet_final.rds"), "\n")
cat("============================================\n")

