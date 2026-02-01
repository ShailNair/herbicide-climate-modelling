# ============================================================================
# herbicide discharge - erosivity relationship analysis
# ============================================================================

rm(list = ls())
set.seed(123)

library(tidyverse)
library(data.table)
library(sf)
library(spdep)
library(blockCV)
library(xgboost)
library(car)
library(FNN)
library(zoo)

# --- Configuration ---
analysis_config <- list(
  cv_folds = 5,
  spatial_block_range_m = 200000,
  random_seed = 123,
  vif_threshold = 5,
  min_years_per_location = 5,
  outlier_percentiles = c(0.01, 0.99),
  moran_k_neighbors = 8
)

set.seed(analysis_config$random_seed)

# ============================================================================
# Data loading and processing
# ============================================================================

discharge_raw <- read_csv("DIS_annual.csv", show_col_types = FALSE)

discharge <- discharge_raw %>%
  mutate(
    quality_flag = case_when(
      sum_discharge <= 0 | is.na(sum_discharge) ~ "invalid",
      sum_discharge < quantile(sum_discharge, analysis_config$outlier_percentiles[1], na.rm = TRUE) ~ "low_outlier",
      sum_discharge > quantile(sum_discharge, analysis_config$outlier_percentiles[2], na.rm = TRUE) ~ "high_outlier",
      TRUE ~ "valid"
    )
  ) %>%
  filter(quality_flag == "valid") %>%
  mutate(
    raw_discharge = sum_discharge,
    log_discharge = log10(sum_discharge),
    sqrt_discharge = sqrt(sum_discharge),
    boxcox_discharge = (sum_discharge^0.5 - 1) / 0.5
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

log_discharge <- discharge_raw %>%
  mutate(
    log_sum_discharge = if_else(sum_discharge == 0 | is.na(sum_discharge), 
                                log10(0.001), log10(sum_discharge)),
    log_mean_discharge = if_else(mean_discharge == 0 | is.na(mean_discharge), 
                                 log10(0.001), log10(mean_discharge))
  ) %>%
  filter(
    log_sum_discharge > quantile(log_sum_discharge, 0.01, na.rm = TRUE),
    log_sum_discharge < quantile(log_sum_discharge, 0.99, na.rm = TRUE)
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Load covariates
apr <- fread("annual_apr.csv") %>%
  rename(apr = apr) %>%
  filter(year %in% unique(discharge$year)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rnf <- fread("annual_rnf.csv") %>%
  rename(rnf = rnf) %>%
  filter(year %in% unique(discharge$year)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

ERA5Land <- fread("ERA5Land_2001_2021_erosivity.csv") %>%
  rename(erosivity = data) %>%
  filter(year %in% unique(discharge$year)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rcp45_change <- fread("2010-2050-rcp4.5_final.csv") %>%
  rename(erosivity_change_2050 = data) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# ============================================================================
# Spatial aggregation (50 Km buffer)
# ============================================================================

join_and_aggregate_covariates_safe <- function(discharge_sf, covariate_sf, covariate_name, buffer_size_m) {
  discharge_proj <- st_transform(discharge_sf, crs = 3857)
  buffers <- st_buffer(discharge_proj, dist = buffer_size_m)
  covariate_proj <- st_transform(covariate_sf, crs = 3857)
  joined_data <- st_join(buffers, covariate_proj, join = st_intersects)
  
  aggregated_data <- joined_data %>%
    st_drop_geometry() %>%
    group_by(lon.x, lat.x, year.x) %>%
    summarise(
      !!paste0(covariate_name, "_mean") := ifelse(all(is.na(!!sym(covariate_name))), NA_real_, 
                                                   mean(!!sym(covariate_name), na.rm = TRUE)),
      !!paste0(covariate_name, "_median") := ifelse(all(is.na(!!sym(covariate_name))), NA_real_,
                                                    median(!!sym(covariate_name), na.rm = TRUE)),
      !!paste0(covariate_name, "_sd") := ifelse(all(is.na(!!sym(covariate_name))), NA_real_,
                                                sd(!!sym(covariate_name), na.rm = TRUE)),
      !!paste0(covariate_name, "_count") := sum(!is.na(!!sym(covariate_name))),
      .groups = "drop"
    ) %>%
    rename(lon = lon.x, lat = lat.x, year = year.x) %>%
    filter(!is.na(!!sym(paste0(covariate_name, "_mean"))))
  
  return(aggregated_data)
}

BUFFER_M <- 50000

apr_buffered <- join_and_aggregate_covariates_safe(discharge, apr, "apr", BUFFER_M)
rnf_buffered <- join_and_aggregate_covariates_safe(discharge, rnf, "rnf", BUFFER_M)
erosivity_buffered <- join_and_aggregate_covariates_safe(discharge, ERA5Land, "erosivity", BUFFER_M)

data_merged <- discharge %>%
  left_join(apr_buffered, by = c("lon", "lat", "year")) %>%
  left_join(rnf_buffered, by = c("lon", "lat", "year")) %>%
  left_join(erosivity_buffered, by = c("lon", "lat", "year"))

data_50km_clean <- data_merged %>%
  filter(
    !is.na(log_sum_discharge),
    !is.na(apr_mean),
    !is.na(rnf_mean),
    !is.na(erosivity_mean)
  )

# ============================================================================
# Feature engineering
# Interactions, lagged variables (1-year), and transformations
# ============================================================================

data_engineered <- data_50km_clean %>%
  arrange(lon, lat, year) %>%
  group_by(lon, lat) %>%
  mutate(
    erosivity_apr_interaction = erosivity_mean * apr_mean,
    erosivity_rnf_interaction = erosivity_mean * rnf_mean,
    apr_rnf_interaction = apr_mean * rnf_mean,
    
    apr_lag1 = lag(apr_mean, 1),
    erosivity_lag1 = lag(erosivity_mean, 1),
    discharge_lag1 = lag(log_sum_discharge, 1),
    
    apr_ma3 = zoo::rollmean(apr_mean, k = 3, fill = NA, align = "right"),
    erosivity_ma3 = zoo::rollmean(erosivity_mean, k = 3, fill = NA, align = "right"),
    rnf_ma3 = zoo::rollmean(rnf_mean, k = 3, fill = NA, align = "right"),
    
    cumulative_erosivity = cumsum(erosivity_mean),
    cumulative_apr = cumsum(apr_mean),
    
    erosivity_to_apr_ratio = erosivity_mean / (apr_mean + 0.001),
    normalized_rnf = rnf_mean / (apr_mean + 0.001),
    
    log_apr = log10(apr_mean + 0.001),
    log_erosivity = log10(erosivity_mean + 0.001),
    log_rnf = log10(rnf_mean + 0.001),
    
    erosivity_squared = erosivity_mean^2,
    apr_squared = apr_mean^2,
    
    year_scaled = scale(year)[,1]
  ) %>%
  ungroup()

# Regional stratification
data_engineered <- data_engineered %>%
  mutate(
    climate_zone = case_when(
      abs(lat) > 60 ~ "Polar",
      abs(lat) > 40 ~ "Temperate",
      abs(lat) > 23.5 ~ "Subtropical",
      TRUE ~ "Tropical"
    ),
    hemisphere = if_else(lat >= 0, "Northern", "Southern"),
    coastal_region = case_when(
      lon < -100 & lat > 30 ~ "North America West",
      lon > -100 & lon < -50 & lat > 30 ~ "North America East",
      lon > -20 & lon < 40 & lat > 35 ~ "Europe",
      lon > 100 & lat > 20 ~ "East Asia",
      lon > 100 & lat < 20 & lat > -10 ~ "Southeast Asia",
      TRUE ~ "Other"
    )
  )

# Quality filtering and minimum temporal coverage
data_clean <- data_engineered %>%
  filter(
    !is.na(log_sum_discharge),
    !is.na(apr_mean),
    !is.na(rnf_mean),
    !is.na(erosivity_mean),
    is.finite(log_apr),
    is.finite(log_erosivity),
    is.finite(log_rnf)
  ) %>%
  group_by(lon, lat) %>%
  filter(row_number() > 1) %>%
  ungroup() %>%
  group_by(lon, lat) %>%
  mutate(n_years = n()) %>%
  filter(n_years >= 5) %>%
  ungroup() %>%
  dplyr::select(-n_years)

# ============================================================================
# Feature selection - VIF < 5
# ============================================================================

model_data <- data_clean %>%
  st_drop_geometry() %>%
  filter(complete.cases(.))

# Feature set: erosivity, application rates, runoff, interactions, and 1-year lags
selected_features <- c(
  "erosivity_mean", "apr_mean", "rnf_mean",
  "erosivity_apr_interaction", "erosivity_rnf_interaction",
  "log_erosivity", "log_apr", "log_rnf",
  "apr_lag1", "erosivity_lag1",
  "erosivity_ma3", "apr_ma3"
)

# VIF assessment for multicollinearity (threshold < 5)
feature_data <- model_data[, selected_features]
vif_values <- car::vif(lm(log_sum_discharge ~ ., data = cbind(log_sum_discharge = model_data$log_sum_discharge, feature_data)))

selected_features <- names(vif_values[vif_values < analysis_config$vif_threshold])

# ============================================================================
# Spatial cross-validation (200 Km blocks, 5-fold)
# ============================================================================

model_data_sf <- st_as_sf(model_data, coords = c("lon", "lat"), crs = 4326)
model_data_sf_proj <- st_transform(model_data_sf, crs = 3857)

spatial_folds <- spatialBlock(
  speciesData = model_data_sf_proj,
  species = "log_sum_discharge",
  k = analysis_config$cv_folds,
  selection = "random",
  iteration = 100,
  biomod2Format = FALSE,
  xOffset = 0,
  yOffset = 0,
  r = NULL
)

# Model performance evaluation function
evaluate_model <- function(observed, predicted, model_name) {
  valid_idx <- !is.na(predicted) & !is.na(observed)
  obs <- observed[valid_idx]
  pred <- predicted[valid_idx]
  
  if (length(obs) < 2) {
    return(list(
      Model = model_name, R2 = NA, RMSE = NA, MAE = NA, 
      MAPE = NA, Bias = NA, n = length(obs)
    ))
  }
  
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  mape <- mean(abs((obs - pred) / obs)) * 100
  bias <- mean(pred - obs)
  
  return(list(
    Model = model_name,
    R2 = r2,
    RMSE = rmse,
    MAE = mae,
    MAPE = mape,
    Bias = bias,
    n = length(obs)
  ))
}

# ============================================================================
# Xgboost model
# ============================================================================
param_grid <- expand.grid(
  max_depth = c(4, 6, 8),
  eta = c(0.01, 0.05, 0.1),
  subsample = c(0.7, 0.8, 0.9),
  colsample_bytree = c(0.7, 0.8, 0.9),
  min_child_weight = c(1, 3, 5),
  gamma = c(0, 0.1, 0.2)
)

param_sample <- param_grid[sample(nrow(param_grid), min(20, nrow(param_grid))), ]

evaluate_xgb_params <- function(params, folds, data, features, target_col) {
  cv_scores <- numeric(length(folds))
  
  for (i in 1:length(folds)) {
    fold_indices <- folds[[i]]
    train_idx <- fold_indices[[1]]
    test_idx <- fold_indices[[2]]
    
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    
    dtrain <- xgb.DMatrix(
      data = as.matrix(train_data[, features]),
      label = train_data[[target_col]]
    )
    
    dtest <- xgb.DMatrix(
      data = as.matrix(test_data[, features]),
      label = test_data[[target_col]]
    )
    
    # XGBoost with L1 (alpha) and L2 (lambda) regularization
    xgb_model <- xgb.train(
      params = list(
        objective = "reg:squarederror",
        max_depth = params$max_depth,
        eta = params$eta,
        subsample = params$subsample,
        colsample_bytree = params$colsample_bytree,
        min_child_weight = params$min_child_weight,
        gamma = params$gamma,
        lambda = 1,
        alpha = 0.1
      ),
      data = dtrain,
      nrounds = 500,
      early_stopping_rounds = 50,
      watchlist = list(test = dtest),
      verbose = 0
    )
    
    preds <- predict(xgb_model, dtest)
    cv_scores[i] <- sqrt(mean((test_data[[target_col]] - preds)^2))
  }
  
  return(mean(cv_scores))
}

best_rmse <- Inf
best_params <- NULL

for (i in 1:nrow(param_sample)) {
  params <- param_sample[i, ]
  
  rmse <- evaluate_xgb_params(
    params = params,
    folds = spatial_folds$folds,
    data = model_data,
    features = selected_features,
    target_col = "log_sum_discharge"
  )
  
  if (rmse < best_rmse) {
    best_rmse <- rmse
    best_params <- params
  }
}

xgb_train_preds_full <- rep(NA, nrow(model_data))
xgb_test_preds_full <- rep(NA, nrow(model_data))
xgb_models <- list()

for (i in 1:length(spatial_folds$folds)) {
  fold_indices <- spatial_folds$folds[[i]]
  train_idx <- fold_indices[[1]]
  test_idx <- fold_indices[[2]]
  
  train_data <- model_data[train_idx, ]
  test_data <- model_data[test_idx, ]
  
  dtrain <- xgb.DMatrix(
    data = as.matrix(train_data[, selected_features]),
    label = train_data$log_sum_discharge
  )
  
  dtest <- xgb.DMatrix(
    data = as.matrix(test_data[, selected_features]),
    label = test_data$log_sum_discharge
  )
  
  xgb_fold <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      max_depth = best_params$max_depth,
      eta = best_params$eta,
      subsample = best_params$subsample,
      colsample_bytree = best_params$colsample_bytree,
      min_child_weight = best_params$min_child_weight,
      gamma = best_params$gamma,
      lambda = 1,
      alpha = 0.1,
      eval_metric = "rmse"
    ),
    data = dtrain,
    nrounds = 1000,
    early_stopping_rounds = 50,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )
  
  xgb_models[[i]] <- xgb_fold
  xgb_train_preds_full[train_idx] <- predict(xgb_fold, dtrain)
  xgb_test_preds_full[test_idx] <- predict(xgb_fold, dtest)
}

# Evaluate model performance
xgb_train_metrics <- evaluate_model(
  model_data$log_sum_discharge, 
  xgb_train_preds_full, 
  "XGBoost_Full"
)

xgb_test_metrics <- evaluate_model(
  model_data$log_sum_discharge, 
  xgb_test_preds_full, 
  "XGBoost_Full"
)

# ============================================================================
# Feature importance (Gain-based)
# ============================================================================

importance_list <- lapply(xgb_models, function(model) {
  imp <- xgb.importance(model = model)
  return(imp)
})

all_features <- unique(unlist(lapply(importance_list, function(x) x$Feature)))
avg_importance <- sapply(all_features, function(feat) {
  gains <- sapply(importance_list, function(imp) {
    idx <- which(imp$Feature == feat)
    if (length(idx) > 0) imp$Gain[idx] else 0
  })
  mean(gains)
})

importance_df <- data.frame(
  Feature = all_features,
  Importance = avg_importance
) %>%
  arrange(desc(Importance)) %>%
  mutate(
    Cumulative_Importance = cumsum(Importance) / sum(Importance),
    Feature_Clean = case_when(
      Feature == "erosivity_mean" ~ "Erosivity",
      Feature == "apr_mean" ~ "Application Rate",
      Feature == "rnf_mean" ~ "Runoff",
      Feature == "erosivity_apr_interaction" ~ "Erosivity × Application",
      Feature == "erosivity_rnf_interaction" ~ "Erosivity × Runoff",
      Feature == "log_erosivity" ~ "Log(Erosivity)",
      Feature == "log_apr" ~ "Log(Application)",
      Feature == "log_rnf" ~ "Log(Runoff)",
      TRUE ~ Feature
    )
  )

# ============================================================================
# Spatial autocorrelation assessment
# Moran's I with k-nearest neighbors (k=8)
# ============================================================================

xgb_residuals <- model_data$log_sum_discharge - xgb_test_preds_full

valid_idx <- !is.na(xgb_residuals)
xgb_residuals_clean <- xgb_residuals[valid_idx]
coords_clean <- model_data[valid_idx, c("lon", "lat")]

coords_check <- coords_clean %>%
  group_by(lon, lat) %>%
  summarise(n = n(), .groups = "drop")

if (any(coords_check$n > 1)) {
  residual_by_location <- data.frame(
    lon = coords_clean$lon,
    lat = coords_clean$lat,
    residuals = xgb_residuals_clean
  ) %>%
    group_by(lon, lat) %>%
    summarise(residuals = mean(residuals, na.rm = TRUE), .groups = "drop")
  
  coords_clean <- residual_by_location[, c("lon", "lat")]
  xgb_residuals_clean <- residual_by_location$residuals
}

coords_sf <- st_as_sf(coords_clean, coords = c("lon", "lat"), crs = 4326)
coords_proj <- st_transform(coords_sf, crs = 3857)
coords_mat <- st_coordinates(coords_proj)

if (any(duplicated(coords_mat))) {
  dup_idx <- !duplicated(coords_mat)
  coords_mat <- coords_mat[dup_idx, ]
  xgb_residuals_clean <- xgb_residuals_clean[dup_idx]
}

# Use k-nearest neighbors (k=8) as per methods
knn <- knn2nb(knearneigh(coords_mat, k = analysis_config$moran_k_neighbors))
listw <- nb2listw(knn, style = "W", zero.policy = TRUE)

moran_xgb <- tryCatch({
  moran.test(xgb_residuals_clean, listw, zero.policy = TRUE, alternative = "two.sided")
}, error = function(e) {
  tryCatch({
    spatial_lag <- lag.listw(listw, xgb_residuals_clean, zero.policy = TRUE)
    simple_cor <- cor(xgb_residuals_clean, spatial_lag, use = "complete.obs")
    list(
      statistic = simple_cor,
      p.value = NA,
      method = "Simple spatial correlation (Moran's I approximation)"
    )
  }, error = function(e2) {
    list(statistic = NA, p.value = NA, method = "Failed")
  })
})

if ("statistic" %in% names(moran_xgb)) {
  morans_i_value <- as.numeric(moran_xgb$statistic)
} else if ("observed" %in% names(moran_xgb)) {
  morans_i_value <- as.numeric(moran_xgb$observed)
} else {
  morans_i_value <- NA
}

# ============================================================================
# Final model training
# ============================================================================

dtrain_full <- xgb.DMatrix(
  data = as.matrix(model_data[, selected_features]),
  label = model_data$log_sum_discharge
)

xgb_final <- xgb.train(
  params = list(
    objective = "reg:squarederror",
    max_depth = best_params$max_depth,
    eta = best_params$eta,
    subsample = best_params$subsample,
    colsample_bytree = best_params$colsample_bytree,
    min_child_weight = best_params$min_child_weight,
    gamma = best_params$gamma,
    lambda = 1,
    alpha = 0.1,
    eval_metric = "rmse"
  ),
  data = dtrain_full,
  nrounds = 1000,
  early_stopping_rounds = 50,
  watchlist = list(train = dtrain_full),
  verbose = 0
)

xgb_package <- list(
  model = xgb_final,
  features = selected_features,
  best_params = best_params,
  train_metrics = xgb_train_metrics,
  test_metrics = xgb_test_metrics,
  importance = importance_df
)

# ============================================================================
# Future projections(2050)
# ============================================================================

baseline_summary <- model_data %>%
  group_by(lon, lat) %>%
  summarise(
    baseline_discharge = mean(log_sum_discharge, na.rm = TRUE),
    baseline_erosivity = mean(erosivity_mean, na.rm = TRUE),
    apr_mean = mean(apr_mean, na.rm = TRUE),
    rnf_mean = mean(rnf_mean, na.rm = TRUE),
    apr_lag1 = mean(apr_lag1, na.rm = TRUE),
    erosivity_lag1 = mean(erosivity_lag1, na.rm = TRUE),
    erosivity_ma3 = mean(erosivity_ma3, na.rm = TRUE),
    apr_ma3 = mean(apr_ma3, na.rm = TRUE),
    .groups = "drop"
  )

rcp45_summary <- rcp45_change %>%
  st_drop_geometry() %>%
  rename(erosivity_change_2050_raw = erosivity_change_2050)

baseline_sf <- st_as_sf(baseline_summary, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
rcp45_sf <- st_as_sf(rcp45_summary, coords = c("lon", "lat"), crs = 4326)

baseline_proj <- st_transform(baseline_sf, crs = 3857)
rcp45_proj <- st_transform(rcp45_sf, crs = 3857)

buffers <- st_buffer(baseline_proj, dist = 50000)

joined <- st_join(buffers, rcp45_proj, join = st_intersects)

erosivity_2050_aggregated <- joined %>%
  st_drop_geometry() %>%
  group_by(lon, lat) %>%
  summarise(
    erosivity_change_2050 = mean(erosivity_change_2050_raw, na.rm = TRUE),
    .groups = "drop"
  )

prediction_2050 <- baseline_summary %>%
  left_join(erosivity_2050_aggregated, by = c("lon", "lat")) %>%
  filter(!is.na(erosivity_change_2050))

prediction_2050 <- prediction_2050 %>%
  mutate(
    erosivity_2050 = baseline_erosivity * (1 + erosivity_change_2050 / 100),
    erosivity_change_pct = erosivity_change_2050,
    
    erosivity_mean = erosivity_2050,
    erosivity_apr_interaction = erosivity_2050 * apr_mean,
    erosivity_rnf_interaction = erosivity_2050 * rnf_mean,
    apr_rnf_interaction = apr_mean * rnf_mean,
    
    erosivity_to_apr_ratio = erosivity_2050 / (apr_mean + 0.001),
    normalized_rnf = rnf_mean / (apr_mean + 0.001),
    
    log_apr = log10(apr_mean + 0.001),
    log_erosivity = log10(erosivity_2050 + 0.001),
    log_rnf = log10(rnf_mean + 0.001),
    
    erosivity_squared = erosivity_2050^2,
    apr_squared = apr_mean^2
  )

dpredict_2050 <- xgb.DMatrix(data = as.matrix(prediction_2050[, selected_features]))

prediction_2050$predicted_discharge_2050 <- predict(xgb_package$model, dpredict_2050)

prediction_2050 <- prediction_2050 %>%
  mutate(
    discharge_change = 10^predicted_discharge_2050 - 10^baseline_discharge,
    discharge_change_pct = (10^predicted_discharge_2050 - 10^baseline_discharge) / 10^baseline_discharge * 100,
    
    risk_baseline = case_when(
      10^baseline_discharge < quantile(10^baseline_discharge, 0.33, na.rm = TRUE) ~ "Low",
      10^baseline_discharge < quantile(10^baseline_discharge, 0.67, na.rm = TRUE) ~ "Medium",
      TRUE ~ "High"
    ),
    
    risk_2050 = case_when(
      10^predicted_discharge_2050 < quantile(10^baseline_discharge, 0.33, na.rm = TRUE) ~ "Low",
      10^predicted_discharge_2050 < quantile(10^baseline_discharge, 0.67, na.rm = TRUE) ~ "Medium",
      TRUE ~ "High"
    ),
    
    region = if_else(abs(lat) > 23.5, "Temperate", "Tropical")
  )

summary_stats <- prediction_2050 %>%
  summarise(
    mean_discharge_change_pct = mean(discharge_change_pct, na.rm = TRUE),
    median_discharge_change_pct = median(discharge_change_pct, na.rm = TRUE),
    mean_erosivity_change_pct = mean(erosivity_change_pct, na.rm = TRUE),
    median_erosivity_change_pct = median(erosivity_change_pct, na.rm = TRUE),
    sites_with_increased_risk = sum(discharge_change_pct > 0, na.rm = TRUE),
    total_sites = n(),
    sites_with_increased_risk_pct = (sites_with_increased_risk / total_sites) * 100,
    mean_discharge_change = mean(discharge_change, na.rm = TRUE),
    median_discharge_change = median(discharge_change, na.rm = TRUE),
    .groups = "drop"
  )

transition_matrix <- table(prediction_2050$risk_baseline, prediction_2050$risk_2050)

# ============================================================================
# Monte Carlo uncertainity (1000 simulations)
# ============================================================================

monte_carlo_sim <- function(prediction_data, xgb_model_package, num_sims = 1000) {
  sim_results <- replicate(num_sims, {
    noise_erosivity <- rnorm(nrow(prediction_data), mean = 0, 
                            sd = sd(prediction_data$erosivity_2050, na.rm = TRUE) * 0.1)
    sim_erosivity <- pmax(prediction_data$erosivity_2050 + noise_erosivity, 0)
    
    sim_data <- prediction_data %>%
      mutate(erosivity_mean = sim_erosivity,
             erosivity_apr_interaction = sim_erosivity * apr_mean,
             erosivity_rnf_interaction = sim_erosivity * rnf_mean)
    
    dpredict_sim <- xgb.DMatrix(data = as.matrix(sim_data[, selected_features]))
    sim_prediction <- predict(xgb_model_package$model, dpredict_sim)
    return(sim_prediction)
  })
  
  sim_df <- t(apply(sim_results, 1, function(x) {
    c(mean_pred = mean(x),
      lower_ci = quantile(x, 0.025),
      upper_ci = quantile(x, 0.975))
  }))
  
  return(as.data.frame(sim_df))
}

mc_results <- monte_carlo_sim(prediction_2050, xgb_package, num_sims = 1000)

if (nrow(mc_results) == nrow(prediction_2050)) {
  prediction_2050 <- prediction_2050 %>%
    mutate(
      predicted_discharge_2050_mc_mean = mc_results$mean_pred,
      predicted_discharge_2050_lower_ci = mc_results$lower_ci,
      predicted_discharge_2050_upper_ci = mc_results$upper_ci,
      discharge_change_pct_mc = (10^predicted_discharge_2050_mc_mean - 10^baseline_discharge) / 10^baseline_discharge * 100
    )
}

mc_summary <- prediction_2050 %>%
  summarise(
    mean_discharge_change_pct_mc = mean(discharge_change_pct_mc, na.rm = TRUE),
    sd_discharge_change_pct_mc = sd(discharge_change_pct_mc, na.rm = TRUE),
    median_discharge_change_pct_mc = median(discharge_change_pct_mc, na.rm = TRUE),
    lower_ci_95 = quantile(discharge_change_pct_mc, 0.025, na.rm = TRUE),
    upper_ci_95 = quantile(discharge_change_pct_mc, 0.975, na.rm = TRUE),
    n_observations = n()
  )

# ============================================================================
# Xgboost vs linear regression
# ============================================================================

formula_linear <- as.formula(paste("log_sum_discharge ~", 
                                   paste(selected_features, collapse = " + ")))
linear_final <- lm(formula_linear, data = model_data)

prediction_2050$predicted_discharge_2050_linear <- predict(
  linear_final, 
  newdata = prediction_2050
)

prediction_2050 <- prediction_2050 %>%
  mutate(
    discharge_change_pct_linear = (10^predicted_discharge_2050_linear - 10^baseline_discharge) / 
                                  10^baseline_discharge * 100,
    model_agreement = sign(discharge_change_pct) == sign(discharge_change_pct_linear)
  )

model_comparison_2050 <- data.frame(
  Model = c("XGBoost_Full", "Linear"),
  Mean_Change_Pct = c(
    mean(prediction_2050$discharge_change_pct, na.rm = TRUE),
    mean(prediction_2050$discharge_change_pct_linear, na.rm = TRUE)
  ),
  Median_Change_Pct = c(
    median(prediction_2050$discharge_change_pct, na.rm = TRUE),
    median(prediction_2050$discharge_change_pct_linear, na.rm = TRUE)
  ),
  Sites_Increased_Risk = c(
    sum(prediction_2050$discharge_change_pct > 0, na.rm = TRUE),
    sum(prediction_2050$discharge_change_pct_linear > 0, na.rm = TRUE))
  )
)