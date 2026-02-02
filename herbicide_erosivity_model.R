# ============================================================================
# XGBoost machine learning framework with erosivity-application interactions
# ============================================================================
library(tidyverse)
library(data.table)
library(sf)
library(spdep)
library(blockCV)
library(xgboost)
library(SHAPforxgboost)
library(car)
library(zoo)
library(parallel)

# --- Configuration ---
config <- list(
  spatial_block_size_m = 200000,
  cv_folds = 5,
  temporal_validation_years = 3,
  monte_carlo_sims = 1000,
  vif_threshold = 5,
  min_years_per_location = 5,
  outlier_percentiles = c(0.01, 0.99),
  random_seed = 123,
  rcp45_erosivity_uncertainty_sd = 0.20,
  baseline_variability_factor = 0.15,
  apr_uncertainty_sd = 0.15,
  rnf_uncertainty_sd = 0.15
)

set.seed(config$random_seed)

# ============================================================================
# Data loading and quality Control
# ============================================================================

discharge_raw <- read_csv("DIS_annual.csv", 
                         show_col_types = FALSE)

discharge <- discharge_raw %>%
  filter(
    !is.na(sum_discharge),
    sum_discharge > 0,
    sum_discharge >= quantile(sum_discharge, config$outlier_percentiles[1], na.rm = TRUE),
    sum_discharge <= quantile(sum_discharge, config$outlier_percentiles[2], na.rm = TRUE)
  ) %>%
  mutate(
    log_discharge = log10(sum_discharge),
    location_id = paste(round(lon, 2), round(lat, 2), sep = "_")
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# --- Load covariates ---
apr <- fread("annual_apr.csv") %>%
  rename(apr = apr) %>%
  filter(year %in% unique(discharge$year)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rnf <- fread("annual_rnf.csv") %>%
  rename(rnf = rnf) %>%
  filter(year %in% unique(discharge$year)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

erosivity <- fread("ERA5Land_2001_2021_erosivity.csv") %>%
  rename(erosivity = data) %>%
  filter(year %in% unique(discharge$year)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

rcp45_change <- fread("2010-2050-rcp4.5_final.csv") %>%
  rename(erosivity_change_2050 = data) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# ============================================================================
# Spatial Aggregation (50 km Buffer)
# ============================================================================

join_aggregate_covariate <- function(discharge_sf, covariate_sf, 
                                    covariate_name, buffer_m = 50000) {
  discharge_proj <- st_transform(discharge_sf, crs = 3857)
  buffers <- st_buffer(discharge_proj, dist = buffer_m)
  covariate_proj <- st_transform(covariate_sf, crs = 3857)
  
  joined <- st_join(buffers, covariate_proj, join = st_intersects)
  
  aggregated <- joined %>%
    st_drop_geometry() %>%
    group_by(lon.x, lat.x, year.x, location_id) %>%
    summarise(
      !!paste0(covariate_name, "_mean") := mean(!!sym(covariate_name), na.rm = TRUE),
      !!paste0(covariate_name, "_sd") := sd(!!sym(covariate_name), na.rm = TRUE),
      !!paste0(covariate_name, "_count") := sum(!is.na(!!sym(covariate_name))),
      .groups = "drop"
    ) %>%
    rename(lon = lon.x, lat = lat.x, year = year.x) %>%
    filter(!is.na(!!sym(paste0(covariate_name, "_mean"))))
  
  return(aggregated)
}

apr_agg <- join_aggregate_covariate(discharge, apr, "apr")
rnf_agg <- join_aggregate_covariate(discharge, rnf, "rnf")
erosivity_agg <- join_aggregate_covariate(discharge, erosivity, "erosivity")

data_merged <- discharge %>%
  left_join(apr_agg, by = c("lon", "lat", "year", "location_id")) %>%
  left_join(rnf_agg, by = c("lon", "lat", "year", "location_id")) %>%
  left_join(erosivity_agg, by = c("lon", "lat", "year", "location_id")) %>%
  filter(
    !is.na(log_discharge),
    !is.na(apr_mean),
    !is.na(rnf_mean),
    !is.na(erosivity_mean)
  )

# ============================================================================
# Feature engineering
# ============================================================================

data_features <- data_merged %>%
  arrange(location_id, year) %>%
  group_by(location_id) %>%
  mutate(
    apr_lag1 = lag(apr_mean, 1),
    erosivity_lag1 = lag(erosivity_mean, 1),
    
    apr_ma3 = zoo::rollmean(apr_mean, k = 3, fill = NA, align = "right"),
    erosivity_ma3 = zoo::rollmean(erosivity_mean, k = 3, fill = NA, align = "right"),
    
    erosivity_apr_interaction = erosivity_mean * apr_mean,
    erosivity_rnf_interaction = erosivity_mean * rnf_mean,
    
    log_apr = log10(apr_mean + 0.001),
    log_erosivity = log10(erosivity_mean + 0.001),
    log_rnf = log10(rnf_mean + 0.001),
    
    erosivity_squared = erosivity_mean^2,
    apr_squared = apr_mean^2,
    
    erosivity_to_apr_ratio = erosivity_mean / (apr_mean + 0.001),
    
    year_scaled = scale(year)[,1],
    
    n_years = n()
  ) %>%
  ungroup() %>%
  mutate(
    climate_zone = case_when(
      abs(lat) > 60 ~ "Polar",
      abs(lat) > 40 ~ "Temperate",
      abs(lat) > 23.5 ~ "Subtropical",
      TRUE ~ "Tropical"
    ),
    region = if_else(abs(lat) > 23.5, "Temperate", "Tropical")
  )

data_clean <- data_features %>%
  filter(
    n_years >= config$min_years_per_location,
    !is.na(apr_lag1),
    is.finite(log_apr),
    is.finite(log_erosivity),
    is.finite(log_rnf)
  ) %>%
  dplyr::select(-n_years)

# ============================================================================
# Spatiotemporal cross-validation (200 km spatial blocks, 5-fold)
# ============================================================================

years_all <- sort(unique(data_clean$year))
years_train <- years_all[1:(length(years_all) - config$temporal_validation_years)]
years_temporal_test <- years_all[(length(years_all) - config$temporal_validation_years + 1):length(years_all)]

data_train <- data_clean %>% filter(year %in% years_train)
data_temporal_test <- data_clean %>% filter(year %in% years_temporal_test)

data_train_sf <- st_as_sf(data_train, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 3857)

spatial_blocks <- spatialBlock(
  speciesData = data_train_sf,
  species = "log_discharge",
  theRange = config$spatial_block_size_m,
  k = config$cv_folds,
  selection = "random",
  iteration = 100,
  biomod2Format = FALSE,
  verbose = FALSE
)

# ============================================================================
# Feature selection - VIF < 5
# ============================================================================

candidate_features <- c(
  "erosivity_mean", "apr_mean", "rnf_mean",
  "erosivity_apr_interaction", "erosivity_rnf_interaction",
  "apr_lag1", "erosivity_lag1",
  "erosivity_ma3", "apr_ma3",
  "log_erosivity", "log_apr", "log_rnf",
  "erosivity_squared", "apr_squared",
  "year_scaled"
)

# VIF-based feature selection
X_vif <- data_train %>%
  st_drop_geometry() %>%
  dplyr::select(all_of(candidate_features)) %>%
  as.data.frame()

vif_values <- car::vif(lm(log_discharge ~ ., 
                         data = cbind(log_discharge = data_train$log_discharge, X_vif)))

final_features <- names(vif_values[vif_values < config$vif_threshold])

# ============================================================================
# XGBoost hyperparameter tuning with L1/L2 regularization
# ============================================================================

param_grid <- expand.grid(
  max_depth = c(4, 6, 8),
  eta = c(0.01, 0.05, 0.1),
  subsample = c(0.7, 0.8, 0.9),
  colsample_bytree = c(0.7, 0.8, 0.9),
  min_child_weight = c(1, 3, 5),
  gamma = c(0, 0.1, 0.2)
)

evaluate_xgb_params <- function(params, folds, data, features, target) {
  cv_scores <- numeric(length(folds))
  
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]][[1]]
    test_idx <- folds[[i]][[2]]
    
    dtrain <- xgb.DMatrix(
      data = as.matrix(data[train_idx, features, drop = FALSE]),
      label = data[[target]][train_idx]
    )
    dtest <- xgb.DMatrix(
      data = as.matrix(data[test_idx, features, drop = FALSE]),
      label = data[[target]][test_idx]
    )
    
    model <- xgb.train(
      params = c(
        objective = "reg:squarederror",
        as.list(params),
        lambda = 1,
        alpha = 0.1,
        eval_metric = "rmse"
      ),
      data = dtrain,
      nrounds = 50,
      early_stopping_rounds = 10,
      watchlist = list(val = dtest),
      verbose = 0
    )
    
    preds <- predict(model, dtest)
    cv_scores[i] <- sqrt(mean((data[[target]][test_idx] - preds)^2))
  }
  
  return(mean(cv_scores))
}

# Parallel hyperparameter search
n_cores <- min(detectCores() - 1, 12)
data_train_df <- data_train %>% st_drop_geometry()

cl <- makeCluster(n_cores)
clusterExport(cl, c("evaluate_xgb_params", "spatial_blocks", "data_train_df", 
                   "final_features", "param_grid"))
clusterEvalQ(cl, {
  library(xgboost)
  library(dplyr)
})

param_results <- parLapply(cl, 1:nrow(param_grid), function(i) {
  params <- as.list(param_grid[i, ])
  rmse <- evaluate_xgb_params(
    params = params,
    folds = spatial_blocks$folds,
    data = data_train_df,
    features = final_features,
    target = "log_discharge"
  )
  return(list(params = params, rmse = rmse))
})

stopCluster(cl)

rmse_values <- sapply(param_results, function(x) x$rmse)
best_idx <- which.min(rmse_values)
best_params <- param_results[[best_idx]]$params

# ============================================================================
# Train XGBoost models on each spatial fold
# ============================================================================

data_train_clean <- data_train %>%
  st_drop_geometry() %>%
  dplyr::select(all_of(c("log_discharge", final_features))) %>%
  as.data.frame()

xgb_models <- list()
xgb_predictions_train <- numeric(nrow(data_train_clean))
xgb_predictions_test <- numeric(nrow(data_train_clean))

for(i in 1:length(spatial_blocks$folds)) {
  train_idx <- spatial_blocks$folds[[i]][[1]]
  test_idx <- spatial_blocks$folds[[i]][[2]]
  
  dtrain <- xgb.DMatrix(
    data = as.matrix(data_train_clean[train_idx, final_features, drop = FALSE]),
    label = data_train_clean$log_discharge[train_idx]
  )
  dtest <- xgb.DMatrix(
    data = as.matrix(data_train_clean[test_idx, final_features, drop = FALSE]),
    label = data_train_clean$log_discharge[test_idx]
  )
  
  model <- xgb.train(
    params = c(
      objective = "reg:squarederror",
      as.list(best_params),
      lambda = 1,
      alpha = 0.1,
      eval_metric = "rmse"
    ),
    data = dtrain,
    nrounds = 200,
    early_stopping_rounds = 20,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )
  
  xgb_models[[i]] <- model
  xgb_predictions_train[train_idx] <- predict(model, dtrain)
  xgb_predictions_test[test_idx] <- predict(model, dtest)
}

# --- Performance metrics ---
calculate_metrics <- function(actual, predicted) {
  valid_idx <- !is.na(actual) & !is.na(predicted)
  actual <- actual[valid_idx]
  predicted <- predicted[valid_idx]
  
  list(
    R2 = cor(actual, predicted)^2,
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAE = mean(abs(actual - predicted)),
    Bias = mean(predicted - actual)
  )
}

metrics_train_cv <- calculate_metrics(data_train_clean$log_discharge, xgb_predictions_train)
metrics_test_cv <- calculate_metrics(data_train_clean$log_discharge, xgb_predictions_test)

cv_performance <- data.frame(
  Dataset = c("Train", "Test"),
  R2 = c(metrics_train_cv$R2, metrics_test_cv$R2),
  RMSE = c(metrics_train_cv$RMSE, metrics_test_cv$RMSE),
  MAE = c(metrics_train_cv$MAE, metrics_test_cv$MAE),
  Bias = c(metrics_train_cv$Bias, metrics_test_cv$Bias),
  Overfitting_Ratio = c(NA, metrics_test_cv$R2 / metrics_train_cv$R2)
)

# ============================================================================
# Train final model on all training Data
# ============================================================================

final_nrounds <- round(median(sapply(xgb_models, function(m) m$best_iteration)))

dtrain_full <- xgb.DMatrix(
  data  = as.matrix(data_train_clean[, final_features, drop = FALSE]),
  label = data_train_clean$log_discharge
)

xgb_final <- xgb.train(
  params = c(
    objective = "reg:squarederror",
    as.list(best_params),
    lambda = 1,
    alpha  = 0.1,
    eval_metric = "rmse"
  ),
  data    = dtrain_full,
  nrounds = final_nrounds,
  verbose = 0
)

# ============================================================================
# Temporal validation
# ============================================================================

dtest_temporal <- xgb.DMatrix(
  data = as.matrix(data_temporal_test %>% st_drop_geometry() %>% 
                  dplyr::select(all_of(final_features)))
)

preds_temporal <- predict(xgb_final, dtest_temporal)
metrics_temporal <- calculate_metrics(data_temporal_test$log_discharge, preds_temporal)

temporal_performance <- data.frame(
  R2 = metrics_temporal$R2,
  RMSE = metrics_temporal$RMSE,
  MAE = metrics_temporal$MAE,
  Bias = metrics_temporal$Bias
)

# ============================================================================
# Spatial autocorrelation assessment (Moran's I, k=8)
# ============================================================================

residuals_cv <- data_train_clean$log_discharge - xgb_predictions_test

residual_spatial <- data_train %>%
  st_drop_geometry() %>%
  mutate(residual = residuals_cv) %>%
  group_by(location_id, lon, lat) %>%
  summarise(
    mean_residual = mean(residual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_residual))

coords <- as.matrix(data.frame(
  lon = as.numeric(residual_spatial$lon),
  lat = as.numeric(residual_spatial$lat)
))

valid_coords <- complete.cases(coords)
coords <- coords[valid_coords, ]
mean_residuals <- residual_spatial$mean_residual[valid_coords]

knn <- knearneigh(coords, k = 8)
nb <- knn2nb(knn)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_result <- tryCatch({
  moran.mc(mean_residuals, listw, nsim = 999, zero.policy = TRUE)
}, error = function(e) {
  spatial_lag <- lag.listw(listw, mean_residuals, zero.policy = TRUE)
  list(
    statistic = cor(mean_residuals, spatial_lag, use = "complete.obs"),
    p.value = NA
  )
})

spatial_autocorr_results <- data.frame(
  Metric = "Moran's I",
  Value = moran_result$statistic,
  P_Value = if(!is.na(moran_result$p.value)) moran_result$p.value else NA,
  N_Locations = nrow(coords)
)

# ============================================================================
# Feature importance (SHAP Values)
# ============================================================================

gain_list <- lapply(xgb_models, function(m) {
  imp <- xgb.importance(model = m)
  imp$Gain
})

gain_df <- data.frame(
  Feature = xgb.importance(model = xgb_models[[1]])$Feature,
  Gain_Mean = rowMeans(do.call(cbind, gain_list)),
  Gain_SD   = apply(do.call(cbind, gain_list), 1, sd)
)

X_shap <- as.matrix(data_train_clean[, final_features, drop = FALSE])

shap_vals <- shap.values(
  xgb_model = xgb_final,
  X_train   = X_shap
)

shap_long <- shap.prep(
  shap_contrib = shap_vals$shap_score,
  X_train      = X_shap
)

shap_df <- shap_long %>%
  group_by(variable) %>%
  summarise(
    SHAP_MeanAbs = mean(abs(value)),
    SHAP_Mean    = mean(value),
    SHAP_PosFrac = mean(value > 0),
    .groups = "drop"
  )

importance_final <- gain_df %>%
  left_join(shap_df, by = c("Feature" = "variable")) %>%
  arrange(desc(SHAP_MeanAbs))

# ============================================================================
# 2050 Projections (RCP4.5)
# ============================================================================

baseline_years <- years_all[years_all >= 2001]

baseline_data <- data_clean %>%
  filter(year %in% baseline_years) %>%
  group_by(location_id, lon, lat) %>%
  summarise(
    baseline_apr        = mean(apr_mean, na.rm = TRUE),
    baseline_rnf        = mean(rnf_mean, na.rm = TRUE),
    baseline_erosivity  = mean(erosivity_mean, na.rm = TRUE),
    baseline_apr_lag1   = mean(apr_lag1, na.rm = TRUE),
    baseline_erosivity_lag1 = mean(erosivity_lag1, na.rm = TRUE),
    region        = first(region),
    climate_zone  = first(climate_zone),
    .groups = "drop"
  )

coastal_sf <- baseline_data %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

nearest_idx <- st_nearest_feature(
  st_transform(coastal_sf, 3857),
  st_transform(rcp45_change, 3857)
)

rcp45_matched <- rcp45_change[nearest_idx, ] %>%
  st_drop_geometry()

prediction_2050 <- baseline_data %>%
  mutate(
    erosivity_mean = pmax(baseline_erosivity + rcp45_matched$erosivity_change_2050, 0.1),
    apr_mean = baseline_apr,
    rnf_mean = baseline_rnf,
    apr_lag1 = baseline_apr_lag1,
    erosivity_lag1 = baseline_erosivity_lag1,
    log_apr = log10(apr_mean + 0.001),
    log_erosivity = log10(erosivity_mean + 0.001),
    log_rnf = log10(rnf_mean + 0.001),
    erosivity_apr_interaction = erosivity_mean * apr_mean,
    year_scaled = 0
  )

pred_check <- prediction_2050 %>%
  st_drop_geometry() %>%
  as.data.frame()

dpredict_2050 <- xgb.DMatrix(
  data = as.matrix(pred_check[, final_features, drop = FALSE])
)

prediction_2050$predicted_log_discharge_2050 <-
  predict(xgb_final, dpredict_2050)

baseline_pred_data <- baseline_data %>%
  mutate(
    apr_mean = baseline_apr,
    rnf_mean = baseline_rnf,
    erosivity_mean = baseline_erosivity,
    apr_lag1 = baseline_apr_lag1,
    erosivity_lag1 = baseline_erosivity_lag1,
    log_apr = log10(apr_mean + 0.001),
    log_rnf = log10(rnf_mean + 0.001),
    log_erosivity = log10(erosivity_mean + 0.001),
    erosivity_apr_interaction = erosivity_mean * apr_mean,
    year_scaled = 0
  ) %>%
  st_drop_geometry() %>%
  as.data.frame()

dpredict_baseline <- xgb.DMatrix(
  data = as.matrix(baseline_pred_data[, final_features, drop = FALSE])
)

prediction_2050$predicted_log_discharge_baseline <-
  predict(xgb_final, dpredict_baseline)

prediction_2050 <- prediction_2050 %>%
  mutate(
    discharge_change_log =
      predicted_log_discharge_2050 - predicted_log_discharge_baseline,
    discharge_change_pct =
      (10^discharge_change_log - 1) * 100
  )

# ============================================================================
# Monte Carlo uncertainty quantification (1000 iterations)
# ============================================================================

monte_carlo_robust <- function(prediction_data_clean, cv_models,
                               config, num_sims = 1000) {
  n_sites <- nrow(prediction_data_clean)
  sim_results <- matrix(NA, nrow = n_sites, ncol = num_sims)
  
  for (i in seq_len(num_sims)) {
    noise_rcp <- rnorm(n_sites, 1, config$rcp45_erosivity_uncertainty_sd)
    noise_apr <- rnorm(n_sites, 1, config$apr_uncertainty_sd)
    noise_rnf <- rnorm(n_sites, 1, config$rnf_uncertainty_sd)
    
    sim_data <- prediction_data_clean %>%
      mutate(
        apr_sim = pmax(baseline_apr * noise_apr, 0.01),
        rnf_sim = pmax(baseline_rnf * noise_rnf, 0.01),
        erosivity_sim = pmax(baseline_erosivity * noise_rcp, 0.1),
        log_apr = log10(apr_sim + 0.001),
        log_rnf = log10(rnf_sim + 0.001),
        log_erosivity = log10(erosivity_sim + 0.001),
        erosivity_apr_interaction = erosivity_sim * apr_sim
      )
    
    model <- cv_models[[sample(seq_along(cv_models), 1)]]
    
    dpredict_sim <- xgb.DMatrix(
      data = as.matrix(sim_data[, final_features, drop = FALSE])
    )
    
    sim_results[, i] <- predict(model, dpredict_sim)
  }
  
  uncertainty <- data.frame(
    mean_pred = rowMeans(sim_results),
    lower_ci = apply(sim_results, 1, quantile, 0.025),
    upper_ci = apply(sim_results, 1, quantile, 0.975),
    sd_pred = apply(sim_results, 1, sd)
  )
  
  return(uncertainty)
}

mc_uncertainty <- monte_carlo_robust(
  prediction_data_clean = prediction_2050 %>% st_drop_geometry(),
  cv_models = xgb_models,
  config = config,
  num_sims = config$monte_carlo_sims
)

prediction_2050 <- prediction_2050 %>%
  bind_cols(mc_uncertainty) %>%
  mutate(
    discharge_change_pct_lower = (10^(lower_ci - predicted_log_discharge_baseline) - 1) * 100,
    discharge_change_pct_upper = (10^(upper_ci - predicted_log_discharge_baseline) - 1) * 100
  )

# --- Regional summary ---
regional_summary <- prediction_2050 %>%
  st_drop_geometry() %>%
  group_by(region) %>%
  summarise(
    n_sites = n(),
    mean_change_pct = mean(discharge_change_pct, na.rm = TRUE),
    median_change_pct = median(discharge_change_pct, na.rm = TRUE),
    sites_increasing = sum(discharge_change_pct > 0, na.rm = TRUE),
    pct_increasing = (sites_increasing / n()) * 100,
    ci_lower = mean(discharge_change_pct_lower, na.rm = TRUE),
    ci_upper = mean(discharge_change_pct_upper, na.rm = TRUE),
    .groups = "drop"
  )

# --- Global summary ---
global_summary <- prediction_2050 %>%
  st_drop_geometry() %>%
  summarise(
    n_sites = n(),
    mean_change_pct = mean(discharge_change_pct, na.rm = TRUE),
    median_change_pct = median(discharge_change_pct, na.rm = TRUE),
    sites_increasing = sum(discharge_change_pct > 0, na.rm = TRUE),
    pct_increasing = (sites_increasing / n()) * 100,
    ci_lower = mean(discharge_change_pct_lower, na.rm = TRUE),
    ci_upper = mean(discharge_change_pct_upper, na.rm = TRUE)
  )

# ============================================================================
# Linear model comparison
# ============================================================================

linear_final <- lm(
  log_discharge ~ .,
  data = data_train_clean[, c("log_discharge", final_features)]
)

prediction_2050$predicted_discharge_2050_linear <- predict(
  linear_final,
  newdata = pred_check
)

prediction_2050$predicted_discharge_baseline_linear <- predict(
  linear_final,
  newdata = baseline_pred_data
)

prediction_2050 <- prediction_2050 %>%
  mutate(
    discharge_change_pct_linear = 
      (10^(predicted_discharge_2050_linear - predicted_discharge_baseline_linear) - 1) * 100,
    model_agreement = sign(discharge_change_pct) == sign(discharge_change_pct_linear)
  )

model_comparison <- data.frame(
  Model = c("XGBoost", "Linear"),
  Mean_Change_Pct = c(
    mean(prediction_2050$discharge_change_pct, na.rm = TRUE),
    mean(prediction_2050$discharge_change_pct_linear, na.rm = TRUE)
  ),
  Median_Change_Pct = c(
    median(prediction_2050$discharge_change_pct, na.rm = TRUE),
    median(prediction_2050$discharge_change_pct_linear, na.rm = TRUE)
  ),
  Sites_Increased = c(
    sum(prediction_2050$discharge_change_pct > 0, na.rm = TRUE),
    sum(prediction_2050$discharge_change_pct_linear > 0, na.rm = TRUE)
  ),
  Directional_Agreement_Pct = c(
    NA,
    mean(prediction_2050$model_agreement, na.rm = TRUE) * 100
  )
)
