# ============================================================================
# Phytoplankton-herbicide-SST interaction modeling
# ============================================================================
rm(list = ls())
set.seed(123)

library(tidyverse)
library(data.table)
library(sf)
library(spdep)
library(blockCV)
library(mgcv)
library(car)
library(FNN)
library(parallel)
library(future.apply)

plan(multisession, workers = 30)

# --- Configuration ---
config <- list(
  common_years = 2003:2015,
  polar_threshold = 60,
  buffer_km = 50,
  k_spatial = 100,
  cv_folds = 5,
  spatial_block_size_m = 200000,
  moran_k = 8,
  monte_carlo_sims = 1000,
  sst_increase_2050 = 1.5,
  random_seed = 123,
  vif_threshold = 5
)

set.seed(config$random_seed)

# ============================================================================
# Data loading and preprocessing
# ============================================================================

# --- Herbicide discharge ---
discharge_raw <- read_csv("DIS_annual.csv", show_col_types = FALSE) %>%
  mutate(
    log_sum_discharge = if_else(sum_discharge == 0 | is.na(sum_discharge), 
                                log10(0.001), log10(sum_discharge)),
    log_mean_discharge = if_else(mean_discharge == 0 | is.na(mean_discharge), 
                                 log10(0.001), log10(mean_discharge))
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# --- Phytoplankton groups ---
phyto_groups <- c(
  "chlorophyll_a" = "tmp/chlorophyll_a_coastal_chlorophyll_annual.csv",
  "diatoms" = "diatoms_coastal_chlorophyll_annual.csv",
  "dinoflagellates" = "dinoflagellates_coastal_chlorophyll_annual.csv",
  "green_algae" = "green_algae_coastal_chlorophyll_annual.csv",
  "haptophytes" = "haptophytes_coastal_chlorophyll_annual.csv",
  "prochlorococcus" = "prochlorococcus_coastal_chlorophyll_annual.csv",
  "microphytoplankton" = "microphytoplankton_coastal_chlorophyll_annual.csv",
  "nanophytoplankton" = "nanophytoplankton_coastal_chlorophyll_annual.csv",
  "picophytoplankton" = "picophytoplankton_coastal_chlorophyll_annual.csv"
)

chlorophyll_raw_list <- list()
for (group_name in names(phyto_groups)) {
  file_path <- phyto_groups[[group_name]]
  if (file.exists(file_path)) {
    chlo_raw <- read_csv(file_path, show_col_types = FALSE) %>%
      rename(chlorophyll = variable) %>%
      filter(year %in% config$common_years) %>%
      mutate(chlorophyll = ifelse(chlorophyll <= 0, NA, chlorophyll)) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
      mutate(phyto_group = group_name)
    chlorophyll_raw_list[[group_name]] <- chlo_raw
  }
}

chlorophyll_raw_combined <- bind_rows(chlorophyll_raw_list)

# --- SST data ---
sst_stack <- terra::rast("sst_stack.tif")
satellite_years_sst <- 2003:2020
nodata_threshold <- -3e38
sst_stack[sst_stack < nodata_threshold] <- NA

common_sst_years_idx <- which(satellite_years_sst %in% config$common_years)
sst_stack_common <- sst_stack[[common_sst_years_idx]]
names(sst_stack_common) <- paste0("sst_", config$common_years[config$common_years %in% satellite_years_sst])

# --- Herbicide budget for conversion ---
river_budget <- readxl::read_excel("tmp/data/Table_PAS92_BUDGET.xlsx", sheet = "RIVER budget") %>%
  filter(`Type (main)` == "HBC") %>%
  dplyr::select(`PAS`, `Discharge avg [%]`) %>%
  rename(discharge_pct = `Discharge avg [%]`)

hbc_discharge_pct <- river_budget %>%
  summarise(avg_discharge_pct = mean(discharge_pct, na.rm = TRUE)) %>%
  pull(avg_discharge_pct)


# Nutrient data

nutrients <- list()

# Nitrate
if (file.exists("nitrate_merged.csv")) {
  nutrients$nitrate <- fread("tmp/nutrients/nitrate_merged.csv") %>%
    filter(year %in% config$common_years) %>%
    rename(nitrate = variable) %>%
    mutate(
      lon = round(lon, 4),
      lat = round(lat, 4)
    )
}

# Phosphate
if (file.exists("phosphate_merged.csv")) {
  nutrients$phosphate <- fread("tmp/nutrients/phosphate_merged.csv") %>%
    filter(year %in% config$common_years) %>%
    rename(phosphate = variable) %>%
    mutate(
      lon = round(lon, 4),
      lat = round(lat, 4)
    )
}

# Silicate
if (file.exists("silicate_merged.csv")) {
  nutrients$silicate <- fread("tmp/nutrients/silicate_merged.csv") %>%
    filter(year %in% config$common_years) %>%
    rename(silicate = variable) %>%
    mutate(
      lon = round(lon, 4),
      lat = round(lat, 4)
    )
}

# ============================================================================
# Nutrient data availability check
# ============================================================================

if (length(nutrients) > 0) {
  # Merge all nutrient data
  nutrient_combined <- nutrients$nitrate %>%
    full_join(nutrients$phosphate, by = c("lon", "lat", "year")) %>%
    full_join(nutrients$silicate, by = c("lon", "lat", "year")) %>%
    filter(!is.na(nitrate) & !is.na(phosphate) & !is.na(silicate))
  
  # For each phytoplankton group, check nutrient overlap
  nutrient_coverage_summary <- list()
  
  for (group_name in names(model_data_list)) {
    group_data <- model_data_list[[group_name]]
    
    # Find locations with both phytoplankton and complete nutrient data
    group_with_nutrients <- group_data %>%
      inner_join(nutrient_combined, by = c("lon", "lat", "year"))
    
    total_locations <- group_data %>%
      distinct(lon, lat) %>%
      nrow()
    
    locations_with_nutrients <- group_with_nutrients %>%
      distinct(lon, lat) %>%
      nrow()
    
    pct_coverage <- (locations_with_nutrients / total_locations) * 100
    
    nutrient_coverage_summary[[group_name]] <- data.frame(
      phyto_group = group_name,
      total_observations = nrow(group_data),
      total_locations = total_locations,
      nutrient_observations = nrow(group_with_nutrients),
      nutrient_locations = locations_with_nutrients,
      percent_coverage = pct_coverage
    )
  }
  
  nutrient_coverage_df <- bind_rows(nutrient_coverage_summary)
  
  # Store for later use
  model_data_with_nutrients <- list()
  for (group_name in names(model_data_list)) {
    group_data <- model_data_list[[group_name]]
    group_with_nutrients <- group_data %>%
      inner_join(nutrient_combined, by = c("lon", "lat", "year"))
    
    if (nrow(group_with_nutrients) >= 30) {
      model_data_with_nutrients[[group_name]] <- group_with_nutrients
    }
  }
} else {
  nutrient_coverage_df <- NULL
  model_data_with_nutrients <- NULL
}

# --- 2050 projections from erosivity model ---
if(file.exists("herb_erosivity_log.RData")) {
  load("herb_erosivity_log.RData")
  herbicide_2050_predictions <- prediction_2050_sf %>%
    mutate(
      lon = st_coordinates(.)[,1],
      lat = st_coordinates(.)[,2]
    ) %>%
    dplyr::select(lon, lat, predicted_discharge_2050, discharge_change_pct) %>%
    st_drop_geometry()
}

# ============================================================================
# Spatial filtering and aggregation (50 km buffer)
# ============================================================================

filter_polar <- function(data_sf) {
  data_sf %>% filter(abs(lat) < config$polar_threshold)
}

discharge_filtered <- filter_polar(discharge_raw) %>%
  filter(year %in% config$common_years)

# --- Spatial aggregation function ---
aggregate_to_discharge_buffer <- function(chlo_data, discharge_data, buffer_km) {
  discharge_proj <- st_transform(discharge_data, crs = 3857)
  buffers <- st_buffer(discharge_proj, dist = buffer_km * 1000)
  
  chlo_proj <- st_transform(chlo_data, crs = 3857)
  
  joined <- st_join(buffers, chlo_proj, join = st_intersects)
  
  aggregated <- joined %>%
    st_drop_geometry() %>%
    filter(!is.na(chlorophyll)) %>%
    group_by(lon.x, lat.x, year.x, phyto_group) %>%
    summarise(
      chlorophyll_mean = mean(chlorophyll, na.rm = TRUE),
      chlorophyll_sd = sd(chlorophyll, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    rename(lon = lon.x, lat = lat.x, year = year.x)
  
  return(aggregated)
}

# --- Apply spatial aggregation for each group ---
chlorophyll_aggregated_list <- list()
for (group_name in names(chlorophyll_raw_list)) {
  chlo_agg <- aggregate_to_discharge_buffer(
    chlorophyll_raw_list[[group_name]],
    discharge_filtered,
    config$buffer_km
  )
  chlorophyll_aggregated_list[[group_name]] <- chlo_agg
}

# ============================================================================
# Extract SST and merge datasets
# ============================================================================

extract_sst_for_group <- function(aggregated_data, sst_stack_layers) {
  coords_sf <- st_as_sf(aggregated_data, coords = c("lon", "lat"), crs = 4326)
  
  sst_values <- terra::extract(sst_stack_layers, coords_sf, ID = FALSE)
  
  aggregated_with_sst <- bind_cols(aggregated_data, sst_values)
  
  return(aggregated_with_sst)
}

# --- Merge discharge, SST, and chlorophyll ---
merge_all_data <- function(discharge_data, chlorophyll_agg, sst_extracted, hbc_pct) {
  discharge_clean <- discharge_data %>%
    st_drop_geometry() %>%
    dplyr::select(lon, lat, year, sum_discharge, log_sum_discharge)
  
  merged <- chlorophyll_agg %>%
    left_join(discharge_clean, by = c("lon", "lat", "year")) %>%
    mutate(
      herbicide_discharge = sum_discharge * (hbc_pct / 100),
      log_herbicide_discharge = log10(herbicide_discharge + 0.001)
    )
  
  sst_year_col <- paste0("sst_", merged$year)
  merged$sst <- NA_real_
  
  for (i in 1:nrow(merged)) {
    year_col <- paste0("sst_", merged$year[i])
    if (year_col %in% names(sst_extracted)) {
      merged$sst[i] <- sst_extracted[[year_col]][i]
    }
  }
  
  merged <- merged %>%
    filter(!is.na(sst), !is.na(chlorophyll_mean), !is.na(log_herbicide_discharge)) %>%
    mutate(
      log_chlorophyll = log10(chlorophyll_mean),
      lat_band = cut(lat, breaks = seq(-60, 60, by = 10), labels = FALSE)
    )
  
  return(merged)
}

# --- Process all groups ---
model_data_list <- list()
for (group_name in names(chlorophyll_aggregated_list)) {
  chlo_with_sst <- extract_sst_for_group(
    chlorophyll_aggregated_list[[group_name]],
    sst_stack_common
  )
  
  model_data <- merge_all_data(
    discharge_filtered,
    chlo_with_sst,
    chlo_with_sst,
    hbc_discharge_pct
  )
  
  model_data_list[[group_name]] <- model_data
}

# ============================================================================
# Spatial cross-validation setup (200 km blocks, 5-fold)
# ============================================================================

create_spatial_folds <- function(data, k_folds, block_size_m) {
  data_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326)
  data_proj <- st_transform(data_sf, crs = 3857)
  
  spatial_folds <- spatialBlock(
    speciesData = data_proj,
    species = "log_chlorophyll",
    k = k_folds,
    selection = "random",
    iteration = 100,
    biomod2Format = FALSE,
    xOffset = 0,
    yOffset = 0,
    r = NULL
  )
  
  return(spatial_folds)
}

# ============================================================================
# GAM fodel fitting
# ============================================================================

fit_gam_models <- function(data, k_spatial = 100) {
  data <- data %>%
    mutate(
      lon_scaled = scale(lon)[,1],
      lat_scaled = scale(lat)[,1]
    )
  
  # M1: Herbicide + Spatial
  formula_herbicide <- as.formula(
    paste0("log_chlorophyll ~ s(log_herbicide_discharge, k=10, bs='tp') + ",
           "lat_band + s(lon_scaled, lat_scaled, k=", k_spatial, ", bs='gp', m=c(2, 0.1))")
  )
  
  # M2: SST + Spatial
  formula_sst <- as.formula(
    paste0("log_chlorophyll ~ s(sst, k=10, bs='tp') + ",
           "lat_band + s(lon_scaled, lat_scaled, k=", k_spatial, ", bs='gp', m=c(2, 0.1))")
  )
  
  # M3: Additive (Herbicide + SST + Spatial)
  formula_additive <- as.formula(
    paste0("log_chlorophyll ~ s(log_herbicide_discharge, k=10, bs='tp') + ",
           "s(sst, k=10, bs='tp') + lat_band + ",
           "s(lon_scaled, lat_scaled, k=", k_spatial, ", bs='gp', m=c(2, 0.1))")
  )
  
  # M4: Interaction (Herbicide × SST tensor product + Spatial)
  formula_interaction <- as.formula(
    paste0("log_chlorophyll ~ te(log_herbicide_discharge, sst, k=10, bs='tp') + ",
           "lat_band + s(lon_scaled, lat_scaled, k=", k_spatial, ", bs='gp', m=c(2, 0.1))")
  )
  
  # Fit models using fast REML
  mod_herbicide <- bam(
    formula_herbicide, 
    data = data, 
    method = "fREML", 
    discrete = TRUE,
    nthreads = 4,
    control = gam.control(trace = FALSE)
  )
  
  mod_sst <- bam(
    formula_sst, 
    data = data, 
    method = "fREML", 
    discrete = TRUE,
    nthreads = 4,
    control = gam.control(trace = FALSE)
  )
  
  mod_additive <- bam(
    formula_additive, 
    data = data, 
    method = "fREML", 
    discrete = TRUE,
    nthreads = 4,
    control = gam.control(trace = FALSE)
  )
  
  mod_interaction <- bam(
    formula_interaction, 
    data = data, 
    method = "fREML", 
    discrete = TRUE,
    nthreads = 4,
    control = gam.control(trace = FALSE)
  )
  
  return(list(
    herbicide = mod_herbicide,
    sst = mod_sst,
    additive = mod_additive,
    interaction = mod_interaction
  ))
}

# ============================================================================
# Residual spatial autocorrelation assessment (Moran's I, k=8)
# ============================================================================

assess_residual_autocorrelation <- function(model, data, k = 8) {
  resid <- residuals(model)
  coords <- data %>% dplyr::select(lon, lat) %>% as.matrix()
  
  if (nrow(data) > 2000) {
    set.seed(123)
    sample_idx <- sample(1:nrow(data), 2000)
    resid_sample <- resid[sample_idx]
    coords_sample <- coords[sample_idx, ]
  } else {
    resid_sample <- resid
    coords_sample <- coords
  }
  
  nb_resid <- knn2nb(knearneigh(coords_sample, k = k))
  lw_resid <- nb2listw(nb_resid, style = "W", zero.policy = TRUE)
  
  moran_result <- tryCatch({
    moran.test(resid_sample, lw_resid, zero.policy = TRUE)
  }, error = function(e) {
    list(estimate = c("Moran I statistic" = NA), p.value = NA)
  })
  
  return(list(
    moran_i = moran_result$estimate[1],
    p_value = moran_result$p.value
  ))
}

# ============================================================================
# Model comparison and selection (AIC)
# ============================================================================

compare_models <- function(models, data) {
  comparison <- data.frame(
    Model = c("M1_Herbicide+Spatial", "M2_SST+Spatial", 
              "M3_Additive+Spatial", "M4_Interaction+Spatial"),
    AIC = c(AIC(models$herbicide), AIC(models$sst), 
            AIC(models$additive), AIC(models$interaction)),
    R2_adj = c(summary(models$herbicide)$r.sq, summary(models$sst)$r.sq, 
               summary(models$additive)$r.sq, summary(models$interaction)$r.sq),
    Dev_expl = c(summary(models$herbicide)$dev.expl, summary(models$sst)$dev.expl,
                 summary(models$additive)$dev.expl, summary(models$interaction)$dev.expl)
  )
  
  comparison$Delta_AIC <- comparison$AIC - min(comparison$AIC)
  
  best_idx <- which.min(comparison$AIC)
  best_model <- list(models$herbicide, models$sst, 
                    models$additive, models$interaction)[[best_idx]]
  
  moran_result <- assess_residual_autocorrelation(best_model, data, k = config$moran_k)
  
  return(list(
    comparison = comparison,
    best_model = best_model,
    best_name = comparison$Model[best_idx],
    moran_i = moran_result$moran_i,
    moran_p = moran_result$p_value
  ))
}

# ============================================================================
# Hierarchical variance partitioning
# ============================================================================

partition_variance <- function(models) {
  dev_herb <- summary(models$herbicide)$dev.expl
  dev_sst <- summary(models$sst)$dev.expl
  dev_add <- summary(models$additive)$dev.expl
  dev_int <- summary(models$interaction)$dev.expl
  
  herb_unique <- max(0, dev_add - dev_sst)
  sst_unique <- max(0, dev_add - dev_herb)
  shared <- max(0, dev_add - herb_unique - sst_unique)
  interaction_term <- dev_int - dev_add
  spatial_resid <- 1 - dev_int
  
  partitioning <- data.frame(
    Component = c("Herbicide_Unique", "SST_Unique", "Shared", 
                  "Interaction", "Spatial_Residual"),
    Deviance_Explained = c(herb_unique, sst_unique, shared, 
                          interaction_term, spatial_resid),
    Percent_of_Total = c(herb_unique, sst_unique, shared, 
                         interaction_term, spatial_resid) * 100,
    Percent_of_Explained = c(herb_unique/dev_int, sst_unique/dev_int, 
                            shared/dev_int, interaction_term/dev_int, 
                            spatial_resid/dev_int) * 100
  )
  
  # Extract F-statistic and p-value for interaction term
  sm <- summary(models$interaction)
  te_row <- grep("^te\\(", rownames(sm$s.table))
  
  if (length(te_row) > 0) {
    f_val <- sm$s.table[te_row[1], "F"]
    p_val <- sm$s.table[te_row[1], "p-value"]
    edf <- sm$s.table[te_row[1], "edf"]
  } else {
    f_val <- NA
    p_val <- NA
    edf <- NA
  }
  
  partitioning$F_stat <- NA
  partitioning$P_value <- NA
  partitioning$EDF <- NA
  partitioning$F_stat[partitioning$Component == "Interaction"] <- f_val
  partitioning$P_value[partitioning$Component == "Interaction"] <- p_val
  partitioning$EDF[partitioning$Component == "Interaction"] <- edf
  
  return(partitioning)
}

# ============================================================================
# Cross-validation (5-fold spatial blocking)
# ============================================================================

spatial_cv_gam <- function(data, k_folds, block_size_m, k_spatial = 100) {
  spatial_folds <- create_spatial_folds(data, k_folds, block_size_m)
  
  cv_results <- list()
  
  for (i in 1:length(spatial_folds$folds)) {
    fold <- spatial_folds$folds[[i]]
    train_idx <- fold[[1]]
    test_idx <- fold[[2]]
    
    train_data <- data[train_idx, ]
    test_data <- data[test_idx, ]
    
    models <- fit_gam_models(train_data, k_spatial)
    
    # Predictions on test set
    test_data_scaled <- test_data %>%
      mutate(
        lon_scaled = (lon - mean(train_data$lon)) / sd(train_data$lon),
        lat_scaled = (lat - mean(train_data$lat)) / sd(train_data$lat)
      )
    
    preds_interaction <- predict(models$interaction, newdata = test_data_scaled)
    
    obs <- test_data$log_chlorophyll
    ss_res <- sum((obs - preds_interaction)^2, na.rm = TRUE)
    ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - (ss_res / ss_tot)
    rmse <- sqrt(mean((obs - preds_interaction)^2, na.rm = TRUE))
    
    cv_results[[i]] <- list(
      fold = i,
      r2 = r2,
      rmse = rmse,
      n_train = nrow(train_data),
      n_test = nrow(test_data)
    )
  }
  
  cv_summary <- data.frame(
    median_r2 = median(sapply(cv_results, function(x) x$r2)),
    mean_r2 = mean(sapply(cv_results, function(x) x$r2)),
    median_rmse = median(sapply(cv_results, function(x) x$rmse)),
    mean_rmse = mean(sapply(cv_results, function(x) x$rmse))
  )
  
  return(list(
    results = cv_results,
    summary = cv_summary
  ))
}

# ============================================================================
# 2050 Projections with Monte Carlo uncertainty (1000 iterations)
# ============================================================================

project_2050_monte_carlo <- function(historical_data, herbicide_2050_pred, 
                                    best_model, n_sims = 1000) {
  baseline_summary <- historical_data %>%
    group_by(lon, lat) %>%
    summarise(
      sst_baseline = mean(sst, na.rm = TRUE),
      chlorophyll_baseline = mean(chlorophyll_mean, na.rm = TRUE),
      log_chlorophyll_baseline = mean(log_chlorophyll, na.rm = TRUE),
      herbicide_baseline = mean(herbicide_discharge, na.rm = TRUE),
      .groups = "drop"
    )
  
  baseline_sf <- st_as_sf(baseline_summary, coords = c("lon", "lat"), crs = 4326)
  herb_2050_sf <- st_as_sf(herbicide_2050_pred, coords = c("lon", "lat"), crs = 4326)
  
  joined <- st_join(baseline_sf, herb_2050_sf, join = st_nearest_feature)
  
  projection_data <- joined %>%
    st_drop_geometry() %>%
    mutate(
      herbicide_2050 = predicted_discharge_2050 * (hbc_discharge_pct / 100),
      log_herbicide_2050 = log10(herbicide_2050 + 0.001),
      sst_2050 = sst_baseline + config$sst_increase_2050
    )
  
  # Monte Carlo simulation
  mc_results <- replicate(n_sims, {
    noise_sst <- rnorm(nrow(projection_data), mean = 0, sd = 0.1)
    noise_herb <- rnorm(nrow(projection_data), mean = 0, sd = 0.1)
    
    sim_data <- projection_data %>%
      mutate(
        log_herbicide_discharge = log_herbicide_2050 + noise_herb,
        sst = sst_2050 + noise_sst,
        lon_scaled = scale(lon)[,1],
        lat_scaled = scale(lat)[,1]
      )
    
    pred <- predict(best_model, newdata = sim_data, type = "response")
    return(pred)
  })
  
  projection_data$predicted_2050_mean <- rowMeans(mc_results, na.rm = TRUE)
  projection_data$predicted_2050_lower <- apply(mc_results, 1, quantile, 0.025, na.rm = TRUE)
  projection_data$predicted_2050_upper <- apply(mc_results, 1, quantile, 0.975, na.rm = TRUE)
  
  projection_data <- projection_data %>%
    mutate(
      chlorophyll_change_pct = ((10^predicted_2050_mean - 10^log_chlorophyll_baseline) / 
                                10^log_chlorophyll_baseline) * 100,
      region = if_else(abs(lat) > 23.5, "Temperate", "Tropical")
    )
  
  return(projection_data)
}

# ============================================================================
# Fit models for all phytoplankton groups
# ============================================================================

all_results <- list()

for (group_name in names(model_data_list)) {
  data <- model_data_list[[group_name]]
  
  if (nrow(data) < 100) next
  
  models <- fit_gam_models(data, k_spatial = config$k_spatial)
  
  model_comparison <- compare_models(models, data)
  
  variance_partition <- partition_variance(models)
  
  cv_results <- spatial_cv_gam(
    data, 
    k_folds = config$cv_folds,
    block_size_m = config$spatial_block_size_m,
    k_spatial = config$k_spatial
  )
  
  projections_2050 <- project_2050_monte_carlo(
    data,
    herbicide_2050_predictions,
    model_comparison$best_model,
    n_sims = config$monte_carlo_sims
  )
  
  all_results[[group_name]] <- list(
    models = models,
    comparison = model_comparison$comparison,
    best_model = model_comparison$best_model,
    best_name = model_comparison$best_name,
    moran_i = model_comparison$moran_i,
    moran_p = model_comparison$moran_p,
    variance_partition = variance_partition,
    cv_results = cv_results,
    projections_2050 = projections_2050
  )
}

# ============================================================================
# Buffer sensitivity analysis (25, 50, 75 km)
# ============================================================================

buffer_sensitivity_analysis <- function(group_data, discharge_data, buffers = c(25, 50, 75)) {
  sensitivity_results <- list()
  
  for (buffer_km in buffers) {
    chlo_agg <- aggregate_to_discharge_buffer(group_data, discharge_data, buffer_km)
    
    if (nrow(chlo_agg) < 100) next
    
    chlo_with_sst <- extract_sst_for_group(chlo_agg, sst_stack_common)
    model_data <- merge_all_data(discharge_data, chlo_with_sst, chlo_with_sst, hbc_discharge_pct)
    
    models <- fit_gam_models(model_data, k_spatial = config$k_spatial)
    comparison <- compare_models(models, model_data)
    
    sensitivity_results[[paste0("buffer_", buffer_km, "km")]] <- list(
      buffer_km = buffer_km,
      n_obs = nrow(model_data),
      r2 = comparison$comparison$R2_adj[comparison$comparison$Model == "M4_Interaction+Spatial"],
      aic = comparison$comparison$AIC[comparison$comparison$Model == "M4_Interaction+Spatial"]
    )
  }
  
  return(sensitivity_results)
}

# ============================================================================
# Nutrient confounding analysis
# ============================================================================

nutrient_confounding_analysis <- function(model_data_with_nutrients) {
  # VIF assessment
  vif_check <- car::vif(lm(log_chlorophyll ~ log_herbicide_discharge + sst + 
                            nitrate + phosphate + silicate, 
                          data = model_data_with_nutrients))
  
  # Nested model comparison
  mod_nutrient_only <- bam(
    log_chlorophyll ~ s(nitrate, k=10) + s(phosphate, k=10) + s(silicate, k=10) +
      lat_band + s(lon_scaled, lat_scaled, k=50, bs='gp', m=c(2, 0.1)),
    data = model_data_with_nutrients,
    method = "fREML"
  )
  
  mod_herbicide_temp <- bam(
    log_chlorophyll ~ s(log_herbicide_discharge, k=10) + s(sst, k=10) +
      lat_band + s(lon_scaled, lat_scaled, k=50, bs='gp', m=c(2, 0.1)),
    data = model_data_with_nutrients,
    method = "fREML"
  )
  
  mod_all_additive <- bam(
    log_chlorophyll ~ s(log_herbicide_discharge, k=10) + s(sst, k=10) +
      s(nitrate, k=10) + s(phosphate, k=10) + s(silicate, k=10) +
      lat_band + s(lon_scaled, lat_scaled, k=50, bs='gp', m=c(2, 0.1)),
    data = model_data_with_nutrients,
    method = "fREML"
  )
  
  # Akaike weights
  aic_vals <- c(AIC(mod_nutrient_only), AIC(mod_herbicide_temp), AIC(mod_all_additive))
  delta_aic <- aic_vals - min(aic_vals)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  
  evidence_ratio <- weights[2] / weights[1]
  
  return(list(
    vif = vif_check,
    akaike_weights = weights,
    evidence_ratio = evidence_ratio,
    delta_r2 = summary(mod_all_additive)$r.sq - summary(mod_nutrient_only)$r.sq
  ))
}

