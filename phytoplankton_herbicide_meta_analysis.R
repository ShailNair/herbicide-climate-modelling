# ============================================================================
# Phytoplankton Herbicide Meta-Analysis
# ============================================================================

rm(list = ls())

library(readxl)
library(dplyr)
library(metafor)
library(robumeta)
library(brms)

# --- Configuration ---
config <- list(
  concentration_range = c(0.05, 100),
  min_sample_size = 3,
  max_cv = 1.0,
  max_lnrr = 3,
  min_se = 0.01,
  max_se = 2,
  variance_quantiles = c(0.01, 0.99),
  random_seed = 123,
  bayesian_iterations = 2000,
  bayesian_warmup = 1000,
  bayesian_chains = 4
)

set.seed(config$random_seed)

# ============================================================================
# Data loading and preprocessing
# ============================================================================

herbicide_data <- read_excel("Herbicide_toxicity_data.xlsx", sheet = "Sheet2")

herbicide_clean <- herbicide_data %>%
  rename(
    exposure_ugL = `Exposure conc. (ug/L)`,
    control_mean_cells = `Control (Mean), cells/ml*10^4`,
    control_sd_cells = `Control (SD)...25`,
    control_n_cells = `Control (n)...26`,
    treatment_mean_cells = `Treatment_Mean...27`,
    treatment_sd_cells = `Treatment_SD...28`,
    treatment_n_cells = `Treatment (n)...29`
  ) %>%
  mutate(
    control_mean = case_when(
      !is.na(control_mean_cells) ~ as.numeric(control_mean_cells),
      !is.na(`Control (Mean)`) ~ as.numeric(`Control (Mean)`),
      TRUE ~ NA_real_
    ),
    control_sd = case_when(
      !is.na(control_sd_cells) ~ as.numeric(control_sd_cells),
      !is.na(`Control (SD)...19`) ~ as.numeric(`Control (SD)...19`),
      TRUE ~ NA_real_
    ),
    control_n = case_when(
      !is.na(control_n_cells) ~ as.numeric(control_n_cells),
      !is.na(`Control (n)...20`) ~ as.numeric(`Control (n)...20`),
      TRUE ~ NA_real_
    ),
    treatment_mean = case_when(
      !is.na(treatment_mean_cells) ~ as.numeric(treatment_mean_cells),
      !is.na(`Treatment_Mean...21`) ~ as.numeric(`Treatment_Mean...21`),
      TRUE ~ NA_real_
    ),
    treatment_sd = case_when(
      !is.na(treatment_sd_cells) ~ as.numeric(treatment_sd_cells),
      !is.na(`Treatment_SD...22`) ~ as.numeric(`Treatment_SD...22`),
      TRUE ~ NA_real_
    ),
    treatment_n = case_when(
      !is.na(treatment_n_cells) ~ as.numeric(treatment_n_cells),
      !is.na(`Treatment (n)...23`) ~ as.numeric(`Treatment (n)...23`),
      TRUE ~ NA_real_
    ),
    exposure_ugL = as.numeric(exposure_ugL)
  ) %>%
  filter(
    exposure_ugL >= config$concentration_range[1],
    exposure_ugL <= config$concentration_range[2],
    !is.na(control_mean) & !is.na(treatment_mean),
    control_mean > 0 & treatment_mean > 0,
    !is.na(control_sd) & !is.na(treatment_sd),
    !is.na(control_n) & !is.na(treatment_n),
    control_n >= config$min_sample_size,
    treatment_n >= config$min_sample_size
  )

# ============================================================================
# Effect size calculation (log Response Ratio)
# ============================================================================

herbicide_es <- herbicide_clean %>%
  mutate(
    lnRR = log(treatment_mean / control_mean),
    var_lnRR = (treatment_sd^2) / (treatment_n * treatment_mean^2) + 
               (control_sd^2) / (control_n * control_mean^2),
    se_lnRR = sqrt(var_lnRR),
    
    hedges_g_pooled_sd = sqrt(((control_n - 1) * control_sd^2 + 
                               (treatment_n - 1) * treatment_sd^2) / 
                              (control_n + treatment_n - 2)),
    cohens_d = (treatment_mean - control_mean) / hedges_g_pooled_sd,
    hedges_g = cohens_d * (1 - 3 / (4 * (control_n + treatment_n) - 9)),
    var_hedges_g = ((control_n + treatment_n) / (control_n * treatment_n)) + 
                   (hedges_g^2 / (2 * (control_n + treatment_n))),
    
    percent_change = ((treatment_mean - control_mean) / control_mean) * 100,
    
    log_conc = log10(exposure_ugL),
    
    conc_category = case_when(
      exposure_ugL < 1 ~ "Background (<1 µg/L)",
      exposure_ugL >= 1 & exposure_ugL < 5 ~ "Low (1-5 µg/L)",
      exposure_ugL >= 5 & exposure_ugL < 20 ~ "Moderate (5-20 µg/L)",
      exposure_ugL >= 20 & exposure_ugL < 50 ~ "High (20-50 µg/L)",
      exposure_ugL >= 50 ~ "Very High (≥50 µg/L)"
    ),
    
    total_n = control_n + treatment_n,
    cv_control = control_sd / control_mean,
    cv_treatment = treatment_sd / treatment_mean,
    
    weight = 1 / var_lnRR
  ) %>%
  filter(
    !is.na(lnRR) & !is.na(var_lnRR),
    is.finite(lnRR) & is.finite(var_lnRR),
    var_lnRR > 0
  )

# ============================================================================
# Quality filtering
# ============================================================================
herbicide_es_robust <- herbicide_es %>%
  filter(
    var_lnRR > quantile(var_lnRR, config$variance_quantiles[1], na.rm = TRUE),
    var_lnRR < quantile(var_lnRR, config$variance_quantiles[2], na.rm = TRUE),
    
    control_n >= config$min_sample_size,
    treatment_n >= config$min_sample_size,
    
    abs(lnRR) < config$max_lnrr,
    
    cv_control < config$max_cv,
    cv_treatment < config$max_cv,
    
    se_lnRR < config$max_se,
    se_lnRR > config$min_se
  ) %>%
  mutate(
    var_lnRR_robust = pmax(var_lnRR, 0.001),
    se_lnRR_robust = sqrt(var_lnRR_robust),
    weight_robust = 1 / var_lnRR_robust,
    
    Study = paste(Author, Year, sep = "_")
  )

# ============================================================================
# Primary random-effects meta-analysis (REML)
# ============================================================================
meta_overall <- rma(
  yi = lnRR,
  vi = var_lnRR_robust,
  data = herbicide_es_robust,
  method = "REML"
)

overall_results <- data.frame(
  estimate_lnRR = meta_overall$beta[1],
  se = meta_overall$se,
  ci_lower = meta_overall$ci.lb,
  ci_upper = meta_overall$ci.ub,
  percent_change = (exp(meta_overall$beta[1]) - 1) * 100,
  percent_lower = (exp(meta_overall$ci.lb) - 1) * 100,
  percent_upper = (exp(meta_overall$ci.ub) - 1) * 100,
  tau2 = meta_overall$tau2,
  I2 = meta_overall$I2,
  H2 = meta_overall$H2,
  Q = meta_overall$QE,
  Q_p = meta_overall$QEp,
  k = meta_overall$k
)

# ============================================================================
# Three-level hierarchical model
# ============================================================================
meta_3level <- rma.mv(
  yi = lnRR,
  V = var_lnRR_robust,
  random = list(~1 | Study, ~1 | Species),
  data = herbicide_es_robust,
  method = "REML"
)

# ============================================================================
# Meta-regression (concentration-dependent response)
# ============================================================================
meta_conc <- rma(
  yi = lnRR,
  vi = var_lnRR_robust,
  mods = ~ log_conc,
  data = herbicide_es_robust,
  method = "REML"
)

conc_regression_results <- data.frame(
  intercept = meta_conc$beta[1],
  slope = meta_conc$beta[2],
  slope_se = meta_conc$se[2],
  slope_p = meta_conc$pval[2],
  R2 = meta_conc$R2,
  tau2 = meta_conc$tau2
)

# ============================================================================
# Subgroup analysis by taxonomy
# ============================================================================

taxonomy_groups <- herbicide_es_robust %>%
  filter(!is.na(Group)) %>%
  group_by(Group) %>%
  filter(n() >= 3) %>%
  ungroup()

meta_by_taxonomy <- taxonomy_groups %>%
  group_by(Group) %>%
  nest() %>%
  mutate(
    meta = map(data, ~rma(yi = lnRR, vi = var_lnRR_robust, data = .x, method = "REML")),
    k = map_dbl(data, nrow),
    estimate = map_dbl(meta, ~.x$beta[1]),
    se = map_dbl(meta, ~.x$se),
    ci_lower = map_dbl(meta, ~.x$ci.lb),
    ci_upper = map_dbl(meta, ~.x$ci.ub),
    percent_change = (exp(estimate) - 1) * 100,
    percent_lower = (exp(ci_lower) - 1) * 100,
    percent_upper = (exp(ci_upper) - 1) * 100,
    tau2 = map_dbl(meta, ~.x$tau2),
    I2 = map_dbl(meta, ~.x$I2),
    Q = map_dbl(meta, ~.x$QE),
    Q_p = map_dbl(meta, ~.x$QEp)
  ) %>%
  select(-data, -meta)

# ============================================================================
# Subgroup Analysis by herbicide class
# ============================================================================

herbicide_class_groups <- herbicide_es_robust %>%
  filter(!is.na(Herbicide_class)) %>%
  group_by(Herbicide_class) %>%
  filter(n() >= 3) %>%
  ungroup()

meta_by_herbicide <- herbicide_class_groups %>%
  group_by(Herbicide_class) %>%
  nest() %>%
  mutate(
    meta = map(data, ~rma(yi = lnRR, vi = var_lnRR_robust, data = .x, method = "REML")),
    k = map_dbl(data, nrow),
    estimate = map_dbl(meta, ~.x$beta[1]),
    se = map_dbl(meta, ~.x$se),
    ci_lower = map_dbl(meta, ~.x$ci.lb),
    ci_upper = map_dbl(meta, ~.x$ci.ub),
    percent_change = (exp(estimate) - 1) * 100,
    tau2 = map_dbl(meta, ~.x$tau2),
    I2 = map_dbl(meta, ~.x$I2)
  ) %>%
  select(-data, -meta)

# ============================================================================
# Subgroup analysis: low concentration (<1 µg/L)
# ============================================================================

low_conc_data <- herbicide_es_robust %>%
  filter(exposure_ugL < 1)

if (nrow(low_conc_data) >= 3) {
  meta_low_conc <- rma(
    yi = lnRR,
    vi = var_lnRR_robust,
    data = low_conc_data,
    method = "REML"
  )
  
  low_conc_results <- data.frame(
    estimate_lnRR = meta_low_conc$beta[1],
    ci_lower = meta_low_conc$ci.lb,
    ci_upper = meta_low_conc$ci.ub,
    percent_change = (exp(meta_low_conc$beta[1]) - 1) * 100,
    percent_lower = (exp(meta_low_conc$ci.lb) - 1) * 100,
    percent_upper = (exp(meta_low_conc$ci.ub) - 1) * 100,
    k = meta_low_conc$k,
    tau2 = meta_low_conc$tau2,
    I2 = meta_low_conc$I2
  )
}

# ============================================================================
# Publication bias assessment
# ============================================================================

# Egger's regression test
egger_test <- regtest(meta_overall, model = "lm")

# Trim and fill analysis
trimfill_analysis <- trimfill(meta_overall)

# PET (Precision Effect Test)
pet_model <- rma(
  yi = lnRR,
  vi = var_lnRR_robust,
  mods = ~ se_lnRR_robust,
  method = "REML",
  data = herbicide_es_robust
)

# PEESE (Precision effect estimate with standard Error)
peese_model <- rma(
  yi = lnRR,
  vi = var_lnRR_robust,
  mods = ~ I(se_lnRR_robust^2),
  method = "REML",
  data = herbicide_es_robust
)

# Selection model
selection_model <- selmodel(
  meta_overall,
  type = "stepfun",
  steps = c(0.025, 0.05, 0.5, 1)
)

publication_bias_results <- data.frame(
  test = c("Egger", "PET_intercept", "PEESE_intercept", "TrimFill"),
  estimate = c(
    egger_test$fit$estimate[1],
    pet_model$beta[1],
    peese_model$beta[1],
    trimfill_analysis$beta[1]
  ),
  p_value = c(
    egger_test$pval,
    pet_model$pval[1],
    peese_model$pval[1],
    NA
  ),
  k_trimmed = c(NA, NA, NA, trimfill_analysis$k0)
)

# ============================================================================
# Sensitivity analysis
# ============================================================================
# leave-one-out
loo_results <- leave1out(meta_overall)

loo_summary <- data.frame(
  estimate = loo_results$estimate,
  se = loo_results$se,
  ci_lb = loo_results$ci.lb,
  ci_ub = loo_results$ci.ub,
  percent_change = (exp(loo_results$estimate) - 1) * 100,
  tau2 = loo_results$tau2,
  I2 = loo_results$I2,
  Q = loo_results$Q,
  Qp = loo_results$Qp
)

loo_stats <- data.frame(
  min_estimate = min(loo_summary$estimate),
  max_estimate = max(loo_summary$estimate),
  min_percent = min(loo_summary$percent_change),
  max_percent = max(loo_summary$percent_change),
  range_percent = max(loo_summary$percent_change) - min(loo_summary$percent_change),
  all_significant = all(loo_summary$Qp < 0.05)
)


# Robust variance estimation

rve_model <- robu(
  formula = lnRR ~ log_conc,
  data = herbicide_es_robust,
  studynum = Study,
  var.eff.size = var_lnRR_robust,
  modelweights = "HIER"
)

rve_results <- data.frame(
  estimate = rve_model$reg_table$b.r[1],
  se = rve_model$reg_table$SE[1],
  ci_lower = rve_model$reg_table$CI.L[1],
  ci_upper = rve_model$reg_table$CI.U[1],
  p_value = rve_model$reg_table$prob[1],
  tau2 = rve_model$mod_info$tau.sq
)


# Alternative effect size (Hedges' g)

meta_hedges <- rma(
  yi = hedges_g,
  vi = var_hedges_g,
  data = herbicide_es_robust,
  method = "REML"
)

hedges_comparison <- data.frame(
  Metric = c("lnRR", "Hedges_g"),
  Estimate = c(meta_overall$beta[1], meta_hedges$beta[1]),
  CI_lower = c(meta_overall$ci.lb, meta_hedges$ci.lb),
  CI_upper = c(meta_overall$ci.ub, meta_hedges$ci.ub),
  I2 = c(meta_overall$I2, meta_hedges$I2)
)


# Bayesian hierarchical model

bayesian_model <- brm(
  lnRR | se(se_lnRR_robust) ~ log_conc * Group + (1 + log_conc | Group),
  data = herbicide_es_robust,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(normal(0, 0.5), class = b),
    prior(exponential(1), class = sd),
    prior(lkj(2), class = cor)
  ),
  chains = config$bayesian_chains,
  iter = config$bayesian_iterations,
  warmup = config$bayesian_warmup,
  cores = config$bayesian_chains,
  seed = config$random_seed
)

# ============================================================================
# Summary stats
# ============================================================================

summary_stats <- data.frame(
  n_effect_sizes = nrow(herbicide_es_robust),
  n_species = length(unique(herbicide_es_robust$Species)),
  n_studies = length(unique(herbicide_es_robust$Study)),
  conc_min = min(herbicide_es_robust$exposure_ugL),
  conc_max = max(herbicide_es_robust$exposure_ugL),
  conc_median = median(herbicide_es_robust$exposure_ugL)
)
