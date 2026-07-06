# Reanalysis of the European earwig egg size-number data
# =============================================================================
# This script is designed to replace the covariance and repeatability sections
# of the rejected-manuscript analysis. It treats the two broods as the complete
# biologically relevant reproductive sequence and asks three complementary
# questions:
#   1. Is egg size-number covariance expressed among or within females?
#   2. Does the association differ between first and second broods?
#   3. Do females compensate for changes in clutch size by changing egg size?
#
# Required input: df.rds (one row per female x brood)

#rm(list = ls())

library(tidyverse)
library(brms)
library(posterior)
library(tidybayes)
library(patchwork)
library(flextable)

options(mc.cores = min(4, parallel::detectCores()))
# Note: this seed applies only to the cluster bootstrap functions below.
# Each brm() call uses its own seed argument and is unaffected by set.seed().
set.seed(123)

# =============================================================================
# Helper functions
# =============================================================================

fmt <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits)
}

summarise_draw <- function(x, term) {
  # Coerce posterior columns to ordinary numeric vectors before summarising.
  # This avoids harmless draws_df metadata warnings when summaries are combined.
  x <- as.numeric(x)

  tibble(
    term = term,
    estimate = median(x),
    lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975)),
    p_negative = mean(x < 0),
    p_positive = mean(x > 0)
  )
}

est <- function(x, term, digits = 2) {
  fmt(x$estimate[match(term, x$term)], digits)
}

lcl <- function(x, term, digits = 2) {
  fmt(x$lower[match(term, x$term)], digits)
}

ucl <- function(x, term, digits = 2) {
  fmt(x$upper[match(term, x$term)], digits)
}

pneg <- function(x, term, digits = 2) {
  fmt(x$p_negative[match(term, x$term)], digits)
}

cluster_boot_cor <- function(data, x, y, cluster, R = 5000, seed = 123) {
  set.seed(seed)
  ids <- unique(data[[cluster]])
  n   <- length(ids)

  r_obs <- cor(data[[x]], data[[y]], use = "complete.obs")

  # Build index list once, then sample by position (avoids repeated data-frame
  # construction inside the loop)
  idx_by_id <- lapply(ids, function(i) which(data[[cluster]] == i))
  xv <- data[[x]]
  yv <- data[[y]]

  r_boot <- replicate(R, {
    rows <- unlist(idx_by_id[sample(n, n, replace = TRUE)], use.names = FALSE)
    cor(xv[rows], yv[rows], use = "complete.obs")
  })

  tibble(
    estimate = r_obs,
    lower = unname(quantile(r_boot, 0.025, na.rm = TRUE)),
    upper = unname(quantile(r_boot, 0.975, na.rm = TRUE))
  )
}

# =============================================================================
# Load and prepare data
# =============================================================================

data_path <- case_when(
  file.exists("data/processed/df.rds") ~ "data/processed/df.rds",
  file.exists("df.rds") ~ "df.rds",
  file.exists("/mnt/data/df.rds") ~ "/mnt/data/df.rds",
  TRUE ~ NA_character_
)

if (is.na(data_path)) {
  stop("Could not find df.rds. Put it in data/processed/ or the working directory.")
}

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

df <- readRDS(data_path) %>%
  mutate(
    id = factor(id),
    brood = factor(brood, levels = c("one", "two")),
    brood_label = factor(
      brood,
      levels = c("one", "two"),
      labels = c("First brood", "Second brood")
    ),
    log_egg_number = log(egg.number),
    log_egg_size = log(mean.egg.size),
    mean_egg_size_z = as.numeric(scale(mean.egg.size))
  )

# Pooled scaling is used in the primary brood-specific model so that a one-SD
# difference in log egg number has the same meaning in both reproductive bouts.
mean_log_egg_number <- mean(df$log_egg_number, na.rm = TRUE)
sd_log_egg_number <- sd(df$log_egg_number, na.rm = TRUE)
mean_log_egg_size <- mean(df$log_egg_size, na.rm = TRUE)
sd_log_egg_size <- sd(df$log_egg_size, na.rm = TRUE)

df <- df %>%
  mutate(
    log_egg_number_z =
      (log_egg_number - mean_log_egg_number) / sd_log_egg_number,
    log_egg_size_z =
      (log_egg_size - mean_log_egg_size) / sd_log_egg_size
  )

# Brood-specific scaling is retained only for a sensitivity analysis requested
# during manuscript review. It is not the scale used in the primary model.
brood_scaling <- df %>%
  group_by(brood) %>%
  summarise(
    mean_log_n = mean(log_egg_number, na.rm = TRUE),
    sd_log_n = sd(log_egg_number, na.rm = TRUE),
    .groups = "drop"
  )

df <- df %>%
  left_join(brood_scaling, by = "brood") %>%
  mutate(
    log_egg_number_brood_z =
      (log_egg_number - mean_log_n) / sd_log_n
  )

# Female pronotum length is a fixed structural trait in adult earwigs. It is
# standardized once at the female level and then joined to both brood records.
size_by_female <- df %>%
  distinct(id, mean.pro) %>%
  mutate(
    log_pronotum_z = as.numeric(scale(log(mean.pro)))
  )

df <- df %>%
  select(-any_of("log_pronotum_z")) %>%
  left_join(
    size_by_female %>% select(id, log_pronotum_z),
    by = "id"
  )

# Within-between decomposition on both the absolute and proportional scales.
# Standardization is performed after decomposition so that coefficients are in
# SD units while retaining the exact within- and among-female separation.
df <- df %>%
  group_by(id) %>%
  mutate(
    egg_number_between = mean(egg.number),
    egg_number_within = egg.number - egg_number_between,
    log_egg_number_between = mean(log_egg_number),
    log_egg_number_within = log_egg_number - log_egg_number_between
  ) %>%
  ungroup()

raw_between_values <- df %>%
  distinct(id, egg_number_between) %>%
  pull(egg_number_between)

log_between_values <- df %>%
  distinct(id, log_egg_number_between) %>%
  pull(log_egg_number_between)

df <- df %>%
  mutate(
    egg_number_within_z = egg_number_within / sd(egg_number_within),
    egg_number_between_z =
      (egg_number_between - mean(raw_between_values)) /
      sd(raw_between_values),
    log_egg_number_within_z =
      log_egg_number_within / sd(log_egg_number_within),
    log_egg_number_between_z =
      (log_egg_number_between - mean(log_between_values)) /
      sd(log_between_values)
  )

# Change scores. These are transparent here because two broods are the
# biologically relevant repeated reproductive sequence rather than a sparse
# sample from a long time series.
change_dat <- df %>%
  select(
    id,
    brood,
    egg.number,
    mean.egg.size,
    log_egg_number,
    log_egg_size,
    log_pronotum_z
  ) %>%
  pivot_wider(
    names_from = brood,
    values_from = c(
      egg.number,
      mean.egg.size,
      log_egg_number,
      log_egg_size,
      log_pronotum_z
    )
  ) %>%
  mutate(
    # The two joined pronotum values should be identical; retain one copy.
    log_pronotum_z = log_pronotum_z_one,
    delta_egg_number = egg.number_two - egg.number_one,
    delta_egg_size = mean.egg.size_two - mean.egg.size_one,
    delta_log_egg_number = log_egg_number_two - log_egg_number_one,
    delta_log_egg_size = log_egg_size_two - log_egg_size_one,
    delta_egg_number_z = as.numeric(scale(delta_egg_number)),
    delta_egg_size_z = as.numeric(scale(delta_egg_size)),
    delta_log_egg_number_z = as.numeric(scale(delta_log_egg_number)),
    delta_log_egg_size_z = as.numeric(scale(delta_log_egg_size))
  ) %>%
  select(-log_pronotum_z_one, -log_pronotum_z_two)

# =============================================================================
# Descriptive statistics and raw correlations
# =============================================================================

summary_stat <- df %>%
  group_by(brood, brood_label) %>%
  summarise(
    n = n(),
    mean_number = mean(egg.number),
    se_number = sd(egg.number) / sqrt(n()),
    mean_size = mean(mean.egg.size),
    se_size = sd(mean.egg.size) / sqrt(n()),
    .groups = "drop"
  )


change_direction_summary <- change_dat %>%
  summarise(
    n_females = n(),
    n_number_decline = sum(delta_egg_number < 0, na.rm = TRUE),
    n_size_decline = sum(delta_egg_size < 0, na.rm = TRUE),
    n_both_decline = sum(
      delta_egg_number < 0 & delta_egg_size < 0,
      na.rm = TRUE
    ),
    prop_number_decline = n_number_decline / n_females,
    prop_size_decline = n_size_decline / n_females,
    prop_both_decline = n_both_decline / n_females
  )

raw_cor_summary <- bind_rows(
  # Cluster resampling retains the paired brood structure in the pooled analysis.
  cluster_boot_cor(df, "egg.number", "mean.egg.size", "id") %>%
    mutate(term = "pooled"),
  # For brood-specific subsets each female contributes exactly one observation,
  # so cluster resampling is equivalent to a standard bootstrap here.
  cluster_boot_cor(filter(df, brood == "one"),
                   "egg.number", "mean.egg.size", "id") %>%
    mutate(term = "first"),
  cluster_boot_cor(filter(df, brood == "two"),
                   "egg.number", "mean.egg.size", "id") %>%
    mutate(term = "second")
) %>%
  select(term, everything())

# =============================================================================
# Model 1: bivariate Gaussian covariance decomposition
# =============================================================================
# Both traits are log transformed and standardized. The female-level random
# correlation estimates among-female covariance; the residual correlation
# estimates within-female covariance after accounting for the mean brood shift.

bf_number <- bf(
  log_egg_number_z ~ brood + log_pronotum_z + (1 | p | id),
  family = gaussian()
)

bf_size <- bf(
  log_egg_size_z ~ brood + log_pronotum_z + (1 | p | id),
  family = gaussian()
)

bivariate_formula <- bf_number + bf_size + set_rescor(TRUE)

prior_bivariate <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  
  # Female-level standard deviations for each response
  prior(
    student_t(3, 0, 1),
    class = "sd",
    group = "id",
    resp = "logeggnumberz"
  ),
  prior(
    student_t(3, 0, 1),
    class = "sd",
    group = "id",
    resp = "logeggsizez"
  ),
  
  # Residual standard deviations for each response
  prior(
    student_t(3, 0, 1),
    class = "sigma",
    resp = "logeggnumberz"
  ),
  prior(
    student_t(3, 0, 1),
    class = "sigma",
    resp = "logeggsizez"
  ),
  
  # Among-female and residual correlation matrices
  prior(lkj(2), class = "cor", group = "id"),
  prior(lkj(2), class = "rescor")
)

fit_bivariate <- brm(
  bf_number + bf_size + set_rescor(TRUE),
  data = df,
  prior = prior_bivariate,
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 6000,
  warmup = 2000,
  backend = "cmdstanr",
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  file = "data/processed/fit_bivariate_revised",
  file_refit = "on_change"
)

biv_draws <- as_draws_df(fit_bivariate)

among_cor_name <- grep(
  "^cor_id__.*Intercept__.*Intercept$",
  names(biv_draws),
  value = TRUE
)
within_cor_name <- grep("^rescor__", names(biv_draws), value = TRUE)

if (length(among_cor_name) != 1 || length(within_cor_name) != 1) {
  stop("Could not uniquely identify the bivariate correlation parameters.")
}

biv_cor_summary <- bind_rows(
  summarise_draw(biv_draws[[among_cor_name]], "among_female"),
  summarise_draw(biv_draws[[within_cor_name]], "within_female")
)

# Brood effects in each response
number_brood_name <- grep("^b_logeggnumberz_broodtwo$", names(biv_draws), value = TRUE)
size_brood_name <- grep("^b_logeggsizez_broodtwo$", names(biv_draws), value = TRUE)

biv_brood_summary <- bind_rows(
  summarise_draw(biv_draws[[number_brood_name]], "egg_number_brood"),
  summarise_draw(biv_draws[[size_brood_name]], "egg_size_brood")
)

# =============================================================================
# Model 2: brood-specific size-number slopes
# =============================================================================
# The primary model uses one pooled standardization of log egg number. The
# predictor therefore has the same scale in both broods, and the interaction is
# directly interpretable as a difference between brood-specific slopes.

prior_context <- c(
  prior(normal(0, 0.7), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(student_t(3, 0, 1), class = "sigma")
)

fit_context <- brm(
  log_egg_size_z ~
    brood * log_egg_number_z +
    log_pronotum_z +
    (1 | id),
  data = df,
  family = gaussian(),
  prior = prior_context,
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 5000,
  warmup = 1500,
  backend = "cmdstanr",
  seed = 124,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  file = "data/processed/fit_context_pooled_scaling",
  file_refit = "on_change"
)

context_draws <- as_draws_df(fit_context)
first_slope <- context_draws$b_log_egg_number_z
slope_difference <- context_draws[["b_broodtwo:log_egg_number_z"]]
second_slope <- first_slope + slope_difference

context_slope_summary <- bind_rows(
  summarise_draw(first_slope, "first"),
  summarise_draw(second_slope, "second"),
  summarise_draw(slope_difference, "second_minus_first")
)

# Sensitivity analysis using separate within-brood standardization. This model
# asks about one brood-specific SD of log egg number rather than one common SD.
fit_context_within_brood_scaling <- brm(
  log_egg_size_z ~
    brood * log_egg_number_brood_z +
    log_pronotum_z +
    (1 | id),
  data = df,
  family = gaussian(),
  prior = prior_context,
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 5000,
  warmup = 1500,
  backend = "cmdstanr",
  seed = 1241,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  file = "data/processed/fit_context_within_brood_scaling",
  file_refit = "on_change"
)

context_wbscale_draws <- as_draws_df(fit_context_within_brood_scaling)
first_slope_wbscale <-
  context_wbscale_draws$b_log_egg_number_brood_z
slope_difference_wbscale <-
  context_wbscale_draws[["b_broodtwo:log_egg_number_brood_z"]]
second_slope_wbscale <-
  first_slope_wbscale + slope_difference_wbscale

context_scaling_sensitivity_summary <- bind_rows(
  summarise_draw(first_slope_wbscale, "first"),
  summarise_draw(second_slope_wbscale, "second"),
  summarise_draw(
    slope_difference_wbscale,
    "second_minus_first"
  )
)

# Posterior predictions in original units for an interpretable second-brood
# contrast. Values of 20 and 40 eggs are well within the observed range.
second_brood_contrast_grid <- tibble(
  brood = factor(c("two", "two"), levels = levels(df$brood)),
  egg.number = c(20, 40),
  log_egg_number_z =
    (log(egg.number) - mean_log_egg_number) / sd_log_egg_number,
  log_pronotum_z = 0
)

second_brood_contrast_draws <- second_brood_contrast_grid %>%
  add_epred_draws(
    fit_context,
    re_formula = NA
  ) %>%
  mutate(
    predicted_log_size =
      .epred * sd_log_egg_size + mean_log_egg_size,
    predicted_size = exp(predicted_log_size)
  )

# Summarise the predicted perimeter separately at 20 and 40 eggs
second_brood_level_summary <- second_brood_contrast_draws %>%
  ungroup() %>%
  group_by(egg.number) %>%
  summarise(
    estimate = median(predicted_size, na.rm = TRUE),
    lower = quantile(
      predicted_size,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      predicted_size,
      probs = 0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    term = paste0("eggs_", as.integer(egg.number))
  ) %>%
  select(
    term,
    estimate,
    lower,
    upper
  )

# Pair the predictions at 20 and 40 eggs within each posterior draw
second_brood_paired_draws <- second_brood_contrast_draws %>%
  ungroup() %>%
  transmute(
    .draw,
    term = paste0("eggs_", as.integer(egg.number)),
    predicted_size
  ) %>%
  pivot_wider(
    id_cols = .draw,
    names_from = term,
    values_from = predicted_size,
    values_fn = mean
  ) %>%
  filter(
    !is.na(eggs_20),
    !is.na(eggs_40)
  ) %>%
  mutate(
    difference = eggs_40 - eggs_20
  )

# Summarise the posterior difference between 40 and 20 eggs
second_brood_difference_summary <- second_brood_paired_draws %>%
  summarise(
    estimate = median(difference, na.rm = TRUE),
    lower = quantile(
      difference,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      difference,
      probs = 0.975,
      na.rm = TRUE
    )
  ) %>%
  mutate(
    term = "difference_40_minus_20"
  ) %>%
  select(
    term,
    estimate,
    lower,
    upper
  )

# Combine the two predicted values and their difference
second_brood_contrast_summary <- bind_rows(
  second_brood_level_summary,
  second_brood_difference_summary
)

second_brood_contrast_summary

# =============================================================================
# Model 3: within-between regression sensitivity analyses
# =============================================================================
# These models preserve the original biological interpretation but avoid the
# redundant joint response model. The raw-count version addresses absolute
# changes in clutch size; the log-count version addresses proportional changes.

prior_wb <- c(
  prior(normal(0, 0.7), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(student_t(3, 0, 1), class = "sigma")
)

fit_wb_raw <- brm(
  mean_egg_size_z ~ brood + egg_number_within_z + egg_number_between_z +
    log_pronotum_z + (1 | id),
  data = df,
  family = gaussian(),
  prior = prior_wb,
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 5000,
  warmup = 1500,
  backend = "cmdstanr",
  seed = 125,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  file = "data/processed/fit_wb_raw_pronotum",
  file_refit = "on_change"
)

fit_wb_log <- brm(
  log_egg_size_z ~ brood + log_egg_number_within_z +
    log_egg_number_between_z + log_pronotum_z + (1 | id),
  data = df,
  family = gaussian(),
  prior = prior_wb,
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 5000,
  warmup = 1500,
  backend = "cmdstanr",
  seed = 126,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  file = "data/processed/fit_wb_log_pronotum",
  file_refit = "on_change"
)

wb_raw_draws <- as_draws_df(fit_wb_raw)
wb_log_draws <- as_draws_df(fit_wb_log)

wb_summary <- bind_rows(
  summarise_draw(wb_raw_draws$b_egg_number_within_z, "raw_within"),
  summarise_draw(wb_raw_draws$b_egg_number_between_z, "raw_between"),
  summarise_draw(wb_log_draws$b_log_egg_number_within_z, "log_within"),
  summarise_draw(wb_log_draws$b_log_egg_number_between_z, "log_between")
)

# =============================================================================
# Model 4: change-score sensitivity analyses
# =============================================================================

fit_change_raw <- brm(
  delta_egg_size_z ~ delta_egg_number_z,
  data = change_dat,
  family = gaussian(),
  prior = c(
    prior(normal(0, 0.7), class = "b"),
    prior(normal(0, 1), class = "Intercept"),
    prior(student_t(3, 0, 1), class = "sigma")
  ),
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 5000,
  warmup = 1500,
  backend = "cmdstanr",
  seed = 127,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "data/processed/fit_change_raw_revised",
  file_refit = "on_change"
)

fit_change_log <- brm(
  delta_log_egg_size_z ~ delta_log_egg_number_z,
  data = change_dat,
  family = gaussian(),
  prior = c(
    prior(normal(0, 0.7), class = "b"),
    prior(normal(0, 1), class = "Intercept"),
    prior(student_t(3, 0, 1), class = "sigma")
  ),
  chains = 4,
  cores = min(4, parallel::detectCores()),
  iter = 5000,
  warmup = 1500,
  backend = "cmdstanr",
  seed = 128,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "data/processed/fit_change_log_revised",
  file_refit = "on_change"
)

change_raw_draws <- as_draws_df(fit_change_raw)
change_log_draws <- as_draws_df(fit_change_log)

change_summary <- bind_rows(
  summarise_draw(change_raw_draws$b_delta_egg_number_z, "raw"),
  summarise_draw(change_log_draws$b_delta_log_egg_number_z, "log")
)

# =============================================================================
# Table 1: compact summary of focal posterior estimates
# =============================================================================

extract_table_rows <- function(summary_object, terms, analysis, parameter, scale) {
  idx <- match(terms, summary_object$term)

  tibble(
    Analysis = analysis,
    Parameter = parameter,
    Scale = scale,
    estimate = summary_object$estimate[idx],
    lower = summary_object$lower[idx],
    upper = summary_object$upper[idx],
    p_negative = summary_object$p_negative[idx]
  )
}

format_table_number <- function(x, digits = 2) {
  # Prevent very small estimates from printing as negative zero.
  threshold <- 0.5 * 10^(-digits)
  x[abs(x) < threshold] <- 0
  fmt(x, digits)
}

format_table_probability <- function(x) {
  case_when(
    x < 0.005 ~ "<0.01",
    x > 0.995 ~ ">0.99",
    TRUE ~ fmt(x, 2)
  )
}

table_1_data <- bind_rows(
  extract_table_rows(
    biv_cor_summary,
    terms = c("among_female", "within_female"),
    analysis = "Bivariate covariance",
    parameter = c(
      "Persistent among-female correlation",
      "Within-female residual correlation"
    ),
    scale = c("Proportional (log)", "Proportional (log)")
  ),
  extract_table_rows(
    context_slope_summary,
    terms = c("first", "second", "second_minus_first"),
    analysis = "Brood-specific slope",
    parameter = c(
      "First brood",
      "Second brood",
      "Second minus first"
    ),
    scale = c(
      "Proportional (pooled SD)",
      "Proportional (pooled SD)",
      "Proportional (pooled SD)"
    )
  ),
  extract_table_rows(
    context_scaling_sensitivity_summary,
    terms = c("first", "second", "second_minus_first"),
    analysis = "Brood-specific slope",
    parameter = c(
      "First brood",
      "Second brood",
      "Second minus first"
    ),
    scale = c(
      "Proportional (within-brood SD)",
      "Proportional (within-brood SD)",
      "Proportional (within-brood SD)"
    )
  ),
  extract_table_rows(
    change_summary,
    terms = c("raw", "log"),
    analysis = "Change score",
    parameter = c(
      "Within-female change",
      "Within-female change"
    ),
    scale = c("Absolute", "Proportional (log)")
  ),
  extract_table_rows(
    wb_summary,
    terms = c(
      "raw_within",
      "log_within",
      "raw_between",
      "log_between"
    ),
    analysis = "Within-between regression",
    parameter = c(
      "Within-female slope",
      "Within-female slope",
      "Among-female slope",
      "Among-female slope"
    ),
    scale = c(
      "Absolute",
      "Proportional (log)",
      "Absolute",
      "Proportional (log)"
    )
  )
) %>%
  mutate(
    Estimate = format_table_number(estimate),
    `95% CrI` = paste0(
      "[",
      format_table_number(lower),
      ", ",
      format_table_number(upper),
      "]"
    ),
    `P(<0)` = format_table_probability(p_negative)
  ) %>%
  select(Analysis, Parameter, Scale, Estimate, `95% CrI`, `P(<0)`)

table_1_flex <- table_1_data %>%
  flextable() %>%
  theme_booktabs() %>%
  merge_v(j = "Analysis") %>%
  valign(j = "Analysis", valign = "top", part = "body") %>%
  align(
    j = c("Estimate", "95% CrI", "P(<0)"),
    align = "center",
    part = "all"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9.5, part = "all") %>%
  padding(padding = 2, part = "all") %>%
  autofit() %>%
  set_table_properties(layout = "autofit", width = 1)

# =============================================================================
# Diagnostics
# =============================================================================

fits <- list(
  bivariate = fit_bivariate,
  context = fit_context,
  context_within_brood_scaling = fit_context_within_brood_scaling,
  within_between_raw = fit_wb_raw,
  within_between_log = fit_wb_log,
  change_raw = fit_change_raw,
  change_log = fit_change_log
)

model_diagnostics <- imap_dfr(fits, function(fit, model_name) {
  s <- posterior::summarise_draws(as_draws_array(fit))
  np <- brms::nuts_params(fit)

  tibble(
    model = model_name,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_bulk_ess = min(s$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(s$ess_tail, na.rm = TRUE),
    n_divergent = sum(
      np$Parameter == "divergent__" & np$Value == 1,
      na.rm = TRUE
    ),
    n_max_treedepth = sum(
      np$Parameter == "treedepth__" & np$Value >= 15,
      na.rm = TRUE
    )
  )
})

print(model_diagnostics)

# Posterior predictive checks for the two primary models. These are saved for
# inspection but need not be included in the manuscript unless a problem is
# detected.
dir.create(
  "figures/diagnostics",
  recursive = TRUE,
  showWarnings = FALSE
)

pp_biv_number <- pp_check(
  fit_bivariate,
  resp = "logeggnumberz",
  ndraws = 100
)
pp_biv_size <- pp_check(
  fit_bivariate,
  resp = "logeggsizez",
  ndraws = 100
)
pp_context <- pp_check(fit_context, ndraws = 100)

ggsave(
  "figures/diagnostics/pp_bivariate_egg_number.png",
  pp_biv_number,
  width = 6,
  height = 4,
  dpi = 300
)
ggsave(
  "figures/diagnostics/pp_bivariate_egg_size.png",
  pp_biv_size,
  width = 6,
  height = 4,
  dpi = 300
)
ggsave(
  "figures/diagnostics/pp_context.png",
  pp_context,
  width = 6,
  height = 4,
  dpi = 300
)

# ==========================================================
# Figure 1: brood-specific egg size-number relationships
# ==========================================================

# Construct a prediction sequence over the observed egg-number range within
# each brood. The model uses a common pooled standardization of log egg number.
prediction_grid <- df %>%
  filter(
    !is.na(egg.number),
    !is.na(log_egg_number_z),
    !is.na(log_pronotum_z)
  ) %>%
  group_by(brood, brood_label) %>%
  summarise(
    min_egg_number = min(egg.number, na.rm = TRUE),
    max_egg_number = max(egg.number, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    egg.number = list(
      seq(
        from = min_egg_number,
        to = max_egg_number,
        length.out = 120
      )
    )
  ) %>%
  ungroup() %>%
  unnest(egg.number) %>%
  mutate(
    log_egg_number = log(egg.number),
    log_egg_number_z =
      (log_egg_number - mean_log_egg_number) / sd_log_egg_number,
    # Predictions are conditional on average female structural size.
    log_pronotum_z = 0
  ) %>%
  select(
    brood,
    brood_label,
    egg.number,
    log_egg_number_z,
    log_pronotum_z
  )

context_predictions <- prediction_grid %>%
  add_epred_draws(
    fit_context,
    re_formula = NA,
    ndraws = 500
  ) %>%
  mutate(
    predicted_log_size =
      .epred * sd_log_egg_size + mean_log_egg_size,
    predicted_size = exp(predicted_log_size)
  ) %>%
  group_by(brood, brood_label, egg.number) %>%
  summarise(
    estimate = median(predicted_size),
    lower = quantile(predicted_size, 0.025),
    upper = quantile(predicted_size, 0.975),
    .groups = "drop"
  )

figure_1 <- ggplot(
  data = df,
  aes(x = egg.number, y = mean.egg.size)
) +
  geom_ribbon(
    data = context_predictions,
    aes(x = egg.number, ymin = lower, ymax = upper),
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = context_predictions,
    aes(x = egg.number, y = estimate),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(size = 2.2, alpha = 0.75) +
  facet_wrap(~ brood_label, scales = "fixed") +
  labs(
    x = "Egg number",
    y = "Mean egg perimeter (mm)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

figure_1

ggsave(
  filename = "figures/figure_1_context.png",
  plot = figure_1,
  width = 8.2,
  height = 4.2,
  units = "in",
  dpi = 400
)

# ==========================================================
# Figure 2: within-female changes between reproductive bouts
# ==========================================================

# Retain complete observations used in the figure
change_plot_dat <- change_dat %>%
  filter(
    !is.na(delta_egg_number),
    !is.na(delta_egg_number_z),
    !is.na(delta_egg_size)
  )

# Values needed to return standardized model predictions
# to the original measurement scales
mean_delta_number <- mean(
  change_plot_dat$delta_egg_number,
  na.rm = TRUE
)

sd_delta_number <- sd(
  change_plot_dat$delta_egg_number,
  na.rm = TRUE
)

mean_delta_size <- mean(
  change_plot_dat$delta_egg_size,
  na.rm = TRUE
)

sd_delta_size <- sd(
  change_plot_dat$delta_egg_size,
  na.rm = TRUE
)

# Construct a prediction sequence spanning the observed range
# of standardized changes in egg number
change_grid <- tibble(
  delta_egg_number_z = seq(
    from = min(
      change_plot_dat$delta_egg_number_z,
      na.rm = TRUE
    ),
    to = max(
      change_plot_dat$delta_egg_number_z,
      na.rm = TRUE
    ),
    length.out = 120
  )
)

# Generate posterior expected predictions
change_predictions <- change_grid %>%
  add_epred_draws(
    fit_change_raw,
    re_formula = NA,
    ndraws = 500
  ) %>%
  mutate(
    # Return the predictor to the original egg-number scale
    delta_number =
      delta_egg_number_z * sd_delta_number +
      mean_delta_number,
    
    # Return the predicted response to the original
    # egg-perimeter scale
    delta_size =
      .epred * sd_delta_size +
      mean_delta_size
  ) %>%
  group_by(
    delta_egg_number_z,
    delta_number
  ) %>%
  summarise(
    estimate = median(delta_size, na.rm = TRUE),
    lower = quantile(
      delta_size,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      delta_size,
      probs = 0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Construct Figure 2
figure_2 <- ggplot(
  data = change_plot_dat,
  aes(
    x = delta_egg_number,
    y = delta_egg_size
  )
) +
  geom_ribbon(
    data = change_predictions,
    aes(
      x = delta_number,
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_line(
    data = change_predictions,
    aes(
      x = delta_number,
      y = estimate
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    size = 2.2,
    alpha = 0.75
  ) +
  labs(
    x = expression(
      Delta * " egg number (second - first)"
    ),
    y = expression(
      Delta * " mean egg perimeter (mm; second - first)"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank()
  )

# Display Figure 2
figure_2

# Create the output directory if necessary
dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

# Save Figure 2
ggsave(
  filename = "figures/figure_2_change.png",
  plot = figure_2,
  width = 6.2,
  height = 5.2,
  units = "in",
  dpi = 400
)
# =============================================================================
# Save compact result tables used by the revised manuscript
# =============================================================================

write_csv(summary_stat, "data/processed/revised_summary_statistics.csv")
write_csv(change_direction_summary, "data/processed/revised_change_directions.csv")
write_csv(raw_cor_summary, "data/processed/revised_raw_correlations.csv")
write_csv(biv_cor_summary, "data/processed/revised_bivariate_correlations.csv")
write_csv(context_slope_summary, "data/processed/revised_context_slopes.csv")
write_csv(
  context_scaling_sensitivity_summary,
  "data/processed/revised_context_scaling_sensitivity.csv"
)
write_csv(
  second_brood_contrast_summary,
  "data/processed/revised_second_brood_original_units.csv"
)
write_csv(wb_summary, "data/processed/revised_within_between_sensitivity.csv")
write_csv(change_summary, "data/processed/revised_change_score_sensitivity.csv")
write_csv(model_diagnostics, "data/processed/revised_model_diagnostics.csv")
write_csv(table_1_data, "data/processed/revised_focal_results_table.csv")

cat("\n\n===== DESCRIPTIVE STATISTICS =====\n")
print(summary_stat)

cat("\n\n===== DIRECTION OF CHANGES BETWEEN BROODS =====\n")
print(change_direction_summary)

cat("\n\n===== RAW CORRELATIONS =====\n")
print(raw_cor_summary)

cat("\n\n===== BROOD-SPECIFIC SCALING VALUES =====\n")
print(brood_scaling)

cat("\n\n===== BIVARIATE MODEL: BROOD EFFECTS =====\n")
print(biv_brood_summary)

cat("\n\n===== BIVARIATE MODEL: CORRELATIONS =====\n")
print(biv_cor_summary)

cat("\n\n===== PRIMARY BROOD-SPECIFIC SLOPES: POOLED SCALING =====\n")
print(context_slope_summary)

cat("\n\n===== SENSITIVITY ANALYSIS: WITHIN-BROOD SCALING =====\n")
print(context_scaling_sensitivity_summary)

cat("\n\n===== SECOND-BROOD ORIGINAL-UNIT CONTRAST =====\n")
print(second_brood_contrast_summary)

cat("\n\n===== CHANGE-SCORE MODELS =====\n")
print(change_summary)

cat("\n\n===== WITHIN-BETWEEN MODELS =====\n")
print(wb_summary)

cat("\n\n===== MODEL DIAGNOSTICS =====\n")
print(model_diagnostics)

cat("\n\n===== FOCAL RESULTS TABLE =====\n")
print(table_1_data)
