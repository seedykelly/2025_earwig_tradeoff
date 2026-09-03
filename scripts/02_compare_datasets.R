# Compare original df.rds (42 females) with complete ingested data (45 females)
# This script helps you decide whether to include the three excluded females

library(tidyverse)
library(tidybayes)

# Load both versions
df_original <- readRDS("data/processed/df.rds")
df_complete <- readRDS("data/processed/df_complete_ingested.rds")

# Identify excluded females
excluded_ids <- setdiff(unique(df_complete$id), unique(df_original$id))

cat("===== DATASET COMPARISON =====\n")
cat("Original dataset (df.rds):", n_distinct(df_original$id), "females,", 
    nrow(df_original), "observations\n")
cat("Complete dataset (ingested):", n_distinct(df_complete$id), "females,", 
    nrow(df_complete), "observations\n")
cat("Excluded females:", paste(excluded_ids, collapse = ", "), "\n\n")

# Descriptive comparison
cat("===== SAMPLE CHARACTERISTICS =====\n")
cat("Included females (n=42):\n")
print(
  df_complete %>%
    filter(!id %in% excluded_ids) %>%
    summarise(
      Mean_eggs = round(mean(egg.number), 1),
      SD_eggs = round(sd(egg.number), 1),
      Mean_size = round(mean(mean.egg.size), 2),
      SD_size = round(sd(mean.egg.size), 2),
      Mean_pronotum = round(mean(mean.pro), 2),
      SD_pronotum = round(sd(mean.pro), 2)
    )
)

cat("\nExcluded females (n=3):\n")
print(
  df_complete %>%
    filter(id %in% excluded_ids) %>%
    summarise(
      Mean_eggs = round(mean(egg.number), 1),
      SD_eggs = round(sd(egg.number), 1),
      Mean_size = round(mean(mean.egg.size), 2),
      SD_size = round(sd(mean.egg.size), 2),
      Mean_pronotum = round(mean(mean.pro), 2),
      SD_pronotum = round(sd(mean.pro), 2)
    )
)

cat("\n===== OUTLIER ANALYSIS =====\n")
cat("Are excluded females statistical outliers on key variables?\n\n")

# Egg number outlier detection
egg_mean <- mean(df_complete$egg.number)
egg_sd <- sd(df_complete$egg.number)
excluded_eggs <- df_complete %>%
  filter(id %in% excluded_ids) %>%
  pull(egg.number)

cat("Egg number distribution:\n")
cat("  Mean:", round(egg_mean, 1), "± SD", round(egg_sd, 1), "\n")
cat("  Excluded females' values:", paste(excluded_eggs, collapse = ", "), "\n")
cat("  Z-scores:", paste(round((excluded_eggs - egg_mean) / egg_sd, 2), collapse = ", "), "\n")

# Egg size outlier detection
size_mean <- mean(df_complete$mean.egg.size)
size_sd <- sd(df_complete$mean.egg.size)
excluded_sizes <- df_complete %>%
  filter(id %in% excluded_ids) %>%
  pull(mean.egg.size)

cat("\nEgg size distribution:\n")
cat("  Mean:", round(size_mean, 2), "± SD", round(size_sd, 2), "\n")
cat("  Excluded females' mean sizes:", paste(round(excluded_sizes, 2), collapse = ", "), "\n")
cat("  Z-scores:", paste(round((excluded_sizes - size_mean) / size_sd, 2), collapse = ", "), "\n")

cat("\n===== RECOMMENDATION =====\n")
cat("If you want to include the three excluded females:\n")
cat("  1. Update line 109 of earwig_tradeoff_reanalysis.R to load df_complete:\n")
cat("     df <- readRDS('data/processed/df_complete_ingested.rds') %>% ...\n")
cat("  2. Note this change in your methods: 'We analyzed data from 45 females'\n")
cat("  3. Consider noting in supplementary materials why these three were initially excluded\n\n")

cat("To run a sensitivity analysis (both with and without):\n")
cat("  1. Save original results from df.rds (42 females)\n")
cat("  2. Re-run analysis with df_complete (45 females)\n")
cat("  3. Compare effect sizes, confidence intervals, and conclusions\n")
