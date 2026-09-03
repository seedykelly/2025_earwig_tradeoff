# Ingest raw egg measurement data and aggregate to female × brood level
# ==============================================================================
# This script reads individual egg measurements from the raw CSV files and
# produces a summary dataframe with one row per female per brood, containing:
#   - egg.number: count of eggs measured
#   - mean.egg.size: mean perimeter across eggs (average of perim, perim2, perim3)
#   - length: mean of the length column (interpreted as female structural size)

library(tidyverse)

# Create output directory
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Find all raw egg data files
raw_files <- list.files(
  "data/raw/egg_data_files",
  pattern = "\\.csv$",
  full.names = TRUE
)

cat("Found", length(raw_files), "raw egg data files\n\n")

# Read all files and combine
raw_data <- raw_files %>%
  map_df(read_csv, show_col_types = FALSE) %>%
  # Extract female ID and brood from filename if needed, but use file contents
  mutate(
    # Standardize brood to "one"/"two" format
    brood_label = ifelse(brood == 1, "one", "two"),
    # Average the three perimeter measurements per egg
    mean_perim = (perim + perim2 + perim3) / 3
  )

cat("Raw data dimensions:", nrow(raw_data), "eggs x", ncol(raw_data), "columns\n")
cat("Unique females:", n_distinct(raw_data$id), "\n")
cat("Brood values in raw data:", unique(raw_data$brood), "\n\n")

# Aggregate to female × brood level (eggs only)
df_ingested <- raw_data %>%
  group_by(id, brood, brood_label) %>%
  summarise(
    egg.number = n(),
    mean.egg.size = round(mean(mean_perim, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(id, brood)

# Load female body morphology measurements
# Note: These are separate from the egg measurement files
# ID mapping fixes for typos in body measurements file:
#   "ODF8T" (letter O) -> "0DF8T" (number 0)
#   "L39Q9" (typo) -> "LE9Q9"
#   "NZA18" (typo) -> "NZAI8"
body_measurements <- read_csv("data/raw/earwig_body_measurements.csv", show_col_types = FALSE) %>%
  mutate(
    id = case_when(
      id == "ODF8T" ~ "0DF8T",
      id == "L39Q9" ~ "LE9Q9",
      id == "NZA18" ~ "NZAI8",
      TRUE ~ id
    ),
    mean.pro = round((pronotum_length_1 + pronotum_length_2 + pronotum_length_3) / 3, 3),
    mean.mass = round((body_mass_1 + body_mass_2 + body_mass_3) / 3, 1)
  ) %>%
  select(id, mean.pro, mean.mass)

# Join body measurements to egg data (once per female, so use distinct)
df_ingested <- df_ingested %>%
  left_join(body_measurements, by = "id")

cat("Aggregated dataframe:\n")
print(df_ingested, n = Inf)

# Identify any data quality issues
cat("\n\n===== DATA QUALITY CHECK =====\n")
cat("Eggs per female × brood (range):\n")
print(
  df_ingested %>%
    summarise(
      min_eggs = min(egg.number),
      max_eggs = max(egg.number),
      mean_eggs = round(mean(egg.number), 1)
    )
)

cat("\nFemales with <5 eggs in any brood:\n")
print(
  df_ingested %>%
    filter(egg.number < 5) %>%
    select(id, brood_label, egg.number)
)

cat("\nFemales with missing broods:\n")
broods_by_id <- df_ingested %>%
  group_by(id) %>%
  summarise(n_broods = n(), .groups = "drop")

ids_with_both_broods <- broods_by_id %>%
  filter(n_broods == 2) %>%
  pull(id)

ids_missing_brood <- broods_by_id %>%
  filter(n_broods < 2) %>%
  pull(id)

if (length(ids_missing_brood) > 0) {
  print(
    df_ingested %>%
      filter(id %in% ids_missing_brood)
  )
} else {
  cat("(All females have both broods)\n")
}

cat("\nFemales in current df.rds:", length(unique(df_original$id)), "\n")
cat("Females in ingested data:", n_distinct(df_ingested$id), "\n")

# Compare
excluded_ids <- setdiff(unique(df_ingested$id), unique(df_original$id))
if (length(excluded_ids) > 0) {
  cat("\nFemales in raw data but NOT in df.rds:\n")
  print(
    df_ingested %>%
      filter(id %in% excluded_ids) %>%
      arrange(id, brood)
  )
}

# Add columns that exist in original df.rds but not in raw data
# pc1 requires PCA which needs all morphological variables; mean.sci is unclear
df_ingested <- df_ingested %>%
  mutate(
    pc1 = NA_real_,           # Principal component (computed from original analysis)
    mean.sci = NA_real_       # Reproductive output (source unclear)
  )

# Save the complete ingested dataframe
saveRDS(df_ingested, "data/processed/df_complete_ingested.rds")
cat("\n\nSaved complete ingested dataframe to data/processed/df_complete_ingested.rds\n")
cat("Note: pc1, mean.mass, and mean.sci columns added as NA (not in raw egg files)\n")
