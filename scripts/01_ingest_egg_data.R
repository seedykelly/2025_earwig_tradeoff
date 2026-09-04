# Ingest raw egg measurement data and aggregate to female × brood level
# ==============================================================================
# This script reads individual egg measurements from the raw CSV files and
# produces a summary dataframe with one row per female per brood, containing:
#   - egg.number: count of eggs measured
#   - mean.egg.size: mean perimeter across eggs (average of perim, perim2, perim3)
#   - mean.pro: mean female pronotum length from the body-measurement file
#   - mean.mass: mean female body mass from the body-measurement file

library(tidyverse)

# Create output directory
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Find all raw egg data files
raw_files <- list.files(
  "data/raw/egg_data_files",
  pattern = "\\.csv$",
  full.names = TRUE
)

expected_females <- 45L
expected_broods <- 2L
expected_files <- expected_females * expected_broods

if (length(raw_files) != expected_files) {
  stop(
    "Expected ", expected_files, " raw egg files (", expected_females,
    " females x ", expected_broods, " broods), but found ",
    length(raw_files), "."
  )
}

cat("Found", length(raw_files), "raw egg data files\n\n")

# Read all files and combine
raw_data <- raw_files %>%
  map_dfr(
    ~ read_csv(.x, show_col_types = FALSE) %>%
      mutate(source_file = basename(.x))
  ) %>%
  mutate(
    file_id = str_remove(source_file, "_[[:digit:]]{2}\\.csv$"),
    file_brood = as.integer(
      str_extract(source_file, "(?<=_)[[:digit:]]{2}(?=\\.csv$)")
    ),
    # Standardize brood to "one"/"two" format
    brood_label = ifelse(brood == 1, "one", "two"),
    # Average the three perimeter measurements per egg
    mean_perim = (perim + perim2 + perim3) / 3
  )

id_or_brood_mismatches <- raw_data %>%
  filter(id != file_id | brood != file_brood)

if (nrow(id_or_brood_mismatches) > 0) {
  print(
    id_or_brood_mismatches %>%
      distinct(source_file, id, brood, file_id, file_brood)
  )
  stop("Egg-file names disagree with IDs or brood numbers inside the files.")
}

if (
  !all(raw_data$brood %in% c(1L, 2L)) ||
  !all(raw_data$file_brood %in% c(1L, 2L))
) {
  stop("Brood codes must be exactly 1 (initial) or 2 (replacement).")
}

if (anyNA(raw_data[, c("num", "id", "brood", "perim", "perim2", "perim3")])) {
  stop("Missing egg number, ID, brood, or replicate perimeter measurement in raw egg data.")
}

duplicate_egg_rows <- raw_data %>%
  count(source_file, num, name = "n_rows") %>%
  filter(n_rows != 1L)

if (nrow(duplicate_egg_rows) > 0) {
  print(duplicate_egg_rows)
  warning(
    "Repeated egg-number labels were found within one or more raw files. ",
    "Egg number is calculated from measurement rows; confirm that each row ",
    "represents a distinct egg."
  )
}

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

duplicated_body_ids <- body_measurements %>%
  count(id) %>%
  filter(n != 1)

if (nrow(duplicated_body_ids) > 0) {
  print(duplicated_body_ids)
  stop("Body-measurement IDs are not unique after ID harmonization.")
}

# Join body measurements to egg data (once per female, so use distinct)
df_ingested <- df_ingested %>%
  left_join(body_measurements, by = "id", relationship = "many-to-one")

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

cat("Females in ingested data:", n_distinct(df_ingested$id), "\n")

if (file.exists("data/processed/df.rds")) {
  df_original <- readRDS("data/processed/df.rds")
  cat("Females in previous df.rds:", n_distinct(df_original$id), "\n")

  previously_unmatched_ids <- setdiff(
    unique(df_ingested$id),
    unique(df_original$id)
  )

  if (length(previously_unmatched_ids) > 0) {
    cat("\nFemales in raw data but not in previous df.rds:\n")
    print(
      df_ingested %>%
        filter(id %in% previously_unmatched_ids) %>%
        arrange(id, brood)
    )
  }
}

# Add columns that exist in original df.rds but not in raw data
# pc1 requires PCA which needs all morphological variables; mean.sci is unclear
df_ingested <- df_ingested %>%
  mutate(
    pc1 = NA_real_,           # Principal component (computed from original analysis)
    mean.sci = NA_real_       # Reproductive output (source unclear)
  )

validation_by_id <- df_ingested %>%
  count(id, brood, name = "n_rows") %>%
  group_by(id) %>%
  summarise(
    n_broods = n(),
    rows_per_brood_ok = all(n_rows == 1),
    .groups = "drop"
  )

if (
  nrow(df_ingested) != expected_files ||
  n_distinct(df_ingested$id) != expected_females ||
  any(validation_by_id$n_broods != expected_broods) ||
  !all(validation_by_id$rows_per_brood_ok) ||
  anyNA(df_ingested[, c("egg.number", "mean.egg.size", "mean.pro")]) ||
  any(df_ingested$egg.number <= 0) ||
  any(df_ingested$mean.egg.size <= 0) ||
  any(df_ingested$mean.pro <= 0)
) {
  stop("Final dataset failed structural or completeness validation.")
}

cat(
  "Validation passed:", n_distinct(df_ingested$id),
  "females,", nrow(df_ingested),
  "female x brood observations, complete focal variables.\n"
)

# Save the complete ingested dataframe
saveRDS(df_ingested, "data/processed/df_complete_ingested.rds")
cat("\n\nSaved complete ingested dataframe to data/processed/df_complete_ingested.rds\n")
cat("Note: pc1 and mean.sci are retained as NA for compatibility; mean.mass comes from the body-measurement file.\n")
