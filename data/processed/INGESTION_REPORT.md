# Egg Data Ingestion Report

## Summary

Successfully ingested all 90 raw egg measurement CSV files from `data/raw/egg_data_files/`, aggregating individual egg measurements to female × brood level.

### Key Findings

| Metric | Value |
|--------|-------|
| **Total females** | 45 |
| **Females in original df.rds** | 42 |
| **Females in ingested data** | 45 |
| **Excluded females** | 3 (NZAI8, LE9Q9, 0DF8T) |

## Excluded Females

### Data Summary

| Female ID | Brood | Egg Count | Mean Perimeter (mm) | Pronotum Length (mm) |
|-----------|-------|-----------|---------------------|----------------------|
| **NZAI8** | 1 (First) | 89 | 3.44 | 1.15 |
| **NZAI8** | 2 (Second) | 44 | 3.10 | 0.97 |
| **LE9Q9** | 1 (First) | 40 | 3.92 | 1.35 |
| **LE9Q9** | 2 (Second) | 32 | 3.53 | 1.18 |
| **0DF8T** | 1 (First) | 30 | 3.58 | 1.20 |
| **0DF8T** | 2 (Second) | 47 | 3.22 | 1.09 |

### Quality Assessment

#### NZAI8
- **Observation**: Notably HIGH first-brood egg count (89 eggs)
- **Assessment**: Data quality appears sound. No missing values, complete measurements, and all replicate perimeters recorded. Egg size (3.44) is within normal range.
- **Note**: High fecundity in first brood is biologically plausible and not necessarily problematic.

#### LE9Q9
- **Observation**: Slightly LARGER mean egg size (3.92 mm) compared to sample average (3.48 mm)
- **Assessment**: Data quality appears sound. Complete measurements in both broods, normal egg counts.

#### 0DF8T
- **Observation**: Lower first-brood egg count (30 eggs) but compensated with high second-brood count (47 eggs)
- **Assessment**: Data quality appears sound. Well-measured eggs, normal sizes, complete data.

## Comparison: Included vs. Excluded Females

| Group | Min Eggs | Max Eggs | Mean Eggs | SD Eggs | Mean Size |
|-------|----------|----------|-----------|---------|-----------|
| Included (n=42) | 1 | 86 | 42.8 | 18.7 | 3.48 |
| Excluded (n=3) | 30 | 89 | 47.0 | 21.6 | 3.46 |

## Conclusion

**The three excluded females do NOT show obvious data quality problems.** They have:
- Complete measurements for both broods
- Reasonable egg counts (well within biological range)
- Normal egg sizes (no extreme outliers)
- All three perimeter replicates recorded
- No patterns suggesting systematic measurement error

**Likely reasons for original exclusion** (speculative, as not documented in original code):
1. **Chronological exclusion**: May have been collected or processed at a different time than the 42 included females
2. **Study design decision**: Possibly excluded to reach a target sample size
3. **Undocumented quality criterion**: Some criterion beyond raw data quality (e.g., handling time, labeling issues, incomplete metadata in original field notes)
4. **Late addition**: May have been collected after the main analysis was underway

## Files Generated

- `data/processed/df_complete_ingested.rds` - Complete dataframe with all 45 females (both broods)
- `scripts/01_ingest_egg_data.R` - Ingestion script that reads raw CSVs and aggregates to summary level

## Next Steps

To include the three excluded females in your analysis:

1. **Load the complete dataframe**:
   ```r
   df_complete <- readRDS("data/processed/df_complete_ingested.rds")
   ```

2. **Review original exclusion reasoning** - Check lab notebooks, field notes, or protocol documents to understand why these three were originally excluded

3. **Decide on inclusion**: 
   - If the original exclusion reason no longer applies, you can use `df_complete`
   - If you want sensitivity analysis, run models both with and without these three females
   - Document your decision in the methods/supplementary materials

