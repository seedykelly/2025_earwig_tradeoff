# Analysis of Three Excluded Female Earwigs

## Executive Summary

Three female earwigs (NZAI8, LE9Q9, 0DF8T) are present in the raw data files but were excluded from the original `df.rds` dataset used in the reanalysis. A complete ingestion of all 90 raw CSV files reveals:

**The three excluded females do NOT have obvious data quality problems and their exclusion reason is unknown.**

## Investigation Results

### Data Quality Assessment

All three excluded females have:
- ✓ Complete measurements for both broods  
- ✓ Multiple eggs measured per brood (30–89 eggs)
- ✓ Three replicate perimeter measurements per egg
- ✓ Reasonable structural size measurements (pronotum length)
- ✓ No systematic measurement errors or missing values

### Statistical Outlier Analysis

| Metric | Mean ± SD | Excluded Females | Z-scores |
|--------|-----------|------------------|----------|
| **Egg count** | 43.1 ± 18.8 | 30, 32, 40, 44, 47, 89 | −0.70 to 2.44 |
| **Egg size (mm)** | 3.48 ± 0.22 | 3.1, 3.22, 3.44, 3.53, 3.58, 3.92 | −1.69 to 2.00 |

**Finding**: NZAI8 brood 1 (89 eggs, Z = 2.44) is ~2.4 SD above the mean in egg count—moderately high but not extreme enough to warrant automatic exclusion. LE9Q9 brood 1 (3.92 mm eggs) is ~2.0 SD above the mean in size. Most other values are well within ±2 SD.

### Comparison with Included Females

| Group | N Females | Mean Eggs | SD Eggs | Mean Size | SD Size |
|-------|-----------|-----------|---------|-----------|---------|
| **Included (df.rds)** | 42 | 42.8 | 18.7 | 3.48 | 0.22 |
| **Excluded (new)** | 3 | 47.0 | 21.6 | 3.46 | 0.29 |

The excluded females are **not systematically different** from included females in their central tendency or variance.

## Data Files Generated

1. **`data/processed/df_complete_ingested.rds`**  
   Complete dataframe with all 45 females (90 observations, both broods)  
   Variables: `id`, `brood`, `brood_label`, `egg.number`, `mean.egg.size`, `mean.pro`

2. **`scripts/01_ingest_egg_data.R`**  
   Script that reads all 90 raw CSVs and aggregates to female × brood level  
   Includes data quality checks and comparison summaries

3. **`scripts/02_compare_datasets.R`**  
   Script to compare df.rds (42 females) with df_complete_ingested.rds (45 females)  
   Provides outlier detection and recommendations

4. **`data/processed/INGESTION_REPORT.md`**  
   Detailed report with egg counts, size measurements, and quality assessment

## Possible Reasons for Original Exclusion

Given the data quality appears sound, the exclusion may have been due to:

1. **Chronological or batch processing issue**: The three females may have been collected or processed at a different time
2. **Metadata or labeling discrepancy**: Possible issues in linking raw files to phenotypic/morphological data (not visible in egg CSVs alone)
3. **Subjective exclusion criterion**: A quality threshold applied during original analysis that isn't documented in the data itself
4. **Sample size target**: May have been excluded to reach a specific target sample size (42 is sometimes preferred for balanced designs)
5. **Late discovery or availability**: May have been added to the dataset after the original analysis was finalized

## Recommendations

### Option 1: Use Complete Dataset (45 females)
If you have no reason to exclude these three:

```r
# Update line 109 in earwig_tradeoff_reanalysis.R
df <- readRDS("data/processed/df_complete_ingested.rds") %>%
  mutate(
    id = factor(id),
    brood = factor(brood, levels = c(1, 2)),
    # ... rest of processing
  )
```

Then update your Methods: _"We analyzed data from **45 female European earwigs** producing **90 observations** (45 females × 2 broods)..."_

### Option 2: Run Sensitivity Analysis
Compare results with both datasets:

1. Save all results/figures from df.rds (42 females)
2. Rerun all analyses with df_complete (45 females)
3. Compare effect sizes, credible intervals, and conclusions
4. Document in supplementary materials: "Results are robust to inclusion/exclusion of three females originally excluded from analysis"

### Option 3: Investigate Original Reason
Before deciding, check:
- Lab notebooks or field notes for these three female IDs
- Original data entry sheets or morphological measurements
- Any protocol documents mentioning sample selection criteria
- Comments in version control history (if available)

## Conclusion

**The exclusion of NZAI8, LE9Q9, and 0DF8T cannot be justified on data quality grounds alone.** The decision to include or exclude them should be based on:

1. **Biological or experimental rationale** (not visible in egg measurement data)
2. **Original protocol documentation** (not present in current files)
3. **Sensitivity analysis** (recommended: show results are/are not sensitive to this choice)

All necessary code and data for working with the complete 45-female dataset are now available.
