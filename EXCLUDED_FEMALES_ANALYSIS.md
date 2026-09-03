# Analysis of Three Excluded Female Earwigs (Updated)

## Executive Summary

Three female earwigs (NZAI8, LE9Q9, 0DF8T) were present in the raw data files but absent from the original `df.rds` dataset (42 females). A complete data ingestion integrating raw egg measurements with female body morphology reveals:

**The three excluded females have complete, high-quality data with no apparent quality defects. Their exclusion from the original analysis appears to be due to ID encoding mismatches or data collection timing, not biological or measurement issues.**

## Investigation Results

### Data Quality Assessment

All three excluded females have:
- ✓ Complete measurements for both broods  
- ✓ Multiple eggs measured per brood (30–89 eggs)
- ✓ Three replicate perimeter measurements per egg
- ✓ Female pronotum length measurements (1.56–2.08 mm)
- ✓ Female body mass measurements (62.8–71.0 mg)
- ✓ No missing values or systematic measurement errors

### The Three Excluded Females

| Female | Brood | Eggs | Egg Size (mm) | Pronotum (mm) | Mass (mg) |
|--------|-------|------|---|---|---|
| 0DF8T | 1 | 30 | 3.58 | 1.93 | 64.4 |
| 0DF8T | 2 | 47 | 3.22 | 1.93 | 64.4 |
| LE9Q9 | 1 | 40 | 3.92 | 1.56 | 62.8 |
| LE9Q9 | 2 | 32 | 3.53 | 1.56 | 62.8 |
| NZAI8 | 1 | 89 | 3.44 | 2.08 | 71.0 |
| NZAI8 | 2 | 44 | 3.10 | 2.08 | 71.0 |

### Statistical Outlier Analysis

#### Egg Count Distribution

```
All 45 females: Mean = 43.1 ± 18.8 eggs per observation

Excluded females' counts: 30, 32, 40, 44, 47, 89
Z-scores:                  -0.72, -0.59, -0.17, 0.05, 0.21, 2.44
```

**Finding**: NZAI8's first brood (89 eggs) is 2.44 SD above the mean—moderately elevated but within acceptable range (not >3 SD). The high fecundity is biologically plausible and consistent with observed variation in earwig reproductive output.

#### Egg Size Distribution

```
All 45 females: Mean = 3.48 ± 0.22 mm

Excluded females' sizes: 3.10, 3.22, 3.44, 3.53, 3.58, 3.92
Z-scores:               -1.69, -1.15, -0.18, 0.25, 0.45, 2.00
```

**Finding**: LE9Q9's first brood (3.92 mm) is 2.00 SD above the mean—moderately large but again within acceptable bounds. The slightly larger eggs reflect individual variation in reproductive strategy.

### Comparison: Included vs. Excluded

| Statistic | Included (n=42) | Excluded (n=3) | Difference |
|-----------|---|---|---|
| **Mean eggs** | 42.8 ± 18.7 | 47.0 ± 21.6 | Not significant |
| **Mean size** | 3.48 ± 0.22 | 3.46 ± 0.29 | Not significant |
| **Mean pronotum** | — | 1.86 ± 0.26 | — |
| **Mean mass** | — | 65.4 ± 4.2 | — |

The three excluded females are **statistically indistinguishable** from the included 42 on key reproductive traits.

## Root Cause: ID Encoding Mismatches

The exclusion of these three females appears to stem from **systematic ID encoding problems** across data files rather than biological or measurement issues:

### 0DF8T
- **Egg files**: ID recorded as "0DF8T" (number zero)
- **Body file**: ID recorded as "ODF8T" (letter O)
- **Status**: Would not match in automated join, causing loss of morphological data

### LE9Q9
- **Egg files**: ID recorded as "LE9Q9" 
- **Body file**: ID recorded as "L39Q9" (OCR misread E as 3, Q as 9)
- **Status**: Would not match in automated join

### NZAI8
- **Egg files**: ID recorded as "NZAI8"
- **Body file**: ID recorded as "NZA18" (OCR misread I as 1)
- **Status**: Would not match in automated join

All three mismatches are **automatically corrected** in the updated ingestion script.

## Likely Reasons for Original Exclusion

Given that data quality is sound, the exclusion probably resulted from:

1. **ID mismatch prevention**: If the original analysis required successful morphological data joins, unmatched IDs would naturally be excluded
2. **Chronological collection difference**: May have been collected/measured in a separate batch with inconsistent ID coding
3. **Late inclusion**: May have been added to the study after the main analysis pipeline was established
4. **Manual curation**: Possible manual removal due to some criterion not visible in the raw data (e.g., incomplete collection notes, labeling issues)

## How the Three Were Recovered

The complete data ingestion process:

1. Read all 90 raw egg CSV files (accounting for all 45 females)
2. Discovered three ID mismatches when attempting to join with body measurements
3. Implemented automatic ID correction (O↔0, OCR typos)
4. Successfully integrated all 45 females with both broods and morphology

## Recommendation: Include or Sensitivity Analyze

**Option A: Analyze with all 45 females**
- Use `df_complete_ingested.rds` in the reanalysis script
- Update methods to state: "45 female earwigs" instead of "42"
- Document the ID corrections in supplementary materials

**Option B: Sensitivity analysis**
- Run analyses both with 42 and 45 females
- Compare results: effect sizes, credible intervals, statistical significance
- Report: "Results were/were not sensitive to inclusion of three additional females with corrected IDs"

**Option C: Investigate original exclusion criterion**
- Check lab notebooks, field notes, or protocol documents
- If a specific exclusion reason exists, document it and decide based on that criterion
- If no documentation exists, Option A or B is appropriate

## Technical Integration

The updated analysis pipeline now:

1. **Loads**: `data/processed/df_complete_ingested.rds` (45 females by default)
2. **Corrects IDs**: Automatically during ingestion (`scripts/01_ingest_egg_data.R`)
3. **Integrates morphology**: Joins female body measurements correctly despite ID mismatches
4. **Transforms data**: All downstream analyses in `scripts/earwig_tradeoff_reanalysis.R` work correctly

All transformation steps (scaling, standardization, within-between decomposition, change scores) execute without errors with the complete 45-female dataset.

## Conclusion

The three excluded females (NZAI8, LE9Q9, 0DF8T) represent **recovered data** resulting from ID encoding issues, not quality problems. They should be **included in the analysis** unless the original exclusion criterion can be specifically documented and justified.

The integration of egg measurements with female morphology is now complete and consistent across all 45 females.
