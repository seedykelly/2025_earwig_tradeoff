# Egg Data Ingestion Report (Updated)

## Summary

Successfully ingested all 90 raw egg measurement CSV files from `data/raw/egg_data_files/`, aggregating individual egg measurements to female × brood level. Integrated female body morphology data from `earwig_body_measurements.csv`.

### Key Findings

| Metric | Value |
|--------|-------|
| **Total females** | 45 |
| **Observations** | 90 (45 females × 2 broods) |
| **Complete broods** | All females have both broods |
| **Eggs per female per brood** | Range: 1–89; Mean: 43.1 |

## Data Sources

### Primary: Egg Measurements
- **Location**: `data/raw/egg_data_files/*.csv` (90 files)
- **Structure**: One row per egg per female per brood
- **Measured variables**: 
  - Egg perimeter (3 replicate measurements)
  - Individual egg length
  - Brood assignment (1 = first, 2 = second)

### Secondary: Female Morphology
- **Location**: `data/raw/earwig_body_measurements.csv`
- **Structure**: One row per female
- **Measured variables**:
  - Pronotum length (3 replicate measurements)
  - Body mass (3 replicate measurements)

### Aggregation
```
Individual eggs (per female × brood)
       ↓
       Aggregate by female × brood
       ├─ Count eggs (egg.number)
       ├─ Mean perimeter (mean.egg.size)
       └─ Join female body measurements
           ├─ Mean pronotum length (mean.pro)
           └─ Mean body mass (mean.mass)
```

## Three Previously Excluded Females

Now included in complete dataset with full morphological data:

### 0DF8T (ID mismatch: "ODF8T" in body file—letter O vs number 0)
- **Brood 1**: 30 eggs, 3.58 mm mean size
- **Brood 2**: 47 eggs, 3.22 mm mean size
- **Pronotum**: 1.93 mm | **Mass**: 64.4 mg
- **Status**: Complete data, no quality issues

### LE9Q9 (ID mismatch: "L39Q9" in body file—OCR typo)
- **Brood 1**: 40 eggs, 3.92 mm mean size (notably large)
- **Brood 2**: 32 eggs, 3.53 mm mean size
- **Pronotum**: 1.56 mm | **Mass**: 62.8 mg
- **Status**: Complete data, slightly larger eggs than average

### NZAI8 (ID mismatch: "NZA18" in body file—OCR typo)
- **Brood 1**: 89 eggs, 3.44 mm mean size (notably high fecundity)
- **Brood 2**: 44 eggs, 3.10 mm mean size
- **Pronotum**: 2.08 mm | **Mass**: 71.0 mg (larger female)
- **Status**: Complete data, high first-brood fecundity is biologically plausible

## Data Quality: Comparison Across Groups

| Group | N | Mean Eggs ± SD | Mean Size ± SD | Mean Pronotum ± SD |
|-------|---|---|---|---|
| **All 45 females** | 45 | 43.1 ± 18.8 | 3.48 ± 0.22 | 1.86 ± 0.18 |
| **Original 42** | 42 | 42.8 ± 18.7 | 3.48 ± 0.22 | — |
| **Three excluded** | 3 | 47.0 ± 21.6 | 3.46 ± 0.29 | 1.86 ± 0.26 |

The excluded females are **not systematically different** from the original 42. Their inclusion slightly increases the standard deviation in egg counts but does not meaningfully alter the sample characteristics.

## Data Completeness

All 45 females have:
- ✓ Both broods (90 observations total)
- ✓ Egg counts and mean sizes per brood
- ✓ Female pronotum length measurements
- ✓ Female body mass measurements
- ✓ Corrected ID mappings

## ID Corrections Applied

The ingestion script automatically corrects three ID mismatches between data files:

| Issue | Body File | Egg Files | Type |
|-------|-----------|-----------|------|
| 0DF8T | ODF8T | 0DF8T | Letter O vs number 0 |
| LE9Q9 | L39Q9 | LE9Q9 | OCR reading error |
| NZAI8 | NZA18 | NZAI8 | OCR reading error |

## Technical Notes

### Pronotum Length
- Stored as **constant per female across both broods** (as it should be—it's a structural trait measured post-mortem)
- Calculated as mean of three replicate measurements from body measurements file
- Varies across females (range: 1.43–2.15 mm) reflecting individual size variation

### Body Mass
- Similarly constant within females across broods
- Calculated as mean of three replicate measurements
- Ranges from 36.1–82.2 mg, reflecting substantial variation in female condition

### Egg Size
- Varies by brood within females, reflecting reproductive trade-off
- Mean across both broods and all females: 3.48 ± 0.22 mm
- Calculated as mean of three replicate perimeter measurements per egg

## Files Generated/Updated

- `data/processed/df_complete_ingested.rds` - Complete dataframe (45 females, 90 observations)
- `scripts/01_ingest_egg_data.R` - Ingestion script with ID corrections and morphology integration
- `scripts/02_compare_datasets.R` - Comparison and sensitivity analysis helper
- `scripts/earwig_tradeoff_reanalysis.R` - Analysis script (updated to load complete data)

## Conclusion

All 45 females have complete, high-quality egg and morphological measurements. The three previously excluded females (0DF8T, LE9Q9, NZAI8) show no data quality issues and are now fully integrated into the complete dataset ready for analysis.
