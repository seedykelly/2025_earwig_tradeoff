# Data Integration Fixes

## Summary

Fixed multiple issues that prevented the ingested data from working with the reanalysis script:

### Issue 1: Brood Column Data Type Mismatch
**Problem**: Ingested data had `brood` as numeric (1, 2) but script expected character ("one", "two")
**Fix**: Updated script line 114-116 to use `brood_label` as the source for factor conversion
**Status**: ✓ Fixed

### Issue 2: Pronotum Length (mean.pro) Variation Within Females
**Problem**: Ingestion script was averaging egg `length` measurements per brood, which vary because they represent individual egg dimensions, not female structural traits
**Solution**: Integrated female body morphology data from `earwig_body_measurements.csv` 
**Result**: Now correctly has constant pronotum length per female across both broods

### Issue 3: ID Mismatches Between Data Files
Three different ID inconsistencies were found and corrected:

| Issue | Source | Target | Fix |
|-------|--------|--------|-----|
| Letter/Number | `ODF8T` (body file) | `0DF8T` (egg files) | Fixed in ingestion script |
| OCR Typo | `L39Q9` (body file) | `LE9Q9` (egg files) | Fixed in ingestion script |
| OCR Typo | `NZA18` (body file) | `NZAI8` (egg files) | Fixed in ingestion script |

**Status**: ✓ All corrected

### Issue 4: Missing Morphological Variables
**Problem**: Original df.rds had columns (pc1, mean.mass, mean.sci) not in raw egg files
**Solution**: 
- Added `mean.mass` (calculated from body measurements file)
- Left `pc1` and `mean.sci` as NA (require original analysis context)

**Status**: ✓ Resolved

## Files Updated

1. **scripts/earwig_tradeoff_reanalysis.R**
   - Line 96: Now checks for `df_complete_ingested.rds` first
   - Line 109: Added informative message showing which dataset is being loaded
   - Lines 114-116: Fixed brood factor conversion to use `brood_label`

2. **scripts/01_ingest_egg_data.R**
   - Integrated female body morphology data from `earwig_body_measurements.csv`
   - Added ID mismatch corrections (ODF8T, L39Q9, NZA18)
   - Proper calculation of `mean.pro` and `mean.mass` from body measurements

## Verification

✓ Data transformation pipeline runs successfully with all 45 females
✓ All three previously excluded females (0DF8T, LE9Q9, NZAI8) now included
✓ Pronotum length consistent within females across both broods
✓ Change score calculations work without errors
✓ All delta variables computed correctly

## Result

The reanalysis script [scripts/earwig_tradeoff_reanalysis.R](scripts/earwig_tradeoff_reanalysis.R) now successfully loads and analyzes all 45 female earwigs with both broods (90 observations total).
