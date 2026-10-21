---
name: eco-data-intake
description: Inspect and validate a researcher's raw ecological data (detection histories, counts, telemetry, capture histories, genotypes, length-age data) before analysis. Use whenever new data files arrive.
---

# eco-data-intake — look before you model

1. **Never modify the researcher's original files.** Copy to
   `<project>/data_raw/`; write cleaned versions to `<project>/data/` with
   a script (`clean.R`) so every change is reproducible.
2. **Describe the data back** to the researcher in one short paragraph:
   rows, columns, sites, surveys, date range, number of detections, NAs.
   Ask them to confirm. Most analysis errors are misunderstandings here.
3. **Run the module's `check_data()`** and fix what it reports with the
   researcher (never silently drop rows).
4. Common traps:
   - Detection histories: NA (not surveyed) vs 0 (surveyed, not detected).
   - Counts: zero-inflation from unsurveyed sites coded as 0.
   - Covariates: units, factor levels misspelled ("Forest" vs "forest"),
     site IDs that don't match between files, covariates measured after
     the response (post-treatment).
   - Telemetry: duplicated timestamps, time zones, projected vs lon/lat
     (steps must be in metres), irregular fix intervals, fix failures.
   - Capture-recapture: individuals with identical IDs across sessions,
     trap operation (effort) matrix missing.
   - Genotypes: allele coding, missing data per locus/individual (drop
     loci with > 20% missing, individuals > 30%, with the researcher's
     agreement), duplicate samples (same individual genotyped twice), HWE
     and linkage checks before structure/Ne analysis.
   - Fisheries: length units (mm vs cm), ages from different readers,
     gear selectivity (catch curves use only fully recruited ages).
5. Record the final dimensions; they set the design for eco-planted-truth.
