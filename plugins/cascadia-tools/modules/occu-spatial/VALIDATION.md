# occu-spatial: validation

Status: **VERIFIED WITH LIMITS** for spPGOcc (default design), msPGOcc community means and
intPGOcc. **Not gate-passing:** sfMsPGOcc at fast settings (convergence), community SDs (bias),
spPGOcc at small or low-detection designs (bias), and the spatial range. spIntPGOcc is not validated.

All recovery runs use `fast_settings` (3 chains x 15000 iterations, burn-in 5000, thin 10, 10
neighbours, `converge_on = "core"`). They are gated with the unchanged `core/R/recovery.R` thresholds:
|bias| <= max(10% of truth, 3 MC SE); coverage >= 0.95 - 3 MC SE (0.885 at nsim 100); >= 90% of reps converged.

**Gate caveat (core bug, not fixed here):** `summarise_recovery()` sets `n_ok` to the number of
*converged* reps, so `n_converged / (n_ok + n_failed)` is always 1 unless a rep errored. The
">= 90% converged" clause therefore never fails for non-convergence. The convergence fractions
below were counted by hand from the reps files.

## 1. Builder's recovery runs (nsim 100 each; numbers re-read from results/occu-spatial/recover_*/)

Commands: `bash dev/occu-spatial/run_recover.sh <cfg>` (config in `dev/occu-spatial/<cfg>.yaml`).

### single_default: spPGOcc, 200 sites x 4 surveys, psi 0.5, p 0.4, sigma^2 1.5, range 6 km
100/100 converged, median 118 s/rep.

| param | truth | mean est | rel. bias | coverage | gate |
|---|---|---|---|---|---|
| occupancy_prob | 0.5 | 0.509 | +1.7% | 0.96 | PASS |
| psi_beta_forest | 0.8 | 0.819 | +2.4% | 0.98 | PASS |
| detection_prob | 0.4 | 0.406 | +1.5% | 0.97 | PASS |
| p_beta_effort | 0.5 | 0.517 | +3.3% | 0.96 | PASS |
| mean_occupancy | 0.5 | 0.503 | +0.9% | 0.94 | PASS |
| prop_sites_occupied | 0.504 | 0.503 | -0.1% | 0.94 | PASS |
| spatial_variance | 1.5 | 1.67 | +11.0% | 0.97 | PASS (descriptive only) |
| effective_range_km | 6 | 7.0 | +16.7% | 1.00 | FAIL bias (weakly identified; descriptive only) |

Power for psi_beta_forest is 0.96.

### single_hard: spPGOcc, 100 sites x 3 surveys, p 0.25
100/100 converged. occupancy_prob +13.3% (FAIL bias), mean_occupancy +11.0% (FAIL),
prop_sites_occupied +12.1% (FAIL), psi_beta_forest +9.7% (PASS, power 0.40), p_beta_effort +6.1% (PASS),
detection_prob -5.2% (PASS). Coverage is 0.95-0.98 for every estimand. The range (+10.8%) and
variance (-17.6%) fail on bias.

### strong_spatial: spPGOcc, sigma^2 3, range 10 km
97/100 converged. Covariates, occupancy and detection all PASS (psi_beta_forest -2.3%, coverage 0.97).

### negctrl_nonspatial: PGOcc fitted to the strong-spatial data (negative control)
psi_beta_forest -27.5%, coverage 0.74 (**FAIL**); occupancy_prob coverage 0.69 (**FAIL**).
The gate can fail when spatial autocorrelation is ignored.

### fpr_zero: spPGOcc with both effects = 0
psi_beta_forest rejection rate 0.07 and p_beta_effort 0.05 (100 reps).

### multi_ms: msPGOcc, 15 species, 150 sites
100/100 converged. Community means and mean_site_richness PASS (comm_psi_beta_forest +4.9%,
richness +0.5%). comm_sd_p_beta_effort is +16.7% (**FAIL bias**); comm_sd_detection_logit is +9.5% (PASS).

### multi_sf: sfMsPGOcc, 15 species, 150 sites, 1 factor (range 6 km)
**Only 76/100 converged** (R-hat < 1.05 on community-level estimands). That is below the 90% clause,
so the whole design FAILS. The core bug above hid this. Among the converged reps:
- PASS: community means, comm_sd_occupancy_logit, comm_sd_psi_beta_forest and richness;
- FAIL bias: comm_sd_detection_logit +11.3%, comm_sd_p_beta_effort +19.5%, effective_range_km_f1 +24.3%.

### integrated: intPGOcc, 150 point counts x 3 + 100 ARU x 5 (50 shared)
100/100 converged. All 8 estimands PASS (psi_beta_forest +6.7%, coverage 0.96).

### power over n_sites (50/100/200/400; nsim 50)
Power for psi_beta_forest is 0.34, 0.74, 0.90 and 1.00. At 50 sites p_beta_effort is biased +15.5%.

### Cross-checks
- **Simulator vs spOccupancy::simOcc** (400 datasets each, identical coordinates):
  - naive occupancy 0.4408 vs 0.4416 (KS p 0.98);
  - detection frequency 0.4628 vs 0.4636 (KS p 0.91);
  - Moran's I 0.0523 vs 0.0525 (KS p 0.64);
  - the covariance function is identical to `mkSpCov` (0.2030029 both).
- **intPGOcc vs hand-coded marginal ML** (48 datasets): mean |difference| 0.0008-0.004 for
  probabilities and detection slopes, 0.020 for psi_beta_forest; correlation 1.00.
- **msPGOcc vs a hand-written Nimble community model** (same priors, 3 datasets): mean |difference|
  is <= 0.008 for every community estimand and 0.006 for richness.
- The builder's chinook jobs for "spPGOcc vs PGOcc vs unmarked on non-spatial data" and for
  "spPGOcc fitted to simOcc data" **hit the SLURM time limit with no output**. The
  non-spatial engine comparison was rerun by the verifier on njobvu (section 2).
- **Nimble exact-GP vs spPGOcc:** still running on chinook at the time of verification (no result).

## 2. Independent verification (2026-10-19, separate verifier agent)

### What was checked
- **Literature and package source.** I checked spOccupancy 0.8.2 (installed source):
  - coefficient priors N(0, var 2.72) and sigma^2 ~ IG(2, 1) match the module;
  - phi ~ U(3/max d, 3/min d) is the package default, while the module uses the 5% distance quantile
    for the lower range bound (deliberate and documented);
  - the exponential covariance is sigma^2 exp(-phi d), so effective range = 3/phi (correlation 0.05);
  - sfMsPGOcc fixes factor variances at 1 and loadings lower-triangular with unit diagonal, and
    sim.R matches;
  - msPGOcc gives species coefficients N(mu, tau^2) with tau^2 ~ IG(0.1, 0.1), and the module reports sqrt(tau^2).
- **Shared code.** simulate() and fit() share only `d_or()` and the term-naming convention. My
  independent simulator shares nothing.

### Independent planted truth
`dev/occu-spatial/verify/vsim.R` is my own generator: eigen-decomposition GP draws, jittered-grid
sites in UTM metres, ragged 3-6 visits, raw-unit covariates. It writes researcher-style CSVs and
reads them back through the module's `read_data()` → `check_data()` → `fit()`, so io,
standardisation and estimand naming are also tested. Designs and parameters differ from the builder's.
Seeds are 7000 + rep. Results are in results/occu-spatial/verify/<scenario>/summary.csv.

| scenario | design | nsim | converged | result |
|---|---|---|---|---|
| sp | spPGOcc, 120 sites, p 0.3, psi 0.35, sigma^2 1, range 4 km, 2 covariates (elev partly spatial) | 150 | 150 | **FAIL**: psi_beta_elev -0.866 vs -0.70 (-24%), psi_beta_shrub 0.513 vs 0.40 (+28%), occupancy_prob +10.5%; p_beta_temp coverage 0.88 (bound 0.897). Coverage of both slopes 0.95-0.97. detection_prob, mean_occupancy, prop_occupied PASS |
| spbig | same generator, 480 sites | 100 | 95 | **all gated estimands PASS**: psi_beta_elev -2.8%, psi_beta_shrub +3.4%, occupancy_prob +3.8%, detection_prob -2.6%, p_beta_temp -2.2% (coverage 0.95), mean_occupancy +2.6%; coverage 0.93-0.98. spatial_variance -7% (PASS); effective_range_km +85% (FAIL, descriptive only) |
| fpr | sp with shrub effect 0 | 200 | 200 | false-positive rate for shrub **0.07** (14/200); psi_beta_elev -14.7% (FAIL bias) |
| neg | **PGOcc** on strongly spatial data (sigma^2 3, range 10 km) | 100 | 100 | **negative control fails, as it should**: psi_beta_elev coverage 0.86 (bound 0.885), psi_beta_shrub -20%, occupancy_prob coverage 0.63 |
| negsp | spPGOcc on the same data | 100 | 100 | shrub -0.9%, coverage 0.98 (PASS); elev (90% spatially smooth) -20%, coverage 0.90 (spatial confounding) |
| ms | msPGOcc, 10 species, 100 sites, ragged visits | 100 | 100 | community means and richness PASS (comm_psi_beta_conifer -0.4%); comm_sd_p_beta_effort_h +30%, comm_sd_psi_beta_conifer +18% (FAIL bias, coverage 0.95-0.96) |
| sf | sfMsPGOcc, same with 1 factor | 100 | 97 | community means and richness PASS; comm_sd_occupancy_logit +13%, comm_sd_p_beta +27%, comm_sd_psi_beta +16% (FAIL bias) |
| int | intPGOcc, transects 100 x 2 + cameras 100 x 6 (40 shared), lognormal covariate | 100 | 100 | all 8 PASS (psi_beta_wet +10.8% within 3 MC SE; coverage 0.93-0.97) |

**Re-running the builder's commands.** `recover_integrated` reproduced every number exactly
(results/occu-spatial/verify/rerun_integrated). `recover_single_default` rerun with the patched module:
identical table (100/100 converged; every mean, coverage and gate outcome matches section 1;
results/occu-spatial/verify/rerun_single_default). The `fit` on `example/config.yaml` with the
patched module also reproduced the builder's estimates exactly (results/occu-spatial/verify/example_fit).

**Interpretation of the small-design failure.** Same generator, parameters and seeds, but with 4x the sites (480 vs 120):
- psi_beta_elev bias went from -24% to -2.8%;
- psi_beta_shrub went from +28% to +3.4%;
- occupancy_prob went from +10.5% to +3.8%.

So the 120-site failure is small-sample behaviour of the posterior mean (logit slopes pulled away
from zero with about 40 occupied sites, p = 0.3 and a spatial random effect), not a module or io bug.
The builder's single_hard design (100 sites, p 0.25) shows the same pattern.

Coverage of the slopes stayed at 0.95-0.97 throughout, so the intervals remain honest; only the
point estimates lean high. This is recorded as a design limit in README Pitfalls and in the report's
"Validated limits" sentence. Users must run `recover` at their own design.

### Engine cross-check on non-spatial data (spec item; builder's job had timed out)
`dev/occu-spatial/verify/xc_nonspatial_v.R` uses 64 datasets from the default design with sigma^2 = 0. It fits
spPGOcc, PGOcc and `unmarked::occu` to each dataset
(results/occu-spatial/verify/xc_nonspatial_engines.csv; log logs/occu-spatial_verify_xc_nonspatial.log).

**PGOcc vs unmarked (independent likelihood engine):**
- per-dataset correlation >= 0.999 for all four parameters;
- mean |difference|: 0.002 detection_prob, 0.006 occupancy_prob, 0.002 p_beta_effort, 0.026 psi_beta_forest
  (the posterior mean is slightly larger than the MLE, as expected).

**spPGOcc vs PGOcc on the same non-spatial data:**
- detection agrees (mean |difference| 0.002-0.003);
- **psi_beta_forest is inflated:** relative bias +20.8% for spPGOcc vs +8.2% for PGOcc and +5.0% for
  unmarked; mean |sp - ns| = 0.10; coverage 0.92 vs 0.94.

The IG(2,1) prior keeps sigma^2 away from 0 (median posterior mean 0.66), so the slopes become
conditional on a spurious spatial effect. **Limit:** when the spatial variance posterior merely
restates its prior, also fit PGOcc and report both. The report now says this in that case.

### Usability (messy copies of the fisher example; `dev/occu-spatial/verify/messy.R`)
Before the fixes, these problems were handled silently or with unclear messages:
- NA detection values were silently treated as surveys not done;
- duplicated site x survey rows were silently accepted;
- **lat/long coordinates were silently accepted** (every distance would be ~0 km);
- a misspelled covariate column gave `object 'elevation_m' not found`.

All of these now give actionable messages. Missing site or survey covariates, a misspelled site id,
non-0/1 values and duplicate coordinates were already caught. The site-id error now suggests the
closest match ("'CAM-14' (did you mean 'CAM-014'?)").

Report check (fisher example): every number in the interpretation was recomputed from the
estimates table. Examples:
- OR per SD exp(1.010) = 2.75;
- elevation per 200 m exp(-0.551 x 200/246) = 0.64;
- detection per day exp(0.397/1.01) = 1.48;
- the IG(2,1) prior 95% interval is 0.18-4.13.

The report uses association language, keeps detection and habitat effects apart, and gives
effects on the probability scale. The conditional-effects curves (occupancy vs canopy; detection
vs days active) were inspected: they pass through the reported values at the covariate means and
include rug and 50/95% bands.

### Changes made by the verifier
- `R/io.R`:
  - missing-column errors list the file's actual columns and the config key to change;
  - formula terms are checked against the file;
  - empty detections are dropped with a note;
  - duplicate site x survey rows and duplicate site ids are flagged;
  - degree-like coordinates are flagged;
  - a misspelled site id gets a closest-match suggestion;
  - missing survey covariates report a count and example sites.
- `R/fit.R`: data notes are passed into the diagnostics. With `converge_on = "all"`, a note is added
  when R-hat >= 1.01 or bulk ESS < 400 (the real-data standard).
- `R/report.R`:
  - the per-unit odds ratio now uses exp(posterior-mean slope), consistent with the per-SD odds
    ratio (previously the mean of exp(), e.g. 2.26 vs 2.19);
  - the spatial paragraph now compares the spatial-variance interval with its prior interval and
    fixes the awkward "8.8 (95% CrI ...) km" wording;
  - trace panels have readable labels;
  - the "Validated limits" sentence now states the measured limits instead of claiming unbiased
    recovery at all designs.
- Estimation code is unchanged. Originals are kept in `dev/occu-spatial/verify/*.R.orig`.
- `R/report.R` (later): when the spatial-variance interval matches its prior, the report now
  recommends also fitting PGOcc (see the engine cross-check above).
- Note: these edits were made while the builder's `xc_nimble.R spatial` job (chinook, started
  23:07) was still running. That R process had already loaded the module and only uses `fit()`,
  which is unchanged, so its results are unaffected.

### Verifier verdict: VERIFIED WITH LIMITS
- **Pass (gated estimands: covariate effects, occupancy, detection, mean occupancy):**
  - spPGOcc at the default design (reproduced exactly) and at 480 sites on independent data;
  - msPGOcc community means and richness (builder's design and independent design);
  - intPGOcc (builder's design reproduced exactly; independent design passes).
- **Controls:**
  - the negative control (PGOcc on spatial data) fails the gate;
  - false-positive rates are 0.07 (independent, 200 reps), 0.07 and 0.05 (builder).
- **Limits (gate fails, documented):**
  - spPGOcc at about 100-120 sites with p 0.25-0.3: slope and baseline bias 10-28%, coverage OK;
  - spatial range is never gated;
  - community SDs +15-30%;
  - sfMsPGOcc reached only 76% convergence at fast settings (builder's design); the independent
    design reached 97%;
  - spPGOcc on non-spatial data inflates slopes about 10% relative to PGOcc;
  - spIntPGOcc is not validated.
- **Open item:** Nimble exact-GP vs spPGOcc result still pending (chinook job).
