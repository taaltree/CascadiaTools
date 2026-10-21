---
id: movement-ssf
title: Habitat selection from telemetry - integrated step selection (iSSA), RSF and home ranges
when_to_use: GPS (or other regular-interval) telemetry from several animals plus habitat rasters; the question is which habitats animals select or avoid while moving, how strongly (relative selection strength), how much individuals differ, what their movement kernel looks like (step lengths, turning) and how large their home ranges are.
not_for: Behavioural states / switching between encamped and travelling (use movement-hmm); very irregular or sparse fixes (> ~12 h apart, VHF triangulation) where steps lose meaning (use an RSF/home-range analysis alone); presence-only or survey data without individual tracks (use sdm-presence-only); categorical land-cover rasters (convert to proportion-of-class rasters, e.g. % forest within 100 m, first); migration or dispersal tracks that are not range-resident (iSSA still works with home_attraction false and home_range false, but home ranges are meaningless).
data_needed: GPS fix CSV (animal id, timestamp, projected x/y in metres, failed fixes allowed as empty rows) and one GeoTIFF per continuous habitat covariate in the same projected CRS
engine: amt (track regularisation, steps, available steps), glmmTMB (Poisson formulation of the mixed conditional logistic model, Muff et al. 2020), ctmm (AKDE), amt (KDE, MCP)
---

# Habitat selection from telemetry (iSSA)

## Question it answers

"Do my collared elk select forest and avoid high elevation, how strongly,
and do individuals differ?" An integrated step-selection analysis (iSSA;
Avgar et al. 2016) compares each observed step (fix to fix) with steps the
animal *could* have taken from the same start point, so availability changes
with where the animal is. It returns:

| Output | Meaning |
|---|---|
| `sel_<cov>_mean` | population-mean selection coefficient, log relative selection strength per 1 SD of the covariate (link scale; power is computed for this) |
| `sel_<cov>_sd` | among-animal SD of that coefficient (random slope SD) |
| `rss_<cov>_per_sd`, `rss_<cov>_per_unit` | relative selection strength (Avgar et al. 2017): how many times more likely a typical animal is to choose a location 1 SD (or 1 meaningful unit, e.g. 10% forest) higher in the covariate than an equally reachable one |
| `step_shape`, `step_scale_m`, `step_mean_m` | selection-free gamma step-length distribution (what movement would look like in homogeneous habitat), updated from the tentative kernel |
| `turn_kappa` | selection-free von Mises concentration of turning angles (0 = random turning, > 0 = directional persistence) |
| `home_attraction_per_km2` | range-residency nuisance term (on by default; overestimated, do not interpret its size; see Pitfalls) |
| `rsf_<cov>_mean`, `rsf_<cov>_sd` | use-availability RSF within each animal's MCP, for comparison only |
| `hr_akde_km2_<id>` | 95% AKDE home-range area per animal with CI (KDE and MCP areas in the figures and fit object) |

## Model

For animal *i* at step *t* (start point fixed, end point *x*), the iSSF is

    f(x | start, previous heading) ∝ phi(sl, ta; theta) * exp( sum_k beta_ik * z_k(x) )

* `z_k(x)` covariate *k* at the step end, standardised by the mean and SD of
  the raster across the study landscape (so beta is per 1 SD of the
  landscape; per-unit RSS is `exp(beta * unit / SD)`).
* `phi` is the selection-free movement kernel: step length ~ Gamma(shape,
  scale), turn angle ~ von Mises(0, kappa).
* Available steps are drawn from a *tentative* gamma / von Mises kernel fitted
  to all observed steps (pooled across animals, so the update is a
  population kernel). Including `sl`, `log(sl)` and `cos(ta)` as covariates
  corrects the tentative kernel (Avgar et al. 2016; Fieberg et al. 2021):
  shape = shape_tent + beta_log_sl; scale = 1 / (1/scale_tent - beta_sl);
  kappa = kappa_tent + beta_cos_ta (identical to `amt::update_gamma` and
  `amt::update_vonmises`). Intervals by the delta method (log scale for
  shape, scale and mean step length).
* Individual variation: beta_ik = mu_k + u_ik, u_ik ~ Normal(0, sd_k^2),
  independent across covariates (Muff et al. 2020 also assume independent
  slopes). The conditional logistic likelihood is fitted as a Poisson GLMM
  with one intercept per step (stratum) whose variance is fixed at 10^6, which
  is equivalent to conditioning on the stratum (Muff et al. 2020), in
  glmmTMB (Laplace approximation). Variance components by REML (default
  `reml: true`): with 12 animals ML gave the among-animal SDs a relative bias
  of -15% to -33%, REML -2% to -3% at the default design (VALIDATION.md).
* Relative selection strength (Avgar et al. 2017):
  RSS(x1, x2) = exp(beta * (z(x1) - z(x2))) for two locations equally
  accessible from the start point; `rss_<cov>` figures plot RSS of each
  covariate value against a location at the landscape mean.
* RSF (comparison): use (fixes, weight 1) vs available points (10 per fix,
  uniform in each animal's 100% MCP, weight 5000; Fieberg et al. 2021),
  binomial glmmTMB with a random intercept of fixed large variance and random
  slopes (Muff et al. 2020). It ignores movement and treats fixes as
  independent.
* Home ranges: per animal, ctmm fits continuous-time movement models (IID, OU,
  OUF, isotropic/anisotropic), selects by AICc, and computes the
  autocorrelated KDE with CI on area (Fleming et al. 2015; Calabrese et al.
  2016). KDE (Gaussian reference bandwidth) and MCP from amt are shown for
  comparison.
* Model check: used-habitat calibration (UHC; Fieberg et al. 2018): 20% of
  steps are held out, a population clogit is fitted to the rest, and the
  model's predicted distribution of each covariate at chosen steps (200
  parameter draws) is compared with the observed held-out used steps.

## Assumptions (and what breaks if they fail)

* **Fix interval is regular** (after `track_resample` with a tolerance). Steps
  spanning missed fixes are dropped (bursts). Mixed schedules: give
  `interval_hours` for the one you want; others are discarded.
* **The covariate at the step end drives selection.** If animals respond to
  what they cross along the step, or to habitat at a larger scale, the
  coefficient is a proxy at this scale.
* **Movement kernel shape is gamma x von Mises and the same everywhere.**
  Habitat-dependent speed (e.g. faster in open habitat) is not modelled here;
  it would partly load onto the selection coefficients.
* **Independent random slopes, normal among animals.** Needs several animals;
  with 4-8 animals the SD was biased low by 20-26% and its interval can be degenerate (see
  VALIDATION.md).
* **GPS location error is small relative to step length** and fix success is
  not habitat-dependent (habitat-biased fix loss under dense canopy biases
  selection for canopy downward).
* **Associations, not causes.** Covariates correlated with unmapped features
  (forage, predators, roads) carry their effects.

## Data format

`fixes` CSV, one row per fix attempt (extra columns are ignored):

| column (config key) | content |
|---|---|
| `id` (default `animal_id`) | animal / collar deployment id |
| `time` (default `timestamp`) | date-time; set `time_format` and `timezone` if not `YYYY-mm-dd HH:MM:SS` UTC |
| `x`, `y` | projected coordinates in metres (UTM etc.); rows with empty x/y are failed fixes and are dropped. Longitude/latitude: set `lonlat: true` and they are projected to the raster CRS |

```
animal_id,timestamp,utm_e,utm_n,dop,fix_type
E01,2024-05-01 00:00:23,511014.9,5007670.1,4,3D
E01,2024-05-01 02:00:16,511255.5,5007538.4,2.3,3D
E01,2024-05-01 04:01:33,,,,no fix
```

Rasters: one GeoTIFF per continuous covariate (`rasters: {name: file.tif}`)
in the same projected CRS as the fixes, covering the animals' ranges plus at
least a few step lengths of margin. Rasters with different grids are
resampled to the first one. Rasters are standardised with their own mean and
SD over the whole file, then cropped to the fixes + `buffer_m` (default
5000 m).

## Config

```yaml
data:
  fixes: elk_gps_2024.csv
  id: animal_id
  time: timestamp
  timezone: UTC
  time_format: "%Y-%m-%d %H:%M:%S"
  x: utm_e
  "y": utm_n               # quote "y": unquoted, YAML reads it as TRUE (handled, but quoting is safer)
  crs: "EPSG:32611"          # default: the rasters' CRS
  lonlat: false              # true if x/y are longitude/latitude
  rasters: {forest: forest_cover.tif, elevation: elevation.tif}
  covariates: [forest, elevation]   # default: all rasters
  labels: {forest: "Forest cover (%)", elevation: "Elevation (m)"}   # axis labels
  names: {forest: "forest cover", elevation: "elevation"}            # words in sentences
  units: {forest: 10, elevation: 100}          # meaningful change for per-unit RSS
  unit_text: {forest: "10 percentage points of forest cover", elevation: "100 m of elevation"}
  animal_label: elk          # animal_plural: if the plural is not label + "s"
  buffer_m: 5000
settings:
  n_control: 20              # available steps per used step (>= 10; more is better)
  interval_hours: 2          # default: median fix interval
  tolerance_minutes: 12      # default: 10% of the interval
  random_slopes: true        # animal-level random slopes for every covariate
  reml: true                 # REML for the among-animal SDs (ML biases them low)
  home_attraction: true      # squared distance (km^2) to an estimated centre per animal (range residency)
  home_centre: estimated     # or 'mean' (plug-in mean location; not recommended)
  min_steps: 20              # animals with fewer usable steps are left out of the iSSA
  rsf: true                  # also fit the RSF
  rsf_n_avail: 10
  rsf_weight: 5000
  home_range: true           # AKDE + KDE + MCP
  hr_level: 0.95
  uhc: true                  # used-habitat calibration check
  uhc_draws: 200
  seed: 1
  method: issa               # "rsf" = RSF only, reported under the sel_* names (validation negative control)
figures: {...}               # see docs/FIGURES.md
```

## Outputs

`<id>_estimates.csv` (all estimands above), `<id>_report.html` (methods,
diagnostics including the UHC check, plain-language interpretation, figures),
`figures/` with `rss_<cov>`, `selection_coefficients`, `movement_kernel`,
`used_available`, `uhc`, `home_ranges`, `hr_area`. `fit$issa$indiv` holds
animal-level coefficients, `fit$home_ranges$table` the AKDE/KDE/MCP areas and
effective sample sizes.

## Figures

`core/run.R fit` exports these to `figures/` (PDF + PNG, gallery.html with
captions); `core/run.R figures movement-ssf --fit <dir>/movement-ssf_fit.rds
--config <cfg>` redraws them from a saved fit after editing the `figures:` block.

| id | what it shows |
|---|---|
| `rss_<cov>` | Relative selection strength (Avgar et al. 2017) vs the covariate in original units, relative to the landscape mean, log y scale (so the line is log-RSS); population-level estimate with `ci` confidence band, dashed RSS = 1, rug of used step ends |
| `selection_coefficients` | Forest plot of log-RSS per 1 SD: population mean (diamond, `ci` interval) and each animal's coefficient (mean + random slope), one panel per covariate |
| `movement_kernel` | A: step length, B: turn angle; observed-step histograms with the tentative (dashed) and selection-free updated (solid) gamma / von Mises distributions |
| `used_available` | Kernel densities of each covariate at used (solid) vs available (dashed) step ends, original units |
| `uhc` | Used-habitat calibration on held-out steps: model-predicted densities (grey), observed used, available |
| `home_ranges` | One panel per animal: AKDE home-range contour (solid) with its 95% CI contours (dotted) and the GPS fixes, km, shared scale |
| `hr_area` | Home-range area per animal: AKDE with 95% CI, KDE (href) and MCP |

Colour roles are fixed across figures: used = vermillion, available = blue,
population mean = black, individuals = grey; AKDE/KDE/MCP = Okabe-Ito 1-3
(each also has its own linetype or shape, so greyscale works). Animals are
facets, not colours, because collar studies often exceed 8 animals.
Relabel anything through `figures: labels:` (covariate keys, `rss`,
`sel_axis`, `step_length`, `turn_angle`, `density`, `used`, `available`,
`population_mean`, `individual`, `kernel_tentative`, `kernel_updated`,
`akde`, `kde`, `mcp`, `easting`, `northing`, `hr_area`, `animal`, and animal
ids). AKDE intervals and CI contours are ctmm's 95% regardless of `ci`.

## Simulation and power

`simulate()` builds a 16 x 16 km landscape (50-m cells) with smooth forest
cover (%) and elevation (m) fields, places animals' range centres in the
middle 6 km, and moves each animal step by step: 200 candidate steps from
Gamma(shape, scale) x von Mises(0, kappa), one chosen with probability
proportional to exp(beta_i' z(end) - attraction * d^2), where d is distance
(km) from the animal's true range centre. Fixes are then lost at random
(`fix_rate`) and timestamps jittered.

* params: `sel_mean`, `sel_sd` (named lists per covariate, per SD), `step_shape`,
  `step_scale` (m), `turn_kappa`, `attraction` (per km^2; keeps animals
  range-resident).
* design: `n_animals`, `n_steps` (per animal), `interval_hours`, `fix_rate`,
  `covariates` (`forest`, `elevation`, or other names = generic smooth fields),
  `field_range_m` (spatial scale of habitat patches), `n_candidates`,
  landscape size/cell, CRS.
* Sensible knobs for power: `n_animals` (precision of the population mean
  and SD is driven by the number of animals far more than by steps per
  animal once each animal has a few hundred steps) and `n_steps`.

## Pitfalls

* **Number of animals.** The population-level interval reflects variation
  among animals: with 5 animals it is wide, and the random-slope SD is biased
  low. Report individual estimates and treat the SD as rough.
* **Home-range attraction term (on by default).** Range-resident animals are
  pulled toward their range. Leaving this out when it is real is a
  mis-specification: in validation (simulated attraction 0.25 per km^2,
  fit without the term) the elevation coefficient's interval covered the
  truth in only 83% of replicates and the turning-concentration interval in 80%.
  With the term, the centre is estimated per animal through animal-specific
  linear terms in the coordinates (`-a|x - c|^2` expands to `-a(hx^2 + hy^2)
  + 2a c_x hx + 2a c_y hy + const`). The attraction coefficient itself is a
  nuisance parameter: with the centre estimated from the same track it is
  overestimated (about +24% with 400 steps per animal, +6% with 1500 steps,
  unbiased when the true centre is supplied; the familiar small-sample bias
  of a mean-reversion rate when the mean is estimated). Do not interpret its
  size. Set `home_attraction: false` for migrating or dispersing animals.
* **Too few available steps** gives noisier coefficients; use >= 10, ideally
  20-100 (Fieberg et al. 2021: "the more the merrier").
* **RSF on GPS data**: coefficients mean something different (home-range scale
  selection, conflated with movement constraints) and naive intervals are too
  narrow because fixes are autocorrelated. Do not use the RSF intervals for
  inference.
* **KDE/MCP on autocorrelated data** underestimate home-range area; use AKDE and
  check its effective sample size (DOF area); < ~10 means few range crossings.
* **Rasters must cover the available steps.** Available steps falling outside
  the rasters are dropped and reported; a large share biases availability.
* **Scale of covariates.** beta is per 1 SD of the raster; RSS per meaningful
  unit is printed too. Do not compare per-SD coefficients across studies with
  different landscapes.

## Validation summary (details in VALIDATION.md)

* Default design (12 animals x 400 two-hour steps), 100 replicates: all selection means, SDs,
  RSS and movement-kernel parameters pass the gate (|rel. bias| <= 3.1%, coverage 0.91-0.98).
  False-positive rate for a zero mean effect 0.06.
* 5 animals x 150 steps: forest SD biased -23% (fails); means pass.
* Negative controls fail as required: an RSF on the same tracks (elevation bias +40%, coverage
  0.52) and omitting range residency (coverage 0.80-0.85 for elevation, turning, step mean).
* Not recovered: the attraction coefficient (nuisance, +25%); AKDE interval coverage 0.85 at
  ~13 effective range crossings.
* Cross-checks: glmmTMB vs survival::clogit vs hand-coded likelihood agree to 5.5e-4; movement
  updates identical to amt::update_gamma / update_vonmises; random-slope mean/SD agree with a
  two-step meta-analysis.

## References

* Avgar T, Potts JR, Lewis MA, Boyce MS (2016) Integrated step selection analysis: bridging the gap between resource selection and animal movement. Methods in Ecology and Evolution 7:619-630.
* Avgar T, Lele SR, Keim JL, Boyce MS (2017) Relative selection strength: quantifying effect size in habitat- and step-selection inference. Ecology and Evolution 7:5322-5330.
* Fieberg J, Signer J, Smith B, Avgar T (2021) A 'How to' guide for interpreting parameters in habitat-selection analyses. Journal of Animal Ecology 90:1027-1043.
* Muff S, Signer J, Fieberg J (2020) Accounting for individual-specific variation in habitat-selection studies: efficient estimation of mixed-effects models using Bayesian or frequentist computation. Journal of Animal Ecology 89:80-92.
* Signer J, Fieberg J, Avgar T (2019) Animal movement tools (amt): R package for managing tracking data and conducting habitat selection analyses. Ecology and Evolution 9:880-890.
* Signer J, Fieberg J, Avgar T (2017) Estimating utilization distributions from fitted step-selection functions. Ecosphere 8:e01771.
* Fieberg JR, Forester JD, Street GM, Johnson DH, ArchMiller AA, Matthiopoulos J (2018) Used-habitat calibration plots. Ecography 41:737-752.
* Fleming CH, Fagan WF, Mueller T, Olson KA, Leimgruber P, Calabrese JM (2015) Rigorous home range estimation with movement data: a new autocorrelated kernel density estimator. Ecology 96:1182-1188.
* Calabrese JM, Fleming CH, Gurarie E (2016) ctmm: an R package for analyzing animal relocation data as a continuous-time stochastic process. Methods in Ecology and Evolution 7:1124-1132.
* Brooks ME et al. (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9:378-400.
* Fortin D et al. (2005) Wolves influence elk movements: behavior shapes a trophic cascade in Yellowstone National Park. Ecology 86:1320-1330.
