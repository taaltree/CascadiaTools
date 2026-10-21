---
id: scr-search
title: Search-encounter SCR (scat-detection dogs, snow-tracking, transect searches for genetic samples)
when_to_use: Individuals identified from noninvasive genetic samples (scats, hair, feathers, urine in snow) found along recorded search routes (GPS tracks of detection-dog teams, snow-tracking, walked transects) rather than at fixed traps; question is density / abundance, how density relates to habitat, and how search conditions affect detection, with uneven search effort.
not_for: Fixed traps, cameras or hair snares (use scr); unidentified detections (use occupancy or distance sampling); open populations across years (use scr-open); searches without GPS/effort records (effort cannot be reconstructed, results will be biased where searching was concentrated).
data_needed: GPS track log of every search (track_id, x, y; projected metres) with an occasion per track, genotyped sample table (individual ID + coordinates, ideally track ID), habitat grid with density covariates (optional), track-level detection covariates (team, snow) (optional)
engine: nimble (Bayesian marginal likelihood with local evaluation; default); secr (ML: tracks as transect detectors with exact locations, and the same cell model) as independent checks
---

# Search-encounter spatial capture-recapture

## Question it answers

How many animals live in the study area and at what density, when individuals
are identified from genetic samples found **along search routes**? Is density
higher in some habitats (e.g. forest cover)? Do some search conditions (an
experienced dog team, fresh snow) find more samples per km searched? The model
also gives the spatial scale of space use (sigma, and a 95% home-range area)
and the sample-encounter rate per km of search.

Typical studies: scat-detection-dog surveys of fishers, martens, bears, wolves,
kit foxes (Thompson et al. 2012); hair/tissue collected while snow-tracking
mountain lions or wolverines (Russell et al. 2012); large-scale noninvasive
genetic sampling of wolverines and bears by searchers whose GPS tracks are
recorded (Milleret et al. 2019; Moqanaki et al. 2021); walked area searches of
plots (Efford 2011; Royle et al. 2011).

## Model

Each animal i has a fixed activity centre s_i in the state space (habitat mask
= the searched area plus a buffer). Activity centres follow a Poisson point
process with intensity

    D(x) = exp(beta0 + beta' z(x))        z = standardised density covariates

The recorded search tracks are cut into square **detector cells** (side
`cell_size`). A *search unit* u is one track inside one cell; its effort E_u is
the length of track inside the cell (km). The expected number of samples of
animal i found in unit u is

    lambda_iu = lambda0 * E_u * exp(alpha' w_u) * exp(-d(s_i, c_u)^2 / (2 sigma^2))

with c_u the detector location of the cell: by default (`detector_location: centroid`) the
length-weighted centroid of the track searched inside the cell, which approximates the hazard
integrated along the track more closely than the geometric cell centre (`detector_location:
centre`; with 1 km cells = 0.67 sigma the centre inflated sigma by 3-5% and the 95% home-range
area by up to 10%, VALIDATION.md "Independent verification"), lambda0 the number of samples found per km searched
right at an activity centre, and w_u search-condition covariates (track-level:
dog team, snow; or cell-level from the habitat grid). Observation model:

- `observation: poisson` (default): number of samples y_iu ~ Poisson(lambda_iu)
- `observation: bernoulli`: sampled at least once in u ~ Bernoulli(1 - exp(-lambda_iu))

This is the search-encounter model of Royle et al. (2011) and Thompson et al.
(2012) with search paths discretised into cells as in Thompson et al. (2012),
Milleret et al. (2019) and Moqanaki et al. (2021); it is the discrete form of
Efford's (2011) transect/area-search likelihood, where the expected count is
the hazard integrated along the searched line.

**Engines.**

| engine | method | role |
|---|---|---|
| `nimble` (default) | Bayesian; data augmentation (M, psi ~ U(0,1)) with z and s **integrated out over the mask** (marginal likelihood); custom `dSearchMarg` with **local evaluation** (pixel-cell pairs beyond radius R = max(6 x sigma guess, ...) dropped; each identified animal's centre only evaluated within R of all its sample cells: Milleret et al. 2019); joint RW_block sampler on the 3-5 top-level parameters | main fit |
| `secr_transect` | maximum likelihood in secr; every GPS track is a `transect` detector on its own occasion, samples enter with their **exact** positions (projected onto the track), Poisson counts (`binomN = 0`), hazard half-normal (HHN) | independent route: no cells, no local evaluation, different software. Slow (~25 min at the default design) |
| `secr_cells` | maximum likelihood in secr of the same cell model (one `count` detector per cell, each track an occasion, usage = km searched in the cell; track covariates as `timecov`) | likelihood check of the Nimble engine (~30 s) |

Priors (nimble): psi ~ U(0,1); lambda0 ~ U(0, 50 per km); sigma ~ U(0, half
the state-space width); density and detection coefficients ~ Normal(0, SD 2)
on the standardised scale. Abundance per draw: n + Binomial(M - n, P(z = 1 |
never sampled)); density = psi M / area of the state space. If N presses
against M (99.5% quantile > 0.9 M) the fit is repeated with 2M; if sigma is
large relative to R the fit is repeated with a larger radius (both noted in
the report).

Effort models (`settings: effort_model`):

- `track_length` (default, validated): E_u = km of that track inside the cell.
- `visits`: each track passing through a cell counts as one standard visit
  (length inside the cell ignored). Sensitivity analysis only.
- `uniform`: search effort ignored, every cell of the searched rectangle is
  treated as equally searched. Shown in VALIDATION.md to bias the density
  covariate effect when search effort is correlated with habitat. Never use
  for inference.

## Assumptions (and what breaks if they fail)

1. **Closure** during the survey (no births, deaths, permanent emigration).
   Violations inflate N. Keep the survey within one season.
2. **Fixed activity centres** and circular half-normal space use, the same
   sigma for all animals. Sex differences in sigma (males range more) cause
   heterogeneity; with strong differences fit sexes separately.
3. **All samples come from recorded searches.** An opportunistic scat found
   off-track has no effort; it inflates density where it is. Drop such samples
   or record the route that found them.
4. **The GPS track is the searched line.** A dog team searches a strip either
   side of the walked line; that strip is absorbed into lambda0 (samples per
   km), so it must be similar across tracks or modelled with a detection
   covariate (team, terrain, snow).
5. **Correct genotypes.** Allelic dropout / false alleles create "ghost"
   individuals with one sample each, inflating N. Use strict QC (re-amplify
   singletons).
6. **Independence** of samples given the activity centre. Latrines (many scats
   at one spot) produce over-dispersion: the posterior predictive check
   `disp_samples` flags it. Use `observation: bernoulli` (one detection per
   search unit) as a robust alternative.
7. **Cell size <= sigma.** Cells coarser than sigma blur where samples were
   found and bias sigma upward; the report warns when cell_size > sigma.
8. **Buffer >= 3-4 sigma** around the searched area; otherwise animals living
   outside the state space are sampled and density is inflated.

## Data format

**Search tracks** (CSV, one row per GPS fix, in walking order or with a time
stamp; or GPX files with `data: crs`):

| track_id | gps_time | utm_e | utm_n |
|---|---|---|---|
| T001 | 2025-06-12T07:34:35Z | 568304.8 | 4883730.0 |
| T001 | 2025-06-12T07:35:35Z | 568266.6 | 4883762.3 |

**Track field sheet** (optional `track_info` CSV; one row per track; any
column can be a detection covariate):

| track_id | date | session | dog_team | experienced_team | km_walked |
|---|---|---|---|---|---|
| T001 | 2025-06-12 | 1 | Pepin & Sara | 1 | 9.2 |

**Samples** (CSV, one row per collected sample; a blank / `fail` individual
ID = failed genotype, dropped):

| scat_id | track_id | date | utm_e | utm_n | fisher_id | genotype_qc |
|---|---|---|---|---|---|---|
| PEPE-25-0001 | T001 | 2025-06-12 | 564731.7 | 4885101.8 | F013 | pass |

If the sample table has no track ID, each sample is attributed to the nearest
track on its occasion (or the nearest track overall).

**Habitat grid** (optional CSV of regular pixel centres covering the searched
area + buffer, with density covariates; or `data: raster` (GeoTIFF) whose
layers are extracted at an automatically built mask). Without a grid the mask
is built from the tracks with a buffer of 4 x a rough sigma guess.

## Config

```yaml
data:
  tracks: search_tracks_gps.csv     # CSV of GPS fixes, or a list of .gpx files (then give crs:)
  crs: "EPSG:32610"                 # only for GPX (projects lon/lat to metres)
  track_id_col: track_id            # default track_id
  x_col: utm_e                      # default x (also used for samples unless sample_x_col)
  y_col: utm_n                      # default y
  order_col: gps_time               # optional: sort fixes within each track
  track_info: track_field_sheet.csv # optional per-track table joined by track_id
  occasion_col: session             # occasion (survey session) column in tracks or track_info (default occasion)
  n_occasions: 3                    # optional (default = max occasion)
  samples: fisher_scats.csv
  sample_id_col: scat_id            # optional
  individual_col: fisher_id         # genotype ID; blank/NA/"fail" = failed sample
  sample_x_col: utm_e               # optional (defaults to x_col)
  sample_y_col: utm_n
  sample_track_col: track_id        # optional: route on which each sample was found
  sample_occasion_col: session      # optional
  mask: habitat_grid.csv            # optional habitat grid CSV
  mask_x_col: utm_e                 # default x
  mask_y_col: utm_n                 # default y
  habitat_col: habitat              # optional 0/1 habitat flag (0 = not habitat, removed)
  raster: covariates.tif            # optional alternative to mask (layers named as covariates)
  mask_spacing: 700                 # optional (inferred from the grid)
  buffer: 6000                      # optional (used when the mask is built automatically)
  units: m                          # m or km (or omit if unknown: estimates then per map unit)
  cell_size: 1000                   # REQUIRED in practice: detector cell side (about sigma/2 to sigma)
  density_covariates: [forest_pct]  # columns of the habitat grid / raster
  detection_covariates: [experienced_team]   # track columns (team, snow) or grid columns
  labels: {forest_pct: "Forest cover (%)"}
settings:
  engine: nimble                    # nimble | secr_transect | secr_cells
  observation: poisson              # poisson (number of samples) | bernoulli (sampled or not per unit)
  effort_model: track_length        # track_length | visits | uniform (sensitivity only)
  detector_location: centroid       # centroid (track-length-weighted, default) | centre (geometric cell centre)
  niter: 40000                      # MCMC (per chain), nburnin: 10000, nchains: 3, thin: 1
  M: null                           # data augmentation size (default max(4n, n + 150); doubled automatically if hit)
  local_radius: null                # local evaluation radius in map units (default 6 x sigma guess); no_local: true to disable
  sigma_max: null                   # upper bound of the sigma prior (default half the state-space width)
  compare_engines: [secr_cells]     # also fit with these engines and draw engine_comparison
  surface_draws: 150                # posterior draws used for density / detection maps
  ppc_draws: 200                    # posterior predictive check draws
```

Categorical detection covariates with two levels are coded 0/1 (alphabetically
second level = 1, noted in the report); use 0/1 columns for clarity.
Continuous covariates are standardised (mean 0, SD 1); figures show original
units.

## Outputs

Estimates (`fit()$estimates$param`):

| param | meaning |
|---|---|
| `density_per_km2` | animals per km² in the state space (mean density) |
| `abundance` | realised number of activity centres in the state space (searched area + buffer) |
| `lambda0_per_km` | samples of an animal found per km searched at its activity centre (reference level of detection covariates) |
| `sigma` | spatial scale of space use (map units, m) |
| `home_range_area_95` | pi (2.448 sigma)² in km²: 95% area of a circular bivariate-normal home range |
| `samples_per_track_through_centre` | lambda0 sqrt(2 pi) sigma: expected samples from a long straight search route through the activity centre (a design-friendly summary of detectability) |
| `density_beta_<cov>` | log density ratio per 1 SD of the covariate (exp = rate ratio) |
| `detection_beta_<cov>` | log ratio of samples per km for level 1 vs 0 (or per 1 SD) |

Figures (`figures()`): `effort_map` (km searched per cell, tracks, samples),
`density_surface` (realised density + estimated activity centres),
`density_expected_surface`, `detection_map` (P(sampled at least once)),
`density_<cov>` (conditional effect, animals per km²), `detection_function`
(samples per km vs distance, by team), `detection_<cov>`,
`abundance_posterior`, `captures_per_individual` (observed vs posterior
predictive), `gof` (posterior predictive checks), `coefficients`,
`engine_comparison` (when `compare_engines`), `trace`, `posterior`.

## Simulation and power

`design` knobs: `search_side` (side of the square searched region, m; 16000),
`buffer` (6000), `mask_spacing` (700), `cell_size` (1000), `n_occasions` (3),
`tracks_per_occasion` (12), `track_km` (c(4, 10): length range of each search
bout), `n_roads` (3), `road_bias` (share of tracks that follow a road; 0.3),
`cov_road_corr` (correlation between closeness to roads and the density
covariate; 0), `covariate_range` (spatial range of the habitat field; 3000),
`density_covariate`, `detection_covariate` (TRUE/FALSE).
`params`: `density` (0.10 per km²), `lambda0_per_km` (0.20),
`sigma` (1500 m), `density_beta_forest` (0.5), `detection_beta_experienced_team` (0.6).

The simulator deposits samples **continuously** along the walked GPS paths
(50 m pieces) at a rate lambda0 x length x half-normal(distance), with 10 m
location jitter, so the cell discretisation used by the fitted model is
tested, not assumed. Searches follow winding roads (share `road_bias`) or are
correlated random walks off-road.

Sensible ranges: aim for >= 20 identified individuals and >= 10 spatial
recaptures; total effort such that `samples_per_track_through_centre` x
(number of tracks passing near a typical home range) is >= 2. Power for the
dog-team effect and the forest effect: VALIDATION.md (tracks per occasion).

## Pitfalls

- **Not recording effort.** Without GPS tracks of every search, effort must be
  guessed; uneven searching (along roads, in good habitat) then looks like
  higher density there. VALIDATION.md shows ignoring effort biases the habitat
  effect.
- **Unprojected coordinates.** Use projected metres (UTM). GPX files need `crs`.
- **Merged tracks.** Two search days under one track ID create a straight
  "jump"; `check_data` flags jumps > 5 cells.
- **Too-large cells** (> sigma) or too small a buffer (< 3 sigma): flagged.
- **Abundance refers to the whole state space** (searched area + buffer);
  quote density, or N for a defined region, not N alone.
- **Few spatial recaptures** (< 5): sigma, and therefore density, is weakly
  identified; the posterior of N can be very right-skewed.
- **Habitat effect far from the searches.** Density in unsearched parts of
  the state space is extrapolated from the habitat model (see `detection_map`).
  When searching is concentrated in one habitat, the posterior of N for the
  whole state space can have a long right tail (VALIDATION.md, Known limits 1-2).
  Report the posterior median alongside the mean and cross-check with
  `compare_engines: [secr_cells]`.
- **Posterior mean of density / N at low precision.** In independent simulations
  with 20-35 identified animals (posterior CV of density 0.23-0.28) the posterior
  MEAN of density was 9-13% above the truth (ML within 4%; interval coverage stayed
  0.90-0.94). The fit adds a note with the posterior medians when CV > 0.2.

## References

- Royle JA, Kéry M, Guélat J (2011) Spatial capture-recapture models for search-encounter data. Methods in Ecology and Evolution 2:602-611. doi:10.1111/j.2041-210X.2011.00116.x
- Thompson CM, Royle JA, Garner JD (2012) A framework for inference about carnivore density from unstructured spatial sampling of scat using detector dogs. Journal of Wildlife Management 76:863-871. doi:10.1002/jwmg.317
- Russell RE, Royle JA, Desimone R, Schwartz MK, Edwards VL, Pilgrim KP, McKelvey KS (2012) Estimating abundance of mountain lions from unstructured spatial sampling. Journal of Wildlife Management 76:1551-1561. doi:10.1002/jwmg.412
- Efford MG (2011) Estimation of population density by spatially explicit capture-recapture analysis of data from area searches. Ecology 92:2202-2207. doi:10.1890/11-0332.1
- Efford MG. secr: Spatially explicit capture-recapture models (R package, v5.5); detector types polygon, transect, polygonX, transectX ("secr-polygondetectors" vignette); SECRbook Appendix: area and transect search.
- Milleret C, Dupont P, Bonenfant C, Brøseth H, Flagstad Ø, Sutherland C, Bischof R (2019) A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363. doi:10.1002/ece3.4751
- Moqanaki EM, Milleret C, Tourani M, Dupont P, Bischof R (2021) Consequences of ignoring variable and spatially autocorrelated detection probability in spatial capture-recapture. Landscape Ecology 36:2879-2895. doi:10.1007/s10980-021-01283-x
- Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press (Ch. 15, search-encounter models).
- Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152:S521-S537.
- Borchers DL, Efford MG (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64:377-385.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26:403-413.
