---
id: sdm-presence-only
title: Species distribution from presence-only records (point process, Maxent, integrated with surveys)
when_to_use: You have presence-only records of a species (museum/herbarium specimens, iNaturalist/eBird-style incidental sightings, agency incidental reports, trapper/hunter reports) with coordinates, plus habitat covariates covering the study region (a raster stack or a regular covariate grid), and the question is which habitat the species is associated with and where in the region it is relatively most likely to occur. Optionally a sampling-bias covariate (distance to road/town/trail, access, observer density) and optionally a set of planned survey sites with detections AND non-detections (camera traps, plots, repeat visits).
not_for: Data with real absences at every site and no presence-only records (use occu-single or glm); counts at sites (nmixture/glm); telemetry locations of marked animals (movement-ssf, rsf); estimating prevalence, occupancy or abundance from presence-only records alone (impossible without survey data, see below); records with positional error much larger than the covariate grid cell; spatially clustered records with strong residual clustering (an LGCP / spatial random effect is not implemented here).
data_needed: presence CSV (x, y in a projected metric CRS, optional species column); covariate grid CSV (x, y, covariates on a regular grid) or GeoTIFF rasters; names of habitat covariates and (strongly recommended) sampling-bias covariates; optional survey CSV (site, x, y, visit, detected 0/1) with the area each survey covers
engine: ppm = inhomogeneous Poisson point process by down-weighted Poisson regression (stats::glm.fit) on grid quadrature; integrated = joint maximum likelihood of the presence-only IPP and a cloglog occupancy model of survey data with detection probability (optim, BFGS with analytic gradient); maxnet = Maxent via glmnet with percentile-bootstrap intervals; survey_only = survey part alone (comparison)
---

# Species distribution from presence-only records

## Question it answers

"Which habitat is the species associated with, and where in my region is it
most likely to be?" when all you have are places where it **was** recorded.

The module treats the records as a spatial point pattern: an
**inhomogeneous Poisson point process (IPP)** whose intensity (expected
number of records per km²) depends on habitat, multiplied by a **thinning
(recording) rate** that depends on where people look (Warton & Shepherd
2010; Renner et al. 2015; Fithian et al. 2015). This is the same model that
Maxent fits (Renner & Warton 2013; Fithian & Hastie 2013; Phillips et al.
2017), written so that the coefficients have a clear meaning and honest
intervals.

What you get:

| Output | Presence-only (ppm / maxnet) | + survey data (integrated) |
|---|---|---|
| Habitat effects on log intensity, per SD and per natural unit, with 95% CI | yes | yes (better protected against bias) |
| Sampling-bias (observer-effort) effect | yes, if a bias covariate is given | yes |
| Map | RELATIVE intensity (1 = landscape average) | absolute density (individuals per km²) |
| Density at mean habitat, expected total abundance, mean plot occupancy, detection probability | **not identifiable** | yes |
| Conditional-effect (response) curves | relative intensity | density per km² |
| Spatial block cross-validation (AUC, continuous Boyce index), quadrature-convergence check, block-count goodness of fit | yes | yes (+ AUC on held-out survey sites) |

**What presence-only data cannot tell you.** The overall number of records
mixes how common the species is with how hard people looked. Prevalence,
occupancy probability and abundance are therefore not identifiable from
presence-only data (Ward et al. 2009; Fithian & Hastie 2013; Hastie &
Fithian 2013). The map is relative. A Maxent "probability of presence" from
the logistic/cloglog output is a rescaled relative intensity, not a
probability. Adding survey data with non-detections fixes the level
(integrated engine).

**Sampling bias confounds habitat.** If people look mostly near roads, and
roads run along valleys, then valley habitat looks preferred. This
happens even when the species is indifferent to it. Include the variable
that drives effort (distance to road/town/trail, human footprint) as a
`bias_covariate`. It enters the model for the records but is held at a
common value for prediction (Warton, Renner & Ramp 2013). This only works
if the bias covariate captures the effort pattern and is not itself a
habitat driver for the species. The VALIDATION negative control shows that
ignoring a habitat-correlated bias makes the habitat effects wrong (bias
-73%, coverage 4.5%), with confident-looking intervals.

## Model

Covariates are standardised over the study-region grid (mean 0, SD 1), so
every coefficient is "per 1 SD of the covariate over the region".

Species intensity (individuals per km²):

  log λ(s) = α + Σ_k β_k x_k(s)            (+ β_kk x_k(s)² for `quadratic` covariates)

Presence-only records: an IPP with intensity λ(s)·b(s), where the recording
(thinning) rate is

  log b(s) = γ + Σ_m δ_m z_m(s)              (z = bias covariates)

The PO log-likelihood is Σ_records log[λ(s_i) b(s_i)] − ∫_D λ(s) b(s) ds.
The integral is approximated by quadrature over the grid cells (weight = cell
area). Only α + γ (`po_intercept`) is estimable from PO data, never α alone.

- **ppm engine**: down-weighted Poisson regression (Berman & Turner 1992;
  Renner et al. 2015). Presence points get a tiny weight w = 1e-6 and
  response 1/w, and quadrature points get the cell area and response 0.
  A Poisson GLM then maximises exactly the IPP likelihood above. Wald 95%
  confidence intervals.
- **integrated engine** (Fithian et al. 2015; Koshkina et al. 2017;
  Dorazio 2014): adds survey sites of area A_i. Abundance in a plot is
  Poisson(λ(s_i)A_i), so the plot is occupied with probability
  ψ_i = 1 − exp(−λ(s_i)A_i) (complementary log-log link). With repeated
  visits, the species is detected on each visit to an occupied plot with
  probability p. The two likelihoods share α (survey part) and the habitat
  β's (both parts), so α, p and the β's become separately identifiable.
  Joint ML by BFGS with an analytic gradient, and Wald CIs on link-scale
  parameters. Derived quantities (density at mean covariates, expected
  total abundance, mean plot occupancy, detection probability) get 95% CIs
  from 2000 draws of the asymptotic normal distribution.
- **maxnet engine** (Phillips et al. 2017): Maxent fitted by penalised
  (lasso) Poisson/logistic regression via glmnet with all grid cells as
  background. The default uses linear features (`l`) and adds quadratic
  features for `quadratic` covariates. Bias covariates enter as linear
  features and are set to their mean for prediction. Intervals are 95%
  percentile bootstrap intervals over resamples of the presence records.
  Regularisation shrinks effects toward 0.
- **survey_only**: the survey likelihood alone, for comparison.

## Assumptions (and what breaks)

1. **Records are independent given the covariates** (Poisson). Clustered
   records (several records of the same animal/den, observers revisiting
   the same spot, unmodelled habitat) make intervals too narrow. The block
   goodness-of-fit test flags this. Remedies: thin duplicate records to one
   per cell per year, add the missing covariate, or treat intervals as a
   lower bound.
2. **Sampling bias is captured by the bias covariates.** Unmeasured effort
   that correlates with habitat biases the habitat effects, and nothing in
   PO data alone can detect it. Survey data (integrated) give partial
   protection.
3. **The bias covariate does not also affect the species.** If martens
   avoid roads, the "bias" coefficient absorbs that avoidance, and the
   habitat map will not show it.
4. **The study region (grid) is the region the records could have come
   from.** Records outside it are dropped. Do not include ocean, other
   countries, or areas never accessible to recorders unless a bias covariate
   describes that.
5. **Grid cells are small relative to habitat variation**, and
   record-location error is smaller than a cell. The quadrature check
   refits on every 2nd/4th/8th cell. If estimates move by more than about
   0.5 SE, the grid is too coarse.
6. **Integrated model:** survey sites and PO records respond to habitat in
   the same way. Survey sites are placed independently of the species. Each
   survey samples a known area A. Detection p is constant across sites and
   visits. Abundance-related outputs depend on the Poisson abundance
   assumption and are best read as an index.

## Data format

All coordinates must be in the **same projected, metric CRS** (UTM metres or
km). Longitude/latitude is rejected.

`presence_records.csv`: one row per record. Extra columns are ignored.

```
record_id,species,date,source,easting,northing,coord_uncertainty_m
MACA-0104,Martes caurina,1995-02-21,iNaturalist,590978,4886164,100
```

`covariate_grid.csv`: one row per cell centre of a **regular** grid
covering the study region. Missing values drop the cell.

```
easting,northing,forest_cover_pct,elevation_m,dist_road_km
560500,4850500,58.3,1654,29.983
```

Or give rasters (`data: rasters: {forest_cover_pct: forest.tif, ...}`).
They must share a CRS, extent and resolution, and are converted to the grid.

Optional `camera_survey.csv` (long format, one row per site x visit actually
done; visits not done are simply absent):

```
station,easting,northing,occasion,detected
CAM001,586500,4874500,1,0
```

A wide format is also accepted (`survey_format: wide`,
`survey_visit_cols: [v1, v2, v3]`).

## Config

```yaml
data:
  presence: presence_records.csv     # required
  presence_x: easting                # default x
  presence_y: northing               # default y
  species_col: species               # optional filter...
  species: "Martes caurina"          # ...to these species
  covariate_grid: covariate_grid.csv # OR rasters: {name: file.tif, ...}
  grid_x: easting                    # default x
  grid_y: northing                   # default y
  coord_unit: m                      # m (default) or km; sets cell area in km2
  habitat_covariates: [forest_cover_pct, elevation_m]
  bias_covariates: [dist_road_km]    # recommended; [] if none
  quadratic: [elevation_m]           # optional hump-shaped responses
  survey: camera_survey.csv          # optional
  survey_site_col: station           # default site
  survey_x: easting                  # default = presence_x
  survey_y: northing
  survey_visit_col: occasion         # default visit
  survey_detection_col: detected     # default detected
  survey_area_km2: 1                 # area each survey samples (or survey_area_col)
  labels: {forest_cover_pct: "Forest cover (%)", ...}   # axis labels with units
settings:
  engine: auto            # auto (= integrated if survey given, else ppm) | ppm | integrated | maxnet | survey_only
  bias_correction: true   # false only to demonstrate what ignoring bias does
  perfect_detection: false  # true: treat survey 0/1 as true absence/presence (single visits)
  draws: 2000             # parametric draws for derived quantities and bands
  evaluate: true          # spatial block cross-validation
  cv_folds: 5
  cv_block_size: 12000    # block side in coordinate units (default: 1/5 of the extent)
  quadrature_check: true
  gof: true
  gof_reps: 1000
  maxnet_classes: l       # maxnet features: l, lq, lqh, lqp... (h/p terms: read the curves, not the table)
  maxnet_regmult: 1
  maxnet_boot: 100        # bootstrap resamples for maxnet intervals
  compare_engines: []     # e.g. [ppm, maxnet] to fit and plot side by side
  seed: 1
```

## Outputs

Estimands (`fit()$estimates$param`):

- `intensity_beta_<cov>`: habitat slope on log intensity per 1 SD
  (`exp()` = multiplicative change in intensity/density per SD).
  `intensity_beta_<cov>_sq` for quadratic terms.
- `bias_beta_<cov>`: slope of log recording rate per 1 SD of the bias
  covariate. It describes the observers, not the species.
- Integrated only: `density_at_mean_covariates` (individuals per km²),
  `expected_total_abundance` (whole region), `mean_occupancy_prob` (mean
  probability that a survey plot of the median survey area is occupied),
  `detection_prob` (per visit).

Report: plain-language interpretation with effects per SD and per a natural
unit, plus the across-range fold change. It also shows what the habitat
effects would have been without the bias correction, the cross-validation
scores, the goodness of fit, and a methods paragraph.

Figures (`figures/`):
- `intensity_map`: relative intensity (or density) with records.
- `bias_map`: estimated recording effort.
- `effect_<cov>`: conditional effect on the response scale with 95%/50% bands
  and a rug of covariate values at records.
- `bias_effect`.
- `response_curves_panel`.
- `coefficients`: rate ratios per SD.
- `evaluation`: block-CV AUC/Boyce and Boyce P/E curves.
- `quadrature_convergence`.
- `gof`.
- `engine_comparison` if `compare_engines` is set.

## Simulation and power

`simulate()` builds an `nx` x `ny` grid of `cell_km` cells. It has two
smooth habitat fields (forest cover %, elevation m) and `n_roads` straight
roads. Elevation correlates with distance to road
(`bias_habitat_cor`; roads follow valleys). Records are an IPP thinned by
exp(δ · z(dist_road)), scaled so that E[records] = `n_presence`. Survey
sites (`survey_sites` > 0) are random cells with `survey_visits` visits
of detection probability `detection_prob`.

| knob | default | sensible range |
|---|---|---|
| `n_presence` | 200 | 25-1000 (power: 50 records detect a 0.5 SD effect ~68% of the time; see VALIDATION) |
| `survey_sites` | 0 (100 for the integrated validation) | 60+ for density; < 40 sites with p < 0.4 do not identify density |
| `survey_visits` | 3 | 2-6 |
| `bias_habitat_cor` | 0.6 | 0-0.9 (how strongly effort and habitat are confounded) |
| params `intensity_beta_*`, `bias_beta_dist_road_km`, `density_at_mean_covariates`, `detection_prob` | 0.8, -0.5, -1, 0.5, 0.5 | |

`core/run.R power sdm-presence-only --vary n_presence=25,50,100,200,400`.

## Pitfalls

- Using Maxent's logistic/cloglog output as "probability of occurrence":
  it is relative.
- Leaving out a bias covariate because "the records look spread out".
  Effort is almost never uniform.
- A study region much larger than the area searched (e.g. a whole state
  for records from two parks) makes everything look "avoided" except the
  parks.
- Duplicate records (same animal, same site, many years) give
  overconfident intervals. Thin to one per cell per season first.
- Strongly collinear habitat covariates (|r| > 0.85) are flagged by
  `check_data()` and cannot be separated; drop one. A habitat covariate
  strongly correlated with a BIAS covariate is not blocked: the fit stays
  unbiased with a wide interval (independent verification, VALIDATION.md),
  and the report notes it. Never drop the bias covariate to "fix" this.
- Survey files: `check_data()` flags detections other than 0/1 (naming site
  and visit), a site id at several locations, two site ids at identical
  coordinates (usually a typo), and duplicated site x visit rows. Presence
  rows with a missing coordinate are dropped and counted in the report notes.
- Hinge/product Maxent features (`lqh`, `lqhp`) give flexible maps, but
  the coefficient table then shows only part of each effect. Use the
  response curves.
- Integrated density with few survey sites (< ~40) and low detection is
  weakly identified. The fit is then flagged as not converged rather than
  returning nonsense with tiny intervals.

## References

- Warton DI, Shepherd LC (2010) Poisson point process models solve the "pseudo-absence problem" for presence-only data in ecology. Annals of Applied Statistics 4:1383-1402.
- Renner IW, Elith J, Baddeley A, Fithian W, Hastie T, Phillips SJ, Popovic G, Warton DI (2015) Point process models for presence-only analysis. Methods in Ecology and Evolution 6:366-379. doi:10.1111/2041-210X.12352
- Fithian W, Hastie T (2013) Finite-sample equivalence in statistical models for presence-only data. Annals of Applied Statistics 7:1917-1939.
- Renner IW, Warton DI (2013) Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics 69:274-281.
- Phillips SJ, Anderson RP, Dudík M, Schapire RE, Blair ME (2017) Opening the black box: an open-source release of Maxent. Ecography 40:887-893. doi:10.1111/ecog.03049
- Warton DI, Renner IW, Ramp D (2013) Model-based control of observer bias for the analysis of presence-only data in ecology. PLoS ONE 8:e79168.
- Fithian W, Elith J, Hastie T, Keith DA (2015) Bias correction in species distribution models: pooling survey and collection data for multiple species. Methods in Ecology and Evolution 6:424-438. doi:10.1111/2041-210X.12242
- Koshkina V, Wang Y, Gordon A, Dorazio RM, White M, Stone L (2017) Integrated species distribution models: combining presence-background data and site-occupancy data with imperfect detection. Methods in Ecology and Evolution 8:420-430. doi:10.1111/2041-210X.12738
- Dorazio RM (2014) Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. Global Ecology and Biogeography 23:1472-1484. doi:10.1111/geb.12216
- Ward G, Hastie T, Barry S, Elith J, Leathwick JR (2009) Presence-only data and the EM algorithm. Biometrics 65:554-563.
- Berman M, Turner TR (1992) Approximating point process likelihoods with GLIM. Applied Statistics 41:31-38.
- Hirzel AH, Le Lay G, Helfer V, Randin C, Guisan A (2006) Evaluating the ability of habitat suitability models to predict species presences. Ecological Modelling 199:142-152.
- Roberts DR et al. (2017) Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography 40:913-929.
- Baddeley A, Rubak E, Turner R (2015) Spatial Point Patterns: Methodology and Applications with R. CRC (spatstat, used for the cross-check).
