---
id: distance
title: Distance sampling (conventional, multiple-covariate and hierarchical)
when_to_use: Observers walked line transects or stood at points and recorded the distance (perpendicular for lines, radial for points) to every animal or group they detected, with detection falling off with distance; the question is how many animals there are per km2 and in the study area, and (optionally) whether density is associated with habitat measured per transect. Exact or binned distances, single or repeated visits, groups (clusters), strata.
not_for: Counts without distances (use nmixture), detection/non-detection (occu-single), individually identified animals (scr or capture-recapture), camera-trap distance sampling with snapshot moments (not implemented), surveys where animals on the line are missed (g(0) < 1, e.g. diving whales or cryptic birds: needs double-observer mrds, not implemented), spatial density surface models over a full study area (DSM; only covariate-based prediction to a grid is offered).
data_needed: observations CSV (transect id, distance, optional group size, optional visit, optional observer/weather), transects CSV (id, line length or number of visits, stratum, habitat covariates; every surveyed transect, including those with no detections), study-area size (km2, per stratum) for abundance, optional prediction grid
engine: Distance::ds (CDS/MCDS, design-based density); unmarked::distsamp / gdistsamp (hierarchical, covariates on density, availability); hand-written likelihoods (cross-check)
---

# Distance sampling

## Question it answers

"How many deer (birds, primates, tortoises...) are there per square kilometre and in my study
area, and is density higher in some habitats?" Observers miss more animals the farther away
they are. Distance sampling models how detection falls with distance (the **detection
function**), estimates the share of animals within the truncation distance that were seen,
and corrects the counts for the ones that were missed.

| estimand | meaning |
|---|---|
| `density_per_km2` | individuals per km2 in the study area (design-based, Distance::ds; for `engine: hds/gds/rlik` from that engine) |
| `abundance` | individuals in the study area (needs `study_area_km2`) |
| `density_per_km2_<stratum>`, `abundance_<stratum>` | the same per stratum |
| `group_density_per_km2`, `mean_group_size` | groups per km2 and mean size of detected groups (when sizes are recorded) |
| `detection_prob` | probability that an animal within the truncation distance is detected (P_a) |
| `esw_m` / `edr_m` | effective strip half-width (lines) / effective detection radius (points) |
| `sigma_m` | half-normal scale (m) at the reference level of detection covariates (reported when the half-normal key without adjustments is selected); `hr_scale_m`, `hr_shape` for the hazard-rate |
| `detection_beta_<cov>` | MCDS: log-scale effect of a covariate on the detection scale (e.g. observer) |
| `encounter_rate_per_km` / `_per_point` | detections per km of line / per point visit |
| `density_at_average_site` | HDS: expected density (individuals/km2) at a transect with average covariates |
| `density_beta_<cov>`, `rate_ratio_<cov>` | HDS: log-scale effect of 1 SD of a covariate on density, and the multiplicative change exp(beta) |
| `density_per_km2_hds`, `abundance_hds` | HDS density/abundance over the study area (prediction grid if given, else mean of transects per stratum) |
| `availability_prob` | gdistsamp: probability an animal of the superpopulation is present/available on a visit |

## Model

**Detection functions** (distance x in m, truncation w):
half-normal g(x) = exp(-x^2 / 2 sigma^2); hazard-rate g(x) = 1 - exp(-(x/sigma)^-b) (Hayes &
Buckland 1983); optional cosine or Hermite adjustment terms added by forward AIC selection
(0, 1, 2 terms; Buckland et al. 2001). Line transects: perpendicular distances have pdf
f(x) = g(x)/int_0^w g. Point transects: f(r) = r g(r)/int_0^w r g. P_a = int_0^w g / w (lines),
2 int_0^w r g / w^2 (points). Multiple-covariate distance sampling (Marques et al. 2007):
log sigma_i = alpha_0 + alpha' z_i (observer, weather, cluster size; no adjustments).

**Design-based density** (Distance::dht2): D = n * E[s] / (a * P_a), a = covered area
(2wL for lines, k pi w^2 for points). Var(D) combines the encounter-rate variance (R2 estimator
for lines, P2 for points; Fewster et al. 2009), the detection-function variance (delta method)
and the group-size variance; 95% log-normal intervals with Satterthwaite df (Buckland et al. 2001).
Abundance = D x study area; strata are estimated separately and combined area-weighted.
With `keys: [hn, hr]` the best model of each key family (adjustments chosen by AIC within the family)
is fitted and density, abundance, P_a, ESW/EDR and encounter rate are **model-averaged with Akaike
weights**, with the unconditional SE of Buckland, Burnham & Augustin (1997); key-specific parameters
(sigma, shape, covariate effects) come from the AIC-best model. `variance: bootstrap` instead
resamples transects within strata `nboot` times, repeating the whole detection-function selection in
every resample (Buckland et al. 2001 sec. 3.6.4; the approach of `Distance::bootdht`) and reports the
bootstrap SD and 2.5/97.5% percentile interval. It is slow (each resample refits all candidates:
roughly 10-15 s per resample with adjustments, i.e. 30-50 min for 200 resamples), so it is optional.

**Hierarchical distance sampling** (Royle et al. 2004; unmarked::distsamp): counts y_ij of
detections on transect i in distance band j (exact distances are binned into `hds_bins`
equal bands) are Poisson(D_i a_i pi_ij), log D_i = beta_0 + beta' x_i (numeric covariates
standardised), pi_ij the band probability under g. With group sizes, density of individuals =
group density x mean observed group size. **Repeated visits** (Chandler et al. 2011;
gdistsamp): M_i ~ Poisson(D_i a_i), available on visit t N_it ~ Binomial(M_i, phi), detections
multinomial. Derived quantities and their intervals come from 2000 multivariate-normal draws
of the MLEs.
The HDS key (half-normal or hazard-rate) is chosen by AIC on the binned counts after a stability
screen: a key is dropped (when another key fits stably) if it has a link-scale SE > 3, a singular
Hessian, or a CV of the average detection probability > 50%, the same rule applied to the Distance::ds
candidates. Without the CV rule a hazard-rate fitted as a near-spike to six point-transect bands
occasionally won by AIC and gave densities 100-1000x too high (VALIDATION.md).

`engine: auto` uses ds when there are no density covariates, ds + hds when there are (ds for
density/abundance/detection, hds for the covariate effects), and gds for repeated visits.
`rlik` refits both likelihoods with hand-written code (no package) as an independent check.

## Assumptions (and what breaks)

1. **g(0) = 1**: everything on the line/point is detected. If not, density is biased low by the
   factor g(0); needs double observers (not implemented).
2. **Detection at the initial location.** Animals that move away before detection (evasive
   movement) produce a dip near zero and **density biased low**; attraction biases it high.
   Simulation (VALIDATION.md): 75 m evasive movement -> density -27%, interval coverage 0.45.
3. **Accurate distances.** Heaping at 0 or at round numbers distorts the shape. Rounding small
   distances to 0 (common with estimated, not measured, distances) inflates the spike at zero
   and **biases density high** (VALIDATION.md). Check the Q-Q plot; bin the data if heaping.
4. **Transects placed independently of animals** (random/systematic design), so the encounter
   rate represents the study area. Roads/trails as transects violate this.
5. **Independent detections** (groups counted as one detection with a size).
6. HDS: Poisson abundance per transect given covariates; density is constant within a transect.
7. gds: availability is constant and animals move around their home range between visits, so
   distances are independent across visits (Chandler et al. 2011). Availability is identified only
   from the visit-to-visit variation in counts and is weakly identified (see VALIDATION.md, Known
   limits). If animals are sedentary (same distance on every visit) availability is over-estimated
   and density under-estimated: pool the visits as effort with the single-visit engines instead.

## Data format

`observations.csv` (one row per detection):

```
transect_id,date,observer,perp_distance_m,group_size
MD01,2025-01-19,JK,61,4
MD01,2025-01-19,JK,83,2
```

`transects.csv` (one row per transect; transects with zero detections MUST be listed):

```
transect_id,stratum,length_km,shrub_cover_pct,road_dist_km,x_km,y_km,observer
MD01,south,4.06,24.9,5.32,14.5,4.5,JK
```

Points: no length column; one row per point, `visits_col` gives the number of visits.
Binned data: `distance_type: binned` with `cutpoints_m` and either `distbegin_col`/`distend_col`
or a distance column whose values are assigned to bins. Optional `prediction_grid` CSV:
`x, y, area_km2` and the density covariates (for abundance by prediction and the density map).

## Config

```yaml
data:
  survey: line                 # line | point
  species: mule deer
  group_word: group            # word used in the report ("animal", "group", "covey", "pod")
  observations: deer_observations.csv
  transects: deer_transects.csv
  transect_col: transect_id
  distance_col: perp_distance_m
  distance_units: m            # m | km | ft | yd
  distance_type: exact         # exact | binned
  cutpoints_m: [0, 50, 100, 200, 400]   # binned only
  distbegin_col: null          # binned only, optional
  distend_col: null
  size_col: group_size         # optional: cluster (group) sizes
  visit_col: null              # optional: visit number (repeated visits -> gdistsamp)
  visits_col: null             # optional: number of visits per transect/point (effort)
  length_col: length_km        # lines
  length_units: km             # km | m
  stratum_col: stratum         # optional
  study_area_km2: {north: 180, south: 180}   # or one number; omit for density only
  truncation_m: 400            # or "10%" (drop farthest 10%); decide BEFORE looking at estimates
  prediction_grid: study_area_grid.csv       # optional
formulas:
  density: "~ shrub_cover_pct + road_dist_km"   # transect covariates on density (HDS)
  detection: "~ 1"                              # e.g. "~ observer" (MCDS)
settings:
  engine: auto                 # auto | ds | hds | ds+hds | gds | rlik
  keys: [hn, hr]               # candidate key functions (AIC)
  adjustment: cos              # cos | herm | none
  max_adjustments: 2
  hds_bins: 10                 # bands used by HDS when distances are exact
  n_draws: 2000                # MVN draws for HDS derived quantities
  gof_nboot: 200               # parametric bootstrap for the HDS GOF (0 = skip)
  K: null                      # gdistsamp upper bound (auto)
  rlik_key: hn
  model_average: true          # AIC-weight average over key families (false = AIC-best model only)
  variance: analytic           # analytic (delta method + R2/P2) | bootstrap (transects resampled)
  nboot: 200                   # bootstrap resamples when variance: bootstrap
figures: {labels: {shrub_cover_pct: "Shrub cover (%)"}}
```

`example/config.yaml`: mule deer winter-range survey (36 lines, 145 detected groups, two
strata, prediction grid). `example/config_mcds.yaml`: the same with `detection: "~ observer"`.

## Outputs

Estimates table (above), an AIC table of candidate detection functions, interpretation
paragraphs, and figures: `detection_function` (histogram scaled to detection probability with
the fitted curve and 95%/50% bands), `qq_gof` (exact distances; Cramer-von Mises) or `gof_bins`
+ `gof_ppc` (binned/HDS; Freeman-Tukey parametric bootstrap), `density_<cov>` (conditional
effects on density, response scale, original units, per-transect estimates as points),
`effects_panel`, `coefficients` (rate ratios), `density_by_stratum`, `density_map` (prediction
grid: mean and SD).

## Simulation and power

`simulate()` places groups by an inhomogeneous Poisson process (log density linear in a
standardised habitat covariate that varies between survey cells), each cell with one line or
point; group sizes 1 + Poisson(mean - 1); detection by the chosen key, optional observer effect
on sigma; distances simulated to 1.5 w so truncation is realistic. Violations for negative
controls: `responsive_move_m` (evasive movement), `heap_zero_m` (distances below it recorded as
0), `round_m`; `fixed_positions` (repeated visits: animals keep their distance on every visit).
Truth: realised density/abundance in the study area, parameter values otherwise.

params: `group_density_at_average_site` (4 /km2), `density_beta_habitat` (0.4), `sigma_m` (150),
`hr_shape`, `detection_key` (hn/hr), `mean_group_size` (1 = no groups), `availability_prob`,
`detection_beta_observer`. design: `survey`, `n_transects` (30), `transect_length_km` (3),
`cell_width_km`, `point_cell_km2`, `truncation_m` (400), `distance_type`, `n_bins`, `n_visits`,
`n_strata`, `observers`, `responsive_move_m`, `heap_zero_m`, `round_m`.
Useful power knobs: `n_transects`, `transect_length_km`, `density_beta_habitat`.

## Validation summary (VALIDATION.md)

Default design (30 lines, ~145 detections), nsim = 200: all 8 estimands pass (density -1.6%, coverage 0.985;
habitat slope -1.2%, 0.92). MCDS + groups + strata, binned half-normal points, hazard-rate points (200 points)
and gdistsamp at 120 transects x 6 visits pass. Fails (documented limits): 8 lines / ~38 detections (HDS
intercept coverage 0.89), hazard-rate binned points with ~38 detections (density +14%), gdistsamp at 60 x 4
(density +23%). False-positive rate for a zero habitat effect 0.055. Cross-check: Distance::ds = hand-coded
likelihood to < 1e-6; unmarked::distsamp = hand-coded HDS to < 1e-3 on the slope.

## Pitfalls

- Transects with no detections must be in the transects file, or density is biased high.
- Choose the truncation distance before looking at density (5-10% of the farthest distances is usual).
- Fewer than ~60 detections: detection function shape and AIC choice are unstable.
- Systematic designs with strong trends: R2/P2 overstate variance (conservative; Fewster et al. 2009).
- HDS covariate effects are associations; a habitat model fitted to non-randomly placed transects
  should not be extrapolated outside the surveyed covariate range.
- Distance's reserved column names (e.g. `observer`) are handled internally; you can use them.

## References

Buckland et al. (2001, 2015); Miller et al. (2019) J. Stat. Softw. 89(1); Marques et al. (2007)
Auk 124:1229; Fewster et al. (2009) Biometrics 65:225; Innes et al. (2002) NAMMCO Sci. Publ. 4:169;
Royle, Dawson & Bates (2004) Ecology 85:1591; Chandler, Royle & King (2011) Ecology 92:1429;
Fiske & Chandler (2011) J. Stat. Softw. 43(10); Buckland, Burnham & Augustin (1997) Biometrics 53:603; Hayes & Buckland (1983) Biometrics 39:29;
Kery & Royle (2016) AHM vol. 1 ch. 8. Full citations in `module.R` (`references`).
