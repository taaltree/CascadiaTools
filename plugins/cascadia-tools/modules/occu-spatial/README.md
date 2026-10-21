---
id: occu-spatial
title: Spatial, multi-species and integrated occupancy (spOccupancy)
when_to_use: Detection/non-detection surveys repeated at sites with known coordinates, when (a) occupancy is likely spatially clustered beyond what the covariates explain (single species, spPGOcc), (b) several species were recorded at the same sites and community-level responses or species richness are wanted (msPGOcc, sfMsPGOcc), or (c) the same species was surveyed with two or more methods (point counts + ARUs, cameras + sign surveys) that should be analysed together (intPGOcc).
not_for: Single-species non-spatial analyses with likelihood output or dynamic (multi-season) occupancy (use occu-single); species interactions (occu-rota / occu-waddle); counts or abundance (N-mixture); estimating the spatial range itself precisely (weakly identified, see Known limits); presence-only data (sdm-presence-only); false-positive detections.
data_needed: long detection table (site, survey, detected 0/1 [, species] [, method/source] [, survey covariates]); site table with projected x/y coordinates (m or km) and site covariates; optional prediction grid (x, y, same covariates)
engine: spOccupancy 0.8.2 (Polya-Gamma Gibbs samplers; NNGP spatial random effects)
---

# Spatial, multi-species and integrated occupancy

## Question it answers

*Where does the species (or each species of a community) occur, which habitat
covariates are associated with occurrence, and is there residual spatial
clustering, after accounting for imperfect detection?*

Variants (choose with `settings: model:`; `auto` picks the bold one for the data type):

| data | `model` | what it adds |
|---|---|---|
| one species | **spPGOcc** | spatial random effect (Gaussian process, exponential covariance, NNGP approximation) |
| one species | PGOcc | non-spatial occupancy (same priors); used for comparison |
| several species (`species_col`) | **sfMsPGOcc** | community model + spatial factor model (residual species co-occurrence and spatial structure) |
| several species | msPGOcc | community model, non-spatial (Dorazio & Royle 2005 style) |
| several methods (`source_col`) | **intPGOcc** | one occupancy process observed by several methods, each with its own detection model |
| several methods | spIntPGOcc | integrated + spatial effect (**not validated**; experimental) |

Ignoring spatial autocorrelation does not usually bias a covariate effect much,
but it makes its interval too narrow. Our negative control shows that: see
VALIDATION.md. The spatial model is the safer default when sites are clustered
in space, but the price is longer run times.

## Model

Single species, site j (coordinates s_j, km), survey k:

    z_j ~ Bernoulli(psi_j),            logit(psi_j) = b0 + x_j' b + w(s_j)
    y_jk | z_j ~ Bernoulli(z_j p_jk),  logit(p_jk)  = a0 + v_jk' a
    w ~ MVN(0, sigma2 * exp(-phi * d_jj'))         (spPGOcc; w = 0 for PGOcc)

The effective range is 3/phi (km): beyond this distance the residual spatial
correlation is below exp(-3) = 0.05. spOccupancy replaces the full GP with a
nearest-neighbour GP (Datta et al. 2016), using `n_neighbors` neighbours
(default 15).

Multi-species (species i): each species has its own coefficients drawn from a
community distribution, `b_i ~ Normal(mu_b, tau_b^2)` and `a_i ~ Normal(mu_a, tau_a^2)`.
sfMsPGOcc adds `w_i(s) = lambda_i' f(s)`, where f holds q independent
unit-variance spatial factors (exponential GP each, own phi_r). The loadings Lambda
have their upper triangle fixed at 0 and diagonal at 1 (Doser, Finley & Banerjee 2023).
Species richness at site j = sum_i z_ij, over the modelled species only.

Integrated (method r): z_j is shared, and `y_r,jk ~ Bernoulli(z_j p_r,jk)` with
`logit(p_r,jk) = a0_r + v' a_r`. Sites may be surveyed by one method or by
several.

Priors (as in spOccupancy defaults unless stated): Normal(0, var 2.72) on all
coefficients and community means; inverse-gamma(0.1, 0.1) on among-species
variances; inverse-gamma(2, 1) on sigma2; phi ~ Uniform(3/max_dist, 3/d_05). The
lower bound for the range is the **5% quantile** of inter-site distances, not the
package default of the minimum distance (see Pitfalls). Override it with
`settings: range_prior_km: [lo, hi]`.

Estimands (names in the estimates table):

- `occupancy_prob`: occupancy probability at average covariates (and w = 0).
- `psi_beta_<term>`: occupancy slope (logit scale, per 1 SD of the term).
- `detection_prob[_<source>]`, `p_beta_<term>[_<source>]`.
- `mean_occupancy`: mean of psi_j over surveyed sites (includes w).
- `prop_sites_occupied`: finite-sample proportion of surveyed sites occupied (mean z).
- spatial: `spatial_variance` (sigma2), `spatial_decay_per_km` (phi), `effective_range_km` (3/phi);
  `effective_range_km_f<r>` per factor.
- community: `comm_occupancy_prob`, `comm_psi_beta_<term>`, `comm_sd_occupancy_logit`,
  `comm_sd_psi_beta_<term>`, and the same with `detection` / `p_beta`.
  `mean_site_richness`. Species-level rows end in `_<species>`.

Numeric model terms are standardised (mean 0, SD 1) after formula transforms, so
`log(area_ha)` works. Factors become 0/1 dummies against the first level.

## Assumptions (and what breaks if they fail)

- **Closure** within the season: no change in occupancy between surveys. If it
  fails, occupancy is read as "use" and detection is biased low.
- **No false positives.** Misidentifications inflate occupancy.
- **Independent detections given z.** Surveys on consecutive days at one camera
  can be correlated; pool them into longer occasions.
- **Spatial effect is a nuisance smooth.** If a covariate is itself spatially
  smooth at the same scale as w, the spatial effect can absorb part of its
  effect ("spatial confounding"). The covariate estimate is then less precise
  and can shift toward 0. Compare with PGOcc.
- **Multi-species:** species responses are exchangeable draws from one
  community. Very rare species are shrunk toward the community mean, and
  species never detected are not in the model.
- **Integrated:** all methods observe the same latent state at a site (same
  species, same season, same spatial unit). An ARU that hears birds 300 m
  away observes a larger area than a 50 m point count. The model does not fix
  that mismatch.

## Data format

`detections.csv`, one row per survey done (per species for multi-species):

| station | week | days_active | detected |
|---|---|---|---|
| CAM-001 | 1 | 7 | 0 |
| CAM-001 | 2 | 5 | 1 |

Multi-species: add `species`. A surveyed site-week with no row for a species is
treated as "not detected", and the report notes how many. Integrated: add a
method column (`source_col`); each method can have a different number of visits
and its own survey covariates (leave other methods' columns empty).

`sites.csv`: `station, utm_e, utm_n, canopy_cover_pct, ...`. Coordinates must be
projected (UTM, Albers...), not latitude/longitude. Set `coord_units: m` or `km`.
Each site needs distinct coordinates.

`prediction` (optional): grid with `x, y` columns (same names as `x_col`/`y_col`)
and the same site covariates. It is used for occupancy, richness and spatial-effect
maps. Keep it to 10,000 cells or fewer.

## Config

```yaml
data:
  detections: fisher_detections.csv
  sites: camera_stations.csv
  site_col: station          # default site
  survey_col: week           # default survey (ordered within site)
  detection_col: detected    # default detected
  species_col: species       # -> multi-species
  source_col: method         # -> integrated
  x_col: utm_e
  y_col: utm_n
  coord_units: m             # m (default) or km
  prediction: grid.csv       # optional
formulas:
  psi: "~ canopy_cover_pct + elevation_m"
  p: "~ days_active"         # integrated: a map {method: formula}
settings:
  model: spPGOcc             # auto | spPGOcc | PGOcc | msPGOcc | sfMsPGOcc | intPGOcc | spIntPGOcc
  n_iter: 50000              # iterations per chain (spatial: run in batches of batch_length)
  n_burn: 20000
  n_thin: 20
  n_chains: 3
  n_neighbors: 15            # NNGP neighbours
  n_factors: 1               # sfMsPGOcc spatial factors
  range_prior_km: [lo, hi]   # optional bounds on the effective range
  gof: true                  # posterior predictive check (Freeman-Tukey)
  converge_on: all           # all | core (excludes range/variance and species-level rows)
  seed: 1
figures:
  labels: {canopy_cover_pct: "Canopy cover (%)"}
  journal: esa
  include: [occupancy_map, effect_canopy_cover_pct, coefficients]
```

Examples: `example/config.yaml` (fisher, spPGOcc, with prediction grid),
`example/config_multispecies.yaml` (12 mesocarnivores, sfMsPGOcc),
`example/config_integrated.yaml` (hermit warbler point counts + ARUs, intPGOcc).

## Outputs

- An estimates table with R-hat and bulk ESS for every estimand.
- Plain-language interpretation: odds ratios per 1 SD and per round unit,
  predicted occupancy at plus and minus 1 SD, and community summaries.
- Publication figures (see the Figures section below).

## Figures

`core/run.R fit` exports these to `<out>/figures/` (PDF + PNG, plus `gallery.html`);
`core/run.R figures occu-spatial --fit <out>/occu-spatial_fit.rds --config <cfg>` redraws
them without refitting. Bands are `ci` (default 95%) and `inner_ci` (50%) credible
intervals. Covariate axes are in original units (centre/scale recorded in `data$scaling`).

| id | variants | what it shows |
|---|---|---|
| `effect_<cov>` | all | Occupancy probability vs each occupancy covariate (original units), other covariates at their means and spatial effect at 0, with 95%/50% bands and a rug of surveyed sites. Multi-species: community-average curve with bands, plus one thin grey posterior-mean curve per species. |
| `detection_<cov>[_<source>]` | single, integrated | Per-survey detection probability vs each detection covariate, with bands and a rug of surveys. |
| `coefficients` | all | Forest plot of occupancy and detection slopes (logit scale, per 1 SD); thin = 95%, thick = 50% interval. Community models show community means. |
| `occupancy_map` | single, integrated | Two panels: posterior mean occupancy (A) and posterior SD (B), on the prediction grid if `data: prediction:` is given, otherwise at surveyed sites; points = surveyed sites. |
| `spatial_effect_map` | spPGOcc, sfMsPGOcc | Posterior mean spatial random effect (logit scale, diverging scale centred on 0). For sfMsPGOcc it shows the first latent spatial factor at surveyed sites. |
| `community_effects` | multi | Species-level occupancy slopes with the community mean (red diamond), 95% intervals. With several occupancy covariates: `community_effects_<cov>`. |
| `richness_map` | multi | Estimated species richness (posterior mean, A) and its SD (B). |
| `range` | spatial | Posterior of the effective spatial range (km) against the prior implied by the phi prior (dashed). One panel per spatial factor in sfMsPGOcc. |
| `gof` | `gof: true` | Posterior predictive check (Freeman-Tukey), replicated vs observed discrepancy, with the Bayesian p-value in the caption. |
| `trace` | all | MCMC traces per chain for the main estimands. |

`community_effects` and `richness_map` exist only for multi-species models, and
`spatial_effect_map` and `range` only for spatial ones. Relabel anything through
`figures: labels:`. The keys are covariate names, species codes, parameter names
such as `psi_beta_<cov>`, and the axis and legend keys `occupancy`, `detection`,
`easting`, `northing`, `occupancy_map_fill`, `occupancy_sd_fill`, `spatial_effect`,
`richness`, `richness_sd`, `range_x`, `coefficients_x`, `gof_x`, `gof_y`,
`community_mean`, `species`, `occupancy_process` and `detection_process`.

## Simulation and power

`simulate()` implements the equations above directly. The GP is drawn by
Cholesky decomposition of the exact exponential covariance, not with
spOccupancy::simOcc, which is used only as a cross-check.

- Design knobs: `type` (single, multi, integrated), `n_sites`, `n_surveys`,
  `extent_km` (square study area), `prop_missing`, `cov_spatial_frac` and
  `cov_range_km` (how spatially smooth the simulated covariate is), `n_species`,
  `n_factors`, `loading_sd`.
- Integrated: `sources`, `n_sites_source`, `n_overlap`, `n_surveys_source`.
- Parameter knobs: `occupancy_prob`, `psi_beta_<cov>`, `detection_prob`,
  `p_beta_<cov>`, `spatial_variance`, `effective_range_km`; `comm_*` for
  communities, `factor_range_km`; `detection_prob_<source>` and
  `p_beta_<cov>_<source>`.

Power for an occupancy slope of 0.8 per SD with 4 surveys and p = 0.4: see
VALIDATION.md (`power` over `n_sites`).

## Pitfalls

- **Spatial range and variance are weakly identified** from 0/1 data. Their
  intervals are wide and depend on the prior. Do not interpret a range
  estimate as a biological dispersal or home-range scale. The `range` figure
  shows how much the posterior moved from the prior.
- **Package-default phi prior.** spOccupancy's default lower bound uses the
  *minimum* inter-site distance. With irregular sites, that puts most of the
  prior mass on ranges shorter than the site spacing, and in our test the fitted
  range collapsed to 0.1 km when the truth was 6 km. This module uses the 5%
  distance quantile instead.
- Spatial fits take minutes. For 200 sites, 50,000 iterations x 3 chains takes
  about 4-5 min, and sfMsPGOcc with 15 species is longer. Increase `n_iter` if
  R-hat > 1.01.
- `occupancy_prob` mixes slowly in spatial models, because the intercept trades
  off against the mean of w. Check its ESS.
- In sfMsPGOcc, the first `n_factors` species (alphabetical) carry the fixed
  diagonal loadings. Mixing is better when those are common, widespread species.
  Reorder species labels if chains mix poorly.
- Richness counts only the species in the data. It is not total community size,
  because no data augmentation is done.
- **Small or low-detection surveys overstate effects** (VALIDATION.md).
  - With about 100-120 sites and per-survey detection of 0.25-0.3, posterior
    means of occupancy slopes and baseline occupancy came out 10-28% too far
    from zero. The 95% intervals still covered the truth at about the nominal rate.
  - At 200 sites x 4 surveys with p = 0.4, bias was under 4%.
  - Run `recover` at the researcher's own design before reporting effect sizes.
- **The spatial model on data without residual clustering** gave occupancy
  slopes about 10% larger than PGOcc (bias +21% vs +8%). The IG(2,1) prior
  keeps the spatial variance away from zero. If the spatial variance
  interval matches its prior (the report says so), also fit `model: PGOcc`
  and report both.
- **Community SDs** (between-species variation in a coefficient) were
  overestimated by 15-30% with 10-15 species. Community means were unbiased.
- **sfMsPGOcc convergence.** At fast settings (15,000 iterations), 24% of
  simulated fits did not reach R-hat < 1.05. Use the real-data settings
  (50,000 iterations) or more, and check R-hat.

## References

- Doser JW, Finley AO, Kery M, Zipkin EF (2022) spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. *Methods in Ecology and Evolution* 13:1670-1678. doi:10.1111/2041-210X.13897
- Doser JW, Finley AO, Banerjee S (2023) Joint species distribution models with imperfect detection for high-dimensional spatial data. *Ecology* 104(9):e4137.
- Doser JW, Leuenberger W, Sillett TS, Hallworth MT, Zipkin EF (2022) Integrated community occupancy models: a framework to assess occurrence and biodiversity dynamics using multiple data sources. *Methods in Ecology and Evolution* 13. doi:10.1111/2041-210X.13811
- Datta A, Banerjee S, Finley AO, Gelfand AE (2016) Hierarchical nearest-neighbor Gaussian process models for large geostatistical datasets. *JASA* 111:800-812.
- Polson NG, Scott JG, Windle J (2013) Bayesian inference for logistic models using Polya-Gamma latent variables. *JASA* 108:1339-1349.
- Dorazio RM, Royle JA (2005) Estimating size and composition of biological communities by modeling the occurrence of species. *JASA* 100:389-398.
- MacKenzie DI et al. (2002) Estimating site occupancy rates when detection probabilities are less than one. *Ecology* 83:2248-2255.
- Banerjee S, Carlin BP, Gelfand AE (2014) *Hierarchical Modeling and Analysis for Spatial Data*, 2nd ed. CRC Press.
- spOccupancy vignettes (model fitting; factor models), https://doserlab.com/files/spoccupancy-web/
