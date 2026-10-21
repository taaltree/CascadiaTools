---
id: occu-single
title: Single-species occupancy (static single-season and dynamic multi-season)
when_to_use: Detection/non-detection surveys repeated at sites (camera traps, call surveys, eDNA replicates, point counts reduced to 0/1); question is where a species occurs and what is associated with occurrence, detection, colonisation or local extinction, accounting for imperfect detection. Static = one season; dynamic = several seasons/years at the same sites.
not_for: Counts or abundance (use an N-mixture module), several interacting species (occu-rota / occu-waddle), spatially autocorrelated occupancy (spOccupancy), false-positive detections, a single survey per site (detection cannot be separated from occupancy), designs with very low detection and few surveys (see Known limits in VALIDATION.md).
data_needed: detection history (long table site[, season], survey, detected 0/1; or wide site x survey), site covariates, optional survey covariates and site-by-season covariates
engine: unmarked (occu, colext; maximum likelihood, default) or nimble (Bayesian, settings$engine = nimble)
---

# Single-species occupancy

## Question it answers

*Where does the species occur, and which habitat or survey conditions are
associated with occurrence and with detecting it?* A species can be present at
a site but missed on every visit; the naive "proportion of sites where it was
seen" is then biased low, and covariates that affect detection can look like
habitat effects. Repeat surveys let the model separate the two.

- **static** (MacKenzie et al. 2002): occupancy probability (psi) and per-survey
  detection probability (p), both with covariates; finite-sample proportion of
  the surveyed sites occupied.
- **dynamic** (MacKenzie et al. 2003): initial occupancy (psi1), colonisation
  (gamma: unoccupied -> occupied between seasons), local extinction (epsilon:
  occupied -> unoccupied), detection (p); proportion of surveyed sites occupied
  in each season and its change from first to last season.

## Model

Static, site i, survey j:

    z_i ~ Bernoulli(psi_i),           logit(psi_i) = a0 + a' x_i
    y_ij | z_i ~ Bernoulli(z_i p_ij),  logit(p_ij)  = b0 + b' w_ij

Dynamic, season t = 1..T:

    z_i1 ~ Bernoulli(psi_i)
    z_i,t+1 | z_it ~ Bernoulli(z_it (1 - eps_it) + (1 - z_it) gamma_it)
    logit(gamma_it) = c0 + c' v_it,  logit(eps_it) = e0 + e' v_it   (v_it: season-t value of a site or site-by-season covariate)
    y_ijt | z_it ~ Bernoulli(z_it p_ijt)

Unsurveyed occasions (NA) contribute nothing. Each numeric model term is
standardised (mean 0, SD 1 over the rows it applies to) after formula
transforms, so:

- the intercept is reported on the probability scale at average covariate
  values (factors at their reference level): `occupancy_prob`,
  `initial_occupancy_prob`, `colonisation_prob`, `extinction_prob`,
  `detection_prob`;
- slopes are on the logit scale **per 1 SD** of the term: `psi_beta_<term>`,
  `p_beta_<term>`, `col_beta_<term>`, `ext_beta_<term>` (e.g.
  `psi_beta_forest_cover_pct`, `p_beta_observerJK` for factor level JK vs the
  reference level, `psi_beta_log_wetland_area_ha` for `log(wetland_area_ha)`).

Finite-sample quantities (Royle & Kery 2007; MacKenzie et al. 2018):
`prop_sites_occupied` (static), `prop_sites_occupied_t1..tT` and
`prop_occupied_change` = season T minus season 1 (dynamic). These describe the
sites actually surveyed, not the wider landscape.

Engines:

- **unmarked** (default): `occu()` / `colext()` (Fiske & Chandler 2011),
  Wald 95% CIs on the logit scale, probabilities back-transformed. The
  finite-sample proportions use each site's conditional probability of
  occupancy given its detection history (forward-backward), with parameter
  uncertainty propagated by drawing from the asymptotic normal of the MLEs.
- **nimble**: same model with latent z sampled (Royle & Kery 2007); priors
  Logistic(0, 1) on intercepts (flat on the probability scale) and
  Normal(0, SD = `prior_sd` = 2) on slopes; slice samplers on coefficients.
- **naive**: logistic regression of "ever detected"; a negative control only,
  never for inference.

Goodness of fit: MacKenzie & Bailey (2004) Pearson chi-square over detection
histories, sites grouped into cohorts by missing-survey pattern; dynamic models
compute it within each season using model-implied occupancy and sum the seasons
(as `AICcmodavg::mb.gof.test` does). unmarked: parametric bootstrap
(`gof_nboot` refits), c-hat = observed / bootstrap mean. nimble: posterior
predictive p-value of the same statistic.

## Assumptions (and what breaks if they fail)

1. **Closure within a season**: a site's occupancy does not change between
   surveys of the same season. Violations (animals moving in and out, surveys
   spread over migration) make psi mean "use" rather than occupancy and
   underestimate p.
2. **No false positives**: a detection means the species was there.
   Misidentifications inflate psi; use a false-positive model instead.
3. **Independence** between sites and surveys given covariates. Clustered or
   spatially autocorrelated sites make intervals too narrow (GOF c-hat > 1 is
   one symptom).
4. **Detection heterogeneity is modelled**: unmodelled variation in p among
   sites (abundance, observer) biases psi low. It can show up as c-hat > 1,
   but the MacKenzie-Bailey test has low power for it (in VALIDATION.md it
   rejected in only ~10% of data sets where psi was biased -23%), so a
   passing GOF test does not rule it out.
5. **Covariate effects are linear on the logit scale** unless you add terms
   (`I(x^2)`, `log(x)`).
6. Dynamic: transitions are first-order Markov; the covariate value in season
   t drives the change from t to t+1.

## Data format

**Detections, long format (recommended)**: one row per survey actually done.
Missing visits are simply absent (or `detected` blank).

| site | survey | julian_day | observer | detected |
|------|--------|------------|----------|----------|
| STAND-001 | 1 | 127 | MR | 1 |
| STAND-001 | 2 | 138 | AB | 1 |
| STAND-002 | 1 | 128 | AB | 0 |

- `site`: site ID (text). `survey`: visit number 1, 2, 3... within the site
  (and within the season for dynamic). Dates are fine as extra covariate
  columns but visits must be numbered.
- `detected`: 1/0 (counts must be converted to 1 = detected). Text such as
  `yes`/`no` or `Y` is reported as a data problem, not guessed; blank = survey not done.
- Survey covariates: any extra column listed in `obs_covs` (numeric or text;
  text becomes a factor). Must be filled in for every survey done.
- Blank cells are missing values in every file and column (text columns
  included), so a blank observer or habitat is reported as missing, never
  treated as a category. Numbers stored with units, `%` or thousands commas
  (`39.8%`, `1,200`) and text covariates with very many distinct values (IDs,
  dates) are reported as data problems rather than fitted as factors.
- Dynamic: add a `season` column. Numbers (years) are ordered numerically;
  labels with one number in a common pattern (`Y1`...`Y12`, `season 3`) by that
  number; any other labels (`spring`, `2021-22`) need `data: season_order:` in
  the config (earliest first), otherwise reading stops with that instruction.

**Detections, wide format**: one row per site, detection columns in survey
order (season-major for dynamic, with `surveys_per_season`), blank = not
surveyed. Survey covariates as matching wide files in `obs_cov_files`.

**Site covariates** (`site_covs`): one row per site, `site` + covariate columns.

**Site-by-season covariates** (dynamic, `season_covs`): `site`, `season`,
covariates. The season-t value is used for the change from t to t+1 (the last
season's value is not used).

## Config

```yaml
data:
  model: static | dynamic
  detections: detections.csv        # path relative to the config file
  format: long | wide               # default long
  site_col: site                    # default "site"
  survey_col: survey                # long: visit number column (default "survey")
  season_col: season                # long, dynamic: season column (default "season")
  season_order: [spring, summer, autumn]   # dynamic: time order of text season labels (see above)
  detection_col: detected           # long: 0/1 column (default "detected")
  obs_covs: [effort, observer]      # long: survey-level covariate columns
  detection_cols: [v1, v2, v3]      # wide: detection columns (default: all but site)
  surveys_per_season: 3             # wide, dynamic
  obs_cov_files: {effort: effort_wide.csv}   # wide: survey covariates
  site_covs: sites.csv
  site_cov_cols: [forest, elev]     # optional subset
  season_covs: site_seasons.csv     # dynamic
formulas:                           # one-sided R formulas; default "~ 1"
  psi: "~ forest + elev"            # (initial) occupancy: site covariates
  p:   "~ effort + observer"        # detection: site, season or survey covariates
  col: "~ forest"                   # dynamic: site or site-by-season covariates
  ext: "~ drought"                  # dynamic
settings:
  engine: unmarked                  # unmarked | nimble | naive (control only)
  gof_nboot: 200                    # unmarked: bootstrap refits for MacKenzie-Bailey (0 = skip). nimble: >0 turns on the PPC
  gof_ndraws: 300                   # nimble: posterior draws used in the PPC
  n_draws: 2000                     # unmarked: normal draws for intervals of derived quantities / plots
  n_prop_draws: 1000                # unmarked: draws for finite-sample proportion occupied
  standardize: true                 # standardise numeric model terms (recommended)
  niter: 12000                      # nimble: iterations per chain
  nburnin: 2000
  nchains: 3
  thin: 1
  prior_sd: 2                       # nimble: SD of Normal prior on slopes (logit per SD)
  slice: true                       # nimble: slice samplers on coefficients
  seed: 1
```

Avoid `poly()` in formulas (use `I(x^2)`); formulas must keep the intercept.

## Outputs

Estimates table (95% intervals): probabilities at average covariates, logit
slopes per 1 SD, finite-sample proportions occupied (and change, dynamic).
Diagnostics: convergence (unmarked: optimiser + invertible Hessian + no
boundary estimate, i.e. no |logit coefficient| or SE > 8; nimble: R-hat and
bulk ESS), MacKenzie-Bailey GOF. Report plots: conditional effects of every
covariate on occupancy, detection, colonisation and extinction on the
probability scale (other covariates at their means, factors at reference,
95% band, rug of observed values); dynamic: proportion of sites occupied by
season with naive values. Interpretation: plain-language effects from the
10th to the 90th percentile of each covariate, odds ratios per 1 SD and per a
round unit, cumulative detection, GOF verdict.

## Figures

`core/run.R fit` exports the publication suite to `<out>/figures/` (PDF + PNG,
`gallery.html`, `figures_index.csv` with captions); `core/run.R figures
occu-single --fit <out>/occu-single_fit.rds --config <cfg>` redraws without
refitting. Code: `R/figures.R` (`occ_figures`). All axes are probabilities on
the 0-1 scale and covariates are in the units of your files.

| id | what it shows |
|---|---|
| `occupancy_<cov>` | (initial) occupancy probability vs each occupancy covariate; the other covariates of that submodel held at their mean / reference level, with the values stated in the caption; band at `ci` (plus `inner_ci` for nimble); rug of observed values; factors as points with intervals |
| `detection_<cov>` | per-survey detection probability vs each detection covariate (same layout) |
| `col_ext_<cov>` | dynamic: colonisation and/or extinction probability vs a covariate; one colour + line type per process when the covariate is in both |
| `coefficients` | odds ratios (log scale, per 1 SD; factor levels vs reference) of every slope, one panel per process, dashed line at 1; caption lists the SDs |
| `detection_by_survey` | per-survey detection probability by visit number (one panel per season for dynamic), averaged over the sites surveyed on that visit |
| `occupancy_map` | only when the site file has coordinates (`x/y`, `lon/lat`, `easting/northing`, or `figures: coords: [col1, col2]`): P(site occupied given its detection history), circles = detected, triangles = never detected; one panel per season for dynamic |
| `trajectory` | dynamic: proportion of surveyed sites occupied per season with interval, naive proportion as points |
| `gof` | unmarked: histogram of the MacKenzie-Bailey bootstrap chi-square with the observed value; nimble: replicated minus observed statistic over posterior draws, line at 0 |
| `trace`, `posterior` | nimble only: MCMC traces per chain and posterior densities of the reported parameters |

Colours follow the process (occupancy, detection, colonisation, extinction)
in every figure (with `greyscale: true`, processes are told apart by panel, line type or axis, never by grey shade alone). Captions use the axis titles as drawn, including `per_figure: axis:` overrides. Rename anything via `figures: labels:` keys: covariate
column names (`forest_cover_pct`), factor levels (`observer_JK`), parameter
names (`psi_beta_forest_cover_pct`, `occupancy_prob`), processes
(`process_psi`, `process_p`, `process_col`, `process_ext`), axes
(`occupancy`, `detection`, `colonisation`, `extinction`, `col_ext`, `survey`,
`season`, `trajectory`, `coefficients_x`, `gof_x`, `occupancy_map_fill`).
Use `units:` only for covariates that were already standardised in your files.
The naive engine (negative control) produces no figures. Figures for fits saved
before this suite existed have no `gof` (the replicate statistics were not
stored); refit to get it.

## Simulation and power

`simulate()` knobs (`design`): `model` (static/dynamic), `n_sites`,
`n_surveys` (per season), `n_seasons` (dynamic), `prop_missing` (share of
surveys not done, at random), `p_site_sd` (SD of unmodelled site-level
detection heterogeneity on the logit scale; 0 = none; use it to test GOF and
robustness). `params`: `occupancy_prob`, `detection_prob`,
`psi_beta_forest` (site covariate), `p_beta_effort` (survey covariate);
dynamic: `initial_occupancy_prob`, `colonisation_prob`, `extinction_prob`,
`col_beta_forest` (site covariate), `ext_beta_drought` (site-by-season
covariate). Covariates are standard normal. Power parameters:
`psi_beta_forest`, `p_beta_effort`, `col_beta_forest`, `ext_beta_drought`,
`prop_occupied_change`.

Sensible ranges: detection 0.2-0.8 per survey, 2-8 surveys, 50-500 sites.
Check cumulative detection p* = 1 - (1 - p)^J: below about 0.5 occupancy
estimates become unstable (VALIDATION.md).

    scripts/Rscript core/run.R recover occu-single --config my_design.yaml --nsim 200 --cores 16
    scripts/Rscript core/run.R power occu-single --vary n_sites=50,100,200,400 --nsim 100 --cores 16

## Pitfalls

- **Few surveys and low detection**: estimates hit the boundary (psi = 1,
  SE huge) or are badly biased. The module flags boundary fits as not
  converged. With 2 surveys and p = 0.2 recovery failed (VALIDATION.md).
- **Small samples inflate slopes**: maximum-likelihood logit slopes are
  biased away from zero at ~100 sites (see VALIDATION.md for size).
- **Detection covariates are not habitat effects.** A survey covariate in p
  says how easily the species is found, not where it lives.
- **Same covariate in psi and p** is allowed but needs good data; a strong
  correlation between occupancy and detection covariates weakens both.
- **Finite-sample vs population occupancy**: `prop_sites_occupied` is about
  the surveyed sites; `occupancy_prob` is the model's probability at average
  covariates.
- **GOF c-hat > 1** (overdispersion): add covariates, or inflate SEs by
  sqrt(c-hat). c-hat >> 4 suggests the model structure is wrong.
- Naive engine is a control: never report it.

## References

- MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA, Langtimm CA (2002) Estimating site occupancy rates when detection probabilities are less than one. Ecology 83:2248-2255.
- MacKenzie DI, Nichols JD, Hines JE, Knutson MG, Franklin AB (2003) Estimating site occupancy, colonization, and local extinction when a species is detected imperfectly. Ecology 84:2200-2207.
- MacKenzie DI, Bailey LL (2004) Assessing the fit of site-occupancy models. JABES 9:300-318.
- Fiske I, Chandler R (2011) unmarked: an R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23.
- Royle JA, Kery M (2007) A Bayesian state-space formulation of dynamic occupancy models. Ecology 88:1813-1823.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26:403-413.
- Guillera-Arroita G, Ridout MS, Morgan BJT (2010) Design of occupancy studies with imperfect detection. Methods in Ecology and Evolution 1:131-139.
- MacKenzie DI, Nichols JD, Royle JA, Pollock KH, Bailey LL, Hines JE (2018) Occupancy Estimation and Modeling, 2nd ed. Academic Press.
- unmarked colext vignette (CRAN): parameterisation and projected vs smoothed occupancy.
- AICcmodavg `mb.gof.test` documentation: cohorts by missing-value pattern; per-season chi-square summed for dynamic models.
