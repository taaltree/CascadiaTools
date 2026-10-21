---
id: glm
title: GLM / GLMM with biological interpretation
when_to_use: One response measured on sampling units (sites, plots, nests, animals, quadrats, surveys) and the question is how it changes with habitat, environmental or treatment predictors, possibly with repeated measures within sites/individuals/years (random intercepts). Counts, presence/absence (0/1), successes out of trials, proportions in (0,1), positive continuous measurements (mass, length, biomass) and normal measurements.
not_for: Imperfect detection that must be separated from occurrence or abundance (use occupancy / N-mixture modules); spatial or temporal autocorrelation that must be modelled explicitly; random slopes, nested (1|a/b) syntax, splines/GAMs, ordinal or multinomial responses; causal effect estimation (use the causal/SEM modules); proportions containing exact 0s or 1s (transform first or use binomial with trials).
data_needed: one CSV, one row per observation, with a response column, predictor columns, optional grouping column(s) (site, individual, year), optional effort column (offset) and, for binomial counts, a trials column
engine: glmmTMB (REML for models with random intercepts, ML otherwise; Laplace approximation); stats::glm / MASS::glm.nb fallback for fixed-effects-only models
---

# GLM / GLMM made easy

## Question it answers

"How does *Y* change with *X*?", where *Y* is something you measured on each
sampling unit:

| Your response | Family | What effects are reported as |
|---|---|---|
| Measurement roughly symmetric around a mean (length, temperature, log-mass) | `gaussian` | difference in mean per SD / per unit |
| Presence/absence, survived/died, used/unused (0/1) | `binomial` | odds ratio; probability on plots |
| Successes out of a known number of trials (eggs hatched of eggs laid) | `binomial` + `trials:` | odds ratio; probability |
| Counts (animals, nests, seedlings) | `poisson` (only if no extra spread) or `nbinom2` (usual choice) | rate ratio = multiplicative change in expected count |
| Counts with many more zeros than expected | `zipoisson`, `zinb` | rate ratio + share of excess zeros |
| Proportion strictly between 0 and 1 (percent cover / 100, diet proportion) | `beta` | odds ratio of the mean proportion |
| Positive continuous, spread grows with the mean (body mass, biomass, distance moved) | `Gamma` (log link) | ratio of means |

Every effect is reported (i) on the link scale (`beta_*`), (ii) as an
odds / rate ratio per 1 SD of the predictor and per a biologically
meaningful unit you choose (`units:`), and (iii) on the response scale in
plain-language sentences and conditional-effects plots.

## Model

For observation *i* in group *j*:

    g(mu_ij) = b0 + b1 * x1_ij + ... + u_j + log(effort_ij)          u_j ~ Normal(0, sd_group^2)
    y_ij ~ Family(mu_ij, dispersion)

* `g` is the link: identity (gaussian), logit (binomial, beta), log (counts, Gamma).
* Continuous predictors are standardised by default, x' = (x - mean) / SD, so
  `beta_x` is the change per 1 SD and `beta_intercept` is the link-scale value
  at the mean of every continuous predictor and the reference level of every
  factor (Schielzeth 2010). Per-unit effects are rescaled exactly:
  beta_per_unit = beta_per_SD x unit / SD.
* Random intercepts `u_j` (one per `(1|group)` term) capture extra similarity
  among observations from the same site, individual, year, etc. (Bolker et al. 2009).
* Variance functions (glmmTMB parameterisation; Brooks et al. 2017):
  gaussian Var = sigma^2; Poisson Var = mu; NB2 Var = mu + mu^2 / theta;
  beta Var = mu(1 - mu) / (1 + phi) (Ferrari & Cribari-Neto 2004);
  Gamma Var = (CV x mu)^2; zero-inflated: y = 0 with probability pi, otherwise from the count model,
  so E[y] = (1 - pi) mu.
* Estimation: models with random intercepts are fitted by REML (restricted
  maximum likelihood: the fixed effects are integrated out together with the
  random effects by the Laplace approximation, started from the ML solution);
  models without random effects by maximum likelihood (TMB). For non-Gaussian
  families the REML fit supplies only the random-effect SDs and dispersion; the
  fixed effects are then re-estimated by maximum likelihood with those held at
  their REML values (the GLMM analogue of LMM REML). The REML mode of the
  fixed effects itself is biased for non-Gaussian data (NB GLMM intercept
  about -0.05 at site SD 0.7, not shrinking with more sites; found and fixed in
  independent verification, see VALIDATION.md). REML removes most
  of the downward bias of random-effect SDs that ML shows with ~30 groups
  (ML: -10 to -18% in validation; REML: within +/-3%) and with it the
  slightly too narrow intervals for group-level predictors. Set
  `settings: {REML: false}` for ML (needed only for likelihood-ratio tests
  or AIC comparisons of models with different fixed effects). Intervals: Wald, on the link scale (coefficients) or
  log scale (SDs, dispersion), then back-transformed. Conditional-effect bands
  use 2000 draws from the multivariate-normal sampling distribution of the estimates.
* Checks: DHARMa simulated quantile residuals (Dunn & Smyth 1996; Hartig 2026)
  - uniformity (KS), dispersion, outliers and (counts) zero-inflation tests;
  R2 of Nakagawa et al. (2017) for mixed models (marginal = predictors only,
  conditional = predictors + groups) via `performance`.

## Assumptions (and what breaks)

1. **Right family / mean-variance relation.** Poisson counts that are really
   overdispersed give intervals that are far too narrow (validated below:
   coverage of effects falls to ~0.6). Always look at the DHARMa dispersion
   test; switch to `nbinom2` when it flags.
2. **Independence given the random effects.** Repeated visits to a site,
   repeated measures on an animal, quadrats on a transect are NOT independent:
   include `(1|site)` etc. Ignoring it makes site-level covariate intervals too
   narrow (validated below). Spatial/temporal autocorrelation beyond the
   grouping is not modelled; check residuals against space/time yourself.
3. **Linearity on the link scale.** Effects are straight lines on the
   log/logit scale (curves on the response scale). Add `I(x^2)` or `poly(x, 2)`
   if the relationship should peak; the plots will show the fitted shape.
4. **Random intercepts need enough groups.** With < 5 levels the SD cannot be
   estimated; use the grouping as a fixed factor. With 10-15 groups and few
   observations per group the SD is poorly estimated and intervals for
   group-level predictors can be too narrow (see VALIDATION.md, hard designs).
5. **No measurement error in predictors; no unmeasured confounding** if you
   want to read effects causally. Results are associations.

## Data format

One CSV, one row per observation (survey, nest, quadrat, animal-capture...):

```
station,visit,date,observer,survey_min,forest_type,canopy_cover_pct,elevation_m,wood_thrush
WT-01,1,2025-05-20,JL,5,mixed,66,725,0
WT-01,2,2025-05-29,JL,10,mixed,66,725,3
WT-02,1,2025-05-21,AK,10,conifer,48,312,0
```

* Response: numeric. 0/1 for presence; whole numbers for counts; strictly
  between 0 and 1 for beta (divide percentages by 100); > 0 for Gamma.
* Predictors: numeric columns are continuous; text columns are factors (the
  first level alphabetically is the reference unless you set `reference:`).
* Grouping columns (site, individual, year) can be text or numbers.
* Extra columns (dates, observers, notes) are ignored unless used in the formula.
* Missing values in used columns stop the analysis with a list of rows;
  set `drop_missing: true` to drop them instead (rows are reported).
* Column names with spaces/symbols are converted by R (spaces become dots).

## Config (every key)

```yaml
data:
  file: wood_thrush_counts.csv        # path relative to this config file
  formula: wood_thrush ~ canopy_cover_pct + elevation_m + forest_type + (1|station)
  family: nbinom2                     # gaussian | binomial | poisson | nbinom2 | beta | Gamma | zipoisson | zinb
  link: log                           # optional; defaults identity / logit / log as in the table above
  offset: survey_min                  # optional effort column (> 0); log(effort) enters as an offset (log link only)
  offset_per: 10                      # optional; report expected values per this much effort (default 1)
  trials: eggs_laid                   # binomial successes/trials only: column with the number of trials
  scale: all                          # all (default) | none | [list of columns] to standardise
  units: {canopy_cover_pct: 10, elevation_m: 100}   # also report effects per these changes
  unit_text: {canopy_cover_pct: percentage points, elevation_m: m}
  labels: {canopy_cover_pct: canopy cover (%)}      # names used in sentences and plot axes
  reference: {forest_type: deciduous}               # reference level for factors
  response_label: Wood Thrush count                 # used in sentences ("expected Wood Thrush count")
  zi: "~1"                            # zero-inflation formula for zipoisson/zinb (default ~1 = constant share)
  drop_missing: false
settings:
  diagnostics: true                   # DHARMa residual checks + R2 (turn off only for simulation)
  REML: true                          # REML for models with random effects (false = ML)
  engine: glmmTMB                     # or glm (fixed-effects-only gaussian/binomial/poisson/nbinom2/Gamma)
  n_sim_resid: 250                    # DHARMa simulations
  n_draws: 2000                       # draws for plot bands / zero-inflated mean
  seed: 1
```

Formula rules: one response column on the left; `+` for main effects, `a:b`
or `a*b` for interactions, `I(x^2)` / `poly(x,2)` for curvature,
`(1|group)` for each random intercept (several allowed: `(1|site) + (1|year)`).
Not supported: random slopes `(x|site)`, nested shorthand `(1|a/b)` (make one
column per level and use `(1|a) + (1|b)`), `offset()` inside the formula
(use `offset:`), `cbind()` on the left (use `trials:`).

## Worked examples (one per family)

Run any of these with
`scripts/Rscript core/run.R fit glm --config <config.yaml> --out projects/<name>/results`.

**Negative binomial counts with effort and sites** (`example/config.yaml`, data
`example/wood_thrush_counts.csv`): Wood Thrush counts at 45 stations x 4
visits, survey length 5 or 10 min. Output reads like: "Each 10-percentage-point
increase in canopy cover is associated with a 39% higher expected Wood Thrush
count per 10-min survey (rate ratio 1.39, 95% CI 1.19 to 1.62)".

**Binomial successes out of trials, with colonies** (`example/config_binomial.yaml`,
`example/nest_hatching.csv`): eggs hatched out of eggs laid; "Per 50 m further
from the forest edge the odds of an egg hatching are multiplied by 1.66 (1.31 to 2.10)".

**Beta proportions, with transects** (`example/config_beta.yaml`,
`example/seagrass_cover.csv`): seagrass cover proportion vs depth and turbidity.

**Presence/absence (0/1)**:
```yaml
data: {file: detections.csv, formula: present ~ shrub_cover + road_dist_km + (1|site), family: binomial,
       units: {shrub_cover: 10, road_dist_km: 1}, response_label: hare presence}
```
**Poisson counts without groups** (use only if DHARMa shows no overdispersion):
```yaml
data: {file: nests.csv, formula: n_nests ~ cliff_height_m + aspect, family: poisson, offset: area_km2}
```
**Zero-inflated counts** (many sites where the species can never occur):
```yaml
data: {file: mussels.csv, formula: mussels ~ substrate + flow_cm_s + (1|river), family: zinb, offset: quadrats}
```
**Gamma (positive, right-skewed)**:
```yaml
data: {file: fish_mass.csv, formula: mass_g ~ length_class + temp_c + (1|lake), family: Gamma, units: {temp_c: 1}}
```
**Gaussian (normal) mixed model**:
```yaml
data: {file: wing.csv, formula: wing_mm ~ age + sex + elevation_m + (1|population), family: gaussian, units: {elevation_m: 500}}
```
**Interaction** (does the canopy effect differ by forest type?):
```yaml
data: {file: wood_thrush_counts.csv, formula: wood_thrush ~ canopy_cover_pct * forest_type + (1|station), family: nbinom2, offset: survey_min}
```
The canopy plot then shows one curve per forest type; the main effect of
canopy is its effect in the reference forest type.

## Outputs

* `glm_estimates.csv`: one row per estimand with estimate, SE (`sd`), 95% CI,
  the `scale` and what it is `per`:
  * `beta_<term>` link-scale coefficients (`beta_intercept`, `beta_canopy_cover_pct`, `beta_forest_type_conifer`, `beta_a_x_b` for interactions);
  * `odds_ratio_<term>` (logit), `rate_ratio_<term>` (counts), `ratio_of_means_<term>` (Gamma / log-link gaussian), per 1 SD (continuous) or vs reference (factor levels);
  * `<ratio>_<term>_per_<u>` (or `beta_<term>_per_<u>` for identity link) per user unit;
  * `mean_response_at_reference`: expected response at the means / reference levels, per `offset_per` units of effort, for a typical group (random effect 0), including zero-inflation;
  * dispersion: `residual_sd` (gaussian), `nb_size_theta` (NB), `beta_precision_phi`, `gamma_cv`;
  * `zero_inflation_prob` (or `zi_beta_*` for a zi formula);
  * `re_sd_<group>`: SD of the random intercepts on the link scale.
* `glm_report.html`: methods paragraph, diagnostics (convergence, DHARMa tests, R2,
  collinearity and small-sample warnings), estimates, plain-language sentences
  for every term, figures: one conditional-effects plot per predictor
  (response scale; others at mean / reference; factors as point + 95% CI;
  interactions as several curves), an effect-size forest plot, and DHARMa
  residual plots.

## Figures

`core/run.R fit` exports the publication figure suite to `<out>/figures/`
(PDF + PNG by default, `gallery.html`, `figures_index.csv` with the
captions). Restyle without refitting:
`scripts/Rscript core/run.R figures glm --fit <out>/glm_fit.rds --config <config.yaml>`
after editing the `figures:` block (see `example/config.yaml` and
`skills/eco-figures`). Code: `figures.R` (`glm_figures`).

| id | what it shows |
|---|---|
| `effect_<column>` | One per predictor in the formula. Continuous: expected response (response scale, per `offset_per` effort, zero-inflation included) across the observed range, x axis in original units, band at `ci` from 2000 draws of the fixed effects; others at mean / reference level, random effects 0. Factors: point + interval per level (`fig_cat_effect`). Interactions: one curve per level of the moderating factor, or at the 10th / 50th / 90th percentile of a continuous moderator (colour + line type, so greyscale works). |
| `effects_panel` | All `effect_*` panels composed with A/B/C tags (double width). |
| `coefficients` | Forest plot of odds ratios (logit), rate ratios (counts) or ratios of means (log link), on a log axis with the null line at 1; differences on the identity scale. Continuous predictors per 1 SD (SD in original units shown in the row label), factor levels vs reference, interaction terms as multiplicative modifiers. |
| `residuals` | DHARMa scaled residuals: (A) QQ plot vs uniform with the KS p-value; (B) residuals vs rank-transformed prediction with binned 0.25/0.5/0.75 quantiles (should be flat on the dashed lines). Uses the stored DHARMa simulation, or runs one if `diagnostics: false`. |
| `random_effects` | Mixed models only: caterpillar plot of the random intercepts (conditional modes +/- `ci` interval from their conditional SD), as multipliers of the odds / expected value vs a typical group for logit / log links. Labels are hidden above 50 groups. |

Data shown (`show_data: true`): a rug of observed predictor values on every
continuous effect; raw points as well for proportions (binomial with trials,
beta) and continuous responses (gaussian, Gamma). Raw counts and 0/1 data are
not drawn as points: individual counts spread far beyond the expected value
and would flatten the curve (counts and Gamma means use a y axis that starts
at 0). Raw points are not adjusted for the other predictors or random effects.

Relabelling keys for `labels:`: any column name (`canopy_cover_pct`), a factor
level (`conifer` or `forest_type_conifer`), a grouping column (`station`) or a
group id (`station_WT-01`), `response` (y axis of the effect plots), a
coefficient row (`beta_<term>`, or `rate_ratio_<term>` / `odds_ratio_<term>`, e.g. `rate_ratio_forest_type_conifer`; a column-name key renames the variable inside the row text),
`coefficients_axis`, `random_effects_axis`, `qq_expected`, `qq_observed`,
`resid_predicted`. `units: {<column>: {center, scale}}` overrides the
back-transformation of a standardised predictor; by default the mean / SD
used in fitting (`fit$prep$scaling`) are applied, so axes are always in the
data's own units. Intervals are frequentist (Wald / sampling-distribution
draws), so `inner_ci` is not used.

## Simulation and power

`simulate(params, design, seed)` is generic:

* `design`: `family`, `link`, `n_groups` and `n_per_group` (or `n_groups: 0`
  and `n`), `group_name`, `covariates` (named list; each `type: normal`
  with `mean`, `sd`; `type: uniform` with `min`, `max`; or `type: factor` with
  `levels`; `level: group` makes it constant within groups), `interactions`
  (e.g. `["forest:temp"]`), `offset` (`name`, `min`, `max`; or `null`),
  `trials` (binomial: a number or `[min, max]`), `units`, `response`.
* `params`: `beta_<term>` per 1 SD of each continuous covariate (exactly as
  fitted, since the simulator standardises with the sample mean/SD), factor
  contrasts `beta_<factor>_<level>`, `re_sd`, `dispersion` (residual SD,
  NB theta, beta phi or Gamma CV according to family), `zi_prob`.
* Defaults: NB GLMM, 30 sites x 8 visits, site-level forest cover and habitat
  (3 levels), visit-level temperature, effort offset 0.5-2.

Examples:
`core/run.R power glm --vary n_groups=10,20,30,50,80 --nsim 100` (sites needed),
`core/run.R power glm --what params --vary beta_forest=0,0.1,0.2,0.3` (detectable effect),
and `--config` with `design: {family: binomial, offset: null}` to switch family.
For a site-level covariate, power is driven by the number of sites far more
than visits per site (see VALIDATION.md).

## Pitfalls

* Poisson on overdispersed counts: intervals too narrow. Use nbinom2 by default for ecological counts.
* A site-level covariate with repeated visits and no `(1|site)`: pseudoreplication, intervals too narrow.
* Proportions from counts (3 of 10 eggs) should be binomial with `trials:`, not beta.
* Percent cover with exact 0 or 100: beta cannot take 0/1; transform (Smithson & Verkuilen 2006) or use a different model.
* Correlated predictors (|r| > 0.7) are flagged: their separate effects are poorly identified.
* Rate ratios multiply: a rate ratio of 0.5 per SD over 3 SD is 0.125, not -150%.
* `mean_response_at_reference` is for a *typical* group (random effect 0); with a
  log link the average over groups is larger by exp(sd^2/2).
* Logistic coefficients beyond about +/-8 signal complete separation; intervals are then meaningless.

## References

* Brooks ME, Kristensen K, van Benthem KJ, Magnusson A, Berg CW, Nielsen A, Skaug HJ, Maechler M, Bolker BM (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9(2):378-400. doi:10.32614/RJ-2017-066
* Bolker BM, Brooks ME, Clark CJ, Geange SW, Poulsen JR, Stevens MHH, White J-SS (2009). Generalized linear mixed models: a practical guide for ecology and evolution. Trends in Ecology and Evolution 24(3):127-135.
* McCullagh P, Nelder JA (1989). Generalized Linear Models, 2nd ed. Chapman & Hall.
* Hardin JW, Hilbe JM (2007). Generalized Linear Models and Extensions, 2nd ed. Stata Press.
* Ferrari SLP, Cribari-Neto F (2004). Beta regression for modelling rates and proportions. Journal of Applied Statistics 31(7):799-815.
* Smithson M, Verkuilen J (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological Methods 11(1):54-71.
* Dunn PK, Smyth GK (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5(3):236-244.
* Hartig F (2026). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package 0.5.0. doi:10.32614/CRAN.package.DHARMa
* Nakagawa S, Johnson PCD, Schielzeth H (2017). The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface 14:20170213.
* Luedecke D, Ben-Shachar MS, Patil I, Waggoner P, Makowski D (2021). performance: an R package for assessment, comparison and testing of statistical models. Journal of Open Source Software 6(60):3139.
* Schielzeth H (2010). Simple means to improve the interpretability of regression coefficients. Methods in Ecology and Evolution 1(2):103-113.
