---
id: fisheries-growth
title: Fish growth (von Bertalanffy) and length-weight condition
when_to_use: Individual age-length data (otolith/scale ages) for one or more groups (populations, lakes, sexes, years); question is how fast fish grow, how large they get, at what age they approach maximum size, whether growth or condition differs among groups, or whether growth varies with a group-level covariate (temperature, density, year). Optional weights give the length-weight relationship and condition (Wr or Kn).
not_for: Length-frequency data without ages (needs ELEFAN/mixture methods), tag-recapture growth increments (Fabens model), back-calculated lengths from annuli (repeated measures per fish), stock assessment / mortality (use fisheries-stock), growth models other than von Bertalanffy (Gompertz, logistic, Richards are not implemented).
data_needed: one row per fish with age (years), length, group (optional), weight (optional); optional group-covariate table (one row per group)
engine: stats::nls with FSA starting values + Kimura (1980) likelihood-ratio tests (default), or nimble (Bayesian, fixed or hierarchical random group effects, covariates on log Linf and log K); stats::lm for length-weight
---

# fisheries-growth

## Question it answers

- How large do fish of this population get (asymptotic length Linf), how
  fast do they approach it (K), and at what age do they reach 90% of Linf?
- Do two or more groups (lakes, sexes, years, management units) grow
  differently, and in which parameter? (Kimura's likelihood-ratio tests, or
  posterior differences.)
- How much does growth vary among many groups (lakes or years), and is the
  variation associated with a group covariate such as temperature or year?
  (Hierarchical Bayesian model.)
- What is the length-weight relationship W = aL^b, are fish isometric
  (b = 3), and are fish in one group heavier at a given length (condition)?

## Model

**Growth (typical von Bertalanffy, von Bertalanffy 1938; Beverton & Holt 1957):**

    L(age) = Linf * (1 - exp(-K * (age - t0)))

- `error: multiplicative` (default): log L_i ~ Normal(log L(age_i), sigma).
  Spread grows with length, which is what most fish data show. The fitted
  curve is the *median* length at age.
- `error: additive`: L_i ~ Normal(L(age_i), sigma), Kimura's (1980) "method a"
  (constant variance). The curve is the *mean* length at age.

**Engines (`settings: engine`)**

1. `nls` (default): all three parameters group-specific (Kimura's general
   model), fitted by nonlinear least squares (port algorithm) on the scale
   (log Linf, log K, t0), with starting values from `FSA::findGrowthStarts`
   (the successor of `FSA::vbStarts`; Ogle 2016). Wald 95% intervals are
   back-transformed for Linf and K (so they stay positive); differences
   between groups and age at 90% Linf (t0 + ln(10)/K) use the delta method.
   Kimura (1980) likelihood-ratio tests compare the general model with
   models having a common Linf, common K, common t0, or all common:
   -2 ln(Lambda) = N ln(RSS_H / RSS_general), chi-square with
   (groups - 1) df per constrained parameter.
2. `nimble` with `group_effects: fixed`: the same general model, Bayesian;
   groups compared by posterior differences.
3. `nimble` with `group_effects: random` (hierarchical; Helser & Lai 2004):

       log Linf_g ~ Normal(mu_logLinf + Linf_beta * z_g, sd_log_Linf)
       log K_g    ~ Normal(mu_logK    + K_beta    * z_g, sd_log_K)
       t0 common to all groups

   z_g is the group covariate standardised by its mean and SD across groups
   (optional). `Linf_typical = exp(mu_logLinf)` is the Linf of a typical
   group at the average covariate value; `Linf_beta` is the change in log
   Linf per 1 SD of the covariate, so `exp(beta) - 1` is the proportional
   change.

Priors (Bayesian engines): log Linf (or its mean) ~ N(log max observed
length, 1); log K ~ N(log 0.3, 1.5); t0 ~ N(0, 2) truncated below the
youngest age (keeps predicted lengths positive); among-group SDs ~
half-normal(0, 0.5); betas ~ N(0, 1); sigma ~ U(0, 1) on the log scale
(or U(0, max length) for additive error). Samplers: automated-factor slice
on each group's (log Linf, log K) [and t0], plus block samplers that move
t0 jointly with all K_g.

**Length-weight and condition (all engines, least squares):**

    log W_i = log a + b log L_i + log(condition_group) + e_i,  e ~ N(0, s)

common slope with group intercepts (ANCOVA). The interaction model
(group-specific b) is compared with an F test and reported. Predicted mean
weights are back-transformed with the bias correction exp(s^2/2) (Sprugel
1983). Condition is reported as

- `condition_ratio_<g>_vs_<ref>`: weight of a group-g fish relative to a
  reference-group fish of the same length (1.05 = 5% heavier);
- relative weight Wr = 100 W / Ws (Wege & Anderson 1978) when a
  standard-weight equation is given (`standard_weight:` in config; the FSA
  `WSlit` database, e.g. walleye from Murphy et al. 1990), or relative
  condition Kn = W / predicted W (Le Cren 1951) otherwise.

## Estimands (names in the estimates table)

| name | meaning | scale |
|---|---|---|
| `Linf_<g>`, `K_<g>`, `t0_<g>` | VB parameters for group g (t0 only for fixed groups) | length units; per year; years |
| `age_at_90pct_Linf_<g>` | age at which the median fish reaches 90% of Linf | years |
| `Linf_diff_<g>_vs_<ref>`, `K_diff_...`, `t0_diff_...` | difference from the reference group (null 0) | as above |
| `Linf_typical`, `K_typical`, `t0` | hierarchical: typical group at average covariate | |
| `sd_log_Linf`, `sd_log_K` | hierarchical: among-group SD on the log scale (0.1 = about 10% variation); point estimate is the posterior median (the posterior is right-skewed with few groups) | |
| `Linf_beta_<cov>`, `K_beta_<cov>` | hierarchical: effect on log Linf / log K per 1 SD of the covariate (null 0) | link (log) |
| `residual_sd` | sigma (log scale for multiplicative, length units for additive) | |
| `lw_b`, `lw_a`, `lw_residual_sd` | W = a L^b for the reference group, residual SD on log scale | |
| `mean_weight_at_<L>_<ref>` | bias-corrected mean weight of a reference-group fish of length L | weight units |
| `condition_ratio_<g>_vs_<ref>` | relative weight at the same length (null 1) | ratio |

## Assumptions (and what breaks if they fail)

- **Ages are correct.** Ageing error is not modelled. Under-ageing old fish
  flattens the curve (Linf too low, K too high).
- **Each age class is sampled at random with respect to length.** Gear that
  misses small fish (gill-net mesh, minimum-size regulation) catches only
  the fastest-growing young fish. Our negative control shows this biases K
  down by about 14% and t0 by about -0.36 years, with interval coverage as
  low as 0.58 (VALIDATION.md, section 4). If your gear is size-selective, drop the youngest
  partially-selected ages or model the selectivity outside this module.
- **Enough old fish.** Linf is identified only if fish approach the
  asymptote. With few fish beyond the age at 90% Linf, Linf (and K, which
  is negatively correlated with it) are extrapolations. The module prints a
  note when fewer than 5 fish are beyond that age or the oldest fish is
  younger than it.
- **One curve per group.** Fish within a group share Linf, K, t0
  (individual variation goes into sigma). Growth of cohorts that
  experienced different conditions is averaged.
- **Independent fish.** Fish from the same net set are treated as
  independent.
- **Hierarchical model:** group effects are exchangeable (normal), t0 is
  shared, and Linf and K effects are independent (no Linf-K correlation).
  With fewer than 10 groups the among-group SD of log Linf is biased high
  (posterior median +15% at 8 groups x 70 fish, +24.5% at 6 groups x 40 fish;
  both fail the recovery gate): report it through its 95% interval, not as a
  point value. Group-level Linf, K, the typical values and the covariate
  effects are still recovered at 8 groups.
- **Covariate effects are associations among groups**, not causal effects.

## Data format

One row per fish (CSV). Column names are set in the config.

| column | required | notes |
|---|---|---|
| fish id | no | unique per row; used in messages |
| age | yes | years (integer otolith age is fine; add `age_offset` for capture date) |
| length | yes | total or fork length, one unit throughout |
| group | no | population, lake, sex, year...; omit for a single group |
| weight | no | leave blank where not weighed |

```
fish_id,year,location,sex,age,tl_mm,wt_g
WAE2014-00001,2014,Toledo-Huron,female,2,412,679
WAE2014-00002,2014,Toledo-Huron,male,3,455,
```

Optional group-covariate CSV (hierarchical engine only): one row per group,
a column matching the group values and a numeric covariate column.

## Config

```yaml
data:
  file: walleye_lake_erie.csv       # relative to the config file
  fish_id: fish_id
  age: age
  length: tl_mm
  weight: wt_g                      # optional
  group: sex                        # optional
  group_levels: [female, male]      # optional order; first = reference unless reference_group given
  reference_group: female           # differences and condition ratios are vs this group
  filter: {year: 2014}              # optional row filter: column -> allowed values
  age_offset: 0.5                   # optional, added to every age (fraction of year at capture)
  length_units: mm
  weight_units: g
  reference_length: 500             # length for mean_weight_at_<L> (default: median length)
  standard_weight: {species: Walleye, group: overall, units: metric}   # FSA::wsVal; or
  # standard_weight: {log10_a: -5.453, b: 3.18, min_length: 150, source: "Murphy et al. 1990"}
  group_covariates: {file: lake_covariates.csv, group: lake, covariate: temperature}  # hierarchical only
settings:
  engine: nls                       # nls | nimble
  error: multiplicative             # multiplicative | additive
  group_effects: random             # nimble only: random (>= 3 groups, 8+ recommended) | fixed
  kimura: true                      # nls only: run the 4 likelihood-ratio tests
  niter: 20000                      # nimble only (hierarchical models need >= 20000)
  nburnin: 5000
  nchains: 3
figures:
  labels: {female: "Female", male: "Male", calendar_year: "Year"}
```

## Outputs

- Estimates table (above), diagnostics (R-hat, bulk ESS for Bayesian
  fits), Kimura LRT table (`fit$object$growth$kimura`), length-weight slope
  test and Wr/Kn by group (`fit$object$lw`).
- Interpretation paragraphs: Linf, K and age at 90% Linf per group with
  95% intervals; group differences; Kimura tests; hierarchical variation and
  covariate effects as % change per 1 SD; b with an isometry statement;
  condition differences; caveats and data notes.
- Figures: see the next section.

## Figures

`core/run.R fit` exports these to `<out>/figures/` (PDF + PNG and a
`gallery.html`); `core/run.R figures fisheries-growth --fit <out>/fisheries-growth_fit.rds --config cfg.yaml`
redraws them from the saved fit after you edit the `figures:` block.
Each group keeps the same colour, symbol and line type in every figure
(Okabe-Ito order; shape/linetype keep groups apart in greyscale). With
more than 8 groups the figures switch to facets or single-colour layouts.

| id | what it shows |
|---|---|
| `growth_curves` | Length at age: individual fish (ages jittered), fitted von Bertalanffy curve per group drawn over that group's observed ages, `ci` confidence (nls) or credible band (nimble, plus the `inner_ci` band) for the curve, and the `ci` prediction interval for single fish (dotted). |
| `parameter_comparison` | Linf, K, t0 (fixed-effects models) and age at 90% Linf per group, with units in the panel titles: estimate with `ci` Wald interval (nls) or posterior median with `ci` and `inner_ci` credible intervals (nimble). |
| `residuals` | A: residuals against fitted length with a loess smooth (curve-shape check). B: residuals by age as box plots (error-model check). Log-scale residuals for the multiplicative error model. |
| `length_weight` | Log-log length-weight plot: weighed fish and the bias-corrected mean weight at length per group with a `ci` confidence band (common slope b). Only when weights are present. |
| `condition` | Relative weight Wr (with a standard-weight equation; dashed line at 100) or relative condition Kn (dashed line at 1) by group: individual fish and the group mean with `ci` confidence interval. |
| `Linf_<cov>`, `K_<cov>` | Hierarchical model with a group covariate only: Linf or K of a typical group against the covariate in its original units (posterior median, `ci` and `inner_ci` bands), with each group's posterior mean and 95% interval. |
| `trace`, `posterior` | Bayesian fits only: MCMC traces per chain and posterior densities for key parameters. |

Relabel anything with `figures: labels:`. Keys are group names (`female`,
`2003`), axes (`age`, `length`, `weight`, `fitted_length`, `residual`,
`condition`), parameters (`Linf`, `K`, `t0`, `age_at_90pct_Linf`, and in
traces `Linf_<group>`, `K_typical`, `sd_log_Linf` and so on) and the
covariate name (`calendar_year`). `units: {<cov>: {center, scale}}`
overrides the covariate back-transform that the hierarchical model records.

## Simulation and power

`simulate()` draws ages from a catch-at-age distribution
P(age = a) ∝ exp(-Z (a - age_min)), a = age_min..max_age, lengths from VB
with the chosen error, and weights from W = a_g L^b exp(e).

Design knobs: `n_groups`, `n_per_group` (scalar or per group), `Z` (total
mortality: higher = fewer old fish), `age_min`, `max_age`, `age_offset`,
`prop_weighed`, `min_length` (length-selective gear; for negative controls),
`reference_length`, `hierarchical`, `cov_mean`, `cov_sd`, `cov_name`.

Parameters: fixed groups `Linf`, `K`, `t0` (vectors, one per group),
`sigma`, `error`; hierarchical `Linf_typical`, `K_typical`, `t0`,
`sd_log_Linf`, `sd_log_K`, `Linf_beta`, `K_beta`; length-weight `lw_a`,
`lw_b`, `lw_sigma`, `condition` (per group, relative to group 1).

Simulated groups are named A, B, C...; power parameters are
`Linf_diff_B_vs_A`, `K_diff_B_vs_A`, `t0_diff_*`, `condition_ratio_*`,
`Linf_beta_temperature`, `K_beta_temperature`. Example:
`core/run.R power fisheries-growth --vary n_per_group=25,50,100,200 --nsim 100`.
At the default design (Linf 600 vs 650 mm, K 0.30 vs 0.25, Z = 0.35) the
single-parameter Wald tests are weak because Linf and K are strongly
correlated: power for the Linf difference is 0.13, 0.30, 0.43 and 0.59 with
60, 100, 150 and 250 aged fish per group (K difference 0.12-0.48), whereas
Kimura's joint test that the whole curve differs rejects in 99.6% of data sets
at 150 fish per group, and a 5% condition difference is detected with power
0.94 at 60 fish per group (VALIDATION.md section 6). Use the whole-curve test
(or age at 90% Linf) to ask "do these groups grow differently?" and expect
to need several hundred fish per group to attribute the difference to Linf
or K specifically.

## Pitfalls

- Do not compare groups whose age ranges differ a lot: the curves then
  differ in where they are extrapolated, not in growth.
- Linf is not "maximum length". Some fish are longer than Linf (sigma) and
  with few old fish Linf can be far above the largest fish.
- The parameters are strongly correlated. A difference in Linf often shows
  up as a compensating difference in K. Compare whole curves (figure) and
  ages at 90% Linf, not just one parameter.
- Kimura tests with very large samples detect tiny differences; read the
  size of the difference (estimates table), not only the p-value.
- Relative weight depends on the standard-weight equation and its minimum
  length (fish below it get no Wr). Use the equation for your species and
  length type (TL vs FL).
- The example data (Lake Erie walleye, FSAdata::WalleyeErie2, Ohio DNR) are
  for educational use only. Do not use them for other purposes without
  permission from the source.

## References

- von Bertalanffy L (1938) A quantitative theory of organic growth. Human Biology 10:181-213.
- Beverton RJH, Holt SJ (1957) On the Dynamics of Exploited Fish Populations. Fishery Investigations Series II 19. HMSO, London.
- Kimura DK (1980) Likelihood methods for the von Bertalanffy growth curve. Fishery Bulletin 77:765-776.
- Helser TE, Lai HL (2004) A Bayesian hierarchical meta-analysis of fish growth: with an example for North American largemouth bass, Micropterus salmoides. Ecological Modelling 178:399-416.
- Ogle DH (2016) Introductory Fisheries Analyses with R. Chapman & Hall/CRC.
- Ogle DH, Doll JC, Wheeler AP, Dinno A. FSA: Simple Fisheries Stock Assessment Methods. R package 0.10.1.
- Nelson GA. fishmethods: Fishery Science Methods and Models. R package (growthlrt, used for the cross-check).
- Froese R (2006) Cube law, condition factor and weight-length relationships: history, meta-analysis and recommendations. Journal of Applied Ichthyology 22:241-253.
- Sprugel DG (1983) Correcting for bias in log-transformed allometric equations. Ecology 64:209-210.
- Le Cren ED (1951) The length-weight relationship and seasonal cycle in gonad weight and condition in the perch (Perca fluviatilis). Journal of Animal Ecology 20:201-219.
- Wege GJ, Anderson RO (1978) Relative weight (Wr): a new index of condition for largemouth bass. Pages 79-91 in New Approaches to the Management of Small Impoundments. AFS North Central Division Special Publication 5.
- Murphy BR, Brown ML, Springer TA (1990) Evaluation of the relative weight (Wr) index, with new applications to walleye. North American Journal of Fisheries Management 10:85-97.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.
- Pinheiro JC, Bates DM (2000) Mixed-Effects Models in S and S-PLUS. Springer (nlme, used for the cross-check).
