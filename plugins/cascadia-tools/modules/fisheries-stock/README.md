---
id: fisheries-stock
title: Fish stock assessment - catch curves, stock-recruitment and surplus production
when_to_use: Fisheries questions about mortality, recruitment or stock status. (1) catch_curve - aged fish from a survey or fishery sample (one year or several years); question is total mortality Z / annual survival of fully recruited fish. (2) stock_recruit - a brood table of spawners and the recruits they produced; question is productivity, density dependence, S_MSY, U_MSY, steepness. (3) surplus_production - an annual catch series plus an abundance index (CPUE or survey); question is MSY, B/B_MSY, F/F_MSY (Kobe status).
not_for: Length-only data (convert to age with an age-length key first, or use a length-based method), separating fishing from natural mortality (needs an M estimate or tagging), growth curves (use fisheries-growth), age-structured assessment (stock synthesis, SAM), catch-only methods without an index, stocks whose index is not proportional to biomass (hyperstability), recruitment series with errors-in-spawners problems, one-way-trip catch/CPUE histories for production models (weakly identified, see Known limits).
data_needed: catch_curve - one row per aged fish (age) or a catch-at-age table (age, count); optionally a year column for age samples from several years (multi-year cohort catch curve). stock_recruit - brood year, spawners, recruits. surplus_production - year, catch (every year), abundance index (may have gaps).
engine: catch_curve - nimble Bayesian Poisson-lognormal catch curve (default) + glmmTMB Poisson-lognormal GLMM (single year) or cohort GLMM (several years) + FSA::catchCurve (weighted regression) + Chapman-Robson (hand code, checked against FSA); stock_recruit - least squares (lm / optim) or nimble (Bayesian, AR(1) residuals); surplus_production - nimble (Bayesian state-space, Meyer & Millar 1999)
---

# Fish stock assessment: catch curves, stock-recruitment and surplus production

Choose the analysis with `settings: model:`.

| `model` | Question | Main outputs |
|---|---|---|
| `catch_curve` | How fast are adult fish dying? | Total instantaneous mortality Z, annual survival S = exp(-Z), annual mortality A = 1 - S |
| `stock_recruit` | How many recruits does a spawning stock produce, and what escapement gives the maximum sustainable yield? | ln(alpha) (productivity), beta (density dependence), S_MSY, U_MSY, equilibrium spawners, steepness, recruitment SD, AR(1) |
| `surplus_production` | How large and productive is the stock, and is it overfished? | r, K, q, MSY, B_MSY, F_MSY, B/B_MSY and F/F_MSY in the final year, Kobe plot |

## 1. Catch curves (`catch_curve`)

**Model.** If recruitment and total mortality are constant and every age from the
first age used is fully vulnerable, numbers at age decline geometrically:
N_a = N_0 exp(-Z a), so log(catch) is linear in age with slope -Z.

**Real populations have strong and weak year classes.** Recruitment is never
constant: some spawning years produce two or three times the average number of
young, others a fraction of it, and those year classes stay strong or weak for
their whole life. In a single age sample this shows up as points scattered above
and below the straight catch-curve line (a strong year class at age 5 sits above
the line, a weak one at age 7 below it). This scatter is not sampling noise that
more fish would remove: each year class is one draw from a variable process. An
estimator that treats the ages as independent Poisson or normal observations
with only sampling error gives intervals that are too narrow. The estimators
below differ mainly in how they handle this.

| estimator (`primary:`) | how it treats year-class variability | 95% interval covered the true Z (single-year designs, VALIDATION.md section 5) | use when |
|---|---|---|---|
| `bayes` (**default**) | Year-class random effect with its SD estimated *and its uncertainty carried into the interval* (Bayesian, Nimble) | 0.95-0.98 (incl. two independent simulators); bias +1% to +10%, always upward | Default for one age sample of randomly sampled fish. Slower (Nimble compiles for ~1 min). With fish from a few net sets/hauls (cluster sampling) its Z was biased +8% to +15% (see Pitfalls). |
| `glmm` | Year-class random effect (Millar 2015), REML fit and a t interval (`glmm_interval: reml_t`, default) | 0.89-0.96 | Fast cross-check. Passed the builder's three designs and all six designs of the second independent verification (0.927-0.957, incl. cluster sampling) but covered only 0.89 in the first independent simulator (Z 0.8, 250 fish). Least biased under cluster sampling (+2% to +8%). |
| `glmm` with `glmm_interval: wald` | Same model, ML fit, Wald interval (the module's original GLMM) | 0.86-0.90 (0.80-0.86 without the zero tail) | Not recommended; kept for comparison. |
| `chapman_robson` | Inflates the SE by sqrt(c-hat), a single overdispersion factor | 0.77-0.92 | Constant recruitment; comparison with published work. |
| `weighted_regression` | Residual variance from the regression; weights favour young ages | 0.84-0.93 | Comparison with published work (Smith et al. 2012). |
| `cohort_glmm` (multi-year data, default there) | Each year class is seen in several years, so its strength is estimated directly | 0.93-0.97 (incl. a missing survey year and cluster sampling) | Age samples from 2 or more years. |
| `cohort_bayes` (multi-year) | As `cohort_glmm`, Bayesian | 0.93-0.97 | Cross-check of `cohort_glmm`. |

Estimators:

- **Bayesian Poisson-lognormal catch curve** (`primary: bayes`, default for one
  sample). catch_a ~ Poisson(mu_a), log mu_a = b0 - Z (a - mean age) + sigma u_a,
  u_a ~ Normal(0, 1) (one year class per age). **Priors:** Z ~ Uniform(0, 5);
  sigma (year-class SD, log scale) ~ half-Normal(0, 1) (`prior_yc_sd`; typical
  recruitment SDs of fish are 0.3-1); b0 ~ Normal(0, 10). Fitted in Nimble, 3 chains
  x 40,000 iterations (10,000 burn-in, thin 10; joint automated-factor slice sampler on
  Z, sigma and the intercepts). Estimate = posterior median, interval = 95%
  equal-tailed credible interval. Ages from the Chapman-Robson first age onwards
  plus `zero_tail` (default 5) ages of zero catch beyond the oldest fish (see below).
  It reports R-hat and bulk ESS; the fit flags R-hat >= 1.05.
- **Poisson-lognormal GLMM** (Millar 2015, `glmm`). The same model by likelihood
  with glmmTMB. `glmm_interval:` chooses the interval: `reml_t` (default: REML
  fit, Z +- t SE on (number of year classes - 2) df), `wald` (ML, Z +- 1.96 SE),
  `profile` (profile likelihood), `boot_pct` / `boot_t` (parametric bootstrap,
  percentile / studentized, `n_boot`, default 999; 299 in recovery runs).
  The ML estimate of the year-class SD is biased low with ~10-15 ages (mean 0.18
  when the truth is 0.30), which is why Wald, profile and bootstrap intervals
  (which all condition on that estimate) stay too narrow; REML with a t quantile
  and the Bayesian posterior both allow for it.
- **Zero tail** (`zero_tail: 5`). The oldest fish in the sample is, by
  definition, a non-empty age. Fitting a Poisson model only up to that age
  conditions on a positive last count and flattens the curve (the old GLMM was
  -4% to -12% biased because of this). The GLMM and Bayesian fits therefore add
  zero catches for 5 ages beyond the oldest fish: the sample did not contain
  older fish, and that is information. Set `zero_tail: 0` to switch it off.
- **Weighted linear regression** (`weighted_regression`). log(catch) on age from
  the modal (peak) age onwards, ages with zero catch excluded, weights equal to the
  predicted log catch of an unweighted fit (Maceina & Bettoli 1998), using
  `FSA::catchCurve(weighted = TRUE)`. 95% t interval. Nearly unbiased but the
  interval is too narrow under year-class variability (treat it as ~90%).
- **Chapman-Robson** (Chapman & Robson 1960). Ages recoded x = 0, 1, ... from the
  first age used; S = T / (n + T - 1), T = sum(x * catch), n = number of fish;
  Z = -log(S) - (n-1)(n-2) / (n (T+1)(n+T-1)) (Hoenig et al. 1983 bias
  correction). SE of Z = (1 - e^-Z) / sqrt(n e^-Z) (Smith et al. 2012 eq. 2, as
  in FSA) multiplied by sqrt(c-hat), c-hat = Pearson chi-square of the geometric
  age distribution / df (floored at 1); t interval on those df. First age = peak
  age + 1 ("peak plus", Smith et al. 2012).

Every fit reports all single-year estimators side by side (and the Bayesian one
only when it is the primary or `bayes: true`); `primary:` chooses which feeds
`total_mortality_Z`, `annual_survival` and `annual_mortality`.

**Multi-year (cohort) catch curve.** When the data hold age samples from two or
more years (`data: year_col`), the module fits

    log E[catch(a, t)] = alpha_t + c_(t-a) - Z a

with alpha_t a fixed effect for each sampling year (it absorbs the number of fish
aged and the sampling effort that year), c_(t-a) ~ Normal(0, s^2) a random effect
for the year class born in year t-a, **shared by that year class in every year it
was sampled**, and a common Z. Because a strong year class is seen again one age
older the next year, its strength is estimated rather than mistaken for a change
in mortality, so this model gives honest intervals even when year classes vary a
lot (coverage 0.95-0.97 with 5 or 10 years, VALIDATION.md section 5). It is
fitted with glmmTMB (`cohort_glmm`, default, REML + t interval) and, with
`bayes: true` or `primary: cohort_bayes`, in Nimble with the priors above
(year intercepts ~ Normal(0, 10)). Chapman-Robson and the weighted regression are
also reported on the age composition **pooled** over years; pooling does not
remove year-class effects, and their intervals covered only 0.58-0.85 in the
multi-year simulations. First age = peak of the pooled age composition + 1;
unobserved year-age cells are zeros.
- `z_trend: true` adds a linear change in Z over the sampling years: Z_y = Z +
  dZ (y - mid year), Z before the first sample = Z in the first year, and each
  year class accumulates the Z of the years it lived through (an exposure
  covariate, not a year-by-age interaction). It reports `Z_trend` (change in Z
  per year) and Z in the middle year. **Weakly identified:** with 10 years of 250
  fish, Z_trend intervals covered the truth (0.94-0.95) and the false-positive
  rate was 0.052, but a true trend of +0.03/yr was estimated at +0.019 (bias
  -36%) and switching the trend on widens the Z interval about tenfold. Use it to
  ask "is there evidence Z has changed?", not to measure the change.

**Assumptions.** No trend in recruitment over the year classes sampled (year
classes vary at random around a constant mean), constant Z over the ages used,
full and equal selectivity from the first age used, representative sampling and
accurate ages. The random year-class effects are assumed independent and
lognormal.
- **Recruitment trend.** If recruitment has been increasing, Z is over-estimated;
  if it has been decreasing, Z is under-estimated. No estimator here protects
  against this (negative controls: +25 to +29% bias, VALIDATION.md).
- **Dome-shaped selectivity.** Z is over-estimated.
- **Ageing error.** Old fish are typically under-aged, which flattens the curve.
- **Autocorrelated recruitment** (runs of good and poor year classes) is not modelled.
- **Cluster sampling.** The models treat every aged fish as an independent random
  draw. Fish aged from a few net sets, hauls or sites, each catching its own mix of
  ages, violate this (Nelson 2019). With 20 hauls the single-year Bayesian Z was
  biased +8% (intra-haul correlation 0.1) to +15% (0.2, fails the gate), the GLMM
  +2% to +8% and Chapman-Robson -2% to -3%; the multi-year cohort GLMM still passed
  (+2.5%, coverage 0.93). VALIDATION.md section 6.
- **Counts must be fish actually aged.** Catch-at-age expanded with an age-length
  key or scaled to the total catch (non-integer or inflated counts) is rejected by
  `check_data`; supply the aged-fish counts.

**Data.** Single year: one row per fish with an `age` column (other columns are
ignored, rows with missing age are dropped), or a catch-at-age table with
`count_col`. Several years: add `year_col`; either one row per fish (year, age) or
a long table with one row per year x age and a count column. Missing ages
and year-age combinations are filled with 0. Ages must be whole years (text such as
'four' or '6+' and fractional ages such as 3.5 stop with the rows to fix), and counts
must be whole numbers of aged fish.

```
fish_id,sample_date,basin,gear,total_length_mm,age
WAE25-0001,2025-09-22,North Basin,experimental gill net,412,3
```
```
survey_year,age,n_fish
2019,3,77
2019,4,64
```

## 2. Stock-recruitment (`stock_recruit`)

**Models** (Ricker 1954; Beverton & Holt 1957; Hilborn & Walters 1992), each with
lognormal errors:

    Ricker:         R_t = alpha S_t exp(-beta S_t) exp(e_t)
    Beverton-Holt:  R_t = alpha S_t / (1 + beta S_t) exp(e_t)

- **`sr_engine: ls`** (default). Least squares on the log scale. For the Ricker
  model this is a linear regression of log(R/S) on S; for Beverton-Holt it is
  nonlinear least squares in (ln alpha, log beta). Intervals for derived quantities
  come from 4000 draws of the multivariate-t sampling distribution of the
  parameters. Residuals are assumed independent. If the lag-1 autocorrelation of
  the residuals exceeds 2/sqrt(n), the fit notes warn you.
- **`sr_engine: nimble`** with `ar1: true`. A Bayesian model with AR(1) residuals:
  e_t = phi e_{t-1} + nu_t, nu_t ~ N(0, sigma_nu^2), stationary start, marginal SD
  sigma = sigma_nu / sqrt(1 - phi^2). Priors are ln(alpha) ~ N(0, 5^2),
  beta ~ U(0, 20/median(S)) (set with `beta_max`), sigma_nu ~ U(0, 3) and
  phi ~ U(-0.99, 0.99). Sampling uses joint automated-factor slice samplers.
  **Use this when recruitment comes in runs of good and poor years.** Least squares
  then gives intervals that are far too narrow (see the negative control).

**Reference points.** These assume recruits are measured in spawner units, i.e.
adult returns, as in salmon brood tables. The replacement line is R = S.
- **Ricker:** S_MSY = (1 - W(e^(1 - ln alpha))) / beta and U_MSY = 1 - W(e^(1 - ln alpha)),
  where W is the Lambert W function (exact solution, Scheuerell 2016; checked
  numerically). The Hilborn (1985) approximation S_MSY ~ (ln alpha / beta)(0.5 - 0.07 ln alpha)
  is not used. The equilibrium is S_eq = ln(alpha) / beta, and the spawner
  abundance giving maximum recruitment is 1/beta.
- **Beverton-Holt:** S_MSY = (sqrt(alpha) - 1) / beta, U_MSY = 1 - 1/sqrt(alpha),
  S_eq = (alpha - 1) / beta and R_max = alpha / beta.
- **Steepness** (Mace & Doonan 1988) is h = R(0.2 S0) / R0 = 0.2 (alpha phi0)^0.8
  (Ricker; this can exceed 1) or alpha phi0 / (4 + alpha phi0) (Beverton-Holt).
  phi0 = unfished spawners per recruit (`spr0`, default 1 under the returns
  convention). If your recruits are age-0 or age-1 fish, give `spr0` from a
  spawner-per-recruit calculation. S_MSY / U_MSY are then not yield-based MSY
  quantities; use steepness only.
- All reference points are on the **median** recruitment scale. No lognormal bias
  correction is applied.

**Assumptions.**
- Stationarity: no regime shift.
- Spawners measured with small error. Measurement error in S biases beta and
  alpha (errors-in-variables).
- Enough contrast in S. The fit warns if max(S) < 0.5/beta.
- Simulations here treat spawners as exogenous. In a real closed loop (S_t depends
  on earlier recruitment), short series overestimate productivity (time-series
  bias, Walters 1985). That case is not covered by the validation.

**Data.**

```
brood_year,escapement_thousands,total_recruits_thousands
1985,213.4,61.6
```

## 3. Surplus production (`surplus_production`)

**Model.** A Bayesian state-space Schaefer model (Schaefer 1954; Meyer & Millar
1999) in Nimble:

    P_t = B_t / K
    P_1 ~ lognormal(log p1, s_p1)
    P_t ~ lognormal(log(max(P_{t-1} + r P_{t-1}(1 - P_{t-1}) - C_{t-1}/K, 0.001)), sigma_proc)
    I_t ~ lognormal(log(q K P_t), sigma_obs)
    MSY = rK/4,  B_MSY = K/2,  F_MSY = r/2 (exploitation rate),
    B/B_MSY = 2 P_T,  F/F_MSY = (C_T / B_T) / (r/2)

**Priors** (defaults; set your own in `settings`):
- r ~ lognormal(median 0.3, log-SD 0.6) (`prior_r: [median, sdlog]`). Use a
  resilience-based prior, e.g. from FishBase.
- K ~ lognormal(median 10 x mean catch, log-SD 1) (`prior_K: [median, sdlog]`).
- log q ~ uniform over 30 log units. This is equivalent to the vague 1/q prior of
  Meyer & Millar.
- 1/sigma_proc^2 ~ Gamma(3.785518, 0.010223) and 1/sigma_obs^2 ~
  Gamma(1.708603, 0.008613854). These are the Meyer & Millar (1999) priors, taken
  from their published WinBUGS code (`prior_sigma_proc`, `prior_sigma_obs`). They
  put process-error SD near 0.05-0.1 and observation-error SD near 0.1-0.2.
- P_1 ~ lognormal(median 1, log-SD 0.2) (`prior_P1`). Set the median below 1 if
  the series starts after substantial fishing.

**Sampler.** The default run is 3 chains x 100,000 iterations (30,000 burn-in, thin
20). Joint AF_slice and RW_block samplers act on (log r, log K, log q). Biomass
states use random-walk samplers.

**Assumptions.**
- Catch is known without error.
- The index is proportional to biomass (constant q). Hyperstable CPUE biases
  status optimistic.
- Production follows a symmetric Schaefer curve (B_MSY = K/2).
- There is enough contrast. Hilborn (1979) and Polacheck et al. (1993) showed that
  a "one-way trip" (catch rising while the index falls, with no recovery) cannot
  separate a small productive stock from a large unproductive one. The fit reports
  the posterior/prior SD ratio for r and K and warns when the data add little
  information (ratio > 0.7).

**Data.** One row per year; catch is required in every year, and the index may be
missing in some years.

```
year,landings_t,cpue_kg_per_hook
1990,410,
1992,745,1.808
```

## Config (every key)

```yaml
settings:
  model: catch_curve | stock_recruit | surplus_production
  # catch_curve
  primary: bayes                  # one year: bayes (default) | glmm | chapman_robson | weighted_regression
                                  # several years: cohort_glmm (default) | cohort_bayes | chapman_robson | weighted_regression
  glmm_interval: reml_t           # GLMM interval: reml_t (default) | wald | profile | boot_pct | boot_t
  n_boot: 999                     # bootstrap refits (boot_pct / boot_t)
  bayes: false                    # true: also fit the Bayesian model when it is not the primary
  prior_yc_sd: 1                  # half-Normal scale of the year-class SD prior (Bayes)
  zero_tail: 5                    # zero-catch ages added beyond the oldest fish (GLMM, Bayes)
  z_trend: false                  # several years: linear change in Z over the sampling years
  age_full: peak_plus             # first age for CR/GLMM/Bayes: peak_plus | peak | <age>
  age_full_regression: peak       # first age for the regression
  age_last: null                  # oldest age used (default: oldest in data); below the oldest fish, older ages are dropped and no zero tail is added
  overdispersion: true            # c-hat inflation of the Chapman-Robson SE
  # MCMC (bayes): niter (40000), nburnin (10000), thin (10), nchains (3), seed
  # stock_recruit
  sr_form: ricker | beverton_holt
  sr_engine: ls | nimble
  ar1: true                       # nimble only
  spr0: 1                         # unfished spawners per recruit (steepness)
  beta_max: null                  # upper bound of the beta prior (nimble)
  n_draws: 4000                   # LS: draws for derived-quantity intervals
  # MCMC (nimble): niter, nburnin, thin, nchains (3), seed
  # surplus_production priors
  prior_r: [0.3, 0.6]
  prior_K: [<10 x mean catch>, 1]
  prior_sigma_proc: [3.785518, 0.010223]
  prior_sigma_obs: [1.708603, 0.008613854]
  prior_P1: [1, 0.2]
data:
  file: <csv>
  # catch_curve: age_col (age), count_col (optional), year_col (optional: several sampling years), species, gear
  # stock_recruit: year_col (brood_year), spawners_col, recruits_col, species, spawner_units, recruit_units
  # surplus_production: year_col (year), catch_col (catch), index_col (cpue), species, catch_units, index_name
figures: {...}                    # see skills/eco-figures
```

Examples are in `example/`:
- `config.yaml` (walleye catch curve, one year; Bayesian default)
- `config_multiyear.yaml` (smallmouth bass, 6 years of ages; cohort catch curve) and `recovery_multiyear.yaml`
- `config_stock_recruit.yaml` (sockeye brood table, Bayesian Ricker with AR(1))
- `config_surplus_production.yaml` (longline catch and CPUE)

## Outputs

- **Estimates table** (contract form, with rhat / ess_bulk for Bayesian fits).
  For the Nimble engines the `estimate` column is the **posterior median** (with
  95% equal-tailed credible interval), as in Meyer & Millar (1999) and as plotted
  in the figures; posterior means of K, steepness and the SDs are pulled up by
  their right tails (independent verification, VALIDATION.md).
- **Diagnostics notes:**
  - catch_curve: ages used, c-hat, year-class SD, interval method, estimator disagreement; R-hat and bulk ESS for the Bayesian fit.
  - stock_recruit: residual autocorrelation, contrast.
  - surplus_production: prior-vs-posterior information, R-hat.
- **Interpretation** paragraphs on the biological scale.
- **Figures** (`docs/FIGURES.md`): see the Figures section below.

## Figures

`core/run.R fit` exports the suite for the chosen `model` to `<out>/figures/`
(PDF + PNG, `gallery.html`, `figures_index.csv` with captions). Redraw without
refitting after editing the `figures:` block:
`scripts/Rscript core/run.R figures fisheries-stock --fit <out>/fisheries-stock_fit.rds --config <config.yaml>`.
Each analysis can only draw the figures its data support (a catch-at-age table has no
stock-recruit pairs or biomass series), so the six required ids are split across the three models.

| Model | Figure id | What it shows |
|---|---|---|
| catch_curve | `catch_curve` | Number caught at each age on a log axis (fish, not log units; pooled over years for multi-year data); filled circles = ages used by Chapman-Robson and the regression, open squares = regression only, crosses = not fully recruited (excluded). Chapman-Robson expected catch (solid), weighted regression (dashed) with its confidence band, and the primary GLMM / Bayesian line for an average year class. |
| catch_curve | `mortality_estimates` | Z (per year) with 95% intervals for every estimator fitted (Chapman-Robson, weighted regression, GLMM, Bayesian; cohort GLMM / Bayes for several years); dashed line = primary estimate. Estimator colours match `catch_curve`. |
| catch_curve (several years) | `catch_by_year` | Catch at age in each sampling year (one panel per year, log scale) with the cohort GLMM's expected catch with and without each year class's strength. |
| catch_curve (several years) | `yearclass_strength` | Relative strength of each year class (exp(random effect), 1 = average) with intervals. |
| catch_curve (Bayes) | `trace` | MCMC traces of Z and the year-class SD (and dZ with `z_trend`). |
| stock_recruit | `stock_recruitment` | Fitted Ricker / Beverton-Holt curve with its band (95% and 50% credible bands for nimble; 95% confidence band for least squares), brood-year points labelled by year (ggrepel), replacement line R = S (when spr0 = 1) and S_MSY. |
| stock_recruit | `recruitment_residuals` | log(R / expected R) by brood year with intervals from parameter uncertainty; dashed line = no deviation; lag-1 autocorrelation in the caption. |
| stock_recruit (nimble) | `trace`, `posterior` | MCMC traces per chain and posterior densities of ln(alpha), beta, sigma (and phi, sigma_nu with AR(1)). |
| surplus_production | `biomass_trajectory` | Biomass (posterior median, 95% and 50% credible bands), index / median q (open points), catch (grey bars, same tonnes axis), B_MSY line. |
| surplus_production | `kobe` | B/B_MSY vs F/F_MSY: median trajectory with first and last years labelled, final-year posterior cloud, conventional quadrant colours (greys with `greyscale: true`). |
| surplus_production | `surplus_production` | Schaefer production curve rB(1 - B/K) with 95% and 50% bands and the annual surplus production implied by the median biomass (points shaded by year). |
| surplus_production | `index_fit` | Observed index (points) and predicted qB_t with 95% and 50% bands. |
| surplus_production | `trace`, `posterior` | Traces of r, K, q; posterior vs prior densities of log10 r, K and MSY. |

Customising (all through the `figures:` block, no code):
- `ci` / `inner_ci` set the band levels.
- `greyscale: true` switches to greys. Estimators keep distinct linetypes and the age classes keep distinct shapes.
- Estimator colours come from `palette` in a fixed order (Chapman-Robson, weighted regression, GLMM, Bayesian, cohort GLMM, cohort Bayesian), so they stay the same across figures.
- `labels:` keys:
  - axes: `age`, `log_catch`, `Z`, `spawners`, `recruits`, `year`, `resid`, `biomass`, `biomass_catch`, `surplus`, `index`, `B_over_BMSY`, `F_over_FMSY`
  - estimators and age classes: `chapman_robson`, `weighted_regression`, `glmm`, `bayes`, `cohort_glmm`, `cohort_bayes`, `used_cr`, `used_rg`, `not_used`, `estimator_title`, `ages_title`
  - multi-year: `fit_with_yc`, `fit_average_yc`, `year_class`, `yc_strength`
  - annotations: `replacement`, `S_MSY_lab`, `B_MSY_lab` (plotmath), `kobe_red`, `kobe_green`, `kobe_yellow`, `kobe_orange`
  - parameters: `ln_alpha`, `beta`, `phi`, `sigma`, `sigma_nu`, `r`, `K`, `q`, `log10_r`, `log10_K`, `log10_MSY`
- `per_figure: {<id>: {width, aspect, height}}` overrides the default sizes. The module sets a height for each figure from its aspect, so the user's values take precedence.
- Units come from `data:` (`spawner_units`, `recruit_units`, `catch_units`, `index_name`). The module does not standardise covariates, so `data$scaling` is not needed.

## Simulation and power

`design: model:` picks the simulator.

**catch_curve**
- Params: `Z`, `recruit_sd_cc` (lognormal year-class SD), `recruit_ar1_cc`,
  `recruit_trend_cc` (log change per cohort), `sel_age50`, `sel_slope`,
  `Z_age_slope`, `sel_dome_age` / `sel_dome_rate`, `Z_trend` (multi-year only:
  change in Z per year during the sampling years).
- Design: `n_fish` (aged fish per sample), `max_age`, `cc_years` (number of annual
  age samples; > 1 simulates a multi-year survey, with the same cohorts followed
  through the years; called `cc_years` because `n_years` is the stock-recruitment
  series length), `first_year`.

**stock_recruit**
- Params: `sr_ln_alpha`, `sr_beta`, `sr_sigma` (marginal SD), `sr_phi`.
- Design: `sr_form`, `n_years`, `S_mean`, `S_cv`, `S_phi` (spawners: lognormal
  AR(1), exogenous), `spr0`, `truth_ar1`.

**surplus_production**
- Params: `r`, `K`, `q`, `sigma_proc`, `sigma_obs`, `P1`.
- Design: `sp_years`, `catch_pattern` (`contrast`: exploitation rises to `F_peak`
  x F_MSY then falls to 0.5 F_MSY, so the stock declines and recovers; `one_way`:
  exploitation rises steadily to `F_end` x F_MSY), `F_cv`, `index_missing`
  (fraction of years without an index).

**Power.** `power_params` are `density_dependence_beta` (null 0) and
`B_over_BMSY_final` (null 1). Use, for example:

```
core/run.R power fisheries-stock --config <cfg> --vary n_years=10,20,30
```

## Pitfalls

- A catch curve gives **Z**, not F. F = Z - M needs an independent M.
- Strong and weak year classes make single-sample catch curves noisy. If you have
  (or can get) age samples from several years, use them together: the
  multi-year cohort catch curve follows each year class and gives honest, much
  narrower intervals (mean width 0.14 with 5 years of 250 fish vs 0.25 from one
  sample of 400 fish, VALIDATION.md section 5).
- Fish aged from a few net sets or hauls (cluster sampling) push the single-year
  Bayesian Z up (+8% to +15% in verification). With such data compare the Bayesian,
  GLMM and Chapman-Robson estimates; prefer several years of samples (cohort GLMM).
- The first age matters. Look at the catch-curve plot: points left of the line
  are not fully recruited.
- Least-squares stock-recruitment intervals are **too narrow** when recruitment is
  autocorrelated. Check the residual plot and the lag-1 note.
- Surplus-production results depend on priors when data lack contrast. Always
  compare posterior and prior (the `posterior` figure) and report the r-K
  correlation.

## References

- Beverton RJH, Holt SJ (1957) On the Dynamics of Exploited Fish Populations. Fishery Investigations Series II, 19.
- Chapman DG, Robson DS (1960) The analysis of a catch curve. Biometrics 16:354-368.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. J Comput Graph Stat 26:403-413.
- Hilborn R (1979) Comparison of fisheries control systems that utilize catch and effort data. J Fish Res Board Can 36:1477-1489.
- Hilborn R, Walters CJ (1992) Quantitative Fisheries Stock Assessment: Choice, Dynamics and Uncertainty. Chapman & Hall.
- Hoenig JM, Lawing WD, Hoenig NA (1983) Using mean age, mean length and median length data to estimate the total mortality rate. ICES CM 1983/D:23.
- Mace PM, Doonan IJ (1988) A generalised bioeconomic simulation model for fish population dynamics. NZ Fish Assess Res Doc 88/4.
- Maceina MJ, Bettoli PW (1998) Variation in largemouth bass recruitment in four mainstream impoundments on the Tennessee River. N Am J Fish Manage 18:998-1003.
- Meyer R, Millar RB (1999) BUGS in Bayesian stock assessments. Can J Fish Aquat Sci 56:1078-1086. (WinBUGS code: stat.auckland.ac.nz/~millar/Bayesian/Surtuna.bugs)
- Brooks ME et al. (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. R Journal 9:378-400.
- Millar RB (2015) A better estimator of mortality rate from age-frequency data. Can J Fish Aquat Sci 72:364-375.
- Nelson GA (2019) Bias in common catch-curve methods applied to age frequency data from fish surveys. ICES J Mar Sci 76(7):2090, doi:10.1093/icesjms/fsz085 (with simple random sampling the random-intercept Poisson GLMM had bias of only -1.6% to +0.8%; with multi-stage cluster sampling (fish aged from a few hauls or net sets, or subsamples expanded to the catch) its bias was mostly positive and up to +139%, the largest of six methods, and the standard errors of every method were too small. See the cluster-sampling pitfall below and VALIDATION.md section 6).
- Ogle DH, Doll JC, Wheeler AP, Dinno A (2025) FSA: Simple Fisheries Stock Assessment Methods. R package 0.10.1.
- Polacheck T, Hilborn R, Punt AE (1993) Fitting surplus production models: comparing methods and measuring uncertainty. Can J Fish Aquat Sci 50:2597-2607.
- Ricker WE (1954) Stock and recruitment. J Fish Res Board Can 11:559-623.
- Schaefer MB (1954) Some aspects of the dynamics of populations important to the management of the commercial marine fisheries. IATTC Bull 1:27-56.
- Scheuerell MD (2016) An explicit solution for calculating optimum spawning stock size from Ricker's stock recruitment model. PeerJ 4:e1623.
- Smith MW, Then AY, Wor C, Ralph G, Pollock KH, Hoenig JM (2012) Recommendations for catch-curve analysis. N Am J Fish Manage 32:956-967.
- Walters CJ (1985) Bias in the estimation of functional relationships from time series data. Can J Fish Aquat Sci 42:147-149.
