---
id: fisheries-stock
title: Fish stock assessment - catch curves, stock-recruitment and surplus production
when_to_use: Fisheries questions about mortality, recruitment or stock status. (1) catch_curve - aged fish from a survey or fishery sample; question is total mortality Z / annual survival of fully recruited fish. (2) stock_recruit - a brood table of spawners and the recruits they produced; question is productivity, density dependence, S_MSY, U_MSY, steepness. (3) surplus_production - an annual catch series plus an abundance index (CPUE or survey); question is MSY, B/B_MSY, F/F_MSY (Kobe status).
not_for: Length-only data (convert to age with an age-length key first, or use a length-based method), separating fishing from natural mortality (needs an M estimate or tagging), growth curves (use fisheries-growth), age-structured assessment (stock synthesis, SAM), catch-only methods without an index, stocks whose index is not proportional to biomass (hyperstability), recruitment series with errors-in-spawners problems, one-way-trip catch/CPUE histories for production models (weakly identified, see Known limits).
data_needed: catch_curve - one row per aged fish (age) or a catch-at-age table (age, count). stock_recruit - brood year, spawners, recruits. surplus_production - year, catch (every year), abundance index (may have gaps).
engine: catch_curve - FSA::catchCurve (weighted regression) + Chapman-Robson (hand code, checked against FSA) + glmmTMB Poisson-lognormal GLMM; stock_recruit - least squares (lm / optim) or nimble (Bayesian, AR(1) residuals); surplus_production - nimble (Bayesian state-space, Meyer & Millar 1999)
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
N_a = N_0 exp(-Z a), so log(catch) is linear in age with slope -Z. Three estimators
are reported:

- **Weighted linear regression** (default primary estimate). log(catch) on age from
  the modal (peak) age onwards, ages with zero catch excluded, with weights equal to
  the predicted log catch of an unweighted fit (Maceina & Bettoli 1998), using
  `FSA::catchCurve(weighted = TRUE)`. 95% t interval. Smith et al. (2012) found
  unweighted regression should no longer be used. Weighted regression is the
  primary estimate because it has the smallest bias and the best (though still
  below-nominal) interval coverage when year-class strength varies: about 0.91-0.93
  at the default design, 0.84-0.88 in an independent age-structured simulator (VALIDATION.md).
  All three catch-curve intervals are somewhat too narrow under year-class
  variability; treat the 95% interval as roughly a 90% interval.
- **Chapman-Robson** (Chapman & Robson 1960). Ages are recoded x = 0, 1, ... from
  the first age used. Then S = T / (n + T - 1), with T = sum(x * catch) and
  n = number of fish. Z = -log(S) - (n-1)(n-2) / (n (T+1)(n+T-1)) (the Hoenig et
  al. 1983 bias correction). The SE of Z is (1 - e^-Z) / sqrt(n e^-Z) (Smith et
  al. 2012 eq. 2, as in FSA), multiplied by sqrt(c-hat). c-hat is the Pearson
  chi-square of the geometric age distribution divided by its degrees of freedom,
  floored at 1. The interval uses a t-quantile on those degrees of freedom. The
  first age is peak age + 1 (the "peak plus" criterion, Smith et al. 2012).
- **Poisson-lognormal GLMM** (Millar 2015). catch_a ~ Poisson(exp(b0 - Z a + e_a)),
  with e_a ~ N(0, s^2) a year-class random effect, fitted on the same ages as
  Chapman-Robson with glmmTMB. Wald interval.

`settings: primary:` chooses which estimator feeds `total_mortality_Z`,
`annual_survival` and `annual_mortality`. The options are `weighted_regression`,
`chapman_robson` and `glmm`.

**Assumptions.** These are constant Z over the ages used, no trend in recruitment
over the cohorts sampled, full and equal selectivity from the first age used,
representative sampling and accurate ages.
- **Recruitment trend.** If recruitment has been increasing, Z is over-estimated;
  if it has been decreasing, Z is under-estimated. The negative control in
  VALIDATION.md shows +29% bias.
- **Dome-shaped selectivity.** Z is over-estimated.
- **Ageing error.** Old fish are typically under-aged, which flattens the curve.

**Data.** Give either one row per fish with an `age` column (other columns are
ignored, and rows with missing age are dropped), or a catch-at-age table with
`count_col`. Ages missing between the youngest and oldest are filled with 0.

```
fish_id,sample_date,basin,gear,total_length_mm,age
WAE25-0001,2025-09-22,North Basin,experimental gill net,412,3
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
  primary: weighted_regression    # | chapman_robson | glmm
  age_full: peak_plus             # first age for CR/GLMM: peak_plus | peak | <age>
  age_full_regression: peak       # first age for the regression
  age_last: null                  # oldest age used (default: oldest in data)
  overdispersion: true            # c-hat inflation of the Chapman-Robson SE
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
  # catch_curve: age_col (age), count_col (optional), species, gear
  # stock_recruit: year_col (brood_year), spawners_col, recruits_col, species, spawner_units, recruit_units
  # surplus_production: year_col (year), catch_col (catch), index_col (cpue), species, catch_units, index_name
figures: {...}                    # see skills/eco-figures
```

Examples are in `example/`:
- `config.yaml` (walleye catch curve)
- `config_stock_recruit.yaml` (sockeye brood table, Bayesian Ricker with AR(1))
- `config_surplus_production.yaml` (longline catch and CPUE)

## Outputs

- **Estimates table** (contract form, with rhat / ess_bulk for Bayesian fits).
  For the Nimble engines the `estimate` column is the **posterior median** (with
  95% equal-tailed credible interval), as in Meyer & Millar (1999) and as plotted
  in the figures; posterior means of K, steepness and the SDs are pulled up by
  their right tails (independent verification, VALIDATION.md).
- **Diagnostics notes:**
  - catch_curve: ages used, c-hat, estimator disagreement.
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
| catch_curve | `catch_curve` | Number caught at each age on a log axis (fish, not log units); filled circles = ages used by Chapman-Robson and the regression, open squares = regression only, crosses = not fully recruited (excluded). Chapman-Robson expected catch (solid) and weighted regression (dashed) with its confidence band. |
| catch_curve | `mortality_estimates` | Z (per year) with 95% CI for Chapman-Robson, weighted regression and the Poisson-lognormal GLMM; dashed line = primary estimate. Estimator colours match `catch_curve`. |
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
- Estimator colours come from `palette` in a fixed order (Chapman-Robson, weighted regression, GLMM), so they stay the same across figures.
- `labels:` keys:
  - axes: `age`, `log_catch`, `Z`, `spawners`, `recruits`, `year`, `resid`, `biomass`, `biomass_catch`, `surplus`, `index`, `B_over_BMSY`, `F_over_FMSY`
  - estimators and age classes: `chapman_robson`, `weighted_regression`, `glmm`, `used_cr`, `used_rg`, `not_used`, `estimator_title`, `ages_title`
  - annotations: `replacement`, `S_MSY_lab`, `B_MSY_lab` (plotmath), `kobe_red`, `kobe_green`, `kobe_yellow`, `kobe_orange`
  - parameters: `ln_alpha`, `beta`, `phi`, `sigma`, `sigma_nu`, `r`, `K`, `q`, `log10_r`, `log10_K`, `log10_MSY`
- `per_figure: {<id>: {width, aspect, height}}` overrides the default sizes. The module sets a height for each figure from its aspect, so the user's values take precedence.
- Units come from `data:` (`spawner_units`, `recruit_units`, `catch_units`, `index_name`). The module does not standardise covariates, so `data$scaling` is not needed.

## Simulation and power

`design: model:` picks the simulator.

**catch_curve**
- Params: `Z`, `recruit_sd_cc` (lognormal year-class SD), `recruit_ar1_cc`,
  `recruit_trend_cc` (log change per cohort), `sel_age50`, `sel_slope`,
  `Z_age_slope`, `sel_dome_age` / `sel_dome_rate`.
- Design: `n_fish`, `max_age`.

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
- Millar RB (2015) A better estimator of mortality rate from age-frequency data. Can J Fish Aquat Sci 72:364-375.
- Ogle DH, Doll JC, Wheeler AP, Dinno A (2025) FSA: Simple Fisheries Stock Assessment Methods. R package 0.10.1.
- Polacheck T, Hilborn R, Punt AE (1993) Fitting surplus production models: comparing methods and measuring uncertainty. Can J Fish Aquat Sci 50:2597-2607.
- Ricker WE (1954) Stock and recruitment. J Fish Res Board Can 11:559-623.
- Schaefer MB (1954) Some aspects of the dynamics of populations important to the management of the commercial marine fisheries. IATTC Bull 1:27-56.
- Scheuerell MD (2016) An explicit solution for calculating optimum spawning stock size from Ricker's stock recruitment model. PeerJ 4:e1623.
- Smith MW, Then AY, Wor C, Ralph G, Pollock KH, Hoenig JM (2012) Recommendations for catch-curve analysis. N Am J Fish Manage 32:956-967.
- Walters CJ (1985) Bias in the estimation of functional relationships from time series data. Can J Fish Aquat Sci 42:147-149.
