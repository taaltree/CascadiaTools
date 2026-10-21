---
id: timeseries
title: Population time series - trends, density dependence and drivers of growth (state-space models)
when_to_use: One count or abundance index per year (aerial/ground counts, nest counts, redd counts, CPUE, BBS-style route totals) for one population or several sub-populations, 10+ years; the question is whether the population is increasing or declining and by how much per year, how large it really is once counting error is removed, whether growth is density dependent, whether a covariate (winter severity, harvest, drought) is associated with good and bad years, and how many years of monitoring are needed to detect a decline.
not_for: Repeated counts within a season at many sites (use nmixture, which estimates detection); marked animals (cjs / jolly-seber / ipm give vital rates); detection/non-detection (occu-*); age- or stage-structured projections and extinction risk (pva); series shorter than ~8 years (no model can separate process from observation error).
data_needed: CSV with one row per year (per site): year, count or index, optional site/herd column, optional yearly covariates; years not surveyed may be blank or absent
engine: hand-written REML / Kalman-filter maximum likelihood (exponential, Gompertz, multi-site default); MARSS package (other multi-site structures); nimble MCMC (Bayesian with covariates and missing years); glmmTMB (log-linear trend with AR(1), comparison)
---

# Population time series: trends, density dependence and drivers of growth

## Question it answers

"Is my population declining, how fast, and what drives good and bad years?"
Annual counts bounce for two different reasons: the population really changes
(births, deaths, a hard winter: **process noise**) and the count misses or
double-counts animals (**observation error**). State-space models separate the
two. That matters because

* the trend's confidence interval is too narrow if all bounce is treated as
  counting noise around a line (log-linear regression), and too wide/biased if
  all of it is treated as real change (Humbert et al. 2009);
* density dependence estimated while ignoring counting error is **spuriously
  strong**: counting error makes a high count tend to be followed by a lower
  one (regression to the mean) even when there is no regulation
  (Dennis et al. 2006; Knape 2008);
* the estimated true trajectory (with CIs) is what managers should look at,
  not the raw counts.

Outputs (names in the estimates table):

| estimand | meaning | models |
|---|---|---|
| `trend_pct_per_year` | average % change per year, 100(exp(mu) - 1) | exponential, marss, loglinear |
| `growth_rate_mu` | the same on the log scale (mean of log(N[t]/N[t-1])); what power is computed for | exponential, marss, loglinear |
| `trend_pct_per_year_<site>`, `growth_rate_mu_<site>` | per-site trends (`trend: separate`) | marss |
| `realised_trend_pct_per_year` | average annual change of the estimated true population from the first to the last year | bayes |
| `growth_rate_r0` | mean log growth rate when there is no density dependence | bayes (`density_dependence: false`) |
| `N_final` | estimated true abundance (on the count scale) in the last year | exponential, gompertz, bayes |
| `sd_process` | SD of real year-to-year changes in log abundance | all state-space models |
| `sd_observation` | SD of counting error on the log scale (0.15 ~ a count typically 15% off) | all state-space models |
| `density_dependence_b` | Gompertz density dependence: 0 = none, -1 = full return to equilibrium each year | gompertz, bayes |
| `equilibrium_N` | abundance at which growth averages zero (carrying capacity in the Gompertz sense) | gompertz, bayes |
| `growth_beta_<cov>` | change in log growth rate per 1 SD of the covariate (link scale) | bayes |
| `growth_effect_pct_<cov>` | % change in the annual growth multiplier N[t]/N[t-1] per 1 SD | bayes |
| `count_beta_<cov>`, `count_effect_pct_<cov>` | effect of a covariate on the *level* of counts in a year (not on growth) | loglinear |
| `process_corr` | correlation of process noise between sites | marss (`process_cov: equalvarcov`) |
| `ar1_rho`, `sd_ar1` | AR(1) residual correlation and SD (no intervals) | loglinear |

Plus the estimated true trajectory with 95% bands (figure `trajectory`) and a
forecast (figure `forecast`).

## Models (settings: model)

All state-space models work on x_t = log N_t and y_t = log(count_t).

**exponential** (EGSS; Dennis et al. 2006; Humbert et al. 2009; Staples et al. 2004):

    x_t = x_{t-1} + mu + e_t,     e_t ~ N(0, sd_process^2)
    y_t = x_t + f_t,              f_t ~ N(0, sd_observation^2)

Default `estimation: reml`: restricted maximum likelihood on the first
differences of the log counts (Staples et al. 2004; Humbert et al. 2009), which
do not depend on the unknown starting abundance. The trend's interval
integrates the restricted likelihood over both variances (flat prior on the
SDs) so it stays honest when one variance is estimated at 0; SD intervals are
profile-REML intervals; the true trajectory comes from a diffuse-start Kalman
smoother with parameter uncertainty propagated by 1000 draws.
`estimation: ml` gives exact ML with a Kalman filter. Use REML: ML trend
intervals cover only 68% at 25 years (VALIDATION.md). REML trend intervals are
conservative (coverage 0.97-1.0, false-positive rate 0.025), which lowers power.

**gompertz** (GSS; Dennis et al. 2006; Knape 2008):

    x_t = x_{t-1} + b (x_{t-1} - log K) + e_t
    y_t = x_t + f_t

Exact ML (Kalman filter, initial state free one year before the first count),
profile-likelihood intervals for the parameters. Abundance (trajectory,
`N_final`) intervals propagate parameter uncertainty by adaptive
sampling-importance-resampling from the likelihood (flat priors), so they stay
valid when a variance sits on its boundary (coverage 0.92 at 30 years;
VALIDATION.md). **Testing density dependence:** the fit also runs
a parametric-bootstrap likelihood-ratio test of exponential growth with drift
(H0) against the Gompertz model (Dennis & Taper 1994; Dennis et al. 2006;
`dd_test_boot: 199` simulations, about 1-2 min). Report its p-value. Do not use
"the interval for b excludes 0": under b = 0 that happens in 47% of
30-year series, whereas the bootstrap test rejected in 4% (30 yr) and 7%
(15 yr) of them (VALIDATION.md). Its power is low: about 6% at 15 years for
b = -0.3 with realistic noise. Earlier tests that ignored sampling error were
shown to be invalid when counts are noisy (Shenk et al. 1998). This test
includes observation error in both models, and its false-positive rate is
checked in VALIDATION.md.

`settings: observation_error: false` fits either model without observation
error (the classic diffusion/AR(1) regression). It exists as the negative
control: it overestimates density dependence.

**marss** (Holmes, Ward & Wills 2012; MARSS User Guide): several sites,

    x_t = x_{t-1} + u + w_t,  w_t ~ MVN(0, Q);   y_t = Z x_t + a + v_t,  v_t ~ MVN(0, R)

`trend: shared` (one u) or `separate` (u per site); `states: separate` (each
site its own population, Z = I) or `single` (all sites are replicate surveys
of one population; Z = column of ones, site offsets a); `process_cov:`
`diagonal and equal` (default), `diagonal and unequal`, `equalvarcov`
(correlated sites, estimates `process_corr`). R is diagonal and equal.
Estimation: with the default structure (`states: separate`, `process_cov:
diagonal and equal`) sites are independent given the two variances, so the
module uses pooled REML over each site's differenced log counts, the
multi-site version of the exponential model, with the same integrated
intervals (`estimation: reml`, default). All other structures, or
`estimation: ml`, use ML in the MARSS package (BFGS; EM for `equalvarcov`)
with Hessian intervals. These under-cover (trend 0.88, process SD 0.79 at 4
sites x 25 years; VALIDATION.md).

**bayes** (Kery & Schaub 2012 ch. 5; NIMBLE, de Valpine et al. 2017):

    x_t = x_{t-1} + r0 + b (x_{t-1} - xbar) + sum_k beta_k z_{k,t} + e_t
    y_t ~ N(x_t, sd_observation^2)   for surveyed years only

xbar = mean log count (centring). z are standardised covariates; the value in
year t's row acts on growth from t-1 to t (e.g. snow depth in the winter that
ends with the year-t count). Missing years are latent states.
`density_dependence: false` drops b (exponential growth with covariates).
Equilibrium (at average covariates) log K = xbar - r0/b.
Priors: r0 ~ N(0, 1); b ~ U(-2, 0); beta_k ~ N(0, 1) per SD;
sd_process ~ half-normal(0, 1), sd_observation ~ half-normal(0, 1)
(`prior_sd_process`, `prior_sd_observation` = c(mean, sd) of a normal
truncated at 0: an informative observation-error prior from replicate counts or
sightability trials is the standard cure for the Knape problem); x_1 ~ N(first
log count, sd 2). Samplers: AF_slice on all top-level parameters + RW_block on
the two SDs; latent states conjugate. Defaults 4 chains x 40,000 iterations
(10,000 burn-in, thin 10).

**loglinear** (comparison; Brooks et al. 2017):

    count_{s,t} ~ NegBin(lambda, theta);  log lambda = alpha + alpha_s + mu (t - tbar) + sum beta_k z_k,t + AR(1)_s,t

`family: nbinom2` (counts; default), `poisson`, `lognormal` (positive indices,
default when `response: index`). `ar1: false` drops the AR(1) term. Covariates
shift the level of counts in that year; they are not effects on growth.

## Assumptions (and what breaks)

* **Constant detection/sightability.** Abundance is on the scale of the count.
  If the fraction counted changes over time (new aircraft, observer, method),
  the trend is biased. Add a covariate for the change or split the series.
* **Lognormal counting error, independent between years.** Very low counts
  (< ~5) and zeros violate this; use `loglinear` (negative binomial) or
  `zero_offset: 0.5` (biased for small counts).
* **Gaussian process noise** with constant variance; catastrophes (one crash)
  inflate sd_process and can pull the trend.
* **Constant growth rate** (exponential, marss) or constant equilibrium
  (gompertz). A long series with a change in trend needs splitting.
* **Density dependence acts on log abundance (Gompertz).** b is weakly
  identified when observation error is also estimated (Knape 2008); see
  VALIDATION.md before interpreting b.
* **Covariates (bayes)** are associations with growth, not causes.

## Data format

`example/herd_counts.csv` (single population, covariate, 3 winters not flown):

    year,survey,animals_counted,snow_depth_cm,notes
    1994,aerial winter count,590,83,
    ...
    2001,aerial winter count,,62,not flown (weather)

`example/herds_multi.csv` (several sub-populations):

    year,herd,count
    2000,North Fork,412

* one row per year (per site); years with no survey: leave the count blank or
  omit the row (the model inserts the gap);
* counts must be > 0 for the log-scale models (0 -> use loglinear);
* `model: bayes` with covariates needs covariate values in every year after the
  first, including years without a count;
* seasonal/replicate counts must first be combined to one value per year.

## Config

```yaml
data:
  file: herd_counts.csv         # relative to the config
  year: year                    # column names
  count: animals_counted
  site: herd                    # optional; omit for one population
  covariates: [snow_depth_cm]   # optional; bayes (growth) or loglinear (level)
  response: count               # count | index
  count_label: "Animals counted"
  covariate_labels: {snow_depth_cm: "Mean winter snow depth (cm)"}
settings:
  model: bayes                  # exponential | gompertz | marss | bayes | loglinear
  estimation: reml              # exponential & default marss: reml (default) | ml
  observation_error: true       # false = negative-control model without counting error
  zero_offset: 0                # added to counts before logging
  # marss
  trend: shared                 # shared | separate
  states: separate              # separate | single
  process_cov: diagonal and equal   # diagonal and unequal | equalvarcov
  # gompertz
  dd_test_boot: 199             # bootstrap simulations for the density-dependence test (0 = skip)
  # bayes
  density_dependence: true
  prior_sd_process: [0, 1]      # half-normal (mean, sd)
  prior_sd_observation: [0, 1]
  niter: 40000; nburnin: 10000; thin: 10; nchains: 4
  # loglinear
  family: nbinom2               # nbinom2 | poisson | lognormal
  ar1: true
figures:
  power_csv: results/timeseries/power_5pct/timeseries_power.csv   # adds figure power_years
  forecast_years: 10
```

Other example configs: `config_exponential.yaml`, `config_gompertz.yaml`,
`config_marss.yaml`, `config_loglinear.yaml`.

## Model checking (built in)

* exponential, gompertz, default marss: standardised one-step-ahead
  prediction errors (Kalman innovations) at the estimates. Ljung-Box test for
  autocorrelation and Shapiro-Wilk test for normality, in the diagnostics
  (`gof_ljung_box_p`, `gof_shapiro_p`), the interpretation and figure
  `residuals`.
* bayes: posterior predictive check of the SD and the lag-1 autocorrelation
  of year-to-year changes in log counts (`ppc_p_sd_change`,
  `ppc_p_acf1_change`; figure `ppc`).
* loglinear: none built in. Run DHARMa on `fit$object$glmm`.

## Outputs

Estimates table (above), plain-language interpretation, and figures
`trajectory`, `trend`, `density_dependence` (gompertz/bayes),
`growth_<cov>` (growth multiplier vs covariate in original units, bayes),
`count_<cov>` (loglinear), `coefficients`, `variances`, `forecast`,
`power_years` (when `power_csv` is set), `residuals` (state-space ML/REML),
`ppc`, `trace`, `posterior` (bayes).

## Simulation and power

`design: model` chooses the simulator. Knobs: `n_years` (default 25),
`frac_missing` (share of interior years without a count), `start_N`,
`n_sites` (marss, loglinear), `density_dependence` (bayes),
`n_covariates`, `process_cov` (marss), `first_year`.
Params: `growth_rate_mu` (-0.03 = -3%/yr), `sd_process` (0.1),
`sd_observation` (0.15), `density_dependence_b` (-0.3), `equilibrium_N` (500),
`growth_beta` (-0.15 per SD), `count_beta`, `process_corr`, `sd_site`,
`ar1_rho`, `sd_ar1`, `mean_count`, `nb_size`.

Key manager question - years of monitoring needed to detect a decline:

    scripts/Rscript core/run.R power timeseries --config my_design.yaml --vary n_years=8,10,15,20,25,30 --nsim 100

with `my_design.yaml` holding `design: {model: exponential}` and
`params: {growth_rate_mu: <log(1 + decline)>, sd_process: ..., sd_observation: ...}`
set from a pilot fit. Power = share of simulated series whose 95% interval for
the trend excludes 0. At sd_process 0.1 and sd_observation 0.15 a single
series declining 5%/yr is detected with power 0.29 after 20 years, 0.49 after
30 and 0.71 after 40; four sub-populations with a shared trend reach 0.74 for a
3%/yr decline in 25 years (single site: 0.15). Full tables are in VALIDATION.md.

## Validation status (VALIDATION.md)

**OPEN, estimand-specific.**
- **Recovered with honest intervals:**
  - trend (`exponential` REML at 12-25 years, `marss`, `loglinear` at 5 x 20);
  - true abundance (`N_final`);
  - covariate effects on growth and equilibrium abundance (`bayes`).
- **Not reliably recovered from one count series:**
  - process SD;
  - observation SD (Gompertz, Bayes);
  - density dependence `b` at realistic noise ratios (it was recovered by
    Gompertz ML only when process noise exceeded counting error and b = -0.5,
    40 yr).

Cross-checks agree with MARSS, brute-force integration, glm.nb and Nimble to
numerical precision.

## Pitfalls

* A variance estimated at **zero** (ML boundary) is common at < 20 years. It
  means the model could not split the noise, not that there is no counting
  error; trend and abundance intervals are then too narrow. The fit flags it.
  `model: bayes` averages over both explanations.
* **Density dependence**: do not report b from a single short series as
  evidence for or against regulation. In simulations b was overestimated by
  30-140% (Gompertz ML and Bayes alike; Knape 2008; Auger-Methe et al. 2016).
  In `bayes` the U(-2, 0) prior means its interval never contains 0. Use the
  bootstrap test of `model: gompertz`, and an informative observation-error
  prior if one exists.
* **"No clear trend" from one series is weak evidence.** Default REML trend
  intervals stay valid but are conservative. When counting error dominates
  (20 yr, sd_obs 0.25, sd_proc 0.06), a +4%/yr trend was never detected,
  whereas known-variance GLS would detect it 60% of the time (VALIDATION.md,
  Independent verification). Always run a power analysis at the study's design.
* **loglinear under real population wandering:** if the population follows
  a random walk with drift, its trend interval covers only ~83%. The report
  warns when AR(1) rho > 0.9.
* **Variance components** from a single series are biased (process SD -17%
  with REML at 25 years, +31% in `bayes`); trend and abundance intervals stay
  valid.
* A strong trend in the log-linear model with `pdHess` false: the AR(1)
  variance collapsed to 0 (competing with NB overdispersion); fixed-effect
  intervals are still usable but check against `exponential`.
* The realised trend (bayes) describes the past; with density dependence it is
  not a rate to project forward.
* Trend over the last few years is much less certain than over the series.

## References

* Dennis B, Taper ML (1994) Density dependence in time series observations of natural populations: estimation and testing. Ecological Monographs 64:205-224.
* Dennis B, Ponciano JM, Lele SR, Taper ML, Staples DF (2006) Estimating density dependence, process noise, and observation error. Ecological Monographs 76:323-341.
* Shenk TM, White GC, Burnham KP (1998) Sampling-variance effects on detecting density dependence from temporal trends in natural populations. Ecological Monographs 68:445-463.
* Knape J (2008) Estimability of density dependence in models of time series data. Ecology 89:2994-3000.
* Auger-Methe M, Field C, Albertsen CM, Derocher AE, Lewis MA, Jonsen ID, Mills Flemming J (2016) State-space models' dirty little secrets: even simple linear Gaussian models can have estimation problems. Scientific Reports 6:26677.
* Humbert J-Y, Mills LS, Horne JS, Dennis B (2009) A better way to estimate population trends. Oikos 118:1940-1946.
* Staples DF, Taper ML, Dennis B (2004) Estimating population trend and process variation for PVA in the presence of sampling error. Ecology 85:923-929.
* Holmes EE, Ward EJ, Wills K (2012) MARSS: multivariate autoregressive state-space models for analyzing time-series data. The R Journal 4(1):11-19.
* Holmes EE, Scheuerell MD, Ward EJ. Analysis of multivariate time series using the MARSS package (MARSS User Guide). NOAA Fisheries, Seattle.
* Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS. Academic Press, ch. 5.
* Auger-Methe M, Newman K, Cole D, et al. (2021) A guide to state-space modeling of ecological time series. Ecological Monographs 91:e01470.
* Brooks ME, Kristensen K, van Benthem KJ, et al. (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9:378-400.
* de Valpine P, et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. J. Computational and Graphical Statistics 26:403-413.
* Harvey AC (1989) Forecasting, Structural Time Series Models and the Kalman Filter. Cambridge University Press.
