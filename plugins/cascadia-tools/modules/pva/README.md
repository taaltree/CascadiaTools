---
id: pva
title: Population viability analysis - count-based (diffusion / state-space) and matrix-model PVA
when_to_use: The question is extinction or quasi-extinction RISK - "how likely is this population to fall below N animals within 10, 25 or 50 years, and which vital rates or management actions change that?" Use the count model when you have one count or index per year (10+ years, ideally 20+). Use the matrix model when you have stage- or age-specific vital rates (survival, fecundity) with standard errors and year-to-year variation, e.g. from a CJS, IPM or nest-monitoring analysis, or raw annual monitoring counts of animals at risk and survivors / young produced.
not_for: Estimating trend, density dependence or covariate effects on growth from counts (use timeseries, which this module complements); estimating survival from capture histories (cjs / multistate / ipm first, then feed their estimates here); spatially explicit or metapopulation PVA (not supported); fewer than ~10 counts (risk is unestimable); populations with strong density dependence near carrying capacity in the count model (see Pitfalls).
data_needed: count model - CSV with one row per year (year, count; blank = not surveyed). matrix model - a matrix definition (stages, entries as formulas in named vital rates) in the config plus EITHER a vital-rate CSV (rate, type, mean, se, sd_env, n_years) OR an annual monitoring CSV (year, rate, n, events); initial abundance per stage; a quasi-extinction threshold.
engine: hand-written R - Dennis et al. diffusion estimators and inverse-Gaussian first-passage risk with a discrete-census correction; Bayesian exponential-growth state-space (EGSS) model computed exactly on a grid over the REML likelihood (no MCMC); vectorised stochastic matrix projection with demographic and environmental stochasticity and parameter uncertainty. popbio, MARSS and a Nimble MCMC of the same state-space model used as independent cross-checks.
---

# Population viability analysis (PVA)

## Question it answers

"How likely is this population to drop to a level we would consider
functionally extinct (the **quasi-extinction threshold**) within the next 10,
25 or 50 years; how fast is it growing or declining in the long run; and
which vital rates, or which management actions, would most reduce the risk?"

Two analyses, chosen with `settings: model:`

| model | data | what it estimates |
|---|---|---|
| `count` (`method: diffusion`) | annual counts treated as exact | long-run growth rate mu (log scale) and lambda = exp(mu), process variance sigma^2, quasi-extinction probability by each horizon (with bootstrap CIs) and the full risk curve |
| `count` (`method: state_space`) | annual counts with counting error, missing years allowed | the same plus observation-error variance tau^2 and the true current abundance |
| `matrix` | stage-structured vital rates (estimates with SEs, or raw monitoring) | deterministic and stochastic lambda (simulation and Tuljapurkar approximation), quasi-extinction curves with uncertainty, vital-rate and matrix-entry elasticities, management-scenario comparisons with paired differences |

### Estimands (names in the estimates table)

| name | meaning | model |
|---|---|---|
| `growth_rate_mu` | mean annual log growth rate, mu (power is computed for this; null 0) | count |
| `lambda_stochastic` | long-run (geometric-mean) growth rate per year: exp(mu) for counts; simulated stochastic lambda for matrices (null 1) | both |
| `process_variance_sigma2` | variance of the annual log growth rate caused by real environmental variation | count |
| `obs_variance_tau2` | counting-error variance on the log scale | count, state_space |
| `N_current` | (state-space) estimated true abundance in the last year; (diffusion) the last count, taken as known | count |
| `p_qext_<h>yr` | probability that abundance is at or below the threshold in at least one annual census within h years | both |
| `<rate>`, `sd_env_<rate>` | mean vital rate and its environmental (year-to-year) SD | matrix |
| `lambda_deterministic` | dominant eigenvalue of the mean matrix | matrix |
| `lambda_stochastic_tuljapurkar` | small-variance approximation to stochastic lambda (Tuljapurkar 1982/1990) | matrix |
| `elasticity_<rate>` | proportional change in deterministic lambda per proportional change in the rate | matrix |
| `scenario_<name>_lambda_stochastic`, `scenario_<name>_p_qext_<h>yr`, `scenario_<name>_change_p_qext_<h>yr` | the same under a management scenario, and the paired change in risk relative to baseline | matrix |

## Model

### Count-based PVA (Dennis, Munholland & Scott 1991; Morris & Doak 2002 ch. 3-5)

True log abundance follows a random walk with drift:
x(t+s) = x(t) + mu*s + e, e ~ Normal(0, sigma^2 * s).

* **diffusion**: counts are exact, y = x. With q intervals of lengths s_i and
  log ratios r_i, mu-hat = sum(r_i)/sum(s_i), sigma^2-hat from the regression
  of r_i/sqrt(s_i) on sqrt(s_i) through the origin (q-1 df). CI for mu is t-based,
  CI for sigma^2 chi-square based (Dennis et al. 1991).
* **state_space** (EGSS; Dennis et al. 2006; Humbert et al. 2009):
  y(t) = x(t) + f, f ~ Normal(0, tau^2), unsurveyed years simply absent.
  Bayesian, with flat priors on the initial log abundance and mu and uniform
  priors on the SDs sigma and tau (Gelman 2006). Flat priors on the two fixed
  effects make the marginal likelihood of (sigma^2, tau^2) exactly the REML
  likelihood (Harville 1974; the same likelihood Staples, Taper & Dennis 2004
  maximise), so p(sigma, tau | y) is computed on a 60 x 60 grid (after a
  coarse pass that finds the mass), and (x1, mu) and all true log abundances
  are drawn from their exact Gaussian conditionals. Estimates are posterior
  medians with 95% credible intervals; risk is computed per draw, so its
  interval includes uncertainty in the current true abundance. Why not REML +
  bootstrap (the previous version of this module): in ~10% of 25-year series
  REML puts sigma^2 at 0, a bootstrap from that fit has no process noise and
  the intervals collapse (coverage of mu 0.90, of risk 0.87-0.88 at nsim 200;
  VALIDATION.md). The REML estimates are still computed and shown in the
  methods text for reference.

**Quasi-extinction probability** within t years from current log distance
d = log(N_current / N_threshold) is the inverse-Gaussian first-passage CDF
(Dennis et al. 1991 eq. 16; Holmes 2004):

G(t) = Phi((-d - mu t)/(sigma sqrt t)) + exp(-2 mu d / sigma^2) Phi((-d + mu t)/(sigma sqrt t)).

That formula assumes the population is watched continuously and overstates
the risk of falling below the threshold at an annual census. By default
(`correction: discrete`) the barrier is moved out by 0.5826*sigma (Broadie,
Glasserman & Kou 1997 continuity correction), which we checked against
brute-force simulation of annual censuses (VALIDATION.md); `correction: none`
gives the classic formula, which is what popbio::countCDFxt computes.
Uncertainty is propagated by recomputing G(t) for parameter sets drawn from
the sampling distributions of (mu, sigma^2) (and of the current state for
state_space): diffusion uses the exact chi-square / normal sampling
distributions (as in Morris & Doak 2002 Box 3.4 and popbio::countCDFxt), state-space
uses the posterior draws.

### Matrix-model PVA (Caswell 2001; Morris & Doak 2002 ch. 7-9)

n(t+1) = A(t) n(t), entries of A written as formulas in named vital rates
(e.g. `fec * s_juv`). Each year:

* **environmental stochasticity**: every vital rate is redrawn - survival and
  other probabilities from a beta with the estimated mean and environmental
  SD; fecundities from a lognormal (default) or stretched beta on [min, max]
  (Morris & Doak 2002 ch. 8). Optional within-year correlation among rates
  (one number, or a correlation-matrix CSV) via a Gaussian copula.
* **demographic stochasticity** (`demographic: true`, default): survivors and
  stage transitions binomial (sequentially multinomial when a stage has several
  fates), offspring Poisson.
* **density dependence** (optional): a ceiling on the counted stages; above it
  every stage is thinned proportionally (binomially).
* **parameter uncertainty**: `n_param_draws` sets of mean rates are drawn from
  their sampling distributions (beta / lognormal with the given SE) and the
  total (environmental + sampling) variance from the chi-square pivot,
  s^2 (n_years-1) / chi^2(n_years-1) (a scaled inverse chi-square; sampling
  variance then removed when known); `n_traj_per_draw` trajectories each. Intervals
  for risk and lambda are percentiles across parameter draws; point estimates
  use the estimated rates with `n_traj_point` trajectories.
* **stochastic lambda**: (i) simulation - mean annual log growth of an
  environment-only projection after a 10-year burn-in; (ii) Tuljapurkar's
  approximation log lambda_s ~ log lambda_1 - tau^2/(2 lambda_1^2), tau^2 =
  sum of sensitivities x covariances of matrix entries (Tuljapurkar 1990;
  Caswell 2001 eq. 14.70).
* **elasticities**: matrix-entry elasticities (Caswell 2001) and lower-level
  vital-rate elasticities by finite differences (agree with popbio::vitalsens).
* **scenarios**: `multiply` or `add` to mean vital rates; scenarios reuse the
  baseline random numbers, so the change in risk is a paired difference with
  its own interval.

When raw monitoring data are given (`demography:`), each rate's annual
estimate is events/n; the mean is the average over years, its SE the SD of
annual estimates / sqrt(years) (inflated by t(0.975, years-1)/1.96 so that
intervals are t-type in short series), and the environmental variance is the variance
of the annual estimates minus their mean binomial (survival) or Poisson
(fecundity) sampling variance (moment estimator; Kendall 1998; Morris & Doak
2002 ch. 8).

## Assumptions (and what breaks)

* **Constant conditions**: mean growth and variability in the future equal
  those in the data. Trends in habitat or harvest make every number wrong.
* **No density dependence** (count model; matrix model unless `ceiling`).
  Near carrying capacity, density-independent PVA tends to overstate long-term
  risk (Sabo et al. 2004). The count model reports a density-dependence check.
* **Independent years** (no autocorrelation in growth). Positive
  autocorrelation (runs of bad years) makes true risk higher than estimated.
* **Counts are exact** (diffusion). Counting error inflates sigma^2 and hence
  the risk; use `method: state_space` (our negative control shows the effect).
* **No catastrophes beyond those in the data.**
* **Count model ignores stage structure and demographic stochasticity**; at
  very small N (tens of animals) prefer the matrix model with `demographic: true`.
* **Matrix**: vital rates vary around constant means; female-only
  (or one-sex) model; the matrix definition (pre- vs post-breeding census)
  is correct. Elasticities are for deterministic lambda.

## Data format

### Count model
```
year,adult_females
1959,44
1960,47
1961,46
```
Column names are set with `data: year_col`, `count_col`. Blank cells = year not
surveyed (handled exactly by state_space; diffusion uses the longer interval).
Zeros cannot be used on the log scale: remove them or reconsider count PVA.

### Matrix model, option A: vital-rate estimates (e.g. from cjs / ipm output)
```
rate,type,mean,se,sd_env,n_years,description
s_juv,survival,0.32,0.035,0.08,12,first-year survival (CJS)
s_ad,survival,0.84,0.018,0.035,12,adult survival (CJS)
fec,fecundity,0.95,0.07,0.28,12,female fledglings per breeding female
```
`type`: survival (any probability: survival, breeding probability,
transition) or fecundity. `se`: standard error of the mean (parameter
uncertainty; 0 = known). `sd_env`: year-to-year (process) SD, with sampling
variance already removed (e.g. a random-effect SD from a CJS/IPM fit).
`n_years`: years behind sd_env (gives the uncertainty of sd_env; blank = treat
as known). Optional: `sampling_var`, `dist` (beta, lognormal, stretched_beta),
`min`, `max` (for stretched beta).

### Matrix model, option B: raw annual monitoring
```
year,rate,n,events
2009,s_juv,60,25      # 60 marked fledglings, 25 alive next year
2009,s_ad,80,61
2009,fec,40,72        # 40 breeding females produced 72 female young
```
Survival rows: events = survivors among n at risk (assumes detection = 1, e.g.
known-fate radio-tracking or resighting with perfect detection; otherwise fit
cjs first and use option A). Fecundity rows: events = offspring of n breeders.

## Config

```yaml
settings:
  model: count | matrix
  threshold: 20                # quasi-extinction threshold (same units as counts / counted stages)
  horizons: [10, 25, 50]       # years ahead at which risk is reported
  # count
  method: diffusion | state_space
  correction: discrete | none  # annual-census correction (default) or classic continuous formula
  nboot: 4000                  # bootstrap parameter sets (diffusion)
  ndraw_ss: 4000               # posterior draws (state_space; ~1 s per fit)
  # matrix
  demographic: true
  n_param_draws: 400
  n_traj_per_draw: 250
  n_traj_point: 10000
  ignore_parameter_uncertainty: false
  scenarios:
    - name: Juvenile survival +10%
      multiply: {s_juv: 1.10}
    - name: Adult survival +0.02
      add: {s_ad: 0.02}
  seed: 1
data:
  # count
  file: counts.csv
  year_col: year
  count_col: count
  # matrix
  vital_rates: vital_rates.csv     # or demography: monitoring.csv + rate_types: {s_juv: survival, fec: fecundity}
  stages: [juvenile, subadult, adult]
  matrix:                          # rows = stage next year, columns = stage this year
    - ["0", "b_sub * fec * s_juv", "fec * s_juv"]
    - ["s_sub", "0", "0"]
    - ["0", "s_ad", "s_ad"]
  initial_abundance: {juvenile: 40, subadult: 30, adult: 110}
  count_stages: [adult]            # stages counted against threshold and ceiling (default all)
  ceiling: 300                     # optional ceiling density dependence
  env_correlation: 0.3             # or a CSV correlation matrix with rate names as row/col names
  species: the population
  units: adult females
figures:
  journal: esa
```

Examples: `example/config.yaml` (Yellowstone grizzly, diffusion),
`example/config_state_space.yaml` (shorebird survey with counting error and
missing years), `example/config_matrix.yaml` (3-stage bird, vital-rate table,
three management scenarios), `example/config_matrix_monitoring.yaml` (rates
estimated from raw monitoring).

## Outputs

Estimates table (above), HTML report with plain-language interpretation and
methods paragraph, and figures (docs/FIGURES.md): `quasi_extinction`
(cumulative risk curve with 50%/95% bands), `trajectories` (fan chart of
projected abundance; counts and state-space estimates for the count model),
`lambda_distribution`, `growth_checks` (count model: growth-rate time series
and density-dependence check), `elasticity`, `elasticity_matrix`, `scenarios`
(matrix model).

## Simulation and power

`design: model: count` simulates a random walk (params `mu`, `sigma2`, `tau2`)
for `n_years` from `N0`, with `p_missing` share of skipped years; the true
risk for each replicate is computed from the TRUE parameters and the TRUE
final state by an exact numerical first-passage recursion for annual censuses
(`pva_qext_exact`, agrees with 1e6-path brute force to 4.8e-4; set
`truth_method: simulation` for brute force with `truth_paths` paths). `design: model: matrix` simulates
`n_mon_years` of monitoring of a 2-stage model (params `s_juv, s_sub, s_ad,
fec, sd_*`; design `n_marked_*`, `n_breeders`, `n0`, `threshold`, `ceiling`),
with risk truth from 1e5 projections at the true rates. Useful power
questions: how many years of counts are needed to tell a 3%/yr decline from
stability (`--vary n_years=...`), and how many marked animals to estimate
stochastic lambda precisely.

## Pitfalls

* **Long-horizon risk is imprecise.** With 20-40 years of data, risk beyond
  10-20% of the series length is usually unestimable: intervals span most of
  0-1 (Ludwig 1999; Fieberg & Ellner 2000; Ellner & Holmes 2008). The module
  flags horizons beyond that and the report tells the reader to read the
  interval, not the point. Use PVA to compare options, not to predict an
  extinction date.
* **Point estimates of tiny risks are biased and skewed.** A risk of 0.1% with
  interval 0% to 22% (grizzly, 50 yr) means "low if the past continues, but not
  ruled out".
* **Counting error** makes the diffusion model overstate both sigma^2 and risk.
  If counts are surveys rather than censuses, use `method: state_space`.
* **Threshold choice** drives everything; report results for more than one.
* **Matrix elasticities** say where a given proportional change matters most,
  not what management can achieve.
* **Survival from raw monitoring** assumes perfect detection; otherwise
  estimate survival with cjs and supply estimates (option A).

## Validation status (details in VALIDATION.md)

* **PASS** (nsim = 200, default designs): growth rate mu / lambda, process and
  observation variances, current true abundance (state-space), all matrix vital
  rates, environmental SDs, deterministic and stochastic lambda, elasticities.
* **OPEN: quasi-extinction point estimates.** Their 95% intervals are calibrated
  against the true risk (coverage 0.94-1.00 at default designs), and the point
  estimates are median-unbiased. However, they are more extreme than the truth,
  so they fail the mean-bias gate (e.g. 25-yr risk: mean predicted 0.40 vs true
  0.34 for counts, 0.22 vs 0.13 for the matrix model). Report the interval.
* Negative controls fail as they should: diffusion fitted to data with counting
  error inflates sigma^2 by +209% (coverage 0.015); ignoring parameter
  uncertainty gives coverage below 0.06.
* False-positive rate for a declining trend: 0.055 (diffusion), 0.04 (state-space), 0.00 (matrix).
* Independent verification (own simulators and truth, new designs: irregular
  survey intervals, 20-yr state-space with tau2 > sigma2, 3-stage vital-rate-table
  input with stretched-beta fecundity): the same pattern held. Two caveats:
  - At a 5-yr horizon, diffusion risk-interval coverage was 0.88-0.89, almost
    entirely from misses at true risks < 0.3% (reported as "<0.1%").
  - With a vital-rate table the module trusts the supplied sd_env and SEs.
    Small environmental SDs estimated from about 12 years are biased low
    (-15 to -18%), and z-type SEs give vital-rate coverage of only 0.91-0.93
    (inflate SEs by t(0.975, years-1)/1.96 for short series).
* Cross-checks:
  - popbio (countCDFxt, vitalsens, stoch.growth.rate, stoch.quasi.ext);
  - MARSS;
  - a Nimble MCMC of the state-space model (posteriors agree within 0.03 posterior SD median, 0.28 max);
  - Morris & Doak (2002) grizzly estimates.

## References

* Broadie M, Glasserman P, Kou SG (1997) A continuity correction for discrete barrier options. Mathematical Finance 7:325-349.
* Caswell H (2001) Matrix Population Models, 2nd ed. Sinauer.
* Dennis B, Munholland PL, Scott JM (1991) Estimation of growth and extinction parameters for endangered species. Ecological Monographs 61:115-143.
* Dennis B, Ponciano JM, Lele SR, Taper ML, Staples DF (2006) Estimating density dependence, process noise, and observation error. Ecological Monographs 76:323-341.
* Ellner SP, Holmes EE (2008) Resolving the debate on when extinction risk is predictable. Ecology Letters 11:E1-E5.
* Fieberg J, Ellner SP (2000) When is it meaningful to estimate an extinction probability? Ecology 81:2040-2047.
* Holmes EE (2001) Estimating risks in declining populations with poor data. PNAS 98:5072-5077.
* Holmes EE (2004) Beyond theory to application and evaluation: diffusion approximations for population viability analysis. Ecological Applications 14:1272-1293.
* Humbert J-Y, Mills LS, Horne JS, Dennis B (2009) A better way to estimate population trends. Oikos 118:1940-1946.
* Kendall BE (1998) Estimating the magnitude of environmental stochasticity in survivorship data. Ecological Applications 8:184-193.
* Ludwig D (1999) Is it meaningful to estimate a probability of extinction? Ecology 80:298-310.
* Morris WF, Doak DF (2002) Quantitative Conservation Biology: Theory and Practice of Population Viability Analysis. Sinauer.
* Sabo JL, Holmes EE, Kareiva P (2004) Efficacy of simple viability models in ecological risk assessment: does density dependence matter? Ecology 85:328-341.
* Staples DF, Taper ML, Dennis B (2004) Estimating population trend and process variation for PVA in the presence of sampling error. Ecology 85:923-929.
* Stubben CJ, Milligan BG (2007) Estimating and analyzing demographic models using the popbio package in R. Journal of Statistical Software 22(11).
* Tuljapurkar S (1990) Population Dynamics in Variable Environments. Springer.
* Gelman A (2006) Prior distributions for variance parameters in hierarchical models. Bayesian Analysis 1:515-533.
* Harville DA (1974) Bayesian inference for variance components using only error contrasts. Biometrika 61:383-385.
* Eberhardt LL, Knight RR, Blanchard BM (1986) Monitoring grizzly bear population trends. Journal of Wildlife Management 50:613-618 (grizzly example data, via Morris & Doak 2002 Table 3.1 and popbio::grizzly).
