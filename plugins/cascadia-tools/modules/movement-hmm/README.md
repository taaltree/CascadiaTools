---
id: movement-hmm
title: Hidden Markov models for animal movement (behavioural states)
when_to_use: GPS / Argos-quality-filtered telemetry with fixes at a (roughly) regular interval, and the question is which movement modes an animal switches between (e.g. encamped/foraging vs travelling), how much time it spends in each (activity budget), and how time of day or habitat at its location is associated with switching between modes.
not_for: Habitat selection / resource selection (use movement-ssf); home-range size (movement-ssf home ranges / ctmm); fixes so irregular or error-prone that a regular grid is meaningless (fit a continuous-time model first, e.g. crawl multiple imputation); accelerometer or dive data (possible with momentuHMM but not built in here); individual random effects on state switching; hidden semi-Markov (non-geometric dwell times).
data_needed: one CSV with one row per fix - animal id, timestamp, coordinates (lon/lat or projected x/y in metres) - plus optional covariate columns measured at each fix (e.g. % forest cover, distance to road)
engine: momentuHMM (maximum likelihood via nlm, multiple random starts); moveHMM and a hand-written forward-algorithm likelihood used for cross-checks
---

# Hidden Markov models for animal movement

## Question it answers

"When is my animal resting/foraging and when is it travelling, how much
time does it spend in each mode, and is switching between modes associated
with time of day and habitat?"

Each regular time step of a track is summarised by a **step length**
(distance between consecutive fixes) and a **turning angle** (change of
direction). An HMM assumes the animal is in one of *N* unobserved
behavioural states, each with its own distribution of step lengths and
turning angles, and that it switches between states according to a Markov
chain whose transition probabilities may depend on covariates (Langrock et
al. 2012; McClintock et al. 2020). The module reports:

| Output | Name in the estimates table | Scale |
|---|---|---|
| Mean and SD of step length in each state | `step_mean_m_state<k>`, `step_sd_m_state<k>` | metres per time step |
| Directional persistence in each state | `turn_concentration_state<k>` | von Mises concentration (0 = no persistence, > 3 = nearly straight) |
| Probability of a zero-length step (only if zeros occur) | `zero_step_prob_state<k>` | probability |
| Activity budget: long-run share of time in each state, averaged over the 24-h cycle at average habitat | `activity_budget_state<k>` | proportion |
| Transition coefficients | `tpm_<i>to<j>_intercept`, `tpm_<i>to<j>_tod_cos`, `tpm_<i>to<j>_tod_sin`, `tpm_<i>to<j>_<covariate>` | log-odds of switching i -> j vs staying in i; covariate slopes **per original unit** (e.g. per 1% forest); intercept at the covariate mean and the daily-mean time of day |

Conditional-effects plots show transition probabilities and stationary
state probabilities (share of time in each state if conditions were held
fixed) against time of day and against each covariate, on the probability
scale with 95% and 50% confidence bands.

## Model

For animal *a* at regular time step *t* (row *t* holds location x_t, the step
from x_t to x_{t+1} and the turning angle at x_t):

    S_t | S_{t-1} = i  ~ Categorical(gamma_i1(t), ..., gamma_iN(t))
    gamma_ij(t) = exp(eta_ij(t)) / sum_k exp(eta_ik(t)),   eta_ii = 0 (staying is the reference)
    eta_ij(t)  = b0_ij + b_cos_ij cos(2 pi h_t / 24) + b_sin_ij sin(2 pi h_t / 24) + sum_l b_lij z_l(t)
    step_t  | S_t = k ~ Gamma(mean mu_k, sd sigma_k)   [ x (1 - pz_k), point mass pz_k at 0 if zero steps occur ]
    angle_t | S_t = k ~ von Mises(0, kappa_k)

* h_t is the local clock hour; z_l(t) are covariates standardised before fitting
  (z = (x - mean)/SD). Slopes are reported per original unit (b / SD), which is
  an exact rescaling.
* The off-diagonal coefficient order is (1,2), (1,3), (2,1), (2,3), (3,1), (3,2)
  as in moveHMM/momentuHMM (Michelot et al. 2016).
* The initial state distribution is estimated (a nuisance parameter).
* Stationary distribution delta(z) solves delta = delta Gamma(z), sum = 1. With
  time of day on the transitions the chain is periodic (period 24 h) and has no
  single stationary distribution. The module uses the **periodically stationary
  distribution** (Koslik et al. 2025): delta_h, the probability of each state at
  a fix with clock hour h, solves delta_h = delta_{h-dt} Gamma(h) around the daily
  cycle. The `stationary_tod` figure shows delta_h; the **activity budget** is
  delta_h averaged over the fixes of one day (the true long-run share of time),
  at mean habitat; habitat curves (`stationary_<cov>`) are the same long-run share
  with the habitat covariate held at each value. (Averaging the hour-specific
  stationary distributions of moveHMM's `plotStationary` instead is not the
  long-run share: it was off by ~2 percentage points for the 3-state design.)
  If 24 h is not a whole number of fix intervals, the hour-specific average is
  used as an approximation.
* Likelihood: forward algorithm (Zucchini et al. 2016), maximised by
  `momentuHMM::fitHMM` (nlm). Missing fixes on the regular grid contribute
  no observation (NA step/angle) but the state process continues.
* **Local optima**: the fit is run from `nstarts` starting values (one
  data-driven: step means at step-length quantiles; the rest randomly
  perturbed on the log scale) and the best likelihood is kept; the report
  says how many starts reached it. `retryFits` additionally uses momentuHMM's
  own perturbation refits.
* **Label switching**: states are relabelled so that state 1 has the shortest
  mean step (refit from the permuted MLE, which reproduces the same optimum).
* **Uncertainty**: 95% confidence intervals for every reported quantity are
  quantiles of `ndraws` parameter vectors drawn from the asymptotic normal
  distribution of the working-scale MLE (log scale for means, SDs and
  concentrations, so intervals stay positive; equal to Wald intervals for
  the coefficients). If the initial-state distribution sits on a boundary
  (all animals start in the same state) the full Hessian is singular; the
  intervals then condition on that nuisance parameter (reported in the notes).
  Likewise, a zero-step probability estimated at ~0 (e.g. a single 0-m step
  switches the zero mass on for every state) is held fixed, with a note, so the
  other intervals remain available (fix from the 2026-10-19 verification;
  before it, ~8% of fits with all animals starting in one state wrongly
  reported "Hessian not invertible").
* **Decoding and checks**: Viterbi most-likely state sequence; pseudo-residuals
  (normal QQ, KS test, autocorrelation) for step lengths and turning angles.

## Assumptions (and what breaks if they fail)

1. **Regular sampling, negligible location error** (Langrock et al. 2012). Fixes are
   snapped to a regular grid (tolerance 25% of the interval); missed fixes are
   treated as missing. `regularise: crawl` (single best crawl prediction) exists
   as a last resort, but on a strongly irregular test set it shortened steps by
   10-18% and doubled turning concentrations (VALIDATION.md section 9); for truly
   irregular data use momentuHMM's `MIfitHMM` multiple imputation outside this
   module. Location error comparable to
   the short steps of an encamped state inflates its step mean and turning.
2. **The number of states is right.** Too few states merge behaviours and bias
   every state's parameters and the activity budget (negative control in
   VALIDATION.md); information criteria tend to pick too many states for real data
   (Pohle et al. 2017). Choose N from biology + pseudo-residuals, not AIC alone.
3. **Conditional independence of steps and angles given the state**, and geometric
   dwell times. Strong residual autocorrelation in pseudo-residuals signals a violation.
4. **All animals share the same parameters** (pooled model). Individual differences
   appear as residual autocorrelation; fit animals separately to compare.
5. **Covariates are measured at the fix and are exogenous.** Habitat effects on
   switching are associations, not causal effects: animals choose where to be.
6. Covariate values at missed fixes are carried forward from the last fix.

## Data format

One CSV, one row per fix (failed fixes can be absent or have empty coordinates):

```
animal_id,timestamp,longitude,latitude,forest_cover_pct
elk_F01,2024-06-01 12:00:33,-123.2671,44.53311,51.4
elk_F01,2024-06-01 12:58:34,-123.2683,44.53602,52.5
elk_F01,2024-06-01 13:59:49,-123.2717,44.54307,45.8
```

* `id`: animal or collar id (text).
* Fix interval: estimated as the most common gap between fixes (snapped to 5, 10, 15, 30 min, 1, 2, 4 h ... when within 10%); the report states the value used. Set `interval_h` if your schedule is unusual.
* `time`: timestamp; set `time_format` (R `strptime` codes) and **`tz` = local
  time zone** so "hour of day" means local clock time (sunrise matters, not UTC).
  If `time_format` does not match, ISO 8601 (`2024-06-01T13:00:00`, optional
  fractional seconds) and `YYYY-MM-DD HH:MM` are tried; a trailing `Z` marks UTC
  and is converted to `tz`. Timestamps that would only be read as dates (clock
  time lost) stop with an error rather than being accepted.
* `lon`/`lat` (WGS84) are projected with sf to the UTM zone of the data (or
  `crs_out`); or give `x`/`y` already in metres.
* Covariate columns: numeric values at each fix (extract rasters beforehand). 0/1
  columns work for categories.

## Config

```yaml
data:
  file: elk_gps.csv
  id: animal_id
  time: timestamp
  time_format: "%Y-%m-%d %H:%M:%S"
  tz: America/Los_Angeles
  lon: longitude          # or x: / y: in metres (+ crs: for documentation)
  lat: latitude
  crs_out: 32610          # optional
  covariates: [forest_cover_pct]
  time_of_day: true       # cos/sin of the hour on transition probabilities
  interval_h: 1           # optional; default = most common gap, snapped to a standard schedule
settings:
  nstates: 2              # 2 or 3 (any N works but is not validated beyond 3)
  nstarts: 10
  retryFits: 0
  ndraws: 2000
  regularise: grid        # or crawl
  tolerance: 0.25
  checks: true            # Viterbi + pseudo-residuals
  state_names: ["Encamped", "Travelling"]
figures:
  labels: {forest_cover_pct: "Forest cover (%)", state1: "Encamped"}
  ids: [elk_F01, elk_F02] # animals in the decoded-track figure (default: first 6)
```

## Outputs

Estimates table (above); diagnostics (optimiser code, Hessian, how many starts
reached the best likelihood, missing/dropped fixes, pseudo-residual KS tests and
lag-1 autocorrelation); plain-language interpretation; figures (docs/FIGURES.md):
`state_distributions`, `transition_tod`, `transition_<cov>`, `stationary_tod`,
`stationary_<cov>`, `decoded_track`, `pseudo_residuals`, `coefficients`.

## Figures

Exported by `core/run.R fit` to `<out>/figures/` (PDF + PNG, `gallery.html`,
`figures_index.csv` with captions); redraw without refitting with
`scripts/Rscript core/run.R figures movement-hmm --fit <out>/movement-hmm_fit.rds --config config.yaml`.
Colour follows the state in every figure (Okabe-Ito order = states ordered by
mean step length); transitions take the colour of the state the animal leaves.
Linetype (curves) and point shape (tracks, coefficients) repeat the state
encoding so `greyscale: true` still separates them. Bands are `ci` (default
95%) and `inner_ci` (50%) confidence intervals from draws of the asymptotic
sampling distribution of the MLE.

| id | shows |
|---|---|
| `state_distributions` | (A) gamma step-length and (B) von Mises turning-angle densities per state, weighted by the activity budget, with bands, the mixture and histograms of the observed steps/turns |
| `transition_tod` | probability of switching between states per time step vs hour of day (local clock, data `tz`), habitat at its mean |
| `transition_<cov>` | the same vs each habitat covariate in original units (rug = values at the fixes), averaged over the daily cycle; with 3+ states one panel per state left |
| `stationary_tod` | probability of being in each state at each clock hour (periodically stationary distribution; points at the fix times) |
| `stationary_<cov>` | long-run share of time in each state vs each covariate in original units |
| `decoded_track` | tracks (km from the first fix) coloured and shaped by the Viterbi state, one panel per animal (`ids:` picks animals; default first 6) |
| `pseudo_residuals` | QQ plots of step and angle pseudo-residuals and the ACF of step residuals |
| `coefficients` | transition coefficients (log-odds of switching vs staying; per SD for covariates) with `ci` and `inner_ci` intervals |

Relabel keys (`figures: labels:`): `state1`, `state2`, ... (state names; else
`settings: state_names`), each covariate name, `hour`, `step`, `angle`,
`density`, `transition_prob`, `stationary_prob`, `state_prob_hour`,
`track_x`, `track_y`, `normal_quantile`, `step_residual`, `angle_residual`,
`acf`, `acf_lag`, `tod_cos`, `tod_sin`, `coef_axis`. Covariates are drawn in the
units of the data file; `units: {<cov>: {center, scale}}` is only needed if the
file itself holds a standardised covariate.

## Simulation and power

`simulate()` generates correlated-random-walk tracks on a smooth 20 x 20 km
habitat surface (% forest), with state switching driven by time of day and the
forest value at the animal's location, 5% failed fixes and +/-2 min fix-time jitter.

* params: `step_mean`, `step_sd`, `angle_kappa` (length N), `tpm_intercept`,
  `tpm_tod_cos`, `tpm_tod_sin`, `tpm_forest` (length N(N-1), order above; forest
  slope per %; intercept at `forest_ref` = 50%), optional `zero_mass`.
* design: `n_animals` (6), `n_fixes` per animal (600), `fix_interval_h` (1),
  `missing_prop` (0.05), `jitter_min` (2), `landscape_km` (20).
* Power knobs: `n_fixes`, `n_animals`, effect size `tpm_forest`, state separation.
  At the default effect (-0.025 / +0.025 per % forest), see VALIDATION.md for the
  power curve over track length.

## Validation (summary; details in VALIDATION.md)

* Default 2-state design (6 animals x 600 hourly fixes): all 16 estimands pass the recovery gate
  over 200 replicates (|relative bias| <= 2.8%, coverage 0.92-0.99).
* 3-state model: 36/36 pass after the boundary-Hessian fix (rare-transition coefficient at coverage
  0.889, just above the bound 0.884); rare-transition coefficients remain poorly identified.
* Independent simulator (other author; 2-h fixes, local time with DST, elevation covariate):
  15/16 pass (one time-of-day coefficient at coverage 0.889 vs bound 0.896); with zero-length
  steps all 18 pass, including the zero-step probabilities.
* Harder designs (short tracks, weakly separated states): 15/16 pass each.
* Negative control (2-state model on 3-state data) fails 12/16, as it must.
* Second independent simulator (30-min fixes, two habitat covariates, lon/lat CSV read through
  read_data): all 18 pass, 120/120 converged; null version all 18 pass.
* False-positive rate for a null habitat effect: 0.07-0.08 (nominal 0.05); 0.069 pooled over four
  null studies (1596 tests). Treat effects whose interval only just excludes zero as weak evidence.
* Agrees with moveHMM and a hand-written likelihood within 0.014 SE on every parameter.
* Power: about 400 hourly fixes per animal (6 animals) for 80% power at the default habitat effect.

## Pitfalls

* Calling states "foraging" or "resting" is an interpretation: an HMM sees only
  step lengths and angles. Validate with observations or accelerometers.
* An HMM always finds N states, even in data without behavioural structure.
* Steps of exactly 0 m (stationary collars, rounding) trigger a zero-mass
  component automatically; many zeros from GPS rounding may deserve cleaning first.
* Different fix schedules between animals (1 h vs 4 h) must not be pooled:
  step-length distributions depend on the interval.
* Very short encamped steps are dominated by GPS error (~5-30 m).
* Transition coefficients are hard to read directly; use the probability-scale
  plots and the activity budget.

## References

* Langrock R, King R, Matthiopoulos J, Thomas L, Fortin D, Morales JM (2012) Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions. Ecology 93(11):2336-2342. doi:10.1890/11-2241.1
* McClintock BT, Michelot T (2018) momentuHMM: R package for generalized hidden Markov models of animal movement. Methods in Ecology and Evolution 9(6):1518-1530. doi:10.1111/2041-210X.12995
* Michelot T, Langrock R, Patterson TA (2016) moveHMM: an R package for the statistical modelling of animal movement data using hidden Markov models. Methods in Ecology and Evolution 7(11):1308-1315. doi:10.1111/2041-210X.12578
* McClintock BT, Langrock R, Gimenez O, Cam E, Borchers DL, Glennie R, Patterson TA (2020) Uncovering ecological state dynamics with hidden Markov models. Ecology Letters 23(12):1878-1903. doi:10.1111/ele.13610
* Pohle J, Langrock R, van Beest FM, Schmidt NM (2017) Selecting the number of states in hidden Markov models: pragmatic solutions illustrated using animal movement. Journal of Agricultural, Biological and Environmental Statistics 22(3):270-293. doi:10.1007/s13253-017-0283-8
* Johnson DS, London JM, Lea MA, Durban JW (2008) Continuous-time correlated random walk model for animal telemetry data. Ecology 89(5):1208-1215. doi:10.1890/07-1032.1
* Koslik J-O, Feldmann CC, Mews S, Michels R, Langrock R (2025) Inference on the state process of periodically inhomogeneous hidden Markov models for animal behavior. Annals of Applied Statistics 19(4):2724-2737. arXiv:2312.14583
* Zucchini W, MacDonald IL, Langrock R (2016) Hidden Markov Models for Time Series: An Introduction Using R, 2nd edn. Chapman & Hall/CRC.
