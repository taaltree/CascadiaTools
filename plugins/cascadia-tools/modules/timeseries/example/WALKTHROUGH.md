# Worked example: is a mule deer herd declining, and do deep-snow winters matter?

Module: `timeseries` (population state-space models). This walkthrough shows
what a researcher who does not code gives the agent, what the agent does at
each step of `skills/eco-analyze`, and every product. All numbers below come
from real runs of the commands in section (g).

**The data are simulated** (see `README_data.md`). We therefore know the true
answer and can check the analysis against it in section (e).

---

## (a) The question, in the researcher's words

> "We have flown a winter count of our mule deer herd every year since 1994,
> except three winters we couldn't fly. The counts bounce around a lot. Is the
> herd going down? How many deer are really out there? Is the herd limiting
> itself? And do the deep-snow winters knock it back?"

## (b) The files the researcher has

One CSV, `herd_counts.csv` (31 rows, one per winter). The first 5 lines:

```
"year","survey","animals_counted","snow_depth_cm","notes"
1994,"aerial winter count",590,83,""
1995,"aerial winter count",625,48,""
1996,"aerial winter count",629,96,""
1997,"aerial winter count",421,33,""
```

A winter without a flight looks like this (count left blank):

```
2001,"aerial winter count",,43,"not flown (weather)"
```

## (c) What the researcher types to Claude

> I've attached herd_counts.csv: aerial winter counts of our mule deer herd,
> 1994-2024, with mean winter snow depth from the snow course. Three winters
> weren't flown. Please tell me whether the herd is declining and by how much
> per year, how many deer there really are, whether growth is density
> dependent, and whether deep-snow winters reduce herd growth. I don't code,
> so please give me a report with figures I can use in a paper.

## (d) What the agent does, step by step

### Step 1-2. Project folder and data intake (eco-data-intake)

The agent makes a project folder, copies the CSV in unchanged, and describes
the data back:

- 31 winters (1994-2024). 28 have a count. 3 interior winters (2001, 2009,
  2016) were not flown.
- Counts run from 316 to 1,125 deer (median 532). No zeros, so log-scale
  models are fine.
- Snow depth is present in every year, including the 3 unflown winters
  (5 to 96 cm; mean 54.1 cm, SD 20.0 cm). This matters, because the model
  needs the covariate in every year to carry the herd through the gaps.
- One count per winter, one herd, one survey method throughout. The agent
  asks: "Did the aircraft, observers or survey area change?" (a change in how
  many deer are seen would look like a population trend). Here: no.

### Step 3. Choose the module

`catalog.md` points to `timeseries`: one count per year, 10+ years, questions
about trend, true abundance, density dependence and a yearly driver. Not
`nmixture` (no repeat counts within a winter) and not `pva` (no age data).
Within the module the agent picks `model: bayes`, the only option that
puts a covariate on growth and handles the missing winters as unknown years.

### Step 4. Analysis plan (eco-analysis-plan)

Written before fitting:

- **Estimands:** realised trend (% per year, 1994-2024); true herd size in
  2024 (`N_final`); equilibrium herd size; effect of 1 SD (20 cm) more snow on
  the annual growth multiplier; density dependence `b`; process and
  observation SDs.
- **Model:** Bayesian Gompertz state-space model. Log herd size changes each
  year by a growth rate that depends on herd size (density dependence) and on
  that winter's snow depth, plus random year effects. The count equals the
  true herd size plus lognormal counting error.
- **Priors:** module defaults (r0 ~ N(0,1); b ~ U(-2,0); beta ~ N(0,1) per SD;
  both SDs half-normal(0,1)). No independent estimate of counting error exists.
- **Adequate fit:** R-hat < 1.01, bulk ESS >= 400 for every reported number;
  posterior predictive p-values between 0.05 and 0.95.
- **Declared in advance:** `VALIDATION.md` says density dependence and the two
  SDs are not reliably recovered from one count series. The plan states that
  they will be reported only if the planted-truth check passes for them.

### Step 5. Planted truth at this design (eco-planted-truth)

The agent simulates 100 herds with exactly this design: 31 winters, 3 interior
winters not flown, one snow covariate, and realistic values (b = -0.3,
equilibrium 500 deer, process SD 0.10, counting error SD 0.15, snow effect
-0.15 per SD). It fits each one and checks whether the truth is recovered.
Config: `example/recovery.yaml`. Each fit takes about 74 s; 93 of 100
replicates met the convergence rule (the gate needs 90%).

| quantity | truth | mean estimate | relative bias | 95% coverage | gate |
|---|---|---|---|---|---|
| Snow effect on growth (log scale, per SD) | -0.15 | -0.156 | -4.2% | 0.97 | **PASS** |
| Snow effect (% change in growth multiplier per SD) | -13.9% | -14.4% | -3.5% | 0.97 | **PASS** |
| True herd size in final year | 517 (mean over replicates) | 527 | +1.6% | 0.96 | **PASS** |
| Equilibrium herd size | 500 | 504 | +0.8% | 0.96 | **PASS** |
| Realised trend (% per year) | +0.12 | +0.12 | bias +0.004 points (MC SE 0.067) | 0.94 | **PASS** |
| Density dependence b | -0.30 | -0.48 | **-59%** | 0.89 | **FAIL** (bias) |
| Process SD | 0.10 | 0.123 | **+23%** | 0.98 | **FAIL** (bias) |
| Observation SD | 0.15 | 0.130 | **-13%** | 0.99 | **FAIL** (bias) |

The snow effect was detected (95% interval excluding 0) in 96.8% of
replicates. The realised trend's relative bias looks huge (-156%) only
because its true value is close to 0; its absolute bias is 0.004 percentage
points, well within Monte Carlo error.

**Gate result, in one sentence:** at this design the model recovers the
trend, the true herd size, the equilibrium herd size and the snow effect with
honest intervals, but it overstates density dependence by about 59% and
misallocates noise between real change and counting error, so `b` and the two
SDs must not be interpreted.

(The recovery used the module's fast MCMC settings: 3 chains x 12,000
iterations. The real fit below uses 4 chains x 40,000.)

### Step 6. Fit

`example/config.yaml` points at the CSV, names the columns, sets
`model: bayes` with density dependence, and gives the figure labels. The
agent runs the fit as a cluster job (section g). Run time: about 2 minutes.

### Step 7. Convergence and model checking (eco-convergence, eco-model-checking)

- **Convergence:** all reported quantities have R-hat <= 1.006 (max) and bulk
  ESS >= 684 (min; process SD). Chains overlap in the trace plot. PASS.
- **Goodness of fit:** posterior predictive check on year-to-year changes in
  log counts. Variability p = 0.61, lag-1 autocorrelation p = 0.75, both
  inside 0.05-0.95. No sign of misfit. PASS.

### Step 8. Interpretation (from the report)

Quoted from `results/timeseries/si/fit/timeseries_report.html`:

> "The estimated true population changed by +0.1% per year on average between
> 1994 and 2024 (95% credible interval -1.4% to +1.3%). This realised trend
> describes what happened; with density dependence it is not a constant
> growth rate to project forward."

> "Estimated true abundance (on the scale of the animals counted) in 2024: 614
> (95% credible interval 499 to 768). If counts miss a constant fraction of
> animals this is an index of abundance, not a census; trends are unaffected
> by a constant detection rate."

> "An increase of one standard deviation (20) in Mean winter snow depth (cm)
> is associated with a -15.3% change in the annual growth multiplier (95%
> credible interval -22.2% to -7.0%; P(effect > 0) = 0.00). At average
> abundance the population is expected to change by +19.9% in a year with
> Mean winter snow depth (cm) at its 10th percentile (30.9) vs -16.4% at its
> 90th percentile (74.5). This is an association, not proof that Mean winter
> snow depth (cm) causes the change."

> "Because the prior allows only b between -2 and 0, the credible interval for
> b can never include 0, so it is not evidence that density dependence
> exists."

What the agent tells the researcher, in plain words:

- **No clear decline.** Over 30 years the herd changed by +0.1% per year
  (95% CI -1.4% to +1.3%). Compounded, that range runs from about a 35% loss
  to a 48% gain over the period, so "no clear trend" is not proof of
  stability.
- **About 614 deer in 2024** (95% CI 499 to 768), on the count scale. The
  2022 count of 1,125 was mostly counting error: the model puts the true
  herd at about 900 that winter.
- **The herd fluctuates around about 536 deer** (95% CI 478 to 599) at
  average snow.
- **Deep snow is associated with poorer growth.** Each extra 20 cm of mean
  winter snow goes with a 15% lower growth multiplier (95% CI -22% to -7%).
- **Density dependence cannot be judged from these counts.** The estimate
  (b = -0.59) is known to be inflated at this design (planted truth: -59%
  bias). The agent does not report it as evidence of regulation. As a
  second check it runs the valid test (`model: gompertz`, parametric-bootstrap
  likelihood ratio, 199 simulations): p = 0.585. The counts are consistent
  with no density dependence, but this test has very little power (about
  2-6% at b = -0.3; VALIDATION.md), so this is not evidence that the herd is
  unregulated either.

### Figures (in `results/timeseries/si/fit/figures/`, PNG and PDF)

| id | caption |
|---|---|
| `trajectory` | Estimated true herd size (line) with 95% (light) and 50% (dark) credible bands; points are the raw counts. The gap between points and line is counting error. |
| `trend` | Realised trend over 1994-2024 in % per year with 95% credible interval; dashed line = no change. |
| `density_dependence` | Expected annual growth rate vs herd size (line, 95% band) at average snow; points are naive growth from consecutive raw counts. Shown for completeness; b is not recovered at this design. |
| `growth_snow_depth_cm` | Expected growth multiplier (1 = stable) vs mean winter snow depth at average herd size; 95% and 50% bands; rug = observed winters. |
| `coefficients` | Snow effect per 1 SD (20 cm) on the growth multiplier, % with 95% CI. |
| `variances` | Process SD (real change) and observation SD (counting error) with 95% CIs. Not reliably separated at this design. |
| `forecast` | Estimated herd size and 10-year forecast (dashed) at average snow, with 95% predictive band. |
| `power_years` | Power to detect a 5%/yr decline vs years of monitoring (single series, process SD 0.1, counting SD 0.15; 200 simulations per point). |
| `ppc` | Posterior predictive check: (A) variability and (B) lag-1 autocorrelation of year-to-year changes in log counts; line = observed (p = 0.61, 0.75). |
| `trace` | MCMC traces for r0, b, sigma_p, sigma_o and the snow coefficient, 4 chains. |
| `posterior` | Posterior densities of the reported quantities. |

The single most informative figure is **`trajectory`**: it shows the answer to
the trend question, the 2024 herd size, and how much of the bounce in the raw
counts is counting error.

### Step 9. Acceptance gate (eco-acceptance-gate)

| gate step | result |
|---|---|
| 0 Declared before fitting | Yes (plan above; DD and SDs flagged in advance) |
| 1 Planted truth at this design | PASS for trend, N_final, equilibrium_N, snow effect. FAIL (bias) for b, process SD, observation SD |
| 2 Convergence | PASS (max R-hat 1.006, min ESS 684) |
| 3 Model checking | PASS (PPC p = 0.61, 0.75) |
| 4 Independent route | Module cross-checks (VALIDATION.md section 5): Nimble vs Kalman ML on 60-year series agree on N_final (574 vs 567) but not on b (-0.735 vs -0.230), the same non-identifiability |
| 5 Negative control | Ignoring counting error (`observation_error: false`) fails recovery as it should: b -162%, N_final coverage 0.000 (VALIDATION.md section 2) |
| 6 Reproducible | Commands, configs, seed below |

**Verdict: OPEN (estimand-specific).** The trend, the 2024 herd size, the
equilibrium herd size and the snow effect pass every step and can be
reported. The density-dependence question is not answered: b and the split
between real change and counting error fail planted truth at this design.
To answer it, the herd needs an independent estimate of counting error
(replicate flights or sightability trials, supplied as
`prior_sd_observation`), or the Gompertz bootstrap test
(`model: gompertz`), which keeps its false-positive rate but has little power.

### Step 10. Methods paragraph (from the report)

> We fitted a Bayesian state-space population model (Kery & Schaub 2012) in
> NIMBLE (de Valpine et al. 2017): x_t = x_{t-1} + r0 + b (x_{t-1} - mean log
> count) + beta z_t + e_t, e_t ~ N(0, sigma_p^2), log count_t ~ N(x_t,
> sigma_o^2) for surveyed years, where z_t is standardised snow depth acting
> on growth into year t. Missing years were handled as unobserved latent
> states. Priors: r0 ~ N(0, 1), b ~ U(-2, 0), beta ~ N(0, 1), sigma_p ~
> half-N(0, 1), sigma_o ~ half-N(0, 1). MCMC: 4 chains of 40000 iterations
> (10000 burn-in, thinning 10).

## (e) Estimates versus the known truth

| quantity | estimate (95% CI) | truth | inside CI? |
|---|---|---|---|
| Realised trend, % per year | +0.07 (-1.44 to +1.31) | +0.10 | yes |
| True herd size 2024 | 614 (499 to 768) | 643 | yes |
| Equilibrium herd size | 536 (478 to 599) | 500 | yes |
| Snow effect, log scale per SD | -0.165 (-0.251 to -0.073) | -0.150 | yes |
| Snow effect, % per SD | -15.3% (-22.2% to -7.0%) | -13.9% | yes |
| Density dependence b | -0.59 (-1.11 to -0.29) | -0.30 | yes, barely (point estimate ~2x too strong) |
| Process SD | 0.097 (0.004 to 0.265) | 0.10 | yes (very wide) |
| Observation SD | 0.184 (0.033 to 0.278) | 0.15 | yes (very wide) |

Every interval contains the truth. The point estimate of b is about twice the
true strength, exactly the bias the planted-truth check predicted.

## (f) Limits from VALIDATION.md that matter here

1. **Density dependence is not recoverable from one count series** (Known
   limit 2). In `bayes` the U(-2, 0) prior means the interval for b can never
   include 0. Here b was biased by -59% in planted truth.
2. **Process and observation SDs are biased** (Known limit 1): +23% and -13%
   here. Trend and abundance intervals stay valid.
3. **Constant detection is assumed** (Known limit 7). Herd size is on the
   count scale (an index), and a change in sightability would look like a
   trend.
4. **"No clear trend" is weak evidence.** A single series declining 5%/yr is
   detected with power 0.29 after 20 years and 0.49 after 30 (figure
   `power_years`).
5. **Bayes DD fits can mix slowly** (Known limit 9). Here min ESS was 684,
   above the 400 threshold; in planted truth 7% of fast-setting fits missed
   the convergence rule.
6. **Covariate effects are associations**, not proof that snow causes poor
   growth.

## (g) Reproduce

From `.`:

```bash
# 0. (optional) regenerate the example data and check it matches the CSV exactly
scripts/ejob -c 1 -m 4 -t 0:15:00 -l logs/timeseries_si_repro.log -- scripts/Rscript dev/timeseries/si/check_example_repro.R

# 1. planted truth at the example's design (100 replicates, ~10 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/timeseries_si_recover.log -- scripts/Rscript core/run.R recover timeseries --config modules/timeseries/example/recovery.yaml --nsim 100 --cores 16 --out results/timeseries/si/recovery
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/timeseries_si_recovery_extra.log -- scripts/Rscript dev/timeseries/si/recovery_extra.R

# 2. fit the real (example) data: estimates, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/timeseries_si_fit.log -- scripts/Rscript core/run.R fit timeseries --config modules/timeseries/example/config.yaml --out results/timeseries/si/fit

# 2b. valid test for density dependence (Gompertz bootstrap LRT, ~2 min)
scripts/ejob -c 2 -m 8 -t 0:30:00 -l logs/timeseries_si_fit_gompertz.log -- scripts/Rscript core/run.R fit timeseries --config modules/timeseries/example/config_gompertz.yaml --out results/timeseries/si/fit_gompertz

# 3. (optional) redraw figures after editing the config's figures: block
scripts/ejob -c 1 -m 8 -t 0:20:00 -l logs/timeseries_si_figures.log -- scripts/Rscript core/run.R figures timeseries --config modules/timeseries/example/config.yaml --fit results/timeseries/si/fit/timeseries_fit.rds
```

Outputs: `results/timeseries/si/recovery/` (gate table, report),
`results/timeseries/si/fit/` (estimates CSV, HTML report, figures),
`results/timeseries/si/si_summary.json` (machine-readable summary).
