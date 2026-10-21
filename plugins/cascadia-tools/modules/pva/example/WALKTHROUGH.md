# Worked example: is a shorebird population at risk? (pva module)

This walkthrough shows what happens when a researcher who does not code asks
CascadiaTools for a population viability analysis. The data are **simulated**,
so we can compare every estimate with the truth. Background and true values:
`README_data.md`.

## (a) The question, in the researcher's words

> "We've counted breeding pairs on our beaches every spring since 1994, but we
> missed a few years. Counts have dropped a lot. Is the population really
> declining, or is it just noisy counts? And how likely is it to fall to 40
> pairs in the next 10 or 25 years?"

## (b) The files the researcher has

One file, `shorebird_breeding_pairs.csv` (30 rows, one per year):

```
"survey_year","breeding_pairs","observer_notes"
1994,412,""
1995,350,""
1996,430,""
1997,,"no survey"
```

A blank count means the beaches were not surveyed that year.

## (c) The prompt they type

> I have 30 years of spring counts of breeding pairs for a shorebird
> (shorebird_breeding_pairs.csv; blank = no survey). The counts are surveys, so
> they're not exact. Please tell me whether the population is declining and how
> likely it is to drop to 40 breeding pairs or fewer within 10 and 25 years.
> I don't code, so please do the analysis and give me a report and figures I
> can use.

## (d) What the agent does, step by step (eco-analyze)

**1. Project folder.** The agent copies the CSV into `projects/<name>/data_raw/`
and never edits it. (For this example, outputs go to `results/pva/si/`.)

**2. Data intake.** The agent describes the data back:
- 30 years (1994-2023), 25 surveyed, 5 blank (1997, 2003, 2004, 2011, 2020).
- Counts range from 79 to 528 pairs; no zeros (needed for the log scale).
- First counts about 410; last three counts 97, 83, 79.
- One large drop between 2008 (401) and 2009 (155).

**3. Choose the module.** The question is extinction risk from one count per
year, so `pva` with `model: count`. The counts are surveys with error, and some
years are missing, so the agent picks `method: state_space` instead of the
classic diffusion model. (Diffusion treats counts as exact; counting error then
inflates the variance and the risk. The module's negative control shows +209%
bias in the variance.)

**4. Plan.** Estimands: long-run growth rate (lambda), year-to-year
variability, counting error, the true current number of pairs, and the
probability of reaching 40 or fewer pairs at an annual census within 10 and
25 years. No covariates. Priors: flat on mu and starting size, uniform on the
two standard deviations. The plan warns in advance that 10- and 25-year risks
are beyond 20% of the 29-year data span, so their intervals will be wide.

**5. Planted truth (eco-planted-truth).** Before touching the real counts, the
agent simulates 200 data sets at this exact design (30 years, about 18% of
years missing, 420 pairs at the start, threshold 40) and fits each one.
Result (`results/pva/si/recovery/`, nsim = 200, 0 failures):

| quantity | truth* | mean estimate | relative bias | 95% coverage | gate |
|---|---|---|---|---|---|
| growth rate mu | -0.040 | -0.042 | -3.8% | 0.970 | PASS |
| lambda | 0.961 | 0.960 | -0.1% | 0.970 | PASS |
| process variance sigma2 | 0.012 | 0.014 | +16.5% | 0.985 | PASS (within Monte Carlo error) |
| counting-error variance tau2 | 0.030 | 0.029 | -4.1% | 0.975 | PASS |
| true current pairs | 149 (mean) | 150 | +0.5% | 0.940 | PASS |
| risk within 10 yr | 0.162 (mean) | 0.238 | bias +0.076 | 0.985 | **FAIL (bias)** |
| risk within 25 yr | 0.463 (mean) | 0.492 | bias +0.029 | 0.985 | PASS |

*Current abundance and risk differ in each simulated data set, so the mean
truth is shown.

In one sentence: the model recovers growth, variability, counting error and
current size at this design, and the risk intervals contain the true risk 98.5%
of the time, but the 10-year risk **point estimate** is too extreme on average
(Known limit 1). So the report must quote the risk interval, not the point.

Power: with this much noise, a 30-year series detects a true 3.9%/yr decline
(95% interval for mu excludes 0) in only 22.5% of data sets.

**6. Fit.** The agent writes `config_state_space.yaml` (file, column names,
`method: state_space`, threshold 40, horizons 10 and 25 years) and runs the
fit. It takes about a second: the posterior is computed exactly on a grid,
with 4000 draws and no MCMC.

**7. Checks (eco-convergence, eco-model-checking).**
- Convergence: not an issue. No MCMC; the grid posterior is exact (R-hat and
  ESS do not apply). The engine was cross-checked against a Nimble MCMC of the
  same model on this very data set (posterior quantiles agree within a median
  0.03 posterior SD).
- Autocorrelation of annual log growth: lag-1 r = -0.31. Negative values are
  expected when counts have error, so no sign of runs of bad years.
- Density dependence: slope of growth on log abundance -0.012 (p = 0.87). No
  sign of it.
- Normality: the log growth increments of the raw counts fail Shapiro-Wilk
  (p = 0.003). This is driven by one value, 2008 to 2009 (401 to 155 pairs).
  The agent flags it and would ask the researcher whether something happened
  that year (storm, method change). Here the data are simulated with normal
  errors, so it is a false alarm: a high 2008 count (counting error) followed by
  a low 2009 count. Raw-count increments also mix in counting error, so the
  check is conservative.

**8. Interpretation.** Excerpts from the report (`results/pva/si/fit/pva_report.html`):

> "the population was not clearly increasing or decreasing: the long-run
> (median) growth rate is lambda = 0.942 per year (95% CI 0.879 to 1.013), i.e.
> a typical change of -5.8% per year (95% CI -12.1% to +1.3%)."

> "Counting (observation) error variance on the log scale is 0.0256 (95% CI
> 0.0003 to 0.0727), equivalent to a coefficient of variation of roughly 16% in
> a single count. The model's estimate of the true current abundance is 79
> breeding pairs (95% CI 61 to 100)."

> "Probability that the population falls to or below the quasi-extinction
> threshold of 40 breeding pairs in at least one annual census within 10 years:
> 49% (95% CI 6% to 92%); 25 years: 90% (95% CI 31% to >99.9%)."

The agent adds context for the researcher: 96% of posterior draws have lambda
below 1, so a decline is very likely, even though the 95% interval just
reaches 1.0. "Not clearly" means "not certain", not "stable".

Figures (`results/pva/si/fit/figures/`):

| figure | caption |
|---|---|
| `quasi_extinction.png` (key figure) | Cumulative probability that the number of breeding pairs falls to 40 or fewer by each year from now. Line: estimate; dark band: 50% and light band: 95% credible interval. |
| `trajectories.png` | Counts (points), state-space estimate of the true number of pairs with 95% band (grey), and projected numbers to 2048: median line, 50% and 95% bands. Dashed red line: threshold of 40 pairs. |
| `lambda_distribution.png` | Posterior distribution of the long-run growth rate. Solid line: estimate (0.942); dashed line: lambda = 1 (stable). 96% of draws are below 1. |
| `growth_checks.png` | Model checks. A: annual log growth rate by year with the estimated mean; B: growth rate against log abundance (no density dependence: slope -0.012, p = 0.87). The 2009 outlier is the 401-to-155 drop. |

**9. Acceptance gate: OPEN.** Growth rate, variances and current abundance are
validated at this design and the fit is healthy. But the main estimand includes
quasi-extinction risk, and its point estimate fails the planted-truth bias gate
at this design (10-year horizon). The risk intervals are calibrated (coverage
0.985), so the result is usable as an interval. An honest OPEN: "use the
intervals, do not quote 49% as a forecast".

**10. Delivery (what the researcher receives).** The counts point to a decline
of about 6% per year (lambda 0.94, 95% CI 0.88 to 1.01; 96% probability that
the population is shrinking). Allowing for counting error, there are probably
about 79 breeding pairs now (61 to 100). If the past 30 years continue, the
chance of reaching 40 pairs or fewer within 10 years is very uncertain (6% to
92%), and within 25 years it is likely (31% to >99.9%). These forecasts assume
no change in conditions or management; use them to compare options, not to
predict a date.

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% interval) | truth inside? |
|---|---|---|---|
| growth rate mu | -0.040 | -0.059 (-0.129 to 0.013) | yes |
| lambda | 0.961 | 0.942 (0.879 to 1.013) | yes |
| process variance sigma2 | 0.012 | 0.019 (0.003 to 0.097) | yes |
| counting-error variance tau2 | 0.030 | 0.026 (0.0003 to 0.073) | yes |
| true pairs in 2023 | 66.0 | 78.9 (61.5 to 100.0) | yes |
| risk within 10 yr | 0.446 | 0.495 (0.063 to 0.924) | yes |
| risk within 25 yr | 0.882 | 0.904 (0.310 to 1.000) | yes |

All seven intervals contain the truth. The decline is slightly overestimated
and the current size slightly overestimated (the 2023 count of 79 was above the
true 66), both well within the intervals.

## (f) Limits from VALIDATION.md that matter here

1. Quasi-extinction point estimates are median-unbiased but too extreme on
   average; they fail the mean-bias gate (here at 10 years). Read the interval.
2. Horizons beyond about 20% of the data span (here about 6 years) give very
   wide risk intervals (Fieberg & Ellner 2000). Both 10 and 25 years are beyond it.
3. With this much noise, power to detect a real decline of this size is low
   (22.5% in our simulations), so "interval includes 1" is not evidence of stability.
4. The counting-error variance can be biased low (-18%) when counting error is
   larger than process noise in short series with gaps; here (30 years) the bias
   was -4%.
5. Not tested by recovery because the simulator obeys the model: density
   dependence, autocorrelated environments, catastrophes, trends in conditions.

## (g) Commands to reproduce

Run from `.` (R runs through `scripts/ejob`):

```
# example data (already in this folder) and the true values
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/pva_si_truth.log -- scripts/Rscript dev/pva/si_truth.R
# planted-truth recovery at this design
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/pva_si_recover.log -- scripts/Rscript core/run.R recover pva --config modules/pva/example/recovery.yaml --nsim 200 --cores 16 --out results/pva/si/recovery
# fit, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/pva_si_fit.log -- scripts/Rscript core/run.R fit pva --config modules/pva/example/config_state_space.yaml --out results/pva/si/fit
```

(`dev/pva/make_examples.R` wrote the CSV originally; `si_truth.R` checks that
the module's own simulator reproduces it exactly.)
