# Worked example: how fast are adult walleye dying? (catch curve)

Module: `fisheries-stock` (`model: catch_curve`). The data are **simulated**
(see `README_data.md`), so we can compare the answer with the truth.
All numbers below come from runs made on 2026-10-21, after the catch-curve
estimators were revised (VALIDATION.md section 5). The previous version of this
walkthrough (weighted regression primary, verdict OPEN) is kept in
`dev/fisheries-stock/backup_pre_cc2/WALKTHROUGH.md`.

## (a) The question, in the researcher's words

> "We aged every walleye from our 2025 fall gill-net survey. What is the total
> annual mortality of adult walleye in the lake? What fraction survive each year?"

## (b) The files the researcher has

One file, `walleye_fall_gillnet_ages.csv`, with one row per fish (350 rows):

```
"fish_id","sample_date","basin","gear","total_length_mm","age"
"WAE25-0001",2025-09-22,"South Basin","experimental gill net",426,2
"WAE25-0002",2025-09-22,"South Basin","experimental gill net",546,6
"WAE25-0003",2025-09-22,"North Basin","experimental gill net",431,3
"WAE25-0004",2025-09-22,"North Basin","experimental gill net",396,2
```

## (c) What they type to Claude

> "Here are otolith ages from our fall walleye gill-net survey
> (walleye_fall_gillnet_ages.csv). I want to know the total annual mortality and
> survival of adult walleye. Please do the analysis and give me a report with
> figures I can use."

## (d) What the agent does (eco-analyze steps)

**1-2. Project folder and data intake.** The agent reads the file and describes it
back. There are 350 fish from two basins, 22 Sep to 10 Oct 2025. Ages range from 1 to 13.
Two fish have no age (unreadable otoliths), so they are dropped, which leaves 348.
Catch at age is 11, 66, 129, 60, 24, 31, 4, 11, 4, 2, 5, 0, 1 for ages 1-13.
Catch peaks at age 3, so younger fish are not yet fully caught by the nets.
The counts jump around (31 at age 6 but 4 at age 7): strong and weak year classes,
as in every real fish population.

**3. Choose the module.** Ages of individual fish, question = total mortality:
`fisheries-stock` with `model: catch_curve`. Length is not needed. There is one
year of ages, so the single-year estimators apply (with several years of ages the
agent would use the multi-year cohort catch curve, see section (h)).

**4. Plan.**
- Estimand: total instantaneous mortality Z for fully recruited fish. Annual survival is S = exp(-Z) and annual mortality is A = 1 - S.
- Primary estimator: the module default, a Bayesian Poisson-lognormal catch curve (Nimble) with a random effect for each year class, from age 4 (peak age 3 + 1) on, plus 5 zero-catch ages beyond the oldest fish. Priors: Z ~ Uniform(0, 5), year-class SD ~ half-Normal(0, 1), intercept ~ Normal(0, 10).
- Cross-checks: Poisson-lognormal GLMM (glmmTMB, REML + t interval), Chapman-Robson, weighted regression.
- Assumptions to tell the researcher: year-class strength varies at random with no trend, constant Z over ages 4+, full selectivity from age 4, accurate ages.

**5. Planted truth** (`recovery.yaml`, 200 simulated surveys at this design:
348 fish, Z = 0.55, year-class SD 0.4, the same gill-net selectivity; gate bound for
coverage 0.904).

| parameter | truth | mean estimate | rel. bias | 95% coverage | gate |
|---|---|---|---|---|---|
| **Z (Bayesian, primary)** | 0.550 | 0.556 | +1.0% | 0.960 | **PASS** |
| **annual survival** | 0.577 | 0.576 | -0.2% | 0.960 | **PASS** |
| **annual mortality** | 0.423 | 0.424 | +0.3% | 0.960 | **PASS** |
| Z (GLMM, REML + t; cross-check) | 0.550 | 0.537 | -2.3% | 0.930 | PASS |
| Z (Chapman-Robson; cross-check) | 0.550 | 0.525 | -4.6% | 0.850 | FAIL (coverage) |
| Z (weighted regression; cross-check) | 0.550 | 0.547 | -0.5% | 0.900 | FAIL (coverage) |

All 200 fits converged (R-hat < 1.05). The Bayesian interval is honest: the
average reported posterior SD of Z was 0.084 and the estimates actually spread
with SD 0.083. For comparison, the weighted regression reported SD 0.070 against a
real spread of 0.090, which is why its interval covered only 90%.
In one sentence to the researcher: *at your design the default method finds the
right answer on average (bias +1%) and its 95% interval contains the truth 96% of
the time.*

**6. Fit** (`config.yaml`; about 2 minutes, most of it Nimble compiling).

**7. Checks.**
- Convergence (Bayesian fit, 3 chains x 40,000 iterations, 10,000 burn-in, thin 10): R-hat 1.001, minimum bulk ESS 2,901 for Z and the year-class SD; the `trace.png` chains overlap with no trend. Passes R-hat < 1.01 and ESS >= 400.
- Goodness of fit: the age distribution is far more variable than a smooth geometric decline (Chapman-Robson c-hat = 2.97, Pearson chi-square 23.8 on 8 df, P = 0.002). This is exactly the year-class variability the Bayesian model includes: its year-class SD is 0.54 (95% interval 0.16-1.24) on the log scale.
- Independent route: the GLMM fitted by likelihood gives Z = 0.50 (0.35-0.66), within 0.3 posterior SD of the Bayesian 0.53.
- The oldest ages got weight 0 in the weighted regression (predicted log catch <= 0); normal under the Maceina & Bettoli weighting and irrelevant to the primary estimate.

**8. Interpretation** (quoted from the report):

> "Total instantaneous mortality of fully recruited walleye (ages 4 and older) is
> Z = 0.53 per year (95% credible interval 0.39 to 0.74; Bayesian Poisson-lognormal
> catch curve). Equivalently, about 59% (48% to 68%) of fully recruited fish
> survive each year and 41% (32% to 52%) die from all causes combined (fishing
> plus natural mortality)."

> "For comparison, the other estimators give Z = 0.50 (0.33 to 0.66;
> Chapman-Robson); Z = 0.53 (0.33 to 0.72; weighted regression); Z = 0.50 (0.35 to
> 0.66; Poisson-lognormal GLMM). The intervals of Chapman-Robson, weighted
> regression and Poisson-lognormal GLMM are known to be too narrow when year
> classes vary (VALIDATION.md), so they are shown only as a cross-check of the
> estimate. All estimators agree within their uncertainty."

> "Real fish populations have strong and weak year classes, so the counts at age
> scatter around the straight line. Here the scatter is much larger than what
> constant recruitment would give (Chapman-Robson c-hat = 2.97); the estimated
> year-class SD is 0.54 on the log scale, i.e. a typical year class is about 72%
> stronger or weaker than average. The Bayesian Poisson-lognormal catch curve
> treats each year class as a random draw and carries that variability into its
> interval; in this module's simulations with variable year classes its 95%
> intervals covered the true Z in 95-98% of samples, as they should."

**Figures** (`results/fisheries-stock/si/fit/figures/`):
1. `catch_curve.png`: number caught at each age on a log scale. Filled circles = ages used by Chapman-Robson, the GLMM and the Bayesian model; open square = regression only (age 3); crosses = ages 1-2, not fully recruited. Lines: Chapman-Robson (solid), weighted regression (dashed, with its 95% band) and the Bayesian estimate for an average year class (dotted).
2. `mortality_estimates.png`: Z with its 95% interval from each of the four estimators; dashed line = the primary (Bayesian) estimate, 0.53.
3. `trace.png`: MCMC traces of Z and the year-class SD for the three chains.

**9. Gate.** The verdict is **DONE**. Every interpreted estimand (Z, annual
survival, annual mortality from the Bayesian catch curve) passes planted-truth
recovery at this design; the fit converged; the goodness-of-fit signal (c-hat 2.97)
is the year-class variability the model includes; the likelihood GLMM agrees;
negative controls fail as they should (VALIDATION.md section 5: unmodelled
recruitment trend +32% bias; Bayes with year-class SD forced to 0 covers 0.65).
Chapman-Robson and the weighted regression still fail coverage, so the report
quotes them only as cross-checks.

**10. Deliver.** The agent's summary for the researcher:
> Adult walleye (age 4+) die at a rate of Z = 0.53 per year (95% credible interval
> 0.39-0.74). That means about 59% survive each year and 41% die from fishing and
> natural causes combined. Your age data show strong and weak year classes (a
> typical year class is ~70% above or below average); the method allows for that,
> which is why the interval is wider than a simple catch-curve regression would
> suggest, and simulations at your sample size show the interval is trustworthy.
> Splitting the mortality into fishing and natural deaths needs a separate
> estimate of natural mortality.

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% interval) | truth inside? |
|---|---|---|---|
| Z (Bayesian, primary) | 0.55 | 0.53 (0.39-0.74) | yes |
| annual survival | 0.577 | 0.59 (0.48-0.68) | yes |
| annual mortality | 0.423 | 0.41 (0.32-0.52) | yes |
| year-class SD (log scale) | 0.40 | 0.54 (0.16-1.24) | yes |
| Z (GLMM, REML + t) | 0.55 | 0.50 (0.35-0.66) | yes |
| Z (Chapman-Robson) | 0.55 | 0.50 (0.33-0.66) | yes |
| Z (weighted regression) | 0.55 | 0.53 (0.33-0.72) | yes |

## (f) Limits that matter at this design (VALIDATION.md)

- One age sample can only average over the year classes it happens to contain; the interval is honest but wide (0.39-0.74). More years of ages (section h) narrow it much more than more fish in one year.
- A trend in recruitment (year classes getting steadily stronger or weaker) biases Z for every estimator: +32% in the negative control.
- Dome-shaped gill-net selectivity and under-ageing of old fish (otoliths) also bias Z.
- At harder designs (Z 0.8, 150 fish) the Bayesian Z is biased about +10% (still within the gate).
- These example data were simulated as a simple random sample of fish. Real gill-net ages come from a few net sets, each catching its own mix of ages (cluster sampling). In the independent verification (VALIDATION.md section 6) that biased the Bayesian Z by +8% (intra-haul correlation 0.1) to +15% (0.2, which fails the gate), against +2% to +8% for the GLMM and -2% to -3% for Chapman-Robson. With real net-set data, compare the three; if the Bayesian Z is clearly the highest, report that it may be too high (Nelson 2019).
- Z includes both fishing and natural mortality. F = Z - M needs an independent estimate of M.

## (g) Commands to reproduce (from `ecotools/`)

```
# regenerate the example data (seed 42)
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/fisheries-stock_si_make.log -- scripts/Rscript dev/fisheries-stock/make_examples.R
# planted-truth recovery at this design (about 10 min on 16 cores; Nimble compiles once per worker and number of ages)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/fisheries-stock_si_recover.log -- scripts/Rscript core/run.R recover fisheries-stock --config modules/fisheries-stock/example/recovery.yaml --nsim 200 --cores 16 --out results/fisheries-stock/si/recovery
# fit, report and figures
scripts/ejob -c 4 -m 16 -t 1:00:00 -l logs/fisheries-stock_si_fit.log -- scripts/Rscript core/run.R fit fisheries-stock --config modules/fisheries-stock/example/config.yaml --out results/fisheries-stock/si/fit
```

## (h) Several years of ages: the multi-year cohort catch curve

`multiyear_ages.csv` (simulated smallmouth bass, 6 fall surveys 2019-2024, ~200 aged
fish each; true Z = 0.6, year-class SD 0.5) has one row per year x age with a count.
`config_multiyear.yaml` sets `year_col: survey_year` and `count_col: n_fish`, which
switches on the cohort model: log E[catch] = year effect + year-class effect - Z age,
with each year class's strength shared across the years it was caught.

| estimator | Z (95% interval) | planted-truth recovery at this design (nsim 200) |
|---|---|---|
| cohort GLMM (default, REML + t) | 0.54 (0.44-0.64) | bias +0.2%, coverage 0.980, PASS |
| cohort Bayesian (Nimble; R-hat 1.002, ESS 906) | 0.55 (0.46-0.66) | bias +2.0%, coverage 0.960, PASS |
| Chapman-Robson, pooled over years | 0.57 (0.52-0.61): misses the truth | coverage 0.835, FAIL |
| weighted regression, pooled | 0.54 (0.49-0.59): misses the truth | coverage 0.825, FAIL |

The `catch_by_year` figure shows the weak 2017 year class moving along the ages
(low at age 3 in 2020, age 4 in 2021, ...), and `yearclass_strength` estimates it
at about 0.2 x average. Pooling the years hides this (pooled c-hat = 1.00), which is
why the pooled intervals are far too narrow.

```
scripts/ejob -c 4 -m 16 -t 1:00:00 -l logs/fisheries-stock_my_fit.log -- scripts/Rscript core/run.R fit fisheries-stock --config modules/fisheries-stock/example/config_multiyear.yaml --out results/fisheries-stock/example_multiyear/fit
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/fisheries-stock_my_recover.log -- scripts/Rscript core/run.R recover fisheries-stock --config modules/fisheries-stock/example/recovery_multiyear.yaml --nsim 200 --cores 16 --out results/fisheries-stock/example_multiyear/recovery
```
