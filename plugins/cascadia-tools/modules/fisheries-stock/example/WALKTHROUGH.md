# Worked example: how fast are adult walleye dying? (catch curve)

Module: `fisheries-stock` (`model: catch_curve`). The data are **simulated**
(see `README_data.md`), so we can compare the answer with the truth.
All numbers below come from runs made on 2026-10-20.

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

**3. Choose the module.** The data are ages of individual fish and the question is
total mortality, so the agent uses `fisheries-stock` with `model: catch_curve`.
Length is not needed.

**4. Plan.**
- Estimand: total instantaneous mortality Z for fully recruited fish. Annual survival is S = exp(-Z) and annual mortality is A = 1 - S.
- Primary estimator: weighted regression of log catch on age, from the peak age (age 3) on.
- Cross-checks: Chapman-Robson (ages 4 and older) and a Poisson-lognormal GLMM.
- Assumptions to tell the researcher: constant recruitment with no trend, constant Z, full selectivity from the first age used, and accurate ages.

**5. Planted truth** (`recovery.yaml`, 200 simulated surveys at this design:
348 fish, Z = 0.55, year-class SD 0.4, the same gill-net selectivity).

| parameter | truth | mean estimate | rel. bias | 95% CI coverage | gate |
|---|---|---|---|---|---|
| Z (weighted regression, primary) | 0.550 | 0.547 | -0.5% | 0.900 | FAIL (coverage; bound 0.904) |
| annual survival | 0.577 | 0.581 | +0.7% | 0.900 | FAIL (coverage) |
| annual mortality | 0.423 | 0.419 | -0.9% | 0.900 | FAIL (coverage) |
| Z (Chapman-Robson) | 0.550 | 0.525 | -4.6% | 0.850 | FAIL (coverage) |
| Z (GLMM) | 0.550 | 0.513 | -6.8% | 0.785 | FAIL (coverage) |

In one sentence, the agent tells the researcher: *the method finds the right answer
on average (bias under 1%), but at this design its 95% interval contains the truth
only 90% of the time.* All 200 fits ran without errors.
The intervals are too narrow because year-class strength varies. The average reported SE
for Z is 0.070, but the estimates actually spread with SD 0.090.
This is the known limit recorded in VALIDATION.md. The gate is failed on coverage only.
The agent continues, but it reports the interval as roughly a 90% interval and gives the
verdict OPEN.

**6. Fit** (`config.yaml`, about 1 minute in total).

**7. Checks.**
- Convergence: this is not MCMC. All three estimators returned finite estimates, so R-hat and ESS do not apply. The report's diagnostics line prints "max R-hat NA" for this reason.
- Goodness of fit: the age distribution is overdispersed. Chapman-Robson c-hat = 2.97 (Pearson chi-square 23.8 on 8 df, P = 0.002), which points to variable year-class strength. The Chapman-Robson SE was inflated by sqrt(c-hat).
- The oldest ages got weight 0 in the weighted regression, because their predicted log catch was 0 or below. This is normal under the Maceina & Bettoli weighting.

**8. Interpretation** (quoted from the report):

> "Total instantaneous mortality of fully recruited walleye (ages 3 and older) is
> Z = 0.53 per year (95% CI 0.33 to 0.72; weighted regression). Equivalently, about
> 59% (48% to 72%) of fully recruited fish survive each year and 41% (28% to 52%)
> die from all causes combined (fishing plus natural mortality)."

> "The other estimators give Z = 0.50 (0.33 to 0.66; Chapman-Robson) and Z = 0.47
> (0.33 to 0.61; Poisson-lognormal GLMM). All estimators agree within their uncertainty."

> "...treat the interval as closer to a 90% interval."

**Figures** (`results/fisheries-stock/si/fit/figures/`):
1. `catch_curve.png`: number caught at each age on a log scale. The legend marks which ages each method used (filled circles = used by both, open square = regression only, crosses = ages 1-2, not fully recruited). The lines are the Chapman-Robson fit (solid) and the weighted regression fit (dashed) with its 95% band.
2. `mortality_estimates.png`: Z with its 95% interval from each of the three estimators. The dashed line marks the primary estimate (0.53).

**9. Gate.** The verdict is **OPEN**. Bias is small, the fit is fine and the
figures are correct. But the planted-truth coverage (0.900) is just below the gate
bound (0.904), and this is a documented limit of every catch-curve estimator when
year classes vary.

**10. Deliver.** The agent's summary for the researcher:
> Adult walleye (age 3+) die at a rate of Z = 0.53 per year (95% CI 0.33-0.72).
> That means about 59% survive each year and 41% die from fishing and natural
> causes combined. Year-class strength varies a lot, so treat the interval as
> roughly 90%, not 95%. Splitting the mortality into fishing and natural
> deaths needs a separate estimate of natural mortality.

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% CI) | truth inside CI? |
|---|---|---|---|
| Z (weighted regression) | 0.55 | 0.53 (0.33-0.72) | yes |
| annual survival | 0.577 | 0.59 (0.48-0.72) | yes |
| annual mortality | 0.423 | 0.41 (0.28-0.52) | yes |
| Z (Chapman-Robson) | 0.55 | 0.50 (0.33-0.66) | yes |
| Z (GLMM) | 0.55 | 0.47 (0.33-0.61) | yes |

## (f) Limits that matter at this design (VALIDATION.md)

- Every catch-curve interval is too narrow when year classes vary: the 95% interval covers the truth only 84-93% of the time (90% here). Read it as about a 90% interval.
- Chapman-Robson and the GLMM read Z low by 5-7% at this design. The weighted regression bias is under 1%.
- A trend in recruitment biases Z. In the negative control, Z was biased +29%.
- Dome-shaped gill-net selectivity and under-ageing of old fish (otoliths) also bias Z.
- Z includes both fishing and natural mortality. F = Z - M needs an independent estimate of M.

## (g) Commands to reproduce (from `ecotools/`)

```
# regenerate the example data (seed 42)
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/fisheries-stock_si_make.log -- scripts/Rscript dev/fisheries-stock/make_examples.R
# planted-truth recovery at this design (about 1 min)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/fisheries-stock_si_recover.log -- scripts/Rscript core/run.R recover fisheries-stock --config modules/fisheries-stock/example/recovery.yaml --nsim 200 --cores 16 --out results/fisheries-stock/si/recovery
# fit, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/fisheries-stock_si_fit.log -- scripts/Rscript core/run.R fit fisheries-stock --config modules/fisheries-stock/example/config.yaml --out results/fisheries-stock/si/fit
```
