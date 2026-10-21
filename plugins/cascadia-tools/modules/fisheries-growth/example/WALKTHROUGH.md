# Worked example: do female and male walleye grow differently?

Module: `fisheries-growth`. Appendix S1 and website walkthrough.
**The data are simulated** (see `README_data.md`), so we know the true answer and can check it.

## (a) The question, in the researcher's words

> "We caught 1,000 walleye in our October 2014 gill-net survey. Do females and males differ
> in how big they get and how fast they grow? And are males or females in better condition?"

## (b) The files the researcher has

`sim_walleye_2014.csv`: one row per fish (1,000 rows).

```
"fish_id","year","sex","age","tl_mm","wt_g"
"SIM14-F-0001",2014,"female",1,360,454.3
"SIM14-F-0002",2014,"female",2,431,606.5
"SIM14-F-0003",2014,"female",5,567,2065.4
"SIM14-F-0004",2014,"female",1,344,416.4
```

That is the only data file. Ages are otolith ages in whole years. Length is total length (mm). Weight is wet weight (g).

## (c) What the researcher types to Claude

> "I have age, length, weight and sex for 1,000 walleye from an October 2014 gill-net
> survey (sim_walleye_2014.csv). Do females and males differ in growth and in condition?
> I'd like a report and figures I can put in a paper (CJFAS format)."

## (d) What the agent does, step by step (eco-analyze)

**1. Project folder.** The agent keeps all work outside the raw data. In this example the files sit in
`modules/fisheries-growth/example/`, and outputs go to `results/fisheries-growth/si/`.

**2. Data intake.** The agent describes the data back:
- 1,000 fish: 530 female, 470 male. One year (2014).
- Otolith ages 1 to 16. Most fish are young: 424 are age 1, 249 age 2, 136 age 3.
- Only 19 females and 54 males are older than about the age at 90% of their maximum size.
  That is enough to estimate maximum length, but only just for females.
- Lengths 312 to 741 mm. Every fish was weighed. No missing values.
- Fish were caught in October, so the agent adds 0.5 year to otolith age (`age_offset: 0.5`).

**3. Choose the module.** `fisheries-growth` fits a von Bertalanffy growth curve
(L = Linf(1 - e^(-K(age - t0)))) for each sex and compares the curves. It also fits
the length-weight relationship (W = aL^b) and condition (relative weight, Wr). A mixed
model is not needed because there are only two groups, so the agent uses the fast
likelihood engine (`nls`) with Kimura's tests.

**4. Analysis plan.**
- Estimands: Linf, K, t0 and age at 90% Linf for each sex.
- Also: the male minus female difference in each growth parameter, the length-weight slope b,
  the male:female weight ratio at the same length, and Wr by sex.
- Error: multiplicative (lognormal), because spread in length grows with size.
- Adequate fit means: the optimiser converged, all standard errors are finite, residuals show
  no trend against fitted length or age, and residual spread is even.

**5. Planted-truth check (at this exact design).** The agent simulates 200 new datasets with
the same design (530 + 470 fish, same age structure) and the true values from the module's
example, then fits each one. Each fit takes about 0.2 s.
- **Result: all 17 quantities PASS the gate** (|bias| at most 4.2% of truth; 95% coverage
  0.925 to 0.965; 200 of 200 fits converged).

| Quantity | Truth | Mean estimate | Rel. bias | Coverage | Power |
|---|---|---|---|---|---|
| Linf female (mm) | 680 | 681.0 | +0.2% | 0.95 | |
| Linf male (mm) | 560 | 559.5 | -0.1% | 0.945 | |
| K female | 0.31 | 0.309 | -0.2% | 0.93 | |
| K male | 0.39 | 0.392 | +0.5% | 0.965 | |
| t0 female | -1.0 | -1.006 | -0.6% | 0.925 | |
| t0 male | -1.1 | -1.102 | -0.2% | 0.96 | |
| Age at 90% Linf, female | 6.43 | 6.46 | +0.6% | 0.94 | |
| Age at 90% Linf, male | 4.80 | 4.80 | -0.1% | 0.94 | |
| Linf difference (M - F) | -120 | -121.5 | +1.3% | 0.955 | 1.00 |
| K difference (M - F) | 0.08 | 0.083 | +3.2% | 0.935 | 0.81 |
| t0 difference (M - F) | -0.10 | -0.096 | -4.2% | 0.93 | 0.10 |
| Condition ratio (M / F) | 1.02 | 1.020 | -0.03% | 0.945 | 0.93 |
| Length-weight b | 3.34 | 3.341 | +0.03% | 0.965 | |

  The other quantities (a, both residual SDs, mean weight of a 500 mm female) also pass.
  "Power" is how often the 95% interval excludes zero (or 1 for the ratio). At this design a
  120 mm difference in Linf is always detected. A 0.08 difference in K is detected 81% of the
  time. A 0.1 year difference in t0 is almost never detected.

**6. Fit.** The agent writes `config_si.yaml` (file, column names, `reference_group: female`,
`age_offset: 0.5`, walleye standard weight, CJFAS figure settings) and runs the fit.

**7. Checks.**
- Convergence: the optimiser converged ("relative convergence"). All standard errors are finite.
  No estimate is on a boundary. (This is a likelihood fit, so R-hat and ESS do not apply.)
- Growth residuals: mean 0, SD 0.050. A smooth line through residuals against fitted length
  never moves more than 0.005 from zero, so the curve shape fits. Shapiro-Wilk p = 0.16.
  Box plots by age are centred on zero with similar spread.
- Length-weight: one common slope fits both sexes (F = 0.00, df 1 and 996, p = 0.97).
  Residuals are close to normal (Shapiro-Wilk p = 0.12).
- One fish (SIM14-F-0433) has a weight residual of 4.3 SD. The report flags it as a possible
  data-entry error. In simulated data this is just a rare chance value, so it is kept.

**8. Interpretation.** Excerpt from the report:

> "female: fish approach an asymptotic length (Linf) of 685 mm (95% CI 665 to 706) and grow
> toward it at K = 0.302 (95% CI 0.272 to 0.335) per year. They reach 90% of Linf (about
> 617 mm) at age 6.6 (95% CI 6.0 to 7.2) years."
>
> "male: fish approach an asymptotic length (Linf) of 561 mm (95% CI 547 to 575) and grow
> toward it at K = 0.385 (95% CI 0.339 to 0.437) per year. They reach 90% of Linf (about
> 505 mm) at age 4.8 (95% CI 4.3 to 5.4) years."
>
> "Linf in male is smaller than in female by 124 mm (95% CI 99 to 149). K in male is larger
> than in female by 0.084 per year (95% CI 0.025 to 0.142)."
>
> "Kimura (1980) likelihood-ratio tests ... common linf: chi-square = 81.7, df = 1,
> p < 0.0001; common k: chi-square = 8.3, df = 1, p = 0.0041; common t0: chi-square = 0.4,
> df = 1, p = 0.51."
>
> "At the same length, fish in male weigh 2.1% more than fish in female (95% CI +0.9% to +3.3%)."

In plain words: females grow larger, males level off sooner. At the same length, males are
about 2% heavier. Both sexes average Wr of about 90 (female 90.2, male 91.1). That is
typical, because the walleye standard is set at the 75th percentile of populations.

**Figures** (`results/fisheries-growth/si/fit/figures/`, PNG and PDF, CJFAS size):

| Figure | Caption (short) |
|---|---|
| `growth_curves.png` | Length at age for each sex, with fitted von Bertalanffy curve, 95% confidence band and 95% prediction interval; points are fish. **Key figure.** |
| `parameter_comparison.png` | Linf, K, t0 and age at 90% Linf by sex, with 95% intervals. |
| `residuals.png` | Growth residuals against fitted length (A) and by age (B). |
| `length_weight.png` | Length-weight relationship on log-log axes, common slope b = 3.340, with 95% bands. |
| `condition.png` | Relative weight (Wr) by sex: mean with 95% CI over individual fish; dashed line Wr = 100. |

**9. Acceptance gate: DONE.**
- Planted truth passed at this design (17 of 17).
- The fit converged and passed the residual checks.
- An independent route already exists for this module: Kimura statistics and estimates match
  `fishmethods::growthlrt` exactly (VALIDATION.md).
- A negative control already exists: fitting with length-selective gear fails the gate (VALIDATION.md).
- Every number here comes from the commands in (g).

**10. Deliver.** The report is `results/fisheries-growth/si/fit/fisheries-growth_report.html`.
It contains methods, estimates, interpretation, figures and references.

## (e) Estimates versus the known truth

Every true value lies inside its 95% interval.

| Quantity | Truth | Estimate (95% CI) |
|---|---|---|
| Linf female (mm) | 680 | 685 (665 to 706) |
| Linf male (mm) | 560 | 561 (547 to 575) |
| K female (per year) | 0.31 | 0.302 (0.272 to 0.335) |
| K male (per year) | 0.39 | 0.385 (0.339 to 0.437) |
| t0 female (years) | -1.0 | -1.04 (-1.21 to -0.86) |
| t0 male (years) | -1.1 | -1.14 (-1.38 to -0.89) |
| Age at 90% Linf, female | 6.43 | 6.60 (5.96 to 7.24) |
| Age at 90% Linf, male | 4.80 | 4.84 (4.32 to 5.37) |
| Linf difference (M - F, mm) | -120 | -124 (-149 to -99) |
| K difference (M - F) | 0.08 | 0.084 (0.025 to 0.142) |
| t0 difference (M - F) | -0.10 | -0.10 (-0.40 to 0.20) |
| Residual SD (log length) | 0.05 | 0.0504 (0.0483 to 0.0528) |
| Length-weight b | 3.34 | 3.340 (3.308 to 3.372) |
| Length-weight a | 1.2e-6 | 1.19e-6 (0.98e-6 to 1.45e-6) |
| Residual SD (log weight) | 0.09 | 0.0901 (0.0863 to 0.0942) |
| Mean weight, 500 mm female (g) | 1,245.9 | 1,239.7 (1,229.1 to 1,250.4) |
| Weight ratio male / female | 1.02 | 1.021 (1.009 to 1.033) |

## (f) Limits that matter at this design (from VALIDATION.md)

- **Linf needs old fish.** Here 19 females are older than about the age at 90% of Linf. That is enough,
  but the female Linf interval (41 mm wide) is wider than the male one (28 mm). With fewer
  than 5 old fish, the module warns that Linf is unreliable.
- **Length-selective gear biases growth.** Gill nets catch the faster-growing young fish.
  In tests this biased K by -14% and t0 by -0.36 years, and coverage fell to 0.58. These
  simulated data have no selectivity. With real gill-net data, drop ages that are only
  partly selected.
- **Single-parameter differences have low power.** At this design a 0.08 difference in K
  is detected 81% of the time. A 0.1 year difference in t0 is detected about 10% of the time. "No
  difference in t0" here means "not detectable", not "equal".
- **Not modelled:** ageing error, individual growth variation beyond the residual, and
  growth models other than von Bertalanffy.

## (g) Commands to reproduce (run from the `ecotools/` folder)

```
# 1. Make the simulated data (seed 20141015)
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/fisheries-growth_si_make.log -- \
  scripts/Rscript dev/fisheries-growth/make_si_example.R

# 2. Planted-truth recovery at this design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/fisheries-growth_si_recover.log -- \
  scripts/Rscript core/run.R recover fisheries-growth \
  --config modules/fisheries-growth/example/recovery.yaml --nsim 200 --cores 16 \
  --out results/fisheries-growth/si/recovery

# 3. Fit the example data (report + figures)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/fisheries-growth_si_fit.log -- \
  scripts/Rscript core/run.R fit fisheries-growth \
  --config modules/fisheries-growth/example/config_si.yaml \
  --out results/fisheries-growth/si/fit
```

The same analysis on the real Lake Erie data uses `--config modules/fisheries-growth/example/config.yaml`.
