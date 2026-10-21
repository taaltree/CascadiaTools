# Worked example: red-backed shrike integrated population model (module `ipm`)

Appendix S1 worked example and website "try an example" walkthrough.
**The data are simulated with known truth** (see `README_data.md`), so every
estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "We have followed the red-backed shrikes on our study plot for 15 years, 2006 to 2020.
> Every June we count the breeding females. We ring nestlings and breeding females, and
> we follow about half of the nests to see how many young fledge. The number of females
> went up and down a lot: from about 100 to over 200 and back. Is that because young
> birds survive their first year better in some years, because adults survive better,
> because breeding success changes, or because birds move in from elsewhere? How many of
> our breeders are immigrants? Would the plot hold its own without them? And do young
> birds survive better after wet years in the Sahel, where they pass on migration?"

## (b) The files the researcher has

`shrike_counts.csv`: one row per year (15 rows; 2013 is empty, no survey):
```
"year","females_counted"
2006,106
2007,106
2008,130
2009,139
```

`shrike_captures.csv`: one row per ring per breeding season (1,309 rows, 980 rings):
```
"ring","year","age","sex"
"AF11066",2006,"adult","F"
"AF25872",2006,"adult","F"
"AF27115",2006,"adult","F"
"AF32180",2006,"adult","F"
```

`shrike_broods.csv`: one row per monitored nest (1,020 rows):
```
"year","nest_id","female_age","fledglings"
2006,"N2006-001","yearling",1
2006,"N2006-002","yearling",5
2006,"N2006-003","yearling",3
2006,"N2006-004","yearling",4
```

`sahel_rainfall.csv`: one row per year (15 rows):
```
"year","sahel_rain_mm"
2006,499
2007,462
2008,420
2009,413
```

## (c) What the researcher types to Claude

> I have 15 years (2006-2020) of data on red-backed shrikes from my study plot, in four
> files: shrike_counts.csv (breeding females counted each June; no count in 2013),
> shrike_captures.csv (every ringed bird each year it was caught, ringed as a nestling or
> as a breeding female), shrike_broods.csv (fledglings from each monitored nest, with the
> mother's age, yearling or older) and sahel_rainfall.csv (July-September rainfall in the
> Sahel each year). I want to know how the number of breeding females changed, which
> vital rates drove the ups and downs, how many breeders are immigrants, and whether
> first-year survival is higher after wet years in the Sahel. I only have data on
> females. Please analyse this and give me a report with figures I can use in a paper.

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

**1. Project folder.** The data are copied to a project folder. The raw files are never edited.

**2. Data intake (eco-data-intake).** The agent reads the four files and describes them back:
- 15 breeding seasons, 2006-2020. Counts in 14 of them (99 to 209 females). 2013 has no
  count; the model fills it in from the other data.
- 980 rings: 700 nestlings (50 per year, 2006-2019) and 280 breeding females (20 new
  females per year). 229 rings were caught again in a later year; 127 of them had been
  ringed as nestlings and came back to breed.
- 1,020 monitored nests: 549 with yearling mothers (mean 2.91 fledglings) and 471 with older
  mothers (mean 3.62). 47 nests failed (0 fledglings) and are included, as they must be.
- Sahel rainfall for all 15 years: 191-522 mm (mean 420, SD 90). The value for year t acts
  on survival from t to t+1, so the 2020 value (191 mm, the lowest) is not used.
- Sex is `U` for nestlings and `F` for adults. The analysis is female-based, which matches.
- No duplicate counts, no impossible years, no unrecognised age codes.

**3. Choose the module.** The study has counts, capture-recapture and productivity data on
the same population, at least 8 years, and the question is about population change,
vital rates and immigration. That is the `ipm` module (integrated population model).
Counts alone (`timeseries`) or rings alone (`cjs`) could not estimate immigration or
explain the changes in numbers.

**4. Analysis plan (eco-analysis-plan), written before fitting.**
- Estimands: number of breeding females per year; mean first-year (juvenile) and adult
  apparent survival; fecundity of yearling and older mothers; immigration rate and number
  of immigrants; mean growth rate; the rainfall effect on first-year survival; which vital
  rate contributes most to the year-to-year variation in growth (LTRE ranking).
- Model: female-based, two age classes plus immigrants, random year effects on every
  rate, Poisson counts, one recapture probability, rainfall on juvenile survival (logit
  scale, per SD). Fitted in Nimble: 3 chains x 150,000 iterations (30,000 burn-in, thin 30).
  A separate analysis (rings + nests only) is fitted too, to show what the counts add.
- "Adequate fit": posterior predictive p-values between 0.025 and 0.975 for each data set;
  R-hat < 1.01 and bulk ESS >= 400 for every reported quantity.
- Stated up front, from VALIDATION.md: the year-to-year SDs and the sizes of the LTRE
  contributions are descriptive only, and the immigration contribution in the LTRE is
  overestimated. Only the ranking of the rates with direct data is interpreted.

**5. Planted truth at this design (eco-planted-truth).** Before touching the real data,
the agent simulates 100 data sets that look like this study: 15 years, about 90 females at
the start, 50 nestlings and 20 adults ringed per year, half of the nests monitored,
immigration about 12% per year, and a rainfall effect of +0.4 (logit, per SD). One fit
takes about 7 minutes (median 408 s with the faster recovery settings), so 100 replicates
were run (`recovery.yaml`). 96/100 fits converged (R-hat < 1.05); 0 errors. The coverage
bound is 0.883. **All the estimands that will be interpreted PASS:**

| estimand | truth | mean estimate | relative bias | 95% coverage | gate |
|---|---|---|---|---|---|
| juvenile (first-year) survival | 0.25 | 0.254 | +1.7% | 0.98 | PASS |
| adult survival | 0.48 | 0.483 | +0.5% | 0.97 | PASS |
| fecundity, yearling mothers | 2.8 | 2.795 | -0.2% | 0.91 | PASS |
| fecundity, older mothers | 3.6 | 3.622 | +0.6% | 0.95 | PASS |
| recapture probability | 0.6 | 0.596 | -0.6% | 0.94 | PASS |
| rainfall effect on juvenile survival (logit/SD) | 0.40 | 0.398 | -0.5% | 0.98 | PASS (power 0.875) |
| immigration rate (study-period mean) | 0.127 | 0.118 | -7.1% | 0.96 | PASS |
| immigrants per year | 13.3 | 12.4 | -9.6% | 0.98 | PASS |
| females in first year (N_first) | 89.8 | 92.6 | +3.2% | 0.98 | PASS |
| females in last year (N_final) | 112.6 | 114.1 | +1.3% | 0.98 | PASS |
| mean number of females | 105.2 | 105.3 | -0.0% | 0.94 | PASS |
| mean growth rate | 1.008 | 1.007 | -0.1% | 0.97 | PASS |
| juvenile survival SD (logit) | 0.30 | 0.307 | +2.3% | 0.97 | PASS |
| fecundity SDs (log) | 0.15 | 0.161 / 0.154 | +7.5% / +2.8% | 0.94 / 0.96 | PASS |
| LTRE: juvenile survival | 0.0213 | 0.0231 | +12.7% | 0.95 | PASS (within 3 MC SE) |
| adult survival SD (logit) | 0.20 | 0.233 | +16.5% | 0.97 | FAIL (bias) |
| immigration rate, median (not interpreted) | 0.12 | 0.103 | -14.2% | 0.97 | FAIL (bias) |
| immigration SD (log) | 0.30 | 0.751 | +150% | 1.00 | FAIL (bias) |
| LTRE: immigration | 0.0016 | 0.0073 | +927% | 0.94 | FAIL (bias) |
| LTRE: population structure | 0.00006 | 0.00026 | +168% | 0.98 | FAIL (bias) |
| LTRE: total explained variance | 0.0272 | 0.0357 | +39% | 0.95 | FAIL (bias) |

(The small LTRE terms for adult survival and fecundity pass only because 3 MC SE is larger
than 10% of a tiny truth; they are not interpreted.) **LTRE ranking:** juvenile survival was
the true top contributor in all 96 converged replicates. The fit ranked it first in 90/96
(94%); the 6 misses put immigration first. Among the rates with direct data (immigration
excluded) the ranking was right in 96/96.

What this tells the researcher. Population size, survival, fecundity, immigration rate and
the rainfall effect come out without meaningful bias and with honest intervals at this
design, and a rainfall effect of this size is detected 88% of the time. The failures are
the known limits: the year-to-year SDs of adult survival and immigration, and the sizes of
the LTRE contributions (especially immigration), are inflated. They are reported but not
interpreted.

**6. Fit.** The agent writes `config.yaml` (file names, column names, which rate the
rainfall acts on, figure labels) and runs the fit as a cluster job. With the separate
analysis and posterior predictive checks it takes about 30 minutes.

**7. Checks (eco-convergence, eco-model-checking).**
- Convergence: max R-hat 1.004 over all reported quantities; min bulk ESS 630 (the
  lowest values are for the 2013 population size and immigration, the year without a
  count). Trace plots: the 3 chains overlap with no trends (`trace` figure).
- Goodness of fit (posterior predictive checks, Freeman-Tukey discrepancy, 600 draws):
  Bayesian p = 0.56 for the counts, 0.46 for productivity and 0.37 for capture-recapture.
  All are far from 0 and 1: no sign of lack of fit for any data set.

**8. Interpretation (eco-interpretation).** Excerpts from the report
(`results/ipm/si/fit/ipm_report.html`):

> Population size. The model estimates 108 (95% CrI 91-127) breeding females in 2006 and
> 135 (95% CrI 116-157) in 2020 (counts were 106 and 135). Over the 15 years the population
> changed on average by a factor of 1.016 (95% CrI 1.000-1.033) per year, i.e. +1.6% per
> year (95% CrI +0.0% to +3.3%).

> Survival. Apparent survival from one breeding season to the next (at average conditions)
> was 0.22 (95% CrI 0.18-0.27) for first-year birds (from fledging to their first breeding
> season) and 0.45 (95% CrI 0.39-0.50) for adults. 'Apparent' means that animals which
> emigrate permanently are counted as dead [...]

> Fecundity. Females breeding for the first time (yearlings, and immigrants, which are
> assumed to breed like yearlings) raised 2.89 (95% CrI 2.63-3.17) fledglings (both sexes)
> per year; older females raised 3.66 (95% CrI 3.32-4.03).

> Immigration. [...] The estimated immigration rate, averaged over the study years, is 0.20
> (95% CrI 0.12-0.30) new female breeders per female already present, i.e. on average 26.8
> (95% CrI 15.8-38.1) immigrant breeding females per year, making up 17% (95% CrI 10-25%)
> of the breeders each year. Immigration is the least certain part of the model [...]

> Local recruitment alone. With the estimated survival and fecundity but no immigrants, the
> population would change by a factor of 0.81 per year (95% CrI 0.72-0.90; a decline: the
> local population behaves as a sink that is maintained by immigration).

> Covariate. Sahel rainfall (mm, Jul-Sep) was associated with juvenile survival: odds ratio
> 1.36 (95% CrI 1.02-1.80) per 1 SD of Sahel rainfall (1 SD = 90 [mm]); P(effect > 0) = 0.98.
> Across the observed range of Sahel rainfall (307 to 522), the expected juvenile survival
> goes from about 0.17 to 0.29 (other conditions average). This is an association across
> 14 years, not proof of cause.

> Caution on immigration in the LTRE. The immigration contribution (0.0422 (95% CrI
> 0.0032-0.1318)) and therefore the total explained variance are NOT reliable [...] Treat it
> as an upper bound. Among the rates with direct data, the largest contribution came from
> juvenile survival (0.0075 (95% CrI 0.0003-0.0187)).

> Benefit of the joint analysis. Compared with analysing the capture-recapture and
> productivity data separately, the IPM changed the posterior SD by: juvenile survival -7%;
> adult survival -22%; fecundity yearling -0%; fecundity adult -1%; recapture prob -2%.

The agent checked each sentence against `ipm_estimates.csv` (for example exp(0.305) = 1.36 and
exp(0.022) to exp(0.585) = 1.02 to 1.80; N_final 135 matches the 2020 count of 135). It adds one sentence for this project:
the LTRE in this data set ranks immigration first, but the planted-truth run shows that the
immigration contribution is inflated about tenfold at this design; the ranking among the rates
with direct data (juvenile survival first) is the supported result.

**Figures** (`results/ipm/si/fit/figures/`, PNG + PDF; gallery in `gallery.html`):

| figure | caption |
|---|---|
| `population_trajectory` | Estimated number of breeding females per year (line: posterior mean; bands: 95% and 50% credible intervals) with the observed counts (points). The wide band in 2013 is the year without a count. |
| `population_composition` | Composition of the breeding population each year (posterior means): local recruits (born on the plot the year before), surviving adults and immigrants. |
| `growth_rate` | Realised annual growth rate N[t+1]/N[t] (points between the two years; thick and thin lines: 50% and 95% CrI); dashed line = no change. |
| `vital_rates` | Annual juvenile survival, adult survival, fecundity of yearling and older mothers and immigration rate (posterior means, 50% and 95% CrI); dashed lines = mean rate. |
| `ltre` | Life-table response experiment: contribution of each vital rate and of population structure to the variance of the growth rate. Read the ranking, not the sizes; the immigration bar is overestimated. |
| `transient_ltre` | Contributions of each vital rate and of population structure to the change in growth rate from one year to the next; white points = estimated change. |
| `ipm_vs_separate` | Precision of the IPM versus separate analyses of the ring and nest data (95% CrI and coefficient of variation). The rainfall effect's CV drops from 63% to 48%; adult survival from 8% to 6%. |
| `effect_sahel_rain_mm` | Expected first-year survival versus Sahel rainfall (line: posterior mean; bands: 95% and 50% CrI); points: annual estimates. |
| `gof` | Posterior predictive checks for each data set (P = 0.37 capture-recapture, 0.56 counts, 0.46 productivity). |
| `trace` | MCMC traces of the main estimands (3 chains). |
| `posterior` | Posterior densities of the main estimands. |

The key figure is `population_composition`. It answers the main question in one picture:
how many breeders each year were local recruits, surviving adults and immigrants, and how
the 2013-2015 peak was built.

Figure fixes made by the agent (config `figures:` block only): in `population_composition`
the legend was moved below the plot and the figure widened, because the year labels
overlapped; in `growth_rate` the y-axis title was shortened and the x range set to
2005.5-2020.5, because the title and the 2020 label were cut off. The redrawn files are in
`figures/`. (The copies embedded in the HTML report use the module's default layout.)

**9. Acceptance gate (eco-acceptance-gate).** The plan, with the estimands to be
interpreted, was set before fitting. The planted-truth check passed at this design for
every interpreted estimand (96/100 converged). The failing estimands (adult survival SD,
immigration SD and median, LTRE sizes) are reported only as descriptive, as the plan said.
The fit converged (R-hat <= 1.004, ESS >= 630). All three posterior predictive checks are
fine. The module has independent routes (JAGS book code; `marked` for the CJS part) and a
negative control that fails as it should (VALIDATION.md). Every number comes from the
commands in (g). **Verdict: DONE**, with the immigration and LTRE limits stated.

**10. Deliver.** Summary given to the researcher:
> The number of breeding females on the plot grew from about 108 (95% CrI 91-127) in 2006
> to 135 (116-157) in 2020, about +1.6% per year (0.0 to +3.3%), with a peak of about 200 in
> 2015. First-year birds survived to breed with probability 0.22 (0.18-0.27) and adults
> 0.45 (0.39-0.50); yearling mothers raised 2.9 (2.6-3.2) fledglings and older mothers 3.7
> (3.3-4.0). The plot depends on immigrants: about 27 (16-38) immigrant females arrived per
> year, roughly a sixth of the breeders, and without them the local population would shrink
> by about 19% a year (growth 0.81, 0.72-0.90). Of the rates we measured directly, first-year
> survival drove most of the year-to-year changes, and it was higher after wet years in the
> Sahel (odds ratio 1.36 per 90 mm, 1.02-1.80; 0.17 at 307 mm vs 0.29 at 522 mm). This is
> an association, not proof of cause. Survival is "apparent": birds that settle elsewhere
> count as lost. The immigration estimate is the least certain number, and the model
> overstates immigration's role in the year-to-year changes.

## (e) Estimates versus the known truth

| parameter | truth | estimate (95% CrI) | truth inside CrI? |
|---|---|---|---|
| juvenile survival | 0.25 | 0.225 (0.184-0.268) | yes |
| adult survival | 0.48 | 0.448 (0.394-0.502) | yes |
| fecundity, yearling | 2.8 | 2.89 (2.63-3.17) | yes |
| fecundity, older | 3.6 | 3.66 (3.32-4.03) | yes |
| recapture probability | 0.6 | 0.654 (0.585-0.722) | yes |
| rainfall effect (logit per SD) | 0.40 | 0.305 (0.022-0.585) | yes |
| immigration rate (study mean) | 0.127 | 0.204 (0.115-0.296) | yes |
| immigrants per year | 17.5 | 26.8 (15.8-38.1) | yes |
| N 2006 (N_first) | 86 | 108 (91-127) | **no** |
| N 2020 (N_final) | 136 | 135 (116-157) | yes |
| mean N | 137.3 | 136.4 (128.9-144.8) | yes |
| mean growth rate | 1.0333 | 1.0165 (1.0000-1.0333) | at the edge (upper limit 1.03326 vs truth 1.03328) |
| juvenile survival SD | 0.30 | 0.146 (0.010-0.418) | yes |
| adult survival SD | 0.20 | 0.217 (0.017-0.529) | yes |
| fecundity SDs | 0.15 / 0.15 | 0.145 (0.068-0.252) / 0.156 (0.082-0.263) | yes |
| immigration SD | 0.30 | 1.29 (0.36-2.70) | **no** (prior-driven; known limit) |
| LTRE juvenile survival | 0.0191 | 0.0075 (0.0003-0.0187) | **no** (just above the upper limit) |
| LTRE immigration | 0.00003 | 0.042 (0.003-0.132) | **no** (known limit) |
| annual N, 15 years | see README_data.md | | 14 of 15 (2006 missed) |

The mean rates, the rainfall effect, recapture and the last-year and mean N are all
inside their intervals. Three things are worth explaining:

- **2006 is too high.** The first count (106) was well above the true 86 females (Poisson
  count error of +20). The first year has no earlier data to correct it, so the model
  follows the count. This pulls the mean growth rate down, so its interval only just reaches
  the truth. It is one data set: in the 100-replicate check, N_first was unbiased (+3.2%,
  coverage 0.98).
- **Immigration is high here** (0.20 vs 0.127) but its interval still contains the truth.
  The 100-replicate check shows no systematic bias at this design (-7%). Immigration and
  first-year survival trade off: this fit has juvenile survival a little low and immigration
  a little high.
- **The LTRE sizes are off, as expected.** The fit puts immigration first (0.042) when the
  truth is juvenile survival (0.019 vs 0.00003). This is the documented overestimation of the
  immigration contribution (Paquet et al. 2021). Among the rates with direct data, the fit
  correctly ranks juvenile survival first. This is why only that ranking is interpreted.

## (f) Limits that matter at this design (from `VALIDATION.md`)

- **The LTRE immigration contribution is overestimated severalfold** (+927% here, +937% to
  +2565% in VALIDATION.md). Read the ranking among the rates with direct data; treat the
  immigration bar as an upper bound. In this example the fitted LTRE puts immigration first,
  wrongly.
- **LTRE sizes are biased upward** (total explained variance +39% at this design). Only the
  ranking is supported (top contributor right in 94% of replicates, 100% excluding immigration).
- **Year-to-year SDs** are descriptive only with 15 years (adult survival SD +17%,
  immigration SD +150%, prior-driven).
- **Immigration is the least certain estimate.** It is unbiased at about 12% immigration
  (this design), but biased about +19% upward when immigration is low (about 8%). It also
  absorbs any error elsewhere: counts that miss females, adult emigration, ring loss.
- **Leaving immigration out would bias juvenile survival up by about 19%** (negative control,
  coverage 0.53). It must stay in for an open population like this one.
- **Apparent survival.** Emigration counts as death, so first-year "survival" is local
  recruitment.
- **Model assumptions:** counts are unbiased for breeding females; one recapture probability
  for nestling-ringed and adult-ringed birds; immigrants breed like yearlings; no density
  dependence; the three data sets are treated as independent.
- **The rainfall effect is an association** across 14 winters. Power to detect an effect of
  this size is 0.875 at this design.

## (g) Commands to reproduce (run from `ecotools/`)

```bash
# 0. regenerate the example data (seed 2024) - optional, the CSVs are already in the folder.
#    Note: this script also overwrites TRUTH.txt and drops the 2013 count on purpose.
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/ipm_make_example.log -- scripts/Rscript dev/ipm/make_example.R
# 1. planted-truth recovery at this design (100 replicates, ~7 min per fit, ~52 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/ipm_si_recover.log -- scripts/Rscript core/run.R recover ipm --config modules/ipm/example/recovery.yaml --nsim 100 --cores 16 --out results/ipm/si/recovery
# 2. fit the example data: report + figures (~30 min)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/ipm_si_fit.log -- scripts/Rscript core/run.R fit ipm --config modules/ipm/example/config.yaml --out results/ipm/si/fit
# 3. (optional) redraw figures after editing the figures: block, without refitting
scripts/ejob -c 1 -m 8 -t 0:20:00 -l logs/ipm_si_figures.log -- scripts/Rscript core/run.R figures ipm --config modules/ipm/example/config.yaml --out results/ipm/si/fit
# 4. machine-readable summary for the SI / website
python3 dev/ipm/make_si_summary.py
```

Outputs: `results/ipm/si/recovery/ipm_recover_summary.csv` (gate table),
`results/ipm/si/recovery/ipm_recover_reps.csv` (per-replicate estimates; the power for the
rainfall effect is computed from it, because the module's built-in power column refers to the
default covariate name), `results/ipm/si/fit/ipm_report.html` (report),
`results/ipm/si/fit/ipm_estimates.csv`, `results/ipm/si/fit/figures/`, and
`results/ipm/si/si_summary.json`.
