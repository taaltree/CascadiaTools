# Worked example: western pond turtle abundance, recruitment and survival (module `jolly-seber`)

Appendix S1 worked example and website "try an example" walkthrough.
**The data are simulated with known truth** (see `README_data.md`), so every
estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "We have hoop-trapped western pond turtles at the same ponds every summer
> since 2015. Every turtle gets a PIT tag. I want to know how many turtles
> live there each year, and how many turtles in total used the site. How many
> new turtles join each year? What share survive from one year to the next,
> and do males survive less well than females? Is the population stable or
> declining? We trapped harder in some years than others, and a few turtles
> could not be sexed."

## (b) The files the researcher has

`turtle_captures.csv`: one row per capture (459 rows):
```
"turtle_id","pit_tag","capture_date","season","trap","sex","carapace_mm","mass_g"
"AM-089","982.149470056","2015-05-21",2015,"HT02","M",164,585
"AM-150","982.176990615","2015-05-23",2015,"HT04","F",145,759
"AM-134","982.676811193","2015-05-23",2015,"HT10","M",127,814
"AM-082","982.146166475","2015-05-26",2015,"HT10","F",145,746
```

`turtle_seasons.csv`: one row per trapping season (8 rows):
```
"season","start_date","end_date","trap_nights","observers"
2015,"2015-05-20","2015-08-03",93,2
2016,"2016-05-20","2016-08-03",80,2
2017,"2017-05-20","2017-08-03",142,3
2018,"2018-05-20","2018-08-03",163,3
```

## (c) What the researcher types to Claude

> I have eight summers (2015 to 2022) of hoop-trapping data for western pond
> turtles at one pond complex, in two files: turtle_captures.csv (every capture,
> with the turtle's ID, the season and its sex; U means we could not sex it) and
> turtle_seasons.csv (trap-nights of effort in each season). I want to know how
> many turtles were there each year and in total, how many new turtles joined
> each year, what share survive from year to year and whether males survive
> less well than females, and whether the population is declining. Please
> analyse this and give me a report with figures I can use in a paper.

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

**1. Project folder.** The data are copied to a project folder. The raw files are never edited.

**2. Data intake (eco-data-intake).** The agent reads both files and describes them back:
- 175 turtles: 100 female, 65 male, 10 of unknown sex (6%).
- 459 captures over 8 seasons (2015-2022). Repeat captures in the same summer
  are merged. That leaves 322 turtle-seasons.
- Turtles caught per season: 34, 23, 47, 55, 38, 44, 46, 35.
- 88 turtles (50%) were caught in only one season; 87 were caught in two or more.
- Effort ranged from 80 to 163 trap-nights (mean 120). Every season has an effort value.
- No missing IDs, no unknown season codes, no turtle with two different sexes.

**3. Choose the module.** The turtles are individually marked and caught over
several years in an open population. The questions are abundance, recruitment,
survival and growth. So the module is `jolly-seber`. (`cjs` would give survival
only. `scr` would need trap locations and gives density.)

**4. Analysis plan (eco-analysis-plan), written before fitting.**
- Estimands: turtles present each year (N_t); total turtles that used the site
  (superpopulation); new turtles per year (entrants) and per-capita recruitment;
  apparent annual survival for females and males; the sex effect on survival;
  realised growth rates and mean growth rate; capture probability and the effect
  of trap-nights on it.
- Model: Bayesian Jolly-Seber (superpopulation form, data augmentation).
  Survival constant over years with a sex effect. Capture probability depends on
  trap-nights (standardised). Turtles of unknown sex are a mixture of the two sexes.
- Engine: Nimble, marginal (forward-algorithm) likelihood; 3 chains x 20,000
  iterations. Maximum-likelihood POPAN (`marked`) is fitted alongside as a cross-check.
- "Adequate fit": R-hat < 1.01 and bulk ESS >= 400; posterior predictive
  p-values between 0.025 and 0.975; no R2ucare test with p < 0.05.
- Caveat stated up front: "apparent survival" = survival x staying at the site.

**5. Planted truth at this design (eco-planted-truth).** Before touching the
real data, the agent simulates 200 turtle studies like this one (8 years,
about 260 turtles in total, 45% present at the start, survival 0.82 for
females and 0.75 for males, 45% males, 6% unknown sex, capture 0.40,
effort effect +0.45 per SD). It fits the same model to each one. Each fit
takes about 11 seconds. All 200 fits converged (largest R-hat 1.016).
Result: **42 of 45 estimands PASS** (coverage must be at least 0.904):

| estimand | truth | mean estimate | relative bias | 95% coverage | power | gate |
|---|---|---|---|---|---|---|
| superpopulation size | 260 | 263.0 | +1.1% | 0.92 | | PASS |
| N 2015 ... N 2022 (8 values) | 102-118 | 104-116 | -2.5% to +2.8% | 0.925-0.965 | | PASS (all 8) |
| female survival | 0.82 | 0.815 | -0.6% | 0.955 | | PASS |
| male survival | 0.753 | 0.754 | +0.1% | 0.95 | | PASS |
| sex effect on survival (logit) | -0.40 | -0.366 | +8.4% (toward 0) | 0.975 | 0.265 | PASS |
| share of males | 0.451 | 0.450 | -0.3% | 0.985 | | PASS |
| capture probability | 0.40 | 0.401 | +0.2% | 0.925 | | PASS |
| effort effect on capture (logit/SD) | 0.45 | 0.458 | +1.8% | 0.975 | 0.99 | PASS |
| mean growth rate | 0.979 | 0.988 | +0.8% | 0.925 | | PASS |
| entry probability, 2015 | 0.45 | 0.427 | -5.1% | 0.955 | | PASS |
| entrants 2017-2022 (6 values) | 19.6-20.5 | 19.8-21.5 | -1.5% to +6.7% | 0.95-0.985 | | PASS (all 6) |
| growth rates, recruitment rates 2-7, entry probabilities 3-8 | | | -2.0% to +10.5% | 0.94-0.995 | | PASS (all 18) |
| **entrants 2016** | 20.8 | 25.6 | **+26%** | 0.975 | | **FAIL** (bias) |
| **entry probability, 2016** | 0.079 | 0.099 | **+26%** | 0.975 | | **FAIL** (bias) |
| **per-capita recruitment, 2015-16** | 0.178 | 0.253 | **+47%** | 0.97 | | **FAIL** (bias) |

What this tells the researcher:
- Abundance, total turtles, survival of each sex, the sex ratio, capture
  probability and the mean growth rate come out without bias and with honest
  intervals at this design.
- **The number of new turtles in the first interval (2015-2016) is not
  reliable.** It is biased about 25-50% high, although its intervals still
  cover. In the first interval, new arrivals are hard to tell apart from
  turtles that were already there but not caught in 2015. This is Known limit
  3 in `VALIDATION.md`. The agent therefore does not interpret 2015-2016
  recruitment. It reports the later years and the abundance trajectory instead.
- **Power for the sex effect is low (27%).** If males really survive only a
  little less well (logit -0.4), a study like this detects it about one time
  in four. An interval that includes zero would be inconclusive.

**6. Fit.** The agent writes `config.yaml` (file names, columns, sex coding,
effort covariate, figure labels) and runs the fit. It takes about 2 minutes,
including compilation and the `marked` cross-check.

**7. Checks (eco-convergence, eco-model-checking).**
- Convergence: max R-hat **1.001**, min bulk ESS **2,583** (thresholds 1.01
  and 400). The trace plots show the three chains overlapping.
- Data augmentation: M = 875 potential turtles. The 99.5% posterior quantile
  of the superpopulation is 328, well below 0.9 M, so M does not cap the estimate.
- Posterior predictive checks: Bayesian p = **0.52** (captures per turtle) and
  **0.62** (turtles caught per season). Neither is extreme, so there is no sign
  of capture heterogeneity.
- R2ucare tests by sex: transients (3.SR) p = 0.67 (female), 0.90 (male),
  0.63 (unknown); trap dependence (2.CT) p = 0.85, 0.73, 1.00; overall CJS
  test p = 0.88, 0.97, 0.53. No test flags lack of fit. (These tests have low
  power with 8 occasions.)
- Cross-check: `marked` ML POPAN on the same data gives superpopulation 264
  (Nimble 268), N 2015 114 (108), N 2022 111 (112), capture 0.375 (0.374) and
  effort effect 0.388 (0.376). marked dropped the sex effect because 10
  turtles have unknown sex, so its survival (0.79) is pooled over both sexes.

**8. Interpretation (eco-interpretation).** Excerpts from the report
(`results/jolly-seber/si/fit/jolly-seber_report.html`):

> 175 western pond turtles were caught at least once over 8 occasions (2015 to
> 2022); 87 (50%) were caught on more than one occasion, and 40.2 were caught
> per occasion on average.

> Estimated abundance was 108 (95% CI 76 to 148) at the first occasion (2015)
> and 112 (95% CI 87 to 145) at the last (2022). Across the study the average
> growth rate was 1.006 (95% CI 0.960 to 1.055) per year, i.e. a change of
> +0.6% per year; the posterior probability that the population declined
> (growth rate < 1) is 0.39.

> In total 268 (95% CI 236 to 310) western pond turtles are estimated to have
> used the study population at some time during the study (the
> superpopulation), of which 175 (65%) were ever caught.

> Apparent survival per year was 0.85 (95% CI 0.78 to 0.90) for female and
> 0.68 (95% CI 0.58 to 0.78) for male. The odds of surviving were 0.39 times
> as high for male as for female (95% CI 0.21 to 0.71); posterior probability
> that male animals survive less well than female animals > 0.99. 'Apparent'
> survival means that death and permanent emigration cannot be told apart.

> Capture probability of an animal present was 0.37 (95% CI 0.31 to 0.44) at
> average sampling effort. Capture probability rose from 0.29 to 0.47 between
> trapping effort (trap-nights) = 90 and 150 [...]. This is a sampling
> (detection) effect, not a change in the population.

The agent checked each sentence against the estimates table. It adds one
caveat for the researcher. The report's average per-capita recruitment (0.23,
95% CI 0.16-0.30) includes the 2015-2016 interval, which failed planted truth.
Read recruitment from the 2016-2017 to 2021-2022 intervals (entrants 11-35
turtles per year) and from the abundance trajectory.

**Figures** (`results/jolly-seber/si/fit/figures/`, PNG + PDF):

| figure | caption |
|---|---|
| `abundance_trajectory` | Turtles present each year (line: posterior mean; bands: 95% and 50% credible intervals). Black points: turtles actually caught that year (a minimum count). |
| `superpopulation` | Posterior of the total number of turtles that used the site 2015-2022; dashed line = 175 turtles caught. M = 875 does not constrain it. |
| `survival_sex` | Apparent annual survival by sex (points: posterior means; bars: 95% and 50% credible intervals). |
| `coefficients` | Odds ratios: survival of males vs females; capture per +1 SD of trap-nights. Dashed line: no effect. |
| `recruitment` | (A) New turtles entering at each season; (B) per-capita recruitment per interval. The 2016 / 2015-16 values are not reliable at this design (planted truth). |
| `entry_probabilities` | Share of all turtles that entered at each season; the 2015 value is the share already present at the start. |
| `growth_rate` | Realised growth rate N(t+1)/N(t) for each interval; dashed line = no change; green line and band = mean growth rate per year with its 95% interval. |
| `detection_time` | Capture probability of a turtle present in each season, as predicted from trap-nights. |
| `detection_trap_nights` | Capture probability against trapping effort (trap-nights), with 95% and 50% bands; points = turtles caught / estimated abundance in each season. |
| `gof` | Posterior predictive checks (replicate minus observed Freeman-Tukey discrepancy); p = 0.52 and 0.62. |
| `trace`, `posterior` | MCMC traces and posterior densities of key quantities, 3 chains. |
| `engine_comparison` | Nimble vs `marked` ML POPAN on the same data. The marked survival is pooled over sexes (the panel title says "female" because that is the Nimble reference group). |

Key figure: `abundance_trajectory`.

**9. Gate (eco-acceptance-gate): DONE**, with one estimand group excluded.
Planted truth passes for every estimand that is interpreted. The fit converged.
The posterior predictive checks and R2ucare tests show no misfit. The module was
cross-checked against `marked` and a latent-state Nimble model, and its negative
controls fail as they should (`VALIDATION.md`). 2015-2016 recruitment (entrants 2016,
entry probability 2016, per-capita recruitment 2015-16) failed the bias criterion
and is not interpreted.

**10. Deliver.** Summary for the researcher:
> About 100-115 turtles used the ponds in any one year (2015: 108, 95% CI
> 76-148; 2022: 112, 87-145). About 268 turtles (236-310) used the site at
> some point over the eight years. The population looks stable: mean growth
> was +0.6% per year (95% CI -4.0% to +5.5%). Females have higher apparent
> annual survival (0.85, 0.78-0.90) than males (0.68, 0.58-0.78). Roughly
> 11-35 new turtles joined each year from 2017 on. The 2016 figure is not
> reliable with this design. Capture probability was about 0.37 per season
> and rose with trap-nights, so effort must stay in the model.

## (e) Estimates versus the known truth

| quantity | truth (simulated) | estimate (95% CI) | inside CI? |
|---|---|---|---|
| superpopulation | 260 | 268 (236-310) | yes |
| N 2015 | 97 | 108 (76-148) | yes |
| N 2016 | 100 | 101 (70-138) | yes |
| N 2017 | 101 | 114 (92-140) | yes |
| N 2018 | 109 | 110 (92-134) | yes |
| N 2019 | 102 | 112 (93-135) | yes |
| N 2020 | 107 | 99 (81-120) | yes |
| N 2021 | 105 | 105 (85-130) | yes |
| N 2022 | 107 | 112 (87-145) | yes |
| entrants 2017-2022 | 16, 31, 17, 28, 21, 27 | 35, 21, 24, 11, 25, 29 | yes (all 6) |
| entrants 2016 (not interpreted) | 23 | 16 (0-46) | yes |
| mean growth rate | 1.014 | 1.006 (0.960-1.055) | yes |
| female survival | 0.82 | 0.845 (0.783-0.902) | yes |
| male survival | 0.753 | 0.685 (0.584-0.782) | yes |
| sex effect on survival (logit) | -0.40 | -0.93 (-1.54 to -0.34) | yes (near the upper edge) |
| share of males | 0.45 | 0.432 (0.378-0.488) | yes |
| capture probability | 0.40 | 0.374 (0.309-0.444) | yes |
| effort effect (logit/SD) | 0.45 | 0.376 (0.142-0.613) | yes |
| entry probability 2015 | 0.45 | 0.395 (0.284-0.523) | yes |

Every true value lies inside its 95% interval. The sex effect estimate
(-0.93) is about twice the true value (-0.40). Its interval still covers the
truth, but only just. This is what low power looks like: in this
particular data set, males happened to disappear faster than average. The
recovery shows the method is unbiased on average (-0.37 vs -0.40), so the
researcher should quote the interval, not the point estimate.

## (f) Limits from VALIDATION.md that matter at this design

1. **First-interval recruitment is weakly identified** (Known limit 3). Here,
   2015-2016 entrants, entry probability and per-capita recruitment were biased
   +26% to +47% in planted truth. They are not interpreted.
2. **Low power for the sex effect on survival** (27% at logit -0.4 with 8
   occasions and ~260 turtles). A null result would be inconclusive.
3. **Temporary emigration is not modelled** (Known limit 4). If turtles skip
   years (for example, overwintering or nesting away from the ponds), N_super
   is inflated (+22% in the negative control) and capture probability is
   biased low. A robust design (several trapping sessions within each summer)
   would be needed to test for it.
4. **No individual capture heterogeneity** (Known limit 5). Trap-shy or
   trap-happy turtles bias abundance low (-19% to -24% in the verifier's
   negative control). Here the posterior predictive checks and R2ucare 2.CT
   did not flag it.
5. **Apparent survival** includes permanent emigration. Abundance counts the
   turtles using the trapped ponds, not a density.
6. **The `marked` cross-check** cannot use turtles of unknown sex, so its
   survival is pooled over sexes, and its N_t is expected (not realised)
   abundance (Known limit 6). Use it for point-estimate comparison only.

## (g) Reproduce

```
# 1. (optional) regenerate the example data and TRUTH.txt (seed 2024 + set.seed(11))
scripts/ejob -c 1 -m 4 -t 0:10:00 -- scripts/Rscript dev/jolly-seber/make_example.R

# 2. planted-truth recovery at this design (200 replicates, about 5 minutes on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/jolly-seber_si_recover.log -- scripts/Rscript core/run.R recover jolly-seber --config modules/jolly-seber/example/recovery.yaml --nsim 200 --cores 16 --out results/jolly-seber/si/recovery

# 3. fit the example data: report, estimates and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/jolly-seber_si_fit.log -- scripts/Rscript core/run.R fit jolly-seber --config modules/jolly-seber/example/config.yaml --out results/jolly-seber/si/fit

# 4. (optional) redraw figures after editing the config's figures: block
scripts/Rscript core/run.R figures jolly-seber --fit results/jolly-seber/si/fit/jolly-seber_fit.rds --config modules/jolly-seber/example/config.yaml
```

Outputs: `results/jolly-seber/si/recovery/` (`*_summary.csv`, `*_reps.csv`,
HTML report) and `results/jolly-seber/si/fit/` (`jolly-seber_report.html`,
`jolly-seber_estimates.csv`, `figures/`).
