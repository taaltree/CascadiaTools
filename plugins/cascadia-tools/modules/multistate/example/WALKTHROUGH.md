# Worked example: chytrid infection, survival and clearance in marked frogs (module `multistate`)

This walkthrough shows what happens when a researcher who does not code brings
capture-recapture data with infection status to the agent. Every number below
comes from the real run in `results/multistate/si/`. **The data are simulated**
with known true values (see `README_data.md`), so you can check the answers
against the truth in section (e).

## (a) The question, in the researcher's words

"We PIT-tag and swab frogs at our breeding pond every year. Some frogs are
infected with chytrid (Bd) and some are not, and frogs can lose or pick up the
infection between years. I want to know three things. Do infected frogs survive
worse? How often do frogs get infected, and how often do they clear it? Does
water temperature affect infection or clearance? Some swabs were lost or gave no
result, so for those captures we don't know the frog's status."

## (b) The files the researcher has

Three CSV files in `modules/multistate/example/`.

`frog_captures.csv`: one row per frog per year it was caught (516 rows). `U` = swab lost or inconclusive.
```
"frog_id","year","bd_status"
"RL-0001","2015","U"
"RL-0001","2016","uninfected"
"RL-0001","2021","infected"
"RL-0002","2015","infected"
```

`frog_individuals.csv`: one row per frog (280 rows), body length at first capture.
```
"frog_id","svl_mm"
"RL-0001",57.9
"RL-0002",54.7
"RL-0003",51.2
"RL-0004",50.6
```

`survey_years.csv`: one row per survey year (8 rows).
```
"year","water_temp_C","survey_hours"
"2015",18.26,6.49
"2016",16.97,8.54
"2017",20.38,3.94
"2018",16.92,6.07
```

## (c) What the researcher types to Claude

> I have eight years (2015-2022) of capture data for a PIT-tagged frog
> population in three files: frog_captures.csv (each capture with the frog's id,
> year and Bd swab result: uninfected, infected, or U when the swab was lost or
> inconclusive), frog_individuals.csv (snout-vent length of each frog) and
> survey_years.csv (mean water temperature and survey effort in person-hours for
> each year). I want to know whether infected frogs survive less well, how often
> frogs become infected or clear the infection each year, whether water
> temperature is related to infection and clearance, and whether bigger frogs
> survive better. Please analyse this and give me a report with figures I can
> use in a paper.

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

### Step 1. Project folder
The agent makes `projects/<name>/` with `data_raw/` (untouched copies), `data/`,
`analysis_plan.md`, `config.yaml`, `recovery/`, `results/` and `log.md`. For this
SI example the files live in `modules/multistate/example/` and the outputs in
`results/multistate/si/`.

### Step 2. Data intake (`eco-data-intake`)
The agent reads the files and describes them back:

- 280 frogs, each tagged in one year from 2015 to 2021 (40 new frogs per year). No new frogs in 2022.
- 516 captures over 8 annual surveys: 258 uninfected, 144 infected and **114 unknown (22%)**.
- 141 of 280 frogs (50%) were recaptured in a later year. Captures per frog range from 1 to 6.
- Snout-vent length: mean 52.0 mm (SD 6.0, range 31.8-66.5). Every frog has a value.
- Water temperature: 13.9-20.4 deg C (SD 2.0). Effort: 3.9-8.5 person-hours.
- No duplicate ids. No conflicting statuses within a year. Every year has a covariate row.
- `check_data()` reports no problems.

The agent explains one point in plain terms: the 2022 water temperature cannot
affect anything, because there is no later survey to show what happened next.

### Step 3. Choose the module (`catalog.md`)
The agent picks `multistate` because each capture records a **state** (infection
status), and the question is about survival in each state and switching between
states. `cjs` is ruled out because it has no states. Because some captures have
an unknown state, the module runs the **multievent** version (Pradel 2005;
`model: auto`). That version uses those captures instead of dropping them or
guessing their status.

### Step 4. Analysis plan (`eco-analysis-plan`)
Plan summary, confirmed with the researcher:

- **Estimands:**
  - annual survival of uninfected and infected frogs, and the difference;
  - annual infection probability (uninfected to infected, among survivors);
  - annual clearance probability (infected to uninfected, among survivors);
  - the effect of water temperature on infection and on clearance;
  - the effect of body length on survival.
- **Nuisance parameters:**
  - detection per state, with survey effort as a covariate;
  - the probability that a swab gives a result ("assignment");
  - infection status at first capture.
- **Model:** multievent Arnason-Schwarz, fitted in NIMBLE with the forward algorithm.
  - Priors: logistic(0,1) on intercepts (uniform on probabilities); Normal(0, SD 2) on slopes per SD of the covariate.
  - MCMC: 3 chains x 6000 iterations, 2000 burn-in.
- **Checks:**
  - R-hat < 1.01 and bulk ESS >= 400;
  - three posterior predictive checks;
  - the U-CARE tests.
- **Caveat agreed up front:** the temperature effects rest on only 7 yearly values. They are associations, not causes.

### Step 5. Planted truth at this design (`eco-planted-truth`)
The agent simulates 100 datasets with the module's own simulator. Each matches
this study: 2 states, 8 annual surveys, 40 new frogs per year, multievent, and
the same three covariates. The true values are realistic (`recovery.yaml`). The
agent fits each dataset (3 chains x 2000 iterations; median 5.5 min per fit,
16 in parallel) and applies the acceptance gate:
|bias| <= max(10% of truth, 3 MC SE), coverage >= 0.885, and >= 90% converged.

**Result: PASS 16/16. All 100 replicates converged.**

| estimand | truth | mean estimate | rel. bias | 95% CI coverage | power | gate |
|---|---|---|---|---|---|---|
| survival_uninfected | 0.75 | 0.759 | +1.2% | 0.93 | | PASS |
| survival_infected | 0.55 | 0.551 | +0.1% | 0.94 | | PASS |
| transition_uninfected_to_infected | 0.30 | 0.321 | +7.1% | 0.94 | | PASS |
| transition_infected_to_uninfected | 0.35 | 0.347 | -0.9% | 0.97 | | PASS |
| detection_uninfected | 0.60 | 0.615 | +2.5% | 0.93 | | PASS |
| detection_infected | 0.50 | 0.521 | +4.2% | 0.96 | | PASS |
| assignment_uninfected | 0.75 | 0.758 | +1.1% | 0.94 | | PASS |
| assignment_infected | 0.85 | 0.836 | -1.7% | 0.95 | | PASS |
| initial_uninfected | 0.65 | 0.641 | -1.4% | 0.93 | | PASS |
| initial_infected | 0.35 | 0.359 | +2.5% | 0.93 | | PASS |
| stationary_uninfected | 0.538 | 0.517 | -4.1% | 0.96 | | PASS |
| stationary_infected | 0.462 | 0.483 | +4.7% | 0.96 | | PASS |
| survival_beta_svl_mm | +0.4 | 0.436 | +9.1% | 0.91 | 0.97 | PASS |
| transition_beta_water_temp_C_uninfected_to_infected | -0.7 | -0.674 | -3.7% (toward 0) | 0.94 | 0.76 | PASS |
| transition_beta_water_temp_C_infected_to_uninfected | +0.5 | 0.548 | +9.5% | 0.97 | 0.26 | PASS |
| detection_beta_survey_hours | +0.4 | 0.425 | +6.2% | 0.96 | 0.86 | PASS |

"Power" is the share of replicates whose 95% interval excludes 0. The planted-truth
step gives the researcher a warning **before** they see their own results. At this
design a temperature effect on **clearance** of this size (+0.5 per SD) is found
only 26% of the time. An interval that excludes 0 will usually come with an
estimate larger than the true effect (see step 8).

In one sentence to the researcher: *"With 280 frogs over 8 years, the model
recovers survival, infection and clearance without bias. It reliably detects a
body-size effect on survival, and usually detects a moderate temperature effect
on infection. It is weak for a temperature effect on clearance."*

### Step 6. Fit
The agent writes `config.yaml`: files, column names, states, the `U` code, which
covariate goes on which process, MCMC length and plain-language figure labels.
It then runs the fit as a SLURM job (command in section (g)). The fit takes
about 5 minutes.

### Step 7. Convergence and goodness of fit (`eco-convergence`, `eco-model-checking`)
- **Convergence:** max R-hat 1.003 and min bulk ESS 788 (the infection probability), so the real-data standard is met. The trace plot shows well-mixed chains.
- **Posterior predictive checks** (300 draws), with misfit flagged below 0.05 or above 0.95:
  - sightings per year and status (Freeman-Tukey): p = 0.43;
  - recorded status switches: observed 49 vs replicated median 53, p = 0.64;
  - frogs never seen again: observed 139 vs replicated median 141, p = 0.60.
- **U-CARE overall test** (Arnason-Schwarz JMV, unknown sightings treated as not seen, so approximate): chi-square 24.1 on 52 df, p = 1.00.
- No sign of transients, trap dependence or memory in switching.

### Step 8. Interpretation (`eco-interpretation`)
The report's "What this means biologically" section, quoted (agent verified each number against the estimates table):

> Survival over one year for an animal with average covariate values was 0.84
> (95% credible interval 0.75-0.96) in state 'uninfected'; 0.50 (95% credible
> interval 0.37-0.63) in state 'infected'.
>
> The difference in survival (uninfected minus infected) was 0.34 (95% interval
> 0.18 to 0.50); posterior probability that survival is higher in 'uninfected' = 1.00.
>
> Transitions ... from 'uninfected' to 'infected' 0.33 (95% credible interval
> 0.19-0.52); from 'infected' to 'uninfected' 0.34 (95% credible interval 0.18-0.54).
>
> Expected number of years a surviving animal stays in a state before switching
> ...: 'uninfected' 3.1 (1.9-5.4); 'infected' 3.0 (1.8-5.5).
>
> Survival and svl_mm: odds ratio per 1 SD (= 6) 1.77 (95% interval 1.34-2.50);
> P(effect > 0) = 1.00. Going from the 10th to the 90th percentile of svl_mm
> (44.5 to 59), predicted survival changes in 'uninfected' from 0.72 to 0.91, in
> 'infected' from 0.33 to 0.65.
>
> Transitions and water_temp_C (from 14.6 to 19.1, 10th to 90th percentile):
> 'uninfected' to 'infected': 0.48 at the 10th percentile vs 0.23 at the 90th
> (difference -0.24, 95% interval -0.41 to -0.06 ...); 'infected' to
> 'uninfected': 0.16 at the 10th percentile vs 0.56 at the 90th (difference 0.40,
> 95% interval 0.10 to 0.69 ...). Associations only.
>
> Detection and survey_hours: odds ratio per 1 SD (= 1.5) 1.54 (95% interval
> 1.21-2.01). This is about how easily animals were caught or seen, not about
> their survival or movement.

The agent adds context for the researcher:

- Infected frogs survive about 0.34 less per year. Infection itself is not "sticky": about a third of survivors switch each year in either direction.
- Cooler years go with more new infections and less clearance. That fits Bd growing best in cool water. It rests on 7 years, so it is an association, and other year-to-year differences could produce it.
- The clearance-temperature estimate (+0.92 per SD, 95% CrI 0.21-1.74) is the least reliable number. At this design its power is only 0.26. When such an effect is detected, the estimate tends to be too large.
- Detection is lower for infected frogs (0.39 vs 0.61). The model accounts for this. It means infected frogs are harder to catch, not that they are rarer.
- The stationary split (0.50 / 0.50) describes switching only. It is not a forecast of prevalence.

**Figures** (in `results/multistate/si/fit/figures/`, PNG + PDF):

| figure | caption (short) |
|---|---|
| `vital_rates_by_state` | Survival, detection, swab-result probability and share at first capture for uninfected and infected frogs (mean, 95% CrI). |
| `transitions` | Heatmap of the annual transition matrix for survivors (mean, 95% CrI). |
| `transition_estimates` | Infection and clearance probabilities per year (mean, 95% CrI). |
| `stationary` | Long-run share of time in each state implied by the transitions (not a prevalence forecast). |
| `survival_svl_mm` | Annual survival against snout-vent length, one curve per state, 50% and 95% bands, rug of frogs. |
| `transition_water_temp_C` | Infection and clearance probability against mean water temperature, 50% and 95% bands, rug of years. |
| `detection_survey_hours` | Detection against survey effort, one curve per state. |
| `detection_by_occasion` | Detection in each year and state. |
| `coefficients` | All covariate effects as odds ratios per 1 SD (log scale). |
| `captures_by_occasion` | Raw data: frogs caught per year by recorded status (including unknown). |
| `gof` | Three posterior predictive checks (p = 0.43, 0.64, 0.60). |
| `trace` | MCMC traces of survival and transition probabilities. |
| `posterior` | Posterior densities of the 14 reported quantities other than the stationary distribution. |

The most informative figure for the question is `transition_water_temp_C`. It
shows infection falling and clearance rising with water temperature, and how
wide the uncertainty is.

The agent's figure review (presentation only) found one fix. Panel titles in
`trace` and `posterior` showed raw parameter names ("transition beta water temp
C ..."). The agent added plain-language labels in the `figures: labels:` block
of `config.yaml` and redrew the figures. One cosmetic point was left as is:
`transition_water_temp_C` colours its two curves by transition, not by state,
so orange there means "infected to uninfected". The legend states this.

### Step 9. Acceptance gate (`eco-acceptance-gate`)
| criterion | result |
|---|---|
| planted truth at this design | PASS 16/16 (100/100 converged) |
| convergence | max R-hat 1.003, min ESS 788 (meets R-hat < 1.01, ESS >= 400) |
| goodness of fit | PPC p = 0.43 / 0.64 / 0.60; U-CARE p = 1.00 |
| interpretation checked against numbers | yes; associational wording; detection kept separate from ecology |
| limits stated | yes (section f) |

**Verdict: DONE.** It carries one stated caveat: the clearance-temperature effect is
low-power and its size is probably overstated.

### Step 10. Deliver
Summary for the researcher:

> Infected frogs had lower annual survival than uninfected frogs: 0.50 (95% CrI
> 0.37-0.63) vs 0.84 (0.75-0.96), a difference of 0.34 (0.18-0.50). Among
> survivors, about a third of uninfected frogs became infected each year (0.33,
> 0.19-0.52), and a similar share of infected frogs cleared the infection (0.34,
> 0.18-0.54). Cooler years were associated with more new infections: 0.48 at
> 14.6 deg C vs 0.23 at 19.1 deg C. They were also associated with less
> clearance, but that second association is weakly supported at this design and
> its size is probably overstated. Larger frogs survived better (odds ratio 1.77
> per 6 mm, 1.34-2.50). The model fits the data well. The temperature results
> are associations based on 7 years.

The report (`multistate_report.html`) includes a methods paragraph and references.

## (e) Estimates versus the known truth

| parameter | truth | estimate | 95% CrI | truth inside? |
|---|---|---|---|---|
| survival_uninfected | 0.75 | 0.838 | 0.747-0.957 | yes (near the lower edge) |
| survival_infected | 0.55 | 0.498 | 0.370-0.627 | yes |
| transition_uninfected_to_infected | 0.30 | 0.334 | 0.186-0.517 | yes |
| transition_infected_to_uninfected | 0.35 | 0.338 | 0.182-0.541 | yes |
| detection_uninfected | 0.60 | 0.608 | 0.491-0.745 | yes |
| detection_infected | 0.50 | 0.389 | 0.229-0.633 | yes |
| assignment_uninfected | 0.75 | 0.737 | 0.673-0.806 | yes |
| assignment_infected | 0.85 | 0.865 | 0.734-0.983 | yes |
| initial_uninfected | 0.65 | 0.677 | 0.599-0.744 | yes |
| initial_infected | 0.35 | 0.323 | 0.256-0.401 | yes |
| stationary_uninfected | 0.538 | 0.502 | 0.310-0.709 | yes |
| stationary_infected | 0.462 | 0.498 | 0.291-0.690 | yes |
| survival_beta_svl_mm | +0.4 | 0.569 | 0.291-0.917 | yes |
| transition_beta_water_temp_C_uninfected_to_infected | -0.7 | -0.510 | -0.936 to -0.124 | yes |
| transition_beta_water_temp_C_infected_to_uninfected | +0.5 | 0.922 | 0.206-1.741 | yes (overestimated) |
| detection_beta_survey_hours | +0.4 | 0.434 | 0.187-0.696 | yes |

All 16 true values lie inside their 95% credible intervals. The clearance slope
is overestimated (0.92 vs 0.50), as the low power in step 5 predicted. With only
7 temperature values, a detected weak effect tends to be inflated.

## (f) Limits from `VALIDATION.md` that matter at this design

1. **Transition slopes on a yearly covariate are weak** (Known limit 2). With 8 surveys there are only 7 temperature values.
   - Power here was 0.76 for infection and 0.26 for clearance.
   - Slopes tend to be estimated away from zero.
   - Treat them as associations; any other year effect is confounded with temperature.
2. **False positives for transition slopes** run slightly above nominal (0.10 for one slope in the null run; Known limit 7). A lone transition "effect" whose interval barely excludes 0 needs caution.
3. **Misclassification is not modelled** (Known limit 4). A false-negative qPCR result recorded as "uninfected" would inflate the infection and clearance probabilities. Here the model assumes every recorded result is correct and only lost/inconclusive swabs are uncertain.
4. **Assignment is constant over time** (Assumption 4). If swab failures changed between years (a new lab, for example), that is not captured.
5. **Survival is apparent survival.** Death and permanent emigration from the pond are confounded.
6. **The U-CARE tests are approximate for multievent data** (Known limit 6). Unknown sightings are treated as not seen.
7. **Body length is measured once, at first capture.** It is treated as a fixed individual covariate, and its slope is shared by both states.
8. Known limits 1 and 3 (prior pull on rare transitions and under high state uncertainty) **do not apply** here: transitions are about 0.3, and 74-86% of swabs give a result.

## (g) Commands to reproduce (run from `ecotools/`)

```bash
# 0. regenerate the example data (seed 2024) - optional, the CSVs are already in the folder
scripts/ejob -c 1 -m 4 -- scripts/Rscript dev/multistate/make_example.R
# 1. planted-truth recovery at this design (100 replicates; ~5.5 min per fit, ~55 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/multistate_si_recover.log -- scripts/Rscript core/run.R recover multistate --config modules/multistate/example/recovery.yaml --nsim 100 --cores 16 --out results/multistate/si/recovery
# 2. fit the example data: report + figures (~5 min)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/multistate_si_fit.log -- scripts/Rscript core/run.R fit multistate --config modules/multistate/example/config.yaml --out results/multistate/si/fit
# 3. (optional) redraw figures from the saved fit after changing the figures: block
scripts/ejob -c 1 -m 8 -t 0:30:00 -- scripts/Rscript core/run.R figures multistate --fit results/multistate/si/fit/multistate_fit.rds --config modules/multistate/example/config.yaml --out results/multistate/si/fit
```

Outputs:
- `results/multistate/si/recovery/`: gate table (`multistate_recover_summary.csv`), replicate estimates and HTML report.
- `results/multistate/si/fit/`: `multistate_report.html`, `multistate_estimates.csv` and `figures/`.
- `results/multistate/si/si_summary.json`: machine-readable summary.
