# Worked example: does forest cover affect where pileated woodpeckers occur?

Module: `occu-single` (single-species occupancy, static). Appendix S1 and the website's
"try an example" walkthrough.

**The data are simulated** (see `README_data.md`), so we know the true answer and can check
the analysis against it. Every number below comes from the runs listed in section (g).

---

## (a) The question, in the researcher's words

> "I surveyed 120 forest stands for pileated woodpeckers with call playback, up to four
> visits each. I want to know if woodpeckers are more likely to be in stands with more forest
> cover, and whether elevation matters. I know we missed birds sometimes, especially late in
> the season and with some observers."

## (b) The files she has

`detections.csv`: one row per survey done (405 rows).

```
"site","survey","date","julian_day","observer","detected"
"STAND-001",1,"2025-05-07",127,"MR",1
"STAND-001",2,"2025-05-18",138,"AB",1
"STAND-002",1,"2025-05-08",128,"AB",1
"STAND-002",2,"2025-05-19",139,"MR",0
```

`sites.csv`: one row per stand (120 rows).

```
"site","forest_cover_pct","elevation_m","stand_age_class"
"STAND-001",39.8,693,"young"
"STAND-002",83.1,362,"mature"
"STAND-003",51.4,1037,"old"
"STAND-004",58.1,1283,"mature"
```

## (c) What she types to Claude

> I have pileated woodpecker call-playback surveys from 120 forest stands, up to 4 visits
> each in May-June 2025. The files are detections.csv (one row per visit: stand, visit
> number, date, day of year, observer, detected 1/0) and sites.csv (forest cover %,
> elevation in m, stand age). Does forest cover or elevation affect whether woodpeckers
> occupy a stand? Please account for birds we missed, and allow detection to change with
> date and observer. I need results and figures for a paper.

## (d) What the agent does, step by step (skill `eco-analyze`)

**1. Project folder.** The agent keeps the raw files untouched and works in a project
folder. For this example the files live in `modules/occu-single/example/` and the outputs in
`results/occu-single/si/`.

**2. Data intake (`eco-data-intake`).** The agent reads both files and describes them back:

- 120 stands, 405 surveys. 72 stands had 4 visits, 21 had 3, 27 had 2. So 75 of the 480
  possible visits (15.6%) were not done. Missing visits have no row; they are not zeros.
- Woodpeckers were detected on 52 surveys at 34 stands. Naive occupancy is 34/120 = 0.28.
- Surveys ran from day 121 (1 May) to day 170 (19 June). Observers: AB 168 surveys,
  JK 148, MR 89.
- Forest cover 4 to 98%; elevation 185 to 1395 m. No missing covariate values.
- `stand_age_class` is not part of her question, so it is left out.

**3. Choose the module.** Detection/non-detection with repeat visits in one season, and
she asks where the species occurs. That is a static single-species occupancy model
(MacKenzie et al. 2002): `occu-single`.

**4. Analysis plan (`eco-analysis-plan`).** Written before any fit:

- Estimands: occupancy probability; the effects of forest cover and elevation on occupancy;
  the number of the 120 stands that are occupied.
- Occupancy: `~ forest_cover_pct + elevation_m`. Detection per survey:
  `~ julian_day + observer`.
- Engine: maximum likelihood (unmarked). Numeric covariates are scaled to 1 SD, so effects
  are "per 1 SD".
- Checks: planted-truth recovery at this design; optimiser convergence and no boundary
  estimates; MacKenzie-Bailey goodness of fit with 500 bootstrap replicates.

**5. Planted truth (`eco-planted-truth`).** Before touching her data, the agent simulates
200 data sets that look like hers and checks that the model gets the right answer.

Run A, the module's standard recovery (`recovery.yaml`): 120 sites, 4 visits, 15.6% of
visits missing, with true values set to the example's truth. The module's simulator has
one occupancy covariate and one detection covariate, so elevation and observer are not in
this run. 200 of 200 fits converged.

| quantity | truth | mean estimate | relative bias | coverage of 95% CI | power | gate |
|---|---|---|---|---|---|---|
| occupancy probability | 0.406 | 0.402 | -0.9% | 0.945 | - | PASS |
| detection probability | 0.402 | 0.407 | +1.3% | 0.950 | - | PASS |
| forest effect (per SD) | 1.05 | 1.11 | +6.0% | 0.960 | 0.99 | PASS |
| date effect (per SD) | -0.532 | -0.546 | +2.7% (more negative) | 0.960 | 0.865 | PASS |
| proportion of sites occupied | 0.420 | 0.422 | +0.5% | 0.955 | - | PASS |

Run B, the full example model (`dev/occu-single/si_full_recovery.R`): the same generating
code that made her data, with new random seeds, written as CSVs and fitted with her exact
formulas (forest + elevation; date + observer). 200 of 200 fits converged, 0 errors.

| quantity | truth (mean) | mean estimate | relative bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| occupancy probability | 0.398 | 0.406 | +1.9% | 0.965 | - | PASS |
| detection probability | 0.397 | 0.393 | -0.8% | 0.965 | - | PASS |
| forest effect (per SD) | 1.02 | 1.11 | +9.6% | 0.960 | 0.975 | PASS |
| elevation effect (per SD) | -0.513 | -0.544 | +6.0% (more negative) | 0.970 | 0.43 | PASS |
| date effect (per SD) | -0.532 | -0.548 | +2.9% (more negative) | 0.970 | 0.88 | PASS |
| observer JK vs AB | -0.600 | -0.635 | +5.8% | 0.935 | 0.315 | PASS |
| observer MR vs AB | 0.400 | 0.428 | +6.9% | 0.980 | 0.155 | PASS |
| proportion of sites occupied | 0.421 | 0.433 | +2.8% | 0.975 | - | PASS |

What the agent tells her, in one sentence: *"At your design the model recovers the true
values: every quantity passed (bias under 10%, 95% intervals contain the truth 94-98% of the
time), but the forest effect tends to come out about 10% too strong, and there is only about
a 43% chance of detecting an elevation effect of the size we simulated."*

**6. Fit.** The agent runs `core/run.R fit` with `example/config.yaml` (section g).

**7. Checks (`eco-convergence`, `eco-model-checking`).**

- The optimiser converged; the Hessian was invertible; no estimate is on a boundary
  (largest logit SE 0.50). R-hat and ESS do not apply to a likelihood fit.
- Goodness of fit, MacKenzie-Bailey test with 500 bootstrap replicates: chi-square 28.3,
  c-hat 1.31, p = 0.138. No evidence of lack of fit. The test has low power against
  hidden differences in detection among stands (known limit 3), so a pass rules out gross
  misfit, not all bias.

**8. Interpretation (`eco-interpretation`).** Excerpt from the report, checked sentence by
sentence against the estimates table:

> "The species was detected at least once at 34 of 120 sites (naive occupancy 0.28).
> Allowing for missed detections, the estimated probability that a site with average
> covariate values is occupied is 0.33 (95% CI 0.20 to 0.49)."
>
> "Of the 120 surveyed sites, an estimated 47 (95% CI 37 to 58) were actually occupied ...
> this finite-sample estimate applies to these sites only, not to the wider landscape."
>
> "With the typical 4 surveys per site, an occupied site is detected at least once with
> probability 0.87 (95% CI 0.70 to 0.96), so about 13% of occupied sites go undetected
> with this survey effort."
>
> "Occupancy and forest_cover_pct: going from the 10th to the 90th percentile ... (22 to
> 85.3), occupancy changes from 0.09 (95% CI 0.03 to 0.28) to 0.62 (95% CI 0.36 to 0.81)
> ... each 10-unit increase [multiplies the odds] by 1.54."
>
> "Occupancy and elevation_m: going from the 10th to the 90th percentile ... (361 to
> 1250), occupancy changes from 0.69 (95% CI 0.41 to 0.88) to 0.12 (95% CI 0.04 to 0.27)."
>
> "Detection per survey by observer: AB 0.39 (95% CI 0.26 to 0.55); JK 0.13 (95% CI 0.05
> to 0.26); MR 0.52 (95% CI 0.30 to 0.71)."
>
> "These are associations from an observational design; they do not by themselves show
> that a covariate causes the change."

In plain words: woodpeckers were more likely to occupy stands with more forest and at lower
elevation. About 47 of the 120 stands were occupied, not the 34 where birds were heard.
Detection fell through the season and was lowest for observer JK.

Figures (all in `results/occu-single/si/fit/figures/`, PNG and PDF):

| figure | caption (short) |
|---|---|
| `occupancy_forest_cover_pct.png` | Occupancy probability vs forest cover, with 95% band, elevation at its mean (834 m); rug = observed values. **Key figure.** |
| `occupancy_elevation_m.png` | Occupancy probability vs elevation, forest cover at its mean (58.2%). |
| `detection_julian_day.png` | Detection per survey vs survey date, observer AB. An effect on finding birds, not on where they live. |
| `detection_observer.png` | Detection per survey for each observer (points and 95% intervals), date at its mean (day 142). |
| `coefficients.png` | Odds ratios (log scale) per 1 SD for every effect, grouped into occupancy and detection; dashed line = no effect. |
| `detection_by_survey.png` | Average detection per survey by visit number (visits 1-4); falls from about 0.44 to 0.17 because later visits are later in the season. |
| `gof.png` | Goodness of fit: bootstrap distribution of the MacKenzie-Bailey chi-square; line = observed (28.3, p = 0.138). |

All seven figures were opened and checked: labels are in plain units, probabilities are on
a 0-1 axis, bands and rugs are present, and the factor (observer) is drawn as points with
intervals. No changes to the `figures:` block were needed.

**9. Gate (`eco-acceptance-gate`): DONE.** Recovery passed at this design (both runs), the
fit converged with no boundary estimates, goodness of fit showed no lack of fit, the report
sentences match the table, and every number was made by a script listed below.

**10. Deliver.** Summary for the researcher: *"Woodpeckers were more likely to occupy
stands with more forest cover: occupancy rose from 0.09 (95% CI 0.03-0.28) at 22% forest to
0.62 (0.36-0.81) at 85% forest. They were less likely to occupy high stands: 0.69
(0.41-0.88) at 361 m vs 0.12 (0.04-0.27) at 1250 m. About 47 of your 120 stands (95% CI
37-58) were occupied, compared with 34 where you heard birds. Four visits found birds at
87% of occupied stands. Treat the size of the forest effect as possibly a little
overstated (simulation showed about +10% at 120 stands); its direction is solid."*
The report is `results/occu-single/si/fit/occu-single_report.html`; its Methods section is
the methods paragraph for the paper.

## (e) Estimates versus the known truth

Because the data were simulated, we can check the answer. Every 95% interval contains the
true value.

| quantity | truth | estimate (95% CI) | interval contains truth? |
|---|---|---|---|
| occupancy probability (average stand) | 0.406 | 0.332 (0.202 to 0.493) | yes |
| forest effect, logit per SD | 1.05 | 1.03 (0.31 to 1.76) | yes |
| elevation effect, logit per SD | -0.475 | -1.07 (-1.72 to -0.42) | yes (near the edge) |
| detection per survey (day 142, AB) | 0.402 | 0.394 (0.262 to 0.543) | yes |
| date effect, logit per SD | -0.532 | -0.49 (-0.90 to -0.09) | yes |
| observer JK vs AB | -0.6 | -1.49 (-2.48 to -0.51) | yes |
| observer MR vs AB | 0.4 | 0.50 (-0.43 to 1.44) | yes |
| proportion of stands occupied | 0.375 (45 stands) | 0.389, 47 stands (0.308 to 0.483; 37 to 58) | yes |

The naive answer (0.28, 34 stands) is well below the truth (45 stands). The elevation and
JK effects came out about twice as strong as the truth. That is chance in one data set: in
Run B the average elevation estimate was only 6% too strong. Wide intervals are the honest
signal that 120 stands cannot pin these effects down precisely.

## (f) Limits that matter at this design (from `VALIDATION.md`)

1. **Small samples inflate occupancy slopes** (known limit 1). Below about 150 sites the
   forest effect is biased away from zero (+11.5% at 100 sites, failing the gate). At 120
   stands it was +6.0% (Run A) and +9.6% (Run B), just inside the 10% limit. Trust the
   direction; treat the size as possibly a little high.
2. **Elevation and observer effects have low power here.** Run B detected a true elevation
   effect of this size only 43% of the time, and the observer differences 16-32% of the
   time. A "no effect" result for these would not have meant no effect.
3. **Hidden differences in detection among stands bias occupancy low** (known limit 3),
   about -23% at a logit SD of 1.5, and the goodness-of-fit test catches this only about 1
   time in 10. The model includes date and observer to reduce this risk.
4. **Few surveys with low detection fail** (known limit 2). Not a problem here: cumulative
   detection over 4 visits is 0.87, well above the 0.5 warning level. But 27 stands had
   only 2 visits.
5. No false-positive detections, no spatial autocorrelation, and the model assumes a stand's
   occupancy did not change between 1 May and 19 June (closure).

## (g) Reproduce

From `.` (R runs go through `scripts/ejob` on the cluster):

```bash
# 0. (optional) regenerate the example CSVs and the true values
scripts/Rscript dev/occu-single/make_example.R
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/occu-single_si_truth.log -- scripts/Rscript dev/occu-single/si_truth.R

# 1. planted truth, Run A (module recovery at the example design)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/occu-single_si_recover.log -- scripts/Rscript core/run.R recover occu-single --config modules/occu-single/example/recovery.yaml --nsim 200 --cores 16 --out results/occu-single/si/recovery

# 2. planted truth, Run B (full example model, researcher-style CSVs)
scripts/ejob -c 16 -m 32 -t 1:00:00 -l logs/occu-single_si_fullrecover.log -- scripts/Rscript dev/occu-single/si_full_recovery.R 200 16

# 3. fit the example data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/occu-single_si_fit.log -- scripts/Rscript core/run.R fit occu-single --config modules/occu-single/example/config.yaml --out results/occu-single/si/fit
```

Outputs: `results/occu-single/si/recovery/` (Run A), `results/occu-single/si/recovery_full/`
(Run B), `results/occu-single/si/fit/` (report, estimates, figures),
`results/occu-single/si/si_summary.json` (machine-readable summary of this walkthrough).
