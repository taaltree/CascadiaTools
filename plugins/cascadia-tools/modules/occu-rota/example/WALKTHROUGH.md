# Worked example: do gray foxes avoid coyotes, and does housing density change that?

Module: `occu-rota` (multispecies occupancy with species interactions, Rota et al. 2016).
Appendix S1 and the website's "try an example" walkthrough.

**The data are simulated** (see `README_data.md`), so we know the true answer and can check
the analysis against it. Every number below comes from the runs listed in section (g).

---

## (a) The question, in the researcher's words

> "I ran 220 camera traps for six weeks last summer, in a landscape that goes from forest to
> farms to suburbs. I want to know whether gray foxes are less likely to be at a camera
> where coyotes are, and whether that avoidance is weaker near houses. I know the cameras
> miss animals, and some cameras were down for part of some weeks."

## (b) The files she has

`camera_detections.csv`: one row per camera per week (1308 rows). 1 = photographed that
week, 0 = not photographed, blank = camera not working (`days_active` = 0).

```
"camera_id","week","start_date","days_active","coyote","gray_fox","bobcat"
"CAM-001",1,"2024-06-03",7,1,0,0
"CAM-001",2,"2024-06-10",5,1,0,0
"CAM-001",3,"2024-06-17",7,0,0,0
"CAM-001",4,"2024-06-24",7,0,0,0
```

`camera_sites.csv`: one row per camera (220 rows).

```
"camera_id","forest_pct","housing_per_km2","on_trail","lat","lon"
"CAM-001",70.9,12.1,"no",38.81867,-77.69299
"CAM-002",47.8,26,"yes",38.62226,-77.31816
"CAM-003",51.3,29.3,"no",38.83137,-77.02135
"CAM-004",49.9,6,"yes",38.84026,-77.3148
```

## (c) What she types to Claude

> I have camera-trap data from 220 stations, 6 weekly occasions in June-July 2024. The files
> are camera_detections.csv (one row per camera and week: days the camera worked, and 1/0
> for coyote, gray fox and bobcat) and camera_sites.csv (forest cover within 500 m, housing
> units per km2, whether the camera is on a trail). Are gray foxes less likely to occupy
> sites where coyotes are present, and does that change with housing density? Please
> account for animals the cameras missed and for cameras that were down. Detection probably
> depends on days active and on whether the camera faces a trail. I need results and
> figures for a paper.

## (d) What the agent does, step by step (skill `eco-analyze`)

**1. Project folder.** Raw files are never edited. For this example the files live in
`modules/occu-rota/example/` and the outputs in `results/occu-rota/si/`.

**2. Data intake (`eco-data-intake`).** The agent reads both files and describes them back:

- 220 cameras, 6 weeks. 1308 camera-week rows. 12 rows are missing because 6 cameras were
  pulled after week 4. 81 rows have `days_active` = 0 (camera down) and blank detections.
  So 1227 of 1320 possible camera-weeks (93%) were surveyed.
- Surveyed weeks per camera: 6 weeks at 144 cameras, 5 at 63, 4 at 10, 3 at 2, 2 at 1.
- Camera effort: 1021 surveyed weeks had 7 active days; the rest 2-6 days (mean 6.5, SD 1.2).
- Detections: coyote 242 camera-weeks at 113 cameras (naive occupancy 0.51); gray fox 101
  at 53 cameras (0.24); bobcat 90 at 52 cameras. Coyote and gray fox were both photographed
  at 19 cameras, and in the same week at 13 camera-weeks.
- Forest cover 5.4-94.8%. Housing 2.3-482 units/km2, strongly right-skewed (so the model
  uses its logarithm). 97 cameras on trails, 123 off. No missing covariates.
- Bobcat is not part of her question. The main analysis uses coyote and gray fox;
  `config_3species.yaml` adds bobcat if she wants it later.

**3. Choose the module.** Two species, detection/non-detection with repeat weeks, and a
question about whether one species occurs less often where the other is present, along a
gradient. That is the Rota et al. (2016) multispecies occupancy model: `occu-rota`.

**4. Analysis plan (`eco-analysis-plan`).** Written before any fit:

- Estimands: probability that a gray fox occupies a site where coyotes are present vs
  absent (and the reverse); the coyote-gray fox interaction (log odds ratio of
  co-occurrence) at average conditions; its change with housing density; marginal
  occupancy of each species.
- Occupancy (each species): `~ forest_pct + log(housing_per_km2)`. Interaction:
  `~ log(housing_per_km2)`. Detection (each species): `~ days_active + on_trail`.
  Weeks with 0 active days are "not surveyed".
- Engine: maximum likelihood (unmarked `occuMulti`). Numeric covariates are scaled to
  1 SD after the log, so effects are "per 1 SD".
- Checks: planted-truth recovery at this design; optimiser convergence and no boundary
  estimates; parametric-bootstrap goodness of fit with 200 replicates.
- Caveat stated in advance: co-occurrence is not proof of avoidance or competition.

**5. Planted truth (`eco-planted-truth`).** Before touching her data, the agent simulates
data sets that look like hers and checks that the model gets the right answer.

Run A, the module's standard recovery (`recovery.yaml`): 2 species, 220 sites, 6 surveys,
7% of surveys missing, true values set to the example's truth. The module's simulator has
one site covariate (it stands in for log housing, which carries the interaction slope) and
one survey covariate (it stands in for days active). 200 of 200 fits converged, 0 failed.
**17 of 17 quantities PASS.**

| quantity | truth | mean estimate | relative bias | coverage of 95% CI | power | gate |
|---|---|---|---|---|---|---|
| coyote occupancy | 0.563 | 0.568 | +0.8% | 0.950 | - | PASS |
| gray fox occupancy | 0.258 | 0.258 | +0.1% | 0.950 | - | PASS |
| both present | 0.103 | 0.102 | -0.4% | 0.975 | - | PASS |
| P(gray fox, coyote present) | 0.182 | 0.180 | -1.3% | 0.980 | - | PASS |
| P(gray fox, coyote absent) | 0.354 | 0.362 | +2.1% | 0.960 | - | PASS |
| P(coyote, gray fox present) | 0.399 | 0.395 | -0.9% | 0.995 | - | PASS |
| P(coyote, gray fox absent) | 0.620 | 0.627 | +1.2% | 0.945 | - | PASS |
| interaction (log odds ratio) | -0.90 | -0.97 | -8.0% (stronger) | 0.980 | 0.535 | PASS |
| interaction slope (per SD log housing) | 0.60 | 0.66 | +10.3% (MC SE 0.035; within 3 MC SE) | 0.950 | 0.255 | PASS |
| coyote intercept / housing slope | 0.49 / 0.17 | 0.53 / 0.17 | +8.4% / +0.6% | 0.950 / 0.955 | - / 0.09 | PASS |
| gray fox intercept / housing slope | -0.60 / 0.30 | -0.58 / 0.30 | +3.1% / -0.3% | 0.960 / 0.985 | - / 0.165 | PASS |
| detection per week, coyote / gray fox | 0.25 / 0.29 | 0.248 / 0.289 | -0.7% / -0.3% | 0.930 / 0.950 | - | PASS |
| detection effort slope, coyote / gray fox | 0.31 / 0.31 | 0.319 / 0.316 | +2.9% / +2.0% | 0.945 / 0.935 | 0.915 / 0.695 | PASS |

(Power = share of replicates whose 95% interval excludes 0. The CLI summary prints NA for
power when species are given names in `recovery.yaml`; these values were computed from the
saved replicates, `results/occu-rota/si/recovery/occu-rota_recover_recovery.rds`.)

Run B, the full example model (`dev/occu-rota/si/si_full_recovery.R`): the same generating
code that made her data (three species, bobcat generated and then left out; forest and
log housing on occupancy; days active and trail on detection; cameras down and pulled),
with her 220 cameras held fixed and new random seeds. Each replicate is written as CSVs
and fitted through `read_data()` with her exact `config.yaml` formulas. Two batches of 200
(seeds 2026001-200 and 2027001-200); 400 of 400 fits converged, 0 errors.

| quantity | truth | mean estimate (400 reps) | relative bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| P(gray fox, coyote present) | 0.182 | 0.181 | -0.7% | 0.958 | - | PASS |
| P(gray fox, coyote absent) | 0.354 | 0.359 | +1.4% | 0.965 | - | PASS |
| P(coyote, gray fox present) | 0.399 | 0.398 | -0.2% | 0.960 | - | PASS |
| P(coyote, gray fox absent) | 0.620 | 0.627 | +1.1% | 0.975 | - | PASS |
| coyote / gray fox occupancy | 0.563 / 0.258 | 0.568 / 0.258 | +0.8% / +0.2% | 0.965 / 0.962 | - | PASS |
| interaction (log odds ratio) | -0.90 | -0.955 | -6.1% (stronger) | 0.950 | 0.54 | PASS |
| interaction slope (per SD log housing) | 0.60 | 0.623 | +3.9% | 0.960 | 0.265 | PASS |
| gray fox forest / housing (per SD) | 0.70 / 0.30 | 0.721 / 0.323 | +3.0% / +7.7% | 0.945 / 0.950 | 0.93 / 0.215 | PASS |
| coyote forest / housing (per SD)* | -0.24 / 0.17 | -0.258 / 0.168 | -5.8% / -0.2% | 0.958 / 0.965 | 0.30 / 0.10 | PASS |
| detection days active, coyote / gray fox | 0.30 / 0.30 | 0.302 / 0.327 | +1.0% / +9.4% | 0.958 / 0.945 | 0.86 / 0.54 | PASS |
| detection on trail, coyote / gray fox | 0.80 / -0.20 | 0.784 / -0.225 | -2.0% / -12.6% (within 3 MC SE) | 0.942 / 0.938 | 0.98 / 0.10 | PASS |

\* Approximate truth: once bobcat is summed out, the coyote terms are not exactly linear
(`README_data.md`). All 21 quantities PASS in the pooled 400 replicates. In the first batch
alone, one quantity failed: the gray fox days-active detection slope (+12.0%, just over
both the 10% and 3 MC SE limits); in the second batch it was +6.8% and everything passed.
This is a nuisance detection parameter, not part of her question.

What the agent tells her, in one sentence: *"At your design the model gets the fox-coyote
pattern right (conditional occupancy within 2% of the truth, 95% intervals contain the
truth 95-98% of the time), but there is only about a 54% chance of detecting an avoidance
of the size we simulated, and only about a 27% chance of detecting that it changes with
housing density, so a 'no change' result would not mean no change."*

**6. Fit.** The agent runs `core/run.R fit` with `example/config.yaml` (section g). The fit
takes about 1 second; the 200-replicate goodness-of-fit bootstrap takes a few minutes.

**7. Checks (`eco-convergence`, `eco-model-checking`).**

- The optimiser converged; the Hessian was invertible; no estimate is on a boundary (the
  module flags |estimate| > 8 or SE > 10; the largest logit SE here is 0.50, for the
  interaction). Log-likelihood -825.7, AIC 1679. R-hat and ESS do not apply to a
  likelihood fit.
- Independent route (from `VALIDATION.md`, same data): a hand-written likelihood gives the
  same log-likelihood (-825.6936) and coefficients to within 1.3e-4 SE; Nimble posterior
  means are within 0.22 SE of these estimates.
- Goodness of fit, parametric bootstrap with 200 replicates (chi-square on the number of
  cameras with 0, 1, 2, ... detections of each species, plus co-detections): statistic
  42.1, c-hat 1.12, p = 0.25. No evidence of lack of fit.

**8. Interpretation (`eco-interpretation`).** Excerpt from the report, checked sentence by
sentence against the estimates table:

> "Allowing for missed detections, the estimated probability that each species occupies a
> site with average covariate values is: coyote 0.60 (95% CI 0.50 to 0.69); gray fox 0.25
> (95% CI 0.18 to 0.34)."
>
> "At average covariate values, coyote and gray fox occur together less often than expected
> if they occurred independently (negative co-occurrence; this pattern alone does not show
> that one species avoids the other). Where coyote is present, gray fox occupies a site
> with probability 0.15 (95% CI 0.08 to 0.27); where coyote is absent, 0.41 (95% CI 0.27 to
> 0.57). ... The co-occurrence odds ratio is 0.25 (95% CI 0.09 to 0.66); a value of 1 would
> mean no association."
>
> "The coyote-gray fox co-occurrence odds ratio is multiplied by 2.06 (95% CI 0.85 to 4.96)
> per 1 SD (0.984) increase in log(housing_per_km2), i.e. by 1.66 (95% CI 0.89 to 3.09) per
> doubling of housing_per_km2 (the interval includes 1, so the data do not show that the
> association changes along this gradient). At low housing_per_km2 (10th percentile),
> P(gray fox | coyote present) = 0.04 (95% CI 0.01 to 0.14) vs coyote absent 0.27 (95% CI
> 0.13 to 0.49); at high housing_per_km2 (90th percentile), 0.42 (95% CI 0.25 to 0.61) vs
> 0.55 (95% CI 0.32 to 0.77)."
>
> "Across the observed range of forest_pct (10th percentile 15.8 to 90th percentile 87),
> marginal occupancy changes: coyote from 0.71 (95% CI 0.57 to 0.83) to 0.46 (95% CI 0.34
> to 0.59); gray fox from 0.12 (95% CI 0.06 to 0.22) to 0.46 (95% CI 0.32 to 0.59). This is
> an association, not necessarily a causal effect."
>
> "Per-survey detection probability where the species is present (average survey
> conditions): coyote 0.23 (95% CI 0.18 to 0.29); gray fox 0.30 (95% CI 0.22 to 0.39). With
> 6 surveys, the chance of detecting a species at least once at an occupied site is coyote
> 0.80, gray fox 0.88."
>
> "Caution: co-occurrence is not proof of a species interaction. A positive or negative
> association can arise from shared or opposite responses to habitat that the model does
> not include, ... or from one species affecting the other's detectability rather than its
> occupancy (Blanchet et al. 2020)."

"Average conditions" means forest 52.4%, housing 39.4 units/km2 (geometric mean), 6.5
active days, camera off trail.

In plain words: gray foxes were found much less often where coyotes were (0.15 vs 0.41).
The pattern looks weaker in the suburbs, but this survey cannot confirm that. Foxes were
more common in forest; coyotes in open land.

Figures (all in `results/occu-rota/si/fit/figures/`, PNG and PDF):

| figure | caption (short) |
|---|---|
| `conditional_occupancy.png` | Two panels against housing density (log axis): (A) P(coyote) where gray fox is present vs absent; (B) P(gray fox) where coyote is present vs absent; 95% bands; forest at its mean (52.4%); rug = camera values. **Key figure.** |
| `interaction_housing_per_km2.png` | Co-occurrence odds ratio (log axis) against housing density, 95% band; dashed line = 1 (independent). |
| `conditional_occupancy_average.png` | Marginal and conditional occupancy of each species at average conditions, 95% intervals. |
| `co_occurrence_matrix.png` | Change in each species' occupancy where the other is present vs absent: coyote -0.34 (-0.54 to -0.09), gray fox -0.26 (-0.45 to -0.07). |
| `marginal_occupancy_forest_pct.png` | Marginal occupancy of each species against forest cover, housing at its geometric mean. |
| `marginal_occupancy_housing_per_km2.png` | Marginal occupancy of each species against housing density, forest at its mean. |
| `conditional_gray_fox_given_coyote_housing_per_km2.png` | P(gray fox) where coyote present vs absent, against housing density (panel B of the key figure, alone). |
| `conditional_coyote_given_gray_fox_housing_per_km2.png` | P(coyote) where gray fox present vs absent, against housing density. |
| `conditional_gray_fox_given_coyote_forest_pct.png` | P(gray fox) where coyote present vs absent, against forest cover. |
| `conditional_coyote_given_gray_fox_forest_pct.png` | P(coyote) where gray fox present vs absent, against forest cover. |
| `detection.png` | Detection per week against (A) days active and (B) trail vs off trail, both species. |
| `detection_days_active.png` | Detection per week against days active (off trail). |
| `detection_on_trail.png` | Detection per week on vs off trail (6.5 active days): coyotes are easier, foxes harder, to detect on trails. |
| `coefficients.png` | Odds ratios (log scale) for every effect, grouped into occupancy, interaction and detection; dashed line = no effect. |
| `gof.png` | Goodness of fit: bootstrap distribution of the chi-square statistic; line = observed (42.1, p = 0.25). |

All 15 PNGs were opened and checked: probabilities on a 0-1 axis, housing on a log axis
with readable tick labels, plain-language axis and legend labels (set in the `figures:`
block of `config.yaml`, including short names for the coefficient plot), and the
two-panel key figure has its legends under the panels so the curves get the full width.

**9. Gate (`eco-acceptance-gate`): DONE.** Plan declared before fitting; planted truth
passed at this design (Run A 17/17; Run B 21/21 pooled over 400 replicates); the fit
converged with no boundary estimates; goodness of fit showed no lack of fit; independent
routes (hand-written likelihood, Nimble) agree; the negative control (ignoring detection)
fails the gate as it should (`VALIDATION.md`); every number here was made by a script
listed below.

**10. Deliver.** Summary for the researcher: *"Gray foxes occurred less often where coyotes
were present: at average forest and housing, the chance a fox used a site was 0.15 (95% CI
0.08-0.27) with coyotes vs 0.41 (0.27-0.57) without (co-occurrence odds ratio 0.25,
0.09-0.66). The estimated avoidance weakens toward the suburbs (odds ratio x1.66 per
doubling of housing, 95% CI 0.89-3.09), but the interval includes no change, and at 220
cameras there is only about a 1-in-4 chance of detecting a change of this size. Report the
gradient as 'suggestive, not shown'. This is a co-occurrence pattern; it fits avoidance of
coyotes but could also come from habitat differences the model does not include."*
The report is `results/occu-rota/si/fit/occu-rota_report.html`; its Methods section is the
methods paragraph for the paper.

## (e) Estimates versus the known truth

Because the data were simulated, we can check the answer. **All 21 95% intervals contain
the true value.**

| quantity | truth | estimate (95% CI) | contains truth? |
|---|---|---|---|
| P(gray fox, coyote present) | 0.182 | 0.147 (0.079 to 0.272) | yes |
| P(gray fox, coyote absent) | 0.354 | 0.412 (0.274 to 0.569) | yes |
| P(coyote, gray fox present) | 0.399 | 0.347 (0.184 to 0.553) | yes |
| P(coyote, gray fox absent) | 0.620 | 0.683 (0.558 to 0.782) | yes |
| interaction, log odds ratio (average site) | -0.90 (OR 0.41) | -1.40 (-2.38 to -0.42) (OR 0.25) | yes |
| interaction slope per SD log housing | 0.60 | 0.72 (-0.16 to 1.60) | yes |
| coyote occupancy | 0.563 | 0.597 (0.502 to 0.686) | yes |
| gray fox occupancy | 0.258 | 0.254 (0.183 to 0.343) | yes |
| both present | 0.103 | 0.088 (0.045 to 0.166) | yes |
| gray fox: intercept / forest / log housing | -0.60 / 0.70 / 0.30 | -0.35 (-0.99 to 0.28) / 0.61 (0.20 to 1.03) / 0.48 (-0.08 to 1.03) | yes / yes / yes |
| coyote: intercept / forest / log housing (approx. truth) | 0.49 / -0.24 / 0.17 | 0.77 (0.25 to 1.28) / -0.25 (-0.61 to 0.11) / 0.12 (-0.32 to 0.56) | yes / yes / yes |
| detection per week, coyote / gray fox | 0.250 / 0.294 | 0.234 (0.184 to 0.294) / 0.299 (0.223 to 0.388) | yes / yes |
| detection: days active per SD, coyote / gray fox | 0.31 / 0.31 | 0.40 (0.20 to 0.61) / 0.47 (0.06 to 0.89) | yes / yes |
| detection: on trail, coyote / gray fox | 0.80 / -0.20 | 0.88 (0.50 to 1.26) / -0.58 (-1.14 to -0.03) | yes / yes |

The naive answer (fox at 0.24 of cameras, coyote at 0.51) is below the truth (0.26 and
0.56 in expectation; 67 and 133 cameras truly occupied in this landscape). The avoidance in
this one data set came out stronger than the truth (odds ratio 0.25 vs 0.41). That is
chance in one data set: Run B shows the average estimate is only 6% too strong. The
interaction slope was estimated close to the truth (0.72 vs 0.60), but its interval
includes 0. This is exactly the low-power situation Run A and Run B warned about.

## (f) Limits that matter at this design (from `VALIDATION.md`)

1. **Interaction slopes are inflated at a few hundred sites** (known limit 1): +10-13% at
   400 sites in the validation; +10.3% (Run A) and +3.9% (Run B) here. Trust the sign and
   the conditional-probability curves; treat the size of the slope cautiously.
2. **Low power for the interaction at 220 cameras.** The average-site interaction was
   detected in 54% of replicates and the housing slope in only 26-27%. A non-significant
   slope here is not evidence that avoidance is the same everywhere. Good power on the
   slope needs a larger survey; `eco-power` can find how many cameras (not run here).
3. **Co-occurrence is not interaction** (known limit 7): nothing in the validation tests
   whether a negative association reflects avoidance, predation or opposite habitat
   preferences that the model does not include. Here the forest effects go in opposite
   directions for the two species, which is included; anything unmeasured is not.
4. **Rare species** (known limit 2): gray fox marginal occupancy is about 0.25, above the
   ~0.12 level where the interaction becomes biased at 400 sites. If she adds a rarer
   species (e.g. a three-species model with bobcat), recheck recovery.
5. **Boundary fits** (known limit 3): not a problem here (largest SE 0.50), but a pair that
   is almost never photographed together would push the MLE to a boundary.
6. **Goodness of fit**: the bootstrap test is calibrated for this ML engine but has limited
   power against hidden differences in detection among cameras.
7. Three-way interactions and factor covariates on interactions are not validated (known
   limits 4-5); this example uses neither.

## (g) Reproduce

From `.` (R runs go through `scripts/ejob` on the cluster):

```bash
# 0. (optional) regenerate the example CSVs and derive the true values (checks the CSVs are byte-identical)
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/occu-rota_si_truth.log -- scripts/Rscript dev/occu-rota/si/si_truth.R

# 1. planted truth, Run A (module recovery at the example design)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/occu-rota_si_recover.log -- scripts/Rscript core/run.R recover occu-rota --config modules/occu-rota/example/recovery.yaml --nsim 200 --cores 16 --out results/occu-rota/si/recovery

# 2. planted truth, Run B (full example model, researcher-style CSVs), two batches + pooled summary
scripts/ejob -c 16 -m 32 -t 1:00:00 -l logs/occu-rota_si_fullrecover.log -- scripts/Rscript dev/occu-rota/si/si_full_recovery.R 200 16
scripts/ejob -c 16 -m 32 -t 1:00:00 -l logs/occu-rota_si_fullrecover_s2027000.log -- scripts/Rscript dev/occu-rota/si/si_full_recovery.R 200 16 2027000
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/occu-rota_si_fullrecover_pool.log -- scripts/Rscript dev/occu-rota/si/pool_full.R

# 3. fit the example data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/occu-rota_si_fit.log -- scripts/Rscript core/run.R fit occu-rota --config modules/occu-rota/example/config.yaml --out results/occu-rota/si/fit
```

Outputs: `results/occu-rota/si/recovery/` (Run A), `results/occu-rota/si/recovery_full/`
(Run B), `results/occu-rota/si/fit/` (report, estimates, figures),
`results/occu-rota/si/si_summary.json` (machine-readable summary of this walkthrough).
