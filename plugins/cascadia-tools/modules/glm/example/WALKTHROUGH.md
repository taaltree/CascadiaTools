# Worked example (Appendix S1): Wood Thrush counts with the `glm` module

This walkthrough shows what happens when a researcher who does not code brings
count data and a question to the agent. Every number below comes from the runs
listed in section (g). The data are **simulated** (see `README_data.md`), so we
can compare the estimates with the true values.

## (a) The question, in the researcher's words

> "We did point counts for Wood Thrush at 45 stations, four visits each. Are
> there more thrushes where the canopy is denser? Does elevation or forest type
> matter? Some counts were only 5 minutes instead of 10."

## (b) The files they have

One CSV, `wood_thrush_counts.csv` (180 rows, one per station visit):

```
"station","visit","date","observer","survey_min","forest_type","canopy_cover_pct","elevation_m","wood_thrush"
"WT-01",1,"2025-05-20","JL",5,"mixed",66,725,0
"WT-01",2,"2025-05-29","JL",10,"mixed",66,725,3
"WT-01",3,"2025-06-09","MR",10,"mixed",66,725,1
"WT-01",4,"2025-06-18","JL",10,"mixed",66,725,1
```

Nothing else is needed. The agent writes the config (`config.yaml`) for them.

## (c) What they type to Claude

> I have Wood Thrush point counts in wood_thrush_counts.csv: 45 stations, each
> counted 4 times in May-June 2025. Most counts were 10 minutes, some were 5.
> I want to know whether the number of thrushes changes with canopy cover,
> elevation and forest type. Please analyse it and give me a report with
> figures I can use in a paper.

## (d) What the agent does, step by step (skill `eco-analyze`)

**1. Project folder.** It copies the CSV unchanged into `data_raw/` and keeps a
log of every command.

**2. Data intake** (`eco-data-intake`). The agent describes the data back:

* 180 counts from 45 stations, exactly 4 visits each; no missing values.
* Counts range 0-10, mean 1.13; half of all counts are 0.
* Effort differs: 118 ten-minute and 62 five-minute counts. So survey length
  must enter as effort (an offset), not be ignored.
* Canopy cover (42-98 %), elevation (184-941 m) and forest type (15 stations
  each of deciduous, mixed, conifer) are measured once per station.
* `date`, `visit` and `observer` are present but not part of the question.

**3. Choose the module.** Counts on repeated visits to stations, with
habitat predictors, points to `glm`. The agent picks a negative-binomial GLMM
(`nbinom2`), because bird counts usually vary more than Poisson allows. It adds a
random intercept for station, because four visits to one station are not
independent. It notes the alternative: an N-mixture model would separate
abundance from detection, but this question is about counts, and the report
says so (see the detection caveat below).

**4. Analysis plan** (`eco-analysis-plan`), written before fitting:

* Estimands: rate ratios for canopy cover (per 10 percentage points), elevation
  (per 100 m), and mixed and conifer forest versus deciduous.
* Model: `wood_thrush ~ canopy_cover_pct + elevation_m + forest_type + (1|station)`,
  NB2, log link, offset log(survey minutes), results per 10-minute count.
* Station SD and NB dispersion are nuisance parameters, needed for honest
  intervals but not part of the question.
* Adequate fit: optimiser converged, Hessian positive definite, all DHARMa tests p > 0.05.

**5. Planted truth** (`eco-planted-truth`). Before touching the real counts,
the agent simulates 200 datasets at this exact design and checks that the model
recovers the true values. The design is 45 stations x 4 visits, station-level
canopy, elevation and forest type, and 5-10 min effort. The values are canopy
+0.39 per SD, elevation -0.28 per SD, mixed -0.25, conifer -0.90, station SD
0.45 and theta 3 (`recovery.yaml`). A selftest showed each fit takes about
0.6 s, so all 200 replicates ran in about 3 minutes on 16 cores.
199 of 200 fits converged.

| estimand | truth | mean estimate | rel. bias | 95% coverage | gate |
|---|---|---|---|---|---|
| canopy (log scale, per SD) | 0.390 | 0.386 | -0.9% | 0.940 | PASS |
| elevation (per SD) | -0.280 | -0.275 | -1.8% | 0.905 | PASS |
| mixed vs deciduous | -0.250 | -0.245 | -2.0% | 0.930 | PASS |
| conifer vs deciduous | -0.900 | -0.890 | -1.1% | 0.940 | PASS |
| intercept | -2.080 | -2.101 | +1.0% | 0.915 | PASS |
| rate ratio canopy (per SD) | 1.477 | 1.484 | +0.4% | 0.940 | PASS |
| rate ratio elevation (per SD) | 0.756 | 0.767 | +1.5% | 0.905 | PASS |
| rate ratio mixed | 0.779 | 0.817 | +4.9% | 0.930 | PASS |
| rate ratio conifer | 0.407 | 0.436 | +7.2% | 0.940 | PASS |
| mean count at reference (per min) | 0.125 | 0.125 | +0.2% | 0.915 | PASS |
| station SD | 0.450 | 0.371 | **-17.6%** | 0.940 | **FAIL** (bias) |
| NB theta | 3 | undefined (7 fits at theta > 100) | median +7% | - | not gradable |

(Coverage bound for 200 replicates: 0.904. Bias bound: max(10% of truth, 3 Monte Carlo SE).)

What the agent tells the researcher: *"At your design, the model recovers
every effect you asked about: bias of 7% or less and coverage of 0.90-0.94.
It does not recover the spread among stations. In 45 of 199 simulated datasets
the station SD collapsed to zero, because with four visits and about one bird
per count, extra-Poisson noise and station differences are hard to tell apart.
So I will not interpret the station SD or the dispersion. Both only make the
intervals honest."*

**6. Fit.** `core/run.R fit glm --config config.yaml` (via `scripts/ejob`).

**7. Checks** (`eco-convergence`, `eco-model-checking`):

* Converged (glmmTMB, REML for the variance parameters); 180 observations, 45 stations. No warnings.
* DHARMa residuals (250 simulations): uniformity p = 0.67; dispersion ratio 1.1 (p = 0.52);
  outliers p = 1; zero-inflation p = 0.98. No sign of a poor fit or excess zeros.
* R2 (Nakagawa): marginal 0.165 (predictors), conditional 0.275 (plus stations).

**8. Interpretation.** Excerpts quoted from `glm_report.html`:

> "Per 10 percentage points of canopy cover (%) this is associated with a 39%
> higher expected Wood Thrush count (rate ratio 1.39, 95% CI 1.18 to 1.64).
> Across the observed range of canopy cover (%) (42 to 98) ... the expected
> Wood Thrush count goes from 0.432 to 2.78."

> "Each 1-SD increase in elevation (m) (1 SD = 234) is associated with a 12% lower
> expected Wood Thrush count (rate ratio 0.884, 95% CI 0.692 to 1.13). Per 100 m of
> elevation (m) this is associated with a 5.1% lower expected Wood Thrush count
> (rate ratio 0.949, 95% CI 0.854 to 1.05). ... The 95% interval includes no
> effect: the data are compatible with effects from 31% lower to 13% higher, including none."

> "Compared with forest type = deciduous (reference), forest type = conifer is
> associated with a 52% lower expected Wood Thrush count (rate ratio 0.482, 95% CI
> 0.248 to 0.935): 0.528 vs 1.1."

> "These are effects on the expected count as recorded, which combines how
> many/whether animals are present with how likely they are to be detected or
> recorded. If detectability changes with a predictor (e.g. dense canopy hides
> birds ...), its effect on counts is not purely an effect on abundance."

The agent checked every sentence against `glm_estimates.csv`.

Figures (`results/glm/si/fit/figures/`, PDF + PNG, ESA single-column size):

| id | caption (short) |
|---|---|
| `effect_canopy_cover_pct` | Expected Wood Thrush count per 10-min survey vs canopy cover, 95% band; elevation at its mean, deciduous forest, typical station; rug = observed values. |
| `effect_elevation_m` | Same, vs elevation. |
| `effect_forest_type` | Expected count per 10-min survey in each forest type, with 95% intervals. |
| `effects_panel` | Panels A-C combining the three effect plots. |
| `coefficients` | Rate ratios with 95% CIs on a log axis (per SD for continuous predictors, vs deciduous for forest types); dashed line = no effect. |
| `residuals` | DHARMa QQ plot (KS p = 0.67) and residuals vs prediction (dispersion p = 0.52). |
| `random_effects` | Station intercepts as multipliers vs a typical station, with 95% intervals. |

The most informative single figure is `effects_panel`: it answers the
question on the scale the researcher counts on (birds per 10-minute count).

**9. Gate** (`eco-acceptance-gate`). Verdict **DONE** for the question asked:
planted truth passes for every interpreted estimand at this design. The fit
converged and DHARMa found no problem. The module was cross-checked against
lme4/MASS and passed negative controls (`VALIDATION.md` sections 4-5). The run
is reproducible from the commands below. The station SD and theta are
explicitly not interpreted because they failed or could not be graded.

**10. Deliver.** Summary for the researcher: *"Wood Thrush counts rose with
canopy cover. Each 10 more percentage points of canopy meant 39% more birds per
10-minute count (95% CI 18-64%). Counts in conifer stands were about half those in
deciduous stands (rate ratio 0.48, 0.25-0.94). Mixed stands did not differ
clearly from deciduous (1.15, 0.65-2.01). Elevation had no clear effect
(0.95 per 100 m, 0.85-1.05). These are counts, not abundance corrected for
detection. If dense canopy changes how easily thrushes are heard, part of the
canopy effect could be detection."*

## (e) Estimates versus the known truth

| quantity | estimate (95% CI) | true value | inside CI? |
|---|---|---|---|
| canopy, rate ratio per 10 points | 1.39 (1.18-1.64) | 1.28 | yes |
| elevation, rate ratio per 100 m | 0.949 (0.854-1.054) | 0.887 | yes |
| conifer vs deciduous | 0.48 (0.25-0.94) | 0.41 | yes |
| mixed vs deciduous | 1.15 (0.65-2.01) | 0.78 | yes |
| expected count per 10 min, deciduous, at mean canopy/elevation | 1.10 (0.73-1.65) | 1.25 | yes |
| station SD (log scale) | 0.56 (0.34-0.93) | 0.45 | yes (not interpreted) |
| NB theta | 5.9 (1.1-30.8) | 3 | yes (not interpreted) |

All seven true values lie inside their intervals. Mixed forest is the weakest
case: the truth is 22% fewer birds than deciduous, but this dataset cannot
tell that apart from no difference. Recovery showed this contrast has wide
intervals (mean CI width 1.12 on the log scale).

## (f) Limits that matter at this design (from `VALIDATION.md` and this run)

* **Station SD and dispersion are poorly identified with 4 visits and low counts.** At this
  design the station SD is biased -18% in the mean, and 23% of simulated fits put it at zero.
  NB theta is sometimes unbounded. Do not interpret either. (This finding comes from this
  example's planted-truth run; it is consistent with Known limit 2 on few observations per group.)
* **Rate ratios are biased upward in the mean when the contrast is imprecise** (Known limit 1).
  Here the conifer ratio was +7.2% and the mixed ratio +4.9%. Both are within the gate, and the
  log-scale coefficients are nearly unbiased. Report ratios with their CIs.
* **Wald intervals only** (Known limit 4). Coverage was 0.90-0.94 in recovery.
* **The expected count is for a typical station** (random effect 0), not the landscape average (Known limit 5).
* **Counts mix abundance and detection.** Repeat-visit data like these could also go to the N-mixture module.
* **Spatial autocorrelation among stations is not modelled** (Known limit 7).
* **Associations, not causes.** Canopy, elevation and forest type were not manipulated.

## (g) Commands to reproduce (from the `ecotools/` folder)

```bash
# 0. (already done) how the example data and truth were made
scripts/ejob -c 1 -m 2 -t 0:10:00 -- scripts/Rscript dev/glm/make_examples.R
scripts/ejob -c 1 -m 2 -t 0:10:00 -l logs/glm_si_truth.log -- scripts/Rscript dev/glm/si_truth.R
# 1. time a few fits
scripts/ejob -c 4 -m 8 -t 0:20:00 -l logs/glm_si_selftest.log -- scripts/Rscript core/run.R selftest glm \
  --config modules/glm/example/recovery.yaml --nsim 8 --cores 4 --out dev/glm/si_selftest
# 2. planted truth at this design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/glm_si_recover.log -- scripts/Rscript core/run.R recover glm \
  --config modules/glm/example/recovery.yaml --nsim 200 --cores 16 --out results/glm/si/recovery
scripts/ejob -c 1 -m 2 -t 0:10:00 -l logs/glm_si_reccheck.log -- scripts/Rscript dev/glm/si_recovery_check.R
# 3. fit the field data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/glm_si_fit.log -- scripts/Rscript core/run.R fit glm \
  --config modules/glm/example/config.yaml --out results/glm/si/fit
# 4. restyle figures without refitting (after editing the figures: block)
scripts/ejob -c 2 -m 8 -t 0:20:00 -l logs/glm_si_figs.log -- scripts/Rscript core/run.R figures glm \
  --fit results/glm/si/fit/glm_fit.rds --config modules/glm/example/config.yaml
```

Outputs: `results/glm/si/fit/glm_report.html`, `glm_estimates.csv`, `figures/`;
`results/glm/si/recovery/glm_recover_summary.csv` and `glm_recover_report.html`;
`results/glm/si/si_summary.json`.
