# Worked example: bobcat density from camera traps (scr module)

This walkthrough follows one researcher from question to report. The data
are **simulated**, so we know the true answer and can check the analysis
against it (see `README_data.md`). Every number below comes from the runs
listed in section (g).

## (a) The question

> "I ran 48 trail cameras for 8 weeks last summer and identified individual
> bobcats from their spots. How many bobcats are there per square
> kilometre, and are there more bobcats where the forest canopy is denser?"

## (b) The files the researcher has

`cameras.csv` (48 rows, one per camera):

```
"camera_id","utm_e","utm_n","on_trail","weeks_active","habitat_note"
"BC01",478360,4921185,0,"11111111","off-trail game path"
"BC02",479899,4921576,0,"11111111","off-trail game path"
"BC03",481487,4921585,0,"11111111","off-trail game path"
"BC04",483093,4921284,0,"11100000","off-trail game path"
```

`weeks_active` has one digit per week: camera BC04 stopped after week 3.

`bobcat_detections.csv` (114 rows, one per photo-detection of an identified cat):

```
"cat_id","week","camera_id","date","photo"
"LYRU-001",5,"BC47","2025-07-06","BC47_20250706_5920.JPG"
"LYRU-002",1,"BC16","2025-06-05","BC16_20250605_3332.JPG"
"LYRU-002",2,"BC08","2025-06-11","BC08_20250611_6561.JPG"
"LYRU-002",2,"BC16","2025-06-09","BC16_20250609_0543.JPG"
```

`habitat_grid.csv` (2,255 rows, one per 400 m pixel within 5 km of a camera):

```
"x","y","canopy_cover_pct","elevation_m","habitat"
477354,4916385,87.7,596,1
477754,4916385,87.4,600,1
478154,4916385,86.9,600,1
478554,4916385,86.1,601,1
```

`habitat` = 0 marks a reservoir, which is not bobcat habitat.

## (c) What the researcher types to Claude

> I have bobcat camera-trap data in modules/scr/example: a camera file with
> UTM coordinates, whether each camera was on a trail and which weeks it
> worked; a file of photo-detections of individually identified bobcats by
> week and camera; and a 400 m habitat grid with canopy cover and a flag for
> the reservoir. Please estimate bobcat density and how many bobcats live in
> the area, and tell me whether density is higher where canopy cover is
> higher. Check that the method works for a design like mine before you
> trust it, and give me a report and figures I can put in a paper.

## (d) What the agent does, step by step (skills/eco-analyze)

**Step 1-2. Project folder and data intake.** The agent reads the three
files and describes them back:

- 48 cameras, about 1.6 km apart, 18 on trails. 8 weekly occasions.
  354 of 384 camera-weeks were working (8 cameras lost 1-6 weeks; their
  missing weeks are set to "not working", not "no bobcat").
- 114 detections of 33 individual bobcats. 7 cats were seen once; one was
  seen 12 times. 25 of the 33 cats were photographed at more than one
  camera (53 "spatial recaptures"). These are what tell the model how far
  bobcats move.
- Habitat grid: 2,189 habitat pixels = 350 km2, reaching 5 km beyond the
  cameras. Canopy cover 19-91% (mean 54%, SD 14.7).
- Coordinates are projected (UTM, metres). No missing values. Every camera
  in the detections file is in the camera file.

**Step 3. Choose the module.** Individually identified animals at fixed
detectors in one short season, question about density and a habitat
covariate: `scr` (closed spatial capture-recapture). Not `scr-smr` (all
cats are identifiable) and not `scr-open` (one season only).

**Step 4. Plan** (written before fitting):

- Estimands: density (bobcats per km2), abundance in the 350 km2 of
  habitat, and the canopy-cover effect on density (per SD of canopy).
- Model: SCR with log-linear density in canopy cover over the habitat grid;
  half-normal detection; baseline detection differs for trail cameras;
  binomial detections summed over the weeks each camera worked.
- Engine: Bayesian (NIMBLE, marginal likelihood), 3 chains x 20,000
  iterations, 5,000 burn-in. Cross-check with maximum likelihood (secr).
- Priors: flat on psi, p0 and sigma; Normal(0, SD 2) on coefficients.
- Adequate fit: R-hat < 1.01, bulk ESS >= 400, posterior predictive
  p-values between 0.025 and 0.975, Nimble and secr agree.
- Elevation is in the grid but not in the question, so it is left out.

**Step 5. Planted-truth check at this design.** Before touching the real
data, the agent simulates 100 data sets that look like this study (49
cameras 1,600 m apart, 8 weeks, 15% camera failures, 5 km buffer, a lake,
canopy effect 0.45 per SD, the true values of `README_data.md`) and fits
each one (`recovery.yaml`; about 3-4 minutes per fit).

Result: **PASS for all 7 quantities**; 96 of 100 fits converged (R-hat <
1.05 with the shorter recovery chains). Gate: |bias| <= 10% of the truth
(or 3 Monte Carlo SE), 95% interval coverage >= 0.883, >= 90% converged.

| quantity | truth | mean estimate | relative bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| density (per km2) | 0.25 | 0.261 | +4.2% | 0.95 | | PASS |
| abundance | 96.7 | 98.4 | +1.8% | 0.95 | | PASS |
| p0 (off-trail) | 0.12 | 0.127 | +5.7% | 0.94 | | PASS |
| sigma (m) | 1,300 | 1,291 | -0.7% | 0.95 | | PASS |
| 95% home range (km2) | 31.8 | 31.6 | -0.8% | 0.95 | | PASS |
| canopy effect (per SD) | 0.45 | 0.465 | +3.4% | 0.95 | 0.67 | PASS |
| on-trail effect | 1.10 | 1.06 | -3.7% | 0.97 | 1.00 | PASS |

(Abundance truth is larger than 82 because the simulator's state space is
a full rectangle around a square array.)

The agent notices one difference. These simulated arrays sat in average
habitat and detected about 46 bobcats each. The real array sits in thinner
canopy (44.5% near the cameras vs 53.8% over the whole area) and detected
only 33. Fewer animals means more small-sample bias, so the agent reruns
the check with density lowered to 0.17, which detects about 32 bobcats per
data set, like the real study (`recovery_matched_n.yaml`). Result: **PASS
for all 7 quantities** again, 96 of 100 converged. Density +5.2% (coverage
0.95), abundance +2.1% (0.97), p0 +5.8% (0.95), sigma -0.2% (0.95), canopy
effect -1.3% (coverage 0.885, just above the 0.883 bound), on-trail effect
-2.2% (0.97). Power to detect the canopy effect drops to 0.51: at this
sample size a real canopy effect of this size is found only about half the
time.

In one sentence to the researcher: "At a design like yours, the method
recovered density, home range and the habitat effect within 6% on average,
and its 95% intervals contained the truth about 95% of the time."

**Step 6. Fit.** `core/run.R fit scr --config modules/scr/example/config.yaml`
(about 5 minutes on 4 cores, including two compilations).

**Step 7. Check.**

- Convergence: max R-hat 1.004, min bulk ESS 992 (needs < 1.01 and >= 400).
  Traces of the three chains overlap with no trends.
- Goodness of fit, posterior predictive p-values (extreme if < 0.025 or
  > 0.975): individuals detected 0.52, total detections 0.50, cameras per
  individual 0.41, spread of detections per individual 0.47, Freeman-Tukey
  of camera totals 0.92. No sign of misfit.
- Design: the state space reaches 3.8 sigma beyond the cameras (need
  >= 3-4). Mean nearest-camera distance 1,350 m = 1.0 sigma (need <= 2).
- Independent route: maximum likelihood (secr) on the same data gives
  density 0.230 (0.151-0.350) vs Nimble 0.237 (0.157-0.341): a 2.8%
  difference, 0.14 posterior SD. sigma 1,324 vs 1,326 m; canopy effect
  0.41 vs 0.41; on-trail effect 1.19 vs 1.16.

**Step 8. Interpret.** Quoted from the report (`scr_report.html`):

> "Density was estimated at 0.24 (95% CrI 0.16 to 0.34) animals per km²,
> i.e. 23.7 per 100 km2. Over the 350 km² of habitat in the state space
> (the detector array plus a buffer), that is 82.7 (95% CrI 58.0 to 118)
> animals, of whom 33 (40%) were detected."

> "The spatial scale of detection, sigma, was 1326 (95% CrI 1165 to 1510)
> m. ... a home-range area of about 33 (95% CrI 26 to 43) km². Real home
> ranges are rarely circular, so treat this as an index of ranging scale."

> "Detection (detector-level): 'Camera on trail' (1 vs 0) multiplied the
> odds of detection by 3.2 (95% interval 1.9 to 5.4). Baseline detection
> was 0.30 (0.21 to 0.41) where Camera on trail = 1 versus 0.12 where it =
> 0. This is a detection effect (how easily animals are recorded), not a
> difference in where animals live."

> "Density and Canopy cover (%): each 1 SD (15 units) higher Canopy cover
> (%) is associated with 51% higher density (rate ratio 1.5, 95% interval
> 1.1 to 2.1). Expected density is 0.13 per km² (0.070 to 0.20) where
> Canopy cover (%) = 34.9 (10th percentile of the habitat) and 0.37 (0.19
> to 0.62) where it = 72.6 (90th percentile) ... Posterior probability that
> the association is positive: 0.99. This is a correlation between density
> and habitat across the state space, not proof that the habitat causes
> higher density."

The agent adds one caution: the lower end of the canopy interval (0.06 on
the log scale) is close to zero, and VALIDATION.md found that ~10% of
intervals exclude zero when there is no true effect at a similar number of
animals. The canopy association is supported, but its size is uncertain.

**Figures** (`results/scr/si/fit/figures/`, PNG and PDF):

| figure | caption |
|---|---|
| `density_surface` | Posterior mean realised bobcat density (per km2) over the habitat; white dots = cameras; orange crosses = most likely activity centre of each of the 33 detected cats; blank = reservoir. |
| `density_expected_surface` | Expected density from the canopy model alone; denser canopy in the south-west predicts more bobcats. |
| `density_canopy_cover_pct` | Expected density vs canopy cover with 50% and 95% credible bands; rug = canopy values in the habitat. |
| `detection_function` | Weekly detection probability vs distance from a cat's activity centre, for off-trail and on-trail cameras; dashed line = camera spacing. |
| `detection_on_trail` | Baseline detection p0 for off-trail (0.12) and on-trail (0.30) cameras with 95% intervals. |
| `abundance_posterior` | Posterior of the number of bobcats in the 350 km2 of habitat; dashed = 95% interval, dotted = 33 cats detected. |
| `captures_per_individual` | Observed numbers of cats detected 1, 2, ... times (bars) vs model predictions (points, 95% intervals). |
| `gof` | Posterior predictive checks (five statistics). |
| `coefficients` | Canopy effect on density (log scale, per SD) and trail effect on detection (logit scale), 95% intervals. |
| `engine_comparison` | Nimble (Bayesian) vs secr (maximum likelihood) estimates on the same data. |
| `trace` | MCMC traces, three chains. |
| `posterior` | Posterior densities of all reported quantities. |

**Step 9. Gate.** Verdict **DONE**. Planted truth at this design passed
(7/7, twice: at the true density and at a matched number of detected
animals). The real fit converged (R-hat 1.004, ESS >= 992) and passed all
five predictive checks. An independent engine (secr) agreed within 0.14
posterior SD. Negative controls for this module (too small a buffer,
unmodelled heterogeneity, ignored effort) fail the gate as they should
(VALIDATION.md sections 3 and 12). Every number comes from the commands in
(g). Caveat carried into the report: the canopy association is supported
(posterior probability 0.99) but its size is uncertain.

Note on figures: the `figures:` block of `config.yaml` moves the
detection-function legend below the plot, shortens two axis titles that
were clipped, and gives readable coefficient names. The trace and
posterior panels still use internal parameter names (a module default).

**Step 10. Deliver.** Summary for the researcher:

> About 24 bobcats per 100 km2 (95% interval 16-34) used the study area in
> summer 2025, or about 83 cats (58-118) in the 350 km2 of habitat around
> your cameras. You photographed 33 of them. Density was higher where the
> canopy was denser: each 15-point increase in canopy cover went with about
> 50% more bobcats (rate ratio 1.5, 1.1-2.1). This is an association, not
> proof of cause. Trail cameras were about three times as likely to record
> a bobcat (odds ratio 3.2, 1.9-5.4), which the model accounts for.

## (e) Estimates versus the known truth

| quantity | estimate (95% interval) | truth | inside interval? |
|---|---|---|---|
| density (per km2) | 0.237 (0.157-0.341) | 0.234 realised (0.25 expected) | yes |
| abundance | 82.7 (58-118) | 82 | yes |
| sigma (m) | 1,326 (1,165-1,510) | 1,300 | yes |
| 95% home range (km2) | 33.3 (25.6-42.9) | 31.8 | yes |
| p0, off-trail | 0.119 (0.082-0.164) | 0.12 | yes |
| on-trail effect (logit) | 1.16 (0.65-1.68) | 1.10 | yes |
| canopy effect (log, per SD) | 0.41 (0.06-0.76) | 0.45 | yes |

## (f) Limits that matter at this design (from VALIDATION.md)

- **Small-sample upward bias.** Posterior means of density and abundance
  lean high when few animals are detected (+12% with ~18 detected; +4-8%
  with ~42). With 33 detected cats, report the interval, not just the mean.
- **Density-covariate false positives.** With ~33-42 detected animals and no
  true effect, ~10% of 95% intervals for a habitat effect excluded zero
  (5% expected). An interval that only just clears zero is weak evidence.
- **Buffer.** Density is inflated if the state space is too small (+47% at
  1 sigma). Here it is 3.8 sigma, inside the safe range.
- **Model scope.** Half-normal detection, Euclidean distance, a closed
  population, correct photo IDs and the same ranging for males and females
  are assumed. Sex was not recorded here; if males range much farther,
  fit sex-specific detection.
- **Density covariates are associations,** not causal effects.

## (g) Commands to reproduce

```bash
cd ecotools
# 0. (already done) generate the simulated example data, seed 20250601
scripts/ejob -c 1 -m 8 -t 0:30:00 -- scripts/Rscript dev/scr/make_example.R
# 1. planted-truth recovery at this design (100 simulated data sets)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/scr_si_recover.log -- \
  scripts/Rscript core/run.R recover scr --config modules/scr/example/recovery.yaml \
  --nsim 100 --cores 16 --out results/scr/si/recovery
# 1b. same check with the number of detected animals matched to the example (~32)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/scr_si_recover_matched.log -- \
  scripts/Rscript core/run.R recover scr --config modules/scr/example/recovery_matched_n.yaml \
  --nsim 100 --cores 16 --out results/scr/si/recovery_matched_n
# 2. fit the example data: report, estimates, figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/scr_si_fit.log -- \
  scripts/Rscript core/run.R fit scr --config modules/scr/example/config.yaml \
  --out results/scr/si/fit
# 3. (optional) redraw figures after editing the figures: block
scripts/ejob -c 1 -m 8 -t 0:30:00 -- scripts/Rscript core/run.R figures scr \
  --fit results/scr/si/fit/scr_fit.rds --config modules/scr/example/config.yaml
```

Outputs: `results/scr/si/recovery/` (gate table, per-replicate results,
HTML report), `results/scr/si/fit/scr_report.html`,
`results/scr/si/fit/scr_estimates.csv`, `results/scr/si/fit/figures/`,
`results/scr/si/recovery_matched_n/`, and `results/scr/si/si_summary.json`
(built by `../.venv/bin/python dev/scr/make_si_summary.py`).
