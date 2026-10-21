# Worked example: do invasive Cuban treefrogs push out native green treefrogs?

Module: `occu-waddle` (two-species conditional occupancy; Richmond et al.
2010, Waddle et al. 2010). Appendix S1 and the website's "try an example".

> The data in this example are **simulated** (`README_data.md`), so we know
> the true answer and can check the analysis against it. Every number below
> comes from the runs listed in section (g).

## (a) The question, in the researcher's words

"We check PVC pipes for treefrogs at 180 sites in the Everglades. Are green
treefrogs found at fewer sites where Cuban treefrogs live? And when Cuban
treefrogs are in the pipes, are green treefrogs harder to find? Distance
to roads and canopy cover might matter too."

## (b) The files the researcher has

`pipe_surveys.csv`: one row per site and visit (900 rows). A visit that
did not happen has `pipes_checked = 0` and empty detection cells.

```
site_id,visit,date,observer,air_temp_c,pipes_checked,cuban_treefrog,green_treefrog
WCA-001,1,2023-06-05,JHW,22.9,10,1,0
WCA-001,2,2023-06-29,MS,29.6,10,1,0
WCA-001,3,2023-07-20,JHW,21.3,10,0,1
WCA-001,4,2023-08-09,MS,23.7,10,0,0
```

`sites.csv`: one row per site (180 rows).

```
site_id,dist_road_m,canopy_pct,lat,lon
WCA-001,370,60,25.72798,-80.4748
WCA-002,597,15,25.9306,-80.60891
WCA-003,1483,42,25.92173,-80.47058
WCA-004,727,41,25.97727,-80.84056
```

## (c) What the researcher types to Claude

> I have pipe-survey data for Cuban treefrogs and green treefrogs at 180
> sites, 5 visits each, in `pipe_surveys.csv`, and site information
> (distance to road, canopy cover) in `sites.csv`. I want to know whether
> green treefrogs occupy fewer sites where Cuban treefrogs are present, and
> whether green treefrogs are harder to detect on visits when Cuban
> treefrogs are found. Cuban treefrog detection probably depends on air
> temperature. Please run the analysis and give me a report with figures I
> can use in a paper.

## (d) What the agent does at each step (eco-analyze)

**1. Project folder.** The agent makes a project folder and keeps a log of
every command. (For this example the files live in
`modules/occu-waddle/example/` and `results/occu-waddle/si/`.)

**2. Data intake (eco-data-intake).** The agent reads both files and
describes them back:

- 180 sites, 5 planned visits, 5 June to 31 August 2023.
- 840 visits done, 60 missed (6.7%). 126 sites have 5 visits, 48 have 4, 6 have 3.
- Cuban treefrogs were found on 226 visits at 90 sites. Green treefrogs on
  209 visits at 81 sites. Both were found at 23 sites.
- Distance to road 29-5409 m (median 584; strongly skewed, so it enters as
  a log). Canopy cover 3-89%. Air temperature 17.8-32.9 C.
- Every site in `pipe_surveys.csv` is in `sites.csv`; no detections are
  missing on done visits; `check_data()` reports no problems.

**3. Choose the module.** Two species, repeated 0/1 surveys, and a clear
dominant species (the large invasive predator). The catalog points to
`occu-waddle`. The symmetric alternative, `occu-rota`, is for when no
species can be called dominant. The agent recommends `occu-waddle` and says
that choosing Cuban treefrog as dominant is an assumption, not a result.

**4. Plan (eco-analysis-plan).** Agreed with the researcher before fitting:

- Estimands: green treefrog occupancy where Cuban treefrogs are present
  (psiBA) vs absent (psiBa); their log-odds ratio and the species
  interaction factor (SIF); green treefrog detection on visits when Cuban
  treefrogs were found (rBA) vs missed (rBa).
- Covariates: Cuban occupancy ~ log(distance to road); green occupancy
  (both conditions) ~ canopy cover; detection of each species alone ~ air
  temperature; detection where both occur ~ 1.
- Engine: Bayesian (Nimble), 3 chains x 6000 iterations (1000 burn-in),
  weakly informative priors (Logistic(0, 1) intercepts, Normal(0, 1)
  slopes per SD).
- Checks: R-hat < 1.01, bulk ESS >= 400, two posterior predictive checks.
- Not claimed: that Cuban treefrogs *cause* the pattern.

**5. Planted truth (eco-planted-truth).** Before touching the real data,
the agent simulates 200 data sets at the same design (180 sites x 5
visits, 6.7% missed visits) with values like the ones expected, fits each
one and checks that the truth is recovered (`recovery.yaml`). A selftest
(4 fits) took about 100 s per fit including compilation, so 200 replicates
on 16 cores were affordable (median 36 s per fit).

Result: **all 18 quantities PASS** the gate. 200 of 200 fits converged.
Relative bias was -4.5% to +4.3%. 95% interval coverage was 0.94-0.97
(needed: at least 0.904).

| quantity | truth | mean estimate | rel. bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| Cuban occupancy (psiA) | 0.525 | 0.527 | +0.4% | 0.94 | | PASS |
| ... slope (road distance) | -0.90 | -0.94 | -4.5% | 0.95 | 1.00 | PASS |
| Green occupancy, Cuban present (psiBA) | 0.250 | 0.256 | +2.3% | 0.96 | | PASS |
| ... slope (canopy) | 0.60 | 0.61 | +1.3% | 0.94 | 0.53 | PASS |
| Green occupancy, Cuban absent (psiBa) | 0.690 | 0.698 | +1.2% | 0.96 | | PASS |
| ... slope (canopy) | 0.50 | 0.51 | +1.4% | 0.94 | 0.39 | PASS |
| Cuban detection alone (pA) | 0.550 | 0.549 | -0.2% | 0.96 | | PASS |
| ... slope (temperature) | 0.378 | 0.394 | +4.3% | 0.95 | 0.93 | PASS |
| Green detection alone (pB) | 0.600 | 0.606 | +1.0% | 0.94 | | PASS |
| ... slope (temperature) | 0.270 | 0.278 | +2.9% | 0.94 | 0.57 | PASS |
| Cuban detection, both present (rA) | 0.550 | 0.544 | -1.1% | 0.96 | | PASS |
| Green detection, Cuban found (rBA) | 0.300 | 0.305 | +1.7% | 0.94 | | PASS |
| Green detection, Cuban missed (rBa) | 0.550 | 0.540 | -1.8% | 0.94 | | PASS |
| Green occupancy overall (psiB) | 0.459 | 0.465 | +1.3% | 0.95 | | PASS |
| P(both) | 0.131 | 0.134 | +2.2% | 0.94 | | PASS |
| SIF | 0.544 | 0.550 | +1.0% | 0.97 | 1.00 | PASS |
| Occupancy log-odds ratio | -1.90 | -1.95 | -2.5% | 0.96 | 1.00 | PASS |
| Detection log-odds ratio | -1.05 | -1.01 | +3.2% | 0.94 | 0.61 | PASS |

In one sentence for the researcher: "At your design the method gets the
right answer on average for every quantity, and it finds an occupancy
effect of this size almost every time, but it finds a detection effect of
this size only 6 times in 10, and canopy effects of this size only 4-5
times in 10."

The simulator has one shared habitat covariate for all three occupancy
probabilities, whereas the real design has road distance for Cuban and
canopy for green treefrogs. Intercepts and slopes match the example
exactly; only that sharing differs.

**6. Fit.** `core/run.R fit` with `config.yaml` (section g). It takes a few
minutes on 4 cores (3 chains in parallel, plus the report and figures).

**7. Checks (eco-convergence, eco-model-checking).**

- Convergence: max R-hat 1.000, min bulk ESS 11,939 (need < 1.01 and >= 400). Trace plots overlap with no trend.
- Goodness of fit: co-detection pattern p = 0.48; per-site detection
  counts p = 0.65. Both are far from 0 and 1, so the model reproduces
  these features of the data.

**8. Interpretation (eco-interpretation).** Excerpts from the report
(`results/occu-waddle/si/fit/occu-waddle_report.html`):

> Where Cuban treefrog is present, green treefrog occupies 0.31 (95% CrI
> 0.21 to 0.43) of sites; where Cuban treefrog is absent, 0.67 (95% CrI
> 0.56 to 0.78). The odds of green treefrog occupying a site are 0.21 (95%
> CrI 0.10 to 0.44) times as high where Cuban treefrog is present.

> The species interaction factor ... is 0.64 (95% CrI 0.46 to 0.81): at
> average conditions the two species co-occur less often than expected if
> they were distributed independently (a segregation pattern).

> Occupancy of green treefrog where Cuban treefrog is absent (psiBa): the
> odds are multiplied by 2.13 (95% CrI 1.35 to 3.63) per 1 SD of canopy_pct
> (1 SD = 23.1); i.e. it increases with canopy_pct (from 0.33 at canopy_pct
> = 3 to 0.89 at canopy_pct = 89).

> Where both occur, green treefrog is detected on 0.23 (95% CrI 0.13 to
> 0.36) of surveys that detected Cuban treefrog (rBA) and on 0.50 (95% CrI
> 0.36 to 0.64) of surveys that missed Cuban treefrog (rBa) ... it is a
> detection effect, not an occupancy effect.

> A difference between psiBA and psiBa, or an SIF away from 1, is a
> co-occurrence pattern. It is consistent with competition or predation
> but equally with both species responding to habitat features not in the
> model (Blanchet et al. 2020).

The agent checked every number in the report against the estimates table
(odds ratios = exp of the log-scale values).

Figures (`results/occu-waddle/si/fit/figures/`, PDF + PNG):

| figure | caption (short) |
|---|---|
| `conditional_occupancy` | **Key figure.** (A) Occupancy at average conditions: Cuban, green overall, green with Cuban absent / present. (B) Green treefrog occupancy along canopy cover, Cuban present vs absent. |
| `sif` | Species interaction factor (1 = independent): (A) at average conditions, (B) along distance to road. |
| `detection_interaction` | Per-visit detection: pA, rA, pB, rBa, rBA. rBA vs rBa is the detection interaction. |
| `coefficients` | Every covariate effect as an odds ratio per 1 SD. |
| `occupancy_dist_road_m` | Occupancy of both species along distance to road (log axis). |
| `occupancy_canopy_pct` | Green treefrog occupancy along canopy cover. |
| `detection_air_temp_c` | Detection of each species along air temperature. |
| `gof` | Posterior predictive check of the co-detection classes (p = 0.48). |
| `trace` | MCMC traces of all 13 coefficients, 3 chains. |
| `posterior` | Posterior densities of psiBA, psiBa, SIF, the two log-odds ratios, rBA and rBa, with null lines. |

All figure PNGs were opened and checked. No presentation problems needed
fixing; the `figures:` block of `config.yaml` was left as is.

**9. Gate (eco-acceptance-gate).** Verdict **DONE**: the planted-truth
gate passed at this design (18/18), convergence and fit checks passed, the
report's numbers were verified, and the limits below are stated.

**10. Deliver.** Summary for the researcher:

> Green treefrogs occupied fewer sites where Cuban treefrogs were present
> (0.31, 95% CrI 0.21-0.43) than where they were absent (0.67, 0.56-0.78);
> the odds were 0.21 (0.10-0.44) times as high. The two species occur
> together less often than chance (SIF 0.64, 0.46-0.81). Where both live,
> green treefrogs were found on fewer visits when Cuban treefrogs were also
> in the pipes (0.23 vs 0.50; log-odds ratio -1.18, -2.01 to -0.40). Cuban
> treefrogs were more common near roads, and green treefrogs (where Cuban
> treefrogs were absent) more common under closed canopy. This is a
> co-occurrence pattern: it fits predation or competition, but shared
> habitat preferences not in the model could also produce it.

## (e) Estimates versus the known truth

Every one of the 18 95% credible intervals contains the true value.

| quantity | truth | estimate (95% CrI) |
|---|---|---|
| Cuban occupancy (psiA) | 0.525 | 0.52 (0.44-0.61) |
| ... slope, log road distance | -0.90 | -0.92 (-1.34 to -0.56) |
| Green occupancy, Cuban present (psiBA) | 0.250 | 0.31 (0.21-0.43) |
| ... slope, canopy | 0.60 | 0.13 (-0.40 to 0.68) |
| Green occupancy, Cuban absent (psiBa) | 0.690 | 0.67 (0.56-0.78) |
| ... slope, canopy | 0.50 | 0.76 (0.30-1.29) |
| Cuban detection alone (pA) | 0.550 | 0.56 (0.49-0.63) |
| ... slope, temperature | 0.378 | 0.59 (0.34-0.87) |
| Green detection alone (pB) | 0.600 | 0.60 (0.54-0.66) |
| ... slope, temperature | 0.270 | 0.19 (-0.06 to 0.45) |
| Cuban detection, both present (rA) | 0.55 | 0.46 (0.35-0.57) |
| Green detection, Cuban found (rBA) | 0.30 | 0.23 (0.13-0.36) |
| Green detection, Cuban missed (rBa) | 0.55 | 0.50 (0.36-0.64) |
| Green occupancy overall (psiB) | 0.459 | 0.48 (0.40-0.57) |
| P(both) | 0.131 | 0.16 (0.10-0.23) |
| SIF | 0.544 | 0.64 (0.46-0.81) |
| Occupancy log-odds ratio | -1.90 | -1.55 (-2.31 to -0.82) |
| Detection log-odds ratio | -1.05 | -1.18 (-2.01 to -0.40) |

Two lessons a reader can see here:

- **One data set is not the average.** The canopy slope where Cuban
  treefrogs are present was truly +0.60 but estimated at +0.13; the
  interval (-0.40 to 0.68) still contains the truth. The recovery run
  predicted this: power for that slope is only 0.53. "No clear effect" is
  not "no effect".
- The occupancy interaction was estimated slightly weaker than the truth
  (SIF 0.64 vs 0.544), in the direction VALIDATION.md warns about, but
  well within the interval. The direction of every interaction is right.

## (f) Limits from VALIDATION.md that matter at this design

- **Power for the detection interaction is modest.** Here 0.61 for a
  log-odds ratio of -1.05 at 180 x 5. A null result for rBA vs rBa would
  not mean "no effect".
- **Canopy effects on green treefrog occupancy have low power** (0.39 and
  0.53 here). Their sizes are uncertain by about +-15% even on average
  (psiBa slope; VALIDATION.md sections 1 and 8).
- **Dominance is assumed.** Results are conditional on calling the Cuban
  treefrog dominant; the model does not test the direction.
- **Co-occurrence is not interaction.** Shared, unmeasured habitat drivers
  give the same pattern (Blanchet et al. 2020). The gate checks
  estimation, not causation.
- **A detection interaction can come from shared survey conditions.** Air
  temperature is in the alone-detection formulas, but other day-level
  conditions that affect both species could mimic rBA < rBa.
- With 180 sites the design is close to the validated 200 x 4; the
  small-sample biases seen at 60 sites (psiBA, P(both), SIF pulled toward
  "no interaction") were not seen here (all within +-2.5%).
- Not validated: low detection (< 0.3), factor covariates, covariates on
  rBA / rBa. None is used here.

## (g) Commands to reproduce

Run from `.`. Always use `scripts/ejob`
(plain `srun` hangs on the interactive node).

```bash
# 0. (optional) regenerate the example data and the truth table
scripts/ejob -c 1 -m 4 -t 0:15:00 -l logs/occu-waddle_si_truth.log -- \
  scripts/Rscript dev/occu-waddle/si_truth.R
# 1. time a short selftest at this design
scripts/ejob -c 4 -m 8 -t 0:30:00 -l logs/occu-waddle_si_selftest.log -- \
  scripts/Rscript core/run.R selftest occu-waddle --config modules/occu-waddle/example/recovery.yaml \
  --nsim 4 --cores 4 --out dev/occu-waddle/si_selftest
# 2. planted-truth recovery, 200 replicates
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/occu-waddle_si_recover.log -- \
  scripts/Rscript core/run.R recover occu-waddle --config modules/occu-waddle/example/recovery.yaml \
  --nsim 200 --cores 16 --out results/occu-waddle/si/recovery
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/occu-waddle_si_mdtable.log -- \
  scripts/Rscript dev/occu-waddle/mdtable.R results/occu-waddle/si/recovery
# 3. fit the example data (report + figures)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/occu-waddle_si_fit.log -- \
  scripts/Rscript core/run.R fit occu-waddle --config modules/occu-waddle/example/config.yaml \
  --out results/occu-waddle/si/fit
# 4. redraw figures after editing the figures: block (no refit)
scripts/Rscript core/run.R figures occu-waddle --fit results/occu-waddle/si/fit/occu-waddle_fit.rds \
  --config modules/occu-waddle/example/config.yaml
```

Outputs: `results/occu-waddle/si/` (`truth.csv`, `recovery/`, `fit/`,
`si_summary.json`).
