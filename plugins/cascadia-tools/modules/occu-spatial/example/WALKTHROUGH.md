# Worked example: where do fishers occur, and does canopy cover matter?

Module: `occu-spatial` (spatial single-species occupancy, spOccupancy `spPGOcc`). This is
Appendix S1 and the website's "try an example" walkthrough.

**The data are simulated** (see `README_data.md`). We know the true answer, so we can check
the analysis against it. Every number below comes from the runs listed in section (g).

---

## (a) The question, in the researcher's words

> "We ran 160 baited camera stations for fishers across a 30 x 30 km block of the western
> Cascades last summer, up to six weeks each. I want to know whether fishers are more likely
> to use areas with more canopy cover, and whether elevation matters. We certainly missed
> fishers at some stations. I also suspect fishers are clumped in space for reasons our
> covariates don't capture. I'd like a map of where they probably are."

## (b) The files she has

`fisher_detections.csv`: one row per station-week that the camera was working (903 rows).

```
"station","week","start_date","days_active","detected"
"CAM-001",1,"2024-06-12",7,0
"CAM-001",2,"2024-06-19",7,0
"CAM-001",3,"2024-06-26",7,0
"CAM-001",4,"2024-07-03",7,0
```

`camera_stations.csv`: one row per station (160 rows), UTM zone 10N coordinates in metres.

```
"station","utm_e","utm_n","canopy_cover_pct","elevation_m","ownership"
"CAM-001",498484,4932345,69.1,922,"federal"
"CAM-002",495420,4941570,61,1368,"private industrial"
"CAM-003",487599,4920094,57.1,830,"federal"
"CAM-004",483436,4948255,36.1,1243,"federal"
```

`prediction_grid_1km.csv`: 900 cells of 1 km^2 covering the block, with the same two
covariates. It is used for the maps.

```
"utm_e","utm_n","canopy_cover_pct","elevation_m"
470500,4920500,77.7,615
471500,4920500,66.3,708
472500,4920500,60.7,797
473500,4920500,63.5,855
```

## (c) What she types to Claude

> I have fisher camera-trap data from 160 stations in the western Cascades, summer 2024, up
> to 6 weekly occasions per camera (some cameras died early). fisher_detections.csv has one
> row per station-week (station, week, start date, days the camera worked, detected 1/0).
> camera_stations.csv has UTM coordinates, canopy cover (%), elevation (m) and ownership.
> prediction_grid_1km.csv is a 1-km grid with the same covariates. Is fisher occupancy related
> to canopy cover and elevation? Please account for fishers we missed and for spatial
> clustering, and make a map of predicted occupancy. I need results and figures for a paper.

## (d) What the agent does, step by step (skill `eco-analyze`)

**1. Project folder.** The raw files stay untouched. For this example they live in
`modules/occu-spatial/example/` and the outputs go to `results/occu-spatial/si/`.

**2. Data intake (`eco-data-intake`).** The agent reads the three files and describes them
back:

- 160 stations, 903 station-weeks. 138 stations ran all 6 weeks. 22 cameras failed early and
  have 2-5 weeks. So 57 of 960 possible weeks (5.9%) are missing. Missing weeks have no row;
  they are not zeros.
- Fishers were photographed in 140 station-weeks at 70 stations. Naive occupancy is
  70/160 = 0.44.
- Weeks started from 3 June to 28 July 2024. 819 of 903 weeks (91%) were full 7-day weeks.
  The rest ran 2-6 days.
- Canopy cover is 24.7-88.7%. Elevation is 538-1678 m. The two are weakly correlated
  (r = -0.22). No covariate values are missing.
- Coordinates are in metres. Stations are 1.0-40.3 km apart (5% quantile 4.2 km).
- `ownership` is not part of her question, so it is left out.

**3. Choose the module.** The data are detection/non-detection records with repeat weeks in
one season, at mapped locations, and she asks about habitat and spatial clustering. That
calls for a spatial occupancy model with a nearest-neighbour Gaussian process (Doser et al.
2022): `occu-spatial`, engine `spPGOcc`.

**4. Analysis plan (`eco-analysis-plan`).** The plan is written before any fit:

- Estimands:
  - occupancy at average habitat;
  - the canopy and elevation effects on occupancy;
  - the share of the 160 stations that are occupied;
  - weekly detection, and the effect of camera days on detection;
  - an occupancy map.
- Occupancy: `~ canopy_cover_pct + elevation_m` plus a spatial random effect. Detection per
  week: `~ days_active`.
- Covariates are scaled to mean 0 and SD 1, so effects are "per 1 SD". 1 SD is 12.9% canopy,
  246 m elevation and 1.0 camera-day.
- Priors are the spOccupancy defaults: Normal(0, var 2.72) on coefficients and IG(2, 1) on the
  spatial variance. The range prior is uniform between the 5% quantile and the maximum distance.
- Checks:
  - planted-truth recovery at this design;
  - R-hat < 1.01 and bulk ESS >= 400;
  - a Freeman-Tukey posterior predictive check.
- The spatial variance and range are described, not interpreted. The module cannot recover
  them reliably (VALIDATION.md).

**5. Planted truth (`eco-planted-truth`).** Before touching her data, the agent simulates
data sets that look like hers and checks that the model gets the right answer. Settings
(`recovery.yaml`):
- the design is 160 sites x 6 weeks in a 30-km block, with 6% of weeks missing;
- the true values are those of the example.

One fit takes about 77 s, so the agent ran 100 replicates (16 cores). All 100 converged, and
none failed.

| quantity | truth | mean estimate | relative bias | coverage of 95% CrI | power | gate |
|---|---|---|---|---|---|---|
| occupancy probability (average habitat) | 0.426 | 0.448 | +5.1% | 0.93 | - | PASS |
| canopy effect (per SD) | 0.90 | 0.974 | +8.2% | 0.99 | 0.97 | PASS |
| elevation effect (per SD) | -0.40 | -0.446 | +11.5% (more negative; within 3 MC SE) | 0.95 | 0.33 | PASS |
| detection probability (weekly) | 0.285 | 0.281 | -1.3% | 0.96 | - | PASS |
| camera-days effect on detection (per SD) | 0.38 | 0.380 | +0.1% | 0.90 | 0.91 | PASS |
| mean occupancy over sites | 0.445 | 0.461 | +3.6% | 0.95 | - | PASS |
| proportion of sites occupied | 0.452 | 0.461 | +2.0% | 0.97 | - | PASS |
| spatial variance | 1.0 | 1.48 | +48% | 1.00 | - | FAIL (bias) |
| effective spatial range (km) | 8 | 9.70 | +21% | 1.00 | - | FAIL (bias) |

Power is the share of replicates whose 95% interval excluded zero. Truths for mean occupancy
and proportion occupied are averages over the replicates' own realised values. The elevation
bias (-0.046) is above 10% of the truth but below 3 Monte Carlo SE (0.086), so it passes.

What the agent tells her, in plain words: *"At your design the model recovers everything you
asked about. Bias is under about 10%, and the 95% intervals contain the truth 90-99% of the
time. The canopy effect tends to come out about 8% too strong. An elevation effect as weak as
the one we simulated is detected only one time in three. The amount and scale of spatial
clustering are not recoverable from these data (+48% and +21% bias), so we will describe them,
not interpret them."*

**6. Fit.** The agent runs `core/run.R fit` with `example/config.yaml` (section g). That is
3 chains x 50,000 iterations, burn-in 20,000, thin 20, 15 neighbours, and it took 283 s.

**7. Checks (`eco-convergence`, `eco-model-checking`).**

- Convergence: max R-hat 1.006, min bulk ESS 672 (spatial range and variance). Every
  covariate effect and probability has R-hat <= 1.002 and ESS >= 2,190. The traces of all
  3 chains overlap with no trend (`trace.png`). The spatial-variance chain has occasional
  spikes (up to about 30), which is typical of a weakly identified variance.
- Goodness of fit: Freeman-Tukey posterior predictive check on site-level detection counts,
  Bayesian p = 0.55. No sign of lack of fit (values near 0 or 1 would show it).
- Prior check: the spatial variance posterior (95% CrI 0.17-3.89) is almost the same as its
  prior (0.18-4.13). The range posterior also lies on top of its prior (`range.png`). The
  report therefore recommends also fitting the non-spatial model. The agent did that
  (`model: PGOcc`; section e).

**8. Interpretation (`eco-interpretation`).** These are excerpts from the report. Each
number was checked against the estimates table.

> "Occupancy probability at average covariate values and zero spatial effect is 0.50 (95% CrI
> 0.36 to 0.65). Averaged over the 160 surveyed sites, mean occupancy probability is 0.51 (95%
> CrI 0.43 to 0.59), and the estimated proportion of surveyed sites that are occupied is 0.51
> (95% CrI 0.46 to 0.57), compared with 0.44 of sites where the species was detected at least
> once (the difference is the correction for missed detections)."
>
> "Occupancy and Canopy cover (%): each 1-SD increase (1 SD = 12.9 in original units)
> multiplies the odds of occupancy by 2.75 (95% CrI 1.65 to 5.26); predicted occupancy goes
> from 0.28 at 1 SD below the mean to 0.73 at 1 SD above it ... Per 10-unit increase the odds
> of occupancy are multiplied by 2.19 (95% CrI 1.48 to 3.64)."
>
> "Occupancy and Elevation (m): each 1-SD increase (1 SD = 246 in original units) multiplies
> the odds of occupancy by 0.58 (95% CrI 0.32 to 0.97) ... Per 200-unit increase the odds of
> occupancy are multiplied by 0.64 (95% CrI 0.40 to 0.98)."
>
> "Detection: when the species is present, a single survey detects it with probability 0.30
> (95% CrI 0.25 to 0.35) ... Per 1-unit increase [in camera days] the odds of detection are
> multiplied by 1.48 (95% CrI 1.12 to 2.09). This is a detection effect (how easily the
> species is found), not a habitat effect."
>
> "The residual spatial variance is 0.95 (95% CrI 0.17 to 3.89) on the logit scale; its prior
> alone gives a 95% interval of 0.18 to 4.13, so the data add little information: these
> numbers mostly restate the prior and do not show whether residual clustering is strong or
> weak."
>
> "All effects are associations from observational data, not causal effects."

**Figures** (in `results/occu-spatial/si/fit/figures/`, as PNG and PDF, ESA column widths):

| id | caption |
|---|---|
| `effect_canopy_cover_pct` | Occupancy probability vs canopy cover (%) (posterior mean, 95% and 50% credible bands), other covariates at their means and spatial effect at zero; rug = surveyed sites. **Key figure.** |
| `effect_elevation_m` | Occupancy probability vs elevation (m), as above. |
| `detection_days_active` | Per-survey detection probability vs days the camera was active in the week (posterior mean, 95% and 50% bands); rug = surveys. |
| `coefficients` | Occupancy (circles) and detection (triangles) coefficients on the logit scale per 1 SD: posterior means with 95% (thin) and 50% (thick) intervals; dashed line = no association. |
| `occupancy_map` | Predicted occupancy probability (A, posterior mean) and its uncertainty (B, posterior SD) on the 1-km grid; points = stations. |
| `spatial_effect_map` | Posterior mean of the spatial random effect (logit scale) on the grid: red = more occupied than the covariates predict, blue = less. |
| `range` | Posterior of the effective spatial range (filled) with its prior (dashed). The overlap shows how little the data say about the range. |
| `gof` | Posterior predictive check: Freeman-Tukey discrepancy, replicated vs observed, one point per draw; 1:1 line; Bayesian p = 0.55. |
| `trace` | MCMC traces of the 3 chains for the main estimands. |

The agent opened every PNG. Labels, units, bands and rugs are readable at print size. The
occupancy map uses the full page width (`per_figure: occupancy_map: {width: double}` in the
config), so its two colour bars do not crowd. No other figure changes were needed.

**9. Gate (`eco-acceptance-gate`).** Verdict: **DONE**, for the estimands that answer her
question.
- Every interpreted quantity passed planted truth at her design. These are the canopy and
  elevation effects, occupancy, detection, the camera-days effect, mean occupancy and the
  proportion occupied.
- The real-data fit converged (R-hat <= 1.006, ESS >= 672) and passed the predictive check
  (p = 0.55).
- The module-level independent routes and negative control are in VALIDATION.md:
  - PGOcc agrees with `unmarked` (r >= 0.999);
  - the simulator agrees with `spOccupancy::simOcc`;
  - the non-spatial negative control fails the gate, as it should.
- Spatial variance and range failed on bias. They are reported as descriptive only, not
  interpreted.

**10. Deliver.** This is the summary she receives:

> Fishers were more likely to occupy stations with more canopy cover. Each extra 10% canopy
> multiplied the odds of occupancy by 2.2 (95% CrI 1.5 to 3.6), and predicted occupancy rose
> from 0.28 to 0.73 between 1 SD below and 1 SD above average canopy. Occupancy was lower at
> higher elevation: the odds fell by about a third per 200 m (odds ratio 0.64, 95% CrI 0.40 to
> 0.98). That interval only just excludes no effect, and our simulations show only a 1-in-3
> chance of detecting an effect this weak at your design. About 51% of your stations were
> occupied (95% CrI 46 to 57%), against 44% where a fisher was photographed. A camera finds a
> fisher at an occupied station in a given full week with probability about 0.30. Your data
> cannot say how strongly fishers cluster beyond what canopy and elevation explain. A
> non-spatial model gives the same answer with a somewhat smaller canopy effect, so the
> conclusions don't depend on that choice.

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% CrI), spatial model | inside CrI? | non-spatial PGOcc (95% CrI) |
|---|---|---|---|---|
| occupancy probability (average habitat) | 0.426 | 0.50 (0.36-0.65) | yes | 0.50 (0.40-0.61) |
| canopy effect (per SD) | 0.90 | 1.01 (0.50-1.66) | yes | 0.87 (0.43-1.37) |
| elevation effect (per SD) | -0.40 | -0.55 (-1.14 to -0.03) | yes | -0.45 (-0.92 to -0.05) |
| weekly detection probability | 0.285 | 0.300 (0.249-0.353) | yes | 0.301 (0.249-0.354) |
| camera-days effect (per SD) | 0.38 | 0.40 (0.12-0.74) | yes | 0.40 (0.12-0.76) |
| mean occupancy over stations | 0.483 | 0.51 (0.43-0.59) | yes | 0.51 (0.43-0.60) |
| proportion of stations occupied | 0.5125 (82/160) | 0.508 (0.463-0.566) | yes | 0.508 (0.463-0.569) |
| spatial variance | 1.0 | 0.95 (0.17-3.89) | yes (but prior-dominated) | - |
| effective range (km) | 8 | 8.8 (4.3-28.5) | yes (but prior-dominated) | - |

All nine 95% intervals contain the truth. The spatial model's slopes are a little stronger
than the truth (+12% canopy, +38% elevation). That matches the upward lean seen in
planted truth (+8% and +11.5% on average) and the slope inflation that VALIDATION.md
documents for spPGOcc. The non-spatial fit is closer to the truth here (canopy 0.87,
elevation -0.45), but it would have narrow intervals if real clustering were strong. So the
report gives both.

## (f) Limits from VALIDATION.md that matter at this design

1. **Spatial variance and range are not recoverable.** In planted truth they were +48% and
   +21% biased. Their posteriors here restate the priors. Read the spatial maps as smoothing
   of residuals, not as a measured clustering or home-range scale.
2. **Slopes lean away from zero at modest sample sizes.** VALIDATION.md found:
   - 10-28% bias at 100-120 sites with p 0.25-0.3;
   - about 10% inflation when the spatial model is fitted to non-spatial data.
   At 160 sites the lean was +8% (canopy) and +11.5% (elevation), within the gate. Trust the
   direction, and treat the effect size as possibly a little high. Report the PGOcc fit beside it.
3. **Low power for the elevation effect.** An effect of -0.40 per SD was detected in only 33%
   of replicates. A null elevation result in another season would not mean "no effect".
4. **Spatial confounding.** When a covariate is itself smooth in space (elevation here has a
   15-km range), the spatial effect can absorb part of its effect (VALIDATION.md: -20% for a
   90%-spatial covariate). The recovery simulator used less spatial covariates than the
   example data, so this is a reason for caution on elevation.
5. **Simulator differences.** In recovery:
   - the covariates were independent rather than correlated (r = -0.22);
   - missing weeks were scattered at random rather than as early camera failures;
   - `days_active` was a standard-normal covariate rather than 91% full weeks.
   These differences are documented in `recovery.yaml`.
6. **Assumptions not tested here.** The model assumes closure over the 8-week season, no false
   positives (misidentified photos), and that detection does not vary among stations beyond
   camera days. It cannot separate fisher home ranges that span several stations from true
   site occupancy. Stations 1 km apart may share individuals, so "occupancy" is better read as
   "use".
7. **Known core issue.** The gate's ">= 90% converged" clause cannot fail for non-convergence
   (VALIDATION.md). The agent counted by hand: 100/100 replicates converged.

## (g) Commands to reproduce

Run from `.`. Always use `scripts/ejob`, because plain `srun`
hangs on the interactive node.

```bash
# 1. Example data (already in example/; set.seed(20251))
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/occu-spatial_make_example.log -- scripts/Rscript dev/occu-spatial/make_example.R
# 2. True values on the module's scale; confirms the CSVs are reproduced byte for byte
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/occu-spatial_si_truth.log -- scripts/Rscript dev/occu-spatial/si_truth.R
# 3. Planted-truth recovery at this design (100 reps, ~77 s each, about 10 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/occu-spatial_si_recover.log -- scripts/Rscript core/run.R recover occu-spatial --config modules/occu-spatial/example/recovery.yaml --nsim 100 --cores 16 --out results/occu-spatial/si/recovery
# 4. Real-data fit, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/occu-spatial_si_fit.log -- scripts/Rscript core/run.R fit occu-spatial --config modules/occu-spatial/example/config.yaml --out results/occu-spatial/si/fit
# 5. Non-spatial comparison fit recommended by the report (same config, model: PGOcc)
scripts/ejob -c 4 -m 16 -t 1:00:00 -l logs/occu-spatial_si_fit_pgocc.log -- scripts/Rscript core/run.R fit occu-spatial --config dev/occu-spatial/si_config_pgocc.yaml --out results/occu-spatial/si/fit_pgocc
```

Outputs:
- `results/occu-spatial/si/recovery/occu-spatial_recover_summary.csv` (gate table);
- `results/occu-spatial/si/fit/occu-spatial_report.html` (report);
- `results/occu-spatial/si/fit/figures/` (figures);
- `results/occu-spatial/si/si_summary.json` (machine-readable summary of this appendix).
