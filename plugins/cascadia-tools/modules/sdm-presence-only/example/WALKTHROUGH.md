# Worked example (Appendix S1): Pacific marten presence-only records with the `sdm-presence-only` module

This walkthrough shows what happens when a researcher who does not code brings
presence-only records and a question to the agent. Every number below comes
from the runs listed in section (g). The data are **simulated** (see
`README_data.md`), so we can compare the estimates with the true values.

## (a) The question, in the researcher's words

> "I have about 180 marten sightings from the last 30 years in my part of the
> Cascades: iNaturalist, ODFW reports, trappers, museum specimens. Which
> habitat are martens associated with: more forest, higher elevation? Where
> are they most likely to be? Most sightings are near roads, so I don't want
> that to fool the map."

## (b) The files they have

`presence_records.csv` (192 rows: 180 marten, 12 fisher):

```
"record_id","species","date","source","easting","northing","coord_uncertainty_m"
"MACA-0104","Martes caurina","1995-02-21","iNaturalist",590978,4886164,100
"MACA-0082","Martes caurina","1995-04-01","trapper report",574705,4882671,30
"MACA-0122","Martes caurina","1995-06-09","iNaturalist",590115,4889744,30
"MACA-0093","Martes caurina","1995-06-27","museum specimen",590498,4884870,10
```

`covariate_grid.csv` (3,600 rows, one per 1 km² cell):

```
"easting","northing","forest_cover_pct","elevation_m","dist_road_km"
560500,4850500,58.3,1654,29.983
561500,4850500,53.6,1675,28.983
562500,4850500,49.1,1696,27.983
563500,4850500,44.9,1714,26.983
```

(`camera_survey.csv`, 45 camera stations, is used only in the optional
integrated analysis; see the end of section (f).)

## (c) What they type to Claude

> I have marten sightings in presence_records.csv (UTM 10N, metres; the file
> also has a few fisher records) and a 1 km habitat grid in
> covariate_grid.csv with forest cover, elevation and distance to road.
> Which habitat are martens associated with, and where are they most likely
> to occur? Most records come from near roads, so please correct for that.
> I'd like a report and figures for a paper.

## (d) What the agent does, step by step (skill `eco-analyze`)

**1. Project folder.** Raw files are copied unchanged into `data_raw/`; every
command goes into `log.md`.

**2. Data intake** (`eco-data-intake`). The agent describes the data back:

* 192 records: 180 *Martes caurina* and 12 *Pekania pennanti*. The fisher
  records are filtered out with `species: "Martes caurina"`.
* Marten records span 1995-02-21 to 2024-12-22. Sources: iNaturalist 85,
  ODFW incidental 44, trapper 34, museum 17.
* Location error is 10-500 m, smaller than the 1 km cells. No missing values.
* The grid is regular, 60 x 60 cells of 1 km² (3,600 km²). Forest cover
  7.6-98.3 %, elevation 718-2135 m, distance to road 0-30 km.
* Warning raised to the researcher: elevation and distance to road are
  correlated (r = 0.67), because roads follow valleys. So road-side recording
  can masquerade as a preference for low ground. Distance to road must be in
  the model as a sampling-bias covariate, not dropped.

**3. Choose the module.** Presence-only coordinates plus a habitat grid point
to `sdm-presence-only`. With no survey data in this question, the engine is
`ppm` (inhomogeneous Poisson point process). The agent tells the researcher
that this gives a *relative* map and habitat effects, not how many martens
there are. Absolute density needs surveys with non-detections (the
integrated engine).

**4. Analysis plan** (`eco-analysis-plan`), written before fitting:

* Estimands: habitat effects on the relative intensity of marten occurrence,
  per SD and per natural unit (per 10 percentage points of forest, per 100 m
  of elevation), with 95% CIs.
* Model: log λ(s) = α + β_forest·forest + β_elev·elevation (standardised);
  recording rate log b(s) = γ + δ·distance-to-road. The bias term is held
  constant for maps and effects.
* Not estimable and not reported: prevalence, occupancy, abundance.
* Checks: quadrature convergence, block-count goodness of fit, spatial
  block cross-validation (5 folds, 12 km blocks).

**5. Planted truth** (`eco-planted-truth`). Before fitting the real records,
the agent simulates 200 datasets at this design and checks that the model
recovers the truth. Design: 60 x 60 km of 1 km² cells, 4 roads, elevation-road
correlation 0.6, about 180 records per dataset. Truth: forest +0.9, elevation
+0.4, distance to road -1.1 per SD (`recovery.yaml`). One fit takes well under
a second, so all 200 ran in under a minute on 16 cores. 200 of 200 converged.

| estimand (per SD) | truth | mean estimate | rel. bias | RMSE | 95% coverage | CI width | power | gate |
|---|---|---|---|---|---|---|---|---|
| forest cover | 0.9 | 0.908 | +0.9% | 0.102 | 0.955 | 0.40 | 1.00 | PASS |
| elevation | 0.4 | 0.396 | -1.0% | 0.099 | 0.965 | 0.39 | 0.97 | PASS |
| distance to road (bias) | -1.1 | -1.095 | +0.5% (toward 0) | 0.130 | 0.970 | 0.54 | 1.00 | PASS |

(Coverage bound for 200 replicates: 0.904. Bias bound: max(10% of truth, 3 Monte Carlo SE).)

What the agent tells the researcher: *"At your design (about 180 records, roads
along valleys), the model recovers both habitat effects and the road effect
with bias of 1% or less and 95% intervals that contain the truth 96-97% of
the time. You have 97% power to detect an elevation effect of this size."*

**6. Fit.** `core/run.R fit sdm-presence-only --config config.yaml` (via
`scripts/ejob`; the fit takes seconds).

**7. Checks** (`eco-convergence`, `eco-model-checking`):

* Model converged (Poisson GLM by down-weighted Poisson regression; 180
  records, 3,600 quadrature points).
* Goodness of fit: record counts in 25 spatial blocks vs 1,000 simulated
  datasets. Pearson dispersion 0.94 (1 = Poisson), P = 0.69. No sign that
  records are more clustered than the model allows.
* Quadrature: using every 2nd cell changes coefficients by at most 0.17 SE
  (every 4th 0.41 SE, every 8th 1.00 SE). The full grid is adequate.
* Spatial block cross-validation: habitat map AUC 0.67, Boyce 0.68; records
  map (habitat x effort) AUC 0.76, Boyce 0.81.

**8. Interpretation.** Excerpts quoted from `sdm-presence-only_report.html`:

> "Forest cover (%): holding the other covariates at their landscape means,
> each 1 SD increase (24.6 percentage points) is associated with a change in
> relative intensity of occurrence ... of x2.99 (95% CI 2.42 to 3.70; +199%);
> per 10 percentage points, x1.56 (95% CI 1.43 to 1.70). Across the observed
> range (7.6% to 98.3%) this is a 56.4-fold difference."

> "Elevation (m): ... each 1 SD increase (280 m) is associated with a change
> in relative intensity of occurrence ... of x1.48 (95% CI 1.22 to 1.80; +48%);
> per 100 m, x1.15 (95% CI 1.07 to 1.23)."

> "Sampling bias, Distance to road (km): the rate at which records were made
> changed by x0.36 (95% CI 0.28 to 0.47; -64%) per 1 SD (6.34 km) and x0.45
> (95% CI 0.36 to 0.55) per 5 km. This describes where people looked, not
> where the species lives. ... For comparison, ignoring the bias covariate
> would give log-intensity slopes per SD of Forest cover (%) 0.98 instead of
> 1.10; Elevation (m) -0.23 instead of 0.39."

> "The map is RELATIVE: a value of 2 means twice the landscape-average
> intensity of occurrence, not a 2% or 20% chance of presence."

The agent checked each number against `sdm-presence-only_estimates.csv`
(for example exp(1.0955) = 2.99; exp(1.0955 x 10/24.65) = 1.56;
exp(0.3923 x 100/280) = 1.15).

The key biological message: once road-side recording is accounted for,
martens are associated with *higher* elevation. A naive analysis would have
said the opposite (-0.23 per SD), because records cluster along valley roads.

**Figures** (`results/sdm-presence-only/si/fit/figures/`, PDF + PNG):

| id | caption (short) |
|---|---|
| `intensity_map` | Relative intensity of marten occurrence (1 = landscape average, log colour scale), road bias removed; white points = records. |
| `bias_map` | Estimated relative recording effort from distance to road, with the records. |
| `effect_forest_cover_pct` | Relative intensity vs forest cover, elevation at its mean, bias removed; 95% and 50% bands; rug = forest cover at records. |
| `effect_elevation_m` | Same, vs elevation. |
| `bias_effect` | Relative recording rate vs distance to road (a property of the observers). |
| `response_curves_panel` | Panels A-B: forest and elevation response curves side by side. |
| `coefficients` | Multiplicative change per 1 SD (log axis) for habitat (orange) and sampling bias (blue), 95% CIs. |
| `evaluation` | Spatial block CV: (A) AUC and Boyce index per fold for the habitat map and records map; (B) predicted/expected curves. |
| `quadrature_convergence` | Coefficients as the quadrature grid is refined; band = full-grid estimate ± 0.5 SE. |
| `gof` | Observed block-count Pearson statistic vs 1,000 simulated datasets (P = 0.69). |

One presentation fix was made in the `figures:` block of `config.yaml`: the
"2200" tick label on the elevation axis of `effect_elevation_m` was clipped,
so the x-range was widened (`xlim: [700, 2300]`). No model setting changed.

The most informative single figure is `response_curves_panel`: it answers the
question (which habitat?) directly on the relative-intensity scale, with
uncertainty.

**9. Gate** (`eco-acceptance-gate`): **DONE.** Recovery passes at this design,
the fit converged, the goodness-of-fit and quadrature checks are clean, and
every report sentence matches the estimates. The answer is limited to relative
intensity, and the report says so.

**10. Deliver.** Summary for the researcher:

> After allowing for the fact that most records come from near roads, marten
> records are strongly associated with forest cover: about 1.6 times as many
> per 10 extra percentage points of forest (95% CI 1.4-1.7). They are also more
> frequent at higher elevation, about 1.15 times per 100 m (1.07-1.23). Without
> the road correction the elevation effect would have looked negative. The map
> shows where martens are relatively most likely, not how many there are; for
> that, use the camera data in the integrated model.

## (e) Estimates versus the known truth

| parameter (per SD) | truth | estimate | 95% CI | truth inside CI? |
|---|---|---|---|---|
| forest cover | 0.90 | 1.10 | 0.88 to 1.31 | yes (near the lower end) |
| elevation | 0.40 | 0.39 | 0.20 to 0.59 | yes |
| distance to road (bias) | -1.10 | -1.02 | -1.29 to -0.76 | yes |

On the natural scale: forest x1.56 per 10 points (truth x1.44), elevation x1.15
per 100 m (truth x1.15). The forest estimate is about 1.8 SE above the truth.
That is ordinary sampling variation for one dataset; the 200-replicate
recovery shows no bias (+0.9%).

## (f) Limits from VALIDATION.md that matter at this design

1. **The map is relative.** Prevalence, occupancy and abundance cannot be
   estimated from presence-only records. The camera data
   (`config_integrated.yaml`) are needed for those.
2. **Sampling bias must be measured.** The correction works only if distance
   to road captures where people looked. If effort also follows something else
   that correlates with habitat (e.g. trailheads at high elevation), habitat
   effects can be badly wrong with confident intervals (negative control:
   -73% bias, coverage 4.5%).
3. **Distance to road must not affect martens themselves.** If martens avoid
   roads, that avoidance is absorbed into the "bias" effect and missing from
   the habitat map.
4. **Independence.** Repeat records of the same animal or den make intervals too
   narrow. Here the block goodness-of-fit test shows no clustering (P = 0.69),
   but thinning duplicates to one per cell per season is still good practice.
   LGCP / spatial random effects are not implemented.
5. **Grid resolution.** Location errors (≤ 500 m) are below the 1 km cell, and
   the quadrature check is clean. Positional error larger than a cell was
   not validated.
6. If the camera data are added (integrated engine), absolute density needs
   about 60-100+ survey sites. With 45 stations, density and abundance are
   only weakly identified and should be read as an index (Known limit 2).

## (g) Commands to reproduce

From `.`:

```bash
# (optional) regenerate the example data and TRUTH.txt (seed 2025)
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/sdm-presence-only_make_example.log -- scripts/Rscript dev/sdm-presence-only/make_example.R

# planted-truth recovery at the example design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/sdm-presence-only_si_recover.log -- \
  scripts/Rscript core/run.R recover sdm-presence-only \
  --config modules/sdm-presence-only/example/recovery.yaml --nsim 200 --cores 16 \
  --out results/sdm-presence-only/si/recovery

# fit the example records (report + figures)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/sdm-presence-only_si_fit.log -- \
  scripts/Rscript core/run.R fit sdm-presence-only \
  --config modules/sdm-presence-only/example/config.yaml \
  --out results/sdm-presence-only/si/fit
```

Outputs: `results/sdm-presence-only/si/recovery/sdm-presence-only_recover_summary.csv`,
`results/sdm-presence-only/si/fit/sdm-presence-only_report.html`,
`results/sdm-presence-only/si/fit/figures/`, and
`results/sdm-presence-only/si/si_summary.json`.
