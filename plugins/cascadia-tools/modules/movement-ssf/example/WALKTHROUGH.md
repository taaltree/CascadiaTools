# Worked example: do elk choose forest and avoid high ground? (step selection, iSSA)

This example shows the full CascadiaTools workflow for GPS telemetry. The data are
**simulated**, so we know the true answer and can check every number against it.
The data, the story and the true values are in `README_data.md`.

## (a) The question, in the researcher's words

> "I collared 10 cow elk on their summer range and got a GPS fix every 2 hours from
> 1 May to 15 June. Do elk pick forest and stay off high ground when they move? How
> strongly? Do some elk do it more than others? And how big are their summer home ranges?"

## (b) The files the researcher has

`elk_gps_2024.csv`: one row per fix attempt, as the collar company delivers it (5,410 rows).
Failed fixes are rows with empty coordinates.

```
"animal_id","timestamp","utm_e","utm_n","dop","fix_type"
"E01","2024-05-01 00:00:23",511014.9,5007670.1,4,"3D"
"E01","2024-05-01 02:00:16",511255.5,5007538.4,2.3,"3D"
"E01","2024-05-01 04:01:33",511324.1,5007658.1,1.4,"3D"
"E01","2024-05-03 04:01:55",,,,"no fix"
```

`forest_cover.tif` and `elevation.tif`: two habitat maps (GeoTIFF rasters, 50-m cells,
UTM zone 11N). A raster has no rows to show. Forest cover runs from 0 to 100 % (mean 55.9,
SD 24.7). Elevation has a mean of 1,250 m (SD 220 m).

## (c) What the researcher types to Claude

> I have GPS collar data from 10 elk (elk_gps_2024.csv, one row per fix attempt, fixes
> every 2 hours, coordinates in UTM zone 11N) and two habitat maps of my study area:
> forest_cover.tif (% canopy) and elevation.tif (m). Please tell me whether elk select
> forest and avoid high elevation as they move, how strong that is, whether individual
> elk differ, and how large their summer home ranges are.

## (d) What the agent does, step by step (`skills/eco-analyze`)

**1. Project folder.** The raw files are copied and never edited. All work, with every
command, is logged in one folder.

**2. Data intake.** The agent reads the files and describes them back:
- 10 elk (E01-E10), 5,410 fix attempts from 1 May 00:00 to 15 June 00:02 UTC.
- 5,000 fixes have coordinates (494-507 per elk). 410 attempts (7.6%) failed.
  59 fixes are 2D and 4,941 are 3D; DOP median 2.0 (range 0.5-7.8).
- Fixes are about 2 h apart, with 0-2 min of jitter. No fix is duplicated.
- All fixes fall inside both rasters (easting 502.6-514.1 km, northing 5,002.6-5,014.2 km).
- `dop` and `fix_type` are recorded but are not used in the model.

**3. Choose the module.** Regular GPS fixes from several animals plus continuous habitat
maps. The question is about selection while moving, individual differences and home
ranges. That is the `movement-ssf` module (integrated step-selection analysis, iSSA).
`movement-hmm` would be the choice for behavioural states; it is not the question here.

**4. Plan.** Written before any fit:
- Estimands: population-mean selection for forest and elevation (log relative selection
  strength per 1 SD), and as relative selection strength (RSS) per 10 percentage points
  of forest and per 100 m of elevation. Also the among-elk SD of each, the selection-free
  step length and turning, and 95% home-range areas per elk (AKDE).
- Model: 20 available steps per used step. Random slopes per elk. A range-residency
  (home-attraction) term is included but is a nuisance and is not interpreted.
- Checks: model convergence, used-habitat calibration on 20% held-out steps, and an
  RSF and KDE/MCP comparison only for context.

**5. Planted truth (the gate).** Before touching the real data, the agent simulates 200
datasets at this exact design (10 elk, 540 two-hour steps, 8% failed fixes) with known
values, fits each one, and checks the answers. One fit takes about 37 s; all 200
converged. Coverage needs to be at least 0.904 and bias within 10% (or 3 Monte Carlo SE).

| estimand | truth | mean estimate | rel. bias | 95% CI coverage | power | gate |
|---|---|---|---|---|---|---|
| forest selection, mean | 0.80 | 0.791 | -1.1% | 0.920 | 1.00 | PASS |
| elevation selection, mean | -0.50 | -0.500 | 0.0% | 0.935 | 0.955 | PASS |
| forest selection, SD among elk | 0.40 | 0.384 | -4.1% | 0.975 | | PASS |
| elevation selection, SD among elk | 0.30 | 0.294 | -1.9% | 0.970 | | PASS |
| RSS forest per SD | 2.23 | 2.23 | 0.0% | 0.920 | | PASS |
| RSS elevation per SD | 0.607 | 0.612 | +0.9% | 0.935 | | PASS |
| step-length shape | 2.00 | 1.999 | 0.0% | 0.930 | | PASS |
| step-length scale (m) | 120 | 120.3 | +0.2% | 0.930 | | PASS |
| mean step (m per 2 h) | 240 | 240.3 | +0.1% | 0.955 | | PASS |
| turning concentration | 1.00 | 1.001 | +0.1% | 0.955 | | PASS |
| home attraction (nuisance) | 0.25 | 0.294 | +17.6% | 0.720 | | FAIL |

What the agent tells the researcher: *"At your design the method recovers all ten
quantities we will interpret: bias at most 4% and interval coverage 92-98%. With 10 elk it
detects an elevation effect of this size 96% of the time. The range-residency term is
overestimated, as the module's validation documents, so we will not interpret its size."*

**6. Fit.** The agent writes `config.yaml` (columns, CRS, raster files, meaningful units:
10 percentage points of forest, 100 m of elevation) and runs the fit on the cluster.

**7. Check.**
- Convergence: glmmTMB converged; the Hessian is positive definite and all SEs are finite.
  (This is a likelihood fit, so there is no R-hat.)
- Data used: 5,000 located fixes became 3,918 two-hour steps, each with 20 available steps.
- Goodness of fit (used-habitat calibration on 784 held-out steps): the observed mean
  forest cover at used steps is 66.9% vs the model's 67.2% (95% interval 66.6 to 67.9).
  Observed mean elevation is 1,142.6 m vs 1,142.2 m (1,138.6 to 1,146.4). Both are
  inside the interval, so the model reproduces where elk actually went.

**8. Interpret.** Quoted from the report (`results/movement-ssf/si/fit/movement-ssf_report.html`):

> "Elk selected locations with higher forest cover at the scale of one 2-h step
> (selection coefficient 0.78 per SD, 95% CI 0.42 to 1.13). [...] Per 10 percentage
> points of forest cover, RSS = 1.37 (1.19 to 1.58)."

> "Elk avoided locations with higher elevation at the scale of one 2-h step (selection
> coefficient -0.50 per SD, 95% CI -0.81 to -0.19). [...] Per 100 m of elevation,
> RSS = 0.80 (0.69 to 0.92)."

> "Individuals differ: the among-animal SD of the selection coefficient is 0.53 (95% CI
> 0.31 to 0.92) [...] so an estimated 93% of elk in this population would show positive
> selection (log-RSS > 0) for forest cover."

> "Home ranges (95% AKDE, which accounts for autocorrelation): median 16.9 km^2 across
> 10 elk (range 4.4 to 49.5 km^2)."

> "These are associations between where animals stepped and the mapped covariates,
> conditional on what was reachable; they are not proof that the covariate causes the
> choice."

In plain words: of two places an elk could equally reach, the one with 10 points more
forest is about 1.4 times as likely to be chosen. The one 100 m higher is about 0.8 times
as likely. Most elk prefer forest; only about 12% would be expected to prefer high ground.

Figures (in `results/movement-ssf/si/fit/figures/`, PNG and PDF):

| figure | caption (short) |
|---|---|
| `rss_forest.png` | RSS for forest cover vs the landscape average (55.9%), typical elk, 95% band, log scale |
| `rss_elevation.png` | RSS for elevation vs the landscape average (1,250 m), typical elk, 95% band |
| `selection_coefficients.png` | **Key figure.** Population mean and each elk's selection coefficient, 95% intervals |
| `movement_kernel.png` | Step lengths and turn angles: observed, tentative and selection-free distributions |
| `used_available.png` | Forest and elevation at used vs available step ends (the raw selection signal) |
| `uhc.png` | Used-habitat calibration on held-out steps: observed use falls inside the model envelope |
| `home_ranges.png` | 95% AKDE home range of each elk with its confidence contours and GPS fixes |
| `hr_area.png` | 95% home-range area per elk: AKDE with CI, compared with KDE and MCP |

**9. Gate (`eco-acceptance-gate`).** Verdict: **DONE**. The planted-truth check passes for
every interpreted quantity at this design, the fit converged, the model check passes, and
each report sentence matches the estimates table. The one failing row (home attraction)
is a nuisance term that the report does not interpret.

**10. Deliver.** A short summary, the report, the figures and a methods paragraph (the
report's Methods section).

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% CI) | truth inside? |
|---|---|---|---|
| forest selection, mean (per SD) | 0.80 | 0.78 (0.42 to 1.13) | yes |
| elevation selection, mean (per SD) | -0.50 | -0.50 (-0.81 to -0.19) | yes |
| forest selection, SD among elk | 0.40 | 0.53 (0.31 to 0.92) | yes |
| elevation selection, SD among elk | 0.30 | 0.42 (0.22 to 0.80) | yes |
| RSS forest per 10 points | 1.38 | 1.37 (1.19 to 1.58) | yes |
| RSS elevation per 100 m | 0.80 | 0.80 (0.69 to 0.92) | yes |
| step-length shape | 2.00 | 2.10 (2.01 to 2.19) | no (just) |
| step-length scale (m) | 120 | 117 (111 to 123) | yes |
| mean step (m per 2 h) | 240 | 245 (239 to 251) | yes |
| turning concentration | 1.00 | 1.02 (0.96 to 1.07) | yes |
| home attraction (nuisance) | 0.25 | 0.31 (0.24 to 0.37) | yes |

11 of 12 intervals contain the truth. With 95% intervals, about one miss in 12 is
expected by chance. The step-shape miss is small (2.10 vs 2.00), and the planted-truth
runs show no bias for it (coverage 0.93). Home ranges have no single true value in this
simulation (see `README_data.md`).

## (f) Limits that matter at this design (from `VALIDATION.md`)

- The home-attraction (range-residency) coefficient is overestimated (+18% here, coverage
  0.72). It must be in the model, but its size should not be interpreted.
- With 10 elk the among-elk SDs are imprecise (intervals span about a factor of 3). The
  Wald interval for an SD breaks down when the estimate is near zero.
- AKDE intervals were too narrow in validation (coverage 0.85 with about 13 effective range
  crossings; here the median is 15). Treat home-range intervals as somewhat optimistic.
- Fix loss is assumed unrelated to habitat. If collars fail more under dense canopy,
  forest selection is underestimated. Check fix success by habitat with real collars.
- Only continuous, smooth covariates were validated. Categorical maps, habitat-dependent
  movement and habitat-biased fix loss were not.
- The RSF in the report is for comparison only. It answers a different (home-range scale)
  question and its intervals are too narrow for GPS data.
- A weaker effect (for example -0.25 per SD of elevation) would need about 30 or more
  elk for 80% power.

## (g) Commands to reproduce

Run from `ecotools/`. `scripts/ejob` submits to SLURM and waits.

```
# example data (already in example/; seed 2024, fix loss seed 99)
scripts/ejob -c 1 -m 8 -- scripts/Rscript dev/movement-ssf/make_example.R

# planted-truth recovery at the example design (about 20 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/movement-ssf_si_recover.log -- scripts/Rscript core/run.R recover movement-ssf --config modules/movement-ssf/example/recovery.yaml --nsim 200 --cores 16 --out results/movement-ssf/si/recovery

# fit, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/movement-ssf_si_fit.log -- scripts/Rscript core/run.R fit movement-ssf --config modules/movement-ssf/example/config.yaml --out results/movement-ssf/si/fit

# redraw figures from the saved fit after editing the figures: block
scripts/Rscript core/run.R figures movement-ssf --fit results/movement-ssf/si/fit/movement-ssf_fit.rds --config modules/movement-ssf/example/config.yaml
```

Outputs: `results/movement-ssf/si/recovery/` (gate table and report),
`results/movement-ssf/si/fit/` (report, estimates CSV, figures),
`results/movement-ssf/si/si_summary.json` (machine-readable summary).
