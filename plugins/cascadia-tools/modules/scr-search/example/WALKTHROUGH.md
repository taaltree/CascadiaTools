# Worked example: how many fishers, and do they favour forest? (search-encounter SCR)

This example shows the full CascadiaTools workflow for scat-detection-dog data. The data are
**simulated**, so we know the true answer and can check every number against it.
`README_data.md` describes the files and lists the true values.

## (a) The question, in the researcher's words

> "Our dog teams walked 36 routes in the western Cascades last summer and picked up fisher scats.
> The lab has identified which fisher left each scat. How many fishers live in the area, and at
> what density? Are there more fishers where there is more forest? And do our experienced dog
> teams really find more scats than the newer teams?"

## (b) The files the researcher has

`search_tracks_gps.csv`: the GPS log, one row per fix (4,715 rows).

```
"track_id","gps_time","utm_e","utm_n"
"T001","2025-06-12T07:34:35Z",568304.8,4883730
"T001","2025-06-12T07:35:35Z",568266.6,4883762.3
"T001","2025-06-12T07:36:35Z",568216.6,4883760.2
"T001","2025-06-12T07:37:35Z",568182.5,4883796.7
```

`track_field_sheet.csv`: one row per route (36 rows).

```
"track_id","date","session","dog_team","experienced_team","km_walked"
"T001","2025-06-12",1,"Pepin & Sara",1,9.2
"T002","2025-06-11",1,"Rue & Lena",0,6.3
"T003","2025-06-17",1,"Rue & Lena",0,8.65
"T004","2025-06-18",1,"Odin & Kim",0,6.35
```

`fisher_scats.csv`: one row per scat collected (116 rows; blank `fisher_id` = failed genotype).

```
"scat_id","track_id","date","utm_e","utm_n","fisher_id","genotype_qc"
"PEPE-25-0001","T001","2025-06-12",564731.7,4885101.8,"F013","pass"
"PEPE-25-0002","T001","2025-06-12",565451.3,4884393.7,"F013","pass"
"PEPE-25-0003","T001","2025-06-12",565593.6,4884254.2,"F014","pass"
"PEPE-25-0004","T001","2025-06-12",566055.2,4884190.9,"F016","pass"
```

`habitat_grid.csv`: 700 m habitat pixels over the searched area plus 6 km (1,600 rows).

```
"utm_e","utm_n","forest_pct"
555350,4866350,58.6
556050,4866350,53.3
556750,4866350,48.6
557450,4866350,44.6
```

(The snippets are the first five lines of each file.)

## (c) What the researcher types to Claude

> I ran a scat-detection-dog survey for fishers in the western Oregon Cascades in summer 2025:
> 36 search routes in three monthly sessions, each logged by GPS (search_tracks_gps.csv), with a
> field sheet per route saying which dog team searched and whether it was one of our experienced
> teams (track_field_sheet.csv). The scats and their genotypes are in fisher_scats.csv (blank
> fisher_id = failed genotype). I also have forest cover on a 700 m grid (habitat_grid.csv).
> Please estimate fisher density and abundance, tell me whether density is related to forest
> cover, and whether experienced teams find more scats per km. I'd like a report and figures
> for a paper.

## (d) What the agent does, step by step (`skills/eco-analyze`)

**1. Project folder.** Raw files are copied, never edited. All work goes in one folder, with a
log of every command.

**2. Data intake.** The agent reads the four files and describes them back:
- 36 routes in 3 sessions (12 each), searched 10 June to 14 August 2025. 4,715 GPS fixes, one per
  minute. Routes are 4.0-10.7 km; 246.7 km in total (85.8, 82.2 and 78.8 km per session).
- Four dog teams. Experienced teams: 15 routes, 100.5 km. Newer teams: 21 routes, 146.2 km.
- 116 scats. 35 (30%) failed genotyping and are dropped. They were still part of the search
  effort, which is counted from the GPS tracks.
- 81 genotyped scats from 22 fishers. Scats per fisher: 1 to 12 (median 3). 20 of 22 fishers
  were found in more than one 1-km cell (spatial recaptures); 14 in more than one session.
- Every scat has a route ID that is on the field sheet; the module's data check raised no warnings
  (no GPS jumps, no scat far from its route, no ID typos).
- Habitat grid: 40 x 40 pixels of 700 m (784 km2), forest cover 32.6-85.5% (mean 55.0, SD 12.0).

**3. Choose the module.** Individuals identified from scats found along recorded search routes,
not at fixed traps: the `scr-search` module (search-encounter spatial capture-recapture). Fixed
traps (`scr`) would ignore where and how much the dogs searched.

**4. Plan.** Written before any fit:
- Estimands: density (fishers per km2) and abundance in the state space; spatial scale sigma and
  95% home-range area; scats found per km searched (lambda0); the forest effect on density
  (rate ratio per 1 SD of forest cover); the experienced-team effect on scats per km.
- Model: GPS tracks cut into 1 km cells; effort = km searched in each cell; Poisson number of
  scats per cell and route; half-normal decline with distance from each activity centre.
  Density log-linear in forest cover; detection depends on experienced team (0/1).
- Engine: Bayesian (Nimble), 3 chains x 40,000 iterations, 10,000 burn-in. Maximum-likelihood
  check of the same model with secr (`compare_engines: [secr_cells]`).
- Checks: R-hat < 1.01 and bulk ESS >= 400; posterior predictive checks; cell size vs sigma;
  buffer vs sigma.

**5. Planted truth at this design.** The agent simulated 100 data sets with the same design
(16 x 16 km searched, 6 km buffer, 36 routes of 4-10 km over 3 sessions, 40% of routes on roads
that run moderately through forest, 1 km cells) and the true values of this example. It fitted
each one and checked the answers (`recovery.yaml`). One fit takes about 4 minutes, so 100
replicates were used (the rule: 100 when a fit takes more than a minute).
Result: **all 8 estimands PASS**; 100 of 100 fits ran and converged (largest R-hat 1.008).

| quantity | truth | mean estimate | relative bias | 95% CI coverage | power |
|---|---|---|---|---|---|
| density (per km2) | 0.100 | 0.105 | +4.8% | 0.91 | - |
| abundance (state space) | 78.6 | 82.1 | +4.3% | 0.94 | - |
| lambda0 (scats per km at centre) | 0.200 | 0.206 | +2.9% | 0.95 | - |
| sigma (m) | 1500 | 1524 | +1.6% | 0.96 | - |
| 95% home range (km2) | 42.4 | 44.1 | +4.2% | 0.96 | - |
| scats from a route through the centre | 0.752 | 0.775 | +3.1% | 0.97 | - |
| forest effect (log, per SD) | 0.50 | 0.511 | +2.2% | 0.97 | 0.57 |
| experienced-team effect (log) | 0.60 | 0.544 | -9.4% | 0.93 | 0.66 |

The pass rule is |bias| <= 10% of truth (or 3 Monte Carlo SE), coverage >= 0.885 at 100
replicates, and >= 90% converged. Abundance truth (78.6) is the average over replicates; the
one example population has 63 activity centres. Across replicates, N-hat / N had median 1.02
(5-95%: 0.71-1.42). Power is the share of 95% intervals that exclude zero: at this design the
forest effect is found about 6 times in 10 and the team effect 2 times in 3. The team effect
passed but is the closest to the bias limit (-9.4%; limit 10%).

**6. Fit.** `core/run.R fit` with `config.yaml` (about 4 minutes for the Nimble fit, plus the
secr check).

**7. Checks.**
- Convergence: max R-hat 1.001; min bulk ESS 3,250 (needs < 1.01 and >= 400). Traces mix well.
- Posterior predictive checks (200 draws; extreme if < 0.025 or > 0.975): fishers identified
  0.49, total scats 0.45, cells per fisher 0.45, Freeman-Tukey per cell 0.81, dispersion of
  scats per fisher 0.30. No lack of fit; no sign of latrine clustering.
- Design: cells are 0.64 x sigma (should be <= 1); the buffer is 3.8 x sigma (should be >= 3-4).
- Precision note from the fit: the posterior SD of density is 20% of its mean. Posterior
  medians (density 0.0697, N 54) are reported next to the means (0.0704, 54.7).
- Second engine: maximum likelihood (secr, same cell model) gives density 0.068 (0.044-0.104),
  N 53.3 (39.7-77.2), sigma 1,554 m, forest effect 0.94, team effect 0.53. Agreement is close.

**8. Interpretation.** Excerpt from the report (`results/scr-search/si/fit/scr-search_report.html`):

> "Density was estimated at 0.070 (95% CrI 0.045 to 0.10) animals per km², i.e. 7.04 per
> 100 km2. Over the 784 km² state space (the searched area plus a buffer), that is 54.7 (95% CrI
> 38.0 to 76.0) animals, of which 22 (40%) were identified from samples."
>
> "The spatial scale of space use, sigma, was 1561 (95% CrI 1345 to 1828) m. ... a home-range
> area of about 46 (95% CrI 34 to 63) km²; treat this as an index of ranging scale."
>
> "Detection: 'Experienced dog team' (1 vs 0) multiplied the number of samples found per km
> searched by 1.65 (95% interval 1.04 to 2.58) ... This is a detection (search efficiency)
> effect, not a difference in where animals live."
>
> "Density and Forest cover (%): each 1 SD (12 units) higher Forest cover (%) is associated with
> 155% higher density (rate ratio 2.55, 95% interval 1.60 to 4.20). Expected density is 0.015 per
> km² (0.0036 to 0.037) where Forest cover (%) = 39.6 (10th percentile of the state space) and
> 0.16 (0.10 to 0.23) where it = 71.6 (90th percentile). ... This is a correlation between
> density and habitat across the state space, not proof that the habitat causes higher density.
> Because search effort was modelled explicitly (km searched per cell), more samples in heavily
> searched areas are not mistaken for higher density."

The agent checked each sentence against the estimates table: exp(0.937) = 2.55; exp(0.501) =
1.65; 0.18 x 1.65 = 0.30 scats per km for experienced teams; 22 / 54.7 = 40%.

Figures (`results/scr-search/si/fit/figures/`, PDF and PNG):

| figure | caption |
|---|---|
| `effort_map` | km searched in each 1 km cell, GPS tracks, and the 81 genotyped scats. Unsearched ground has no detector. |
| `density_surface` | Posterior mean realised density (fishers per km2), with search tracks and the estimated activity centres of the 22 identified fishers. |
| `density_expected_surface` | Expected density surface implied by the forest-cover model. |
| `detection_map` | Probability that a fisher living in each pixel left at least one scat that was found. Low values: density there rests on the habitat model. |
| `density_forest_pct` | Expected density vs forest cover, with 95% and 50% bands. |
| `detection_function` | Scats found per km searched vs distance from the activity centre, by team; dashed line = cell size. |
| `detection_experienced_team` | Scats per km at the activity centre: newer vs experienced teams (95% intervals). |
| `abundance_posterior` | Posterior of the number of fishers in the 784 km2 state space; dashed = 95% interval; dotted = 22 identified. |
| `captures_per_individual` | Observed number of fishers with 1, 2, ... scats vs posterior predictive intervals. |
| `gof` | Posterior predictive checks (panels A-E). |
| `coefficients` | Forest and team effects on the log scale (dashed line = no effect). |
| `engine_comparison` | Nimble vs maximum likelihood (secr) on the same data, every estimate. |
| `trace` | MCMC traces, one colour per chain. |
| `posterior` | Posterior densities of the reported quantities. |

The single most informative figure is `density_forest_pct`: it answers the habitat question in
fishers per km2 across the forest cover found in the study area.

Presentation fixes (config `figures:` block only): `density_expected_surface` is drawn
one-and-a-half columns wide, because at single width the long legend title squeezed the map and
the x tick labels overlapped; `coefficients` uses x limits -0.1 to 1.6 so the last tick label
is not cut off.

**9. Gate verdict: DONE.** The planted-truth check passed at this design (8/8), the fit
converged, the fit checks passed, a second engine agrees, and every number comes from a
rerunnable command.

## (e) Estimates versus the known truth

| quantity | estimate (95% CrI) | truth | inside CrI? |
|---|---|---|---|
| density (per km2) | 0.070 (0.045 to 0.101) | 0.10 expected; 0.080 realised (63 / 784) | yes (both) |
| abundance (state space) | 54.7 (38 to 76) | 63 | yes |
| sigma (m) | 1561 (1345 to 1828) | 1500 | yes |
| 95% home range (km2) | 46.2 (34.1 to 62.9) | 42.4 | yes |
| lambda0, newer team (scats per km) | 0.180 (0.113 to 0.274) | 0.20 | yes |
| scats from a route through the centre | 0.70 (0.47 to 1.01) | 0.752 | yes |
| forest effect (log, per SD) | 0.94 (0.47 to 1.43) | 0.50 | yes (near the lower end) |
| team effect (log) | 0.50 (0.04 to 0.95) | 0.60 | yes |

All eight truths are inside the 95% intervals. The forest effect in this one data set is about
twice the truth, but the truth is still inside the interval. In the 100 simulated data sets the
forest estimate was unbiased (+2.2%), so this is chance in this realisation, not a flaw.
Density is under the expected 0.10 partly because this population happened to be small (63
centres = 0.080 per km2, versus 78 on average).

## (f) Limits from VALIDATION.md that matter at this design

- **Sparse data and the posterior mean.** With 20-35 fishers identified (here 22) and a density
  CV near 0.2, the posterior MEAN of density and N ran 9-13% high in some independent
  simulations, while intervals kept nominal coverage. Here the bias was +4.8%, but quote the
  posterior median and the secr estimate alongside the mean.
- **Abundance is for the whole state space** (searched area + 6 km buffer, 784 km2). Away from
  the routes, density comes from the forest model (see `detection_map`). Quote density, not N
  alone.
- **Effort must be recorded.** Ignoring effort more than doubled the apparent forest effect in
  simulations where routes followed roads through forest. Here effort is modelled from the GPS
  tracks.
- **Discretisation.** 1 km cells (0.64 sigma) inflate sigma by about 3%. Small here.
- **Power.** At this effort (~250 km) the forest effect (true rate ratio 1.65) is detected in
  57% of simulated surveys and the team effect in 66%. About 400 km of searching is needed for
  ~75% power on the forest effect.
- **Assumptions not tested by the data**: a closed population over June-August, one sigma for
  both sexes, independent scats (no latrines), correct genotypes.

## (g) Commands to reproduce

```bash
# how the example data were generated (seed 2025); not needed, the CSVs are already in example/
scripts/ejob -c 1 -m 4 -t 0:20:00 -- scripts/Rscript dev/scr-search/make_example.R
# planted-truth recovery at this design (about 45 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/scr-search_si_recover.log -- scripts/Rscript core/run.R recover scr-search \
  --config modules/scr-search/example/recovery.yaml --nsim 100 --cores 16 --out results/scr-search/si/recovery
# fit, report and figures (about 6 min)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/scr-search_si_fit.log -- scripts/Rscript core/run.R fit scr-search \
  --config modules/scr-search/example/config.yaml --out results/scr-search/si/fit
# redraw figures only, after editing the figures: block
scripts/ejob -c 2 -m 8 -t 0:20:00 -- scripts/Rscript core/run.R figures scr-search \
  --fit results/scr-search/si/fit/scr-search_fit.rds --config modules/scr-search/example/config.yaml
```
