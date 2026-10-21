# Worked example: mule deer density from line transects (distance module)

This walkthrough follows one researcher from question to report. The data are
**simulated**, so the true answer is known and the analysis can be checked against it
(see `README_data.md`). Every number below comes from the runs listed in section (g).

## (a) The question

> "We walked 36 transects on the deer winter range in January and measured how far
> each deer group was from the line. How many mule deer are on the winter range, in each
> management unit, and are there more deer where there is more shrub forage and farther
> from the highway?"

## (b) The files the researcher has

`deer_observations.csv` (146 rows, one per deer group seen):

```
"transect_id","date","observer","perp_distance_m","group_size"
"MD01","2025-01-19","JK",61,4
"MD01","2025-01-19","JK",83,2
"MD01","2025-01-19","JK",88,1
"MD01","2025-01-19","JK",96,5
```

`deer_transects.csv` (36 rows, one per transect, including MD18 where no deer were seen):

```
"transect_id","stratum","length_km","shrub_cover_pct","road_dist_km","x_km","y_km","observer"
"MD01","south",4.06,24.9,5.32,14.5,4.5,"JK"
"MD02","south",2.46,28.8,2.23,6.5,8.5,"AB"
"MD03","north",3.96,21.4,6.17,5.5,12.5,"JK"
"MD04","south",3.03,17.3,4.29,17.5,7.5,"JK"
```

`study_area_grid.csv` (360 rows, one per 1 km2 cell of the range):

```
"x","y","area_km2","stratum","shrub_cover_pct","road_dist_km"
0.5,0.5,1,"south",33.1,1.54
1.5,0.5,1,"south",36,2.06
2.5,0.5,1,"south",38.9,2.57
3.5,0.5,1,"south",41.4,3.09
```

## (c) What the researcher types to Claude

> I have a mule deer line-transect survey in modules/distance/example: a file of deer
> groups seen with the perpendicular distance from the line (m), group size, observer and
> date; a file of the 36 transects with their length, management unit (north or south,
> 180 km2 each), shrub cover and distance to the highway; and a 1 km2 grid of the whole
> 360 km2 winter range with the same two habitat values. Please estimate how many deer
> there are per km2 and in total, for each unit and overall, and tell me whether deer are
> more common where shrub cover is higher and farther from the road. Check that the method
> works for a survey like mine before you trust it, and give me a report and figures I can
> put in a paper.

## (d) What the agent does, step by step (eco-analyze)

**1. Project folder.** For a real researcher the agent would copy the files into
`projects/<name>/data_raw/` and keep a `log.md`. In this example the files stay in
`modules/distance/example/` and the outputs go to `results/distance/si/`.

**2. Data intake (eco-data-intake).** The agent reads the files and describes them back:
- 146 deer groups (486 deer) on 36 lines, 1.9-4.1 km long, 114.9 km of line in total.
- 97 groups were in the south unit and 49 in the north unit.
- One line (MD18) had no deer. It is in the transects file, as it must be.
- Distances run from 0 to 461 m. Only one is beyond 400 m.
- Group sizes are 1-9 (mean 3.33).
- The surveyed lines cover shrub cover 13.9-45.4 % and road distance 0-10.8 km.
  The grid has shrub cover 3-55 %, so map predictions at the extremes are extrapolations.
- The histogram of distances has no pile-up at zero and no heaping at round numbers.

**3. Choose the module.** `distance` fits: walked lines, measured distances, groups,
strata, habitat per transect. It is not `nmixture` (that needs counts without distances)
and not `scr` (that needs individually identified animals).

**4. Plan (eco-analysis-plan).**
- Estimands: deer density (per km2) and abundance, overall and per unit; detection
  probability; the change in density with shrub cover and with road distance.
- Detection function: half-normal and hazard-rate keys with up to 2 cosine adjustment
  terms, chosen by AIC and model-averaged (Distance::ds).
- Truncation at 400 m, fixed **before** fitting.
- Density and abundance: design-based (R2 encounter-rate variance), separately per unit.
- Habitat effects: hierarchical distance sampling (unmarked::distsamp), with log density
  linear in shrub cover and road distance, both standardised.
- Checks: Cramer-von Mises test and Q-Q plot for the distances; a Freeman-Tukey parametric
  bootstrap (200 sets) for the hierarchical model.

**5. Planted truth (eco-planted-truth).** The agent simulated 200 surveys like this one
(`recovery.yaml`) and fitted each one the same way:
- 36 lines of 3.19 km in a 360 km2 area with two units, truncation 400 m.
- Groups of mean size 3.4, half-normal detection (sigma 146 m).
- Group density 3.05 per km2 at an average transect, and a habitat effect of 0.49 per SD
  (the true shrub effect).

Each fit took about 9 s. All 200 converged, with a median of about 145 detections, like the real survey.
The gate: |bias| at most 10 % (or 3 MC SE), and coverage at least 0.904.

| estimand | truth | mean estimate | rel. bias | 95% coverage | power | gate |
|---|---|---|---|---|---|---|
| density_per_km2 (deer) | 11.67 (realised) | 11.54 | -1.0 % | 0.975 | | PASS |
| abundance | 4,193 | 4,149 | -1.0 % | 0.975 | | PASS |
| density per unit (S1 / S2) | 11.62 / 11.71 | 11.47 / 11.61 | -1.5 % / -0.5 % | 0.970 / 0.985 | | PASS |
| group_density_per_km2 | 3.43 | 3.40 | -0.9 % | 0.995 | | PASS |
| mean_group_size | 3.40 | 3.39 | -0.2 % | 0.925 | | PASS |
| detection_prob (P_a) | 0.455 | 0.465 | +2.4 % | 0.955 | | PASS |
| esw_m | 181.9 | 186.2 | +2.4 % | 0.955 | | PASS |
| sigma_m (142 fits that chose a plain half-normal) | 146 | 145.8 | -0.1 % | 0.972 | | PASS |
| density_beta_habitat | 0.49 | 0.488 | -0.4 % | 0.950 | 1.00 | PASS |
| rate_ratio_habitat | 1.632 | 1.635 | +0.2 % | 0.950 | 1.00 | PASS |
| density_at_average_site (deer, HDS) | 10.37 | 10.32 | -0.5 % | **0.885** | | **FAIL** |
| group_density_at_average_site (HDS) | 3.05 | 3.04 | -0.4 % | **0.890** | | **FAIL** |

In one sentence: **the method recovers density, abundance, detection and the habitat
effect at this design (bias within 2.4 %, coverage 0.95-0.995, power 1.00 for the shrub
effect), but the interval for "density at an average transect" from the hierarchical model is
a little too narrow (coverage 0.885 against a bound of 0.904), although its estimate is unbiased.**
This is Known limit 9 in VALIDATION.md. That interval uses only the detection key with the best AIC.
In the 142 replicates where a plain half-normal model was chosen, coverage was 0.92. In the other
58 it was 0.81. The agent tells the researcher this. The design-based density
and the habitat effect are not affected.

**6. Fit.** `core/run.R fit distance --config modules/distance/example/config.yaml`
(engine `auto`: Distance::ds for density and detection, unmarked::distsamp for the habitat
effects). It runs in about 10 minutes, most of it the 200-set bootstrap for the goodness-of-fit check.

**7. Check (eco-convergence, eco-model-checking).**
- Both optimisers converged. All standard errors are finite and no estimate is on a boundary.
  The fit raised no warnings.
- Detection function: half-normal was chosen by AIC. Hazard-rate was 1.9 AIC units worse
  (Akaike weight of half-normal 0.72), and the two were model-averaged.
- Cramer-von Mises test of the distances: W = 0.075, p = 0.72. The Q-Q plot follows the line.
- Hierarchical model: half-normal chosen (hazard-rate +1.4 AIC). Freeman-Tukey parametric
  bootstrap p = 0.32 (0.31 in an earlier identical run; the bootstrap is random). No evidence of lack of fit.
- The hierarchical density over the grid (12.2 per km2) agrees with the design-based estimate (12.0).

**8. Interpret (eco-interpretation).** Excerpt from the report (`distance_report.html`), with each
sentence checked against the estimates table:

> A group within 400 m of the line had probability 0.44 (95% CI 0.38 to 0.50) of being
> detected. Equivalently, the effective strip half-width is 175 m.
>
> Density: 12.0 (95% CI 9.24 to 15.5) individuals per km2, which over the 360 km2 study area
> is 4,314 (95% CI 3,328 to 5,593) individuals. This is 3.59 (95% CI 2.82 to 4.58) groups per
> km2 times a mean group size of 3.34 (95% CI 3.12 to 3.56).
>
> By stratum: north 7.97 (95% CI 5.68 to 11.2) per km2; south 16.0 (95% CI 11.6 to 22.0) per km2.
>
> Each increase of 10 percentage points in shrub cover multiplies expected density by 1.89
> (95% CI 1.56 to 2.30) ... This is an association from an observational survey, not proof
> that shrub cover causes the difference.
>
> Each increase of 3 km in distance from road multiplies expected density by 1.40 (95% CI
> 1.12 to 1.74).

Check: exp(0.624 x 10 / 9.78) = 1.89, and 3.59 x 3.34 x 360 = 4,314.

Figures (`results/distance/si/fit/figures/`, PNG + PDF):

| id | caption |
|---|---|
| `detection_function` | Fitted half-normal detection function with 95% and 50% bands over the histogram of the 145 distances, scaled to detection probability. |
| `qq_gof` | Q-Q plot of the detection function. Cramer-von Mises W = 0.075, p = 0.72. |
| `density_shrub_cover_pct` | Expected deer density versus shrub cover (hierarchical model; road distance at its mean), with 95%/50% bands. Points are rough per-transect densities. |
| `density_road_dist_km` | Expected deer density versus distance from the road (shrub cover at its mean). |
| `effects_panel` | Both conditional effects side by side (A shrub, B road). |
| `coefficients` | Density ratio per 1 SD of each covariate (shrub 9.8 %, road 2.7 km) with 95% CIs. |
| `density_by_stratum` | Design-based density in the north unit, south unit and whole area with 95% CIs. |
| `density_map` | A: predicted deer density per km2 over the 360 km2 grid. B: its SD. White dots are the transects. |

The only presentation change was in the config's `figures:` block. The coefficient plot's axis title was
cut off, so it now gives the SD units and the plot is wider.

**9. Gate (eco-acceptance-gate): OPEN.**
- Planted truth passes for every estimand that answers the question: density, abundance,
  density per unit, detection and the habitat effect.
- Convergence and goodness of fit pass.
- Independent cross-checks (hand-coded likelihoods) and negative controls (evasive movement,
  heaping) are in VALIDATION.md.
- One reported quantity fails coverage: the hierarchical model's density at an average transect
  (0.885). The bands on the two conditional-effect curves and the map are built from the same intercept.
- So the verdict is OPEN, not DONE. Report the design-based density and abundance, and the
  rate ratios, as validated. Treat the bands on the curves and the map as slightly too narrow.
- In this data set the hierarchical model chose the half-normal, where coverage in the
  simulation was 0.92. The prepared fix (VALIDATION.md, limit 9) would close the gap.

**10. Deliver.** Summary for the researcher:

> About 12.0 mule deer per km2 used the winter range (95% CI 9.2-15.5), or about 4,300 deer
> (3,300-5,600) on the 360 km2. Density was about twice as high in the south unit (16.0 per
> km2, 11.6-22.0) as in the north (8.0, 5.7-11.2). Deer were more common where shrub cover was
> higher: 1.89 times more per 10 extra percentage points of shrub cover (1.56-2.30). They were also
> more common farther from the highway: 1.40 times more per 3 km (1.12-1.74). These are associations,
> not proof of cause. Observers saw about 44 % of the groups within 400 m of the line, and this has been
> corrected for. The estimates assume that every group on the line was seen and that deer did not
> move away before they were spotted.

## (e) Estimates versus the known truth

| quantity | estimate (95% CI) | true value | inside CI? |
|---|---|---|---|
| deer density, whole area (per km2) | 12.0 (9.24-15.5) | 11.31 | yes |
| abundance | 4,314 (3,328-5,593) | 4,071 | yes |
| density north / south | 7.97 (5.68-11.2) / 16.0 (11.6-22.0) | 8.43 / 14.19 | yes / yes |
| groups per km2 | 3.59 (2.82-4.58) | 3.33 | yes |
| mean group size | 3.34 (3.12-3.56) | 3.4 | yes |
| detection probability within 400 m | 0.437 (0.378-0.499) | 0.456 | yes |
| effective strip half-width (m) | 175 (151-200) | 182 | yes |
| half-normal sigma (m) | 138 (122-156) | 135 (JK) / 165 (AB); 146 pooled | yes (pooled) |
| shrub effect per 10 % points | 1.89 (1.56-2.30) | 1.65 | yes |
| road effect per 3 km | 1.40 (1.12-1.74) | 1.28 | yes |
| density at an average transect (HDS) | 10.3 (8.12-13.1) | 10.37 | yes |
| HDS density over the grid | 12.2 (9.70-15.3) | 11.31 | yes |

Every interval contains the truth. The two habitat effects are estimated about 10-15 %
too strong, well inside their intervals. The two observers had different
detection scales (135 and 165 m). The pooled detection function still gave an
accurate detection probability (0.437 vs 0.456), as "pooling robustness" predicts.

## (f) Limits from VALIDATION.md that matter at this design

1. **Known limit 9 (HDS intercept):** "density at an average transect" uses only the
   AIC-best key. Coverage here was 0.885 (0.92 when the half-normal is chosen, 0.81 otherwise).
   The design-based density and the habitat slopes are unaffected.
2. **Design-based intervals are conservative** when density follows habitat across
   transects (coverage 0.975 here, limit 3). The density interval is wider than it needs to be.
3. **g(0) = 1, no evasive movement, accurate distances** cannot be tested from these data.
   75 m of evasive movement biased density by -27 % in simulation, and heaping at zero biased it
   upward (limit 4). Deer in winter may flush before being seen, so field protocol matters.
4. **About 145 detections** is above the ~60 needed for a stable detection function
   (limit 1 applies only below ~40).
5. Recovery simulated one habitat covariate and one detection scale. The example has two
   covariates and two observers. The road effect and the observer difference were therefore
   checked only on the example itself (where both are recovered), not by 200 replicates.
6. The map extrapolates to shrub cover of 3-14 % and 45-55 %, outside the surveyed range.

## (g) Commands to reproduce

```bash
cd .
# true values (re-runs the generator, checks it reproduces the CSVs)
scripts/Rscript -e 'source("core/R/env.R"); source("dev/distance/si/example_truth.R")'
# planted truth at this design (200 simulated surveys, ~12 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/distance_si_recover.log -- \
  scripts/Rscript core/run.R recover distance --config modules/distance/example/recovery.yaml \
  --nsim 200 --cores 16 --out results/distance/si/recovery
# fit the example data: report + figures (~10 min)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/distance_si_fit.log -- \
  scripts/Rscript core/run.R fit distance --config modules/distance/example/config.yaml \
  --out results/distance/si/fit
# redraw figures only, after editing the figures: block
scripts/ejob -c 1 -m 8 -t 0:20:00 -- scripts/Rscript core/run.R figures distance \
  --fit results/distance/si/fit/distance_fit.rds --config modules/distance/example/config.yaml
```
