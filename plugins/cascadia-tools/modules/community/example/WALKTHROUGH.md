# Worked example: forest bird communities along elevation and stand type

Module: `community` (ordination, RDA and variation partitioning, PERMANOVA with
a dispersion test, beta-diversity partitioning, indicator species).
Appendix S1 example and website "try an example" walkthrough.

**The data are simulated** (see `README_data.md`), so we know the right
answers and can check the analysis against them.

## (a) The question, in the researcher's words

> "I counted breeding birds at 60 stations in old growth, thinned stands and
> young plantations in the Oregon Cascades. How much do elevation and canopy
> cover shape which birds live where? Do the three stand types have different
> bird communities? Is it species swapping out or just poorer sites missing
> species? And which birds go with which stand type?"

## (b) The files they have

`cascades_bird_counts.csv`: one row per station, one column per species
(first 10 of 30 species shown):

```
"station","PSFL","PAWR","HEWA","CBCH","BRCR","GCKI","VATH","SWTH","WIWA","MGWA"
"CAS-01",0,3,1,2,0,17,2,0,0,1
"CAS-02",0,3,2,0,1,0,1,0,0,6
"CAS-03",0,0,0,5,2,2,3,0,0,0
"CAS-04",1,8,0,0,0,0,0,1,1,6
```

`cascades_stations.csv`: one row per station:

```
"station","stand_type","elevation_m","canopy_cover_pct","survey_date"
"CAS-01","young_plantation",753,14,"2025-05-24"
"CAS-02","thinned",833,66,"2025-05-26"
"CAS-03","thinned",724,39,"2025-05-22"
"CAS-04","old_growth",1212,76,"2025-05-25"
```

## (c) What they type to Claude

> I have bird point counts from 60 forest stations (cascades_bird_counts.csv,
> one column per species) and a station table with stand type, elevation and
> canopy cover (cascades_stations.csv). How much of the variation in the bird
> community goes with elevation versus canopy cover, do old growth, thinned
> and young plantation stands hold different communities, is it turnover or
> nestedness, and which species are indicators of each stand type? Please
> give me a report and figures I can use in a paper.

## (d) What the agent does, step by step (skills/eco-analyze)

**1. Project folder.** The agent makes `projects/<name>/` and copies the raw
files into `data_raw/` untouched. (In this appendix the files stay in
`modules/community/example/`, and outputs go to `results/community/si/`.)

**2. Data intake (eco-data-intake).** It reads both files back to the
researcher:

- 60 stations, 30 species, 1,774 birds in total.
- 3 to 17 species per station (median 10).
- The most widespread species is Pacific Wren (47 stations); the rarest is
  Hutton's Vireo (5 stations).
- Stand types: old growth 22, thinned 20, young plantation 18.
- Elevation 407-1582 m; canopy cover 14-88%. The two are correlated
  (r = 0.33), so part of their effect cannot be told apart.
- Station ids match in both files. No empty stations. No missing values.
- One count per station, the same effort everywhere. `survey_date` is not
  used.

**3. Choose the module.** `community` fits: a station-by-species table, two
gradients and a grouping factor, and questions about whole-community
composition. The agent notes what it does not do: it does not correct for
birds that were present but missed (that needs repeat visits and an
occupancy model such as `occu-spatial` or `jsdm`).

**4. Analysis plan (eco-analysis-plan).** Written before fitting:

- Estimands: the adjusted R2 for elevation, for canopy cover, their unique
  and shared parts, and the total (RDA on Hellinger-transformed counts).
  Also the PERMANOVA adjusted R2 for stand type (Bray-Curtis), always read
  with the PERMDISP dispersion test. Also mean pairwise Sorensen
  dissimilarity split into turnover and nestedness, and IndVal.g indicator
  species.
- Intervals: delete-one-station jackknife 95% CIs. Tests: 999 permutations.
- Checks: NMDS stress, DCA gradient length, agreement of CCA and db-RDA with
  RDA, and PERMDISP.

**5. Planted truth at this design (eco-planted-truth).** The agent simulates
200 datasets with the same design: 60 stations; groups of 22/20/18; two
gradients correlated 0.4; stand-type effects SD 0.9; young plantations 1.4
times as variable. Each dataset has a new species pool. It then checks that
the analysis recovers the known values
(`modules/community/example/recovery.yaml`, 16 cores, about 2 minutes; about 6 s per dataset).

| Estimand | mean truth | mean estimate | rel. bias | 95% CI coverage | gate |
|---|---|---|---|---|---|
| R2adj elevation | 0.293 | 0.291 | -0.8% | 0.955 | PASS |
| R2adj canopy cover | 0.110 | 0.109 | -0.3% | 0.935 | PASS |
| unique elevation | 0.247 | 0.245 | -0.6% | 0.940 | PASS |
| unique canopy cover | 0.064 | 0.064 | -0.4% | 0.955 | PASS |
| **shared** | 0.046 | 0.046 | +2.8% | **0.890** | **FAIL (coverage)** |
| total (both gradients) | 0.356 | 0.354 | -0.6% | 0.935 | PASS |
| stand type, PERMANOVA R2adj | 0.097 | 0.100 | +3.0% | 0.970 | PASS |
| pairwise Sorensen | 0.524 | 0.525 | +0.3% | 0.965 | PASS |
| pairwise turnover | 0.442 | 0.444 | +0.6% | 0.965 | PASS |
| pairwise nestedness | 0.082 | 0.081 | -1.0% | 0.925 | PASS |
| turnover share | 0.843 | 0.844 | +0.2% | 0.955 | PASS |

All 200 datasets ran. The coverage threshold at 200 replicates is 0.904.

In one sentence to the researcher: *"In 200 simulated surveys like yours,
10 of 11 quantities were recovered with less than 3% bias and correct
intervals. The shared elevation-canopy part was unbiased, but its interval
was too narrow (it contained the truth 89% of the time instead of 95%), so
we report that number without trusting its interval."* This matches
VALIDATION.md known limits 1-2: the shared fraction's jackknife interval runs
narrow when one group is more variable.

**6. Fit.** `core/run.R fit community --config modules/community/example/config.yaml`
(4 cores, under 1 minute including 999-permutation tests, NMDS, DCA, CCA,
db-RDA, indicator species and response curves).

**7. Checks (eco-convergence, eco-model-checking).** There is no MCMC. These
are least-squares and permutation methods, so "converged" is always true.
The checks are:

- **NMDS stress 0.163.** Usable (below 0.2). Do not over-read small
  distances. NMDS did not find the same best solution twice in 50 random
  starts, so the picture is approximate. Inference does not depend on it.
- **DCA axis 1 = 4.6 SD.** This is a long gradient: most species are
  replaced along it. Raw-count RDA would be wrong. Hellinger-transformed RDA
  is designed for this case.
- **Other methods agree.** Total explained: RDA 29.3%, db-RDA on Bray-Curtis
  30.0%, CCA 23.9%. CCA uses a different distance that weights rare species
  more, so some difference is expected.
- **Dispersion test (PERMDISP) p = 0.78.** Mean distance to the group centre
  is 0.525 for old growth, 0.527 for thinned and 0.509 for young plantations
  (largest/smallest 1.04). So the PERMANOVA result is not driven by a
  difference in spread. Note that young plantations *were* simulated 1.4
  times as variable. PERMDISP missed this modest difference, which is
  consistent with its limited power.

**8. Interpretation (eco-interpretation).** Estimates with 95% intervals:

| Quantity | Estimate (95% CI) | p (999 perms) |
|---|---|---|
| Elevation alone | 24.6% (20.0 to 29.1) | 0.001 |
| Canopy cover alone | 7.7% (2.7 to 12.6) | 0.001 |
| Elevation after canopy | 21.7% (15.7 to 27.6) | 0.001 |
| Canopy after elevation | 4.8% (2.4 to 7.2) | 0.001 |
| Shared | 2.9% (-1.9 to 7.8), descriptive | no test |
| Both together | 29.3% (24.6 to 34.0) | 0.001 |
| Stand type (PERMANOVA) | 3.9% (0.9 to 6.9) | 0.011 |
| Pairwise Sorensen dissimilarity | 0.63 (0.61 to 0.66) | |
| Turnover / nestedness | 0.54 (0.50 to 0.58) / 0.10 (0.08 to 0.12) | |
| Turnover share | 85% (81 to 88) | |

Excerpts from the report's "What this means biologically" section (each was
checked against the table):

> "Together the variables accounted for 29.3% (95% CI 24.6% to 34.0%) of the
> among-site variation (adjusted R2, p = 0.001), leaving 70.7% unexplained."

> "These percentages are shares of variation in community composition, not
> effect sizes for any single species ... The association is correlational:
> unmeasured variables that co-vary with elevation (m) could produce the same
> pattern."

> "Stand type groups ... differed in composition: group membership accounted
> for 3.9% (95% CI 0.9% to 6.9%) of the variation in Bray-Curtis
> dissimilarities (PERMANOVA adjusted R2, p = 0.011). ... Because PERMDISP
> found no dispersion difference, the PERMANOVA result is most consistent with
> a shift in average species composition (location). PERMDISP has limited
> power for modest differences in spread, though, so a smaller contribution
> from dispersion cannot be ruled out."

> "Pairwise PERMANOVA (Holm-adjusted): old growth vs thinned R2 = 0.038,
> p = 0.124; old growth vs young plantation R2 = 0.059, p = 0.040; thinned vs
> young plantation R2 = 0.070, p = 0.027."

> "Accounting for the environmental variables, stand type still explained
> 6.9% of the dissimilarity (marginal PERMANOVA, unadjusted R2, p = 0.001)."

> "1 species was an indicator of one stand type after false-discovery-rate
> correction ...: Golden-crowned Kinglet for young plantation (IndVal 0.78,
> A 0.83, B 0.72, p = 0.030)."

> "Of this, 85% (95% CI 81% to 88%) was turnover (species replaced by other
> species) and the rest nestedness ... Composition changes mainly by
> replacement."

The agent adds one project-specific caveat. The shared fraction's interval
is not trusted at this design (step 5). The report's text says it "has no
test", which is correct; the agent also tells the researcher not to quote its
interval.

Figures (`results/community/si/fit/figures/`, PDF and 600-dpi PNG; every PNG
was opened and checked):

| Figure | Caption (short) |
|---|---|
| `ordination` | NMDS of Bray-Curtis dissimilarities (stress 0.163), stations coloured by stand type, with fitted elevation (r2 = 0.94) and canopy (r2 = 0.56) arrows and 8 species. |
| `nmds_shepard` | Shepard diagram: ordination distance against observed dissimilarity. |
| `rda_triplot` | RDA triplot: stations, gradient arrows and the 8 best-represented species (constrained axes: adjusted R2 29.3%). |
| `variance_partition` | **Key figure.** Unique, shared, total and unexplained shares of variation with 95% jackknife CIs. |
| `effect_sizes` | Adjusted R2 with 95% CIs for every gradient fraction and for stand type (PERMANOVA). |
| `beta_components` | Pairwise Sorensen split into turnover and nestedness over all 1,770 station pairs; turnover share 85%. |
| `distance_decay` | Turnover rises from about 0.4 to about 1 as the elevation difference grows to 1,200 m; nestedness stays low. |
| `dispersion` | Distance of each station to its stand-type centroid (PERMDISP p = 0.776). |
| `indicator_species` | IndVal.g: Golden-crowned Kinglet (young plantation) is significant after FDR correction; 7 others only at uncorrected p < 0.05. |
| `species_response_elevation_m` | Expected count against elevation for the 6 species most tied to it (e.g. Swainson's Thrush peaking near 1,500 m). |
| `species_response_canopy_cover_pct` | Expected count against canopy cover for the 6 species most tied to it. |

One presentation fix was made in `config.yaml` (`figures: per_figure`).
Hutton's Vireo was seen at only 5 stations, and its confidence band ran to
about 300 birds at low canopy cover. That stretched the whole panel. The
y-axis is now capped at 15 birds (the band is clipped, not removed).

**9. Acceptance gate (eco-acceptance-gate): OPEN.**

- Planted truth: 10 of 11 estimands pass. The shared fraction fails on
  coverage (0.890 vs 0.904).
- Fit checks: pass (stress usable, methods agree, dispersion test clean).
- Independent route: the module's estimates equal vegan::varpart, adonis2
  and betapart to 3e-16 (VALIDATION.md section 5).
- Negative controls: VALIDATION.md section 4.
- Reproducible: commands below.

Everything except the shared fraction's interval can be reported. The
verdict stays OPEN because the report prints an interval for an estimand
that failed the gate at this design. The researcher should give the shared
fraction as a descriptive value (2.9%) without its CI.

**10. Delivery: summary for the researcher.**

> Elevation was the main correlate of which birds occurred where: alone it
> accounted for 24.6% (95% CI 20.0-29.1%) of the variation in community
> composition, and 21.7% (15.7-27.6%) after allowing for canopy cover. Canopy
> cover added a smaller, real part: 4.8% (2.4-7.2%) beyond elevation. Stand
> types differed modestly in composition (3.9%, 0.9-6.9%; p = 0.011), mostly
> young plantations versus the other two. The dispersion test found no
> difference in variability (p = 0.78), so this looks like a shift in which
> species are present. Stations differ mostly by species replacing each
> other (85% of dissimilarity is turnover). Golden-crowned Kinglet was the
> one clear indicator (of young plantations). These are associations, not
> causes, and they describe birds as counted on one visit.

## (e) Estimates versus the known truth

The truths are for this dataset's own species pool (`README_data.md`).

| Estimand | Truth | Estimate (95% CI) | Inside CI? |
|---|---|---|---|
| R2adj elevation | 0.252 | 0.246 (0.200 to 0.291) | yes |
| R2adj canopy cover | 0.106 | 0.077 (0.027 to 0.126) | yes |
| unique elevation | 0.212 | 0.217 (0.157 to 0.276) | yes |
| unique canopy cover | 0.067 | 0.048 (0.024 to 0.072) | yes |
| shared | 0.040 | 0.029 (-0.019 to 0.078) | yes |
| total | 0.319 | 0.293 (0.246 to 0.340) | yes |
| stand type PERMANOVA R2adj | 0.050 | 0.039 (0.009 to 0.069) | yes |
| pairwise Sorensen | 0.619 | 0.634 (0.605 to 0.663) | yes |
| pairwise turnover | 0.517 | 0.537 (0.499 to 0.576) | yes |
| pairwise nestedness | 0.102 | 0.097 (0.075 to 0.118) | yes |
| turnover share | 0.835 | 0.847 (0.812 to 0.883) | yes |

All 11 intervals contain the truth. The young plantations' extra variability
(1.4 times) was real but was not detected (PERMDISP p = 0.78). That is a
correct "no detectable difference", not proof of equal spread.

## (f) Limits from VALIDATION.md that matter here

- **Shared fraction.** Its jackknife interval is too narrow when groups
  differ in variability or with fewer than about 30 sites. Here coverage was
  0.890. Report it without an interval. It also has no test.
- **Jackknife CIs are not tests.** Use the permutation p-values to decide
  whether an effect exists. The PERMANOVA CI is conservative.
- **PERMANOVA mixes location and dispersion** in unbalanced designs. Groups
  here are nearly balanced (22/20/18) and PERMDISP was clean. PERMDISP has
  limited power, though, and is slightly liberal (7% false positives).
- **No detection model.** Results describe birds as counted on one 10-minute
  visit. Missed birds inflate apparent turnover.
- **Sites are treated as independent.** Clustered stations (several per
  stand) would need restricted permutations, which are not built in.
- **Correlation, not causation.** Elevation stands for everything that
  changes with it (temperature, snow, forest type).
- **Descriptive outputs** (NMDS, DCA, CCA, db-RDA, multi-site beta,
  indicator values) are checked by tests and inspection, not by recovery.

## (g) Commands to reproduce

From `.`:

```bash
# 0. (optional) regenerate the example data and truth (seed 2026), then check them
scripts/ejob -c 1 -m 4 -- scripts/Rscript dev/community/make_example.R
scripts/ejob -c 1 -m 4 -l logs/community_si_check.log -- scripts/Rscript dev/community/si/check_example.R

# 1. planted-truth recovery at the example design (200 datasets, ~2 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/community_si_recover.log -- \
  scripts/Rscript core/run.R recover community --config modules/community/example/recovery.yaml \
  --nsim 200 --cores 16 --out results/community/si/recovery

# 2. fit the example data: report, estimates and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/community_si_fit.log -- \
  scripts/Rscript core/run.R fit community --config modules/community/example/config.yaml \
  --out results/community/si/fit

# 3. (optional) redraw figures only, after editing the figures: block
scripts/Rscript core/run.R figures community --config modules/community/example/config.yaml \
  --fit results/community/si/fit/community_fit.rds
```

Outputs: `results/community/si/fit/community_report.html`,
`results/community/si/fit/community_estimates.csv`,
`results/community/si/fit/figures/`,
`results/community/si/recovery/community_recover_summary.csv`,
`results/community/si/si_summary.json`.
