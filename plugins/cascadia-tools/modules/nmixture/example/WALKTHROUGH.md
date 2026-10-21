# Worked example: how many Hermit warblers, and where? (N-mixture)

This example shows the full CascadiaTools workflow for repeated counts. The data are
**simulated**, so we know the true answer and can check every number against it.
The data and the true values are described in `README_data.md`.

## (a) The question, in the researcher's words

> "We counted Hermit warblers at 120 points in the Coast Range, up to three times each this
> spring. How many warblers are actually out there, given that we surely missed some? And are
> there more warblers where the canopy is denser, or higher up?"

## (b) The files the researcher has

`point_counts.csv`: one row per visit (346 rows).

```
"point","visit","date","julian_day","start_time","wind_beaufort","observer","count"
"PT-001",1,"2025-05-13",133,"06:15",1,"SAT",1
"PT-001",2,"2025-05-26",146,"07:00",1,"SAT",0
"PT-001",3,"2025-06-11",162,"05:30",2,"MLR",0
"PT-002",1,"2025-05-10",130,"05:45",3,"SAT",2
```

`points.csv`: one row per point (120 rows).

```
"point","canopy_cover_pct","elevation_m","stand_type"
"PT-001",50,251,"old"
"PT-002",26,722,"old"
"PT-003",58,715,"young"
"PT-004",51,232,"mature"
```

## (c) What the researcher types to Claude

> I have Hermit warbler point counts from 120 points, each visited up to 3 times between
> mid-May and mid-June. The counts are in point_counts.csv (one row per visit, with wind,
> day of year and observer). The habitat for each point is in points.csv (canopy cover and
> elevation). Please estimate how many warblers there are, allowing for birds we missed,
> and tell me whether abundance is related to canopy cover and elevation.

## (d) What the agent does, step by step (`skills/eco-analyze`)

**1. Project folder.** Raw files are copied, never edited. All work goes in one folder.

**2. Data intake.** The agent reads both files and describes them back:
- 120 points, 3 planned visits each; 346 of 360 visits were done (14 third visits missing, 3.9%).
- Counts range from 0 to 6 per visit. The sum over points of the largest count is 217.
- Visit covariates: wind 0-4 Beaufort (mean 1.4), day of year 130-164, three observers.
- Point covariates: canopy cover 16-98 %, elevation 152-1098 m. `stand_type` and `start_time`
  are recorded but not in the question, so they are left out.
- No duplicated visits; every point in the counts file appears in the points file.

**3. Choose the module.** Repeated counts of unmarked animals at fixed sites, within one
season: the `nmixture` module, static (closed) model.

**4. Plan.** Written before any fit:
- Estimands: total warblers at the 120 points; expected abundance per point; effects of
  canopy cover and elevation on abundance (rate ratios per 1 SD); detection per bird per visit.
- Model: abundance ~ canopy + elevation (Poisson, log link); detection ~ wind + day + observer
  (logit link). Poisson is compared with negative binomial and zero-inflated Poisson by AIC.
- Engine: maximum likelihood (`unmarked::pcount`).
- Checks: parametric-bootstrap goodness of fit (300 replicates); sensitivity to the bound K.

**5. Planted truth at this design.** The agent simulates 200 data sets with the same design
(120 points, 3 visits, 3.9% visits missing, Poisson abundance) and realistic values (the
generating values of this example). It fits each one and checks the answers
(`recovery.yaml`). Result: **all 8 estimands PASS**, 0 of 200 fits failed.

| quantity | truth | mean estimate | relative bias | 95% CI coverage | power |
|---|---|---|---|---|---|
| expected abundance, average point | 2.10 | 2.12 | +0.8% | 0.950 | - |
| total abundance (120 points) | 277 | 282 | +1.9% | 0.985 | - |
| mean abundance per point | 2.31 | 2.35 | +1.9% | 0.985 | - |
| detection per visit | 0.553 | 0.549 | -0.7% | 0.960 | - |
| canopy effect (log, per SD) | 0.434 | 0.433 | -0.3% | 0.965 | 1.00 |
| canopy rate ratio per SD | 1.54 | 1.54 | +0.1% | 0.965 | 1.00 |
| wind effect on detection (logit, per SD) | -0.496 | -0.501 | +1.0% | 0.970 | 1.00 |
| wind odds ratio per SD | 0.609 | 0.608 | -0.1% | 0.970 | - |

The pass rule is |bias| <= 10% of truth (or 3 Monte Carlo SE) and coverage >= 0.904 at 200
replicates. (Total abundance truth here is the average over replicates, 277, not the 273 of
the one example population.) The simulator has one point covariate and one visit covariate,
so elevation, day of year and observer were not in this check; they are checked against the
known truth in step (e) below.

**6. Fit.** `core/run.R fit` with `config.yaml`. Takes under a minute.

**7. Checks.**
- Convergence: the optimiser converged; all standard errors finite; no estimate on a boundary.
  (Maximum likelihood, so there is no R-hat or ESS.)
- Abundance distribution by AIC: Poisson 836.5, NB 838.5, ZIP 838.5. Poisson is best;
  all three give 2.00 birds at an average point.
- Goodness of fit (300 bootstrap replicates): chi-square 397.3, c-hat 1.17, p = 0.060;
  Freeman-Tukey p = 0.143. No lack of fit at the 5% level, but borderline.
- K sensitivity: total abundance changes by 0.00% when K goes from 36 to 54 to 72.

**8. Interpretation.** Excerpt from the report (`results/nmixture/si/fit/nmixture_report.html`):

> "Summed over the 120 surveyed sites, estimated abundance is 272 (95% CI 244 to 308) (an
> average of 2.26 (95% CI 2.03 to 2.57) per site). The sum of the largest count at each site
> was 217, so the model infers that about 20% of the animals present were missed even on the
> visit with the highest count."
>
> "Each increase of 1 SD (19.74 units) in Canopy cover is associated with expected abundance
> being multiplied by 1.59 (95% CI 1.37 to 1.84) ... Across the observed range (16.00 to 98.00)
> expected abundance per site goes from 0.71 to 4.90."
>
> "Each increase of 1 SD (293.79 units) in Elevation is associated with expected abundance being
> multiplied by 1.19 (95% CI 1.05 to 1.36)."
>
> "Detection: on a visit with average covariate values (Observer 'JDK') each animal present had a
> 0.51 (95% CI 0.41 to 0.61) probability of being counted; over 3 such visits the probability of
> counting an animal at least once is 0.89 (0.80 to 0.94)."
>
> "Relative differences among sites and covariate effects (rate ratios) are more robust than the
> absolute totals. All covariate relationships are associations, not demonstrated causes."

Figures (`results/nmixture/si/fit/figures/`, PDF and PNG):

| figure | caption |
|---|---|
| `abundance_canopy_cover_pct` | Expected warblers per point vs canopy cover (95% and 50% bands; other covariates at their mean). |
| `abundance_elevation_m` | Expected warblers per point vs elevation. |
| `abundance_effects` | Both abundance effects as panels A and B. |
| `detection_wind_beaufort` | Detection per visit vs wind. A counting effect, not a habitat effect. |
| `detection_julian_day` | Detection per visit vs day of year. |
| `detection_observer` | Detection per visit by observer. |
| `detection_effects` | All three detection effects as panels A-C. |
| `coefficients` | Rate ratios and odds ratios per 1 SD (log scale; dashed line = no effect). |
| `abundance_by_site` | Estimated warblers at each point (95% interval) and the largest count there (crosses). |
| `total_abundance` | Uncertainty in the total over 120 points; dot-dash line = sum of largest counts (217). |
| `gof` | Bootstrap goodness of fit: observed chi-square vs simulated (p = 0.06). |
| `K_sensitivity` | Total abundance at K = 36, 54, 72 (flat = K is large enough). |

The single most informative figure is `abundance_canopy_cover_pct`: it answers the main
habitat question on the scale the researcher thinks in (birds per point).

**9. Gate verdict: DONE.** Planted truth passed at this design, the fit converged, the fit
checks passed (borderline GOF stated), and every number comes from a rerunnable command.

## (e) Estimates versus the known truth

| quantity | estimate (95% CI) | truth | inside CI? |
|---|---|---|---|
| total warblers at 120 points | 272 (244 to 308) | 273 | yes |
| mean per point | 2.26 (2.03 to 2.57) | 2.275 | yes |
| expected per average point | 2.00 (1.69 to 2.38) | 2.10 | yes |
| canopy rate ratio per SD | 1.59 (1.37 to 1.84) | 1.54 | yes |
| elevation rate ratio per SD | 1.19 (1.05 to 1.36) | 1.19 | yes |
| detection per visit (JDK, average conditions) | 0.51 (0.41 to 0.61) | 0.55 | yes |
| wind odds ratio per SD | 0.49 (0.40 to 0.60) | 0.61 | **no** |
| day-of-year odds ratio per SD | 0.83 (0.71 to 0.97) | 0.79 | yes |
| observer MLR vs JDK odds ratio | 1.57 (1.03 to 2.40) | 1.00 | **no** |
| observer SAT vs JDK odds ratio | 0.75 (0.49 to 1.15) | 0.67 | yes |

All abundance quantities, which answer the question, contain the truth. Two of the detection
intervals miss it. Each is about 2 SE off. With 95% intervals, about 1 in 20 misses by chance,
and the detection terms share the same counts, so their errors move together. The recovery
study (step 5) shows the wind effect has 97% coverage over 200 data sets, so this is chance in
this one data set, not a flaw. It is a useful reminder: one data set can mislead on a single
detection covariate, even when the abundance answers are right.

## (f) Limits that matter at this design (from `VALIDATION.md`)

- Absolute N relies on closure, no double counting and independent detection. Double counting
  and closure violations were not tested, and the GOF test rarely detects them. Treat the total
  (272) as conditional on these assumptions; the canopy and elevation effects are more robust.
- Unmodelled detection heterogeneity biases N downward with intervals that are too narrow.
  Strong site-level heterogeneity was flagged by the GOF test 95% of the time; individual-level
  heterogeneity was not tested. Here GOF p = 0.06 is borderline, so it is worth stating.
- Overdispersed abundance fitted as Poisson biases N by about -20%. Here AIC favoured Poisson and
  NB/ZIP gave the same answer, so this is unlikely to matter.
- With low detection (p about 0.2) and few visits, absolute N can be badly biased. Here p is about
  0.5 with 3 visits, a design where recovery passed.
- The recovery study used fast settings (no GOF, K-sensitivity or mixture comparison) and one
  covariate per submodel. Elevation, day and observer effects were checked only against the
  known truth of this one example.

## (g) Commands to reproduce

Run from `.`. Use `scripts/ejob` (it submits a SLURM job).

```
# example data (seed 20250512) and the true values on the report scale
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/nmixture_make_example.log -- scripts/Rscript dev/nmixture/make_example.R
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/nmixture_si_truth.log -- scripts/Rscript dev/nmixture/si_truth.R

# planted truth at this design (200 replicates, about 4 minutes on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/nmixture_si_recover.log -- scripts/Rscript core/run.R recover nmixture \
  --config modules/nmixture/example/recovery.yaml --nsim 200 --cores 16 --out results/nmixture/si/recovery

# fit, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/nmixture_si_fit.log -- scripts/Rscript core/run.R fit nmixture \
  --config modules/nmixture/example/config.yaml --out results/nmixture/si/fit
```
