# Worked example: Pacific marten numbers, survival and habitat from five years of hair snares (module `scr-open`)

Appendix S1 worked example and website "try an example" walkthrough.
**The data are simulated with known truth** (see `README_data.md`), so every
estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "We run 49 hair snares for Pacific marten every summer in the Oregon Cascades,
> and the lab gives us an ID for every hair sample. We have five summers now,
> 2020 to 2024. How many martens are out there each year? Is the population
> going down? How many survive from one year to the next, and how many new ones
> show up? And do martens really live at higher density under thick canopy?"

## (b) The files the researcher has

`hair_snares.csv`: one row per snare (49 rows). `ops_<year>` says which of the 5 weekly checks the snare was working.
```
"snare_id","utm_e","utm_n","elevation_m","ops_2020","ops_2021","ops_2022","ops_2023","ops_2024"
"HS01",574000,4864000,1351,"11111","11111","11111","11111","11111"
"HS02",575000,4864000,1212,"11111","11111","11111","11111","11111"
"HS03",576000,4864000,1232,"11111","11111","11111","11111","11111"
"HS04",577000,4864000,1244,"11111","11111","11111","11111","11111"
```

`marten_genotypes.csv`: one row per genotyped hair sample (598 rows).
```
"sample_id","year","check","snare_id","marten_id","n_hairs"
"S42482","2020",1,"HS01","MACA-476",8
"S75084","2020",1,"HS01","MACA-463",3
"S77658","2020",1,"HS01","MACA-724",16
"S56664","2020",1,"HS02","MACA-669",24
```

`habitat_grid.csv`: one row per 400 m habitat pixel, covering the snares plus 3 km (900 rows).
```
"x","y","canopy_cover_pct"
571200,4861200,49.2
571600,4861200,49.6
572000,4861200,50.3
572400,4861200,51.5
```

## (c) What the researcher types to Claude

> I have five summers (2020-2024) of marten hair-snare data from 49 snares on a
> 1-km grid, each checked 5 times a summer. hair_snares.csv has the snare
> locations (UTM) and which checks each snare was working; marten_genotypes.csv
> has one row per genotyped sample with year, check, snare and marten ID; and
> habitat_grid.csv has canopy cover on a 400 m grid around the snares. I want to
> know how many martens there are each year, whether they are declining, annual
> survival and recruitment, and whether density is higher under denser canopy.
> Please give me a report with figures I can use in a paper.

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

**1. Project folder.** The three files are copied to `projects/<name>/data_raw/`.
The raw files are never edited. Every command is logged.

**2. Data intake.** The agent reads the files and describes them back:
- 49 snares, 1,000 m apart, 7 x 7 grid. Coordinates are UTM metres.
- 5 summers x 5 weekly checks. A few snares failed part-way: 2, 0, 2, 3 and 5 snares in 2020-2024.
- 95 different martens; 32, 38, 39, 41 and 48 detected per year.
- 48 martens were detected in two or more years. These recaptures carry the information on survival.
- 598 detections (one per marten, snare and check).
- Canopy cover in the grid runs from 24% to 76% (mean 50%).
`check_data()` found no problems (all snare IDs match, all checks are 1-5, no
detections at closed snares).

**3. Choose the module.** Individually identified animals at fixed detectors,
over several years, with questions about yearly numbers, survival and
recruitment: this is open-population spatial capture-recapture (`scr-open`).
Closed SCR (`scr`) would treat each year separately and give no survival.
Non-spatial Jolly-Seber (`jolly-seber`) would ignore where martens were caught,
so it could not give density or the canopy effect.

**4. Plan.** Estimands: martens alive in the 144 km2 study area each year,
density, superpopulation, recruits and per-capita recruitment, growth rate,
apparent survival, and the canopy effect on density. Model: Jolly-Seber SCR
with constant survival and detection, static home-range centres, half-normal
detection, density proportional to exp(beta x canopy). Engine: Nimble,
marginal likelihood (3 chains x 30,000 iterations, 10,000 burn-in). Cross-checks:
NumPyro (a second Bayesian engine) and openCR (maximum likelihood).

**5. Planted truth at this design** (`recovery.yaml`, 100 simulated datasets,
each fitted exactly as above but with the module's shorter settings).
The design and values match the example: 7 x 7 snares 1,000 m apart, 5 years x 5 checks,
3 km buffer, density 0.8 per km2, p0 = 0.30, sigma = 800 m, survival 0.80,
per-capita recruitment 0.25, canopy effect 1.0.
**Gate result: 28 of 28 quantities pass.** All 100 fits converged. Median fit
time was 7 min (max 11 min) on one core.

| quantity | truth (mean over datasets) | mean estimate | relative bias | 95% coverage |
|---|---|---|---|---|
| abundance 2020 | 114.4 | 115.3 | +0.6% | 0.98 |
| abundance 2024 | 139.0 | 140.9 | +1.6% | 0.98 |
| superpopulation | 236.5 | 239.5 | +1.4% | 0.93 |
| recruits 2021-2024 | 27-33 | | -0.7% to +3.7% | 0.94-0.98 |
| per-capita recruitment 2021-2024 | 0.24-0.26 | | -0.5% to +5.3% | 0.93-0.99 |
| growth rate 2021-2024 | 1.04-1.06 | | +0.1% to +1.3% | 0.94-0.98 |
| survival | 0.80 | 0.798 | -0.3% | 0.98 |
| baseline detection p0 | 0.30 | 0.301 | +0.2% | 0.91 |
| sigma (m) | 800 | 800.5 | +0.1% | 0.96 |
| canopy effect | 1.00 | 1.01 | +1.0% | 0.98 (power 1.00) |

In one sentence: at this design the model recovers every quantity to within
about 5% and its 95% intervals contain the truth 91-99% of the time.

**6. Fit.** `config.yaml` names the files and columns. The fit ran on 4 cores:
Nimble took 8.7 min; with the NumPyro (CPU) and openCR cross-checks the whole
job took 74 min.

**7. Check.**
- Convergence: max R-hat 1.005, min bulk ESS 1,035 (the standard is R-hat < 1.01 and ESS >= 400). The traces of all three chains overlap.
- Goodness of fit: posterior predictive p-values were all between 0.39 and 0.64:
  - individuals detected per year 0.47 / 0.49 / 0.55 / 0.64 / 0.52, and in total 0.53;
  - detected in 2+ years 0.61; missed years within histories 0.39;
  - total detections 0.49; snares per marten-year 0.46.
  Values near 0 or 1 would flag a problem. None did.
- Design notes: snare spacing is 1.2 sigma, and the buffer is 3.7 sigma. Both are adequate.
- **Warning.** The snares sit in below-average canopy (0.45 SD below the area mean).
  Under the fitted model, 72% of martens live more than 2 sigma from any snare.
  Average density over the whole area is 1.94 x the density around the snares.
  So total abundance depends on extrapolating the canopy effect into dense
  forest that the snares did not sample.
- Second engine: NumPyro gives the same answers. Abundance in 2020 is 121 (87-166) vs Nimble's 119 (85-161).
  Superpopulation 279 vs 277. Survival 0.839 in both. Canopy effect 1.32 in both.
  openCR (maximum likelihood, *homogeneous* density) gives very similar survival (0.852),
  sigma (788 m) and per-capita recruitment (0.34 in 2021).
  Its abundance is lower (72 in 2020, 108 in 2024), as the warning predicts.

**8. Interpret.** Excerpt from the report (`results/scr-open/si/fit/scr-open_report.html`):

> Abundance in the 144 km² state space (traps plus buffer) was 119 (95% CrI 85.0 to 161) in period 2020 and
> 186 (95% CrI 144 to 238) in period 2024 (density 0.83 and 1.3 animals per km²). Over the study the
> population changed by a factor of 1.58 (95% CrI 1.15 to 2.16); posterior probability of a decline from the
> first to the last period: < 0.01.
>
> Apparent survival between consecutive periods was 0.84 (95% CrI 0.77 to 0.91): an animal alive in one period
> had an estimated 84% chance of still being alive and resident in the next. 'Apparent' because death and
> permanent emigration from the study area cannot be told apart.
>
> In total 277 (95% CrI 222 to 341) different animals were alive in the state space at some point during the
> study (superpopulation), of whom 95 (34%) were ever detected.
>
> Density and Canopy cover (%): each 1 SD (12) higher Canopy cover (%) is associated with 273% higher density
> (rate ratio 3.7, 95% CrI 2.9 to 5.0). [...] This describes where animals' activity centres are, not proof that
> Canopy cover (%) causes higher density.

The agent checked each sentence against the estimates table. For example,
exp(1.32) = 3.7, and 95% of activity within 2.448 x 803 = 1,966 m of the centre.

Figures (`results/scr-open/si/fit/figures/`, PNG + PDF):

| figure | caption |
|---|---|
| `abundance_trajectory` | Martens alive each year (line, 50% and 95% bands) and martens detected (black points, a minimum count). |
| `recruitment` | (A) New martens entering each year; (B) per-capita recruitment. |
| `growth_rate` | Year-to-year growth rate N_t / N_t-1; dashed line = stable. |
| `survival` | Posterior of apparent survival between years, with 95% interval. |
| `density_surface` | Map of realised density in each year, including undetected martens; snares as white circles. |
| `density_canopy_cover_pct` | Expected density vs canopy cover (%) in 2020 and 2024, with 95% bands. |
| `detection_function` | Chance a snare detects a marten on one check vs distance from its home-range centre. |
| `detection_by_period` | Martens detected each year (bars) vs the number the model predicts (points, 95% intervals). |
| `gof` | Posterior predictive checks for five summary statistics. |
| `coefficients` | Canopy effect on log density per SD, with 95% interval. |
| `engine_comparison` | Same data fitted by Nimble, NumPyro and openCR. |
| `trace`, `posterior` | MCMC traces and posterior densities of the key quantities. |

Two labels were fixed in the `figures:` block of `config.yaml`: the coefficient
plot showed the raw name `density_beta_canopy_cover_pct` and a clipped axis title,
and the survival axis said "periods" instead of "years". The figures were then
redrawn from the saved fit (no refit).

**9. Gate.** Planted truth passes (28/28). Convergence passes. Posterior predictive checks pass.
A second engine agrees. The module's negative control (state space too small) fails the gate as it should
(0/27, `VALIDATION.md`). Every number comes from the commands below.
**Verdict: DONE**, with the habitat-extrapolation caveat stated.

**10. Deliver** (the summary the researcher receives):

> We estimate 119 martens (95% CrI 85-161) in the 144 km2 study area in 2020 and 186 (144-238) in 2024.
> There is no sign of a decline: the probability of a decline is below 1%.
> Apparent annual survival was 0.84 (0.77-0.91), and 32-50 new martens arrived each year
> (per-capita recruitment 0.23-0.36). Density was strongly associated with canopy cover: 3.7 times higher
> per 12% more canopy (2.9-5.0).
> One caveat matters: your snares sample thinner canopy than the area as a whole. Most of the estimated
> martens live in dense forest away from the snares. So the total numbers depend on the canopy model.
> The trend, survival and recruitment depend much less on it. In future, put some snares in dense canopy.

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% CrI) | interval covers truth? |
|---|---|---|---|
| abundance 2020 | 127 | 119 (85-161) | yes |
| abundance 2021 | 131 | 143 (110-187) | yes |
| abundance 2022 | 133 | 151 (117-195) | yes |
| abundance 2023 | 145 | 162 (125-207) | yes |
| abundance 2024 | 158 | 186 (144-238) | yes |
| superpopulation | 262 | 277 (222-341) | yes |
| recruits 2021 / 2022 / 2023 / 2024 | 26 / 31 / 35 / 43 | 41 (19-72) / 32 (12-60) / 35 (14-62) / 50 (25-82) | yes (all) |
| per-capita recruitment 2021 / 2024 | 0.20 / 0.30 | 0.36 (0.14-0.68) / 0.31 (0.15-0.52) | yes |
| growth rate 2021 / 2024 | 1.03 / 1.09 | 1.21 (0.98-1.55) / 1.15 (0.98-1.37) | yes |
| survival | 0.80 | 0.839 (0.771-0.905) | yes |
| baseline detection p0 | 0.30 | 0.317 (0.282-0.354) | yes |
| sigma (m) | 800 | 803 (764-846) | yes |
| home range area (km2) | 12.0 | 12.1 (11.0-13.5) | yes |
| canopy effect (log, per SD) | 1.00 | 1.32 (1.06-1.61) | **no** |

Every interval but one contains the truth. The canopy effect misses. This is a
property of this one simulated population, not a model failure:
- Fitting the density model to the *true* home-range centres gives 1.18 (SE 0.08). So these martens really do favour canopy more than the simulation value.
- The snares sample thin canopy, so the steepness of the effect is partly extrapolated.
- The 100-dataset recovery shows the canopy effect is unbiased at this design (+1.0%, coverage 0.98).

The abundance estimates for 2021-2024 are 9-18% above the truth. That is
consistent with the canopy effect being estimated too steep.

## (f) Limits from VALIDATION.md that matter at this design

- **Snares in unrepresentative habitat.** Total abundance rests on the canopy model (extrapolation factor 1.94; the report warns above 1.5). A homogeneous-density fit (openCR) gives about 1.7-1.9 x lower abundance.
- **Static home-range centres.** Martens that shift ranges between years look like deaths plus recruits. Survival would then be biased low.
- **Constant survival and detection.** Sex/age differences, trap response and year-to-year detection changes are not modelled. The posterior predictive checks found no sign of them here, but they were not validated against such violations.
- **Apparent survival.** Death and permanent emigration cannot be told apart.
- **Recruitment.** Recruitment in the first year cannot be estimated. The last interval's survival and recruitment are less precise.
- **Equal intervals.** Survival is per year. A skipped summer would break this.
- **Sample size.** At ~40 detected animals (about 20 per year), abundance is 4-8% high and per-capita recruitment 12-20% high. This study has 95 animals (32-48 per year), where recovery is within ~5%.
- **Engines.** Recovery was validated with Nimble marginal only. NumPyro agrees on this dataset.

## (g) Commands to reproduce

Run from `ecotools/`:

```bash
# 0. (Already done) regenerate the example data from the simulator, seed 2025
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/scr-open_make_example.log -- scripts/Rscript dev/scr-open/make_example.R

# 1. Planted-truth recovery at the example design (100 datasets, ~45 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/scr-open_si_recover.log -- \
  scripts/Rscript core/run.R recover scr-open --config modules/scr-open/example/recovery.yaml \
  --nsim 100 --cores 16 --out results/scr-open/si/recovery

# 2. Fit the example data, with report, figures and engine cross-checks (~75 min on 4 cores)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/scr-open_si_fit.log -- \
  scripts/Rscript core/run.R fit scr-open --config modules/scr-open/example/config.yaml \
  --out results/scr-open/si/fit

# 3. Redraw the figures after editing the figures: block (no refit)
scripts/ejob -c 1 -m 8 -t 0:30:00 -l logs/scr-open_si_figures.log -- \
  scripts/Rscript core/run.R figures scr-open --config modules/scr-open/example/config.yaml \
  --fit results/scr-open/si/fit/scr-open_fit.rds
```

Outputs: `results/scr-open/si/recovery/` (gate table, per-replicate CSV, HTML),
`results/scr-open/si/fit/` (report, estimates CSV, saved fit, figures),
`results/scr-open/si/si_summary.json` (machine-readable summary).
