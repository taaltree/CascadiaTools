# Worked example: which forest birds respond to temperature and canopy, and who occurs with whom?

Module: `jsdm` (joint species distribution model, HMSC). Appendix S1 worked example.

> The data in this example are **simulated** (see `README_data.md`). The true values are known,
> so every estimate below can be checked against the truth (section 5).

All numbers below come from these runs (2026-10-20, R 4.5.3, Hmsc 3.3.7):
- fit: `results/jsdm/si/fit/` (report `jsdm_report.html`, figures in `figures/`)
- planted-truth recovery: `results/jsdm/si/recovery/`
- truth and comparison: `results/jsdm/si/truth_full.csv`, `results/jsdm/si/estimates_vs_truth.csv`
- machine-readable summary: `results/jsdm/si/si_summary.json`

---

## 1. The question, in the researcher's words

> "We did point counts for breeding birds at 120 stations in the western Cascades, three visits
> each. I want to know how each species responds to summer temperature and canopy cover. Do big
> birds respond differently from small ones? And which species turn up together more, or less,
> often than you'd expect from habitat alone? I don't code."

## 2. The files the researcher has

Three CSV files, in the format a field crew keeps. Nothing was reshaped by hand.

`bird_detections.csv`: one row per species detected on a visit (1,446 rows)

```
"station","visit","date","species","count"
"CAS-001",1,"2025-05-22","Setophaga occidentalis",1
"CAS-001",1,"2025-05-22","Troglodytes pacificus",1
"CAS-001",1,"2025-05-22","Certhia americana",1
"CAS-001",1,"2025-05-22","Ixoreus naevius",2
```

`stations.csv`: one row per station (120 rows)

```
"station","utm_e","utm_n","elevation_m","july_temp_c","canopy_cover_pct","survey_year"
"CAS-001",576563,4880362,118,22.5,71,2025
"CAS-002",576764,4881624,503,21.8,48,2025
"CAS-003",576156,4880325,510,20.5,28,2025
"CAS-004",577719,4880443,226,23.4,63,2025
```

`bird_traits.csv`: one row per species (15 rows)

```
"species","common","family","order","body_mass_g"
"Setophaga occidentalis","Hermit Warbler","Parulidae","Passeriformes",9.3
"Cardellina pusilla","Wilson's Warbler","Parulidae","Passeriformes",7.5
"Geothlypis tolmiei","MacGillivray's Warbler","Parulidae","Passeriformes",10
"Troglodytes pacificus","Pacific Wren","Troglodytidae","Passeriformes",9
```

## 3. What the researcher types to Claude

> I have bird point-count data in `modules/jsdm/example/`: `bird_detections.csv` (every detection,
> 3 visits per station), `stations.csv` (coordinates, July temperature, canopy cover) and
> `bird_traits.csv` (body mass and family for each species). How does each species respond to
> July temperature and canopy cover? Does body size explain who responds how? Which species occur
> together more or less than the habitat predicts? Please give me a report with figures I can put
> in a paper.

## 4. What the agent does, step by step (skill `eco-analyze`)

### Step 1. Project folder
The agent keeps every file and command in one place. For this appendix the outputs are in
`results/jsdm/si/` and the scripts in `dev/jsdm/`. The raw CSVs are never edited.

### Step 2. Data intake (`eco-data-intake`)
The agent reads the three files and describes them back:

- 120 stations, all with temperature and canopy values. No missing values.
- 15 species. Every species in the detections is in the traits table.
- 1,446 detection rows over 3 visits. A species seen on several visits at a station is counted once
  (detected / not detected).
- Species were recorded at 14% (Chestnut-backed Chickadee, 17 stations) to 74% (Varied Thrush,
  89 stations) of stations, median 32%.
- July temperature 12.8 to 25.0 °C; canopy cover 5 to 98%. The two are only weakly correlated
  (r = -0.27), so their effects can be separated.
- Elevation is in the file but is not used: temperature is derived from it, and using both would
  duplicate the same gradient.

### Step 3. Choose the module
From `catalog.md`, the question ("how do many species respond to the same gradients, do traits
explain it, and who co-occurs") is exactly the `jsdm` module's `when_to_use`.
The agent notes one alternative. With 3 visits per station, a multi-species occupancy model
(`occu-spatial`) could separate "absent" from "present but missed". The researcher's questions are
about traits and co-occurrence, which only `jsdm` answers. The agent therefore chooses `jsdm` and
records the detection caveat (section 6).

### Step 4. Analysis plan (`eco-analysis-plan`)
Declared before fitting:

| item | choice |
|---|---|
| response | detected / not detected per station and species (probit link) |
| covariates | July mean max. temperature, canopy cover (both standardised: effects are per 1 SD) |
| trait | log body mass (standardised) |
| shared ancestry | taxonomic tree from order and family |
| residual associations | up to 5 latent site factors |
| estimands interpreted | population-averaged species responses, community averages, trait effects (with Gamma as the test for the wider species pool), latent-scale residual correlations, expected richness |
| adequate fit | R-hat < 1.01 and bulk ESS >= 400 for every reported quantity; posterior predictive P between 0.05 and 0.95 |
| MCMC | 4 chains x 1,000 samples, thin 50, burn-in 25,000 (module defaults) |

### Step 5. Planted truth at this design (`eco-planted-truth`)
Before touching the real data, the agent simulates 100 datasets that look like this study:
120 stations, 15 species, 2 covariates with correlation -0.27, 1 trait, 2 true latent factors,
presence/absence. Values are plausible community values (`example/recovery.yaml`). Each simulated
dataset is fitted exactly as the real data will be, and the gate checks bias and 95% interval
coverage for every interpreted quantity.

One fit takes about 7-10 minutes (module selftest: about 7 min per replicate), so the agent runs
100 replicates instead of 200, on 16 cores (about 75 min in total):

```
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/jsdm_si_recover.log -- scripts/Rscript core/run.R recover jsdm \
  --config modules/jsdm/example/recovery.yaml --nsim 100 --cores 16 --out results/jsdm/si/recovery
```

Result: 0 errors, 98 of 100 replicates converged (max R-hat < 1.05; median 522 s per replicate).

| estimand class | rows | gate pass | mean coverage (min) | median abs. rel. bias |
|---|---|---|---|---|
| population-averaged species responses `popavg_*` (interpreted) | 30 | **30/30** | 0.942 (0.91) | 20%* |
| community averages `comm_popavg_*` (interpreted) | 2 | **2/2** | 0.934 (0.92) | 0.8% |
| trait effects `popavg_trait_*` (interpreted) | 3 | **3/3** | 0.963 (0.95) | 1.8% |
| latent residual correlations `latcor_*` (interpreted) | 105 | **102/105** | 0.954 (0.88) | 32%* |
| prevalences `prev_*` (interpreted) | 15 | **15/15** | 0.952 (0.92) | 1.5% |
| richness at average conditions (interpreted) | 1 | **1/1** | 0.959 | 0.0% |
| Gamma community means `comm_mean_*` (test for the species pool) | 3 | 3/3 | 0.959 | 9.4% |
| Gamma trait effects `trait_*` (test for the species pool) | 3 | 1/3 | 0.946 | 15.8% |
| Hmsc conditional Beta (reported, not interpreted) | 45 | 32/45 | 0.953 | 44%* |
| Hmsc `rescor_*` (reported, not interpreted) | 105 | 0/105 | 0.760 | 20% |

\* Relative bias against per-species truths that are averaged over random communities and often
near 0. In those cases the gate uses the 3-Monte-Carlo-SE rule.

Community-level rows (truth / mean estimate / relative bias / coverage / power):

| param | truth | mean est. | rel. bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| comm_popavg_temperature | 0.373 | 0.378 | +0.8% | 0.92 | 1.00 | PASS |
| comm_popavg_forest (= canopy) | 0.301 | 0.300 | -0.8% | 0.95 | 0.97 | PASS |
| popavg_trait on intercept | 0.264 | 0.276 | -1.8%** | 0.97 | - | PASS |
| popavg_trait on temperature | -0.356 | -0.362 | -1.2%** | 0.97 | 0.99 | PASS |
| popavg_trait on canopy | 0.220 | 0.224 | -3.1%** | 0.95 | 0.87 | PASS |
| richness_mean_env | 5.93 | 5.93 | 0.0% | 0.96 | - | PASS |
| Gamma comm_mean_temperature | 0.50 | 0.544 | +8.8% | 0.96 | 0.89 | PASS |
| Gamma comm_mean_forest (= canopy) | 0.40 | 0.438 | +9.4% | 0.96 | 0.77 | PASS |
| Gamma trait on intercept | 0.35 | 0.405 | +15.8% | 0.94 | 0.60 | FAIL (bias) |
| Gamma trait on temperature | -0.45 | -0.531 | +18.0% | 0.95 | 0.82 | FAIL (bias) |
| Gamma trait on canopy | 0.30 | 0.321 | +7.1% | 0.95 | 0.45 | PASS |

\*\* For the realised (per-community) quantities the truth changes in every replicate. The truth column is
its average, and relative bias is averaged per replicate. It therefore need not equal (mean est. - truth) / truth.

Sign recovery (all converged replicates pooled): among population-averaged responses whose 95% CrI
excludes 0, the sign was correct 99.7% of the time, and 93% of responses with |truth| > 0.3 were
detected. Among supported residual correlations, the sign was correct 99.8% of the time. 80% of pairs
with |r| > 0.3 were detected, and only 2.2% of pairs with |r| < 0.05 were flagged. Strong
correlations are shrunk toward 0 (slope 0.85).

**Gate result in one sentence (what the agent tells the researcher):** "At your design, the model
recovered every species response, community average, trait effect, prevalence and richness value
without bias and with honest 95% intervals. 3 of the 105 species-pair residual correlations failed
the check, all of them pairs with true correlation near zero, so I will not interpret individual
weak associations."

The Gamma bias failures (+16% and +18%) are the known inflation of Hmsc's conditional scale
(VALIDATION.md limit 2). Their coverage is nominal, so the report uses Gamma only for "does the
effect extend to the species pool", never for effect size.

### Step 6. Fit
```
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/jsdm_si_fit.log -- \
  scripts/Rscript core/run.R fit jsdm --config modules/jsdm/example/config.yaml --out results/jsdm/si/fit
```
MCMC time 952 s (about 16 minutes) on 4 cores.

### Step 7. Checks (`eco-convergence`, `eco-model-checking`)
- **Convergence:** max R-hat 1.004 and minimum bulk ESS 1,300 over every reported quantity
  (targets < 1.01 and >= 400). Trace plots show the 4 chains overlapping with no trends
  (`trace.png`). Hmsc kept all 5 latent factors.
- **Goodness of fit (posterior predictive check at new sites):** the observed variance of site
  richness (4.52) falls at P = 0.71 of the replicated distribution; mean pairwise co-occurrence at
  P = 0.67. Both are inside 0.05-0.95, so the model reproduces the community structure.
- **Explanatory power:** mean Tjur R² 0.46 (range 0.15-0.66), mean AUC 0.96. This is in-sample and
  includes the latent factors, so it is optimistic.

### Step 8. Interpretation (`eco-interpretation`)
The report (`results/jsdm/si/fit/jsdm_report.html`) writes its own interpretation. The agent
checked each sentence against the estimates table. Quoted excerpts:

> July mean max. temperature (°C): 8 of 15 species respond positively and 0 negatively with at
> least 95% posterior support. The average response across these species is 0.26 per 1 SD of July
> mean max. temperature (°C) (1 SD = 2.47 in the units of july_temp_c) on the probit
> (population-averaged) scale (95% CrI 0.19 to 0.33). For a species with the average occurrence
> probability (0.35 at average conditions), that is a change to 0.45 (95% CrI 0.43 to 0.48) per SD
> increase. [...] the hierarchical community mean (Gamma [...]) is 0.37 (95% CrI -0.06 to 0.78),
> which does not exclude zero: a new species from the same pool could respond either way.

> Canopy cover (%): 9 of 15 species respond positively and 1 negatively with at least 95% posterior
> support. The average response across these species is 0.30 per 1 SD of Canopy cover (%) (1 SD = 23
> in the units of canopy_cover_pct) [...] (95% CrI 0.23 to 0.38).

> Expected species richness (number of these species present) is 4.2 (95% CrI 3.7-4.7) at July mean
> max. temperature (°C) = 15.5 (10th percentile of sites) and 7.1 (6.6-7.7) at 21.8 (90th
> percentile), other covariates at their mean. This is an association, not evidence that July mean
> max. temperature (°C) causes the change.

> Traits (across the modelled species): species with higher Body mass (log g) have a higher baseline
> prevalence (0.33 [...], 95% CrI 0.26 to 0.39; pool-level trait effect Gamma 0.50, 95% CrI 0.12 to
> 0.88, which also excludes zero); species with higher Body mass (log g) have a more negative
> response to July mean max. temperature (°C) (-0.20 [...], 95% CrI -0.27 to -0.13; pool-level trait
> effect Gamma -0.26, 95% CrI -0.55 to 0.04, which does not exclude zero, so this pattern may not
> hold beyond these species); species with higher Body mass (log g) have a more positive response to
> Canopy cover (%) (0.25 [...], 95% CrI 0.18 to 0.32; pool-level trait effect Gamma 0.39, 95% CrI
> 0.06 to 0.75, which also excludes zero).

> Residual associations: of 105 species pairs, 39 co-occur more and 37 less often than expected from
> the covariates, with at least 95% posterior support. Strongest: Hairy Woodpecker & Spotted Towhee
> (r = -0.71, 95% CrI -0.86 to -0.51); Spotted Towhee & Red-breasted Nuthatch (r = 0.69, 95% CrI
> 0.51 to 0.84) [...]. These residual correlations [...] can arise from unmeasured habitat variables,
> shared history or sampling, as well as from biotic interactions, so they are not evidence of
> interactions (Blanchet et al. 2020).

> Variance partitioning (average over species, of the explained variance on the linear-predictor
> scale): July mean max. temperature (°C) 17%; Canopy cover (%) 20%; Residual associations (site
> latent factors) 63%.

Also reported: body mass explains 34% (temperature) and 41% (canopy) of the among-species variation
in responses. Phylogenetic (taxonomic) signal rho = 0.31 (95% CrI 0.00 to 0.96), i.e. the data say
nothing about it.

**In plain words for the researcher.** Most species are recorded more often at warmer and more
closed-canopy stations. A typical species goes from 35% to about 45% of stations per 1 SD
(2.5 °C) of warming, and to about 47% per 1 SD (23 percentage points) of extra canopy. Golden-crowned
Kinglet is the only species with a supported negative response (to canopy). Larger species are more
widespread and respond more strongly to canopy. For these 15 species, larger species also respond
less to temperature, but that may not hold for other species. Many species pairs co-occur more or less
than habitat predicts. Most of that variation (63%) is not explained by the two covariates, which
points to unmeasured habitat structure rather than to interactions.

**Figures** (`results/jsdm/si/fit/figures/`, PDF + PNG; all inspected)

| id | caption (short) |
|---|---|
| `species_responses` | Heat map of each species' population-averaged response (per SD) to temperature and canopy; coloured = >= 95% support. |
| `residual_associations` | Species x species latent-scale residual correlations; coloured = >= 95% support. Not evidence of interactions. |
| `variance_partitioning` | Share of each species' explained variance due to temperature, canopy and latent factors. |
| `trait_effects` | Community-average responses and body-mass effects on responses, with 50% and 95% CrIs. |
| `richness_july_temp_c` | Expected richness along July temperature, with observed richness per station. |
| `richness_canopy_cover_pct` | Expected richness along canopy cover, with observed richness per station. |
| `occurrence_july_temp_c` | Probability of occurrence of each species along temperature (15 panels). |
| `occurrence_canopy_cover_pct` | Probability of occurrence of each species along canopy cover (15 panels). |
| `model_fit` | Tjur R² per species against observed prevalence. |
| `ppc` | Posterior predictive check: replicated vs observed variance of site richness. |
| `trace` | MCMC traces of the community-level parameters, 4 chains. |
| `posterior` | Posterior densities of the community-level parameters. |

One presentation fix was made, in the config's `figures:` block only: in
`occurrence_canopy_cover_pct` the "100" tick label ran into the neighbouring panel, so the x range
was widened (`per_figure: occurrence_canopy_cover_pct: {xlim: [0, 110]}`) and the figures were
redrawn from the saved fit with `core/run.R figures`.

### Step 9. Acceptance gate (`eco-acceptance-gate`)

| gate step | result |
|---|---|
| 0. declared before fitting | analysis plan above (estimands, covariates, priors, adequacy targets) |
| 1. planted truth at this design | 153 of 156 interpreted rows pass; 98/100 converged. **3 near-zero `latcor` rows fail** |
| 2. convergence on the data | max R-hat 1.004, min bulk ESS 1,300: PASS |
| 3. model checking | PPC P = 0.71 (richness variance), 0.67 (co-occurrence): PASS |
| 4. independent route | module cross-check against glmmTMB (VALIDATION.md section 7: popavg r = 0.986, latcor r = 0.918) |
| 5. negative control | module run omitting a covariate fails the gate (VALIDATION.md section 6) |
| 6. reproducible | config, recovery config, data and all commands in this folder and section 7 |

**Verdict: OPEN (narrowly).** The gate requires every interpreted estimand to pass at the study's
design, and its thresholds are never widened. Three of 105 residual-correlation rows failed. Two were
on coverage (0.878 against a threshold of 0.884 for 98 replicates). One was on bias (+0.040 against a
truth of -0.045, 3.2 Monte Carlo SE). All three belong to pairs whose true correlation is near zero. If
coverage were exactly nominal, about 0.4 such failures would be expected by chance across 105 pairs,
so this is a small real under-coverage for weak pairs, not just noise.

What this means for the researcher:
- The species responses, community averages, body-mass effects, prevalences and richness gradients
  passed every check and can be reported as they are.
- The strong supported associations can be reported as residual co-occurrence patterns, with the
  usual caveat that they are not interactions. In recovery, supported signs were correct 99.8% of the time.
- Do not interpret individual weak associations (|r| below about 0.2), even when they are flagged.
- To close the gate: re-run the recovery with more replicates (nsim 200, to tell real under-coverage
  apart from chance), or interpret associations only above an |r| threshold that is declared in advance.

### Step 10. Deliver
The researcher gets the HTML report, the 12 figures (PDF for the journal, PNG for slides), the
estimates table (`jsdm_estimates.csv`) and the methods paragraph printed in the report.

## 5. Estimates versus the known truth

Because the data are simulated, the truth is known (`README_data.md`; rebuilt from seed 20261019 by
`dev/jsdm/si_truth.R`; comparison by `dev/jsdm/si_compare.R`).

| quantity | truth | estimate (95% CrI) | covered? |
|---|---|---|---|
| community average response, temperature (population-averaged) | 0.227 | 0.264 (0.193 to 0.334) | yes |
| community average response, canopy (population-averaged) | 0.257 | 0.302 (0.225 to 0.380) | yes |
| body mass on prevalence (realised, population-averaged) | 0.327 | 0.330 (0.265 to 0.393) | yes |
| body mass on temperature response (realised) | -0.223 | -0.199 (-0.267 to -0.132) | yes |
| body mass on canopy response (realised) | 0.245 | 0.248 (0.178 to 0.317) | yes |
| Gamma: community mean, temperature | 0.50 | 0.37 (-0.06 to 0.78) | yes |
| Gamma: community mean, canopy | 0.40 | 0.43 (-0.06 to 0.88) | yes |
| Gamma: body mass on intercept | 0.35 | 0.50 (0.12 to 0.88) | yes |
| Gamma: body mass on temperature | -0.45 | -0.26 (-0.55 to 0.04) | yes |
| Gamma: body mass on canopy | 0.30 | 0.39 (0.06 to 0.75) | yes |
| expected richness at average conditions | 5.65 | 5.27 (4.94 to 5.62) | no (see below) |

By class:
- **Species responses (population-averaged):** 29 of 30 intervals cover the truth (97%);
  correlation with the truth 0.94. All 16 supported responses have the correct sign.
- **Residual correlations (latent scale):** 99 of 105 intervals cover the truth (94%);
  correlation with the truth 0.95. 64 of 65 supported pairs have the correct sign. 56 of the
  62 pairs with a true |r| > 0.3 were detected with the correct sign. 2 of 15 pairs with a true
  |r| < 0.1 were flagged. Strong correlations were shrunk toward zero (slope 0.81), as
  VALIDATION.md documents.
- **Prevalences:** 14 of 15 intervals cover. Western Tanager (truth 0.47, estimate 0.34, CrI
  0.23-0.44) is the exception. Its occurrence switches steeply with canopy, so the value at
  average canopy is sensitive.
- **Richness:** the truth (5.65) lies just above the interval (4.94 to 5.62). The truth is true
  occupancy. The data are detections, and 2% of occupied station x species combinations were
  never detected over 3 visits. Expected *recorded* richness is about 5.65 x 0.973 = 5.50, inside
  the interval. This is the detection caveat in action, not a model failure.
- **Hmsc's conditional Beta** (not interpreted): 43 of 45 cover.

## 6. Limits that matter at this design (from VALIDATION.md)

1. **Detection is not modelled.** The model describes where birds were *recorded*. Here 2% of
   occupied combinations were missed, which is enough to put true richness just outside its interval.
   With fewer visits, or with canopy affecting detectability, "responses" would mix habitat and detection.
2. **Residual associations are not interactions.** Leaving out one important covariate produced
   spurious supported "associations" in 57% of species pairs in simulation. Here 63% of explained
   variance is in the latent factors, so unmeasured habitat probably drives most of them.
3. **Strong residual correlations are shrunk toward zero** by about 15-20% (slope 0.81 here).
4. **Gamma (pool-level) has low power with 15 species.** Neither community-mean Gamma excludes zero
   here, although the realised averages for these 15 species clearly do. Adding sites does not fix this;
   adding species does.
5. **Hmsc's conditional Beta and Gamma are inflated by 10-25%.** Report effect sizes on the
   population-averaged scale, as the report does.
6. **The taxonomic tree (rho) and spatial random effects have no recovery validation.** Here rho is
   uninformative (0.31, CrI 0.00-0.96). The recovery runs had no tree, because the truth has no phylogenetic signal.
7. **In-sample Tjur R² / AUC are optimistic.** No cross-validation is implemented.

## 7. Reproduce everything

From `.`:

```
# 0. (optional) regenerate the example data (seed 20261019)
scripts/ejob -c 1 -m 4 -t 0:15:00 -- scripts/Rscript modules/jsdm/example/make_example.R
# 1. full truth (popavg, latcor, prevalence, richness) from the same seed
scripts/ejob -c 1 -m 4 -t 0:15:00 -l logs/jsdm_si_truth.log -- scripts/Rscript dev/jsdm/si_truth.R
# 2. planted-truth recovery at this design (100 replicates, about 65 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/jsdm_si_recover.log -- scripts/Rscript core/run.R recover jsdm \
  --config modules/jsdm/example/recovery.yaml --nsim 100 --cores 16 --out results/jsdm/si/recovery
scripts/ejob -c 1 -m 4 -t 0:15:00 -l logs/jsdm_si_recsumm.log -- scripts/Rscript dev/jsdm/summ.R results/jsdm/si/recovery recover
# 3. fit the example data (report + figures, about 16 min)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/jsdm_si_fit.log -- scripts/Rscript core/run.R fit jsdm \
  --config modules/jsdm/example/config.yaml --out results/jsdm/si/fit
# 4. redraw figures from the saved fit (after editing the config's figures: block)
scripts/ejob -c 2 -m 8 -t 0:30:00 -l logs/jsdm_si_figs.log -- scripts/Rscript core/run.R figures jsdm \
  --fit results/jsdm/si/fit/jsdm_fit.rds --config modules/jsdm/example/config.yaml
# 5. estimates vs truth, and the SI summary JSON
scripts/ejob -c 1 -m 4 -t 0:15:00 -l logs/jsdm_si_compare.log -- scripts/Rscript dev/jsdm/si_compare.R
scripts/ejob -c 1 -m 4 -t 0:15:00 -l logs/jsdm_si_summary.log -- scripts/Rscript dev/jsdm/si_summary.R
```

Do not run the R commands directly on the login node; use `scripts/ejob` as shown (plain `srun`
hangs inside the interactive allocation).
