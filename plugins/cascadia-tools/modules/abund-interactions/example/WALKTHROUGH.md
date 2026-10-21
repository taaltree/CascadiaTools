# Worked example: do more pine martens mean fewer red squirrels?

Module: `abund-interactions` (abundance-mediated species interactions).
This example uses **simulated** data, so we know the true answer and can check
the tool against it. See `README_data.md` for the full story and the true values.

---

## (a) The question, in the researcher's words

> "We put a camera trap at 260 woodland sites in Northern Ireland last spring.
> Pine martens are back, and we want to know if red squirrels are less likely
> to be at sites where there are more martens. We also think red squirrels
> cope better in native broadleaf woods than in conifer plantations. Does the
> marten effect depend on how much broadleaf there is?"

## (b) The files they have

Three CSV files, as exported from their camera software and GIS.

`deployments.csv`: when each camera was working (279 rows).
```
"site","camera_id","start_date","end_date"
"NI001","C001","2024-03-02","2024-03-29"
"NI002","C002","2024-03-05","2024-04-01"
"NI003","C003","2024-03-05","2024-04-01"
"NI004","C004","2024-03-07","2024-04-03"
```

`records.csv`: one row per photo of an animal (2,377 rows, five species).
```
"site","date","species","count"
"NI001","2024-03-03 11:12:07","Red squirrel",1
"NI001","2024-03-03 18:29:43","Red squirrel",1
"NI001","2024-03-03 08:48:36","Badger",1
"NI001","2024-03-05 08:35:35","Red squirrel",1
```

`sites.csv`: one row per site, with habitat measured in GIS (260 rows).
```
"site","easting","northing","woodland_pct","broadleaf_pct","dist_settlement_km"
"NI001",330600,341600,40.7,22.3,0.29
"NI002",237800,328100,38.8,30,0.77
"NI003",302500,391900,54.5,23.2,0.75
"NI004",305700,413800,31.4,0,0.5
```

## (c) What they type to Claude

> I have camera-trap data from 260 woodland sites in Northern Ireland:
> deployments.csv (camera start and end dates), records.csv (every photo, with
> species) and sites.csv (woodland cover and broadleaf share within 1 km).
> Is red squirrel occupancy lower where there are more pine martens, and does
> that depend on how much of the woodland is native broadleaf? Please give me
> a report with figures I can use in a paper.

## (d) What the agent does, step by step (skills/eco-analyze)

**1. Project folder.** The agent copies the files into a project folder and
keeps a log of every command. The raw files are never changed.

**2. Data intake** (eco-data-intake). The agent reads the files and describes
them back:
- 260 sites, 279 camera deployments, 6,939 camera-days between 1 March and
  17 April 2024. 19 sites had a memory-card swap with a gap.
- Camera-days are grouped into 4 weekly occasions of 7 days. 8.9% of
  site-weeks have lost days; 3.0% have none.
- Pine marten photographed at 142 sites (55%). Red squirrel at 130 sites (50%).
- Badger, red fox and grey squirrel records are present but not needed.
- Woodland cover: 5-95% (mean 46%). Broadleaf share: 0-100% (mean 39%).
- `dist_settlement_km` is in the file but the researcher did not ask about it,
  so it is left out.

**3. Choose the module.** Martens cannot be counted directly, but how often a
camera sees *any* marten rises with how many martens use the site. The
`abund-interactions` module (Twining et al. 2025) estimates marten abundance
at each site this way (a Royle-Nichols model). It then uses that abundance,
with its uncertainty, to predict red squirrel occupancy. A plain two-species
occupancy model would only use marten *presence*, not *how many*.

**4. Plan** (eco-analysis-plan). Written down before any fitting:
- Marten abundance depends on woodland cover.
- Red squirrel occupancy depends on woodland cover, broadleaf share and the
  number of martens.
- The marten effect is allowed to change with broadleaf share.
- Detection is constant per day for each species.
- Main quantities: gamma0 (the marten effect at average broadleaf share) and
  gamma1 (how much that effect changes per SD of broadleaf share).
- Adequate fit: R-hat < 1.01, effective sample size (ESS) >= 400, and a
  posterior predictive p-value between 0.05 and 0.95 for both species.

**5. Planted-truth check** (eco-planted-truth). Before touching the real
data, the agent simulates 100 data sets with the same design (260 sites,
4 weeks of 7 days, 9% lost effort). It plants realistic values, fits each one
and checks the answers. Each fit takes about 2.7 minutes, so 100 replicates
were used rather than 200.

| parameter | truth | mean estimate | rel. bias | coverage | power | gate |
|---|---|---|---|---|---|---|
| marten abundance (average woodland) | 1.0 | 0.99 | -1.4% | 0.96 | | PASS |
| woodland -> marten abundance | 0.5 | 0.51 | +1.5% | 0.94 | 1.00 | PASS |
| marten daily detection (per animal) | 0.08 | 0.081 | +1.7% | 0.97 | | PASS |
| squirrel occupancy with no martens | 0.6 | 0.61 | +1.6% | 0.92 | | PASS |
| woodland -> squirrel occupancy | 0.5 | 0.55 | +10.2% | 0.92 | 0.86 | PASS |
| broadleaf -> squirrel occupancy | 0.3 | 0.37 | +21.8% | 0.95 | 0.34 | PASS* |
| **gamma0: per marten, average habitat** | **-0.6** | **-0.67** | **-11.0%** | **0.95** | **0.97** | **FAIL** |
| gamma1: per marten x broadleaf | 0.6 | 0.65 | +9.1% | 0.96 | 0.92 | PASS |
| squirrel daily detection | 0.15 | 0.150 | +0.2% | 0.95 | | PASS |
| total marten-site abundance | 293.6 | 292.8 | -0.3% | 0.98 | | PASS |
| share of sites with martens | 0.63 | 0.63 | -0.2% | 0.92 | | PASS |
| share of sites with squirrels | 0.46 | 0.46 | +0.2% | 1.00 | | PASS |

\* Passes because its bias is within 3 Monte Carlo SE. With power 0.34 it is
weakly estimated at this design.

All 100 replicates converged. **The gate passed 11 of 12 quantities.** It
failed for gamma0: the marten effect at average habitat was overstated by
11% (limit 10%; bias 3.2 times its Monte Carlo SE). Its 95% intervals still
contained the truth 95% of the time. This matches a known limit of the module
(fewer than about 600 sites with a habitat-dependent interaction). The agent
tells the researcher plainly: *the interval is trustworthy; the point estimate
may make the marten effect look about 10% stronger than it is.*

**6. Fit.** The agent writes `config.yaml` and fits the model: 3 chains of
20,000 iterations (5,000 burn-in, thin 5), about 8 minutes on 4 cores.

**7. Check** (eco-convergence, eco-model-checking).
- Converged: largest R-hat 1.001, smallest bulk ESS 3,279 (needs >= 400).
  The trace plots show well-mixed chains.
- Goodness of fit: posterior predictive p = 0.51 (pine marten) and 0.32
  (red squirrel). Both are far from 0 and 1, so there is no sign of lack of fit.

**8. Interpret** (eco-interpretation). Excerpts from the report:

> "Interaction at average habitat: each additional Pine marten individual at a
> site is associated with 0.59 (95% CI 0.39 to 0.85) times the odds that Red
> squirrel occupies the site (odds ratio per individual) ... The posterior
> probability that the association is negative is > 0.99."

> "Predicted Red squirrel occupancy at average habitat: 0.62 (0.52-0.72) at a
> site with no Pine marten, 0.50 (0.42-0.57) with 1, 0.37 (0.25-0.49) with 2
> and 0.26 (0.12-0.44) with 3 individuals."

> "Habitat dependence: the interaction changes by 0.69 (95% CI 0.30 to 1.14)
> logit units per individual for each 1-SD (20.8 units) increase in
> broadleaf_pct; the 95% interval excludes zero, so the data support an
> interaction whose strength depends on habitat."

> "Where broadleaf_pct = 12.9 (10th percentile): odds ratio per Pine marten
> individual 0.25 (95% CI 0.12 to 0.50) (negative). Where broadleaf_pct = 65.2
> (90th percentile): odds ratio per Pine marten individual 1.42 (95% CI 0.84
> to 2.53) (not clearly different from zero)."

> "Caution: this is an association between the dominant's estimated abundance
> and the subordinate's occupancy, not proof of predation or competition."

In plain words: in conifer-dominated woods, each extra marten cuts the odds
of finding red squirrels by about three-quarters. In broadleaf-rich woods
there is no clear marten effect.

**Figures** (in `results/abund-interactions/si/fit/figures/`, PNG and PDF):

| figure | caption (short) |
|---|---|
| `interaction_by_habitat` | Red squirrel occupancy against marten numbers at the 10th, 50th and 90th percentiles of broadleaf share. The slopes differ: steeply down in conifer woods, flat to rising in broadleaf woods. **Key figure.** |
| `subordinate_vs_dominant` | Red squirrel occupancy against marten numbers at average habitat (median, 50% and 95% bands). |
| `interaction_strength_broadleaf_pct` | Marten effect (logit per marten) across broadleaf share; below zero = fewer squirrels with more martens. The estimate crosses zero near 55% broadleaf (truth: 59%). |
| `dominant_abundance_woodland_pct` | Expected martens per site against woodland cover. |
| `subordinate_occupancy_woodland_pct` | Red squirrel occupancy against woodland cover with 0 and 2 martens. |
| `dominant_abundance_map` | Estimated martens per site on a map. |
| `coefficients` | All regression coefficients with 95% intervals. |
| `joint_vs_twostep` | Joint model against the "estimate martens first, then plug in" shortcut. |
| `gof` | Posterior predictive check for each species (p = 0.51 and 0.32). |
| `trace` | MCMC chains for the key parameters. |
| `posterior` | Posterior densities for the key parameters. |

**9. Gate** (eco-acceptance-gate). **Verdict: OPEN.** Convergence, model
checking, cross-checks (Nimble vs own likelihood code vs `unmarked`), the
negative control and reproducibility all pass. The planted-truth check at
this design fails for gamma0 (bias -11%). The analysis can be reported, but
the size of the marten effect at average habitat should be read from its
interval, not its point estimate. Doubling the survey to about 600 sites
would remove the problem.

**10. Deliver.** A short summary, the HTML report
(`results/abund-interactions/si/fit/abund-interactions_report.html`) and the
methods paragraph from the report.

## (e) Estimates versus the known truth

| quantity | truth | estimate (95% CrI) | truth inside? |
|---|---|---|---|
| martens per site, average woodland | 1.0 | 0.85 (0.70-1.03) | yes |
| woodland -> marten abundance (log, per SD) | 0.5 | 0.53 (0.39-0.67) | yes |
| marten daily detection per animal | 0.08 | 0.084 (0.071-0.098) | yes |
| squirrel occupancy, no martens | 0.6 | 0.62 (0.53-0.72) | yes |
| woodland -> squirrel occupancy (logit, per SD) | 0.5 | 0.49 (0.16-0.83) | yes |
| broadleaf -> squirrel occupancy (logit, per SD) | 0.3 | 0.02 (-0.37 to 0.40) | yes |
| gamma0: per marten, average habitat | -0.6 | -0.53 (-0.93 to -0.17) | yes |
| gamma1: per marten x broadleaf (per SD) | 0.6 | 0.69 (0.30-1.14) | yes |
| per marten at 12.9% broadleaf | -1.34 | -1.39 (-2.16 to -0.70) | yes |
| per marten at 65.2% broadleaf | +0.16 | +0.35 (-0.17 to 0.93) | yes |
| squirrel daily detection | 0.15 | 0.152 (0.140-0.165) | yes |
| total marten-site abundance (realised) | 274 | 256 (224-298) | yes |
| share of sites with martens (realised) | 0.592 | 0.582 (0.558-0.612) | yes |
| share of sites with squirrels (realised) | 0.504 | 0.513 (0.500-0.527) | yes |

Every true value lies inside its 95% interval. The broadleaf main effect on
squirrels is the weakest: the data cannot tell it from zero at this design
(planted-truth power was 0.34).

## (f) Limits that matter at this design (from VALIDATION.md)

- **260 sites is below the validated 600** for a habitat-dependent
  interaction. Here gamma0 is overstated by about 11%; intervals keep their
  coverage. Read the interval, not the point estimate.
- **Association, not causation.** Habitat that both species respond to, and
  that is not in the model, will show up as an "interaction".
- **Royle-Nichols assumptions**: abundance is Poisson across sites, martens
  are detected independently, and the population is closed for the 4 weeks.
  Martens moving in pairs or families would bias abundance and so the
  interaction. The posterior predictive check is the only guard.
- **Abundance means "martens using the camera site"**, not density. One
  marten can use several sites.
- **The two-step shortcut** (estimate martens, then plug in) is shown only
  for comparison. It fails the gate (pulled towards zero, intervals too narrow).

## (g) Reproduce everything

Run from `.`:

```bash
# 0. (optional) regenerate the example data (seed 2024)
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/abund-interactions_make_example.log -- \
  scripts/Rscript dev/abund-interactions/make_example.R

# 1. planted-truth check at the example design (100 replicates, ~20 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/abund-interactions_si_recover.log -- \
  scripts/Rscript core/run.R recover abund-interactions \
  --config modules/abund-interactions/example/recovery.yaml --nsim 100 --cores 16 \
  --out results/abund-interactions/si/recovery

# 2. fit the example data: report + figures (~8 min on 4 cores)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/abund-interactions_si_fit.log -- \
  scripts/Rscript core/run.R fit abund-interactions \
  --config modules/abund-interactions/example/config.yaml \
  --out results/abund-interactions/si/fit

# 3. (optional) redraw figures after editing the figures: block, without refitting
scripts/ejob -c 1 -m 8 -t 0:30:00 -l logs/abund-interactions_si_figures.log -- \
  scripts/Rscript core/run.R figures abund-interactions \
  --fit results/abund-interactions/si/fit/abund-interactions_fit.rds \
  --config modules/abund-interactions/example/config.yaml

# 4. machine-readable summary for the SI / website
../.venv/bin/python dev/abund-interactions/si/make_si_summary.py
```

The MCMC seed is fixed (`seed = 1`), so the fit reproduces the same numbers.

## References
- Twining JP, Augustine BC, Royle JA, Fuller AK (2025) Abundance-mediated species interactions. Ecology 106:e4468.
- Twining JP, Sutherland C, Reid N, Tosh DG (2022) Habitat mediates coevolved but not novel species interactions. Proc R Soc B 289:20212338.
- Royle JA, Nichols JD (2003) Estimating abundance from repeated presence-absence data or point counts. Ecology 84:777-790.
