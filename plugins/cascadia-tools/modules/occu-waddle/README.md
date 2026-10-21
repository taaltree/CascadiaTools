---
id: occu-waddle
title: Two-species conditional occupancy (dominant/subordinate; Richmond et al. 2010, Waddle et al. 2010)
when_to_use: Repeated detection/non-detection surveys of TWO species at the same sites, where one species can reasonably be treated as dominant (predator, invasive competitor, larger species) and the question is whether the other (subordinate) species occupies sites less (or more) often where the dominant is present, and whether detection of one depends on the other; covariates on every occupancy and detection probability.
not_for: Three or more species, or no defensible dominant/subordinate ordering (use occu-rota, which is symmetric); counts or abundance of the dominant (abund-interactions); single species (occu-single); multi-season colonisation/extinction dynamics.
data_needed: detection history of both species (site x survey 0/1, both recorded on every survey), site covariates, optional survey covariates
engine: nimble (Bayesian, marginal likelihood; default) or hand-coded likelihood + optim (maximum likelihood)
---

# occu-waddle: two-species conditional occupancy

## Question it answers

"Is the subordinate species (B) found at fewer sites where the dominant
species (A) occurs than where it is absent, and is that pattern stronger in
some habitats?", "Is B harder to detect on surveys when A is detected?",
after allowing for both species being missed on some surveys.

Typical uses: native vs invasive frogs (Waddle et al. 2010: green and
squirrel treefrogs vs the invasive Cuban treefrog), small vs large rails
(Richmond et al. 2010: California black rail vs Virginia rail), prey vs
predator at camera stations, native vs introduced fish in ponds.

## Model

Sites i = 1..M, surveys j = 1..J. Species A is dominant, B subordinate.

State (latent, summed out of the likelihood):

- z_A,i ~ Bernoulli(psiA_i)
- z_B,i | z_A,i = 1 ~ Bernoulli(psiBA_i); z_B,i | z_A,i = 0 ~ Bernoulli(psiBa_i)

Detection on survey j:

| sites where | detection of A | detection of B |
|---|---|---|
| only A occurs | pA_ij | - |
| only B occurs | - | pB_ij |
| both occur | rA_ij | rBA_ij if A was detected on that survey, rBa_ij if not |

Each of the eight probabilities is logit-linear in covariates:
logit(theta) = b0 + b1 x1 + ... Occupancy probabilities take site
covariates; detection probabilities take site or survey covariates.
Numeric covariates (columns of the model matrix, so log(x) or x^2 terms
too) are standardised to mean 0, SD 1, so the intercept is the probability
at average conditions and slopes are per 1 SD.

Site likelihood (Richmond et al. 2010; the same likelihood as wiqid::occ2sps):

L_i = psiA psiBA prod_j [rA^yA (1-rA)^(1-yA) q^yB (1-q)^(1-yB)], q = rBA if yA = 1 else rBa
    + psiA (1 - psiBA) prod_j pA^yA (1-pA)^(1-yA) * I(B never detected)
    + (1 - psiA) psiBa prod_j pB^yB (1-pB)^(1-yB) * I(A never detected)
    + (1 - psiA)(1 - psiBa) * I(neither detected)

Derived quantities (at average covariate values):

- occ_subordinate: psiB = psiA psiBA + (1 - psiA) psiBa
- occ_both: psiA psiBA
- species_interaction_factor (SIF; MacKenzie et al. 2004): psiA psiBA / (psiA psiB) = psiBA / psiB.
  1 = independent occurrence, < 1 = the species co-occur less often than expected by chance, > 1 more often.
- occ_log_odds_ratio_sub: logit(psiBA) - logit(psiBa), the change in the log-odds of B occupancy when A is present (0 = no effect).
- det_log_odds_ratio_sub: logit(rBA) - logit(rBa), the detection interaction (0 = detecting A on a survey does not change the chance of detecting B).

Unlike the SIF parameterisation of MacKenzie et al. (2004), this
conditional parameterisation stays numerically stable when covariates are
added (Richmond et al. 2010; Waddle et al. 2010), but it is asymmetric:
the dominant species must be chosen in advance.

**Constraints (ties).** A parameter with no formula is constrained equal to
another (the wiqid / PRESENCE convention): psiBA -> psiBa (no occupancy
interaction), rA -> pA, rBa -> pB, rBA -> rBa. psiA, psiBa, pA and pB
default to "~ 1". Omitting rA, rBA and rBa gives the model with no
detection interaction; omitting psiBA gives the null model of independent
occupancy. Compare such models with AIC from the `mle` engine if needed.

**Priors (nimble engine).** Intercepts ~ Logistic(0, 1) (uniform on the
probability scale at average covariates; Northrup & Gerber 2018); slopes ~
Normal(0, SD = 1) per 1 SD of the covariate (`prior_sd`). This is weakly
informative: 95% prior mass on odds ratios 0.14-7.1 per SD of a covariate,
in the spirit of Gelman et al. (2008). It is deliberately not vaguer: with
SD = 2 the posterior median of the slope of subordinate occupancy where the
dominant is absent was exaggerated by 27% at 200 sites (recovery gate FAIL,
VALIDATION.md section 1b); with SD = 1 every estimand passes. The cost is
shrinkage of large effects toward 0: a psiBa slope of -1.5 per SD was
recovered as -1.33 (-11%) at 200 x 4 (VALIDATION.md section 8). Set `prior_sd: 2`
as a sensitivity check (eco-convergence / eco-model-checking).

## Assumptions (and what breaks)

- **Closure** within the season for both species. Movement in/out between
  surveys makes "occupancy" mean "use" and inflates psi.
- **Dominance is one-directional.** B's occupancy may depend on A, not the
  reverse. This is an assumption. With intercept-only occupancy the occupancy
  part describes the same joint distribution either way round, but the
  detection part (rBA / rBa) is not symmetric, and with covariates the two
  orderings are different models that can give different answers.
- **No false positives** (misidentification). A few false detections of A
  at B-only sites are read as co-occurrence.
- **Independence among sites** and among surveys given the states and
  covariates. Unmodelled detection heterogeneity (e.g. observer, weather)
  can masquerade as a detection interaction if it affects both species on
  the same survey: rBA != rBa can then mean "good survey days", not
  behaviour. Put shared survey conditions in the detection formulas.
- **Co-occurrence is not interaction.** psiBA < psiBa can arise from shared
  habitat preferences that are not in the model (Blanchet et al. 2020).
  Include the habitat covariates that plausibly drive both species.

## Data format

One CSV with one row per site x survey (the field sheet), and one CSV
with one row per site.

pipe_surveys.csv:

```
site_id,visit,date,observer,air_temp_c,pipes_checked,cuban_treefrog,green_treefrog
WCA-001,1,2023-06-05,JHW,22.9,10,1,0
WCA-001,2,2023-06-29,MS,29.6,10,1,0
WCA-001,3,2023-07-20,JHW,21.3,10,0,1
WCA-002,1,2023-06-06,KR,,0,,
```

sites.csv:

```
site_id,dist_road_m,canopy_pct,lat,lon
WCA-001,370,60,25.72798,-80.4748
```

- detections are 0/1 (yes/no and true/false are also read; convert counts to presence); both species must be
  recorded (0 or 1) on every survey that was done; a survey that was not
  done has an empty cell or 0 in `surveyed_col`.
- Sites can have different numbers of surveys (unequal effort is handled
  as missing surveys).
- A "long" file (one row per site x survey x species) is also accepted.

## Config

```yaml
data:
  detections: pipe_surveys.csv
  format: wide                # wide (one column per species) | long
  site_col: site_id
  survey_col: visit
  dominant_col: cuban_treefrog      # wide: detection column of species A
  subordinate_col: green_treefrog   # wide: detection column of species B
  # long format instead: species_col, detection_col, dominant: <name>, subordinate: <name>
  dominant: Cuban treefrog          # display names (optional for wide)
  subordinate: green treefrog
  surveyed_col: pipes_checked       # optional; 0/NA = survey not done
  obs_covs: [air_temp_c]            # survey-level covariate columns
  site_covs: sites.csv
formulas:                           # omit a key to tie it (see Constraints)
  psiA:  "~ log(dist_road_m)"
  psiBA: "~ canopy_pct"
  psiBa: "~ canopy_pct"
  pA:    "~ air_temp_c"
  pB:    "~ air_temp_c"
  rA:    "~ 1"
  rBA:   "~ 1"
  rBa:   "~ 1"
settings:
  engine: nimble      # nimble (default) | mle | naive (negative control only)
  niter: 6000         # nimble: iterations per chain
  nburnin: 1000
  nchains: 3
  prior_sd: 1         # SD of Normal prior on standardised slopes
  sampler: slice      # slice (per coefficient) | af_slice (joint) | default
  mle_starts: 3       # mle: random restarts
  n_draws: 2000       # mle: draws from the asymptotic normal for derived quantities
  gof_ndraws: 500     # predictive-check draws (0 = skip)
figures:              # see docs/FIGURES.md; e.g. labels: {canopy_pct: "Canopy cover (%)"}
```

## Outputs

Estimates table (`occu-waddle_estimates.csv`):

| param | meaning |
|---|---|
| occ_dominant | psiA at average covariates |
| occ_sub_given_dom_present / _absent | psiBA / psiBa at average covariates |
| det_dominant_alone / det_sub_alone | pA / pB |
| det_dominant_with_sub | rA |
| det_sub_with_dom_detected / _missed | rBA / rBa |
| `<above>_beta_<term>` | logit-scale slope per 1 SD of the term (or level contrast) |
| occ_subordinate, occ_both | psiB, P(both) |
| species_interaction_factor | SIF (null = 1) |
| occ_log_odds_ratio_sub | logit(psiBA) - logit(psiBa) (null = 0) |
| det_log_odds_ratio_sub | logit(rBA) - logit(rBa) (null = 0) |

Tied parameters are not reported separately; derived quantities fixed by a
tie are dropped.

Diagnostics: R-hat, bulk ESS (nimble), notes on sparse data (e.g. few
sites with both species detected), and two predictive checks: number of
sites in each co-detection class (neither / A only / B only / both), and
per-site detection counts. p-values near 0 or 1 indicate lack of fit.

## Figures

`core/run.R fit` exports these to `<out>/figures/` (PDF + PNG, plus
`gallery.html` and `figures_index.csv` with captions). Redraw without
refitting: `scripts/Rscript core/run.R figures occu-waddle --fit <out>/occu-waddle_fit.rds --config config.yaml`.

| id | what it shows |
|---|---|
| `conditional_occupancy` | (A) psiA, psiB, psiBa and psiBA at average covariates (posterior median, 95% and 50% intervals); (B) psiBA vs psiBa along the first psiBA/psiBa covariate, with bands and a rug of sites. Only panel A when psiBA is tied to psiBa. |
| `conditional_occupancy_<cov>` | psiBA vs psiBa along each further psiBA/psiBa covariate |
| `sif` | species interaction factor with the null line at 1: (A) at average covariates, (B) along the first numeric occupancy covariate. Omitted when psiBA is tied to psiBa. |
| `detection_interaction` | pA, rA, pB, rBa, rBA at average survey conditions; rBA vs rBa is the detection interaction |
| `coefficients` | every slope as an odds ratio per 1 SD (or vs the reference level), coloured by species; the caption gives each SD in original units |
| `occupancy_<cov>` | psiA, psiB (overall), psiBA, psiBa along each occupancy covariate |
| `detection_<cov>` | detection parameters whose formula includes the covariate |
| `gof` | predictive check of the co-detection classes (only when it ran) |
| `trace` | MCMC traces of all coefficients (nimble) |
| `posterior` | posterior densities of psiBA, psiBa, SIF, the occupancy and detection log-odds ratios, rBA and rBa, with null lines (nimble; truth lines on simulated data) |

Conventions: x axes are in the covariate's original units (covariates
entered as `log(x)` are drawn on a log axis labelled in x). Colour follows the
quantity in every figure: dominant species = orange, subordinate overall =
sky blue, subordinate with the dominant absent / missed = bluish green, present /
detected = vermillion. Linetype and shape repeat this for greyscale. You can
rename any label with `figures: labels:`, using a covariate name, `dominant`,
`subordinate`, a quantity (`psiA psiB psiBA psiBa sif pA pB rA rBA rBa`), an
estimates-table name (for coefficient rows), or an axis key (`occupancy_axis`,
`sub_occupancy_axis`, `detection_axis`, `sif_axis`, `sif_short_axis`,
`odds_ratio_axis`, `gof_axis`). The naive engine produces no figures.

## Simulation and power

`params` (probabilities at average covariates and logit slopes):
`occ_dominant`, `occ_dominant_beta_habitat`, `occ_sub_given_dom_present`,
`occ_sub_given_dom_present_beta_habitat`, `occ_sub_given_dom_absent`,
`occ_sub_given_dom_absent_beta_habitat`, `det_dominant_alone`,
`det_dominant_alone_beta_effort`, `det_sub_alone`,
`det_sub_alone_beta_effort`, `det_dominant_with_sub`,
`det_sub_with_dom_detected`, `det_sub_with_dom_missed`.

`design`: `n_sites`, `n_surveys`, `prop_missing`, `det_interaction`
(FALSE: rA, rBA, rBa equal pA, pB and are not estimated),
`occ_interaction` (FALSE: psiBA = psiBa and is not estimated).

Power parameters: `occ_log_odds_ratio_sub`, `species_interaction_factor`
(null 1), `det_log_odds_ratio_sub`, and every slope.

Validated (VALIDATION.md): at 200 sites x 4 surveys all 18 estimands pass
the recovery gate (nsim 100) and the false-positive rate for an interaction
is 0.045-0.055 (nsim 200). At 60 sites psiBA, P(both) and SIF are biased
toward "no interaction" (+12% to +21%; gate FAIL) and power for SIF = 0.6 is
about 0.5. Power for that occupancy interaction is 0.86 at 150 sites and
0.98 at 200; a detection interaction of -0.8 on the log-odds scale reaches
power 0.7 only at 300 sites. Waddle et al. (2010) likewise found occupancy
estimates unbiased with 50-100 sites x 3 surveys only when detection >= 0.5,
and > 200 sites needed for harder-to-detect species.

## Pitfalls

- **Choose the dominant species before fitting**, on biological grounds.
  Do not try both orders and report the one with the more striking result.
- **Few sites with both species detected** (< 10) leave psiBA, rA, rBA and
  rBa weakly identified: intervals are wide and the mle engine may hit a
  boundary (reported in Diagnostics).
- A difference between rBA and rBa is a *detection* effect (calling
  suppression, hiding in refugia, shared survey conditions), not a habitat
  effect.
- SIF depends on psiA: the same psiBA/psiBa contrast gives a SIF closer to
  1 where A is rare. Report the log-odds ratio alongside SIF.
- Real-data standard: R-hat < 1.01 and bulk ESS >= 400 for every
  parameter; raise `niter` if the report says otherwise.

## References

- Richmond OMW, Hines JE, Beissinger SR (2010) Two-species occupancy models: a new parameterization applied to co-occurrence of secretive rails. Ecological Applications 20:2036-2046. doi:10.1890/09-0470.1
- Waddle JH, Dorazio RM, Walls SC, Rice KG, Beauchamp J, Schuman MJ, Mazzotti FJ (2010) A new parameterization for estimating co-occurrence of interacting species. Ecological Applications 20:1467-1475. doi:10.1890/09-0850.1
- MacKenzie DI, Bailey LL, Nichols JD (2004) Investigating species co-occurrence patterns when species are detected imperfectly. Journal of Animal Ecology 73:546-555. doi:10.1111/j.0021-8790.2004.00828.x
- MacKenzie DI, Lombardi JV, Tewes ME (2021) A note on investigating co-occurrence patterns and dynamics for many species, with imperfect detection and a log-linear modeling parameterization. Ecology and Evolution 11:8507-8515. doi:10.1002/ece3.7604
- Meredith M (2022) wiqid: Quick and dirty estimates for wildlife populations. R package 0.3.3, function occ2sps (cross-check).
- Rota CT et al. (2016) A multispecies occupancy model for two or more interacting species. Methods in Ecology and Evolution 7:1164-1173 (unmarked::occuMulti; cross-check).
- Blanchet FG, Cazelles K, Gravel D (2020) Co-occurrence is not evidence of ecological interactions. Ecology Letters 23:1050-1063. doi:10.1111/ele.13525
- Northrup JM, Gerber BD (2018) A comment on priors for Bayesian occupancy models. PLoS ONE 13:e0192819. doi:10.1371/journal.pone.0192819
- Gelman A, Jakulin A, Pittau MG, Su Y-S (2008) A weakly informative default prior distribution for logistic and other regression models. Annals of Applied Statistics 2:1360-1383. doi:10.1214/08-AOAS191
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.
