---
id: occu-rota
title: Multispecies occupancy with species interactions (Rota et al. 2016)
when_to_use: Detection/non-detection surveys of 2 or 3 species at the same sites (camera traps, call surveys, eDNA, point counts reduced to 0/1), repeated so detection can be estimated; question is whether the species occur together more or less than expected under independence, how one species' occupancy differs where another is present vs absent, and whether that association changes along a habitat or disturbance gradient.
not_for: One species (use occu-single); more than 3 species or whole communities (joint species distribution models / community occupancy); asymmetric dominant-subordinate questions where one species' occupancy is assumed to depend on the other but not vice versa (occu-waddle); counts or abundance; single survey per site; multi-season co-occurrence dynamics; spatially autocorrelated occupancy; proving competition or predation (co-occurrence is not interaction).
data_needed: detection table (site x survey rows with one 0/1 column per species, or long site x survey x species), site covariate table, optional survey covariates (effort, days active, observer...)
engine: unmarked::occuMulti (maximum likelihood, default) or nimble (Bayesian, categorical latent state; settings engine = nimble)
---

# Multispecies occupancy with species interactions

## Question it answers

*Do two or three species occur at the same sites more, or less, often than we
would expect if they occurred independently, and does this change along a
habitat or disturbance gradient, after allowing for each species being missed
on some surveys?*

Outputs on the scale a biologist uses:

- **marginal occupancy** of each species (`occupancy_prob_<sp>`),
- **conditional occupancy**: the probability species A occupies a site where
  species B is present vs where B is absent (`prob_<A>_given_<B>_present`,
  `prob_<A>_given_<B>_absent`) - the classic Rota et al. (2016) comparison,
- **co-occurrence probability** (both present, `cooccurrence_prob_<A>_<B>`),
- the **interaction parameter** (`interaction_<A>_<B>`: log odds ratio of
  co-occurrence) and how it changes with covariates
  (`interaction_<A>_<B>_beta_<covariate>`),
- detection probability per survey per species and detection covariate effects,
- conditional-effects curves of all of these against each covariate.

**Co-occurrence is not interaction.** A positive or negative interaction
parameter means "the species occur together more/less than expected given the
covariates in the model". Shared or opposite habitat preferences that are not
in the model, different spatial scales of use, or one species changing the
other's detectability all produce the same signal (Blanchet et al. 2020). The
module's interpretation text always says so.

## Model

Rota et al. (2016). For S species, the latent occupancy state of site i is
the vector Z_i = (z_i1, ..., z_iS) of 0/1 values, a multivariate Bernoulli
random variable with 2^S possible states. Its probabilities are a
multinomial-logit function of *natural parameters*:

    P(Z_i = z) = exp( sum_k f_ik * prod_{s in set k} z_s ) / sum over all 2^S states of the same

with one natural parameter for every species set k up to `max_order`:

- first order f_A (one per species): log odds of A present when all other
  species are absent,
- second order f_AB (one per pair): the log odds ratio of co-occurrence;
  for two species f_12 = log( psi_11 psi_00 / (psi_10 psi_01) ),
- third order f_ABC (3 species with `max_order: 3`): change in the pairwise
  log odds ratios when the third species is present.

Each natural parameter is linear in site covariates, f_ik = x_ik' beta_k, with
its own formula; a formula of `"0"` fixes that term at 0 (species assumed
independent).

Detection, survey j: y_ijs | z_is ~ Bernoulli(z_is p_ijs),
logit(p_ijs) = w_ijs' alpha_s, one formula per species. Surveys not done
(NA) contribute nothing.

Derived quantities: marginal occupancy P(z_A = 1) = sum of state
probabilities with A present; conditional occupancy P(z_A = 1 | z_B = 1) =
P(both) / P(z_B = 1) and P(z_A = 1 | z_B = 0) = P(A, not B) / (1 - P(z_B = 1)),
summing over any third species. They are computed at **average conditions**
and along each covariate. Average conditions = the covariate value at which
each standardised model term is at its mean, so they describe the same site
as the intercepts: the arithmetic mean for a covariate used as is, the
geometric mean for one used as `log(x)` (other single transforms are solved
numerically), the reference level for factors. The report states these
values ("forest_pct = 52.4 (mean), housing_per_km2 = 39.4 (geometric mean
...)"). Curves along one covariate hold the others at these values.

Numeric covariates are standardised (mean 0, SD 1) after any transformation
(e.g. `log(housing_per_km2)`), so intercepts refer to average conditions and
slopes are **per 1 SD** of the term. The interpretation also gives effects per
a round number of raw units where possible.

Engines (identical model and design matrices):

- `unmarked` (default): `occuMulti` maximum likelihood (Fiske & Chandler 2011;
  Kellner et al. 2023). Wald 95% intervals for link-scale coefficients;
  intervals for probabilities from 2000 draws of the asymptotic multivariate
  normal of the estimates. Optional penalized likelihood (`penalty`, Clipp et
  al. 2021) for sparse data; its standard errors are bootstrapped. In our
  validation the penalty shrank a true interaction of 1.0 to 0.59 with 71%
  interval coverage (rare-species design), so use it for prediction, not
  for testing interactions.
- `nimble`: Bayesian, categorical latent-state formulation
  (z_i ~ Categorical(psi_i[1:2^S]), y_ijs ~ Bernoulli(Zmat[z_i, s] p_ijs)).
  Priors Normal(0, SD 2) on every natural-parameter and logit-detection
  coefficient (`prior_sd`). Slice samplers on coefficients (default; random
  walk mixed 5-10x worse), categorical sampler on z_i.
- `naive`: negative control only. Treats "detected at least once" as the true
  state (no detection model). Never report it.

## Assumptions (and what breaks)

1. **Closure** within the season: no species colonises or leaves a site
   between surveys. Violations inflate apparent detection heterogeneity and
   bias occupancy; choose a season short relative to turnover.
2. **No false positives** (misidentified species). False detections of a
   species where it is absent inflate its occupancy and can create spurious
   co-occurrence.
3. **Detection independent across species given occupancy.** If one species
   changes the other's detectability (e.g. a fox avoids cameras when coyotes
   pass), that appears as an occupancy interaction. Add a detection covariate
   for the other species' detection if this is plausible.
4. **Detection heterogeneity** among sites not explained by covariates biases
   occupancy low; use detection covariates (effort, trail, observer).
5. **Sites independent** and the model's covariates capture the habitat
   drivers of both species. Omitted shared habitat drivers are the main reason
   co-occurrence is not evidence of interaction.
6. **Enough data per pair**: interactions need many sites where each species
   is detected and both are detected. With a rare species, MLE interaction
   estimates are biased away from 0 (see VALIDATION.md, Known limits).

## Data format

`camera_detections.csv` (wide_species format; one row per camera x week):

| camera_id | week | start_date | days_active | coyote | gray_fox | bobcat |
|---|---|---|---|---|---|---|
| CAM-001 | 1 | 2024-06-03 | 7 | 1 | 0 | 0 |
| CAM-001 | 2 | 2024-06-10 | 5 | 1 | 0 | 0 |
| CAM-002 | 1 | 2024-06-03 | 0 |  |  |  |

- species columns: 1 = detected, 0 = not detected, blank = not surveyed;
- `surveyed_col` (optional): a week with 0 or blank effort counts as not
  surveyed, whatever the species columns say;
- missing rows (e.g. camera removed early) are treated as not surveyed;
- a survey with a blank for any species is treated as not done for all.

Long format is also accepted (`format: long`): columns site, survey, species,
detected (one row per site x survey x species).

`camera_sites.csv`: one row per site: `camera_id, forest_pct,
housing_per_km2, on_trail, lat, lon`. Character columns are factors (first
level alphabetically is the reference). Site covariates can be used in both
occupancy and detection formulas; survey covariates only in detection.

## Config

```yaml
data:
  detections: camera_detections.csv   # path relative to this config
  format: wide_species                # or long
  site_col: camera_id
  survey_col: week
  species_cols: [coyote, gray_fox]    # wide: 2 or 3 species columns (names are cleaned to letters/digits/_)
  # long format instead: species_col: species, detection_col: detected, species: [coyote, gray_fox]
  surveyed_col: days_active           # optional effort column (0/blank = survey not done)
  obs_covs: [days_active]             # survey-level covariates
  site_covs: camera_sites.csv
  max_order: 2                        # 2 = pairwise; 3 adds the three-way term (3 species)
formulas:
  occupancy: "~ forest_pct + log(housing_per_km2)"   # one string for all species, or per species:
  #   occupancy: {coyote: "~ forest_pct", gray_fox: "~ forest_pct + log(housing_per_km2)"}
  interaction: "~ log(housing_per_km2)"              # one string for all pairs, or per pair:
  #   interaction: {"coyote:gray_fox": "~ log(housing_per_km2)", "coyote:bobcat": "~ 1", "gray_fox:bobcat": "0"}
  detection: "~ days_active + on_trail"              # one string or per species
settings:
  engine: unmarked      # unmarked | nimble | naive (control only)
  standardize: true     # standardise numeric model terms
  n_draws: 2000         # unmarked: parameter draws for derived-quantity intervals
  penalty: 0            # unmarked: Clipp et al. (2021) penalty (e.g. 0.5-2) for sparse data; 0 = plain ML
  penalty_boot: 30      # bootstrap samples for SEs when penalty > 0
  gof_nboot: 100        # unmarked: parametric-bootstrap GOF replicates (0 = skip)
  gof_ndraws: 300       # nimble: posterior predictive check draws
  niter: 6000           # nimble iterations per chain
  nburnin: 1000
  nchains: 3
  thin: 1
  prior_sd: 2           # nimble: Normal(0, prior_sd) priors on all coefficients
  sampler: slice        # nimble: slice (default) | block (RW_block) | default (nimble's RW)
  seed: 1
```

Defaults: interaction formulas default to `~ 1` for pairs and `0` for the
three-way term; occupancy and detection default to `~ 1`.

## Outputs

Estimates table (`param`, `estimate`, `sd`, `lower`, `upper`; Nimble adds
`rhat`, `ess_bulk`):

| param | meaning | scale |
|---|---|---|
| `occupancy_prob_<sp>` | marginal occupancy at average covariates | probability |
| `prob_<A>_given_<B>_present` / `_absent` | occupancy of A where B is present / absent | probability |
| `cooccurrence_prob_<A>_<B>` | both present | probability |
| `interaction_<A>_<B>` | f_AB intercept: log odds ratio of co-occurrence at average covariates (3 species: where the third is absent) | logit (exp = odds ratio) |
| `interaction_<A>_<B>_beta_<cov>` | change in that log odds ratio per 1 SD of cov | logit |
| `interaction_<A>_<B>_<C>` | three-way natural parameter (max_order 3) | logit |
| `occ_<sp>_intercept` | f_sp intercept: logit P(sp present | others absent) at average covariates | logit |
| `occ_<sp>_beta_<cov>` | effect of cov on f_sp per 1 SD | logit |
| `detection_prob_<sp>` | detection per survey at average survey conditions | probability |
| `p_<sp>_beta_<cov>` | detection effect per 1 SD (factor level vs reference) | logit |

Report: methods paragraph, diagnostics (convergence, boundary estimates,
sparse-pair warnings, goodness of fit), plain-language interpretation
(including the average-conditions values used), and figures.

Publication figure suite (`figures/`, ids usable in the config `figures:`
block `include` / `exclude` / `per_figure` / `labels`):

| id | content |
|---|---|
| `conditional_occupancy` | headline: for every ordered pair, P(A present) where B is present vs absent along the covariate on the interaction (else the first numeric covariate), one panel per pair |
| `conditional_<A>_given_<B>_<cov>` | the same for every pair x numeric site covariate |
| `conditional_occupancy_average` | marginal and conditional occupancy at average conditions (dot plot) |
| `marginal_occupancy_<cov>` | marginal occupancy of every species vs each site covariate (factor: by level) |
| `interaction_<cov>` | co-occurrence odds ratio exp(f_AB) vs the covariate (log axis, null line at 1), for pairs whose interaction formula uses it |
| `co_occurrence_matrix` | P(row | column present) - P(row | column absent) at average conditions, with intervals, diverging colour centred on 0 |
| `detection`, `detection_<cov>` | detection probability per survey vs each detection covariate |
| `coefficients` | all slopes and interaction terms as odds ratios (log axis), faceted by process |
| `gof` | bootstrap / posterior predictive distribution of the fit statistic with the observed value |
| `trace`, `posterior` | Nimble engine only |

Covariates are drawn in original units; covariates entered as `log(x)` get
a log x axis. Species and covariate names can be relabelled with
`figures: labels: {gray_fox: "Gray fox", forest_pct: "Forest cover (%)"}`.
Axis and legend words are relabellable too: `occupancy_axis`,
`marginal_axis`, `conditional_axis` (`%s` = species), `interaction_axis`,
`detection_axis`, `coef_axis` (`%d` = interval level), `gof_axis`,
`matrix_x`, `matrix_y`, `matrix_fill`, `marginal`, `other_present`,
`other_absent`, `present`, `absent`, `process_occupancy`,
`process_interaction`, `process_detection`; any coefficient by its
parameter name (e.g. `occ_coyote_beta_forest_pct`). Colours follow the
entity: each species keeps its Okabe-Ito colour in every figure; the state
of the other species (present / absent) is drawn in ink / grey with solid /
dashed lines, and species pairs use Okabe-Ito 7, 5, 6 (reddish purple
first, the hue furthest from the species colours), so state and pair
curves are never mistaken for a species. `units:` back-transforms a
covariate column that was already standardised before it reached the
module (the module itself predicts on raw-unit grids).
Covariates that enter a formula as `log(x)` are drawn on a log axis with a
log-spaced prediction grid. Captions state the value at which every other
covariate is held (the mean; the geometric mean for log-transformed
covariates, which is where the log term is at its mean; factors at the
reference level). The two composites (`conditional_occupancy`, `detection`)
stack their panels in one column, with legends underneath, when their width
is under 120 mm (e.g. `per_figure: {conditional_occupancy: {width: single}}`);
a per-figure `aspect` then does not apply, because stacked panels need their
own height.
Bands are 95% confidence (unmarked; from the MVN draws) or credible
(Nimble, with an inner 50% band).

Goodness of fit (calibration in VALIDATION.md, "Independent verification"):
statistic = chi-square on the frequency of detection counts (for each species,
the number of sites with 0, 1, 2, ... detections among the surveys done,
observed vs expected under the fitted model) + chi-square on the number of
sites where each pair was detected together. Parametric bootstrap with refits
(unmarked) or posterior predictive check (Nimble). p < 0.05 or c-hat well
above 1 means the model does not reproduce the data (typically detection
varies among sites in ways the detection covariates do not capture: too many
sites with no detections and too many with many). In simulation (300 sites,
5 surveys) the bootstrap p-value was uniform under the true model (4% below
0.05) and rejected in 66% / 96% of datasets with unmodelled site-level
detection heterogeneity of SD 1 / 1.5 on the logit scale. The power of the
co-detection part against a misspecified interaction structure was not
tested, so a good fit does not validate the interaction formulas. The posterior predictive version
(Nimble) was not calibrated and, like all posterior predictive p-values, is
conservative.

## Simulation and power

`simulate()` generates S = 2 or 3 species with one standardised site
covariate `forest` (on every first-order and, if
`interaction_covariate: true`, every pairwise natural parameter) and one
survey covariate `effort` (on detection). Species are named sp1, sp2, sp3.

Params (vectors indexed by species, or by pair in the order 1:2, 1:3, 2:3;
a scalar is recycled):

| param | default |
|---|---|
| `occ_intercept` | -0.3, -0.5, -0.6 |
| `occ_beta_forest` | 0.8, -0.6, 0.4 |
| `interaction` | 1.0, -0.8, 0.5 |
| `interaction_beta_forest` | -0.7, 0.5, 0 |
| `interaction3` | 0 |
| `detection_prob` | 0.4, 0.35, 0.45 |
| `p_beta_effort` | 0.5, -0.4, 0.3 |

Design: `n_species` (2), `n_sites` (400), `n_surveys` (4), `max_order` (2),
`prop_missing` (0.1, share of surveys not done), `interaction_covariate`
(true), `species` (optional names).

Power: `power_params` are the interaction terms and slopes, occupancy and
detection slopes. Useful knobs: `--vary n_sites=...`, `--vary n_surveys=...`,
`--vary interaction=0,0.5,1,1.5 --what params`. The interaction intercept
needs roughly 400 sites at 4 surveys for 80% power at a log odds ratio of 1
(odds ratio 2.7) with both species at occupancy ~0.5; interaction slopes and
rare species need far more (VALIDATION.md).

Example recovery at a researcher's design:

```
scripts/ejob -c 16 -m 32 -t 1:00:00 -- scripts/Rscript core/run.R recover occu-rota --config my_design.yaml --nsim 200 --cores 16
```

## Pitfalls

- **Rare species / few co-detections**: interaction estimates are biased away
  from 0 (MLE +24% at marginal occupancy 0.12) and have huge intervals; the
  fit may hit a boundary (estimate > 8 or SE > 10 is flagged as not
  converged). Use simpler interaction formulas (`~ 1`), fix interactions you
  cannot estimate at `"0"`, and report the conditional probabilities rather
  than the raw interaction coefficient (they are recovered well even when the
  coefficient is not). `penalty` (Clipp et al. 2021) stabilises predictions
  but its intervals under-cover the interaction (VALIDATION.md). The Nimble
  engine (N(0, 2) priors) recovered the rare-species interaction with +5%
  bias and nominal coverage at 400 sites, so prefer `engine: nimble` here.
- **Interaction slopes are data-hungry**: at 400 sites they had power 0.56-0.65
  for a slope of -0.7 and MLE slopes are inflated by ~10% (3 species: 13%,
  failing the bias gate; VALIDATION.md).
- **Do not put many covariates on interactions**; each pair term costs as
  much information as a whole species model.
- **Three-way terms** (`max_order: 3`) are rarely estimable; leave at 0 unless
  you have many hundreds of sites.
- **Detection interactions**: if one species changes the other's detection,
  the model attributes it to occupancy.
- **Species naming**: the order of `species_cols` determines parameter names
  (`interaction_coyote_gray_fox`); conditional probabilities are given in both
  directions.
- **unmarked parameter draws** are an asymptotic approximation; with sparse
  data prefer the Nimble engine and compare.

## Validation summary

Full numbers in VALIDATION.md. Default design (2 species, 400 sites, 4
surveys): every estimand passes the recovery gate with both engines
(unmarked nsim 200, Nimble nsim 100; coverage 0.91-0.99). False-positive
rate for an interaction of 0: 0.036 (intercept) and 0.052 (slope) at nsim
500. Ignoring detection fails 11 of 13 estimands. Derived probabilities match
`unmarked::predict(species, cond)` to < 1e-14; the hand-written likelihood
matches occuMulti's to 4 decimals; unmarked and Nimble agree per dataset
(r >= 0.994). Limits: interaction slopes inflated ~10-13% at 400 sites
(3 species: fails the bias gate; fine at 800), rare species (occupancy
~0.12) ML interaction biased +24% at 400 sites (fine at 800; Nimble +5%,
passes), three-way terms untested. An independent verification with a
separately written simulator (different species, covariates incl. a log
transform and a factor, 300 and 500 sites, 3 species with a pair fixed at 0)
passed every estimand with both engines; see VALIDATION.md.

## References

- Rota CT, Ferreira MAR, Kays RW, Forrester TD, Kalies EL, McShea WJ, Parsons AW, Millspaugh JJ (2016) A multispecies occupancy model for two or more interacting species. Methods in Ecology and Evolution 7:1164-1173. doi:10.1111/2041-210X.12587
- Fiske I, Chandler R (2011) unmarked: An R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23. doi:10.18637/jss.v043.i10
- Kellner KF, Smith AD, Royle JA, Kery M, Belant JL, Chandler RB (2023) The unmarked R package: Twelve years of advances in occurrence and abundance modelling in ecology. Methods in Ecology and Evolution 14:1408-1415. doi:10.1111/2041-210X.14123
- Clipp HL, Evans AL, Kessinger BE, Kellner K, Rota CT (2021) A penalized likelihood for multispecies occupancy models improves predictions of species interactions. Ecology 102:e03520. doi:10.1002/ecy.3520
- MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA, Langtimm CA (2002) Estimating site occupancy rates when detection probabilities are less than one. Ecology 83:2248-2255.
- Blanchet FG, Cazelles K, Gravel D (2020) Co-occurrence is not evidence of ecological interactions. Ecology Letters 23:1050-1063. doi:10.1111/ele.13525
- de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.
- unmarked package vignette "Multispecies occupancy models with occuMulti", CRAN: https://cran.r-project.org/web/packages/unmarked/vignettes/occuMulti.html (natural-parameter definitions, formula ordering, `predict(species, cond)`, penalty interface).
