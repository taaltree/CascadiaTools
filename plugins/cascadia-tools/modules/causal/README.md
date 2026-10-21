---
id: causal
title: Causal effects from observational data (DAG + g-computation)
when_to_use: The question is causal ("does canopy cover increase trout?", "what would happen to native plants if pastures were not grazed?") and the data are observational (field surveys, monitoring). The researcher can draw a DAG of how their variables cause each other. Also panels (sites measured repeatedly) with unmeasured site-level or year-level confounders, via fixed effects.
not_for: Pure prediction or description (use glm); experiments with randomised treatment (use glm; no adjustment needed); hypotheses about a whole network of paths with latent constructs (use sem); imperfect detection of the outcome (occupancy/N-mixture modules; their covariates can be chosen with this module's DAG); feedback over time (past outcome -> future exposure), mediation with exposure-mediator interaction, or effects identified only by the front-door criterion or instruments (not implemented).
data_needed: one CSV, one row per site (or per site x time for panels); a column for the exposure, the outcome and each measured variable in the DAG; a DAG written in dagitty syntax in the config
engine: dagitty (identification, local tests) + stats::glm / MASS::glm.nb (outcome model) + marginaleffects (g-computation) + sandwich (cluster-robust SE)
---

# causal: causal effects from observational data with DAGs

## Question it answers

*What would happen to the outcome if we changed the exposure?* For example:
"If riparian canopy cover increased by 10 percentage points at every reach,
how many more juvenile trout would there be per 100 m?" The answer is an
**average causal effect** on the scale the biologist measures (count,
probability, mean), plus a ratio (rate or risk ratio) and the link-scale
coefficient.

A regression coefficient is only a causal effect when the model adjusts for
the right variables. Adding variables "to be safe" can *create* bias:
adjusting for a **mediator** removes part of the effect you want, and
adjusting for a **collider** (a variable caused by both exposure and
outcome) or an outcome descendant manufactures a spurious association. Leaving out a **confounder** biases
the estimate the other way. Model-selection tools (AIC, stepwise) cannot
tell these apart (Arif & MacNeil 2023; Laubach et al. 2021). The DAG can.

## Model

1. **DAG**: you write your causal assumptions as a directed acyclic graph,
   with one arrow per direct cause. Nodes that are not columns of the data are
   unmeasured (latent).
2. **Identification** (dagitty; Textor et al. 2016): the module lists every
   *minimal adjustment set*, the smallest sets of measured variables that
   block every back-door (non-causal) path between exposure X and outcome Y
   without opening new ones (Pearl 2009, adjustment criterion). By default
   it uses the smallest one; `adjust: optimal` uses the optimal set of
   Henckel et al. (2022), which gives the smallest variance among valid sets
   when it is fully measured. If no set of measured variables works, the
   module stops and says which unmeasured variable is the problem.
3. **Outcome model**: a GLM `g(E[Y]) = a + b X + c'Z`, where Z is the adjustment
   set. Families: gaussian (identity), poisson (log), negbin (log; MASS),
   binomial (logit; 0/1 or successes/trials), Gamma (log). Optional effort
   offset.
4. **Average causal effect by g-computation** (standardisation;
   marginaleffects, Arel-Bundock et al. 2024): predict every site's outcome
   with its exposure set to the contrast (continuous: X + delta vs X;
   binary: 1 vs 0), keep each site's own Z, and average.
   `ate = mean(mu_hi) - mean(mu_lo)`; `ratio = mean(mu_hi) / mean(mu_lo)`.
   Standard errors: delta method plus an influence-function term for the
   sampling of sites (`ci_method: unconditional`, default), so intervals
   refer to the population of sites. `delta` (conditional on these sites) and
   `bootstrap` (percentile, whole units for panels) are also available.
5. **Controlled direct effect** (`effect: direct` or `both`): a second GLM
   with the set dagitty gives for the direct effect, which holds the
   mediators fixed.
6. **Panels with unmeasured confounders** (Byrnes & Dee 2025): if an
   unmeasured variable is constant within each site (`unit_constant`), site
   fixed effects (`panel: fixed`) remove it. If it is shared by all sites in
   a year (`time_constant`), year fixed effects remove it (`panel: twoway`
   for both). The DAG treats these variables as adjusted for. Standard
   errors are cluster-robust by site (sandwich::vcovCL, HC1).
7. **Checks**: the DAG's testable implications. One conditional-independence
   test per non-adjacent pair of measured variables (dagitty::localTests,
   partial correlations; `cis.pillai` if factors are present), combined with
   Fisher's C (Shipley 2000), with Holm-adjusted individual tests. The
   outcome model is checked with DHARMa simulated residuals. An E-value
   (VanderWeele & Ding 2017) for unmeasured confounding, and positivity/extrapolation checks.
8. **What-if table** (educational negative controls): the same effect
   re-estimated with no adjustment, omitting each confounder, adding each
   mediator, collider or outcome descendant, other valid sets, and "all
   measured variables" (Cinelli et al. 2022). Each row is labelled valid or
   biased according to the DAG. Valid rows should agree with the main
   estimate; biased rows show how much the answer would move.

## Assumptions (and what breaks if they fail)

- **The DAG is right and complete**: no unmeasured common cause of X and Y
  beyond those drawn, and every missing arrow is a real "no direct effect".
  If this fails, the estimate is biased by an unknown amount. Local tests
  can reveal some missing arrows, but they *cannot* detect an unmeasured
  confounder of X and Y. The E-value tells you how strong one would need to be.
- **Positivity**: every site could plausibly have had each exposure value.
  The module flags extrapolation (continuous X) and extreme propensities
  (binary X).
- **Consistency / no interference**: the contrast is a well-defined
  intervention, and one site's exposure does not change another site's
  outcome.
- **Model form**: the effect is linear on the link scale and does not
  interact with the adjustment variables. For logit links the conditional
  odds ratio is non-collapsible, so causal statements use the probability
  difference.
- **Measurement**: confounders are measured without serious error. Error
  leaves residual confounding.
- **Fixed effects**: unmeasured confounders are constant within a site
  (one-way FE), or shared across sites within a year (two-way FE). Past
  outcomes do not change future exposure, and effects are homogeneous.
  Time-varying site-specific confounders are not removed.

## Data format

One CSV. Column names are the DAG node names (or map them with `columns:`).

```
reach_id,stream,survey_date,canopy_pct,trout_per_reach,elevation_m,agriculture_pct,fine_sediment_pct,wetted_width_m,water_temp_c,periphyton_chla
R001,Mill Cr,2025-08-01,66,15,780,43.8,14.1,4.7,11.2,57.8
R002,Beaver Cr,2025-07-04,43,3,678,39.2,25.7,8,17.7,63
```

Panels: one row per site x time, plus the site and time columns
(`example/snail_panel.csv`: `site,year,sst_c,upwelling_index,snail_count`).
Text columns are treated as categorical adjusters. A two-level text exposure
is coded 0/1 (`exposure_reference:` sets the 0 level). Rows with missing
values in DAG variables are dropped, with a note.

## Config (`data:` block)

| key | meaning |
|---|---|
| `file` | CSV path (relative to the config) |
| `dag` | DAG in dagitty syntax (`dag { a -> b; c [latent] ... }`); or `edges:` list of `"a -> b"` |
| `exposure`, `outcome` | node names |
| `family`, `link` | gaussian, poisson, negbin, binomial, Gamma; link optional |
| `trials` | binomial: column of trials (default 1 = 0/1 outcome) |
| `offset`, `offset_per` | effort column (log offset); effects reported per `offset_per` units of effort (default 1) |
| `contrast` | continuous exposure: the change to evaluate (default 1 SD, rounded); report it in biological units (10 % canopy, 1 deg C) |
| `exposure_reference` | 2-level text exposure: the reference level |
| `effect` | `total` (default), `direct`, `both` |
| `adjust` | `minimal` (default), `optimal`, or an explicit list. An explicit set that the DAG says is invalid is rejected by check_data |
| `panel` | `none` (default), `fixed` (unit FE), `twoway` (unit + time FE) |
| `unit`, `time` | panel id columns |
| `unit_constant`, `time_constant` | latent DAG nodes absorbed by the unit / time fixed effects |
| `columns` | map `node: column` if names differ |
| `labels`, `unit_text`, `response_label` | text for reports and figures |

`settings:` `ci_method` (unconditional | delta | bootstrap), `n_boot`,
`vcov` (auto = cluster-robust by unit for panels, sandwich HC3 for poisson,
model-based otherwise | model | robust), `whatif`, `local_tests`, `ci_test_type`,
`diagnostics`. `force_adjust` and `force_panel` exist only for negative
controls.

## Outputs

Estimates table (`*_estimates.csv`, report). Each row has `adjusted_for`,
`role`, `valid_by_dag` and `strategy`.

| param | meaning |
|---|---|
| `ate_<x>` | average causal effect of the contrast on the response scale (count, probability, mean) |
| `rate_ratio_<x>` / `risk_ratio_<x>` / `mean_ratio_<x>` | ratio of the two average potential outcomes (non-identity links) |
| `beta_<x>` | link-scale coefficient per 1 unit of x (logit: conditional log odds ratio, non-collapsible) |
| `mean_<y>_if_<x>_0/1` | potential-outcome means (binary exposure) |
| `direct_*` | the same for the controlled direct effect |
| `whatif_<strategy>__ate_<x>` | the effect under an alternative adjustment strategy (what-if table) |

## Figures

`core/run.R fit` exports these to `<out>/figures/` (PDF + PNG, with a
`gallery.html`). Redraw without refitting with
`scripts/Rscript core/run.R figures causal --fit <out>/causal_fit.rds --config config.yaml`.

| id | what it shows |
|---|---|
| `dag` | The DAG as drawn in config. Node fill shows the role: exposure (orange), outcome (sky blue), adjustment set used (bluish green), absorbed by fixed effects (pink), unmeasured (grey), not adjusted (white). Exposure and outcome are bold and unmeasured nodes are italic, so roles survive greyscale. Solid thick arrows lie on causal paths from exposure to outcome; dashed arrows are the other assumed effects. Double width. |
| `effect` | Causal dose-response: the average predicted outcome (response scale: count, probability or mean) if every site were set to each exposure value, by g-computation over the observed sites, with an `ci` confidence band and a rug of observed exposure values. The x axis is in the exposure's own units over its 2nd-98th percentiles; the y axis includes 0 for counts and positive responses. For a binary exposure it shows the two potential-outcome means with intervals. |
| `adjustment_comparison` | The what-if forest plot: the same effect under the DAG's set (green circle and vertical line), other valid sets (blue triangles) and biased choices (orange squares: crude, omitting a confounder, adding a mediator or collider, all measured variables, pooled instead of fixed effects). Intervals are redrawn at `ci` (Wald, from the stored SE). Variable names are shown with their labels. |
| `implied_independencies` | The DAG's testable implications (`dagitty::localTests`): each partial correlation the DAG says is zero, with its interval (Fisher z, redrawn at `ci`), coloured and shaped by whether it is contradicted after Holm correction. Only drawn when local tests ran on numeric data. |
| `residuals` | DHARMa simulated quantile residuals of the outcome model against uniform quantiles (1:1 line). |

All figures are frequentist (confidence intervals), so `inner_ci` does not
apply. Rename anything with `figures: labels:`. Use DAG node names
(`canopy_pct: "Riparian canopy cover (%)"`; when a node is relabelled here, the
`unit_text` suffix is dropped from the effect axis, so put units in the label).
Legend entries: `role_exposure`, `role_outcome`, `role_adjusted`,
`role_absorbed`, `role_unmeasured`, `role_not_adjusted`, `path_causal`,
`path_other`, `whatif_main`, `whatif_valid`, `whatif_biased`,
`test_consistent`, `test_contradicted`. Axes: `effect_x`, `effect_y`,
`comparison_x`, `tests_x`, `tests_y`, `resid_x`, `resid_y`. For a binary
exposure, use `<exposure>_<level>`. Legends go below the plot by default; set
`legend:` globally or in `per_figure:` to move them. `greyscale: true`
switches the role and status colours to greys, and the second encodings
(bold/italic nodes, point shapes) keep them distinguishable. Covariate units
come from the data. A `units:` entry for the exposure back-transforms its axis
if you supplied it standardised.

## Simulation and power

`simulate()` is a structural causal model built from the DAG. Every
non-outcome node is linear-Gaussian and standardised, and the
`b_<from>__<to>` params are standardised path coefficients. The exposure can be
`exposure_type: binary` (logistic in its parents). The outcome uses the
family's noise (`intercept`, `sigma`, `theta`, `shape`, design `trials`).
Panels: `n_units`, `n_times`, `unit_level`/`time_level` nodes drawn once per
unit/time. Truths are computed by intervention on the SCM (Monte Carlo,
4e5 draws): the total effect (sum of path products) on the link scale,
`ate`, and ratios. `whatif_truth: true` adds truths for the what-if rows.
Knobs for power: `n_sites` (or `n_units`, `n_times`), any path coefficient
(e.g. `b_canopy__trout`), `contrast`. At the default (standardised effects
0.2 direct + 0.16 via temperature, NB theta 3), 80% power for the
rate ratio needs about 60 sites (VALIDATION.md).

## Pitfalls

- Drawing the DAG after looking at the data, or pruning arrows until local
  tests pass. Write the DAG from biology first; report any revisions.
- "Controlling for everything": the what-if table shows why that is wrong.
- Reading the coefficients of the *adjustment variables* as causal effects
  (the "Table 2 fallacy"). Each needs its own DAG analysis with itself as the exposure.
- Logit coefficients change when variables are added even without
  confounding (non-collapsibility). Use `ate`/`risk_ratio`.
- With fewer than ~60 sites, response-scale effects on counts are skewed:
  `ate` was biased +18% at 40 sites (VALIDATION.md); the ratio and `beta` were fine.
- A total-effect model leaves the mediators and other causes of the
  outcome out of the regression, so counts are overdispersed *given the
  adjustment set* even when the process is Poisson. Prefer `negbin`. With
  `poisson`, the module uses sandwich SEs by default, because model-based
  SEs under-cover (VALIDATION.md).
- Fixed effects only remove *time-invariant* site confounding, and they discard
  between-site information (wider intervals).

## References

- Arif S, MacNeil MA (2023) Applying the structural causal model framework for observational causal inference in ecology. Ecological Monographs 93(1): e1554. doi:10.1002/ecm.1554
- Laubach ZM, Murray EJ, Hoke KL, Safran RJ, Perng W (2021) A biologist's guide to model selection and causal inference. Proc R Soc B 288: 20202815. doi:10.1098/rspb.2020.2815
- Byrnes JEK, Dee LE (2025) Causal inference with observational data and unobserved confounding variables. Ecology Letters 28(1): e70023. doi:10.1111/ele.70023
- Pearl J (2009) Causality, 2nd edn. Cambridge University Press.
- Textor J et al. (2016) Robust causal inference using directed acyclic graphs: the R package 'dagitty'. Int J Epidemiol 45: 1887-1894. doi:10.1093/ije/dyw341
- Cinelli C, Forney A, Pearl J (2022) A crash course in good and bad controls. Sociological Methods & Research 53(3): 1071-1104. doi:10.1177/00491241221099552
- Henckel L, Perkovic E, Maathuis MH (2022) Graphical criteria for efficient total effect estimation via adjustment in causal linear models. JRSS B 84: 579-599. doi:10.1111/rssb.12451
- Arel-Bundock V, Greifer N, Heiss A (2024) How to interpret statistical models using marginaleffects for R and Python. J Stat Softw 111(9). doi:10.18637/jss.v111.i09
- Shipley B (2000) A new inferential test for path models based on directed acyclic graphs. Structural Equation Modeling 7: 206-218.
- VanderWeele TJ, Ding P (2017) Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med 167: 268-274. doi:10.7326/M16-2607
- Mundlak Y (1978) On the pooling of time series and cross section data. Econometrica 46: 69-85.
