---
id: sem
title: Structural equation models (piecewise SEM and lavaan)
when_to_use: Several measured variables on the same sampling units (plots, sites, lakes, animals) are thought to influence each other as a network (e.g. fire -> canopy -> shrubs -> birds), and the question is which links the data support, how strong each direct path is, and how much of a driver's association with an outcome runs indirectly through mediators. Piecewise engine for count / binary / Gaussian responses and random site intercepts; lavaan engine for latent variables measured by several indicators, multigroup comparisons and correlated errors.
not_for: A single response with several predictors and no mediators (use glm); estimating one causal effect of an exposure with confounder adjustment (use causal); imperfect detection (occupancy / N-mixture modules; an SEM on raw counts treats detection as part of abundance); feedback loops (A -> B -> A) or time-series dynamics; interactions or nonlinear (spline) paths; very small samples (< ~30 units, or fewer than ~10 units per estimated coefficient in a component).
data_needed: one CSV, one row per sampling unit, with a numeric column for every variable in the diagram (counts as whole numbers, presence as 0/1), optional grouping column (site / burn unit) for random intercepts or for a multigroup comparison
engine: piecewiseSEM 2.3 (component lm / glm / MASS::glm.nb / lme4 models + Shipley's d-separation test) or lavaan 0.7 (ML / MLR)
---

# Structural equation models (SEM)

## Question it answers

"How does *X* affect *Y*, and how much of that goes through *M*?" when you
have a diagram of hypothesised causal links among measured variables. For
example: fire severity reduces canopy cover; less canopy lets shrubs grow;
shrubs support more birds; fire may also reduce birds directly. An SEM fits
every arrow at once, tests whether the arrows you left out are really
unnecessary, and splits the association of fire with birds into a
**direct** part and an **indirect** part (through canopy and shrubs).

Two engines:

| | `piecewise` (default) | `lavaan` |
|---|---|---|
| Responses | Gaussian, Poisson, negative binomial, binary (logit / probit) | Gaussian (continuous) only |
| Random site intercepts `(1\|site)` | yes | no |
| Latent variables (e.g. "soil fertility" measured by N, C, organic matter) | no | yes |
| Multigroup (does a path differ between burned and unburned?) | no | yes (`group:`) |
| Global test | Shipley's d-separation test (Fisher's C) | chi-square, CFI, TLI, RMSEA, SRMR |
| Missing links reported as | d-sep claims with P < 0.05 | modification indices > 3.84 |

## Model

**Piecewise** (Shipley 2009; Lefcheck 2016). Each endogenous variable *y_k*
(one with at least one arrow into it) gets its own generalised linear
(mixed) model with its parents as predictors:

    g_k(E[y_k]) = a_k + sum_{j in parents(k)} b_jk x_j  (+ u_site)

The global model is checked with the test of directed separation. For every
pair of variables with no arrow between them (the basis set), the omitted
link is added to the child's model, conditioning on both variables' parents;
its P value tests "these two are independent given their causes". The *k*
P values are combined as Fisher's C = -2 sum ln(P_i), which follows a
chi-square distribution with 2k df when the model is correct (Shipley 2000).
Each claim is tested in the model of the causally later variable (Shipley's
basis set; `dsep_direction: causal`, the default). piecewiseSEM's
alternative for claims involving non-Gaussian variables, testing both
directions and keeping the smaller P value (`dsep_direction: conserve`), is
available. In every design validated here the two options gave the same
test (piecewiseSEM duplicated a claim only for two unlinked non-Gaussian
responses, and there both held 5.6-6.1% rejection of a true model;
VALIDATION.md).

**lavaan** (Rosseel 2012). All equations are fitted together by maximum
likelihood on the covariance matrix; latent variables are defined by
`=~` (first loading fixed to 1, so the latent variable is in the units of
its first indicator).

**Standardised coefficients** (Grace et al. 2018; Lefcheck's SEM book,
chapter "Coefficients"): b_std = b * SD(x) / SD(y), where SD(y) is the
observed SD for Gaussian responses, sqrt(var(eta) + pi^2/3) (logit) or
sqrt(var(eta) + 1) (probit) for binary responses (latent-theoretic), and
SD(eta) / cor(y, mu-hat) for counts (observation-empirical). These match
`piecewiseSEM::coefs()` exactly (VALIDATION.md). lavaan uses its `std.all`
solution.

**Indirect and total effects.** For a pair (X, Y), every directed route is
found; each route's effect is the product of its path coefficients; the
indirect effect is the sum over routes with at least one mediator, and the
total effect = direct + indirect. Unstandardised indirect effects are on the
**link scale of Y** (log expected count, log odds, or Y's units for a
Gaussian Y). They are exact effects when every mediator is Gaussian (identity
link); when a mediator is a count or binary variable, the unstandardised
product mixes scales and the module warns that only the standardised value
is comparable. Intervals: Monte Carlo (default; Preacher & Selig 2012:
4,000 draws from the joint normal sampling distribution of the coefficients,
block-diagonal across components for piecewise and the full lavaan
covariance otherwise), delta method, or nonparametric bootstrap (rows, or
whole sites when there is a random intercept).

## Assumptions (and what breaks if they fail)

1. **The diagram is causally correct and complete.** No unmeasured variable
   drives two measured ones. If one does, path estimates and indirect
   effects are biased; the d-sep / chi-square test detects this only when
   the omission implies an independence that the data contradict. Passing
   the test does not prove the diagram. Read paths as associations unless
   the design (experiment) or strong prior knowledge justifies a causal
   reading.
2. **Acyclic (piecewise).** No feedback loops.
3. **Each component model is right**: correct family (overdispersed counts
   -> `nbinom`; the module warns when Poisson dispersion > 1.5), linear on
   the link scale, independent rows apart from the random intercept.
4. **lavaan: multivariate normality** of the endogenous variables. Use
   `estimator: MLR` for robust SEs and a scaled chi-square with skewed data.
   Binary or count outcomes belong in the piecewise engine.
5. **Enough data.** Rules of thumb: >= 10 units per estimated coefficient in
   each component; latent variables with >= 3 indicators. Validation shows
   where the methods hold (VALIDATION.md).

## Data format

One CSV, one row per sampling unit:

```
plot_id,burn_unit,fire_severity_dnbr,canopy_cover_pct,soil_moisture_pct,shrub_cover_pct,bird_detections
P001,BU01,732,28.1,19.2,47.9,2
P002,BU01,633,33.4,23.3,48.1,7
```

* Every variable in the diagram must be a numeric column. Counts must be
  whole numbers >= 0; binary responses 0/1. Text columns may be used as a
  fixed-factor predictor (direct paths only; no indirect effects through a
  factor) or as a random-intercept grouping variable `(1|burn_unit)`.
* Missing values: every component must use the same rows. The module stops
  and lists them; `data: drop_missing: true` drops incomplete rows.
* Transformations (log area, rescaling) go in `data: transform:`.

Examples: `example/burned_forest_plots.csv` (piecewise, 120 plots in 12
burn units) and `example/grassland_plots.csv` (lavaan latent variable, 150
plots).

## Config

```yaml
data:
  file: burned_forest_plots.csv     # path relative to the config
  drop_missing: false               # true = drop rows with NA in model variables
  transform:                        # optional new columns (R expressions of existing columns)
    dnbr_100: fire_severity_dnbr / 100
engine: piecewise                   # or lavaan
model:                              # piecewise: one entry per response
  - formula: canopy_cover_pct ~ fire_severity_dnbr
    family: gaussian                # gaussian | poisson | nbinom | binomial
  - formula: bird_detections ~ shrub_cover_pct + fire_severity_dnbr + (1|burn_unit)
    family: poisson
    link: log                       # optional; binomial may use probit
lavaan: |                           # lavaan engine: model syntax instead of model:
  fertility =~ soil_om_pct + soil_c_pct + soil_n_x10
  biomass_100g ~ fertility + grazing_au_ha
group: pasture                      # lavaan only: multigroup model by this column
estimator: ML                       # lavaan only: ML or MLR
effects:                            # pairs for indirect/total effects (default: every pair joined by a route with a mediator)
  - {from: fire_severity_dnbr, to: bird_detections}
units: {fire_severity_dnbr: 100}    # meaningful step for the text ("per 100 dNBR")
labels: {fire_severity_dnbr: "Fire severity (dNBR)"}   # used in text and figures
settings:
  ci_method: montecarlo             # montecarlo (default) | delta | bootstrap
  n_mc: 4000                        # Monte Carlo draws
  n_boot: 500                       # bootstrap resamples
  dsep_direction: causal            # d-sep claim direction: causal (default) | conserve (piecewiseSEM: both directions, min P)
figures: {...}                      # see skills/eco-figures
```

## Outputs

Estimates table (`param` names):

| name | meaning |
|---|---|
| `path_<x>_to_<y>` | unstandardised coefficient on the link scale of y (factor levels: `path_<x>_<level>_to_<y>`) |
| `std_<x>_to_<y>` | standardised coefficient |
| `rate_ratio_<x>_to_<y>` / `odds_ratio_<x>_to_<y>` | exp(path) per unit x, for log / logit links |
| `direct_`, `indirect_`, `total_<x>_to_<y>` | effects for each pair (link scale of y) |
| `std_direct_`, `std_indirect_`, `std_total_<x>_to_<y>` | standardised versions |
| `rate_ratio_total_<x>_to_<y>` | exp(total) per unit x (count outcome, Gaussian mediators) |
| `loading_<latent>_<indicator>`, `std_loading_...` | lavaan factor loadings |
| `cov_<a>_<b>` | lavaan (residual) covariances; user-defined `:=` parameters keep their names |
| lavaan multigroup | every name gets the suffix `_<group>` |

Diagnostics notes: d-sep result and every missing link with P < 0.05
(piecewise) or chi-square/CFI/TLI/RMSEA/SRMR and the largest modification
indices (lavaan); R2 per response; convergence, singular random effects,
overdispersion, inadmissible lavaan solutions.

Figures: see the next section.

## Figures

`core/run.R fit` exports the suite to `<out>/figures/` (PDF + PNG and
`gallery.html`); `core/run.R figures sem --fit <out>/sem_fit.rds --config
<config>` redraws it after you edit the `figures:` block (no refit). Ids:

| id | engine | what it shows |
|---|---|---|
| `path_diagram` | both | The diagram with standardised coefficients. Arrow width is proportional to \|coefficient\|, colour gives the sign (blue positive, orange negative with the default palette), and arrows are dashed when the `ci` interval includes 0. Arrows that skip a level are drawn as arcs so they do not cross intermediate boxes. Dotted arrows labelled "missing? P = ..." are omitted links flagged by the d-sep test (piecewise) or by modification indices > 3.84 (lavaan). Shaded boxes are latent variables. |
| `indirect_effects` | both | Forest plot of standardised direct, indirect (summed over mediated routes) and total effects for each pair in `effects:`, with intervals at `ci`. Colour and shape both mark the effect type, so it works in greyscale. |
| `dsep` | piecewise | P value (log scale) of each independence claim in the basis set. Claims with P < 0.05 are highlighted triangles (the same colour as the dotted "missing" arrows); Fisher's C, df and P are in the caption. |
| `missing_paths` | lavaan | The largest modification indices (shown only if any exist). |
| `loadings` | lavaan | Standardised factor loadings with intervals; a dotted reference line at 0.7. |
| `partial_<x>_to_<y>` | both | Partial-residual plot for each numeric path x -> y, on the response scale, with a `ci` band and other parents at their means (random effects at zero). Points: partial residuals (Gaussian); working residuals added to the conditional linear predictor and back-transformed (counts, log link; the y axis starts at 0); a rug of x (binary). lavaan: line through the means with the unstandardised slope; latent variables as centred factor scores. |
| `effects_panel` | both | All `partial_` plots composed with A/B/C tags. |

Labels: every variable name can be renamed with `figures: labels:` (the
top-level `labels:` block is used too). Other keys: `sign_positive`,
`sign_negative`, `ci_excludes_0`, `ci_includes_0`, `missing_link`,
`effect_direct`, `effect_indirect`, `effect_total`, `dsep_flagged`,
`dsep_consistent`, `factor_score`, and the axis titles `axis_std_effect`,
`axis_dsep_p`, `axis_dsep_claim`, `axis_mi`, `axis_std_loading`.

Units: the module does not standardise covariates, so the axes are in the
units of your CSV columns (after any `data: transform:`). If a column is
already a z-score, `figures: units: {<column>: {center: ..., scale: ...}}`
puts the axis back into original units.

Intervals: the stored intervals are 95% (Monte Carlo by default). With
`figures: ci:` set to another level, coefficient and effect intervals use
estimate +/- z x SE, and the captions say so. There is no inner band,
because these are not Bayesian fits.

## Simulation and power

`simulate()` generates data from any recursive diagram given in
`design$model` (same format as the config), with parameters named
`<x>_to_<y>` (path), `int_<y>`, `sd_<y>` (Gaussian residual SD),
`theta_<y>` (NB size), `re_sd_<y>`, `mean_<x>`, `sd_<x>` (exogenous).
`design$fit_model` fits a different diagram from the one that generated the
data (to study omitted paths); `design$engine: lavaan` + `design$lavaan` +
`design$latent` simulate and fit latent-variable models. Truth for
standardised coefficients is computed from a 200,000-row simulation of the
same model. Knobs: `n_units` (sampling units; write `n_units`, not `n`, in YAML
because YAML reads a bare `n` key as FALSE; on the command line
`--vary n=...` also works when `n_units` is not set), `n_groups`, any
path coefficient.
Power for the default design (fire -> canopy -> shrub -> birds) is in
VALIDATION.md: e.g. `power sem --vary n=30,50,...`.

## Pitfalls

* **A non-rejected model is not a proven model.** With small n the d-sep
  test has little power (VALIDATION.md: an omitted path with standardised
  size ~0.3 is detected in only part of the datasets at n = 40).
* **Omitted paths bias indirect effects.** If X also affects M2 directly
  but the model omits it, the indirect effect of X through M1 -> M2 is
  wrong (VALIDATION.md negative control: about -40% bias).
* **Small samples with a binary outcome**: at n = 40, logit path
  coefficients are inflated by about 11-17% (the known small-sample bias of
  logistic maximum likelihood); the standardised values were not
  (VALIDATION.md). Report standardised effects, or collect more units.
* **Indirect effects through count or binary mediators** are not
  well-defined on the unstandardised scale; compare standardised values only.
  Both mediator types were recovery-tested at n = 150 (Poisson and logit
  mediators, PASS 24/24 each; VALIDATION.md "Mediator designs"), and an
  independent check added NB-count and probit mediators (PASS). The raw
  product can even have the opposite sign to the standardised total (seen on
  the example data with shrub cover coded present/absent), so the report text
  describes these routes with standardised values only. A yes/no mediator must
  be coded 1/0.
* **Delta-method intervals** for products are symmetric, while the sampling
  distribution of a product is skewed when the indirect effect is small
  relative to its SE (Preacher & Selig 2012). At the validated designs delta
  coverage was still 0.94-0.95, but Monte Carlo (default) or bootstrap is the
  safer choice at small n.
* **Collinearity** between parents of the same response makes paths
  imprecise; standardised coefficients near or above 1 are a warning sign.
* **Scale differences** of >1000x among variables can stop lavaan from
  converging: rescale with `transform:`.
* **Many independence claims**: individual claims with P < 0.05 can occur
  by chance; judge the model by Fisher's C first, then add a missing link
  only if it makes biological sense.
* Interactions and nonlinear paths are not supported (indirect effects
  through them are not defined by products of coefficients).

## References

* Grace JB (2006). Structural Equation Modeling and Natural Systems. Cambridge University Press.
* Grace JB, Anderson TM, Olff H, Scheiner SM (2010). On the specification of structural equation models for ecological systems. Ecological Monographs 80(1):67-87.
* Shipley B (2000). A new inferential test for path models based on directed acyclic graphs. Structural Equation Modeling 7(2):206-218.
* Shipley B (2009). Confirmatory path analysis in a generalized multilevel context. Ecology 90(2):363-368.
* Shipley B, Douma JC (2020). Generalized AIC and chi-squared statistics for path models consistent with directed acyclic graphs. Ecology 101(3):e02960.
* Lefcheck JS (2016). piecewiseSEM: Piecewise structural equation modelling in R for ecology, evolution, and systematics. Methods in Ecology and Evolution 7(5):573-579.
* Lefcheck JS. Structural Equation Modeling in R for Ecology and Evolution. https://jslefche.github.io/sem_book/ (chapters "Local estimation" and "Coefficients").
* Grace JB, Johnson DJ, Lefcheck JS, Byrnes JEK (2018). Quantifying relative importance: computing standardized effects in models with binary outcomes. Ecosphere 9(6):e02283.
* Rosseel Y (2012). lavaan: An R package for structural equation modeling. Journal of Statistical Software 48(2):1-36. doi:10.18637/jss.v048.i02
* Preacher KJ, Selig JP (2012). Advantages of Monte Carlo confidence intervals for indirect effects. Communication Methods and Measures 6(2):77-98. doi:10.1080/19312458.2012.679848
* Hu L, Bentler PM (1999). Cutoff criteria for fit indexes in covariance structure analysis. Structural Equation Modeling 6(1):1-55.
* Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1):1-48.
