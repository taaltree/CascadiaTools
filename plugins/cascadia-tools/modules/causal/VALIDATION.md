# VALIDATION: causal (DAG adjustment + GLM g-computation, fixed-effects panels)

Date: 2026-10-18. R 4.5.3, dagitty 0.3.4, marginaleffects 1.0.0, MASS, sandwich.
Code state for every recovery/power number below: `dev/causal/code_md5_final_suite.txt`.
After that, only `output.R` changed: DAG-figure edge drawing, then (2026-10-19)
extra data-intake checks in `causal_read_data()`/`causal_check_data()`.
Recovery never calls either. Final hashes: `dev/causal/code_md5_after_input_checks.txt`.
An independent re-run of default, hard_small, null_effect and force_collider
(`dev/causal/verify/run_all.sh`, `logs/causal_verify.log`) reproduced every number below.

Gate: `recovery_gate()` with its fixed thresholds: |bias| <= max(10% of truth,
3 MC SE); coverage >= 0.95 - 3 MC SE (0.904 at nsim 200, 0.921 at 500);
>= 90% converged. Recovery uses `fast_settings` (real-data settings minus the
what-if table, local tests and DHARMa, which do not change the estimates).
Truths come from the structural causal model by intervention (Monte Carlo,
4e5 draws, common random numbers; MC error about 0.2%).

Reproduce:
- `bash dev/causal/final_suite.sh`: every config in `modules/causal/validation/` plus two power runs; outputs in `results/causal/validation/<run>/`, `results/causal/power*/`.
- `scripts/Rscript dev/causal/crosscheck.R`
- `scripts/Rscript dev/causal/calib_checks.R`
- `scripts/Rscript dev/causal/tables.R` (the tables below)

All 22 runs took about 6 minutes on 16 cores. A fit takes 0.3-0.9 s; 2-3.5 s with the what-if table.

## Default design (`default`): PASS on every estimand

Stream DAG (12 arrows: elevation, agriculture -> sediment, width as confounders;
water temperature as mediator; periphyton as collider caused by canopy and trout).
Standardised path coefficients; total effect of canopy on log E[trout] 0.36
(0.20 direct + 0.16 via temperature). Negative binomial trout counts (theta 3,
mean about 6.7), 120 reaches, contrast +1 SD of canopy. 200/200 converged.

| estimand | truth | mean est | rel. bias % | bias MCSE | coverage | CI width | power | gate |
|---|---|---|---|---|---|---|---|---|
| ate_canopy | 2.884 | 2.929 | 1.6 | 0.0674 | 0.955 | 3.66 | 0.965 | PASS |
| beta_canopy | 0.36 | 0.359 | -0.3 | 0.00625 | 0.945 | 0.34 | 0.985 | PASS |
| rate_ratio_canopy | 1.433 | 1.438 | 0.3 | 0.00913 | 0.945 | 0.492 | 0.985 | PASS |

## Other families, effects and estimators (nsim 200 each): all PASS

| run | design | result |
|---|---|---|
| gaussian_direct | Gaussian outcome, 100 sites, total **and controlled direct** effect | PASS all 4: total 0.3626 vs 0.36 (+0.7%, cov 0.96); direct 0.1976 vs 0.20 (-1.2%, cov 0.96) |
| binary_binomial | grazing DAG (10 arrows, collider "erosion", mediator "litter"), **binary exposure**, 0/1 outcome, 300 pastures | PASS all 4: ATE -0.188 vs -0.194 (+2.9%, cov 0.965); risk ratio 0.673 vs 0.664 (cov 0.96); potential-outcome means 0.572/0.384 vs 0.577/0.383 (cov 0.96/0.95) |
| poisson_default | Poisson counts, sandwich SEs (default for poisson) | PASS all 3: bias <= 1.1%, coverage 0.945-0.95 with HC3 (default since 2026-10-19; HC0 gave 0.935) |
| poisson_optimal | `adjust: optimal` (Henckel et al. 2022) | identical to poisson_default: in this DAG the optimal set equals the smallest minimal set {elevation, sediment, width} |
| panel_twoway | Byrnes & Dee-style panel: 30 sites x 6 years, unmeasured site confounder U and year shock, observed time-varying confounder; Poisson; **two-way FE**, cluster-robust SE | PASS all 3: beta 0.2546 vs 0.25 (+1.9%, cov 0.97); ATE 2.067 vs 2.054 (+0.6%, cov 0.94); rate ratio cov 0.97 |
| panel_fixed_gaussian | Gaussian, 20 sites x 5 years, unit FE | PASS: 0.2552 vs 0.25 (+2.1%, cov 0.93) |
| ci_delta / ci_bootstrap | default design with delta-method (conditional) or bootstrap (200 resamples) intervals | both PASS; coverage ATE 0.955 / 0.945, beta 0.945 / 0.940; mean ATE CI width 3.55 (delta) vs 3.66 (default unconditional) vs 3.74 (bootstrap) |

## Harder design (`hard_small`, 40 sites): FAIL on the response-scale ATE

| estimand | truth | mean est | rel. bias % | bias MCSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|
| ate_canopy | 2.884 | 3.40 | **+17.9** | 0.141 | 0.935 | 0.43 | **FAIL (bias)** |
| beta_canopy | 0.36 | 0.381 | 5.7 | 0.0124 | 0.915 | 0.68 | PASS |
| rate_ratio_canopy | 1.433 | 1.486 | 3.7 | 0.0185 | 0.915 | 0.68 | PASS |

With 40 sites the count-scale ATE is right-skewed and its mean is biased upward,
though coverage is fine (0.935). The power curve shows it shrinking with n:
+19.7% at 40 sites, -3.5% at 60, +5.8% at 90, +3.7% at 120. See Known limits.

## Negative controls: the gate CAN fail

### Forced wrong adjustment on the main estimate (nsim 200)

| run | what was done | ate_canopy | beta_canopy | rate_ratio_canopy |
|---|---|---|---|---|
| force_mediator | also adjust for water_temp (mediator) | -45.9%, cov 0.57 **FAIL** | -43.0%, cov 0.625 **FAIL** | -13.9%, cov 0.625 **FAIL** |
| force_omit_confounder | drop width (confounder) | -34.0%, cov 0.67 **FAIL** | -31.8%, cov 0.70 **FAIL** | -10.4%, cov 0.70 **FAIL** |
| force_collider | also adjust for periphyton (collider) | +18.2%, cov 0.97 **FAIL (bias)** | +13.0%, cov 0.925 **FAIL (bias)** | +5.1%, cov 0.925 PASS |
| panel_pooled | panel data, ignore site and year (pooled GLM) | +154.6%, cov 0.205 (HC3; 0.18 with HC0) **FAIL** | +115.1%, cov 0.15 **FAIL** | +35.0% **FAIL** |
| panel_oneway | site FE only, year shocks ignored | -131.0%, cov 0.065 **FAIL** | -138.1%, cov 0.10 **FAIL** | -28.4% **FAIL** |

The collider row shows a real limitation of the gate. At this collider
strength the bias on the ratio scale (+5.1%) is inside the 10% tolerance and
coverage stays at 0.925. A moderate collider can therefore slip past a
ratio-only check, while the same error fails on the difference and the
coefficient. Report all three estimands.

### Built-in what-if table scored against the truth (`whatif_controls`, `binary_whatif`, `panel_whatif`; nsim 200)

Every what-if row targets the same true causal effect, and the DAG labels each
strategy valid or biased in advance. **In every case the DAG's verdict matched
the gate's.** The one DAG-valid alternative strategy (and the 3 main DAG-set
estimates) passed, and all 19 DAG-biased strategies failed (7 stream, 6 grazing, 6 panel).

Stream DAG (truth ATE 2.884):

| strategy | DAG says | rel. bias % | coverage | gate |
|---|---|---|---|---|
| DAG set {elevation, sediment, width} | valid | +1.6 | 0.955 | PASS |
| alternative minimal set {agriculture, elevation, width} | valid | +3.6 | 0.935 | PASS |
| no adjustment | biased | +47.0 | 0.885 | FAIL |
| omit elevation | biased | +41.4 | 0.905 | FAIL (bias) |
| omit sediment | biased | +58.0 | 0.865 | FAIL |
| omit width | biased | -34.0 | 0.670 | FAIL |
| + mediator water_temp | biased | -45.9 | 0.570 | FAIL |
| + collider periphyton | biased | +18.2 | 0.970 | FAIL (bias) |
| all measured variables | biased | -24.9 | 0.805 | FAIL |

Grazing DAG, binary exposure (truth ATE -0.194): DAG set +2.9% PASS. Biased:
crude +29.8%, omit soil_moisture +53.2%, omit slope -18.3%, + mediator litter
+37.0%, + collider erosion +91.2% (coverage 0.09), all measured +115.3% (sign
flipped, coverage 0.015). All FAIL.

Panel DAG (truth ATE 2.054): two-way FE +0.6% PASS. Biased: no FE (pooled)
+154.6%, omit site FE (year FE only) +318.1% (coverage 0.000), omit year FE
-131.0%, omit upwelling -100.4%, crude +80.8%, all measured (no FE) +154.6%.
All FAIL.

### False-positive rate (`null_effect`, nsim 500): canopy has no effect (direct and via temperature both 0)

| estimand | rejection rate | coverage |
|---|---|---|
| ate_canopy | 0.046 | 0.954 |
| beta_canopy | 0.050 | 0.950 |
| rate_ratio_canopy | 0.050 | 0.950 |

All are consistent with 0.05 (MC SE 0.0097).

### Poisson standard errors

The total-effect model leaves the mediator and the noise in the outcome's
other causes out of the regression, so counts are overdispersed given
(X, Z) even when the outcome is Poisson given all its parents. On default-design
Poisson data the plain-GLM Pearson dispersion was 1.27-1.39, and DHARMa flagged
dispersion in 53% of correctly specified fits (calibration below).
Model-based SEs (`poisson_modelvcov`) gave ATE coverage 0.920 (gate bound 0.904:
PASS, but low); sandwich SEs, now the default for poisson (`poisson_default`),
gave 0.935 with HC0 and 0.950 with HC3 (the default since the independent verification, below). negbin with model-based SEs gave 0.955.

## Calibration of the model checks (`dev/causal/calib_checks.R`, 300 datasets each)

| scenario | Fisher's C p < 0.05 | any Holm-adjusted local test p < 0.05 | DHARMa uniformity / dispersion / outliers flagged |
|---|---|---|---|
| correct DAG, negbin | 0.060 | 0.043 | 0.000 / 0.013 / 0.003 |
| correct DAG, poisson | 0.050 | 0.053 | 0.003 / **0.533** / 0.080 |
| correct DAG, gaussian | 0.060 | 0.050 | 0.000 / 0.000 / 0.003 |
| WRONG DAG (data have elevation -> trout and agriculture -> elevation; analysed without them) | **0.950** | **0.960** | 0.000 / 0.007 / 0.000 |
| Poisson model fitted to negbin data | - | - | 0.990 / **1.000** / - |

The local tests hold their nominal size and detect two missing arrows 95% of
the time at 120 sites. DHARMa is conservative when the family is right (dispersion
flagged 0-1.3%) and flags a Poisson fit to overdispersed counts 100% of the time.

## Cross-checks against independent implementations (`dev/causal/crosscheck.R`)

1. **Adjustment sets, dagitty vs a from-scratch implementation** of the
   adjustment criterion (Bayes-ball d-separation on the proper back-door graph,
   `adjustment_criterion()` in dag.R). The test covered every subset of measured
   covariates for 5 named DAGs (stream, grazing, M-bias, butterfly, latent
   front-door) and 92 random 8-node DAGs with 0-2 latent nodes. Result:
   2,592/2,592 candidate sets agree on validity, and minimal-set lists are identical in 92/92 DAGs.
2. **Regression with the adjustment set vs a full lavaan path model** (Gaussian
   SCM, 50 datasets of 500 sites, truth 0.36 total / 0.20 direct). Mean total effect:
   module 0.3626, lavaan (product of paths) 0.3621; mean |difference| 0.011
   (max 0.036), correlation 0.94. Direct effect: identical to 1e-4. Mean SE: module
   0.0460 vs lavaan 0.0432; the full SEM is more efficient because it uses
   every arrow, but it is valid only if the whole DAG is right.
3. **g-computation by hand vs marginaleffects** (example data, negbin):
   ATE 1.181192 (SE 0.467750) in both; rate ratio 1.187395 (SE 0.072162) in both.
4. **Fixed effects**: the unit-dummy GLM equals the hand within-transformation
   (beta 0.213901 both), and sandwich::vcovCL HC1 equals a hand CR1 cluster SE
   (0.107669 both).
5. **Local tests**: dagitty partial correlations equal hand residual-on-residual
   correlations, max |difference| 8.6e-16 over 16 tests.
6. **Interval methods on the example data**: SE of ATE is 0.468 (delta), 0.471
   (default unconditional) and 0.646 (2000-resample bootstrap); for beta, 0.0061
   vs 0.0080 (bootstrap). On this dataset the bootstrap is 30-40% wider (DHARMa
   dispersion p = 0.016 here: a heavy-tailed draw). Across 200 simulated
   datasets the three methods have the same coverage (0.945-0.955), so the
   delta/unconditional default is calibrated. The bootstrap remains available
   (`ci_method: bootstrap`) for small or skewed samples.

## Power (`core/run.R power`, nsim 100 per point)

Over n_sites, at the default effects (`results/causal/power/`):

| n_sites | 20 | 40 | 60 | 90 | 120 | 200 |
|---|---|---|---|---|---|---|
| power, rate ratio / beta | 0.41 | 0.70 | 0.79 | 0.94 | 0.98 | 1.00 |
| power, ATE | 0.13 | 0.42 | 0.65 | 0.87 | 0.96 | 1.00 |
| coverage (beta) | 0.84 | 0.90 | 0.90 | 0.94 | 0.92 | 0.94 |

About 60 sites gives 80% power for the rate ratio and about 90 for the ATE.
Below 60 sites coverage drops under 0.90; do not trust intervals there.

Over the direct path b_canopy__trout at 120 sites (`results/causal/power_effect/`;
the indirect path via temperature stays at 0.16, so 0 is not a null):
power for beta 0.50, 0.88, 0.99, 0.99 at 0, 0.1, 0.2, 0.3.

## Example data (simulated, known truth: `dev/causal/make_examples.R`)

- `example/config.yaml`: 140 reaches, negbin, +10% canopy. Truth vs estimate (95% CI):
  - ATE 1.95 vs 1.18 (0.26-2.10)
  - rate ratio 1.293 vs 1.187 (1.054-1.338)
  - beta 0.0257 vs 0.0172 (0.0053-0.0291)
  - direct ATE 1.02 vs 0.57 (-0.30-1.43)

  All intervals cover the truth. The 16 local tests have Fisher's C p = 0.50.
  The interpretation sentences in `results/causal/example/causal_report.html`
  were checked against the estimates table.
- `example/config_panel.yaml`: 24 sites x 6 years, two-way FE. Truth vs estimate:
  - beta 0.227 vs 0.219 (0.089-0.348)
  - ATE 3.71 vs 4.74 (1.28-8.21)

  Pooled GLM: 13.6. The DAG has no testable implications among measured variables, which the report states.

## Data-intake checks on messy files (`dev/causal/verify/run_messy.sh`, `logs/causal_messy.log`)

Each broken copy of the example data stops before fitting with a specific message:

| file problem | message (abridged) |
|---|---|
| decimal comma / stray text in a covariate (`12,5`, `16?`) | 'fine_sediment_pct' is numeric except for value(s) '12,5'; fix or blank them |
| column renamed (`wetted_width` vs DAG node `wetted_width_m`) | node would be treated as UNMEASURED; column 'wetted_width' looks similar |
| exposure misspelt in config | not a node of the DAG / not a column (lists both) |
| non-integer or negative counts with negbin | needs non-negative whole-number counts |
| user `adjust:` that misses a confounder / includes a mediator | not a valid set; lists the valid minimal sets |
| rows with a different number of fields | header has k fields but line(s) ... differ |

Before 2026-10-19 the first two cases gave an unhelpful "text with 93 distinct
values" message and a misleading "not identifiable" message; both are fixed.
Missing values in DAG variables are dropped (complete cases) and the count is
stated in the report.

## Known limits

1. **Small samples, count outcomes**: the response-scale `ate` is biased upward
   (+18% at 40 sites, FAIL), and coverage falls below 0.90 under about 60 sites.
   `beta` and the ratio pass at 40. The report warns when a count dataset has
   fewer than 80 rows. Use at least 60-90 sites, or report the ratio.
2. **A moderate collider can pass on the ratio scale** (+5.1% bias, inside
   tolerance) while failing on the ATE and beta. The DAG, not the data, must exclude it.
3. **The DAG is untestable where it matters most**: local tests detect missing
   arrows between measured variables (95% power here), but not an unmeasured
   confounder of exposure and outcome. The E-value is reported for that reason.
   Local tests on panel rows ignore clustering, so their p-values are too small.
4. **Poisson total-effect models are overdispersed by construction** (DHARMa
   flags 53% of correct fits). The module uses sandwich SEs for poisson; negbin is preferred.
5. **Not implemented**:
   - front-door and instrumental-variable estimators;
   - exposure x covariate interactions (effect heterogeneity);
   - natural direct/indirect effects (only controlled direct effects, and the
     "share through mediators" sentence is a rough log-scale split);
   - categorical exposures with more than 2 levels;
   - time-varying site-specific confounding and dynamic feedback (past outcome -> future exposure);
   - Mundlak / correlated random effects;
   - binomial outcomes with unit fixed effects. They run, but are not validated;
     the incidental-parameter bias of logit FE with few repeats is expected.
6. **Simulation scope**: the SCM simulator makes every non-outcome node
   linear-Gaussian, with a binary exposure as the only alternative. Recovery
   with categorical or skewed confounders and non-linear confounding was not
   tested. A misspecified functional form for a confounder leaves residual bias.
7. `selftest` (nsim 4) shows coverage 0.75 because one of four intervals
   missed; that is Monte Carlo noise. nsim 200 is the evidence.

## Independent verification (2026-10-19, separate verifier)

The verifier wrote its own simulator (`dev/causal/verify/verify_sim.R`), which uses no module code to
generate data or compute truths. Truths come from 2e6-draw interventions in the verifier's own SCMs.
Every dataset was fitted with the module's `fit()` (fast settings) and scored with core
`summarise_recovery()` / `recovery_gate()`. Logs: `logs/causal_verify.log`, `logs/causal_verify_final.log`,
`logs/causal_verify_hc.log`, `logs/causal_verify_hc3.log`.

**Literature check.**
- `causal_optimal_set()` is exactly Henckel et al. (2022) Def. 3.12: O = pa(cn) \ forb, with forb = de(cn) ∪ {X}.
- E-value: RR + sqrt(RR(RR-1)), using the CI limit nearest 1; continuous outcomes use RR ≈ exp(0.91 d) (VanderWeele & Ding 2017).
- Fisher's C: -2 Σ ln p on 2k df (Shipley 2000).
- The unconditional SE adds the influence-function variance of the g-computation mean, summed by unit for panels.
- The simulator's truth is computed by intervention on the SCM, not from fit code. The only code shared with `fit()` is
  `causal_identify()`, which is used to label the what-if truths.

**Independent planted truth (designs differ from the builder's; 200 datasets each, 2000 for the null).**

| scenario | design | result |
|---|---|---|
| A_poisson_offset | camera-trap Poisson with effort offset (`offset_per` 30 d), natural-unit exposure (forest %), mediator + collider + precision covariate, 200 sites | PASS: ATE +4.0% (cov 0.96), beta +1.6% (0.945), rate ratio +0.7% (0.945) |
| B_gamma_binary | binary exposure (logged), Gamma biomass, 3-level text confounder, lognormal confounder, 250 stands | PASS: ATE +4.9% (cov 0.93), beta +5.0% (0.925), mean ratio +1.3% (0.925), PO means ±0.5% (0.94/0.95) |
| F_gauss_small | Gaussian, 60 sites, contrast 0.5 | PASS: -0.6%, cov 0.955 |
| C_null (nsim 2000) | scenario A with zero effect | PASS. Rejection rate 0.053-0.054 with HC3. **With the original HC0 it was 0.062-0.063** (MCSE 0.005; SE/empirical SD 0.955) |
| D_force_mediator | A's data, mediator forced in | FAIL as required: ATE -37%, beta -36%, ratio cov 0.855 |
| E_omit_confounder | A's data, precip omitted | FAIL as required: ATE +63% (cov 0.58), beta +50% |

Reruns of the builder's own commands reproduce the numbers above exactly, for the same seeds:
- default: ATE +1.6%, cov 0.955, PASS
- hard_small: ATE +17.9%, FAIL
- null_effect (nsim 500): rejection 0.046-0.050
- force_collider: ATE +18.2% and beta +13.0% FAIL; ratio +5.1% PASS

**What the verifier changed.**
1. *Poisson sandwich SEs HC0 -> HC3* (`engine.R`). HC0 was mildly anti-conservative at n = 200: the null rejection
   rate was 0.063 with 2000 replicates. HC3 gives 0.053, and SE/empirical SD goes to 0.99. Affected validations were rerun at nsim 200:
   - poisson_default / poisson_optimal: PASS, coverage 0.945-0.95
   - panel_pooled: still FAIL, +155%
   - verifier scenarios A, D and E: unchanged verdicts

   Panel (cluster HC1) and negbin/Gaussian/Gamma SEs are unaffected.
2. *Panel id checks* (`check_data`): flags unit ids that differ only in case or spaces, rare ids one or two characters from
   an existing id (a likely typo; the closest ids are named), and duplicated unit x time rows. Before this change, a misspelt site
   id silently became a new one-row site. Tested on `dev/causal/verify/messy2/panel_typo.csv`.
3. *Report text* (`output.R`):
   - a "Detection versus ecology" paragraph for count and binomial outcomes;
   - the simulation-limits sentence corrected to "bias within 5%, coverage 92-97%" so that it covers the independent scenarios;
   - exposure units added to the effect-plot x-axis and the what-if axis ("+10 percentage points of ...").

**Usability.**
- Both examples were re-run on the final code. I checked every number in the interpretation sentences against the estimates CSV,
  by hand: ATE, CI, PO means, rate ratio, % change, exp(beta), the 1 - direct/total split, what-if % differences, and E-values
  (1.66 / 1.29 for the stream example). All are correct.
- Causal language is conditional on the stated assumptions, and effects are given on the response scale first.
- Messy-data checks stop with actionable messages. These cases were tested:
  - decimal comma or stray text
  - renamed column
  - misspelt exposure
  - non-integer counts
  - invalid user adjustment set
  - panel id typo, case/space clash and duplicated row
- Missing values (blank or `NA` in a confounder and in the outcome: 24 of 140 rows) are dropped as complete cases,
  and the report states the count.
- The dose-response figure is a g-computation average-prediction curve over the 2-98% exposure range, with a rug of the observed values. It is correct.

**Remaining limits (in addition to the list above).** The nsim 4 `selftest` still reports FAIL on coverage (0.75; one of four).
That is expected Monte Carlo noise, not evidence. The methods paragraph still uses raw column names.
