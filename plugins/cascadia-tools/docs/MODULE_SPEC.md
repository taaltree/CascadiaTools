# Writing an ecotools module

A module is one analysis family in `modules/<id>/`. Users never call your
code directly: an agent driving ecotools (see `CLAUDE.md`) reads your
`README.md` to decide whether your module fits the researcher's question,
builds a config, and runs `core/run.R`. Write for that reader.

## Files

```
modules/<id>/
  module.R        assigns `module <- ecomodule(...)`  (contract: core/R/contract.R)
  README.md       frontmatter + user-facing guide (format below)
  model.R / *.py  engine code (Nimble code, NumPyro model, helpers) sourced via MODULE_DIR
  example/        a small example dataset in the REAL-DATA format (CSV) + config.yaml
  VALIDATION.md   written by you: recovery results, controls, known limits
```

Look at `modules/_template/module.R` first: it is the smallest complete
module. Core helpers available inside module.R (core is loaded first):
`eco_require()`, `nimble_run()`, `summarise_draws_eco()`, `stack_draws()`,
`derived_rows()`, `wald_rows()`, `merge_list()`, `eco_python()`, `MODULE_DIR`.

## Contract essentials

- `simulate(params, design, seed)` returns `list(data, truth)`. `data` has
  exactly the structure `read_data(config)` returns, so the same `fit()`
  serves simulated and real data. `truth` is a named numeric vector; the
  names are the estimand names in `fit()$estimates$param`.
- `fit(data, settings)` returns `estimates` (param, estimate, sd, lower,
  upper; 95% intervals), `diagnostics` (converged, rhat_max, ess_min, notes)
  and optionally `object`.
- Estimands should be on the scale biologists care about: occupancy
  probability, abundance, density per km2, survival, selection coefficient,
  Fst... and also the regression coefficients on the link scale (these are
  what power is computed for). Name coefficients descriptively
  (`psi_beta_forest`, not `b[2]`).
- `power_params` lists effect parameters; `null_values` gives the no-effect
  value when it is not 0 (e.g. Fst = 0, rho = 0, ratio = 1).
- `fast_settings`: the cheapest settings that still give honest intervals
  (fewer MCMC iterations, fewer chains); recovery uses them.
- `read_data(config)`: read the researcher's own files. Accept the formats
  ecologists actually have (long-format detection CSVs, site covariate CSVs,
  GPS track CSVs, genotype tables/VCF...). `config$.dir` is the folder of the
  config file; resolve relative paths against it. Document every column.
- `check_data(data)`: catch real-world mistakes with messages a non-coder can
  act on ("site S14 has 7 survey columns but the others have 6").
- `interpret(fit, data)`: plain-language paragraphs with numbers and units
  ("each 10% increase in forest cover raises occupancy from 0.42 to 0.49").
  Never claim causation unless the module is a causal module whose
  assumptions are stated.
- `figures(fit, data, opts)`: the full publication figure suite for your
  analysis type, exactly as listed in `docs/FIGURES.md` (read it; use the
  `fig_*` constructors in core/R/figures.R; honour opts). Export it from
  the example fit and LOOK at every PNG (Read tool) before declaring done.
- `plots(fit, data)`: (legacy) named list of ggplot objects. Effects of covariates
  are shown as conditional-effects curves on the response scale with 95%
  bands by default, other covariates at their mean.
- `methods_text()`: a paragraph that could go into a paper's methods.
- Bayesian modules: use Nimble via `nimble_run()` (pass covariates as
  `data`, dimensions as `constants`, give a `cache_key` built from all
  dimensions so recovery reuses the compiled model). Report R-hat and bulk
  ESS. Default settings for real data must reach R-hat < 1.01 and bulk ESS
  >= 400 on the example data. Priors: state them in README and methods_text.
- Engines must be deterministic given `seed`.

## README.md format

```
---
id: occu-single
title: Single-species occupancy
when_to_use: Detection/non-detection surveys repeated at sites; question is where a species occurs and what drives it, accounting for imperfect detection.
not_for: Counts (use nmixture), multiple interacting species (occu-rota / occu-waddle).
data_needed: detection history (site x survey 0/1), site covariates, optional survey covariates
engine: unmarked (likelihood) or nimble (Bayesian)
---
```
Then sections: Question it answers; Model (equations in plain text/LaTeX);
Assumptions (and what breaks if they fail); Data format (columns, with a
tiny example); Config (every key); Outputs; Simulation and power (which
design/params knobs exist and what is sensible); Pitfalls; References.

## Validation you must do before declaring done (protocols/acceptance-gate.md)

Use `scripts/Rscript` (R 4.5.3), never bare `Rscript`. `/scratch` is node-local,
so put throwaway scripts in `dev/<id>/`, not scratch. Run everything heavy
through SLURM, never on the login node (it has one core):
`scripts/ejob -c 16 -m 32 -t 4:00:00 -- scripts/Rscript core/run.R ...`

1. `selftest` (4 reps) passes the contract.
2. `recover` with nsim >= 100 at the default design: every estimand passes
   `recovery_gate()` (bias within max(10%, 3 MC SE), coverage >= bound,
   >= 90% converged). Record the table in VALIDATION.md.
3. A second, harder design (smaller sample / rarer species / stronger
   effect). Report honestly where it fails; that is the module's limit.
4. Negative control: fit data simulated with an effect = 0 and confirm
   power ~ 0.05 (false-positive rate); and show the gate fails if you fit a
   deliberately mis-specified model (e.g. ignore detection) - this proves
   the check can fail.
5. Where an established package exists (unmarked, spOccupancy, secr,
   lavaan, hierfstat...), cross-check your fit against it on one dataset
   (independent route). Where you use the package as the engine, cross-check
   against a hand-written likelihood or Nimble model instead.
6. `fit` on `example/config.yaml` produces a report; open the HTML and check
   that the interpretation sentences are correct for the numbers shown.
7. `power` over one design knob produces a sensible curve.

Write results with numbers, not adjectives. If a gate fails and you cannot
fix it, say so in VALIDATION.md under "Known limits": an honest OPEN beats a
false PASS.
