---
name: eco-analyze
description: Run a complete, validated ecological / wildlife / fisheries / conservation-genetics analysis for a researcher who may not code. Use when someone brings data and a question ("does forest cover affect occupancy?", "estimate density from my camera traps", "how many loci do I need?").
---

> **Toolkit:** `${CLAUDE_PLUGIN_ROOT}` (scripts in `${CLAUDE_PLUGIN_ROOT}/scripts`, module index in `${CLAUDE_PLUGIN_ROOT}/catalog.md`). Work and outputs go in the researcher's project folder, never in the toolkit. If R packages are missing, run the `cascadia-setup` skill first; if the data are sensitive, offer `cascadia-private`.

# eco-analyze — from question to validated report

The researcher should never need to write code. You do the coding through
CascadiaTools modules, and you hold yourself to the protocol skills. Talk to the
researcher in biological terms; keep statistical jargon for the report's
methods section.

## Steps

1. **Project folder.** `projects/<short-name>/` with `data_raw/`, `data/`,
   `analysis_plan.md`, `config.yaml`, `recovery/`, `results/`, `log.md`
   (what you did, in order, with commands).
2. **Intake** (eco-data-intake): inspect files, describe back, confirm.
3. **Choose the module** from `${CLAUDE_PLUGIN_ROOT}/catalog.md` (read each candidate's
   `${CLAUDE_PLUGIN_ROOT}/modules/<id>/README.md`, especially `when_to_use` / `not_for`). If two
   fit, explain the trade-off in one or two sentences and recommend one. If
   none fits, say so; propose a new module (eco-new-module) rather than
   forcing a poor fit.
4. **Plan** (eco-analysis-plan): write `analysis_plan.md`; confirm with the
   researcher the estimand and covariates.
5. **Planted truth** (eco-planted-truth): recovery at their design. Report
   the gate result in one sentence. If it fails, stop and discuss.
6. **Fit**: write `config.yaml` (data paths, columns, settings) and run
   `"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/core/run.R" fit <module> --config config.yaml --out results/`
   (via `"${CLAUDE_PLUGIN_ROOT}/scripts/ejob"` for anything Bayesian or large).
7. **Check** (eco-convergence, eco-model-checking).
8. **Interpret** (eco-interpretation): the HTML report in `results/`
   contains estimates, conditional-effects plots and plain-language
   interpretation. Read it, verify every sentence against the numbers,
   and add project-specific context.
9. **Gate** (eco-acceptance-gate): give the verdict, DONE or OPEN, with
   reasons.
10. **Deliver**: summary for the researcher (3-6 sentences, numbers with
    intervals, caveats), path to the report, and a methods paragraph for
    their paper.

## Compute

Quick fits (most GLMs, likelihood occupancy and N-mixture models) run directly with
`"${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/core/run.R" ...`. Anything Bayesian, any recovery or power study, and
anything longer than a few minutes goes through `"${CLAUDE_PLUGIN_ROOT}/scripts/ejob" [-c cores] [-m GB] [-t time] -- <command>`.
On a SLURM cluster it submits the job and waits; on a laptop or workstation it runs the
command locally with the same logging. Measure a `selftest` first, then size the run.
For runs longer than about 8 minutes, start them in the background and check back,
telling the researcher roughly how long to wait.
