---
name: eco-analysis-plan
description: Turn a researcher's question into a written analysis plan (estimand, model, covariates, priors, checks) before any model is fit. Use at the start of every analysis and whenever the question changes.
---

# eco-analysis-plan — declare before fitting

Write `analysis_plan.md` in the project folder **before** fitting real data.
Fill it in by asking the researcher, in plain language, only what you
cannot infer from their data.

## Template

1. **Question** in the researcher's words, then restated as an estimand:
   "the change in occupancy probability per 10% forest cover, accounting
   for detection".
2. **Data**: what one row is, sampling design, dimensions, what was
   recorded, known problems (missed surveys, effort varies, GPS fix
   failures).
3. **Module and model**: which CascadiaTools module (see `${CLAUDE_PLUGIN_ROOT}/catalog.md`) and why
   it fits; which alternatives were rejected and why.
4. **Covariates**: each one with a biological reason. Scaled (mean 0, SD 1)
   unless there is a reason not to; collinearity check (|r| > 0.7 flagged).
   For causal questions: draw the DAG and use the `causal` module's
   adjustment-set output. Do not adjust for mediators or colliders.
5. **Priors** (Bayesian): module defaults unless the researcher has
   information; state them.
6. **Checks**: recovery design (eco-planted-truth), convergence thresholds,
   goodness-of-fit test, sensitivity analyses.
7. **Model selection**: decided now. Prefer one pre-specified model plus a
   small set of biologically motivated alternatives. No stepwise
   selection, no dredging all subsets then interpreting the winner's
   p-values as if pre-specified.

## Rules

- If the question is causal ("does logging reduce..."), say what design or
  assumptions would justify a causal reading; otherwise use associational
  language.
- Changes after seeing results are allowed but logged in the plan with the
  reason ("added quadratic elevation after DHARMa showed curvature").
