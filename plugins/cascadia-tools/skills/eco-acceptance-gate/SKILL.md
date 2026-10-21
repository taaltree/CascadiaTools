---
name: eco-acceptance-gate
description: The standard an ecological analysis must meet before it is called done. Use whenever an estimate, a fitted model, a power result or a new module is about to be reported to a researcher or declared finished.
---

# eco-acceptance-gate — what "done" means for an ecological analysis

**An analysis is done when the method has been shown to recover known truth
at this study's design, the fit to the real data has converged and passes
model checks, and every number reported was produced by a script that can be
rerun. It is never done merely because the model ran without error.**

Ecological models fail silently. An occupancy model with too few surveys
returns psi = 1.00 with a tiny SE. An N-mixture model on data that violate
closure returns abundance 10x too high with good-looking intervals. A
Nimble chain stuck in one mode reports a crisp posterior. An SSF without
enough available steps returns biased selection coefficients. None of
these announce themselves. The researcher, who often cannot read the code,
is relying on the gate.

## The gate (run in order; declare 0 before seeing results)

**0. Declare before fitting** (see eco-analysis-plan). Estimands, model,
covariates, priors, and what would count as an adequate fit, written into
`analysis_plan.md` in the project folder before the real data are fitted.

**1. Planted truth at this design** (see eco-planted-truth). Simulate data
that match the researcher's design (number of sites, surveys, individuals,
traps, loci, detection level) with plausible parameter values, and run
`${CLAUDE_PLUGIN_ROOT}/core/run.R recover` with nsim >= 100 (>= 50 if each fit takes minutes).
`recovery_gate()` must PASS for every estimand that will be interpreted:
|bias| <= max(10% of truth, 3 Monte Carlo SE); 95% coverage >= 0.95 - 3 MC
SE (floor 0.80); >= 90% of replicates converged. Thresholds are fixed in
`${CLAUDE_PLUGIN_ROOT}/core/R/recovery.R`; never widen them to make a case pass.

**2. Convergence on the real data** (see eco-convergence). Bayesian:
R-hat < 1.01 and bulk ESS >= 400 for every reported quantity, no stuck
chains, traces inspected. Likelihood: optimiser converged, Hessian
invertible, no estimates on a boundary (probabilities at 0/1, SE NaN or
huge).

**3. Model checking** (see eco-model-checking). A goodness-of-fit check
appropriate to the model (DHARMa residuals, MacKenzie-Bailey chi-square,
posterior predictive checks, d-separation for SEM). A failed check is
reported, not hidden; the model is revised or the limitation stated.

**4. Independent route for anything new.** A new module, or a custom model
written for one project, is cross-checked against an established
implementation (unmarked, spOccupancy, secr, lavaan, hierfstat...) or a
second engine (Nimble vs NumPyro, likelihood vs Bayesian) on the same data.
Estimates must agree within Monte Carlo/numerical error, stated as a
number.

**5. Negative control.** The check must be able to fail. For recovery:
fitting a deliberately mis-specified model (ignore detection, ignore
spatial autocorrelation, omit the confounder) must fail the gate. For
power: an effect of zero must give rejection rate near 0.05.

**6. Reproducible.** The report states the command and config that produced
it; data, config and outputs are saved in the project folder.

**7. Verdict.** DONE only if every step passed; otherwise OPEN with what
failed and what it means for interpretation. Say this plainly to the
researcher: "The model cannot separate occupancy from detection with 2
surveys per site; recovery showed psi biased by +25%. Add surveys or treat
psi as an index."

## Failure modes this gate exists to catch

- **Running without error mistaken for success.** Check, don't assume.
- **Estimates on a boundary.** psi = 1.000 (SE 0.000) means non-identified,
  not "present everywhere".
- **Recovery at a convenient design, not the user's.** A model that works
  with 200 sites and 5 surveys may fail at the user's 40 sites and 3 surveys.
  Simulate at their design.
- **Threshold drift.** Moving a threshold until the case passes.
- **Adjectives instead of numbers.** "Recovered well" -> "relative bias
  -2.1% (MC SE 1.5%), coverage 0.94 over 200 replicates".
- **Single-chain confidence.** One chain cannot show non-convergence.
- **Causal language from correlational models** (see eco-interpretation).
