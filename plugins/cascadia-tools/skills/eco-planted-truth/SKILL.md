---
name: eco-planted-truth
description: Simulate data with known parameters at the researcher's own design and confirm the model recovers them before fitting real data. Use before every real-data fit, when choosing a sample size, and whenever a model or its settings change.
---

# eco-planted-truth — recover a known answer before touching real data

Every CascadiaTools module can simulate its own data (`module$simulate`). Use it.

## Procedure

1. **Match the design.** Read the researcher's data first (dimensions only:
   sites, surveys, years, individuals, traps, loci, fixes per animal). Put
   those into `design:` in a recovery config. Missing surveys, unequal
   effort and covariate distributions matter: if the module's design knobs
   allow it, mimic them.
2. **Choose plausible truths.** Use values from the literature, a pilot,
   or a first fit of the real data rounded to one significant figure. Also
   try a pessimistic scenario (lower detection, weaker effect).
3. **Run recovery** via SLURM, never the login node:
   ```
   "${CLAUDE_PLUGIN_ROOT}/scripts/ejob" -c 16 -m 32 -t 4:00:00 -- "${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/core/run.R" recover <module> --config recovery.yaml --nsim 100 --cores 16 --out <project>/recovery
   ```
4. **Read the gate.** Every estimand that will be interpreted must PASS.
   If one fails, that is a finding: report it and change the design,
   simplify the model, or do not interpret that estimand.
5. **Measure time first** (timing discipline): run `selftest` (4 reps) and
   extrapolate before launching 100+ replicates. If a projection exceeds a
   few hours, reduce MCMC iterations in `settings` only if R-hat stays
   < 1.05 in recovery, or use the module's faster engine.

## Rules

- Same `fit()` and settings for simulated and real data. Recovering with a
  different model than the one used on real data proves nothing.
- Seeds are recorded; any replicate can be reproduced exactly.
- Don't simulate from the model fitted to the real data and present the
  recovery as validation of those same estimates. It checks the method,
  not the answer. Say so.
