---
name: eco-convergence
description: Convergence and numerical-health checks for MCMC (Nimble, NumPyro) and likelihood fits (unmarked, glmmTMB, lavaan, momentuHMM). Use after every fit before reading any estimate.
---

# eco-convergence

## Bayesian (Nimble / NumPyro)
- >= 3 chains from dispersed inits. Thresholds for reporting: R-hat < 1.01,
  bulk and tail ESS >= 400 for every reported parameter and derived
  quantity. (Recovery replicates may use R-hat < 1.05 to save time; final
  real-data fits may not.)
- If not met: run longer (double iterations), then reparameterise (centre
  covariates, non-centred random effects), then use block samplers
  (Nimble RW_block / AF_slice on correlated parameters). Report what was
  done.
- NumPyro: zero divergences after warmup; if any, raise target_accept to
  0.95 and reparameterise.
- Data augmentation (SCR, N-mixture with M): posterior of psi or N must
  not pile up near the upper bound M; if it does, increase M.
- Prior sensitivity for weakly identified parameters: refit with a wider
  prior; report if conclusions change.

## Likelihood
- Optimiser converged (code 0); Hessian invertible; no NaN SEs.
- Boundary estimates (p or psi at 0/1, variance at 0) mean the parameter is
  not estimable from these data. Report this and simplify.
- Multimodality: refit from 3+ starting values (especially HMMs: use
  momentuHMM's retryFits or random starts); keep the best likelihood and
  report whether the others agreed.
