---
name: eco-model-checking
description: Goodness-of-fit and assumption checks for ecological models (GLMs, occupancy, N-mixture, SCR, HMMs, SSFs, SEM, popgen). Use after convergence and before interpretation.
---

# eco-model-checking

| Model | Check | Failing means |
|---|---|---|
| GLM/GLMM | DHARMa simulated residuals: uniformity, dispersion, zero-inflation, residuals vs each covariate | wrong family, missing nonlinearity or covariate |
| Occupancy | MacKenzie-Bailey chi-square parametric bootstrap (c-hat); Bayesian PPC on detection-frequency | heterogeneity in detection; inflate SEs by sqrt(c-hat) if 1 < c-hat < ~3, otherwise revise |
| N-mixture | PPC/chi-square, compare Poisson vs NB vs ZIP; check that abundance isn't at extreme values; sensitivity to K | non-closure, double counting, unmodelled heterogeneity; abundance not trustworthy as absolute |
| SCR | Bayesian p-values on individual and trap encounter frequencies; sigma vs trap spacing (spacing < 2 sigma); M not limiting | behavioural response, transience, mask too small |
| Movement HMM | pseudo-residuals (QQ, ACF); state-dependent distributions overlay | missing states or covariates, irregular sampling |
| SSF | number of available steps (>= 10 per used), k-fold CV (UHC plots) | habitat availability misrepresented |
| SEM | Fisher's C / d-sep tests (piecewiseSEM), chi-square/RMSEA/CFI (lavaan) | missing paths; do not interpret a rejected SEM |
| Popgen | HWE per locus and population, linkage disequilibrium, null alleles, sample size per pop | Wahlund effect, null alleles, family structure |
| Fisheries | residual plots of growth/S-R fits, catch-curve linearity over the ages used | wrong age range, gear selectivity, recruitment variation |

Report each check with its statistic and a sentence on what it means.
A failed check is information; never re-run with different seeds until it
passes.
