# fisheries-stock: validation

Gate rules (core/R/recovery.R, unchanged): |bias| <= max(10% of truth, 3 MC SE);
95% coverage >= 0.95 - 3 MC SE (floor 0.80); >= 90% of replicates converged.

**Overall status: VERIFIED_WITH_LIMITS for stock_recruit (Ricker) and
surplus_production with a contrasting catch history; catch_curve FAILS the
coverage part of the gate at every design tried (bias passes). See "Known limits".**

All numbers below come from scripts that can be rerun:
- builder runs: `dev/fisheries-stock/recover_all.sh` with `dev/fisheries-stock/*.yaml`
  (logs `logs/fisheries-stock_recover_{fast,mcmc}.log`);
- verification: `dev/fisheries-stock/verify/` (logs `logs/fisheries-stock_verify_*.log`,
  results `results/fisheries-stock/verify/`).
Recovery runs use `fast_settings` (Nimble SR 10k iter, SP 30k iter, 3 chains).

## 1. Builder's planted-truth recovery (`core/run.R recover`)

Rows marked (v2) were rerun after the verification change to posterior medians
(section 3); other rows did not involve Nimble and are unaffected.

### Catch curve (true Z = 0.5, year-class SD 0.3, 400 aged fish, logistic selectivity)

| design | nsim | estimator | rel. bias | coverage | gate |
|---|---|---|---|---|---|
| cc_default | 500 | weighted regression (primary) | +0.1% | 0.900 | FAIL (coverage) |
| cc_default | 500 | Chapman-Robson | -2.4% | 0.916 | FAIL (coverage) |
| cc_default | 500 | GLMM | -4.5% | 0.858 | FAIL (coverage) |
| cc_default, seed 77 | 2000 | weighted regression | -0.3% | 0.915 (bound 0.935) | FAIL (coverage) |
| cc_default, seed 77 | 2000 | Chapman-Robson | -2.4% | 0.914 | FAIL |
| cc_default, seed 77 | 2000 | GLMM | -4.4% | 0.865 | FAIL |
| cc_hard (Z=1, SD 0.6, 100 fish) | 500 | WR / CR / GLMM | -1% / -11% / -19% | 0.895 / 0.899 / 0.765 | FAIL; 45/500 fits returned errors |
| neg. control: CR without c-hat, SD 0.5 | 500 | Chapman-Robson | -5.1% | **0.558** | FAIL (as intended) |
| neg. control: recruitment trend 0.15/yr | 300 | all | +25 to +29% | 0.55-0.69 | FAIL (as intended) |

### Stock-recruitment

| design | nsim | result |
|---|---|---|
| sr_default (Ricker LS, 30 yr, sigma 0.4) | 500 | **PASS all 7** (|rel bias| <= 0.9%, coverage 0.950-0.964) |
| sr_bh (Beverton-Holt LS, builder defaults) | 500 | FAIL: beta +35.5% (coverage 0.998), S_eq +11% (cov 0.910), S_MSY +12.7%; ln alpha, U_MSY, steepness, sigma PASS |
| sr_hard (15 yr, sigma 0.6) | 500 | ln alpha, beta, U_MSY, sigma, steepness PASS; S_MSY and S_eq undefined (Inf) in reps with beta_hat <= 0, so mean is Inf and `recovery_report` errors |
| sr_fpr (beta = 0) | 500 | rejection rate of beta = 0: **0.050** |
| neg. control: LS on AR(1) residuals (phi 0.6) | 500 | coverage 0.76-0.88: **FAIL (as intended)** |
| sr_nimble_ar1 (Nimble AR(1), phi 0.6) (v2) | 100 | **PASS all 8** (before v2: recruitment_sd +14.9% FAIL) |
| sr_nimble_iid (Nimble AR(1) on iid data) (v2) | 100 | **PASS all 8** (before v2: recruitment_sd +10.3% FAIL) |
| power for beta, low-contrast spawners | 200/size | power 0.09 (12 yr), 0.19 (16), 0.25 (20), 0.33 (30), 0.43 (40) |

### Surplus production (r 0.4, K 1000, q 0.002, sigma_proc 0.05, sigma_obs 0.15, 30 yr)

| design | nsim | converged (R-hat < 1.05) | result |
|---|---|---|---|
| sp_contrast (v2) | 200 | 189/200 = 0.945 | **PASS all 7** (r -6.2%, K +8.6%, q -8.5%, MSY +0.6%, B/B_MSY +1.3%, F/F_MSY -0.5%; coverage 0.97-0.98). Before v2 (posterior means): K +12.6% FAIL |
| sp_oneway (v2) | 200 | 152/200 = **0.76** | FAIL: convergence; r -34%, K +44%, q -32% (prior-driven); MSY, B/B_MSY, F/F_MSY, depletion within 10% |

Note: `summarise_recovery()` computes `n_ok` from converged replicates only, so the
gate's `converged_ok` column is always TRUE; convergence fractions above were
computed by hand from the replicate files (see "Core issue" below).

### Cross-checks (`dev/fisheries-stock/crosscheck.R`, `logs/fisheries-stock_crosscheck.log`)
- Chapman-Robson Z and raw SE vs FSA::chapmanRobson: identical (max |diff| 0).
- Weighted regression vs hand-written weighted lm: max |diff| 4e-16.
- GLMM Z, glmmTMB vs lme4::glmer: max |diff| 1.3e-4.
- Ricker LS vs Nimble (flat priors): ln alpha within 0.008, U_MSY within 0.01 (5 datasets).
- Surplus production, Nimble vs JAGS running the Meyer & Millar (1999) BUGS model
  with the same priors: posterior medians within 2.4% for r, K, q, MSY, B/B_MSY, F/F_MSY.
- Beverton-Holt Nimble vs LS disagree strongly (posterior beta 2-4x the LS value
  under the U(0, beta_max) prior): the Nimble BH engine is NOT validated.

## 2. Known limits

- **Catch curves under-cover.** With lognormal year-class variability, every
  estimator's nominal 95% interval covers the truth only 84-93% (weighted
  regression), 77-92% (Chapman-Robson with c-hat) or 79-91% (GLMM) of the time.
  HC3/HC2 sandwich SEs do not fix it (section 3). Bias of the weighted regression
  is small (< 1%) at the default design but -5% (Z = 0.8, 250 fish) to -12%
  (no year-class variability, Z = 0.8; zero-catch ages dropped) in harder designs.
  The report says so and tells the reader to treat the interval as ~90%.
- **Catch-curve assumptions.** Recruitment trend (negative control: +29% bias),
  dome-shaped selectivity and ageing error bias Z; Z cannot be split into F and M
  without an external M.
- **Stock-recruitment.** LS ignores autocorrelated residuals (negative control
  fails); use `sr_engine: nimble, ar1: true`. Beverton-Holt beta and S_eq/S_MSY are
  biased (+11 to +35%) when spawners do not reach the flat part of the curve
  (builder default design); BH passes when the spawner range spans ~0.2-6 x 1/beta
  (section 3). Closed-loop (spawners from past recruits) data give beta +26%
  (time-series bias, Walters 1985). Nimble BH is untested. S_MSY/U_MSY assume
  recruits in spawner units. Steepness posterior is long-tailed (Ricker h can exceed 1).
- **Surplus production.** One-way-trip data: r, K and q are prior-driven (bias
  -34% to +93%; posterior/prior SD ratio for r 0.87) and ~25% of fits do not reach
  R-hat < 1.05 at fast settings; the fit flags this. K and q are only borderline
  identified even with contrast (bias -10% to +9%; q fails at -10.1% in one
  design); interpret MSY and relative status (B/B_MSY, F/F_MSY), not K. Process SD
  is prior-dominated (bias +17% at truth 0.05, -42% at truth 0.10, coverage ~1.0)
  and is a nuisance parameter. Catch assumed known; index assumed proportional to
  biomass; symmetric Schaefer production.

## 3. Independent verification (2026-10-19, separate agent)

**What was checked.** Code read line by line against Chapman & Robson (1960) /
Smith et al. (2012) / FSA (S = T/(n+T-1), Hoenig correction, SE (1-e^-Z)/sqrt(n e^-Z),
c-hat inflation), Maceina & Bettoli (1998) weighting (FSA::catchCurve), Millar
(2015) GLMM, Scheuerell (2016) Lambert-W Ricker S_MSY/U_MSY, Beverton-Holt
closed forms, Mace & Doonan steepness (0.2(alpha phi0)^0.8; alpha phi0/(4+alpha phi0)),
and Meyer & Millar (1999) state-space Schaefer equations and precision priors
(Gamma(3.785518, 0.010223), Gamma(1.708603, 0.008613854)). All correct. The
module's simulators and fitters share only `sr_mean`/`sr_refpoints`; the truth in my
own simulator was computed by numerical optimisation instead, and agrees.

**Independent simulator** (`dev/fisheries-stock/verify/vsim.R`, no module code):
catch curve from a year-by-year age-structured population (M + F, lognormal
recruitment, Poisson survey catch with logistic gear selectivity); stock-recruit
with AR(1) residuals (marginal SD) and optional closed-loop spawners; Schaefer
population with an effort-driven exploitation history and the first 3 index years
missing. Different parameters from the builder: Z 0.8, 250 fish; Ricker ln alpha 1,
beta 5e-4, sigma 0.6, 25 yr; BH ln alpha 1.6, sigma 0.35, 30 yr; SP r 0.6, K 5000,
q 1e-4, 25 yr. Fitted with the module's `fit()`, gate = core `recovery_gate()`
(`dev/fisheries-stock/verify/vrun.R`).

| scenario | nsim | result |
|---|---|---|
| cc_v1 (Z 0.8, YC SD 0.4, 250 fish) | 400 | WR bias -3.2%, coverage 0.875; CR -6.0%, 0.900; GLMM -8.9%, 0.812: **FAIL coverage** |
| cc_v_noyc (no YC variability, 300 fish) | 400 | CR -2.1%, 0.965 PASS; GLMM -3.9%, 0.942 PASS; WR -11.7%, 0.795 FAIL |
| cc_v_neg (CR, no c-hat, YC SD 0.5) | 400 | coverage **0.458**: negative control FAILS as intended |
| sr_v_ricker (LS) | 500 | **PASS all 7** (coverage 0.946-0.962) |
| sr_v_bh (LS, wide spawner range) | 500 | **PASS all 7** |
| sr_v_endog (closed-loop spawners) | 500 | beta +26% FAIL (time-series bias, as expected); others PASS |
| sr_v_fpr (beta = 0) | 1000 | **false-positive rate 0.044** (MC SE 0.006) |
| sr_v_neg_ar (LS on AR(1) rho 0.7) | 500 | coverage 0.67-0.82: **FAIL as intended** |
| sr_v_nim_ar (Nimble AR(1), rho 0.7) | 160 | before fix: recruitment_sd +17.7%, steepness +11.2% FAIL. **After fix: PASS all 8** (sigma +3.9%, steepness +1.6%, phi -8.1%, coverage 0.95-0.97) |
| sp_v_contrast | 200 | converged 0.97. After fix: r -9.0%, K +9.4%, MSY -1.3%, B/B_MSY +0.3%, F/F_MSY +2.2% PASS; q -10.1% FAIL (borderline); sigma_proc +17% FAIL (nuisance). Before fix: K +13.5%, q -10.6% and sigma_proc +27% FAIL |
| sp_v_noisy (sigma_proc 0.10, sigma_obs 0.25) | 200 | converged **0.755 FAIL**; B/B_MSY, F/F_MSY, MSY PASS; r -13%, K +15.5%, q -14% FAIL |
| sp_v_oneway | 100 | converged 0.75; r -53%, K +93% (prior-driven; posterior/prior SD ratio 0.87, warning issued); B/B_MSY +9% PASS |

Catch-curve interval alternatives (`cc_hc3.R`, 2000 reps per design): weighted
regression coverage 0.919-0.933 on builder-type designs and 0.841/0.881 on my
simulator (Z 0.8 / Z 0.5); HC3 SE 0.909-0.929 and 0.869/0.891; HC2 worse. No
estimator or SE reached the bound, so the catch-curve under-coverage is reported
as a limit, not patched.

**Usability.** `fit` on all three example configs works. Interpretation sentences
checked number by number against the estimate tables (alpha = exp(ln alpha), the
factor exp(2 sigma), percentages, Kobe probabilities): all correct, hedged, no
causal claims. Figures checked: catch curve, Kobe, S-R curve, biomass trajectory
are correct. Messy-data tests (`verify/messy/test.R`, `test2.R`): wrong column name,
missing model, zero spawners, missing catch, year gap and index <= 0 already gave
clear messages.

**Changes made by the verifier** (originals in `dev/fisheries-stock/verify/backup_pre_verify/`):
1. Nimble engines (stock_recruit nimble, surplus_production) now report the
   **posterior median** as `estimate` (`fs_use_medians` in R/common.R). Means of
   the right-skewed posteriors of recruitment SD, steepness and K were biased +10
   to +18%. Medians are the Meyer & Millar (1999) convention and are what the
   figures already plotted, so tables and figures now agree. Methods text says so.
2. Non-numeric data values (e.g. age "4+", recruits "1,234.5") used to fail with
   "'from' must be a finite number" or produced a silent text column. They now
   stop with the column, rows and values and a hint (`fs_num`). Dropped
   incomplete brood years and missing ages are reported. Duplicated brood years,
   index <= 0 years and year gaps are now named in check_data messages.
3. The catch-curve interpretation claimed the weighted-regression interval "stays
   close to its nominal 95%"; it now gives the measured coverage ranges. README
   updated likewise (and the trend negative control is +29%, not +31%).

**Core issue (not fixed; core not edited).** `summarise_recovery()` sets `n_ok`
to the number of *converged* replicates, so `recovery_gate()`'s
`converged_ok = n_converged / (n_ok + n_failed)` is always 1 and the ">= 90%
converged" rule is never enforced for any module. Also `recovery_report()` errors
when a gate row is NA (Inf estimates, sr_hard).

**Verdict: VERIFIED_WITH_LIMITS** for stock_recruit (Ricker LS and Nimble AR(1);
BH LS only with wide spawner contrast) and surplus_production with contrast
(management quantities). **catch_curve remains OPEN on interval coverage**
(bias passes; intervals ~3-11 points short of 95%).

## 4. Follow-up (2026-10-19, builder continuation)

| design | nsim | result |
|---|---|---|
| sr_bh_wide (BH LS, builder simulator, S_mean 100, S_cv 1.2: S spans ~0.2-8 x 1/beta) | 500 | ln alpha +1.9%, U_MSY +0.5%, S_MSY +2.1%, S_eq +2.2%, steepness +1.2%, sigma -0.4% PASS (coverage 0.946-0.964); **beta +12.7% FAIL on bias** (coverage 0.970; right-skewed LS estimate of beta, so the mean overshoots although intervals cover). Interpret BH reference points, not beta itself. |
| sr_nimble_bh_wide, sr_nimble_bh_default (Nimble BH) | 100 each | submitted (`logs/fisheries-stock_recover_bh.log`); results not yet recorded here. Until they are, the Nimble Beverton-Holt engine remains NOT validated. |

Core issue still open: a minimal fix to `recovery_gate()` (use all returned replicates,
not only converged ones, as the convergence denominator) was attempted but the edit
to core/R/recovery.R was refused by the permission system; core is unchanged.
