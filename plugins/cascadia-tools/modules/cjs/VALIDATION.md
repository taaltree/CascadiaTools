# cjs: validation

All recovery runs: `scripts/Rscript core/run.R recover cjs --config dev/cjs/cfg/<label>.yaml --nsim N --cores 16 --out results/cjs/recover_<label>`
on njobvu (SLURM via `scripts/ejob`, 16 cores), module `fast_settings` (marked: ML via TMB;
Nimble: 3 chains x 3,000 iterations, 1,000 burn-in, no GOF). Drivers: `dev/cjs/run_marked.sh`,
`dev/cjs/run_nimble.sh`. Raw results: `results/cjs/recover_<label>/` (summary CSV,
per-replicate CSV, HTML report). The tables below come straight from the summary CSVs
(`dev/cjs/md_tables.R`). The marked runs and the Nimble default run were repeated on 2026-10-20
with the final module code and reproduced the earlier numbers exactly (same seeds).

Gate (`core/R/recovery.R::recovery_gate`, thresholds NOT changed): |bias| <= max(10% of
truth, 3 MC SE) (for truth = 0: max(0.1 x empirical SD, 3 MC SE)); 95% interval coverage
>= 0.95 - 3 MC SE (floor 0.80); converged in >= 90% of ALL replicates. "Converged" means
max R-hat < 1.05 for Nimble. For marked it means a finite Hessian, no probability estimate
outside 0.001-0.999, and no logit coefficient with |estimate| or SE > 10 (boundary or
non-identifiable fits are flagged).

**Default design** (module defaults): 8 annual occasions; 50 newly marked per occasion at
occasions 1-7 (350 animals); half male; survival 0.60 for females at average covariates;
male +0.4 logit; winter severity (occasion covariate) -0.5 per SD; body mass (individual
covariate) +0.3 per SD; constant recapture p = 0.5. Fitted: phi ~ sex + winter_severity + body_mass, p ~ 1.

**Timing** (median per replicate): marked 0.5 s; Nimble marginal (dCJS_vv) 44 s;
Nimble latent state 163 s; Nimble with `(1 | time)` at 12 occasions 292 s. The selftest
(4 reps, marked) takes 3.3 s per replicate including GOF.

## 1. Recovery at the default design: PASS (all three engines)

### marked (default engine), nsim 200: 200/200 converged, coverage bound 0.904

| estimand | truth | mean est | rel. bias | bias MC SE | RMSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|---|
| detection | 0.500 | 0.500 | +0.1% | 0.0026 | 0.036 | 0.940 |  | PASS |
| phi_beta_body_mass | 0.300 | 0.317 | +5.5% | 0.0076 | 0.109 | 0.960 | 0.89 | PASS |
| phi_beta_sex_male | 0.400 | 0.409 | +2.3% | 0.0138 | 0.195 | 0.955 | 0.56 | PASS |
| phi_beta_winter_severity | -0.500 | -0.516 | -3.3% | 0.0126 | 0.178 | 0.980 | 0.89 | PASS |
| survival_female | 0.600 | 0.599 | -0.2% | 0.0027 | 0.039 | 0.945 |  | PASS |
| survival_male | 0.691 | 0.691 | +0.0% | 0.0025 | 0.035 | 0.930 |  | PASS |

### Nimble marginal (nimbleEcology dCJS_vv), nsim 100: 100/100 R-hat < 1.05 (max 1.003), bound 0.885

| estimand | truth | mean est | rel. bias | bias MC SE | RMSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|---|
| detection | 0.500 | 0.496 | -0.9% | 0.0036 | 0.037 | 0.930 |  | PASS |
| phi_beta_body_mass | 0.300 | 0.322 | +7.4% | 0.0115 | 0.116 | 0.940 | 0.88 | PASS |
| phi_beta_sex_male | 0.400 | 0.437 | +9.3% | 0.0192 | 0.195 | 0.960 | 0.60 | PASS |
| phi_beta_winter_severity | -0.500 | -0.539 | -7.9% | 0.0201 | 0.204 | 0.980 | 0.90 | PASS |
| survival_female | 0.600 | 0.602 | +0.3% | 0.0038 | 0.038 | 0.950 |  | PASS |
| survival_male | 0.691 | 0.699 | +1.1% | 0.0037 | 0.037 | 0.900 |  | PASS |

### Nimble latent state (z sampled), nsim 100: 100/100 R-hat < 1.05 (max 1.016), bound 0.885

| estimand | truth | mean est | rel. bias | bias MC SE | RMSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|---|
| detection | 0.500 | 0.495 | -1.0% | 0.0037 | 0.037 | 0.930 |  | PASS |
| phi_beta_body_mass | 0.300 | 0.322 | +7.5% | 0.0114 | 0.116 | 0.940 | 0.88 | PASS |
| phi_beta_sex_male | 0.400 | 0.438 | +9.4% | 0.0193 | 0.195 | 0.960 | 0.60 | PASS |
| phi_beta_winter_severity | -0.500 | -0.540 | -8.1% | 0.0199 | 0.202 | 0.970 | 0.91 | PASS |
| survival_female | 0.600 | 0.602 | +0.3% | 0.0038 | 0.038 | 0.960 |  | PASS |
| survival_male | 0.691 | 0.699 | +1.2% | 0.0037 | 0.037 | 0.910 |  | PASS |

The two Nimble formulations use the same seeds and give nearly identical estimates, which
confirms that the marginal (forward algorithm) and latent-state likelihoods agree. The
posterior-mean slopes are 7-9% further from zero than the ML estimates, the usual
small-sample skew of logit-scale posteriors under N(0, 2) priors. This is within tolerance.

## 2. Model variants (marked unless stated): all PASS

| variant (config) | nsim | converged | result |
|---|---|---|---|
| Fully time-dependent phi(t) p(t) with random year variation in the truth (`timedep`) | 200 | 189/200 (94.5%) | 13/13 PASS: per-interval survival bias -3.3% to +4.3%, coverage 0.937-0.989; last phi x p product -0.1%, coverage 0.968 |
| same, Nimble (`timedep_nimble`) | 100 | 100/100 | 13/13 PASS: bias -1.6% to +5.3%, coverage 0.93-0.99 |
| 40% juveniles (age_class), trap-happy recapture (td, +0.8 logit), sex effect on p (`agetd`) | 200 | 200/200 | 13/13 PASS: juvenile effect -1.033 (truth -1.0, coverage 0.945); td 0.787 (truth 0.8, coverage 0.95); survival in 4 sex x age classes within 1.2% |
| 30% transients, tsm fitted (crossed with the other terms, the default) (`transient_tsm`) | 200 | 200/200 | 8/8 PASS: resident survival -0.9 / -0.7%; proportion transient 0.282 / 0.279 (truth 0.30, coverage 0.97 / 0.955); covariate effects within 4% |
| Random year effects on survival `(1 | time)`, 12 occasions, Nimble (`randomtime_nimble`) | 100 | 100/100 | 5/6: survival, detection and slopes PASS (within 2.2%); **survival_temporal_sd FAILS on bias** (0.569 vs 0.5, +13.8%, MC SE 0.020; coverage 0.94). See Known limits |

The 11 replicates of `timedep` that were not converged were boundary fits (some
interval's survival or p estimated at 0 or 1), which the module flags as not converged.

## 3. Harder design: 5 occasions, 30 newly marked per occasion (120 animals), p = 0.25

### marked (`hard`), nsim 200: 190/200 converged (95%), bound 0.903. All 6 PASS, but see note

| estimand | truth | mean est | rel. bias | bias MC SE | RMSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|---|
| detection | 0.250 | 0.267 | +6.8% | 0.0053 | 0.075 | 0.911 |  | PASS |
| phi_beta_body_mass | 0.300 | 0.423 | +40.9% | 0.0423 | 0.595 | 0.979 | 0.08 | PASS |
| phi_beta_sex_male | 0.400 | 0.477 | +19.2% | 0.0630 | 0.870 | 0.968 | 0.09 | PASS |
| phi_beta_winter_severity | -0.500 | -0.575 | -15.0% | 0.0744 | 1.026 | 0.979 | 0.08 | PASS |
| survival_female | 0.600 | 0.600 | +0.0% | 0.0116 | 0.160 | 0.953 |  | PASS |
| survival_male | 0.691 | 0.694 | +0.4% | 0.0111 | 0.153 | 0.947 |  | PASS |

Survival and p are recovered without bias and with nominal coverage, but the intervals are
wide (RMSE about 0.16). The covariate effects pass only because their MC SE is large. Power is
8-9%, the RMSE (0.6-1.0 logit) is larger than the effects, and the large relative bias comes
from a few replicates with near-separation. **At this design, survival can be estimated but
effects on survival cannot be detected.** The report says so through the wide intervals.

### Nimble (`hard_nimble`), nsim 100: 100/100 R-hat < 1.05. **3/6 FAIL (bias)**

| estimand | truth | mean est | rel. bias | bias MC SE | RMSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|---|
| detection | 0.250 | 0.243 | -2.9% | 0.0065 | 0.066 | 0.920 |  | PASS |
| phi_beta_body_mass | 0.300 | 0.511 | +70.2% | 0.0475 | 0.518 | 0.910 | 0.23 | **FAIL** (bias) |
| phi_beta_sex_male | 0.400 | 0.703 | +75.7% | 0.0889 | 0.935 | 0.910 | 0.14 | **FAIL** (bias) |
| phi_beta_winter_severity | -0.500 | -0.704 | -40.7% | 0.0704 | 0.730 | 0.960 | 0.21 | PASS |
| survival_female | 0.600 | 0.668 | +11.3% | 0.0125 | 0.141 | 0.940 |  | **FAIL** (bias) |
| survival_male | 0.691 | 0.756 | +9.4% | 0.0134 | 0.148 | 0.900 |  | PASS |

With about 120 animals and p = 0.25, the posterior mean is pulled by the priors
(Logistic(0,1) on the intercept, N(0, 2) on slopes) and by the skew of the posterior. Coverage
stays at 0.90-0.96, so the credible intervals remain honest, but point estimates are biased.
This is a documented limit (below). At such designs, use the marked engine, or report
intervals rather than point estimates.

## 4. Negative control: 30% transients ignored (`transient_ignored`). The gate FAILS, as it should

| estimand | truth | mean est | rel. bias | bias MC SE | RMSE | coverage | gate |
|---|---|---|---|---|---|---|---|
| detection | 0.500 | 0.444 | -11.2% | 0.0031 | 0.071 | 0.710 | **FAIL** |
| phi_beta_body_mass | 0.300 | 0.273 | -9.0% | 0.0082 | 0.119 | 0.925 | PASS |
| phi_beta_sex_male | 0.400 | 0.368 | -7.9% | 0.0165 | 0.234 | 0.920 | PASS |
| phi_beta_winter_severity | -0.500 | -0.409 | +18.3% | 0.0138 | 0.215 | 0.900 | **FAIL** |
| survival_female | 0.600 | 0.518 | **-13.7%** | 0.0035 | 0.096 | **0.500** | **FAIL** |
| survival_male | 0.691 | 0.608 | **-12.1%** | 0.0030 | 0.094 | **0.450** | **FAIL** |

When transients (animals that leave after marking) are not modelled, apparent survival is
biased low by 12-14%, and the 95% intervals contain the truth only half of the time. So the
recovery check catches this misspecification. In real data, R2ucare TEST 3.SR detects it: in
the example data set, 3.SR p = 0.027 and the report recommends `tsm`. Fitting `tsm`
(section 2) removes the bias.

A second, related failure: `tsm` fitted **additively** on the logit scale
(`transient_tsm_additive`, `tsm_interact: false`). Survival passes (-0.8% / -2.8%), but
covariate effects are attenuated by 16-25%. The body-mass, sex and winter effects and
proportion_transient_male FAIL. Transients scale first-interval survival multiplicatively,
so an additive logit term mixes them into the covariate slopes. For this reason the module
crosses `tsm` with the other survival terms by default (MARK's phi(a2 x ...)), and that
version passes (section 2).

## 5. False-positive rate with all survival effects = 0 (`null`, marked, nsim 200): PASS

| estimand | truth | mean est | bias | coverage | rejection rate (alpha 0.05) |
|---|---|---|---|---|---|
| phi_beta_body_mass | 0 | 0.001 | +0.001 | 0.920 | **0.080** |
| phi_beta_sex_male | 0 | -0.003 | -0.003 | 0.925 | **0.075** |
| phi_beta_winter_severity | 0 | 0.011 | +0.011 | 0.970 | **0.030** |
| survival_female / male | 0.6 | 0.597 / 0.596 | -0.5% / -0.6% | 0.935 / 0.960 | |

The false-positive rates are 3-8% (MC SE about 1.5%). Individual-level effects are slightly
above nominal (Wald test, 350 animals), and the winter effect (8 years) is conservative.
Treat a lone individual-covariate effect with p just under 0.05 cautiously.

## 6. Cross-checks against independent implementations (`dev/cjs/crosscheck.R`, log `logs/cjs_crosscheck2.log`)

**European dipper** (Lebreton et al. 1992; 294 birds, 7 occasions), read from a long CSV by
the module's own reader. Its histories are identical to `marked::dipper` (all 294).

| model | quantity | module marked | hand-coded forward likelihood | wiqid::survCJS (checked against MARK by its author) | marked::crm direct | module Nimble (posterior mean) |
|---|---|---|---|---|---|---|
| phi(.) p(.) | survival | 0.5602 | (logit 0.2421) | 0.5602 | (logit 0.242) | 0.5615 |
| | recapture | 0.9027 | (logit 2.2263) | 0.9026 | (logit 2.227) | 0.8957 |
| | -2 lnL | 666.84 | 666.838 | 666.838 | 666.838 | |
| phi(flood) p(.) | survival, normal years | 0.6071 | | 0.6071 | | 0.6079 |
| | survival, flood years (1982-83, 1983-84) | 0.4688 | | 0.4688 | | |
| | recapture | 0.8998 | | 0.8998 | | 0.8924 |
| | -2 lnL | 660.10 | 660.1028 | 660.1028 | | |
| phi(sex) p(.) | logit sex effect | 0.07929 | 0.07929 | (no individual covariates in wiqid) | | 0.0804 |

R2ucare on the dipper females gives overall chi-square 10.276 on 12 df (p = 0.592). This is
the value printed in the R2ucare vignette (Gimenez et al. 2018).

**Simulated data, default design (seed 11).** The module marked fit matches the hand-coded
likelihood to 5 decimals for every coefficient and SE (e.g. sex 0.48591, SE 0.2084).
Nimble marginal and latent agree with ML within 0.03 logit units and have posterior SDs
within 3% of the ML SEs. The fit with unequal intervals (one occasion dropped, so a 2-year
gap) agrees with the hand likelihood to 4 decimals (intercept 0.4647, sex 0.2765, mass 0.1550)
and with Nimble within 0.006 logit units.

**phi(t) p(t), simulated (seed 21).** The module marked fit equals wiqid::survCJS for all six
identifiable interval survivals (e.g. 0.7959, 0.7151, 0.6367) and recaptures (0.5277,
0.4983, ...). wiqid's confounded last pair gives phi7 x p7 = 0.8287 x 0.5501 = 0.4559, which
equals the module's `survival_x_detection_2021` = 0.4559. Nimble gives 0.4468.

## 7. Power (`core/run.R power cjs --vary n_marked_per_occasion`, marked, nsim 100 per level, `results/cjs/power/`)

Power at alpha = 0.05 (default effects; 8 occasions, p = 0.5):

| newly marked per occasion (total animals) | body mass (+0.3/SD) | sex (+0.4) | winter severity (-0.5/SD) | survival coverage F / M |
|---|---|---|---|---|
| 10 (70) | 0.25 | 0.11 | 0.21 | 0.98 / 0.95 |
| 20 (140) | 0.51 | 0.22 | 0.43 | 0.97 / 0.95 |
| 40 (280) | 0.80 | 0.41 | 0.77 | 0.96 / 0.94 |
| 80 (560) | 0.99 | 0.78 | 0.98 | 0.94 / 0.94 |
| 160 (1120) | 1.00 | 0.96 | 1.00 | 0.95 / 0.97 |

Coverage holds at every sample size (0.90-0.99 for all estimands). To have 80% power for a
male-female difference of +0.4 logit (about +9 percentage points of survival), plan for
roughly 80-100 newly marked animals per year over 8 years at p = 0.5.

## 8. Example fits (`modules/cjs/example/`)

- `config.yaml`: song sparrow banding station, simulated with known truth (`TRUTH.txt`), dates
  format, transients, sex, winter temperature, body mass. Marked estimates: resident female survival
  0.544 (truth 0.55); male logit effect 0.30 (0.3); winter 0.27 (0.4); mass 0.16 (0.25);
  proportion transient 0.26 / 0.32 (0.25); female p 0.41 (0.45); male p effect 0.63 (0.5).
  All truths are inside the 95% CIs. 3.SR flags transience (p = 0.027).
- `config_nimble.yaml`: same data, Bayesian, with `(1 | time)`. Estimates are within 0.04 of
  marked, R-hat <= 1.001, bulk ESS >= 3,600.
- `config_dipper.yaml`: the real dipper data with phi(flood). The output matches section 6.

## Known limits

1. **Low-information designs, Bayesian engine.** With about 5 occasions, about 120 animals and
   p of about 0.25, Nimble posterior means of survival (+11%) and of slopes (+40-76%) are
   biased (the gate FAILS on bias; coverage is still 0.90-0.96). marked passes at that design,
   but its effect estimates are uninformative (power below 10%). With few occasions or low p,
   prefer marked, and report intervals rather than point estimates.
2. **SD of random year effects** (`(1 | time)`, Nimble): biased upward by about 14% (0.569
   vs 0.5) with 11 intervals under the half-Normal(0, 1) prior. The gate FAILS on bias, while
   coverage is 0.94. Treat the temporal SD as weakly identified with fewer than about 15
   intervals. The other parameters in that model recover well.
3. **Transients must be modelled.** If they are ignored, survival is about 13% too low and the
   gate fails (section 4). The module flags 3.SR. An additive `tsm` attenuates covariate
   effects, so keep the default `tsm_interact: true`.
4. **False-positive rate** for individual covariates is 7.5-8% at nominal 5% (Wald, 350
   animals).
5. **ML boundary fits.** About 5% of replicates in the hard and phi(t)p(t) designs hit a
   boundary (a probability of 0 or 1). These fits are flagged as not converged and not
   interpreted.
6. **Not modelled:** individual heterogeneity in p beyond covariates (no mixtures), tag loss,
   temporary emigration, losses on capture, and multistate movement. Apparent survival is
   always true survival x site fidelity.
7. `proportion_transient` = 1 - phi_first/phi_later is an unconstrained ratio (per unit time when intervals are unequal; not validated with tsm and unequal intervals). Its table
   interval can extend below 0 (no evidence of transients). The report text truncates it to
   0-1 for display.
8. Random time effects are available only with `engine: nimble`. Random effects on p are not
   implemented (`p_temporal_sd` is used only in simulation).

## Independent verification (2026-10-20, separate verifier agent)

**Literature / code check.** Checked line by line against the CJS likelihood (Lebreton et al. 1992),
the latent-state formulation (Kery & Schaub 2012 ch. 7: z[f]=1, z[t] ~ Bern(phi z[t-1]),
y ~ Bern(p z)), the installed `nimbleEcology::dCJS_vv` source (forward algorithm, x[1] must be 1,
probCapture[1] unused: the module's padding with y = 1, phi = 1, p = 1 before first capture adds
log 1 = 0 and is correct), Pradel et al. 1997 (proportion transient = 1 - phi_first/phi_resident)
and Pradel 1993 (td = caught at t-1). Priors: Logistic(0,1) intercepts (= U(0,1)), N(0, 2) slopes,
half-N(0,1) on the random-year SD (non-centred). Unequal intervals enter as phi^dt in both
engines. Simulator and fit share only `cjs_scaling()`/`cjs_prepare()` (used to define truth on the
standardised scale), so the independent simulator below was written without any module code.

**Independent planted truth** (`dev/cjs/verify/verify_sim.R`, own simulator; logs
`logs/cjs_verify.log`, `logs/cjs_verify_nimble.log`; results `results/cjs/verify/`). Design differs
from the builder's: 10 occasions with **unequal intervals** (one 2-unit gap), 35 marked per
occasion at 1-9 (315 animals), 45% male, group levels "F"/"M", survival 0.72 (F), male -0.5,
wing length +0.4/SD (individual), rainfall -0.35/SD (occasion, gamma-distributed), p = 0.35 (F)
with male +0.6. Fitted phi ~ sex + wing + rain, p ~ sex.

| run | nsim | converged | result |
|---|---|---|---|
| main, marked | 200 | 200/200 | **8/8 PASS**: survival F -0.2% (cov 0.96), M -1.4% (0.945); p F -0.3%, M +1.8%; slopes: sex -6.7%, wing +1.0%, rain +5.0%, p sex +7.6%; coverage 0.93-0.98 |
| main, Nimble marginal | 100 | 100/100 (max R-hat 1.003) | **8/8 PASS**: survival -0.5% / -0.5%; slopes -0.4% to +6.9%; coverage 0.93-0.99 |
| null (all effects 0), marked | 200 | 200/200 | 8/8 PASS; **rejection rates 0.04 (p sex), 0.05 (rain), 0.06 (sex), 0.04 (wing)** |
| 35% transients, tsm fitted, marked | 200 | 200/200 | **10/10 PASS**: proportion transient 0.352 / 0.329 (truth 0.35, cov 0.965 / 0.955); resident survival +0.8% / -0.3% |
| negative control: 35% transients ignored | 200 | 200/200 | **FAILS as it should**: survival -10.7% / -14.4%, coverage 0.58 / 0.585; p -13 to -17%; rain effect +24% |
| negative control: unmodelled p heterogeneity (p 0.15 / 0.65 mixture) | 100 | 100/100 | **FAILS**: survival -8.7% / -10.1%, coverage 0.47 / 0.72 |

**Builder's commands rerun** (`dev/cjs/verify/run_marked_verify.sh`, `results/cjs/verify/rerun_*`):
`default`, `null` and `transient_ignored` reproduced sections 1, 4 and 5 exactly (e.g. survival
-0.2% / 0.0%, coverage 0.945 / 0.93; FPR 0.080 / 0.075 / 0.030; transients ignored -13.7% /
-12.1%, coverage 0.50 / 0.45).

**Usability.** The three examples were refitted (`dev/cjs/verify/run_examples_verify.sh`); every
interpretation sentence of the song sparrow report was recomputed from the fit (odds ratios,
CIs, 10th/90th-percentile survival 0.45 / 0.61 and 0.49 / 0.59, interval range) and is correct;
dipper flood survival 0.607 / 0.469 matches section 6. Figures checked by eye (conditional
effects in original units with rug, per-interval survival, odds-ratio forest plot, GOF, captures,
Nimble trace/posterior). Messy-data tests (`dev/cjs/verify/messy_test.R`): NA covariate, NA
covariate year, missing year row, misspelled band id, duplicated id, bad date all gave clear
messages already.

**Changes made** (backups `dev/cjs/verify/*.preverify`; fit code untouched, so no recovery result changes):
- `R/io.R`: a renamed/missing individual-covariate, group or age column was silently dropped
  (check_data reported nothing and the fit failed later); now `read_data` stops and lists the
  available columns. A group value not in `group_levels` (e.g. "f" vs "F") was silently turned
  into NA and reported as "missing"; now it names the bad value and the allowed levels. Missing
  `date_col`/`occasion_col`/`occasion_cov_col` messages list the available columns. Blank ids in
  the captures file are reported with row numbers. `check_data` now checks that every formula
  variable exists (a misspelled covariate in the formula was not caught before fitting).
- `R/report.R`: when survival varies over time only through occasion covariates (no `time`,
  `Time` or `(1 | time)`), the report said "Survival varied between intervals" and the
  `survival_time` caption implied year-specific estimates. Both now say the values are predicted
  from the named covariate only (no separate year effect).

**Additional limit noted.** With `tsm` and an unequal first interval (dt != 1), the first-interval
survival is modelled as (phi_first)^dt, so proportion_transient is per unit time, not per
release; not validated (the transient runs above use equal intervals).

**Verdict: VERIFIED_WITH_LIMITS.** Passes at the default design and at an independent design
with unequal intervals on both engines; controls fail as they should and the FPR is about 0.05;
the documented limits (Bayesian bias at low-information designs, random-year SD bias) stand.
