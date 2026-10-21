# abund-interactions: validation

All recovery runs: `scripts/Rscript core/run.R recover abund-interactions --config dev/abund-interactions/<cfg>.yaml --nsim 100 --cores 16 --out results/abund-interactions/recover_<cfg>`
(configs in `dev/abund-interactions/`; tables printed by `dev/abund-interactions/val_tables.R` straight from
the saved `*_recovery.rds`, so no number below was typed by hand). Runs were made 2026-10-18/19 with R 4.5.3,
nimble 1.4.3 and unmarked (ecotools `Rlib/`), on njobvu with 16 workers.

Recovery uses `fast_settings`: Nimble 3 chains x 5000 iterations, 1500 burn-in, thin 2 (5250 draws), no GOF and no
two-step side fit. Replicate r uses seed 100000 + r. The Nimble, ML and two-step runs of a design therefore
fit **identical data sets** (paired comparisons below).

Gate (core `recovery_gate`, unchanged): |bias| <= max(10% of |truth|, 3 MC SE), or <= max(0.1 x empirical SD, 3 MC SE)
when truth = 0; coverage of the 95% interval >= 0.95 - 3 MC SE (0.885 at nsim 100); >= 90% of replicates
converged. Nimble "converged" means max R-hat < 1.05 over all monitored quantities. The real-data
standard printed in reports (R-hat < 1.01, bulk ESS >= 400) is checked separately on every fit.

**Default design** (module defaults): 600 camera sites; 4 occasions of 7 camera-days each, with 10% of
occasions losing a random number of days. Dominant species: mean 1 individual per site at average
woodland, log-slope on woodland 0.5, per-individual daily detection r = 0.10 (weekly ~0.52 for one
individual). Subordinate: occupancy 0.6 where the dominant is absent, logit slopes woodland 0.5 and
broadleaf 0.3, **interaction gamma0 = -0.7 per dominant individual at average broadleaf share** and
**gamma1 = +0.5 per SD of broadleaf share** (the interaction weakens, then reverses, in broadleaf-rich
woodland, as in Twining et al. 2022 Proc R Soc B), daily detection 0.15. Simulated with the model in
README.md (Royle-Nichols with K camera-days per occasion), as in Twining et al. (2025) simulation I.

## Verdict

| design / check | engine | result |
|---|---|---|
| **default (600 sites)** | nimble (joint) | **PASS** all 12 estimands; 100/100 R-hat < 1.05 (max 1.012) |
| default | ml (exact marginal likelihood, independent code) | **PASS** all 12 |
| **harder: low dominant detection** (r = 0.03/day, ~0.19/week) | nimble | **PASS** all 12; 99/100 converged |
| **false-positive rate** (gamma0 = gamma1 = 0) | nimble | **PASS**; P(95% CrI excludes 0) = **0.08** (gamma0), **0.03** (gamma1) |
| **negative control**: occupancy-mediated model (psi depends on 1[N_A > 0]) fitted to abundance-mediated data | ml | **FAIL (as required)**: gamma0 -49%, gamma1 +52%, psi_beta_woodland -27%; coverage 0.64 / 0.82 / 0.70 |
| comparison: two-step plug-in (occuRN -> E[N\|y] -> occu) | twostep | **FAIL**: gamma0 attenuated +11.6% (coverage 0.84), gamma1 -16.9% (0.87); intervals 13-14% too narrow |
| smaller: 300 sites | nimble | FAIL gamma1 (+19.4% bias) and r coverage 0.86 (others pass) |
| smaller: 300 sites | ml | FAIL gamma1 (+11.3% bias) |
| cross-check vs unmarked occuRN + occu-single module (one data set) | nimble vs ml vs unmarked | agree within 0.05 SD-units on all shared parameters (section 5) |
| power over n_sites (150 / 300 / 600; nsim 50) | nimble | power for gamma1 0.70 / 0.88 / 1.00 (section 6) |

The module is **DONE at the default design** (600 sites, 4 weekly occasions, daily records). Smaller
designs fail the bias rule for the habitat-by-interaction term. That is reported as a limit, and no
threshold was changed.

Honesty note on the default design: the first full runs used 300 sites as the default (folders
`results/abund-interactions/n300_recover_*`). There gamma1 was biased +19% (Nimble) and +11% (ML on the same
data), so the joint model failed the gate. The ML engine shows the same direction of bias, so it is
finite-sample bias of the estimator at that design, not an MCMC problem. The default was raised to 600
sites, the smallest design in the power run with |bias| < 10% for both interaction terms. Twining et al.
(2025) also report that a habitat-varying interaction needs > 600 sites when occasions are single 0/1
records. The 300-site results stay in this file.

## 1. Recovery at the default design (600 sites)

### nimble, joint model (nsim 100, 0 errors, 100/100 R-hat < 1.05, max R-hat 1.012, median 349 s per replicate on 1 core)

| param | truth | mean est | rel. bias | bias MCSE | coverage (bound) | CI width | mean post. SD | emp. SD | power | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| dominant_detection_per_individual |  0.1 | 0.101 | +0.6% | 0.000496 | 0.93 (0.885) | 0.0188 | 0.0048 | 0.00496 | - | PASS |
| dominant_mean_abundance |    1 | 0.997 | -0.3% | 0.00561 | 0.97 (0.885) | 0.236 | 0.0603 | 0.0561 | - | PASS |
| dominant_prop_occupied | 0.632 | 0.632 | -0.0% | 0.000668 | 0.96 (0.885) | 0.0267 | 0.00688 | 0.0189 | - | PASS |
| dominant_total_abundance |  678 |  676 | -0.2% | 3.16 | 0.95 (0.885) |  117 |   30 | 35.7 | - | PASS |
| lambda_beta_woodland |  0.5 | 0.492 | -1.6% | 0.00435 | 0.98 (0.885) | 0.173 | 0.0442 | 0.0435 |    1 | PASS |
| psi_beta_broadleaf |  0.3 | 0.29 | -3.3% | 0.0126 | 0.98 (0.885) | 0.547 | 0.139 | 0.126 | 0.61 | PASS |
| **psi_beta_dominant_abundance** | -0.7 | -0.709 | -1.3% | 0.0134 | 0.96 (0.885) | 0.51 | 0.13 | 0.134 |    1 | PASS |
| **psi_beta_dominant_abundance_x_broadleaf** |  0.5 | 0.519 | +3.8% | 0.0134 | 0.95 (0.885) | 0.511 | 0.131 | 0.134 |    1 | PASS |
| psi_beta_woodland |  0.5 | 0.514 | +2.9% | 0.0132 | 0.92 (0.885) | 0.452 | 0.115 | 0.132 | 0.98 | PASS |
| subordinate_detection | 0.15 | 0.15 | -0.2% | 0.000435 | 0.94 (0.885) | 0.0173 | 0.00442 | 0.00435 | - | PASS |
| subordinate_occupancy_no_dominant |  0.6 | 0.601 | +0.2% | 0.00379 | 0.94 (0.885) | 0.137 | 0.0349 | 0.0379 | - | PASS |
| subordinate_prop_occupied | 0.443 | 0.443 | +0.0% | 0.000358 | 0.96 (0.885) | 0.0137 | 0.0036 | 0.0215 | - | PASS |

Finite-sample quantities (`dominant_total_abundance`, `*_prop_occupied`) are compared with the realised
latent states of each data set. Their posterior SD is (correctly) smaller than the between-data-set SD.

### ml engine, same 100 data sets (exact marginal likelihood: N_i in 0..60 and z_i summed out; Wald intervals; median 62 s/rep)

| param | truth | mean est | rel. bias | bias MCSE | coverage (bound) | CI width | mean SE | emp. SD | power | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| dominant_detection_per_individual |  0.1 | 0.0993 | -0.7% | 0.000502 | 0.93 (0.885) | 0.0189 | 0.00497 | 0.00502 | - | PASS |
| dominant_mean_abundance |    1 | 1.01 | +0.8% | 0.00576 | 0.97 (0.885) | 0.241 | 0.0634 | 0.0576 | - | PASS |
| dominant_prop_occupied | 0.632 | 0.633 | +0.1% | 0.000681 | 0.96 (0.885) | 0.0276 | 0.00718 | 0.019 | - | PASS |
| dominant_total_abundance |  678 |  684 | +1.0% | 3.26 | 0.94 (0.885) |  124 | 31.8 | 36.9 | - | PASS |
| lambda_beta_woodland |  0.5 | 0.494 | -1.2% | 0.00437 | 0.98 (0.885) | 0.173 | 0.0441 | 0.0437 |    1 | PASS |
| psi_beta_broadleaf |  0.3 | 0.29 | -3.2% | 0.0125 | 0.98 (0.885) | 0.542 | 0.138 | 0.125 |  0.6 | PASS |
| psi_beta_dominant_abundance | -0.7 | -0.691 | +1.2% | 0.0131 | 0.95 (0.885) | 0.503 | 0.128 | 0.131 |    1 | PASS |
| psi_beta_dominant_abundance_x_broadleaf |  0.5 | 0.498 | -0.4% | 0.0132 | 0.95 (0.885) | 0.501 | 0.128 | 0.132 |    1 | PASS |
| psi_beta_woodland |  0.5 | 0.511 | +2.2% | 0.0132 | 0.93 (0.885) | 0.451 | 0.115 | 0.132 | 0.98 | PASS |
| subordinate_detection | 0.15 | 0.15 | -0.2% | 0.000434 | 0.94 (0.885) | 0.0173 | 0.00456 | 0.00434 | - | PASS |
| subordinate_occupancy_no_dominant |  0.6 | 0.601 | +0.1% | 0.00381 | 0.94 (0.885) | 0.137 | 0.0345 | 0.0381 | - | PASS |
| subordinate_prop_occupied | 0.443 | 0.443 | +0.0% | 0.000362 | 0.94 (0.885) | 0.0137 | 0.00362 | 0.0215 | - | PASS |

## 2. Harder design: low detection of the dominant (r = 0.03 per individual per day)

Config `dev/abund-interactions/hard_lowdet.yaml`; everything else at the default. One marten is then
detected in a given week with probability ~0.19 instead of ~0.52, so N_i is much less certain.

nimble (nsim 100, 0 errors, 99/100 R-hat < 1.05; the exception had max R-hat 1.058; median R-hat_max 1.011; median 325 s/rep)

| param | truth | mean est | rel. bias | bias MCSE | coverage (bound) | CI width | mean post. SD | emp. SD | power | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| dominant_detection_per_individual | 0.03 | 0.0313 | +4.3% | 0.000348 | 0.939 (0.884) | 0.0136 | 0.00347 | 0.00346 | - | PASS |
| dominant_mean_abundance |    1 | 0.982 | -1.8% | 0.0122 | 0.919 (0.884) | 0.446 | 0.114 | 0.121 | - | PASS |
| dominant_prop_occupied | 0.632 | 0.623 | -1.4% | 0.00327 | 0.929 (0.884) | 0.122 | 0.0312 | 0.0393 | - | PASS |
| dominant_total_abundance |  678 |  666 | -1.7% | 7.99 | 0.919 (0.884) |  282 | 71.9 |   82 | - | PASS |
| lambda_beta_woodland |  0.5 | 0.491 | -1.8% | 0.00569 | 0.97 (0.884) | 0.217 | 0.0555 | 0.0566 |    1 | PASS |
| psi_beta_broadleaf |  0.3 | 0.294 | -1.9% | 0.0148 | 0.97 (0.884) | 0.614 | 0.157 | 0.147 | 0.556 | PASS |
| psi_beta_dominant_abundance | -0.7 | -0.714 | -2.0% | 0.0159 |    1 (0.884) | 0.663 | 0.169 | 0.158 |    1 | PASS |
| psi_beta_dominant_abundance_x_broadleaf |  0.5 | 0.527 | +5.3% | 0.0153 | 0.96 (0.884) | 0.632 | 0.161 | 0.152 | 0.96 | PASS |
| psi_beta_woodland |  0.5 | 0.511 | +2.2% | 0.0159 | 0.889 (0.884) | 0.512 | 0.131 | 0.158 | 0.98 | PASS |
| subordinate_detection | 0.15 | 0.15 | -0.2% | 0.000431 | 0.929 (0.884) | 0.0173 | 0.00442 | 0.00428 | - | PASS |
| subordinate_occupancy_no_dominant |  0.6 | 0.599 | -0.2% | 0.00438 | 0.96 (0.884) | 0.17 | 0.0434 | 0.0436 | - | PASS |
| subordinate_prop_occupied | 0.443 | 0.443 | +0.0% | 0.000359 | 0.96 (0.884) | 0.0137 | 0.00359 | 0.0216 | - | PASS |

What low detection costs: the CrI for the dominant's total abundance is 2.4x wider (282 vs 117). The
interaction intervals are 24-30% wider (gamma0 0.663 vs 0.51; gamma1 0.632 vs 0.511). Power for gamma1 falls
from 1.00 to 0.96. The joint model **widens** its intervals (gamma0 coverage 1.00) instead of becoming biased.
psi_beta_woodland coverage, 0.889, sits just above its bound of 0.884.

## 3. False-positive rate: no interaction (gamma0 = gamma1 = 0)

Config `dev/abund-interactions/null_effect.yaml`. nimble (nsim 100, 0 errors, 100/100 converged, median 350 s/rep)

| param | truth | mean est | bias | bias MCSE | coverage (bound) | CI width | power (= FPR) | gate |
|---|---|---|---|---|---|---|---|---|
| psi_beta_dominant_abundance | 0 | 0.011 | +0.011 | 0.0102 | 0.92 (0.885) | 0.394 | **0.08** | PASS |
| psi_beta_dominant_abundance_x_broadleaf | 0 | 0.0137 | +0.014 | 0.00895 | 0.97 (0.885) | 0.36 | **0.03** | PASS |

All other estimands pass as well (full table in `results/abund-interactions/recover_null_effect/`). The
FPR of 0.08 for gamma0 lies within Monte Carlo error of 0.05: the binomial SE at nsim 100 is 0.022, so the 95% range is 0.01-0.09.
Pooled over both terms, 11 of 200 tests rejected (0.055).

## 4. Negative control and the two-step shortcut (the gate CAN fail)

### 4a. Occupancy-mediated model fitted to abundance-mediated data (Twining et al. 2025, their central result)

Config `occmediated.yaml`: `engine: ml, mediator: occupancy`. The subordinate's occupancy then depends
only on whether the dominant is present (1[N_i > 0]), which is the assumption behind Waddle et al. (2010)
and Rota et al. (2016) style interaction models. The data are the same 100 default data sets.

| param | truth | mean est | rel. bias | coverage (bound) | gate |
|---|---|---|---|---|---|
| psi_beta_dominant_abundance | -0.7 | -1.05 | -49.4% | 0.64 (0.885) | **FAIL** (bias, coverage) |
| psi_beta_dominant_abundance_x_broadleaf | 0.5 | 0.761 | +52.2% | 0.82 (0.885) | **FAIL** (bias, coverage) |
| psi_beta_woodland | 0.5 | 0.364 | -27.2% | 0.70 (0.885) | **FAIL** (bias, coverage) |
| psi_beta_broadleaf | 0.3 | 0.277 | -7.7% | 0.96 | PASS |
| dominant / detection parameters | | | <= 2.2% | >= 0.93 | PASS |

The gate fails, as required. A presence-only interaction misstates the per-individual effect: the
presence effect, -1.05, absorbs the effect of 1.6 individuals on average where martens are present. It
also pulls the subordinate's woodland response towards zero, because woodland drives marten abundance.
The same control at 300 sites fails in the same way (gamma0 -48.9%, gamma1 +67.1%, woodland -30.5%).

### 4b. Joint model vs two-step plug-in: does the joint model propagate the uncertainty in N_A?

The `twostep` engine fits unmarked::occuRN to the dominant and takes the empirical-Bayes posterior mean
E[N_i | y]. It then fits unmarked::occu to the subordinate with that as a fixed covariate. This is the
common shortcut. It runs on the same 100 default data sets as section 1 (median 4 s/rep):

| param | joint nimble: mean est / coverage / CI width | two-step: mean est / coverage / CI width | two-step gate |
|---|---|---|---|
| psi_beta_dominant_abundance (truth -0.7) | -0.709 / 0.96 / 0.510 | **-0.619** / **0.84** / 0.443 | FAIL (bias +11.6%, coverage) |
| psi_beta_dominant_abundance_x_broadleaf (truth 0.5) | 0.519 / 0.95 / 0.511 | **0.416** / **0.87** / 0.437 | FAIL (bias -16.9%, coverage) |
| psi_beta_woodland (0.5) | 0.514 / 0.92 / 0.452 | 0.506 / 0.93 / 0.446 | PASS |

Paired on identical data (n = 100): the two-step interval is narrower in the median data set, by 0.869x
the joint width for gamma0 and 0.858x for gamma1. The two-step point estimates correlate 0.92-0.94 with
the joint ones but are shrunk towards zero. This is classical errors-in-variables attenuation: E[N|y]
is a noisy, shrunk proxy for N. Ignoring that noise makes the shortcut both **biased towards "no
interaction"** and **overconfident**. The joint model gets both right. At 300 sites the two-step gamma0 is
again attenuated (+10.6%, coverage 0.87). Every report shows the two-step fit beside the joint fit
(`compare_twostep: TRUE`, figure `joint_vs_twostep`) but never uses it for inference.

## 5. Cross-checks against independent implementations (single data sets)

Script `dev/abund-interactions/crosscheck.R`, log `logs/abund-interactions_crosscheck.log`, tables in
`results/abund-interactions/crosscheck/`. It uses K = 1 (6 plain 0/1 occasions, 300 sites), so unmarked's
occuRN and occu fit exactly the same observation model.

**A. Interaction present** (seed 7; r = 0.3, p = 0.4). Columns give the estimate (95% interval). Nimble:
R-hat max 1.001, min ESS 2602; posterior predictive Bayesian p = 0.53 for the dominant and 0.39 for the
subordinate.

| param | truth | nimble joint | ml (own code) | two-step | unmarked occuRN |
|---|---|---|---|---|---|
| dominant_mean_abundance | 1 | 1.10 (0.91, 1.32) | 1.12 (0.92, 1.35) | 1.12 | 1.12 (0.92, 1.35) |
| lambda_beta_woodland | 0.5 | 0.469 (0.339, 0.597) | 0.473 (0.344, 0.601) | 0.465 | 0.465 (0.336, 0.594) |
| dominant_detection_per_individual | 0.3 | 0.282 (0.239, 0.326) | 0.277 (0.235, 0.323) | 0.278 | 0.278 (0.235, 0.325) |
| psi_beta_dominant_abundance | -0.7 | -0.662 (-1.08, -0.30) | -0.622 (-1.00, -0.25) | -0.613 (-0.95, -0.28) | - |
| psi_beta_dominant_abundance_x_broadleaf | 0.5 | 0.604 (0.196, 1.09) | 0.542 (0.112, 0.972) | 0.502 (0.132, 0.871) | - |
| interaction_at_low_broadleaf | | -1.39 (-2.21, -0.70) | -1.28 (-2.01, -0.55) | -1.22 (-1.84, -0.60) | - |
| interaction_at_high_broadleaf | | 0.104 (-0.48, 0.70) | 0.066 (-0.49, 0.62) | 0.023 (-0.48, 0.53) | - |
| psi_beta_woodland | 0.5 | 0.694 | 0.672 | 0.673 | - |
| subordinate_detection | 0.4 | 0.400 | 0.401 | 0.401 | - |
| subordinate_occupancy_no_dominant | 0.6 | 0.580 | 0.576 | 0.583 | - |

The dominant-species part of the ML engine reproduces unmarked::occuRN to 2-3 significant figures. Nimble
posterior means differ from the MLEs by at most 0.06 on the logit scale (< 0.3 posterior SD), as expected for
skewed posteriors under weak priors. On this one data set the two-step interval for gamma1 is 0.74 wide,
against 0.89 for the joint model.

**B. Interaction = 0: the reduced model equals single-species occupancy** (seed 8; data simulated with
gamma = 0). Nimble with `interaction: false` (gamma fixed at 0) vs the **occu-single** module (unmarked
engine) on the subordinate's detection histories, and vs unmarked::occuRN on the dominant's:

| param | truth | reduced nimble (95% CrI) | occu-single / occuRN (95% CI) | full joint (gamma free) |
|---|---|---|---|---|
| subordinate_occupancy_no_dominant | 0.6 | 0.586 (0.520, 0.651) | 0.584 (0.518, 0.648) | 0.556 |
| psi_beta_woodland | 0.5 | 0.894 (0.595, 1.23) | 0.880 (0.573, 1.19) | 0.849 |
| psi_beta_broadleaf | 0.3 | 0.315 (0.046, 0.601) | 0.309 (0.036, 0.581) | 0.308 |
| subordinate_detection | 0.4 | 0.423 (0.391, 0.456) | 0.424 (0.391, 0.458) | 0.423 |
| dominant_mean_abundance | 1 | 1.05 | 1.06 | 1.04 |
| lambda_beta_woodland | 0.5 | 0.425 | 0.427 | 0.426 |
| dominant_detection_per_individual | 0.3 | 0.296 | 0.293 | 0.298 |

The reduced model matches occu-single and occuRN to within 0.015 on every parameter. The full joint model fitted
to the same null data estimates gamma0 = 0.116 (95% CrI -0.197, 0.440), which correctly includes 0.

## 6. Power over number of sites (`core/run.R power --vary n_sites=150,300,600 --nsim 50`)

Output: `results/abund-interactions/power_sites/` (Nimble, fast settings, all other values at the default).

| n_sites | param | rel. bias | coverage | CI width | power | converged |
|---|---|---|---|---|---|---|
| 150 | psi_beta_dominant_abundance | -18.2% | 0.94 | 1.10 | 0.92 | 50/50 |
| 150 | psi_beta_dominant_abundance_x_broadleaf | +23.1% | 0.98 | 1.13 | 0.70 | 50/50 |
| 150 | dominant_mean_abundance | +0.6% | 0.96 | 0.477 | - | 50/50 |
| 150 | psi_beta_broadleaf | +2.7% | 1.00 | 1.16 | 0.18 | 50/50 |
| 300 | psi_beta_dominant_abundance | -13.4% | 0.94 | 0.762 | 1.00 | 50/50 |
| 300 | psi_beta_dominant_abundance_x_broadleaf | +15.9% | 0.94 | 0.773 | 0.88 | 50/50 |
| 300 | dominant_mean_abundance | -2.0% | 0.96 | 0.327 | - | 50/50 |
| 300 | psi_beta_broadleaf | +5.0% | 0.92 | 0.796 | 0.36 | 50/50 |
| 600 | psi_beta_dominant_abundance | -4.0% | 0.98 | 0.513 | 1.00 | 50/50 |
| 600 | psi_beta_dominant_abundance_x_broadleaf | +9.5% | 0.96 | 0.523 | 1.00 | 50/50 |
| 600 | dominant_mean_abundance | -1.9% | 0.92 | 0.232 | - | 50/50 |
| 600 | psi_beta_broadleaf | -7.4% | 0.96 | 0.545 | 0.56 | 50/50 |

Reading: with gamma0 = -0.7 per individual, the constant part of the interaction is detected (power >= 0.92)
even at 150 sites. The habitat-dependence (gamma1 = 0.5 per SD) needs about 300 sites for 80% power and 600
for 100%. Below 600 sites both interaction terms are biased away from zero, by 13-23%: magnitudes are
overstated, an effect that falls as n grows (the usual small-sample inflation of logistic coefficients,
amplified by the latent covariate). Coverage stays >= 0.94 throughout, so intervals remain honest even
where the point estimate is biased. With nsim 50, power and coverage have MC SE of about 0.03-0.07 and
bias about 3-5%.

## 7. Example data set

`example/` is a simulated pine marten - red squirrel camera survey in the formats camera-trap software
exports: `deployments.csv` (with gaps and re-deployments), `records.csv` (timestamped detections of
several species, including non-target ones) and `sites.csv`. It has 260 sites, up to 28 camera-days
each, and truth in `example/TRUTH.md`. `core/run.R fit` on `example/config.yaml` (default settings: 3 x 20 000 iterations, GOF, two-step comparison)
took 7 min 41 s on 4 cores. It gives `results/abund-interactions/example/abund-interactions_report.html`:
max R-hat 1.001, min bulk ESS 3279, posterior predictive p 0.51 (marten) and 0.32 (squirrel). Estimates
against the simulation truth:

| quantity | truth | estimate (95% CrI) |
|---|---|---|
| gamma0, logit per marten at mean broadleaf | -0.6 | -0.53 (-0.93, -0.17) |
| gamma1, per SD broadleaf share | 0.6 | 0.69 (0.31, 1.14) |
| marten daily detection per individual | 0.08 | 0.084 (0.071, 0.098) |
| marten log-slope on woodland | 0.5 | 0.53 (0.39, 0.67) |
| squirrel occupancy with no marten | 0.6 | 0.62 (0.53, 0.72) |
| total marten-site abundance | 274 | 256 (224, 298) |
| proportion of sites with squirrel | 0.504 | 0.513 (0.500, 0.527) |

I read the report and every interpretation sentence matches the table. For example, odds ratio per marten
exp(-0.53) = 0.59. At the 10th / 90th percentile of broadleaf share the odds ratio is 0.25 / 1.42, the
second with a CrI overlapping 1, worded "not clearly different from zero". The two-step interval is 86% as
wide as the joint one. The figures render on the response scale with 50% and 95% posterior bands:
squirrel occupancy vs marten abundance, interaction curves at the 10th/50th/90th percentile of
broadleaf, interaction strength across broadleaf, abundance and occupancy vs woodland, abundance map,
coefficients, joint vs two-step, GOF and trace. The example has 260 sites, below the validated 600, so
the report carries the SMALL DESIGN note.

## Known limits

1. **Fewer than ~600 sites with a habitat-varying interaction.** At 300 sites gamma1 is biased +19%
   (Nimble) / +11% (ML on the same data), and the gate fails. At 150 sites the biases are -18% (gamma0)
   and +23% (gamma1). Interval coverage stays near nominal, so CrIs can still be read. The point
   estimate overstates the strength of the interaction. A constant interaction (no habitat term) is
   better determined, but it was not separately gate-tested below 600 sites.
2. **Dominant per-individual detection at the default (r = 0.10) at 300 sites**: coverage of r was 0.86
   (bound 0.885). At 600 sites it is 0.93.
3. **Low dominant detection (r = 0.03/day)** passes at 600 sites. Lower values and fewer sites were not
   tested. Expect interactions to be weakly identified there, and look at the width of the CrI.
4. **Royle-Nichols assumptions** (Poisson N, independent and identical detection of individuals,
   closure) were not stress-tested. Overdispersed or group-living dominants will bias N_i and so gamma.
   The posterior predictive check is the only guard, and its power was not measured.
5. **Direction and causality**: the model assumes A affects S's occupancy and that S's detection does not
   depend on N_A. Shared unmodelled habitat drivers load onto gamma. Reports state "associated with",
   not "caused by".
6. **Two-step engine** fails the gate (attenuated, overconfident). It is shown for comparison only.
7. **Run time**: about 6 min per Nimble fit at 600 sites with fast settings (1 core per replicate); the full
   default settings took 7.7 min at 260 sites (4 cores, incl. GOF and two-step). Scale roughly linearly in sites.
8. Not implemented: abundance (N-mixture) model for the subordinate; time-varying (dynamic)
   interactions; more than two species; spatially explicit N (SCR) for the dominant; covariates on
   the subordinate's detection that involve N_A.

## Independent verification (2026-10-19, separate verifier agent)

Done by an agent that did not build the module. Scripts are in `dev/abund-interactions/verify/`, saved results in
`results/abund-interactions/verify/`, and logs in `logs/abund-interactions_verify*.log`. The gate is core `recovery_gate` with no
threshold changes.

**Literature.** The model was checked line by line against Twining, Augustine, Royle & Fuller (2025, Ecology 106:e4468; PMC11725697):
- N_i^D ~ Poisson(lambda_i), with log(lambda_i) linear.
- y^D_ij ~ Binomial(K, 1 - (1 - r_ij)^N_i).
- z^S_i ~ Bernoulli(psi_i), with logit(psi_i) = beta X + gamma0 N_i + gamma1 N_i x cov_i.
- y^S_ij ~ Bernoulli(p z_i). The module generalises this to Binomial(K_ij, z_i p_ij) over camera-days.

The Nimble code (`abi_nimble_code`), the exact-likelihood engine (`abi_site_loglik`) and the simulator all implement this model.
- gamma0 is the slope at the standardised modifier = 0, which is average habitat.
- The paper's sample-size findings match the README. A constant interaction was unbiased at 300 sites x 4 occasions. A covariate-varying interaction was biased below 600 sites x 10 single 0/1 occasions and unbiased at 1000 sites.
- The paper's priors are "vague" and given only in its Appendix S1. The module's Logistic(0,1) and N(0, SD 2) priors are a reasonable weakly informative choice.

The simulator and fit share only the formula and standardisation conventions. The planted-truth check below uses a separately written simulator, so a shared bug would have shown up there.

**Independent planted truth** (`verify_sim.R`, written without module code; nsim 100 per run).
- Design:
  - 700 sites, 5 occasions x 5 camera-days.
  - 8% of occasions entirely missing (NA) and 12% with partial effort.
  - Raw covariates in percent, and the two covariates are correlated.
  - The modifier (`conifer`) is not in the psi formula, so the module has to add its main effect.
  - An occasion covariate (rain) on subordinate detection.
- Parameters:
  - Dominant abundance: lambda = 1.5 at the mean.
  - Dominant detection r = 0.12 per day.
  - Subordinate: psi0 = 0.5, p = 0.2.
  - A weaker interaction with the opposite sign of habitat dependence: gamma0 = -0.4, gamma1 = -0.35.

| run | engine | result | gamma0 rel. bias (cov.) | gamma1 rel. bias (cov.) |
|---|---|---|---|---|
| 700 sites | nimble | **FAIL 1 of 13**: gamma0 bias -11.6% (limit 10%, MCSE 0.009); all else passes, 100/100 R-hat < 1.05 (max 1.014) | -11.6% (0.94) | -3.8% (0.92) |
| 700 sites, same data | ml | PASS all 13 | -7.0% (0.95) | +0.3% (0.95) |
| 1400 sites | nimble | PASS all 13 (max R-hat 1.015) | -3.2% (0.94) | -5.4% (0.98) |
| 1400 sites, same data | ml | PASS all 13 | -0.4% (0.94) | -2.7% (0.98) |
| null (gamma0 = gamma1 = 0), 700 sites | nimble | PASS all 13; **FPR 0.01 (gamma0), 0.05 (gamma1)** | bias -0.0002 | bias -0.009 |

All other estimands passed at 700 sites, with |bias| <= 3.6% and coverage 0.92-0.98. They include dominant abundance, total abundance, r, p, the rain effect on p, woodland and conifer slopes, and the proportions occupied.

The 700-site failure is the same finite-sample inflation of the interaction's magnitude that the builder found at 300 sites. It is not a coding error:
- The exact ML fit on the same data sets also overstates gamma0 (-7%).
- The bias disappears at 1400 sites.
- Coverage stays nominal.
- The posterior median is no better than the mean (-10.6% vs -11.2% on the first 48 replicates; `median_diag.R`).

The pooled false-positive rate over both null terms is 6/200 = 0.03.

**Builder's numbers reproduced.** These reruns used the builder's own command: `core/run.R recover ... --nsim 100` with default.yaml and occmediated.yaml.
- Default design, nimble: all 12 estimands pass. The table matches section 1 to every printed digit (identical seeds).
- Negative control (occupancy-mediated): fails as reported. gamma0 coverage 0.64, gamma1 coverage 0.82, woodland coverage 0.70.
- The gate can therefore fail. The builder's own null run gave an FPR of 0.08 / 0.03, and this independent null gave 0.01 / 0.05. Both are consistent with 0.05.

**Usability.**
- `core/run.R fit` on `example/config.yaml` was rerun with the final code. It took about 8 min on 4 cores, with max R-hat 1.001, min ESS 3279 and estimates identical to section 7.
- Every interpretation sentence was checked by hand against the estimate table. Examples:
  - Odds ratio exp(-0.53) = 0.59 (0.39-0.85).
  - Per-10-unit rate ratio exp(0.531 x 10/17.4) = 1.36.
  - Low/high-broadleaf odds ratios 0.25 / 1.42.
  - Two-step interval width 86% of the joint interval.
  - Predicted occupancy 0.62/0.50/0.37/0.26 for 0-3 martens.
- The wording is associational throughout ("is associated with"). Detection is labelled as a detection process. Results are on the response scale: odds ratios, rate ratios and probabilities.
- Every figure required by docs/FIGURES.md is present: subordinate_vs_dominant, interaction_by_habitat, dominant_abundance_<cov>, dominant_abundance_map and coefficients. So are interaction_strength, subordinate_occupancy, joint_vs_twostep, gof and trace. The curves agree with `fit$estimates`, for example psi(N=0) = 0.62 and psi(N=1) = 0.50 at average habitat. The axes use original units.
- Messy data (`messy_test.R`) gives clear, actionable messages for each problem tested:
  - NA in a covariate.
  - A misspelled site id: the run stops with "did you mean 'NI010'?".
  - A wrong column name: the message lists the columns that are present.
  - An empty date.
  - A dd/mm/yyyy date format.
  - Species-name case.
  - Trailing spaces: these are trimmed, with a warning.
  - End date before start date.
  - A misspelled modifier.
  - An empty site id.

**Changes made by the verifier** (only `R/io.R` and `R/report.R` were changed; none affect the likelihood or the samplers, so the recovery results above stand):
1. `read_data`:
   - `fm$p` partial-matched `fm$psi`, so leaving out `p:` silently gave detection the occupancy covariates. Fixed by using `[[ ]]`.
   - Unknown keys under `formulas:` now stop with an error.
   - A missing `formulas:` section now gives a warning.
   - Site/species whitespace is trimmed.
   - Empty dates or site ids, duplicate covariate rows and site-id typos are caught.
2. Factor habitat modifier: the report said gamma0 applied "at average habitat". It now says "where habitat = <reference level>", and each contrast is worded with its level and reference.
3. The detection unit is now "camera-day" or "survey occasion" as appropriate. "Multiplies the odds" was reworded to "is associated with ... times the odds".
4. Figures:
   - The coefficient forest-plot caption now says the slopes are per 1 SD (the axis labels carry raw units) and uses "credible" or "confidence" according to the engine.
   - The x-axis labels of the abundance map overlapped at 6-digit eastings. Fixed with 3 breaks and `check.overlap`.
5. The SMALL DESIGN note and the report caution now cite the 700-site independent result. The README not_for and Simulation sections say the same.

**Verdict: VERIFIED_WITH_LIMITS.**
- The module passes the gate at its default design (600 sites, builder) and at an independent 1400-site design.
- The negative control fails, and the FPR is about 0.05.
- At an independent 700-site design with a weaker interaction, the Nimble posterior mean overstates the interaction at average habitat by 11.6%, which is just outside the 10% rule. Coverage stays nominal.
- For a habitat-varying interaction with fewer than about 1000 sites, read the credible interval rather than the point estimate, and run a planted-truth check at the user's own design.
