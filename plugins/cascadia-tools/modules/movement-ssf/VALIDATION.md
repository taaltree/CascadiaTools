# movement-ssf: validation

All runs 2026-10-18/19, R 4.5.3 (amt, glmmTMB, survival, ctmm, adehabitatHR from `Rlib/`) on
njobvu via `scripts/ejob` with 16 cores. Reproduce with `bash dev/movement-ssf/run_validation.sh`
(configs in `modules/movement-ssf/validation/`, outputs in `results/movement-ssf/validation/<run>/`,
logs in `logs/movement-ssf_*.log`). Recovery uses `fast_settings` (no RSF, home ranges or UHC; the
iSSA fit is identical to the full fit). Replicate r uses seed 100000 + r. Median 201 s per
replicate at the default design.

Gate (core `recovery_gate`, unchanged): |bias| <= max(10% of |truth|, 3 MC SE); 95% coverage >=
0.95 - 3 MC SE (floor 0.80); >= 90% of replicates converged. "Converged" = glmmTMB convergence
code 0, positive-definite Hessian and finite SEs. This is a likelihood module; there is no R-hat.

Simulation (`simulate.R`): 16 x 16 km landscape (50-m cells) with smooth forest-cover (%) and
elevation (m) fields (forest partly negatively correlated with elevation). Each animal moves
step by step: 200 candidate steps from Gamma(shape 2, scale 120 m) x von Mises(0, kappa 1)
relative to the current heading; one is chosen with probability proportional to
exp(beta_i' z(end) - 0.25 d^2), d = distance (km) to the animal's true range centre.
beta_i ~ Normal(mu, sd) per covariate (mu = 0.8 forest, -0.5 elevation per landscape SD;
sd = 0.4, 0.3). 2-h fixes, 5% failed fixes at random, 0-2 min time jitter, so the amt
regularisation/burst code runs in every replicate. Default design: 12 animals x 400 steps.

## Verdict

| design | reps (conv.) | result |
|---|---|---|
| **default** 12 animals x 400 steps (REML, estimated range centre) | 100 (100) | **PASS all 10 interpreted estimands**; `home_attraction_per_km2` (nuisance, not interpreted) FAILS |
| false-positive control: elevation mean = 0 (SD 0.3) | 100 (100) | **PASS all 10**; rejection rate for elevation 0.06 |
| harder: 5 animals x 150 steps | 100 (99) | 9/10 PASS; `sel_forest_sd` bias -23% (MC SE 6%) |
| negative control: RSF (use vs MCP availability) on the same autocorrelated tracks | 100 (100) | **FAILS as required**: elevation mean bias +40%, coverage 0.52; elevation SD coverage 0.65 |
| negative control: omit range residency (fit without the attraction term) | 100 (100) | **FAILS as required**: elevation mean coverage 0.83, turning kappa coverage 0.80, mean step coverage 0.85 |
| 1 available step per used step | 100 (100) | passes (unbiased) but RMSE +10% to +85% and elevation power 0.81 vs 0.93 -- not a failing control, see below |
| power over n_animals (elevation mean -0.25), 50 reps each | 4 x 50 | power 0.29 / 0.48 / 0.46 / 0.64 at 4 / 8 / 12 / 20 animals |
| independent simulator (3 covariates, 16 animals x 300 steps, CSV read path; section 8) | 100 (99) | **PASS all 13 interpreted estimands**; null FP rate 0.05; RSF control fails 6/6 |
| home ranges: AKDE on simulated OUF tracks (45 days, 2-h fixes) | 100 | AKDE bias -9%, coverage 0.85 (below 0.885); KDE -24%, MCP -40% |

The iSSA part of the module is **DONE at the default design** for the selection means and SDs,
RSS and the movement-kernel parameters. Not validated to the gate: the attraction coefficient
(nuisance), the among-animal SD with ~5 animals, and AKDE interval coverage with short tracks
(Known limits). No threshold was changed.

## 1. Default design (`results/movement-ssf/validation/default/`)

| estimand | truth | mean est | rel. bias | bias MC SE | coverage | power | gate |
|---|---|---|---|---|---|---|---|
| sel_forest_mean | 0.80 | 0.788 | -1.5% | 0.012 | 0.92 | 1.00 | PASS |
| sel_elevation_mean | -0.50 | -0.516 | -3.1% | 0.014 | 0.91 | 0.93 | PASS |
| sel_forest_sd | 0.40 | 0.391 | -2.3% | 0.011 | 0.98 | | PASS |
| sel_elevation_sd | 0.30 | 0.291 | -2.9% | 0.014 | 0.96 | | PASS |
| rss_forest_per_sd | 2.226 | 2.216 | -0.4% | 0.027 | 0.92 | | PASS |
| rss_elevation_per_sd | 0.607 | 0.603 | -0.5% | 0.009 | 0.91 | | PASS |
| step_shape | 2.00 | 2.007 | +0.3% | 0.004 | 0.97 | | PASS |
| step_scale_m | 120 | 119.8 | -0.2% | 0.30 | 0.97 | | PASS |
| step_mean_m | 240 | 240.3 | +0.1% | 0.30 | 0.95 | | PASS |
| turn_kappa | 1.00 | 0.999 | -0.1% | 0.003 | 0.98 | | PASS |
| home_attraction_per_km2 (nuisance) | 0.25 | 0.312 | +24.8% | 0.003 | 0.50 | | FAIL |

**Why REML is the default.** The first full run used ML (glmmTMB default; archived
`dev/movement-ssf/recover_default_ML_estcentre.log`): the means and movement parameters passed,
but `sel_elevation_sd` was 0.200 (rel. bias -33%) and `sel_forest_sd` -9%. Diagnostics
(32 reps each, `dev/movement-ssf/diag_*.log`): ML without any attraction in simulation or fit
still gave -15% (elevation SD), so this is the usual downward bias of ML variance components with
12 groups plus many fixed effects (the per-animal centre terms), and REML (`reml: true`) removed
it (-11% on 32 reps, -2.9% on the 100-rep run above).

**Wald intervals for the SDs.** The interval is built on the log-SD scale. When the SD estimate
is near zero the upper limit explodes (in 9% of default replicates the elevation-SD upper limit
exceeded 10; mean interval widths in the summary CSV are therefore astronomically large; median
width 0.68). These degenerate intervals count as "covering" in the coverage figure. Read an SD
whose interval runs to very large values as "not estimable from these animals".

## 2. False-positive rate (`null_elevation/`)

Elevation population mean 0 (individual SD 0.3 retained, so some animals select and some avoid).
The 95% interval for `sel_elevation_mean` excluded 0 in 6/100 replicates (nominal 5%); mean
estimate -0.015 (MC SE 0.013); coverage 0.94. All other estimands passed.

## 3. Harder design: 5 animals x 150 steps (`hard_few/`)

99/100 converged. Means: forest -2.5% bias, coverage 0.90; elevation +4.7%, coverage 0.93
(power only 0.28 for -0.5). `sel_forest_sd` 0.306 vs 0.40 (-23%, MC SE 0.022) **FAIL**;
`sel_elevation_sd` passed (+5%) but its RMSE (0.39) exceeds the truth. Attraction +87%.
With ~5 animals, report individual coefficients and treat the SD as rough (fit() writes a note).

## 4. Negative controls

* **RSF on autocorrelated GPS data** (`negctrl_rsf/`, `method: rsf`): use-availability logistic
  model (10 available points per fix in each animal's MCP, weight 5000, random slopes) reported
  under the iSSA names. Elevation mean -0.298 vs -0.5 (bias +40%), coverage 0.52; elevation SD
  0.417 vs 0.3, coverage 0.65. Forest passed by chance (bias -5%, coverage 0.91): the RSF estimates
  home-range-scale selection, which here differs from step-scale selection mainly for the
  broad-scale covariate (elevation). The gate fails, as it must.
* **Omitting range residency** (`negctrl_no_attraction/`, `home_attraction: false` on data with
  attraction): elevation mean coverage 0.83, turning kappa -3% with coverage 0.80, mean step
  coverage 0.85 -> 4 estimands FAIL.
* **One available step per used step** (`negctrl_one_control/`): the conditional logit is still
  consistent with one control, so this does not fail the gate (all 10 interpreted estimands
  pass); it costs efficiency: RMSE of the elevation mean 0.185 vs 0.144, forest mean 0.137 vs
  0.125, mean step 5.6 vs 3.0 m, and power for the elevation mean falls from 0.93 to 0.81.
  Recorded honestly: "too few available steps" is a precision problem, not a bias problem, in
  this model.
* ML-era plug-in centre (archived `dev/movement-ssf/recover_plugin_ML.log`): using each animal's
  mean location as a fixed centre (`home_centre: mean`) biased the elevation mean by +18%
  (coverage 0.73), because the mean location is pulled toward preferred habitat. This is why
  the centre is estimated (`home_centre: estimated`, default).

## 5. Power (`power/`, elevation mean -0.25, forest 0.8, 50 reps per design, 400 steps)

| animals | power, elevation (-0.25) | coverage | power, forest (0.8) | coverage | elevation-SD rel. bias |
|---|---|---|---|---|---|
| 4 | 0.29 | 0.94 | 0.96 | 0.86 | -20% |
| 8 | 0.48 | 0.88 | 1.00 | 0.98 | -26% |
| 12 | 0.46 | 0.90 | 1.00 | 0.90 | -13% |
| 20 | 0.64 | 0.96 | 1.00 | 0.92 | +1% |

With individual SD 0.3 comparable to a weak mean effect (-0.25), power is driven by the number
of animals, not steps; ~30+ animals would be needed for 80% power. Coverage of the means stays
near nominal (MC SE ~0.03 at 50 reps).

## 6. Home ranges (`hr/`, `dev/movement-ssf/hr_recovery.R`)

100 isotropic OUF tracks (ctmm::simulate; sigma^2 = 8e5 m^2, position autocorrelation 3 days,
velocity 2 h), 45 days of 2-h fixes; true 95% area 15.06 km^2. AKDE (ctmm.select + akde):
mean 13.67 km^2 (rel. bias -9.2%), 95% CI coverage **0.85** (below the 0.885 bound),
median effective sample size 13.5, 0 failures. KDE (reference bandwidth, amt) -24%, MCP -40%.
AKDE is far better than KDE/MCP but its interval under-covers at ~13 effective range crossings.

## 7. Cross-checks (`dev/movement-ssf/crosscheck.R`, one dataset, seed 7)

* **Likelihood engine.** Fixed-effects glmmTMB Poisson formulation vs `survival::clogit` (exact)
  vs a hand-written conditional-logit likelihood maximised with `optim`: max |glmmTMB - clogit|
  = 5.5e-4 (coefficients ~0.03-0.75, SEs agree to 4 digits); max |optim - clogit| = 2.1e-6.
* **Random slopes.** Population mean / SD from glmmTMB (ML) vs two-step approach (per-animal
  clogit + DerSimonian-Laird random-effects meta-analysis): forest mean 0.881 vs 0.885, SD 0.439
  vs 0.397; elevation mean -0.221 vs -0.208, SD 0.302 vs 0.334 (agree within one SE ~0.13).
* **Movement updates.** Module vs `amt::update_gamma` / `amt::update_vonmises`: shape 2.02812 vs
  2.02812, scale 113.8449 vs 113.8449 m, kappa 0.91898 vs 0.91898 (identical).
* **KDE area.** amt `hr_kde` vs adehabitatHR `kernelUD(h = "href")`, 95% area of one animal:
  16.11 vs 16.32 km^2 (1.3%, different grid/extent).
* **Attraction term with the true centre** (`dev/movement-ssf/att_test.R`, 8 datasets, 12 animals
  x 1500 steps, no habitat selection): clogit with the true squared distance 0.250 (truth 0.25);
  with the estimated centre 0.266; with the mean-location centre 0.265. The implementation is
  right; the bias comes from estimating the centre.

## 8. Independent verification (2026-10-19, `dev/movement-ssf/verify/`)

A second simulator written without any module code (`verify/indep_sim.R`): terra landscapes
(40-m cells) with THREE covariates (shrub %, slope deg, NDVI), an inverse-CDF von Mises sampler,
Gamma(1.5, 200 m) steps, kappa 0.5, attraction 0.15 /km^2, 16 animals x 300 steps, 3-h schedule,
10% failed fixes kept as empty rows, 0-5 min jitter, collar-style CSV with different column names
and a dd/mm/yyyy time format read back through the module's `read_data()` (the real-data path).

| run | reps (ok) | result |
|---|---|---|
| default (mu shrub 0.6, slope -0.4, NDVI 0.3; SD 0.3/0.2/0.25) | 100 (99 converged; 1 failed the Hessian/convergence check) | **all 13 interpreted estimands PASS** (bias -4.3% to +1.5%, coverage 0.90-0.97; slope SD -13.7%, within 3 MC SE); attraction +42%, coverage 0.29 FAIL (Known limit 1) |
| null (shrub mean 0) | 100 (99) | all 13 PASS; **false-positive rate for shrub 5/99 = 0.05** |
| RSF negative control (`method: rsf`) | 100 (100) | **all 6 FAIL** (SDs +35% to +87%, slope SD coverage 0.18; means coverage 0.87-0.88) |
| rerun of the default validation config (`verify/rerun_default`) | 100 (100) | identical to section 1 to all printed digits (deterministic seeds) |

Messy-data intake (`verify/messy_test.R`, `verify/messy_applied.R`): eight corrupted copies of
the example CSV. The verifier found that ids with different spellings ("E3" vs "E03",
"e05 " vs "E05"), rows with no id, and coordinates written as "511,014.9" were silently
dropped or treated as separate animals. The fix (applied 2026-10-19; pre-fix copies in
`dev/movement-ssf/pre_verify_backup/`) makes `check_data()` flag all of them, and it names the
animal whose coordinates are off the raster (wrong CRS for one collar). Wrong column names and
unreadable timestamps stop with a message that says what to set. Duplicated fixes were
already flagged. The same patch added wording to the report: a fix-success / habitat-biased
fix-loss caveat, the AKDE coverage caveat from section 6, and "higher" instead of "more" for
non-percentage covariates. It also themes the movement-kernel figure. None of these
changes touch the fitting code (`engine.R` and `simulate.R` are unchanged), so the recovery
results above still apply.

**Continuation (2026-10-19, after a session restart).** These checks were added to the work above.
* The logs (`logs/movement-ssf_verify_indep_{default,null,rsf}.log`) were re-read, and the numbers
  in the table above match them. The patched files in `verify/patched/` are byte-identical to the
  installed module. `engine.R` (01:32) predates every recovery run.
* Literature check of `engine.R` line by line:
  - The stratum intercept SD is fixed at 1e3 (variance 1e6), with independent `(0 + z | id)`
    slopes (Muff et al. 2020).
  - Gamma update: shape + beta_log(sl); rate 1/scale - beta_sl, where beta_sl is converted
    from per-km to per-m.
  - von Mises update: kappa + beta_cos(ta) (Avgar et al. 2016; Fieberg et al. 2021, Appendix).
  - Delta-method gradients for scale and mean.
  - RSS: exp(beta Delta z), per unit via unit/SD (Avgar et al. 2017).
  - The attraction expansion -a|x-c|^2 into home_d2 plus id:hx and id:hy terms.
  - RSF weight 5000.
  - UHC draws.
  All are correct. `simulate()` and `fit()` share only `ssf_lookup()`. The independent terra
  simulator (`cellFromXY`, GeoTIFFs read through `read_data()`) rules out cancelling lookup bugs.
* Example report checked against the estimates table:
  - Per-unit RSS: exp(0.777 x 10/24.7) = 1.37 and exp(-0.501 x 100/220) = 0.80.
  - Share of animals selecting: Phi(0.777/0.533) = 93% and Phi(-0.501/0.420) = 12%.
  - Fix loss: 410/5410 = 8%.
  - Median AKDE: (16.71 + 17.01)/2 = 16.9 km^2.
  - The RSS curve is consistent with the coefficient: RSS of about 3.2 at 93% forest and about
    0.25 at 11%.
  - Wording is associational, with the fix-loss caveat that separates detection from behaviour.
* Fix: the RSS curves are relative to the landscape-average covariate value, but core
  `report.R` does not embed figure captions, so the main report never said what the reference
  was. That was only in `figures/gallery.html`. `output.R` now adds one sentence per covariate
  that gives the reference value (e.g. 55.9% forest, 1250 m) and says the curve is for a typical
  animal on a log scale. The example was refit (`logs/movement-ssf_verify2_example.log`) and the
  estimates are unchanged. No fitting code was touched.
* Not fixed (core, not blocking):
  - The report's "Reproduce" line says `Rscript core/run.R ...`. On this cluster it should be
    `scripts/Rscript`.
  - Figure captions are not shown in the main HTML report for any module.

## Known limits

1. **Attraction coefficient is not recovered** (`home_attraction_per_km2`): +25% at 400
   steps/animal (coverage 0.50), +6% at 1500 steps, +87% at 150 steps, but unbiased when the
   true centre is supplied, and not fixed by 2000 candidate steps in simulation (`diag_nc2000`)
   or REML. This is the small-sample bias of a mean-reversion rate when the mean (centre) is
   estimated from the same autocorrelated track. It is a nuisance term: the report gives it only
   qualitatively and says not to interpret its size. Including it is still necessary (omitting it
   fails the gate, section 4).
2. **Among-animal SD with few animals**: biased low (-20% to -26% with 4-8 animals in the power
   runs; -23% for forest with 5 animals x 150 steps). Its Wald interval (log scale) is
   degenerate when the estimate is near zero.
3. **AKDE interval coverage 0.85** at ~13 effective range crossings (bias -9%). Treat AKDE CIs as
   somewhat too narrow for short deployments.
4. **"Too few available steps" is not a failing control here**: one available step per used step
   is unbiased (less precise). The negative controls that do fail are the RSF and the omitted
   range-residency term.
5. **Power for weak effects with individual variation** is low (0.46 at 12 animals for -0.25 per
   SD with SD 0.3); `run.R power` at the researcher's design is advised.
6. Simulated landscapes only use continuous, smooth covariates; categorical rasters,
   habitat-dependent movement (speed or turning varying with habitat), interactions between
   movement and habitat, and habitat-biased fix loss are not validated.
7. The RSF and UHC parts are comparison/diagnostic outputs; the RSF is shown to be biased for
   step-scale selection (section 4) and is not gated as an estimator of anything else.
