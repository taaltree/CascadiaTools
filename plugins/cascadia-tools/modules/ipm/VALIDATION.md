# ipm: validation

Gate (`core/R/recovery.R`, unchanged): |bias| <= max(10% of truth, 3 MC SE);
95% coverage >= 0.95 - 3 MC SE (about 0.884 at nsim = 100); >= 90% of
replicates converged (R-hat < 1.05 under the recovery settings
`fast_settings`: 3 chains x 35,000 iterations, 5,000 burn-in, thin 6).
Full tables (every estimand: truth, mean, bias, MC SE, RMSE, coverage, CI
width, power) were produced by `dev/ipm/summarise_validation.R` and are in
`dev/ipm/verify/validation_tables_verified.md`. Reps and summaries are in
`results/ipm/recover_<design>/` and `results/ipm/verify/<run>/`.

**Verdict: VERIFIED WITH LIMITS.** Population size, growth rate, mean
survival, fecundity, recapture and the covariate effect pass at both the
default design and an independent design. The temporal SDs, the LTRE
contribution sizes and immigration (at low immigration) do not pass and are
documented limits below.

## 1. Default design (15 years, ~100 females, 50 juveniles + 20 adults marked per year, p = 0.6, Poisson counts)

`core/run.R recover ipm --nsim 100 --cores 16 --out results/ipm/recover_default`
(rerun 2026-10-20 with the current estimand definitions). 97/100 converged,
0 errors, median 558 s per replicate.

| estimand | truth | rel. bias | coverage | gate |
|---|---|---|---|---|
| juvenile_survival | 0.25 | +1.6% | 0.97 | PASS |
| adult_survival | 0.48 | +0.1% | 0.95 | PASS |
| fecundity_yearling | 2.8 | +0.0% | 0.92 | PASS |
| fecundity_adult | 3.6 | +0.6% | 0.95 | PASS |
| recapture_prob | 0.6 | -0.8% | 0.96 | PASS |
| juvenile_survival_beta_winter_severity | -0.4 | -5.1% | 0.97 | PASS (power 0.80) |
| juvenile_survival_sd | 0.3 | +7.6% | 0.99 | PASS |
| fecundity_yearling_sd / fecundity_adult_sd | 0.15 | +6.2% / +3.0% | 0.98 / 0.98 | PASS |
| adult_survival_sd | 0.2 | +18.6% | 0.99 | FAIL (bias) |
| immigration_rate (study-period mean of the annual rates) | 0.127 (realised) | -5.3% | 0.93 | PASS |
| immigration_rate_median (exp(mu); not interpreted) | 0.12 | -12.3% | 0.93 | FAIL (bias) |
| immigration_sd | 0.3 | +145% | 1.00 | FAIL (bias) |
| immigrants_mean | 14.2 | -6.9% | 0.91 | PASS |
| N_first / N_final / N_mean | | +3.0% / +0.2% / +0.5% | 0.97 / 0.96 / 0.94 | PASS |
| lambda_geomean | 1.011 | -0.2% | 0.97 | PASS |
| ltre_juvenile_survival | 0.0221 | +25.4% | 0.93 | FAIL (bias) |
| ltre_immigration | 0.0012 | +937% | 0.93 | FAIL (bias) |
| var_lambda_explained | 0.0282 | +51.8% | 0.91 | FAIL (bias) |
| ltre_adult_survival, ltre_fecundity_*, ltre_population_structure | small | +50% to +425% | 0.93-0.98 | PASS only because 3 MC SE > 10% of a tiny truth |

LTRE ranking: the largest true contributor (juvenile survival in all 100
replicates) was ranked first in 96/100 fits; the 4 misses ranked immigration
first.

## 2. Harder designs

**Short series (6 years), nsim = 100, 100/100 converged.** Mean rates,
recapture, beta (-6.2%, but power only 0.07), N and lambda pass. Every temporal
SD fails (+38% to +293%), as do immigration_rate_median (+56%) and all large
LTRE terms. immigration_rate passes (+12.6%, within 3 MC SE). The LTRE top
contributor is right in 77/100 (97/100 if immigration is excluded).

**Sparse capture-recapture (20 juveniles + 8 adults marked per year), nsim = 100, 97/100 converged.**
Mean rates, N, lambda and immigration_rate (-10.0%, coverage 1.00) pass;
power for beta drops to 0.70. adult_survival_sd fails (+51%).

## 3. IPM vs separate analyses (same simulated data sets, same seeds)

Separate = the same CJS m-array + productivity models without counts
(`components: separate`). Ratio of mean 95% interval widths, IPM / separate:

| estimand | default design | sparse capture-recapture |
|---|---|---|
| juvenile_survival | 0.92 | 0.79 |
| adult_survival | 0.94 | 0.87 |
| juvenile_survival_beta | 0.87 | 0.71 |
| juvenile_survival_sd | 0.86 | 0.70 |
| fecundities | 1.00-1.01 | 1.00-1.01 |
| recapture_prob | 0.98 | 0.95 |

RMSE is lower in the IPM in each case (juvenile survival 0.027 vs 0.029, and
0.030 vs 0.038 with sparse marking). Coverage is similar. As expected, the
gain is largest for the rates with the weakest direct data and is zero for
fecundity, which the productivity data already pin down.

## 4. Controls

- **Negative control (data with immigration of 0.12, fitted without immigration),** nsim = 100:
  juvenile_survival +18.7% (coverage 0.53), adult_survival +4.8% (coverage 0.88),
  recapture_prob -7.0% (coverage 0.71). Gate FAILS as intended. Fecundity is not biased
  (+0.8%, +1.4%): the productivity data hold it in place.
- **False-positive rate (beta = 0)**, nsim = 100: 3/98 = 0.031.
- **Power** over `marked_juveniles` 15 / 30 / 60 (nsim = 50): 0.67 / 0.88 / 0.82 for beta = -0.4
  (`results/ipm/power`). The dip at 60 is within MC error (SE about 0.055).

## 5. Independent routes

- **Full IPM in JAGS (book-style code, `dev/ipm/crosscheck.R`) vs Nimble on the example data:**
  every posterior mean agrees within 0.05 posterior SD; SD ratios 0.97-1.07. In
  `results/ipm/crosscheck/ipm_nimble_vs_jags.csv` the row `immigration_rate` is the
  median, exp(mu) (the file predates the change of estimand).
- **CJS component vs `marked` (maximum likelihood, individual histories):** differences of
  0.03-0.11 SE for juvenile survival, adult survival and recapture. The custom
  `dMarr` log-likelihood equals `dmultinom` to 2e-13.

## 6. Example data (`example/config.yaml`, truth in `example/TRUTH.txt`)

Default settings converge (max R-hat 1.004, min bulk ESS 630 over all reported
quantities). The posterior predictive p-values are 0.56 (counts), 0.46
(productivity) and 0.37 (capture-recapture). The truth lies inside the 95% CrI
for every mean rate, beta, immigration_rate (0.20, CrI 0.12-0.30, truth
0.127), immigrants_mean and N_final. N_first (108, CrI 91-127, truth 86) is
outside: the first count was 106 for 86 true females. In this data set the
LTRE puts immigration first (0.042), while the truth has juvenile survival
first (0.019 vs immigration 3e-5). This is the documented overestimation of the
immigration contribution, and the report says so.

## Known limits

- **Temporal SDs** (`*_sd`) are biased upward with 15 or fewer years (adult_survival_sd
  +19% in the IPM, +31% in the separate analysis; +59% at the 12-year independent design).
  They are descriptive only.
- **Immigration.** immigration_sd is essentially prior-driven (+145% to +293%). The
  study-period immigration_rate passes at about 12% immigration, but at low immigration
  (about 8%, independent design below) it is biased upward by +19% (immigrants_mean +14%),
  with good coverage (0.98-0.99). The skew comes from the overestimated SD. Treat
  low immigration estimates as upper values. immigration_rate_median is not interpreted.
- **LTRE contributions:** absolute sizes are biased upward (the posterior variance of the
  annual rates adds to their temporal variance), and the immigration contribution is
  overestimated severalfold (+937% to +2565%; Paquet et al. 2021). Only the ranking is
  supported: the top contributor was correct in 96% of replicates (default), 92%
  (independent design) and 77% (6 years).
- **Short series (< ~8 years):** report N and mean rates only.
- **Model assumptions:** one recapture probability for both age classes (the book's
  hoopoe model has a separate pj), immigrants breed like yearlings, no density
  dependence, independent data sets.

## Independent verification (2026-10-20, separate agent)

**Literature and code.** I checked the process model, priors, m-array cell
probabilities and LTRE line by line against the book code: Schaub & Kéry (2022),
`IPM_09/IPM_09.2.R` (hoopoe IPM) and `IPM_09.4.R` (fixed-design and transient
LTRE). The code is in `dev/ipm/ref/IPM_code-main`.

- `N1[t+1] ~ Pois(rho phij (f1 (N1+N3) + f2 N2))`, `N2 ~ Bin(Ntot, phia)` and
  `N3 ~ Pois(omega Ntot)` are identical to the book with rho = 0.5.
- The priors match (U(0,1) for means and p, U(0,3) for SDs and mean omega). The
  fecundity bound is U(0,20) here vs U(0,10) in the book.
- `dMarr` cell probabilities equal the book's `pr.j`/`pr.a` when pj = p.
- The lambda expression, sensitivities at hyper-means and mean stage proportions,
  and `rowSums(V * outer(s, s))` are the book's 9.4 code.
- The transient LTRE uses stage proportions, where the book code differences
  absolute N. Proportions are the consistent choice (Koons et al. 2016).
- The truth for the LTRE comes from `ipm_ltre_draws`, the same function as the
  estimate (shared code). I recomputed it with analytic sensitivities in my own
  simulator, and the two agree to machine precision on three data sets.
- Paquet et al. (2021, MEE 12: 1899-1910) citation confirmed.

**Independent planted truth** (`dev/ipm/verify/vsim.R`, `vrun.R`, written without
module code). My simulator writes raw researcher CSVs, so `read_data`, m-array
construction and covariate standardisation are exercised:

- a long capture file with age codes `pullus`/`ad`;
- per-nest productivity with mother age `first-time`/`older`;
- non-default column names;
- the covariate in original units (°C).

The design and parameters differ from the defaults:

- 12 years, about 80 females, 40 juveniles + 25 adults marked per year, 60% of broods monitored;
- phij 0.30, phia 0.55, f 2.2/3.0, p 0.5, median omega 0.08;
- SDs 0.25/0.15/0.12/0.12/0.25, beta 0.5 per SD.

Recruits are drawn per female, thinned by sex and survival. nsim = 100, 99/100
converged, 0 errors. Results are in `results/ipm/verify/main`.

| estimand | rel. bias | coverage | gate |
|---|---|---|---|
| juvenile_survival | -1.2% | 0.97 | PASS |
| adult_survival | -0.1% | 0.96 | PASS |
| fecundity_yearling / adult | -0.5% / -0.1% | 0.96 / 0.97 | PASS |
| recapture_prob | +1.0% | 0.95 | PASS |
| juvenile_survival_beta_temp_C | +7.6% | 1.00 | PASS (power 0.90) |
| N_first / N_final / N_mean | +2.5% / +0.4% / +0.5% | 0.92 / 0.97 / 0.95 | PASS |
| lambda_geomean | -0.1% | 0.95 | PASS |
| immigration_rate_median | +6.3% | 1.00 | PASS |
| immigration_rate (study mean) | +18.8% | 0.98 | FAIL (bias) |
| immigrants_mean | +14.1% | 0.99 | FAIL (bias) |
| immigration_sd | +258% | 0.98 | FAIL |
| juvenile_survival_sd / adult_survival_sd | +36% / +59% | 0.98 / 0.99 | FAIL |
| ltre_juvenile_survival / ltre_immigration / var_lambda_explained | +26% / +2565% / +61% | 0.94-0.97 | FAIL |

LTRE top contributor correct in 92% of replicates (99% excluding immigration).

- **Null effect (beta = 0), independent design:** false-positive rate 3/99 = 0.030
  (`results/ipm/verify/null`). In that run immigration_rate was +5.6% (PASS),
  so the immigration bias varies between replicate sets.
- **Negative control, independent design (immigration fitted out),** nsim = 48:
  juvenile_survival +11.7% FAIL (bias, just above the 10% bound; coverage 0.92), adult_survival +3.5%. The gate fails, but only barely, because immigration is low (8%) in this design; at 12% immigration the failure is clear (+18.7%, coverage 0.53).
- **Builder's recovery reruns** (default, short, sparse_cr, null_effect) checked: the numbers above
  were regenerated from those runs, and every claim in the builder's summary was confirmed.
  The run with the new immigration_rate estimand passes at the default design (-5.3%).
- **Usability.** I ran a messy copy of the example data (`dev/ipm/verify/messy`) with a wrong
  count column name, an NA count, an NA fledgling number, a misspelled ring (`AF7339l`), age
  codes `juv.` and `fledgeling`, and a misspelled covariate name. The wide capture format
  gives the same m-arrays as the long format.

**Changes made in this verification** (`io.R`, `report.R`; the model, simulate and fit are unchanged,
so no recovery result is affected):

1. A missing-column error now names the real config key (e.g. `count_col:`); it used to suggest a
   non-existent `counts_*_col`.
2. Age codes with trailing dots (`juv.`) are recognised. Unrecognised age codes are listed by value
   and id in the data notes (they used to be counted only).
3. Broods with an empty fledgling number are reported (they used to be dropped silently), with
   advice to enter 0 for failed broods.
4. Ring ids with an unusual letter/digit format are flagged with the closest common id
   (`'AF7339l' (close to 'AF73391')`). There are no false alarms on the example or simulated data.
5. A misspelled covariate name now lists the available covariate columns.
6. Interpretation: the caveat said omitting immigration biases "juvenile survival and
   fecundity". Simulation shows juvenile and adult survival, not fecundity, so I corrected it.
   I added a caution that low immigration estimates are biased upward. Covariate sentences
   use the config figure label (e.g. "Sahel rainfall (mm, Jul-Sep)"), and missing first/last
   counts print as "not taken" instead of "NA".
7. Figures: the immigration-covariate effect plot looked up a non-existent annual column
   (`immigration_rate_median[year]`), so its annual points were missing. It now uses `immigration_rate[year]`.

## Independent verification, round 2 (2026-10-20, separate agent)

This round continued the interrupted second pass. It reused its finished SLURM runs
(`results/ipm/verify2/`) after checking them, and applied the fixes that pass had
drafted but never applied.

**Second independent simulator** (`dev/ipm/verify/v2/sim2.R`, no module code):

- individual-based: per-female Poisson fledglings, Binomial sex, per-bird survival and recapture;
- 18 years, about 150 females, 45 juveniles + 25 adults marked per year, p = 0.7;
- phij 0.30, phia 0.50, f 2.2/2.8, median omega 0.15;
- **normal** count error (SD 8) with one missing count year;
- **wide** capture format and **aggregated** productivity;
- the covariate (snow days, raw units) acts on **adult** survival (beta 0.3 per SD);
- the LTRE truth uses hand-derived sensitivities.

nsim = 100 per run, fitted with `fast_settings` through `read_data()` → `check_data()` → `fit()`.

| run | converged | result |
|---|---|---|
| main (beta 0.3) | 97/100 | PASS: all mean rates (bias <= 2%, coverage 0.93-0.98), all temporal SDs except immigration_sd (-3% to +8%), beta (+2.3%, power 0.89), immigration_rate (-4.4%), immigration_rate_median (-7.1%), immigrants_mean (-5.6%), N_first/N_final/N_mean, lambda_geomean, ltre_adult_survival, ltre_juvenile_survival. FAIL: immigration_sd +81%, count_sd +12% (coverage 1.00), ltre_immigration +495%, var_lambda_explained +19%, ltre_fecundity_yearling (tiny truth) |
| null (beta 0) | 92/100 | false-positive rate for the covariate 3/92 = **0.03**; the other estimands are as in main |
| negative control (immigration present, fitted without it) | 99/100 | gate FAILS as intended: juvenile_survival **+23.8%, coverage 0.14**; adult_survival +6.4%, coverage 0.76; recapture_prob -6.4%, coverage 0.58; fecundity unbiased |

- LTRE top vital-rate contributor correct in 73% (main) and 57% (null) of replicates.
  It is lower here than at the default design because the true top two contributors
  were within a factor of 2 of each other in about 60% of replicates (correct 72% then,
  86% when the factor is >= 2). The ranking is reliable only when one rate clearly dominates.
- With 18 years and normal counts, the temporal SDs of survival and fecundity were
  recovered without bias. Their upward bias is a short-series effect (<= 15 years).
- **Builder's default recovery rerun** (`results/ipm/verify2/rerun_default`, same seeds): the
  estimates are identical to the builder's run (max |difference| = 0), so every number in
  section 1 is reproduced.

**Changes in round 2** (`io.R`, `report.R`, README; model, simulate and fit unchanged, so no
recovery result is affected; backups are in `dev/ipm/verify/module_backup_v2/`):

1. A count file with two rows for the same year now stops with both values and advice. Before,
   it was accepted silently.
2. A mistyped year (e.g. 2106 in the capture file) now stops, naming the file and the year.
   Before, it silently stretched the study to 101 years.
3. Broods with an unrecognised mother age (e.g. `3yr`) are now reported. Before, they were
   dropped from the productivity data without a note.
4. Report: covariate sentences give the size of 1 SD in the covariate's own units. The
   transient-LTRE components are ordered (and coloured) like the vital-rate series.
5. README: count_sd is documented as a biased nuisance parameter.

Messy-data test (`dev/ipm/verify/v2/messy_test2.R`, output `messy_test2.out`): each of the
three problems above, plus a misspelled age (`juvenille`), a misread ring (`NX5607l`) and an
empty fledgling number, gives a specific, actionable message. The clean example gives none.

Example refit with the patched module (`core/run.R fit ipm --config modules/ipm/example/config.yaml
--out results/ipm/verify3_example`) gave the same estimates as section 6 (max R-hat 1.002). I checked
every report sentence against `ipm_estimates.csv`, e.g. odds ratio exp(0.305) = 1.36 (1.02-1.80),
survival 0.17-0.29 over 307-522 mm, and 26.8 immigrants per year = 17% of breeders. Survival is
called "apparent" and the covariate an "association". All 11 figures listed in the README (with one effect plot for the one covariate) are
present; the effect and transient-LTRE PNGs were inspected visually.

**Round-2 verdict: VERIFIED WITH LIMITS** (unchanged limits: temporal SDs at <= 15 years,
immigration SD and low-immigration bias, LTRE sizes and the immigration contribution, count_sd).
