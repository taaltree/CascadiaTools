---
id: nmixture
title: N-mixture models for repeated counts (closed and open)
when_to_use: The same sites (points, plots, transects, stream reaches) were counted two or more times within a period when the number of animals there did not change; the question is how many animals there are and what habitat or survey conditions are associated with abundance, allowing for animals that were present but not counted. With several seasons/years of such counts, the dynamic (Dail-Madsen) version estimates recruitment, apparent survival and the abundance trajectory.
not_for: Detection/non-detection only (use occu-single); individually identified animals (use scr, cjs or closed capture-recapture, which give far more reliable absolute abundance); single counts per site (no repeat visits, so detection cannot be separated from abundance; use a Poisson/NB GLM (glm module) for a relative index); counts where animals move freely in and out between visits (closure fails); distance-sampling data (use distance).
data_needed: counts per site per visit (long CSV, one row per visit, or wide site x visit table), site covariates (one row per site), optional visit covariates (wind, date, observer, effort); for the dynamic model also a season/year column and optional site-by-season covariates
engine: unmarked::pcount / pcountOpen (maximum likelihood, default); nimble + nimbleEcology (Bayesian, closed model only); hand-written likelihood (cross-check)
---

# N-mixture models for repeated counts

## Question it answers

"How many animals are at my sites, how does abundance change with habitat, and
how much do my counts underestimate it?" From repeated counts at many sites the
binomial N-mixture model (Royle 2004) separates two things the raw counts mix
together: the number of animals present at each site (abundance, N) and the
chance that each animal is counted on a visit (detection, p).

Outputs (names in the estimates table):

| estimand | meaning |
|---|---|
| `abundance_at_average_site` | expected number of animals at a site with average covariate values (mean of numeric covariates, reference level of factors) |
| `total_abundance` | estimated number of animals summed over the surveyed sites (realised N, with parameter uncertainty propagated) |
| `mean_abundance_per_site` | `total_abundance` / number of surveyed sites |
| `detection_prob` | probability that one animal present is counted on an average visit |
| `lambda_beta_<cov>` | log-scale effect of 1 SD of a site covariate on expected abundance |
| `rate_ratio_<cov>` | exp(lambda_beta): multiplicative change in expected abundance per 1 SD |
| `p_beta_<cov>`, `odds_ratio_detection_<cov>` | effect of a covariate on detection (logit scale / odds ratio) |
| `nb_size` | negative binomial size (mixture NB): variance = mean + mean^2/size |
| `zero_inflation_prob` | share of sites that are structurally empty (mixture ZIP) |
| dynamic: `initial_abundance_per_site`, `recruitment_per_site`, `apparent_survival`, `total_abundance_t1..tT`, `abundance_ratio_last_first` | first-season expected abundance, recruits per site per season, probability an animal survives and stays between seasons, total N per season, and ratio of last to first season |

The report also gives per-site abundance (caterpillar plot), conditional-effect
curves on the response scale, a goodness-of-fit test, a K-sensitivity check,
and (static) an AIC comparison of the Poisson, NB and ZIP abundance models.

## Model

Closed ("static", Royle 2004; Kéry et al. 2005):

    N_i ~ Poisson(lambda_i)  |  NegBin(mean lambda_i, size)  |  ZIP(lambda_i, zero-inflation psi)
    log(lambda_i) = beta0 + beta' x_i                (site covariates)
    y_ij | N_i ~ Binomial(N_i, p_ij)                  (count at site i, visit j)
    logit(p_ij) = alpha0 + alpha' w_ij                (site or visit covariates)

The likelihood sums over N_i = max(y_i.), ..., K. K is chosen automatically (start
at max count + max(30, 2 x max count) and doubled until the fitted abundance
distribution puts < 1e-4 probability above K) and then checked by refitting at
1.5K and 2K (Dennis et al. 2015 show that a too-small K biases N downward).

Open ("dynamic", Dail & Madsen 2011; unmarked::pcountOpen, "constant" dynamics):

    N_i1 ~ Poisson/NB/ZIP(lambda_i)
    N_it = S_it + G_it,  S_it ~ Binomial(N_i,t-1, omega_it),  G_it ~ Poisson(gamma_it)
    y_ijt | N_it ~ Binomial(N_it, p_ijt)

omega is apparent survival (death and permanent emigration are not separated),
gamma is recruits per site per season (births plus immigration). `dynamics: notrend`
fixes gamma = lambda (1 - omega) (stable expected abundance).

Derived abundance: for maximum-likelihood fits, `n_draws` parameter vectors are
drawn from the asymptotic normal distribution of the estimates; for each, every
site's N is drawn from its exact conditional distribution given its counts
(closed: P(N_i = k | y_i); open: forward filtering, backward sampling). Totals
therefore include both parameter and latent-state uncertainty. The Bayesian engine
does the same at each posterior draw.

## Assumptions (read these; they matter more than for most models)

1. **Closure within a season.** No births, deaths, immigration or emigration
   between repeat visits, and no animals moving between sites. Temporary
   movement in and out of a site between visits makes the model estimate the
   number of animals that *use* the site (a superpopulation), and detection
   becomes "available and detected". Strong violations bias N upward
   (Link et al. 2018). Keep the visits within a short window.
2. **No double counting and no false positives.** Each animal is counted at most
   once per visit. Double counting inflates N badly, and goodness-of-fit tests
   usually do not detect it (Link et al. 2018).
3. **Independent detection of individuals with the modelled p.** All animals at
   a site on a visit share the same p, apart from the covariates you include.
   Unmodelled heterogeneity in detection, among individuals, sites or visits,
   biases N (Duarte et al. 2018; Barker et al. 2018). Our negative control shows
   that site-level heterogeneity in p (logit SD 1) makes total N about 12% too low,
   with 95% intervals covering the truth in only 26% of simulations. The module's
   bootstrap GOF test flagged this strong site-level heterogeneity in 95% of
   simulated datasets (VALIDATION.md). Weaker or individual-level heterogeneity
   is much harder to detect: in the independent verification, individual-level
   heterogeneity (logit SD 1.2) made total N 16% too low (interval coverage 0.37),
   yet the GOF test flagged only 12.5% of those datasets.
4. **Sites are independent** and abundance follows the chosen distribution
   (Poisson, NB or ZIP) given the covariates. Overdispersed abundance fitted as
   Poisson biases N and gives over-confident intervals (Knape et al. 2018). In our
   negative control, NB(size 1) abundance fitted as Poisson underestimated total N by
   20%, and the intervals never covered the truth. The GOF test flagged it every time
   (c-hat about 3). Check the
   goodness-of-fit result and the AIC comparison.
5. **Dynamic model:** closure within seasons; between seasons, survival and
   recruitment follow the stated model; equal spacing between seasons. If a site
   was not surveyed at all in the first season, unmarked's `pcountOpen` starts that
   site's initial-abundance distribution at its first surveyed season, whereas the
   module's derived per-site and total N assume it applies to season 1 (the formal
   Dail-Madsen model). The two differ slightly for such sites, so survey every site
   in the first season if possible.

**Absolute N is fragile.** Barker et al. (2018) show that, without marked animals,
counts alone carry little information to separate N from p. Plausible alternative
models with no identifiable N produce data that look the same as N-mixture data,
especially when counts are sparse. Link et al. (2018) show that small assumption
violations cause large biases that GOF tests miss. Kéry (2018) found
identifiability problems in about 25% of real bird datasets under the NB mixture,
and almost none under Poisson/ZIP. **Report relative abundance (rate ratios, and
differences among sites or years) with more confidence than absolute totals.**
If absolute N matters (harvest quotas, population size for listing), collect
auxiliary data that inform p: marked animals, double-observer or removal
protocols, or distance sampling.

## Data format

Long format (default): one row per visit actually done.

```
point,visit,date,julian_day,wind_beaufort,observer,count
PT-001,1,2025-05-13,133,1,SAT,1
PT-001,2,2025-05-26,146,1,SAT,0
PT-001,3,2025-06-11,162,2,MLR,0
```

- site column (`site_col`, default `site`): site identifier.
- visit column (`visit_col`, default `visit`): 1, 2, 3, ... within site (and season).
- count column (`count_col`, default `count`): whole number >= 0. A visit not
  done is simply absent, or has a blank count.
- visit covariates: any further columns listed in `obs_covs` (numeric or text;
  text becomes a factor such as observer).
- dynamic: a season column (`season_col`, default `season`, e.g. year).

Site covariates (`site_covs`): one row per site, same site column. Season
covariates (`season_covs`, dynamic): site, season, covariates; the value in
season t drives the change from t to t+1.

Wide format (`format: wide`): one row per site with the count columns in visit
order (`count_cols`; for dynamic, all visits of season 1, then season 2, ...
with `visits_per_season`). Visit covariates come as wide files listed in
`obs_cov_files: {name: file.csv}`.

## Config

```yaml
data:
  model: static            # static | dynamic
  species: Hermit warblers # used in sentences and axis labels
  counts: point_counts.csv
  format: long             # long | wide
  site_col: point
  visit_col: visit
  count_col: count
  season_col: year         # dynamic only
  obs_covs: [julian_day, wind_beaufort, observer]
  site_covs: points.csv
  season_covs: reach_years.csv   # dynamic only, optional
formulas:
  lambda: "~ canopy_cover_pct + elevation_m"   # abundance (initial abundance if dynamic)
  p: "~ wind_beaufort + julian_day + observer" # detection
  gamma: "~ 1"             # dynamic: recruitment
  omega: "~ 1"             # dynamic: apparent survival
settings:
  engine: unmarked         # unmarked | nimble (static only) | rlik
  mixture: P               # P | NB | ZIP
  dynamics: constant       # constant | notrend
  standardize: true        # numeric covariates scaled to mean 0, SD 1 (effects per SD)
  K: null                  # upper bound on N per site; null = automatic + checked
  k_sensitivity: true
  compare_mixtures: true   # AIC of P / NB / ZIP (static, likelihood engines)
  gof_nboot: 200           # parametric bootstrap (ML) or posterior predictive replicates (nimble); 0 = skip
  n_draws: 1000            # parameter draws for derived abundance (ML engines)
  niter: 6000              # nimble only: iterations, burn-in, chains, thinning
  nburnin: 1000
  nchains: 3
  thin: 5
  prior_sd_lambda_int: 3   # nimble priors (normal SDs on link scale)
  prior_sd_p_int: 1.5
  prior_sd_beta: 2
  prior_sd_log_size: 2
```

Priors (nimble): log expected-abundance intercept ~ N(0, 3^2); logit detection
intercept ~ N(0, 1.5^2), roughly flat on the probability scale; coefficients of
standardised covariates ~ N(0, 2^2); log NB size ~ N(0, 2^2); zero-inflation ~ U(0, 1).

## Outputs

`<id>_estimates.csv`, `<id>_report.html` (methods, diagnostics, estimates,
interpretation, figures), and `figures/` (see Figures below). Diagnostics notes
give the GOF result, K sensitivity and the mixture AIC table.

## Figures

`core/run.R fit` exports these to `<out>/figures/` (PDF + PNG, `gallery.html`,
`figures_index.csv` with captions). Redraw from a saved fit without refitting:
`scripts/Rscript core/run.R figures nmixture --fit <out>/nmixture_fit.rds --config <cfg>`.

| id | what it shows |
|---|---|
| `abundance_<cov>` | expected abundance per site vs a site covariate (original units), other covariates at mean / reference level; line, `ci` and `inner_ci` bands, rug of observed values; categorical covariates as point + two intervals |
| `detection_<cov>` | detection probability per individual per visit vs a site or visit covariate (0-1 axis) |
| `recruitment_<cov>`, `survival_<cov>` | dynamic model only, when gamma / omega have covariates |
| `abundance_effects`, `detection_effects` (also `recruitment_effects`, `survival_effects`) | all effects on one process composed into one A/B/C figure (only when the process has two or more covariates) |
| `coefficients` | forest plot of rate ratios (abundance, recruitment) and odds ratios (detection, survival) per 1 SD or vs the reference level, log scale, null line at 1, `ci` and `inner_ci` intervals |
| `abundance_by_site` | caterpillar of estimated N per site (95% interval, fixed at fit time) with the largest count at each site (crosses); last season for the dynamic model |
| `total_abundance` | distribution of total N over surveyed sites (posterior or sampling-distribution draws), estimate, `ci` and `inner_ci` intervals, and the sum of maximum counts as a minimum |
| `trajectory` | dynamic only: total N per season with bands, sum of maximum counts as points |
| `gof` | parametric bootstrap (unmarked) or posterior predictive (nimble) chi-square check |
| `K_sensitivity` | plug-in total N refitted at increasing upper bounds K (axis spans +/- 5%) |

Label keys for `figures: labels:` — raw covariate column names (`canopy_cover_pct`),
factor levels (`JDK`), process axes (`abundance`, `detection`, `recruitment`,
`survival`), `coefficients` (x axis), any estimate name such as
`rate_ratio_canopy_cover_pct` (row of the forest plot), `abundance_by_site`,
`site_rank`, `site_estimate`, `site_max_count`, `total_abundance`, `draws`,
`min_count`, `trajectory`, `season`, `gof`, `replicates`, `K`, `K_sensitivity`.
Covariates are standardised internally; axes are back-transformed with the
recorded centre/scale, or with `figures: units:` if you give one.
Renamed factor levels (`JDK: "Observer 1"`) are used on the categorical effect
axis and in coefficient rows (not prefixed again with the variable name if the new
label already contains it). `greyscale: true` turns every module-drawn accent
(bands, largest-count crosses, minimum-count line) into print greys; shape and
line type remain the second cue. Long y-axis labels wrap onto two lines when
they would not fit the panel height after `per_figure` size overrides. Rugs of
covariates with few distinct or integer values (wind force, day of year) are
jittered (integers by at most +/- 0.4 units) so they are not hidden by axis
ticks; captions say so.

## Simulation and power

`design`: `model`, `n_sites`, `n_visits`, `n_seasons`, `prop_missing` (share of
visits not done), `abundance_dist` (P/NB/ZIP), and violations for negative
controls: `p_site_sd` (site-level logit-normal detection heterogeneity),
`p_visit_sd` (site-by-visit heterogeneity), `double_count_prob`.
`params`: `abundance_at_average_site` (or `initial_abundance_per_site`),
`detection_prob`, `lambda_beta_forest`, `p_beta_wind`, `nb_size`,
`zero_inflation_prob`, `recruitment_per_site`, `apparent_survival`.
Sensible knobs for power: `n_sites` (the dominant one for covariate effects),
`n_visits` (dominant for separating N from p), `detection_prob`. Validated at
150 sites x 3 visits; with p = 0.2 and only 2 visits absolute abundance is
unreliable (see VALIDATION.md).

## Pitfalls

- **Huge N with huge SE, or N rising with K:** N and p are not separately
  identified (low p, few visits, NB mixture). The fit is flagged as not
  converged. Treat abundance as an index, or add visits or auxiliary data.
- **Lack of fit (GOF p < 0.05):** try `mixture: NB`, add detection covariates
  (observer, wind, time), and check closure. An NB fit that sends N toward K is a
  known failure mode (Kéry 2018).
- **Small dynamic designs:** with about 70 sites x 4 seasons, a covariate effect on apparent survival
  (omega) was overestimated by about 14% on average (intervals still covered the truth 95% of the time); at
  200 sites the bias was 3% (VALIDATION.md, Independent verification). `nb_size` is also biased upward
  (+7-12%) in moderate samples; it is a nuisance parameter.
- **Visit-level covariates cannot go in the abundance formula.**
- **Effects are per 1 SD** of the covariate when `standardize: true`. The report
  also translates them to a round number of original units.
- **Factor covariates:** `abundance_at_average_site` refers to the reference
  (first alphabetical) level.
- **Bayesian engine:** closed model only; dynamic models use unmarked.
- **Run time (dynamic):** one `pcountOpen` fit takes 0.5-2 min at 60-100 sites x 4-5 seasons, and the
  K-sensitivity refits (at 1.5K and 2K) and the bootstrap GOF (`gof_nboot` refits, run in parallel over
  the job's cores) multiply that. Run it as a SLURM job with 16 cores (`scripts/ejob -c 16 -t 3:00:00`).
  For a quick first look, set `gof_nboot: 0` and `k_sensitivity: false`.

## References

- Royle JA (2004) N-mixture models for estimating population size from spatially replicated counts. Biometrics 60:108-115. doi:10.1111/j.0006-341X.2004.00142.x
- Kéry M, Royle JA, Schmid H (2005) Modeling avian abundance from replicated counts using binomial mixture models. Ecological Applications 15:1450-1461.
- Dail D, Madsen L (2011) Models for estimating abundance from repeated counts of an open metapopulation. Biometrics 67:577-587.
- Hostetler JA, Chandler RB (2015) Improved state-space models for inference about spatial and temporal variation in abundance from count data. Ecology 96:1713-1723. doi:10.1890/14-1487.1
- Dennis EB, Morgan BJT, Ridout MS (2015) Computational aspects of N-mixture models. Biometrics 71:237-246. doi:10.1111/biom.12246
- Barker RJ, Schofield MR, Link WA, Sauer JR (2018) On the reliability of N-mixture models for count data. Biometrics 74:369-377. doi:10.1111/biom.12734
- Link WA, Schofield MR, Barker RJ, Sauer JR (2018) On the robustness of N-mixture models. Ecology 99:1547-1551. doi:10.1002/ecy.2362
- Kéry M (2018) Identifiability in N-mixture models: a large-scale screening test with bird data. Ecology 99:281-288.
- Knape J et al. (2018) Sensitivity of binomial N-mixture models to overdispersion: the importance of assessing model fit. Methods in Ecology and Evolution 9:2102-2114. doi:10.1111/2041-210X.13062
- Duarte A, Adams MJ, Peterson JT (2018) Fitting N-mixture models to count data with unmodeled heterogeneity: bias, diagnostics, and alternative approaches. Ecological Modelling 374:51-59.
- Fiske I, Chandler R (2011) unmarked: an R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10).
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26:403-413.
- Goldstein BR, de Valpine P (2022) Comparing N-mixture models and GLMMs for relative abundance estimation in a citizen science dataset. Scientific Reports 12:12276 (nimbleEcology N-mixture distributions).
- Kéry M, Royle JA (2016) Applied Hierarchical Modeling in Ecology, Vol. 1, ch. 6. Academic Press.
