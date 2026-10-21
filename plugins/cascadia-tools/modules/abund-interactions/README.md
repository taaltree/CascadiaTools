---
id: abund-interactions
title: Abundance-mediated species interactions (dominant abundance -> subordinate occupancy)
when_to_use: Camera-trap (or other repeated detection/non-detection) surveys of two species where one (dominant, e.g. a predator or competitor) may affect where the other (subordinate) occurs, and the question is whether the subordinate's occupancy is associated with HOW MANY dominant individuals use a site, possibly differently in different habitats.
not_for: Interactions that act through presence/absence only and symmetric co-occurrence (use occu-rota or occu-waddle); fewer than ~300 sites (a habitat-dependent interaction passed the gate at 600 sites in the default design but its size was overstated by ~12% at 700 sites in an independent design; Twining et al. 2025 recommend ~1000 sites; see VALIDATION.md); very rarely detected dominant species (per-individual daily detection well below 0.03, untested); counts of individually identified animals (use scr); temporal avoidance at the scale of hours (activity-overlap analyses).
data_needed: camera deployments (site, start, end) + detection records (site, date, species), or a summarised site x occasion table of detection days and camera-days; site covariates (habitat), optional coordinates
engine: nimble (joint Bayesian, default); ml (exact marginal likelihood, cross-check); twostep (unmarked occuRN + occu plug-in, comparison only)
---

# Abundance-mediated species interactions

## Question it answers
Is the probability that a subordinate species (red squirrel) occupies a site associated with the
**number** of individuals of a dominant species (pine marten) using that site, and does the
strength or even the sign of that association depend on habitat (native broadleaf vs conifer
plantation)? It also estimates the dominant species' site abundance and its habitat drivers.

This is the model of Twining, Augustine, Royle & Fuller (2025, *Ecology* 106:e4468), who showed
that when an interaction acts through the dominant's abundance, occupancy-based interaction
models (Waddle et al. 2010; Rota et al. 2016) give biased interaction estimates. The habitat
modifier follows Twining, Sutherland, Reid & Tosh (2022, *Proc R Soc B*), who found pine martens
were associated with more red squirrels in broadleaf woodland but fewer in conifer plantations
(there with a Rota occupancy model).

## Model
Sites i = 1..M, occasions j = 1..J (e.g. weeks), each with K_ij camera-days.

Dominant species A (Royle & Nichols 2003):
- N_i ~ Poisson(lambda_i), log(lambda_i) = X_i beta_lambda
- yA_ij ~ Binomial(K_ij, 1 - (1 - r_ij)^N_i), logit(r_ij) = W_ij beta_r
  (yA_ij = number of days in occasion j with a detection; r = daily detection probability of
  one individual; K_ij = 1 everywhere gives the classic Royle-Nichols model.)

Subordinate species S:
- z_i ~ Bernoulli(psi_i)
- logit(psi_i) = X_i beta_psi + N_i (gamma0 + gamma1 H_i)   (H = habitat modifier; gamma1 optional)
- yS_ij ~ Binomial(K_ij, z_i p_ij), logit(p_ij) = W_ij beta_p

N_i is latent: the joint model uses all of its posterior uncertainty when estimating gamma.
gamma0 is the change in log-odds of subordinate occupancy per additional dominant individual at
average habitat (exp(gamma0) = odds ratio per individual); gamma1 is how much that changes per
1 SD (or level) of the habitat modifier.

Priors: logit intercepts of r, psi, p ~ Logistic(0,1) (uniform on the probability scale);
log-abundance intercept, all slopes and gamma ~ Normal(0, SD 2) (settings$prior_sd). Continuous
covariates are standardised.

Engines:
- `nimble` (default): the joint model above, 3 chains, AF_slice block samplers for
  (beta_lambda, beta_r) and (beta_psi, gamma); slice samplers for N_i.
- `ml`: same model, N_i (0..Nmax) and z_i summed out exactly; optim + Hessian. Used as an
  independent cross-check. `mediator: occupancy` makes psi depend on 1[N_i > 0] instead (the
  occupancy-mediated model; negative control only).
- `twostep`: unmarked::occuRN for A, then unmarked::occu for S with E[N_i | y] as a fixed
  covariate. Common shortcut; shown in every report for comparison (`compare_twostep`). Do not
  use it for inference (see VALIDATION.md: biased towards zero, intervals too narrow).
- `interaction: false` (nimble): reduced model with gamma fixed at 0; the subordinate part is
  then exactly single-species occupancy (cross-checked against occu-single).

## Assumptions (and what breaks)
- Closure of both species over the survey period (keep surveys short, e.g. 3-6 weeks).
- Royle-Nichols: each individual is detected independently with the same r; N_i Poisson among
  sites given covariates. Group-living or overdispersed abundance inflates N and distorts gamma.
  The posterior predictive check (Freeman-Tukey) flags gross violations.
- Sites far enough apart that the dominant's N_i are not shared (an individual can still use
  several sites; N_i is "individuals exposed to the camera", not density).
- The interaction acts in one direction (A on S occupancy). Detection of S is assumed not to
  depend on N_A (not modelled here).
- No unmodelled shared habitat drivers: if both species respond to an omitted covariate, gamma
  absorbs it. Include the main habitat drivers of both species.
- Association, not causation.

## Data format
`format: records` (default, camera-trap software exports):

deployments.csv
```
site,camera_id,start_date,end_date
NI001,C001,2024-03-02,2024-03-29
NI002,C002,2024-03-05,2024-03-17     # second row for NI002 if the camera was re-deployed
```
records.csv (one row per photo / independent event; only the day matters)
```
site,date,species,count
NI001,2024-03-03 11:12:07,Red squirrel,1
NI001,2024-03-04 02:41:10,Pine marten,1
```
sites.csv
```
site,easting,northing,woodland_pct,broadleaf_pct,dist_settlement_km
NI001,330600,341600,40.7,22.3,0.29
```
Camera-days are cut into occasions of `occasion_days` days from each site's first day; K_ij =
days the camera ran in the occasion, y_ij = days with >= 1 record of the species. Records on days
with no deployment are ignored (with a warning).

`format: occasions` (already summarised): one row per site x occasion with columns `site`,
`occasion` (1, 2, ...), `effort_days`, and one column per species holding the number of days with
detections (use 0/1 and effort_days = 1 for plain detection/non-detection occasions); optional
occasion covariates listed under `obs_covs`.

## Config
```yaml
data:
  format: records | occasions
  deployments, records          # records format
  occasions                     # occasions format (+ occasion_col, effort_col, dominant_col, subordinate_col, obs_covs)
  site_col, start_col, end_col, date_col, species_col, date_format, occasion_days (default 7), n_occasions (optional cap)
  dominant, subordinate         # species names exactly as in the data
  dominant_label, subordinate_label   # optional display names
  site_covs, site_cov_cols, coords: [x, y]
  interaction_modifier          # site covariate (numeric or factor) that modifies the interaction; omit for a constant one
formulas: {lambda, r, psi, p}   # one-sided R formulas; N_A (and N_A x modifier) are added to psi automatically
settings:
  engine: nimble | ml | twostep
  niter 20000, nburnin 5000, thin 5, nchains 3, prior_sd 2, interaction TRUE
  Nmax 60 (ml/twostep truncation), gof TRUE, compare_twostep TRUE
figures: (see docs/FIGURES.md; labels for covariates, `dominant`, `subordinate`)
```

## Outputs
Estimates (95% CrI): `dominant_mean_abundance` (individuals per site at average covariates),
`lambda_beta_<cov>` (log scale), `dominant_detection_per_individual` (per day),
`subordinate_occupancy_no_dominant` (psi where N_A = 0, average covariates), `psi_beta_<cov>`,
**`psi_beta_dominant_abundance`** (gamma0, logit change per dominant individual at average
habitat), **`psi_beta_dominant_abundance_x_<habitat>`** (gamma1),
`interaction_at_low_/high_<habitat>` (gamma at the 10th/90th percentile of habitat, or
`interaction_in_<level>` for a factor), `subordinate_detection`, `dominant_total_abundance`
(sum of N_i over surveyed sites), `dominant_prop_occupied`, `subordinate_prop_occupied`.

Figures: see the Figures section below.

## Figures
Exported by `core/run.R fit` (and redrawn without refitting by `core/run.R figures`) to
`<out>/figures/` as PDF + PNG with a `gallery.html` and `figures_index.csv` of captions.
All curves are on the response scale with covariates in original units; bands are the
`ci` (default 95%) and `inner_ci` (50%) credible intervals. Rename anything with
`figures: labels:` (covariate names, `dominant`, `subordinate`, `x`, `y`, coefficient names,
trace parameter names such as `gamma0`).

| id | shows |
|---|---|
| `subordinate_vs_dominant` | Subordinate occupancy probability vs dominant abundance at a site (0 to ~99th percentile of site N), full posterior band; rug = posterior mean N per site |
| `interaction_by_habitat` | Same curve in one panel per habitat level (10th/50th/90th percentile of a numeric modifier, or each factor level); only with `interaction_modifier` |
| `interaction_strength_<habitat>` | Interaction coefficient (logit change per dominant individual) vs the numeric modifier, null line at 0; only with a numeric modifier |
| `dominant_abundance_<cov>` | Expected dominant individuals per site vs each numeric abundance covariate |
| `subordinate_occupancy_<cov>` | Subordinate occupancy vs each numeric occupancy covariate at N_A = 0 and at the 90th-percentile site N_A (colour + linetype) |
| `dominant_abundance_map` | Posterior mean N_A at each surveyed site (viridis); metre coordinates shown in km; only if `coords` given |
| `coefficients` | Forest plot of all regression coefficients (link scale, per SD of standardised covariates), grouped by process (colour + shape) |
| `joint_vs_twostep` | Interaction estimates from the joint model vs the two-step plug-in; when `compare_twostep` |
| `gof` | Posterior predictive check (Freeman-Tukey) for each species, with Bayesian p-values |
| `trace` | MCMC traces of the intercepts and interaction coefficients (nimble engine) |
| `posterior` | Posterior densities of the same parameters (nimble engine) |

Entity colours follow the core palette order, so each habitat level / abundance scenario keeps
its colour wherever it appears. The `ml`/`twostep` engines give the same figures except
`trace`/`posterior`, with confidence instead of credible bands.

## Simulation and power
`params`: dominant_abundance (lambda at mean), lambda_beta_woodland,
dominant_detection_per_individual (daily r), subordinate_occupancy_no_dominant,
psi_beta_woodland, psi_beta_broadleaf, psi_beta_dominant_abundance (gamma0),
psi_beta_dominant_abundance_x_broadleaf (gamma1), subordinate_detection (daily p).
`design`: n_sites (default 600), n_occasions (4), days_per_occasion (7), prop_effort_loss (0.1:
share of occasions with camera failure), habitat_interaction (TRUE).
Sensible knobs for power: `n_sites`, `n_occasions`, `days_per_occasion`, and the effect sizes.
Twining et al. (2025) found a constant interaction recoverable with 300 sites x 4 occasions but a
habitat-varying one needed > 600 sites with 10 occasions when occasions are single
detection/non-detection. With daily camera records (K = 7 days per occasion) this module passes the
recovery gate at 600 sites x 4 occasions; at 300 sites the habitat x interaction term is overstated by
~11-19%, and the constant term is detected with power >= 0.92 even at 150 sites (VALIDATION.md).
An independent verification at 700 sites x 5 occasions (higher dominant abundance, weaker interaction)
found the interaction at average habitat overstated by 12% with Nimble (7% with the ML engine), with
nominal coverage. Treat 600-1000 sites as borderline for a habitat-varying interaction: run a
planted-truth check at your own design, and read the credible interval rather than the point estimate.

## Validation summary (details and all numbers in VALIDATION.md)
- Default design (600 sites, 4 x 7-day occasions), nimble, nsim 100: all 12 estimands pass the gate;
  interaction gamma0 bias -1.3% (coverage 0.96), gamma1 +3.8% (0.95); 100/100 R-hat < 1.05.
- Exact-likelihood engine on the same data: passes; agrees with unmarked::occuRN on the dominant part.
- Low dominant detection (r = 0.03/day): passes (99/100 converged); wider intervals, not bias.
- No interaction (gamma = 0): false-positive rate 0.08 (gamma0) and 0.03 (gamma1).
- Negative control (occupancy-mediated model on abundance-mediated data): fails as it should
  (gamma0 -49%, coverage 0.64).
- Two-step plug-in: fails (gamma0 attenuated 12%, gamma1 17%; intervals ~14% too narrow).
- Reduced model (interaction fixed at 0) reproduces the occu-single module to <= 0.015.
- Independent verification (own simulator, 700 sites x 5 occasions, weaker interaction): all pass except gamma0 overstated by 11.6% (Nimble; ML -7%, passes); passes at 1400 sites; false-positive rates 0.01 / 0.05. Verdict VERIFIED_WITH_LIMITS.

## Pitfalls
- psi_beta_dominant_abundance is per individual: with many individuals per site the implied
  change is large; read the `subordinate_vs_dominant` figure, not only the coefficient.
- `subordinate_occupancy_no_dominant` is an extrapolation if almost every site has the dominant.
- Low daily detection of the dominant makes N_i poorly known: at r = 0.03/day (600 sites) recovery still
  passes, but the interaction CrIs are ~25-30% wider and total dominant abundance ~2.4x wider
  (VALIDATION.md, hard design). Lower r was not tested.
- With fewer than ~1000 sites and a habitat modifier the size of the interaction tends to be
  overstated (by up to ~12% at 600-700 sites, 10-20% at 300 sites); reports add a SMALL DESIGN note.
- Do not report the two-step estimate as the result; it is shown only to illustrate how much
  uncertainty the shortcut hides.
- Royle-Nichols abundance is sensitive to the Poisson and independence assumptions: interpret
  relative differences in abundance more confidently than absolute numbers.

## References
- Twining JP, Augustine BC, Royle JA, Fuller AK (2025) Abundance-mediated species interactions. Ecology 106(1):e4468. doi:10.1002/ecy.4468
- Twining JP, Sutherland C, Reid N, Tosh DG (2022) Habitat mediates coevolved but not novel species interactions. Proc R Soc B 289:20212338. doi:10.1098/rspb.2021.2338
- Twining JP, Lawton C, White A, Sheehy E, Hobson K, Montgomery WI, Lambin X (2022) Restoring vertebrate predator populations can provide landscape-scale biological control of established invasive vertebrates: insights from pine marten recovery in Europe. Global Change Biology 28:5368-5384. doi:10.1111/gcb.16236
- Royle JA, Nichols JD (2003) Estimating abundance from repeated presence-absence data or point counts. Ecology 84:777-790.
- Waddle JH et al. (2010) A new parameterization for estimating co-occurrence of interacting species. Ecological Applications 20:1467-1475.
- Rota CT et al. (2016) A multispecies occupancy model for two or more interacting species. Methods Ecol Evol 7:1164-1173.
- Fiske I, Chandler R (2011) unmarked. J Stat Softw 43(10). de Valpine P et al. (2017) NIMBLE. J Comput Graph Stat 26:403-413.
