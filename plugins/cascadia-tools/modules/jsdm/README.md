---
id: jsdm
title: Joint species distribution model (HMSC)
when_to_use: Presence/absence or counts of many species (about 5-60) recorded at the same sampling units (plots, point-count stations, traps, lakes), when the question is how species respond to environmental covariates, whether species traits (and phylogeny/taxonomy) explain which species respond how, which species co-occur more or less often than the environment predicts (residual associations), how much of each species' variation is due to each covariate (variance partitioning), or how expected species richness changes along a gradient.
not_for: Imperfect detection that matters for the question (use occu-spatial msPGOcc / sfMsPGOcc, which separates occupancy from detection); one species (glm, occu-single); pairwise interaction hypotheses between two named species (occu-rota, occu-waddle, abund-interactions); ordination / dissimilarity questions (community); cover, biomass or proportion data; inferring biotic interactions from co-occurrence (not possible with any JSDM; see Pitfalls).
data_needed: community records (long: site, species, value; or wide: site x species), a sites table with covariates (optionally x/y coordinates), optional species traits table (and optional taxonomy columns or a Newick phylogeny)
engine: Hmsc 3.x (Bayesian Gibbs MCMC; probit (validated) or lognormal Poisson (not validated: slow mixing); latent factors with multiplicative gamma process shrinkage); glmmTMB reduced-rank model used as an independent cross-check
---

# Joint species distribution model (HMSC)

## Question it answers

*How do the species of a community respond to environmental gradients, do
species traits explain those responses, which species co-occur more (or less)
often than their shared environment predicts, and how does expected species
richness change along the gradients?*

The module fits the Hierarchical Model of Species Communities (HMSC;
Ovaskainen et al. 2017; Ovaskainen & Abrego 2020) with the R package Hmsc
(Tikhonov et al. 2020). All species are fitted in one model: each species
gets its own responses to the covariates, those responses are drawn from a
community distribution whose mean depends on species traits, and the
residual co-variation among species is captured by a few latent "site
factors".

## Model

Site i = 1..n, species j = 1..S, covariates x_i (continuous ones
standardised to mean 0, SD 1), traits t_j (standardised).

    L_ij = x_i' beta_j + eta_i' lambda_j                       (linear predictor)
    presence/absence:  y_ij = 1(L_ij + e_ij > 0),      e_ij ~ N(0, 1)          (probit)
    counts:            y_ij ~ Poisson(exp(L_ij + e_ij)), e_ij ~ N(0, sigma_j^2)  (lognormal Poisson)
    beta_j ~ MVN(Gamma t_j, V)                                  (traits -> responses)
      with phylogeny/taxonomy: vec(Beta) ~ MVN(vec(Gamma T'), V (x) [rho C + (1 - rho) I])
    eta_i ~ N(0, I_H) (non-spatial) or a Gaussian process over site coordinates (spatial: NNGP / full)
    lambda_hj ~ N(0, 1 / (phi_hj tau_h)), tau_h = prod_{l<=h} delta_l   (multiplicative gamma process)

- **Beta** (`beta_<covariate>_<species>`): Hmsc's species-specific responses
  on the probit (presence) or log (count) scale, per 1 SD of the covariate.
- **Gamma**: `comm_mean_<covariate>` = response of a species with average
  traits; `trait_<trait>_on_<covariate>` = change in that response per 1 SD of
  the trait (a trait-environment interaction). `trait_<trait>_on_intercept`
  means the trait shifts baseline prevalence/abundance.
- **Population-averaged responses (presence/absence; these are the ones
  interpreted).** With a probit link, Beta and Gamma are *conditional on the
  latent site factors*: a species' liability has variance 1 + Omega_jj, and
  binary data (one observation per site and species) cannot tell
  species-specific factor variance apart from the probit error. So the scale
  of Beta is only weakly identified, and in simulation Hmsc's Beta and Gamma
  are inflated by 10-30% (VALIDATION.md). The identified quantities are the
  population-averaged slopes, latent factors integrated out:
  - `popavg_<covariate>_<species>` = beta_kj / sqrt(1 + Omega_jj) (probit
    scale; P(occurrence) = Phi(popavg intercept + sum_k popavg_kj x_k) at a
    typical site);
  - `comm_popavg_<covariate>` = average of these over species;
  - `popavg_trait_<trait>_on_<covariate>` = least-squares slope of the
    population-averaged responses on the standardised trait (per posterior
    draw); `..._on_intercept` = effect of the trait on probit(prevalence).
    This describes the modelled species only: with a true Gamma of 0 it
    excludes zero in about half of simulated datasets. The report therefore
    pairs every supported trait effect with Gamma's interval, which is the
    test for the wider species pool.
  For counts (lognormal Poisson, log link) the log-scale Beta is already the
  population-averaged slope (E[y] = exp(x'beta + (Omega_jj + sigma_j^2)/2)),
  so `beta_*`, `comm_mean_*` and `trait_*` are interpreted directly.
  Gamma remains the right parameter for the question "is there a
  community-level/trait effect at all" (it is a hyperparameter; its sign and
  its false-positive rate are fine).
- **Residual associations.** Omega = Lambda'Lambda is the residual covariance
  among species. Two correlation summaries are reported:
  - `latcor_<a>__<b>` (**used for inference**): the latent-scale residual
    correlation cov2cor(Omega + D), with D = I (probit) or diag(sigma^2)
    (lognormal Poisson). It is the correlation of the two species' liabilities
    at a site after the covariates, i.e. how strongly their *occurrences* co-vary
    beyond the environment. It is well identified (see VALIDATION.md).
  - `rescor_<a>__<b>`: Hmsc's conventional `computeAssociations` correlation,
    cov2cor(Omega). It ignores the species-specific residual and is close to
    +/-1 whenever only a couple of factors are needed, and it is poorly
    identified for species with small loadings. It is reported for
    comparability with published HMSC analyses but **not interpreted**.
  - Support = posterior probability that the correlation is > 0; pairs with
    support >= 0.95 or <= 0.05 are flagged.
- **Derived**: `prev_<species>` (probability of occurrence at mean covariate
  values (binary 0/1 covariates at 0, their reference level), latent factors integrated out: Phi(beta_0j / sqrt(1 + Omega_jj))),
  `richness_mean_env` (expected species richness at mean covariates); for counts
  `abund_<species>` (expected count) and `total_abund_mean_env`;
  `phylo_signal_rho` when a phylogeny/taxonomy is supplied.
- **Variance partitioning** (`computeVariancePartitioning`): share of each
  species' linear-predictor variance due to each covariate and to the latent
  factors; R2 of Beta on traits.
- **Model fit**: Tjur R2 and AUC (presence) or pseudo-R2 (counts), in sample.
- **Posterior predictive check**: variance of site richness and mean pairwise
  co-occurrence in data replicated at new sites.

### Priors

Hmsc defaults for Beta, Gamma, V (inverse Wishart), sigma and rho, **except**
the latent-factor shrinkage prior: ecotools uses `a1 = 2, b1 = 1, a2 = 3,
b2 = 1` (nu = 3) instead of the Hmsc default `a1 = a2 = 50`. With the Hmsc
default the first factor's loadings have prior SD of about 0.14, which in our
simulations shrank residual correlations by about half and, because the probit
scale is tied to the latent variance, attenuated Beta by about 25% (slopes of
estimate on truth 0.51 and 0.73; VALIDATION.md section 0). a1 = 2, a2 = 3 (rate 1)
still shrinks later factors more than earlier ones: E(1/tau_h) = 1/{(a1 - 1)(a2 - 1)^(h-1)}
decreases in h because a2 > 2, and a2 > a1 widens the range over which the factor precisions are
stochastically ordered (Durante 2017, Lemma 2 and Table 2); it was chosen empirically in VALIDATION.md section 0. Set `settings: factor_prior:
{a1: 50, b1: 1, a2: 50, b2: 1}` to reproduce Hmsc defaults.

## Assumptions (and what breaks)

- **Detection is perfect or constant within species.** HMSC models where a
  species was *recorded*. Covariates that affect detectability (canopy,
  observer, noise) show up as "responses". Use occu-spatial (msPGOcc) when
  detection matters and repeat visits exist.
- **The covariates in the model are the relevant ones.** Omitting an important
  driver turns its effect into spurious residual "associations" between all
  species that respond to it (negative control in VALIDATION.md).
- **Residual associations are not interactions** (Blanchet et al. 2020; Zurell
  et al. 2018). They can come from missing covariates, shared dispersal,
  sampling, or interactions; the model cannot tell which.
- Sampling units are independent given covariates and latent factors. If sites
  are clustered in space, use `settings: spatial: nngp` with coordinates.
- Linear (on the link scale) responses unless you add quadratic terms in the
  formula (`formula: "~ temp + I(temp^2) + canopy"`).

## Data format

`community` (long; one row per record; site x species combinations not
listed are absences/zeros at every site in `sites`; repeated rows such as
visits are summed):

    station,visit,date,species,count
    CAS-001,1,2025-05-22,Setophaga occidentalis,1
    CAS-001,1,2025-05-22,Troglodytes pacificus,1

or `format: wide` (one row per site, one column per species).

`sites` (one row per sampling unit, every surveyed site including those where
nothing was found):

    station,utm_e,utm_n,elevation_m,july_temp_c,canopy_cover_pct
    CAS-001,576563,4880362,118,22.5,71

`traits` (optional; one row per species):

    species,common,family,order,body_mass_g
    Setophaga occidentalis,Hermit Warbler,Parulidae,Passeriformes,9.3

`phylogeny` (optional): Newick tree with species as tips (spaces may be
underscores), or `taxonomy: [order, family]` to build a taxonomic tree from
trait columns.

## Config

```yaml
data:
  community: bird_detections.csv  # file with species records
  format: long                    # long | wide (default: long if species_col exists)
  site_col: station               # site id column (community and sites files)
  species_col: species            # long format: species column
  value_col: count                # long format: count or 0/1 column
  species_cols: [...]             # wide format: species columns (default: all but site_col)
  response: presence              # presence (probit) | count (lognormal Poisson); default: presence if all 0/1
  sites: stations.csv
  covariates: [july_temp_c, canopy_cover_pct]   # or formula: "~ july_temp_c + I(july_temp_c^2) + habitat"
  traits: bird_traits.csv         # optional
  trait_species_col: species      # default species_col
  trait_formula: "~ log(body_mass_g)"   # or trait_cols: [...]
  taxonomy: [order, family]       # optional taxonomic tree from trait columns
  phylogeny: tree.nwk             # optional Newick phylogeny (instead of taxonomy)
  label_col: common               # optional: species labels used in text/figures
  x_col: utm_e                    # optional coordinates (needed for spatial: nngp/full)
  y_col: utm_n
  coord_units: m                  # m (converted to km) or km
  min_occurrences: 3              # drop species recorded at fewer sites (default 1)
settings:                         # defaults shown in module.R
  samples: 1000                   # posterior samples per chain
  thin: 50
  transient: 25000                # burn-in iterations (example: ~15 min with 4 parallel chains)
  n_chains: 4
  n_parallel: 4                   # chains run in parallel (cores)
  max_factors: 5                  # maximum number of latent factors
  spatial: none                   # none | nngp | full (needs x_col/y_col)
  factor_prior: {a1: 2, b1: 1, a2: 3, b2: 1}
  omit_covariates: []             # sensitivity: refit without these covariates
figures:
  labels: {july_temp_c: "July mean max. temperature (°C)"}
```

## Outputs

Estimates table (CSV) with the parameters above; HTML report with
interpretation, methods paragraph and figures; publication figures
(`species_responses`, `residual_associations`, `variance_partitioning`,
`trait_effects`, `richness_<cov>`, `occurrence_<cov>`, `model_fit`, `ppc`,
`trace`, `posterior`).

## Simulation and power

`simulate()` draws data from the HMSC model. The "community" (traits, Beta,
loadings) is fixed by `design$community_seed`; replicates differ in sites,
covariates, latent site scores and noise.

- params: `comm_mean_intercept`, `comm_mean_temperature`, `comm_mean_forest`,
  `trait_log_body_mass_on_{intercept,temperature,forest}`, `species_sd`
  (SD of species' deviations from the trait model), `loading_sd1`,
  `loading_sd2` (strength of residual associations), `residual_sd` (counts).
- design: `n_sites`, `n_species`, `n_factors` (0 = no residual associations),
  `response` (presence | count), `phylo_rho`, `spatial_range_km`,
  `cov_correlation`.
- `design$community`: `random` (default; each replicate draws new species
  deviations and loadings from the trait model, so Gamma is judged as a
  hyperparameter) or `fixed` (same species in every replicate).
- Power: `power_params` are the community-level effects: Gamma
  (`comm_mean_*`, `trait_*`) and, for presence data, their population-averaged
  counterparts (`comm_popavg_*`, `popavg_trait_*`). Use Gamma for "is there an
  effect" (false-positive rate at Gamma = 0 is in VALIDATION.md).
  `scripts/Rscript core/run.R power jsdm --vary n_sites=50,100,200 --nsim 50`.

## Validation status (VALIDATION.md)

- **Presence/absence (probit): validated** at 100 sites x 12 species and 100 sites x 20 mostly rare
  species. Population-averaged responses, latent residual correlations, prevalences and richness pass
  the recovery gate. Gamma false-positive rate 0.04-0.05; spurious supported associations 0.1% of
  pairs when none exist. Hmsc's conditional Beta/Gamma are inflated by 10-25% (coverage nominal).
  Cross-checked against a glmmTMB reduced-rank probit JSDM (popavg r = 0.986, latcor r = 0.918).
- **Counts (lognormal Poisson): NOT validated.** Species-level parameters did not converge (R-hat
  1.3-2.0) even at 2.5x the recovery budget. Community-level Gamma is recovered. Prefer presence/absence,
  or run much longer chains and check the report's convergence line.
- Spatial random levels and phylogenetic rho run but have no recovery validation.

## Pitfalls

- Rare species (< ~5 occurrences) carry almost no information on their own;
  their Beta is pulled to the trait/community mean (that is intended
  shrinkage, but don't interpret their individual responses). Drop them with
  `min_occurrences` or report them as uncertain.
- Residual correlations need many sites: below ~100 sampling units most
  pairs will lack support, even when associations exist.
- Mixing: Beta and the latent factors trade off on the probit scale; check
  R-hat < 1.01 and ESS >= 400 before reporting and increase `thin` /
  `transient` if needed (the report says so).
- In-sample Tjur R2 includes the latent factors and is optimistic.

## References

- Ovaskainen O, Tikhonov G, Norberg A, et al. (2017) How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters 20:561-576. doi:10.1111/ele.12757
- Ovaskainen O, Abrego N (2020) Joint Species Distribution Modelling: With Applications in R. Cambridge University Press. doi:10.1017/9781108591720
- Tikhonov G, Opedal ØH, Abrego N, Lehikoinen A, de Jonge MMJ, Oksanen J, Ovaskainen O (2020) Joint species distribution modelling with the R-package Hmsc. Methods in Ecology and Evolution 11:442-447. doi:10.1111/2041-210X.13345
- Bhattacharya A, Dunson DB (2011) Sparse Bayesian infinite factor models. Biometrika 98:291-306.
- Durante D (2017) A note on the multiplicative gamma process. Statistics & Probability Letters 122:198-204.
- Blanchet FG, Cazelles K, Gravel D (2020) Co-occurrence is not evidence of ecological interactions. Ecology Letters 23:1050-1063.
- Zurell D, Pollock LJ, Thuiller W (2018) Do joint species distribution models reliably detect interspecific interactions from co-occurrence data in homogenous environments? Ecography 41:1812-1819.
- Warton DI, Blanchet FG, O'Hara RB, et al. (2015) So many variables: joint modeling in community ecology. Trends in Ecology & Evolution 30:766-779.
- McGillycuddy M, Popovic G, Bolker BM, Warton DI (2025) Parsimoniously fitting large multivariate random effects in glmmTMB. Journal of Statistical Software 112(1). doi:10.18637/jss.v112.i01
- Tjur T (2009) Coefficients of determination in logistic regression models. The American Statistician 63:366-372.
- Vehtari A, Gelman A, Simpson D, Carpenter B, Bürkner P-C (2021) Rank-normalization, folding, and localization: an improved R-hat. Bayesian Analysis 16:667-718.
