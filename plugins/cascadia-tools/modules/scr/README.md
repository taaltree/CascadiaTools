---
id: scr
title: Spatial capture-recapture (closed population density)
when_to_use: Individually identified animals detected at fixed detectors over repeated occasions within one closed season (camera traps of marked/spotted species, DNA from hair snares, scat stations or rub objects, live traps); question is population density (per km2), abundance, ranging scale, sex ratio or sex differences in ranging, or whether density is associated with habitat covariates. Handles 0/1 detections (binomial or per-occasion Bernoulli) and counts of genetic samples (Poisson), uneven detector effort, habitat masks from a raster, polygon or grid, detector / occasion / sex covariates on detection and a trap-specific behavioural response.
not_for: Unmarked animals (use nmixture / occupancy / distance sampling), partially marked populations (scr-smr), open populations across seasons or years with survival and recruitment (scr-open / jolly-seber), samples found along search tracks rather than at fixed detectors (scr-search), telemetry-only data (movement modules), genotypes with unresolved identity errors (fix them in the lab first; see NGS section).
data_needed: detector CSV (id, x, y [, detector covariates]); detection or genetic-sample CSV (individual, occasion, detector [, sex, QC status]); optional effort table (detector x occasion), occasion covariate table, individual sex table, habitat grid CSV / GeoTIFF raster / polygon (shapefile, GeoPackage)
engine: nimble (Bayesian marginal likelihood, all variants [default]; latent data augmentation for the basic binomial model), numpyro (Bayesian NUTS, CPU/GPU; binomial/Poisson with mask, density and detector covariates), secr (maximum likelihood; independent cross-check of every variant)
---

# scr: spatial capture-recapture

## Question it answers

*How many animals live here, at what density, and is density associated
with habitat? Do males range farther than females? Are animals trap-happy
at baited sites?* SCR uses **where** each identified individual was
detected across an array of detectors. Animals living near many detectors
are detected often and at several detectors; animals near the edge,
rarely. The model turns that spatial pattern into:

| estimand | meaning |
|---|---|
| `density_per_km2` | expected animals per km2 of habitat, averaged over the state space (`density_per_unit2` if map units are unknown) |
| `abundance` | realised number of activity centres in the habitat of the state space (detector array plus buffer) |
| `baseline_detection` | p0: probability that a detector placed at an animal's activity centre records it on one occasion (binomial / Bernoulli models; at the reference level of detection covariates) |
| `baseline_encounter_rate` | lambda0: expected number of samples a detector at the activity centre collects per occasion with full effort (Poisson model) |
| `sigma` | spatial scale of the half-normal detection function, in map units (m) |
| `home_range_area_95` | pi (2.448 sigma)^2, the 95% area of a circular bivariate-normal home range (km2) |
| `*_female`, `*_male` | the three quantities above per sex when sex-specific detection is fitted |
| `proportion_female` | sex ratio of the population (including undetected animals) |
| `density_beta_<cov>` | log rate ratio of density per 1 SD of a habitat covariate (log-linear density; several allowed) |
| `detection_beta_<cov>` | effect of a detector- or occasion-level covariate on p0 (logit) or lambda0 (log); 0/1 covariates: level 1 vs 0, continuous: per 1 SD |
| `sigma_beta_<cov>` | log ratio of sigma per unit of a detector-level covariate (e.g. a scent lure that draws animals from farther) |
| `detection_beta_male`, `sigma_beta_male` | male vs female on the logit/log scale of p0/lambda0 and on log sigma |
| `detection_beta_trap_response` | trap-specific behavioural response: change in logit p0 (log lambda0) at a detector after the animal's first detection there (> 0 trap-happy, < 0 trap-shy) |

## Model

Activity centres s_i are an inhomogeneous Poisson point process over the
state space S, represented by a habitat mask of pixel centres g (only
habitat pixels; non-habitat such as lakes or farmland is excluded):

    log D(g) = beta0 + sum_k beta_k x_k(g)              (x standardised over the habitat)

Detection of individual i at detector j on occasion k, with effort u_jk
(0 = not operating; days active / nominal days for hair snares):

    binomial   y_ij  ~ Binomial(K_j, p_ij)              (occasions aggregated; 0/1 usage)
    Bernoulli  y_ijk ~ Bernoulli(1 - (1 - p_ijk)^u_jk)  (per occasion; partial effort)
    Poisson    y_ijk ~ Poisson(u_jk lambda_ijk)         (counts of samples per occasion)
    p_ijk = p0_ijk exp(-d(s_i, j)^2 / (2 sigma_ij^2))    lambda_ijk = lambda0_ijk exp(...)   (half-normal)
    link(p0_ijk) = link(p0) + w_j' alpha + v_k' gamma + male_i beta_m + b_ijk beta_b     (logit for p0, log for lambda0)
    log sigma_ij = log sigma + w_j' delta + male_i delta_m
    male_i ~ Bernoulli(1 - proportion_female)            (unknown sex summed out)
    b_ijk = 1 if i was detected at j before occasion k   (trap-specific learned response, secr "bk")

Effort enters on the hazard scale, as in secr (Efford, Borchers & Mowat
2013). The binomial and Bernoulli models give the same likelihood when
nothing varies by occasion; the module aggregates over occasions whenever
it can (fast) and switches to the per-occasion form automatically when
there are occasion covariates, a behavioural response or partial effort.
The Poisson model is the standard model for noninvasive genetic sampling
where one animal can leave several samples at one detector in one session
(Royle et al. 2014 ch. 9; Efford secr "count" detectors with Poisson
counts; Gardner et al. 2009, 2010 for black bear hair snares). With the
hazard parameterisation lambda0 is the expected number of samples, not a
probability.

Sex: individuals of known sex contribute with their sex; those whose sex
could not be determined contribute a mixture over both sexes weighted by
the sex ratio, and undetected animals are a mixture too (sex as a
partially observed individual covariate under data augmentation, Royle et
al. 2014; the same model as the secr hybrid mixture `hcov`, a finite
mixture with known class for sexed animals, Pledger 2000). This assumes sex is equally likely
to be missing for males and females.

This is the closed SCR model of Efford (2004), Borchers & Efford (2008) and
Royle & Young (2008); Royle et al. (2014) is the textbook.

Engines (`settings$engine`):

- **nimble, approach `marginal` (default; every variant).** Data
  augmentation to M potential individuals (Royle & Dorazio 2012), with
  z_i, s_i and unknown sex summed out over the mask in one compiled custom
  distribution (`dSCRgen`, `R/nimble_engines.R`). Only psi, p0/lambda0,
  sigma, the proportion female and the coefficients are sampled (adaptive
  block random walk). Abundance is drawn per posterior draw as
  n + Binomial(M - n, P(z = 1 | never detected)). The compiled likelihood
  equals an independent plain-R implementation (`scr_loglik`) to 1e-13
  for every variant (VALIDATION.md).
- **nimble, approach `latent`.** Classic data augmentation with z_i and s_i
  sampled; binomial model with detector covariates only (cross-check).
- **numpyro.** The same marginal likelihood with NUTS in JAX
  (`python/scr_numpyro.py`): binomial / Bernoulli (aggregated) and Poisson,
  effort, irregular masks, density and detector covariates. CPU by default;
  `numpyro_gpu: true` on the fw_gpu partition for large arrays.
- **secr.** Maximum likelihood (`secr.fit`) on the identical mask:
  proximity detectors (binomial/Bernoulli, half-normal g0) or count
  detectors with Poisson counts (`binomN = 0`, hazard half-normal lambda0),
  usage = effort, occasion covariates via `timecov`, sex via `hcov`
  (`g0 ~ h2`, `sigma ~ h2`, `pmix`), behaviour via `bk`. This is the
  independent route; use `compare_engines: [secr]`.

Priors (Bayesian engines): psi ~ U(0,1) (equivalent to a flat prior on N
up to M), p0 ~ U(0,1) or lambda0 ~ U(0, 10) (`lambda0_max`), sigma ~ U(0,
half the state-space width), proportion female ~ U(0,1), all coefficients
~ Normal(0, SD 2) (covariates standardised).

## Assumptions (and what breaks)

- **Closure:** no births, deaths, immigration or emigration during the
  survey. Keep surveys short relative to the species' turnover. Violations
  inflate N and sigma.
- **Fixed activity centres, circular half-normal space use.** Strong sex
  differences in ranging (male carnivores and bears range 1.5-3x farther)
  cause heterogeneity: fit `sex_effects: [p0, sigma]` when sex is known for
  at least some animals. Sex must be missing at random with respect to sex.
- **Independent detections** among individuals (family groups violate
  this). A trap-happy/shy response is allowed only in its trap-specific form
  (`behaviour: trap`); ignoring a real trap-happy response biased p0
  upward (+43%) and sigma downward (-6%) in simulation, while density was
  little affected in that design (VALIDATION.md 12b).
- **Correct identities.** Every genotype or photo ID is assumed correct.
  See "Noninvasive genetic sampling" below.
- **Effort recorded correctly.** Detector-occasions when a detector was not
  working must be 0; partial effort must be entered. Ignoring effort biases
  detection and density (VALIDATION.md negative control).
- **State space large enough:** the buffer beyond the outermost detectors
  should be at least 3-4 sigma (of the wider-ranging sex), otherwise animals
  living outside S are detected and squeezed into it, and density is
  overestimated (VALIDATION.md: buffer 1 sigma -> density +47%, coverage
  0.34). The module warns when buffer < 3 sigma. **Check it** with
  `buffer_check: true`, which refits with the buffer enlarged by 2 sigma
  (mask built from detectors, polygon or raster) or, for a fixed grid CSV,
  trimmed by 2 sigma, and reports how much density changed within the region
  both fits share (< 5% and < 0.3 posterior SD = adequate; comparing the
  same region keeps a change in habitat composition from looking like a
  buffer problem when density depends on covariates). In validation it
  flagged 16/16 data sets with a 1.5-sigma buffer and 0/16 with 4 sigma.
- **Mask fine enough:** pixel spacing about sigma/2 or finer (warning if
  coarser than sigma). Aim for 1,000-4,000 pixels.
- **Detector spacing relative to sigma.** Spacing no more than about 2 sigma
  (Sun et al. 2014; ideally 1-2 sigma), so most animals are detected at
  more than one detector. The report states your spacing/sigma ratio.
- **Enough animals and spatial recaptures.** The posterior MEAN of density
  leans high (by roughly CV^2 or more; CV = posterior SD / mean): it fails
  the 10% bias tolerance once density is imprecise, posterior CV above ~0.2.
  That happened with ~18-33 detected animals in the camera designs and with
  ~46 detected animals in the independent sex-specific Poisson hair-snare
  design (two density covariates, 30% unknown sex: posterior-mean density
  +19%, while secr's MLE on the same data was +8%; VALIDATION.md, Independent
  verification round 2). More parameters need more animals. The report warns
  below 25 detected individuals or 10 spatial recaptures, and when the
  posterior CV of density exceeds 0.2; then quote the interval and the secr
  estimate (`compare_engines: [secr]`) alongside the posterior mean.
- **Density-covariate false positives.** With ~33-42 detected animals and no
  true effect, ~10% of 95% intervals for a density covariate excluded 0 (both
  Nimble and secr; 5% expected). Treat an interval that only just excludes 0
  as weak evidence.
- **0/1 data with a trap-specific response.** With Bernoulli/binomial data,
  p0 (before the first detection at a detector) and the response trade off:
  in simulation the Bayesian posterior mean of p0 was ~16% too high and the
  response biased toward trap-shy (centred at -0.38 when the truth was 0),
  while density, abundance and sigma were unbiased; secr's MLE was unbiased
  for both (VALIDATION.md 10f, 12c). Use `compare_engines: [secr]` and
  report secr's p0 and response. Poisson counts do not show this.
- **Density covariates** are associations across space, not causal effects.

## Noninvasive genetic sampling (NGS) workflow

Hair snares, scat stations, rub objects / rub trees and similar passive DNA
detectors are fixed detectors: use this module (samples found along search
routes belong in `scr-search`).

1. **Sample-level QC before analysis.** Genotyping errors (allelic dropout,
   false alleles) create false "new" individuals (ghosts) and inflate N; low
   marker variability merges relatives (shadow effect) and deflates N
   (Taberlet et al. 1999; Paetkau 2003; Lukacs & Burnham 2005). Before using
   this module the lab should: genotype with enough loci for a probability of
   identity among siblings well below 1/N; replicate amplifications (or
   re-analyse) for low-quality samples; drop samples with too few loci,
   mixed samples (more than two alleles) and unresolved mismatches; and
   re-check all pairs of genotypes that differ at only 1-2 loci (error-checking
   protocol of Paetkau 2003). Put the QC result in a column and set `qc_col`
   (and `qc_pass`); only passing samples are analysed and the report states
   how many were dropped.
2. **Counts or 0/1.** If one animal can leave several samples at one detector
   per session, use `obs_model: poisson` (counts of samples). If the lab
   genotypes only one sample per detector per session, or you collapse
   multiple samples, use `bernoulli`/`binomial`. Do not collapse counts and
   then call them Poisson.
3. **Effort.** Record the days each detector was active per session
   (`effort:` table; `effort_scale: 7` turns days into a fraction of a
   session). Destroyed, flooded or unchecked detectors get 0.
4. **Sex** from a sex marker (e.g. amelogenin, ZFX/ZFY) goes in `sex_col`
   (per sample; the most frequent call per individual is used) or an
   `individuals:` table. Empty / "U" = unknown (kept, sex treated as latent).
5. **Baited / lured detectors** often produce a trap-happy response
   (Gardner et al. 2010 found a positive trap-dependence for black bears at
   baited hair snares): consider `behaviour: trap` and compare with the
   model without it.
6. **Known limit: identity errors are NOT modelled.** Models that do
   (latent multinomial, Link et al. 2010; misidentification models of Lukacs
   & Burnham 2005 and Wright et al. 2009) are not implemented. If the lab
   reports a non-negligible per-sample error rate after QC, density will be
   biased (ghosts: upward); say so in the report and prefer stricter QC.

## Data format

Detector file (one row per detector; camera, snare, station):

    snare_id,utm_e,utm_n,lure
    HS01,511873,4655089,1

Detections or samples (one row per photo-detection or per genetic sample;
extra columns ignored). Binomial/Bernoulli: repeats of one animal at one
detector in one occasion count once. Poisson: every row is one sample (or
give `count_col`).

    sample_id,session,snare_id,bear_id,sex,loci_typed,qc_status
    BH24-3624,1,HS02,UA004,F,8,pass
    BH24-1187,1,HS07,,,3,fail_low_loci

Effort (optional), wide (detector + one column per occasion) or long
(detector, occasion, effort):

    snare_id,session_1,session_2,session_3,...
    HS01,2,7,7,...            # days active; with effort_scale: 7 -> 0.29, 1, 1

Instead of an effort table a usage string per detector ("11110000", 1 =
operating) can be a column of the detector file (`usage_col`).

Occasion covariates (optional), one row per occasion:

    occasion,rain_mm,mean_temp_c
    1,12.4,18.0

Individuals (optional; alternative to `sex_col`): `individual,sex` with
F/M/female/male, empty or U for unknown.

Habitat / state space, one of:
- `mask:` a regular grid CSV of pixel centres covering detectors + buffer,
  with density covariates and an optional 0/1 habitat flag (`habitat_col`);
- `habitat_polygon:` a polygon of suitable habitat (shapefile .shp or
  GeoPackage .gpkg, same projected CRS as the detectors; holes = non-habitat);
  the mask is the grid within `buffer` of a detector and inside the polygon;
- `raster:` a GeoTIFF whose layers are the density covariates (names must
  match `density_covariates`), optionally with a 0/1 habitat layer
  (`habitat_layer`; 0 or NA = non-habitat); can be combined with a polygon;
- nothing: the mask is built from the detectors with buffer 4 x the rough
  sigma from repeat detections and pixels of sigma/2 (homogeneous density).

Coordinates must be projected (metres, e.g. UTM), not latitude/longitude.

## Config

```yaml
data:
  traps: hair_snares.csv          # required: detectors
  trap_id_col: snare_id           # default id
  x_col: utm_e                    # default x
  y_col: utm_n                    # default y
  units: m                        # m or km -> density per km2; omit if unknown
  detection_covariates: [lure]    # detector columns affecting p0 / lambda0 (0/1 kept, continuous standardised)
  sigma_covariates: []            # detector columns affecting sigma
  usage_col: weeks_active         # optional 0/1 string per occasion, OR:
  effort: snare_effort.csv        # optional effort table (wide or long)
  effort_trap_col: snare_id       # detector id column of the effort table (default trap_id_col)
  effort_scale: 7                 # divide effort by this (days -> fraction of a session)
  n_occasions: 6                  # required if no usage_col / effort table
  captures: bobcat_detections.csv # detections, OR
  samples: hair_samples.csv       # genetic samples (same columns)
  individual_col: bear_id         # default individual
  occasion_col: session           # default occasion (1..K)
  trap_col: snare_id              # default trap
  count_col: null                 # optional number of samples per row (Poisson)
  qc_col: qc_status               # optional sample QC column; rows not in qc_pass are dropped
  qc_pass: [pass]
  sex_col: sex                    # optional sex per row (F/M; empty = unknown), OR
  individuals: bears.csv          # optional individual table (individual, sex)
  sex_effects: [p0, sigma]        # default when sex is given; [] = ignore sex
  occasions: sessions.csv         # optional occasion table (column occasion + covariates)
  occasion_covariates: [rain_mm]  # occasion-level detection covariates
  obs_model: poisson              # binomial (default) | bernoulli | poisson
  behaviour: trap                 # none (default) | trap (trap-specific learned response)
  mask: habitat_grid.csv          # habitat grid (x, y + covariates), OR
  habitat_col: habitat            #   optional 0/1 flag in the grid
  habitat_polygon: study_area.gpkg   # polygon of suitable habitat, and/or
  raster: habitat.tif             # GeoTIFF with density-covariate layers
  habitat_layer: null             # optional 0/1 habitat layer in the raster
  buffer: 12000                   # mask buffer when built here (map units)
  mask_spacing: 1000              # mask pixel size when built here
  density_covariates: [forest_pct, elevation_m]
  labels: {forest_pct: "Forest cover (%)"}
settings:
  engine: nimble                  # nimble | numpyro | secr
  approach: marginal              # marginal | latent (nimble; latent = basic binomial only)
  obs_model: null                 # override data obs_model
  behaviour: null                 # override (none | trap)
  sex_effects: null               # override, e.g. [sigma] or []
  density_covariates: null        # override (e.g. [] for homogeneous density)
  detection_covariates: null      # override detector-level p0 covariates
  occasion_detection_covariates: null
  sigma_covariates: null
  niter: 20000                    # nimble marginal iterations per chain (run longer for many parameters)
  nburnin: 5000
  nchains: 3
  sampler: RW_block               # or AF_slice
  niter_latent: 30000             # nimble latent
  nburnin_latent: 10000
  thin_latent: 5
  numpyro_warmup: 1000
  numpyro_samples: 2000
  numpyro_gpu: false              # true when run with scripts/ejob -g
  M: null                         # augmentation size (default max(5n, n+150); auto-doubled if N presses on it)
  buffer: null                    # trim the state space to this distance from detectors
  buffer_check: false             # refit with buffer +/- 2 sigma and report the change in density
  sigma_max: null                 # upper limit of the sigma prior (map units)
  lambda0_max: 10                 # upper limit of the lambda0 prior (Poisson)
  compare_engines: [secr]         # also fit with these engines; figure engine_comparison
  ppc_draws: 300                  # posterior predictive checks
  surface_draws: 200              # draws for the density maps
```

Examples: `example/config.yaml` (bobcat camera traps, binomial, grid-CSV
mask with a reservoir, canopy covariate, on-trail detection covariate) and
`example/ngs/config.yaml` (black bear hair snares: Poisson sample counts,
effort table, QC column, sex from amelogenin with missing calls, habitat
polygon + GeoTIFF covariates, scent-lure detection covariate, trap-happy
response, buffer check).

## Outputs

Estimates table (above), diagnostics (R-hat, bulk ESS; the real-data
standard is R-hat < 1.01 and ESS >= 400), design notes (spacing/sigma,
buffer/sigma, spatial recaptures, buffer check), predictive p-values, and
figures: `density_surface` (realised density over the habitat mask with
detectors and estimated activity centres), `density_expected_surface` and
`density_<cov>` (conditional effect of each habitat covariate on density,
response scale), `detection_function` (p or expected samples vs distance,
by sex or by a 0/1 detector covariate, detector spacing marked),
`detection_<cov>` (baseline detection vs each detector/occasion covariate,
response scale), `detection_trap_response`, `sex_differences` (p0/lambda0
and sigma by sex, sex ratio), `detector_effort` (map of effort and samples;
NGS / partial effort), `buffer_check`, `abundance_posterior`,
`captures_per_individual` (observed vs predicted), `gof`, `coefficients`,
`engine_comparison` (when `compare_engines` is set), `trace`, `posterior`.

## Simulation and power

`params`: `density` (per km2), `baseline_detection` (binomial/Bernoulli)
or `baseline_encounter_rate` (Poisson), `sigma` (m), `density_beta_canopy`,
`density_beta_elevation`, `detection_beta_on_trail`, `detection_beta_rain`,
`sigma_beta_lure`, `detection_beta_trap_response`, `proportion_female`,
`baseline_detection_male` / `baseline_encounter_rate_male`, `sigma_male`.
`design`: `trap_grid_side` (square grid, side^2 detectors), `trap_spacing`
(m), `n_occasions`, `buffer` (m), `mask_spacing` (m), `obs_model`
(binomial | bernoulli | poisson), `density_covariate`, `density_covariate2`
(elevation), `nonhabitat` (share of the state space that is a lake),
`detection_covariate`, `prop_on_trail`, `occasion_covariate` (rain),
`sigma_covariate` (lure), `behaviour`, `sex`, `prop_sex_unknown`,
`continuous_effort` (effort 0.5-1.5 per detector-occasion), `prop_inactive`,
`trap_failure` (share of detectors that stop working part-way),
`covariate_range` (m), `covariate_trend` (> 0 puts the array in good habitat).

Useful knobs: `trap_spacing` relative to `sigma` (<= 2 sigma),
`n_occasions`, `trap_grid_side`, `density`, `prop_sex_unknown`:

    scripts/Rscript core/run.R power scr --config my_design.yaml --vary n_occasions=3,4,6,8 --nsim 50

## Pitfalls

- Latitude/longitude coordinates: project first (UTM). Polygons and rasters
  must use the same projected CRS as the detectors.
- Occasions where a detector was not working must be 0 in the usage/effort
  table; otherwise p0 is underestimated and density overestimated.
- Do not pool sexes with very different ranging; fit sex-specific detection.
- Too small a buffer inflates density: use >= 3-4 sigma of the
  wider-ranging sex and run `buffer_check: true` on the final model.
- Density covariates are standardised over the habitat in the state space:
  a different buffer changes the SD and therefore the per-SD coefficient
  (the per-unit rate ratio in the interpretation does not change).
- Abundance refers to the habitat in the state space you chose; report
  density, which does not depend on the buffer once the buffer is adequate.
- Genotyping errors: see the NGS section; the model assumes correct IDs.
- Many parameters (sex x covariates x behaviour) need longer chains: check
  R-hat < 1.01 and ESS >= 400 in the report and raise `niter` if not.

## References

- Efford M (2004) Density estimation in live-trapping studies. Oikos 106:598-610.
- Borchers DL, Efford MG (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64:377-385.
- Royle JA, Young KV (2008) A hierarchical model for spatial capture-recapture data. Ecology 89:2281-2289.
- Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press. (ch. 9 alternative observation models: Bernoulli, Poisson; ch. 11 density covariates.)
- Royle JA, Gardner B (2011) Hierarchical spatial capture-recapture models for estimating density from trapping arrays. In: Camera Traps in Animal Ecology. Springer, pp 163-190.
- Efford MG, Borchers DL, Mowat G (2013) Varying effort in capture-recapture studies. Methods in Ecology and Evolution 4:629-636.
- Gardner B, Royle JA, Wegan MT (2009) Hierarchical models for estimating density from DNA mark-recapture studies. Ecology 90:1106-1115.
- Gardner B, Royle JA, Wegan MT, Rainbolt RE, Curtis PD (2010) Estimating black bear density using DNA data from hair snares. Journal of Wildlife Management 74:318-325.
- Pledger S (2000) Unified maximum likelihood estimates for closed capture-recapture models using mixtures. Biometrics 56:434-442.
- Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152:S521-S537.
- Sollmann R, Gardner B, Belant JL (2012) How does spatial study design influence density estimates from spatial capture-recapture models? PLoS ONE 7:e34575.
- Sun CC, Fuller AK, Royle JA (2014) Trap configuration and spacing influences parameter estimates in spatial capture-recapture models. PLoS ONE 9:e88025.
- Milleret C et al. (2019) A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363.
- Taberlet P, Waits LP, Luikart G (1999) Noninvasive genetic sampling: look before you leap. Trends in Ecology & Evolution 14:323-327.
- Paetkau D (2003) An empirical exploration of data quality in DNA-based population inventories. Molecular Ecology 12:1375-1387.
- Lukacs PM, Burnham KP (2005) Review of capture-recapture methods applicable to noninvasive genetic sampling. Molecular Ecology 14:3909-3919.
- Lukacs PM, Burnham KP (2005) Estimating population size from DNA-based closed capture-recapture data incorporating genotyping error. Journal of Wildlife Management 69:396-403.
- Link WA, Yoshizaki J, Bailey LL, Pollock KH (2010) Uncovering a latent multinomial: analysis of mark-recapture data with misidentification. Biometrics 66:178-185.
- Wright JA, Barker RJ, Schofield MR, Frantz AC, Byrom AE, Gleeson DM (2009) Incorporating genotype uncertainty into mark-recapture-type models for estimating abundance using DNA samples. Biometrics 65:833-840.
- Efford MG secr: Spatially explicit capture-recapture models. R package (5.x used). doi:10.32614/CRAN.package.secr
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26:403-413.
- Phan D, Pradhan N, Jankowiak M (2019) Composable effects for flexible and accelerated probabilistic programming in NumPyro. arXiv:1912.11554.
