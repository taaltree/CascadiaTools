---
id: scr-smr
title: Spatial mark-resight (SMR) and generalized spatial mark-resight (gSMR)
when_to_use: Only part of the population carries individual marks (collars, ear tags, natural marks on a subset) and cameras (or other resighting detectors) photograph both marked and unmarked animals; question is density (per km2), abundance, the marked fraction, ranging scale (sigma) or habitat effects on density. Use gsmr when the animals were marked at their own traps during the study (live-trapping, darting at bait sites), so marked animals are concentrated near the marking sites; use smr when the marked animals are a random sample of the population (e.g. collared earlier across the whole area). Optional telemetry (GPS/VHF fixes) of marked animals.
not_for: Every animal individually identifiable (use scr), no marked animals at all (unmarked SCR / spatial count models; density is weakly identified), open populations across seasons (scr-open), mark loss or marks that change between sessions, animals that cannot be told marked vs unmarked on a photo at all.
data_needed: camera CSV (id, x, y [, usage per occasion, camera covariates]); sightings CSV (camera, occasion, status marked / unmarked / marked-unidentified, animal ID for identified marks); either a list of marked animals (SMR) or the marking traps + marking captures (gSMR); optional telemetry CSV (animal ID, x, y); optional habitat grid CSV or raster with density covariates
engine: nimble (Bayesian data augmentation with latent activity centres; default) and secr (maximum pseudo-likelihood mark-resight; independent cross-check)
---

# scr-smr: spatial mark-resight and generalized spatial mark-resight

## Question it answers

*How many animals live here, at what density, when only some of them can
be recognised?* Cameras photograph marked animals (whose identity is read
from a collar, ear tag or natural mark) and unmarked animals (counted, but
not identified). Spatial mark-resight (Chandler & Royle 2013; Sollmann et
al. 2013) uses where each marked animal was resighted to learn how far
animals range (sigma) and how often a camera records an animal, then uses
the **spatial pattern of the unmarked counts** across cameras to estimate
how many unmarked animals there are and where they live.

Generalized SMR (Whittington, Hebblewhite & Chandler 2018; Jimenez et al.
2019) adds the **marking process** itself: animals were marked at traps,
so animals living near the traps are more likely to carry marks. Ignoring
this makes marked animals look like a random sample when they are not, and
biases density (Whittington et al. found conventional SMR gave 7.4
grizzly bears per 1000 km2 against 12.4 for gSMR). This module fits both,
so you can see the difference on your own design.

| estimand | meaning |
|---|---|
| `density_per_km2` | realised animals per km2 in the state space (`density_per_unit2` if map units are unknown) |
| `abundance` | realised number of activity centres (marked + unmarked) in the state space (camera array plus buffer) |
| `abundance_unmarked` | abundance minus the number of marked animals |
| `marked_fraction` | m / N, the share of the population carrying marks |
| `sigma` | spatial scale of the half-normal detection function, map units (m) |
| `home_range_area_95` | pi (2.448 sigma)^2, 95% area of a circular bivariate-normal home range (km2); an index of ranging scale |
| `resight_rate` | lambda0: expected sightings per occasion of an animal by a camera at its activity centre (Poisson model, reference level of camera covariates). `resight_detection` (p0) for the Bernoulli model |
| `id_probability` | delta: probability that a sighting of a marked animal can be identified to individual (only when marked-but-unidentified sightings are recorded) |
| `marking_detection` | gSMR only: per-occasion capture probability at a marking trap at the animal's activity centre |
| `density_beta_<cov>` | log rate ratio of density per 1 SD of a habitat covariate (log-linear density) |
| `detection_beta_<cov>` | camera covariate effect on lambda0 (log scale, Poisson) or p0 (logit scale, Bernoulli); 0/1 covariates: level 1 vs 0 |

## Model

State space S: habitat pixels (a regular grid; camera and marking traps
plus a buffer). Activity centres s_i of the N animals follow a Poisson
point process with intensity D(x) = exp(b0 + b' X(x)) (homogeneous if
there are no density covariates). Fitted by data augmentation (Royle &
Dorazio 2012): M rows, z_i ~ Bernoulli(psi) says whether row i is a real
animal, N = sum z_i.

Resighting (camera j, operating K_j occasions; distance d_ij from s_i):

    lambda_ij = lambda0_j exp(-d_ij^2 / (2 sigma^2)),   log lambda0_j = log lambda0 + bP' W_j
    identified sightings of marked i at j:  y_ij ~ Poisson(K_j delta lambda_ij)
    unmarked sightings at camera j:         nU_j ~ Poisson(K_j sum_{i unmarked} z_i lambda_ij)
    marked-unidentified sightings at j:     nTm_j ~ Poisson(K_j (1 - delta) sum_{i marked} lambda_ij)

These are exact for the Poisson resight model (Poisson thinning and
superposition: Chandler & Royle 2013, Sollmann et al. 2013, Jimenez et al.
2019). With `resight_model: bernoulli` (detected / not per occasion)
identified sightings are Binomial(K_j, delta p_ij) and the pooled counts
use the Poisson approximation of Chandler & Royle (2013).

Marking:

- **SMR** (`model: smr`): the m marked animals are a known set included
  with certainty; psi applies to the remaining rows. Marked animals are
  assumed to be a spatially random sample of the population, and their
  activity centres follow the same point process.
- **gSMR** (`model: gsmr`): every row may be caught at the L marking traps
  over Km occasions, x_il ~ Binomial(Km_l, pm0 exp(-d_il^2 / (2 sigma^2)))
  (an all-zero history for animals never caught = unmarked animals). Marked
  animals are therefore modelled as the animals that happened to be caught,
  not as a random sample (Whittington et al. 2018).
- **Telemetry** (optional): each fix of marked animal i ~ bivariate normal
  N(s_i, sigma^2 I) (Sollmann et al. 2013), entered through sufficient
  statistics. Fixes must be roughly independent (thin GPS data to one fix
  every day or two).

Priors: psi ~ U(0,1); sigma ~ U(0, half the state-space width); lambda0 ~
U(0, 10) (p0 ~ U(0,1) for Bernoulli); delta ~ U(0,1); pm0 ~ U(0,1); all
covariate coefficients ~ Normal(0, SD 2) (covariates standardised; 0/1
covariates left as 0/1).

MCMC: activity centres and inclusion indicators are updated by a custom
joint sampler that caches each camera's pooled expected count, so a sweep
over all M rows costs O(M J) instead of O(M^2 J). Rows outside the
population (z = 0) get a uniform pseudo-prior for s, which leaves the
posterior unchanged but stops them from slowing the density coefficients.
Top-level parameters: adaptive random-walk Metropolis plus a block sampler.
If the posterior of N approaches the augmentation limit M, the model is
refitted automatically with a larger M.

secr engine (`engine: secr`, or `compare_engines: [secr]`): Efford &
Hunter's (2018) maximum pseudo-likelihood mark-resight model in secr, with
sighting-only occasions (SMR) or marking + sighting occasions (gSMR),
half-normal detection, Poisson counts, pID for unidentified marks, and
variances inflated for the overdispersion of the unmarked counts (c-hat by
simulation, two-step fit). It estimates **expected** N (Bayesian:
realised N), ignores telemetry, uses a logit link for g0 (so camera
covariate coefficients are on a different scale) and supports only the
Poisson resight model.

## Assumptions (and what breaks if they fail)

- **Closure** during resighting (and marking, for gSMR): no births,
  deaths, immigration or emigration. Violations inflate N and sigma. Keep
  the survey short (a few months for mammals).
- **Marked and unmarked animals behave alike**: same sigma and same
  sighting rate. Trap-happy or collar-shy animals bias density. Collars
  sometimes change detectability; check with the posterior predictive plot
  of sightings per marked animal.
- **Every marked animal present is listed** (SMR), including those never
  photographed, and none has left or lost its mark. A marked animal that
  died or left before resighting must be removed (`marked_in_area_col`).
  Missing a marked animal underestimates the marked fraction and N.
- **Marked animals are a spatially random sample (SMR only).** If marking
  was concentrated near some locations, use gSMR with the marking captures.
  The negative control in VALIDATION.md shows the size of the bias.
- **Independent sighting events.** A burst of photos of one animal must be
  one event (use an independence interval, e.g. 30 min). Counting the same
  pass several times inflates both lambda0 and the unmarked counts in
  different ways, which biases N.
- **Unmarked animals are counted, not identified.** Two unmarked animals
  in one photo count as 2 (use `count_col`).
- **Status is scored correctly.** A marked animal whose tag is not visible
  in the photo looks unmarked. If that happens often, the marked fraction
  is underestimated and N inflated. Only score "unmarked" when the body
  part carrying the mark is clearly visible, and record "marked,
  unidentified" when the mark is seen but cannot be read.
- **Half-normal space use around fixed activity centres**; state space
  large enough (buffer >= 3 sigma); cameras spaced <= 2 sigma apart (Sun
  et al. 2014).

## Data format

Cameras (`cameras.csv`), one row per camera:

    camera_id,utm_e,utm_n,on_trail,weeks_active
    FC01,516000,5216000,1,1111111111
    FC02,517500,5216000,0,1111111000

Sightings (`fisher_photos.csv`), one row per independent photo event. The
status column says whether the animal was marked and identified, marked
but unreadable, or unmarked. The animal ID is filled only for identified
marks. Extra columns are ignored.

    camera_id,week,date,tag_status,fisher_id,photo
    FC05,1,2025-06-03,ear-tag read,FI-022,FC05_20250603_0947.JPG
    FC01,1,2025-06-03,no tag,,FC01_20250603_0189.JPG
    FC12,2,2025-06-11,"ear-tag, number unreadable",,FC12_20250611_1203.JPG

SMR, list of all marked animals present during resighting (`marked.csv`):

    animal_id
    M001
    M002

gSMR, marking traps (`box_traps.csv`: trap_id, x, y) and marking captures
(`box_trap_captures.csv`: animal ID, occasion, trap ID; one row per
capture). Every animal in the capture file is a marked animal.

Telemetry (`gps_fixes.csv`, optional): animal ID, x, y (one row per fix;
same coordinate system as the cameras).

Habitat grid (`habitat_grid.csv`, optional; required for density
covariates): regular grid of pixel centres x, y covering the detectors
plus a buffer of >= 3 sigma, with covariate columns and an optional 0/1
habitat column. Aim for 800-2500 pixels (spacing about sigma / 2). A
raster (`raster:`, read with terra) can supply covariates instead.

## Config

```yaml
data:
  species: fishers                 # used in sentences
  units: m                         # m | km (omit if coordinates are not metric)
  cameras: cameras.csv
  camera_id_col: camera_id         # default camera_id
  x_col: utm_e                     # default x
  y_col: utm_n                     # default y
  usage_col: weeks_active          # optional "1101..." per camera; else n_occasions
  n_occasions: 10                  # needed if no usage_col
  detection_covariates: [on_trail] # camera columns
  sightings: fisher_photos.csv
  camera_col: camera_id            # default camera_id
  occasion_col: week               # default occasion (1..K)
  status_col: tag_status           # default status
  marked_values: ["ear-tag read"]  # default marked / M / marked_identified
  unmarked_values: ["no tag"]      # default unmarked / U
  unidentified_values: ["ear-tag, number unreadable"]  # default marked_unidentified / MU
  animal_col: fisher_id            # default animal_id
  count_col: n_animals             # optional: animals per sighting row (group photos)
  unidentified_recorded: true      # default: TRUE if any unidentified rows
  # SMR:
  marked: marked.csv               # every marked animal present
  marked_id_col: animal_id
  marked_in_area_col: alive        # optional 0/1 filter (collar working, animal present)
  # gSMR:
  marking_traps: box_traps.csv
  marking_trap_id_col: trap_id
  marking_occasions: 5             # or marking_usage_col
  marking_captures: box_trap_captures.csv
  marking_animal_col: fisher_id
  marking_occasion_col: night
  marking_trap_col: trap_id
  # optional:
  telemetry: gps_fixes.csv
  telemetry_animal_col: fisher_id
  telemetry_x_col: utm_e
  telemetry_y_col: utm_n
  mask: habitat_grid.csv           # or raster: covariates.tif
  habitat_col: habitat
  density_covariates: [canopy_cover_pct]
  mask_spacing: 500                # optional (inferred); buffer: 4000 when no mask given
  labels: {canopy_cover_pct: "Canopy cover (%)"}
settings:
  model: auto                      # auto (gsmr if marking captures given) | smr | gsmr
  resight_model: poisson           # poisson | bernoulli
  engine: nimble                   # nimble | secr
  compare_engines: [secr]          # also fit secr and show both
  niter: 30000                     # per chain (thin 4, 6000 burn-in, 3 chains)
  nburnin: 6000
  nchains: 3
  thin: 4
  M_aug: null                      # augmentation rows for unmarked animals (auto; refits if too small)
  use_unidentified: true           # use marked-unidentified counts
  use_telemetry: true
  surface_draws: 300               # draws for density surfaces
  ppc_draws: 300                   # draws for posterior predictive checks
  secr_nsim: 50                    # simulations for secr's c-hat
```

## Outputs

`fit` writes the estimates table, an HTML report (methods, interpretation,
design notes, figures) and the figure suite (`figures/`): `density_surface`
(realised density map with cameras, marking traps, marked activity
centres), `density_expected_surface`, `density_<cov>` (conditional effect
on the response scale), `detection_function` (resighting and, for gSMR,
marking, with camera spacing marked), `abundance_posterior`,
`marked_fraction`, `captures_per_individual` (observed vs posterior
predictive), `unmarked_by_camera`, `gof` (posterior predictive checks),
`coefficients`, `engine_comparison` (when secr is also fitted), `trace`,
`posterior`.

Diagnostics notes report the design checks (camera spacing / sigma,
buffer / sigma, number of resighted marked animals), posterior predictive
p-values and convergence against the real-data standard (R-hat < 1.01,
bulk ESS >= 400).

## Simulation and power

`simulate()` knobs (`design:`): `model` (smr | gsmr), `resight_model`,
`n_marked` (SMR), `trap_grid_side`, `trap_spacing`, `n_occasions`,
`buffer`, `mask_spacing`, `density_covariate`, `detection_covariate`,
`prop_on_trail`, `trap_failure`, `covariate_range`, `telemetry_fixes`,
`telemetry_fraction`, `marking_grid_side`, `marking_spacing`,
`marking_occasions`, `marking_offset` (shift of the marking grid from the
camera-grid centre). `params:` `density` (per km2), `resight_rate`,
`sigma` (m), `id_probability`, `marking_detection`,
`density_beta_canopy`, `detection_beta_on_trail`. Sensible power knobs:
`n_marked`, `n_occasions`, `trap_grid_side`, `marking_occasions`.

## Pitfalls

- **Few marked animals.** With 5 marked animals of ~40, sigma and lambda0
  rest on very few animals; see VALIDATION.md (harder design) for how
  intervals and coverage behave. Marking more animals is worth more than
  more camera weeks (Chandler & Royle 2013).
- **Conventional SMR after spatially concentrated marking** is biased.
  Use gSMR (supply marking captures) whenever marking happened at traps.
- **Mis-scored status** (marks not visible) is not modelled; see
  Assumptions.
- **secr vs Bayesian**: secr reports expected N and treats unmarked counts
  as overdispersed independent counts. For conventional SMR the two agreed
  closely in validation (mean ratio 1.02 for N). For **gSMR**, secr's joint
  marking + sighting fit was biased low in our simulations (-26% N, coverage
  0.60) while the Bayesian gSMR was unbiased (VALIDATION.md section 8), so a
  lower secr N in a gSMR comparison is expected and is not evidence against
  the Bayesian fit; compare sigma and detection rates instead.
- **Run length**: the default (30,000 iterations, 3 chains) meets R-hat <
  1.01 / ESS >= 400 on the example. Longer runs are needed with few marked
  animals; the report says so when the standard is not met.

## References

- Chandler RB, Royle JA (2013) Spatially explicit models for inference about density in unmarked or partially marked populations. Annals of Applied Statistics 7:936-954. doi:10.1214/12-AOAS610
- Sollmann R, Gardner B, Parsons AW, Stocking JJ, McClintock BT, Simons TR, Pollock KH, O'Connell AF (2013) A spatial mark-resight model augmented with telemetry data. Ecology 94:553-559. doi:10.1890/12-1256.1
- Whittington J, Hebblewhite M, Chandler RB (2018) Generalized spatial mark-resight models with an application to grizzly bears. Journal of Applied Ecology 55:157-168. doi:10.1111/1365-2664.12954
- Jimenez J, Chandler R, Tobajas J, Descalzo E, Mateo R, Ferreras P (2019) Generalized spatial mark-resight models with incomplete identification: an application to red fox density estimates. Ecology and Evolution 9:4739-4748. doi:10.1002/ece3.5077 (corrigendum: use sum(), not dsum(), in JAGS)
- Efford MG, Hunter CM (2018) Spatial capture-mark-resight estimation of animal population density. Biometrics 74:411-420. doi:10.1111/biom.12766
- Efford MG (2025) Mark-resight in secr 5.2 (package vignette secr-markresight.pdf); secr R package.
- Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press (ch. 18 unmarked populations, ch. 19 spatial mark-resight).
- Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(Suppl 2):S521-S537.
- Sun CC, Fuller AK, Royle JA (2014) Trap configuration and spacing influences parameter estimates in spatial capture-recapture models. PLoS ONE 9:e88025.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26:403-413.
