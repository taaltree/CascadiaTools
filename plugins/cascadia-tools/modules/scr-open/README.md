---
id: scr-open
title: Open-population spatial capture-recapture (Jolly-Seber SCR)
when_to_use: Individually identified animals (hair-snare or scat DNA, camera traps of marked/spotted species, live traps) detected at fixed detectors over repeated occasions in each of several seasons or years (robust design); question is abundance/density in each year, trend, apparent survival, recruitment, or whether density is associated with habitat.
not_for: One closed season only (use scr); non-spatial capture histories without detector locations (use jolly-seber or cjs); survival of marked animals only, with no interest in abundance (cjs); animals whose ranges shift a lot between years (static activity centres assumed; see Pitfalls); unmarked animals (nmixture, occupancy, distance).
data_needed: detector CSV (id, x, y, optional per-year operation columns); detection CSV (individual, year/period, occasion within period, detector); optional habitat grid CSV with density covariates
engine: nimble (Bayesian; marginal [default], hmm or latent data augmentation), numpyro (Bayesian NUTS, CPU/GPU); openCR (maximum likelihood) as an independent cross-check
---

# scr-open: open-population spatial capture-recapture

## Question it answers

*How many animals are there each year, is the population growing or
declining, how many survive from year to year and how many new animals
arrive, and is density higher in some habitats?* The module combines a
Jolly-Seber model (animals enter, survive and die between primary periods)
with spatial capture-recapture (where each animal was detected tells the
model where it lives and how detectable it is). It estimates:

| estimand | meaning |
|---|---|
| `abundance_period_<t>` | realised number of animals alive in the state space (detectors plus buffer) in period t |
| `density_per_km2_period_<t>` | abundance / state-space area (`density_per_unit2_...` if map units are unknown) |
| `superpopulation_size` | number of different animals alive in the state space at any time during the study |
| `recruits_period_<t>` | new animals entering between period t-1 and t (births plus immigrants that settle) |
| `per_capita_recruitment_period_<t>` | recruits_t / abundance_{t-1} |
| `growth_rate_period_<t>` | abundance_t / abundance_{t-1} (realised lambda) |
| `survival` | apparent survival between consecutive periods (death and permanent emigration are confounded) |
| `baseline_detection` | p0: per-occasion detection probability for a detector at an animal's activity centre |
| `sigma` | spatial scale of the half-normal detection function (map units, e.g. m) |
| `home_range_area_95` | pi (2.448 sigma)^2, 95% area of a circular bivariate-normal home range (km2) |
| `entry_prob_period_<t>` | gamma_t, data-augmentation entry probability (technical; not a biological rate) |
| `density_beta_<cov>` | log rate ratio of density per 1 SD of a habitat covariate across the state space |

## Model

Robust design: T primary periods (years), each with K secondary
occasions; the population is closed within periods and open between them.
Data augmentation to M potential individuals (Royle & Dorazio 2012). Each
has a state in each period: 1 = not yet entered, 2 = alive, 3 = dead.

    period 1:     P(alive) = gamma_1
    period t > 1: not entered -> alive with prob gamma_t (else stays not entered)
                  alive -> alive with prob phi, -> dead with prob 1 - phi
                  dead stays dead

(the restricted dynamic-occupancy form `z[i,t] ~ Bern(phi z[i,t-1] +
gamma_t A[i,t])` of Kery & Schaub 2012 ch. 10 and Royle & Dorazio 2008).
Each animal has one activity centre s_i, fixed across periods (Gardner et
al. 2010), distributed over a habitat mask of pixels g with

    P(s_i in g) = exp(sum_k beta_k x_k(g)) / sum_h exp(sum_k beta_k x_k(h))   (x standardised over the mask)

While alive in period t, the number of occasions animal i was detected at
detector j is

    y_ijt ~ Binomial(K_jt, p0 exp(-d(s_i, x_j)^2 / (2 sigma^2)))     K_jt = occasions detector j worked in period t

("proximity" detectors: at most one detection per detector per occasion,
several detectors per occasion allowed).

Engines (`settings$engine`, `settings$approach`):

- **nimble, `marginal` (default).** Both the alive-state chain and the
  activity centre are summed out for every augmented animal in one compiled
  custom distribution (`dOSCRmarg`): states with the forward algorithm
  (Turek et al. 2016), the centre over mask pixels; pixels where a detected
  animal's seen-period likelihood is below exp(-30) of its best pixel are
  skipped (local evaluation, Milleret et al. 2019; exact to < 1e-13). Only
  p0, sigma, phi, gamma_1..T and the betas are sampled (adaptive block
  random walk). Realised N_t, recruits, the superpopulation and density maps
  are then drawn exactly from their conditional posterior for each
  retained draw: each animal's pixel from P(g | history), then its state path
  by forward-filtering backward-sampling.
- **nimble, `hmm`.** Alive states summed out per animal (`dOSCRhmm`, the
  benchmark's `dOpenSCR`), activity centres sampled (block random walk);
  N_t by the same FFBS step.
- **nimble, `latent`.** Everything sampled (z[i,t] binary, s_i block RW);
  the textbook BUGS formulation, slowest; monitors N_t, B_t and the
  superpopulation directly.
- **numpyro.** The same fully marginal likelihood with NUTS in JAX
  (`python/scro_numpyro.py`). CPU by default; set `numpyro_gpu: true` and
  run under `scripts/ejob -g` for big arrays.
- **openCR** (`compare_engines: [openCR]` only): maximum-likelihood spatial
  JSSA (`JSSAsecrf`, Efford & Schofield 2020) on the identical mask. Its
  detection function is the hazard half-normal (HHN), so p0/sigma are on a
  slightly different scale; N_t, phi and per-capita recruitment f_t are
  directly comparable.

All Bayesian engines are ported from the SCR GPU benchmark in the parent
repository (`../python/simulate_open.py`, `../python/fit_numpyro_open.py`,
`../R/fit_nimble_open.R`, cross-validated there), extended with per-period
detector operation, irregular masks and named covariates.

Priors: p0, phi ~ Uniform(0, 1); entry in the superpopulation
parameterisation (Kery & Schaub 2012 sec. 10.3.3): psi ~ Uniform(0, 1) (share of
the M augmented rows ever alive, so Nsuper is uniform on 0..M) and entry
probabilities b ~ Dirichlet(1, ..., 1); gamma_1 = psi b_1, gamma_t = psi b_t /
(1 - psi sum_{k<t} b_k). sigma ~ Uniform(0, half the state-space width); density
coefficients ~ Normal(0, SD 2) on standardised covariates. (Until 2026-10-20 the
gamma_t had independent Uniform(0, 1) priors; that implies a prior on Nsuper / M
with median 0.97 at T = 4 and biased abundance upward when the data are weak;
see VALIDATION.md, "Independent verification".)

## Assumptions (and what breaks)

- **Closure within periods** (keep each season short relative to turnover)
  and **turnover only between periods**.
- **Static activity centres.** An animal keeps one centre for the whole
  study. Animals that shift ranges between years look like deaths plus
  recruits (survival biased low) and inflate sigma (Efford & Schofield
  2020; Glennie et al. 2019). Use openCR with a movement model if shifts are
  large.
- **Homogeneous survival and detection** across animals and periods
  (constant p0, sigma, phi). Sex/age differences, trap-happiness, transients
  and year-to-year detection changes violate this; the posterior predictive
  checks (`detection_by_period`, `gof`) are there to flag it.
- **Apparent survival**: permanent emigration counts as death; recruitment
  includes immigration.
- **State space large enough**: buffer >= 3-4 sigma beyond the detectors
  (VALIDATION.md negative control: buffer = 1 sigma inflates density).
- **Detector spacing** <= about 2 sigma (Sun et al. 2014).
- **Enough data**: the default design (~135 animals detected over 5
  periods, ~60 per period, ~80 in two or more periods, ~1,200
  detections) recovers everything to within 2% (abundance) to 5%
  (per-capita recruitment). With ~40 detected animals (~20 per period)
  abundance and the superpopulation are 4-8% high and per-capita
  recruitment 12-20% high (posterior means; intervals still cover;
  VALIDATION.md, "Independent verification"). With few animals detected in a
  period, that period's abundance leans on the survival/recruitment
  model; the report warns below 10.
- **M large enough**: default M = max(6n, n + 300); refitted with 2M
  automatically if the superpopulation's 99.5% quantile reaches 0.9 M.

## Data format

`hair_snares.csv` (one row per detector; `ops_<year>`: 0/1 per occasion,
or the number of occasions worked):

    snare_id,utm_e,utm_n,elevation_m,ops_2020,ops_2021,ops_2022,ops_2023,ops_2024
    HS01,574000,4864000,1351,11111,11111,11111,11111,11111
    HS07,580000,4864000,1206,11111,11100,11111,11111,11111

`marten_genotypes.csv` (one row per detection/sample; extra columns
ignored; several samples of one animal at one detector on one occasion
count once):

    sample_id,year,check,snare_id,marten_id,n_hairs
    S42482,2020,1,HS01,MACA-476,8

`habitat_grid.csv` (optional; regular grid of pixel centres covering the
detectors plus a buffer, with density covariates and an optional 0/1
habitat flag):

    x,y,canopy_cover_pct
    571200,4861200,49.2

Without a grid the mask is built from the detectors with a buffer of
4 x the rough sigma from repeat detections and pixels of sigma/2
(homogeneous density only). Coordinates must be projected (metres, e.g.
UTM), not latitude/longitude. Occasions are numbered 1..K **within each
period**.

## Config

```yaml
data:
  traps: hair_snares.csv          # required
  trap_id_col: snare_id           # default id
  x_col: utm_e                    # default x
  y_col: utm_n                    # default y
  usage_cols: [ops_2020, ops_2021, ...]   # optional, one per period: 0/1 string or count
  captures: marten_genotypes.csv  # required
  individual_col: marten_id       # default individual
  period_col: year                # default period
  occasion_col: check             # default occasion (1..K within period)
  trap_col: snare_id              # default trap
  periods: [2020, 2021, 2022, 2023, 2024]  # order of periods; REQUIRED if a period had no detections
  n_occasions: 5                  # K, required
  units: m                        # m or km -> density per km2; omit if unknown
  mask: habitat_grid.csv          # optional habitat grid
  mask_x_col: x                   # defaults x, y
  habitat_col: habitat            # optional 0/1 flag
  buffer: 3000                    # if the mask is built here (map units)
  mask_spacing: 400
  density_covariates: [canopy_cover_pct]
  labels: {canopy_cover_pct: "Canopy cover (%)"}
settings:
  engine: nimble                  # nimble | numpyro
  approach: marginal              # marginal | hmm | latent (nimble)
  niter: 30000                    # nimble marginal, per chain
  nburnin: 10000
  thin: 2
  nchains: 3
  top_sampler: RW_block           # or AF_slice (nimble marginal)
  niter_latent: 30000             # nimble hmm / latent
  nburnin_latent: 10000
  thin_latent: 5
  numpyro_warmup: 1000
  numpyro_samples: 1500
  numpyro_gpu: false              # true with scripts/ejob -g
  max_draws_per_chain: 500        # draws per chain used for N_t, recruits (FFBS)
  M: null                         # augmentation size (auto)
  buffer: null                    # trim the state space to this distance from detectors
  sigma_max: null                 # upper limit of the sigma prior (map units)
  density_covariates: null        # override (e.g. [] = homogeneous density)
  compare_engines: [numpyro, openCR]   # also fit with these; figure engine_comparison
  ppc_draws: 300
  surface_draws: 200
```

## Outputs

The estimates table (above), diagnostics (R-hat, bulk ESS; the real-data
standard is R-hat < 1.01 and bulk ESS >= 400 — the report says when it is
not met), design notes (spacing/sigma, buffer/sigma, between-period
recaptures), posterior predictive p-values, and figures:
`abundance_trajectory`, `recruitment` (recruits and per-capita rate),
`growth_rate`, `survival` (posterior of the constant survival),
`density_surface` (realised density map per period), `density_<cov>`
(expected density vs covariate, response scale), `detection_function`,
`detection_by_period` (observed vs predicted individuals detected per
period), `gof`, `coefficients`, `engine_comparison`, `trace`, `posterior`.

## Simulation and power

`params`: `density` (animals per km2 alive in period 1), `baseline_detection`,
`sigma` (m), `survival`, `per_capita_recruitment` (expected recruits per
animal alive in the previous period; simulated as Poisson), `density_beta_canopy`.
`design`: `trap_grid_side` (square grid, side^2 detectors), `trap_spacing`
(m), `n_occasions` (K), `n_periods` (T), `buffer` (m), `mask_spacing` (m),
`density_covariate` (TRUE/FALSE), `covariate_range` (smoothing SD of the
simulated canopy field, m), `trap_failure` (per detector and period, chance
it stops working from a random occasion on).

Defaults: 7 x 7 detectors 1,000 m apart (1.25 sigma), 5 periods x 5
occasions, 3,000 m buffer, 400 m pixels (900), density 0.8 / km2 in
period 1 (~115 animals; ~235 over the study, ~135 detected), p0 = 0.3,
sigma = 800 m, survival 0.8, recruitment 0.25, canopy effect 1.0 per SD.
Useful power knobs: `n_periods`, `trap_grid_side`, `n_occasions`,
`density_beta_canopy`. A 100-replicate recovery takes about an hour on 16
cores (~8.5 min per replicate) (see VALIDATION.md for timings).

## Pitfalls

- Give `periods:` explicitly; a year with zero detections would otherwise
  silently vanish and the model would treat the next year as consecutive.
- Unequal intervals between periods are not modelled: survival is per
  period interval. If you skipped a year, survival across that gap is phi^2
  in reality but phi here; fit only consecutive periods or treat results
  for the gap with care.
- `entry_prob_period_<t>` depends on M; interpret recruits and per-capita
  recruitment instead.
- Period-1 recruitment is not identifiable (everything alive in period 1
  counts as "entered"); recruits start at period 2.
- Survival in the last interval and recruitment in the last period are
  less precise (no later periods to confirm presence).
- Density covariate effects are spatial associations, not causal effects.
- **Traps placed in unrepresentative habitat** (e.g. along roads in open
  forest when animals prefer dense canopy): total abundance then rests on
  extrapolating the covariate effect to habitat the traps never sampled. The
  report gives the share of animals living more than 2 sigma from any trap
  and the ratio of state-space-average to trapped-area density, and warns when
  that ratio is above 1.5 or below 1/1.5. The bundled example triggers this
  warning (ratio 2.1); see VALIDATION.md section 8.

## References

- Gardner B, Reppucci J, Lucherini M, Royle JA (2010) Spatially explicit inference for open populations: estimating demographic parameters from camera-trap studies. Ecology 91:3376-3383. doi:10.1890/09-0804.1
- Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press (ch. 16, open population models).
- Efford MG, Schofield MR (2020) A spatial open-population capture-recapture model. Biometrics 76:392-402. doi:10.1111/biom.13150
- Glennie R, Borchers DL, Murchie M, Harmsen BJ, Foster RJ (2019) Open population maximum likelihood spatial capture-recapture. Biometrics 75:1345-1355. doi:10.1111/biom.13078
- Schwarz CJ, Arnason AN (1996) A general methodology for the analysis of capture-recapture experiments in open populations. Biometrics 52:860-873.
- Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. J. Ornithology 152:S521-S537.
- Royle JA, Dorazio RM (2008) Hierarchical Modeling and Inference in Ecology. Academic Press.
- Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS. Academic Press (ch. 10, Jolly-Seber via data augmentation).
- Borchers DL, Efford MG (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64:377-385.
- Turek D, de Valpine P, Paciorek CJ (2016) Efficient Markov chain Monte Carlo sampling for hierarchical hidden Markov models. Environ. Ecol. Stat. 23:549-564.
- Milleret C et al. (2019) A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363.
- Sun CC, Fuller AK, Royle JA (2014) Trap configuration and spacing influences parameter estimates in spatial capture-recapture models. PLoS ONE 9:e88025.
- Efford MG (2024) openCR: Open population capture-recapture. R package 2.2.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26:403-413.
- Phan D, Pradhan N, Jankowiak M (2019) Composable effects for flexible and accelerated probabilistic programming in NumPyro. arXiv:1912.11554.
