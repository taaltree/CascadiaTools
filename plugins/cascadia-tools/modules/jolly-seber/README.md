---
id: jolly-seber
title: Jolly-Seber open-population capture-recapture (abundance, survival, recruitment)
when_to_use: Individually marked animals (tags, PIT tags, bands, photo-ID, genetic ID) captured on 3+ occasions (usually years or seasons) in an open population, where the question is how many animals were present at each occasion, how many entered (births plus immigrants), how many survived, and whether the population is growing or declining.
not_for: Survival only, conditioning on first capture (use cjs); a single closed season (closed capture-recapture or scr); density per km2 or spatial capture histories (scr / scr-open); unmarked counts (nmixture / timeseries); populations with strong temporary emigration (needs a robust design, not provided here).
data_needed: one row per capture (animal ID + occasion, optional sex/group), or one capture-history string per animal; optional occasion table with interval lengths and occasion covariates for capture probability (effort)
engine: nimble (Bayesian data augmentation, restricted dynamic occupancy formulation; marginal forward-algorithm likelihood by default, latent z sampler as an option); marked (maximum-likelihood POPAN) as an independent cross-check
---

# jolly-seber: abundance, recruitment and survival from capture-recapture

## Question it answers

"How many turtles lived in the pond each year? How many new animals joined
the population (births plus immigrants) between years? What share survived and
stayed? Do males survive less well than females? Is the population declining?"

A Jolly-Seber (JS) model uses the full capture history of every marked animal,
including the occasions **before** its first capture. That lets it estimate
abundance and entry as well as survival. A Cormack-Jolly-Seber model (`cjs`)
conditions on first capture, so it estimates only survival.

## Model

The model is the superpopulation formulation of Crosbie & Manly (1985) and
Schwarz & Arnason (1996), known as POPAN in program MARK. It is fitted by
data augmentation in the restricted dynamic occupancy (RDO) form of Royle &
Dorazio (2008) and Kéry & Schaub (2012, ch. 10):

- **Data augmentation.** The n animals caught are padded with M - n all-zero
  histories (potential animals never caught). Each of the M rows is a member
  of the superpopulation (the animals ever present during the study) with
  probability psi, so N_super ~ Binomial(M, psi).
- **Entry.** b_t is the share of the superpopulation that enters at occasion
  t (sum over t = 1). b_1 is the share already present at the first occasion.
  The hazard of entering at t among rows that have not yet entered is
  gamma_1 = psi b_1 and gamma_t = psi b_t / (1 - psi sum_{k<t} b_k).
- **Latent states.**
  - z[i,1] ~ Bern(gamma_1)
  - z[i,t] ~ Bern(phi_{t-1}^{dt} z[i,t-1] + gamma_t A[i,t-1])
  - A[i,t-1] = 1 if row i has not yet entered by t-1.
  - dt is the interval length in time units, so phi is always survival per
    time unit (e.g. per year), even when occasions are unevenly spaced.
- **Observation.** y[i,t] ~ Bern(p_t z[i,t]).
- **Survival.** logit(phi) = intercept (or one value per interval) + an
  effect of group (e.g. sex).
- **Capture.** logit(p_t) = intercept (or one value per occasion) + slopes on
  standardised occasion covariates (e.g. trap-nights).
- **Groups.** Group membership has prior Bern(pi). Animals with unknown
  group (and all augmented rows) are a mixture over the two groups.

**Priors.**
- psi ~ U(0,1) and b ~ Dirichlet(1,...,1). Together these give a uniform
  prior on N_super over 0..M and a uniform prior on the entry probabilities.
  We do not use independent uniform priors on gamma_t, because those pile
  prior mass on N_super near M (Royle & Dorazio 2008; Dorazio 2020).
- Intercepts ~ Logistic(0,1), which is uniform on the probability scale.
- Slopes and the group effect ~ Normal(0, sd = 2) on the logit scale.
- pi ~ U(0,1).

**Engines** (`settings$engine`, `settings$approach`).
- `nimble` + `marginal` (default). The alive state is summed out with a
  forward algorithm over three states (not yet entered, alive, dead). This
  runs once per unique capture history, and the M - n all-zero histories
  count as one weighted pattern (Saracco & Yackulic 2023). Cost does not
  grow with M. The latent counts (N_t, entrants, N_super, sex ratio) are then
  drawn exactly from their conditional posterior for each retained draw. Each
  animal's history is one of T(T+1)/2 + 1 "paths" (entry occasion, last
  occasion alive). The forward algorithm and the path sum agree to 3e-17
  (dev/jolly-seber/check_paths.R).
- `nimble` + `latent`. The textbook RDO code, with every z[i,t] sampled.
  It gives the same posterior (section 4b of VALIDATION.md) but is much
  slower (about 39 times at the default design) and is used as an in-module
  cross-check.
- `marked`. Maximum-likelihood POPAN with `marked::crm(model = "JS")`
  (Laake et al. 2013).
  - Intervals for derived quantities come from 2000 draws of the
    coefficients from their asymptotic normal distribution.
  - Its N_t is the *expected* abundance, not the realised one.
  - It cannot handle animals of unknown group. If any are present, the group
    effect is dropped in the marked fit.

## Outputs (estimand names)

| name | meaning |
|---|---|
| `superpopulation_size` | animals that used the population at some time during the study |
| `abundance_t` | realised number alive and present at occasion t (N_t) |
| `entrants_t` | animals that entered between occasion t-1 and t (B_t, t >= 2: births + immigrants that survived to t) |
| `entry_prob_t` | b_t, share of the superpopulation entering at t; `entry_prob_1` = share present at the start |
| `recruitment_rate_t` | per-capita recruitment B_{t+1} / N_t for interval t |
| `growth_rate_t` | realised growth rate N_{t+1} / N_t |
| `growth_rate_mean` | geometric mean growth rate per time unit, (N_T / N_1)^(1 / total time) |
| `survival` | apparent survival per time unit (reference group); `survival_1..T-1` when `survival: time` |
| `survival_<level>` | apparent survival of the second group level (e.g. `survival_male`) |
| `survival_beta_<level>` | logit-scale effect of group on survival (log odds ratio); a power parameter |
| `prop_<level>` | realised share of that group in the superpopulation |
| `detection` | capture probability at average covariate values (or `detection_t` with `detection: time`) |
| `detection_beta_<cov>` | logit-scale capture effect per 1 SD of the occasion covariate; a power parameter |

The report gives each estimate with a 95% credible interval, R-hat and bulk
ESS. It also includes posterior predictive checks (Freeman-Tukey discrepancies
of captures per animal and per occasion), R2ucare tests (3.SR transients,
2.CT trap dependence) and a check that the data-augmentation limit M does not
truncate N_super. If the 99.5% posterior quantile of N_super is >= 0.9 M, the
fit is redone with M doubled (up to twice). If it still presses on M, the fit
is marked not converged.

Figures (`docs/FIGURES.md`):
- `abundance_trajectory`: N_t with 95% and 50% bands, plus the number caught
- `recruitment`: entrants and per-capita recruitment
- `entry_probabilities`
- `growth_rate`
- `survival_<group>` or `survival_time`
- `detection_time`
- `detection_<cov>`: capture probability against the covariate in original
  units
- `coefficients`: odds ratios
- `gof`: posterior predictive checks
- `superpopulation`: posterior of N_super, with n caught marked
- `trace`, `posterior`
- `engine_comparison`: when `compare_engines` is set

## Assumptions (and what breaks)

1. **Equal capture probability** for all animals present at an occasion,
   marked or unmarked. Individual heterogeneity (some animals trap-shy or
   always in traps) biases abundance low. Trap-happiness after first capture
   biases it too. Check the PPC and R2ucare 2.CT.
2. **Emigration is permanent.** Survival is *apparent*: death and permanent
   emigration cannot be separated. **Temporary emigration** with animals that
   tend to stay away breaks the model. In validation it inflated
   N_super by +22% (coverage 0.43) and biased capture probability by -33%.
   This needs a robust design (Kendall et al. 1997), which this module does
   not provide.
3. **Marks are permanent and read correctly.** Tag loss looks like death, and
   misidentification creates false "new" animals.
4. **Closure within an occasion.** Each occasion is short relative to the
   intervals between occasions.
5. **Independent fates** among animals.
6. **Transients** (animals passing through once) inflate entry and lower
   survival. Check R2ucare 3.SR.
7. **Entry probabilities are free** for each occasion. Recruitment has no
   covariate model in this version.

## Data format

**captures** (one row per capture; repeat captures within an occasion count once):

```
turtle_id,pit_tag,capture_date,season,trap,sex,carapace_mm,mass_g
AM-089,982.149470056,2015-05-21,2015,HT02,M,164,585
AM-150,982.176990615,2015-05-23,2015,HT04,F,145,759
```
Only the ID, occasion and (optional) group columns are used. Other columns
are ignored. Unknown group values (`U`, `unk`, `?`, blank, `NA`) become
unknown and are modelled as a mixture.

**histories** (alternative): one row per animal (or per unique history with
a `freq` column), e.g. `ch = 01101000`.

**occasions** (recommended): one row per occasion, *including occasions with
no captures*:

```
season,start_date,end_date,trap_nights,observers
2015,2015-05-20,2015-08-03,93,2
```

## Config

```yaml
data:
  species: western pond turtles
  captures: turtle_captures.csv   # or histories: file.csv (+ history_col, freq_col)
  id_col: turtle_id               # default animal_id
  occasion_col: season            # default occasion
  group_col: sex                  # optional; exactly 2 levels
  group_recode: {F: female, M: male}
  group_levels: [female, male]    # reference first
  group_name: sex
  occasions: turtle_seasons.csv
  occasion_key: season            # column matching occasion_col
  interval_col: null              # interval to next occasion (default 1), or
  time_col: null                  # numeric time or dates -> intervals (time_unit)
  label_col: null
  detection_covariates: [trap_nights]
  time_unit: year                 # year | month | week | day (for date time_col)
settings:
  engine: nimble                  # or marked
  approach: marginal              # or latent
  survival: constant              # or time (one per interval)
  detection: auto                 # auto = constant + detection_covariates; or time; or [constant, cov...]
  group_effect: true              # survival differs by group
  M: null                         # default max(400, 5 n); grown automatically if limiting
  prior_sd_slope: 2
  niter: 20000                    # nburnin 5000, thin 5, nchains 3
  compare_engines: [marked]       # also fit and show other engines
figures: {...}                    # see docs/FIGURES.md
```

## Simulation and power

`params`:
- `superpopulation_size` (350)
- `entry_prob_first` (0.35, rest spread evenly), or `entry_probs` (one per
  occasion)
- `survival` (0.75 per time unit), or `survival_by_interval`
- `survival_beta_male` (-0.5), `prop_male` (0.5)
- `detection` (0.45), or `detection_by_occasion`
- `detection_beta_effort` (0.5 per SD)

`design`:
- `n_occasions` (8), `interval` (1)
- `groups` (TRUE), `group_levels`
- `prop_group_unknown` (0)
- `detection_covariate` (TRUE)
- `temp_emigration` (c(leave, return), default 0: the negative control)

Power for the sex effect on survival (log odds ratio -0.5, about 350
animals, p = 0.45) rises with the number of occasions. The gate passes at
every T from 5 to 10:

| n_occasions | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|
| power (survival_beta_male) | 0.18 | 0.22 | 0.50 | 0.64 | 0.74 | 0.70 |
| power (detection_beta_effort, 0.5/SD) | 0.74 | 0.90 | 0.92 | 0.98 | 1.00 | 1.00 |

These come from `power` with nsim = 50 per value; see VALIDATION.md. A
moderate sex difference in survival needs about 10 occasions or more animals.

## Pitfalls

- **Few occasions or low capture probability.** With p = 0.2 and 5
  occasions, posterior means of N_t were biased +11% to +24%, but 95%
  intervals still covered the truth (see VALIDATION.md). Report intervals,
  not just means.
- **Small entry probabilities and per-interval recruitment.** When few
  animals enter in an interval (about 20 or fewer, or a short interval
  between occasions), the posterior means of that `entry_prob_t` and
  `recruitment_rate_t` can be 10-20% too high (the ML fit puts them on 0
  instead), although their 95% intervals still cover. Interpret them
  through their intervals, or through `entrants_t` and the abundance
  trajectory (see VALIDATION.md, "Independent verification").
- **Fully time-varying p.** With `detection: time`, N_1, entry_prob_1 and the
  last survival/capture probability are only weakly identified (the classic
  JS confounding; Schwarz & Arnason 1996). The report warns about this.
- **First entrants are not recruits.** `entry_prob_1` is the share present
  at the start. Per-capita recruitment starts at interval 1-2.
- **Abundance is not density.** It counts the animals using the trapped
  area. Use `scr` / `scr-open` for density.
- **M (augmentation)** is checked automatically. An N_super posterior that
  presses against M is an error, not an estimate.

## Relation to other modules

- **cjs** conditions on first capture and models survival only. It is more
  robust to capture heterogeneity at first capture, but gives no abundance or
  recruitment. JS survival estimates should agree with CJS on the same data.
- **scr** (closed spatial capture-recapture) estimates density in one closed
  season from *where* animals were caught.
- **scr-open** adds the same open-population dynamics as this module to SCR:
  entry, survival, and a forward algorithm over not-entered / alive / dead.
  Its implementation follows the Taal SCR benchmark open models. It gives
  density and population dynamics but needs trap locations and individual
  spatial histories. Use jolly-seber when traps are not spatially resolved,
  or when the question is about the population using a site rather than a
  density.
- **ipm** can combine JS-type capture-recapture with counts and productivity
  data.

## References

- Jolly GM (1965) Explicit estimates from capture-recapture data with both death and immigration-stochastic model. Biometrika 52:225-247.
- Seber GAF (1965) A note on the multiple-recapture census. Biometrika 52:249-259.
- Crosbie SF, Manly BFJ (1985) Parsimonious modelling of capture-mark-recapture studies. Biometrics 41:385-398.
- Schwarz CJ, Arnason AN (1996) A general methodology for the analysis of capture-recapture experiments in open populations. Biometrics 52:860-873.
- Pradel R (1996) Utilization of capture-mark-recapture for the study of recruitment and population growth rate. Biometrics 52:703-709.
- Kendall WL, Nichols JD, Hines JE (1997) Estimating temporary emigration using capture-recapture data with Pollock's robust design. Ecology 78:563-578.
- Royle JA, Dorazio RM (2008) Hierarchical Modeling and Inference in Ecology. Academic Press (Jolly-Seber as restricted dynamic occupancy via data augmentation).
- Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(S2):521-537.
- Kéry M, Schaub M (2012) Bayesian Population Analysis using WinBUGS. Academic Press, ch. 10.
- Dorazio RM (2020) Objective prior distributions for Jolly-Seber models of zero-augmented data. Biometrics 76:1285-1296.
- Saracco JF, Yackulic CB (2023) Fast and flexible Bayesian Jolly Seber models and application to populations with transients. bioRxiv doi:10.1101/2023.05.25.542303.
- Yackulic CB, Dodrill M, Dzul M, Sanderlin JS, Reid JA (2020) A need for speed in Bayesian population models: a practical guide to marginalizing and recovering discrete latent states. Ecological Applications 30:e02112.
- Laake JL, Johnson DS, Conn PB (2013) marked: an R package for maximum likelihood and MCMC analysis of capture-recapture data. Methods in Ecology and Evolution 4:885-890.
- Gimenez O, Lebreton J-D, Choquet R, Pradel R (2018) R2ucare: an R package to perform goodness-of-fit tests for capture-recapture models. Methods in Ecology and Evolution 9:1749-1754.
