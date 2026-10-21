---
id: multistate
title: Multistate and multievent capture-recapture (Arnason-Schwarz / Pradel multievent)
when_to_use: Individually marked animals are caught or resighted on repeated occasions and each capture records a state (site or wetland, breeding vs non-breeding, infected vs uninfected, size class), possibly with some captures where the state could not be determined; the question is survival in each state, how often animals move or switch between states, and what is associated with those transitions.
not_for: Plain survival without states (use cjs); abundance or recruitment (jolly-seber, scr-open); states that can be misclassified as another state (only "state recorded correctly or coded unknown" is modelled); unobservable states (e.g. temporary emigration to an unsampled site, skipped breeders that cannot be seen); dead recoveries; individual random effects; memory (second-order) transitions; animals with continuous position (use scr / movement-hmm).
data_needed: capture records with animal id, occasion and recorded state (long CSV, wide table or capture-history strings; 0 = not seen, a state name, or an "unknown state" code), optional individual covariates (id + columns) and occasion covariates (occasion + columns)
engine: nimble + nimbleEcology (dDHMMo, forward algorithm); option engine 'ml' = hand-written forward-algorithm likelihood maximised with optim (also used as a cross-check); marked (crm 'hmmMSCJS') used as an independent cross-check
---

# Multistate and multievent capture-recapture

## Question it answers

"How likely are my animals to survive in each state, how often do they
move or switch between states, and is switching associated with a
covariate?" Typical examples:

- **Movement between sites**: survival at each of three wetlands and the
  annual probability of moving between them (Arnason 1973; Schwarz et al.
  1993; Kéry & Schaub 2012 ch. 9).
- **Disease dynamics**: survival of infected and uninfected frogs, annual
  infection and recovery probabilities, and whether water temperature is
  associated with them. Swabs that are lost or inconclusive are kept as
  "seen, state unknown" (multievent; Conn & Cooch 2009).
- **Breeding status**: survival cost of breeding and the probability of
  breeding again, when breeding status is sometimes not determined (Pradel
  2005; Lebreton et al. 2009).

## Model

Animals are followed from their first capture. States 1..S plus *dead*.
Between occasions t and t+1 (one *interval*, e.g. a year):

- survival in state r: `phi_r = logit^-1(a_phi[r] + sum_k bphi[k] x_k)`.
  A survival covariate can be per individual (body size) or per occasion
  (winter severity). Its slope is shared by all states.
- transition for animals that survive, r -> s:
  `psi_rs = exp(eta_rs) / sum_u exp(eta_ru)` with `eta_rr = 0` (staying is
  the reference) and `eta_rs = apsi[r,s] + sum_k bpsi[r,s,k] x_k`. That is
  a multinomial logit with one slope per transition and covariate.
- detection of a live animal in state s at an occasion:
  `p_s = logit^-1(a_p[s] + sum_k bp[k] x_k)` with occasion covariates
  (effort, weather), slope shared by states.

Observations at each occasion: not seen; seen in state s; or (multievent
only) seen, state unknown. In the **multievent** model (Pradel 2005) a
detected animal in state s has its state recorded with probability
`delta_s` (**assignment**) and is otherwise coded unknown. This is the
"partial observation" structure of Conn & Cooch (2009) and Kendall (2009):
a recorded state is assumed correct. The state at first capture is
uncertain for animals first seen as "unknown", so the model also estimates
the state distribution at first capture `pi_s` (**initial**). The first
event of every animal contributes `pi_s delta_s` (seen in s) or
`sum_s pi_s (1 - delta_s)` (unknown).

In the **Arnason-Schwarz** model (no unknown records) delta = 1 and the
likelihood is conditioned on the observed state at first capture.

The latent state sequence is summed out with the forward algorithm
(hidden Markov model; Pradel 2005; Zucchini et al. 2016) via
nimbleEcology's `dDHMMo` (Turek et al. 2016; Goldstein et al. 2024), so no
latent states are sampled.

Derived quantities, computed per posterior draw:

- `stationary_<s>`: the stationary distribution of the transition matrix
  (left eigenvector of psi). It is the long-run share of time a surviving
  animal spends in each state if transitions stayed as they are. It
  ignores survival and recruitment differences, so it is **not** the
  predicted population structure.
- expected residence time in state r, `1 / (1 - psi_rr)` intervals. This
  is reported in the interpretation.

**Priors** (Nimble): logistic(0,1) on every logit-scale intercept, which is
uniform(0,1) on survival, detection and assignment at average covariates.
logistic(0,1) on each multinomial-logit transition intercept and initial
state logit. Normal(0, SD = 2) on covariate slopes (per SD of the
covariate; set with `settings: prior_sd`). Covariates are standardised
internally. Figures and interpretation use original units.

## Estimands (names in the estimates table)

| name | meaning | scale |
|---|---|---|
| `survival_<s>` | survival over one interval in state s (average covariates) | probability |
| `transition_<r>_to_<s>` | P(in s at next occasion \| in r now and survives) | probability |
| `detection_<s>` | P(detected \| alive in s) (average covariates) | probability |
| `assignment_<s>` | P(state recorded \| detected in s) (multievent) | probability |
| `initial_<s>` | share of newly marked animals in state s (multievent) | probability |
| `stationary_<s>` | stationary distribution of the transition matrix | proportion |
| `survival_beta_<cov>` | logit slope of survival per 1 SD of cov | link |
| `transition_beta_<cov>_<r>_to_<s>` | multinomial-logit slope (log odds of moving r->s vs staying) per 1 SD | link |
| `detection_beta_<cov>` | logit slope of detection per 1 SD | link |

## Assumptions (and what breaks)

1. **Markov transitions.** Where an animal goes next depends only on
   where it is now. Site fidelity to a previous site ("memory") violates
   this. The U-CARE test WBWA and the PPC "state switches" statistic
   flag it.
2. **Homogeneity.** Animals in the same state share survival, transition
   and detection (apart from covariates). Transients and trap-happiness
   inflate the U-CARE 3G.SR and M.ITEC tests. Unmodelled heterogeneity
   biases survival low.
3. **Recorded states are correct.** Misclassification (infected recorded
   as uninfected) is not modelled. If it happens, transition
   probabilities are inflated.
4. **Assignment does not depend on time.** `delta_s` is constant per
   state, and the same `delta_s` applies at first capture and at
   recaptures (Pradel 2005 allows a separate first-capture event matrix;
   that is not offered here).
5. **All states are observable.** Detection > 0 in every state. An
   unobservable state needs a different model.
6. **Constant parameters over time** apart from occasion covariates.
   There are no time-varying random effects. If survival varies by year
   for reasons not in the covariates, intervals are too narrow.
7. Marks are not lost and are read correctly. Animals are independent.
   Death and emigration are confounded: survival is *apparent* survival.

## Data format

Long format (default), one row per capture (example/frog_captures.csv):

```
frog_id,year,bd_status
RL-0001,2015,U
RL-0001,2016,uninfected
RL-0001,2021,infected
RL-0002,2015,infected
```

- `id_col`: animal id.
- `occasion_col`: the occasion (year, session). Occasions are sorted
  numerically, or by the number inside a label such as "S1".."S12". Use
  `occasion_order` for anything else.
- `state_col`: a state name from `states`, the `unknown_code` (default
  `U`), or a number 1..S. Two different known states for the same animal
  and occasion are an error. A known state and a U on the same occasion
  keep the known state.
- `format: wide` means one row per animal with one column per occasion,
  cells holding 0, a state name or the unknown code. `format: ch` means
  one row per animal with a capture-history string (`ch_col`), one
  character per occasion, for example `0A0BU0`.
- `individuals`: CSV with the id column plus numeric covariates (for
  example `svl_mm`). These can be used on survival or transitions.
- `occasions`: CSV with the occasion column plus numeric covariates (for
  example `water_temp_C`, `survey_hours`). A value at occasion t is used
  for detection at t and for survival and transitions over the interval
  t -> t+1.

`check_data()` reports problems such as fewer than 3 occasions, too few
recaptures, a state that is never recorded, missing covariates, unknown
state codes and duplicated ids.

## Config

```yaml
data:
  captures: frog_captures.csv
  format: long            # long | wide | ch
  id_col: frog_id
  occasion_col: year
  state_col: bd_status
  states: [uninfected, infected]
  unknown_code: U
  model: auto             # auto (multievent if any unknown) | multievent | AS
  interval_unit: year     # wording only
  individuals: frog_individuals.csv
  occasions: survey_years.csv
  survival: [svl_mm]
  transition: [water_temp_C]
  detection: [survey_hours]
settings:
  engine: nimble          # nimble (default) | ml (hand-written likelihood + optim; fast, Wald-type intervals)
  niter: 6000             # per chain
  nburnin: 2000
  nchains: 3
  sampler: slice          # slice (default) | RW_block | AF_slice
  prior_sd: 2             # SD of Normal priors on slopes
  unknown_handling: model # model (multievent) | drop (treat unknown as not seen) | assign (force into naive_state: NEGATIVE CONTROL only)
  ppc_draws: 300          # posterior predictive check draws (0 = skip)
  ucare: true             # R2ucare goodness-of-fit tests
figures:                  # see docs/FIGURES.md
  labels: {uninfected: "Uninfected", water_temp_C: "Mean water temperature (°C)"}
```

## Outputs

- An estimates table of the quantities listed above with posterior
  mean, SD and 95% credible interval. Diagnostics give max R-hat and min
  bulk ESS, with a note when the real-data standard (R-hat < 1.01, ESS
  >= 400) is not met. The notes also report the PPC p-values and the
  U-CARE overall test.
- Interpretation:
  - survival per state, with the difference and P(difference > 0);
  - transition probabilities and residence times;
  - the stationary distribution, with its caveat;
  - detection and assignment;
  - covariate effects as changes in the probability from the 10th to the
    90th percentile of the covariate, odds ratios per SD, and
    P(effect > 0);
  - goodness of fit.
- Figures (`figures/`):

  | figure | content |
  |---|---|
  | `vital_rates_by_state` | survival, detection, assignment and initial state for each state |
  | `transitions` | heatmap of the transition matrix |
  | `transition_estimates` | forest plot of the transition probabilities |
  | `stationary` | the stationary distribution |
  | `survival_<cov>`, `transition_<cov>`, `detection_<cov>` | conditional effects on the probability scale, original units, 50% and 95% bands, rug of observed values |
  | `detection_by_occasion` | detection at each occasion |
  | `coefficients` | odds ratios |
  | `captures_by_occasion` | raw data |
  | `gof` | three posterior predictive checks |
  | `trace`, `posterior` | MCMC traces and posterior densities |

## Simulation and power

`simulate(params, design, seed)`:

- Parameters are named on the biological scale with the state names:
  `survival_A`, `transition_A_to_B`, `detection_B`, `assignment_A`,
  `initial_A`, and slopes `survival_beta_<cov>`,
  `transition_beta_<cov>_<r>_to_<s>`, `detection_beta_<cov>`.
- Design knobs:
  - `states` (any number ≥ 2);
  - `n_occasions`;
  - `n_marked` (newly marked animals at each occasion 1..T-1);
  - `model` (`multievent` or `AS`);
  - `survival_covs` (individual covariates);
  - `transition_covs` and `detection_covs` (occasion covariates).
- Default: 2 states, 8 occasions, 40 new animals per occasion (280
  animals), multievent with assignment 0.7/0.8, one covariate on each
  process.
- Sensible power knobs: `n_marked`, `n_occasions`, and the assignment
  probabilities (how often the state is determined).
- Occasion covariates have only T−1 values, so effects on transitions
  are the weakest-estimated quantities. In the power run (VALIDATION.md)
  power for the A->B temperature slope (0.8 per SD) rose from 0.40 to
  0.98 as `n_marked` went from 10 to 80 per occasion. Power for the
  B->A slope (-0.5 per SD) only reached 0.55 at 80 per occasion. More
  occasions (more covariate values) also help, and they are the only
  remedy for confounding with other year effects.

## Pitfalls

- **Do not force unknown states into a default state.** That is the
  negative control in VALIDATION.md, and it biases transitions strongly.
  Dropping unknown records (`unknown_handling: drop`) keeps estimates
  unbiased if assignment does not depend on time, but it wastes data and
  confounds detection with assignment.
- **Last-interval confounding.** Survival and detection in the final
  interval are only separable because parameters are constant over time.
  Fully time-dependent models are not offered for this reason.
- **The stationary distribution is not a population forecast.**
- **Transitions are conditional on survival.** An "infection probability"
  of 0.3 means 30% of uninfected frogs that survive the year are infected
  at the next survey.
- **Occasion covariates with few occasions** (8 values) give weak
  evidence and are confounded with any other year effect. Treat them as
  associations.
- Animals first captured on the last occasion contribute only to `initial`
  (multievent) and nothing else.

## References

- Arnason AN (1973) The estimation of population size, migration rates and survival in a stratified population. Researches on Population Ecology 15:1-8.
- Schwarz CJ, Schweigert JF, Arnason AN (1993) Estimating migration rates using tag-recovery data. Biometrics 49:177-193.
- Brownie C, Hines JE, Nichols JD, Pollock KH, Hestbeck JB (1993) Capture-recapture studies for multiple strata including non-Markovian transitions. Biometrics 49:1173-1187.
- Lebreton J-D, Nichols JD, Barker RJ, Pradel R, Spendelow JA (2009) Modeling individual animal histories with multistate capture-recapture models. Advances in Ecological Research 41:87-173.
- Pradel R (2005) Multievent: an extension of multistate capture-recapture models to uncertain states. Biometrics 61:442-447. doi:10.1111/j.1541-0420.2005.00318.x
- Conn PB, Cooch EG (2009) Multistate capture-recapture analysis under imperfect state observation: an application to disease models. Journal of Applied Ecology 46:486-492.
- Kendall WL (2009) One size does not fit all: adapting mark-recapture and occupancy models for state uncertainty. In Thomson DL et al. (eds) Modeling Demographic Processes in Marked Populations, Springer, pp 765-780.
- Kéry M, Schaub M (2012) Bayesian Population Analysis using WinBUGS. Academic Press. Chapter 9.
- Turek D, de Valpine P, Paciorek CJ (2016) Efficient Markov chain Monte Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics 23:549-564.
- Goldstein B, Turek D, Ponisio L, de Valpine P (2024) nimbleEcology: Distributions for Ecological Models in nimble. R package.
- de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.
- Zucchini W, MacDonald IL, Langrock R (2016) Hidden Markov Models for Time Series, 2nd ed. CRC Press.
- Pradel R, Wintrebert CMA, Gimenez O (2003) A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model. Biometrics 59:43-53.
- Gimenez O, Lebreton J-D, Choquet R, Pradel R (2018) R2ucare: An R package to perform goodness-of-fit tests for capture-recapture models. Methods in Ecology and Evolution 9:1749-1754.
- Laake JL, Johnson DS, Conn PB (2013) marked: an R package for maximum likelihood and Markov Chain Monte Carlo analysis of capture-recapture data. Methods in Ecology and Evolution 4:885-890.
