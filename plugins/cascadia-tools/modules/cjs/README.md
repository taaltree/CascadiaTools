---
id: cjs
title: Cormack-Jolly-Seber apparent survival
when_to_use: Individually marked animals (bands, PIT tags, photo-ID, genetic tags) recaptured or resighted on repeated occasions; the question is how many survive from one occasion to the next and whether survival differs by sex/age class or varies with individual traits, environmental conditions or time.
not_for: Abundance or recruitment (use jolly-seber / scr-open), unmarked counts (nmixture), dead recoveries (band-recovery models), movement between sites or states (multistate), known-fate telemetry (use a known-fate / survival-time model).
data_needed: capture history per individual (wide 0/1 history, long table of id x occasion, or long table of id x capture date), optional individual table (group, age at marking, covariates) and occasion table (environmental covariates)
engine: marked (maximum likelihood, default) or nimble + nimbleEcology (Bayesian; marginal dCJS likelihood or latent alive-state formulation); R2ucare for goodness of fit
---

# cjs: apparent survival from capture-recapture

## Question it answers

"What proportion of my marked animals survive from one year (season, month) to
the next, does that differ between males and females or juveniles and adults,
and is it associated with body condition or winter weather?" The model
separates survival from imperfect recapture: an animal not seen this year may
be alive but missed. It conditions on first capture, so it estimates
survival and recapture probability, not population size.

**Apparent survival = true survival x site fidelity.** An animal that
emigrates permanently is indistinguishable from one that died. Say so in
every interpretation.

## Model

For individual *i* first caught at occasion *f_i*, alive state *z_it*:

- z_{i,f_i} = 1;  z_{i,t+1} | z_{it} ~ Bernoulli(z_{it} phi_{it}^{Delta_t})
- y_{it} | z_{it} ~ Bernoulli(z_{it} p_{it}),  t > f_i
- logit(phi_{it}) = X_phi beta_phi [+ eps_t, eps_t ~ N(0, sigma^2) with `(1 | time)`]
- logit(p_{it}) = X_p beta_p

Delta_t is the interval length (survival is reported per unit time). The
likelihood engines integrate z out (Lebreton et al. 1992; marked: Laake et
al. 2013; Nimble: nimbleEcology `dCJS_vv`, the forward algorithm of Turek
et al. 2016). The latent formulation (`formulation: latent`) samples z as in
Gimenez et al. (2007) and Kery & Schaub (2012, ch. 7).

Variables the formulas can use:

| variable | meaning |
|---|---|
| `time` | one parameter per interval (phi) / occasion (p) |
| `Time` | linear trend over occasions (0, 1, 2, ...) |
| group column (e.g. `sex`) | grouping factor; first level (or `group_levels[1]`) is the reference |
| `tsm` | time since marking: first interval after marking vs later intervals (transients; Pradel et al. 1997). Survival estimates then refer to residents, and `proportion_transient` = 1 - phi_first/phi_later. `tsm` is automatically crossed with the other survival terms (`~ tsm + sex + x` is fitted as `~ tsm * (sex + x)`, MARK's phi(a2 x ...)), because transients scale first-interval survival multiplicatively and an additive logit model attenuates covariate effects by 16-25% (VALIDATION.md). `phi_beta_*` are then resident effects; `phi_beta_tsm_first_*` are first-interval deviations. `settings: tsm_interact: false` keeps the additive model |
| `age_class` | juvenile in the first interval after marking as a juvenile, adult afterwards (needs `age_col`; phi only) |
| `td` | trap dependence: 1 if caught at the previous occasion (p only; Pradel 1993) |
| individual covariates | numeric covariates are standardised (per 1 SD); 0/1 covariates are used as indicators |
| occasion covariates | value in row t applies to the interval t -> t+1 (phi) and to occasion t (p) |
| `(1 | time)` | random interval effect on logit survival (engine: nimble only) |

Priors (Nimble): Logistic(0, 1) on logit-scale intercepts (= Uniform(0, 1) on
the probability scale), Normal(0, SD 2) on all other logit-scale coefficients
(`prior_sd`), half-Normal(0, 1) on the SD of random time effects.

## Assumptions (and what breaks)

- **Homogeneous survival and recapture within each modelled group/covariate
  combination.** Unmodelled individual heterogeneity in p biases p and can bias
  survival; check with R2ucare (TEST2/TEST3).
- **No transients, or tsm in the model.** If a fraction of newly marked animals
  never return (passing migrants), survival is biased LOW (validated: -12% to
  -14% with 30% transients, coverage 0.45-0.50, gate FAILS). Test 3.SR flags
  this; add `tsm`.
- **No trap dependence, or td in the model.** Test 2.CT flags trap-happiness
  (negative signed z) or trap-shyness.
- **Marks are not lost and are read correctly**; tag loss looks like death.
- **Instantaneous sampling relative to the interval**, independent fates.
- **Emigration is permanent or random**; temporary emigration is not modelled.
- Animals removed at capture (died in the net, kept) must be excluded from
  later occasions; this module does not handle losses on capture: drop their
  later occasions by censoring them (remove the animal or truncate its history).

## Data format

Captures (choose one `format`):

- `long`: one row per capture, `id_col` + `occasion_col` (year or occasion label; numeric labels such as 1981, 1982 are treated as calendar years, so survival is annual and gaps between years become interval lengths)
  ```
  bird_id,year,sex
  D001,1981,Female
  D001,1983,Female
  ```
- `dates`: one row per capture, `id_col` + `date_col`; dates are binned into
  occasions by `occasion_unit` (year/month/week) or `occasion_breaks` (start
  dates). Repeat captures within an occasion are collapsed.
  ```
  band_number,capture_date,net
  2841-37593,2016-05-11,N09
  ```
- `wide`: one row per individual with `ch_col` (history string "01011") or
  `occasion_cols` (one 0/1 column per occasion).

Optional `individuals` CSV (one row per individual: id, group, age at marking,
covariates; covariates must be constant per individual, e.g. mass at first
capture) and `occasion_covs` CSV (one row per occasion: label column +
covariates; the last occasion's value is not used for survival).

## Config

```yaml
data:
  species: adult song sparrows
  format: dates            # wide | long | dates
  captures: captures.csv
  id_col: band_number
  date_col: capture_date   # dates format
  date_format: "%Y-%m-%d"
  occasion_unit: year      # year | month | week, or occasion_breaks: [2019-05-01, ...]
  season_months: [5, 6, 7, 8]   # optional: drop captures outside the sampling season
  occasion_col: year       # long format
  ch_col: ch               # wide format (or occasion_cols: [...])
  occasions: [...]         # optional explicit occasion list (in order)
  time_intervals: [...]    # optional; default 1, or year differences for year labels
  individuals: birds.csv
  group_col: sex
  group_levels: [F, M]     # reference level first
  age_col: age_at_banding  # optional: adult / juvenile
  individual_covs: [mass_g]
  occasion_covs: years.csv
  occasion_cov_col: year
  covariate_units: {mass_g: g}
formulas:
  phi: "~ tsm + sex + winter_min_temp_c + mass_g"    # add "+ (1 | time)" for random years (nimble)
  p: "~ sex"
settings:
  engine: marked           # marked | nimble
  formulation: marginal    # nimble: marginal (dCJS) | latent (z sampled; slower, cross-check)
  niter: 6000              # nimble: iterations per chain (incl. burn-in)
  nburnin: 2000
  nchains: 3
  prior_sd: 2              # SD of Normal priors on logit coefficients
  gof: true                # R2ucare tests and c-hat
  chat_adjust: false       # marked: inflate SEs by sqrt(c-hat of the fitted model) when > 1
  tsm_interact: true       # cross tsm with the other survival terms (resident effects)
figures:
  labels: {mass_g: "Body mass (g)", F: Female, M: Male}   # also used in the report text
```

## Outputs

- `survival[_<group>][_<age>][_<interval>]`: apparent survival per unit time
  at average continuous covariates (reference level of 0/1 indicators;
  residents when `tsm` is in the model; average interval with `(1 | time)`).
- `detection[_<group>][_<occasion>][_caught_previous]`: recapture probability.
- `phi_beta_<term>`, `p_beta_<term>`: logit-scale coefficients (per 1 SD for
  continuous covariates; odds ratio = exp).
- `proportion_transient[_<group>]` (with `tsm`; 1 - phi_first/phi_later, unconstrained, so its table interval can dip below 0 = no evidence of transients; the report text truncates to 0-1), `survival_temporal_sd` (with
  `(1 | time)`), `survival_x_detection_<last>` (fully time-dependent models:
  only the product of the last phi and p is identifiable).
- Diagnostics: R2ucare TEST 3.SR (transience), 3.SM, 2.CT (trap dependence),
  2.CL, overall chi-square and c-hat of the general time-dependent model, plus
  c-hat of the FITTED model (3.SR removed when `tsm` is fitted, 2.CT removed
  when `td` is fitted; Pradel et al. 2005). `chat_adjust: true` uses the latter.
  R-hat/ESS for Nimble. Likelihood fits with a probability estimate < 0.001 or
  > 0.999, or a logit coefficient with |estimate| or SE > 10, are flagged NOT
  converged (boundary / not identifiable).
- Coefficients of the `time` factor are not reported as rows (use
  `survival_<interval>` / `detection_<occasion>`); `phi_beta_Time` (trend) is.
- Figures (`figures/`): `survival_<cov>` and `detection_<cov>` conditional
  effects (response scale, original units), `survival_effects_panel`,
  `survival_time`/`detection_time` (per interval/occasion), `survival_groups`,
  `detection_groups`/`detection`, `coefficients` (odds ratios), `gof`,
  `captures` (newly marked vs recaptured per occasion), `trace`/`posterior`
  (Nimble).

## Simulation and power

`design`: `n_occasions` (8), `n_marked_per_occasion` (50, released at
occasions 1..T-1), `prop_male`, `prop_juvenile`, covariate means/SDs,
`phi_formula`, `p_formula` (the model that is FITTED). `params`: `survival`
(reference group, average covariates), `detection_prob`, `phi_beta_sex_male`,
`phi_beta_winter_severity`, `phi_beta_body_mass`,
`phi_beta_age_class_juvenile`, `p_beta_sex_male`, `p_beta_td`,
`prop_transient`, `temporal_sd` (random year SD on logit survival),
`p_temporal_sd`. The simulator always applies every effect in `params`
(sex, winter severity, body mass, juvenile, td), whatever the fitted
formulas say. Set effects that the fitted formula omits to 0, otherwise the
"truth" (reference group, average covariates) is not what the data average
over (e.g. an omitted winter effect inflates the random-year SD). Useful power knobs: `n_marked_per_occasion`,
`n_occasions`, `detection_prob`, effect sizes. See VALIDATION.md for what
recovers and where it breaks (few occasions with low p: Bayesian posterior means biased, effects undetectable; random-year SD biased up with < 15 intervals).

## Pitfalls

- An occasion with no captures at all (not sampled) must be removed and the
  gap expressed through `time_intervals`, otherwise p for that occasion is 0.
- Fully time-dependent phi(t) p(t): last phi and p are confounded. The module
  fixes recapture at the last occasion to 1, so the last survival parameter
  estimates the product, reported as `survival_x_detection_<interval>`. This is
  exact when both last parameters are free; if the formulas add other terms
  shared across occasions (`~ time + sex`) the constraint also acts through them.
- With few occasions (< 5) and low p (< 0.3), survival intervals under-cover;
  the variance of random year effects is poorly estimated with < 10 intervals.
- c-hat > 3 means a structural problem, not overdispersion: fix the model.
- "Survival" here is apparent survival; dispersal-prone age classes (juveniles)
  look like they die.

## References

- Cormack RM (1964) Biometrika 51:429-438. Jolly GM (1965) Biometrika 52:225-247. Seber GAF (1965) Biometrika 52:249-259.
- Lebreton J-D, Burnham KP, Clobert J, Anderson DR (1992) Modeling survival and testing biological hypotheses using marked animals. Ecological Monographs 62:67-118. doi:10.2307/2937171
- Pradel R (1993) Flexibility in survival analysis from recapture data: handling trap-dependence. In: Marked Individuals in the Study of Bird Population, Birkhauser, 29-37.
- Pradel R, Hines JE, Lebreton J-D, Nichols JD (1997) Capture-recapture survival models taking account of transients. Biometrics 53:60-72. doi:10.2307/2533097
- Pradel R, Gimenez O, Lebreton J-D (2005) Principles and interest of GOF tests for multistate capture-recapture models. Animal Biodiversity and Conservation 28:189-204.
- Gimenez O, Lebreton J-D, Choquet R, Pradel R (2018) R2ucare: An R package to perform goodness-of-fit tests for capture-recapture models. Methods in Ecology and Evolution 9:1749-1754. doi:10.1111/2041-210X.13014
- Choquet R, Lebreton J-D, Gimenez O, Reboulet A-M, Pradel R (2009) U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data. Ecography 32:1071-1074. doi:10.1111/j.1600-0587.2009.05968.x
- Meredith M (2024) wiqid: Quick and Dirty Estimates for Wildlife Populations, R package (survCJS; used only for the cross-check in VALIDATION.md).
- Laake JL, Johnson DS, Conn PB (2013) marked: an R package for maximum likelihood and MCMC analysis of capture-recapture data. Methods in Ecology and Evolution 4:885-890. doi:10.1111/2041-210X.12065
- Gimenez O et al. (2007) State-space modelling of data on marked individuals. Ecological Modelling 206:431-438. doi:10.1016/j.ecolmodel.2007.03.040
- Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS, ch. 7. Academic Press.
- Schaub M, Kery M (2022) Integrated Population Models, ch. 3. Academic Press.
- Turek D, de Valpine P, Paciorek CJ (2016) Efficient MCMC sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics 23:549-564. doi:10.1007/s10651-016-0353-z
- Goldstein BR, Turek D, Ponisio L, de Valpine P (2024) nimbleEcology R package v0.5.0.
- de Valpine P et al. (2017) NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487
- Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference, 2nd ed. Springer.
