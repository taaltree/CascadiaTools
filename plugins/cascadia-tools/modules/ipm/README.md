---
id: ipm
title: Integrated population model (counts + capture-recapture + productivity)
when_to_use: A local population studied for several years (>= 8, ideally 12+) with three kinds of data on the same population (annual counts of breeders, capture-recapture or resighting of individually marked young and adults, and brood/litter productivity), and the question is how population size changed, which vital rates drove the changes, how much immigration contributes, or how precise survival estimates can be made by combining the data.
not_for: Only counts (use timeseries), only capture-recapture (cjs / jolly-seber / multistate), only productivity (glm), projecting extinction risk (pva), populations with more than two age classes that differ strongly in vital rates (needs a custom stage structure), species where males and females must both be modelled, density dependence (not in this module).
data_needed: three CSV files (annual counts; capture histories in long, wide or m-array form with age at marking; productivity per brood or per year with optional mother age) plus an optional year-covariate table
engine: nimble (MCMC); cross-checked against the same model in JAGS (book code) and marked (CJS maximum likelihood)
---

# ipm: integrated population model

## Question it answers

"Our breeding population went from 106 to 135 females in 15 years with big
swings in between. Was that driven by first-year survival, adult survival,
breeding success or immigrants? How many breeders are immigrants? Would the
population persist without them?" An integrated population model (IPM)
answers this by fitting one population model jointly to all the data a
long-term study usually has:

- **counts** of breeders each year (with observation error),
- **capture-recapture** of birds ringed as nestlings and as adults (survival),
- **productivity**: fledglings per monitored brood (fecundity).

Because the three data sets must agree through one population model, the IPM
(i) estimates **annual population size** with less noise than the counts,
(ii) estimates **immigration**, for which there are no direct data, as the
part of the change in counts that local survival and reproduction cannot
explain (Abadi et al. 2010), (iii) usually gives more precise vital rates
than separate analyses, especially the rate with the weakest direct data,
typically juvenile (first-year) survival (Besbeas et al. 2002; Abadi et al.
2010b; Schaub & Kéry 2022, ch. 5), and (iv) decomposes the variation of the
population growth rate among the vital rates with a **life-table response
experiment (LTRE)** (Koons et al. 2016, 2017; Schaub & Kéry 2022, 9.4).

## Model

Female-based, pre-breeding census, two age classes at breeding plus
immigrants (Schaub & Kéry 2022, hoopoe model in 9.2 and woodchat shrike in
ch. 11; parameterisation of Abadi et al. 2010):

```
N1[t+1] ~ Poisson( rho * phij[t] * ( f1[t] * (N1[t] + N3[t]) + f2[t] * N2[t] ) )   local recruits (1-year-olds)
N2[t+1] ~ Binomial( Ntot[t], phia[t] )                                            surviving adults
N3[t+1] ~ Poisson( omega[t] * Ntot[t] )                                           immigrants
Ntot[t] = N1[t] + N2[t] + N3[t]
```

- `rho` = proportion of females among fledglings (`sex_ratio`, default 0.5);
  fecundity `f` = fledglings of both sexes per breeding female.
- Immigrants are assumed to breed like yearlings (first-time breeders).
- `omega[t]` = immigration rate: new immigrant females per female present
  the year before (Schaub & Fletcher 2015 compare parameterisations; the
  rate form is the book's default).
- Annual rates: `logit(phij[t]) = logit(mean_phij) + beta * x[t] + eps[t]`,
  `eps ~ Normal(0, sd)`; likewise for `phia` (logit) and `f1`, `f2`,
  `omega` (log). Year covariates optional on every rate.

Data likelihoods (assumed independent):

- Counts: `C[t] ~ Poisson(Ntot[t])` (default) or `Normal(Ntot[t], sigma)`.
- Productivity: `J1[t] ~ Poisson(B1[t] f1[t])`, `J2[t] ~ Poisson(B2[t] f2[t])`
  (J = fledglings, B = monitored broods of yearling / older mothers).
- Capture-recapture: age-dependent Cormack-Jolly-Seber in m-array form
  (multinomial; Burnham et al. 1987): animals released as juveniles survive
  the first interval with `phij[t]` and then with `phia`; animals released as
  adults with `phia`. One constant recapture probability `p` for both age
  classes (the book's hoopoe model allows a separate `pj` for birds marked
  as young; not implemented here).

Priors (as in the book): mean survival and `p` U(0,1); mean fecundity
U(0, `fmax` = 20); mean immigration rate U(0,3); temporal SDs U(0,3);
coefficients Normal(0, SD 2) on standardised covariates; initial stage sizes
discrete uniform on 1..2 x first count.

**Separate analysis** (`components: separate`): the same CJS and
productivity models without counts and population model, i.e. what you would
get analysing the data sets one by one. With `compare_separate: true` (real-
data default) both are fitted and the gain in precision is reported.

**LTRE** (Koons et al. 2016, 2017; Schaub & Kéry 2022, 9.4). Growth rate as
a function of rates and stage proportions n1, n2, n3:
`lambda = rho*phij*(f1*(n1+n3) + f2*n2) + phia + omega`. Fixed-design LTRE:
contribution of rate k to Var(lambda_t) = sum_l cov(theta_k, theta_l) s_k s_l,
with sensitivities s at the mean rates; computed per posterior draw, so
contributions come with credible intervals. Population-structure
contributions (n1, n2, n3) are summed. A transient (year-to-year) LTRE is drawn
as a figure.

## Assumptions (and what breaks if they fail)

- **One closed demographic unit apart from immigration.** Permanent
  emigration of adults appears as lower adult survival; emigration of
  young as lower juvenile survival (apparent survival). That is consistent
  across data sets, so N is fine, but say "apparent".
- **Independence of the data sets.** The same birds are counted, ringed and
  their broods monitored. Simulations show this has small effects on bias and
  coverage (Abadi et al. 2010, Ecology 91: 7-14; Schaub & Kéry 2022, 5.5).
- **Counts are unbiased for breeding females** (Poisson or normal error
  around the true number). A count that systematically misses 20% of
  females makes immigration and survival wrong, because the IPM reconciles
  everything through N. If counts are an index (not complete), say so and
  interpret N as an index.
- **Immigrants breed like yearlings.**
- **Population model right.** Any process missing from the model (e.g.
  immigration when it exists, skipped breeding, marker loss) is absorbed by
  the other rates and biases them (Riecke et al. 2019; our negative control:
  dropping immigration biases juvenile survival upwards by +19%, coverage 0.53,
  and adult survival upwards; fecundity is held in place by the productivity data).
- **Homogeneous CJS.** No trap-dependence or transience. Check with R2ucare
  (cjs module) if in doubt.
- No density dependence. Year effects are independent normal deviates.

## Data format

Three files (four with covariates). Years in any file define the study
period (contiguous); a year with no count is fine (leave it empty).

**counts** (one row per year):

| year | females_counted |
|---|---|
| 2006 | 106 |
| 2013 |  (empty: no survey) |

**captures**, `capture_format: long` (default): one row per capture or
resighting of an individual per year. `age` = age at THAT capture:
`nestling`/`juvenile`/`chick`/`pullus`/`HY`/`0`/`1` = young of the year;
`adult`/`yearling`/`SY`/`ASY`/`2` = breeder. Only the age at first capture
matters. A bird first caught as a breeding yearling counts as an adult (it
already survived its first year). Extra columns (sex) are ignored. Include
only females, plus young of unknown sex, for a female-based model.

| ring | year | age | sex |
|---|---|---|---|
| NX12345 | 2006 | nestling | U |
| NX12345 | 2008 | adult | F |

`capture_format: wide`: one row per animal, `id`, `age` (at marking), and
one 0/1 column per year named like `2006` or `y2006`.
`capture_format: marray`: an age-dependent m-array, columns `age`
(juvenile/adult), `release_year`, one column per recapture year, `never`.
This is what most published IPMs report (as written by `write_data`).

**productivity**, per brood (default): `year`, `fledglings`, optional
`mother_age` (yearling vs older). Or aggregated per year: `year`, `broods`,
`fledglings` (total), optional `mother_age`. If mother age is known for
>= 80% of broods, fecundity is age-specific (`age_specific_fecundity: auto`).
Failed broods (0 fledglings) MUST be included.

**covariates** (optional), one row per year. A covariate in year t acts on
survival and immigration over the interval t -> t+1, and on fecundity in
year t.

## Config

```yaml
data:
  species: red-backed shrike
  count_unit: breeding females
  counts: shrike_counts.csv          # year_col (year), count_col (count)
  captures: shrike_captures.csv
  capture_format: long               # long | wide | marray
  id_col: ring                       # long/wide
  age_col: age
  capture_year_col: year             # long
  release_year_col: release_year     # marray
  productivity: shrike_broods.csv
  productivity_format: broods        # broods | aggregated (auto-detected from a 'broods' column)
  fledglings_col: fledglings
  broods_col: broods                 # aggregated only
  mother_age_col: female_age
  age_specific_fecundity: auto       # auto | true | false
  covariates: sahel_rainfall.csv     # covariate_year_col (year)
  effects:                           # rate: [covariates]; rates: juvenile_survival, adult_survival, fecundity, immigration
    juvenile_survival: [sahel_rain_mm]
  first_year: 2006                   # optional trimming
  last_year: 2020
settings:
  components: ipm                    # ipm | separate
  immigration: true                  # false only if the population is known to be closed (island)
  temporal_variation: true           # random year effects on all rates (needed for LTRE)
  count_model: poisson               # poisson | normal
  sex_ratio: 0.5
  fmax: 20                           # upper limit of the fecundity prior
  niter: 150000                      # per chain; nburnin 30000; thin 30; nchains 3 (defaults)
  compare_separate: true             # also fit separate analyses and report the precision gain
  ppc_draws: 600                     # posterior predictive checks (0 = off)
```

## Outputs

Estimates (posterior mean, SD, 95% CrI, R-hat, bulk ESS):

- `juvenile_survival`, `adult_survival` (mean apparent survival, at average
  covariate), `fecundity_yearling`, `fecundity_adult` (or `fecundity`),
  `recapture_prob`, `immigration_rate` (immigrant females per female already
  present, averaged over the study years; recovered within 5% at ~12% immigration
  but biased about +19% upward when immigration is low, ~8%: VALIDATION.md),
  `immigration_rate_median` (rate in
  a typical year, exp of the log-scale mean; depends on the poorly identified
  immigration SD and is NOT interpreted on its own).
- Temporal SDs on the link scale: `*_sd` (descriptive only: with <= 15 years
  they are biased upward, see VALIDATION.md).
- `count_sd` (only with `count_model: normal`): count observation SD, biased
  upward by about 12-16% at 18 years (coverage 1.00); a nuisance parameter,
  not interpreted.
- Covariate effects on the link scale: `<rate>_beta_<covariate>` (per 1 SD).
- Population: `N_first`, `N_final`, `N_mean`, `lambda_geomean` (geometric
  mean growth rate), `immigrants_mean`; annual `N_total[year]`,
  `local_recruits[year]`, `surviving_adults[year]`, `immigrants[year]`,
  `lambda[year]`, annual rates `juvenile_survival[year]` etc.
- LTRE: `ltre_juvenile_survival`, `ltre_adult_survival`, `ltre_immigration`,
  `ltre_fecundity_*`, `ltre_population_structure`, `var_lambda_explained`.
  Interpret the ranking only: absolute contributions are biased upward and
  the immigration contribution is overestimated (VALIDATION.md; Paquet et al. 2021).
- Interpretation: immigrant share of breeders, growth rate without
  immigration (sink/source), covariate effects on the response scale, IPM
  vs separate precision.

Figures: `population_trajectory`, `population_composition`, `growth_rate`,
`vital_rates`, `ltre`, `transient_ltre`, `ipm_vs_separate`,
`effect_<covariate>`, `gof` (posterior predictive checks for each data set),
`trace`, `posterior`.

## Simulation and power

`design`: `n_years` (15), `N_initial` (100 females), `marked_juveniles` (50
per year), `marked_adults` (20 per year), `prop_broods_monitored` (0.5),
`count_model` (poisson/normal; `count_sd`), `immigration`,
`age_specific_fecundity`, `temporal_variation`, `covariate` (one standardised
year covariate on juvenile survival). `params`: mean rates
(`juvenile_survival` 0.25, `adult_survival` 0.48, `fecundity_yearling` 2.8,
`fecundity_adult` 3.6, `immigration_rate` 0.12, `recapture_prob` 0.6),
temporal SDs and `juvenile_survival_beta` (-0.4 per SD). These are close to
the hoopoe / woodchat shrike values in Schaub & Kéry (2022). The defaults give a
stable population of ~100 females with ~13% immigrants each year.

Power is computed for `juvenile_survival_beta_winter_severity` (does a year
covariate affect first-year survival). VALIDATION.md has a power run over
`marked_juveniles` (15/30/60 per year) and the short-series (6-year) result.

## Pitfalls

- **Fewer than ~8 years**: temporal SDs, immigration SD and LTRE are
  imprecise (VALIDATION.md, short series). Report N and mean rates only.
- **Omitting immigration in an open population** biases juvenile survival
  (and adult survival) upwards, because the model must create the extra
  recruits from local survival and reproduction (VALIDATION.md, negative
  control: juvenile survival +19%, coverage 0.53).
- **Low immigration is overestimated.** With ~8% immigrants per year, the
  study-period immigration rate came out about 19% too high in simulation
  (immigrants_mean +14%; coverage still 0.98), because the immigration SD is
  prior-driven. Report low immigration estimates as upper values.
- **Juvenile survival with few juvenile recaptures** comes mostly from the
  counts. That is the IPM's strength, but it is then confounded with
  immigration: with few marked young, immigration and juvenile survival
  trade off (negative posterior correlation). Check the `ipm_vs_separate`
  figure.
- Counts of pairs or territories = females only if every territory has one
  breeding female. Polygyny or non-breeding floaters break this.
- Fecundity data must include failed nests; otherwise f is too high and
  immigration too low.
- The "separate" analysis cannot estimate N or immigration. Its rows are
  only for the precision comparison.

## References

- Besbeas P, Freeman SN, Morgan BJT, Catchpole EA (2002) Integrating mark-recapture-recovery and census data to estimate animal abundance and demographic parameters. Biometrics 58: 540-547.
- Abadi F, Gimenez O, Ullrich B, Arlettaz R, Schaub M (2010) Estimation of immigration rate using integrated population models. Journal of Applied Ecology 47: 393-400.
- Abadi F, Gimenez O, Arlettaz R, Schaub M (2010b) An assessment of integrated population models: bias, accuracy, and violation of the assumption of independence. Ecology 91: 7-14.
- Schaub M, Abadi F (2011) Integrated population models: a novel analysis framework for deeper insights into population dynamics. Journal of Ornithology 152 (Suppl 1): S227-S237.
- Schaub M, Fletcher D (2015) Estimating immigration using a Bayesian integrated population model: choice of parametrization and priors. Environmental and Ecological Statistics 22: 535-549.
- Schaub M, Kéry M (2022) Integrated Population Models: Theory and Ecological Applications with R and JAGS. Academic Press. Code: github.com/mikemeredith/IPM_code; R package IPMbook.
- Koons DN, Iles DT, Schaub M, Caswell H (2016) A life history perspective on the demographic drivers of structured population dynamics in changing environments. Ecology Letters 19: 1023-1031.
- Koons DN, Arnold TW, Schaub M (2017) Understanding the demographic drivers of realized population growth rates. Ecological Applications 27: 2102-2115.
- Paquet M, Knape J, Arlt D, et al. (2021) Integrated population models poorly estimate the demographic contribution of immigration. Methods in Ecology and Evolution 12: 1899-1910. doi:10.1111/2041-210X.13667
- Riecke TV, Williams PJ, Behnke TL, Gibson D, Leach AG, Sedinger BS, Street PA, Sedinger JS (2019) Integrated population models: model assumptions and inference. Methods in Ecology and Evolution 10 (2019). (Marker effects and mortality heterogeneity bias supplemental parameters such as immigration and fecundity.)
- Burnham KP et al. (1987) Design and analysis methods for fish survival experiments based on release-recapture. AFS Monograph 5.
- de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. JCGS 26: 403-413.
