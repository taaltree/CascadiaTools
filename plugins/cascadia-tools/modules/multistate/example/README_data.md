# Example data: chytrid infection in a marked frog population, 2015-2022 (SIMULATED)

**These data are simulated.** They were generated with known "true" parameter
values so that a reader can check that the analysis recovers them. They are not
field data. The generating script is `dev/multistate/make_example.R`. It calls the
module's own simulator, `module$simulate(params, design, seed = 2024)`, and then
writes the result in the long CSV format a field crew would keep. The true values
are saved in `truth.rds` and listed below.

## The (simulated) study

A pond-breeding frog population is surveyed once a year, at the peak of the
breeding season, from 2015 to 2022. Each frog caught gets a PIT tag (ids
`RL-0001` to `RL-0280`) and has its snout-vent length measured. Each capture
also gets a skin swab, which qPCR tests for the amphibian chytrid fungus
*Batrachochytrium dendrobatidis* (Bd). The result is `uninfected` or
`infected`. Some swabs are lost, contaminated or inconclusive. Those captures
are coded `U`: the frog was seen, but its infection status is unknown.

Frogs can clear Bd and become uninfected again, and they can be reinfected.
Bd grows best in cool water, so infection might be more common after cool years.

The researcher's questions:
1. Do infected frogs survive less well than uninfected frogs from one year to the next?
2. How often do uninfected frogs become infected, and how often do infected frogs clear the infection?
3. Is water temperature associated with infection and clearance? Do larger frogs survive better?

The lost and inconclusive swabs (`U`) mean this needs a **multievent** model
(Pradel 2005). Dropping those captures, or guessing their status, would bias the
results (see VALIDATION.md, negative control).

## Design

| dimension | value |
|---|---|
| states | 2: `uninfected`, `infected` (plus "dead", which is never observed) |
| occasions (annual surveys) | 8 (2015-2022) |
| survival / transition intervals | 7 (2015-16 to 2021-22), each 1 year |
| frogs tagged | 280: 40 new frogs in each year 2015-2021, none new in 2022 |
| capture rows | 516 (one per frog per year caught) |
| status records | 258 uninfected, 144 infected, 114 unknown (`U`, 22%) |
| frogs recaptured in a later year | 141 of 280 (50%); captures per frog 1-6 |
| snout-vent length | mean 52.0 mm, SD 6.0 mm, range 31.8-66.5 mm |
| water temperature | mean 17.0 deg C, SD 2.0 deg C, range 13.9-20.4 deg C |
| survey effort | mean 6.0 person-hours, SD 1.5, range 3.9-8.5 |

## Files and columns

### `frog_captures.csv`: one row per capture (516 rows)
| column | meaning |
|---|---|
| `frog_id` | PIT-tag id, `RL-0001` to `RL-0280` |
| `year` | survey year, 2015-2022 (one survey per year = one occasion) |
| `bd_status` | `uninfected`, `infected`, or `U` (swab lost or inconclusive: seen, status unknown) |

Frogs not caught in a year have no row for that year.

### `frog_individuals.csv`: one row per frog (280 rows)
| column | meaning |
|---|---|
| `frog_id` | matches `frog_captures.csv` |
| `svl_mm` | snout-vent length at first capture, mm (used as a covariate on survival) |

### `survey_years.csv`: one row per survey year (8 rows)
| column | meaning |
|---|---|
| `year` | survey year |
| `water_temp_C` | mean pond water temperature over the following year, deg C. The value in year t is used for survival and infection/clearance over t to t+1. The 2022 value is not used for those, because there is no later interval. |
| `survey_hours` | survey effort that year, person-hours (used as a covariate on detection in year t) |

### Other files
- `config.yaml`: the analysis settings the agent writes (columns, states, covariates, MCMC length, figure labels).
- `recovery.yaml`: the planted-truth recovery at this design (same truths as below).
- `truth.rds`: the true values as saved by the generating script.

## True parameter values (SIMULATED; used to generate the data)

Probabilities are at average covariate values. Slopes are on the logit scale
(multinomial logit for transitions) per 1 SD of the covariate.

| parameter | true value | meaning |
|---|---|---|
| `survival_uninfected` | 0.75 | annual survival of uninfected frogs |
| `survival_infected` | 0.55 | annual survival of infected frogs |
| `transition_uninfected_to_infected` | 0.30 | P(infected next year \| uninfected now and survives) |
| `transition_infected_to_uninfected` | 0.35 | P(uninfected next year \| infected now and survives): clearance |
| `detection_uninfected` | 0.60 | P(caught \| alive and uninfected) |
| `detection_infected` | 0.50 | P(caught \| alive and infected) |
| `assignment_uninfected` | 0.75 | P(swab gives a result \| caught, uninfected) |
| `assignment_infected` | 0.85 | P(swab gives a result \| caught, infected) |
| `initial_uninfected` | 0.65 | share of newly tagged frogs that are uninfected |
| `initial_infected` | 0.35 | share of newly tagged frogs that are infected |
| `stationary_uninfected` | 0.538 | derived: long-run share of time uninfected (0.35 / (0.30 + 0.35)) |
| `stationary_infected` | 0.462 | derived |
| `survival_beta_svl_mm` | +0.4 | larger frogs survive better (per 1 SD, about 6 mm) |
| `transition_beta_water_temp_C_uninfected_to_infected` | -0.7 | warmer years: fewer new infections (per 1 SD, about 2 deg C) |
| `transition_beta_water_temp_C_infected_to_uninfected` | +0.5 | warmer years: more clearance |
| `detection_beta_survey_hours` | +0.4 | more effort: higher capture probability (per 1 SD, about 1.5 h) |

The simulator draws standardised covariates (mean 0, SD 1). The CSVs hold them
back-transformed to field units: `svl_mm = 52 + 6 z`, `water_temp_C = 17 + 2 z`,
`survey_hours = 6 + 1.5 z`. The module re-standardises them before fitting, so
the slopes are per SD of the covariate as it appears in the files. Those SDs are
within 0.1% of the simulated SDs, so the true slopes above apply directly.

**Seed.** `make_example.R` uses seed 2024. Rerun it with
`scripts/ejob -c 1 -m 4 -- scripts/Rscript dev/multistate/make_example.R` to regenerate the files.
