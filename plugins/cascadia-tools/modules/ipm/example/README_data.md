# Example data: red-backed shrike study plot, 2006-2020 (SIMULATED)

**These data are simulated.** They were generated with known "true" parameter
values so that a reader can check that the analysis recovers them. They are not
field data. The generating script is `dev/ipm/make_example.R`. It calls the
module's own simulator (`ipm_simulate`, the same process as `module$simulate()`)
with `seed = 2024`, and then writes the result in the raw formats a field team
keeps (`set.seed(7)` is used only to draw ring numbers and to split each year's
fledglings among nests). The true values are also in `TRUTH.txt`.

## The (simulated) study

The red-backed shrike (*Lanius collurio*) is a small migratory songbird. It breeds
in hedgerows and extensive farmland in Europe and winters in southern Africa. It
passes through the Sahel on migration. A research group has followed the shrikes on
one study plot for 15 breeding seasons (2006-2020). The analysis is female-based:
every number below is about females.

Each year the team does three things:

1. **Counts.** In June they map all territories and count the breeding females on
   the plot. In 2013 the survey was not done (the count is empty).
2. **Ringing.** They ring 50 nestlings per year at the nest (2006-2019) and catch and
   ring 20 unringed breeding females per year at the nest. A female caught again in a later year is
   recorded as a recapture. Ringed nestlings that come back to breed are caught as
   adults.
3. **Nest monitoring.** They follow about half of all nests to fledging and record
   the number of fledglings, including nests that failed (0 fledglings). They also
   note whether the mother is a yearling (first-time breeder) or older.

A fourth file gives July-September rainfall in the Sahel for each year. It is used
as a possible driver of first-year survival: young birds cross the Sahel on their
first migration.

The researcher's questions:
1. How did the number of breeding females change over the 15 years?
2. Which vital rates (first-year survival, adult survival, breeding success,
   immigration) drive the ups and downs?
3. How many breeders are immigrants? Would the population hold up without them?
4. Is first-year survival higher after wet years in the Sahel?

## Design

| dimension | value |
|---|---|
| years | 15 (2006-2020), one breeding season per year |
| survival intervals | 14 (2006-07 to 2019-20) |
| counts | 14 of 15 years (no survey in 2013); 99-209 females counted |
| nestlings ringed | 50 per year, 2006-2019 (700 in total) |
| adult females ringed (first capture) | 280 (20 new females per year, 2006-2019) |
| capture rows | 1,309 (980 rings; one row per ring per season) |
| rings caught again in a later year | 229 (127 of them ringed as nestlings) |
| nests monitored | 1,020 (549 with yearling mothers, 471 with older mothers); 47 failed (0 fledglings) |
| fledglings per nest | 0-11; mean 2.91 (yearling mothers) and 3.62 (older mothers) |
| Sahel rainfall | 191-522 mm, mean 420 mm, SD 90 mm |

## Files and columns

### `shrike_counts.csv`: one row per year (15 rows)
| column | meaning |
|---|---|
| `year` | breeding season |
| `females_counted` | breeding females counted on the plot in June; empty = no survey (2013) |

### `shrike_captures.csv`: one row per ring per breeding season (1,309 rows)
| column | meaning |
|---|---|
| `ring` | ring number; `NX#####` = ringed as a nestling, `AF#####` = ringed as an adult female |
| `year` | breeding season of the capture |
| `age` | age at THAT capture: `nestling` (young of the year) or `adult` (breeding female). Only the age at first capture is used. |
| `sex` | `U` for nestlings (sex unknown), `F` for breeding females; not used by the model |

### `shrike_broods.csv`: one row per monitored nest (1,020 rows)
| column | meaning |
|---|---|
| `year` | breeding season |
| `nest_id` | nest label, `N<year>-<number>` |
| `female_age` | `yearling` (first-time breeder) or `older` |
| `fledglings` | number of young that left the nest (both sexes); 0 = failed nest |

### `sahel_rainfall.csv`: one row per year (15 rows)
| column | meaning |
|---|---|
| `year` | year |
| `sahel_rain_mm` | July-September rainfall in the Sahel (mm). The value for year t acts on survival from t to t+1, so the 2020 value is not used. |

## How the data were generated (the population model)

Female-based, pre-breeding census, two age classes plus immigrants (the model in
`modules/ipm/README.md`):

- local recruits: `N1[t+1] ~ Poisson(0.5 * phij[t] * (f1[t] * (N1[t] + N3[t]) + f2[t] * N2[t]))`
- surviving adults: `N2[t+1] ~ Binomial(Ntot[t], phia[t])`
- immigrants: `N3[t+1] ~ Poisson(omega[t] * Ntot[t])`
- counts: `C[t] ~ Poisson(Ntot[t])`
- productivity: fledglings ~ Poisson(broods x fecundity), by mother age; half of
  all broods are monitored (Binomial(females, 0.5)).
- capture-recapture: every year 50 nestlings and 20 adult females are marked;
  survival as above; recapture probability 0.6 for every bird alive.
- Annual rates vary at random around their means (normal year effects on the
  logit or log scale). The rainfall series was drawn as a standard normal
  variable and then put into millimetres as 420 + 90 x value.

The three data sets are drawn independently from the same annual rates. This is
the independence assumption of the IPM.

## TRUE parameter values (simulated)

Mean rates (the values the model estimates):

| parameter | true value |
|---|---|
| first-year (juvenile) survival, at mean rainfall | 0.25 |
| adult survival | 0.48 |
| fecundity, yearling mothers (fledglings per female) | 2.8 |
| fecundity, older mothers | 3.6 |
| recapture probability | 0.6 |
| rainfall effect on first-year survival (logit, per SD = 90 mm) | +0.40 (odds ratio 1.49) |
| immigration rate, study-period mean of the 14 annual rates | 0.127 immigrant females per female present |
| immigration rate, median (typical year) | 0.12 |
| year-to-year SD (link scale): juvenile survival / adult survival / fecundities / immigration | 0.30 / 0.20 / 0.15, 0.15 / 0.30 |

Realised population (what actually happened in this simulated data set):

| quantity | true value |
|---|---|
| breeding females 2006 (N_first) | 86 (the 2006 count was 106) |
| breeding females 2020 (N_final) | 136 |
| mean number of females (N_mean) | 137.3 |
| geometric mean growth rate | 1.033 per year |
| immigrants per year (mean, 2007-2020) | 17.5, i.e. 12.8% of breeders |
| realised means of the annual rates | phij 0.255, phia 0.473, f1 2.84, f2 3.76, omega 0.127 |

Number of breeding females each year (true):

| 2006 | 2007 | 2008 | 2009 | 2010 | 2011 | 2012 | 2013 | 2014 | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 86 | 126 | 112 | 138 | 105 | 99 | 130 | 184 | 206 | 175 | 165 | 129 | 118 | 150 | 136 |

True life-table response experiment (LTRE) contributions to the variance of the
growth rate: juvenile survival 0.0191, fecundity of older mothers 0.0032,
fecundity of yearlings 0.0027, adult survival 0.0009, immigration 0.00003,
population structure -0.0002 (total 0.0258). First-year survival is by far the
main driver of the year-to-year changes.
