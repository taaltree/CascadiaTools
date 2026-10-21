# Example data: walleye fall gill-net ages (SIMULATED)

**These data are simulated.** They are not from a real lake. They were made with
the module's own simulator so that the true answer is known.

## The story

A state fisheries agency runs a standard fall experimental gill-net survey on a
two-basin lake, 22 September to 10 October 2025. Every walleye caught is measured
(total length) and aged from otoliths. The biologist wants to know how fast
adult walleye are dying: total annual mortality, all causes together (angling
plus natural deaths).

## Design

| dimension | value |
|---|---|
| fish sampled | 350 |
| fish with a readable age | 348 (2 otoliths unreadable, age = NA) |
| ages observed | 1-13 (no age-12 fish) |
| basins | North Basin (178 fish), South Basin (172) |
| gear | experimental gill net |
| oldest age class in simulated population (`max_age`) | 20 |

Catch at age (aged fish): age 1: 11, 2: 66, 3: 129, 4: 60, 5: 24, 6: 31, 7: 4,
8: 11, 9: 4, 10: 2, 11: 5, 12: 0, 13: 1.

## Columns (`walleye_fall_gillnet_ages.csv`, one row per fish)

| column | meaning |
|---|---|
| `fish_id` | tag, WAE25-0001 ... WAE25-0350 |
| `sample_date` | date of the net set (YYYY-MM-DD) |
| `basin` | North Basin or South Basin (assigned at random; no basin effect simulated) |
| `gear` | always "experimental gill net" |
| `total_length_mm` | total length (mm), from a von Bertalanffy curve (Linf 650 mm, k 0.3) with 7% lognormal noise. Not used by the catch curve |
| `age` | otolith age in years; NA = unreadable |

## TRUE values (simulated)

| parameter | true value | meaning |
|---|---|---|
| Z | 0.55 per year | total instantaneous mortality, constant over ages |
| annual survival S = exp(-Z) | 0.577 | |
| annual mortality A = 1 - S | 0.423 | |
| recruit_sd_cc | 0.4 | lognormal SD of year-class strength (no trend, no autocorrelation) |
| sel_age50 | 2.2 years | gill-net logistic selectivity, age at 50% vulnerability |
| sel_slope | 3 | logistic selectivity slope (fish are nearly fully vulnerable by age 3-4) |

## How it was made

`dev/fisheries-stock/make_examples.R` (section 1): `module$simulate()` with the
values above, `n_fish = 350`, seed 42. Then, with `set.seed(42)`, fish were shuffled
and lengths, basins and dates were added, and the ages of rows 17 and 203 were set
to NA. A rerun on 2026-10-20 reproduced the shipped CSV exactly
(`logs/fisheries-stock_si_check.log`). The true values are also saved in
`dev/fisheries-stock/example_truths.rds` (element `cc`).

# Example data: smallmouth bass ages from six fall surveys (SIMULATED, multi-year)

`multiyear_ages.csv` is **simulated** with the module's own multi-year simulator,
so the true answer is known. Story: a river smallmouth bass population is sampled
by boat electrofishing each fall from 2019 to 2024; about 200 fish a year are aged
from scales/otoliths. The biologist wants total annual mortality of adults.

| column | meaning |
|---|---|
| `survey_year` | year of the fall survey (2019-2024) |
| `age` | age in years |
| `n_fish` | number of fish of that age in that year's aged sample (rows with 0 fish are omitted) |

Design: 6 years x 200 aged fish (1,200 fish), ages 0-16, logistic selectivity
(age at 50% 2.3 years, slope 3), oldest age in the population 16. The same year
classes are followed from year to year (the weak 2017 year class shows up as low
counts at age 3 in 2020, age 4 in 2021, age 5 in 2022 ...).

| parameter | true value |
|---|---|
| Z | 0.6 per year (annual survival 0.549, annual mortality 0.451) |
| recruit_sd_cc | 0.5 (lognormal SD of year-class strength; no trend, no autocorrelation) |
| Z_trend | 0 (mortality constant over the six years) |

Made by `dev/fisheries-stock/make_multiyear_example.R` (`module$simulate()` with
`cc_years = 6`, `n_fish = 200`, `max_age = 16`, `first_year = 2019`, seed 2024).
Fit with `config_multiyear.yaml`; planted-truth check with `recovery_multiyear.yaml`
(results in `results/fisheries-stock/example_multiyear/`).
