# Example data: western pond turtles at one pond complex, 2015-2022 (SIMULATED)

**These data are simulated.** They were generated from known "true" parameter
values so that a reader can check that the analysis recovers them. They are not
field data. The generating script is `dev/jolly-seber/make_example.R`. It calls
the module's own simulator (`js_simulate`, seed 2024) for the capture histories,
then `set.seed(11)` for the field details (tag numbers, dates, traps, body size).
The true values are also in `TRUTH.txt`.

## The (simulated) study

A small team hoop-traps western pond turtles (*Actinemys marmorata*) at one
pond complex every summer, 20 May to 3 August, from 2015 to 2022. Each new
turtle gets a PIT tag and a field ID (`AM-001`, ...). Crews record sex,
carapace length and mass. Some juveniles and a few adults could not be sexed;
their sex is recorded as `U`. Effort varied between years (80-163
trap-nights), because the number of crew members changed.

The researcher's questions:
1. How many turtles used the pond complex each year, and over the whole study?
2. How many new turtles joined the population each year (recruits plus immigrants)?
3. What share of turtles survive (and stay) from one year to the next? Do males survive less well than females?
4. Is the population stable, growing or declining?

The Jolly-Seber model fits this well. Turtles are individually marked. There
are 8 annual occasions. The population is open between years (turtles enter
and leave), and trapping within a summer is short relative to a turtle's life.

## Design

| dimension | value |
|---|---|
| occasions (trapping seasons) | 8 (2015-2022), one per year, all intervals 1 year |
| turtles caught | 175 (100 female, 65 male, 10 unknown sex) |
| capture rows | 459 (1-3 captures per turtle per season; collapsed to 322 turtle-seasons) |
| turtles caught per season | 34, 23, 47, 55, 38, 44, 46, 35 |
| new turtles per season | 34, 17, 35, 27, 16, 15, 16, 15 |
| seasons caught per turtle | 1: 88, 2: 44, 3: 30, 4: 11, 6: 2 |
| trapping effort | 80-163 trap-nights per season (mean 120) |

## Files and columns

### `turtle_captures.csv`: one row per capture (459 rows)
| column | meaning | used by the model? |
|---|---|---|
| `turtle_id` | field ID, `AM-xxx`, unique per turtle | yes (animal ID) |
| `pit_tag` | PIT tag number | no |
| `capture_date` | date of capture, `YYYY-MM-DD` | no (the season column defines the occasion) |
| `season` | trapping season (year) the capture belongs to | yes (occasion) |
| `trap` | hoop-trap ID (`HT01`-`HT14`) | no |
| `sex` | `F`, `M` or `U` (unknown); recoded to female / male; unknown is modelled as a mixture | yes (group) |
| `carapace_mm` | straight carapace length, mm | no |
| `mass_g` | body mass, g | no |

Repeat captures of a turtle within a season count as one capture for that season.

### `turtle_seasons.csv`: one row per season (8 rows)
| column | meaning | used by the model? |
|---|---|---|
| `season` | year; matches `season` in the captures file | yes (occasion key) |
| `start_date`, `end_date` | first and last day of trapping | no |
| `trap_nights` | total trap-nights that season | yes (capture-probability covariate, standardised) |
| `observers` | number of crew | no |

## True values (SIMULATED; from `TRUTH.txt`)

Generating parameters:

| parameter | true value |
|---|---|
| superpopulation size (turtles ever present 2015-2022) | 260 |
| share present at the first season (`entry_prob_1`) | 0.45 |
| share entering at each later season (`entry_prob_2..8`) | 0.0786 each |
| apparent annual survival, females | 0.82 |
| sex effect on survival (logit scale, male vs female) | -0.40 |
| apparent annual survival, males | 0.753 |
| share of males | 0.45 |
| capture probability at mean effort | 0.40 |
| effort effect on capture (logit scale, per SD of trap-nights) | +0.45 |
| share of caught turtles with unknown sex | 0.06 |

Realised values in this simulated population:

| season | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 |
|---|---|---|---|---|---|---|---|---|
| turtles alive and present (N_t) | 97 | 100 | 101 | 109 | 102 | 107 | 105 | 107 |
| new turtles since previous season (B_t) | - | 23 | 16 | 31 | 17 | 28 | 21 | 27 |

Realised mean growth rate (N_2022 / N_2015)^(1/7) = 1.014 per year.

Note that the trap-nights column was drawn by the simulator (mean 120, SD 30);
the capture probabilities used for each season depend on it through the
effort effect above.
