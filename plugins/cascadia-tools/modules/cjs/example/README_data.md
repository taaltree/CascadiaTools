# Example data: song sparrow banding station, 2016-2024 (SIMULATED)

**These data are simulated.** They were generated with known "true" parameter
values so that a reader can check that the analysis recovers them. They are not
field data. The generating script is `dev/cjs/make_example.R` (`set.seed(2031)`).
The true values are also in `TRUTH.txt`.

## The (simulated) study

A constant-effort mist-netting station catches adult song sparrows
(*Melospiza melodia*) each breeding season, May to August, from 2016 to 2024.
Every new adult gets a numbered aluminium band. Banders record sex, body mass
and wing chord at first capture. Many birds are caught once and never again.
Some of them die. Others are passing migrants or floaters ("transients") that
never settle at the site.

The researcher's questions:
1. What is the annual apparent survival of adult song sparrows at the station?
2. Do males and females survive at different rates?
3. Is survival lower after cold winters, and higher for heavier birds?

## Design

| dimension | value |
|---|---|
| occasions (breeding seasons) | 9 (2016-2024), one occasion per year |
| survival intervals | 8 (2016-17 to 2023-24), all 1 year |
| birds banded | 575 (78, 70, 74, 66, 72, 69, 75, 71 new adults in 2016-2023; none needed in 2024) |
| sex ratio | 301 males (52%), 274 females |
| capture rows | 1,283 (1-3 captures per bird per season; collapsed to 813 bird-seasons) |
| birds recaptured in a later year | 158 (27%) |
| realised transients | 147 of 575 (25.6%) |
| body mass | mean 21.9 g, SD 1.66 g, range 15.8-27.7 g |
| winter temperature | mean January minimum -3.1 deg C, SD 2.04 deg C, range -6.5 to -0.4 deg C |

## Files and columns

### `sosp_captures.csv`: one row per capture (1,283 rows)
| column | meaning |
|---|---|
| `band_number` | band number, `2841-xxxxx` (unique per bird) |
| `capture_date` | date of capture, `YYYY-MM-DD`; all dates fall between 10 May and 8 August |
| `net` | mist net ID (`N01`-`N12`); recorded but not used in the model |
| `bander` | initials of the bander (`KL`, `MR`, `JT`); recorded but not used |

The module bins dates into one occasion per year (`occasion_unit: year`) and keeps
only May-August captures (`season_months: [5, 6, 7, 8]`). Repeat captures within a
season count as one.

### `sosp_birds.csv`: one row per bird (575 rows)
| column | meaning |
|---|---|
| `band_number` | matches `sosp_captures.csv` |
| `sex` | `F` or `M` (F is the reference level) |
| `mass_g` | body mass at first capture, grams |
| `wing_chord_mm` | wing chord at first capture, mm; not used in the model and has no effect in the simulation |
| `first_year` | year of banding (informational; the module takes first capture from the captures file) |

### `sosp_years.csv`: one row per year (9 rows)
| column | meaning |
|---|---|
| `year` | breeding season |
| `winter_min_temp_c` | mean January minimum temperature (deg C) of the winter AFTER that season, so the 2016 row applies to survival from 2016 to 2017. The 2024 value is not used. |
| `spring_precip_mm` | spring precipitation, mm; not used in the model and has no effect in the simulation |

## TRUE parameter values (simulated)

| parameter | true value | meaning |
|---|---|---|
| resident female survival | **0.55** | annual apparent survival of a resident female at mean mass and mean winter temperature |
| male effect on survival | **+0.30** logit | so resident male survival = **0.623** |
| winter temperature effect | **+0.40** logit per SD (2.04 deg C) | milder winters, higher survival |
| body mass effect | **+0.25** logit per SD (1.66 g) | heavier birds, higher survival |
| proportion transient | **0.25** (both sexes) | chance that a newly banded bird is a transient and never returns |
| female recapture probability | **0.45** | chance of catching a live resident female in a season |
| male effect on recapture | **+0.50** logit | so male recapture = **0.574** |

Effects per SD use the SD over the 9 years (winter) or over the 575 birds (mass).
In the simulation, transients survive the first interval with probability 0. Residents
survive every interval with the logistic-model probability. Detection does not depend on
year, net or bander.

Seed choice: seed 2031 is the first seed from 2025 upward whose fitted estimates are a
typical realisation, with every key estimate within 1.5 SE of its truth (`dev/cjs/pick_seed.R`).
The generator itself was checked to be unbiased over 100 seeds (`dev/cjs/exgen_check.R`).
