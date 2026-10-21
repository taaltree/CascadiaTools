# Example data: Pacific marten hair-snare survey, 2020-2024 (SIMULATED)

**These data are simulated.** No real martens were sampled. The files were made with the
module's own simulator so that the true values are known and the analysis can be checked
against them. The study area, species and field protocol are realistic, but invented.

## The story

A forest carnivore biologist runs a DNA hair-snare survey for Pacific marten (*Martes caurina*)
in the Oregon Cascades. Each summer from 2020 to 2024 the crew sets 49 baited hair snares on a
1-km grid and checks them 5 times, once a week. Hair caught on the snare's glue strips is sent to
a genetics lab, which identifies the individual marten. Snares are sometimes knocked down by bears
or snow-load, so a few stopped working part-way through a summer.

The biologist wants to know:

1. How many martens live in the study area each summer, and is the population declining?
2. What is the year-to-year (apparent) survival, and how many new martens arrive each year?
3. Is marten density higher under denser forest canopy?

## Design

| dimension | value |
|---|---|
| detectors | 49 hair snares, 7 x 7 grid, 1,000 m apart (6 km x 6 km) |
| primary periods (years) | 5 summers, 2020-2024 |
| occasions per year | 5 weekly checks |
| snare outages | 2 / 0 / 2 / 3 / 5 snares not working on every check in 2020-2024 |
| habitat grid | 900 pixels, 400 m apart, covering the snares plus a 3 km buffer (12 km x 12 km = 144 km2) |
| density covariate | canopy cover (%), 24-76% across the grid, mean 50% |
| martens detected | 95 different animals; 32 / 38 / 39 / 41 / 48 per year; 48 seen in two or more years |
| detections | 598 genotyped samples (one per marten, snare and check) |

Coordinates are UTM zone 10N (metres).

## Files and columns

### `hair_snares.csv` (one row per snare; 49 rows)

| column | meaning |
|---|---|
| `snare_id` | snare label, HS01-HS49 |
| `utm_e`, `utm_n` | snare location, UTM easting and northing (m) |
| `elevation_m` | snare elevation (m); recorded by the crew, not used in the analysis |
| `ops_2020` ... `ops_2024` | which of the 5 checks the snare was working that summer, as a 0/1 string. `11111` = all checks; `11100` = working on checks 1-3, then destroyed |

### `marten_genotypes.csv` (one row per genotyped hair sample; 598 rows)

| column | meaning |
|---|---|
| `sample_id` | lab sample number |
| `year` | summer of the survey (2020-2024) |
| `check` | weekly check within that summer, 1-5 |
| `snare_id` | snare where the hair was collected (matches `hair_snares.csv`) |
| `marten_id` | individual identity from the genotype (e.g. MACA-476) |
| `n_hairs` | number of hairs in the sample; not used in the analysis |

### `habitat_grid.csv` (one row per habitat pixel; 900 rows)

| column | meaning |
|---|---|
| `x`, `y` | pixel centre, UTM easting and northing (m) |
| `canopy_cover_pct` | canopy cover (%) of the pixel |

## How the data were generated

Script: `dev/scr-open/make_example.R`. It calls

    module$simulate(module$defaults$params, module$defaults$design, seed = 2025)

and then rewrites the simulated detectors and captures in field-crew format (UTM offsets
571000 E / 4861000 N, years 2020-2024, random marten and sample IDs with `set.seed(7)`,
made-up `elevation_m` and `n_hairs` columns). The true values are saved in
`dev/scr-open/example_truth.rds`.

## TRUE values (simulated; known only because the data are simulated)

Parameters used to simulate:

| parameter | true value |
|---|---|
| density in 2020 (expected) | 0.8 martens per km2 |
| baseline detection p0 (per check, snare at the centre of the home range) | 0.30 |
| sigma (spatial scale of the home range) | 800 m |
| apparent survival between years | 0.80 |
| per-capita recruitment (expected new martens per marten alive the year before) | 0.25 |
| canopy effect on density (log scale, per 1 SD of canopy cover) | 1.0 |

Realised values in this one simulated population:

| quantity | 2020 | 2021 | 2022 | 2023 | 2024 |
|---|---|---|---|---|---|
| martens alive in the 144 km2 area | 127 | 131 | 133 | 145 | 158 |
| density (per km2) | 0.88 | 0.91 | 0.92 | 1.01 | 1.10 |
| recruits (new martens) | - | 26 | 31 | 35 | 43 |
| per-capita recruitment | - | 0.20 | 0.24 | 0.26 | 0.30 |
| growth rate (N_t / N_t-1) | - | 1.03 | 1.02 | 1.09 | 1.09 |

Superpopulation (different martens alive at any time in 2020-2024): **262**.
95% home-range area: pi (2.448 x 800 m)^2 = **12.0 km2**.

Note: in this particular realisation the martens happen to favour canopy more strongly than
the simulation value. Fitting the density model to the TRUE activity centres gives a canopy
effect of 1.18 (SE 0.08), not 1.0. The snare grid also sits in below-average canopy
(VALIDATION.md, section 8).
