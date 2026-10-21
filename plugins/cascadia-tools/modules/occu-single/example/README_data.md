# Example data: pileated woodpecker call-broadcast surveys (SIMULATED)

**These data are simulated.** No real birds were surveyed. The data were made by
`dev/occu-single/make_example.R` (R, `set.seed(2026)`, static block) so that the true values
are known and the analysis can be checked against them. The study story below is invented
to make the example realistic.

This folder also holds a second, dynamic (multi-season) example
(`dynamic_*.csv`, `config_dynamic*.yaml`; simulated chorus-frog surveys at 80 ponds over
5 years). This file describes only the static woodpecker example used in the worked
example (`WALKTHROUGH.md`).

## The (simulated) study

A forest ecologist wants to know where pileated woodpeckers (*Dryocopus pileatus*) occur in a
managed forest. In May and June 2025 her crew visited 120 forest stands. At each visit an
observer played a recorded call and listened for a reply. A stand got up to 4 visits, about
12 days apart. Some stands got only 2 or 3 visits because of weather and road closures.
Three observers (AB, JK, MR) did the surveys. Woodpeckers call less as the season goes on,
and observers differ in how well they hear replies. So a stand with no detections may still
have had woodpeckers.

Question: **Is woodpecker occupancy linked to forest cover and elevation, once we allow for
birds that were present but missed?**

## Design

| dimension | value |
|---|---|
| sites (forest stands) | 120 |
| visits per stand | up to 4 (72 stands had 4, 21 had 3, 27 had 2) |
| surveys done | 405 of 480 possible (15.6% not done) |
| survey dates | 1 May to 19 June 2025 (day of year 121 to 170) |
| observers | 3 (AB 168 surveys, JK 148, MR 89) |
| detections | 52 of 405 surveys, at 34 stands (naive occupancy 0.28) |
| stands truly occupied (simulated truth) | 45 (0.375) |

## Files and columns

### `detections.csv` (one row per survey done; 405 rows)

| column | meaning |
|---|---|
| `site` | stand ID, `STAND-001` to `STAND-120` |
| `survey` | visit number at that stand (1, 2, 3, 4) |
| `date` | survey date (YYYY-MM-DD); not used in the model (`julian_day` is) |
| `julian_day` | survey date as day of year (121 = 1 May) |
| `observer` | observer initials: AB, JK or MR (AB is the reference level) |
| `detected` | 1 = woodpecker replied or was heard, 0 = not detected |

Surveys that were not done are simply absent (no row). They are not zeros.

### `sites.csv` (one row per stand; 120 rows)

| column | meaning |
|---|---|
| `site` | stand ID (matches `detections.csv`) |
| `forest_cover_pct` | percent forest cover around the stand (4 to 98 in these data) |
| `elevation_m` | elevation in metres (185 to 1395 in these data) |
| `stand_age_class` | young / mature / old. Random noise; it has **no** effect in the simulation and is not in the model |

## TRUE values (simulated)

Generating model (raw covariate units):

- occupancy: `logit(psi) = -0.4 + 1.1 * (forest - 55)/25 - 0.5 * (elevation - 750)/350`
- detection per survey: `logit(p) = -0.3 - 0.6 * (julian_day - 140)/15 + observer effect`,
  with observer effect AB 0, JK -0.6, MR +0.4.

The module reports numeric effects per 1 standard deviation (SD) of the covariate, and the
intercepts as probabilities at the average covariate values (observer at AB). The truth on
that scale (computed by `dev/occu-single/si_truth.R`, saved in
`results/occu-single/si/truth.json`):

| quantity | true value |
|---|---|
| occupancy probability at mean forest cover (58.2%) and mean elevation (834 m) | 0.406 |
| forest cover effect, logit per 1 SD (1 SD = 23.9 percentage points) | +1.05 |
| elevation effect, logit per 1 SD (1 SD = 333 m) | -0.475 |
| detection per survey on the mean date (day 142), observer AB | 0.402 |
| survey date effect, logit per 1 SD (1 SD = 13.3 days) | -0.532 |
| observer JK vs AB (logit) | -0.6 |
| observer MR vs AB (logit) | +0.4 |
| proportion of the 120 stands truly occupied | 0.375 (45 stands) |

To regenerate the CSVs exactly: `scripts/Rscript dev/occu-single/make_example.R` (this also
rewrites the dynamic example). `si_truth.R` re-runs the static block and checks that the
regenerated files are byte-identical to the ones here.
