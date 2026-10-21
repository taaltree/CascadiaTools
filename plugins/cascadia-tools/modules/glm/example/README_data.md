# Example data: Wood Thrush point counts (`wood_thrush_counts.csv`)

**These data are SIMULATED.** No real birds were counted. The values were drawn
by a script with known parameters, so you can check whether the analysis gets
the right answer. Generator: `dev/glm/make_examples.R` (`set.seed(20261018)`;
the same script also makes `nest_hatching.csv` and `seagrass_cover.csv`).

## The (simulated) study

Wood Thrush (*Hylocichla mustelina*) is a forest-interior songbird of eastern
North America. It nests in the understorey of closed-canopy, mostly deciduous
forest, and it has declined over much of its range. A field crew wants to know
where the birds are most numerous, so that managers know which stands to protect.

The crew set up **45 point-count stations** (WT-01 to WT-45) across a forested
landscape. Each station was visited **4 times** between 20 May and 18 June
2025, about 9 days apart. On each visit one of three observers (AK, JL, MR)
counted every Wood Thrush seen or heard. Most counts lasted 10 minutes. Some were
cut to 5 minutes because of weather or time. At each station the crew also
recorded the forest type, the percent canopy cover and the elevation.

Question: **does the number of Wood Thrushes counted change with canopy cover,
elevation and forest type?**

## Design

| dimension | value |
|---|---|
| stations (groups) | 45 (15 deciduous, 15 mixed, 15 conifer) |
| visits per station | 4 |
| observations (rows) | 180 |
| survey length | 5 min (62 counts) or 10 min (118 counts) |
| station-level predictors | canopy cover 42-98 % (mean 70, SD 15.6); elevation 184-941 m (mean 568, SD 234); forest type |
| response | Wood Thrush count, 0-10 (mean 1.13; 50% of counts are zero) |

## Columns

| column | type | meaning |
|---|---|---|
| `station` | text | Station id, WT-01 ... WT-45. Same station on every visit, so it is used as a random intercept `(1|station)`. |
| `visit` | 1-4 | Visit number at that station. Not used in the model. |
| `date` | YYYY-MM-DD | Survey date. Not used in the model. |
| `observer` | text | Observer initials (AK, JL, MR). Assigned at random in the simulation and has no effect. Not used in the model. |
| `survey_min` | 5 or 10 | Survey length in minutes. Used as effort (offset), so the model is for birds per minute, reported per 10-minute count. |
| `forest_type` | deciduous / mixed / conifer | Dominant forest type at the station. Reference level: deciduous. |
| `canopy_cover_pct` | 0-100 | Percent canopy cover at the station (one value per station). |
| `elevation_m` | metres | Station elevation (one value per station). |
| `wood_thrush` | whole number | Number of Wood Thrushes counted on that visit (the response). |

## TRUE values used to simulate the data (simulated, known only because this is an example)

Model: count ~ negative binomial (NB2) with mean
`mu = survey_min x exp(b0 + b_canopy (canopy - 68) + b_elev (elev - 560) + b_mixed [mixed] + b_conifer [conifer] + u_station)`,
`u_station ~ Normal(0, 0.45^2)`, NB size theta = 3.

| parameter | true value (per raw unit) | on the scale the module reports | as a rate ratio |
|---|---|---|---|
| baseline (deciduous, canopy 68 %, 560 m) | 0.12 birds per minute | intercept at the sample means = -2.08 (log birds/min); 1.25 birds per 10-min count | - |
| canopy cover | +0.025 per percentage point (log scale) | +0.390 per 1 SD (15.6 points) | 1.28 per 10 points; 1.48 per SD |
| elevation | -0.0012 per m | -0.281 per 1 SD (234 m) | 0.887 per 100 m; 0.755 per SD |
| mixed vs deciduous | -0.25 | -0.25 | 0.779 |
| conifer vs deciduous | -0.90 | -0.90 | 0.407 |
| SD of station random intercepts | 0.45 (log scale) | 0.45 | - |
| NB size theta | 3 | 3 | - |

The conversion to the reported scale is done by `dev/glm/si_truth.R` (output
`results/glm/si/truth_fitted_scale.csv`). The module standardises continuous
predictors with the sample mean and SD, so per-SD effects are `per-unit x SD`
and the intercept is the value at the sample means (canopy 70.0 %, elevation 568 m).
