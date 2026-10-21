# Example data: coyote, gray fox and bobcat on 220 camera traps

**These data are simulated.** No real animals were recorded. They were made with
known ("true") parameter values, so you can check whether the analysis gets the
right answer. Script: `dev/occu-rota/make_example.R` (`set.seed(2025)`). Running
it again gives byte-identical CSV files (checked by `dev/occu-rota/si/si_truth.R`,
which regenerates both files into `dev/occu-rota/si/regen/` and compares them).

## The story

A mammal ecologist ran a summer camera-trap survey across a landscape of
forest, farmland and suburbs (coordinates are placed in northern Virginia, USA,
for realism only). She put 220 cameras out for 6 weeks, from 3 June to
14 July 2024, and scored each camera-week as "species seen" or "not seen".

Her question: **do gray foxes occur less often where coyotes are present, and
does that depend on housing density?** Coyotes kill and displace foxes, but
near houses foxes may find refuges, so the avoidance might weaken there.
Bobcats were also recorded. The main (two-species) example ignores bobcat;
`config_3species.yaml` adds it.

## Design

| dimension | value |
|---|---|
| species | coyote, gray fox (bobcat in the 3-species variant) |
| sites (camera stations) | 220 |
| surveys (weeks) per camera | 6 |
| camera-weeks not surveyed | 93 of 1320 (7.0%): 81 with 0 active days, 12 rows absent (6 cameras pulled after week 4) |
| site covariates | forest cover (%), housing density (units/km2), camera on trail (yes/no) |
| survey covariate | days the camera was active in that week (2-7) |

## Files and columns

`camera_detections.csv`: one row per camera x week (1308 rows).

| column | meaning |
|---|---|
| `camera_id` | camera station, `CAM-001` ... `CAM-220` |
| `week` | survey week, 1-6 |
| `start_date` | first day of that week (YYYY-MM-DD) |
| `days_active` | days the camera worked that week (0-7). 0 = camera down, so the week is treated as not surveyed |
| `coyote`, `gray_fox`, `bobcat` | 1 = photographed at least once that week, 0 = not photographed, blank = not surveyed |

`camera_sites.csv`: one row per camera (220 rows).

| column | meaning |
|---|---|
| `camera_id` | matches the detections file |
| `forest_pct` | forest cover within 500 m of the camera (%), 5-95 |
| `housing_per_km2` | housing units per km2 (log-normal, geometric mean 39.4, range 2.3-482) |
| `on_trail` | `yes` if the camera faces a trail or road, else `no` |
| `lat`, `lon` | location (decimal degrees, WGS84; for mapping only, not used in the model) |

## The TRUE values (simulated)

Occupancy came from a three-species Rota et al. (2016) model. Covariates were
standardised (mean 0, SD 1) over the 220 cameras: forest cover as is, housing
as log(housing). Natural parameters (logit scale):

| term | intercept | forest (per SD) | log housing (per SD) |
|---|---|---|---|
| coyote | 0.4 | -0.3 | +0.2 |
| gray fox | -0.6 | +0.7 | +0.3 |
| bobcat | -1.0 | +0.9 | -0.5 |
| coyote : gray fox (interaction) | **-0.9** | 0 | **+0.6** |
| coyote : bobcat | 0.3 | 0 | 0 |
| gray fox : bobcat | 0 | 0 | 0 |

So, at average housing, the odds that a fox occupies a site are multiplied by
exp(-0.9) = **0.41** where coyotes are present. Each SD of log housing (housing
multiplied by 2.7) raises the log odds ratio by 0.6, so the
avoidance fades in suburbs: odds ratio 0.41 at the geometric mean (39 houses/km2)
and about 1.0 (no association) at about 170 houses/km2 (1.5 SD above).

Detection per week (logit): qlogis(p0) + trail effect x (on_trail - 0.45) +
0.25 x (days_active - 7), with p0 = 0.35 / 0.30 / 0.22 and on-trail effects
+0.8 / -0.2 / +0.6 (coyote / gray fox / bobcat).

### Truth on the scale the two-species analysis reports

Summing bobcat out of the three-species model gives an exact two-species
model. The fox term and the coyote-fox interaction do not change. The coyote
term gains a small, smooth adjustment (between +0.01 and +0.26 across the
cameras), so its values below are a linear approximation at average conditions
(a least-squares fit across the 220 cameras gives 0.50 / -0.25 / +0.17). All
values come from `dev/occu-rota/si/si_truth.R` and are stored in
`dev/occu-rota/si/truth.json`.

Average conditions = forest 52.4%, housing 39.4 units/km2 (geometric mean),
6.5 active days, camera off trail.

| quantity | true value |
|---|---|
| coyote occupancy | 0.563 |
| gray fox occupancy | 0.258 |
| both present | 0.103 |
| P(fox present, coyote present) | 0.182 |
| P(fox present, coyote absent) | 0.354 |
| P(coyote present, fox present) | 0.399 |
| P(coyote present, fox absent) | 0.620 |
| interaction (log odds ratio) | -0.90 (odds ratio 0.41) |
| interaction slope per SD log housing | +0.60 |
| coyote first-order: intercept / forest / log housing | 0.49 / -0.24 / +0.17 (approximate, see above) |
| gray fox first-order: intercept / forest / log housing | -0.60 / +0.70 / +0.30 |
| weekly detection, coyote / gray fox (off trail, 6.5 days) | 0.250 / 0.294 |
| detection: days_active (per SD = 1.23 days) | +0.31 both species |
| detection: on trail vs off (logit) | coyote +0.80, gray fox -0.20 |

Realised in this one simulated landscape: coyotes truly occupied 133 cameras,
gray foxes 67, both 31. Naive occupancy (detected at least once) was 0.51 for
coyote and 0.24 for gray fox.
