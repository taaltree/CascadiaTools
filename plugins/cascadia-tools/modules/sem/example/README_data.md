# Example data: post-fire forest plots (module `sem`)

**These data are SIMULATED.** No real field survey produced them. They were made
with known "true" values so that every estimate can be checked against the
right answer. Generator: `dev/sem/make_example.R` with `set.seed(2024)`.
Rerunning that script reproduces `burned_forest_plots.csv` byte for byte
(checked 2026-10-20 with `dev/sem/si/make_example_check.R`).

## The story

A forest ecologist surveyed a mixed-conifer forest one to two years after a
large wildfire. The fire burned unevenly. Some stands lost almost all their
canopy and others were barely scorched. She wants to know how fire severity
affects breeding songbirds. Fire may harm birds directly (fewer nest trees,
less food in the canopy). Fire may also help them indirectly: canopy loss
lets light reach the ground, shrubs grow back, and many shrub-nesting birds
increase. Wetter soils also grow more shrubs.

Her hypothesised diagram:

```
fire severity --> canopy cover --> shrub cover --> bird detections
      |                               ^                   ^
      |                    soil moisture                  |
      +---------------------------------------------------+   (direct path)
```

She sampled 120 plots in 12 burn units (10 plots per unit). At each plot she
recorded satellite fire severity (dNBR), canopy cover, soil moisture and shrub
cover, and she ran one 10-minute bird point count between 20 May and 14 June 2025.
Plots in the same burn unit share local conditions, so burn unit is a random
effect in the bird model.

## Design

| dimension | value |
|---|---|
| sampling units (plots) | 120 |
| groups (burn units) | 12, with 10 plots each |
| bird point counts | 1 per plot |
| variables in the diagram | 5 (2 drivers, 2 mediators, 1 count outcome) |
| survey window | 2025-05-20 to 2025-06-14 |

## Columns (`burned_forest_plots.csv`, one row per plot)

| column | meaning | units / values | range in the file |
|---|---|---|---|
| `plot_id` | plot label | P001-P120 | |
| `burn_unit` | burn unit the plot is in (random intercept) | BU01-BU12 | 10 plots each |
| `survey_date` | date of the bird point count | YYYY-MM-DD | 2025-05-20 to 2025-06-14 (not used in the model) |
| `fire_severity_dnbr` | differenced Normalized Burn Ratio from satellite images; higher = more severe | dNBR units | 32-1008 (mean 497, SD 212) |
| `canopy_cover_pct` | tree canopy cover | % | 10.5-89.0 (mean 52.2, SD 15.1) |
| `soil_moisture_pct` | volumetric soil moisture | % | 9.3-36.7 (mean 22.1, SD 5.0) |
| `shrub_cover_pct` | shrub cover | % | 22.6-59.5 (mean 40.7, SD 8.2) |
| `bird_detections` | birds detected in one 10-min point count (all species) | count | 1-17 (mean 5.7; variance/mean 1.29; no zeros) |

No values are missing.

## TRUE values (simulated)

The generating model in `dev/sem/make_example.R`:

| quantity | true value | meaning |
|---|---|---|
| burn-unit mean dNBR | Uniform(150, 750) | each burn unit has its own typical severity |
| dNBR within a unit | Normal(unit mean, SD 150), limited to 0-1100, rounded | |
| canopy cover | 78 - **0.055** x dNBR + Normal(0, SD 9) | -5.5 % canopy per 100 dNBR |
| soil moisture | Normal(22, SD 5) | independent of fire |
| shrub cover | 45 - **0.35** x canopy + **0.6** x moisture + Normal(0, SD 7) | -3.5 % shrub per 10 % canopy; +3 % per 5 % moisture |
| bird detections | Poisson(exp(1.0 + **0.025** x shrub - **0.0006** x dNBR + u)), u ~ Normal(0, SD **0.25**) per burn unit | log link |

Derived true effects (products of the paths, on the scale the module reports):

| effect | true value | in plain words |
|---|---|---|
| indirect fire -> shrub (via canopy) | (-0.055)(-0.35) = **+0.01925** % per dNBR | +1.9 % shrub per 100 dNBR |
| direct fire -> birds | **-0.0006** per dNBR (log scale) | x0.94 birds per 100 dNBR |
| indirect fire -> birds (via canopy, shrub) | (-0.055)(-0.35)(0.025) = **+0.00048** per dNBR | x1.05 birds per 100 dNBR |
| total fire -> birds | **-0.00012** per dNBR | x0.99 birds per 100 dNBR: the two routes nearly cancel |
| indirect moisture -> birds (via shrub) | (0.6)(0.025) = **+0.015** per % moisture | x1.08 birds per 5 % moisture |

There is **no** direct soil moisture -> birds path and **no** direct fire ->
shrub path in the truth. The diagram fitted in `config.yaml` is the true diagram.

Standardised true values (population values from a 200,000-plot simulation of
the recovery design, `example/recovery.yaml`): fire -> canopy -0.81,
canopy -> shrub -0.58, moisture -> shrub 0.32, shrub -> birds 0.45,
fire -> birds (direct) -0.27; indirect fire -> birds 0.21, total fire -> birds
-0.05; indirect moisture -> birds 0.15; indirect fire -> shrub 0.47.
These depend on the spread of each variable, so they apply only to a
population like this one.

## Files in this folder

- `burned_forest_plots.csv`: the data above (piecewise SEM example).
- `config.yaml`: the analysis configuration a researcher would use.
- `recovery.yaml`: planted-truth design matching this dataset.
- `grassland_plots.csv`, `config_lavaan.yaml`: a separate latent-variable (lavaan) example, not part of this walkthrough.
- `WALKTHROUGH.md`: the worked example (Appendix S1).
