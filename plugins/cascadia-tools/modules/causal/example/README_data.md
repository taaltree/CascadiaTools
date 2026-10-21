# Example data: riparian canopy and juvenile trout (SIMULATED)

**These data are simulated.** They were generated with known true values so a
reader can check that the `causal` module gets the right answer. No real
stream was sampled.

## The (simulated) study

A fisheries biologist surveyed 140 stream reaches across seven coastal
streams in summer 2025 (1 July to 15 August). Each reach is 100 m long. She
made one electrofishing pass per reach and counted juvenile trout. At each
reach she also measured riparian canopy cover, elevation, the share of
agricultural land in the catchment, fine sediment on the bed, wetted width,
August water temperature and periphyton chlorophyll-a.

Her question: **does more riparian canopy cover cause more juvenile trout?**
Nobody planted trees at random, so this is an observational study. The
canopy-trout relationship is confounded by elevation, land use and stream
size. The `causal` module uses a causal diagram (DAG) to decide what to adjust
for.

## Design

| dimension | value |
|---|---|
| reaches (rows) | 140 |
| streams | 7 (labels only; not used in the model) |
| surveys per reach | 1 electrofishing pass |
| exposure | `canopy_pct` (continuous) |
| outcome | `trout_per_reach` (count, negative binomial) |
| contrast | +10 percentage points of canopy |
| estimands | total effect and controlled direct effect (holding water temperature fixed) |
| DAG | 8 nodes, 12 arrows (see `config.yaml`) |

## Files

- `stream_trout_reaches.csv`: the data (140 rows).
- `config.yaml`: the analysis configuration a researcher would use.
- `recovery.yaml`: the planted-truth check at this design.
- `config_panel.yaml`, `snail_panel.csv`: a second, unrelated panel example (not covered here).

## Columns of `stream_trout_reaches.csv`

| column | meaning | units | role in the DAG |
|---|---|---|---|
| `reach_id` | reach label R001-R140 | - | identifier |
| `stream` | stream name (Elk, Mill, Beaver, Rock, Fish, Cedar, Alder Cr) | - | label only |
| `survey_date` | date of the electrofishing pass | YYYY-MM-DD | label only |
| `canopy_pct` | riparian canopy cover | % | **exposure** |
| `trout_per_reach` | juvenile trout caught per 100 m reach | count | **outcome** |
| `elevation_m` | reach elevation | m | confounder (-> canopy, -> water temp) |
| `agriculture_pct` | agricultural land in the catchment | % | confounder (-> canopy, -> sediment) |
| `fine_sediment_pct` | fine sediment cover on the bed | % | on a back-door path (-> trout) |
| `wetted_width_m` | wetted channel width | m | confounder (-> canopy, -> trout) |
| `water_temp_c` | August mean water temperature | deg C | **mediator** (canopy -> temp -> trout) |
| `periphyton_chla` | periphyton chlorophyll-a | mg/m2 | **collider** (canopy -> chla <- trout) |

## TRUE values (simulated)

Generator: `dev/causal/make_examples.R` (module simulator `causal_simulate`,
seed 2024, which was the first seed tried in 2024-2100 that gave
agriculture >= 0 %, canopy <= 100 % and width > 0.5 m). Running
`dev/causal/si_check.R` with the settings in `recovery.yaml` and seed 2024
reproduces the CSV exactly (max difference 0). The truth is saved in
`dev/causal/example_truth.rds`.

Structural model: standardised path coefficients (per SD of the parent). Trout
counts are negative binomial with a log link, intercept log(5) (about 5 trout
per reach at average conditions) and size theta = 3. Every other node is
linear-Gaussian, rescaled to natural units:

| arrow | coefficient |
|---|---|
| elevation -> canopy | 0.4 |
| elevation -> water temp | -0.6 |
| agriculture -> canopy | -0.5 |
| agriculture -> fine sediment | 0.6 |
| fine sediment -> trout | -0.4 |
| wetted width -> canopy | -0.3 |
| wetted width -> trout | 0.3 |
| canopy -> water temp | -0.4 |
| water temp -> trout | -0.4 |
| canopy -> trout (direct) | 0.2 |
| canopy -> periphyton | -0.4 |
| trout -> periphyton | 0.5 |

Natural-unit scales (mean, SD): elevation 620, 240 m; agriculture 35, 9 %;
fine sediment 18, 5 %; wetted width 7.5, 2 m; canopy 58, 14 %; water temperature
14.5, 2.2 deg C; periphyton 42, 11 mg/m2.

**True causal effects of +10 percentage points of canopy** (found by intervening
on the structural model, 400,000 draws):

| estimand | truth |
|---|---|
| total effect, log scale per % canopy (`beta_canopy_pct`) | 0.02571 |
| total average causal effect (`ate_canopy_pct`) | +1.954 trout per reach |
| total rate ratio (`rate_ratio_canopy_pct`) | 1.293 (+29 %) |
| direct effect, log scale per % (`direct_beta_canopy_pct`) | 0.01429 |
| controlled direct effect (`direct_ate_canopy_pct`) | +1.023 trout per reach |
| direct rate ratio (`direct_rate_ratio_canopy_pct`) | 1.154 (+15 %) |

## Quick facts about this draw

Canopy 27-97 % (mean 61.4). Trout 0-68 per reach (mean 6.3, variance 50.4;
6.4 % zeros). The one reach with 68 trout makes this draw heavy-tailed. It is
why DHARMa flags dispersion (p = 0.016) on this dataset, even though the model
family is correct.
