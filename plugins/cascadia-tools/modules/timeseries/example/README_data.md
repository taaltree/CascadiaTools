# Example data: `herd_counts.csv` (SIMULATED)

**These data are simulated.** They were made by the module's own simulator so
that the true values are known and the analysis can be checked against them.
No real herd, agency or survey is described.

## The story

A wildlife biologist manages a hypothetical mule deer herd on a
mountain winter range. Each winter since 1994 the herd has been counted from a
fixed-wing aircraft on its winter range. Aerial counts miss some deer and
sometimes count groups twice, so a single count is only roughly right. In
three winters (2001, 2009, 2016) the flight was cancelled because of weather.
A nearby snow course records the mean winter snow depth every year, including
the years without a flight.

The biologist wants to know:

1. Is the herd declining, and by how much per year?
2. How many deer are really there, once counting error is removed?
3. Does the herd regulate itself (density dependence)?
4. Do deep-snow winters reduce herd growth?

## Design

| dimension | value |
|---|---|
| populations (herds) | 1 |
| years | 31 (1994-2024) |
| winters counted | 28 |
| winters not flown | 3 (2001, 2009, 2016; interior years) |
| covariates | 1 (mean winter snow depth, cm; observed every year) |
| counts | 316 to 1,125 deer (median 532) |

## Columns

| column | meaning |
|---|---|
| `year` | Winter of the count (the count is flown late in the winter that ends in this year). |
| `survey` | Survey method. Always "aerial winter count" here. |
| `animals_counted` | Number of deer counted. Blank = not flown. |
| `snow_depth_cm` | Mean winter snow depth (cm) for the winter ending in `year`. In the model it acts on growth from the previous winter's count to this one. |
| `notes` | Free text. "not flown (weather)" marks the missing years. |

## How the data were generated

Script: `dev/timeseries/make_example.R` (part 1). It calls
`module$simulate()` with the module's default parameter values and the
design `{model: bayes, n_years: 31, first_year: 1994, covariate_names:
snow_depth_cm}`, seed **20**. It then blanks the counts for 2001, 2009 and
2016, and rounds snow depth to whole cm (minimum 5 cm).
`dev/timeseries/si/check_example_repro.R` re-runs this and confirms the CSV
is reproduced exactly. The true abundance in every year is saved in
`dev/timeseries/si/example_true_N.csv`, and the truth vector in
`dev/timeseries/example_truth.rds`.

The generating model (log scale, x = log true abundance):

    x[t] = x[t-1] + b * (x[t-1] - log K) + beta * z[t] + e[t],   e ~ Normal(0, sd_process)
    log(count[t]) = x[t] + f[t],                                f ~ Normal(0, sd_observation)

z[t] is snow depth standardised to mean 0 and SD 1 (raw snow depth was drawn
as Normal(60, 20) cm; in this dataset mean 54.1 cm, SD 20.0 cm).

## TRUE values (simulated)

| parameter | true value | meaning |
|---|---|---|
| `density_dependence_b` | -0.30 | Gompertz density dependence (0 = none, -1 = full return to equilibrium each year) |
| `equilibrium_N` | 500 | Herd size at which growth averages zero (at average snow) |
| `sd_process` | 0.10 | Real year-to-year variation in log growth (about 10%) |
| `sd_observation` | 0.15 | Counting error on the log scale (a count is typically about 16% off) |
| `growth_beta_snow_depth_cm` | -0.15 | Change in log growth per 1 SD (20 cm) more snow |
| `growth_effect_pct_snow_depth_cm` | -13.9% | The same as a % change in the annual growth multiplier |
| `N_final` | 643 | True herd size in 2024 |
| `realised_trend_pct_per_year` | +0.10% per year | Average annual change of the true herd, 1994 (625) to 2024 (643) |

So the true answer is: the herd is **not** declining. It fluctuates around an
equilibrium of about 500 deer, deep-snow winters reduce growth, and the 2024
herd is about 643 deer.
