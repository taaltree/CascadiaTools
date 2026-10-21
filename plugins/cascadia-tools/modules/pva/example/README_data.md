# Example data for the `pva` module

This folder has four example datasets. The worked example used in the
manuscript Supporting Information (Appendix S1) and on the website is the
**shorebird breeding-pair survey** (`shorebird_breeding_pairs.csv`, fitted with
`config_state_space.yaml`), because it is simulated and so its true values are
known. The other three are described at the end.

## Worked example: coastal shorebird breeding-pair survey (SIMULATED)

### The story

A coastal shorebird nests on a chain of beaches. Each spring from 1994 to 2023
a field crew walked the beaches and counted breeding pairs. The counts are a
survey, not a census: some pairs are missed and some are double-counted, so a
single count is off by roughly 15-20%. In five years (1997, 2003, 2004, 2011,
2020) there was no survey (funding gaps, a storm closure, the pandemic). The
counts fell from about 400-500 pairs in the 1990s to about 80 pairs now.

The managers ask: **is the population really declining, and how likely is it
to fall to 40 breeding pairs (the quasi-extinction threshold they use for a
"critical" listing) within the next 10 or 25 years?**

The data are simulated. No real species or site is described.

### Design

| dimension | value |
|---|---|
| years spanned | 30 (1994-2023) |
| years surveyed | 25 |
| years not surveyed | 5 (1997, 2003, 2004, 2011, 2020) |
| counts per year | 1 |
| quasi-extinction threshold | 40 breeding pairs |
| risk horizons | 10 and 25 years after 2023 |

### Columns of `shorebird_breeding_pairs.csv`

| column | meaning |
|---|---|
| `survey_year` | year of the spring survey |
| `breeding_pairs` | number of breeding pairs counted; **blank = not surveyed** |
| `observer_notes` | free text; "no survey" in the skipped years (ignored by the analysis) |

### How the data were generated

`dev/pva/make_examples.R`, block 2 (`set.seed(2024)`). True log abundance
follows a random walk with drift starting at 420 pairs in 1994; each count is
the true abundance times lognormal counting error; counts are rounded; the five
years were then blanked. `dev/pva/si_truth.R` checks that
`module$simulate(seed = 2024)` with the values below reproduces every count
exactly, and writes the full truth (including the true abundance in every year)
to `results/pva/si/truth.json`.

### TRUE values (simulated; known only because the data are simulated)

| parameter | true value | meaning |
|---|---|---|
| `growth_rate_mu` | -0.04 | mean annual log growth rate |
| `lambda_stochastic` | 0.961 | exp(mu): about -3.9% per year |
| `process_variance_sigma2` | 0.012 | real year-to-year variance of log growth (SD 0.11) |
| `obs_variance_tau2` | 0.03 | counting-error variance on the log scale (CV about 17%) |
| N in 1994 | 420 | starting true abundance |
| `N_current` | 66.0 | true number of breeding pairs in 2023 |
| `p_qext_10yr` | 0.446 | true probability of N <= 40 at an annual census within 10 years (given mu, sigma2 and the true 2023 abundance; exact first-passage recursion) |
| `p_qext_25yr` | 0.882 | the same within 25 years |

Note: the 2023 count was 79, but the true abundance was 66. This is counting
error, which is why the state-space model is used.

## Other examples in this folder

* `yellowstone_grizzly_females.csv` + `config.yaml`: **real data**. Estimated
  adult female grizzly bears in Yellowstone, 1959-1997 (Morris & Doak 2002
  Table 3.1, via `popbio::grizzly`; Eberhardt et al. 1986). Diffusion PVA,
  threshold 20 females. No truth (real data); used as a cross-check against the
  book's estimates.
* `vital_rates.csv` + `config_matrix.yaml`: illustrative vital-rate table for a
  3-stage long-lived bird (values chosen by hand, as if from a CJS/IPM
  analysis). Matrix PVA with three management scenarios.
* `annual_monitoring.csv` + `config_matrix_monitoring.yaml`: raw annual
  monitoring (at risk / survivors, breeders / young), simulated with
  `module$simulate()` (seed 4242) at the module's default 2-stage design.
