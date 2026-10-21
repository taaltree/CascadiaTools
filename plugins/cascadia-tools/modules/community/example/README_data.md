# Example data: breeding birds in western Oregon Cascades forests (SIMULATED)

**These data are simulated.** No real birds were counted. They were made with
this module's own simulator so that the true answers are known. That lets a
reader check whether the analysis gets the right answer.

## The story

A field crew surveyed breeding birds at 60 point-count stations in conifer
forest on the west slope of the Oregon Cascades in late May 2025. Each
station was counted once: 10 minutes, all birds seen or heard within 50 m.
Stations sit in three stand types:

- old growth (22 stations)
- thinned second growth (20 stations)
- young plantations (18 stations)

At each station the crew recorded elevation and canopy cover. The questions:

1. How much of the difference in bird communities among stations goes with
   elevation, and how much with canopy cover?
2. Do the three stand types hold different bird communities, or are they
   just more or less variable?
3. Do stations differ because species replace each other (turnover), or
   because poor stations hold subsets of rich ones (nestedness)?
4. Which species characterise each stand type?

## Design

| Dimension | Value |
|---|---|
| Stations (sites) | 60 |
| Species in the table | 30 (each recorded at 5 to 47 stations) |
| Stand types | 3: old growth 22, thinned 20, young plantation 18 |
| Visits | 1 count per station (10 min, 50 m radius) |
| Gradients | elevation 407-1582 m; canopy cover 14-88% |
| Correlation of elevation and canopy cover | 0.33 |
| Birds counted | 1,774 in total; 3-17 species per station (median 10) |

## Files

### `cascades_bird_counts.csv` (wide format: one row per station)

| Column | Meaning |
|---|---|
| `station` | Station id, `CAS-01` to `CAS-60` |
| `PSFL` ... `HETH` | Number of birds of that species counted (4-letter AOU codes; the full names are in `config.yaml`, `figures: labels`) |

The 30 species: Pacific-slope Flycatcher (PSFL), Pacific Wren (PAWR), Hermit
Warbler (HEWA), Chestnut-backed Chickadee (CBCH), Brown Creeper (BRCR),
Golden-crowned Kinglet (GCKI), Varied Thrush (VATH), Swainson's Thrush
(SWTH), Wilson's Warbler (WIWA), MacGillivray's Warbler (MGWA),
Orange-crowned Warbler (OCWA), Spotted Towhee (SPTO), Dark-eyed Junco (DEJU),
American Robin (AMRO), Red-breasted Nuthatch (RBNU), Hairy Woodpecker (HAWO),
Pileated Woodpecker (PIWO), Western Tanager (WETA), Black-headed Grosbeak
(BHGR), Evening Grosbeak (EVGR), Hutton's Vireo (HUVI), Cassin's Vireo
(CAVI), Olive-sided Flycatcher (OSFL), Western Wood-Pewee (WEWP), Rufous
Hummingbird (RUHU), Song Sparrow (SOSP), Black-throated Gray Warbler (BTYW),
Townsend's Warbler (TOWA), Steller's Jay (STJA), Hermit Thrush (HETH).

Because the counts are simulated, the species' responses are random draws.
They do not reflect these species' real habitat use. The names make the
example readable; do not read biology into which species came out where.

### `cascades_stations.csv` (one row per station)

| Column | Meaning |
|---|---|
| `station` | Station id, matches the counts file |
| `stand_type` | `old_growth`, `thinned` or `young_plantation` |
| `elevation_m` | Elevation (m), rounded to 1 m |
| `canopy_cover_pct` | Canopy cover (%), rounded to 1% |
| `survey_date` | Date of the count (random within 20 May - 19 June 2025; not used) |

## How the data were made

Script: `dev/community/make_example.R` (seed 2026). It calls
`module$simulate()` with the module defaults plus the changes marked below.
`dev/community/si/check_example.R` re-runs the simulation and confirms that
it reproduces both CSVs and the truth file exactly.

Simulation settings (TRUE, simulated values):

| Setting | Value | Meaning |
|---|---|---|
| `n_species` | 30 | species pool |
| `abund_max` | 4 (lognormal SD 0.8) | median peak expected count per station |
| `tol1` | 0.2 | niche width on elevation, as a fraction of the 400-1600 m range (narrow: strong gradient) |
| `tol2` | 0.4 | niche width on canopy cover, as a fraction of the 10-95% range (wider: weaker gradient) |
| `group_effect_sd` | **0.9** (default 0.5) | SD of species-specific stand-type effects on log abundance |
| `site_noise_sd` | 0.5 | extra station-to-station variation (lognormal SD) |
| `dispersion_mult` | **1.4** (default 1) | young plantations 1.4 times as variable among stations |
| `env_cor` | 0.4 | correlation of the two gradient positions |
| `effort` | 1 | one count per station; Poisson counts, no detection model |

## True values of the estimands (SIMULATED)

Each is the population value of the statistic for this species pool and
design: the value the analysis would give with infinitely many stations.
They were computed from 40,000 extra stations and 60,000 station pairs per
stand-type combination, independent of the 60 analysed stations
(`example_truth.rds`; also `results/community/si/truth.csv`).

| Estimand | Truth |
|---|---|
| Share of composition variation explained by elevation (adjusted R2) | 0.252 |
| ... by canopy cover | 0.106 |
| ... by elevation after canopy cover (unique) | 0.212 |
| ... by canopy cover after elevation (unique) | 0.067 |
| ... shared by the two (correlated) gradients | 0.040 |
| ... by both together | 0.319 |
| Stand type, PERMANOVA adjusted R2 (Bray-Curtis) | 0.050 |
| Mean pairwise Sorensen dissimilarity | 0.619 |
| Mean pairwise turnover (Simpson) | 0.517 |
| Mean pairwise nestedness | 0.102 |
| Turnover share (turnover / Sorensen) | 0.835 |

The young plantations were simulated with 1.4 times the station-to-station
noise. So there is a real, small difference in dispersion. It has no
population "truth" value in the gate (see VALIDATION.md).
