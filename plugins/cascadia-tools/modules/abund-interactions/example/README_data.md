# Example data: pine marten and red squirrel camera-trap survey (SIMULATED)

**These data are simulated.** No real animals were recorded. The study is
invented so that readers can see the whole CascadiaTools workflow on data
whose true answer is known. The generating script is
`dev/abund-interactions/make_example.R` (seed `set.seed(2024)`). It was
written separately from the module's own `simulate()`, so a shared coding
error would not hide itself.

## The story

In Ireland and Britain, the native red squirrel declined after the grey
squirrel was introduced. Where the pine marten (a native predator) has
recovered, grey squirrels have declined and red squirrels have come back.
Martens also eat red squirrels, though. Field studies suggest red squirrels
cope with martens better in native broadleaf woodland, where they share a
long evolutionary history, than in conifer plantations (Twining et al. 2022).

A (fictional) team ran one camera trap at each of 260 woodland sites across
Northern Ireland in March-April 2024. Each camera ran for up to 28 days.
They want to know two things:

1. Is red squirrel occupancy lower at sites used by MORE pine martens?
2. Does that association depend on how much of the woodland is native
   broadleaf?

## Design

| dimension | value |
|---|---|
| sites | 260 (one camera each) |
| deployment length | up to 28 camera-days; 6,939 camera-days in total |
| occasions | 4 occasions of 7 camera-days (set in `config.yaml`) |
| effort loss | 8.9% of site-occasions have fewer than 7 camera-days; 3.0% have none (short deployments, 19 sites with a card-swap gap) |
| survey window | 2024-03-01 to 2024-04-17 |
| species in records | Red squirrel 846, Pine marten 649, Badger 397, Red fox 282, Grey squirrel 203 photos (2,377 rows) |
| naive detections | Pine marten at 142 sites (54.6%, 539 detection-days); red squirrel at 130 sites (50.0%, 546 detection-days) |

Badger, red fox and grey squirrel records are bycatch. The module ignores them.
They are there because real camera exports contain other species.

## Files and columns

### `deployments.csv` (279 rows: one per camera deployment)
| column | meaning |
|---|---|
| `site` | site ID (`NI001`-`NI260`); matches `sites.csv` |
| `camera_id` | camera ID (one camera per site) |
| `start_date` | first full day the camera worked (YYYY-MM-DD) |
| `end_date` | last full day the camera worked (YYYY-MM-DD) |

A site with two rows had its memory card swapped. The days between the rows
are lost effort and are counted as zero camera-days.

### `records.csv` (2,377 rows: one per photo record)
| column | meaning |
|---|---|
| `site` | site ID |
| `date` | date-time of the photo (YYYY-MM-DD HH:MM:SS) |
| `species` | common name; the module uses `Pine marten` and `Red squirrel` |
| `count` | animals in the photo (not used; the model uses detection-days) |

### `sites.csv` (260 rows: one per site)
| column | meaning |
|---|---|
| `site` | site ID |
| `easting`, `northing` | Irish Grid-style coordinates in metres (used only for the map) |
| `woodland_pct` | woodland cover within 1 km of the camera (%), range 5-95, mean 45.8, SD 17.4 |
| `broadleaf_pct` | share of that woodland that is native broadleaf (%), range 0-100, mean 38.7, SD 20.8 |
| `dist_settlement_km` | distance to the nearest settlement (km); recorded but not in the model (a distractor, as in real data) |

## TRUE parameter values (simulated; for checking only)

Covariates were standardised over the 260 sites (z = (x - mean)/SD), the same
way the module does it.

| parameter | symbol | true value | meaning |
|---|---|---|---|
| `dominant_abundance` | exp(lambda0) | 1.0 | pine martens using a site at average woodland cover |
| `lambda_beta_woodland` | lambda1 | 0.5 | log-abundance slope per SD of woodland cover |
| `dominant_detection_per_individual` | r | 0.08 | daily chance that the camera records one given marten |
| `subordinate_occupancy_no_dominant` | plogis(psi0) | 0.6 | red squirrel occupancy with no martens, average habitat |
| `psi_beta_woodland` | psi1 | 0.5 | logit occupancy slope per SD of woodland cover |
| `psi_beta_broadleaf` | psi2 | 0.3 | logit occupancy slope per SD of broadleaf share |
| `psi_beta_dominant_abundance` | gamma0 | -0.6 | change in red squirrel logit occupancy per marten, at average broadleaf share |
| `psi_beta_dominant_abundance_x_broadleaf` | gamma1 | 0.6 | change in gamma per SD of broadleaf share |
| `subordinate_detection` | p | 0.15 | daily detection of red squirrel at an occupied site |

What these imply:
- The interaction per marten is -0.6 + 0.6 z(broadleaf). That is **-1.34** at
  the 10th percentile of broadleaf share (12.9%) and **+0.16** at the 90th
  percentile (65.2%). It crosses zero at about 59% broadleaf.
- Realised in this one data set: 274 marten-site associations in total; 59.2%
  of sites used by at least one marten; 50.4% of sites occupied by red squirrels
  (see `TRUTH.md`).

## Model used to generate the data

    N_i  ~ Poisson(1.0 * exp(0.5 z(woodland_i)))                      marten abundance
    marten detected on day d with prob 1 - (1 - 0.08)^N_i               Royle-Nichols
    logit psi_i = logit(0.6) + 0.5 z(woodland_i) + 0.3 z(broadleaf_i)
                  + N_i * (-0.6 + 0.6 z(broadleaf_i))                   red squirrel occupancy
    z_i ~ Bernoulli(psi_i); squirrel detected on day d with prob 0.15 z_i

## Reference
Twining JP, Sutherland C, Reid N, Tosh DG (2022) Habitat mediates coevolved but
not novel species interactions. Proc R Soc B 289:20212338.
