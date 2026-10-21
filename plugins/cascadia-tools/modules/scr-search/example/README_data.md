# Example data: fisher scat-detection-dog survey (SIMULATED)

**These data are simulated.** They were made by `dev/scr-search/make_example.R`
(`set.seed(2025)`, `module$simulate(..., seed = 2025)`), so the true values are known.
The files are written in the formats a real field project would have.

## The (simulated) study

- **Species:** fisher (*Pekania pennanti*), a forest-dependent mesocarnivore.
- **Place:** a 16 x 16 km block of the western Oregon Cascades (UTM zone 10N), summer 2025.
- **Method:** four dog-handler teams searched for fisher scats with scat-detection dogs.
  Every search route was logged by a handheld GPS (one fix per minute). Each scat found was
  collected, its position recorded, and it was genotyped at 12 microsatellite loci to
  identify the individual fisher.
- **Effort:** 3 monthly sessions (June, July, August). 12 routes per session, 36 routes in
  all, 4-11 km each, 247 km searched. About 40% of the routes follow forest roads; the rest
  go cross-country. The roads run moderately through more-forested ground.
- **Teams:** two experienced dog teams (Pepin & Sara, Jack & Tom; 15 routes, 100 km) and two
  newer teams (Odin & Kim, Rue & Lena; 21 routes, 146 km).
- **Results in the field:** 116 scats collected. 35 (30%) failed genotyping. The other 81 came
  from 22 different fishers (1 to 12 scats each; 20 of 22 had 2 or more scats; 14 were
  found in more than one session).
- **Questions:** How many fishers live in the area, and at what density? Is density higher
  where there is more forest? Do experienced dog teams find more scats per km?

## Design dimensions

| dimension | value |
|---|---|
| searched region | 16 x 16 km (256 km2) |
| buffer around it (habitat grid) | 6 km; state space 28 x 28 km = 784 km2 |
| habitat grid | 1,600 pixels of 700 m |
| sessions (occasions) | 3 |
| search routes | 36 (12 per session) |
| total km searched | 246.7 |
| detector cells used by the model | 1 km x 1 km |
| scats collected / genotyped / failed | 116 / 81 / 35 |
| fishers identified | 22 |

## Files and columns

`search_tracks_gps.csv`: GPS log, one row per fix (4,715 rows).

| column | meaning |
|---|---|
| track_id | route ID (T001-T036) |
| gps_time | time stamp of the fix (UTC, ISO 8601); used to order fixes |
| utm_e, utm_n | position, UTM zone 10N (m) |

`track_field_sheet.csv`: one row per route (36 rows).

| column | meaning |
|---|---|
| track_id | route ID, matches the GPS log |
| date | date searched |
| session | survey session 1-3 (June, July, August) |
| dog_team | dog & handler names |
| experienced_team | 1 = experienced team, 0 = newer team |
| km_walked | route length from the GPS (km) |

`fisher_scats.csv`: one row per scat collected (116 rows).

| column | meaning |
|---|---|
| scat_id | lab sample ID |
| track_id | route on which the scat was found |
| date | date collected |
| utm_e, utm_n | where the scat was found (m) |
| fisher_id | individual from the genotype; blank = failed genotype (dropped by the module) |
| genotype_qc | "pass" or "failed (<6 loci)" |

`habitat_grid.csv`: habitat pixels covering the searched area plus 6 km (1,600 rows).

| column | meaning |
|---|---|
| utm_e, utm_n | pixel centre (m), 700 m spacing |
| forest_pct | forest cover in the pixel (%), range 32.6-85.5, mean 55.0, SD 12.0 |

## TRUE values (simulated)

| quantity | true value |
|---|---|
| expected density | 0.10 fishers per km2 |
| realised number of activity centres in the 784 km2 state space | 63 (realised density 0.080 per km2) |
| sigma (spatial scale of space use) | 1,500 m |
| 95% home-range area | 42.4 km2 |
| lambda0: scats found per km searched at an activity centre, newer team | 0.20 |
| samples per straight route through the centre (lambda0 x sqrt(2 pi) x sigma) | 0.752 |
| forest effect on density (log ratio per 1 SD = 12 % forest) | 0.5 (rate ratio 1.65) |
| experienced-team effect on scats per km (log ratio) | 0.6 (rate ratio 1.82) |

The scat-deposition truth refers to genotyped scats (failed scats were added on top, at
random positions on the routes). Scats were deposited continuously along the walked GPS
paths, so the model's 1 km cells are an approximation that is tested, not assumed. The
expected density (0.10) is the model's `density_per_km2` target; the realised number (63)
is the target for `abundance` in this one data set.
