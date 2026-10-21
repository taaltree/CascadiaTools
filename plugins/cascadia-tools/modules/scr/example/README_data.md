# Example data: bobcat camera-trap survey (SIMULATED)

**These data are simulated.** No real bobcats, cameras or places are
involved. They were generated with known ("true") parameter values so that
readers can see how close the analysis gets to the truth. Generator:
`dev/scr/make_example.R`, random seed `20250601`.

## Species and system

Bobcat (*Lynx rufus*) in a managed conifer forest. Bobcats are solitary and
territorial. Individuals can be told apart from photos by their coat spot
and bar patterns, so each photo can be assigned to a known cat. Density is
expected to be higher where the forest canopy is denser. A reservoir in the
south-east of the area is not bobcat habitat.

## The (simulated) study

A biologist set 48 motion-triggered cameras on a grid about 1.6 km apart
(placement shifted up to 250 m to fit access). The cameras ran for 8 weeks,
2 June to 27 July 2025. Each week is one sampling occasion. Eighteen
cameras were on logging roads or trails, where bobcats travel and are
photographed more often; 30 were on off-trail game paths. Some cameras
failed part-way (dead batteries, theft) and two were set a week late, so
each camera has a record of the weeks it was working.

Photos were sorted, and each bobcat photo was assigned to an individual
(`LYRU-001`, `LYRU-002`, ...). The biologist also has a habitat grid: 400 m
pixels covering the cameras plus 5 km around them, with percent canopy
cover, elevation and a habitat flag (0 = reservoir).

**Question:** how many bobcats live in the area (density per km2), and is
bobcat density higher where canopy cover is higher?

## Design dimensions

| item | value |
|---|---|
| cameras | 48 (8 x 6 grid, ~1,600 m apart); 18 on trails |
| occasions | 8 weeks (2 June - 27 July 2025) |
| camera-weeks working | 354 of 384 (8 cameras lost 1-6 weeks) |
| bobcat detections (photo-detections, one per cat x camera x week) | 114 |
| individuals photographed | 33 |
| individuals photographed at more than one camera | 25 (53 spatial recaptures) |
| cameras that photographed a bobcat | 40 of 48 |
| habitat grid | 2,255 pixels of 400 m x 400 m; 2,189 habitat pixels = 350.2 km2 |
| buffer around the cameras | 5,000 m |

## Files and columns

`cameras.csv` (one row per camera)

| column | meaning |
|---|---|
| `camera_id` | camera code, BC01-BC48 |
| `utm_e`, `utm_n` | UTM easting and northing (m) |
| `on_trail` | 1 = on a logging road or trail, 0 = off-trail game path |
| `weeks_active` | 8 characters, one per week: 1 = camera working, 0 = not working |
| `habitat_note` | free-text note (not used in the analysis) |

`bobcat_detections.csv` (one row per photo-detection)

| column | meaning |
|---|---|
| `cat_id` | individual bobcat, identified from coat pattern |
| `week` | occasion, 1-8 |
| `camera_id` | camera that took the photo (matches `cameras.csv`) |
| `date` | date of the photo |
| `photo` | photo file name (not used in the analysis) |

`habitat_grid.csv` (one row per 400 m pixel centre)

| column | meaning |
|---|---|
| `x`, `y` | UTM easting and northing of the pixel centre (m) |
| `canopy_cover_pct` | forest canopy cover (%) |
| `elevation_m` | elevation (m); has NO true effect on density |
| `habitat` | 1 = bobcat habitat, 0 = reservoir (not habitat) |

`config.yaml` tells ecotools which file and column is which.
`recovery.yaml` holds the same design and true values for the
planted-truth check; `recovery_matched_n.yaml` repeats it with density
lowered so that the number of detected animals (~32) matches this data set.

## TRUE parameter values (simulated)

| parameter | true value |
|---|---|
| expected density | 0.25 bobcats per km2 of habitat (at average canopy cover over the grid) |
| realised abundance | 82 bobcats in the 350.2 km2 of habitat (realised density 0.234 per km2) |
| sigma (spatial scale of movement) | 1,300 m (95% home range = pi x (2.448 x 1300)^2 = 31.8 km2) |
| p0 (weekly detection prob. at the activity centre, off-trail camera) | 0.12 |
| on-trail effect on detection (logit scale) | +1.10 (p0 on trail = 0.29) |
| canopy effect on density (log scale, per SD of canopy cover) | +0.45 (SD = 14.7 percentage points; rate ratio 1.57 per SD) |
| elevation effect | 0 (not used to simulate) |

The realised number of bobcats (82) differs from the expected number
(0.25 x 350.2 = 87.6) because the number of animals is itself random
(Poisson). Abundance estimates should be compared with 82; density
estimates with the realised density 0.234 and the expected 0.25.
