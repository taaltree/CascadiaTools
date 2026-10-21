# Worked example: fisher density from ear tags and trail cameras (module `scr-smr`)

Appendix S1 worked example and website "try an example" walkthrough.
**The data are simulated with known truth** (see `README_data.md`), so every
estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "We live-trapped fishers in May and put coloured, numbered ear tags on every
> one we caught. Nine also got GPS collars. Then we ran 49 baited trail cameras
> all summer. Some photos show a tag we can read. Some show a tag but we can't
> read the number. Most fishers have no tag at all. How many fishers live here?
> Is density higher where the canopy is denser? And do the cameras on skid
> roads get more photos? One worry: our box traps were all in the middle of
> the camera grid, so the tagged fishers are not a random sample."

## (b) The files the researcher has

`cameras.csv`: one row per camera (49 rows):
```
"camera_id","utm_e","utm_n","on_trail","weeks_active","bait","notes"
"FC01",516000,5216000,1,"1111111111","chicken + skunk lure","skid road / trail"
"FC02",517500,5216000,0,"1111111111","chicken + skunk lure","off-trail, bait tree"
"FC03",519000,5216000,0,"1111111111","chicken + skunk lure","off-trail, bait tree"
"FC04",520500,5216000,1,"1111111111","chicken + skunk lure","skid road / trail"
```

`fisher_photos.csv`: one row per independent photo event (360 rows):
```
"camera_id","week","date","tag_status","fisher_id","photo"
"FC01",1,"2025-06-03","no tag",,"FC01_20250603_0189.JPG"
"FC05",1,"2025-06-03","ear-tag read","FI-022","FC05_20250603_0947.JPG"
"FC06",1,"2025-06-03","ear-tag read","FI-023","FC06_20250603_7745.JPG"
"FC12",1,"2025-06-03","ear-tag, number unreadable",,"FC12_20250603_3812.JPG"
```

`box_traps.csv`: the 25 marking traps:
```
"trap_id","utm_e","utm_n","trap_type"
"BX01",517500,5217500,"Tomahawk 108 box trap"
"BX02",519000,5217500,"Tomahawk 108 box trap"
"BX03",520500,5217500,"Tomahawk 108 box trap"
"BX04",522000,5217500,"Tomahawk 108 box trap"
```

`box_trap_captures.csv`: one row per capture (57 rows, 23 fishers):
```
"fisher_id","night","date","trap_id","sex","new_capture","mass_kg"
"FI-001",1,"2025-05-12","BX04","F",TRUE,2.2
"FI-001",2,"2025-05-13","BX04","F",FALSE,2.3
"FI-002",1,"2025-05-12","BX06","M",TRUE,4.2
"FI-002",2,"2025-05-13","BX12","M",FALSE,5.2
```

`gps_fixes.csv`: GPS fixes of the 9 collared fishers (270 rows):
```
"fisher_id","fix_time","utm_e","utm_n","hdop"
"FI-001","2025-06-04 20:25",520816,5216415,1.1
"FI-001","2025-06-06 23:30",521992,5217100,2.8
"FI-001","2025-06-08 20:54",521560,5214867,2.1
"FI-001","2025-06-10 09:12",523686,5216155,2.2
```

`habitat_grid.csv`: 500 m habitat pixels with canopy cover (1,156 rows):
```
"x","y","canopy_cover_pct","elevation_m","habitat"
512250,5212250,26.5,919,1
512750,5212250,30,883,1
513250,5212250,34.1,889,1
513750,5212250,38.2,859,1
```

## (c) What the researcher types to Claude

> "I've put my fisher data in a folder: box-trap captures and trap locations
> from May, trail camera locations and a photo log from June to August, GPS
> fixes from 9 collared fishers, and a habitat grid with canopy cover. Some
> photos show a readable ear tag, some a tag we couldn't read, and most no tag.
> Please estimate fisher density and abundance, tell me whether density
> rises with canopy cover, and whether trail cameras get more photos. The box
> traps were only in the middle of the camera grid. Please check the method
> works for a design like mine before you trust it, and give me a report and
> figures I can use in a paper."

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

**Step 1. Project folder.** The raw files are copied, never edited. All work
happens in the project folder. (In this example the files live in
`modules/scr-smr/example/` and outputs go to `results/scr-smr/si/`.)

**Step 2. Data intake (eco-data-intake).** The agent reads every file and
describes the data back in plain words:

- 49 cameras on a 7 x 7 grid, 1.5 km apart. 24 are on trails, 25 off-trail.
  3 cameras failed part-way (`weeks_active` has zeros).
- 10 weeks of resighting. 360 photo events: 194 with a readable tag,
  54 with a tag that could not be read, 112 with no tag.
- 25 box traps, 5 nights, 57 captures of 23 fishers. All 23 were later
  photographed with a readable tag. 22 of them at more than one camera.
- 270 GPS fixes from 9 collared fishers, about one every 2 days. That is
  sparse enough to treat the fixes as independent.
- 1,156 habitat pixels (289 km2): the cameras plus 4 km all round.
  Canopy cover runs from about 15% to 80% (mean ~50%, SD ~12%).
- No missing coordinates, no unknown tag codes, every tagged ID in the photos
  appears in the capture file.

**Step 3. Choose the module.** Some animals are marked, most are not, and
cameras photograph both. That is spatial mark-resight, module `scr-smr`.
The box traps sat only in the middle of the camera grid. So tagged fishers
are not a random sample, and the marking captures must be modelled too:
**generalized SMR** (Whittington et al. 2018), `settings: model: gsmr`.
The agent reads the module `README.md` and the Known limits in `VALIDATION.md`.

**Step 4. Analysis plan (eco-analysis-plan).**

- Estimands: density (fishers per km2) and abundance N in the 289 km2 state
  space; the marked fraction; sigma (movement scale) and 95% home-range area.
- Density model: log-linear in canopy cover (standardised, per 1 SD).
- Resighting: Poisson counts of photo events per camera and week. Rate
  lambda0 at the activity centre, half-normal fall-off with distance,
  `on_trail` effect on the log scale.
- Unreadable tags: identified with probability delta; the rest are
  "marked, identity unknown" counts.
- Marking: binomial capture at box traps with probability pm0 at the
  activity centre, modelled for every fisher (gSMR).
- Telemetry: GPS fixes share the same sigma.
- Priors: psi, delta, pm0 ~ U(0,1); sigma ~ U(0, 8,500 m); lambda0 ~ U(0, 10);
  coefficients ~ Normal(0, SD 2).
- Checks: R-hat and bulk ESS; posterior predictive Freeman-Tukey checks;
  an independent maximum-likelihood cross-check with `secr`.

**Step 5. Planted truth at this design (eco-planted-truth).** Before touching
the real data, the agent simulates 100 data sets at exactly this design
(`recovery.yaml`: 49 cameras, 25 box traps, 5 nights, 10 weeks, 40% of
tagged fishers collared with 30 fixes) with known values, and fits each one.
One fit took about 2 to 5 minutes, so 100 replicates were used instead of
200. All 100 fits converged. The gate passed for every quantity:

| quantity | truth | mean estimate | relative bias | 95% coverage | power | gate |
|---|---|---|---|---|---|---|
| abundance N | 73.0 | 76.1 | +4.6% | 0.96 | | PASS |
| density (per km2) | 0.253 | 0.263 | +4.6% | 0.96 | | PASS |
| marked fraction | 0.302 | 0.295 | -1.4% | 0.96 | | PASS |
| sigma (m) | 1000 | 1000 | +0.04% | 0.97 | | PASS |
| 95% home range (km2) | 18.8 | 18.9 | +0.2% | 0.97 | | PASS |
| resight rate lambda0 | 0.30 | 0.309 | +3.0% | 0.93 | | PASS |
| P(tag readable) delta | 0.80 | 0.797 | -0.4% | 0.93 | | PASS |
| box-trap capture pm0 | 0.30 | 0.308 | +2.6% | 0.94 | | PASS |
| canopy effect (log, per SD) | 0.50 | 0.509 | +1.8% | 0.97 | 0.77 | PASS |
| on-trail effect (log) | 0.70 | 0.686 | -2.1% | 0.91 | 1.00 | PASS |

(Coverage must be at least 0.885 with 100 replicates.) What this means for
the researcher: at this design the method is nearly unbiased. Abundance runs
about 5% high on average, which is the known small-sample skew. The canopy
effect is detected 77% of the time if it is as strong as assumed, so a
"no clear effect" result would not prove there is none. The trail effect is
always detected.

**Step 6. Fit.** The agent writes `config.yaml` (file names, column names,
which tag codes mean "read", "unreadable" and "no tag", the habitat grid and
covariates, `model: gsmr`) and runs the fit. Nimble ran 3 chains of 24,000
iterations after 6,000 burn-in, thinned by 4 (M = 145 augmented rows).

**Step 7. Checks (eco-convergence, eco-model-checking).**

- Convergence: max R-hat 1.002; min bulk ESS 2,310. Trace plots show
  three well-mixed chains. PASS.
- Posterior predictive p-values (extreme if < 0.025 or > 0.975):
  unmarked photos total 0.65; identified photos total 0.45; marked fishers
  resighted 0.90; unreadable-tag photos 0.50; Freeman-Tukey by camera
  (unmarked) 0.75; Freeman-Tukey by tagged fisher 0.64. No sign of lack of fit.
- Sightings per tagged fisher: the observed bars sit inside the predictive
  intervals, so there is no sign of strong individual differences in
  how often fishers are photographed.
- Design: camera spacing 1,500 m is 1.5 x sigma (should be at most 2);
  the buffer is 4 x sigma (should be at least 3).
- Cross-check with `secr`: sigma (secr ~1,075 m vs 1,010 m), lambda0 and
  pm0 agree. secr's abundance is lower (about 59). That is expected:
  `VALIDATION.md` section 8 shows secr's joint gSMR fit runs ~26% low on
  these designs, so it is used only to check sigma and detection rates.

**Step 8. Interpretation (eco-interpretation).** Excerpt from the report
(`results/scr-smr/si/fit/scr-smr_report.html`):

> "Density was estimated at 0.30 (95% CrI 0.24 to 0.38) fishers per km²
> (30.1 per 100 km2). Over the 289 km² state space (camera array plus buffer)
> that is 87.1 (95% CrI 68.0 to 109) fishers, of which the 23 marked animals
> are a fraction 0.27 (95% CrI 0.21 to 0.34)."

> "Density and Canopy cover (%): each 1 SD (12 units) higher Canopy cover (%)
> is associated with 64% higher density (rate ratio 1.6, 95% interval 1.2 to
> 2.2). Expected density is 0.15 per km² (0.082 to 0.24) where Canopy cover
> (%) = 36.3 (10th percentile of the state space) and 0.59 (0.36 to 0.87)
> where it = 68.8 (90th percentile). ... This is an association between
> density and habitat across the state space, not proof that the habitat
> causes higher density."

> "Detection: 'Camera on trail' (1 vs 0) multiplied the sighting rate by 2.3
> (95% interval 1.8 to 3.0). This is a detection effect (how easily animals
> are photographed), not a difference in where animals live."

> "Marking phase: a marking trap at an animal's activity centre caught it on a
> given occasion with probability 0.28 (95% CrI 0.21 to 0.35). Because marking
> was modelled, marked animals are not assumed to be a random sample."

Figures (`results/scr-smr/si/fit/figures/`, PNG and PDF):

| figure | caption |
|---|---|
| `density_canopy_cover_pct` | **Key figure.** Expected fisher density versus canopy cover, with 95% and 50% credible bands. Rug: canopy values across the state space. |
| `density_surface` | Posterior mean realised density (fishers per km2). Circles: cameras; triangles: box traps; orange crosses: activity centres of the 23 tagged fishers. |
| `density_expected_surface` | Expected density implied by the canopy-cover model. |
| `abundance_posterior` | Posterior of N in the 289 km2 state space; dashed: 95% interval; dotted: the 23 tagged fishers (lower bound). |
| `marked_fraction` | Posterior of the share of fishers carrying tags (23 / N). |
| `detection_function` | Expected photos per week versus distance, for off-trail and on-trail cameras, and nightly box-trap capture probability; dashed: camera spacing. |
| `coefficients` | Canopy effect on density (log scale, per SD) and trail effect on photo rate (log scale), with 95% intervals. |
| `captures_per_individual` | Observed versus predicted number of tagged fishers identified 0, 1, 2, ... times. |
| `unmarked_by_camera` | Untagged photos per camera, observed versus expected, with 95% predictive intervals. |
| `gof` | Posterior predictive checks (six panels). |
| `engine_comparison` | Nimble gSMR versus secr mark-resight estimates. |
| `trace`, `posterior` | MCMC traces and posterior densities. |

Two presentation fixes were made in the `figures:` block of `config.yaml`:
the trail label was shortened (the strip had read "Camera, Camera on trail = 0"),
and the marked-fraction axis was limited to 0-0.6 (the 0.6-1 range was empty).
The figures were then redrawn from the saved fit with `core/run.R figures`.
(The images embedded in the HTML report were drawn at fit time with the
default labels; the PNG/PDF files carry the fixes.)

**Step 9. Acceptance gate (eco-acceptance-gate).** Planted truth PASS at this
design (10/10 quantities); converged; no lack of fit; design within the
recommended spacing and buffer; interpretation states uncertainty and limits.
**Verdict: DONE.**

## (e) Estimates versus the known truth

The realised population in the simulation had 86 fishers (expected 73;
N is a Poisson draw). Every 95% interval contains the truth.

| quantity | truth | estimate (95% CrI) | inside? |
|---|---|---|---|
| abundance N | 86 | 87.1 (68 to 109) | yes |
| density (per km2) | 0.298 | 0.301 (0.235 to 0.377) | yes |
| untagged fishers | 63 | 64.1 (45 to 86) | yes |
| marked fraction | 0.267 | 0.268 (0.211 to 0.338) | yes |
| sigma (m) | 1000 | 1010 (966 to 1057) | yes |
| 95% home range (km2) | 18.8 | 19.2 (17.6 to 21.0) | yes |
| resight rate lambda0 | 0.30 | 0.259 (0.205 to 0.322) | yes |
| P(tag readable) | 0.80 | 0.780 (0.727 to 0.829) | yes |
| box-trap capture pm0 | 0.30 | 0.277 (0.211 to 0.352) | yes |
| canopy effect (log, per SD) | 0.50 | 0.494 (0.183 to 0.795) | yes |
| on-trail effect (log) | 0.70 | 0.833 (0.583 to 1.085) | yes |

## (f) Limits that matter at this design (from `VALIDATION.md`)

- **Marking must be modelled here.** Treating these tagged fishers as a random
  sample (conventional SMR) underestimated N by 31% with 10% coverage in
  validation. Always give the box-trap captures (gSMR).
- **Small upward bias of about +4-7% in N and density** (skewed posterior).
  It was +4.6% at this design, inside the gate. Report the interval.
- **Do not use secr's abundance as the reference for gSMR.** secr's joint fit
  ran about 26% low in validation; use it only for sigma and detection rates.
- **Not modelled:** tag loss, misread tags (a tagged fisher scored "no tag"),
  individual differences in photo rate (e.g. collar effects), double-counted
  photo events, and births, deaths or movement in and out during the 10 weeks.
- **GPS fixes must be thinned** to be roughly independent (one per ~2 days here).
- **With 5 or fewer tagged animals** N runs ~15% high and the gate fails.
  This study has 23, well clear of that.
- Only the Poisson photo-count model was validated at nsim = 100.

## (g) Commands to reproduce (run from `ecotools/`)

```bash
# 0. regenerate the example data (seed 2025; decorative fields seed 7) - optional
scripts/ejob -c 1 -m 8 -t 0:30:00 -l logs/scr-smr_make_example.log -- scripts/Rscript dev/scr-smr/make_example.R
# 1. module selftest (timing; ~2-5 min per fit at this design)
scripts/ejob -c 4 -m 16 -t 1:00:00 -l logs/scr-smr_selftest.log -- scripts/Rscript core/run.R selftest scr-smr --cores 4
# 2. planted-truth recovery at this design (100 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/scr-smr_si_recover.log -- scripts/Rscript core/run.R recover scr-smr --config modules/scr-smr/example/recovery.yaml --nsim 100 --cores 16 --out results/scr-smr/si/recovery
# 3. fit the example data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/scr-smr_si_fit.log -- scripts/Rscript core/run.R fit scr-smr --config modules/scr-smr/example/config.yaml --out results/scr-smr/si/fit
# 4. redraw figures after editing the figures: block (no refit)
scripts/ejob -c 2 -m 8 -t 0:30:00 -l logs/scr-smr_si_figures.log -- scripts/Rscript core/run.R figures scr-smr --config modules/scr-smr/example/config.yaml --out results/scr-smr/si/fit
```

Machine-readable summary: `results/scr-smr/si/si_summary.json`.
