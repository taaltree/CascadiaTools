# Worked example: adult song sparrow survival from banding data (module `cjs`)

Appendix S1 worked example and website "try an example" walkthrough.
**The data are simulated with known truth** (see `README_data.md`), so every
estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "I run a banding station. Every summer since 2016 we mist-net song sparrows and
> band every new adult. I want to know what fraction of adults survive from one
> year to the next. Do males survive better than females? Do birds die more after
> cold winters? Are heavier birds more likely to come back? A lot of birds are only
> ever caught once. I don't know if they died or were just passing through."

## (b) The files the researcher has

`sosp_captures.csv`: one row per capture (1,283 rows):
```
"band_number","capture_date","net","bander"
"2841-22465","2016-05-10","N03","KL"
"2841-80867","2016-05-10","N07","JT"
"2841-42326","2016-05-11","N02","MR"
"2841-59941","2016-05-11","N07","MR"
```

`sosp_birds.csv`: one row per bird (575 rows):
```
"band_number","sex","mass_g","wing_chord_mm","first_year"
"2841-50307","M",21.7,65,2016
"2841-14576","F",18.1,64,2016
"2841-85173","F",21.2,62,2016
"2841-88116","M",21.1,65,2016
```

`sosp_years.csv`: one row per year (9 rows):
```
"year","winter_min_temp_c","spring_precip_mm"
2016,-3.1,212
2017,-1.2,188
2018,-5.8,240
2019,-0.4,199
```

## (c) What the researcher types to Claude

> I have nine years of song sparrow banding data from my station (2016 to 2024) in
> three files: sosp_captures.csv (every capture, with band number and date),
> sosp_birds.csv (sex and body mass of each bird when first banded) and
> sosp_years.csv (mean January minimum temperature for each winter). I want annual
> adult survival, whether it differs between males and females, and whether it is
> related to winter cold and body mass. Many birds are only caught once; I think some
> are migrants passing through. Please analyse this and give me a report with figures
> I can use in a paper.

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

**1. Project folder.** The data are copied to a project folder. The raw files are never edited.

**2. Data intake (eco-data-intake).** The agent reads the three files and describes them back:
- 575 adult birds, banded 2016-2023 (66 to 78 new birds per year). 301 are males (52%).
- 1,283 captures, all between 10 May and 8 August. Repeat captures in the same summer are
  merged, which leaves 813 bird-seasons over 9 annual occasions.
- 158 birds (27%) were caught again in a later year. The other 417 were seen only once.
- Body mass is known for every bird: mean 21.9 g (SD 1.66; range 15.8-27.7 g).
- Winter temperature is given for every year: -6.5 to -0.4 deg C. The value for a year
  applies to the winter after that summer, so the 2024 value is not used.
- No missing values, no duplicate bands, no captures outside the season.

**3. Choose the module.** The data are individually marked animals with repeat captures,
and the question is survival. So the module is `cjs` (Cormack-Jolly-Seber). It
estimates survival and recapture probability separately. A bird that is not caught
may still be alive.

**4. Analysis plan (eco-analysis-plan), written before fitting.**
- Estimands: annual apparent survival of residents (females, males); recapture
  probability (females, males); effects of sex, winter temperature and body mass on
  survival; proportion of new birds that are transients.
- Model: survival `~ tsm + sex + winter_min_temp_c + mass_g`; recapture `~ sex`.
  `tsm` (time since marking) lets the first year after banding differ from later
  years. This is how transients are handled. The module crosses `tsm` with the
  other terms, so the sex, winter and mass effects describe resident birds.
- Engine: maximum likelihood (`marked`). Goodness of fit: U-CARE tests via R2ucare.
- "Adequate fit": an overall c-hat near 1, and any failing component test accounted
  for by the model.
- Caveat stated up front: "apparent survival" = true survival x staying at the site.

**5. Planted truth at this design (eco-planted-truth).** Before touching the real data,
the agent simulates 200 data sets that look like this study: 9 years, 72 new birds per
year, 52% males, 25% transients, plausible survival and recapture values. It fits the
same model to each one. Each fit takes about 4 seconds. Result: **all 10 estimands PASS**
(200/200 fits converged; coverage must be at least 0.904):

| estimand | truth | mean estimate | relative bias | 95% coverage | power | gate |
|---|---|---|---|---|---|---|
| survival, female | 0.550 | 0.548 | -0.3% | 0.94 | | PASS |
| survival, male | 0.623 | 0.621 | -0.2% | 0.96 | | PASS |
| recapture, female | 0.450 | 0.452 | +0.4% | 0.94 | | PASS |
| recapture, male | 0.574 | 0.575 | +0.1% | 0.96 | | PASS |
| sex effect on survival (logit) | 0.30 | 0.302 | +0.7% | 0.95 | 0.17 | PASS |
| winter severity effect (logit/SD) | -0.40 | -0.405 | +1.2% | 0.955 | 0.56 | PASS |
| body mass effect (logit/SD) | 0.25 | 0.274 | +9.7% | 0.965 | 0.51 | PASS |
| sex effect on recapture (logit) | 0.50 | 0.503 | +0.5% | 0.965 | 0.37 | PASS |
| proportion transient, female | 0.25 | 0.228 | -8.7% | 0.96 | | PASS |
| proportion transient, male | 0.25 | 0.248 | -0.8% | 0.96 | | PASS |

What this tells the researcher. The method gives unbiased answers with honest intervals
at this design. But the design has **low power** for the effects: a real sex difference
of this size is detected only 17% of the time, and the winter and mass effects only about
half the time. So "no clear effect" must not be read as "no effect". (The simulator calls
its weather covariate "winter severity", so its sign is flipped: -0.4 per SD of severity
equals +0.4 per SD of milder winters.)

**6. Fit.** The agent writes `config.yaml` (file names, columns, formulas, figure labels)
and runs the fit (under a minute).

**7. Checks (eco-convergence, eco-model-checking).**
- Convergence: the likelihood optimiser converged with an invertible Hessian. No
  probability sits on a boundary (all between 0.28 and 0.70), and every SE is finite and
  small (largest logit SE 0.49). (The report's line "max R-hat NA; min bulk ESS NA" just
  means these MCMC diagnostics do not apply to a maximum-likelihood fit.)
- Goodness of fit (R2ucare, time-dependent CJS): overall chi-square 41.4 on 45 df,
  p = 0.624, **c-hat = 0.92**. Test **3.SR flags transients** (chi-square 25.8, df 14,
  p = 0.027): newly banded birds come back less often than birds banded earlier. The
  model includes `tsm` for exactly this. Without 3.SR, the remaining tests give c-hat 0.50
  (p = 0.990). Trap dependence (2.CT) p = 0.801, 3.SM p = 0.986, 2.CL p = 0.783.
  No sign of remaining lack of fit.

**8. Interpretation (eco-interpretation).** Excerpts from the report
(`results/cjs/si/fit/cjs_report.html`):

> Apparent annual survival at average mean January minimum temperature (°C) and body mass
> (g) of residents (intervals after the first one following marking): Female 0.54 (95% CI
> 0.43-0.65); Male 0.62 (95% CI 0.53-0.70). Apparent survival is the probability of
> surviving AND staying in the study area; animals that emigrate permanently are
> indistinguishable from deaths, so true survival is at least this high.

> Sex: Male vs Female (apparent survival): odds ratio 1.34 (95% CI 0.76-2.40) (logit
> coefficient 0.30). The interval includes 1: the data do not show a clear difference.

> Mean January minimum temperature (°C): each 1 SD (2.04 deg C) increase multiplies the
> odds of apparent survival by 1.30 (95% CI 0.91-1.87). Holding other covariates at their
> means (Female), apparent survival is 0.45 at the 10th percentile (-6.0 deg C) and 0.61 at
> the 90th percentile (-1.0 deg C). The interval includes no effect. This is an
> association, not proof of a causal effect.

> Body mass (g): each 1 SD (1.66 g) increase multiplies the odds of apparent survival by
> 1.17 (95% CI 0.90-1.52). [...] apparent survival is 0.49 at the 10th percentile (19.8 g)
> and 0.59 at the 90th percentile (23.9 g). The interval includes no effect.

> Transients: an estimated 0.26 (95% CI 0.00-0.53) (Female); 0.32 (95% CI 0.11-0.48)
> (Male) of newly marked animals were transients (passing through and never returning).

> Detection (recapture) probability per occasion: Female 0.41 (95% CI 0.28-0.55); Male
> 0.57 (95% CI 0.47-0.65).

**Figures** (`results/cjs/si/fit/figures/`, PNG + PDF; gallery in `gallery.html`):

| figure | caption |
|---|---|
| `survival_winter_min_temp_c` | Apparent survival versus mean January minimum temperature, one curve per sex (line: estimate; bands: 50% and 95% CI); body mass at its mean. Rug: the 8 observed winters. |
| `survival_mass_g` | Apparent survival versus body mass, one curve per sex; winter at its mean. Rug: observed masses (one per bird). |
| `survival_effects_panel` | Both conditional effects side by side (A: winter temperature, B: body mass), for a paper. |
| `survival_time` | Apparent survival per interval (2016-17 to 2023-24) by sex, predicted from that winter's temperature only (no separate year effect). |
| `detection_groups` | Recapture probability by sex (95% CI). |
| `coefficients` | Odds ratios (95% CI, log scale) for sex, winter temperature, body mass, the first interval after marking, and sex on recapture. Dashed line: no effect. |
| `gof` | Goodness-of-fit components (chi-square/df); 3.SR (transience) is the only test with p < 0.05. Overall c-hat 0.92. |
| `captures` | Birds caught per year: newly banded versus recaptured. |

The key figure is `survival_effects_panel`. It answers the winter and body mass questions,
in original units, separately for each sex.

**9. Acceptance gate (eco-acceptance-gate).** The plan was declared before fitting. The
planted-truth check passed at this design (10/10). The fit converged. Goodness of fit was
adequate, and the one flag (transients) is handled. The module has been cross-checked
against independent software (wiqid, which was checked against MARK; a hand-coded
likelihood; Nimble) and a negative control fails as it should (VALIDATION.md). Every number
comes from the commands in (g). **Verdict: DONE.**

**10. Deliver.** Summary given to the researcher:
> Resident adult song sparrows had an apparent annual survival of 0.54 (95% CI 0.43-0.65)
> for females and 0.62 (0.53-0.70) for males, at average body mass and winter temperature.
> About a quarter to a third of newly banded birds were transients that never returned
> (female 0.26, 95% CI 0.00-0.53; male 0.32, 0.11-0.48); without accounting for them,
> survival would have been underestimated by about 13%. Survival tended to be higher after
> milder winters (odds ratio 1.30 per 2 deg C, 95% CI 0.91-1.87) and for heavier birds (1.17
> per 1.7 g, 0.90-1.52), and males tended to survive better (1.34, 0.76-2.40), but none of
> these intervals exclude "no effect". The simulation check shows this study design has only
> 17-56% power for effects of this size, so these are inconclusive results, not evidence of no
> effect. "Survival" includes staying at the site: birds that moved away count as lost.

## (e) Estimates versus the known truth

| parameter | truth | estimate (95% CI) | truth inside CI? |
|---|---|---|---|
| survival, resident female | 0.550 | 0.544 (0.431-0.653) | yes |
| survival, resident male | 0.623 | 0.616 (0.526-0.699) | yes |
| recapture, female | 0.450 | 0.410 (0.284-0.548) | yes |
| recapture, male | 0.574 | 0.565 (0.472-0.654) | yes |
| sex effect on survival (logit) | 0.30 | 0.296 (-0.281 to 0.873) | yes |
| winter temperature effect (logit per SD) | 0.40 | 0.266 (-0.094 to 0.626) | yes |
| body mass effect (logit per SD) | 0.25 | 0.156 (-0.105 to 0.416) | yes |
| sex effect on recapture (logit) | 0.50 | 0.629 (-0.042 to 1.299) | yes |
| proportion transient, female | 0.25 | 0.259 (-0.120 to 0.528) | yes |
| proportion transient, male | 0.25 | 0.319 (0.113-0.481) | yes |

All 10 true values lie inside the 95% intervals. The winter and mass effects are
estimated a little low here, but by less than 1 SE. This is one data set: the 200-replicate
check in step 5 shows no systematic bias. The fact that the winter and mass intervals
include zero, even though the true effects are real, is exactly the low power predicted in step 5.

## (f) Limits that matter at this design (from `VALIDATION.md`)

- **Low power for effects.** With 575 birds over 9 years and recapture about 0.4-0.6,
  power is 17% for the sex effect and 51-56% for the winter and mass effects. Report
  intervals. Do not claim "no effect".
- **Transients must be modelled.** If `tsm` is dropped, survival is biased about 13% low
  and the 95% intervals cover the truth only about half the time (validated negative
  control). Here 3.SR flagged transients (p = 0.027) and `tsm` was kept.
- **`proportion_transient` is an unconstrained ratio.** The female table interval goes below
  0 (-0.12). The report text truncates it to 0-0.53. The female value is poorly estimated.
- **False-positive rate for individual covariates** (body mass) is about 7.5-8% at a nominal
  5%. Treat a lone mass effect with p just under 0.05 with caution.
- **Winter is one value per year.** The winter effect rests on only 8 winters. It is an
  association, not a causal effect.
- **Not modelled:** individual differences in catchability beyond sex, band loss, temporary
  emigration. "Apparent survival" = true survival x site fidelity.

## (g) Commands to reproduce (run from `ecotools/`)

```bash
# 0. regenerate the example data (seed 2031) - optional, the CSVs are already in the folder
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/cjs_make_example.log -- scripts/Rscript dev/cjs/make_example.R
# 1. module selftest (4 replicates at module defaults; used only for timing, ~4 s per fit;
#    4 replicates are too few for the gate, so its PASS/FAIL column is not meaningful)
scripts/ejob -c 4 -m 8 -t 0:30:00 -l logs/cjs_si_selftest.log -- scripts/Rscript core/run.R selftest cjs --out dev/cjs/si_selftest
# 2. planted-truth recovery at this design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/cjs_si_recover.log -- scripts/Rscript core/run.R recover cjs --config modules/cjs/example/recovery.yaml --nsim 200 --cores 16 --out results/cjs/si/recovery
# 3. fit the example data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/cjs_si_fit.log -- scripts/Rscript core/run.R fit cjs --config modules/cjs/example/config.yaml --out results/cjs/si/fit
# 4. machine-readable summary for the SI / website
python3 dev/cjs/make_si_summary.py
```

Outputs: `results/cjs/si/recovery/cjs_recover_summary.csv` (gate table),
`results/cjs/si/fit/cjs_report.html` (report), `results/cjs/si/fit/cjs_estimates.csv`,
`results/cjs/si/fit/figures/`, and `results/cjs/si/si_summary.json`.
