# Worked example: how does wildfire affect songbirds, directly and through vegetation? (module `sem`)

Appendix S1 worked example and website "try an example" walkthrough.
**The data are simulated with known truth** (see `README_data.md`), so every
estimate below can be checked against the real answer.

---

## (a) The question, in the researcher's words

> "Last summer I surveyed 120 plots across a big wildfire, in 12 burn units.
> I know fire kills trees. I think losing the canopy lets shrubs come back, and
> shrubs are good for a lot of songbirds. But severe fire might also be bad for
> birds directly. So does fire severity hurt or help birds overall? How much of
> its effect goes through canopy and shrubs? Wetter plots also seem shrubbier.
> Does soil moisture matter for the birds?"

## (b) The files the researcher has

One file, `burned_forest_plots.csv`, one row per plot (120 rows):
```
"plot_id","burn_unit","survey_date","fire_severity_dnbr","canopy_cover_pct","soil_moisture_pct","shrub_cover_pct","bird_detections"
"P001","BU01","2025-05-24",732,28.1,19.2,47.9,2
"P002","BU01","2025-06-04",633,33.4,23.3,48.1,7
"P003","BU01","2025-06-06",468,47.1,16.4,42.5,7
"P004","BU01","2025-06-05",484,56.1,22.8,34.1,5
```

## (c) What the researcher types to Claude

> I have a CSV of 120 forest plots from 12 burn units after a wildfire
> (burned_forest_plots.csv). For each plot I have fire severity (dNBR), canopy
> cover, soil moisture, shrub cover and the number of birds detected in one point
> count. I think fire reduces canopy, less canopy means more shrubs, more shrubs
> means more birds, wetter soil means more shrubs, and fire might also reduce
> birds directly. Can you test this diagram, tell me how strong each link is, and
> how much of the fire effect on birds goes through the vegetation? Plots in the
> same burn unit are not independent. Please give me a report and figures for a
> paper.

## (d) What the agent does, step by step (`skills/eco-analyze/SKILL.md`)

**1. Project folder.** The CSV is copied to a project folder. The raw file is never edited.

**2. Data intake (eco-data-intake).** The agent reads the file and describes it back:
- 120 plots, 12 burn units with exactly 10 plots each. No missing values.
- Fire severity 32-1008 dNBR (mean 497, SD 212). Canopy 10.5-89.0 % (mean 52.2).
  Soil moisture 9.3-36.7 % (mean 22.1). Shrub cover 22.6-59.5 % (mean 40.7).
- Bird detections are whole numbers, 1-17 per count (mean 5.7). Variance/mean = 1.29,
  so a Poisson count model is reasonable (the module warns above 1.5). No zeros.
- One count per plot, 20 May-14 June 2025. Survey date is not in the diagram.

**3. Choose the module.** Several measured variables influence each other as a network,
and the question asks how much of an effect runs through mediators. So the module is
`sem` (structural equation model), piecewise engine. It allows a Poisson count outcome
and a random intercept for burn unit. (lavaan is not used: it handles only continuous
outcomes and no random effects.)

**4. Analysis plan (eco-analysis-plan), written before fitting.**
- Diagram, one model per response:
  `canopy ~ fire` (Gaussian); `shrub ~ canopy + moisture` (Gaussian);
  `birds ~ shrub + fire + (1|burn_unit)` (Poisson, log link).
- Estimands: each path; direct, indirect and total effects of fire on birds, of
  moisture on birds, and of fire on shrubs. 95% Monte Carlo intervals.
- Global check: Shipley's d-separation test (Fisher's C). "Adequate fit" = P >= 0.05.
- Reporting scale: raw paths and rate ratios ("birds multiplied by ... per 100 dNBR").
  Standardised values are used only to compare the size of routes.
- Caveats stated up front: bird counts mix abundance with detectability; paths are
  associations unless the diagram is causally complete.

**5. Planted truth at this design (eco-planted-truth).** Before touching the real data,
the agent simulates 200 data sets that look like this study: 120 plots in 12 burn units,
the same diagram, Poisson birds with a burn-unit effect (SD 0.25), and realistic values
(`example/recovery.yaml`). It fits the same model to each one. One fit takes under a
second (median 0.8 s). All 200 fits converged. Coverage must be at least 0.904.

**Result: 18 of 24 estimands PASS. 6 FAIL**, all of them *standardised* values that
involve the bird count model:

| estimand | truth | mean estimate | rel. bias | bias MCSE | coverage | power | gate |
|---|---|---|---|---|---|---|---|
| path fire -> canopy | -0.055 | -0.0546 | +0.7% | 0.00025 | 0.955 | 1.00 | PASS |
| path canopy -> shrub | -0.35 | -0.352 | -0.6% | 0.0031 | 0.950 | 1.00 | PASS |
| path moisture -> shrub | 0.6 | 0.602 | +0.3% | 0.0089 | 0.970 | 1.00 | PASS |
| path shrub -> birds | 0.025 | 0.0259 | +3.4% | 0.00037 | 0.935 | 1.00 | PASS |
| path fire -> birds (= direct) | -0.0006 | -0.00062 | -2.6% | 0.000014 | 0.970 | 0.87 | PASS |
| indirect fire -> birds | 0.00048 | 0.00050 | +3.1% | 0.0000085 | 0.940 | 1.00 | PASS |
| total fire -> birds | -0.00012 | -0.00012 | -0.8% | 0.000013 | 0.980 | 0.09 | PASS |
| indirect (= total) moisture -> birds | 0.015 | 0.0156 | +3.7% | 0.00031 | 0.950 | 1.00 | PASS |
| indirect (= total) fire -> shrub | 0.01925 | 0.0192 | -0.2% | 0.00018 | 0.970 | 1.00 | PASS |
| std fire -> canopy | -0.814 | -0.809 | +0.7% | 0.0023 | 1.000 | 1.00 | PASS |
| std canopy -> shrub | -0.578 | -0.581 | -0.4% | 0.0040 | 0.990 | 1.00 | PASS |
| std moisture -> shrub | 0.322 | 0.322 | +0.2% | 0.0046 | 0.965 | 1.00 | PASS |
| std indirect (= total) fire -> shrub | 0.471 | 0.470 | -0.3% | 0.0036 | 0.990 | 1.00 | PASS |
| std total fire -> birds | -0.053 | -0.061 | -15.5% | 0.0064 | 0.980 | 0.09 | PASS (within 3 MCSE) |
| **std shrub -> birds** | 0.453 | 0.532 | **+17.4%** | 0.0064 | **0.880** | 1.00 | **FAIL** |
| **std fire -> birds** (and std direct) | -0.266 | -0.311 | **-17.0%** (17% too large in size) | 0.0066 | 0.945 / 0.935 | 0.87 | **FAIL** (2 rows) |
| **std indirect fire -> birds** | 0.213 | 0.250 | **+17.4%** | 0.0038 | 0.935 | 1.00 | **FAIL** |
| **std indirect (and total) moisture -> birds** | 0.146 | 0.171 | **+17.6%** | 0.0031 | 0.965 | 1.00 | **FAIL** (2 rows) |

(24 rows in all. The table merges rows that are identical by definition: direct = path,
and total = indirect where no direct path exists. Power = share of 95% intervals that exclude 0.
Full table: `results/sem/si/recovery/sem_recover_summary.csv`.)

What this tells the researcher:
- **Every raw path, rate ratio and indirect/total effect is recovered.** Bias is at most
  3.7% and coverage is 0.935-0.98. These are the numbers to interpret.
- **Standardised coefficients for the bird model are about 17% too large in size.**
  The agent checked why (`dev/sem/si/std_count_check.R`, 200 data sets). For a count
  response, the module (like piecewiseSEM) divides by sd(eta)/cor(y, fitted). The fitted
  values include the 12 estimated burn-unit effects. These soak up some Poisson noise, so
  the in-sample correlation is too high (0.66 versus 0.62 with the true means), and every
  standardised bird path is inflated by +17%. Using fixed-effect-only fitted values gives
  +2.2% bias. Using the true means gives +6.0%. This is VALIDATION.md Known limit 3
  ("about +7%" at the module's GLMM design). It is larger here because birds have a
  weak signal (Poisson means around 5, small fixed effects). This is a property of the
  standardisation convention, not a coding error. The module matches
  `piecewiseSEM::coefs()` to 5e-5, so module code was not changed.
- Power: the design detects the direct fire effect 87% of the time and every indirect
  route 100% of the time. The total fire -> birds effect is tiny in truth (the routes
  nearly cancel), and it is detected only 9% of the time.

**6. Fit.** The agent writes `config.yaml` (file, three component models, effects to report,
units for the text, labels, figure settings) and runs the fit. It takes under a minute.

**7. Checks (eco-convergence, eco-model-checking).**
- Convergence: all three component models converged. This is a likelihood fit, so R-hat
  does not apply. No singular random effect was reported. Poisson dispersion was not
  flagged (raw variance/mean 1.29; the module warns above 1.5).
- Goodness of fit: **d-separation test, Fisher's C = 11.06, df = 8, P = 0.198**, over 4
  independence claims. So there is no evidence that the diagram leaves out a needed link.
- One single claim had P < 0.05: soil moisture -> birds (P = 0.023). The other three
  claims were well above 0.05 (roughly 0.4-0.9 on the d-sep figure). The agent does not add this link: the global test passes, and
  about 1 in 5 data sets from a true model show at least one claim with P < 0.05
  (VALIDATION.md). Here we know the truth: there is no such link, so this flag is a false alarm.
- R2: canopy 0.68, shrub 0.39, birds 0.23.

**8. Interpretation (eco-interpretation).** Excerpts from the report
(`results/sem/si/fit/sem_report.html`):

> Fire severity (dNBR) -> Canopy cover (%): [...] an increase of 100 units is associated
> with a difference of -5.87 (95% CI -6.61 to -5.14) in Canopy cover (%).

> Canopy cover (%) -> Shrub cover (%): [...] an increase of 10 units is associated with a
> difference of -2.63 (95% CI -3.41 to -1.85) in Shrub cover (%).

> Soil moisture (%) -> Shrub cover (%): [...] an increase of 5 units is associated with a
> difference of 3.29 (95% CI 2.12 to 4.45) in Shrub cover (%).

> Shrub cover (%) -> Bird detections: [...] an increase of 10 units is associated with the
> expected Bird detections multiplied by 1.33 (95% CI 1.20 to 1.47; +33%).

> Fire severity (dNBR) -> Bird detections: [...] an increase of 100 units is associated with
> the expected Bird detections multiplied by 0.91 (95% CI 0.87 to 0.95; -9%).

> Per 1-SD increase in Fire severity (dNBR) (212 units), the indirect route alone is
> associated with the expected Bird detections multiplied by 1.10 (95% CI 1.05 to 1.15;
> +10%); the total association (direct + indirect) corresponds to the expected Bird
> detections multiplied by 0.90 (95% CI 0.83 to 0.99; -10%). The direct (-0.000936) and
> indirect (0.000442) associations have opposite signs and partly cancel [...]; a weak
> total does not mean the driver has little influence.

> The model has no direct path from Soil moisture (%) to Bird detections, so the whole
> association is indirect. [Per 1 SD (5.02 %) of moisture: birds x1.10, 95% CI 1.05-1.16.]

> Detection versus ecology. [...] a path into them can reflect an effect on detection (e.g.
> dense vegetation hiding animals) as well as on the animals; the SEM cannot separate the two.

The agent checked the arithmetic: exp(0.0286 x 10) = 1.33; exp(-0.000936 x 100) = 0.91;
exp(0.000442 x 212) = 1.10; exp(-0.000494 x 212) = 0.90.

**Figures** (`results/sem/si/fit/figures/`, PNG + PDF; gallery in `gallery.html`):

| figure | caption |
|---|---|
| `path_diagram` | The fitted diagram. Numbers are standardised coefficients; arrow width shows size, blue = positive, orange = negative; all paths have 95% CIs excluding 0. The dotted pink arc is the one omitted link with a d-sep P < 0.05 (moisture -> birds, P = 0.023). Standardised bird paths are inflated by about 17% at this design (step 5). |
| `indirect_effects` | Standardised direct, indirect and total effects with 95% Monte Carlo CIs: fire -> birds (direct negative, indirect positive, total negative), moisture -> birds (all indirect), fire -> shrub (all indirect, via canopy). |
| `dsep` | P values (log scale) of the 4 independence claims. One claim (moisture -> birds) is below 0.05; Fisher's C = 11.06, df 8, P = 0.198. |
| `partial_fire_severity_dnbr_to_canopy_cover_pct` | Canopy cover versus fire severity, with 95% band and partial residuals. |
| `partial_canopy_cover_pct_to_shrub_cover_pct` | Shrub cover versus canopy cover, soil moisture at its mean. |
| `partial_soil_moisture_pct_to_shrub_cover_pct` | Shrub cover versus soil moisture, canopy at its mean. |
| `partial_shrub_cover_pct_to_bird_detections` | Expected birds per point count versus shrub cover (about 3.3 at 23 % shrub to 9.4 at 60 %), fire at its mean, burn-unit effect at zero. |
| `partial_fire_severity_dnbr_to_bird_detections` | Expected birds per point count versus fire severity (direct path only; about 8.5 at low to 3.5 at the highest severity), shrub at its mean. |
| `effects_panel` | All five conditional effects (A-E) in one figure for a paper. |

Presentation fix made while checking figures: the bird label "Bird detections (per point
count)" wrapped onto three lines in the diagram and the d-sep plot. In the example config's
`figures: labels:` it is now "Birds per point count". The estimates did not change.

The key figure is `path_diagram`. It shows the whole answer at once: fire cuts canopy, less
canopy means more shrubs, more shrubs mean more birds, and fire also has a direct negative path to birds.

**9. Acceptance gate (eco-acceptance-gate).** The plan was declared before fitting. The fit
converged and the d-sep test passed. The module has been cross-checked against lavaan and
piecewiseSEM, and a negative control (omitted path) fails as it should (VALIDATION.md).
Every number comes from the commands in (g). But the planted-truth gate **failed for 6
standardised estimands** of the bird model (+17% bias). **Verdict: OPEN.**
What it means in practice: raw paths, rate ratios and the indirect/total effects in
original units passed and can be reported. The standardised bird coefficients (0.57, -0.48,
0.23, 0.23 and the -0.25 total) should not be quoted as accurate sizes. If they are used, use
them only to rank routes, with a note that they run about 17% high. To close it, a
fixed-effects-only standardisation for count GLMMs would need to be added and validated
(a module change, outside this example).

**10. Deliver.** Summary given to the researcher:
> More severe fire meant less canopy: -5.9 % canopy per 100 dNBR (95% CI -6.6 to -5.1).
> Less canopy meant more shrubs (+2.6 % shrub per 10 % less canopy, 1.9-3.4), and wetter
> soil meant more shrubs (+3.3 % per 5 % moisture, 2.1-4.5). Plots with more shrub had more
> birds: x1.33 per 10 % shrub (1.20-1.47). Through this vegetation route, more severe fire
> was associated with more birds (x1.10 per 212 dNBR, 1.05-1.15). But fire also had a direct
> negative association with birds (x0.91 per 100 dNBR, 0.87-0.95). Overall, birds were 10%
> fewer per 212 dNBR (x0.90, 0.83-0.99), because the two routes partly cancel. The
> diagram fits the data (Fisher's C = 11.1, df 8, P = 0.20). These are associations from
> one season of counts. Counts also reflect how easily birds are detected in dense shrubs.

## (e) Estimates versus the known truth

Raw paths and effects (the quantities that passed the gate):

| parameter | truth | estimate (95% CI) | truth inside CI? |
|---|---|---|---|
| fire -> canopy (% per dNBR) | -0.055 | -0.0587 (-0.0661 to -0.0514) | yes |
| canopy -> shrub (% per %) | -0.35 | -0.263 (-0.341 to -0.185) | **no** (just outside) |
| moisture -> shrub (% per %) | 0.60 | 0.657 (0.424 to 0.891) | yes |
| shrub -> birds (log per %) | 0.025 | 0.0286 (0.0185 to 0.0387) | yes |
| fire -> birds, direct (log per dNBR) | -0.00060 | -0.00094 (-0.00134 to -0.00049) | yes |
| indirect fire -> birds | 0.00048 | 0.00044 (0.00025 to 0.00066) | yes |
| total fire -> birds | -0.00012 | -0.00049 (-0.00090 to -0.00007) | **no** |
| indirect moisture -> birds | 0.015 | 0.0188 (0.0101 to 0.0293) | yes |
| indirect fire -> shrub | 0.01925 | 0.0155 (0.0107 to 0.0205) | yes |

7 of 9 intervals contain the truth. With 95% intervals, about 1 miss in 20 is expected,
and 2 misses in 9 is within chance for a single data set (the recovery runs show coverage 0.94-0.98
for these quantities). The total fire -> birds miss follows from the direct path. That
estimate is more negative than the truth (-0.00094 vs -0.0006, though its interval contains it), and the total
is direct + indirect. The true total is almost zero (x0.97 birds per 212 dNBR). The data
suggest a clearer decline (x0.90). This is a single-data-set fluctuation, which the 9% power for this
effect predicts.

Standardised values (truth = population values for the recovery design):
fire -> canopy -0.81 vs estimate -0.82 (-0.93 to -0.72); canopy -> shrub -0.58 vs -0.48
(-0.63 to -0.34); moisture -> shrub 0.32 vs 0.40 (0.26 to 0.55); shrub -> birds 0.45 vs
0.57 (0.37 to 0.77); direct fire -> birds -0.27 vs -0.48 (-0.69 to -0.25, misses);
indirect fire -> birds 0.21 vs 0.23 (0.13 to 0.34); total fire -> birds -0.05 vs -0.25
(-0.46 to -0.04, misses); indirect moisture -> birds 0.15 vs 0.23 (0.12 to 0.36). The bird
values sit high in size, as step 5 predicted.

The d-sep test flagged moisture -> birds (P = 0.023). There is no such link in the truth.
This shows why single claims must not be over-read.

## (f) Limits from VALIDATION.md that matter at this design

1. **Standardised count coefficients with a random intercept are biased** (Known limit 3).
   At this design the bias is +17%, which fails the gate (step 5). Report raw paths, rate ratios and
   effects in original units. Treat standardised bird values only as a ranking.
2. **Omitted paths or unmeasured confounders bias indirect effects strongly** (Known limit 5:
   -42% in the negative control). At n = 100 the d-sep test catches an omitted path of
   standardised size ~0.3 only 73% of the time. A passing test (P = 0.198) does not prove the
   diagram. For example, logging or salvage after the fire could drive both canopy and birds.
3. **Detection.** Counts mix bird abundance with detectability. Dense shrubs may hide birds
   or make them easier to hear. The SEM cannot separate these (README "not_for"; use an
   N-mixture or occupancy module with repeat visits for abundance).
4. **Indirect effects through non-Gaussian mediators** (Known limit 4) do not apply here.
   Both mediators (canopy, shrub) are Gaussian, so the raw indirect effects are exact
   products on the log-count scale.
5. Interactions, nonlinear paths (e.g. birds peaking at intermediate severity) and feedback
   loops are not supported (Known limit 6). Check the partial-residual plots for curvature.
   None is obvious here.
6. Small-sample binary outcome bias (Known limit 1) and lavaan limits (Known limit 2) do not
   apply. The outcome is a count, n = 120.

## (g) Commands to reproduce

All from `.`:

```bash
# 0. regenerate the example data (seed 2024; byte-identical check)
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/sem_si_regen.log -- scripts/Rscript dev/sem/si/make_example_check.R
cmp dev/sem/si/regen_burned_forest_plots.csv modules/sem/example/burned_forest_plots.csv

# 1. planted-truth recovery at this design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/sem_si_recover.log -- scripts/Rscript core/run.R recover sem \
  --config modules/sem/example/recovery.yaml --nsim 200 --cores 16 --out results/sem/si/recovery

# 2. why the standardised bird paths fail (diagnostic, 200 replicates)
scripts/ejob -c 8 -m 16 -t 0:30:00 -l logs/sem_si_stdcheck.log -- scripts/Rscript dev/sem/si/std_count_check.R

# 3. fit the example data: report + figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/sem_si_fit.log -- scripts/Rscript core/run.R fit sem \
  --config modules/sem/example/config.yaml --out results/sem/si/fit

# 4. (optional) redraw figures after editing the figures: block, no refit
scripts/Rscript core/run.R figures sem --fit results/sem/si/fit/sem_fit.rds \
  --config modules/sem/example/config.yaml --out results/sem/si/fit
```

Outputs: `results/sem/si/recovery/` (gate summary, replicates, report),
`results/sem/si/fit/` (report, estimates, figures), `results/sem/si/std_count_check.csv`,
`results/sem/si/si_summary.json`.
