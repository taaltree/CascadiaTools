# Worked example: elk movement modes (movement-hmm)

This walkthrough shows a full CascadiaTools (ecotools) analysis from start to finish. A
researcher brings GPS data and a question. The agent picks the analysis, tests it on
simulated data, fits it, checks it, and writes the report.

**The data are simulated** (see `README_data.md`), so we know the right answer. Section (e)
compares the estimates with that answer.

---

## (a) The question, in the researcher's words

> "I collared five elk in the Coast Range with hourly GPS fixes for a month. They seem to
> alternate between hanging around in one spot and moving long distances. How much time do
> they spend doing each? Do they move more at certain times of day? And does forest cover
> change when they settle down or start moving?"

## (b) The files the researcher has

One file, `elk_gps.csv` (3416 rows, one per successful GPS fix):

```
"animal_id","timestamp","longitude","latitude","forest_cover_pct"
"elk_F01","2024-06-01 12:00:33",-123.267139,44.53311,51.4
"elk_F01","2024-06-01 12:58:34",-123.268262,44.536019,52.5
"elk_F01","2024-06-01 13:59:49",-123.271727,44.543074,45.8
"elk_F01","2024-06-01 14:59:24",-123.271775,44.542839,46.2
```

Timestamps are local Pacific time. Failed fixes are simply missing rows.

## (c) What the researcher types to Claude

> "Here is a CSV of hourly GPS fixes from 5 elk (elk_gps.csv, times are Pacific time). I want to
> know whether they switch between resting and travelling, how much time they spend in each, and
> whether time of day and forest cover affect switching. Please analyse it and give me a report
> with figures I can use in a paper."

## (d) What the agent does, step by step (eco-analyze)

**Step 1. Project folder.** The agent creates a project folder with the raw data, a plan, a
config, recovery results and a log of every command. (For this example the files live in
`modules/movement-hmm/example/` and `results/movement-hmm/si/`.)

**Step 2. Data intake.** The agent describes the data back to the researcher:

- 5 elk, 3416 fixes (677 to 690 per animal), 1 June to 1 July 2024.
- The most common gap between fixes is 1 hour. Fix times wander by up to about 2 minutes.
- On a regular 1-hour grid there are 3600 hourly slots. 184 of them (5.1%) have no fix.
- Forest cover at the fixes ranges from 10.5% to 88.0% (mean 52.6%).
- Coordinates are lon/lat. They are projected to UTM zone 10N (EPSG:32610) in metres.

**Step 3. Choose the module.** The question is about behavioural modes in regular-interval
tracks, so the agent picks `movement-hmm` (a hidden Markov model). It rules out `movement-ssf`,
which answers habitat *selection*, not behavioural switching.

**Step 4. Plan.** The agent writes the plan and confirms it with the researcher:

- Model: 2-state HMM. Step lengths follow a gamma distribution and turning angles a von Mises
  distribution, each with its own values per state.
- State switching depends on time of day (cosine and sine of the hour) and on forest cover.
- Estimands: mean step and turning in each state, the activity budget (share of time in each
  state), and the time-of-day and forest effects on switching.
- Fit settings: 10 random starts (to guard against local optima), 2000 parameter draws for
  intervals.
- Adequate fit: most starts reach the same best likelihood, an invertible Hessian, and
  pseudo-residuals with no clear pattern.
- Two states is a choice. The agent tells the researcher that 2 vs 3 states should be judged
  from biology and residuals, not from AIC alone.

**Step 5. Planted truth (recovery at this design).** Before touching the real data, the agent
simulates 200 datasets with the same design (5 elk x 720 hourly fixes, 5% failed fixes,
+/- 2 min jitter, a forest effect of -0.025 / +0.025 per 1%). It then fits each one. Each fit
took about 45 s (median). All 200 converged.

Result: **15 of 16 estimands pass the gate.**

| Parameter | Truth | Mean estimate | Rel. bias | 95% CI coverage | Power | Pass |
|---|---|---|---|---|---|---|
| step_mean_m_state1 | 60 | 60.15 | +0.2% | 0.920 | | yes |
| step_mean_m_state2 | 450 | 450.87 | +0.2% | 0.935 | | yes |
| step_sd_m_state1 | 50 | 49.99 | -0.0% | 0.955 | | yes |
| step_sd_m_state2 | 300 | 299.74 | -0.1% | 0.960 | | yes |
| turn_concentration_state1 | 0.5 | 0.497 | -0.7% | 0.945 | | yes |
| turn_concentration_state2 | 2.5 | 2.506 | +0.2% | 0.920 | | yes |
| tpm_1to2_intercept | -1.683 | -1.688 | +0.3% | 0.960 | | yes |
| tpm_1to2_tod_cos | -0.9 | -0.895 | +0.6% | 0.955 | | yes |
| tpm_1to2_tod_sin | 0.4 | 0.402 | +0.5% | 0.915 | | yes |
| tpm_1to2_forest | -0.025 | -0.0260 | -4.1% | 0.930 | 0.995 | yes |
| tpm_2to1_intercept | -1.117 | -1.111 | +0.5% | 0.955 | | yes |
| tpm_2to1_tod_cos | 0.7 | 0.709 | +1.3% | 0.955 | | yes |
| tpm_2to1_tod_sin | -0.3 | -0.313 | -4.3% | **0.890** | | **no** |
| tpm_2to1_forest | 0.025 | 0.0251 | +0.5% | 0.965 | 1.000 | yes |
| activity_budget_state1 | 0.592 | 0.593 | +0.2% | 0.935 | | yes |
| activity_budget_state2 | 0.409 | 0.407 | -0.3% | 0.935 | | yes |

(The truth for the intercepts and activity budgets is the average over the 200 simulated
landscapes, because these values depend on mean forest cover. The coverage bound is 0.904,
which is 0.95 minus 3 Monte Carlo SE over 200 replicates.)

In one sentence for the researcher: *"At your design the method recovers the step lengths,
turning, activity budget and forest effects without bias (all within 4.3%), and it detects a
forest effect of your size in more than 99% of simulated studies. One time-of-day coefficient
has 95% intervals that are slightly too narrow (89% coverage), so time-of-day intervals should
be read as a little optimistic."* This matches a known limit of the module (VALIDATION.md, Known
limit 3).

**Step 6. Fit.** The agent writes `config.yaml` (file, column names, time zone, covariate, 2
states, 10 starts) and runs the fit. It took a few minutes on 4 cores.

**Step 7. Check.**

- Convergence: the optimiser converged. **10 of 10 starting values reached the same best
  log-likelihood (-24653.6)**, so a local optimum is unlikely.
- Hessian: the initial-state distribution (a nuisance parameter for the first fix of each elk)
  sits on a boundary. Its full Hessian is singular, so the intervals condition on its estimate.
  This does not affect any reported quantity.
- Goodness of fit (pseudo-residuals): step-length KS test **p = 0.93**, turning-angle KS test
  **p = 0.75**, lag-1 autocorrelation of step residuals **-0.02**. The QQ plots follow the line,
  and no autocorrelation lag exceeds the white-noise bounds except a single small spike (0.043 at
  lag 21). The model fits well.

**Step 8. Interpret.** The report (`results/movement-hmm/si/fit/movement-hmm_report.html`)
says, in its own words:

> "Encamped (short steps): mean step 60 m per 1 h (95% CI 57-63 m), turning concentration 0.50
> (0.43-0.57), i.e. little directional persistence (frequent sharp turns and reversals).
> Travelling (long steps): mean step 454 m per 1 h (95% CI 434-473 m), turning concentration 2.32
> (2.17-2.49), i.e. moderate directional persistence."

> "Activity budget (expected long-run share of time in each state, averaged over the 24-h cycle,
> at the average value of the habitat covariates): Encamped (short steps) 58% (95% CI 55-61%);
> Travelling (long steps) 42% (95% CI 39-45%)."

> "Time of day: the probability that an animal is in Travelling (long steps) is highest around
> 11:00 at 71% (95% CI 67-76%) and lowest around 23:00 at 15% (12-19%)."

> "Each 1-SD increase in forest_cover_pct (14.9 units) multiplies the odds of switching from
> Encamped to Travelling (versus staying) per time step by 0.67 (95% CI 0.58-0.79); switching back
> from Travelling to Encamped changes by a factor 1.48 (1.26-1.76) per SD. [...] the model-predicted
> long-run share of time in Encamped would be 41% (95% CI 36-47%) if forest_cover_pct were held at
> 32.5 (10th percentile of the fixes) and 73% (69-77%) if it were held at 72.8 (90th percentile).
> This is an association with the covariate measured at the animal's location, not proof that
> forest_cover_pct causes the behaviour."

The report also warns that missed fixes are treated as random. If fixes fail more often under
dense canopy, the forest association could be misstated.

Figures (in `results/movement-hmm/si/fit/figures/`, PDF and PNG):

| Figure | Caption (short) |
|---|---|
| `state_distributions` | Step-length (A) and turning-angle (B) distributions of each state, weighted by the activity budget, over histograms of the observed data. |
| `transition_tod` | Hourly probability of switching between states against hour of day (local time). |
| `stationary_tod` | Probability of being in each state at each clock hour (periodically stationary distribution). |
| `transition_forest_cover_pct` | Hourly probability of switching against forest cover, averaged over the day. |
| `stationary_forest_cover_pct` | **Key figure.** Long-run share of time in each state if forest cover were held at each value. |
| `decoded_track` | Each elk's track coloured by its most likely state (Viterbi). |
| `pseudo_residuals` | QQ plots of step and angle pseudo-residuals, and the autocorrelation of step residuals. |
| `coefficients` | Transition coefficients (log-odds of switching vs staying) with 95% and 50% intervals. |

One presentation change was made in the config's `figures:` block: the five one-panel figures
are drawn at 1.5-column width (`per_figure: {width: onehalf, aspect: 0.6}`). At single-column
width the state legend took about 40% of the figure. A bottom legend was tried first, but it
clipped axis titles, so it was not kept.

**Step 9. Gate verdict: OPEN.**

- Passed: the model converged from every start, the residual checks were clean, and the run is
  reproducible. The module was cross-checked against moveHMM and a hand-written likelihood
  (agreement within 0.014 SE), and a negative control (2-state model on 3-state data) fails the
  gate as it should (VALIDATION.md sections 5-6).
- Not passed: one of 16 estimands, `tpm_2to1_tod_sin`, had 0.890 interval coverage at this
  design (bound 0.904). Its bias was small (-4.3%).
- What this means: the main conclusions stand. These are the two states, the activity budget
  (58%/42%) and the forest association. The time-of-day *shape* is well estimated (bias of
  every time-of-day coefficient is at most 4.3%), but intervals on the time-of-day curves may be
  slightly too narrow. Report the 11:00 peak, not its exact interval width.

**Step 10. Deliver.** Summary for the researcher:

> Your elk switched between two movement modes. In the "encamped" mode they moved about 60 m per
> hour (95% CI 57-63 m) with frequent turns. In the "travelling" mode they moved about 454 m per
> hour (434-473 m) in fairly straight lines. They spent about 58% of their time encamped (55-61%).
> Travelling peaked around 11:00 (71% of the time, 67-76%) and was least common around 23:00 (15%,
> 12-19%). Elk in more forested places were less likely to start travelling and more likely to stop:
> the predicted share of time encamped rises from 41% (36-47%) at 32% forest to 73% (69-77%) at 73%
> forest. This is an association, not proof that forest causes the behaviour, and the labels
> "encamped" and "travelling" describe movement, not confirmed activities.

The methods paragraph for the paper is in the report's Methods section.

## (e) Estimates versus the known truth

| Parameter | Truth | Estimate | 95% CI | Inside CI? |
|---|---|---|---|---|
| step_mean_m_state1 | 60 | 59.75 | 56.77 to 62.79 | yes |
| step_mean_m_state2 | 450 | 453.73 | 434.33 to 472.72 | yes |
| step_sd_m_state1 | 50 | 50.35 | 47.15 to 53.56 | yes |
| step_sd_m_state2 | 300 | 293.62 | 277.96 to 309.84 | yes |
| turn_concentration_state1 | 0.5 | 0.495 | 0.428 to 0.570 | yes |
| turn_concentration_state2 | 2.5 | 2.325 | 2.166 to 2.494 | **no** (just above) |
| tpm_1to2_intercept | -1.667 | -1.633 | -1.782 to -1.477 | yes |
| tpm_1to2_tod_cos | -0.9 | -0.921 | -1.135 to -0.696 | yes |
| tpm_1to2_tod_sin | 0.4 | 0.480 | 0.266 to 0.702 | yes |
| tpm_1to2_forest (per 1%) | -0.025 | -0.0265 | -0.0365 to -0.0160 | yes |
| tpm_2to1_intercept | -1.133 | -1.131 | -1.302 to -0.959 | yes |
| tpm_2to1_tod_cos | 0.7 | 0.680 | 0.430 to 0.940 | yes |
| tpm_2to1_tod_sin | -0.3 | -0.413 | -0.631 to -0.204 | yes |
| tpm_2to1_forest (per 1%) | 0.025 | 0.0264 | 0.0154 to 0.0381 | yes |
| activity_budget_state1 | 0.586 | 0.582 | 0.549 to 0.611 | yes |
| activity_budget_state2 | 0.414 | 0.418 | 0.389 to 0.451 | yes |

15 of 16 intervals contain the truth. One miss in 16 is what 95% intervals should give.
The recovery study shows that the travelling turning concentration is unbiased (+0.2%,
coverage 0.92), so this miss is chance, not a flaw.

## (f) Limits from VALIDATION.md that matter at this design

1. **Number of states.** Two states is the analyst's choice. A model with too few states still
   converges cleanly but merges behaviours and biases everything. Compare 2 vs 3 states using
   residuals and biology.
2. **Transition-coefficient intervals are slightly narrow.** With no real effect, intervals excluded
   zero in about 7% of simulated studies, not 5%. At this design one time-of-day coefficient
   covered 0.890. Treat effects whose interval only just excludes zero as weak evidence. The forest
   effects here are far from zero (odds ratios 0.67 and 1.48 per SD), so this does not threaten
   them.
3. **Missed fixes are assumed random.** If GPS fixes fail more often under canopy or while elk
   are bedded, the encamped share and the forest association can be misstated.
4. **Pooled animals.** All five elk share one set of parameters. There are no individual random
   effects.
5. **Association, not cause.** Forest effects are associations with the conditions at the fix.
6. **State labels are interpretations.** The model sees only step lengths and turns. "Encamped"
   is not proof of resting or foraging.

## (g) Commands to reproduce

Run from `.`. On this cluster R runs through `scripts/ejob`.

```bash
# 0. (optional) regenerate the simulated example data and its truth (seed 2024)
scripts/ejob -c 1 -m 4 -t 0:15:00 -- scripts/Rscript dev/movement-hmm/make_example.R

# 1. quick timing check (4 replicates)
scripts/ejob -c 4 -m 16 -t 0:30:00 -l logs/movement-hmm_si_selftest.log -- \
  scripts/Rscript core/run.R selftest movement-hmm --config modules/movement-hmm/example/recovery.yaml \
  --nsim 4 --cores 4 --out dev/movement-hmm/si/selftest

# 2. planted-truth recovery at the example design (200 replicates)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/movement-hmm_si_recover.log -- \
  scripts/Rscript core/run.R recover movement-hmm --config modules/movement-hmm/example/recovery.yaml \
  --nsim 200 --cores 16 --out results/movement-hmm/si/recovery

# 3. fit the example data (report + figures)
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/movement-hmm_si_fit.log -- \
  scripts/Rscript core/run.R fit movement-hmm --config modules/movement-hmm/example/config.yaml \
  --out results/movement-hmm/si/fit

# 4. redraw figures after editing the config's figures: block (no refit)
scripts/ejob -c 2 -m 8 -t 0:30:00 -l logs/movement-hmm_si_figs.log -- \
  scripts/Rscript core/run.R figures movement-hmm --fit results/movement-hmm/si/fit/movement-hmm_fit.rds \
  --config modules/movement-hmm/example/config.yaml
```

Outputs:

- `results/movement-hmm/si/recovery/movement-hmm_recover_report.html`
- `results/movement-hmm/si/fit/movement-hmm_report.html`
- `results/movement-hmm/si/fit/figures/gallery.html`
- `results/movement-hmm/si/si_summary.json`
