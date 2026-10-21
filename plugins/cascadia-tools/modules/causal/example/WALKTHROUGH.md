# Worked example: does riparian canopy cause more juvenile trout?

Module: `causal` (DAG + g-computation). The data are **simulated**, so we
know the right answer and can check the analysis against it.
See `README_data.md` for the full story and the true values.

## (a) The question, in the researcher's words

> "Reaches with more streamside canopy seem to have more young trout. Is that
> because shade itself helps the fish, or just because shady reaches are
> different in other ways, such as higher up, less farmland or narrower? If we
> restored 10 % more canopy, how many more trout should we expect? And how much
> of that comes through cooler water?"

## (b) The files she has

`stream_trout_reaches.csv`: 140 reaches, one electrofishing pass each.

```
"reach_id","stream","survey_date","canopy_pct","trout_per_reach","elevation_m","agriculture_pct","fine_sediment_pct","wetted_width_m","water_temp_c","periphyton_chla"
"R001","Mill Cr","2025-08-01",66,15,780,43.8,14.1,4.7,11.2,57.8
"R002","Beaver Cr","2025-07-04",43,3,678,39.2,25.7,8,17.7,63
"R003","Elk Cr","2025-08-15",60,3,947,34,17.3,6.8,15.2,32.9
"R004","Alder Cr","2025-07-27",85,5,674,33.1,15.6,4,10.4,31.3
```

She also has a sketch of what she believes causes what. That sketch becomes
the DAG in `config.yaml`. Each arrow means "directly causes". Each missing
arrow means "no direct effect".

## (c) What she types to Claude

> I surveyed 140 stream reaches last summer, one electrofishing pass each. The
> file is stream_trout_reaches.csv. I want to know whether riparian canopy cover
> actually causes more juvenile trout, and by how much for a 10 % increase in
> canopy, not just whether they're correlated. Elevation, farmland in the
> catchment and stream width all affect both canopy and trout. I think shade
> works partly by keeping the water cool. I also measured periphyton, which
> depends on both shade and how many trout graze. Please do a proper causal
> analysis and give me a report with figures I can use in a paper.

## (d) What the agent does, step by step (eco-analyze)

**1. Project folder.** Creates `projects/<name>/` and copies the CSV into
`data_raw/` unchanged.

**2. Data intake.** Describes the data back to her: 140 reaches, 11 columns, no
missing values. Canopy ranges from 27 to 97 % (mean 61). Trout range from 0 to
68 per reach (mean 6.3, variance 50). The counts are strongly overdispersed,
6.4 % are zeros, and one reach has 68 trout. Stream and date are labels only.

**3. Choose the module.** `causal`. The question asks for an *effect of an
intervention* from observational data, so a regression with "all the
covariates" is not enough. With one pass per reach, detection cannot be
modelled. The report says so: the count is treated as an index.

**4. Plan.** The agent writes the DAG with her (8 variables, 12 arrows). From the
DAG, dagitty finds two minimal adjustment sets for the total effect:
{elevation, fine sediment, width} (used) and {agriculture, elevation, width}.
- Water temperature is a **mediator**, so it is *not* adjusted for in the total effect.
- Periphyton is a **collider**, so it is never adjusted for.

The plan specifies:
- outcome model: negative binomial GLM;
- estimand: the average causal effect of +10 percentage points of canopy;
- second estimand: the controlled direct effect with water temperature held fixed.

**5. Planted truth (gate step 1).** The agent simulates 200 datasets at exactly
this design: 140 reaches, the same DAG, natural units and rounding, and true
effects matching the example. It fits each one and checks the answers.
**All six estimands pass.** Bias is at most 2.6 %, and coverage is 0.940-0.955
(bound 0.904). All 200 fits succeeded.

| estimand | truth | mean estimate | rel. bias | 95 % coverage | power | gate |
|---|---|---|---|---|---|---|
| total effect, trout per +10 % canopy | 1.954 | 1.989 | +1.8 % | 0.950 | 1.00 | PASS |
| total rate ratio | 1.293 | 1.295 | +0.2 % | 0.945 | 1.00 | PASS |
| total beta (log, per %) | 0.02571 | 0.02569 | -0.1 % | 0.945 | 1.00 | PASS |
| direct effect, trout per +10 % | 1.023 | 1.050 | +2.6 % | 0.940 | 0.51 | PASS |
| direct rate ratio | 1.154 | 1.156 | +0.2 % | 0.955 | 0.60 | PASS |
| direct beta (log, per %) | 0.01429 | 0.01424 | -0.3 % | 0.955 | 0.60 | PASS |

Power is the share of replicates whose 95 % interval excluded no effect. In
plain words: at 140 reaches the method gets the total effect right and almost
always detects it. The direct effect is estimated without bias, but only about
half of such studies would detect it.

**6. Fit.** `core/run.R fit causal` fits a negative binomial GLM on the 140
reaches (MASS::glm.nb). It then computes the average effect by g-computation:
every reach gets +10 % canopy and the predictions are averaged.

**7. Check.**
- *Convergence:* the GLM converged. There is no MCMC here, so R-hat and ESS do not apply.
- *Is the DAG consistent with the data?* The DAG implies 16 testable
  independencies. None is contradicted (Fisher's C = 31.4, df = 32, p = 0.50).
- *Outcome model (DHARMa):* uniformity p = 0.71 and outliers p = 1.00, but
  **dispersion p = 0.016 (flagged)**. The residual plot lies close to the 1:1
  line. With a correct negative binomial, this flag occurs in only 1.3 % of
  datasets, so this draw is unusually heavy-tailed (one reach has 68 trout).
  The agent reports the flag and re-fits with bootstrap intervals as a
  sensitivity check (next step).
- *Unmeasured confounding:* E-value 1.66. To make the 95 % interval include no
  effect, a hidden confounder would need to change both canopy and trout by a
  factor of 1.29.

**8. Interpret.** From the report (`results/causal/si/fit/causal_report.html`):

> "Average causal effect: increasing riparian canopy cover by 10 percentage
> points at every site would change the expected number of juvenile trout per
> 100 m reach by 1.18 (95% CI 0.258 to 2.10), from 6.30 to 7.48 on average. As a
> ratio (rate ratio), 1.19 (1.05 to 1.34): a 19% increase."

> "Controlled direct effect (holding August water temperature fixed): 0.566
> (95% CI -0.296 to 1.43) on the count scale. On the model scale about 50% of the
> total effect runs through the mediator(s)..."

> "Why the adjustment set matters (what-if table): the DAG-based estimate is
> 1.18. ... No adjustment (crude association) -> 1.78 (+51%); Omit elevation_m
> -> 2.14 (+81%); ... Omit wetted_width_m -> 0.530 (-55%); Also adjust for
> mediator water_temp_c -> 0.552 (-53%); Also adjust for periphyton_chla
> (collider) -> 1.59 (+34%)."

Bootstrap sensitivity (2,000 resamples, `ci_method: bootstrap`):
- total effect 1.18 trout (95 % CI 0.03 to 2.52);
- rate ratio 1.19 (1.005 to 1.385).

The direction still holds, but only just. The agent tells her to lead with
the rate ratio, give both interval methods, and avoid over-stating precision.

Figures (in `results/causal/si/fit/figures/`, PNG and PDF):

| id | caption (short) |
|---|---|
| `dag` | The causal diagram. Exposure orange, outcome blue, adjusted variables green; solid arrows are causal paths from canopy to trout. |
| `effect` | Causal dose-response: expected trout per reach if every reach had a given canopy (g-computation), 95 % band, rug of observed canopy. |
| `adjustment_comparison` | What-if: the same effect under 9 adjustment choices. Green = DAG set, blue = other valid set, orange = choices the DAG says are biased. |
| `implied_independencies` | The 16 partial correlations the DAG says are zero, with 95 % CIs; all consistent. |
| `residuals` | DHARMa simulated residuals vs uniform quantiles (dispersion flagged, p = 0.016). |

The most informative figure is `adjustment_comparison`. It shows at a glance
why the DAG matters. The crude estimate and the "omit elevation" estimate are
inflated by 51-81 %. Adjusting for the mediator or omitting width halves the
effect.

**9. Gate.** Every gate step passed:
- planted truth: PASS;
- fit: converged, no boundary estimates;
- model checks: DAG local tests consistent; the one DHARMa flag is reported and a bootstrap sensitivity fit was run;
- independent routes: dagitty vs a from-scratch adjustment criterion, and the module vs lavaan and hand g-computation (VALIDATION.md);
- negative controls: wrong adjustment sets fail the gate (VALIDATION.md);
- the run is reproducible.

**Verdict: DONE**, with the stated caveats.

**10. Deliver.** Summary for the researcher: "Under your DAG, adding 10
percentage points of canopy would raise juvenile trout by about 1.2 per 100 m
reach (95 % CI 0.3-2.1; bootstrap 0.03-2.5). That is a 19 % increase (rate
ratio 1.19, 1.05-1.34). About half of this may run through cooler water, but
the direct part (0.57, -0.30 to 1.43) is too uncertain to separate with 140
reaches. Without adjustment the effect would be overstated by about 50 %. The
result assumes your DAG has no hidden common cause stronger than an E-value of
1.29-1.66."

## (e) Estimates versus the known truth

| estimand | truth | estimate (95 % CI) | covered? |
|---|---|---|---|
| total effect (trout per +10 % canopy) | 1.954 | 1.18 (0.26 to 2.10) | yes |
| total rate ratio | 1.293 | 1.19 (1.05 to 1.34) | yes |
| total beta (log, per %) | 0.0257 | 0.0172 (0.0053 to 0.0291) | yes |
| direct effect (trout per +10 %) | 1.023 | 0.57 (-0.30 to 1.43) | yes |
| direct rate ratio | 1.154 | 1.09 (0.96 to 1.24) | yes |
| direct beta (log, per %) | 0.0143 | 0.0086 (-0.0039 to 0.0211) | yes |

All six intervals contain the truth. This draw is on the low side: the total
effect is about 1.3 SE below the truth. That is ordinary sampling variation;
in the 200 planted-truth replicates the mean estimate was 1.99.

## (f) Limits that matter at this design (from VALIDATION.md)

- **The DAG cannot be tested where it matters most.** The local tests catch
  missing arrows between measured variables (95 % power at 120 sites). They
  cannot detect an unmeasured confounder of canopy and trout. Use the E-value.
- **A moderate collider can pass on the ratio scale.** Adjusting for periphyton
  gave only +5 % bias on the rate ratio in validation, but +18 % on the count
  effect. Keep periphyton out because of the DAG, not because a test allows it.
- **Small samples:** the count-scale effect is biased upward (+18 %) at 40
  sites. At 140 reaches it passed (+1.8 %).
- **Controlled direct effect only.** The "about 50 % through water temperature"
  split is a rough log-scale calculation, not a natural indirect effect.
- **Linear, homogeneous effects.** No canopy x covariate interactions. Confounders
  are assumed to act linearly on the log scale.
- **Detection is not modelled.** With one pass, a change in capture efficiency
  with canopy would be mixed into the estimate.

## (g) Reproduce

```bash
cd .
# regenerate the example data (seed 2024) and the truth
scripts/ejob -c 1 -m 4 -t 0:20:00 -l logs/causal_make_examples.log -- scripts/Rscript dev/causal/make_examples.R
# planted-truth recovery at this design (200 replicates, ~2 min on 16 cores)
scripts/ejob -c 16 -m 32 -t 3:00:00 -l logs/causal_si_recover.log -- scripts/Rscript core/run.R recover causal \
  --config modules/causal/example/recovery.yaml --nsim 200 --cores 16 --out results/causal/si/recovery
# real-data fit, report and figures
scripts/ejob -c 4 -m 16 -t 2:00:00 -l logs/causal_si_fit.log -- scripts/Rscript core/run.R fit causal \
  --config modules/causal/example/config.yaml --out results/causal/si/fit
# bootstrap-interval sensitivity fit
scripts/ejob -c 4 -m 16 -t 1:00:00 -l logs/causal_si_boot.log -- scripts/Rscript core/run.R fit causal \
  --config dev/causal/si_boot/config_bootstrap.yaml --out dev/causal/si_boot/fit
```
