# Larger planted-truth runs for DONE examples near the gate limits

Requested by Taal Levi, 2026-11-05. Written before any of these runs started.
Why: the occu-single forest slope passed at +9.6% with 200 replicates, then failed at +10.8%
with 5,000. Other DONE examples have reported rows near the limit from 100-200 replicates.

## Scope
DONE examples with a passing row at |relative bias| >= 7.5%: cjs, jolly-seber, multistate, occu-rota,
occu-spatial, scr-search, ipm, popgen-structure. Every row in each config is rerun and read, not just
the near-limit rows. Same module, config and settings as the published run, with new seeds.

| Module | Replicates | Replicates in published run |
|---|---|---|
| cjs | 5,000 | 200 |
| occu-rota, module recovery | 5,000 | 200 |
| occu-rota, full example model | 5,000 | 400 |
| jolly-seber | 5,000 | 200 |
| popgen-structure | 2,000 | see summary |
| occu-spatial | 3,000 | 100 |
| multistate | 2,000 | 100 |
| scr-search | 2,000 | 100 |
| ipm | 1,000 | 100 |

## Decision rule
- The gate is unchanged (core/R/recovery.R): |bias| <= max(10% of truth, 3 MC SE), coverage >=
  max(0.80, 0.95 - 3 MC SE), and at least 90% of all attempted replicates converged.
- Rows the published example already fails and doesn't interpret stay not interpreted. Rows the example
  already treats as approximate truth (occu-rota coyote first-order terms) stay out of the verdict.
- Any other row that fails at the larger size makes the example OPEN on that estimand, with the numbers.
- All rows pass: the example stays DONE, and the record is updated with the larger run.
- Near-zero truths: when |truth| < 0.1 x empirical SD of the estimates, relative bias has no useful
  meaning, and a tolerance of 10% of truth becomes impossible as replicates grow. The gate result is still
  reported. Such a row doesn't flip a verdict by itself if |bias| <= 0.1 x empirical SD, which is the rule
  the gate already applies when truth is exactly 0. These rows are listed separately, for review.
- No rerun with new seeds after seeing a result. If a run is lost to a crash or time limit, its completed
  batches are reported as they are.

## Results (added after the runs, 2026-11-05; nothing above was changed)
No batches were lost. The full tables are at the end of each module's VALIDATION.md.

| Example | Result | Verdict |
|---|---|---|
| occu-single (5,000, run earlier under the same rule) | forest +10.8%, elevation 11.1% too strong | DONE to OPEN |
| occu-rota | gray fox days_active slope +10.2% (Monte Carlo error 0.9 points) | DONE to OPEN |
| popgen-structure | global FST coverage 0.922, one pairwise 0.934 (bound 0.935) | DONE to OPEN |
| scr-search | density coverage 0.930 (bound 0.935) | DONE to OPEN |
| cjs | all rows pass | stays DONE |
| jolly-seber | only the 3 rows already excluded fail | stays DONE |
| multistate | all rows pass | stays DONE |
| occu-spatial | only the 2 rows already excluded fail; canopy +9.5% | stays DONE |
| ipm | the 6 rows already excluded fail, plus ltre_adult_survival | stays DONE (judgement call below) |

Near-zero rule: used for popgen-structure Fis_global, which also passes the gate as is.

Judgement call: ipm ltre_adult_survival wasn't on the published list of failing rows. It stays out of the verdict
because the published limits already say LTRE contribution sizes are ranking-only. A stricter reading would make ipm OPEN.

Unresolved: the occu-rota row's Monte Carlo interval (about 8.5 to 11.9%) includes the limit. The gate result stands.
The other 11 DONE examples weren't rerun.
