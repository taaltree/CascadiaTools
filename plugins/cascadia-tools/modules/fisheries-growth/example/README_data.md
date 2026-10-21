# Worked-example data: simulated walleye, autumn 2014 (`sim_walleye_2014.csv`)

**These data are SIMULATED.** No fish were caught. The file was generated with the
module's own simulator so that every true value is known. This lets a reader check that
the analysis recovers the truth. It is the dataset for Appendix S1 and the website's
"try an example" walkthrough.

The other dataset in this folder, `walleye_lake_erie.csv` (used by `config.yaml`), is
real (FSAdata::WalleyeErie2, Ohio DNR; educational use only) and has no known truth.

## The (simulated) study

- **Species / system:** walleye (*Sander vitreus*), a large lake (modelled on western Lake Erie).
- **Story:** in October 2014 a gill-net survey caught 1,000 walleye. Each fish was sexed,
  measured (total length, mm), weighed (g) and aged from otoliths.
- **Question:** do females and males differ in growth (how big they get, how fast they get
  there) and in condition (weight at a given length)?
- **Why it matters:** female walleye usually grow larger. Sex-specific growth affects
  size limits, egg production and stock assessment.

## Design

| Dimension | Value |
|---|---|
| Groups | 2 (female = reference, male) |
| Fish per group | 530 female, 470 male (1,000 total) |
| Ages (otolith, integer years) | 1 to 16; most fish are 1 to 4 years old |
| Age used in the model | otolith age + 0.5 (fish caught in October, about half a year after spring hatch) |
| Age structure | total mortality Z = 0.55 per year, from age 1 |
| Fish weighed | all 1,000 |
| Year | 2014 only |

Fish per otolith age (female / male): 1: 219/205, 2: 142/107, 3: 68/68, 4: 43/36,
5: 22/22, 6: 17/12, 7: 7/10, 8: 4/3, 9: 2/3, 10: 3/1, 12: 0/1, 13: 3/0, 15: 0/1, 16: 0/1.

## Columns

| Column | Meaning |
|---|---|
| `fish_id` | Unique fish tag, `SIM14-<F/M>-<number>` |
| `year` | Survey year (2014 for every fish) |
| `sex` | `female` or `male` (the comparison group) |
| `age` | Otolith age in whole years |
| `tl_mm` | Total length, mm (whole mm) |
| `wt_g` | Wet weight, g (0.1 g); none missing |

## TRUE parameter values (simulated)

Growth: von Bertalanffy, L = Linf (1 - exp(-K (age - t0))), multiplicative (lognormal) error.

| Parameter | Female | Male |
|---|---|---|
| Linf (mm) | 680 | 560 |
| K (per year) | 0.31 | 0.39 |
| t0 (years) | -1.0 | -1.1 |
| Age at 90% of Linf (years) | 6.43 | 4.80 |

| Shared / derived parameter | True value |
|---|---|
| Residual SD of log length | 0.05 |
| Male minus female: Linf | -120 mm |
| Male minus female: K | +0.08 per year |
| Male minus female: t0 | -0.10 years |
| Length-weight: W = a L^b, a | 1.2e-6 |
| Length-weight: b (common slope) | 3.34 |
| Residual SD of log weight | 0.09 |
| Weight ratio male / female at the same length | 1.02 (males 2% heavier) |
| Mean weight of a 500 mm female | 1,245.9 g |

## How the file was made

```
scripts/ejob -c 1 -m 4 -t 0:10:00 -l logs/fisheries-growth_si_make.log -- \
  scripts/Rscript dev/fisheries-growth/make_si_example.R
```

- Generator: `dev/fisheries-growth/make_si_example.R`, which calls `module$simulate()`
  (`modules/fisheries-growth/simulate.R`).
- Seed: **20141015**.
- All true values, the design and the seed are stored in `sim_walleye_2014_truth.json`.
  The module labels the groups A (female) and B (male); the JSON gives both naming schemes.
- The fit config for these data is `config_si.yaml`. The planted-truth config is `recovery.yaml`.
