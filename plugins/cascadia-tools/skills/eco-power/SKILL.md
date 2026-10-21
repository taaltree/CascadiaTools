---
name: eco-power
description: Simulation-based power analysis and study design for any CascadiaTools module (how many sites, surveys, traps, animals, loci, years). Use when a researcher asks about sample size, power, survey design, or "can I detect an effect of X".
---

# eco-power — design a study by simulation

1. **Pin down the effect.** Which parameter, what size is biologically
   meaningful (ask: "what change would matter for management?"), on which
   scale. Convert to the module's parameterisation (e.g. a 20% decline in
   occupancy per SD of road density -> logit-scale slope) and show the
   researcher the conversion.
2. **Set realistic nuisance parameters**: detection probability, baseline
   occupancy/density, movement scale, number of loci... from literature or
   pilot data. Uncertain ones become extra scenarios.
3. **Choose the knob** to vary (n_sites, n_surveys, n_traps, n_years, n
   individuals, n_loci) and a grid of 4-6 values spanning feasible effort.
4. **Run** (SLURM):
   ```
   "${CLAUDE_PLUGIN_ROOT}/scripts/ejob" -c 32 -m 64 -t 12:00:00 -- "${CLAUDE_PLUGIN_ROOT}/scripts/Rscript" "${CLAUDE_PLUGIN_ROOT}/core/run.R" power <module> --vary n_sites=40,80,120,200 \
     --config power.yaml --nsim 200 --cores 32 --out <project>/power
   ```
   Time a `selftest` first and scale nsim to the budget. 200 reps give
   power +/- ~0.03 near 0.8.
5. **Read power only where the estimator works.** At each design value,
   check bias and coverage too; power from a biased estimator is
   meaningless. Small designs often fail coverage before they fail power.
6. **Also run effect = 0** at the chosen design: false-positive rate
   should be ~0.05.
7. **Report**: the curve, the smallest design reaching 80% power, the
   assumptions it rests on, and how sensitive the answer is to the
   nuisance-parameter scenarios. Trade-offs (more sites vs more surveys)
   deserve a two-knob comparison at fixed total effort.
