---
name: eco-interpretation
description: Rules for translating model output into biological statements for researchers and manuscripts. Use whenever writing results, report text or answering "what does this mean?".
---

# eco-interpretation

1. **Response scale first.** Report effects as changes in the quantity the
   biologist measures: probability of occupancy, expected count, density
   per km2, survival, selection strength. Coefficients on the logit/log
   scale go in the table, not the sentence.
2. **Meaningful units.** Per 1 SD of the covariate *and* per a biologically
   meaningful change (10% canopy, 100 m elevation, 1 km from road).
   Odds ratios for logit models; rate ratios (% change) for log links:
   exp(beta) - 1.
3. **Conditional effects plots by default**: response vs covariate, other
   covariates at their means (or reference levels), 95% band, rug of
   observed covariate values; never extrapolate beyond observed range
   without saying so.
4. **Uncertainty always.** Estimate (95% CI/CrI). Avoid "significant";
   describe the size and direction and how certain we are. For Bayesian
   fits, P(effect > 0) is acceptable as a supplement.
5. **Causal language** ("affects", "reduces", "drives") only when the
   causal module's adjustment set was used and its assumptions are stated
   or the design is experimental. Otherwise "is associated with".
6. **Detection vs ecology.** In hierarchical models, say which process an
   effect is in ("cameras near trails detected bobcats more often" is a
   detection effect, not habitat use).
7. **Species interactions** (Rota, Waddle, Twining): co-occurrence patterns
   are not proof of competition or predation; shared unmodelled habitat
   preferences produce the same signal.
8. **Abundance from N-mixture models** is sensitive to assumptions; present
   relative trends with more confidence than absolute N.
9. **State limits** discovered by planted-truth recovery.
