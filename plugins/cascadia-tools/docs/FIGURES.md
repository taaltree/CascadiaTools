# Publication figure suites

Every module implements `figures(fit, data, opts)`. It returns a named list of
ggplots, each tagged with `with_caption()`, built from the constructors in
`core/R/figures.R`. The core applies the user's styling and exports at
journal sizes (`core/run.R fit` does this automatically, and
`core/run.R figures` redraws from a saved fit without refitting).

## Rules for every figure

1. **Response scale and original units.** Effects are drawn as probability,
   count, density per km², survival and so on. X axes show covariates in
   original units: if the module scales covariates, store the centre and
   scale in `data` and use `unscale()`. A link-scale plot is optional and
   never the default.
2. **Uncertainty is always visible.** Draw bands or intervals at `opts$ci`
   (default 95%), plus an inner `opts$inner_ci` band (50%) for Bayesian
   fits. Captions state the interval type (credible or confidence) and
   level.
3. **Show the data when it's meaningful** (`opts$show_data`): raw points,
   a rug, or observed values. Don't plot raw 0/1 detections against
   occupancy probability; use a rug.
4. **Everything can be relabelled.** Pass every axis label, parameter name,
   species, state and group name through `fig_label(key, default, opts)`,
   so users can rename anything from config without touching code. Use
   descriptive defaults ("Forest cover (%)", not "x1").
5. **Colour follows the dataviz rules.**
   - Categorical colours come from Okabe-Ito in a fixed order, and colour
     follows the entity (a species keeps its colour across every figure).
     With more than 8 groups, use facets.
   - Sequential maps use viridis.
   - Signed quantities use the diverging scale centred on the null value.
   - No dual axes.
   - A greyscale option must still work: use linetype or shape as a second
     encoding when there are groups.
6. **Use constructors before custom code.** Use `fig_cond_effect`,
   `fig_cat_effect`, `fig_coef`, `fig_posterior`, `fig_trace`, `fig_map`,
   `fig_timeseries`, `fig_matrix`, `fig_ppc`, `fig_recovery`, `fig_power`
   and `fig_path` where they fit. Custom ggplots are fine for
   analysis-specific figures, but don't hard-code theme, palette or sizes:
   the core applies them.
7. **Multi-panel figures** use `compose_figures()` (patchwork) with A/B/C
   tags, e.g. all covariate effects on one process.
8. **Captions.** Each caption says what is plotted, what the band means,
   and what the points are, in a form that can be pasted into a
   manuscript and edited.
9. **Render and look.** Export every figure from the example data, open
   the PNGs and check for overlapping labels, clipped text, empty panels,
   wrong units and illegible sizes at single-column width. Fix what you
   see.

## Required suite per analysis type

Names in `code` are the figure ids users refer to in `include:`,
`exclude:` and `per_figure:`.

| Module | Required figures |
|---|---|
| glm | `effect_<term>` for every predictor (conditional effects; factors via `fig_cat_effect`; interactions as one curve per level of the moderator), `effects_panel` (all effects composed), `coefficients` (forest plot on the transformed scale: odds/rate ratios), `residuals` (DHARMa QQ + residual vs predicted), `random_effects` (caterpillar plot, if any) |
| occu-single | `occupancy_<cov>` and `detection_<cov>` (conditional effects), `coefficients` (both processes, grouped), `occupancy_map` (if coordinates given), `detection_by_survey` (p over surveys), `gof` (MB c-hat bootstrap via `fig_ppc`); dynamic: `trajectory` (proportion occupied per season with CI), `col_ext_<cov>` |
| occu-rota | `conditional_occupancy` (P(sp B | A present/absent) vs covariate, one curve per state of A), `marginal_occupancy_<cov>`, `interaction_<cov>` (natural parameter f12 vs covariate), `co_occurrence_matrix`, `detection`, `coefficients` |
| occu-waddle | `conditional_occupancy` (psiBA vs psiBa with CI, and vs covariate), `sif` (species interaction factor with null line at 1), `detection_interaction` (rBA vs rBa), `coefficients` |
| occu-spatial | `occupancy_map` (mean and SD maps side by side), `spatial_effect_map`, `effect_<cov>`, `community_effects` (species-level effects with the community mean, forest plot), `richness_map`, `range` (posterior of the effective spatial range) |
| abund-interactions | `subordinate_vs_dominant` (subordinate occupancy vs dominant abundance, full posterior band), `interaction_by_habitat` (one curve per habitat level), `dominant_abundance_<cov>`, `dominant_abundance_map`, `coefficients` |
| nmixture | `abundance_<cov>`, `detection_<cov>`, `abundance_by_site` (caterpillar), `total_abundance` (posterior or profile), `gof` (PPC), `K_sensitivity`; dynamic: `trajectory` |
| scr | `density_surface` (map + traps + detected activity centres optional), `density_<cov>`, `detection_function` (p vs distance with band and trap spacing marked), `abundance_posterior`, `captures_per_individual`, `engine_comparison` (estimates across engines, when run) |
| movement-hmm | `state_distributions` (step and turn densities per state), `transition_<cov>` (transition probabilities vs covariate), `stationary_<cov>` (stationary state probabilities vs covariate), `decoded_track` (track coloured by Viterbi state, one panel per animal or chosen ids), `pseudo_residuals` (QQ and ACF) |
| movement-ssf | `rss_<cov>` (log-RSS vs covariate, Avgar et al. 2017), `selection_coefficients` (population mean plus individual random slopes), `movement_kernel` (updated step and turn distributions), `used_available` (density of used vs available per covariate), `home_ranges` (AKDE polygons per animal with CI contours), `hr_area` (areas with CI) |
| causal | `dag` (`fig_path`, adjustment set highlighted), `effect` (marginal effect or contrast on the response scale), `adjustment_comparison` (estimates under correct vs wrong adjustment sets, the educational "what if"), `implied_independencies` (test results) |
| sem | `path_diagram` (`fig_path`, standardised coefficients, width proportional to effect size, colour = sign, non-significant paths dashed), `indirect_effects` (forest plot of indirect and total effects), `dsep` (tests of missing paths), `partial_<path>` (partial residual plots for key paths) |
| popgen-structure | `diversity` (Ho, He, AR per population with CIs), `fst_matrix` (`fig_matrix`), `ibd` (FST/(1-FST) vs log distance with fit), `pca` (PC1 vs PC2 coloured by population, with % variance on the axes), `dapc_assignment` (STRUCTURE-style bar plot of membership), `bic_k`, `sample_map` (if coordinates) |
| popgen-ne | `ne_estimates` (Ne with CIs by population and method, log scale, reference lines at 50 and 500), `ne_by_maf` (sensitivity to the allele-frequency cutoff), `relatedness_distribution` (with known-relationship reference values), `kin_network` (pairs above threshold), `design_precision` (CV of Ne vs individuals and loci) |
| fisheries-growth | `growth_curves` (length-at-age with fitted VB per group, band, data points), `parameter_comparison` (Linf, K, t0 by group), `residuals`, `length_weight` (log-log with fit and observed points), `condition` (relative weight by group) |
| fisheries-stock | `catch_curve` (log catch vs age, ages used highlighted), `stock_recruitment` (fitted S-R with band and points labelled by year), `recruitment_residuals` (time series), `biomass_trajectory`, `kobe` (B/BMSY vs F/FMSY with uncertainty cloud and quadrants), `surplus_production` |
| cjs / multistate / jolly-seber / scr-open | `survival_<cov or time>`, `detection_<time>`, `transitions` (matrix or alluvial for multistate), `abundance_trajectory` (JS / open SCR), `recruitment`, `gof` (R2ucare / PPC) |
| distance | `detection_function` (histogram of distances with fitted function), `qq_gof`, `density_<cov>`, `density_by_stratum`, `density_map` (DSM / HDS) |
| ipm | `population_trajectory` (N_t with counts overlaid), `vital_rates` (time series of survival, fecundity, immigration), `ltre` (contributions to variance in lambda), `ipm_vs_separate` (precision comparison) |
| sdm-presence-only | `intensity_map` (relative intensity + presence points), `effect_<cov>`, `bias_effect`, `evaluation` (block-CV AUC / Boyce), `response_curves_panel` |
| community | `ordination` (biplot with sites coloured by group, environmental vectors, optional species labels via ggrepel), `variance_partition` (bars or Euler), `beta_components` (turnover vs nestedness), `dispersion` (betadisper distances by group), `indicator_species` |
| jsdm | `species_responses` (heatmap of Beta with support), `residual_associations` (`fig_matrix` of Omega with support mask), `variance_partitioning` (stacked bar per species), `trait_effects`, `richness_<cov>` |
| timeseries | `trajectory` (latent N with observations), `trend` (% change per year), `density_dependence`, `growth_<cov>`, `forecast` (with widening band), `power_years` (from power analysis) |
| pva | `quasi_extinction` (cumulative probability vs time with band), `trajectories` (fan chart of simulated N), `elasticity`, `scenarios` (extinction risk by management scenario), `lambda_distribution` |

Every module also gets these from the core: `recovery` (`fig_recovery`),
`power_curve` (`fig_power`), `trace` and `posterior` for Bayesian engines
(`fig_trace`, `fig_posterior`, with truth lines when fitting simulated
data).

## User customisation (config `figures:` block)

```yaml
figures:
  journal: esa            # generic | esa | besjournals | wiley_jwm | cjfas | plos | nature | slide
  width: single           # single | onehalf | double | number in mm
  aspect: 0.75            # height/width when height isn't given
  formats: [pdf, png, tiff]   # pdf (vector, fonts embedded), svg, eps, png, tiff (LZW)
  dpi: 600
  base_size: 9            # pt; default from journal preset
  font: "Helvetica"       # any installed font family
  style: classic          # classic | minimal
  palette: okabe_ito      # or a list of hex colours, or greyscale: true
  greyscale: false
  legend: right           # right | bottom | top | none
  titles: false           # journals put titles in captions
  ci: 0.95
  inner_ci: 0.5
  show_data: true
  include: [occupancy_forest, coefficients]   # only these
  exclude: [gof]
  labels:                 # rename anything: axes, params, species, states
    forest: "Forest cover (%)"
    psi_beta_forest: "Forest cover"
    sp1: "Pine marten"
  units:                  # back-transform scaled covariates for axes
    forest: {center: 45.2, scale: 18.1}
  per_figure:
    occupancy_forest: {width: double, aspect: 0.5, axis: {x: "Canopy cover (%)"}, ylim: [0, 1]}
```
