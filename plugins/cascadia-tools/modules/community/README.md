---
id: community
title: Community composition - ordination, RDA and variation partitioning, PERMANOVA with dispersion test, beta-diversity partitioning, indicator species
when_to_use: A site x species table (counts, cover, biomass or presence/absence) from many sites plus a site table of environmental variables and/or groups (habitat types, treatments, regions, years); questions like "how much of the variation in species composition is related to elevation vs canopy cover?", "do burned and unburned plots hold different communities?", "is beta diversity driven by species replacement or species loss?", "which species characterise each habitat?", or "draw me an ordination".
not_for: Individual species' responses with uncertainty from imperfect detection (use occu-single / occu-spatial / nmixture, or a joint species distribution model, jsdm); species interactions (occu-rota, occu-waddle); fewer than about 10 sites or 3 species; time series of one community (timeseries); phylogenetic or functional beta diversity (betapart can do it, not wired here); causal claims about what drives composition (use causal / sem with stated assumptions).
data_needed: species CSV, wide (one row per site, one column per species) or long (site, species, count); site CSV with the same site ids, numeric environmental variables and an optional grouping column
engine: hand-written RDA / variation partitioning (adjusted R2), PERMANOVA sums of squares and Baselga (2010) pairwise components with delete-one-site jackknife intervals (identical to vegan::varpart, RsquareAdj, adonis2 and betapart::beta.pair); vegan permutation tests, NMDS, DCA, CCA, db-RDA, envfit, betadisper; betapart (multi-site); indicspecies (IndVal.g)
---

# Community composition

## Question it answers

* **How much of the among-site variation in species composition is
  associated with each environmental gradient?** Redundancy analysis (RDA) of
  Hellinger-transformed data, reported as adjusted R2 with a 95% jackknife
  interval and a permutation test, and **variation partitioning** into the
  part unique to each variable (set), the shared part, and the unexplained
  remainder.
* **Do groups of sites differ in composition?** PERMANOVA (Bray-Curtis by
  default) with adjusted R2 and interval, Holm-adjusted pairwise tests, the
  group effect conditional on the environment, and **always** the
  test of homogeneity of multivariate dispersion (PERMDISP / betadisper),
  because PERMANOVA responds to differences in spread as well as in average
  composition (Warton et al. 2012; Anderson & Walsh 2013).
* **Is beta diversity turnover or nestedness?** Baselga's (2010) partition
  of Sorensen dissimilarity into turnover (Simpson, species replacement) and
  nestedness-resultant components, as the mean over all pairs of sites (with
  intervals) and as the multi-site value.
* **Which species characterise each group?** IndVal.g indicator values with
  permutation p-values and false-discovery-rate correction.
* **Pictures**: NMDS (stress, Shepard diagram, fitted environmental vectors),
  RDA triplot, variance-partition bars, beta components and distance decay,
  dispersion plot, indicator species, and per-species response curves on the
  count/probability scale along each gradient. DCA axis length and CCA/db-RDA
  are reported to check whether linear RDA is adequate.

## Model

Let Y be the n x S site-by-species matrix.

* **Transformation.** Hellinger: y'_ij = sqrt(y_ij / y_i+) (Legendre &
  Gallagher 2001). Euclidean distance between Hellinger rows is the Hellinger
  distance, which does not treat shared absences as similarity, so linear
  RDA/PCA becomes appropriate for community data. Options: `chord`, `log`
  (log1p), `none`.
* **RDA and adjusted R2.** Multivariate least-squares regression of the
  centred Y' on the explanatory matrix X (p columns):
  R2 = ||fitted||^2 / ||Y'_c||^2 (the constrained share of total inertia).
  The adjusted R2 = 1 - (1 - R2)(n - 1)/(n - p - 1) (Ezekiel) is an unbiased
  estimator of the population share (Peres-Neto et al. 2006); unadjusted R2
  grows with every added variable even when nothing is related.
* **Variation partitioning** (Borcard et al. 1992; Peres-Neto et al. 2006)
  for sets X1..Xk: marginal fraction of set j = R2adj(Xj); unique fraction
  = R2adj(all) - R2adj(all but j) (tested by partial RDA, permutations
  of residuals under the reduced model); shared = R2adj(all) - sum of unique
  fractions (with two sets, fraction [b] = R2adj(X1) + R2adj(X2) - R2adj(X1,X2));
  it is not testable and can be negative.
* **PERMANOVA** (Anderson 2001; McArdle & Anderson 2001) on a dissimilarity
  d (Bray-Curtis default): SS_T = (1/n) sum_{i<j} d_ij^2,
  SS_W = sum_g (1/n_g) sum_{i<j in g} d_ij^2, pseudo-F =
  [(SS_T - SS_W)/(G - 1)] / [SS_W/(n - G)], p by permuting sites
  (adonis2). Effect size R2 = 1 - SS_W/SS_T, adjusted with G - 1 df (exactly
  the permutation expectation of R2 under no effect).
* **PERMDISP** (Anderson 2006): distance of each site to its group centroid
  in the principal-coordinate space of d, then ANOVA F with a permutation
  test (vegan::betadisper + permutest).
* **Beta partition** (Baselga 2010): for a site pair with a shared species,
  b and c species unique to each site, beta_sor = (b + c)/(2a + b + c),
  beta_sim = min(b, c)/(a + min(b, c)) (turnover),
  beta_sne = beta_sor - beta_sim (nestedness-resultant). Reported as means
  over all site pairs (turnover share = mean beta_sim / mean beta_sor) and
  as multi-site values (betapart::beta.multi).
* **Indicator value** (Dufrene & Legendre 1997; De Caceres & Legendre 2009):
  IndVal.g = sqrt(A x B), A = group-equalised specificity, B = fidelity;
  permutation test per species (indicspecies::multipatt, single groups);
  Benjamini-Hochberg correction across species.
* **Uncertainty.** Adjusted R2 values, PERMANOVA R2 and pairwise beta means
  get delete-one-site jackknife standard errors and t-intervals
  (theta_hat +/- t_{n-1} SE). Sites are the independent units. Tests are
  permutation tests (999 by default).
* **Species response curves**: per species, quasi-Poisson (counts) or
  quasi-binomial (presence/absence) GLM quadratic in the gradient (the
  Gaussian response model, ter Braak & Looman 1986), other numeric gradients
  at their means; the six species most associated with each gradient are
  drawn.

## Assumptions (and what breaks if they fail)

* **Sites are independent replicates.** Spatially clustered sites (several
  plots per stand) inflate the evidence; analyse at the cluster level or
  restrict permutations to blocks (not wired; see Pitfalls).
* **Equal sampling effort per site.** Dissimilarities and beta diversity
  describe communities *as sampled*. Unequal effort (different plot sizes,
  survey durations, number of visits) changes richness and abundance totals
  and appears as composition differences and nestedness. Standardise effort
  first (rarefy, or use relative abundances via Hellinger, which removes
  total abundance but not detection of rare species).
* **No correction for imperfect detection.** Species missed where present
  inflate turnover. If detection is a concern and you have repeat visits, use
  a multi-species occupancy model (occu-spatial).
* **Linear RDA of raw abundances assumes roughly monotone responses.** With
  long gradients (DCA axis 1 > ~4 SD) species respond unimodally; the
  Hellinger/chord transformations used here make RDA appropriate for such data
  (Legendre & Gallagher 2001), and CCA and db-RDA are reported alongside as a
  check. Validated: recovery passes with narrow niches (turnover share 0.92).
* **PERMANOVA tests location and dispersion together.** In unbalanced
  designs a more variable small group produces false "differences in
  composition" (validated below). Always read the dispersion test.
* **Association, not causation.** Fractions of explained variation are
  correlational.

## Data format

Species file, wide (`format: wide`, default):

```
station,PSFL,PAWR,HEWA,CBCH,...
CAS-01,0,0,0,0,...
CAS-02,1,0,0,0,...
```

or long (`format: long`; missing site x species combinations are zero,
repeated rows are summed):

```
site,species,count
CAS-01,AMRO,3
CAS-02,PSFL,1
```

Site file:

```
station,stand_type,elevation_m,canopy_cover_pct
CAS-01,old_growth,1540,55
```

Values can be counts, percent cover, biomass or 0/1 (detected automatically).
Sites with no species must be removed (Hellinger, Bray-Curtis and Sorensen
are undefined for them); `check_data` lists them. Species never recorded are
dropped automatically.

## Config

```yaml
data:
  species: counts.csv        # required
  format: wide               # wide | long
  site_col: station          # site id column (default: first column)
  species_col: species       # long format only
  count_col: count           # long format only
  drop_cols: [notes]         # non-species columns to ignore in a wide file
  env: sites.csv             # site table (optional, but needed for RDA/PERMANOVA)
  env_site_col: station      # default = site_col
  env_vars: [elevation_m, canopy_cover_pct]   # default: all numeric columns
  group: stand_type          # optional grouping factor for PERMANOVA/indicators
  partition:                 # optional variable sets for variation partitioning
    topography: [elevation_m, slope]
    stand: [canopy_cover_pct, basal_area]
  data_type: auto            # auto | abundance | presence
  min_sites: 0               # drop species found at fewer sites (e.g. 2)
  labels: {elevation_m: "Elevation (m)"}
settings:
  transform: hellinger       # hellinger | chord | log | none (RDA)
  distance: bray             # any vegan::vegdist method (PERMANOVA, NMDS, db-RDA, PERMDISP)
  binary: false              # true = presence/absence analysis of abundance data
  poly: 1                    # 2 adds quadratic terms for numeric variables in RDA
  n_perm: 999
  nmds_k: 2
  nmds_trymax: 50
  dispersion_type: centroid  # or median (spatial median)
  ci: 0.95
  dca, cca, dbrda, indicator, pairwise, multisite, species_curves: true
figures:                     # see skills/eco-figures; n_species_labels: 12
```

## Outputs

Estimates table (`param, estimate, sd, lower, upper, p_value`):

| param | meaning |
|---|---|
| `R2adj_<set>` | share of composition variation explained by the set alone (RDA, adjusted R2) |
| `R2adj_unique_<set>` | share explained by the set after the other sets (partial RDA) |
| `R2adj_shared` | share explained jointly (correlated variables), untestable |
| `R2adj_env_total` | share explained by all sets together |
| `R2adj_<group>_permanova` | PERMANOVA adjusted R2 for the group factor; p from adonis2 |
| `dispersion_ratio_<group>` | largest/smallest group mean distance-to-centroid; p from PERMDISP |
| `beta_sor_pair_mean`, `beta_sim_pair_mean`, `beta_sne_pair_mean` | mean pairwise Sorensen, turnover, nestedness |
| `turnover_share` | mean beta_sim / mean beta_sor |
| `beta_SOR/SIM/SNE_multisite` | Baselga multi-site values (depend on the number of sites) |
| `nmds_stress`, `dca_axis1_length_sd` | ordination quality; gradient length |
| `cca_R2adj_env_total`, `dbrda_R2adj_env_total` | CCA (permutation-adjusted) and db-RDA adjusted R2 |

The fit object also holds pairwise PERMANOVA, the PERMANOVA conditional on
the environment, the indicator table, vegan's varpart, RDA axis tests and
envfit. Figures: `ordination`, `nmds_shepard`, `rda_triplot`,
`variance_partition`, `effect_sizes`, `beta_components`, `distance_decay`,
`dispersion`, `indicator_species`, `species_response_<variable>`.

## Simulation and power

`simulate()` draws a species pool with Gaussian niches on two gradients
(`tol1`, `tol2`: niche widths as a fraction of the sampled gradient; smaller
= stronger gradient effect), peak abundances (`abund_max`, `abund_sd`),
species-specific group effects on log abundance (`group_effect_sd`, 0 = no
location difference), mean-preserving site-level noise (`site_noise_sd`), and
extra noise in the last group (`dispersion_mult` > 1 = dispersion-only
difference). Design knobs: `n_sites`, `n_groups`, `group_sizes` (unbalanced
designs), `n_gradients`, `env_cor` (gradient correlation -> shared fraction),
`effort` (multiplies expected counts: lower = more undersampling),
`data_type` (abundance | presence). Truth for every estimand is the
population value of the same statistic from 40,000 extra sites / 60,000 site
pairs per group combination (see VALIDATION.md, "What truth means").

Useful power runs: `--vary n_sites=20,30,45,60,90` (precision of R2 and
detecting a group effect); `--what params --vary group_effect_sd=0,0.2,0.4`.
Power from the core is the share of jackknife intervals excluding 0; the
PERMANOVA permutation test itself has slightly higher power (VALIDATION.md).

## Validation summary (details in VALIDATION.md)

* Default design (60 sites, nsim = 200): all 11 estimands pass the gate
  (coverage 0.91-0.975, |relative bias| <= 3%). Presence/absence and strong
  gradient designs pass too.
* 24 sites, presence/absence, half effort: shared fraction and nestedness
  intervals under-cover (0.855, 0.885); treat them as descriptive at < ~30 sites.
* Permutation tests hold their level (PERMANOVA 3.7-5.2%, PERMDISP 4.7-6.8%,
  indicator any-FDR 1%). Use p-values, not jackknife CIs, to decide whether an
  effect exists.
* Dispersion-only differences in unbalanced designs bias PERMANOVA. With the
  small group 3x noisier (48/12), R2 is inflated by 36-42% and the gate fails
  on it, as designed. PERMDISP flags the difference in 92-100% of datasets.
  The same heteroscedasticity also narrows some jackknife intervals
  (coverage 0.84-0.86).
* Point estimates equal vegan::varpart, ade4::varipart, adonis2 and
  betapart::beta.pair to < 5e-16.

## Pitfalls

* Reading PERMANOVA without PERMDISP. A significant PERMANOVA with a
  significant dispersion test, especially with unequal group sizes, is not
  evidence of a shift in average composition.
* Unadjusted R2 (vegan `rda()` default summary) overstates explained
  variation, by about p/(n - 1) when nothing is related.
* Treating fractions as effect sizes for species, or the shared fraction as
  an interaction.
* Multi-site beta diversity depends on how many sites you have; compare
  pairwise means or resample to equal numbers of sites.
* NMDS pictures with stress > 0.2 mean little; base inference on tests.
* Nested sampling (plots within stands) needs restricted permutations
  (`permute::how(blocks=)`) - not exposed; ask for a project-specific script.
* Raw counts in PCA/RDA (no transformation) let a few abundant species dominate
  and treat double zeros as similarity.

## References

See `module$references` (printed by `core/run.R describe community`). Key
sources: Legendre & Gallagher 2001; Borcard et al. 1992; Peres-Neto et al.
2006; Anderson 2001, 2006; Anderson & Walsh 2013; Warton et al. 2012;
Baselga 2010; Baselga & Orme 2012; Dufrene & Legendre 1997; De Caceres &
Legendre 2009; ter Braak 1986; Legendre & Legendre 2012; vegan (Oksanen et al.).
