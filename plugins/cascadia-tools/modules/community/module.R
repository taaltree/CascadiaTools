# community: multivariate community ecology (ordination, constrained
# ordination with variation partitioning, PERMANOVA + dispersion, Baselga
# beta-diversity partitioning, indicator species).
#
# Engine: hand-written least squares / sum-of-squares / pairwise estimators
# for every estimand (so delete-one-site jackknife intervals are exact and
# cheap), cross-checked against vegan::varpart, RsquareAdj, adonis2 and
# betapart::beta.pair; vegan for permutation tests, NMDS, DCA, CCA, db-RDA,
# envfit and betadisper; indicspecies for IndVal.g. simulate() builds
# communities from Gaussian niche models along two gradients with known group
# location and dispersion effects; truth is the population value of each
# statistic from large independent Monte Carlo samples (see simulate.R).

source(file.path(MODULE_DIR, "simulate.R"), local = TRUE)
source(file.path(MODULE_DIR, "engine.R"), local = TRUE)
source(file.path(MODULE_DIR, "io.R"), local = TRUE)
source(file.path(MODULE_DIR, "output.R"), local = TRUE)

module <- ecomodule(
  id = "community",
  title = "Community composition: ordination, RDA/variation partitioning, PERMANOVA, beta diversity, indicator species",
  description = paste(
    "How much of the variation in species composition among sites is associated with environmental gradients",
    "(and how much is unique to each gradient or shared), do groups of sites (habitats, treatments, regions)",
    "differ in average composition or only in how variable they are, is beta diversity driven by species",
    "turnover or by nestedness, and which species characterise each group? Site x species matrix",
    "(abundance or presence/absence) plus a site table of environmental variables and groups."),
  defaults = list(
    params = list(n_species = 30, abund_max = 4, abund_sd = 0.8, tol1 = 0.2, tol2 = 0.4,
                  group_effect_sd = 0.5, site_noise_sd = 0.5, dispersion_mult = 1),
    design = list(n_sites = 60, n_groups = 3, group_sizes = NULL, n_gradients = 2, env_cor = 0.4,
                  effort = 1, data_type = "abundance", truth_n = 40000, truth_pairs = 60000)),
  simulate = com_simulate,
  fit = com_fit,
  settings = com_default_settings(),
  fast_settings = merge_list(com_default_settings(), list(
    n_perm = 0, nmds = FALSE, dca = FALSE, cca = FALSE, dbrda = FALSE, indicator = FALSE,
    pairwise = FALSE, multisite = FALSE, species_curves = FALSE)),
  read_data = com_read_data,
  check_data = com_check_data,
  write_data = com_write_data,
  interpret = com_interpret,
  figures = com_figures,
  plots = function(fit, data) com_figures(fit, data, fig_opts()),
  methods_text = com_methods_text,
  power_params = c("R2adj_elevation", "R2adj_soil_moisture", "R2adj_unique_elevation",
                   "R2adj_unique_soil_moisture", "R2adj_habitat_permanova"),
  null_values = c(R2adj_elevation = 0, R2adj_soil_moisture = 0, R2adj_unique_elevation = 0,
                  R2adj_unique_soil_moisture = 0, R2adj_habitat_permanova = 0),
  engine = "R: hand-written RDA/variation partitioning, PERMANOVA sums of squares and Baselga components with delete-one-site jackknife CIs (match vegan/betapart); vegan (permutation tests, NMDS, DCA, CCA, db-RDA, envfit, betadisper); betapart; indicspecies",
  references = c(
    "Legendre P, Gallagher ED (2001) Ecologically meaningful transformations for ordination of species data. Oecologia 129:271-280. doi:10.1007/s004420100716",
    "Borcard D, Legendre P, Drapeau P (1992) Partialling out the spatial component of ecological variation. Ecology 73:1045-1055.",
    "Peres-Neto PR, Legendre P, Dray S, Borcard D (2006) Variation partitioning of species data matrices: estimation and comparison of fractions. Ecology 87:2614-2625. doi:10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2",
    "Anderson MJ (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecology 26:32-46.",
    "McArdle BH, Anderson MJ (2001) Fitting multivariate models to community data: a comment on distance-based redundancy analysis. Ecology 82:290-297.",
    "Anderson MJ (2006) Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62:245-253.",
    "Anderson MJ, Walsh DCI (2013) PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: what null hypothesis are you testing? Ecological Monographs 83:557-574. doi:10.1890/12-2010.1",
    "Warton DI, Wright ST, Wang Y (2012) Distance-based multivariate analyses confound location and dispersion effects. Methods in Ecology and Evolution 3:89-101.",
    "Baselga A (2010) Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143. doi:10.1111/j.1466-8238.2009.00490.x",
    "Baselga A, Orme CDL (2012) betapart: an R package for the study of beta diversity. Methods in Ecology and Evolution 3:808-812.",
    "Dufrene M, Legendre P (1997) Species assemblages and indicator species: the need for a flexible asymmetrical approach. Ecological Monographs 67:345-366.",
    "De Caceres M, Legendre P (2009) Associations between species and groups of sites: indices and statistical inference. Ecology 90:3566-3574. doi:10.1890/08-1823.1",
    "ter Braak CJF (1986) Canonical correspondence analysis: a new eigenvector technique for multivariate direct gradient analysis. Ecology 67:1167-1179.",
    "ter Braak CJF, Looman CWN (1986) Weighted averaging, logistic regression and the Gaussian response model. Vegetatio 65:3-11.",
    "Hill MO, Gauch HG (1980) Detrended correspondence analysis: an improved ordination technique. Vegetatio 42:47-58.",
    "Legendre P, Anderson MJ (1999) Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. Ecological Monographs 69:1-24.",
    "Minchin PR (1987) An evaluation of the relative robustness of techniques for ecological ordination. Vegetatio 69:89-107.",
    "Clarke KR (1993) Non-parametric multivariate analyses of changes in community structure. Australian Journal of Ecology 18:117-143.",
    "Legendre P, Legendre L (2012) Numerical Ecology, 3rd English edition. Elsevier, Amsterdam.",
    "Oksanen J et al. vegan: Community Ecology Package. R package version 2.7. https://CRAN.R-project.org/package=vegan")
)
