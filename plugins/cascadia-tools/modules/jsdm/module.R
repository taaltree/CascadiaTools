# jsdm: joint species distribution models with HMSC (Hmsc R package).
#   presence/absence (probit) or counts (lognormal Poisson); site covariates;
#   species traits (Gamma) and phylogeny (rho) informing species' responses;
#   latent factors for residual species associations at the sampling-unit
#   level (optionally spatial: NNGP / full GP).

for (f in c("sim.R", "fit.R", "io.R", "report.R")) sys.source(file.path(MODULE_DIR, "R", f), envir = environment())

module <- ecomodule(
  id = "jsdm",
  title = "Joint species distribution model (HMSC)",
  description = paste(
    "How do the species of a community respond to environmental gradients, do species traits (and phylogeny) explain",
    "which species respond how, which species co-occur more or less often than the environment predicts (residual",
    "associations), and how does expected species richness change along the gradients? Fits a Hierarchical Model of",
    "Species Communities (HMSC) to presence/absence or count data from many species at the same sampling units."),
  defaults = list(
    params = list(comm_mean_intercept = -0.3, comm_mean_temperature = 0.6, comm_mean_forest = 0.3,
                  trait_log_body_mass_on_intercept = 0.4, trait_log_body_mass_on_temperature = -0.5,
                  trait_log_body_mass_on_forest = 0,
                  species_sd = 0.5, loading_sd1 = 0.8, loading_sd2 = 0.6, residual_sd = 0.5),
    design = list(n_sites = 100, n_species = 12, n_factors = 2, response = "presence",
                  community_seed = 2024, community = "random", phylo_rho = 0, spatial_range_km = 0, extent_km = 20,
                  cov_correlation = 0, fit_spatial = FALSE)),
  simulate = simulate_jsdm,
  fit = fit_jsdm,
  # Latent-factor shrinkage prior a1 = 2, a2 = 3 (Durante 2017) instead of the
  # Hmsc default a1 = a2 = 50, which attenuated Beta and residual correlations
  # in simulation (VALIDATION.md section 0).
  settings = list(samples = 1000, thin = 50, transient = 25000, n_chains = 4, n_parallel = 4,
                  max_factors = 5, spatial = "none", converge_on = "core", seed = 1, full_output = TRUE,
                  factor_prior = list(a1 = 2, b1 = 1, a2 = 3, b2 = 1)),
  fast_settings = list(samples = 250, thin = 80, transient = 20000, n_chains = 3, n_parallel = 1,
                       max_factors = 5, spatial = "none", converge_on = "core", seed = 1, full_output = FALSE,
                       factor_prior = list(a1 = 2, b1 = 1, a2 = 3, b2 = 1)),
  read_data = read_jsdm,
  check_data = check_jsdm,
  interpret = interpret_jsdm,
  figures = figures_jsdm,
  plots = plots_jsdm,
  methods_text = methods_jsdm,
  # Gamma (hyperparameters of the trait model; the right target for "is there
  # an effect at all" / false-positive rate) and their population-averaged,
  # identified counterparts for the realised community (probit only).
  power_params = c("comm_mean_temperature", "comm_mean_forest",
                   "trait_log_body_mass_on_intercept", "trait_log_body_mass_on_temperature",
                   "trait_log_body_mass_on_forest",
                   "comm_popavg_temperature", "comm_popavg_forest",
                   "popavg_trait_log_body_mass_on_temperature", "popavg_trait_log_body_mass_on_forest"),
  engine = "Hmsc (Bayesian Gibbs MCMC; probit / lognormal Poisson; latent factors with multiplicative gamma shrinkage)",
  references = c(
    "Ovaskainen O, Tikhonov G, Norberg A, Blanchet FG, Duan L, Dunson D, Roslin T, Abrego N (2017) How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters 20:561-576. doi:10.1111/ele.12757",
    "Ovaskainen O, Abrego N (2020) Joint Species Distribution Modelling: With Applications in R. Cambridge University Press. doi:10.1017/9781108591720",
    "Tikhonov G, Opedal OH, Abrego N, Lehikoinen A, de Jonge MMJ, Oksanen J, Ovaskainen O (2020) Joint species distribution modelling with the R-package Hmsc. Methods in Ecology and Evolution 11:442-447. doi:10.1111/2041-210X.13345",
    "Bhattacharya A, Dunson DB (2011) Sparse Bayesian infinite factor models. Biometrika 98:291-306. doi:10.1093/biomet/asr013",
    "Durante D (2017) A note on the multiplicative gamma process. Statistics & Probability Letters 122:198-204.",
    "Zurell D, Pollock LJ, Thuiller W (2018) Do joint species distribution models reliably detect interspecific interactions from co-occurrence data in homogenous environments? Ecography 41:1812-1819.",
    "Blanchet FG, Cazelles K, Gravel D (2020) Co-occurrence is not evidence of ecological interactions. Ecology Letters 23:1050-1063. doi:10.1111/ele.13525",
    "Tjur T (2009) Coefficients of determination in logistic regression models - a new proposal: the coefficient of discrimination. The American Statistician 63:366-372.",
    "Pollock LJ, Tingley R, Morris WK, Golding N, O'Hara RB, Parris KM, Vesk PA, McCarthy MA (2014) Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution 5:397-406.",
    "Warton DI, Blanchet FG, O'Hara RB, Ovaskainen O, Taskinen S, Walker SC, Hui FKC (2015) So many variables: joint modeling in community ecology. Trends in Ecology & Evolution 30:766-779.",
    "Vehtari A, Gelman A, Simpson D, Carpenter B, Buerkner P-C (2021) Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. Bayesian Analysis 16:667-718.",
    "McGillycuddy M, Popovic G, Bolker BM, Warton DI (2025) Parsimoniously fitting large multivariate random effects in glmmTMB. Journal of Statistical Software 112(1). doi:10.18637/jss.v112.i01")
)
