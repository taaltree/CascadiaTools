# nmixture: binomial N-mixture models for repeated counts.
#   static  (closed; Royle 2004) with Poisson, negative binomial or zero-inflated
#           Poisson abundance; covariates on abundance (lambda) and detection (p).
#   dynamic (open; Dail & Madsen 2011, unmarked::pcountOpen) with initial abundance,
#           recruitment (gamma) and apparent survival (omega); select with data$model.
# Engines: unmarked (default, maximum likelihood), nimble (Bayesian, static only,
# nimbleEcology marginal N-mixture distributions) and rlik (hand-written likelihood,
# used to cross-check unmarked).

for (f in c("prep.R", "lik.R", "engines.R", "fit.R", "sim.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

module <- ecomodule(
  id = "nmixture",
  title = "N-mixture models for repeated counts (closed and open)",
  description = paste(
    "How many animals are at each site, and what habitat and survey conditions are associated with abundance and with detecting them,",
    "when the same sites are counted several times and some animals present are missed on every count?",
    "Static (closed, Royle 2004): expected abundance per site, abundance at each surveyed site and in total, covariate effects as rate ratios,",
    "and per-visit detection probability. Dynamic (open, Dail & Madsen 2011): initial abundance, recruitment, apparent survival and",
    "total abundance per season. Absolute abundance from counts alone depends strongly on the model assumptions; relative effects are more robust."),
  defaults = list(
    params = list(abundance_at_average_site = 4, detection_prob = 0.5, lambda_beta_forest = 0.5, p_beta_wind = -0.5,
                  nb_size = 2, zero_inflation_prob = 0.3,
                  initial_abundance_per_site = 4, recruitment_per_site = 1.2, apparent_survival = 0.7),
    design = list(model = "static", n_sites = 150, n_visits = 3, n_seasons = 5, prop_missing = 0.1,
                  abundance_dist = "P", p_site_sd = 0, p_visit_sd = 0, double_count_prob = 0)),
  simulate = nm_simulate, fit = nm_fit, settings = nm_settings, fast_settings = nm_fast,
  read_data = nm_read_data, check_data = nm_check_data,
  interpret = nm_interpret, plots = nm_plots, figures = nm_figures, methods_text = nm_methods_text,
  power_params = c("lambda_beta_forest", "p_beta_wind", "rate_ratio_forest", "abundance_ratio_last_first"),
  null_values = c(lambda_beta_forest = 0, p_beta_wind = 0, rate_ratio_forest = 1, abundance_ratio_last_first = 1),
  engine = "unmarked (pcount, pcountOpen), nimble + nimbleEcology, hand-written likelihood",
  references = c(
    "Royle JA (2004) N-mixture models for estimating population size from spatially replicated counts. Biometrics 60:108-115. doi:10.1111/j.0006-341X.2004.00142.x",
    "Kery M, Royle JA, Schmid H (2005) Modeling avian abundance from replicated counts using binomial mixture models. Ecological Applications 15:1450-1461. doi:10.1890/04-1120",
    "Dail D, Madsen L (2011) Models for estimating abundance from repeated counts of an open metapopulation. Biometrics 67:577-587. doi:10.1111/j.1541-0420.2010.01465.x",
    "Hostetler JA, Chandler RB (2015) Improved state-space models for inference about spatial and temporal variation in abundance from count data. Ecology 96:1713-1723. doi:10.1890/14-1487.1",
    "Dennis EB, Morgan BJT, Ridout MS (2015) Computational aspects of N-mixture models. Biometrics 71:237-246. doi:10.1111/biom.12246",
    "Barker RJ, Schofield MR, Link WA, Sauer JR (2018) On the reliability of N-mixture models for count data. Biometrics 74:369-377. doi:10.1111/biom.12734",
    "Link WA, Schofield MR, Barker RJ, Sauer JR (2018) On the robustness of N-mixture models. Ecology 99:1547-1551. doi:10.1002/ecy.2362",
    "Kery M (2018) Identifiability in N-mixture models: a large-scale screening test with bird data. Ecology 99:281-288. doi:10.1002/ecy.2093",
    "Knape J, Arlt D, Barraquand F, et al. (2018) Sensitivity of binomial N-mixture models to overdispersion: the importance of assessing model fit. Methods in Ecology and Evolution 9:2102-2114. doi:10.1111/2041-210X.13062",
    "Duarte A, Adams MJ, Peterson JT (2018) Fitting N-mixture models to count data with unmodeled heterogeneity: bias, diagnostics, and alternative approaches. Ecological Modelling 374:51-59. doi:10.1016/j.ecolmodel.2018.02.007",
    "Fiske I, Chandler R (2011) unmarked: An R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23. doi:10.18637/jss.v043.i10",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487",
    "Goldstein BR, de Valpine P (2022) Comparing N-mixture models and GLMMs for relative abundance estimation in a citizen science dataset. Scientific Reports 12:12276. doi:10.1038/s41598-022-16368-z",
    "Kery M, Royle JA (2016) Applied Hierarchical Modeling in Ecology, Vol. 1, chapter 6. Academic Press."))
