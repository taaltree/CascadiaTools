# ecotools module `fisheries-stock`: mortality, recruitment and stock dynamics.
#
# Three analyses, chosen by settings$model (or data$model; design$model in
# simulation):
#   catch_curve         total mortality Z from catch-at-age (Chapman-Robson,
#                       weighted regression, Poisson-lognormal GLMM, Bayesian
#                       Nimble catch curve; multi-year cohort catch curve)
#   stock_recruit       Ricker / Beverton-Holt stock-recruitment with
#                       reference points; least squares or Bayesian (Nimble)
#                       with AR(1) recruitment residuals
#   surplus_production  Bayesian state-space Schaefer production model
#                       (Meyer & Millar 1999) for catch + abundance index;
#                       MSY, B/B_MSY, F/F_MSY, Kobe plot
#
# Data structure (what simulate()$data and read_data() both return):
#   catch_curve:        list(model, df = data.frame(age, catch), long = NULL or
#                       data.frame(year, age, catch) for multi-year samples, meta)
#   stock_recruit:      list(model, sr_form, spr0, df = data.frame(year, spawners, recruits), meta)
#   surplus_production: list(model, df = data.frame(year, catch, index), meta)

for (f in c("common.R", "catch_curve.R", "stock_recruit.R", "surplus_production.R"))
  sys.source(file.path(MODULE_DIR, "R", f), envir = environment())

FS_MODELS <- c("catch_curve", "stock_recruit", "surplus_production")

fs_model <- function(data, settings = list()) {
  m <- settings$model %||% data$model
  if (is.null(m) || !m %in% FS_MODELS)
    stop("settings: model must be one of ", paste(FS_MODELS, collapse = ", "), call. = FALSE)
  if (!is.null(data$model) && !identical(m, data$model))
    stop(sprintf("settings$model is '%s' but the data were read/simulated for '%s'.", m, data$model), call. = FALSE)
  m
}

fs_simulate <- function(params, design, seed) {
  switch(design$model %||% "catch_curve",
         catch_curve = cc_simulate(params, design, seed),
         stock_recruit = sr_simulate(params, design, seed),
         surplus_production = sp_simulate(params, design, seed),
         stop("design$model must be one of ", paste(FS_MODELS, collapse = ", ")))
}

fs_fit <- function(data, settings) {
  st <- merge_list(FS_SETTINGS, settings)
  switch(fs_model(data, settings),
         catch_curve = cc_fit(data, st),
         stock_recruit = sr_fit(data, st),
         surplus_production = sp_fit(data, st))
}

fs_read_data <- function(config) {
  dc <- config$data
  m <- config$settings$model %||% dc$model
  if (is.null(m)) stop("config: set settings: model: catch_curve | stock_recruit | surplus_production", call. = FALSE)
  switch(m, catch_curve = cc_read_data(config, dc), stock_recruit = sr_read_data(config, dc),
         surplus_production = sp_read_data(config, dc), stop("unknown model ", m))
}

fs_check_data <- function(data) switch(data$model, catch_curve = cc_check_data(data),
                                       stock_recruit = sr_check_data(data), surplus_production = sp_check_data(data))

fs_interpret <- function(fit, data) {
  out <- switch(fit$object$model, catch_curve = cc_interpret(fit, data),
                stock_recruit = sr_interpret(fit, data), surplus_production = sp_interpret(fit, data))
  if (!isTRUE(fit$diagnostics$converged)) out <- c("WARNING: the fit did not converge; do not interpret the numbers below.", out)
  out
}

fs_methods <- function(fit, data, settings) switch(fit$object$model, catch_curve = cc_methods(fit, data, settings),
                                                  stock_recruit = sr_methods(fit, data, settings),
                                                  surplus_production = sp_methods(fit, data, settings))

fs_figures <- function(fit, data, opts = fig_opts()) switch(fit$object$model, catch_curve = cc_figures(fit, data, opts),
                                                           stock_recruit = sr_figures(fit, data, opts),
                                                           surplus_production = sp_figures(fit, data, opts))

fs_plots <- function(fit, data) {
  f <- fs_figures(fit, data, fig_opts())
  f[setdiff(names(f), c("trace"))]
}

FS_SETTINGS <- list(
  model = NULL,
  # catch curve
  # primary = NULL: module default for single-year data (CC_DEFAULT_PRIMARY), cohort_glmm for multi-year data
  age_full = "peak_plus", age_full_regression = "peak", age_last = NULL, overdispersion = TRUE, primary = NULL,
  glmm_interval = NULL, n_boot = NULL, zero_tail = 5, bayes = FALSE, prior_yc_sd = 1, z_trend = FALSE,
  # stock-recruitment
  sr_form = NULL, sr_engine = "ls", ar1 = TRUE, spr0 = NULL, n_draws = 4000, beta_max = NULL,
  # MCMC (stock_recruit with sr_engine = nimble; surplus_production)
  niter = NULL, nburnin = NULL, nchains = 3, thin = NULL, seed = 1, fast = FALSE,
  # surplus-production priors: c(median, sdlog) and gamma(shape, rate) on precisions
  prior_r = c(0.3, 0.6), prior_K = NULL, prior_sigma_proc = c(3.785518, 0.010223),
  prior_sigma_obs = c(1.708603, 0.008613854), prior_P1 = c(1, 0.2))

module <- ecomodule(
  id = "fisheries-stock",
  title = "Fish stock assessment: catch curves, stock-recruitment and surplus production",
  description = paste(
    "How fast are fish dying, how many recruits does a spawning stock produce, and how big is the stock relative to the level that supports maximum sustainable yield?",
    "Catch curves estimate total mortality Z (and annual survival) from catch-at-age; stock-recruitment models (Ricker, Beverton-Holt) estimate productivity, density dependence,",
    "S_MSY, U_MSY and steepness; a Bayesian state-space Schaefer surplus-production model estimates r, K, MSY and stock status (B/B_MSY, F/F_MSY) from catch and an abundance index."),
  defaults = list(
    params = list(
      # catch curve
      Z = 0.5, recruit_sd_cc = 0.3, sel_age50 = 2.5, sel_slope = 4,
      # stock-recruitment (spawners and recruits in thousands of fish)
      sr_ln_alpha = 1.5, sr_beta = 0.02, sr_sigma = 0.4, sr_phi = 0,
      # surplus production (biomass/catch in tonnes, index = q B)
      r = 0.4, K = 1000, q = 0.002, sigma_proc = 0.05, sigma_obs = 0.15, P1 = 1),
    design = list(model = "catch_curve",
                  n_fish = 400, max_age = 20,
                  sr_form = "ricker", n_years = 30, S_mean = 60, S_cv = 0.6, S_phi = 0.5, spr0 = 1,
                  sp_years = 30, catch_pattern = "contrast", F_peak = 2, F_end = 1.6, F_cv = 0.1, index_missing = 0)),
  simulate = fs_simulate, fit = fs_fit, read_data = fs_read_data, check_data = fs_check_data,
  interpret = fs_interpret, figures = fs_figures, plots = fs_plots, methods_text = fs_methods,
  settings = list(),
  fast_settings = list(fast = TRUE, n_draws = 2000),
  power_params = c("density_dependence_beta", "B_over_BMSY_final"),
  null_values = c(density_dependence_beta = 0, B_over_BMSY_final = 1),
  engine = "FSA, glmmTMB and nimble (catch curves); stats lm/optim or nimble (stock-recruitment); nimble (surplus production)",
  references = c(
    "Chapman DG, Robson DS (1960) The analysis of a catch curve. Biometrics 16:354-368.",
    "Hoenig JM, Lawing WD, Hoenig NA (1983) Using mean age, mean length and median length data to estimate the total mortality rate. ICES CM 1983/D:23.",
    "Smith MW, Then AY, Wor C, Ralph G, Pollock KH, Hoenig JM (2012) Recommendations for catch-curve analysis. North American Journal of Fisheries Management 32:956-967.",
    "Millar RB (2015) A better estimator of mortality rate from age-frequency data. Canadian Journal of Fisheries and Aquatic Sciences 72:364-375.",
    "Nelson GA (2019) Bias in common catch-curve methods applied to age frequency data from fish surveys. ICES Journal of Marine Science 76(7):2090. doi:10.1093/icesjms/fsz085",
    "Maceina MJ, Bettoli PW (1998) Variation in largemouth bass recruitment in four mainstream impoundments on the Tennessee River. North American Journal of Fisheries Management 18:998-1003.",
    "Brooks ME, Kristensen K, van Benthem KJ, Magnusson A, Berg CW, Nielsen A, Skaug HJ, Maechler M, Bolker BM (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9:378-400.",
    "Ogle DH, Doll JC, Wheeler AP, Dinno A (2025) FSA: Simple Fisheries Stock Assessment Methods. R package version 0.10.1.",
    "Ricker WE (1954) Stock and recruitment. Journal of the Fisheries Research Board of Canada 11:559-623.",
    "Beverton RJH, Holt SJ (1957) On the Dynamics of Exploited Fish Populations. Fishery Investigations Series II, 19. HMSO, London.",
    "Hilborn R, Walters CJ (1992) Quantitative Fisheries Stock Assessment: Choice, Dynamics and Uncertainty. Chapman & Hall, New York.",
    "Scheuerell MD (2016) An explicit solution for calculating optimum spawning stock size from Ricker's stock recruitment model. PeerJ 4:e1623.",
    "Mace PM, Doonan IJ (1988) A generalised bioeconomic simulation model for fish population dynamics. New Zealand Fisheries Assessment Research Document 88/4.",
    "Walters CJ (1985) Bias in the estimation of functional relationships from time series data. Canadian Journal of Fisheries and Aquatic Sciences 42:147-149.",
    "Schaefer MB (1954) Some aspects of the dynamics of populations important to the management of the commercial marine fisheries. Inter-American Tropical Tuna Commission Bulletin 1:27-56.",
    "Meyer R, Millar RB (1999) BUGS in Bayesian stock assessments. Canadian Journal of Fisheries and Aquatic Sciences 56:1078-1086.",
    "Hilborn R (1979) Comparison of fisheries control systems that utilize catch and effort data. Journal of the Fisheries Research Board of Canada 36:1477-1489.",
    "Polacheck T, Hilborn R, Punt AE (1993) Fitting surplus production models: comparing methods and measuring uncertainty. Canadian Journal of Fisheries and Aquatic Sciences 50:2597-2607.",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.")
)
