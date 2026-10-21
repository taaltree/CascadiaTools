# ecotools module `pva`: population viability analysis.
#
# Two analyses, chosen by settings$model (or data$model; design$model in
# simulation):
#   count   count-based PVA from a time series of abundance counts:
#           diffusion approximation (Dennis et al. 1991; Morris & Doak 2002)
#           or exponential-growth state-space model with observation error
#           (REML; Staples et al. 2004, Humbert et al. 2009). Quasi-extinction
#           probabilities with parametric-bootstrap intervals.
#   matrix  stage-structured stochastic projection: environmental
#           stochasticity (beta / lognormal / stretched-beta vital rates,
#           optional within-year correlation), demographic stochasticity,
#           ceiling density dependence, parameter uncertainty, stochastic
#           lambda (simulation + Tuljapurkar), elasticities, scenarios.
#
# Data structure (what simulate()$data and read_data() both return):
#   count:  list(model, df = data.frame(year, count), threshold, horizons, meta)
#   matrix: list(model, stages, A (k x k character expressions), rates (data.frame
#           rate,type,mean,se,sd_env,n_years,sampling_var,dist,min,max) or NULL,
#           demography (data.frame year,rate,n,events) or NULL, rate_types, n0,
#           threshold, count_stages, ceiling, env_correlation, horizons,
#           scenarios, meta)

`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("count.R", "matrix.R", "report.R"))
  sys.source(file.path(MODULE_DIR, "R", f), envir = environment())

PVA_MODELS <- c("count", "matrix")

PVA_SETTINGS <- list(
  model = NULL, threshold = NULL, horizons = NULL, tmax = NULL, seed = 1,
  # count
  method = "diffusion", correction = "discrete", nboot = 4000, ndraw_ss = 4000,
  # matrix
  n_param_draws = 400, n_traj_per_draw = 250, n_traj_point = 10000, demographic = TRUE,
  scenarios = NULL, ignore_parameter_uncertainty = FALSE)

pva_model <- function(data, settings = list()) {
  m <- settings$model %||% data$model
  if (is.null(m) || !m %in% PVA_MODELS) stop("settings: model must be 'count' or 'matrix'", call. = FALSE)
  if (!is.null(data$model) && !identical(m, data$model))
    stop(sprintf("settings$model is '%s' but the data were read/simulated for '%s'.", m, data$model), call. = FALSE)
  m
}

pva_simulate <- function(params, design, seed) {
  switch(design$model %||% "count",
         count = count_simulate(params, design, seed),
         matrix = matrix_simulate(params, design, seed),
         stop("design$model must be 'count' or 'matrix'"))
}

pva_fit <- function(data, settings) {
  st <- merge_list(PVA_SETTINGS, settings)
  if (pva_model(data, settings) == "count") {
    count_fit(data, st)
  } else matrix_fit(data, st)
}

pva_read_data <- function(config) {
  dc <- config$data
  m <- config$settings$model %||% dc$model
  if (is.null(m)) stop("config: set settings: model: count | matrix", call. = FALSE)
  switch(m, count = count_read_data(config, dc), matrix = matrix_read_data(config, dc), stop(sprintf("settings: model must be 'count' or 'matrix' (got '%s').", m), call. = FALSE))
}

pva_check_data <- function(data) switch(data$model, count = count_check_data(data), matrix = matrix_check_data(data))

pva_interpret <- function(fit, data) {
  out <- switch(fit$object$model, count = count_interpret(fit, data), matrix = matrix_interpret(fit, data))
  if (!isTRUE(fit$diagnostics$converged)) out <- c("WARNING: the fit did not converge; do not interpret the numbers below.", out)
  out
}

pva_methods <- function(fit, data, settings) switch(fit$object$model, count = count_methods(fit, data, settings),
                                                   matrix = matrix_methods(fit, data, settings))

pva_figures <- function(fit, data, opts = fig_opts()) switch(fit$object$model, count = count_figures(fit, data, opts),
                                                            matrix = matrix_figures(fit, data, opts))

pva_plots <- function(fit, data) pva_figures(fit, data, fig_opts())

module <- ecomodule(
  id = "pva",
  title = "Population viability analysis: count-based and matrix-model PVA",
  description = paste(
    "How likely is a population to fall below a critical size within the next 10, 25 or 50 years, how fast is it growing or declining,",
    "and which vital rates or management actions would most reduce that risk? Count-based PVA estimates the long-run growth rate and its",
    "year-to-year variability from a time series of counts (diffusion approximation, or a state-space model that separates counting error);",
    "matrix PVA projects a stage-structured population with environmental and demographic stochasticity, parameter uncertainty and management scenarios."),
  defaults = list(
    params = list(
      # count model (log-scale growth rate, process and observation variance)
      mu = -0.03, sigma2 = 0.03, tau2 = 0,
      # matrix model: first-year, second-year and adult survival, female young per adult female
      s_juv = 0.30, s_sub = 0.60, s_ad = 0.75, fec = 1.4,
      sd_s_juv = 0.08, sd_s_sub = 0.06, sd_s_ad = 0.05, sd_fec = 0.35),
    design = list(model = "count", method = "diffusion",
                  # count
                  n_years = 25, N0 = 500, first_year = 2000, p_missing = 0,
                  threshold = 50, horizons = c(10, 25, 50), truth_paths = 2e5,
                  # matrix (threshold/horizons above are reused; set threshold for the matrix design)
                  n0 = c(15, 45), count_stages = NULL, ceiling = 250, env_correlation = 0,
                  n_mon_years = 15, n_marked_juv = 60, n_marked_sub = 30, n_marked_ad = 80, n_breeders = 40)),
  simulate = pva_simulate, fit = pva_fit, read_data = pva_read_data, check_data = pva_check_data,
  interpret = pva_interpret, figures = pva_figures, plots = pva_plots, methods_text = pva_methods,
  settings = list(),
  fast_settings = list(nboot = 2000, ndraw_ss = 2000),  # matrix PVA recovery uses the full default draws
  power_params = c("growth_rate_mu", "lambda_stochastic"),
  null_values = c(growth_rate_mu = 0, lambda_stochastic = 1),
  engine = "hand-written R (diffusion / EGSS REML / vectorised stochastic matrix projection); popbio for cross-checks and vital-rate sensitivities",
  references = c(
    "Dennis B, Munholland PL, Scott JM (1991) Estimation of growth and extinction parameters for endangered species. Ecological Monographs 61:115-143.",
    "Morris WF, Doak DF (2002) Quantitative Conservation Biology: Theory and Practice of Population Viability Analysis. Sinauer, Sunderland MA.",
    "Holmes EE (2004) Beyond theory to application and evaluation: diffusion approximations for population viability analysis. Ecological Applications 14:1272-1293.",
    "Holmes EE (2001) Estimating risks in declining populations with poor data. PNAS 98:5072-5077.",
    "Staples DF, Taper ML, Dennis B (2004) Estimating population trend and process variation for PVA in the presence of sampling error. Ecology 85:923-929.",
    "Humbert J-Y, Mills LS, Horne JS, Dennis B (2009) A better way to estimate population trends. Oikos 118:1940-1946.",
    "Dennis B, Ponciano JM, Lele SR, Taper ML, Staples DF (2006) Estimating density dependence, process noise, and observation error. Ecological Monographs 76:323-341.",
    "Fieberg J, Ellner SP (2000) When is it meaningful to estimate an extinction probability? Ecology 81:2040-2047.",
    "Ellner SP, Holmes EE (2008) Resolving the debate on when extinction risk is predictable. Ecology Letters 11:E1-E5.",
    "Ludwig D (1999) Is it meaningful to estimate a probability of extinction? Ecology 80:298-310.",
    "Broadie M, Glasserman P, Kou SG (1997) A continuity correction for discrete barrier options. Mathematical Finance 7:325-349.",
    "Tuljapurkar S (1990) Population Dynamics in Variable Environments. Springer, New York.",
    "Caswell H (2001) Matrix Population Models: Construction, Analysis, and Interpretation, 2nd ed. Sinauer, Sunderland MA.",
    "Kendall BE (1998) Estimating the magnitude of environmental stochasticity in survivorship data. Ecological Applications 8:184-193.",
    "Stubben CJ, Milligan BG (2007) Estimating and analyzing demographic models using the popbio package in R. Journal of Statistical Software 22(11).",
    "Sabo JL, Holmes EE, Kareiva P (2004) Efficacy of simple viability models in ecological risk assessment: does density dependence matter? Ecology 85:328-341.",
    "Eberhardt LL, Knight RR, Blanchard BM (1986) Monitoring grizzly bear population trends. Journal of Wildlife Management 50:613-618.")
)
