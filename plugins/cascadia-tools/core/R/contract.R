# The module contract. Every analysis in modules/<id>/module.R ends by
# returning ecomodule(...). The core (recovery, power, reports, CLI) only ever
# talks to modules through these fields, so any module that validates here
# gets simulation-recovery, power analysis and HTML reports for free.
#
# Required fields
#   id, title         short id (matches the folder) and human title
#   description       one paragraph: what question the analysis answers
#   defaults          list(params = <named list>, design = <named list>)
#                     params: true parameter values used by simulate()
#                     design: sampling design (n sites, occasions, ...)
#   simulate(params, design, seed)
#       -> list(data = <same structure read_data() returns>,
#               truth = named numeric vector of estimands, names matching
#                       fit()$estimates$param)
#   fit(data, settings)
#       -> list(estimates = data.frame(param, estimate, sd, lower, upper),
#               diagnostics = list(converged = TRUE/FALSE, rhat_max, ess_min,
#                                  notes = character()),
#               object = <raw engine fit, optional>)
#       `sd` is posterior SD (Bayesian) or SE (likelihood); lower/upper are
#       95% credible or confidence limits.
#
# Optional fields
#   settings          default settings list for fit() on real data
#   fast_settings     cheaper settings for simulation-recovery replicates
#   read_data(config) -> data from user files (config$data paths etc.)
#   check_data(data)  -> character vector of problems (empty = OK)
#   interpret(fit, data) -> character vector of plain-language paragraphs
#   plots(fit, data)  -> named list of ggplot objects (conditional effects...)
#   figures(fit, data, opts) -> the full publication figure suite: named list
#                     of ggplots built with core/R/figures.R constructors, each
#                     tagged with with_caption(); honours opts (fig_opts():
#                     labels, units, ci, show_data...). See docs/FIGURES.md.
#   methods_text(fit, data, settings) -> character; methods paragraph
#   null_values       named numeric: value meaning "no effect" per param
#                     (power = P(CI excludes null)). Default 0 for any param
#                     listed in power_params.
#   power_params      params for which power is meaningful (effects)
#   engine            "nimble", "unmarked", "numpyro", ... (documentation)
#   references        character vector of citations the model follows

ecomodule <- function(id, title, description, defaults, simulate, fit,
                      settings = list(), fast_settings = NULL,
                      read_data = NULL, check_data = NULL, interpret = NULL,
                      plots = NULL, methods_text = NULL, null_values = NULL,
                      power_params = character(), engine = NA_character_,
                      references = character(), ...) {
  m <- list(id = id, title = title, description = description,
            defaults = defaults, simulate = simulate, fit = fit,
            settings = settings,
            fast_settings = if (is.null(fast_settings)) settings else fast_settings,
            read_data = read_data, check_data = check_data,
            interpret = interpret, plots = plots, methods_text = methods_text,
            null_values = null_values, power_params = power_params,
            engine = engine, references = references, ...)
  class(m) <- "ecomodule"
  validate_module(m)
  m
}

validate_module <- function(m) {
  need <- c("id", "title", "description", "defaults", "simulate", "fit")
  miss <- need[vapply(need, function(k) is.null(m[[k]]), logical(1))]
  if (length(miss)) stop("module missing: ", paste(miss, collapse = ", "))
  if (!all(c("params", "design") %in% names(m$defaults)))
    stop("defaults must contain params and design")
  for (f in c("simulate", "fit")) if (!is.function(m[[f]])) stop(f, " must be a function")
  invisible(TRUE)
}

# Load modules/<id>/module.R in its own environment. The file must evaluate to
# an ecomodule as its last expression. MODULE_DIR is visible inside it so
# modules can source helpers (Nimble code, python scripts) relative to self.
load_module <- function(id) {
  dir <- file.path(eco_root(), "modules", id)
  f <- file.path(dir, "module.R")
  if (!file.exists(f)) stop("no module '", id, "' (looked for ", f, ")")
  env <- new.env(parent = globalenv())
  env$MODULE_DIR <- dir
  m <- sys.source(f, envir = env, keep.source = FALSE, toplevel.env = env)
  m <- get_last_module(env, f)
  m$dir <- dir
  m
}

get_last_module <- function(env, f) {
  # sys.source returns invisibly NULL; modules assign `module <- ecomodule(...)`
  if (exists("module", envir = env, inherits = FALSE)) return(get("module", envir = env))
  stop(f, " must assign `module <- ecomodule(...)`")
}

list_modules <- function() {
  d <- list.dirs(file.path(eco_root(), "modules"), recursive = FALSE, full.names = FALSE)
  d[file.exists(file.path(eco_root(), "modules", d, "module.R")) & !startsWith(d, "_")]
}

# Merge user-supplied overrides into defaults (one level deep).
merge_list <- function(base, over) {
  if (is.null(over)) return(base)
  for (k in names(over)) base[[k]] <- over[[k]]
  base
}

# Standard estimates table checker, used by recovery and selftest.
check_estimates <- function(est) {
  need <- c("param", "estimate", "sd", "lower", "upper")
  if (!is.data.frame(est) || !all(need %in% names(est)))
    stop("fit()$estimates must be a data.frame with columns: ", paste(need, collapse = ", "))
  est$param <- as.character(est$param)
  est
}

# Helper for modules: tidy summary of an MCMC sample matrix/list (Bayesian
# engines). Returns the estimates table and diagnostics in contract form.
summarise_draws_eco <- function(samples, params = NULL, prob = 0.95) {
  eco_require(c("posterior"))
  if (inherits(samples, "mcmc.list") || is.list(samples) && !is.data.frame(samples)) {
    dr <- posterior::as_draws_array(lapply(samples, as.matrix) |> as_chain_array())
  } else {
    dr <- posterior::as_draws_matrix(as.matrix(samples))
  }
  if (!is.null(params)) dr <- posterior::subset_draws(dr, variable = params)
  a <- (1 - prob) / 2
  s <- posterior::summarise_draws(
    dr, mean = mean, sd = stats::sd,
    q = ~stats::quantile(.x, c(a, 1 - a), na.rm = TRUE),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  names(s)[4:5] <- c("lower", "upper")
  est <- data.frame(param = s$variable, estimate = s$mean, sd = s$sd,
                    lower = s$lower, upper = s$upper,
                    rhat = s$rhat, ess_bulk = s$ess_bulk)
  diag <- list(rhat_max = max(est$rhat, na.rm = TRUE),
               ess_min = min(est$ess_bulk, na.rm = TRUE))
  diag$converged <- is.finite(diag$rhat_max) && diag$rhat_max < 1.05
  list(estimates = est, diagnostics = diag)
}

# list of per-chain matrices -> iterations x chains x variables array
as_chain_array <- function(mats) {
  mats <- lapply(mats, as.matrix)
  nit <- min(vapply(mats, nrow, 1L))
  v <- colnames(mats[[1]])
  arr <- array(NA_real_, c(nit, length(mats), length(v)),
               dimnames = list(NULL, NULL, v))
  for (k in seq_along(mats)) arr[, k, ] <- mats[[k]][seq_len(nit), v, drop = FALSE]
  arr
}
