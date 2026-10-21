# Nimble helpers shared by Bayesian modules.
#
# nimble_run() builds, compiles and runs a model, returning per-chain sample
# matrices plus the contract-shaped estimates/diagnostics. Compilation costs
# ~30-100 s, so within one R process the compiled model is cached by a key
# (module id + dimension signature). Simulation-recovery replicates with the
# same design then only pay for compilation once per worker: new data and
# inits are pushed into the compiled model and the MCMC is rerun.
#
# For caching to be valid, everything that changes between datasets must be
# passed as `data` (covariates included), and `constants` must hold only
# dimensions. The cache key must capture every dimension.

.eco_nimble_cache <- new.env(parent = emptyenv())

nimble_run <- function(code, constants, data, inits, monitors,
                       niter = 5000, nburnin = 1000, nchains = 3, thin = 1,
                       cache_key = NULL, configure = NULL, summary_params = NULL,
                       seed = NULL, calculate = FALSE) {
  eco_require("nimble")
  suppressPackageStartupMessages(requireNamespace("nimble"))
  if (!is.null(seed)) set.seed(seed)
  init_list <- function(k) if (is.function(inits)) inits() else inits
  t0 <- proc.time()[["elapsed"]]
  ent <- if (!is.null(cache_key)) .eco_nimble_cache[[cache_key]] else NULL
  if (is.null(ent)) {
    m <- nimble::nimbleModel(code, constants = constants, data = data,
                             inits = init_list(1), calculate = calculate, check = FALSE)
    conf <- nimble::configureMCMC(m, monitors = monitors, thin = thin, print = FALSE)
    if (is.function(configure)) configure(conf, m)
    mcmc <- nimble::buildMCMC(conf)
    # One compile directory per process: forked recovery workers share
    # tempdir(), and concurrent compilations there collide on generated file
    # names ("Failed to create the shared library", "not in load table").
    cdir <- file.path(tempdir(), paste0("nimble_pid", Sys.getpid()))
    cm <- nimble::compileNimble(m, showCompilerOutput = FALSE, dirName = cdir)
    cmcmc <- nimble::compileNimble(mcmc, project = m, showCompilerOutput = FALSE, dirName = cdir)
    ent <- list(cm = cm, cmcmc = cmcmc)
    if (!is.null(cache_key)) .eco_nimble_cache[[cache_key]] <- ent
  } else {
    ent$cm$resetData()
    ent$cm$setData(data)
  }
  t_setup <- proc.time()[["elapsed"]] - t0
  samples <- vector("list", nchains)
  for (k in seq_len(nchains)) {
    ini <- init_list(k)
    for (v in names(ini)) ent$cm[[v]] <- ini[[v]]
    ent$cm$calculate()
    ent$cmcmc$run(niter, nburnin = nburnin, thin = thin, reset = TRUE, progressBar = FALSE)
    samples[[k]] <- as.matrix(ent$cmcmc$mvSamples)
  }
  t_total <- proc.time()[["elapsed"]] - t0
  s <- summarise_draws_eco(samples, params = summary_params)
  s$diagnostics$t_setup <- t_setup
  s$diagnostics$t_total <- t_total
  s$samples <- samples
  s
}

# Posterior draws (all chains stacked) as a matrix: handy for derived
# quantities (e.g. N, conditional effects with full uncertainty).
stack_draws <- function(samples) do.call(rbind, lapply(samples, as.matrix))

# Summarise a vector/matrix of derived draws into contract rows.
derived_rows <- function(draws, names = NULL, prob = 0.95) {
  draws <- as.matrix(draws)
  if (!is.null(names)) colnames(draws) <- names
  a <- (1 - prob) / 2
  data.frame(param = colnames(draws),
             estimate = colMeans(draws),
             sd = apply(draws, 2, stats::sd),
             lower = apply(draws, 2, stats::quantile, a, na.rm = TRUE),
             upper = apply(draws, 2, stats::quantile, 1 - a, na.rm = TRUE),
             row.names = NULL)
}

# Contract rows from a likelihood fit: Wald intervals on the link scale.
wald_rows <- function(est, se, names, transform = identity, prob = 0.95) {
  z <- stats::qnorm(1 - (1 - prob) / 2)
  lo <- transform(est - z * se); hi <- transform(est + z * se)
  data.frame(param = names, estimate = transform(est),
             sd = if (identical(transform, identity)) se else NA_real_,
             lower = pmin(lo, hi), upper = pmax(lo, hi), row.names = NULL)
}
