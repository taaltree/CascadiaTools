# scr-open: open-population (Jolly-Seber) spatial capture-recapture
# (Gardner et al. 2010; Royle et al. 2014 ch. 16; Efford & Schofield 2020).
#
# Several primary periods (years/seasons), each with K secondary occasions;
# the population is closed within and open between periods (survival phi,
# entry probabilities gamma_t under data augmentation), activity centres are
# static; detection is half-normal. Density may depend on habitat covariates.
# Engines (settings$engine / settings$approach):
#   nimble  marginal (default; centres and alive states integrated out),
#           hmm (alive states integrated out, centres sampled),
#           latent (everything sampled; the textbook formulation)
#   numpyro marginal NUTS (Python subprocess; CPU, or GPU via scripts/ejob -g)
#   openCR  maximum-likelihood JSSAsecrf cross-check (compare_engines only)
# Ported from the SCR GPU benchmark in ../../ (python/simulate_open.py,
# python/fit_numpyro_open.py, R/fit_nimble_open.R); those files are untouched.

for (f in c("prepare.R", "simulate.R", "nimble_engines.R", "other_engines.R", "post.R", "ppc.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

so_settings <- list(
  engine = "nimble", approach = "marginal", seed = 1,
  niter = 30000, nburnin = 10000, nchains = 3, thin = 2,
  niter_latent = 30000, nburnin_latent = 10000, thin_latent = 5,
  numpyro_warmup = 1000, numpyro_samples = 1500, numpyro_gpu = FALSE, numpyro_float32 = FALSE,
  max_draws_per_chain = 500, surface_draws = 200, ppc_draws = 300,
  compare_engines = character(), M = NULL, buffer = NULL, sigma_max = NULL,
  density_covariates = NULL, max_mask_points = 2500)

so_fast <- merge_list(so_settings, list(niter = 10000, nburnin = 4000, thin = 2, max_draws_per_chain = 200,
                                        surface_draws = 0, ppc_draws = 0,
                                        numpyro_warmup = 500, numpyro_samples = 750,
                                        niter_latent = 15000, nburnin_latent = 5000))

so_run_engine <- function(engine, prep, settings) {
  switch(engine,
    nimble = switch(if (is.null(settings$approach)) "marginal" else settings$approach,
                    marginal = so_engine_nimble_marginal(prep, settings),
                    hmm = so_engine_nimble_s(prep, settings, "hmm"),
                    latent = so_engine_nimble_s(prep, settings, "latent"),
                    stop("approach must be marginal, hmm or latent")),
    numpyro = so_engine_numpyro(prep, settings),
    stop("unknown engine '", engine, "' (nimble or numpyro; openCR only via compare_engines)"))
}

so_bayes <- function(prep, settings, engine) {
  for (attempt in 1:3) {
    eng <- so_run_engine(engine, prep, settings)
    mx <- settings$max_draws_per_chain
    t_post0 <- proc.time()[["elapsed"]]
    rep_chains <- lapply(eng$chains, function(m) so_report_chain(prep, m, so_thin_idx(nrow(m), mx), prep$M))
    Nq <- stats::quantile(unlist(lapply(rep_chains, function(m) m[, "superpopulation_size"])), 0.995)
    if (Nq < 0.9 * prep$M || attempt == 3) break
    prep$notes <- c(prep$notes, sprintf("Data augmentation limit M = %d was too small (99.5%% quantile of the superpopulation = %.0f); refitted with M = %d.",
                                        prep$M, Nq, 2L * prep$M))
    prep$M <- 2L * prep$M
  }
  eng$t_post <- proc.time()[["elapsed"]] - t_post0
  s <- summarise_draws_eco(rep_chains)
  list(eng = eng, prep = prep, rep_chains = rep_chains, estimates = s$estimates, diag = s$diagnostics,
       M_hit = Nq >= 0.9 * prep$M)
}

so_fit <- function(data, settings) {
  settings <- merge_list(so_settings, settings)
  t_fit0 <- proc.time()[["elapsed"]]
  set.seed(settings$seed)
  prep <- so_prepare(data, settings)
  if (prep$n == 0) stop("no individuals were detected: spatial capture-recapture needs captures")
  engine <- settings$engine
  b <- so_bayes(prep, settings, engine)
  prep <- b$prep
  est <- b$estimates
  dg <- b$diag; dg$t_total <- b$eng$t_total; dg$t_post <- b$eng$t_post; dg$attempts_M <- prep$M; dg$notes <- prep$notes
  if (b$M_hit) { dg$converged <- FALSE
    dg$notes <- c(dg$notes, "The superpopulation still presses against the data-augmentation limit M: set settings$M larger.") }
  if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
    dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                    dg$rhat_max, dg$ess_min))
  if (!is.null(b$eng$timing$divergences) && b$eng$timing$divergences > 0)
    dg$notes <- c(dg$notes, sprintf("NumPyro reported %d divergent transitions.", b$eng$timing$divergences))
  pool <- so_pool(b$eng$chains, max(settings$surface_draws, settings$ppc_draws, 1))
  surf <- if (settings$surface_draws > 0) so_surfaces(prep, pool[so_thin_idx(nrow(pool), settings$surface_draws), , drop = FALSE], prep$M) else NULL
  ppc <- if (settings$ppc_draws > 0) so_ppc(prep, pool[so_thin_idx(nrow(pool), settings$ppc_draws), , drop = FALSE], prep$M) else NULL
  if (!is.null(ppc)) dg$notes <- c(dg$notes, sprintf(
    "Posterior predictive checks (%d draws; P(replicate >= observed), extreme if < 0.025 or > 0.975): %s.",
    nrow(ppc$rep), paste(sprintf("%s %.2f", names(ppc$p), ppc$p), collapse = ", ")))
  obj <- list(engine = b$eng$engine, prep = prep, chains = b$eng$chains, rep_chains = b$rep_chains,
              surfaces = surf, ppc = ppc, settings = settings, timing = b$eng$timing)
  dg$notes <- c(dg$notes, so_design_notes(prep, est))
  if (length(settings$compare_engines)) {
    cmp <- list(data.frame(engine = obj$engine, est[, c("param", "estimate", "sd", "lower", "upper")]))
    for (e2 in unlist(settings$compare_engines)) {
      r2 <- tryCatch({
        if (identical(e2, "openCR")) {
          oc <- so_engine_opencr(prep, settings)
          obj$opencr <- oc
          data.frame(engine = oc$engine, so_opencr_rows(prep, oc))
        } else {
          s2 <- settings; s2$compare_engines <- character(); s2$ppc_draws <- 0; s2$surface_draws <- 0; s2$M <- prep$M
          parts <- strsplit(e2, " ")[[1]]; s2$engine <- parts[1]; if (length(parts) > 1) s2$approach <- parts[2]
          if (paste(s2$engine, s2$approach) == obj$engine) NULL else {
            f2 <- so_fit(data, s2)
            data.frame(engine = f2$object$engine, f2$estimates[, c("param", "estimate", "sd", "lower", "upper")])
          }
        }
      }, error = function(e) { message("compare engine ", e2, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r2)) cmp[[length(cmp) + 1]] <- r2
    }
    obj$comparison <- do.call(rbind, cmp)
  }
  dg$t_fit <- proc.time()[["elapsed"]] - t_fit0
  list(estimates = est, diagnostics = dg, object = obj)
}

so_defaults <- list(
  params = list(density = 0.8, baseline_detection = 0.3, sigma = 800, survival = 0.8,
                per_capita_recruitment = 0.25, density_beta_canopy = 1.0),
  design = list(trap_grid_side = 7, trap_spacing = 1000, n_occasions = 5, n_periods = 5,
                buffer = 3000, mask_spacing = 400, density_covariate = TRUE, covariate_range = 2000,
                trap_failure = 0.05))

module <- ecomodule(
  id = "scr-open",
  title = "Open-population spatial capture-recapture (Jolly-Seber SCR)",
  description = paste(
    "How many animals are there each year (or season), how many survive and how many are recruited, and is density associated with habitat?",
    "Open-population spatial capture-recapture uses where individually identified animals (camera traps, hair snares, scat DNA, live traps)",
    "were detected across several primary periods to estimate abundance and density in each period, apparent survival between periods,",
    "recruitment (number of new animals and per-capita recruitment rate), population growth, the superpopulation ever alive,",
    "baseline detection (p0), the spatial scale of movement (sigma) and the effect of a habitat covariate on density,",
    "accounting for animals that were never detected."),
  defaults = so_defaults,
  simulate = function(params, design, seed) { s <- so_simulate(params, design, seed); s[c("data", "truth")] },
  fit = so_fit, settings = so_settings, fast_settings = so_fast,
  read_data = so_read_data, check_data = so_check_data, write_data = so_write_data,
  interpret = so_interpret, plots = function(fit, data) so_figures(fit, data, fig_opts()),
  figures = so_figures, methods_text = so_methods_text,
  power_params = c("density_beta_canopy"),
  null_values = c(density_beta_canopy = 0),
  engine = "nimble (marginal [default], hmm or latent data augmentation), numpyro (marginal NUTS); openCR (ML) cross-check",
  references = so_references)
