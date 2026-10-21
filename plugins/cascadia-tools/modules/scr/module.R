# scr: closed-population spatial capture-recapture (Efford 2004; Borchers &
# Efford 2008; Royle & Young 2008; Royle et al. 2014).
#
# Density (homogeneous or log-linear in several habitat covariates on a
# habitat mask from a CSV grid, raster or polygon), half-normal detection with
# binomial / Bernoulli / Poisson (hazard lambda0) observation models,
# detector-, occasion- and individual-level (sex) covariates on p0/lambda0,
# detector covariates and sex on sigma, sex ratio with unknown sex latent,
# trap-specific behavioural response, effort / trap-operation matrices
# (noninvasive genetic sampling). Engines (settings$engine):
#   nimble  (approach "marginal" [default; all variants] or "latent" [basic binomial])
#   numpyro (marginal, NUTS; binomial/Poisson + mask + density/detector covariates)
#   secr    (maximum likelihood; independent cross-check of every variant)
# Ported from the SCR GPU benchmark in ../../ (python/simulate.py,
# python/fit_numpyro.py, R/fit_nimble.R); those files are not modified.

for (f in c("prepare.R", "simulate.R", "nimble_engines.R", "other_engines.R", "post.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

scr_settings <- list(
  engine = "nimble", approach = "marginal", sampler = "RW_block", seed = 1,
  niter = 20000, nburnin = 5000, nchains = 3, thin = 1,
  niter_latent = 30000, nburnin_latent = 10000, thin_latent = 5,
  numpyro_warmup = 1000, numpyro_samples = 2000, numpyro_gpu = FALSE, numpyro_float32 = FALSE,
  max_draws_per_chain = 2000, surface_draws = 200, ppc_draws = 300,
  compare_engines = character(), M = NULL, buffer = NULL, sigma_max = NULL, lambda0_max = 10,
  obs_model = NULL, density_covariates = NULL, detection_covariates = NULL,
  occasion_detection_covariates = NULL, sigma_covariates = NULL, sex_effects = NULL, behaviour = NULL,
  buffer_check = FALSE, max_mask_points = 4000)

scr_fast <- merge_list(scr_settings, list(niter = 12000, nburnin = 4000, max_draws_per_chain = 1000,
                                          surface_draws = 0, ppc_draws = 0,
                                          numpyro_warmup = 500, numpyro_samples = 1000,
                                          niter_latent = 15000, nburnin_latent = 5000))

scr_run_engine <- function(engine, prep, settings) {
  switch(engine,
    nimble = if (identical(settings$approach, "latent")) scr_engine_nimble_latent(prep, settings)
             else scr_engine_nimble_marginal(prep, settings),
    numpyro = scr_engine_numpyro(prep, settings),
    secr = scr_engine_secr(prep, settings),
    stop("unknown engine '", engine, "' (nimble, numpyro or secr)"))
}

# Bayesian engine -> contract estimates (with augmentation check / auto-refit)
scr_bayes <- function(prep, settings, engine) {
  for (attempt in 1:3) {
    eng <- scr_run_engine(engine, prep, settings)
    prep$M_used <- if (!is.null(eng$M)) eng$M else prep$M
    mx <- settings$max_draws_per_chain
    rep_chains <- lapply(eng$chains, function(m) scr_report_chain(prep, m, scr_thin_idx(nrow(m), mx)))
    Nq <- stats::quantile(unlist(lapply(rep_chains, function(m) m[, "abundance"])), 0.995)
    if (Nq < 0.9 * prep$M_used || attempt == 3) break
    settings$M <- 2 * prep$M_used; prep$M <- 2 * prep$M_used
    prep$notes <- c(prep$notes, sprintf("Data augmentation limit M = %d was too small (99.5%% quantile of N = %.0f); refitted with M = %d.",
                                        prep$M_used, Nq, prep$M))
  }
  s <- summarise_draws_eco(rep_chains)
  list(eng = eng, prep = prep, rep_chains = rep_chains, estimates = s$estimates, diag = s$diagnostics,
       M_hit = Nq >= 0.9 * prep$M_used)
}

scr_fit <- function(data, settings) {
  settings <- merge_list(scr_settings, settings)
  set.seed(settings$seed)
  prep <- scr_prepare(data, settings)
  if (prep$n == 0) stop("no individuals were detected: SCR cannot estimate density from zero captures")
  if (!is.null(data$mask_notes)) prep$notes <- c(data$mask_notes, prep$notes)
  if (!is.null(data$qc)) prep$notes <- c(prep$notes, sprintf("Sample QC: %d of %d genetic samples passed quality control (column '%s') and were analysed; %d were excluded.",
                                                             data$qc$n_pass, data$qc$n_samples, data$qc$column, data$qc$n_fail))
  engine <- settings$engine
  if (engine == "secr") {
    e <- scr_engine_secr(prep, settings)
    est <- e$estimates
    dg <- list(converged = e$converged, notes = prep$notes, t_total = e$t_total)
    if (!isTRUE(e$converged)) dg$notes <- c(dg$notes, "secr: the optimiser did not report convergence or some standard errors are not finite; do not use these estimates.")
    obj <- list(engine = "secr", prep = prep, secr = e$fit, secr_b = e$b, secr_V = e$V, secr_Mbig = e$Mbig, settings = settings)
    nd <- max(settings$surface_draws, settings$ppc_draws)
    if (nd > 0) {
      pd <- scr_pdraws(list(object = obj), nd)
      if (settings$surface_draws > 0) obj$surfaces <- scr_surfaces(prep, scr_with_M(pd[scr_thin_idx(nd, settings$surface_draws), , drop = FALSE], pd))
      if (settings$ppc_draws > 0) obj$ppc <- scr_ppc(prep, scr_with_M(pd[scr_thin_idx(nd, settings$ppc_draws), , drop = FALSE], pd))
      if (!is.null(obj$ppc)) dg$notes <- c(dg$notes, sprintf(
        "Parametric predictive checks (%d draws from the MLE sampling distribution; extreme if < 0.025 or > 0.975): %s.",
        nrow(obj$ppc$rep), paste(sprintf("%s %.2f", names(obj$ppc$p), obj$ppc$p), collapse = ", ")))
    }
  } else {
    b <- scr_bayes(prep, settings, engine)
    prep <- b$prep
    est <- b$estimates
    dg <- b$diag; dg$t_total <- b$eng$t_total
    dg$notes <- prep$notes
    if (b$M_hit) { dg$converged <- FALSE
      dg$notes <- c(dg$notes, "Abundance still presses against the data-augmentation limit M: set settings$M larger.") }
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                      dg$rhat_max, dg$ess_min))
    if (!is.null(b$eng$timing$divergences) && b$eng$timing$divergences > 0)
      dg$notes <- c(dg$notes, sprintf("NumPyro reported %d divergent transitions.", b$eng$timing$divergences))
    pool <- scr_pool(b$eng$chains, max(settings$surface_draws, settings$ppc_draws, 1))
    surf <- if (settings$surface_draws > 0) scr_surfaces(prep, pool[scr_thin_idx(nrow(pool), settings$surface_draws), , drop = FALSE]) else NULL
    ppc <- if (settings$ppc_draws > 0) scr_ppc(prep, pool[scr_thin_idx(nrow(pool), settings$ppc_draws), , drop = FALSE]) else NULL
    if (!is.null(ppc)) dg$notes <- c(dg$notes, sprintf(
      "Posterior predictive checks (%d draws; P(replicate >= observed), extreme if < 0.025 or > 0.975): %s.",
      nrow(ppc$rep), paste(sprintf("%s %.2f", names(ppc$p), ppc$p), collapse = ", ")))
    dg$ppc <- ppc
    obj <- list(engine = trimws(paste(engine, if (engine == "nimble") settings$approach else "")), prep = prep,
                chains = b$eng$chains, rep_chains = b$rep_chains, surfaces = surf, ppc = ppc, settings = settings)
  }
  dg$notes <- c(dg$notes, scr_design_notes(prep, est))
  # optional cross-check engines on the same data
  if (length(settings$compare_engines)) {
    cmp <- list(data.frame(engine = obj$engine, est[, c("param", "estimate", "sd", "lower", "upper")]))
    for (e2 in setdiff(unlist(settings$compare_engines), engine)) {
      s2 <- settings; s2$engine <- e2; s2$compare_engines <- character(); s2$ppc_draws <- 0; s2$surface_draws <- 0; s2$buffer_check <- FALSE
      r2 <- tryCatch(scr_fit(data, s2), error = function(e) { message("compare engine ", e2, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r2)) cmp[[length(cmp) + 1]] <- data.frame(engine = r2$object$engine, r2$estimates[, c("param", "estimate", "sd", "lower", "upper")])
    }
    obj$comparison <- do.call(rbind, cmp)
  }
  if (isTRUE(settings$buffer_check)) {
    bc <- tryCatch(scr_buffer_check(data, settings, list(estimates = est, object = obj)),
                   error = function(e) { message("buffer check failed: ", conditionMessage(e)); NULL })
    if (!is.null(bc)) { obj$buffer_check <- bc$table; dg$notes <- c(dg$notes, bc$note) }
  }
  list(estimates = est, diagnostics = dg, object = obj)
}

# Buffer check: refit with a larger state space (mask rebuilt from traps /
# polygon / raster with buffer + 2 sigma, at least 1.5x) or, when the mask is a
# fixed grid that cannot be enlarged, with the state space trimmed by 2 sigma;
# density should not change if the buffer is adequate (>= 3 sigma).
scr_buffer_check <- function(data, settings, fit) {
  prep <- fit$object$prep; est <- fit$estimates
  sg <- est$estimate[est$param %in% c("sigma", "sigma_female", "sigma_male")]
  sig <- max(sg)
  cur <- prep$buffer
  s2 <- settings; s2$buffer_check <- FALSE; s2$compare_engines <- character(); s2$ppc_draws <- 0; s2$surface_draws <- 0
  d2 <- data
  if (identical(data$mask_kind, "csv")) {
    new <- max(cur - 2 * sig, 2 * sig)
    if (new >= cur * 0.95) return(list(table = NULL, note = "Buffer check skipped: the supplied habitat grid is too small to trim and cannot be enlarged."))
    s2[["buffer"]] <- new
  } else {
    new <- max(cur + 2 * sig, 1.5 * cur)
    d2[["mask"]] <- scr_build_mask(data, new, prep$spacing); d2$buffer <- new; d2$mask_notes <- NULL
    s2[["buffer"]] <- NULL; s2$max_mask_points <- max(settings$max_mask_points, 2.5 * nrow(d2[["mask"]]))
  }
  r2 <- scr_fit(d2, s2)
  # Compare density over the SAME region in both fits: pixels within min(cur, new)
  # of the nearest detector. With density covariates the state-space average
  # legitimately changes when habitat is added or removed, so comparing the
  # overall averages would confuse habitat composition with buffer adequacy.
  rad <- min(cur, new)
  dreg <- function(ob) {
    pr <- ob$prep
    keep <- scr_dist_to_traps(pr$mask, pr$traps) <= rad + 1e-9 * rad
    pd <- scr_pdraws(list(object = ob), 1000); M <- attr(pd, "M")
    v <- vapply(seq_len(nrow(pd)), function(r) pd[r, "psi"] * M * sum(exp(scr_lpi(pr, pd[r, ])[keep])), 1)
    v <- v / (sum(keep) * pr$pix^2) * scr_area_conv(pr)
    data.frame(estimate = mean(v), sd = stats::sd(v), lower = scr_q(v)[1], upper = scr_q(v)[2])
  }
  tab <- rbind(cbind(data.frame(fit = "main", buffer = cur, buffer_sigma = cur / sig), dreg(fit$object)),
               cbind(data.frame(fit = if (new > cur) "larger buffer" else "smaller buffer", buffer = new, buffer_sigma = new / sig), dreg(r2$object)))
  tab$region_buffer <- rad
  ch <- tab$estimate[2] / tab$estimate[1] - 1
  zs <- (tab$estimate[2] - tab$estimate[1]) / tab$sd[1]
  ok <- abs(ch) < 0.05 && abs(zs) < 0.3
  note <- sprintf("Buffer check: refitting with a %s state space (buffer %s = %.1f sigma instead of %s = %.1f sigma) changed density within %s %s of the detectors (the region both fits share) by %+.1f%% (%.2f posterior SD). %s",
                  if (new > cur) "larger" else "smaller", format(signif(new, 3)), new / sig, format(signif(cur, 3)), cur / sig,
                  format(signif(rad, 3)), scr_unit_word(prep$units), 100 * ch, zs,
                  if (ok) "The buffer is adequate: density does not depend on it." else
                    if (new > cur) "WARNING: density changed with the buffer: the original buffer is too small; use the larger one (>= 3-4 sigma)." else
                      "WARNING: density changed when the state space was trimmed by 2 sigma: check that the mask extends at least 3-4 sigma beyond the traps.")
  list(table = tab, note = note)
}

module <- ecomodule(
  id = "scr",
  title = "Spatial capture-recapture (closed population density)",
  description = paste(
    "How many animals live in the study area, at what density, and is density associated with habitat?",
    "Spatial capture-recapture uses where individually identified animals (camera-trapped by markings, DNA from hair",
    "snares, scats or rub objects, live traps) were detected across an array of detectors to estimate population density per km2,",
    "abundance in the habitat, baseline detection (p0, or the expected number of samples lambda0 for count data), the spatial",
    "scale of movement (sigma, from which a 95% home-range area follows), sex ratio and sex differences in detection and ranging,",
    "and the effects of habitat covariates on density and of detector, occasion and behavioural covariates on detection,",
    "accounting for animals that were never detected and for uneven detector effort."),
  defaults = list(
    params = list(density = 0.3, baseline_detection = 0.15, sigma = 1000,
                  density_beta_canopy = 0.5, detection_beta_on_trail = 1.0,
                  baseline_encounter_rate = 0.3, density_beta_elevation = -0.4, detection_beta_rain = -0.5,
                  sigma_beta_lure = 0.3, detection_beta_trap_response = 0.8,
                  proportion_female = 0.55, baseline_detection_male = 0.10, baseline_encounter_rate_male = 0.2, sigma_male = 1500),
    design = list(trap_grid_side = 7, trap_spacing = 1500, n_occasions = 6, buffer = 4000, mask_spacing = 500,
                  density_covariate = TRUE, detection_covariate = TRUE, prop_on_trail = 0.5,
                  trap_failure = 0.1, covariate_range = 2000, covariate_trend = 0,
                  obs_model = "binomial", density_covariate2 = FALSE, nonhabitat = 0, continuous_effort = FALSE,
                  prop_inactive = 0.1, occasion_covariate = FALSE, behaviour = FALSE, sigma_covariate = FALSE,
                  sex = FALSE, prop_sex_unknown = 0.2)),
  simulate = function(params, design, seed) { s <- scr_simulate(params, design, seed); s[c("data", "truth")] },
  fit = scr_fit, settings = scr_settings, fast_settings = scr_fast,
  read_data = scr_read_data, check_data = scr_check_data, write_data = scr_write_data,
  interpret = scr_interpret, plots = function(fit, data) scr_figures(fit, data, fig_opts()),
  figures = scr_figures, methods_text = scr_methods_text,
  power_params = c("density_beta_canopy", "density_beta_elevation", "detection_beta_on_trail", "detection_beta_rain",
                   "sigma_beta_lure", "detection_beta_trap_response", "detection_beta_male", "sigma_beta_male"),
  engine = "nimble (marginal: all variants; latent: basic binomial), numpyro (marginal NUTS: binomial/Poisson, mask, density and detector covariates), secr (maximum likelihood; all variants)",
  references = scr_references)
