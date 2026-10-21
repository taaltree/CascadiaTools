# scr-smr: spatial mark-resight (SMR; Chandler & Royle 2013; Sollmann et al.
# 2013) and generalized spatial mark-resight (gSMR; Whittington et al. 2018;
# Jimenez et al. 2019) for camera-trap resighting of partially marked
# populations, with optional telemetry of marked animals.
#
# Engines (settings$engine):
#   nimble  Bayesian data augmentation with latent activity centres (default)
#   secr    maximum (pseudo)likelihood mark-resight (Efford & Hunter 2018);
#           independent cross-check for the cases secr supports (SMR with
#           known marks; gSMR with marking + sighting occasions)

for (f in c("prepare.R", "simulate.R", "nimble_engine.R", "secr_engine.R", "post.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

smr_settings <- list(
  engine = "nimble", model = "auto", resight_model = "poisson", seed = 1,
  niter = 30000, nburnin = 6000, nchains = 3, thin = 4,
  M_aug = NULL, sigma_max = NULL, lam0_max = 10,
  density_covariates = NULL, detection_covariates = NULL,
  use_unidentified = TRUE, use_telemetry = TRUE,
  max_draws_per_chain = 3000, surface_draws = 300, ppc_draws = 300,
  compare_engines = character(), max_mask_points = 2500, secr_nsim = 50)

smr_fast <- merge_list(smr_settings, list(niter = 9000, nburnin = 3000, thin = 2,
                                          max_draws_per_chain = 1500, surface_draws = 0, ppc_draws = 0))

smr_bayes <- function(prep, settings) {
  st <- settings
  st$.monitor_latent <- max(settings$surface_draws, settings$ppc_draws) > 0
  for (attempt in 1:3) {
    eng <- smr_engine_nimble(prep, st)
    rep_chains <- lapply(eng$chains, function(mm) smr_report_chain(prep, mm, smr_thin_idx(nrow(mm), settings$max_draws_per_chain)))
    Nall <- unlist(lapply(rep_chains, function(mm) mm[, "abundance"]))
    cap <- if (prep$has_marking) eng$M else eng$M - prep$m
    Nq <- stats::quantile(if (prep$has_marking) Nall else Nall - prep$m, 0.995)
    hit <- Nq >= 0.9 * cap
    if (!hit || attempt == 3) break
    st$M_aug <- 2 * eng$naug
    prep$notes <- c(prep$notes, sprintf("Data augmentation limit (M = %d) was too small (99.5%% quantile of N = %.0f); refitted with M = %d.",
                                        eng$M, stats::quantile(Nall, 0.995), eng$MK + st$M_aug))
  }
  s <- summarise_draws_eco(rep_chains)
  list(eng = eng, prep = prep, rep_chains = rep_chains, estimates = s$estimates, diag = s$diagnostics, M_hit = hit)
}

smr_fit <- function(data, settings) {
  settings <- merge_list(smr_settings, settings)
  set.seed(settings$seed)
  prep <- smr_prepare(data, settings)
  if (prep$n_id_sightings + prep$n_unmarked_sightings == 0) stop("no sightings at all: nothing to estimate")
  if (settings$engine == "secr") {
    e <- smr_engine_secr(prep, settings)
    dg <- list(converged = e$converged, notes = c(prep$notes, e$notes), t_total = e$t_total)
    obj <- list(engine = "secr", prep = prep, secr = e$fit, settings = settings)
    est <- e$estimates
  } else {
    b <- smr_bayes(prep, settings)
    prep <- b$prep
    est <- b$estimates
    dg <- b$diag; dg$t_total <- b$eng$t_total; dg$notes <- prep$notes
    if (b$M_hit) { dg$converged <- FALSE
      dg$notes <- c(dg$notes, "Abundance still presses against the data-augmentation limit M: set settings$M_aug larger.") }
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                      dg$rhat_max, dg$ess_min))
    M <- b$eng$M
    surf <- ppc <- NULL
    if (max(settings$surface_draws, settings$ppc_draws) > 0) {
      pool <- smr_pool(b$eng$chains, max(settings$surface_draws, settings$ppc_draws))
      if (settings$surface_draws > 0) surf <- smr_surfaces(prep, pool[smr_thin_idx(nrow(pool), settings$surface_draws), , drop = FALSE], M)
      if (settings$ppc_draws > 0) ppc <- smr_ppc(prep, pool[smr_thin_idx(nrow(pool), settings$ppc_draws), , drop = FALSE], M)
    }
    if (!is.null(ppc)) dg$notes <- c(dg$notes, sprintf(
      "Posterior predictive checks (%d draws; P(replicate >= observed), extreme if < 0.025 or > 0.975): %s.",
      nrow(ppc$rep), paste(sprintf("%s %.2f", names(ppc$p), ppc$p), collapse = ", ")))
    # keep only top-level columns of the raw chains (s/z are large)
    keep <- function(mm) mm[, !grepl("^(s|z)\\[", colnames(mm)), drop = FALSE]
    obj <- list(engine = paste("nimble", prep$model), prep = prep, chains = lapply(b$eng$chains, keep), rep_chains = b$rep_chains,
                surfaces = surf, ppc = ppc, settings = settings, M = M)
  }
  dg$notes <- c(dg$notes, smr_design_notes(prep, est))
  if (length(settings$compare_engines)) {
    cmp <- list(data.frame(engine = obj$engine, est[, c("param", "estimate", "sd", "lower", "upper")]))
    for (e2 in setdiff(unlist(settings$compare_engines), settings$engine)) {
      s2 <- settings; s2$engine <- e2; s2$compare_engines <- character(); s2$ppc_draws <- 0; s2$surface_draws <- 0
      r2 <- tryCatch(smr_fit(data, s2), error = function(e) { message("compare engine ", e2, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r2)) cmp[[length(cmp) + 1]] <- data.frame(engine = r2$object$engine, r2$estimates[, c("param", "estimate", "sd", "lower", "upper")])
    }
    obj$comparison <- do.call(rbind, cmp)
  }
  list(estimates = est, diagnostics = dg, object = obj)
}

module <- ecomodule(
  id = "scr-smr",
  title = "Spatial mark-resight (SMR and generalized SMR)",
  description = paste(
    "How many animals live in the study area, and at what density, when only some of them carry marks (collars, ear tags, natural",
    "marks on a subset) and cameras photograph both marked and unmarked animals? Spatial mark-resight combines where each marked",
    "animal was resighted, the counts of unmarked animals at each camera, sightings of marked animals whose mark could not be read,",
    "optional telemetry of the marked animals and (generalized SMR) the spatial marking process itself, to estimate density per km2,",
    "abundance, the marked fraction, the spatial scale of movement (sigma), resighting and marking rates and habitat effects on density."),
  defaults = list(
    params = list(density = 0.25, resight_rate = 0.3, sigma = 1000, id_probability = 0.8,
                  marking_detection = 0.3, density_beta_canopy = 0.5, detection_beta_on_trail = 0.7),
    design = list(model = "smr", resight_model = "poisson", n_marked = 20,
                  trap_grid_side = 7, trap_spacing = 1500, n_occasions = 10, buffer = 4000, mask_spacing = 500,
                  density_covariate = TRUE, detection_covariate = TRUE, prop_on_trail = 0.5,
                  trap_failure = 0.1, covariate_range = 2000,
                  telemetry_fixes = 0, telemetry_fraction = 1,
                  marking_grid_side = 5, marking_spacing = 1500, marking_occasions = 5, marking_offset = c(0, 0))),
  simulate = function(params, design, seed) { s <- smr_simulate(params, design, seed); s[c("data", "truth")] },
  fit = smr_fit, settings = smr_settings, fast_settings = smr_fast,
  read_data = smr_read_data, check_data = smr_check_data, write_data = smr_write_data,
  interpret = smr_interpret, plots = function(fit, data) smr_figures(fit, data, fig_opts()),
  figures = smr_figures, methods_text = smr_methods_text,
  power_params = c("density_beta_canopy", "detection_beta_on_trail"),
  engine = "nimble (Bayesian data augmentation, latent activity centres); secr (maximum pseudo-likelihood; cross-check)",
  references = smr_references)
