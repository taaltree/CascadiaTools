# scr-search: search-encounter spatial capture-recapture for noninvasive
# genetic sampling along search tracks (scat-detection dogs, snow-tracking,
# transect searches for hair/scat). Royle, Kery & Guelat 2011; Thompson,
# Royle & Garner 2012; Russell et al. 2012; Efford 2011; Milleret et al. 2019;
# Moqanaki et al. 2021.
#
# Search tracks are discretised into square detector cells; effort = length
# of each track inside each cell. Expected samples of an animal centred at s
# in search unit u (cell x track):
#   lambda = lambda0 * E_u * exp(w_u' alpha) * exp(-d(s, cell)^2 / (2 sigma^2))
# Engines (settings$engine):
#   nimble        Bayesian, marginal over the habitat mask, local evaluation (default)
#   secr_transect maximum likelihood, tracks as secr transect detectors with
#                 exact sample locations (independent cross-check)
#   secr_cells    maximum likelihood of the same discretised model (likelihood check)

for (f in c("prepare.R", "simulate.R", "nimble_engine.R", "secr_engine.R", "post.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

ss_settings <- list(
  engine = "nimble", observation = "poisson", effort_model = "track_length", detector_location = "centroid", sampler = "RW_block", seed = 1,
  niter = 40000, nburnin = 10000, nchains = 3, thin = 1,
  max_draws_per_chain = 1500, surface_draws = 150, ppc_draws = 200,
  compare_engines = character(), M = NULL, sigma_max = NULL, local_radius = NULL, no_local = FALSE,
  cell_size = NULL, density_covariates = NULL, detection_covariates = NULL, max_mask_points = 3000,
  secr_method = "Newton-Raphson")

ss_fast <- merge_list(ss_settings, list(niter = 20000, nburnin = 5000, max_draws_per_chain = 800,
                                        surface_draws = 0, ppc_draws = 0))

ss_bayes <- function(prep, settings, data) {
  for (attempt in 1:4) {
    eng <- ss_engine_nimble(prep, settings)
    prep$M_used <- prep$M
    mx <- settings$max_draws_per_chain
    rep_chains <- lapply(eng$chains, function(m) ss_report_chain(prep, m, ss_thin_idx(nrow(m), mx)))
    Nq <- stats::quantile(unlist(lapply(rep_chains, function(m) m[, "abundance"])), 0.995)
    sq <- stats::quantile(unlist(lapply(eng$chains, function(m) m[, "sigma"])), 0.975)
    M_hit <- Nq >= 0.9 * prep$M
    R_hit <- is.finite(prep$R) && 4 * sq > prep$R
    if ((!M_hit && !R_hit) || attempt == 4) break
    if (M_hit) {
      prep$notes <- c(prep$notes, sprintf("Data augmentation limit M = %d was too small (99.5%% quantile of N = %.0f); refitted with M = %d.", prep$M, Nq, 2 * prep$M))
      settings$M <- 2 * prep$M
    }
    if (R_hit) {
      newR <- 1.5 * 4 * sq * prep$cscale
      prep$notes <- c(prep$notes, sprintf("Local-evaluation radius %s %s was less than 4 x the 97.5%% quantile of sigma; refitted with radius %s.",
                                          format(signif(prep$R * prep$cscale, 3)), ss_unit_word(prep$units), format(signif(newR, 3))))
      settings$local_radius <- newR
    }
    old_notes <- prep$notes
    prep <- ss_prepare(data, settings); prep$notes <- unique(c(old_notes, prep$notes))
  }
  s <- summarise_draws_eco(rep_chains)
  list(eng = eng, prep = prep, rep_chains = rep_chains, estimates = s$estimates, diag = s$diagnostics,
       M_hit = M_hit, R_hit = R_hit)
}

ss_fit <- function(data, settings) {
  settings <- merge_list(ss_settings, settings)
  set.seed(settings$seed)
  prep <- ss_prepare(data, settings)
  prep$notes <- c(data$read_notes, prep$notes)
  if (prep$n == 0) stop("no individuals were identified from the samples: density cannot be estimated from zero detections")
  engine <- settings$engine
  if (engine %in% c("secr_transect", "secr_cells")) {
    e <- if (engine == "secr_transect") ss_engine_secr_transect(prep, settings) else ss_engine_secr_cells(prep, settings)
    est <- e$estimates
    dg <- list(converged = e$converged, notes = prep$notes, t_total = e$t_total)
    obj <- list(engine = engine, prep = prep, secr = e$fit, settings = settings)
  } else if (engine == "nimble") {
    b <- ss_bayes(prep, settings, data)
    prep <- b$prep
    est <- b$estimates
    dg <- b$diag; dg$t_total <- b$eng$t_total; dg$notes <- prep$notes
    # posterior skew of density / N (independent verification: posterior MEANS ran
    # +9 to +12% high at designs where the density posterior CV was >= ~0.23)
    dn <- ss_density_name(prep$units); dr <- est[est$param == dn, ]
    if (nrow(dr) && is.finite(dr$sd) && dr$sd / dr$estimate > 0.2) {
      dmed <- stats::median(unlist(lapply(b$rep_chains, function(m) m[, dn])))
      nmed <- stats::median(unlist(lapply(b$rep_chains, function(m) m[, "abundance"])))
      dg$notes <- c(dg$notes, sprintf(paste(
        "NOTE: density is imprecise (posterior SD = %.0f%% of the posterior mean), and its posterior is right-skewed: posterior median density %s and N %s vs posterior means %s and %s.",
        "In recovery simulations with similar precision (posterior CV 0.23-0.28, 20-35 animals identified) the posterior MEAN of density and N was 9-12%% above the truth while secr's ML estimate was within 4%% (modules/scr-search/VALIDATION.md).",
        "Report the 95%% interval, and the posterior median or the secr_cells estimate (settings compare_engines: [secr_cells]) alongside the mean."),
        100 * dr$sd / dr$estimate, format(signif(dmed, 3)), format(round(nmed, 1)), format(signif(dr$estimate, 3)),
        format(round(est$estimate[est$param == "abundance"], 1))))
    }
    if (b$M_hit) { dg$converged <- FALSE
      dg$notes <- c(dg$notes, "Abundance still presses against the data-augmentation limit M: set settings$M larger.") }
    if (b$R_hit) dg$notes <- c(dg$notes, "WARNING: sigma is large relative to the local-evaluation radius even after refitting; set settings$no_local: true.")
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                      dg$rhat_max, dg$ess_min))
    nd <- max(settings$surface_draws, settings$ppc_draws, 1)
    pool <- ss_pool(b$eng$chains, nd)
    surf <- if (settings$surface_draws > 0) ss_surfaces(prep, pool[ss_thin_idx(nrow(pool), settings$surface_draws), , drop = FALSE], prep$M_used) else NULL
    ppc <- if (settings$ppc_draws > 0) ss_ppc(prep, pool[ss_thin_idx(nrow(pool), settings$ppc_draws), , drop = FALSE], prep$M_used) else NULL
    if (!is.null(ppc)) dg$notes <- c(dg$notes, sprintf(
      "Posterior predictive checks (%d draws; P(replicate >= observed), extreme if < 0.025 or > 0.975): %s.",
      nrow(ppc$rep), paste(sprintf("%s %.2f", names(ppc$p), ppc$p), collapse = ", ")))
    dg$ppc_p <- if (!is.null(ppc)) ppc$p
    obj <- list(engine = "nimble", prep = prep, chains = b$eng$chains, rep_chains = b$rep_chains,
                surfaces = surf, ppc = ppc, settings = settings)
  } else stop("unknown engine '", engine, "' (nimble, secr_transect or secr_cells)")
  dg$notes <- c(dg$notes, ss_design_notes(prep, est))
  if (length(settings$compare_engines)) {
    cmp <- list(data.frame(engine = obj$engine, est[, c("param", "estimate", "sd", "lower", "upper")]))
    for (e2 in setdiff(unlist(settings$compare_engines), engine)) {
      s2 <- settings; s2$engine <- e2; s2$compare_engines <- character(); s2$ppc_draws <- 0; s2$surface_draws <- 0
      r2 <- tryCatch(ss_fit(data, s2), error = function(e) { message("compare engine ", e2, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r2)) cmp[[length(cmp) + 1]] <- data.frame(engine = r2$object$engine, r2$estimates[, c("param", "estimate", "sd", "lower", "upper")])
    }
    obj$comparison <- do.call(rbind, cmp)
  }
  list(estimates = est, diagnostics = dg, object = obj)
}

module <- ecomodule(
  id = "scr-search",
  title = "Search-encounter SCR (scat-dog, snow-track and transect searches for genetic samples)",
  description = paste(
    "How many animals live in the study area, at what density, and is density associated with habitat, when individuals",
    "are identified from noninvasive genetic samples (scats, hair) found along search routes (detection-dog teams, snow-tracking,",
    "transect searches) rather than at fixed traps? The search tracks are turned into detector cells whose effort is the",
    "length searched; a spatial capture-recapture model then estimates density per km2, abundance, the spatial scale of",
    "movement (sigma), the sample-encounter rate per km of search (lambda0) and effects of habitat on density and of",
    "search conditions (dog team, snow) on detection, accounting for uneven search effort and for animals never sampled."),
  defaults = list(
    params = list(density = 0.10, lambda0_per_km = 0.20, sigma = 1500,
                  density_beta_forest = 0.5, detection_beta_experienced_team = 0.6),
    design = list(search_side = 16000, buffer = 6000, mask_spacing = 700, cell_size = 1000,
                  n_occasions = 3, tracks_per_occasion = 12, track_km = c(4, 10),
                  n_roads = 3, road_bias = 0.3, cov_road_corr = 0, covariate_range = 3000,
                  density_covariate = TRUE, detection_covariate = TRUE)),
  simulate = function(params, design, seed) { s <- ss_simulate(params, design, seed); s[c("data", "truth")] },
  fit = ss_fit, settings = ss_settings, fast_settings = ss_fast,
  read_data = ss_read_data, check_data = ss_check_data, write_data = ss_write_data,
  interpret = ss_interpret, plots = function(fit, data) ss_figures(fit, data, fig_opts()),
  figures = ss_figures, methods_text = ss_methods_text,
  power_params = c("density_beta_forest", "detection_beta_experienced_team"),
  engine = "nimble (Bayesian marginal likelihood with local evaluation), secr (ML: transect detectors with exact locations; discretised cells)",
  references = ss_references)
