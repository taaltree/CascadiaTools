# abund-interactions: abundance-mediated species interactions (Twining,
# Augustine, Royle & Fuller 2025, Ecology 106:e4468), as used for the pine
# marten - squirrel system (Twining et al. 2022 Proc R Soc B; Twining et al.
# 2022 Global Change Biology).
#
# Dominant species A: Royle-Nichols abundance from repeated detection /
# non-detection (or counts of detection days per occasion).
# Subordinate species S: occupancy, with the dominant's LATENT site abundance
# N_A as a covariate, optionally interacting with a habitat covariate.
# Engines: nimble (default, joint Bayesian model), ml (exact marginal
# likelihood; independent cross-check and occupancy-mediated negative
# control), twostep (plug-in comparison: occuRN then occu with E[N|y]).

for (f in c("prep.R", "engines.R", "gof.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

# ------------------------------------------------------------------ simulate
abi_simulate <- function(params, design, seed) {
  set.seed(seed)
  M <- design$n_sites; J <- design$n_occasions; D <- design$days_per_occasion
  zstd <- function(x) (x - mean(x)) / stats::sd(x)
  woodland <- zstd(stats::rnorm(M)); broadleaf <- zstd(stats::rnorm(M))
  # effort: most occasions run all D days; a share lose some days (camera failure)
  K <- matrix(D, M, J)
  lost <- matrix(stats::runif(M * J) < design$prop_effort_loss, M, J)
  K[lost] <- sample.int(D, sum(lost), replace = TRUE) - 1L
  for (i in which(rowSums(K) == 0)) K[i, 1] <- D
  lam <- params$dominant_abundance * exp(params$lambda_beta_woodland * woodland)
  N <- stats::rpois(M, lam)
  slope <- params$psi_beta_dominant_abundance +
    (if (isTRUE(design$habitat_interaction)) params$psi_beta_dominant_abundance_x_broadleaf * broadleaf else 0)
  psi <- stats::plogis(stats::qlogis(params$subordinate_occupancy_no_dominant) +
                         params$psi_beta_woodland * woodland + params$psi_beta_broadleaf * broadleaf + slope * N)
  z <- stats::rbinom(M, 1, psi)
  pA <- 1 - (1 - params$dominant_detection_per_individual)^N
  yA <- matrix(stats::rbinom(M * J, as.vector(K), rep(pA, J)), M, J)
  yS <- matrix(stats::rbinom(M * J, as.vector(K), rep(z * params$subordinate_detection, J)), M, J)
  yA[K == 0] <- NA; yS[K == 0] <- NA
  hab <- isTRUE(design$habitat_interaction)
  data <- list(yA = yA, yS = yS, K = K, site_covs = data.frame(woodland = woodland, broadleaf = broadleaf),
               obs_covs = list(),
               formulas = list(lambda = "~ woodland", r = "~ 1", psi = "~ woodland + broadleaf", p = "~ 1"),
               modifier = if (hab) "broadleaf" else NULL,
               species = c(dominant = "Dominant species", subordinate = "Subordinate species"),
               site_ids = sprintf("S%03d", seq_len(M)))
  truth <- c(dominant_mean_abundance = params$dominant_abundance,
             lambda_beta_woodland = params$lambda_beta_woodland,
             dominant_detection_per_individual = params$dominant_detection_per_individual,
             subordinate_occupancy_no_dominant = params$subordinate_occupancy_no_dominant,
             psi_beta_woodland = params$psi_beta_woodland, psi_beta_broadleaf = params$psi_beta_broadleaf,
             psi_beta_dominant_abundance = params$psi_beta_dominant_abundance,
             if (hab) c(psi_beta_dominant_abundance_x_broadleaf = params$psi_beta_dominant_abundance_x_broadleaf),
             subordinate_detection = params$subordinate_detection,
             dominant_total_abundance = sum(N), dominant_prop_occupied = mean(N > 0),
             subordinate_prop_occupied = mean(z))
  list(data = data, truth = truth)
}

# ----------------------------------------------------------------- estimates
abi_derived_names <- c("dominant_total_abundance", "dominant_prop_occupied", "subordinate_prop_occupied")

# Interaction slope (log-odds per dominant individual) at chosen habitat values
abi_slope_points <- function(des) {
  if (!length(des$colnames$mod)) return(NULL)
  mod <- des$modifier; raw <- des$site_covs_raw[[mod]]
  if (is.numeric(raw)) {
    q <- stats::quantile(raw, c(0.1, 0.9), na.rm = TRUE)
    sc <- des$scaling[[mod]]; z <- if (is.null(sc)) q else (q - sc$center) / sc$scale
    data.frame(label = paste0("interaction_at_", c("low_", "high_"), abi_clean(mod)),
               desc = sprintf("%s = %s (%s percentile)", mod, signif(q, 3), c("10th", "90th")),
               stringsAsFactors = FALSE) |> cbind(Mod = I(matrix(z, ncol = 1)))
  } else {
    lv <- levels(factor(raw))
    X <- stats::model.matrix(stats::as.formula(paste("~", mod)), data.frame(stats::setNames(list(factor(lv, levels = lv)), mod)))[, -1, drop = FALSE]
    data.frame(label = paste0("interaction_in_", abi_clean(lv)), desc = paste(mod, "=", lv), stringsAsFactors = FALSE) |>
      cbind(Mod = I(X))
  }
}

abi_slope_draws <- function(des, draws) {
  sp <- abi_slope_points(des)
  if (is.null(sp)) return(NULL)
  gm <- draws[, paste0("gamma_mod[", seq_along(des$colnames$mod), "]"), drop = FALSE]
  out <- draws[, "gamma0"] + gm %*% t(unclass(sp$Mod))
  colnames(out) <- sp$label
  out
}

abi_estimates <- function(des, eng) {
  dg <- eng$diag
  if (eng$engine == "nimble") {
    ch <- lapply(eng$chains, function(m) {
      pub <- abi_to_public(des, m)
      der <- cbind(dominant_total_abundance = m[, "N_total"], dominant_prop_occupied = m[, "n_occ_dominant"] / des$M,
                   subordinate_prop_occupied = m[, "n_occ_subordinate"] / des$M)
      cbind(pub, abi_slope_draws(des, m), der)
    })
    s <- summarise_draws_eco(ch)
    est <- s$estimates
    # R-hat is undefined for a quantity constant within chains (e.g. all sites occupied); ignore NA
    dg$rhat_max <- s$diagnostics$rhat_max; dg$ess_min <- s$diagnostics$ess_min
    dg$converged <- s$diagnostics$converged
  } else {
    int <- abi_internal(des); pub <- abi_public(des)
    th <- eng$point[int]; se <- eng$se[int]
    rows <- lapply(seq_along(int), function(k) {
      tr <- switch(pub[k], dominant_mean_abundance = exp,
                   dominant_detection_per_individual = , subordinate_occupancy_no_dominant = ,
                   subordinate_detection = stats::plogis, identity)
      r <- wald_rows(th[k], se[k], pub[k], transform = tr)
      if (!identical(tr, identity)) r$sd <- stats::sd(tr(eng$draws[, int[k]]))
      r
    })
    est <- do.call(rbind, rows)
    sl <- abi_slope_draws(des, eng$draws)
    if (!is.null(sl)) {
      sp <- abi_slope_points(des)
      G <- cbind(1, unclass(sp$Mod)); idx <- c("gamma0", paste0("gamma_mod[", seq_along(des$colnames$mod), "]"))
      V <- eng$vcov[idx, idx, drop = FALSE]
      est <- rbind(est, wald_rows(drop(G %*% eng$point[idx]), sqrt(pmax(0, diag(G %*% V %*% t(G)))), sp$label))
    }
    if (!is.null(eng$N)) est <- rbind(est, derived_rows(cbind(
      dominant_total_abundance = rowSums(eng$N), dominant_prop_occupied = rowMeans(eng$N > 0),
      subordinate_prop_occupied = rowMeans(eng$z))))
    est$rhat <- NA_real_; est$ess_bulk <- NA_real_
  }
  rownames(est) <- NULL
  list(estimates = est, diag = dg)
}

# ----------------------------------------------------------------------- fit
abi_settings <- list(engine = "nimble", interaction = TRUE, standardize = TRUE, seed = 1,
                     niter = 20000, nburnin = 5000, nchains = 3, thin = 5, prior_sd = 2, block = TRUE,
                     Nmax = 60, n_draws = 2000, mediator = "abundance",
                     gof = TRUE, gof_ndraws = 500, compare_twostep = TRUE, keep_latent = 2000)

abi_fast <- merge_list(abi_settings, list(niter = 5000, nburnin = 1500, thin = 2, gof = FALSE,
                                          compare_twostep = FALSE, n_draws = 300, keep_latent = 0))

abi_fit <- function(data, settings) {
  settings <- merge_list(abi_settings, settings)
  set.seed(settings$seed)
  des <- abi_prepare(data, standardize = isTRUE(settings$standardize))
  engine <- match.arg(settings$engine, c("nimble", "ml", "twostep"))
  if (identical(settings$interaction, FALSE)) {
    if (engine != "nimble") stop("settings$interaction = FALSE (reduced model) is available for engine nimble only")
    des$colnames$mod <- character(); des$Mod <- NULL; des$modifier <- NULL
    des$notes <- c(des$notes, "REDUCED MODEL: interaction fixed at 0; the subordinate part is a single-species occupancy model (psi_beta_dominant_abundance is fixed, not estimated).")
  }
  eng <- switch(engine, nimble = abi_engine_nimble(des, settings), ml = abi_engine_ml(des, settings),
                twostep = abi_engine_twostep(des, settings))
  e <- abi_estimates(des, eng)
  dg <- e$diag
  dg$notes <- c(data$warnings, des$notes, dg$notes)
  if (length(des$colnames$mod) && des$M < 1000)
    dg$notes <- c(dg$notes, sprintf("SMALL DESIGN for a habitat-dependent interaction: %d sites. In simulation checks (VALIDATION.md) the interaction terms were recovered without bias at 600 sites in the default design, but in an independent design with 700 sites (more dominant individuals per site, weaker interaction) the interaction at average habitat was overstated by ~12%%, and at 300 sites the habitat x interaction term by ~10-20%%. Intervals kept near-nominal coverage in all of these. Read the credible intervals, treat the size of the interaction as possibly exaggerated, and run a planted-truth check at your own design.", des$M))
  else if (des$M < 300)
    dg$notes <- c(dg$notes, sprintf("SMALL DESIGN: %d sites. Simulation checks (VALIDATION.md) found the interaction overstated by ~15-20%% at 150 sites; treat its size as possibly exaggerated.", des$M))
  if (engine == "ml" && eng$mediator == "occupancy")
    dg$notes <- c(dg$notes, "OCCUPANCY-MEDIATED CONTROL: subordinate occupancy depends only on whether the dominant is present (1[N > 0]); psi_beta_dominant_abundance is then the effect of presence, not of each individual.")
  if (engine == "nimble" && is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
    dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter / nburnin / thin).",
                                    dg$rhat_max, dg$ess_min))
  gof <- NULL
  if (isTRUE(settings$gof) && !is.null(eng$N)) {
    gof <- abi_gof(des, eng, settings$gof_ndraws)
    dg$notes <- c(dg$notes, sprintf("Posterior predictive check (Freeman-Tukey on site detection totals): Bayesian p = %.2f (%s), %.2f (%s); values near 0 or 1 (< 0.05 or > 0.95) indicate lack of fit.",
                                    gof$p[["dominant"]], des$species[["dominant"]], gof$p[["subordinate"]], des$species[["subordinate"]]))
  }
  ts <- NULL
  if (isTRUE(settings$compare_twostep) && engine != "twostep") {
    ts <- tryCatch(abi_estimates(des, abi_engine_twostep(des, settings))$estimates,
                   error = function(e) { dg$notes <<- c(dg$notes, paste("Two-step comparison failed:", conditionMessage(e))); NULL })
  }
  keep <- function(m) if (is.null(m) || settings$keep_latent <= 0) NULL else m[unique(round(seq(1, nrow(m), length.out = min(nrow(m), settings$keep_latent)))), , drop = FALSE]
  dg$gof <- gof
  list(estimates = e$estimates, diagnostics = dg,
       object = list(engine = engine, des = des, draws = keep(eng$draws) %||% eng$draws,
                     N = keep(eng$N), z = keep(eng$z), chains = if (settings$keep_latent > 0) eng$chains,
                     point = eng$point, vcov = eng$vcov, gof = gof, twostep = ts, settings = settings,
                     bup = eng$bup))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

module <- ecomodule(
  id = "abund-interactions",
  title = "Abundance-mediated species interactions (dominant abundance -> subordinate occupancy)",
  description = paste(
    "Is the occurrence of a subordinate species (e.g. red squirrel) associated with how MANY individuals of a dominant species",
    "(e.g. pine marten) use a site, and does that association change with habitat? The dominant's site abundance is estimated",
    "from repeated detection/non-detection (camera days) with a Royle-Nichols model, and enters the subordinate's occupancy model",
    "as a latent covariate, so its uncertainty is carried into the interaction estimate (Twining et al. 2025)."),
  defaults = list(
    params = list(dominant_abundance = 1, lambda_beta_woodland = 0.5, dominant_detection_per_individual = 0.1,
                  subordinate_occupancy_no_dominant = 0.6, psi_beta_woodland = 0.5, psi_beta_broadleaf = 0.3,
                  psi_beta_dominant_abundance = -0.7, psi_beta_dominant_abundance_x_broadleaf = 0.5,
                  subordinate_detection = 0.15),
    design = list(n_sites = 600, n_occasions = 4, days_per_occasion = 7, prop_effort_loss = 0.1,
                  habitat_interaction = TRUE)),
  simulate = abi_simulate, fit = abi_fit, settings = abi_settings, fast_settings = abi_fast,
  read_data = abi_read_data, check_data = abi_check_data,
  interpret = abi_interpret, figures = abi_figures, plots = abi_plots, methods_text = abi_methods_text,
  power_params = c("psi_beta_dominant_abundance", "psi_beta_dominant_abundance_x_broadleaf",
                   "lambda_beta_woodland", "psi_beta_woodland", "psi_beta_broadleaf"),
  engine = "nimble (joint Bayesian); ml (exact marginal likelihood); twostep (unmarked occuRN + occu plug-in, comparison)",
  references = c(
    "Twining JP, Augustine BC, Royle JA, Fuller AK (2025) Abundance-mediated species interactions. Ecology 106(1):e4468. doi:10.1002/ecy.4468",
    "Twining JP, Sutherland C, Reid N, Tosh DG (2022) Habitat mediates coevolved but not novel species interactions. Proceedings of the Royal Society B 289:20212338. doi:10.1098/rspb.2021.2338",
    "Twining JP, Lawton C, White A, Sheehy E, Hobson K, Montgomery WI, Lambin X (2022) Restoring vertebrate predator populations can provide landscape-scale biological control of established invasive vertebrates: insights from pine marten recovery in Europe. Global Change Biology 28:5368-5384. doi:10.1111/gcb.16236",
    "Royle JA, Nichols JD (2003) Estimating abundance from repeated presence-absence data or point counts. Ecology 84:777-790. doi:10.1890/0012-9658(2003)084[0777:EAFRPA]2.0.CO;2",
    "Waddle JH, Dorazio RM, Walls SC, Rice KG, Beauchamp J, Schuman MJ, Mazzotti FJ (2010) A new parameterization for estimating co-occurrence of interacting species. Ecological Applications 20:1467-1475. doi:10.1890/09-0850.1",
    "Rota CT, Ferreira MAR, Kays RW, Forrester TD, Kalies EL, McShea WJ, Parsons AW, Millspaugh JJ (2016) A multispecies occupancy model for two or more interacting species. Methods in Ecology and Evolution 7:1164-1173. doi:10.1111/2041-210X.12587",
    "Fiske I, Chandler R (2011) unmarked: An R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23. doi:10.18637/jss.v043.i10",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487"))
