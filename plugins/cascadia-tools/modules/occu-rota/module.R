# occu-rota: multispecies occupancy with species interactions
# (Rota et al. 2016, Methods Ecol Evol 7:1164-1173), 2 or 3 species.
# Latent occupancy state of the S species at a site is multivariate
# Bernoulli; first-order (f_a) and interaction (f_ab, f_abc) natural
# parameters take covariates; detection is species-specific with covariates.
# Engines: unmarked::occuMulti (maximum likelihood, default), Nimble
# (Bayesian, categorical latent state z_i in 1..2^S), and a "naive" control
# that ignores detection (negative control only).

for (f in c("design.R", "engines.R", "gof.R", "io.R", "report.R", "figures.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

# ------------------------------------------------------------------ simulate
# Parameter vectors are indexed by species (occ_*, detection_*, p_*) or by
# species pair in combn order (1:2, 1:3, 2:3) for interaction_*; a scalar
# is recycled to every species / pair.
rota_pv <- function(x, k) x[min(k, length(x))]

rota_simulate <- function(params, design, seed) {
  set.seed(seed)
  S <- design$n_species; M <- design$n_sites; J <- design$n_surveys
  species <- if (!is.null(design$species)) design$species[seq_len(S)] else paste0("sp", seq_len(S))
  max_order <- min(if (is.null(design$max_order)) 2 else design$max_order, S)
  zstd <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(as.vector(x), na.rm = TRUE)
  w <- matrix(stats::rbinom(M * J, 1, 1 - design$prop_missing), M, J)
  for (i in seq_len(M)) if (sum(w[i, ]) == 0) w[i, sample.int(J, 1)] <- 1
  forest <- zstd(stats::rnorm(M))
  effort <- matrix(stats::rnorm(M * J), M, J); effort[w == 0] <- NA; effort <- zstd(effort)
  terms <- rota_terms(species, max_order)
  Fm <- matrix(0, M, length(terms)); truth <- c()
  for (k in seq_along(terms)) {
    cc <- terms[[k]]
    if (length(cc) == 1) {
      a <- rota_pv(params$occ_intercept, cc); b <- rota_pv(params$occ_beta_forest, cc)
      Fm[, k] <- a + b * forest
      truth[paste0("occ_", species[cc], "_intercept")] <- a
      truth[paste0("occ_", species[cc], "_beta_forest")] <- b
    } else if (length(cc) == 2) {
      pair <- match(paste(cc, collapse = ":"), c("1:2", "1:3", "2:3"))
      a <- rota_pv(params$interaction, pair)
      b <- if (isTRUE(design$interaction_covariate)) rota_pv(params$interaction_beta_forest, pair) else 0
      Fm[, k] <- a + b * forest
      stem <- paste0("interaction_", paste(species[cc], collapse = "_"))
      truth[stem] <- a
      if (isTRUE(design$interaction_covariate)) truth[paste0(stem, "_beta_forest")] <- b
    } else {
      Fm[, k] <- params$interaction3
      truth[paste0("interaction_", paste(species, collapse = "_"))] <- params$interaction3
    }
  }
  Z <- rota_states(species)
  A <- rota_Amat(Z, terms)
  Psi <- rota_psi(Fm, A)
  state <- apply(Psi, 1, function(pr) sample.int(nrow(Z), 1, prob = pr))
  zs <- Z[state, , drop = FALSE]
  y <- array(NA_real_, c(M, J, S))
  for (s in seq_len(S)) {
    p0 <- rota_pv(params$detection_prob, s); bp <- rota_pv(params$p_beta_effort, s)
    p <- stats::plogis(stats::qlogis(p0) + bp * effort)
    ys <- matrix(stats::rbinom(M * J, 1, ifelse(is.na(p), 0, p) * zs[, s]), M, J)
    ys[w == 0] <- NA
    y[, , s] <- ys
    truth[paste0("detection_prob_", species[s])] <- p0
    truth[paste0("p_", species[s], "_beta_effort")] <- bp
  }
  # biological-scale truth at average forest (the standardised covariate = 0)
  F0 <- matrix(vapply(seq_along(terms), function(k) {
    cc <- terms[[k]]
    if (length(cc) == 1) rota_pv(params$occ_intercept, cc)
    else if (length(cc) == 2) rota_pv(params$interaction, match(paste(cc, collapse = ":"), c("1:2", "1:3", "2:3")))
    else params$interaction3
  }, 0), 1)
  der <- rota_derived(rota_psi(F0, A), Z)
  truth <- c(truth, stats::setNames(as.vector(der), colnames(der)))
  int_f <- if (isTRUE(design$interaction_covariate)) "~ forest" else "~ 1"
  data <- list(species = species, y = y, site_covs = data.frame(forest = forest),
               obs_covs = list(effort = effort),
               formulas = list(occupancy = "~ forest", interaction = int_f, detection = "~ effort"),
               max_order = max_order, site_ids = sprintf("S%03d", seq_len(M)))
  list(data = data, truth = truth)
}

# ----------------------------------------------------------------- estimates
# Transform a draws matrix (beta then alpha, link scale) into the reported
# columns: link coefficients (detection intercepts -> probability) and
# biological-scale quantities at average covariate values.
rota_report_cols <- function(des, D, has_det = TRUE) {
  P <- des$P
  B <- D[, seq_len(P), drop = FALSE]
  out <- B
  if (has_det) {
    Al <- D[, P + seq_len(des$Q), drop = FALSE]
    di <- grepl("^detection_prob_", colnames(Al))
    Al[, di] <- stats::plogis(Al[, di])
    out <- cbind(out, Al)
  }
  Fm <- rota_F_from(des, B, rota_Xnew(des, data.frame()))[[1]]
  Fm <- matrix(Fm, nrow(B))
  der <- rota_derived(rota_psi(Fm, des$A), des$Z)
  cbind(out, der)
}

rota_estimates <- function(des, eng) {
  has_det <- eng$engine != "naive"
  if (eng$engine == "nimble") {
    ch <- lapply(eng$chains, function(m) rota_report_cols(des, m, has_det))
    s <- summarise_draws_eco(ch)
    est <- s$estimates
    eng$diag$rhat_max <- s$diagnostics$rhat_max; eng$diag$ess_min <- s$diagnostics$ess_min
    eng$diag$converged <- s$diagnostics$converged
  } else {
    pt <- rota_report_cols(des, matrix(eng$point, 1, dimnames = list(NULL, names(eng$point))), has_det)
    np <- if (has_det) des$P + des$Q else des$P
    link <- wald_rows(eng$point[seq_len(np)], eng$se[seq_len(np)], names(eng$point)[seq_len(np)])
    di <- grepl("^detection_prob_", link$param)
    if (any(di)) {
      b <- eng$point[di]; se <- eng$se[di]
      link$estimate[di] <- stats::plogis(b); link$sd[di] <- stats::plogis(b) * (1 - stats::plogis(b)) * se
      link$lower[di] <- stats::plogis(b - 1.959964 * se); link$upper[di] <- stats::plogis(b + 1.959964 * se)
    }
    dd <- rota_report_cols(des, eng$draws, has_det)
    dn <- rota_derived_names(des$species)
    der <- derived_rows(dd[, dn, drop = FALSE])
    der$estimate <- as.vector(pt[1, dn])
    est <- rbind(link, der)
  }
  rownames(est) <- NULL
  list(estimates = est, diag = eng$diag)
}

# ----------------------------------------------------------------------- fit
# R here links multithreaded OpenBLAS. In forked recovery workers
# (parallel::mclapply) its thread pool oversubscribes the node or deadlocks
# (measured: 0.7 s/fit serial vs 86 s/fit or hang with 4 forks), so fits run
# single-threaded; the models are small and gain nothing from BLAS threads.
rota_single_thread <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
  }
}

rota_fit <- function(data, settings) {
  rota_single_thread()
  settings <- merge_list(rota_settings, settings)
  set.seed(settings$seed)
  des <- rota_prepare(data, standardize = isTRUE(settings$standardize))
  engine <- match.arg(settings$engine, c("unmarked", "nimble", "naive"))
  eng <- switch(engine,
                unmarked = rota_engine_unmarked(des, settings),
                nimble = rota_engine_nimble(des, settings),
                naive = rota_engine_naive(des, settings))
  e <- rota_estimates(des, eng)
  dg <- e$diag
  dg$notes <- c(dg$notes, rota_data_warnings(des))
  if (engine == "nimble" && is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
    dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                    dg$rhat_max, dg$ess_min))
  gof <- NULL
  if (engine != "naive" && isTRUE(settings$gof_nboot > 0) && isTRUE(dg$converged))
    gof <- tryCatch(rota_gof(des, eng, settings), error = function(err) {
      dg$notes <<- c(dg$notes, paste("Goodness-of-fit check failed to run:", conditionMessage(err))); NULL })
  if (!is.null(gof)) dg$notes <- c(dg$notes, sprintf(
    "Goodness of fit (%s, %d replicates): statistic %.1f, c-hat %.2f, p = %.3f%s.",
    gof$method, gof$nboot, gof$stat, gof$c_hat, gof$p_value,
    if (gof$p_value < 0.05) " -- LACK OF FIT: the model does not reproduce the detection data; see the goodness-of-fit section of the README" else ""))
  dg$gof <- gof
  list(estimates = e$estimates, diagnostics = dg,
       object = list(engine = engine, des = des, draws = eng$draws, chains = eng$chains, point = eng$point, vcov = eng$vcov,
                     gof = gof, fm = eng$fm, settings = settings, loglik = eng$diag$loglik))
}

rota_settings <- list(engine = "unmarked", standardize = TRUE, seed = 1,
                      n_draws = 2000, penalty = 0, penalty_boot = 30,
                      gof_nboot = 100, gof_ndraws = 300,
                      niter = 6000, nburnin = 1000, nchains = 3, thin = 1, prior_sd = 2,
                      sampler = "slice")

rota_fast <- merge_list(rota_settings, list(gof_nboot = 0, n_draws = 1000, niter = 3000, nburnin = 1000))

rota_power_params <- c(
  "interaction_sp1_sp2", "interaction_sp1_sp3", "interaction_sp2_sp3",
  "interaction_sp1_sp2_beta_forest", "interaction_sp1_sp3_beta_forest", "interaction_sp2_sp3_beta_forest",
  "interaction_sp1_sp2_sp3",
  "occ_sp1_beta_forest", "occ_sp2_beta_forest", "occ_sp3_beta_forest",
  "p_sp1_beta_effort", "p_sp2_beta_effort", "p_sp3_beta_effort")

module <- ecomodule(
  id = "occu-rota",
  title = "Multispecies occupancy with species interactions (Rota et al. 2016)",
  description = paste(
    "Do two or three species occur together at sites more or less often than expected if they occurred independently,",
    "and does that co-occurrence change along a habitat gradient, after allowing for each species being missed on surveys?",
    "The model (Rota et al. 2016) estimates each species' occupancy, the probability one species occupies a site",
    "given another species is present versus absent, and pairwise interaction (co-occurrence) parameters with covariates.",
    "Co-occurrence is not proof of competition, predation or facilitation: shared unmeasured habitat preferences produce the same pattern."),
  defaults = list(
    params = list(occ_intercept = c(-0.3, -0.5, -0.6), occ_beta_forest = c(0.8, -0.6, 0.4),
                  interaction = c(1.0, -0.8, 0.5), interaction_beta_forest = c(-0.7, 0.5, 0),
                  interaction3 = 0,
                  detection_prob = c(0.4, 0.35, 0.45), p_beta_effort = c(0.5, -0.4, 0.3)),
    design = list(n_species = 2, n_sites = 400, n_surveys = 4, max_order = 2, prop_missing = 0.1,
                  interaction_covariate = TRUE)),
  simulate = rota_simulate, fit = rota_fit, settings = rota_settings, fast_settings = rota_fast,
  read_data = rota_read_data, check_data = rota_check_data,
  interpret = rota_interpret, plots = rota_plots, figures = rota_figures, methods_text = rota_methods_text,
  power_params = rota_power_params,
  engine = "unmarked::occuMulti (maximum likelihood, default) or nimble (Bayesian, categorical latent state)",
  references = c(
    "Rota CT, Ferreira MAR, Kays RW, Forrester TD, Kalies EL, McShea WJ, Parsons AW, Millspaugh JJ (2016) A multispecies occupancy model for two or more interacting species. Methods in Ecology and Evolution 7:1164-1173. doi:10.1111/2041-210X.12587",
    "Fiske I, Chandler R (2011) unmarked: An R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23. doi:10.18637/jss.v043.i10",
    "Kellner KF, Smith AD, Royle JA, Kery M, Belant JL, Chandler RB (2023) The unmarked R package: Twelve years of advances in occurrence and abundance modelling in ecology. Methods in Ecology and Evolution 14:1408-1415. doi:10.1111/2041-210X.14123",
    "Clipp HL, Evans AL, Kessinger BE, Kellner K, Rota CT (2021) A penalized likelihood for multispecies occupancy models improves predictions of species interactions. Ecology 102:e03520. doi:10.1002/ecy.3520",
    "MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA, Langtimm CA (2002) Estimating site occupancy rates when detection probabilities are less than one. Ecology 83:2248-2255. doi:10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2",
    "Blanchet FG, Cazelles K, Gravel D (2020) Co-occurrence is not evidence of ecological interactions. Ecology Letters 23:1050-1063. doi:10.1111/ele.13525",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487"))
