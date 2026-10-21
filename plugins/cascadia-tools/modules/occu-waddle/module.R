# occu-waddle: two-species conditional ("dominant/subordinate") occupancy
# (Richmond, Hines & Beissinger 2010; Waddle et al. 2010).
# Species A (dominant) occupies a site with psiA; species B (subordinate)
# with psiBA where A is present and psiBa where A is absent. Detection of A is
# pA (B absent) or rA (both present); detection of B is pB (A absent), or,
# where both occur, rBA on surveys that detected A and rBa on surveys that
# missed A. Every probability is logit-linear in covariates.
# Engines: Nimble (Bayesian, marginal likelihood as a custom distribution;
# default), hand-written likelihood + optim (maximum likelihood; independent
# cross-check), and "naive" (ignores detection; negative control only).

for (f in c("design.R", "simulate.R", "engines.R", "gof.R", "io.R", "report.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

# ----------------------------------------------------------------- estimates
wd_estimates <- function(des, eng, point = "median") {
  occ_only <- eng$engine == "naive"
  if (eng$engine == "nimble") {
    ch <- lapply(eng$chains, function(m) wd_report_cols(des, m))
    s <- summarise_draws_eco(ch)
    est <- s$estimates
    # point estimate: posterior median (default) or mean (settings$point)
    if (!identical(point, "mean")) {
      D <- do.call(rbind, ch)
      est$estimate <- apply(D[, est$param, drop = FALSE], 2, stats::median)
    }
    eng$diag$rhat_max <- s$diagnostics$rhat_max; eng$diag$ess_min <- s$diagnostics$ess_min
    eng$diag$converged <- s$diagnostics$converged
  } else {
    pt <- wd_report_cols(des, matrix(eng$point, 1, dimnames = list(NULL, names(eng$point))), occ_only)
    dd <- wd_report_cols(des, eng$draws, occ_only)
    z <- stats::qnorm(0.975)
    rows <- lapply(colnames(pt), function(n) {
      k <- match(n, des$coef_names)
      if (!is.na(k)) {                       # link-scale slope: Wald
        return(data.frame(param = n, estimate = eng$point[k], sd = eng$se[k],
                          lower = eng$point[k] - z * eng$se[k], upper = eng$point[k] + z * eng$se[k]))
      }
      k <- match(paste0(n, "_intercept"), des$coef_names)
      if (!is.na(k)) {                       # probability: Wald on logit scale, back-transformed
        b <- eng$point[k]; se <- eng$se[k]
        return(data.frame(param = n, estimate = stats::plogis(b), sd = stats::plogis(b) * (1 - stats::plogis(b)) * se,
                          lower = stats::plogis(b - z * se), upper = stats::plogis(b + z * se)))
      }
      v <- dd[, n]                           # derived: point at MLE, interval from draws
      data.frame(param = n, estimate = pt[1, n], sd = stats::sd(v),
                 lower = stats::quantile(v, 0.025, names = FALSE), upper = stats::quantile(v, 0.975, names = FALSE))
    })
    est <- do.call(rbind, rows)
  }
  rownames(est) <- NULL
  list(estimates = est, diag = eng$diag)
}

# ----------------------------------------------------------------------- fit
wd_single_thread <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
  }
}

wd_fit <- function(data, settings) {
  wd_single_thread()
  settings <- merge_list(wd_settings, settings)
  set.seed(settings$seed)
  des <- wd_prepare(data)
  engine <- match.arg(settings$engine, c("nimble", "mle", "naive"))
  eng <- switch(engine, nimble = wd_engine_nimble(des, settings), mle = wd_engine_mle(des, settings),
                naive = wd_engine_naive(des, settings))
  e <- wd_estimates(des, eng, settings$point)
  dg <- e$diag
  dg$notes <- c(dg$notes, wd_data_warnings(des))
  tied <- wd_pars[des$map[wd_pars] != wd_pars]
  if (length(tied)) dg$notes <- c(dg$notes, sprintf("Constrained equal (no formula given): %s.",
    paste(sprintf("%s = %s", tied, des$map[tied]), collapse = ", ")))
  if (engine == "nimble" && is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
    dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                    dg$rhat_max, dg$ess_min))
  gof <- NULL
  if (engine != "naive" && isTRUE(settings$gof_ndraws > 0) && isTRUE(dg$converged))
    gof <- tryCatch(wd_gof(des, eng, settings), error = function(err) {
      dg$notes <<- c(dg$notes, paste("Goodness-of-fit check failed to run:", conditionMessage(err))); NULL })
  if (!is.null(gof)) dg$notes <- c(dg$notes, sprintf(
    "Goodness of fit (%s, %d draws): co-detection pattern p = %.2f; per-site detection counts p = %.2f%s.",
    gof$method, gof$ndraws, gof$p_pattern, gof$p_counts,
    if (min(gof$p_pattern, gof$p_counts) < 0.05 || max(gof$p_pattern, gof$p_counts) > 0.95)
      " -- LACK OF FIT: the model does not reproduce these data features; see README 'Goodness of fit'" else ""))
  dg$gof <- gof
  list(estimates = e$estimates, diagnostics = dg,
       object = list(engine = engine, des = des, draws = eng$draws, chains = eng$chains, point = eng$point,
                     vcov = eng$vcov, gof = gof, settings = settings, loglik = eng$diag$loglik, aic = eng$diag$aic))
}

wd_settings <- list(engine = "nimble", seed = 1, niter = 6000, nburnin = 1000, nchains = 3, thin = 1,
                    prior_sd = 1, sampler = "slice", point = "median", n_draws = 2000, mle_starts = 3, gof_ndraws = 500)

wd_fast <- merge_list(wd_settings, list(niter = 2500, nburnin = 500, gof_ndraws = 0))

wd_power_params <- c("occ_log_odds_ratio_sub", "species_interaction_factor", "det_log_odds_ratio_sub",
                     "occ_dominant_beta_habitat", "occ_sub_given_dom_present_beta_habitat",
                     "occ_sub_given_dom_absent_beta_habitat", "det_dominant_alone_beta_effort",
                     "det_sub_alone_beta_effort")

module <- ecomodule(
  id = "occu-waddle",
  title = "Two-species conditional occupancy (dominant/subordinate; Richmond et al. 2010, Waddle et al. 2010)",
  description = paste(
    "Does a subordinate species occupy sites less (or more) often where a dominant species is present than where it is absent,",
    "and does detecting one species depend on the other, after allowing for both species being missed on surveys?",
    "The model estimates occupancy of the dominant species, occupancy of the subordinate species conditional on the dominant being present or absent,",
    "and species-specific detection probabilities that may differ where both occur, with covariates on every probability.",
    "The dominant/subordinate direction is an assumption, not a result, and co-occurrence is not proof of competition or predation:",
    "shared unmeasured habitat preferences produce the same pattern."),
  defaults = list(
    params = list(occ_dominant = 0.5, occ_dominant_beta_habitat = 0.8,
                  occ_sub_given_dom_present = 0.3, occ_sub_given_dom_present_beta_habitat = 0.5,
                  occ_sub_given_dom_absent = 0.7, occ_sub_given_dom_absent_beta_habitat = -0.5,
                  det_dominant_alone = 0.5, det_dominant_alone_beta_effort = 0.4,
                  det_sub_alone = 0.5, det_sub_alone_beta_effort = 0.4,
                  det_dominant_with_sub = 0.5, det_sub_with_dom_detected = 0.35, det_sub_with_dom_missed = 0.55),
    design = list(n_sites = 200, n_surveys = 4, prop_missing = 0.05, det_interaction = TRUE, occ_interaction = TRUE)),
  simulate = wd_simulate, fit = wd_fit, settings = wd_settings, fast_settings = wd_fast,
  read_data = wd_read_data, check_data = wd_check_data,
  interpret = wd_interpret, plots = wd_plots, figures = wd_figures, methods_text = wd_methods_text,
  power_params = wd_power_params,
  null_values = c(stats::setNames(rep(0, length(wd_power_params)), wd_power_params)[wd_power_params != "species_interaction_factor"],
                  species_interaction_factor = 1),
  engine = "nimble (Bayesian, marginal likelihood; default) or hand-coded likelihood + optim (maximum likelihood)",
  references = c(
    "Richmond OMW, Hines JE, Beissinger SR (2010) Two-species occupancy models: a new parameterization applied to co-occurrence of secretive rails. Ecological Applications 20:2036-2046. doi:10.1890/09-0470.1",
    "Waddle JH, Dorazio RM, Walls SC, Rice KG, Beauchamp J, Schuman MJ, Mazzotti FJ (2010) A new parameterization for estimating co-occurrence of interacting species. Ecological Applications 20:1467-1475. doi:10.1890/09-0850.1",
    "MacKenzie DI, Bailey LL, Nichols JD (2004) Investigating species co-occurrence patterns when species are detected imperfectly. Journal of Animal Ecology 73:546-555. doi:10.1111/j.0021-8790.2004.00828.x",
    "MacKenzie DI, Lombardi JV, Tewes ME (2021) A note on investigating co-occurrence patterns and dynamics for many species, with imperfect detection and a log-linear modeling parameterization. Ecology and Evolution 11:8507-8515. doi:10.1002/ece3.7604",
    "Meredith M (2022) wiqid: Quick and dirty estimates for wildlife populations. R package version 0.3.3 (function occ2sps).",
    "Blanchet FG, Cazelles K, Gravel D (2020) Co-occurrence is not evidence of ecological interactions. Ecology Letters 23:1050-1063. doi:10.1111/ele.13525",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487"))
