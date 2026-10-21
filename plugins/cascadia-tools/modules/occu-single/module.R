# occu-single: single-species occupancy, single-season (MacKenzie et al. 2002)
# and dynamic multi-season (MacKenzie et al. 2003), selected by data$model
# (config data$model / design$model = "static" | "dynamic").
# Engines: unmarked (default, maximum likelihood) and Nimble (Bayesian),
# plus a "naive" logistic regression that ignores detection (negative control only).

for (f in c("design.R", "engines.R", "gof.R", "io.R", "report.R", "figures.R"))
  source(file.path(MODULE_DIR, "R", f), local = TRUE)

# ------------------------------------------------------------------ simulate
occ_simulate <- function(params, design, seed) {
  set.seed(seed)
  model <- if (is.null(design$model)) "static" else design$model
  M <- design$n_sites; J <- design$n_surveys
  T <- if (model == "dynamic") design$n_seasons else 1
  zstd <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(as.vector(x), na.rm = TRUE)
  # missing surveys at random, at least one survey per site and season
  w <- array(stats::rbinom(M * J * T, 1, 1 - design$prop_missing), c(M, J, T))
  for (i in seq_len(M)) for (t in seq_len(T)) if (sum(w[i, , t]) == 0) w[i, sample.int(J, 1), t] <- 1
  forest <- zstd(stats::rnorm(M))
  effort <- array(stats::rnorm(M * J * T), c(M, J, T)); effort[w == 0] <- NA; effort <- zstd(effort)
  site_re <- if (!is.null(design$p_site_sd) && design$p_site_sd > 0) stats::rnorm(M, 0, design$p_site_sd) else rep(0, M)
  p <- stats::plogis(stats::qlogis(params$detection_prob) + params$p_beta_effort * effort + site_re)
  if (model == "static") {
    psi <- stats::plogis(stats::qlogis(params$occupancy_prob) + params$psi_beta_forest * forest)
    z <- matrix(stats::rbinom(M, 1, psi), M, 1)
  } else {
    drought <- matrix(stats::rnorm(M * T), M, T)   # season t value drives the t -> t+1 transition
    u <- drought[, -T]; drought <- (drought - mean(u)) / stats::sd(as.vector(u))
    psi <- stats::plogis(stats::qlogis(params$initial_occupancy_prob) + params$psi_beta_forest * forest)
    col <- stats::plogis(stats::qlogis(params$colonisation_prob) + params$col_beta_forest * forest)
    z <- matrix(0, M, T); z[, 1] <- stats::rbinom(M, 1, psi)
    for (t in 2:T) {
      ext <- stats::plogis(stats::qlogis(params$extinction_prob) + params$ext_beta_drought * drought[, t - 1])
      z[, t] <- stats::rbinom(M, 1, ifelse(z[, t - 1] == 1, 1 - ext, col))
    }
  }
  zz <- aperm(array(z, c(M, T, J)), c(1, 3, 2))
  pz <- zz * p; pz[w == 0] <- 0
  y <- array(stats::rbinom(M * J * T, 1, pz), c(M, J, T)); y[w == 0] <- NA
  ids <- sprintf("S%03d", seq_len(M))
  if (model == "static") {
    data <- list(model = "static", y = matrix(y, M, J), site_covs = data.frame(forest = forest),
                 obs_covs = list(effort = matrix(effort, M, J)),
                 formulas = list(psi = "~ forest", p = "~ effort"), site_ids = ids)
    truth <- c(occupancy_prob = params$occupancy_prob, detection_prob = params$detection_prob,
               psi_beta_forest = params$psi_beta_forest, p_beta_effort = params$p_beta_effort,
               prop_sites_occupied = mean(z))
  } else {
    data <- list(model = "dynamic", y = y, site_covs = data.frame(forest = forest),
                 obs_covs = list(effort = effort), yearly_covs = list(drought = drought),
                 formulas = list(psi = "~ forest", col = "~ forest", ext = "~ drought", p = "~ effort"),
                 site_ids = ids, season_labels = as.character(seq_len(T)))
    truth <- c(initial_occupancy_prob = params$initial_occupancy_prob,
               colonisation_prob = params$colonisation_prob, extinction_prob = params$extinction_prob,
               detection_prob = params$detection_prob, psi_beta_forest = params$psi_beta_forest,
               col_beta_forest = params$col_beta_forest, ext_beta_drought = params$ext_beta_drought,
               p_beta_effort = params$p_beta_effort,
               stats::setNames(colMeans(z), paste0("prop_sites_occupied_t", seq_len(T))),
               prop_occupied_change = mean(z[, T]) - mean(z[, 1]))
  }
  list(data = data, truth = truth)
}

# ----------------------------------------------------------------- estimates
occ_prop_names <- function(des) if (des$model == "static") "prop_sites_occupied" else
  c(paste0("prop_sites_occupied_t", seq_len(des$T)), "prop_occupied_change")

occ_prop_cols <- function(des, prop) {
  if (des$model == "static") return(stats::setNames(as.data.frame(prop[, 1, drop = FALSE]), "prop_sites_occupied"))
  out <- as.data.frame(prop); names(out) <- paste0("prop_sites_occupied_t", seq_len(des$T))
  out$prop_occupied_change <- prop[, des$T] - prop[, 1]
  out
}

occ_estimates <- function(des, eng) {
  inames <- occ_intercept_names(des$model)
  if (eng$engine == "unmarked") {
    rows <- list()
    for (s in des$subs) {
      nm <- paste0("beta_", s, "[", seq_along(des$colnames[[s]]), "]")
      b <- eng$point[nm]; se <- eng$se[nm]
      r0 <- wald_rows(b[1], se[1], inames[[s]], transform = stats::plogis)
      r0$sd <- stats::plogis(b[1]) * (1 - stats::plogis(b[1])) * se[1]
      rows[[length(rows) + 1]] <- r0
      if (length(b) > 1) rows[[length(rows) + 1]] <- wald_rows(b[-1], se[-1], occ_coef_public(s, des$colnames[[s]][-1]))
    }
    pr <- occ_prop_cols(des, eng$prop)
    rows[[length(rows) + 1]] <- derived_rows(as.matrix(pr))
    est <- do.call(rbind, rows)
  } else if (eng$engine == "nimble") {
    rep_mat <- function(m) {
      out <- list()
      for (s in des$subs) {
        nm <- paste0("beta_", s, "[", seq_along(des$colnames[[s]]), "]")
        out[[inames[[s]]]] <- stats::plogis(m[, nm[1]])
        if (length(nm) > 1) for (k in 2:length(nm)) out[[occ_coef_public(s, des$colnames[[s]][k])]] <- m[, nm[k]]
      }
      no <- m[, grep("^n_occ", colnames(m)), drop = FALSE] / des$M
      cbind(as.data.frame(out, check.names = FALSE), occ_prop_cols(des, no)) |> as.matrix()
    }
    ch <- lapply(eng$chains, rep_mat)
    s <- summarise_draws_eco(ch)
    est <- s$estimates
    eng$diag$rhat_max <- s$diagnostics$rhat_max; eng$diag$ess_min <- s$diagnostics$ess_min
    eng$diag$converged <- s$diagnostics$converged
    # finite-sample change can have rhat NaN when constant; ignore NA
  } else {  # naive
    rows <- list()
    for (s in names(eng$glm)) {
      g <- eng$glm[[s]]
      r0 <- wald_rows(g$coef[1], g$se[1], inames[[s]], transform = stats::plogis)
      r0$sd <- stats::plogis(g$coef[1]) * (1 - stats::plogis(g$coef[1])) * g$se[1]
      rows[[length(rows) + 1]] <- r0
      if (length(g$coef) > 1) rows[[length(rows) + 1]] <- wald_rows(g$coef[-1], g$se[-1], occ_coef_public(s, des$colnames[[s]][-1]))
    }
    pm <- colMeans(eng$ever); se <- sqrt(pm * (1 - pm) / des$M)
    pn <- if (des$model == "static") "prop_sites_occupied" else paste0("prop_sites_occupied_t", seq_len(des$T))
    rows[[length(rows) + 1]] <- data.frame(param = pn, estimate = pm, sd = se, lower = pmax(0, pm - 1.96 * se),
                                           upper = pmin(1, pm + 1.96 * se))
    est <- do.call(rbind, rows)
  }
  rownames(est) <- NULL
  list(estimates = est, diag = eng$diag)
}

# ----------------------------------------------------------------------- fit
occ_fit <- function(data, settings) {
  settings <- merge_list(occ_settings, settings)
  if (!is.null(settings$model) && settings$model != data$model)
    stop(sprintf("settings$model = '%s' but the data are '%s'; set data$model in the config instead", settings$model, data$model))
  set.seed(settings$seed)
  des <- occ_prepare(data, standardize = isTRUE(settings$standardize))
  engine <- match.arg(settings$engine, c("unmarked", "nimble", "naive"))
  eng <- switch(engine,
                unmarked = occ_engine_unmarked(des, settings),
                nimble = occ_engine_nimble(des, settings),
                naive = occ_engine_naive(des, settings))
  if (engine == "unmarked") {
    nprop <- min(settings$n_prop_draws, nrow(eng$draws))
    eng$prop <- occ_prop_draws(des, eng$draws[seq_len(nprop), , drop = FALSE])
    if (nrow(eng$draws) == 1) eng$prop <- eng$prop[rep(1, 2), , drop = FALSE]
  }
  e <- occ_estimates(des, eng)
  dg <- e$diag
  gof <- NULL
  if (engine == "unmarked") gof <- occ_gof_unmarked(des, eng, settings)
  if (engine == "nimble") {
    gof <- if (isTRUE(settings$gof_nboot > 0)) occ_gof_nimble(des, eng, settings) else NULL
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      dg$notes <- c(dg$notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                      dg$rhat_max, dg$ess_min))
  }
  if (!is.null(gof)) dg$notes <- c(dg$notes, if (engine == "nimble")
    sprintf("Goodness of fit (%s, %d draws): mean observed chi-square %.1f, observed/replicated ratio %.2f (not a c-hat), Bayesian p = %.3f.",
            gof$method, gof$nboot, gof$stat, gof$c_hat, gof$p_value)
    else sprintf("Goodness of fit (%s, %d replicates): chi-square %.1f, c-hat %.2f, p = %.3f.",
                 gof$method, gof$nboot, gof$stat, gof$c_hat, gof$p_value))
  dg$gof <- gof
  list(estimates = e$estimates, diagnostics = dg,
       object = list(engine = engine, des = des, draws = eng$draws, point = eng$point, vcov = eng$vcov,
                     prop = eng$prop, gof = gof, fm = eng$fm, glm = eng$glm, settings = settings,
                     loglik = eng$diag$loglik))
}

occ_settings <- list(engine = "unmarked", standardize = TRUE, seed = 1,
                     n_draws = 2000, n_prop_draws = 1000,
                     gof_nboot = 200, gof_ndraws = 300,
                     niter = 12000, nburnin = 2000, nchains = 3, thin = 1, prior_sd = 2, slice = TRUE)

occ_fast <- merge_list(occ_settings, list(gof_nboot = 0, n_draws = 1000, n_prop_draws = 400,
                                          niter = 4000, nburnin = 1000))

module <- ecomodule(
  id = "occu-single",
  title = "Single-species occupancy (static and dynamic)",
  description = paste(
    "Where does a species occur, and what habitat and survey conditions are associated with occupancy and with detecting it,",
    "allowing for the fact that a species can be present but missed? Static (single-season, MacKenzie et al. 2002):",
    "occupancy and per-survey detection with site and survey covariates. Dynamic (multi-season, MacKenzie et al. 2003):",
    "initial occupancy, colonisation and local extinction between seasons, and the proportion of surveyed sites occupied each season."),
  defaults = list(
    params = list(occupancy_prob = 0.5, detection_prob = 0.4, psi_beta_forest = 1, p_beta_effort = 0.5,
                  initial_occupancy_prob = 0.5, colonisation_prob = 0.2, extinction_prob = 0.3,
                  col_beta_forest = 0.75, ext_beta_drought = 0.75),
    design = list(model = "static", n_sites = 200, n_surveys = 4, n_seasons = 5, prop_missing = 0.1, p_site_sd = 0)),
  simulate = occ_simulate, fit = occ_fit, settings = occ_settings, fast_settings = occ_fast,
  read_data = occ_read_data, check_data = occ_check_data,
  interpret = occ_interpret, plots = occ_plots, figures = occ_figures, methods_text = occ_methods_text,
  power_params = c("psi_beta_forest", "p_beta_effort", "col_beta_forest", "ext_beta_drought", "prop_occupied_change"),
  engine = "unmarked (occu, colext) or nimble",
  references = c(
    "MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA, Langtimm CA (2002) Estimating site occupancy rates when detection probabilities are less than one. Ecology 83:2248-2255. doi:10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2",
    "MacKenzie DI, Nichols JD, Hines JE, Knutson MG, Franklin AB (2003) Estimating site occupancy, colonization, and local extinction when a species is detected imperfectly. Ecology 84:2200-2207. doi:10.1890/02-3090",
    "MacKenzie DI, Bailey LL (2004) Assessing the fit of site-occupancy models. Journal of Agricultural, Biological, and Environmental Statistics 9:300-318. doi:10.1198/108571104X3361",
    "Fiske I, Chandler R (2011) unmarked: An R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23. doi:10.18637/jss.v043.i10",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487",
    "Royle JA, Kery M (2007) A Bayesian state-space formulation of dynamic occupancy models. Ecology 88:1813-1823. doi:10.1890/06-0669.1",
    "Guillera-Arroita G, Ridout MS, Morgan BJT (2010) Design of occupancy studies with imperfect detection. Methods in Ecology and Evolution 1:131-139. doi:10.1111/j.2041-210X.2010.00017.x",
    "MacKenzie DI, Nichols JD, Royle JA, Pollock KH, Bailey LL, Hines JE (2018) Occupancy Estimation and Modeling, 2nd ed. Academic Press."))
