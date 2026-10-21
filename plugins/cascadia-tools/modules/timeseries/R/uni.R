# Univariate exponential-growth (EGSS) and Gompertz (GSS) state-space models:
# simulation and fitting via the Kalman-filter ML in kalman.R.

ts_years0 <- function(design) design$first_year %||% 1995

# Random missing years (never the first or last year).
ts_drop_years <- function(n, frac) {
  if (is.null(frac) || frac <= 0 || n < 4) return(integer())
  sample(2:(n - 1), min(n - 2, round(frac * n)))
}

uni_simulate <- function(params, design, seed, form) {
  set.seed(seed)
  n <- design$n_years %||% 25
  sp <- params$sd_process; so <- params$sd_observation
  x <- numeric(n)
  if (form == "exponential") {
    x[1] <- log(design$start_N %||% 500)
    for (t in 2:n) x[t] <- x[t - 1] + params$growth_rate_mu + stats::rnorm(1, 0, sp)
  } else {
    b <- params$density_dependence_b; lK <- log(params$equilibrium_N)
    sd0 <- if (b < 0 && b > -2) sp / sqrt(1 - (1 + b)^2) else 0
    x[1] <- if (!is.null(design$start_N)) log(design$start_N) else lK + stats::rnorm(1, 0, sd0)
    for (t in 2:n) x[t] <- x[t - 1] + b * (x[t - 1] - lK) + stats::rnorm(1, 0, sp)
  }
  y <- exp(x + stats::rnorm(n, 0, so))
  if (isTRUE(design$round_counts %||% TRUE)) y <- pmax(round(y), 1)
  y[ts_drop_years(n, design$frac_missing)] <- NA
  df <- data.frame(year = ts_years0(design) + seq_len(n) - 1, site = "all", count = y)
  truth <- if (form == "exponential")
    c(growth_rate_mu = params$growth_rate_mu, trend_pct_per_year = 100 * (exp(params$growth_rate_mu) - 1))
  else c(density_dependence_b = params$density_dependence_b,
         if (params$density_dependence_b < 0) c(equilibrium_N = params$equilibrium_N))
  truth <- c(truth, sd_process = sp, sd_observation = so, N_final = exp(x[n]))
  list(data = list(model = form, df = df, covariates = character(), response = "count",
                   meta = list(source = "simulated", count_label = "Count", x_true = x)),
       truth = truth)
}

uni_fit <- function(data, settings, form) {
  tm <- ts_matrix(data$df)
  y <- ts_logY(tm$Y[1, ], settings)
  obs_error <- !isFALSE(settings$observation_error)
  est_method <- tolower(settings$estimation %||% "reml")
  if (form == "exponential" && obs_error && est_method == "reml") {
    sm <- egss_reml_fit(y, n_draws = settings$n_draws %||% 1000, seed = settings$seed %||% 1)
    kf <- sm$kf
  } else {
    kf <- kf_fit_uni(y, form, obs_error = obs_error)
    kf$method <- "ML"
    sm <- kf_summaries(kf, n_draws = settings$n_draws %||% 1000, seed = settings$seed %||% 1)
    kf$sir_ess <- sm$sir_ess
  }
  if (form == "gompertz" && (settings$dd_test_boot %||% 199) > 0) {
    tst <- gss_dd_test(y, obs_error, B = settings$dd_test_boot %||% 199, seed = settings$seed %||% 1)
    kf$dd_test <- tst
    kf$notes <- c(kf$notes, sprintf("Test for density dependence (parametric-bootstrap likelihood ratio, Gompertz vs exponential growth, %d simulations): LR = %.2f, p = %.3f. Use this p-value, not whether the interval for b excludes 0.",
                                    tst$B, tst$lr, tst$p))
  }
  ok <- which(!is.na(y)); yrs <- tm$years[min(ok):max(ok)]
  # goodness of fit: standardised innovations at the point estimates
  yy <- y[min(ok):max(ok)]
  innov <- if (!obs_error) {
    # process-only model: innovations are standardised changes between counted years
    th <- kf$theta; sp <- exp(th[["log_sd_process"]])
    cc <- if (form == "exponential") 1 else 1 + th[["b"]]; dd <- if (form == "exponential") th[["mu"]] else -th[["b"]] * th[["logK"]]
    kf_innov(yy, cc, dd, sp, 1e-8, x0 = (yy[1] - dd) / cc)
  } else if (identical(kf$method, "REML")) {
    est <- sm$estimates
    kf_innov(yy, 1, est$estimate[est$param == "growth_rate_mu"], max(est$estimate[est$param == "sd_process"], 1e-4),
             max(est$estimate[est$param == "sd_observation"], 1e-4), diffuse = TRUE)
  } else { p0 <- kf_unpack(kf$theta, form, TRUE); kf_innov(yy, p0$cc, p0$dd, max(p0$sp, 1e-4), p0$so, p0$x0) }
  if (!obs_error) innov[1] <- NA
  gof <- ts_resid_check(innov); gof$innov <- data.frame(year = yrs, z = innov)
  traj <- cbind(data.frame(year = yrs, site = tm$sites[1]), sm$trajectory[, -1])
  # parameter draws for forecasts
  pd <- ts_param_draws_ml(kf, sm$theta_draws, sm$state_draws)
  notes <- kf$notes
  if (!is.null(sm$sir_ess) && sm$sir_ess < 200)
    notes <- c(notes, sprintf("Abundance intervals: the importance sampler that propagates parameter uncertainty had a low effective sample size (%.0f); intervals may be rough.", sm$sir_ess))
  if (!obs_error) notes <- c(notes, "observation_error = FALSE: all year-to-year noise is treated as real population change (no sampling error). Use only if counts are near-complete censuses.")
  list(estimates = sm$estimates,
       diagnostics = list(converged = kf$converged, notes = notes, loglik = kf$loglik, aic = kf$aic,
                          dd_test_p = kf$dd_test$p, dd_test_lr = kf$dd_test$lr,
                          gof_ljung_box_p = gof$ljung_box_p, gof_shapiro_p = gof$shapiro_p),
       object = list(model = form, engine = if (kf$method == "REML") "REML (differenced log counts) + diffuse Kalman smoother" else "Kalman-filter maximum likelihood", kf = kf, method = kf$method, trajectory = traj,
                     years = tm$years, sites = tm$sites, Y = tm$Y, pdraws = pd,
                     obs_error = obs_error, gof = gof))
}

# Parameter draws (asymptotic normal, from kf_summaries) + matching
# final-state draws, used by forecasts and density-dependence curves.
ts_param_draws_ml <- function(kf, D, X) {
  out <- data.frame(sd_process = exp(D[, "log_sd_process"]))
  if (kf$form == "exponential") { out$r0 <- D[, "mu"]; out$b <- 0; out$xc <- 0 }
  else { out$b <- D[, "b"]; out$xc <- D[, "logK"]; out$r0 <- 0 }
  out$xT <- X[, ncol(X)]
  out
}
