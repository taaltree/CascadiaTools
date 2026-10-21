# Multivariate state-space model for several sub-populations (MARSS package;
# Holmes, Ward & Wills 2012; Holmes, Scheuerell & Ward, MARSS User Guide).
#
#   x_t = x_{t-1} + u + w_t,  w_t ~ MVN(0, Q)        (one state per site, B = identity)
#   y_t = Z x_t + a + v_t,    v_t ~ MVN(0, R)        (y = log counts, one row per site)
#
# settings$trend:   "shared" (U = "equal": one growth rate for all sites) or
#                   "separate" (U = "unequal": a growth rate per site)
# settings$states:  "separate" (Z = identity: each site its own population) or
#                   "single" (Z = column of ones: all sites survey one population;
#                   site offsets a estimated, A = "scaling")
# settings$process_cov: Q structure, "diagonal and equal" (default),
#                   "diagonal and unequal", "equalvarcov" (correlated sites)
# R is "diagonal and equal" (same observation SD at every site).
# Fitted by ML with method = "BFGS" ("kem" for equalvarcov); Hessian-based intervals (MARSSparamCIs);
# variances are converted to SDs with log-scale delta-method intervals.

marss_simulate <- function(params, design, seed) {
  set.seed(seed)
  n <- design$n_years %||% 25; S <- design$n_sites %||% 4
  sp <- params$sd_process; so <- params$sd_observation
  rho <- params$process_corr %||% 0
  u <- if (length(params$growth_rate_mu_sites)) rep_len(unlist(params$growth_rate_mu_sites), S) else rep(params$growth_rate_mu, S)
  Sig <- sp^2 * ((1 - rho) * diag(S) + rho)
  L <- t(chol(Sig))
  x <- matrix(NA_real_, S, n)
  x[, 1] <- log(exp(stats::runif(S, log(design$start_N_min %||% 80), log(design$start_N_max %||% 800))))
  for (t in 2:n) x[, t] <- x[, t - 1] + u + as.vector(L %*% stats::rnorm(S))
  y <- exp(x + matrix(stats::rnorm(S * n, 0, so), S, n))
  if (isTRUE(design$round_counts %||% TRUE)) y <- pmax(round(y), 1)
  sites <- sprintf("site%02d", seq_len(S))
  df <- do.call(rbind, lapply(seq_len(S), function(i) {
    yi <- y[i, ]; yi[ts_drop_years(n, design$frac_missing)] <- NA
    data.frame(year = ts_years0(design) + seq_len(n) - 1, site = sites[i], count = yi)
  }))
  truth <- if (length(params$growth_rate_mu_sites)) {
    stats::setNames(c(u, 100 * (exp(u) - 1)), c(paste0("growth_rate_mu_", sites), paste0("trend_pct_per_year_", sites)))
  } else c(growth_rate_mu = u[1], trend_pct_per_year = 100 * (exp(u[1]) - 1))
  truth <- c(truth, sd_process = sp, sd_observation = so)
  if (rho != 0 || identical(design$process_cov, "equalvarcov")) truth <- c(truth, process_corr = rho)
  list(data = list(model = "marss", df = df, covariates = character(), response = "count",
                   meta = list(source = "simulated", count_label = "Count", x_true = x)),
       truth = truth)
}

marss_fit <- function(data, settings) {
  eco_require("MARSS")
  suppressPackageStartupMessages(requireNamespace("MARSS"))
  tm <- ts_matrix(data$df)
  Y <- ts_logY(tm$Y, settings)
  S <- nrow(Y)
  trend <- settings$trend %||% "shared"; states <- settings$states %||% "separate"
  qcov <- settings$process_cov %||% "diagonal and equal"
  if (states == "single") qcov <- "diagonal and equal"
  m <- if (states == "single") 1 else S
  est_method <- tolower(settings$estimation %||% "reml")
  if (est_method == "reml" && states == "separate" && qcov == "diagonal and equal" &&
      identical(settings$obs_cov %||% "diagonal and equal", "diagonal and equal"))
    return(marss_reml_wrap(Y, tm, trend, settings))
  mod <- list(B = "identity", U = if (trend == "shared" || m == 1) "equal" else "unequal",
              Q = if (m == 1) "diagonal and equal" else qcov,
              Z = if (states == "single") matrix(1, S, 1) else "identity",
              A = if (states == "single") "scaling" else "zero",
              R = settings$obs_cov %||% "diagonal and equal",
              x0 = "unequal", tinitx = 0)  # x0 one year before the first count: avoids the unbounded-likelihood spike at R -> 0 (see kalman.R)
  # BFGS cannot estimate shared variances + a shared covariance ("equalvarcov");
  # MARSS's EM algorithm (kem) can.
  meth <- if (identical(mod$Q, "equalvarcov")) "kem" else "BFGS"
  fit <- tryCatch(MARSS::MARSS(Y, model = mod, method = meth, silent = TRUE,
                               control = if (meth == "kem") list(maxit = 5000, conv.test.slope.tol = 0.05) else list(maxit = 5000)),
                  error = function(e) NULL)
  if (is.null(fit)) stop("MARSS fit failed")
  ci <- tryCatch(suppressWarnings(MARSS::MARSSparamCIs(fit)), error = function(e) NULL)
  notes <- character()
  conv <- fit$convergence == 0 && !is.null(ci)
  z <- stats::qnorm(0.975)
  rows <- list()
  U <- fit$par$U; Use <- if (!is.null(ci)) ci$par.se$U else U * NA
  un <- rownames(U)
  site_lab <- if (states == "single") "all" else tm$sites
  if (nrow(U) == 1) {
    rows$mu <- ts_rows("growth_rate_mu", U[1], Use[1], U[1] - z * Use[1], U[1] + z * Use[1])
    rows$tr <- ts_rows("trend_pct_per_year", 100 * (exp(U[1]) - 1), 100 * exp(U[1]) * Use[1],
                       100 * (exp(U[1] - z * Use[1]) - 1), 100 * (exp(U[1] + z * Use[1]) - 1))
  } else {
    rows$mu <- ts_rows(paste0("growth_rate_mu_", site_lab), U[, 1], Use[, 1], U[, 1] - z * Use[, 1], U[, 1] + z * Use[, 1])
    rows$tr <- ts_rows(paste0("trend_pct_per_year_", site_lab), 100 * (exp(U[, 1]) - 1), 100 * exp(U[, 1]) * Use[, 1],
                       100 * (exp(U[, 1] - z * Use[, 1]) - 1), 100 * (exp(U[, 1] + z * Use[, 1]) - 1))
  }
  var_rows <- function(v, se, name, nm2 = NULL) {
    sdv <- sqrt(v); sel <- se / (2 * v)   # SE of log(sd)
    ts_rows(if (is.null(nm2)) name else paste0(name, "_", nm2), sdv, sdv * sel, sdv * exp(-z * sel), sdv * exp(z * sel))
  }
  Qp <- fit$par$Q; Qse <- if (!is.null(ci)) ci$par.se$Q else Qp * NA
  qn <- rownames(Qp)
  if (qcov == "equalvarcov" && m > 1) {
    iv <- which(qn == "diag"); ic <- which(qn == "offdiag")
    rows$sp <- var_rows(Qp[iv, 1], Qse[iv, 1], "sd_process")
    r <- Qp[ic, 1] / Qp[iv, 1]
    # delta method for the correlation; Fisher-z interval
    g <- c(-Qp[ic, 1] / Qp[iv, 1]^2, 1 / Qp[iv, 1])
    rse <- sqrt(sum((g * c(Qse[iv, 1], Qse[ic, 1]))^2))
    zz <- atanh(max(min(r, 0.999), -0.999)); zse <- rse / (1 - r^2)
    rows$rho <- ts_rows("process_corr", r, rse, tanh(zz - z * zse), tanh(zz + z * zse))
    notes <- c(notes, "Interval for process_corr ignores the covariance between the variance and covariance estimates (approximate).")
  } else if (nrow(Qp) == 1) {
    rows$sp <- var_rows(Qp[1, 1], Qse[1, 1], "sd_process")
  } else rows$sp <- var_rows(Qp[, 1], Qse[, 1], "sd_process", site_lab)
  Rp <- fit$par$R; Rse <- if (!is.null(ci)) ci$par.se$R else Rp * NA
  rows$so <- if (nrow(Rp) == 1) var_rows(Rp[1, 1], Rse[1, 1], "sd_observation") else var_rows(Rp[, 1], Rse[, 1], "sd_observation", tm$sites)
  est <- do.call(rbind, rows); rownames(est) <- NULL
  if (any(!is.finite(est$sd))) { conv <- FALSE; notes <- c(notes, "Some standard errors are not finite: a variance is probably estimated at zero (boundary). Intervals unreliable.") }
  if (any(c(Qp[, 1], Rp[, 1]) < 1e-6)) { conv <- FALSE; notes <- c(notes, "A variance was estimated at (near) zero: the data cannot separate process from observation error at this design.") }
  # Smoothed states (conditional on the MLE)
  xs <- fit$states; xse <- fit$states.se
  traj <- do.call(rbind, lapply(seq_len(m), function(i) data.frame(
    year = tm$years, site = site_lab[i], estimate = exp(xs[i, ]),
    lower = exp(xs[i, ] - z * xse[i, ]), upper = exp(xs[i, ] + z * xse[i, ]),
    lower50 = exp(xs[i, ] - 0.674 * xse[i, ]), upper50 = exp(xs[i, ] + 0.674 * xse[i, ]))))
  # parameter draws for forecasts: mu and sd_process from their intervals
  set.seed(settings$seed %||% 1)
  nd <- settings$n_draws %||% 1000
  mu_hat <- est$estimate[grepl("^growth_rate_mu", est$param)]; mu_se <- est$sd[grepl("^growth_rate_mu", est$param)]
  sp_row <- est[grepl("^sd_process", est$param), ][1, ]
  pd <- lapply(seq_len(m), function(i) {
    k <- min(i, length(mu_hat))
    data.frame(site = site_lab[i], r0 = stats::rnorm(nd, mu_hat[k], ifelse(is.finite(mu_se[k]), mu_se[k], 0)), b = 0, xc = 0,
               sd_process = sp_row$estimate * exp(stats::rnorm(nd, 0, ifelse(is.finite(sp_row$sd), sp_row$sd / sp_row$estimate, 0))),
               xT = stats::rnorm(nd, xs[i, ncol(xs)], xse[i, ncol(xse)]))
  })
  list(estimates = est,
       diagnostics = list(converged = conv, notes = notes, loglik = fit$logLik, aic = fit$AIC),
       object = list(model = "marss", engine = paste0("MARSS (", meth, ") maximum likelihood"), method = "ML", marss = fit, trajectory = traj,
                     years = tm$years, sites = tm$sites, Y = tm$Y, pdraws = do.call(rbind, pd),
                     trend = trend, states = states, process_cov = qcov))
}

# model: marss with the default structure, estimated by pooled REML (reml.R).
marss_reml_wrap <- function(Y, tm, trend, settings) {
  shared <- trend == "shared"
  r <- marss_reml_fit(Y, shared = shared, n_draws = settings$n_draws %||% 1000, seed = settings$seed %||% 1)
  S <- nrow(Y); a <- 0.025
  traj <- do.call(rbind, lapply(seq_len(S), function(s) {
    Xs <- r$X[, s, , drop = TRUE]; if (!is.matrix(Xs)) Xs <- matrix(Xs, ncol = ncol(Y))
    q <- function(pp) exp(apply(Xs, 2, stats::quantile, pp, na.rm = TRUE))
    data.frame(year = tm$years, site = tm$sites[s], estimate = q(0.5), lower = q(a), upper = q(1 - a),
               lower50 = q(0.25), upper50 = q(0.75))
  }))
  pd <- do.call(rbind, lapply(seq_len(S), function(s)
    data.frame(site = tm$sites[s], r0 = r$mu_draws[, min(s, ncol(r$mu_draws))], b = 0, xc = 0,
               sd_process = r$sp_draws, xT = r$X[, s, ncol(Y)])))
  e <- r$estimates
  sp <- max(e$estimate[e$param == "sd_process"], 1e-4); so <- max(e$estimate[e$param == "sd_observation"], 1e-4)
  mu <- e$estimate[grepl("^growth_rate_mu", e$param)]
  inn <- do.call(rbind, lapply(seq_len(S), function(s)
    data.frame(year = tm$years, site = tm$sites[s], z = kf_innov(Y[s, ], 1, mu[min(s, length(mu))], sp, so, diffuse = TRUE))))
  gof <- ts_resid_check(inn$z); gof$innov <- inn
  list(estimates = r$estimates,
       diagnostics = list(converged = r$converged, notes = r$notes, loglik_reml = r$loglik_reml,
                          gof_ljung_box_p = gof$ljung_box_p, gof_shapiro_p = gof$shapiro_p),
       object = list(model = "marss", engine = "pooled REML (differenced log counts) + diffuse Kalman smoother",
                     method = "REML", trajectory = traj, years = tm$years, sites = tm$sites, Y = tm$Y, pdraws = pd,
                     trend = trend, states = "separate", process_cov = "diagonal and equal", gof = gof))
}
