# Log-linear trend regression with AR(1) residuals (glmmTMB; Brooks et al. 2017).
#
#   count_{s,t} ~ NegBin(lambda_{s,t}, theta)       (family "nbinom2"; counts)
#   log lambda_{s,t} = alpha + alpha_s + mu * (t - tbar) + sum_k beta_k z_{k,t} + eps_{s,t}
#   alpha_s ~ N(0, sd_site^2) (several sites); eps_{s,.} ~ AR(1)(rho, sd_ar1) per site
#
# family "lognormal" fits log(index) with a Gaussian family instead (positive
# indices such as CPUE). mu is the average log-scale change per year, so the
# trend is 100 (exp(mu) - 1) % per year. Covariates here act on the *level*
# of abundance in a year, not on growth (contrast with model "bayes").
# This is the classic "route regression" style trend (e.g. BBS) with
# autocorrelation; it does not separate process and observation noise.

ll_simulate <- function(params, design, seed) {
  set.seed(seed)
  n <- design$n_years %||% 20; S <- design$n_sites %||% 5
  covs <- design$covariate_names %||% "spring_rain_mm"
  nc <- design$n_covariates %||% 1; covs <- covs[seq_len(nc)]
  beta <- rep_len(unlist(params$count_beta %||% 0), nc)
  raw <- matrix(stats::rnorm(n * nc, 120, 35), n, nc, dimnames = list(NULL, covs))
  zs <- if (nc) scale(raw) else matrix(0, n, 0)
  tt <- seq_len(n) - (n + 1) / 2
  sites <- sprintf("route%02d", seq_len(S))
  a_s <- stats::rnorm(S, 0, params$sd_site %||% 0.5)
  rho <- params$ar1_rho %||% 0.5; sar <- params$sd_ar1 %||% 0.2
  df <- do.call(rbind, lapply(seq_len(S), function(i) {
    e <- numeric(n); e[1] <- stats::rnorm(1, 0, sar)
    for (t in 2:n) e[t] <- rho * e[t - 1] + stats::rnorm(1, 0, sar * sqrt(1 - rho^2))
    eta <- log(params$mean_count %||% 40) + a_s[i] + params$growth_rate_mu * tt + (if (nc) as.vector(zs %*% beta) else 0) + e
    y <- stats::rnbinom(n, mu = exp(eta), size = params$nb_size %||% 10)
    d <- data.frame(year = ts_years0(design) + seq_len(n) - 1, site = sites[i], count = y)
    d$count[ts_drop_years(n, design$frac_missing)] <- NA
    for (k in covs) d[[k]] <- raw[, k]
    d
  }))
  truth <- c(growth_rate_mu = params$growth_rate_mu, trend_pct_per_year = 100 * (exp(params$growth_rate_mu) - 1))
  if (nc) truth <- c(truth, stats::setNames(beta, paste0("count_beta_", covs)),
                     stats::setNames(100 * (exp(beta) - 1), paste0("count_effect_pct_", covs)))
  list(data = list(model = "loglinear", df = df, covariates = covs, response = "count",
                   meta = list(source = "simulated", count_label = "Count")),
       truth = truth)
}

ll_fit <- function(data, settings) {
  eco_require("glmmTMB")
  suppressPackageStartupMessages(requireNamespace("glmmTMB"))
  df <- data$df[!is.na(data$df$count), ]
  covs <- data$covariates %||% character()
  yrs <- seq(min(data$df$year), max(data$df$year))
  cv <- ts_cov_matrix(data$df, covs, yrs)
  tbar <- mean(yrs)
  d <- data.frame(count = df$count, site = factor(df$site), yr = df$year - tbar,
                  yearf = factor(df$year, levels = yrs))
  for (k in covs) d[[k]] <- (df[[k]] - cv$scaling[[k]]$center) / cv$scaling[[k]]$scale
  S <- nlevels(d$site)
  fam <- settings$family %||% (if (identical(data$response, "index")) "lognormal" else "nbinom2")
  ar <- !isFALSE(settings$ar1)
  rhs <- c("yr", covs, if (S > 1) "(1 | site)", if (ar) "ar1(yearf + 0 | site)")
  fml <- stats::as.formula(paste(if (fam == "lognormal") "log(count)" else "count", "~", paste(rhs, collapse = " + ")))
  famobj <- switch(fam, nbinom2 = glmmTMB::nbinom2(), poisson = stats::poisson(), lognormal = stats::gaussian(),
                   stop("settings: family must be nbinom2, poisson or lognormal"))
  m <- suppressWarnings(glmmTMB::glmmTMB(fml, data = d, family = famobj))
  notes <- character()
  pdh <- isTRUE(m$sdr$pdHess)
  if (!pdh) notes <- c(notes, "Hessian not positive definite: a variance component (usually the AR(1) SD, which competes with negative-binomial overdispersion) is estimated at ~0. Fixed-effect intervals are computed with it at zero.")
  co <- summary(m)$coefficients$cond
  z <- stats::qnorm(0.975)
  mu <- co["yr", 1]; se <- co["yr", 2]
  rows <- list(
    ts_rows("growth_rate_mu", mu, se, mu - z * se, mu + z * se),
    ts_rows("trend_pct_per_year", 100 * (exp(mu) - 1), 100 * exp(mu) * se, 100 * (exp(mu - z * se) - 1), 100 * (exp(mu + z * se) - 1)))
  for (k in covs) {
    b <- co[k, 1]; s <- co[k, 2]
    rows[[length(rows) + 1]] <- ts_rows(paste0("count_beta_", k), b, s, b - z * s, b + z * s)
    rows[[length(rows) + 1]] <- ts_rows(paste0("count_effect_pct_", k), 100 * (exp(b) - 1), 100 * exp(b) * s,
                                        100 * (exp(b - z * s) - 1), 100 * (exp(b + z * s) - 1))
  }
  if (ar) {
    vc <- glmmTMB::VarCorr(m)$cond
    arv <- vc[[length(vc)]]
    rows[[length(rows) + 1]] <- ts_rows("ar1_rho", attr(arv, "correlation")[1, 2], NA, NA, NA)
    rows[[length(rows) + 1]] <- ts_rows("sd_ar1", attr(arv, "stddev")[1], NA, NA, NA)
  }
  est <- do.call(rbind, rows)
  # fitted trajectories per site (conditional on random effects) on the count scale
  nd <- expand.grid(yearf = factor(yrs, levels = yrs), site = levels(d$site))
  nd$yr <- as.numeric(as.character(nd$yearf)) - tbar
  for (k in covs) nd[[k]] <- cv$Z[as.character(nd$yearf), k]
  if (anyNA(nd[covs])) for (k in covs) nd[[k]][is.na(nd[[k]])] <- 0
  pr <- tryCatch(stats::predict(m, newdata = nd, se.fit = TRUE, type = "link", allow.new.levels = TRUE), error = function(e) NULL)
  traj <- if (!is.null(pr)) data.frame(year = as.numeric(as.character(nd$yearf)), site = as.character(nd$site),
                                       estimate = exp(pr$fit), lower = exp(pr$fit - z * pr$se.fit), upper = exp(pr$fit + z * pr$se.fit),
                                       lower50 = exp(pr$fit - 0.674 * pr$se.fit), upper50 = exp(pr$fit + 0.674 * pr$se.fit)) else NULL
  fe_ok <- all(is.finite(co[c("yr", covs), 2]))
  conv <- fe_ok && m$fit$convergence == 0
  tm <- ts_matrix(data$df)
  list(estimates = est,
       diagnostics = list(converged = conv, notes = notes, aic = stats::AIC(m)),
       object = list(model = "loglinear", engine = "glmmTMB", glmm = m, trajectory = traj, years = yrs,
                     sites = tm$sites, Y = tm$Y, covariates = covs, scaling = cv$scaling, cov_raw = cv$raw,
                     tbar = tbar, family = fam, d = d))
}
