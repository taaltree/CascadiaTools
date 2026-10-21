# Publication figure suite for the timeseries module (docs/FIGURES.md row
# "timeseries": trajectory, trend, density_dependence, growth_<cov>, forecast,
# power_years) plus variances, coefficients, trace and posterior.

ts_obs_df <- function(obj) {
  Y <- obj$Y
  d <- data.frame(time = rep(as.numeric(colnames(Y)), each = nrow(Y)), site = rep(rownames(Y), ncol(Y)),
                  y = as.vector(Y))
  d[!is.na(d$y), ]
}

# Forecast paths from parameter draws: x_{t+1} = x_t + r0 + b (x_t - xc) + e
ts_forecast_draws <- function(pd, h, seed = 11) {
  set.seed(seed)
  n <- nrow(pd)
  X <- matrix(NA_real_, n, h)
  x <- pd$xT
  for (j in seq_len(h)) {
    x <- x + pd$r0 + pd$b * (x - pd$xc) + stats::rnorm(n, 0, pd$sd_process)
    X[, j] <- x
  }
  X
}

# sparse year breaks so neighbouring facets' tick labels do not collide
ts_year_breaks <- function(l) { b <- pretty(l, n = 3); b[b >= l[1] & b <= l[2]] }

# widen a figure (keeps any height hint the constructor set)
ts_widen <- function(p, width = "onehalf", aspect = NULL) {
  a <- attr(p, "fig_size"); if (is.null(a)) a <- list()
  a$width <- width; if (!is.null(aspect)) a$aspect <- aspect
  attr(p, "fig_size") <- a; p
}

ts_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  obj <- fit$object; est <- fit$estimates
  ci <- opts$ci %||% 0.95; a <- (1 - ci) / 2
  ylab <- fig_label("count", data$meta$count_label %||% "Count", opts)
  ylab_N <- fig_label("abundance", if (obj$model == "loglinear") paste("Expected", tolower(ylab)) else "Estimated abundance", opts)
  xlab <- fig_label("year", "Year", opts)
  multi <- length(unique(obj$trajectory$site)) > 1
  out <- list()
  ci_word <- if (obj$model == "bayes") "credible" else "confidence"
  # ---- trajectory
  tr <- obj$trajectory; tr$time <- tr$year
  obs <- ts_obs_df(obj)
  if (identical(obj$states, "single")) obs$site <- "all"
  p <- fig_timeseries(tr[, c("time", "estimate", "lower", "upper", "site")], ylab_N,
                      obs = if (isTRUE(opts$show_data)) obs else NULL, xlab = xlab)
  if (all(c("lower50", "upper50") %in% names(tr)) && obj$model == "bayes")
    p$layers <- c(ggplot2::geom_ribbon(data = tr, ggplot2::aes(time, ymin = lower50, ymax = upper50), inherit.aes = FALSE,
                                       fill = eco_pal$band, alpha = 0.25), p$layers)
  if (multi || identical(obj$states, "single")) p <- p + ggplot2::facet_wrap(~site, scales = "free_y") + ggplot2::scale_x_continuous(breaks = ts_year_breaks)
  p <- p + ggplot2::scale_y_continuous(labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE))
  if (multi) attr(p, "fig_size") <- list(width = "double", aspect = 0.6)
  cap_tr <- if (obj$model == "loglinear")
    sprintf("Fitted expected count per year (line; log-linear trend plus AR(1) deviations) with %d%% confidence band. Points are the observed %s.",
            round(100 * ci), tolower(ylab))
  else sprintf("Estimated true abundance (line, smoothed latent state) with %d%% %s band%s. Points are the observed %s; the gap between points and line is what the model attributes to observation error.",
               round(100 * ci), ci_word, if (obj$model == "bayes") " (inner band 50%)" else "", tolower(ylab))
  out$trajectory <- with_caption(p, cap_tr)
  # ---- trend
  tr_rows <- est[grepl("trend_pct_per_year", est$param), ]
  if (nrow(tr_rows)) {
    pn <- tr_rows$param
    tr_rows$param <- ifelse(grepl("^realised", pn), "Realised trend over study",
                            ifelse(pn == "trend_pct_per_year", "Average trend", sub("^trend_pct_per_year_", "", pn)))
    out$trend <- with_caption(ts_widen(fig_coef(tr_rows, null_value = 0, xlab = fig_label("trend_axis", "Population change (% per year)", opts), opts = opts)),
                              sprintf("Estimated population trend in %% per year with %d%% %s interval; dashed line = no change.%s",
                                      round(100 * ci), ci_word,
                                      if (obj$model == "bayes") " Realised trend = average annual change of the estimated true population between the first and last year." else ""))
  }
  # ---- density dependence
  if (obj$model %in% c("gompertz", "bayes") && (obj$model == "gompertz" || isTRUE(obj$dd))) {
    pd <- obj$pdraws
    rng <- range(c(tr$lower, tr$upper))
    Ng <- exp(seq(log(rng[1] * 0.7), log(rng[2] * 1.3), length.out = 80))
    G <- sapply(Ng, function(N) pd$r0 + pd$b * (log(N) - pd$xc))
    pred <- data.frame(x = Ng, estimate = apply(G, 2, stats::median), lower = apply(G, 2, stats::quantile, a),
                       upper = apply(G, 2, stats::quantile, 1 - a))
    Y1 <- obj$Y[1, ]; yy <- as.numeric(Y1)
    ob <- data.frame(x = yy[-length(yy)], y = log(yy[-1] / yy[-length(yy)]))
    ob <- ob[is.finite(ob$x) & is.finite(ob$y), ]
    p <- fig_cond_effect(pred, fig_label("dd_x", paste(ylab_N, "in year t (log scale)"), opts),
                         fig_label("dd_y", "Growth rate, log(N[t+1]/N[t])", opts),
                         obs = if (isTRUE(opts$show_data)) ob else NULL) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = eco_pal$null) + ggplot2::scale_x_log10()
    out$density_dependence <- with_caption(p, sprintf(paste(
      "Expected annual growth rate vs abundance (line, %d%% %s band)%s. Points: growth computed naively from consecutive raw counts;",
      "observation error makes these points look more strongly density dependent than the true process (regression-to-the-mean), which is why the state-space fit is flatter.",
      "Where the line crosses 0 is the equilibrium abundance."), round(100 * ci), ci_word,
      if (length(obj$covariates)) ", covariates at their mean" else ""))
  }
  # ---- covariate effects on growth (bayes) / on counts (loglinear)
  if (obj$model == "bayes" && length(obj$covariates)) {
    pd <- obj$pdraws; B <- as.matrix(pd$cov_beta)
    for (k in seq_along(obj$covariates)) {
      cv <- obj$covariates[k]; sc <- obj$scaling[[cv]]
      xr <- range(obj$cov_raw[, cv], na.rm = TRUE)
      xs <- seq(xr[1], xr[2], length.out = 60)
      zz <- (xs - sc$center) / sc$scale
      L <- exp(outer(pd$r0, rep(1, length(xs))) + outer(B[, k], zz))
      pred <- data.frame(x = fig_x(xs, cv, opts), estimate = apply(L, 2, stats::median),
                         lower = apply(L, 2, stats::quantile, a), upper = apply(L, 2, stats::quantile, 1 - a),
                         lower50 = apply(L, 2, stats::quantile, 0.25), upper50 = apply(L, 2, stats::quantile, 0.75))
      lab <- fig_label(cv, data$meta$covariate_labels[[cv]] %||% gsub("_", " ", cv), opts)
      p <- fig_cond_effect(pred, lab, fig_label("lambda", "Growth multiplier, N[t]/N[t-1]", opts),
                           obs = data.frame(x = fig_x(obj$cov_raw[-1, cv], cv, opts))) +
        ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.07)))
      out[[paste0("growth_", cv)]] <- with_caption(p, sprintf(paste(
        "Expected annual growth multiplier (lambda; 1 = stable) vs %s, at average abundance%s;",
        "line = posterior median, bands = %d%% and 50%% credible intervals; rug = observed values. Values below 1 mean decline."),
        lab, if (length(obj$covariates) > 1) " and other covariates at their mean" else "", round(100 * ci)))
    }
  }
  if (obj$model == "loglinear" && length(obj$covariates)) {
    m <- obj$glmm; co <- glmmTMB::fixef(m)$cond; V <- stats::vcov(m)$cond
    for (cv in obj$covariates) {
      sc <- obj$scaling[[cv]]; xr <- range(obj$cov_raw[, cv], na.rm = TRUE)
      xs <- seq(xr[1], xr[2], length.out = 60); zz <- (xs - sc$center) / sc$scale
      X <- cbind(1, zz); idx <- c("(Intercept)", cv)
      eta <- as.vector(X %*% co[idx]); se <- sqrt(rowSums((X %*% V[idx, idx]) * X))
      z <- stats::qnorm(1 - a)
      pred <- data.frame(x = fig_x(xs, cv, opts), estimate = exp(eta), lower = exp(eta - z * se), upper = exp(eta + z * se))
      if (obj$family == "lognormal") pred[, -1] <- pred[, -1]  # geometric-mean index
      lab <- fig_label(cv, data$meta$covariate_labels[[cv]] %||% gsub("_", " ", cv), opts)
      p <- fig_cond_effect(pred, lab, paste(ylab, "(typical site, mid-study year)"),
                           obs = data.frame(x = fig_x(obj$cov_raw[, cv], cv, opts)))
      out[[paste0("count_", cv)]] <- with_caption(p, sprintf("Expected %s at a typical site in the middle year vs %s (line) with %d%% confidence band; rug = observed values.",
                                                            tolower(ylab), lab, round(100 * ci)))
    }
  }
  # ---- coefficients
  cr <- est[grepl("^(growth|count)_effect_pct_", est$param), ]
  if (nrow(cr)) {
    cr$param <- gsub("^(growth|count)_effect_pct_", "", cr$param)
    out$coefficients <- with_caption(ts_widen(fig_coef(cr, null_value = 0, opts = opts,
      xlab = fig_label("coef_axis", if (obj$model == "bayes") "Change in growth per 1 SD (%)" else "Change in count per 1 SD (%)", opts))),
      sprintf("Covariate effects per 1 standard deviation of each covariate with %d%% %s intervals; dashed line = no effect.", round(100 * ci), ci_word))
  }
  # ---- variances
  vr <- est[grepl("^sd_(process|observation)", est$param), ]
  if (nrow(vr)) {
    vr$param <- sub("^sd_process", "Process SD (real change)", sub("^sd_observation", "Observation SD (counting error)", vr$param))
    out$variances <- with_caption(ts_widen(fig_coef(vr, null_value = 0, order = FALSE, opts = opts, xlab = fig_label("sd_axis", "SD on the log scale", opts))),
      sprintf("Process noise (year-to-year variation in true growth) and observation error (counting/sampling noise), both as SDs of log abundance, with %d%% %s intervals. 0.1 on this scale is roughly 10%% variation.",
              round(100 * ci), ci_word))
  }
  # ---- forecast
  if (!is.null(obj$pdraws) && obj$model != "loglinear") {
    h <- opts$forecast_years %||% 10
    pd <- obj$pdraws
    groups <- if ("site" %in% names(pd)) unique(pd$site) else tr$site[1]
    fc <- do.call(rbind, lapply(groups, function(g) {
      pg <- if ("site" %in% names(pd)) pd[pd$site == g, ] else pd
      X <- ts_forecast_draws(pg, h)
      ly <- max(tr$year[tr$site == g])
      data.frame(time = ly + seq_len(h), site = g, estimate = exp(apply(X, 2, stats::median)),
                 lower = exp(apply(X, 2, stats::quantile, a)), upper = exp(apply(X, 2, stats::quantile, 1 - a)))
    }))
    hist <- tr[, c("time", "site", "estimate", "lower", "upper")]
    both <- rbind(cbind(hist, part = "Estimated"), cbind(fc, part = "Forecast"))
    fcd <- both[both$part == "Forecast", ]; hsd <- both[both$part == "Estimated", ]
    lab_d <- do.call(rbind, lapply(split(fcd, fcd$site), function(d) data.frame(site = d$site[1], time = min(d$time), y = max(d$upper))))
    p <- ggplot2::ggplot(both, ggplot2::aes(time, estimate)) +
      ggplot2::geom_ribbon(data = hsd, ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.25, colour = NA) +
      ggplot2::geom_ribbon(data = fcd, ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.12, colour = NA) +
      ggplot2::geom_line(data = hsd, colour = eco_pal$band, linewidth = 0.7) +
      ggplot2::geom_line(data = fcd, colour = eco_pal$band, linewidth = 0.7, linetype = 2) +
      ggplot2::geom_text(data = lab_d, ggplot2::aes(time, y, label = "Forecast"), inherit.aes = FALSE,
                         hjust = 0, vjust = -0.4, size = 2.6, colour = eco_pal$ink) +
      ggplot2::scale_y_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE),
                             expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
      ggplot2::labs(x = xlab, y = paste(ylab_N, "(log scale)"))
    if (isTRUE(opts$show_data)) p <- p + ggplot2::geom_point(data = obs[obs$y > 0, ], ggplot2::aes(time, y), inherit.aes = FALSE, size = 0.9, colour = eco_pal$ink)
    if (length(groups) > 1) { p <- p + ggplot2::facet_wrap(~site, scales = "free_y") + ggplot2::scale_x_continuous(breaks = ts_year_breaks); attr(p, "fig_size") <- list(width = "double", aspect = 0.6) }
    out$forecast <- with_caption(p, sprintf(paste(
      "Estimated abundance (solid) and %d-year forecast (dashed) with %d%% %s band; the band widens because process noise accumulates%s.",
      "Covariates (if any) are held at their average. Forecasts assume the dynamics of the study period continue."),
      h, round(100 * ci), if (obj$model == "bayes") "posterior predictive" else "approximate prediction",
      if (obj$model == "bayes") "" else " (parameter uncertainty propagated by sampling the variance split and growth parameters in proportion to their likelihood)"))
  }
  # ---- power: years of monitoring
  pf <- opts$power_csv
  if (!is.null(pf)) {
    pth <- if (file.exists(pf)) pf else file.path(eco_root(), pf)
    if (file.exists(pth)) {
      pw <- utils::read.csv(pth)
      pw <- pw[pw$param %in% c("growth_rate_mu", "trend_pct_per_year"), ]
      pw <- pw[pw$param == "growth_rate_mu", ]; pw$param <- "Detect the decline"
      pp <- fig_power(pw, xlab = fig_label("power_x", "Years of monitoring", opts))
      pp$mapping$colour <- NULL
      for (i in seq_along(pp$layers)) if (!inherits(pp$layers[[i]]$geom, "GeomHline")) pp$layers[[i]]$aes_params$colour <- eco_pal$band
      out$power_years <- with_caption(pp,
        sprintf("Power to detect the simulated trend (%.1f%% per year) as a function of monitoring length: share of %d simulated datasets per point whose 95%% interval for the trend excludes 0; bars = +/- 2 Monte Carlo SE; dashed line = 80%% power.",
                100 * (exp(pw$truth[1]) - 1), max(pw$n_converged, na.rm = TRUE)))
    }
  }
  # ---- model checks
  if (!is.null(obj$gof$innov)) {
    d <- obj$gof$innov; d <- d[is.finite(d$z), ]
    if (!"site" %in% names(d)) d$site <- "all"
    p <- ggplot2::ggplot(d, ggplot2::aes(year, z)) +
      ggplot2::geom_hline(yintercept = 0, colour = eco_pal$null) +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_line(colour = eco_pal$muted, linewidth = 0.4) +
      ggplot2::geom_point(colour = eco_pal$ink, size = 1.2) +
      ggplot2::labs(x = xlab, y = fig_label("innov", "Standardised prediction error", opts))
    if (length(unique(d$site)) > 1) { p <- p + ggplot2::facet_wrap(~site) + ggplot2::scale_x_continuous(breaks = ts_year_breaks); attr(p, "fig_size") <- list(width = "double", aspect = 0.6) }
    out$residuals <- with_caption(p, sprintf(paste(
      "Standardised one-step-ahead prediction errors (each log count minus its prediction from all earlier counts, divided by its predicted SD).",
      "If the model fits they scatter independently around 0 with ~95%% inside the dashed +/-2 lines.",
      "Ljung-Box test for autocorrelation p = %.2f; Shapiro-Wilk test for normality p = %.2f."), obj$gof$ljung_box_p, obj$gof$shapiro_p))
  }
  if (!is.null(obj$ppc)) {
    pp <- obj$ppc
    p1 <- fig_ppc(pp$T_rep[, 1], pp$T_obs[[1]], fig_label("ppc_sd", "SD of year-to-year change in log count", opts))
    p2 <- fig_ppc(pp$T_rep[, 2], pp$T_obs[[2]], fig_label("ppc_acf", "Lag-1 autocorrelation of changes", opts))
    p <- compose_figures(list(p1 + theme_eco(), p2 + theme_eco()), ncol = 2)
    attr(p, "fig_size") <- list(width = "double", aspect = 0.45)
    out$ppc <- with_caption(p, sprintf(paste(
      "Posterior predictive check: histograms of each statistic in 1000 count series replicated from the fitted model; vertical line = observed counts.",
      "(A) variability of year-to-year changes, (B) their lag-1 autocorrelation (counting error makes it negative).",
      "Bayesian p-values %.2f and %.2f; values below 0.05 or above 0.95 indicate misfit."), pp$p[[1]], pp$p[[2]]))
  }
  # ---- MCMC
  if (obj$model == "bayes") {
    keep <- c("r0", "b", "sigma_p", "sigma_o", if (length(obj$covariates)) paste0("beta[", seq_along(obj$covariates), "]"))
    ptr <- fig_trace(obj$samples, intersect(keep, colnames(obj$samples[[1]]))) +
      ggplot2::scale_x_continuous(n.breaks = 3) + ggplot2::facet_wrap(~ind, scales = "free_y", ncol = 3)
    out$trace <- with_caption(ts_widen(ptr, "double", 0.5),
                              "MCMC traces per chain for the top-level parameters (should overlap like hairy caterpillars).")
    dr <- stack_draws(obj$derived)
    pp <- setdiff(colnames(dr), c("equilibrium_N", grep("^growth_beta_", colnames(dr), value = TRUE)))
    short <- c(density_dependence_b = "Density dependence b", sd_process = "Process SD", sd_observation = "Observation SD",
               N_final = "Final abundance", realised_trend_pct_per_year = "Realised trend (%/yr)", growth_rate_r0 = "Growth rate r0")
    dr2 <- as.data.frame(dr[, pp, drop = FALSE])
    names(dr2) <- ifelse(names(dr2) %in% names(short), short[names(dr2)],
                         gsub("_", " ", sub("^growth_effect_pct_(.*)$", "Effect of \\1 (%/SD)", names(dr2))))
    names(dr2) <- vapply(names(dr2), function(n) paste(strwrap(n, 24), collapse = "\n"), "")
    ppo <- fig_posterior(dr2, names(dr2), ncol = 3) + ggplot2::scale_x_continuous(n.breaks = 3) + ggplot2::scale_y_continuous(n.breaks = 3)
    out$posterior <- with_caption(ts_widen(ppo, "double", 0.6),
                                  "Posterior densities of the reported quantities.")
  }
  out
}
