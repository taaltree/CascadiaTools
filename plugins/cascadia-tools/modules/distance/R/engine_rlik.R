# Hand-written likelihoods: an independent route used to cross-check
# Distance::ds and unmarked::distsamp. No package code is used.
#
# (1) CDS, conditional likelihood of the distances (Buckland et al. 2001 eq. 3.x):
#     line  f(x) = g(x) / int_0^w g,      point f(r) = r g(r) / int_0^w r g
#     binned data: multinomial with bin probabilities pi_j / sum(pi).
#     D = n * E[s] / (a * P_a); var from the R2/P2 encounter rate
#     (Fewster et al. 2009), the delta-method variance of P_a and the
#     group-size variance; log-normal interval.
# (2) HDS, Poisson counts per transect and bin (Royle et al. 2004):
#     y_ij ~ Poisson(D_i a_i pi_ij(sigma_i)), log D_i = X beta, log sigma_i = Z alpha.

ds_cds_nll <- function(par, x, key, w, survey, binned, cut = NULL, counts = NULL) {
  sigma <- exp(par[1]); b <- if (key == "hr") exp(par[2]) else NULL
  if (binned) {
    pj <- ds_bin_probs(key, sigma, cut, survey, b)
    pr <- pj / sum(pj)
    return(-sum(counts * log(pmax(pr, 1e-300))))
  }
  pa <- ds_pa(key, sigma, w, survey, b)
  g <- ds_g(x, key, sigma, b)
  if (survey == "line") -sum(log(pmax(g, 1e-300)) - log(pa * w))
  else -sum(log(pmax(g, 1e-300)) + log(2 * x / w^2) - log(pa))
}

ds_engine_rlik <- function(data, settings) {
  w <- data$truncation_m
  o <- data$obs_w
  binned <- data$distance_type == "binned"
  key <- settings$rlik_key
  cut <- if (binned) data$cutpoints[data$cutpoints <= w + 1e-9] else NULL
  cnt <- if (binned) as.numeric(table(factor(findInterval((o$distbegin + o$distend) / 2, cut, all.inside = TRUE), levels = seq_len(length(cut) - 1)))) else NULL
  x <- if (binned) NULL else pmax(o$distance, 1e-6)
  lo <- log(w / 100); hi <- log(w * 100)
  if (key == "hn") {
    # one parameter: bounded 1-d search (robust: BFGS can run off to sigma -> Inf on flat data)
    oo <- stats::optimize(function(p) ds_cds_nll(p, x = x, key = key, w = w, survey = data$survey, binned = binned, cut = cut, counts = cnt),
                          c(lo, hi), tol = 1e-10)
    op <- list(par = oo$minimum, value = oo$objective, convergence = as.integer(oo$minimum > hi - 0.01))
  } else {
    p0 <- c(log(stats::median(c(o$distance, w / 4))), log(2.5))
    op <- stats::optim(p0, ds_cds_nll, x = x, key = key, w = w, survey = data$survey, binned = binned,
                       cut = cut, counts = cnt, method = "L-BFGS-B", lower = c(lo, log(1.01)), upper = c(hi, log(50)))
  }
  H <- stats::optimHess(op$par, ds_cds_nll, x = x, key = key, w = w, survey = data$survey, binned = binned, cut = cut, counts = cnt)
  Vp <- tryCatch(solve(H), error = function(e) matrix(NA, length(op$par), length(op$par)))
  pa_fun <- function(par) ds_pa(key, exp(par[1]), w, data$survey, if (key == "hr") exp(par[2]))
  pa <- pa_fun(op$par)
  gr <- vapply(seq_along(op$par), function(j) { h <- 1e-5; e <- replace(numeric(length(op$par)), j, h)
    (pa_fun(op$par + e) - pa_fun(op$par - e)) / (2 * h) }, 1)
  pa_se <- sqrt(as.numeric(t(gr) %*% Vp %*% gr))
  # encounter rate
  tr <- data$transects
  n_k <- as.numeric(table(factor(o$transect, levels = tr$transect)))
  l_k <- tr$effort * data$n_visits
  n <- sum(n_k); L <- sum(l_k)
  er <- n / L; er_se <- sqrt(ds_er_var(n_k, l_k))
  gs <- ds_group_size(data)
  a_unit <- if (data$survey == "line") 2 * w / 1000 else pi * (w / 1000)^2   # km2 per km of line / per point visit
  D <- er * gs$mean / (a_unit * pa)
  cv <- sqrt((er_se / er)^2 + (pa_se / pa)^2 + (gs$se / gs$mean)^2)
  rows <- list(ds_lognorm_row("density_per_km2", D, D * cv))
  A <- ds_total_area(data)
  if (is.finite(A)) rows[[length(rows) + 1]] <- ds_lognorm_row("abundance", D * A, D * A * cv)
  rows[[length(rows) + 1]] <- ds_logit_row("detection_prob", pa, pa_se)
  r <- rows[[length(rows)]]
  rows[[length(rows) + 1]] <- if (data$survey == "line") ds_row("esw_m", pa * w, pa_se * w, r$lower * w, r$upper * w)
                              else ds_row("edr_m", w * sqrt(pa), w * pa_se / (2 * sqrt(pa)), w * sqrt(r$lower), w * sqrt(r$upper))
  s1 <- sqrt(Vp[1, 1])
  rows[[length(rows) + 1]] <- ds_row(if (key == "hn") "sigma_m" else "hr_scale_m", exp(op$par[1]), exp(op$par[1]) * s1,
                                     exp(op$par[1] - 1.96 * s1), exp(op$par[1] + 1.96 * s1))
  if (key == "hr") { s2 <- sqrt(Vp[2, 2])
    rows[[length(rows) + 1]] <- ds_row("hr_shape", exp(op$par[2]), exp(op$par[2]) * s2, exp(op$par[2] - 1.96 * s2), exp(op$par[2] + 1.96 * s2)) }
  if (isTRUE(data$clusters)) rows[[length(rows) + 1]] <- ds_lognorm_row("mean_group_size", gs$mean, gs$se)
  cds <- do.call(rbind, rows)

  # HDS Poisson likelihood
  hcut <- ds_hds_cut(data, settings)
  y <- ds_bin_counts(data, hcut)
  Xs <- stats::model.matrix(stats::as.formula(ds_rhs(data$formulas$density)), data$tz)
  Xd <- stats::model.matrix(stats::as.formula(ds_rhs(data$formulas$detection)), data$tz)
  a_i <- ds_cov_area(data, max(hcut))
  nb <- ncol(Xs); na <- ncol(Xd)
  nll <- function(par) {
    D_i <- exp(Xs %*% par[seq_len(nb)]); s_i <- exp(Xd %*% par[nb + seq_len(na)])
    b <- if (key == "hr") exp(par[nb + na + 1]) else NULL
    ll <- 0
    for (i in seq_len(nrow(y))) {
      mu <- D_i[i] * a_i[i] * ds_bin_probs(key, s_i[i], hcut, data$survey, b)
      ll <- ll + sum(stats::dpois(y[i, ], mu, log = TRUE))
    }
    if (!is.finite(ll)) 1e10 else -ll
  }
  st <- c(log(max(sum(y) / sum(a_i) / 0.5, 0.01)), rep(0, nb - 1), op$par[1], rep(0, na - 1), if (key == "hr") op$par[2])
  oh <- stats::optim(st, nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  if (oh$convergence != 0 || abs(oh$par[nb + 1]) > hi) oh <- stats::optim(st, nll, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-12))
  Hh <- stats::optimHess(oh$par, nll)
  Vh <- tryCatch(solve(Hh), error = function(e) matrix(NA, length(st), length(st)))
  sh <- sqrt(pmax(diag(Vh), 0))
  hds <- rbind(
    ds_row("density_at_average_site", exp(oh$par[1]) * gs$mean, exp(oh$par[1]) * gs$mean * sh[1],
           exp(oh$par[1] - 1.96 * sh[1]) * gs$mean, exp(oh$par[1] + 1.96 * sh[1]) * gs$mean),
    ds_row("hds_sigma_m", exp(oh$par[nb + 1]), exp(oh$par[nb + 1]) * sh[nb + 1],
           exp(oh$par[nb + 1] - 1.96 * sh[nb + 1]), exp(oh$par[nb + 1] + 1.96 * sh[nb + 1])))
  for (j in seq_len(nb)[-1]) hds <- rbind(hds,
    wald_rows(oh$par[j], sh[j], paste0("density_beta_", ds_clean(colnames(Xs)[j]))),
    wald_rows(oh$par[j], sh[j], paste0("rate_ratio_", ds_clean(colnames(Xs)[j])), transform = exp))
  for (j in seq_len(na)[-1]) hds <- rbind(hds, wald_rows(oh$par[nb + j], sh[nb + j], paste0("detection_beta_", ds_clean(colnames(Xd)[j]))))
  est <- rbind(cds, hds)
  conv <- op$convergence == 0 && oh$convergence == 0 && all(is.finite(c(pa_se, sh)))
  list(engine = "hand-written likelihood", estimates = est, key = key, converged = conv,
       notes = if (!conv) "hand-written likelihood: optimiser or Hessian problem" else character(),
       cds_par = op$par, hds_par = oh$par, nll_cds = op$value, nll_hds = oh$value)
}
