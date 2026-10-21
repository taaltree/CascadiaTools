# Restricted maximum likelihood (REML) for the exponential-growth state-space
# model (EGSS; Staples, Taper & Dennis 2004; Dennis et al. 2006; Humbert et
# al. 2009). Default estimator for model: exponential.
#
# Why: plain ML for the EGSS underestimates the process variance in short
# series (it is often estimated at exactly 0), and the trend interval then
# collapses (simulation coverage 0.68 for 25 years; VALIDATION.md). REML
# estimates the two variances from the first differences of the log counts,
# which do not depend on the initial state, after integrating the trend out:
#
#   w_i = y_{t_i} - y_{t_{i-1}}   (observed years t_1 < t_2 < ...; gap k_i = t_i - t_{i-1})
#   E[w_i] = mu k_i
#   Var(w_i) = k_i sd_process^2 + 2 sd_observation^2,  Cov(w_i, w_{i+1}) = -sd_observation^2
#
# REML log-likelihood (Patterson & Thompson; Harville 1977):
#   -1/2 [ log|S| + log(k' S^-1 k) + r' S^-1 r ],  r = w - mu_hat k
# mu_hat is the generalised-least-squares trend given the variances, with
# Var(mu_hat) = 1 / (k' S^-1 k). The reported trend and its 95% interval do
# NOT plug in the REML variances: they are the mean and quantiles of the
# mixture of these normals over a grid of the two SDs weighted by the
# restricted likelihood (flat prior; egss_grid / egss_mix_quantile below).
# The true trajectory is the Kalman smoother with a diffuse initial state,
# drawing (mu, variances) from their sampling distribution.

egss_diffs <- function(y) {
  tobs <- which(!is.na(y))
  list(w = diff(y[tobs]), k = diff(tobs), tobs = tobs)
}

egss_S <- function(k, sp, so) {
  m <- length(k)
  S <- diag(k * sp^2 + 2 * so^2, m)
  if (m > 1) { i <- seq_len(m - 1); S[cbind(i, i + 1)] <- -so^2; S[cbind(i + 1, i)] <- -so^2 }
  S
}

egss_gls <- function(w, k, sp, so) {
  S <- egss_S(k, sp, so)
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  Si_k <- backsolve(R, forwardsolve(t(R), k)); Si_w <- backsolve(R, forwardsolve(t(R), w))
  kSk <- sum(k * Si_k)
  mu <- sum(k * Si_w) / kSk
  r <- w - mu * k
  Si_r <- backsolve(R, forwardsolve(t(R), r))
  list(mu = mu, se = sqrt(1 / kSk), logdet = 2 * sum(log(diag(R))), kSk = kSk, q = sum(r * Si_r))
}

egss_negreml <- function(lth, w, k) {
  g <- egss_gls(w, k, exp(lth[1]), exp(lth[2]))
  if (is.null(g)) return(1e10)
  0.5 * (g$logdet + log(g$kSk) + g$q)
}

# Diffuse-start Kalman filter/smoother for the random walk with drift
# (first observed year has a flat prior on its state).
egss_smooth <- function(y, mu, sp, so) {
  n <- length(y); ap <- Pp <- af <- Pf <- numeric(n)
  big <- 1e7
  for (t in seq_len(n)) {
    if (t == 1) { a <- if (is.na(y[1])) 0 else y[1]; P <- big } else { a <- af[t - 1] + mu; P <- Pf[t - 1] + sp^2 }
    ap[t] <- a; Pp[t] <- P
    if (!is.na(y[t])) { Fv <- P + so^2; K <- P / Fv; af[t] <- a + K * (y[t] - a); Pf[t] <- P * (1 - K) }
    else { af[t] <- a; Pf[t] <- P }
  }
  xs <- af; Ps <- Pf
  if (n > 1) for (t in (n - 1):1) {
    J <- if (Pp[t + 1] > 0) Pf[t] / Pp[t + 1] else 0
    xs[t] <- af[t] + J * (xs[t + 1] - ap[t + 1])
    Ps[t] <- Pf[t] + J^2 * (Ps[t + 1] - Pp[t + 1])
  }
  list(xs = xs, Ps = pmax(Ps, 0))
}

egss_reml_fit <- function(y, n_draws = 1000, seed = 1, prob = 0.95) {
  ok <- which(!is.na(y)); y <- y[min(ok):max(ok)]; n <- length(y)
  d <- egss_diffs(y); w <- d$w; k <- d$k; m <- length(w)
  if (m < 4) stop("REML needs at least 5 years with counts")
  lo <- c(-8, -8); hi <- c(3, 3)
  vd <- stats::var(w / sqrt(k)); if (!is.finite(vd) || vd <= 0) vd <- 0.05
  obj <- function(l) egss_negreml(l, w, k)
  best <- NULL
  for (f in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
    s <- pmin(pmax(c(log(sqrt(f * vd)), log(sqrt((1 - f) * vd / 2))), lo + 1e-3), hi - 1e-3)
    o <- tryCatch(stats::nlminb(s, obj, lower = lo, upper = hi), error = function(e) NULL)
    if (!is.null(o) && is.finite(o$objective) && (is.null(best) || o$objective < best$objective)) best <- o
  }
  if (is.null(best)) stop("REML optimisation failed")
  nm <- c("log_sd_process", "log_sd_observation")
  th <- stats::setNames(best$par, nm)
  at_bound <- nm[abs(th - lo) < 1e-3 | abs(th - hi) < 1e-3]
  free <- setdiff(nm, at_bound)
  V <- matrix(0, 2, 2, dimnames = list(nm, nm)); pd <- TRUE
  if (length(free)) {
    H <- tryCatch(stats::optimHess(th[free], function(v) { t2 <- th; t2[free] <- v; obj(t2) }), error = function(e) NULL)
    Vf <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
    pd <- !is.null(Vf) && all(is.finite(Vf)) && all(diag(as.matrix(Vf)) > 0)
    if (pd) V[free, free] <- Vf
  }
  # profile REML interval for each log SD
  crit <- stats::qchisq(prob, 1) / 2; target <- -best$objective - crit
  prof <- function(j, v) {
    o <- stats::optimize(function(u) { l <- numeric(2); l[j] <- v; l[-j] <- u; obj(l) }, c(lo[-j], hi[-j]))
    -o$objective
  }
  pci <- lapply(1:2, function(j) {
    side <- function(dir) {
      bound <- if (dir < 0) lo[j] else hi[j]
      if (prof(j, bound) >= target) return(bound)
      tryCatch(stats::uniroot(function(v) prof(j, v) - target, sort(c(th[[j]], bound)), tol = 1e-4)$root,
               error = function(e) bound)
    }
    c(lower = side(-1), upper = side(1))
  })
  names(pci) <- nm
  sp <- exp(th[[1]]); so <- exp(th[[2]])
  if (abs(th[[1]] - lo[1]) < 1e-3) sp <- 0
  if (abs(th[[2]] - lo[2]) < 1e-3) so <- 0
  notes <- character()
  if (length(at_bound)) notes <- c(notes, sprintf("REML estimate on the boundary for %s (a standard deviation estimated at ~0): the data cannot separate process from observation error well; its interval starts at 0. Trend and abundance intervals still average over all variance splits the data allow (integrated restricted likelihood), so they stay honest.",
                                                 paste(at_bound, collapse = ", ")))
  # ---- integrate the variances out on a grid (see egss_grid) ----
  gp <- egss_grid(w, k, upper = 2.5 * max(sqrt(vd), 0.02))
  a <- (1 - prob) / 2
  mu_q <- egss_mix_quantile(gp, c(a, 1 - a))
  mu_mean <- sum(gp$p * gp$mu); mu_sd <- sqrt(sum(gp$p * (gp$se^2 + gp$mu^2)) - mu_mean^2)
  sdrow <- function(name, p, est) {
    ci <- exp(pci[[p]]); if (abs(pci[[p]][["lower"]] - lo[1]) < 1e-6) ci[1] <- 0
    se <- if (p %in% free && pd) est * sqrt(V[p, p]) else NA_real_
    ts_rows(name, est, se, ci[1], ci[2])
  }
  rows <- list(
    ts_rows("growth_rate_mu", mu_mean, mu_sd, mu_q[1], mu_q[2]),
    ts_rows("trend_pct_per_year", 100 * (exp(mu_mean) - 1), 100 * exp(mu_mean) * mu_sd, 100 * (exp(mu_q[1]) - 1), 100 * (exp(mu_q[2]) - 1)),
    sdrow("sd_process", "log_sd_process", sp), sdrow("sd_observation", "log_sd_observation", so))
  # ---- trajectory: draw (sd_p, sd_o) from the grid, mu | variances, states | all ----
  set.seed(seed)
  idx <- sample.int(length(gp$p), n_draws, replace = TRUE, prob = gp$p)
  X <- matrix(NA_real_, n_draws, n)
  TH <- matrix(NA_real_, n_draws, 3, dimnames = list(NULL, c("mu", nm)))
  for (i in seq_len(n_draws)) {
    j <- idx[i]
    spi <- max(gp$sp[j] + stats::runif(1, -0.5, 0.5) * gp$h[1], 1e-6)
    soi <- max(gp$so[j] + stats::runif(1, -0.5, 0.5) * gp$h[2], 1e-6)
    mui <- gp$mu[j] + gp$se[j] * stats::rnorm(1)
    TH[i, ] <- c(mui, log(spi), log(soi))
    sm <- egss_smooth(y, mui, spi, soi)
    X[i, ] <- sm$xs + sqrt(sm$Ps) * stats::rnorm(n)
  }
  traj <- data.frame(t = seq_len(n), estimate = exp(apply(X, 2, stats::median)),
                     lower = exp(apply(X, 2, stats::quantile, a)), upper = exp(apply(X, 2, stats::quantile, 1 - a)),
                     lower50 = exp(apply(X, 2, stats::quantile, 0.25)), upper50 = exp(apply(X, 2, stats::quantile, 0.75)))
  nf <- exp(X[, n])
  rows$nf <- ts_rows("N_final", stats::median(nf), stats::sd(nf), stats::quantile(nf, a, names = FALSE), stats::quantile(nf, 1 - a, names = FALSE))
  est <- do.call(rbind, rows); rownames(est) <- NULL
  kf <- list(theta = c(mu = mu_mean, th), form = "exponential", obs_error = TRUE, y = y, n = n,
             loglik = -best$objective, converged = best$convergence == 0 && gp$ok, notes = notes,
             at_bound = at_bound, aic = NA_real_, method = "REML",
             grid = list(sp_median = egss_wquant(gp$sp, gp$p, 0.5), so_median = egss_wquant(gp$so, gp$p, 0.5),
                         sp_mean = sum(gp$p * gp$sp), so_mean = sum(gp$p * gp$so),
                         sp_q = egss_wquant(gp$sp, gp$p, c(a, 1 - a)), so_q = egss_wquant(gp$so, gp$p, c(a, 1 - a))))
  if (!gp$ok) kf$notes <- c(kf$notes, "The restricted likelihood is not negligible at the edge of the variance grid; intervals may be too narrow.")
  list(kf = kf, estimates = est, trajectory = traj, state_draws = X, theta_draws = TH)
}

# Integrated restricted likelihood on a grid of (sd_process, sd_observation).
# exp(REML log-likelihood) is the likelihood of the differenced counts with mu
# integrated out under a flat prior; multiplying by a flat prior on the two
# SDs over [0, upper] and normalising gives their joint posterior, and
# mu | SDs is normal(GLS estimate, GLS SE). The marginal for mu is therefore
# an exact mixture of normals. This is what keeps trend intervals honest when
# a variance is estimated at 0: plain GLS/REML plug-in intervals then
# collapse (simulation coverage 0.81 at 25 years, VALIDATION.md).
egss_grid <- function(w, k, upper, ng = 60) {
  h <- upper / ng
  g <- (seq_len(ng) - 0.5) * h
  G <- expand.grid(sp = g, so = g)
  ll <- mu <- se <- numeric(nrow(G))
  for (i in seq_len(nrow(G))) {
    r <- egss_gls(w, k, G$sp[i], G$so[i])
    if (is.null(r)) { ll[i] <- -Inf; next }
    ll[i] <- -0.5 * (r$logdet + log(r$kSk) + r$q); mu[i] <- r$mu; se[i] <- r$se
  }
  p <- exp(ll - max(ll)); p <- p / sum(p)
  edge <- sum(p[G$sp > g[ng - 2] | G$so > g[ng - 2]])
  list(sp = G$sp, so = G$so, p = p, mu = mu, se = se, h = c(h, h), ok = edge < 0.01)
}

egss_mix_quantile <- function(gp, probs) {
  keep <- gp$p > 1e-10
  pm <- gp$p[keep] / sum(gp$p[keep]); m <- gp$mu[keep]; s <- gp$se[keep]
  cdf <- function(x) sum(pm * stats::pnorm(x, m, s))
  lo <- min(m - 8 * s); hi <- max(m + 8 * s)
  vapply(probs, function(q) stats::uniroot(function(x) cdf(x) - q, c(lo, hi), tol = 1e-10)$root, 1)
}

egss_wquant <- function(x, p, probs) {
  o <- order(x); cp <- cumsum(p[o])
  vapply(probs, function(q) x[o][which(cp >= q)[1]], 1)
}

# ---- several sites: pooled REML (model: marss, estimation: reml) ------------
# With one state per site (Z = I), Q and R "diagonal and equal" and B = I,
# the sites are independent given (sd_process, sd_observation), so the
# restricted likelihood of all first differences is the sum over sites of the
# EGSS terms above, with a shared trend (one mu; U "equal") or one mu per site
# (U "unequal"). The variances are integrated out on a grid exactly as for one
# site, so trend intervals average over every split of the noise the data
# allow. Plain ML (MARSS, Hessian intervals) under-covers here: trend 0.88,
# process SD 0.79 at 4 sites x 25 years (VALIDATION.md).

# Per-site GLS pieces at (sp, so): k'S^-1k, k'S^-1w, w'S^-1w, log|S|.
mreml_pieces <- function(D, sp, so) {
  out <- matrix(NA_real_, length(D), 4, dimnames = list(NULL, c("kSk", "kSw", "wSw", "logdet")))
  for (s in seq_along(D)) {
    S <- egss_S(D[[s]]$k, sp, so)
    R <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    zk <- forwardsolve(t(R), D[[s]]$k); zw <- forwardsolve(t(R), D[[s]]$w)
    out[s, ] <- c(sum(zk^2), sum(zk * zw), sum(zw^2), 2 * sum(log(diag(R))))
  }
  out
}

# REML log-likelihood + GLS trend(s) and SE(s).
mreml_eval <- function(D, sp, so, shared) {
  pc <- mreml_pieces(D, sp, so)
  if (is.null(pc)) return(NULL)
  if (shared) {
    kSk <- sum(pc[, "kSk"]); mu <- sum(pc[, "kSw"]) / kSk
    q <- sum(pc[, "wSw"]) - sum(pc[, "kSw"])^2 / kSk
    list(ll = -0.5 * (sum(pc[, "logdet"]) + log(kSk) + q), mu = mu, se = sqrt(1 / kSk))
  } else {
    mu <- pc[, "kSw"] / pc[, "kSk"]
    q <- sum(pc[, "wSw"] - pc[, "kSw"]^2 / pc[, "kSk"])
    list(ll = -0.5 * (sum(pc[, "logdet"]) + sum(log(pc[, "kSk"])) + q), mu = mu, se = sqrt(1 / pc[, "kSk"]))
  }
}

marss_reml_fit <- function(Ylog, shared = TRUE, n_draws = 1000, seed = 1, prob = 0.95, ng = 60) {
  S <- nrow(Ylog); n <- ncol(Ylog)
  D <- lapply(seq_len(S), function(s) egss_diffs(Ylog[s, ]))
  nd <- vapply(D, function(d) length(d$w), 1L)
  if (any(nd < 2)) stop("each site needs at least 3 years with counts for model: marss")
  if (sum(nd) < 6) stop("too few counts for REML")
  wk <- unlist(lapply(D, function(d) d$w / sqrt(d$k)))
  vd <- stats::var(wk); if (!is.finite(vd) || vd <= 0) vd <- 0.05
  lo <- c(-8, -8); hi <- c(3, 3)
  obj <- function(l) { r <- mreml_eval(D, exp(l[1]), exp(l[2]), shared); if (is.null(r)) 1e10 else -r$ll }
  best <- NULL
  for (f in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
    s0 <- pmin(pmax(c(log(sqrt(f * vd)), log(sqrt((1 - f) * vd / 2))), lo + 1e-3), hi - 1e-3)
    o <- tryCatch(stats::nlminb(s0, obj, lower = lo, upper = hi), error = function(e) NULL)
    if (!is.null(o) && is.finite(o$objective) && (is.null(best) || o$objective < best$objective)) best <- o
  }
  if (is.null(best)) stop("REML optimisation failed")
  nm <- c("log_sd_process", "log_sd_observation")
  th <- stats::setNames(best$par, nm)
  at_bound <- nm[abs(th - lo) < 1e-3]
  crit <- stats::qchisq(prob, 1) / 2; target <- -best$objective - crit
  prof <- function(j, v) -stats::optimize(function(u) { l <- numeric(2); l[j] <- v; l[-j] <- u; obj(l) }, c(lo[-j], hi[-j]))$objective
  pci <- lapply(1:2, function(j) {
    side <- function(dir) {
      bound <- if (dir < 0) lo[j] else hi[j]
      if (prof(j, bound) >= target) return(bound)
      tryCatch(stats::uniroot(function(v) prof(j, v) - target, sort(c(th[[j]], bound)), tol = 1e-4)$root,
               error = function(e) bound)
    }
    c(lower = side(-1), upper = side(1))
  })
  names(pci) <- nm
  sp <- if (abs(th[[1]] - lo[1]) < 1e-3) 0 else exp(th[[1]])
  so <- if (abs(th[[2]] - lo[2]) < 1e-3) 0 else exp(th[[2]])
  # grid over the two SDs (flat prior), mu | SDs normal
  upper <- 2.5 * max(sqrt(vd), 0.02); h <- upper / ng; g <- (seq_len(ng) - 0.5) * h
  G <- expand.grid(sp = g, so = g)
  nmu <- if (shared) 1 else S
  ll <- numeric(nrow(G)); MU <- SE <- matrix(NA_real_, nrow(G), nmu)
  for (i in seq_len(nrow(G))) {
    r <- mreml_eval(D, G$sp[i], G$so[i], shared)
    if (is.null(r)) { ll[i] <- -Inf; next }
    ll[i] <- r$ll; MU[i, ] <- r$mu; SE[i, ] <- r$se
  }
  p <- exp(ll - max(ll)); p <- p / sum(p)
  edge_ok <- sum(p[G$sp > g[ng - 2] | G$so > g[ng - 2]]) < 0.01
  a <- (1 - prob) / 2
  mix_q <- function(j) egss_mix_quantile(list(p = p, mu = MU[, j], se = SE[, j]), c(a, 1 - a))
  rows <- list()
  lab <- if (shared) "" else paste0("_", rownames(Ylog))
  for (j in seq_len(nmu)) {
    mm <- sum(p * MU[, j]); ms <- sqrt(sum(p * (SE[, j]^2 + MU[, j]^2)) - mm^2); qq <- mix_q(j)
    rows[[length(rows) + 1]] <- ts_rows(paste0("growth_rate_mu", lab[j]), mm, ms, qq[1], qq[2])
    rows[[length(rows) + 1]] <- ts_rows(paste0("trend_pct_per_year", lab[j]), 100 * (exp(mm) - 1), 100 * exp(mm) * ms,
                                        100 * (exp(qq[1]) - 1), 100 * (exp(qq[2]) - 1))
  }
  sdrow <- function(name, j, est) {
    ci <- exp(pci[[j]]); if (abs(pci[[j]][["lower"]] - lo[1]) < 1e-6) ci[1] <- 0
    ts_rows(name, est, NA_real_, ci[1], ci[2])
  }
  rows$sp <- sdrow("sd_process", 1, sp); rows$so <- sdrow("sd_observation", 2, so)
  est <- do.call(rbind, rows); rownames(est) <- NULL
  # trajectories: draw grid cell, mu | cell, states | all, per site
  set.seed(seed)
  idx <- sample.int(length(p), n_draws, replace = TRUE, prob = p)
  X <- array(NA_real_, c(n_draws, S, n))
  spd <- numeric(n_draws); mud <- matrix(NA_real_, n_draws, nmu)
  for (i in seq_len(n_draws)) {
    j <- idx[i]
    spi <- max(G$sp[j] + stats::runif(1, -0.5, 0.5) * h, 1e-6); soi <- max(G$so[j] + stats::runif(1, -0.5, 0.5) * h, 1e-6)
    mui <- MU[j, ] + SE[j, ] * stats::rnorm(nmu)
    spd[i] <- spi; mud[i, ] <- mui
    for (s in seq_len(S)) {
      sm <- egss_smooth(Ylog[s, ], mui[min(s, nmu)], spi, soi)
      X[i, s, ] <- sm$xs + sqrt(sm$Ps) * stats::rnorm(n)
    }
  }
  notes <- character()
  if (length(at_bound)) notes <- c(notes, sprintf("REML estimate on the boundary for %s (a standard deviation estimated at ~0); its interval starts at 0. Trend intervals still average over all variance splits the data allow.", paste(at_bound, collapse = ", ")))
  if (!edge_ok) notes <- c(notes, "The restricted likelihood is not negligible at the edge of the variance grid; intervals may be too narrow.")
  list(estimates = est, X = X, sp_draws = spd, mu_draws = mud, notes = notes,
       converged = best$convergence == 0 && edge_ok, loglik_reml = -best$objective)
}
