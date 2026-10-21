# Count-based PVA (settings$model = "count").
#
# Two estimators of the stochastic exponential-growth (diffusion) model
#   log N_{t+s} = log N_t + mu s + process noise N(0, sigma2 s)
# from a series of annual (or irregular) counts:
#   diffusion   - Dennis, Munholland & Scott (1991): counts taken as exact;
#                 regression of log ratios on sqrt(interval) through the origin
#                 (Morris & Doak 2002, ch. 3).
#   state_space - exponential-growth state-space (EGSS) model with
#                 observation error tau2 on log counts, REML estimates of
#                 sigma2 and tau2 from second differences (Staples, Taper &
#                 Dennis 2004; Humbert et al. 2009, App. eq. A21-A28), GLS mu.
#
# Quasi-extinction: probability that true abundance falls to or below the
# threshold Nq in at least one annual census within h years. The classic
# inverse-Gaussian first-passage formula (Dennis et al. 1991; Holmes 2004
# eq. 6) assumes continuous monitoring; with annual censuses it over-states
# risk. The default "discrete" correction shifts the barrier by
# 0.5826 * sigma (Broadie, Glasserman & Kou 1997; Siegmund 1979), which we
# checked against brute-force simulation (VALIDATION.md).

# Resolve a data path from the config: absolute paths as given, relative ones against the config folder.
pva_path <- function(config, f) if (grepl("^(/|[A-Za-z]:[/\\\\]|~)", f)) path.expand(f) else file.path(config$.dir %||% ".", f)

PVA_BGK <- 0.5825971579390106   # -zeta(1/2)/sqrt(2*pi)

# P(first passage below log-distance d within time Tt). Vectorised over all
# arguments (recycled). d <= 0 means already at/below the threshold.
pva_G <- function(Tt, mu, s2, d) {
  n <- max(length(Tt), length(mu), length(s2), length(d))
  Tt <- rep_len(Tt, n); mu <- rep_len(mu, n); s2 <- rep_len(s2, n); d <- rep_len(d, n)
  out <- numeric(n)
  hit <- d <= 0
  out[hit] <- 1
  ok <- !hit & s2 > 0 & Tt > 0
  if (any(ok)) {
    s <- sqrt(s2[ok] * Tt[ok])
    a <- stats::pnorm((-d[ok] - mu[ok] * Tt[ok]) / s)
    lb <- -2 * mu[ok] * d[ok] / s2[ok] + stats::pnorm((-d[ok] + mu[ok] * Tt[ok]) / s, log.p = TRUE)
    out[ok] <- pmin(1, a + exp(pmin(lb, 0)))
  }
  # no process noise: deterministic decline
  det <- !hit & s2 <= 0
  out[det] <- as.numeric(mu[det] * Tt[det] <= -d[det])
  out
}

# Quasi-extinction probability by horizon(s) h with optional discrete-census
# correction. d = log(N_current / Nq).
pva_qext <- function(h, mu, s2, d, correction = "discrete") {
  if (identical(correction, "discrete")) d <- ifelse(d > 0, d + PVA_BGK * sqrt(pmax(s2, 0)), d)
  pva_G(h, mu, s2, d)
}

# Brute-force truth: annual random walk on log scale started at x0 (log N),
# P(any of x_1..x_h <= log Nq). Returns probabilities for every horizon in hs.
pva_qext_sim <- function(hs, mu, s2, x0, Nq, npaths = 2e5, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  hmax <- max(hs)
  x <- rep(x0, npaths); dead <- rep(x0 <= log(Nq), npaths)
  out <- numeric(hmax)
  for (t in seq_len(hmax)) {
    x <- x + mu + sqrt(s2) * stats::rnorm(npaths)
    dead <- dead | x <= log(Nq)
    out[t] <- mean(dead)
  }
  out[hs]
}

# Exact (numerical) truth: discrete-time first-passage probability of the
# annual Gaussian random walk, by recursion on the density of surviving paths
# over a fine grid of log N above the threshold. The mass absorbed in each
# year is accumulated directly, so tiny probabilities keep their relative
# precision (a Monte Carlo truth with 2e5 paths cannot resolve risks < 5e-6).
# Agrees with brute-force simulation (dev/pva/crosscheck_exact.R).
pva_qext_exact <- function(hs, mu, s2, x0, Nq, cells_per_sd = 25, max_cells = 4000) {
  b <- log(Nq); H <- max(hs)
  if (x0 <= b) return(rep(1, length(hs)))
  s <- sqrt(s2)
  if (s <= 0) return(as.numeric(x0 + mu * hs <= b))
  U <- max(x0, x0 + mu * H) + 9 * s * sqrt(H)
  n <- min(max_cells, ceiling((U - b) / (s / cells_per_sd)))
  edges <- seq(b, U, length.out = n + 1); xm <- (edges[-1] + edges[-(n + 1)]) / 2
  w <- diff(stats::pnorm((edges - x0 - mu) / s))          # after year 1
  ab <- numeric(H); ab[1] <- stats::pnorm((b - x0 - mu) / s)
  if (H > 1) {
    P <- stats::pnorm(outer(edges, xm, function(e, x) (e - x - mu) / s))
    K <- P[-1, , drop = FALSE] - P[-(n + 1), , drop = FALSE]   # K[i, j]: cell j -> cell i
    a <- P[1, ]                                                 # cell j -> absorbed
    for (t in 2:H) { ab[t] <- sum(w * a); w <- as.vector(K %*% w) }
  }
  pmin(cumsum(ab), 1)[hs]
}

# ---- simulation -------------------------------------------------------------

count_simulate <- function(params, design, seed) {
  set.seed(seed)
  n <- design$n_years
  mu <- params$mu; s2 <- params$sigma2; tau2 <- params$tau2 %||% 0
  x <- log(design$N0) + c(0, cumsum(mu + sqrt(s2) * stats::rnorm(n - 1)))
  y <- x + if (tau2 > 0) sqrt(tau2) * stats::rnorm(n) else 0
  cnt <- pmax(1, round(exp(y)))
  yrs <- design$first_year + seq_len(n) - 1
  keep <- rep(TRUE, n)
  pm <- design$p_missing %||% 0
  if (pm > 0 && n > 4) keep[2:(n - 1)] <- stats::runif(n - 2) > pm
  hs <- design$horizons
  Nq <- design$threshold
  pq <- if (identical(design$truth_method, "simulation"))
    pva_qext_sim(hs, mu, s2, x[n], Nq, npaths = design$truth_paths %||% 2e5, seed = seed + 7919)
  else pva_qext_exact(hs, mu, s2, x[n], Nq)
  truth <- c(growth_rate_mu = mu, lambda_stochastic = exp(mu), process_variance_sigma2 = s2,
             stats::setNames(pq, paste0("p_qext_", hs, "yr")))
  # N_current is an estimand only for the state-space model; the diffusion
  # method takes the last count as known (no interval to check).
  if (tau2 > 0 || identical(design$method, "state_space"))
    truth <- c(truth, N_current = exp(x[n]), obs_variance_tau2 = tau2)
  list(data = list(model = "count", df = data.frame(year = yrs[keep], count = cnt[keep]),
                   threshold = Nq, horizons = hs,
                   meta = list(species = "simulated population", units = "individuals")),
       truth = truth)
}

# ---- Dennis et al. (1991) diffusion estimates -------------------------------

count_dennis <- function(y, tt) {
  s <- diff(tt); q <- length(s)
  xi <- diff(y) / sqrt(s)
  mu <- sum(diff(y)) / sum(s)
  s2 <- sum((xi - mu * sqrt(s))^2) / (q - 1)
  list(mu = mu, s2 = s2, q = q, tq = sum(s), se_mu = sqrt(s2 / sum(s)), resid = xi - mu * sqrt(s))
}

# ---- EGSS (state-space) REML -------------------------------------------------

egss_V <- function(tt, s2, tau2) {
  t0 <- tt - tt[1]
  s2 * outer(t0, t0, pmin) + diag(tau2, length(tt))
}

egss_negll_reml <- function(theta, y, tt) {
  s2 <- exp(theta[1]); tau2 <- exp(theta[2])
  q <- length(y) - 1; s <- diff(tt)
  V <- egss_V(tt, s2, tau2)
  D1 <- cbind(-diag(1 / s, q), 0) + cbind(0, diag(1 / s, q))
  D2 <- cbind(-diag(1, q - 1), 0) + cbind(0, diag(1, q - 1))
  M <- D2 %*% D1
  Phi <- M %*% V %*% t(M)
  u <- as.vector(M %*% y)
  R <- tryCatch(chol(Phi), error = function(e) NULL)
  if (is.null(R)) return(1e10)
  z <- backsolve(R, u, transpose = TRUE)
  ((q - 1) / 2) * log(2 * pi) + sum(log(diag(R))) + 0.5 * sum(z^2)
}

egss_fit <- function(y, tt, starts = NULL) {
  q <- length(y) - 1; s <- diff(tt)
  if (q < 4) stop("The state-space model needs at least 5 counts.", call. = FALSE)
  # EGOE / EGPN starting values (Humbert et al. 2009 appendix code)
  b <- stats::coef(stats::lm(y ~ tt))
  tsq_oe <- sum((y - b[1] - b[2] * tt)^2) / (q - 1)
  dn <- count_dennis(y, tt)
  st <- rbind(c(log(max(dn$s2 / 2, 1e-4)), log(max(tsq_oe / 2, 1e-4))),
              c(log(max(dn$s2, 1e-4)), log(1e-3)),
              c(log(max(dn$s2 / 10, 1e-4)), log(max(tsq_oe, 1e-4))),
              starts)
  best <- NULL
  for (i in seq_len(nrow(st))) {
    o <- tryCatch(stats::optim(st[i, ], egss_negll_reml, y = y, tt = tt, method = "L-BFGS-B",
                               lower = c(-16, -16), upper = c(3, 3)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) stop("State-space REML optimisation failed.", call. = FALSE)
  s2 <- exp(best$par[1]); tau2 <- exp(best$par[2])
  if (s2 < 1e-6) s2 <- 0
  if (tau2 < 1e-6) tau2 <- 0
  # GLS mu from scaled first differences (A26-A27). V is singular when
  # tau2 = 0 (the first state is then known exactly), so floor it for algebra.
  w <- diff(y) / s
  V <- egss_V(tt, max(s2, 1e-10), max(tau2, 1e-8))
  D1 <- cbind(-diag(1 / s, q), 0) + cbind(0, diag(1 / s, q))
  V1 <- D1 %*% V %*% t(D1)
  V1i <- solve(V1)
  j <- rep(1, q)
  vmu <- 1 / sum(V1i)
  mu <- vmu * sum(V1i %*% w)
  # initial condition (A28) and smoothed final state
  t0 <- tt - tt[1]
  Vi <- solve(V)
  jj <- rep(1, q + 1)
  beta <- sum(Vi %*% (y - t0 * mu)) / sum(Vi)
  cvec <- max(s2, 0) * t0
  r <- y - beta - mu * t0
  xT <- beta + mu * t0[q + 1] + sum(cvec * (Vi %*% r))
  PT <- max(max(s2, 0) * t0[q + 1] - sum(cvec * (Vi %*% cvec)), 0)
  # smoothed states at every observation time (for plots)
  C <- max(s2, 0) * outer(t0, t0, pmin)
  xs <- beta + mu * t0 + as.vector(C %*% (Vi %*% r))
  Ps <- pmax(diag(C - C %*% Vi %*% C), 0)
  list(mu = mu, se_mu = sqrt(vmu), s2 = s2, tau2 = tau2, beta = beta, xT = xT, PT = PT,
       xs = xs, Ps = Ps, negll = best$value, convergence = best$convergence)
}

egss_sim <- function(tt, mu, s2, tau2, beta, states = FALSE) {
  t0 <- tt - tt[1]
  inc <- c(0, mu * diff(t0) + sqrt(s2 * diff(t0)) * stats::rnorm(length(tt) - 1))
  x <- beta + cumsum(inc)
  y <- x + sqrt(tau2) * stats::rnorm(length(tt))
  if (states) list(y = y, x = x) else y
}

# ---- EGSS posterior by exact grid integration ------------------------------
#
# y = X b + W + e with X = [1, t - t1], b = (x at first count, mu), W a
# Brownian motion with variance sigma2 per year, e iid N(0, tau2). With flat
# priors on b, the marginal likelihood of (sigma2, tau2) is the REML
# likelihood (Harville 1974), so the joint posterior is
#   p(sigma, tau | y) ∝ L_REML(sigma2, tau2) p(sigma) p(tau)        (2-D grid)
#   b | sigma, tau, y   ~ N(b_GLS, (X'V^-1 X)^-1)                     (exact)
#   x | b, sigma, tau, y ~ N(Xb + C V^-1 (y - Xb), C - C V^-1 C)       (exact)
# with uniform priors on the standard deviations sigma and tau (Gelman 2006).
# No MCMC: the 2-D integral is done on a grid, everything else is Gaussian.
# Unlike a parametric bootstrap from the REML fit, the posterior does not
# collapse when REML puts sigma2 (or tau2) on the boundary at 0, which
# happened in ~10% of 25-year series and gave 0.90 coverage for mu.

egss_cell <- function(y, X, t0, s2, tau2, full = FALSE) {
  n <- length(y)
  C <- s2 * outer(t0, t0, pmin)
  V <- C + diag(tau2, n)
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  Ry <- backsolve(R, y, transpose = TRUE); RX <- backsolve(R, X, transpose = TRUE)
  XtViX <- crossprod(RX); XtViy <- crossprod(RX, Ry)
  Rx <- tryCatch(chol(XtViX), error = function(e) NULL)
  if (is.null(Rx)) return(NULL)
  b <- backsolve(Rx, backsolve(Rx, XtViy, transpose = TRUE))
  rss <- sum(Ry^2) - sum(XtViy * b)
  ll <- -sum(log(diag(R))) - sum(log(diag(Rx))) - 0.5 * rss
  if (!full) return(ll)
  Cb <- chol2inv(Rx)
  RC <- backsolve(R, C, transpose = TRUE)               # R^-T C
  K <- t(backsolve(R, RC))                              # C V^-1
  P <- C - crossprod(RC)                                # C - C V^-1 C
  P <- (P + t(P)) / 2
  Lp <- tryCatch(t(chol(P + diag(1e-12, n), pivot = FALSE)), error = function(e) {
    ev <- eigen(P, symmetric = TRUE); ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), n)
  })
  list(ll = ll, b = as.vector(b), Lb = t(chol(Cb)), K = K, Lp = Lp)
}

egss_posterior <- function(y, tt, ndraw = 4000, ngrid = 60, sd_max = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  if (n < 5) stop("The state-space model needs at least 5 counts.", call. = FALSE)
  t0 <- tt - tt[1]; X <- cbind(1, t0)
  inc <- diff(y) / sqrt(diff(tt))
  U <- sd_max %||% max(4 * stats::sd(inc), 0.05)
  grid_ll <- function(sg, tg) {
    outer(seq_along(sg), seq_along(tg), Vectorize(function(i, j) {
      v <- egss_cell(y, X, t0, sg[i]^2, tg[j]^2); if (is.null(v)) -Inf else v
    }))
  }
  # coarse pass over [0, U]^2 to find where the mass is, then a fine grid
  m0 <- 24
  sg <- (seq_len(m0) - 0.5) * U / m0; tg <- sg
  L0 <- grid_ll(sg, tg)
  keep <- L0 > max(L0) - 14
  ri <- range(which(apply(keep, 1, any))); ci <- range(which(apply(keep, 2, any)))
  slo <- if (ri[1] <= 2) 0 else sg[ri[1] - 1] - U / (2 * m0)
  shi <- min(U, sg[min(ri[2] + 1, m0)] + U / (2 * m0))
  tlo <- if (ci[1] <= 2) 0 else tg[ci[1] - 1] - U / (2 * m0)
  thi <- min(U, tg[min(ci[2] + 1, m0)] + U / (2 * m0))
  ds <- (shi - slo) / ngrid; dtau <- (thi - tlo) / ngrid
  sg <- slo + (seq_len(ngrid) - 0.5) * ds; tg <- tlo + (seq_len(ngrid) - 0.5) * dtau
  L <- grid_ll(sg, tg)
  w <- exp(L - max(L)); w <- w / sum(w)
  edge <- sum(w[ngrid, ]) + sum(w[, ngrid])
  cells <- sample.int(length(w), ndraw, replace = TRUE, prob = as.vector(w))
  ii <- (cells - 1) %% ngrid + 1; jj <- (cells - 1) %/% ngrid + 1
  sig <- sg[ii] + (stats::runif(ndraw) - 0.5) * ds
  tau <- tg[jj] + (stats::runif(ndraw) - 0.5) * dtau
  mu <- b1 <- numeric(ndraw); xs <- matrix(NA_real_, ndraw, n)
  for (cl in unique(cells)) {
    k <- which(cells == cl)
    f <- egss_cell(y, X, t0, sg[ii[k[1]]]^2, tg[jj[k[1]]]^2, full = TRUE)
    bb <- matrix(f$b, length(k), 2, byrow = TRUE) + t(f$Lb %*% matrix(stats::rnorm(2 * length(k)), 2))
    mean_x <- bb %*% t(X) + t(f$K %*% (y - X %*% t(bb)))
    xs[k, ] <- mean_x + t(f$Lp %*% matrix(stats::rnorm(n * length(k)), n))
    mu[k] <- bb[, 2]; b1[k] <- bb[, 1]
  }
  list(draws = data.frame(mu = mu, s2 = sig^2, tau2 = tau^2, xT = xs[, n]), states = xs,
       grid = list(sigma = sg, tau = tg, post = w), edge_mass = edge, U = U)
}

# Bias-corrected (BC) percentile remap of bootstrap replicates (Efron &
# Tibshirani 1993 ch. 14): returns the draws moved to quantile
# pnorm(2 z0 + qnorm(u)), z0 = qnorm(share of replicates below the estimate).
# Corrects for the median bias of REML variance estimates near zero.
pva_bc_remap <- function(b, est) {
  B <- length(b)
  p0 <- mean(b < est) + 0.5 * mean(b == est)
  p0 <- min(max(p0, 1 / (B + 1)), B / (B + 1))
  z0 <- stats::qnorm(p0)
  u <- (rank(b, ties.method = "random") - 0.5) / B
  stats::quantile(b, stats::pnorm(2 * z0 + stats::qnorm(u)), names = FALSE, type = 7)
}

# ---- fit -------------------------------------------------------------------

count_fit <- function(data, st) {
  df <- data$df[order(data$df$year), ]
  y <- log(df$count); tt <- df$year
  Nq <- st$threshold %||% data$threshold
  hs <- sort(unique(as.integer(st$horizons %||% data$horizons)))
  if (is.null(Nq)) stop("settings: threshold (quasi-extinction threshold, in count units) is required.", call. = FALSE)
  method <- st$method %||% "diffusion"
  corr <- st$correction %||% "discrete"
  nb <- st$nboot %||% 2000
  set.seed(st$seed %||% 1)
  tmax <- max(c(hs, st$tmax %||% max(hs)))
  notes <- character()
  if (method == "diffusion") {
    dn <- count_dennis(y, tt)
    mu <- dn$mu; s2 <- dn$s2; q <- dn$q; tq <- dn$tq
    xT <- y[length(y)]; PT <- 0
    tcrit <- stats::qt(0.975, q - 1)
    # parametric bootstrap from the exact sampling distributions
    s2b <- s2 * stats::rchisq(nb, q - 1) / (q - 1)
    mub <- mu + sqrt(s2b / tq) * stats::rnorm(nb)
    xTb <- rep(xT, nb)
    rows <- rbind(
      data.frame(param = "growth_rate_mu", estimate = mu, sd = dn$se_mu,
                 lower = mu - tcrit * dn$se_mu, upper = mu + tcrit * dn$se_mu),
      data.frame(param = "lambda_stochastic", estimate = exp(mu), sd = exp(mu) * dn$se_mu,
                 lower = exp(mu - tcrit * dn$se_mu), upper = exp(mu + tcrit * dn$se_mu)),
      data.frame(param = "process_variance_sigma2", estimate = s2, sd = s2 * sqrt(2 / (q - 1)),
                 lower = (q - 1) * s2 / stats::qchisq(0.975, q - 1), upper = (q - 1) * s2 / stats::qchisq(0.025, q - 1)))
    resid <- dn$resid
    ss <- NULL
  } else if (method == "state_space") {
    # Bayesian EGSS by exact grid integration (egss_posterior); REML kept for reference
    nd <- st$ndraw_ss %||% 4000
    U <- NULL
    for (k in 1:3) {
      pst <- egss_posterior(y, tt, ndraw = nd, sd_max = U, seed = (st$seed %||% 1) + k)
      if (pst$edge_mass < 0.005) break
      U <- 2 * pst$U
    }
    if (pst$edge_mass >= 0.005) notes <- c(notes, sprintf("%.1f%% of the posterior mass of the standard deviations is at the edge of the integration grid; results may be truncated.", 100 * pst$edge_mass))
    reml <- tryCatch(egss_fit(y, tt), error = function(e) NULL)
    D <- pst$draws
    mub <- D$mu; s2b <- D$s2; xTb <- D$xT; tau2b <- D$tau2
    mu <- stats::median(mub); s2 <- stats::median(s2b); xT <- stats::median(xTb); PT <- stats::var(xTb)
    pq <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE)
    rw <- function(p, v, f = identity) data.frame(param = p, estimate = f(stats::median(v)), sd = stats::sd(f(v)),
                                                  lower = f(pq(v)[1]), upper = f(pq(v)[2]))
    rows <- rbind(rw("growth_rate_mu", mub), rw("lambda_stochastic", mub, exp),
                  rw("process_variance_sigma2", s2b), rw("obs_variance_tau2", tau2b), rw("N_current", xTb, exp))
    # Variance point estimates: REML (Staples et al. 2004), which is close to
    # unbiased; the posterior median under uniform-SD priors overestimates
    # sigma2 by ~25% in 25-year series (VALIDATION.md). Intervals stay the
    # posterior credible intervals (they have nominal coverage; REML-based
    # intervals collapse when REML hits 0).
    if (!is.null(reml)) {
      rows$estimate[rows$param == "process_variance_sigma2"] <- reml$s2
      rows$estimate[rows$param == "obs_variance_tau2"] <- reml$tau2
    }
    resid <- NULL
    xs <- pst$states
    ss <- list(xs = apply(xs, 2, stats::median), xs_lo = apply(xs, 2, stats::quantile, 0.025, names = FALSE),
               xs_hi = apply(xs, 2, stats::quantile, 0.975, names = FALSE), Ps = apply(xs, 2, stats::var),
               tau2 = stats::median(tau2b), reml = if (is.null(reml)) NULL else reml[c("mu", "se_mu", "s2", "tau2", "xT")],
               edge_mass = pst$edge_mass, convergence = as.integer(pst$edge_mass >= 0.005), grid = pst$grid)
    if (!is.null(reml)) {
      if (reml$s2 == 0) notes <- c(notes, "The REML point estimate of the process variance is 0 (all year-to-year variation could be counting error); the Bayesian intervals reported here still allow for real process variation, which a REML plug-in would not.")
      if (reml$tau2 == 0) notes <- c(notes, "The REML point estimate of the counting-error variance is 0; the diffusion method (counts taken as exact) is then a reasonable alternative.")
    }
  } else stop("settings: method must be 'diffusion' or 'state_space'", call. = FALSE)

  d <- xT - log(Nq)
  tgrid <- seq_len(tmax)
  cdf_hat <- pva_qext(tgrid, mu, s2, d, corr)
  cdf_b <- vapply(tgrid, function(h) pva_qext(h, mub, s2b, xTb - log(Nq), corr), numeric(length(mub)))
  cdf_b <- matrix(cdf_b, ncol = length(tgrid))
  # state-space: posterior median of the risk (plug-in at the median parameters
  # would ignore the uncertainty in the current abundance)
  if (method == "state_space") cdf_hat <- apply(cdf_b, 2, stats::median)
  qrow <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE)
  for (h in hs) {
    ci <- qrow(cdf_b[, h])
    rows <- rbind(rows, data.frame(param = paste0("p_qext_", h, "yr"), estimate = cdf_hat[h],
                                   sd = stats::sd(cdf_b[, h]), lower = ci[1], upper = ci[2]))
  }
  if (method == "diffusion")
    rows <- rbind(rows, data.frame(param = "N_current", estimate = exp(xT), sd = 0, lower = exp(xT), upper = exp(xT)))
  # Fieberg & Ellner (2000) precision warning
  span <- diff(range(tt))
  long <- hs[hs > 0.2 * span]
  if (length(long))
    notes <- c(notes, sprintf("Horizons %s yr exceed 20%% of the %d-year monitoring span: forecasts this far ahead are typically too imprecise to rank risk (Fieberg & Ellner 2000). Read the interval, not the point estimate.",
                              paste(long, collapse = ", "), round(span)))
  # Simple model checks on the process (diffusion residuals or smoothed states)
  chk <- count_checks(y, tt, if (!is.null(ss)) ss$xs else NULL)
  notes <- c(notes, chk$notes)
  cdf_tab <- data.frame(year = tgrid, estimate = cdf_hat,
                        lower = apply(cdf_b, 2, stats::quantile, 0.025, names = FALSE),
                        upper = apply(cdf_b, 2, stats::quantile, 0.975, names = FALSE),
                        lower50 = apply(cdf_b, 2, stats::quantile, 0.25, names = FALSE),
                        upper50 = apply(cdf_b, 2, stats::quantile, 0.75, names = FALSE))
  # Classic continuous-monitoring diffusion value at each horizon, for reference
  classic <- pva_qext(hs, mu, s2, d, "none")
  rows$rhat <- NA_real_; rows$ess_bulk <- NA_real_
  list(estimates = rows,
       diagnostics = list(converged = is.null(ss) || ss$convergence == 0, notes = notes),
       object = list(model = "count", method = method, correction = corr, mu = mu, s2 = s2,
                     tau2 = if (!is.null(ss)) ss$tau2 else 0, xT = xT, PT = PT, threshold = Nq,
                     horizons = hs, tmax = tmax, cdf = cdf_tab, classic = stats::setNames(classic, hs),
                     boot = data.frame(mu = mub, s2 = s2b, xT = xTb), ss = ss, df = df, checks = chk,
                     nboot = length(mub), span = span))
}

# Checks of the diffusion assumptions: density dependence (slope of growth
# increment on log N), lag-1 autocorrelation of increments, normality.
count_checks <- function(y, tt, xs = NULL) {
  s <- diff(tt)
  inc <- diff(y) / s
  x0 <- if (is.null(xs)) y[-length(y)] else xs[-length(xs)]
  notes <- character()
  dd <- tryCatch(summary(stats::lm(inc ~ x0, weights = s))$coefficients[2, ], error = function(e) rep(NA, 4))
  z <- (diff(y) - s * sum(diff(y)) / sum(s)) / sqrt(s)       # standardised increments (equal variance under the model)
  ac <- if (length(z) > 5) stats::cor(z[-1], z[-length(z)]) else NA
  sw <- if (length(z) >= 8) stats::shapiro.test(z)$p.value else NA
  if (is.finite(dd[4]) && dd[4] < 0.05 && dd[1] < 0)
    notes <- c(notes, sprintf("Check: growth rate declines with abundance (slope %.3f per log-unit, p = %.3f), suggesting density dependence; density-independent PVA then tends to over-state long-term risk for populations near carrying capacity (Sabo et al. 2004). Note the slope test is biased towards finding density dependence when counts contain observation error.", dd[1], dd[4]))
  if (is.finite(ac) && abs(ac) > 0.4)
    notes <- c(notes, sprintf("Check: year-to-year growth rates are autocorrelated (lag-1 r = %.2f); negative values are expected with counting error (use method: state_space), positive values suggest age-structure cycles or trends in the environment.", ac))
  if (is.finite(sw) && sw < 0.01)
    notes <- c(notes, sprintf("Check: log growth increments are not normal (Shapiro-Wilk p = %.3f); look for outlier years (catastrophes, census changes).", sw))
  list(dd_slope = unname(dd[1]), dd_p = unname(dd[4]), lag1 = ac, shapiro_p = sw, notes = notes,
       inc = data.frame(year = tt[-1], increment = inc, logN = x0))
}

# ---- read / check --------------------------------------------------------

count_read_data <- function(config, dc) {
  f <- pva_path(config, dc$file)
  if (!file.exists(f)) stop("Count file not found: ", f, call. = FALSE)
  df <- utils::read.csv(f, check.names = FALSE, strip.white = TRUE)
  yc <- dc$year_col %||% "year"; cc <- dc$count_col %||% "count"
  for (cn in c(yc, cc)) if (!cn %in% names(df))
    stop(sprintf("Column '%s' not in %s. Columns present: %s. Set data: year_col / count_col in the config.",
                 cn, dc$file, paste(names(df), collapse = ", ")), call. = FALSE)
  out <- data.frame(year = suppressWarnings(as.numeric(df[[yc]])), count = suppressWarnings(as.numeric(df[[cc]])))
  raw <- trimws(as.character(df[[cc]]))
  bad <- which(is.na(out$year) | (is.na(out$count) & !is.na(raw) & nzchar(raw) & !toupper(raw) %in% c("NA", "NAN")))
  if (length(bad)) stop(sprintf("Non-numeric year or count in data row(s) %s of %s (%s). Fix the typo; leave the count cell blank (or NA) for a year that was not surveyed.",
                                paste(head(bad, 10), collapse = ", "), dc$file,
                                paste(head(sprintf("%s = '%s', %s = '%s'", yc, trimws(as.character(df[[yc]][bad])), cc, raw[bad]), 10), collapse = "; ")), call. = FALSE)
  out <- out[!is.na(out$count), ]   # empty cells = year not surveyed
  st <- config$settings %||% list()
  list(model = "count", df = out, threshold = st$threshold, horizons = st$horizons %||% c(10, 25, 50),
       meta = list(species = dc$species %||% "the population", units = dc$units %||% "individuals",
                   file = dc$file))
}

count_check_data <- function(data) {
  pr <- character(); df <- data$df
  if (nrow(df) < 5) pr <- c(pr, sprintf("Only %d counts: a count-based PVA needs at least 5 (and realistically 15+) years.", nrow(df)))
  if (any(duplicated(df$year))) pr <- c(pr, sprintf("Year(s) %s appear more than once; give one count per year (sum or average sub-samples first).",
                                                   paste(unique(df$year[duplicated(df$year)]), collapse = ", ")))
  if (any(df$count <= 0)) pr <- c(pr, sprintf("Zero or negative counts in year(s) %s: the log-scale model cannot use them. Remove those years (treat as missing) or, if the population was truly absent, count-based PVA is not appropriate.",
                                              paste(df$year[df$count <= 0], collapse = ", ")))
  if (is.null(data$threshold) || !is.numeric(data$threshold) || data$threshold <= 0)
    pr <- c(pr, "settings: threshold (quasi-extinction threshold, > 0, same units as the counts) is missing.")
  else if (nrow(df) && data$threshold >= df$count[which.max(df$year)])
    pr <- c(pr, sprintf("The quasi-extinction threshold (%g) is at or above the latest count (%g): the population is already quasi-extinct by this definition.",
                        data$threshold, df$count[which.max(df$year)]))
  if (nrow(df) > 1 && any(diff(sort(df$year)) <= 0)) pr <- c(pr, "Years must be distinct.")
  pr
}
