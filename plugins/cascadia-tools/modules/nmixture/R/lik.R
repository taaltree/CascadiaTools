# nmixture: hand-written likelihoods (independent of unmarked), posterior of
# site abundance N given the counts, and forward-filtering backward-sampling
# for the open (Dail-Madsen) model. Used for (1) derived abundance with
# parameter uncertainty propagated, (2) the "rlik" engine that cross-checks
# unmarked, (3) posterior predictive checks.

nm_logsumexp_rows <- function(L) {
  mx <- apply(L, 1, max)
  mx[!is.finite(mx)] <- 0
  mx + log(rowSums(exp(L - mx)))
}

# log P(N = k), k = 0..K, for each site: M x (K+1)
nm_log_prior <- function(lam, K, mixture, size = NA, zi = 0) {
  M <- length(lam)
  kk <- matrix(0:K, M, K + 1, byrow = TRUE)
  lm <- matrix(lam, M, K + 1)
  if (mixture == "P") return(stats::dpois(kk, lm, log = TRUE))
  if (mixture == "NB") return(stats::dnbinom(kk, size = size, mu = lm, log = TRUE))
  if (mixture == "ZIP") {
    lp <- log(1 - zi) + stats::dpois(kk, lm, log = TRUE)
    lp[, 1] <- log(zi + (1 - zi) * exp(-lam))
    return(lp)
  }
  stop("unknown mixture ", mixture)
}

# log P(y_i. | N = k) over visits: M x (K+1). y, p: M x J (NA / p = 0 for visits not done)
nm_log_emis <- function(y, p, K) {
  M <- nrow(y)
  kk <- matrix(0:K, M, K + 1, byrow = TRUE)
  L <- matrix(0, M, K + 1)
  for (j in seq_len(ncol(y))) {
    ok <- !is.na(y[, j])
    if (!any(ok)) next
    l <- stats::dbinom(matrix(y[ok, j], sum(ok), K + 1), kk[ok, , drop = FALSE],
                       matrix(p[ok, j], sum(ok), K + 1), log = TRUE)
    L[ok, ] <- L[ok, ] + l
  }
  L
}

# Static model: per-site log-likelihood and posterior P(N_i = k | y_i)
nm_static_post <- function(u, des, K, mixture) {
  y <- matrix(des$y, des$M, des$J); p <- matrix(u$p, des$M, des$J)
  L <- nm_log_prior(u$lam, K, mixture, u$size, u$zi) + nm_log_emis(y, p, K)
  ll <- nm_logsumexp_rows(L)
  list(loglik = ll, post = exp(L - ll))
}

# Sample one value per row from a probability matrix (values 0..K)
nm_sample_rows <- function(P) {
  cp <- P %*% upper.tri(diag(ncol(P)), diag = TRUE)
  u <- stats::runif(nrow(P)) * cp[, ncol(P)]
  rowSums(cp < u)
}

# Transition matrix for constant/notrend Dail-Madsen dynamics:
# P(N_t = n | N_{t-1} = m) = sum_s Bin(s; m, omega) Pois(n - s; gamma), truncated at K
nm_trans <- function(gam, om, K) {
  B <- outer(0:K, 0:K, function(m, s) stats::dbinom(s, m, om))
  P <- outer(0:K, 0:K, function(s, n) stats::dpois(n - s, gam))
  B %*% P
}

# Forward algorithm (and optional backward sampling) for the open model.
# Returns loglik per site and, if sample = TRUE, an M x T matrix of N draws.
nm_dyn_ffbs <- function(u, des, K, mixture, sample = FALSE) {
  M <- des$M; T <- des$T
  F <- vector("list", T)
  lp0 <- nm_log_prior(u$lam, K, mixture, u$size, u$zi)
  ll <- numeric(M)
  keyfun <- function(t) paste(signif(u$gam[, t], 12), signif(u$om[, t], 12))
  cache <- new.env()
  getT <- function(g, o) {
    k <- paste(signif(g, 12), signif(o, 12))
    if (is.null(cache[[k]])) cache[[k]] <- nm_trans(g, o, K)
    cache[[k]]
  }
  for (t in seq_len(T)) {
    em <- nm_log_emis(matrix(des$y[, , t], M, des$J), matrix(u$p[, , t], M, des$J), K)
    if (t == 1) {
      L <- lp0 + em
      mx <- apply(L, 1, max); f <- exp(L - mx)
    } else {
      keys <- keyfun(t - 1)
      pred <- matrix(0, M, K + 1)
      for (k in unique(keys)) {
        ii <- which(keys == k)
        pred[ii, ] <- F[[t - 1]][ii, , drop = FALSE] %*% getT(u$gam[ii[1], t - 1], u$om[ii[1], t - 1])
      }
      mx <- apply(em, 1, max)
      f <- pred * exp(em - mx)
    }
    cs <- rowSums(f)
    ll <- ll + log(cs) + mx
    F[[t]] <- f / cs
  }
  out <- list(loglik = ll, FT = F[[T]])
  if (sample) {
    N <- matrix(0L, M, T)
    N[, T] <- nm_sample_rows(F[[T]])
    if (T > 1) for (t in (T - 1):1) {
      keys <- keyfun(t)
      W <- matrix(0, M, K + 1)
      for (k in unique(keys)) {
        ii <- which(keys == k)
        Tm <- getT(u$gam[ii[1], t], u$om[ii[1], t])
        W[ii, ] <- F[[t]][ii, , drop = FALSE] * t(Tm[, N[ii, t + 1] + 1, drop = FALSE])
      }
      W <- W / rowSums(W)
      N[, t] <- nm_sample_rows(W)
    }
    out$N <- N
  }
  out
}

# Total log-likelihood for a theta vector (used by the rlik engine)
nm_loglik <- function(th, des, K, mixture, dynamics = "constant") {
  u <- nm_unpack(th, des, mixture, dynamics)
  if (des$model == "static") sum(nm_static_post(u, des, K, mixture)$loglik)
  else sum(nm_dyn_ffbs(u, des, K, mixture)$loglik)
}

# Draws of site abundance N (n_draws x M [x T]) given parameter draws.
nm_N_draws <- function(TH, des, K, mixture, dynamics = "constant") {
  nd <- nrow(TH)
  out <- array(0L, c(nd, des$M, des$T))
  for (d in seq_len(nd)) {
    th <- TH[d, ]
    u <- nm_unpack(th, des, mixture, dynamics)
    if (des$model == "static") {
      pp <- nm_static_post(u, des, K, mixture)$post
      out[d, , 1] <- nm_sample_rows(pp)
    } else out[d, , ] <- nm_dyn_ffbs(u, des, K, mixture, sample = TRUE)$N
  }
  out
}

# Upper-tail probability of the abundance distribution beyond K at the
# largest fitted mean (used to choose / check K).
nm_tail_prob <- function(u, K, mixture, des) {
  mu <- max(u$lam)
  # Dynamic: largest expected N over the observed seasons, E[N_t] = E[N_t-1] omega + gamma
  # (exact Poisson for a Poisson start). The earlier bound, the stationary mean gamma/(1 - omega)
  # at the largest omega, can never be reached in T seasons and inflated K (and run time ~K^3)
  # badly when omega has covariates (verification 2026-10-20).
  if (des$model == "dynamic" && des$T > 1) {
    EN <- u$lam
    for (t in seq_len(des$T - 1)) { EN <- EN * u$om[, t] + u$gam[, t]; mu <- max(mu, max(EN)) }
  }
  if (mixture == "NB") return(stats::pnbinom(K, size = u$size, mu = mu, lower.tail = FALSE))
  stats::ppois(K, mu, lower.tail = FALSE)
}
