# multistate: a hand-written forward-algorithm likelihood in plain R,
# written independently of the Nimble code (different formulation: the joint
# probability of the first event and the state at first capture is the
# starting forward vector, then alpha_{t+1} = (alpha_t G_t) * P(y_{t+1} | state)).
# Used (1) as the maximum-likelihood cross-check of the Nimble engine and
# (2) as the fast "ml" engine (optim + Hessian, Wald intervals on the link
# scale, parametric draws for derived quantities).
#
# theta layout (ms_theta_index): a_phi[S], apsi[M] (pair order), bphi[Kphi],
# bpsi[M x Kpsi] (pair fastest), a_p[S], bp[Kp], a_delta[S] (ME), a_pi[2..S] (ME)

ms_theta_index <- function(des) {
  S <- des$S; M <- nrow(des$pairs); ME <- des$model == "multievent"
  len <- c(a_phi = S, apsi = M, bphi = des$Kphi, bpsi = M * des$Kpsi, a_p = S, bp = des$Kp,
           a_delta = if (ME) S else 0, a_pi = if (ME) S - 1 else 0)
  end <- cumsum(len); start <- end - len + 1
  idx <- lapply(names(len), function(n) if (len[[n]] > 0) start[[n]]:end[[n]] else integer())
  names(idx) <- names(len)
  idx$n <- sum(len)
  idx
}

# Convert theta into the Nimble parameter names (a_phi[1], apsi[1, 2], ...).
ms_theta_to_nimble <- function(theta, des) {
  ix <- ms_theta_index(des); S <- des$S; pr <- des$pairs
  out <- c(stats::setNames(theta[ix$a_phi], sprintf("a_phi[%d]", 1:S)),
           stats::setNames(theta[ix$a_p], sprintf("a_p[%d]", 1:S)),
           stats::setNames(theta[ix$apsi], sprintf("apsi[%d, %d]", pr[, 1], pr[, 2])))
  if (des$Kphi > 0) out <- c(out, stats::setNames(theta[ix$bphi], sprintf("bphi[%d]", seq_len(des$Kphi))))
  if (des$Kpsi > 0) {
    b <- matrix(theta[ix$bpsi], nrow(pr))
    for (k in seq_len(des$Kpsi)) out <- c(out, stats::setNames(b[, k], sprintf("bpsi[%d, %d, %d]", pr[, 1], pr[, 2], k)))
  }
  if (des$Kp > 0) out <- c(out, stats::setNames(theta[ix$bp], sprintf("bp[%d]", seq_len(des$Kp))))
  if (des$model == "multievent")
    out <- c(out, stats::setNames(theta[ix$a_delta], sprintf("a_delta[%d]", 1:S)),
             stats::setNames(theta[ix$a_pi], sprintf("a_pi[%d]", 2:S)))
  out
}

ms_negll <- function(theta, des) {
  S <- des$S; S1 <- S + 1; T <- des$T; N <- des$N; pr <- des$pairs; M <- nrow(pr)
  ix <- ms_theta_index(des); ME <- des$model == "multievent"
  a_phi <- theta[ix$a_phi]; apsi <- theta[ix$apsi]
  bphi <- theta[ix$bphi]; bpsi <- matrix(theta[ix$bpsi], M, des$Kpsi)
  a_p <- theta[ix$a_p]; bp <- theta[ix$bp]
  delta <- if (ME) stats::plogis(theta[ix$a_delta]) else rep(1, S)
  pi0 <- if (ME) { e <- exp(c(0, theta[ix$a_pi])); e / sum(e) } else rep(1 / S, S)
  # covariate rows per individual and interval
  Xrow <- function(X, t, K) {
    if (K == 0) return(matrix(0, N, 0))
    if (des$gam_by_ind) matrix(X[, t, seq_len(K)], N, K) else matrix(X[t, seq_len(K)], N, K, byrow = TRUE)
  }
  # starting forward vector at first capture: P(state s, first event)
  alpha <- matrix(0, N, S1)
  for (s in 1:S) {
    pe <- ifelse(des$e1 == s, delta[s], ifelse(des$e1 == S + 1, 1 - delta[s], 0))
    alpha[, s] <- if (ME) pi0[s] * pe else as.numeric(des$e1 == s)
  }
  ll <- log(rowSums(alpha)); alpha <- alpha / rowSums(alpha)
  for (t in 1:(T - 1)) {
    act <- des$f <= t                       # individuals already marked by t
    if (!any(act)) next
    Xf <- Xrow(des$Xphi, t, des$Kphi); Xs <- Xrow(des$Xpsi, t, des$Kpsi)
    phi <- stats::plogis(outer(rep(1, N), a_phi) + if (des$Kphi) drop(Xf %*% bphi) else 0)  # N x S
    eta <- array(0, c(N, S, S))
    for (m in 1:M) eta[, pr[m, 1], pr[m, 2]] <- apsi[m] + if (des$Kpsi) drop(Xs %*% bpsi[m, ]) else 0
    new <- matrix(0, N, S1)
    for (r in 1:S) {
      e <- exp(eta[, r, , drop = FALSE]); dim(e) <- c(N, S)
      psi <- e / rowSums(e)
      new[, 1:S] <- new[, 1:S] + alpha[, r] * phi[, r] * psi
      new[, S1] <- new[, S1] + alpha[, r] * (1 - phi[, r])
    }
    new[, S1] <- new[, S1] + alpha[, S1]
    # observation at t+1
    lp <- a_p + if (des$Kp) sum(des$Xp[t + 1, seq_len(des$Kp)] * bp) else 0
    p <- stats::plogis(lp)
    yt <- des$y[, t + 1]
    obs <- matrix(0, N, S1)
    for (s in 1:S) obs[, s] <- ifelse(yt == 0, 1 - p[s], ifelse(yt == s, p[s] * delta[s],
                                     ifelse(yt == S + 1, p[s] * (1 - delta[s]), 0)))
    obs[, S1] <- as.numeric(yt == 0)
    new <- new * obs
    sc <- rowSums(new)
    ll[act] <- ll[act] + log(sc[act])
    alpha[act, ] <- new[act, ] / sc[act]
  }
  -sum(ll)
}

ms_fit_ml <- function(des, settings = list(), start = NULL) {
  ix <- ms_theta_index(des)
  if (is.null(start)) {
    start <- numeric(ix$n); start[ix$a_phi] <- 1; start[ix$apsi] <- -1.5
    if (length(ix$a_delta)) start[ix$a_delta] <- 1
  }
  fn <- function(th) { v <- ms_negll(th, des); if (!is.finite(v)) 1e10 else v }
  o <- stats::optim(start, fn, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  o2 <- stats::optim(o$par, fn, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
  if (o2$value < o$value - 1e-8) o <- stats::optim(o2$par, fn, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
  H <- tryCatch(stats::optimHess(o$par, fn), error = function(e) NULL)
  V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
  list(par = o$par, value = o$value, convergence = o$convergence, vcov = V,
       nimble_par = ms_theta_to_nimble(o$par, des))
}
