# multistate: from link-scale draws (Nimble names) to the biological
# quantities reported by the module. Works on a matrix of draws (rows) and is
# also used for ML point estimates / parametric draws.

ms_get <- function(m, nm) {
  if (nm %in% colnames(m)) m[, nm] else rep(0, nrow(m))
}

ms_softmax_rows <- function(E) { E <- exp(E - apply(E, 1, max)); E / rowSums(E) }

# Transition matrix (conditional on survival) for each draw at covariate
# values xs (standardised, length Kpsi): returns array draws x S x S.
ms_psi_draws <- function(m, des, xs = rep(0, des$Kpsi)) {
  S <- des$S; n <- nrow(m)
  out <- array(0, c(n, S, S))
  for (r in 1:S) {
    E <- matrix(0, n, S)
    for (s in 1:S) if (s != r) {
      E[, s] <- ms_get(m, sprintf("apsi[%d, %d]", r, s))
      for (k in seq_len(des$Kpsi)) E[, s] <- E[, s] + ms_get(m, sprintf("bpsi[%d, %d, %d]", r, s, k)) * xs[k]
    }
    out[, r, ] <- ms_softmax_rows(E)
  }
  out
}

ms_phi_draws <- function(m, des, xs = rep(0, des$Kphi)) {
  S <- des$S
  sapply(1:S, function(s) {
    lp <- ms_get(m, sprintf("a_phi[%d]", s))
    for (k in seq_len(des$Kphi)) lp <- lp + ms_get(m, sprintf("bphi[%d]", k)) * xs[k]
    stats::plogis(lp)
  }) |> matrix(nrow = nrow(m))
}

ms_p_draws <- function(m, des, xs = rep(0, des$Kp)) {
  sapply(1:des$S, function(s) {
    lp <- ms_get(m, sprintf("a_p[%d]", s))
    for (k in seq_len(des$Kp)) lp <- lp + ms_get(m, sprintf("bp[%d]", k)) * xs[k]
    stats::plogis(lp)
  }) |> matrix(nrow = nrow(m))
}

ms_stationary_draws <- function(P) {
  n <- dim(P)[1]; S <- dim(P)[2]
  t(vapply(seq_len(n), function(d) {
    v <- tryCatch(ms_stationary(P[d, , ]), error = function(e) rep(NA_real_, S))
    v
  }, numeric(S))) |> matrix(nrow = n)
}

ms_derive <- function(m, des) {
  m <- as.matrix(m); S <- des$S; st <- des$states; n <- nrow(m)
  out <- list()
  phi <- ms_phi_draws(m, des)
  for (s in 1:S) out[[paste0("survival_", st[s])]] <- phi[, s]
  P <- ms_psi_draws(m, des)
  for (j in seq_len(nrow(des$pairs))) {
    r <- des$pairs[j, 1]; s <- des$pairs[j, 2]
    out[[sprintf("transition_%s_to_%s", st[r], st[s])]] <- P[, r, s]
  }
  p <- ms_p_draws(m, des)
  for (s in 1:S) out[[paste0("detection_", st[s])]] <- p[, s]
  if (des$model == "multievent") {
    for (s in 1:S) out[[paste0("assignment_", st[s])]] <- stats::plogis(ms_get(m, sprintf("a_delta[%d]", s)))
    # vapply + explicit matrix: with a single row (ML point estimate) sapply would
    # return a vector and cbind would build the wrong matrix for S >= 3.
    A <- cbind(0, matrix(vapply(seq_len(S - 1) + 1, function(s) ms_get(m, sprintf("a_pi[%d]", s)), numeric(n)), n))
    pi0 <- ms_softmax_rows(A)
    for (s in 1:S) out[[paste0("initial_", st[s])]] <- pi0[, s]
  }
  stn <- ms_stationary_draws(P)
  for (s in 1:S) out[[paste0("stationary_", st[s])]] <- stn[, s]
  for (k in seq_len(des$Kphi)) out[[paste0("survival_beta_", des$covs$survival[k])]] <- ms_get(m, sprintf("bphi[%d]", k))
  for (k in seq_len(des$Kpsi)) for (j in seq_len(nrow(des$pairs))) {
    r <- des$pairs[j, 1]; s <- des$pairs[j, 2]
    out[[sprintf("transition_beta_%s_%s_to_%s", des$covs$transition[k], st[r], st[s])]] <-
      ms_get(m, sprintf("bpsi[%d, %d, %d]", r, s, k))
  }
  for (k in seq_len(des$Kp)) out[[paste0("detection_beta_", des$covs$detection[k])]] <- ms_get(m, sprintf("bp[%d]", k))
  do.call(cbind, out)
}
