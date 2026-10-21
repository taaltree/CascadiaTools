# multistate: simulation of multistate / multievent capture-recapture data.
#
# Parameters are flat, named on the biological scale with the state names, so
# that power analyses can vary any one of them:
#   survival_<s>                 annual (per-interval) survival in state s at mean covariates
#   transition_<r>_to_<s>        P(move r -> s | survives), r != s, at mean covariates
#   detection_<s>                P(detected | alive in s) at mean covariates
#   assignment_<s>               P(state recorded | detected in s)  (multievent only)
#   initial_<s>                  P(in state s at first capture)     (multievent only)
#   survival_beta_<cov>          logit-scale slope (per SD) shared by all states
#   transition_beta_<cov>_<r>_to_<s>   multinomial-logit slope (per SD) of r -> s
#   detection_beta_<cov>         logit-scale slope (per SD) shared by all states
# Missing values fall back to sensible defaults so a 3-state design only needs
# the parameters the user cares about.
#
# Design keys: states (names), n_occasions, n_marked (new individuals marked at
# each occasion 1..T-1), model ("multievent" | "AS"), survival_covs
# (individual covariates), transition_covs (occasion covariates),
# detection_covs (occasion covariates).

ms_param <- function(params, name, default) {
  v <- params[[name]]
  if (is.null(v)) default else as.numeric(v)
}

ms_true_matrices <- function(params, states) {
  S <- length(states)
  phi <- vapply(states, function(s) ms_param(params, paste0("survival_", s), 0.7), 0)
  p <- vapply(states, function(s) ms_param(params, paste0("detection_", s), 0.5), 0)
  delta <- vapply(states, function(s) ms_param(params, paste0("assignment_", s), 1), 0)
  pi0 <- vapply(states, function(s) ms_param(params, paste0("initial_", s), NA_real_), 0)
  if (anyNA(pi0)) {           # unspecified: the rest shared equally
    rest <- max(0, 1 - sum(pi0, na.rm = TRUE)); pi0[is.na(pi0)] <- rest / sum(is.na(pi0))
  }
  pi0 <- pi0 / sum(pi0)
  psi <- matrix(0, S, S, dimnames = list(states, states))
  for (r in seq_len(S)) {
    for (s in seq_len(S)) if (r != s)
      psi[r, s] <- ms_param(params, sprintf("transition_%s_to_%s", states[r], states[s]), 0.2 / (S - 1))
    if (sum(psi[r, ]) >= 1) stop(sprintf("transition probabilities out of state %s sum to >= 1", states[r]))
    psi[r, r] <- 1 - sum(psi[r, ])
  }
  list(phi = phi, p = p, delta = delta, pi = pi0, psi = psi)
}

# Stationary distribution of a row-stochastic matrix (left eigenvector).
ms_stationary <- function(P) {
  S <- nrow(P)
  A <- rbind(t(P) - diag(S), rep(1, S))
  b <- c(rep(0, S), 1)
  as.vector(qr.solve(A, b))
}

ms_simulate <- function(params, design, seed) {
  set.seed(seed)
  states <- as.character(unlist(design$states))
  S <- length(states); T <- design$n_occasions
  model <- if (is.null(design$model)) "multievent" else design$model
  tm <- ms_true_matrices(params, states)
  if (model == "AS") tm$delta[] <- 1
  nm <- design$n_marked
  first <- rep(seq_len(T - 1), each = nm)
  N <- length(first)
  zs <- function(x) as.vector(scale(x))
  ind_covs <- data.frame(id = sprintf("ind%04d", seq_len(N)))
  occ_covs <- data.frame(occasion = seq_len(T))
  scov <- as.character(unlist(design$survival_covs)); tcov <- as.character(unlist(design$transition_covs))
  dcov <- as.character(unlist(design$detection_covs))
  for (k in scov) ind_covs[[k]] <- zs(stats::rnorm(N))
  for (k in unique(c(tcov, dcov))) occ_covs[[k]] <- zs(stats::rnorm(T))
  # linear predictors
  lphi <- matrix(stats::qlogis(tm$phi), N, S, byrow = TRUE)
  ms_lp_add <- function(lp, k, beta) lp + beta * ind_covs[[k]]
  for (k in scov) lphi <- ms_lp_add(lphi, k, ms_param(params, paste0("survival_beta_", k), 0))
  phi_i <- stats::plogis(lphi)                          # N x S (individual covariates only)
  eta0 <- log(tm$psi / diag(tm$psi))                    # row r: log(psi_rs / psi_rr)
  psi_t <- array(0, c(S, S, T - 1))
  for (t in seq_len(T - 1)) {
    eta <- eta0
    for (k in tcov) for (r in seq_len(S)) for (s in seq_len(S)) if (r != s)
      eta[r, s] <- eta[r, s] + ms_param(params, sprintf("transition_beta_%s_%s_to_%s", k, states[r], states[s]), 0) * occ_covs[[k]][t]
    e <- exp(eta); psi_t[, , t] <- e / rowSums(e)
  }
  lp <- matrix(stats::qlogis(tm$p), T, S, byrow = TRUE)
  for (k in dcov) lp <- lp + ms_param(params, paste0("detection_beta_", k), 0) * occ_covs[[k]]
  p_t <- stats::plogis(lp)                              # T x S
  y <- matrix(0L, N, T)
  z <- matrix(NA_integer_, N, T)                        # true state, S+1 = dead
  for (i in seq_len(N)) {
    f <- first[i]
    z[i, f] <- sample.int(S, 1, prob = tm$pi)
    y[i, f] <- if (stats::runif(1) < tm$delta[z[i, f]]) z[i, f] else S + 1L
    if (f < T) for (t in (f + 1):T) {
      zp <- z[i, t - 1]
      if (zp == S + 1L || stats::runif(1) > phi_i[i, zp]) { z[i, t] <- S + 1L; next }
      z[i, t] <- sample.int(S, 1, prob = psi_t[zp, , t - 1])
      if (stats::runif(1) < p_t[t, z[i, t]])
        y[i, t] <- if (stats::runif(1) < tm$delta[z[i, t]]) z[i, t] else S + 1L
    }
  }
  data <- list(y = y, states = states, ids = ind_covs$id, occasion_labels = as.character(seq_len(T)),
               ind_covs = ind_covs, occ_covs = occ_covs, model = model,
               formulas = list(survival = scov, transition = tcov, detection = dcov))
  truth <- ms_truth(params, states, tm, model, scov, tcov, dcov)
  list(data = data, truth = truth)
}

ms_truth <- function(params, states, tm, model, scov, tcov, dcov) {
  S <- length(states)
  tr <- c(stats::setNames(tm$phi, paste0("survival_", states)))
  for (r in seq_len(S)) for (s in seq_len(S)) if (r != s)
    tr[sprintf("transition_%s_to_%s", states[r], states[s])] <- tm$psi[r, s]
  tr <- c(tr, stats::setNames(tm$p, paste0("detection_", states)))
  if (model == "multievent")
    tr <- c(tr, stats::setNames(tm$delta, paste0("assignment_", states)),
            stats::setNames(tm$pi, paste0("initial_", states)))
  tr <- c(tr, stats::setNames(ms_stationary(tm$psi), paste0("stationary_", states)))
  for (k in scov) tr[paste0("survival_beta_", k)] <- ms_param(params, paste0("survival_beta_", k), 0)
  for (k in tcov) for (r in seq_len(S)) for (s in seq_len(S)) if (r != s) {
    nm <- sprintf("transition_beta_%s_%s_to_%s", k, states[r], states[s])
    tr[nm] <- ms_param(params, nm, 0)
  }
  for (k in dcov) tr[paste0("detection_beta_", k)] <- ms_param(params, paste0("detection_beta_", k), 0)
  tr
}
