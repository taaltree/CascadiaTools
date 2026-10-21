# multistate: goodness of fit.
#
# 1. Posterior (or parametric, for ml) predictive check. Replicate capture
#    histories are simulated from each draw, conditional on each animal's
#    first capture (occasion and first event, as in the likelihood). Three
#    discrepancies:
#      - Freeman-Tukey distance between observed and expected numbers of each
#        "seen" event (seen in state s / state unknown) at each occasion;
#      - number of recorded state switches between consecutive sightings with a
#        known state (sensitive to misclassification and to non-Markovian
#        movement);
#      - number of animals never seen again after first capture (sensitive to
#        transients / trap-shyness).
#    Bayesian p-value = P(replicate discrepancy >= observed); values near 0
#    or 1 flag misfit.
# 2. U-CARE tests for the Arnason-Schwarz model (Pradel et al. 2003; R2ucare
#    overall_JMV), on the histories with any "state unknown" sighting
#    treated as not seen.

# Per-draw model matrices at interval t (1..T-1): phi N x S, psi N x S x S.
ms_mats_t <- function(m1, des, t) {
  S <- des$S; N <- des$N
  xrow <- function(X, K) {
    if (K == 0) return(matrix(0, N, 0))
    if (des$gam_by_ind) matrix(X[, t, seq_len(K)], N, K) else matrix(X[t, seq_len(K)], N, K, byrow = TRUE)
  }
  Xf <- xrow(des$Xphi, des$Kphi); Xs <- xrow(des$Xpsi, des$Kpsi)
  phi <- matrix(0, N, S); psi <- array(0, c(N, S, S))
  for (s in 1:S) {
    lp <- m1[[sprintf("a_phi[%d]", s)]] + 0
    for (k in seq_len(des$Kphi)) lp <- lp + m1[[sprintf("bphi[%d]", k)]] * Xf[, k]
    phi[, s] <- stats::plogis(lp)
  }
  for (r in 1:S) {
    E <- matrix(0, N, S)
    for (s in 1:S) if (s != r) {
      E[, s] <- m1[[sprintf("apsi[%d, %d]", r, s)]]
      for (k in seq_len(des$Kpsi)) E[, s] <- E[, s] + m1[[sprintf("bpsi[%d, %d, %d]", r, s, k)]] * Xs[, k]
    }
    psi[, r, ] <- ms_softmax_rows(E)
  }
  list(phi = phi, psi = psi)
}

ms_obs_t <- function(m1, des, t) {   # detection p_s at occasion t, delta_s
  S <- des$S
  p <- vapply(1:S, function(s) {
    lp <- m1[[sprintf("a_p[%d]", s)]]
    for (k in seq_len(des$Kp)) lp <- lp + m1[[sprintf("bp[%d]", k)]] * des$Xp[t, k]
    stats::plogis(lp)
  }, 0)
  delta <- if (des$model == "multievent") stats::plogis(vapply(1:S, function(s) m1[[sprintf("a_delta[%d]", s)]], 0)) else rep(1, S)
  list(p = p, delta = delta)
}

ms_post_first <- function(m1, des) {     # N x S: P(state at first capture | first event)
  S <- des$S
  post <- matrix(0, des$N, S)
  known <- des$e1 <= S
  post[cbind(which(known), des$e1[known])] <- 1
  if (any(!known)) {
    A <- c(0, vapply(seq_len(S - 1) + 1, function(s) m1[[sprintf("a_pi[%d]", s)]], 0))
    pi0 <- exp(A) / sum(exp(A))
    delta <- ms_obs_t(m1, des, 1)$delta
    w <- pi0 * (1 - delta); w <- w / sum(w)
    post[!known, ] <- matrix(w, sum(!known), S, byrow = TRUE)
  }
  post
}

ms_stat_switch <- function(y, S) {
  sw <- 0
  for (i in seq_len(nrow(y))) {
    k <- y[i, y[i, ] >= 1 & y[i, ] <= S]
    if (length(k) > 1) sw <- sw + sum(diff(k) != 0)
  }
  sw
}

ms_stat_never <- function(y, f) sum(vapply(seq_len(nrow(y)), function(i) sum(y[i, -seq_len(f[i])] > 0) == 0, TRUE))

ms_ppc <- function(draws, des, nd = 300) {
  draws <- as.matrix(draws)
  nd <- min(nd, nrow(draws))
  idx <- unique(round(seq(1, nrow(draws), length.out = nd)))
  S <- des$S; S1 <- S + 1; T <- des$T; N <- des$N
  O <- S + as.integer(des$model == "multievent")        # "seen" event codes 1..O
  obs_tab <- sapply(2:T, function(t) tabulate(des$y[, t][des$y[, t] > 0], O))
  obs_tab <- matrix(obs_tab, O)
  sw_obs <- ms_stat_switch(des$y, S); nv_obs <- ms_stat_never(des$y, des$f)
  ft_o <- ft_r <- sw_r <- nv_r <- numeric(length(idx))
  for (j in seq_along(idx)) {
    m1 <- as.list(draws[idx[j], ])
    post <- ms_post_first(m1, des)
    # expected distribution (N x S1) and a replicate
    dist <- cbind(post, 0) * 0
    z <- rep(NA_integer_, N)
    yr <- matrix(0L, N, T)
    for (i in seq_len(N)) yr[i, des$f[i]] <- des$e1[i]
    exp_tab <- matrix(0, O, T - 1); rep_tab <- matrix(0, O, T - 1)
    for (t in 1:(T - 1)) {
      start <- des$f == t
      if (any(start)) {
        dist[start, ] <- cbind(post[start, , drop = FALSE], 0)
        z[start] <- apply(post[start, , drop = FALSE], 1, function(pr) sample.int(S, 1, prob = pr))
      }
      act <- des$f <= t
      mt <- ms_mats_t(m1, des, t)
      nd2 <- matrix(0, N, S1)
      for (r in 1:S) {
        nd2[, 1:S] <- nd2[, 1:S] + dist[, r] * mt$phi[, r] * mt$psi[, r, ]
        nd2[, S1] <- nd2[, S1] + dist[, r] * (1 - mt$phi[, r])
      }
      nd2[, S1] <- nd2[, S1] + dist[, S1]
      dist[act, ] <- nd2[act, ]
      ob <- ms_obs_t(m1, des, t + 1)
      for (s in 1:S) {
        exp_tab[s, t] <- sum(dist[act, s]) * ob$p[s] * ob$delta[s]
        if (O > S) exp_tab[O, t] <- exp_tab[O, t] + sum(dist[act, s]) * ob$p[s] * (1 - ob$delta[s])
      }
      # replicate step
      for (i in which(act & z <= S)) {
        zi <- z[i]
        if (stats::runif(1) > mt$phi[i, zi]) { z[i] <- S1; next }
        z[i] <- sample.int(S, 1, prob = mt$psi[i, zi, ])
        if (stats::runif(1) < ob$p[z[i]]) yr[i, t + 1] <- if (stats::runif(1) < ob$delta[z[i]]) z[i] else S + 1L
      }
      rep_tab[, t] <- tabulate(yr[, t + 1][yr[, t + 1] > 0], O)
    }
    ft_o[j] <- sum((sqrt(obs_tab) - sqrt(exp_tab))^2)
    ft_r[j] <- sum((sqrt(rep_tab) - sqrt(exp_tab))^2)
    sw_r[j] <- ms_stat_switch(yr, S)
    nv_r[j] <- ms_stat_never(yr, des$f)
  }
  list(ft_obs = ft_o, ft_rep = ft_r, p_ft = mean(ft_r >= ft_o),
       switch_obs = sw_obs, switch_rep = sw_r, p_switch = mean(sw_r >= sw_obs),
       never_obs = nv_obs, never_rep = nv_r, p_never = mean(nv_r >= nv_obs), n_draws = length(idx))
}

ms_ppc_note <- function(g) sprintf(paste(
  "Posterior predictive check (%d draws): Freeman-Tukey discrepancy of sightings per occasion and event, Bayesian p = %.2f;",
  "recorded state switches observed %d vs replicated median %.0f (p = %.2f); animals never seen again after first capture",
  "observed %d vs replicated median %.0f (p = %.2f). p-values below 0.05 or above 0.95 indicate misfit."),
  g$n_draws, g$p_ft, g$switch_obs, stats::median(g$switch_rep), g$p_switch, g$never_obs,
  stats::median(g$never_rep), g$p_never)

ms_ucare <- function(des) {
  eco_require("R2ucare")
  y <- des$y; y[y > des$S] <- 0L
  y <- y[rowSums(y > 0) > 0, , drop = FALSE]
  if (des$S < 2) return(NULL)
  ov <- suppressWarnings(R2ucare::overall_JMV(y, rep(1, nrow(y))))
  comps <- list(
    test3Gsr = tryCatch(suppressWarnings(R2ucare::test3Gsr(y, rep(1, nrow(y)))$test3Gsr), error = function(e) NULL),
    test3Gwbwa = tryCatch(suppressWarnings(R2ucare::test3Gwbwa(y, rep(1, nrow(y)))$test3Gwbwa), error = function(e) NULL),
    testM = tryCatch(suppressWarnings(R2ucare::testMitec(y, rep(1, nrow(y)))$testMitec), error = function(e) NULL))
  list(overall = ov, components = comps)
}

ms_ucare_note <- function(u, des) {
  ov <- u$overall
  v <- suppressWarnings(as.numeric(unlist(ov)))
  nm <- names(unlist(ov))
  chi <- v[grep("chi2", nm)][1]; df <- v[grep("degree", nm)][1]; p <- v[grep("p_val", nm)][1]
  sprintf("U-CARE overall goodness-of-fit test of the Arnason-Schwarz (JMV) model%s: chi-square %.1f on %.0f df, p = %.3f%s.",
          if (des$model == "multievent") " on the histories with 'state unknown' sightings treated as not seen" else "",
          chi, df, p, if (is.finite(p) && p < 0.05) " (lack of fit: consider transience, trap-dependence or memory effects; c-hat = chi-square/df could inflate the intervals)" else "")
}
