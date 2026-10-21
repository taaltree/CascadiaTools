# scr: turn engine draws (internal units) into the reported quantities, plus
# the per-draw pieces shared by figures and checks. These functions mirror
# the compiled likelihood dSCRgen (nimble_engines.R) in plain R; scr_loglik()
# is used in dev tests to check the two agree.
#
# Draw columns: psi, det0 (p0 or lambda0), sigma, beta1..betaNB, phi [, N].

log1mexp <- function(lp) ifelse(lp > -0.693, log(-expm1(lp)), log1p(-exp(lp)))
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
scr_beta <- function(prep, th, from, len) if (len > 0) unname(th[paste0("beta", from + seq_len(len))]) else numeric(0)
scr_beta1 <- function(prep, th, idx) if (idx > 0) unname(th[[paste0("beta", idx)]]) else 0

# Per-draw linear predictors. th: named numeric vector (one draw).
scr_lpi <- function(prep, th) {
  PD <- ncol(prep$XD)
  if (!PD) return(rep(-log(prep$G), prep$G))
  lw <- as.vector(prep$XD %*% scr_beta(prep, th, 0, PD))
  lw - lse(lw)
}

# link-scale baseline detection eta (J*Kq x S), log sigma (J x S)
scr_lin <- function(prep, th, behav = 0) {
  b <- prep$bidx
  ld0 <- if (prep$obs == "poisson") log(th[["det0"]]) else stats::qlogis(th[["det0"]])
  eta <- rep(ld0, prep$J * prep$Kq)
  if (b[4] > 0) eta <- eta + as.vector(prep$XPq %*% scr_beta(prep, th, b[3], b[4]))
  if (behav) eta <- eta + scr_beta1(prep, th, b[9])
  ls <- rep(log(th[["sigma"]]), prep$J)
  if (b[6] > 0) ls <- ls + as.vector(prep$XS %*% scr_beta(prep, th, b[5], b[6]))
  E <- cbind(eta); L <- cbind(ls)
  if (prep$S == 2) { E <- cbind(eta, eta + scr_beta1(prep, th, b[7])); L <- cbind(ls, ls + scr_beta1(prep, th, b[8])) }
  list(eta = E, lsig = L)
}

# log p (or log lambda) for occasion column k, sex s: G x J
scr_logp_ks <- function(prep, lin, k, s) {
  e <- lin$eta[(k - 1) * prep$J + seq_len(prep$J), s]
  lq <- if (prep$obs == "poisson") e else stats::plogis(e, log.p = TRUE)
  a <- -0.5 * exp(-2 * lin$lsig[, s])
  sweep(sweep(prep$d2, 2, a, "*"), 2, lq, "+")
}

# base: G x S matrix of log pi_g + log P(not detected anywhere | g, s)
scr_base_all <- function(prep, th, lpi = scr_lpi(prep, th)) {
  lin <- scr_lin(prep, th)
  B <- matrix(lpi, prep$G, prep$S)
  for (s in seq_len(prep$S)) for (k in seq_len(prep$Kq)) {
    lp <- scr_logp_ks(prep, lin, k, s)
    u <- prep$use[, k]
    B[, s] <- B[, s] + if (prep$obs == "poisson") -as.vector(exp(lp) %*% u) else as.vector(log1mexp(lp) %*% u)
  }
  B
}
scr_sexw <- function(prep, th) if (prep$S == 2) c(th[["phi"]], 1 - th[["phi"]]) else 1

# log P(undetected) per animal (mixture over pixels and sex)
scr_lf0 <- function(prep, th, B = scr_base_all(prep, th)) {
  w <- scr_sexw(prep, th)
  lse(as.vector(sweep(B, 2, log(w), "+")))
}

# Per detected individual: n x G x S array of log pi_g + log P(history | g, s) (sex weight not included)
scr_ind_ll <- function(prep, th, B = scr_base_all(prep, th)) {
  n <- prep$n; G <- prep$G; S <- prep$S
  out <- array(0, c(n, G, S))
  lin <- scr_lin(prep, th); linb <- if (prep$behav == "trap") scr_lin(prep, th, behav = 1) else lin
  for (s in seq_len(S)) {
    V <- matrix(B[, s], n, G, byrow = TRUE)
    for (k in seq_len(prep$Kq)) {
      sel <- which(prep$yp[, k] > 0 | prep$bp[, k] > 0)
      if (!length(sel)) next
      lp0 <- scr_logp_ks(prep, lin, k, s)
      lpb <- if (prep$behav == "trap") scr_logp_ks(prep, linb, k, s) else lp0
      for (p in sel) {
        j <- prep$pj[p]; i <- prep$pi[p]; yy <- prep$yp[p, k]; u <- prep$use[j, k]
        l0 <- lp0[, j]; lb <- if (prep$bp[p, k] > 0) lpb[, j] else l0
        add <- switch(prep$obs_code,
          yy * (lb - log1mexp(lb)),
          if (yy > 0) log1mexp(u * log1mexp(lb)) - u * log1mexp(l0) else u * (log1mexp(lb) - log1mexp(l0)),
          yy * (log(u) + lb) - u * exp(lb) + u * exp(l0))
        V[i, ] <- V[i, ] + add
      }
    }
    out[, , s] <- V
  }
  out
}

# Full marginal log-likelihood (up to constants), same as dSCRgen
scr_loglik <- function(prep, th, M = prep$M) {
  B <- scr_base_all(prep, th); A <- scr_ind_ll(prep, th, B); w <- scr_sexw(prep, th)
  tot <- 0
  for (i in seq_len(prep$n)) {
    Li <- apply(A[i, , , drop = FALSE], 3, lse)
    if (prep$S == 1) tot <- tot + Li else {
      si <- prep$sexi[i]
      tot <- tot + if (si == 0) log(w[1]) + Li[1] else if (si == 1) log(w[2]) + Li[2] else lse(log(w) + Li)
    }
  }
  f0 <- exp(scr_lf0(prep, th, B))
  tot + prep$n * log(th[["psi"]]) + (M - prep$n) * log(th[["psi"]] * f0 + 1 - th[["psi"]])
}

scr_hr <- function(sig_internal, prep) pi * (sqrt(stats::qchisq(0.95, 2)) * sig_internal)^2 * scr_hr_conv(prep)

# Reported draws for one chain.
scr_report_chain <- function(prep, m, idx) {
  m <- m[idx, , drop = FALSE]
  S <- nrow(m)
  N <- if ("N" %in% colnames(m)) m[, "N"] else {
    vapply(seq_len(S), function(r) {
      th <- m[r, ]
      f0 <- exp(scr_lf0(prep, th))
      w <- th[["psi"]] * f0 / (th[["psi"]] * f0 + 1 - th[["psi"]])
      prep$n + stats::rbinom(1, prep$M - prep$n, w)
    }, 1)
  }
  Mx <- if (!is.null(prep$M_used)) prep$M_used else prep$M
  dn <- scr_det_name(prep)
  cols <- list()
  cols[[scr_density_name(prep$units)]] <- m[, "psi"] * Mx / prep$area * scr_area_conv(prep)
  cols$abundance <- N
  b <- prep$bidx
  if (prep$S == 1) {
    cols[[dn]] <- m[, "det0"]
    cols$sigma <- m[, "sigma"] * prep$cscale
    cols$home_range_area_95 <- scr_hr(m[, "sigma"], prep)
  } else {
    bsp <- if (b[7] > 0) m[, paste0("beta", b[7])] else 0
    bss <- if (b[8] > 0) m[, paste0("beta", b[8])] else 0
    dm <- if (prep$obs == "poisson") m[, "det0"] * exp(bsp) else stats::plogis(stats::qlogis(m[, "det0"]) + bsp)
    cols[[paste0(dn, "_female")]] <- m[, "det0"]
    cols[[paste0(dn, "_male")]] <- dm
    cols$sigma_female <- m[, "sigma"] * prep$cscale
    cols$sigma_male <- m[, "sigma"] * exp(bss) * prep$cscale
    cols$home_range_area_95_female <- scr_hr(m[, "sigma"], prep)
    cols$home_range_area_95_male <- scr_hr(m[, "sigma"] * exp(bss), prep)
    cols$proportion_female <- m[, "phi"]
  }
  for (k in seq_len(prep$nb)) cols[[prep$bnames[k]]] <- m[, paste0("beta", k)]
  out <- do.call(cbind, cols)
  if (!is.matrix(out)) out <- matrix(out, nrow = S, dimnames = list(NULL, names(cols)))
  colnames(out) <- names(cols)
  out
}

scr_thin_idx <- function(nr, max_draws) if (nr <= max_draws) seq_len(nr) else unique(round(seq(1, nr, length.out = max_draws)))

# Pool a thinned set of draws across chains (internal units)
scr_pool <- function(chains, ndraw) {
  all <- do.call(rbind, chains)
  all[scr_thin_idx(nrow(all), ndraw), , drop = FALSE]
}

# Posterior probability of each detected individual's activity centre over
# pixels (n x G), and the realised density surface (expected animals per
# pixel: detected individuals + expected undetected), averaged over draws.
scr_surfaces <- function(prep, draws) {
  n <- prep$n; G <- prep$G
  Pc <- matrix(0, n, G); real <- numeric(G); expd <- numeric(G)
  Mx <- if (!is.null(attr(draws, "M"))) attr(draws, "M") else if (!is.null(prep$M_used)) prep$M_used else prep$M
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    lpi <- scr_lpi(prep, th)
    B <- scr_base_all(prep, th, lpi); w <- scr_sexw(prep, th)
    if (n) {
      A <- scr_ind_ll(prep, th, B)
      ll <- matrix(-Inf, n, G)
      for (i in seq_len(n)) {
        cols <- lapply(seq_len(prep$S), function(s) {
          if (prep$S == 2 && prep$sexi[i] >= 0 && prep$sexi[i] != s - 1) return(rep(-Inf, G))
          A[i, , s] + log(w[s])
        })
        ll[i, ] <- if (prep$S == 1) cols[[1]] else pmax(cols[[1]], cols[[2]]) + log1p(exp(-abs(cols[[1]] - cols[[2]])))
      }
      ll <- ll - apply(ll, 1, lse)
      pr <- exp(ll)
      Pc <- Pc + pr
      real <- real + colSums(pr)
    }
    b0 <- as.vector(apply(sweep(B, 2, log(w), "+"), 1, lse))
    lf0 <- lse(b0); psi <- th[["psi"]]
    wz <- psi * exp(lf0) / (psi * exp(lf0) + 1 - psi)
    real <- real + (Mx - n) * wz * exp(b0 - lf0)
    expd <- expd + psi * Mx * exp(lpi)
  }
  R <- nrow(draws)
  list(Pc = Pc / R, realised = real / R, expected = expd / R)
}

# Simulate detection data for N animals at pixels g with sex sx (0/1) under
# draw th: returns N x J x Kq array of counts (Kq occasions or aggregated).
scr_sim_det <- function(prep, th, g, sx) {
  N <- length(g); J <- prep$J; Kq <- prep$Kq
  lin <- scr_lin(prep, th); linb <- if (prep$behav == "trap") scr_lin(prep, th, behav = 1) else NULL
  Y <- array(0, c(N, J, Kq))
  if (!N) return(Y)
  caught <- matrix(FALSE, N, J)
  for (k in seq_len(Kq)) {
    for (s in seq_len(prep$S)) {
      ii <- which(sx == s - 1); if (!length(ii)) next
      lp <- scr_logp_ks(prep, lin, k, s)[g[ii], , drop = FALSE]
      if (!is.null(linb)) { lpb <- scr_logp_ks(prep, linb, k, s)[g[ii], , drop = FALSE]; lp[caught[ii, ]] <- lpb[caught[ii, ]] }
      U <- matrix(prep$use[, k], length(ii), J, byrow = TRUE)
      yk <- switch(prep$obs_code,
        stats::rbinom(length(lp), as.vector(U), exp(lp)),
        stats::rbinom(length(lp), 1, 1 - exp(as.vector(U) * log1mexp(lp))) * (as.vector(U) > 0),
        stats::rpois(length(lp), as.vector(U) * exp(lp)))
      Y[ii, , k] <- matrix(yk, length(ii), J)
    }
    caught <- caught | (Y[, , k] > 0)
  }
  Y
}

# Posterior predictive check: simulate a new population and survey per draw
# (N ~ Bin(M, psi), centres ~ pi, sex ~ phi, detections from the observation
# model); compare statistics that react to mis-specified detection/space use.
scr_ppc <- function(prep, draws) {
  stat <- function(y, Ey_trap) {
    det <- rowSums(y) > 0; yd <- y[det, , drop = FALSE]
    tt <- colSums(y)
    c(n_detected = sum(det), total_detections = sum(y),
      traps_per_individual = if (sum(det)) mean(rowSums(yd > 0)) else 0,
      ft_trap_totals = sum((sqrt(tt) - sqrt(Ey_trap))^2),
      disp_detections = if (sum(det) > 1) stats::var(rowSums(yd)) / mean(rowSums(yd)) else NA_real_)
  }
  Mx <- if (!is.null(attr(draws, "M"))) attr(draws, "M") else if (!is.null(prep$M_used)) prep$M_used else prep$M
  res <- matrix(NA_real_, nrow(draws), 5); obs <- matrix(NA_real_, nrow(draws), 5)
  caps_rep <- vector("list", nrow(draws))
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    lpi <- scr_lpi(prep, th); w <- scr_sexw(prep, th)
    N <- stats::rbinom(1, Mx, th[["psi"]])
    g <- sample.int(prep$G, N, replace = TRUE, prob = exp(lpi))
    sx <- if (prep$S == 2) stats::rbinom(N, 1, 1 - th[["phi"]]) else rep(0, N)
    Y <- scr_sim_det(prep, th, g, sx)
    yr <- apply(Y, c(1, 2), sum); if (!is.matrix(yr)) yr <- matrix(yr, N, prep$J)
    # expected trap totals (no behavioural response; used for obs and replicate alike)
    lin <- scr_lin(prep, th); Ey <- numeric(prep$J)
    for (s in seq_len(prep$S)) for (k in seq_len(prep$Kq)) {
      P <- exp(scr_logp_ks(prep, lin, k, s))
      if (prep$obs_code == 2) P <- 1 - exp(sweep(log1mexp(log(P)), 2, prep$use[, k], "*")) else P <- sweep(P, 2, prep$use[, k], "*")
      Ey <- Ey + th[["psi"]] * Mx * w[s] * as.vector(exp(lpi) %*% P)
    }
    res[r, ] <- stat(yr, Ey); obs[r, ] <- stat(prep$y, Ey)
    caps_rep[[r]] <- rowSums(yr)[rowSums(yr) > 0]
  }
  colnames(res) <- colnames(obs) <- names(stat(prep$y, rep(0, prep$J)))
  pval <- vapply(seq_len(5), function(k) mean(res[, k] >= obs[, k], na.rm = TRUE), 1)
  names(pval) <- colnames(res)
  list(rep = res, obs = obs, p = pval, caps_rep = caps_rep)
}
