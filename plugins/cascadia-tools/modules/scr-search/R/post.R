# scr-search: reported quantities from engine draws, realised density
# surfaces and posterior predictive checks. th = one draw (named numeric:
# psi, lam0, sigma, bD1.., bW1..), internal units.

ss_lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

ss_lpi <- function(prep, th) {
  PD <- ncol(prep$XD)
  if (!PD) return(rep(-log(prep$G), prep$G))
  lw <- as.vector(prep$XD %*% th[paste0("bD", seq_len(PD))])
  lw - ss_lse(lw)
}
# log(E_u exp(w_u' alpha)) per unit
ss_lEu <- function(prep, th) {
  PW <- ncol(prep$W)
  if (!PW) return(prep$logE)
  prep$logE + as.vector(prep$W %*% th[paste0("bW", seq_len(PW))])
}
# effective effort per cell
ss_Et <- function(prep, th) as.vector(tapply(exp(ss_lEu(prep, th)), factor(prep$ucell, levels = seq_len(prep$C)), sum))
# expected samples per pixel from all search (local evaluation, as the engine)
ss_Lambda <- function(prep, th) {
  Et <- ss_Et(prep, th)
  K <- exp(-prep$d2gc / (2 * th[["sigma"]]^2))
  if (is.finite(prep$R)) K[prep$d2gc > prep$R^2] <- 0
  th[["lam0"]] * as.vector(K %*% Et)
}

# log-likelihood of each detected individual's records at every pixel (n x G);
# -Inf outside the individual's local window
ss_indiv_ll <- function(prep, th, base) {
  n <- prep$n; G <- prep$G
  ll <- matrix(-Inf, n, G)
  lEu <- ss_lEu(prep, th); a <- -1 / (2 * th[["sigma"]]^2); l0 <- log(th[["lam0"]])
  bern <- prep$observation == "bernoulli"
  for (i in seq_len(n)) {
    s <- prep$sets[[i]]
    v <- base[s]
    rr <- which(prep$rec$ind == i)
    for (r in rr) {
      lr <- l0 + lEu[prep$rec$unit[r]] + a * prep$d2gc[s, prep$rec$cell[r]]
      v <- v + if (bern) ifelse(exp(lr) > 30, exp(lr), log(expm1(exp(lr)))) else prep$rec$count[r] * lr
    }
    ll[i, s] <- v
  }
  ll
}

ss_area_conv <- function(prep) if (!is.na(prep$units) && prep$units %in% c("m", "km")) 1 else 1 / prep$cscale^2
ss_len_conv <- function(prep) if (!is.na(prep$units) && prep$units %in% c("m", "km")) 1 else 1 / prep$cscale
ss_hr_conv <- function(prep) if (!is.na(prep$units) && prep$units %in% c("m", "km")) 1 else prep$cscale^2

# Reported draws for one chain
ss_report_chain <- function(prep, m, idx) {
  m <- m[idx, , drop = FALSE]
  S <- nrow(m)
  Mx <- if (!is.null(prep$M_used)) prep$M_used else prep$M
  N <- vapply(seq_len(S), function(r) {
    th <- m[r, ]
    lf0 <- ss_lse(ss_lpi(prep, th) - ss_Lambda(prep, th))
    w <- th[["psi"]] * exp(lf0) / (th[["psi"]] * exp(lf0) + 1 - th[["psi"]])
    prep$n + stats::rbinom(1, Mx - prep$n, w)
  }, 1)
  lam <- m[, "lam0"] * ss_len_conv(prep)   # per km (or per map unit)
  out <- cbind(m[, "psi"] * Mx / prep$area * ss_area_conv(prep), N, lam, m[, "sigma"] * prep$cscale,
               pi * (sqrt(stats::qchisq(0.95, 2)) * m[, "sigma"])^2 * ss_hr_conv(prep),
               m[, "lam0"] * sqrt(2 * pi) * m[, "sigma"])
  colnames(out) <- c(ss_density_name(prep$units), "abundance", ss_lam_name(prep$units), "sigma",
                     "home_range_area_95", "samples_per_track_through_centre")
  for (k in seq_along(prep$dcov)) { out <- cbind(out, m[, paste0("bD", k)]); colnames(out)[ncol(out)] <- paste0("density_beta_", prep$dcov[k]) }
  for (k in seq_along(prep$wcov)) { out <- cbind(out, m[, paste0("bW", k)]); colnames(out)[ncol(out)] <- paste0("detection_beta_", prep$wcov[k]) }
  out
}

ss_thin_idx <- function(nr, max_draws) if (nr <= max_draws) seq_len(nr) else unique(round(seq(1, nr, length.out = max_draws)))
ss_pool <- function(chains, ndraw) { all <- do.call(rbind, chains); all[ss_thin_idx(nrow(all), ndraw), , drop = FALSE] }

# Posterior activity-centre probabilities of detected animals (n x G),
# realised density (animals per pixel) and expected density surface.
ss_surfaces <- function(prep, draws, M) {
  n <- prep$n; G <- prep$G
  Pc <- matrix(0, n, G); real <- numeric(G); expd <- numeric(G); pdet <- numeric(G)
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    lpi <- ss_lpi(prep, th); Lam <- ss_Lambda(prep, th)
    base <- lpi - Lam
    if (n) {
      ll <- ss_indiv_ll(prep, th, base)
      ll <- ll - apply(ll, 1, ss_lse)
      pr <- exp(ll); Pc <- Pc + pr; real <- real + colSums(pr)
    }
    lf0 <- ss_lse(base); psi <- th[["psi"]]
    w <- psi * exp(lf0) / (psi * exp(lf0) + 1 - psi)
    real <- real + (M - n) * w * exp(base - lf0)
    expd <- expd + psi * M * exp(lpi)
    pdet <- pdet + (1 - exp(-Lam))
  }
  R <- nrow(draws)
  list(Pc = Pc / R, realised = real / R, expected = expd / R, pdetect = pdet / R)
}

# Posterior predictive check: per draw simulate N ~ Bin(M, psi), centres ~ pi,
# samples per unit ~ Poisson (or Bernoulli) and compare statistics.
ss_ppc <- function(prep, draws, M) {
  stat <- function(Y, Ecell) {
    det <- rowSums(Y) > 0; yd <- Y[det, , drop = FALSE]
    ct <- as.vector(tapply(colSums(Y), factor(prep$ucell, levels = seq_len(prep$C)), sum))
    cells_per <- if (sum(det)) mean(apply(yd, 1, function(v) length(unique(prep$ucell[v > 0])))) else 0
    c(n_detected = sum(det), total_samples = sum(Y), cells_per_individual = cells_per,
      ft_cell_totals = sum((sqrt(ct) - sqrt(Ecell))^2),
      disp_samples = if (sum(det) > 1) stats::var(rowSums(yd)) / mean(rowSums(yd)) else NA_real_)
  }
  Yobs <- matrix(0, prep$n, prep$U)
  Yobs[cbind(prep$rec$ind, prep$rec$unit)] <- prep$rec$count
  bern <- prep$observation == "bernoulli"
  res <- obs <- matrix(NA_real_, nrow(draws), 5)
  caps_rep <- vector("list", nrow(draws))
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    lpi <- ss_lpi(prep, th); lEu <- ss_lEu(prep, th)
    K <- exp(-prep$d2gc / (2 * th[["sigma"]]^2)) * (prep$d2gc <= prep$R^2)  # G x C
    lamGU <- th[["lam0"]] * K[, prep$ucell, drop = FALSE] * rep(exp(lEu), each = prep$G)   # G x U
    EGU <- if (bern) 1 - exp(-lamGU) else lamGU
    Eu <- th[["psi"]] * M * as.vector(exp(lpi) %*% EGU)
    Ecell <- as.vector(tapply(Eu, factor(prep$ucell, levels = seq_len(prep$C)), sum))
    N <- stats::rbinom(1, M, th[["psi"]])
    g <- sample.int(prep$G, N, replace = TRUE, prob = exp(lpi))
    L <- lamGU[g, , drop = FALSE]
    Y <- matrix(if (bern) stats::rbinom(length(L), 1, 1 - exp(-L)) else stats::rpois(length(L), L), N, prep$U)
    res[r, ] <- stat(Y, Ecell); obs[r, ] <- stat(Yobs, Ecell)
    caps_rep[[r]] <- rowSums(Y)[rowSums(Y) > 0]
  }
  colnames(res) <- colnames(obs) <- c("n_detected", "total_samples", "cells_per_individual", "ft_cell_totals", "disp_samples")
  pval <- vapply(seq_len(5), function(k) mean(res[, k] >= obs[, k], na.rm = TRUE), 1)
  names(pval) <- colnames(res)
  list(rep = res, obs = obs, p = pval, caps_rep = caps_rep, obs_caps = rowSums(Yobs))
}
