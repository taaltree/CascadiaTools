# scr-smr: engine draws (internal units) -> reported quantities, density
# surfaces and posterior predictive checks.

smr_thin_idx <- function(nr, max_draws) if (nr <= max_draws) seq_len(nr) else unique(round(seq(1, nr, length.out = max_draws)))
smr_q <- function(x, p = c(0.025, 0.975)) stats::quantile(x, p, na.rm = TRUE, names = FALSE)
smr_lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

smr_rate_name <- function(prep) if (prep$pois) "resight_rate" else "resight_detection"

# Reported draws for one chain matrix m (Nimble column names).
smr_report_chain <- function(prep, m, idx = seq_len(nrow(m))) {
  m <- m[idx, , drop = FALSE]
  N <- m[, "N"]
  A <- prep$area * smr_area_conv(prep)
  out <- cbind(N / A, N, prep$m / N, N - prep$m, m[, "sigma"] * prep$cscale,
               pi * (sqrt(stats::qchisq(0.95, 2)) * m[, "sigma"])^2 * smr_area_conv(prep), m[, "lam0"])
  colnames(out) <- c(smr_density_name(prep$units), "abundance", "marked_fraction", "abundance_unmarked", "sigma",
                     "home_range_area_95", smr_rate_name(prep))
  if (prep$has_tm) out <- cbind(out, id_probability = m[, "delta"])
  if (prep$has_marking) out <- cbind(out, marking_detection = m[, "pm0"])
  for (k in seq_along(prep$dcov)) { out <- cbind(out, m[, sprintf("bD[%d]", k)]); colnames(out)[ncol(out)] <- paste0("density_beta_", prep$dcov[k]) }
  for (k in seq_along(prep$pcov)) { out <- cbind(out, m[, sprintf("bP[%d]", k)]); colnames(out)[ncol(out)] <- paste0("detection_beta_", prep$pcov[k]) }
  out
}

# Pool thinned draws across chains (all monitored columns)
smr_pool <- function(chains, ndraw) {
  all <- do.call(rbind, lapply(chains, as.matrix))
  all[smr_thin_idx(nrow(all), ndraw), , drop = FALSE]
}

smr_lpi <- function(prep, th) {
  PD <- ncol(prep$XD)
  if (!PD) return(rep(-log(prep$G), prep$G))
  b <- th[sprintf("bD[%d]", seq_len(PD))]
  lw <- as.vector(prep$XD %*% b)
  lw - smr_lse(lw)
}
smr_lr0 <- function(prep, th) {
  PP <- ncol(prep$XP)
  b <- if (PP) th[sprintf("bP[%d]", seq_len(PP))] else numeric(0)
  eta <- if (PP) as.vector(prep$XP %*% b) else rep(0, prep$J)
  if (prep$pois) log(th[["lam0"]]) + eta else stats::plogis(stats::qlogis(th[["lam0"]]) + eta, log.p = TRUE)
}

# s/z columns of a draw -> M x 2 matrix and z vector
smr_sz <- function(th, M) {
  s <- cbind(th[sprintf("s[%d, 1]", seq_len(M))], th[sprintf("s[%d, 2]", seq_len(M))])
  list(s = s, z = th[sprintf("z[%d]", seq_len(M))])
}

# Realised density surface (expected animals per mask pixel, averaged over
# draws), expected (covariate-model) surface, posterior mean centres of the
# marked animals.
smr_surfaces <- function(prep, draws, M) {
  G <- prep$G
  rc <- smr_rect(prep)
  map <- integer(rc$GR); map[rc$g] <- seq_len(G)
  real <- numeric(G); realU <- numeric(G); expd <- numeric(G)
  cen <- matrix(0, max(prep$m, 1), 2)
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    sz <- smr_sz(th, M)
    on <- sz$z == 1
    ix <- floor(sz$s[on, 1] / prep$pix); iy <- floor(sz$s[on, 2] / prep$pix)
    g <- map[ix + iy * rc$nx + 1]
    g <- g[g > 0]
    real <- real + tabulate(g, G)
    unm <- on & seq_len(M) > prep$m
    ixu <- floor(sz$s[unm, 1] / prep$pix); iyu <- floor(sz$s[unm, 2] / prep$pix)
    gu <- map[ixu + iyu * rc$nx + 1]; gu <- gu[gu > 0]
    realU <- realU + tabulate(gu, G)
    expd <- expd + th[["N"]] * exp(smr_lpi(prep, th))
    if (prep$m) cen <- cen + sz$s[seq_len(prep$m), , drop = FALSE]
  }
  R <- nrow(draws)
  list(realised = real / R, realised_unmarked = realU / R, expected = expd / R, centres = cen / R)
}

# Posterior predictive checks conditional on the latent centres (Royle et al.
# 2014 ch. 8 style Freeman-Tukey discrepancies), plus totals.
smr_ppc <- function(prep, draws, M) {
  J <- prep$J; m <- prep$m; Kj <- prep$Kj
  ft <- function(o, e) sum((sqrt(o) - sqrt(e))^2)
  st_names <- c("unmarked_total", "identified_total", "marked_resighted", "ft_unmarked_cameras", "ft_marked_individuals",
                if (prep$has_tm) "unidentified_total")
  rep <- obs <- matrix(NA_real_, nrow(draws), length(st_names), dimnames = list(NULL, st_names))
  sight_rep <- vector("list", nrow(draws))
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    sz <- smr_sz(th, M)
    sg <- th[["sigma"]]
    lr0 <- smr_lr0(prep, th)
    d2 <- outer(sz$s[, 1], prep$txy[, 1], "-")^2 + outer(sz$s[, 2], prep$txy[, 2], "-")^2
    lam <- exp(sweep(-d2 / (2 * sg^2), 2, lr0, "+"))
    dl <- if (prep$has_tm) th[["delta"]] else 1
    wU <- sz$z * (seq_len(M) > m)
    EU <- Kj * as.vector(crossprod(lam, wU))
    nUr <- stats::rpois(J, EU)
    if (m) {
      lm <- lam[seq_len(m), , drop = FALSE]
      EY <- sweep(dl * lm, 2, Kj, "*")
      Yr <- if (prep$pois) matrix(stats::rpois(m * J, EY), m, J)
            else matrix(stats::rbinom(m * J, rep(Kj, each = m), pmin(dl * lm, 1)), m, J)
      EI <- rowSums(EY)
      yrow <- rowSums(prep$yID); yrrow <- rowSums(Yr)
    } else { Yr <- matrix(0, 0, J); EI <- numeric(0); yrow <- yrrow <- numeric(0) }
    rep[r, 1:5] <- c(sum(nUr), sum(Yr), sum(yrrow > 0), ft(nUr, EU), ft(yrrow, EI))
    obs[r, 1:5] <- c(sum(prep$nU), sum(prep$yID), sum(yrow > 0), ft(prep$nU, EU), ft(yrow, EI))
    if (prep$has_tm) {
      ET <- Kj * (1 - dl) * colSums(if (m) lam[seq_len(m), , drop = FALSE] else matrix(0, 1, J))
      rep[r, 6] <- sum(stats::rpois(J, ET)); obs[r, 6] <- sum(prep$nTm)
    }
    sight_rep[[r]] <- list(per_ind = yrrow, per_cam_unmarked = nUr, EU = EU)
  }
  pval <- vapply(st_names, function(k) mean(rep[, k] >= obs[, k], na.rm = TRUE), 1)
  list(rep = rep, obs = obs, p = pval, sight_rep = sight_rep)
}
