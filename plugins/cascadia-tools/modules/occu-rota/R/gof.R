# occu-rota goodness of fit.
#
# Statistic T(y) = part 1 + part 2:
#  1. Detection-count frequency chi-square: for each species and each number
#     of surveys done n, the number of sites with k = 0..n detections,
#     observed vs expected, sum (O - E)^2 / E. Expected counts mix "absent"
#     (k = 0) with the Poisson-binomial distribution of k given presence
#     (survey-specific p). Checks each species' occupancy x detection model,
#     in particular detection heterogeneity among sites that the detection
#     covariates do not capture (too many 0s and high counts).
#  2. For every species pair, Pearson (O - E)^2 / E on the number of sites
#     where both species were detected, E = sum_i P(both present at i) *
#     P(a detected | present) * P(b detected | present). Checks the
#     dependence part of the model.
# Calibration (independent verification, VALIDATION.md): under the true model
# the bootstrap p-value is uniform (P(p < 0.05) = 0.04); with site-level
# detection heterogeneity SD 1 / 1.5 on the logit scale it rejects in 66% /
# 96% of datasets. Two other statistics are computed for reference only:
# MacKenzie-Bailey chi-square on detection histories (calibrated but weak:
# 6% / 14% power) and the original Freeman-Tukey statistic on per-site
# detection counts (no power, non-uniform p-values).

rota_gof_stat <- function(des, y, Psi, p) {
  y0 <- y; y0[is.na(y0)] <- 0
  S <- des$S; Z <- des$Z
  pres <- Psi %*% Z                                        # M x S marginal presence probs
  # (reference only) MacKenzie-Bailey chi-square on each species' detection histories, by
  #    cohort of sites sharing the same pattern of surveys done. Expected count
  #    of history h = sum over the cohort's sites of P_i(h); all histories never
  #    observed are pooled in one remainder cell (O = 0), which makes the sum
  #    identical to the full chi-square over all 2^J histories.
  cohort <- apply(des$w, 1, paste, collapse = "")
  t_mb <- 0
  for (s in seq_len(S)) {
    lp <- log(pmax(p[, , s], 1e-300)); lq <- log1p(-p[, , s])   # p = 0 on surveys not done
    ys <- y0[, , s]
    for (cc in unique(cohort)) {
      ii <- which(cohort == cc); nc <- length(ii)
      Ys <- ys[ii, , drop = FALSE]
      hk <- apply(Ys, 1, paste, collapse = "")
      H <- Ys[!duplicated(hk), , drop = FALSE]; O <- as.vector(table(factor(hk, levels = hk[!duplicated(hk)])))
      L <- exp(lp[ii, , drop = FALSE] %*% t(H) + lq[ii, , drop = FALSE] %*% t(1 - H))   # nc x nh  P(h | present)
      P <- pres[ii, s] * L + outer(1 - pres[ii, s], as.numeric(rowSums(H) == 0))
      E <- colSums(P)
      t_mb <- t_mb + sum((O - E)^2 / pmax(E, 1e-8)) + max(nc - sum(E), 0)
    }
  }
  # part 2: for every species pair, Pearson (O - E)^2 / E on the number of sites
  #    where both species were detected.
  t2 <- 0
  pdet <- sapply(seq_len(S), function(s) 1 - apply(1 - p[, , s], 1, prod))   # M x S
  for (a in seq_len(S - 1)) for (b in (a + 1):S) {
    O <- sum(apply(y0[, , a], 1, max) * apply(y0[, , b], 1, max))
    E <- sum((Psi %*% (Z[, a] * Z[, b])) * pdet[, a] * pdet[, b])
    t2 <- t2 + (O - E)^2 / max(E, 1e-8)
  }
  # part 1: frequency of detection counts: number of sites with k detections among
  #    the n surveys done, per species and n; expected counts use the
  #    Poisson-binomial distribution of k at each site (survey-specific p).
  ns <- rowSums(des$w); t_cf <- 0
  for (s in seq_len(S)) {
    ps <- p[, , s]; f <- matrix(0, des$M, des$J + 1); f[, 1] <- 1
    for (j in seq_len(des$J)) f <- f * (1 - ps[, j]) + cbind(0, f[, -(des$J + 1), drop = FALSE]) * ps[, j]
    f <- pres[, s] * f; f[, 1] <- f[, 1] + (1 - pres[, s])
    k <- rowSums(y0[, , s])
    for (n in unique(ns)) {
      ii <- which(ns == n)
      E <- colSums(f[ii, seq_len(n + 1), drop = FALSE]); O <- tabulate(k[ii] + 1, n + 1)
      t_cf <- t_cf + sum((O - E)^2 / pmax(E, 1e-8))
    }
  }
  # (legacy, reported for comparison only) Freeman-Tukey on per-site detection counts
  t_ft <- 0
  for (s in seq_len(S)) {
    n <- rowSums(y0[, , s]); E <- pres[, s] * rowSums(p[, , s])
    t_ft <- t_ft + sum((sqrt(n) - sqrt(E))^2)
  }
  c(total = t_cf + t2, count_freq = t_cf, codetect = t2, mb_histories = t_mb, ft_counts = t_ft)
}

rota_sim_y <- function(des, Psi, p) {
  st <- apply(Psi, 1, function(pr) sample.int(des$C, 1, prob = pr))
  zs <- des$Z[st, , drop = FALSE]
  zarr <- as.vector(zs[, rep(seq_len(des$S), each = des$J), drop = FALSE])   # (i, j, s) order
  y <- array(stats::rbinom(des$M * des$J * des$S, 1, as.vector(p) * zarr), c(des$M, des$J, des$S))
  y[rep(des$w == 0, des$S)] <- NA
  y
}

rota_gof <- function(des, eng, settings) {
  set.seed(settings$seed + 7)
  P <- des$P
  stat_at <- function(d, y, par) rota_gof_stat(d, y, rota_psi(rota_F_sites(d, par[seq_len(P)]), d$A),
                                              rota_p_sites(d, par[P + seq_len(d$Q)]))
  if (eng$engine == "unmarked") {
    Psi <- rota_psi(rota_F_sites(des, eng$point[seq_len(P)]), des$A)
    p <- rota_p_sites(des, eng$point[P + seq_len(des$Q)])
    obs <- rota_gof_stat(des, des$y, Psi, p)
    st2 <- merge_list(settings, list(n_draws = 2, penalty = 0))
    # simulate every bootstrap dataset first: the refit reseeds the RNG
    ysim <- lapply(seq_len(settings$gof_nboot), function(b) rota_sim_y(des, Psi, p))
    sims <- t(vapply(seq_len(settings$gof_nboot), function(b) {
      d2 <- des; d2$y <- ysim[[b]]
      e2 <- tryCatch(rota_engine_unmarked(d2, st2), error = function(e) NULL)
      if (is.null(e2)) return(rep(NA_real_, 5))
      stat_at(d2, d2$y, e2$point)
    }, numeric(5)))
    sims <- sims[stats::complete.cases(sims), , drop = FALSE]
    list(method = "parametric bootstrap", nboot = nrow(sims), stat = obs[["total"]],
         components = obs, sims = sims[, 1], sims_components = sims,
         p_value = mean(sims[, 1] >= obs[["total"]]), c_hat = obs[["total"]] / mean(sims[, 1]))
  } else {
    D <- eng$draws
    idx <- unique(round(seq(1, nrow(D), length.out = min(settings$gof_ndraws, nrow(D)))))
    res <- t(vapply(idx, function(r) {
      par <- D[r, ]
      Psi <- rota_psi(rota_F_sites(des, par[seq_len(P)]), des$A)
      p <- rota_p_sites(des, par[P + seq_len(des$Q)])
      c(rota_gof_stat(des, des$y, Psi, p)[["total"]], rota_gof_stat(des, rota_sim_y(des, Psi, p), Psi, p)[["total"]])
    }, numeric(2)))
    list(method = "posterior predictive check", nboot = nrow(res), stat = mean(res[, 1]),
         sims = res[, 2], p_value = mean(res[, 2] >= res[, 1]), c_hat = mean(res[, 1]) / mean(res[, 2]))
  }
}
