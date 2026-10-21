# MacKenzie & Bailey (2004) goodness-of-fit for occupancy models.
#
# Pearson chi-square over detection histories. Sites are grouped into
# cohorts by their pattern of missing surveys; within a cohort every
# possible history h (2^K for K surveys done) gets observed count O_h and
# expected count E_h = sum over sites of Pr(h | site's psi and p). For the
# dynamic model the statistic is computed within each season using the
# model-implied (unconditional) occupancy for that season, and the seasonal
# statistics are summed (as AICcmodavg::mb.gof.test does for colext fits).
# Reference distribution: parametric bootstrap (simulate from the fitted
# model, refit, recompute). c-hat = observed / mean(bootstrap).

occ_mb_chisq <- function(des, psi_mt, p, y = des$y, max_k = 12) {
  chi <- 0
  for (t in seq_len(des$T)) {
    wt <- des$w[, , t, drop = FALSE][, , 1]; wt <- matrix(wt, des$M)
    pat <- apply(wt, 1, paste, collapse = "")
    for (g in unique(pat)) {
      ii <- which(pat == g); jj <- which(wt[ii[1], ] == 1); K <- length(jj)
      if (K == 0 || K > max_k) next
      H <- as.matrix(expand.grid(rep(list(0:1), K)))
      pp <- matrix(p[ii, jj, t], length(ii))
      lp <- log(pp); lq <- log1p(-pp)
      # Pr(h | site) for every history x site
      pr <- exp(H %*% t(lp) + (1 - H) %*% t(lq))
      pr <- sweep(pr, 2, psi_mt[ii, t], "*")
      allzero <- rowSums(H) == 0
      pr[allzero, ] <- pr[allzero, ] + (1 - psi_mt[ii, t])
      E <- rowSums(pr)
      yy <- matrix(y[ii, jj, t], length(ii))
      key_obs <- apply(yy, 1, paste, collapse = "")
      key_H <- apply(H, 1, paste, collapse = "")
      O <- as.vector(table(factor(key_obs, levels = key_H)))
      ok <- E > 0
      chi <- chi + sum((O[ok] - E[ok])^2 / E[ok])
    }
  }
  chi
}

occ_gof_unmarked <- function(des, eng, settings) {
  nb <- settings$gof_nboot
  if (is.null(nb) || nb <= 0 || !isTRUE(eng$diag$converged)) return(NULL)
  stat <- function(dd, cf) {
    pr <- occ_probs(dd, occ_split_coef(cf, dd))
    occ_mb_chisq(dd, occ_marginal_psi(dd, pr), pr$p)
  }
  obs <- stat(des, eng$point)
  pr <- occ_probs(des, occ_split_coef(eng$point, des))
  boot <- rep(NA_real_, nb)
  for (b in seq_len(nb)) {
    yb <- occ_simulate_y(des, pr)$y
    if (sum(yb, na.rm = TRUE) == 0) next
    db <- des; db$y <- yb
    fb <- tryCatch(occ_fit_unmarked(db, settings, y = yb, se = FALSE), error = function(e) NULL)
    if (is.null(fb)) next
    boot[b] <- stat(db, fb$coef)
  }
  boot <- boot[is.finite(boot)]
  list(method = "MacKenzie-Bailey chi-square, parametric bootstrap", stat = obs, nboot = length(boot),
       mean_boot = mean(boot), c_hat = obs / mean(boot), p_value = mean(boot >= obs),
       boot = boot)  # replicate statistics, kept for the gof figure
}

# Bayesian version: posterior predictive check of the same statistic.
occ_gof_nimble <- function(des, eng, settings) {
  nd <- min(settings$gof_ndraws, nrow(eng$draws))
  if (is.null(nd) || nd <= 0) return(NULL)
  idx <- round(seq(1, nrow(eng$draws), length.out = nd))
  to <- tr <- numeric(nd)
  for (k in seq_along(idx)) {
    pr <- occ_probs(des, occ_split_coef(eng$draws[idx[k], ], des))
    ps <- occ_marginal_psi(des, pr)
    to[k] <- occ_mb_chisq(des, ps, pr$p)
    tr[k] <- occ_mb_chisq(des, ps, pr$p, y = occ_simulate_y(des, pr)$y)
  }
  list(method = "MacKenzie-Bailey chi-square, posterior predictive check", stat = mean(to), nboot = nd,
       mean_boot = mean(tr), c_hat = mean(to) / mean(tr), p_value = mean(tr >= to),
       obs = to, rep = tr)  # per-draw statistics, kept for the gof figure
}
