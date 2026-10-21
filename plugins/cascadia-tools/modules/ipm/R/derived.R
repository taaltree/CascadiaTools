# ipm: derived quantities from posterior draws: population totals, growth
# rates, immigration, and the life-table response experiments (LTRE) of
# Schaub & Kery (2022) section 9.4 (after Koons et al. 2016, 2017):
#
#  * fixed-design ("random design") LTRE: contribution of each vital rate
#    and of population structure to the temporal variance of the population
#    growth rate, c_k = sum_l cov(theta_k, theta_l) s_k s_l, with
#    sensitivities s evaluated at the mean rates and mean structure;
#  * transient LTRE: contribution of each rate to the change in growth
#    rate between successive years, (theta_k[t+1] - theta_k[t]) s_k, with
#    s evaluated at the mean of the two years.
#
# Growth rate as a function of rates and stage proportions (book eq. 9.4):
#   lambda = ((n1 + n3) rho phij f1 + n2 rho phij f2 + (phia + om)(n1 + n2 + n3)) / (n1 + n2 + n3)

ipm_lambda_expr <- function(rho = 0.5)
  substitute(((n1 + n3) * RHO * phij * f1 + n2 * RHO * phij * f2 + (phia + om) * (n1 + n2 + n3)) / (n1 + n2 + n3),
             list(RHO = rho))

ipm_ltre_keys <- function(flags) {
  k <- c("phij", "phia")
  if (flags$IMM == 1) k <- c(k, "om")
  k <- c(k, "f1")
  if (flags$AGEF == 1) k <- c(k, "f2", "n1", "n2", if (flags$IMM == 1) "n3")
  k
}

ipm_ltre_names <- function(flags) {
  nm <- c(phij = "ltre_juvenile_survival", phia = "ltre_adult_survival", om = "ltre_immigration",
          f1 = if (flags$AGEF == 1) "ltre_fecundity_yearling" else "ltre_fecundity",
          f2 = "ltre_fecundity_adult")
  nm
}

# d: list of draw matrices (draws x time): phij, phia, omega (T-1 cols),
# f1, f2, N1, N2, N3 (T cols). means: list of vectors (one per draw) of the
# mean rates at which sensitivities are evaluated.
ipm_ltre_draws <- function(d, means, flags, rho = 0.5, transient = FALSE) {
  T1 <- ncol(d$phij)
  nd <- nrow(d$phij)
  if (is.null(d$omega)) d$omega <- d$phij * 0
  if (is.null(d$N3)) d$N3 <- d$N1 * 0
  if (is.null(d$f2)) d$f2 <- d$f1
  if (is.null(means$omega)) means$omega <- rep(0, nd)
  if (is.null(means$f2)) means$f2 <- means$f1
  Ntot <- d$N1 + d$N2 + d$N3
  Ntot[Ntot <= 0] <- NA
  n1 <- (d$N1 / Ntot)[, 1:T1, drop = FALSE]
  n2 <- (d$N2 / Ntot)[, 1:T1, drop = FALSE]
  n3 <- (d$N3 / Ntot)[, 1:T1, drop = FALSE]
  series <- list(phij = d$phij, phia = d$phia, om = d$omega,
                 f1 = d$f1[, 1:T1, drop = FALSE], f2 = d$f2[, 1:T1, drop = FALSE],
                 n1 = n1, n2 = n2, n3 = n3)
  keys <- ipm_ltre_keys(flags)
  lam <- ipm_lambda_expr(rho)
  rm_na <- function(x) rowMeans(x, na.rm = TRUE)
  mu <- list(phij = means$phij, phia = means$phia, om = means$omega, f1 = means$f1, f2 = means$f2,
             n1 = rm_na(n1), n2 = rm_na(n2), n3 = rm_na(n3))
  sens <- sapply(keys, function(k) rep_len(eval(stats::D(lam, k), envir = mu), nd))
  sens <- matrix(sens, nd, length(keys), dimnames = list(NULL, keys))
  cont <- matrix(NA_real_, nd, length(keys), dimnames = list(NULL, keys))
  for (s in seq_len(nd)) {
    dp <- sapply(keys, function(k) series[[k]][s, ])
    dp <- matrix(dp, ncol = length(keys))
    V <- stats::var(dp, na.rm = TRUE)
    cont[s, ] <- rowSums(V * outer(sens[s, ], sens[s, ]))
  }
  nm <- ipm_ltre_names(flags)
  out <- data.frame(row.names = NULL, check.names = FALSE,
                    stats::setNames(as.data.frame(cont[, intersect(keys, names(nm)), drop = FALSE]),
                                    nm[intersect(keys, names(nm))]))
  if (flags$AGEF == 1) out$ltre_population_structure <- rowSums(cont[, intersect(keys, c("n1", "n2", "n3")), drop = FALSE])
  out$var_lambda_explained <- rowSums(cont)
  res <- list(contrib = out, contrib_full = cont)
  if (transient && T1 >= 2) {
    # differences and means between successive years: (T1 - 1) steps
    gd <- function(x) x[, 2:T1, drop = FALSE] - x[, 1:(T1 - 1), drop = FALSE]
    gm <- function(x) (x[, 2:T1, drop = FALSE] + x[, 1:(T1 - 1), drop = FALSE]) / 2
    ms <- lapply(series, gm)
    tc <- array(NA_real_, c(nd, T1 - 1, length(keys)), dimnames = list(NULL, NULL, keys))
    for (k in keys) tc[, , k] <- gd(series[[k]]) * eval(stats::D(lam, k), envir = ms)
    res$transient <- tc
  }
  res
}

# Draw matrices for annual quantities from stacked posterior samples.
ipm_draw_mat <- function(dr, base, n) {
  cn <- paste0(base, "[", seq_len(n), "]")
  if (!all(cn %in% colnames(dr))) return(NULL)
  dr[, cn, drop = FALSE]
}
