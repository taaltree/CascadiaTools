# scr-open: posterior predictive checks. For each posterior draw a new
# population (M augmented rows, alive paths from gamma/phi, centres from pi)
# and survey (Binomial(K_jt, p) per trap and period) are simulated and summary
# statistics compared with the observed ones.

so_ppc_stats <- function(y, seen) {
  # y: n x J x T counts; seen n x T
  det <- rowSums(seen) > 0
  ni <- sum(det)
  pt <- colSums(seen)
  trapsper <- if (sum(seen)) sum(apply(y > 0, c(1, 3), sum)) / sum(seen) else 0
  c(stats::setNames(pt, paste0("detected_period_", seq_along(pt))),
    n_detected = ni,
    detected_in_2plus_periods = sum(rowSums(seen) > 1),
    total_detections = sum(y),
    traps_per_individual_period = trapsper,
    gaps_in_histories = sum(apply(seen, 1, function(s) { w <- which(s); if (length(w) < 2) 0 else sum(!s[min(w):max(w)]) })))
}

so_ppc <- function(prep, draws, M) {
  TT <- prep$TT; J <- prep$J; G <- prep$G
  obs <- so_ppc_stats(prep$y, prep$seen)
  res <- matrix(NA_real_, nrow(draws), length(obs), dimnames = list(NULL, names(obs)))
  nper <- matrix(0, nrow(draws), TT)
  for (r in seq_len(nrow(draws))) {
    th <- draws[r, ]
    gam <- th[paste0("gam", seq_len(TT))]; phi <- th[["phi"]]
    st <- matrix(1L, M, TT)
    st[, 1] <- ifelse(stats::runif(M) < gam[1], 2L, 1L)
    for (t in 2:TT) {
      u <- stats::runif(M); pv <- st[, t - 1]
      st[, t] <- ifelse(pv == 1, ifelse(u < gam[t], 2L, 1L), ifelse(pv == 2, ifelse(u < phi, 2L, 3L), 3L))
    }
    ever <- which(rowSums(st == 2) > 0)
    g <- sample.int(G, length(ever), replace = TRUE, prob = exp(so_lpi(prep, th)))
    p <- th[["p0"]] * exp(-prep$d2[g, , drop = FALSE] / (2 * th[["sigma"]]^2))
    y <- array(0L, c(length(ever), J, TT))
    for (t in seq_len(TT)) {
      al <- st[ever, t] == 2
      if (any(al)) y[al, , t] <- stats::rbinom(sum(al) * J, rep(prep$usage[, t], each = sum(al)), p[al, , drop = FALSE])
    }
    sn <- apply(y, c(1, 3), sum) > 0
    if (length(ever) == 1) sn <- matrix(sn, 1)
    keep <- rowSums(sn) > 0
    res[r, ] <- so_ppc_stats(y[keep, , , drop = FALSE], sn[keep, , drop = FALSE])
  }
  pval <- vapply(names(obs), function(k) mean(res[, k] >= obs[[k]], na.rm = TRUE), 1)
  list(rep = res, obs = obs, p = pval)
}
