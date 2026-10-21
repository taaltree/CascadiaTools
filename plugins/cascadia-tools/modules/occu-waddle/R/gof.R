# occu-waddle: goodness of fit by posterior (Nimble) or parametric
# (maximum likelihood: draws from the asymptotic normal) predictive checks.
# Two discrepancy statistics, both computed against the analytic
# expectation under each parameter draw:
#   T_pattern: chi-square over the 4 site co-detection classes (neither /
#              dominant only / subordinate only / both detected) -- checks the
#              co-occurrence structure the interaction parameters describe.
#   T_counts:  chi-square over per-site detection counts of each species --
#              checks detection heterogeneity.

wd_expect <- function(des, q) {
  M <- des$M; J <- des$J; o <- matrix(des$obs_v, M, J, byrow = TRUE)
  pr <- function(v) { m <- matrix(v, M, J, byrow = TRUE); m[!o] <- 0; m }
  rA <- pr(q$rA); rBA <- pr(q$rBA); rBa <- pr(q$rBa); pA <- pr(q$pA); pB <- pr(q$pB)
  w1 <- q$psiA * q$psiBA; w2 <- q$psiA * (1 - q$psiBA); w3 <- (1 - q$psiA) * q$psiBa
  w4 <- (1 - q$psiA) * (1 - q$psiBa)
  P <- function(m) apply(m, 1, prod)
  noA <- w1 * P(1 - rA) + w2 * P(1 - pA) + w3 + w4
  noB <- w1 * P((1 - rA) * (1 - rBa) + rA * (1 - rBA)) + w2 + w3 * P(1 - pB) + w4
  none <- w1 * P((1 - rA) * (1 - rBa)) + w2 * P(1 - pA) + w3 * P(1 - pB) + w4
  cls <- cbind(none, noB - none, noA - none, 1 - noA - noB + none)
  nA <- w1 * rowSums(rA) + w2 * rowSums(pA)
  nB <- w1 * rowSums(rA * rBA + (1 - rA) * rBa) + w3 * rowSums(pB)
  list(cls = colSums(cls), nA = nA, nB = nB)
}

wd_sim_y <- function(des, q) {
  M <- des$M; J <- des$J
  zA <- stats::rbinom(M, 1, q$psiA); zB <- stats::rbinom(M, 1, ifelse(zA == 1, q$psiBA, q$psiBa))
  zAr <- rep(zA, each = J); zBr <- rep(zB, each = J); both <- zAr * zBr
  ya <- stats::rbinom(M * J, 1, ifelse(both == 1, q$rA, q$pA) * zAr)
  yb <- stats::rbinom(M * J, 1, ifelse(both == 1, ifelse(ya == 1, q$rBA, q$rBa), q$pB) * zBr)
  ya[!des$obs_v] <- 0; yb[!des$obs_v] <- 0
  list(ya = ya, yb = yb)
}

wd_stats <- function(des, ya, yb, ex) {
  M <- des$M; J <- des$J
  na <- rowSums(matrix(ya, M, J, byrow = TRUE)); nb <- rowSums(matrix(yb, M, J, byrow = TRUE))
  oc <- c(sum(na == 0 & nb == 0), sum(na > 0 & nb == 0), sum(na == 0 & nb > 0), sum(na > 0 & nb > 0))
  c(pattern = sum((oc - ex$cls)^2 / pmax(ex$cls, 0.5)),
    counts = sum((na - ex$nA)^2 / pmax(ex$nA, 0.5)) + sum((nb - ex$nB)^2 / pmax(ex$nB, 0.5)))
}

wd_gof <- function(des, eng, settings) {
  D <- eng$draws
  nd <- min(settings$gof_ndraws, nrow(D))
  rows <- unique(round(seq(1, nrow(D), length.out = nd)))
  res <- t(vapply(rows, function(r) {
    q <- wd_probs(des, D[r, ])
    ex <- wd_expect(des, q)
    sim <- wd_sim_y(des, q)
    c(obs = wd_stats(des, des$ya, des$yb, ex), rep = wd_stats(des, sim$ya, sim$yb, ex))
  }, numeric(4)))
  list(method = if (eng$engine == "nimble") "posterior predictive check" else "parametric predictive check (draws from the MLE's normal approximation)",
       ndraws = length(rows),
       obs_pattern = res[, "obs.pattern"], rep_pattern = res[, "rep.pattern"],
       obs_counts = res[, "obs.counts"], rep_counts = res[, "rep.counts"],
       p_pattern = mean(res[, "rep.pattern"] >= res[, "obs.pattern"]),
       p_counts = mean(res[, "rep.counts"] >= res[, "obs.counts"]))
}
