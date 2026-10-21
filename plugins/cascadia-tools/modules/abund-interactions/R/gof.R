# Posterior predictive check for both species (Twining et al. 2025 used a
# Pearson chi-square per species submodel). We use Freeman-Tukey discrepancy
# on site totals of detection days, conditional on latent N_i and z_i:
#   T(y) = sum_i (sqrt(Y_i) - sqrt(E_i))^2,  Y_i = sum_j y_ij,  E_i = sum_j K_ij p_ij
# p-value = P(T(y_rep) >= T(y_obs)); values < 0.05 or > 0.95 flag lack of fit.
abi_gof <- function(des, eng, ndraws = 500) {
  n <- nrow(eng$N)
  idx <- unique(round(seq(1, n, length.out = min(n, ndraws))))
  ft <- function(Y, E) sum((sqrt(Y) - sqrt(E))^2)
  YA <- rowSums(des$yA); YS <- rowSums(des$yS)
  out <- matrix(NA_real_, length(idx), 4, dimnames = list(NULL, c("obsA", "repA", "obsS", "repS")))
  for (k in seq_along(idx)) {
    d <- idx[k]
    lp <- abi_linpred(des, eng$draws[d, ])
    N <- eng$N[d, ]; z <- eng$z[d, ]
    pA <- 1 - (1 - lp$r)^N
    pS <- lp$p * z
    EA <- rowSums(des$K * pA); ES <- rowSums(des$K * pS)
    yrA <- matrix(stats::rbinom(length(pA), as.vector(des$K), as.vector(pA)), des$M)
    yrS <- matrix(stats::rbinom(length(pS), as.vector(des$K), as.vector(pS)), des$M)
    out[k, ] <- c(ft(YA, EA), ft(rowSums(yrA), EA), ft(YS, ES), ft(rowSums(yrS), ES))
  }
  list(stats = out, p = c(dominant = mean(out[, "repA"] >= out[, "obsA"]),
                          subordinate = mean(out[, "repS"] >= out[, "obsS"])), ndraws = length(idx))
}
