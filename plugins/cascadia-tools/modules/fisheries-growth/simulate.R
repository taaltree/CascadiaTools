# Simulation for fisheries-growth.
#
# Ages come from a catch-at-age distribution: fish recruit to the gear at
# age_min and numbers then decline with total mortality Z, so
#   P(age = a) ∝ exp(-Z (a - age_min)),  a = age_min..max_age.
# High Z means few old fish, so Linf is poorly identified.
#
# Lengths follow von Bertalanffy, L = Linf (1 - exp(-K (age - t0))), with
# multiplicative (lognormal, sigma on the log scale) or additive (normal,
# sigma in mm) error.
#
# Weights: W = a_g L^b exp(e), e ~ N(0, lw_sigma). a_g = lw_a * condition_g,
# where condition_g is the group's weight at a given length relative to the
# reference group (first group).
#
# Two kinds of group effects:
#   fixed (design$hierarchical = FALSE): per-group Linf, K, t0 vectors in params.
#   hierarchical (TRUE): log Linf_g = log(Linf_typical) + Linf_beta * z_g + e_g,
#     e_g ~ N(0, sd_log_Linf). log K_g is built the same way. t0 is common to
#     all groups. z_g is the group covariate standardised by its sample mean
#     and SD, exactly as fit() standardises it, so the betas are effects on
#     the log scale per 1 SD of the covariate.

fg_group_names <- function(G) if (G <= 26) LETTERS[seq_len(G)] else sprintf("G%02d", seq_len(G))

fg_vb <- function(age, Linf, K, t0) Linf * (1 - exp(-K * (age - t0)))

fg_rep <- function(x, G) { x <- unlist(x); if (length(x) == 1) rep(x, G) else { stopifnot(length(x) == G); x } }

simulate_fg <- function(params, design, seed) {
  set.seed(seed)
  hier <- isTRUE(design$hierarchical)
  G <- design$n_groups
  gn <- fg_group_names(G)
  covname <- if (is.null(design$cov_name)) "temperature" else design$cov_name
  truth <- c()
  zc <- NULL; gcov <- NULL
  if (hier) {
    raw <- design$cov_mean + design$cov_sd * stats::rnorm(G)
    zc <- (raw - mean(raw)) / stats::sd(raw)
    lL <- log(params$Linf_typical) + params$Linf_beta * zc + params$sd_log_Linf * stats::rnorm(G)
    lK <- log(params$K_typical) + params$K_beta * zc + params$sd_log_K * stats::rnorm(G)
    Linf <- exp(lL); K <- exp(lK); t0 <- rep(params$t0, G)
    gcov <- data.frame(group = gn, cov = round(raw, 3)); names(gcov)[2] <- covname
    truth <- c(Linf_typical = params$Linf_typical, K_typical = params$K_typical,
               t0 = params$t0, sd_log_Linf = params$sd_log_Linf, sd_log_K = params$sd_log_K)
    truth[paste0("Linf_beta_", covname)] <- params$Linf_beta
    truth[paste0("K_beta_", covname)] <- params$K_beta
  } else {
    Linf <- fg_rep(params$Linf, G); K <- fg_rep(params$K, G); t0 <- fg_rep(params$t0, G)
  }
  npg <- fg_rep(design$n_per_group, G)
  ages <- design$age_min:design$max_age
  pa <- exp(-design$Z * (ages - design$age_min))
  cond <- unlist(params$condition)
  # condition is per group; a vector of the wrong length (e.g. the 3-group
  # default in a 10-group hierarchical design) means 'no condition differences'
  cond <- if (is.null(cond) || !(length(cond) %in% c(1, G))) rep(1, G) else fg_rep(cond, G)
  lmin <- if (is.null(design$min_length)) 0 else design$min_length
  dfs <- lapply(seq_len(G), function(g) {
    # Length-selective gear (negative control): fish shorter than
    # design$min_length are never caught, so the youngest ages are
    # represented only by their fastest-growing members.
    nn <- if (lmin > 0) 20 * npg[g] else npg[g]
    a <- ages[sample.int(length(ages), nn, replace = TRUE, prob = pa)] + design$age_offset
    mu <- fg_vb(a, Linf[g], K[g], t0[g])
    L <- if (params$error == "multiplicative") mu * exp(stats::rnorm(nn, 0, params$sigma))
         else mu + stats::rnorm(nn, 0, params$sigma)
    if (lmin > 0) { keep <- which(L >= lmin)[seq_len(npg[g])]; a <- a[keep]; L <- L[keep] }
    L <- pmax(round(L), 1)
    W <- params$lw_a * cond[g] * L^params$lw_b * exp(stats::rnorm(npg[g], 0, params$lw_sigma))
    W <- round(W, 1)
    if (design$prop_weighed < 1) W[stats::runif(npg[g]) > design$prop_weighed] <- NA
    data.frame(fish_id = sprintf("%s-%04d", gn[g], seq_len(npg[g])), age = a, length = L,
               group = gn[g], weight = W)
  })
  df <- do.call(rbind, dfs)
  df$group <- factor(df$group, levels = gn)
  # Per-group truths (realised values in the hierarchical design)
  for (g in seq_len(G)) {
    truth[paste0("Linf_", gn[g])] <- Linf[g]
    truth[paste0("K_", gn[g])] <- K[g]
    if (!hier) truth[paste0("t0_", gn[g])] <- t0[g]
    truth[paste0("age_at_90pct_Linf_", gn[g])] <- t0[g] + log(10) / K[g]
    if (g > 1) {
      truth[paste0("Linf_diff_", gn[g], "_vs_", gn[1])] <- Linf[g] - Linf[1]
      truth[paste0("K_diff_", gn[g], "_vs_", gn[1])] <- K[g] - K[1]
      if (!hier) truth[paste0("t0_diff_", gn[g], "_vs_", gn[1])] <- t0[g] - t0[1]
      truth[paste0("condition_ratio_", gn[g], "_vs_", gn[1])] <- cond[g] / cond[1]
    }
  }
  truth["residual_sd"] <- params$sigma
  truth["lw_b"] <- params$lw_b
  truth["lw_a"] <- params$lw_a
  truth["lw_residual_sd"] <- params$lw_sigma
  lref <- design$reference_length
  truth[paste0("mean_weight_at_", lref, "_", gn[1])] <- params$lw_a * cond[1] * lref^params$lw_b * exp(params$lw_sigma^2 / 2)
  data <- list(df = df, group_cov = gcov, cov_name = if (hier) covname else NULL,
               ref_group = gn[1], reference_length = lref,
               units = list(length = "mm", weight = "g", age = "years"),
               ws = NULL, species = NULL)
  list(data = data, truth = truth)
}
