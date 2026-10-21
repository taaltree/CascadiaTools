# Community simulator: species abundances from Gaussian (unimodal) niche
# models along one or two environmental gradients, plus species-specific
# responses to a site grouping factor (location effect) and group-specific
# site-level heterogeneity (dispersion effect).
#
# Site j, species i:
#   u1_j, u2_j  gradient positions on [0, 1] (Gaussian copula, correlation env_cor),
#               reported in original units: env_k = lo_k + u_k (hi_k - lo_k)
#   log lambda_ij = log(effort * h_i) - (u1_j - o1_i)^2 / (2 t1_i^2)
#                   - (u2_j - o2_i)^2 / (2 t2_i^2) + delta_{i, g(j)} + eps_ij
#   eps_ij ~ N(-s_g^2 / 2, s_g^2)  (mean-preserving site x species noise; s_g is
#                                   site_noise_sd, times dispersion_mult in the LAST group)
#   y_ij ~ Poisson(lambda_ij)       (presence data: 1{y_ij > 0})
#   delta_{i,1} = 0, delta_{i,g} ~ N(0, group_effect_sd^2) for g > 1
# Species: optima o ~ U(-0.1, 1.1), tolerances t1 = tol1 * exp(N(0, 0.3^2)),
# t2 = tol2 * exp(N(0, 0.3^2)), peak abundance h ~ logN(log(abund_max), abund_sd^2).
# Sites with no species at all are redrawn (the population is non-empty sites:
# Hellinger, Bray-Curtis and Sorensen are undefined for empty sites).
#
# TRUTH (all estimands): the population value of the same statistic for an
# infinite set of sites drawn from this process (same species pool, same
# environmental distribution, groups in the design's proportions), computed by
# large Monte Carlo samples that are independent of the dataset:
#   * R2adj_* (RDA on Hellinger data): population R^2 of the linear regression
#     of the Hellinger-transformed community on the environmental set(s),
#     from truth_n (default 40,000) sites. Fractions follow the same algebra
#     as varpart: unique_k = R2(all) - R2(all but k), shared = R2(both) -
#     unique_1 - unique_2. Adjusted R^2 is an (approximately) unbiased
#     estimator of these (Peres-Neto et al. 2006).
#   * R2adj_<group>_permanova (Bray-Curtis): 1 - sum_g pi_g E[d^2 | both in g] /
#     E[d^2 | random pair], the n -> infinity limit of PERMANOVA's
#     SS_within/SS_total, with pi_g = n_g / n; expectations over truth_pairs
#     independent site pairs per group combination.
#   * beta_*_pair_mean: expected Sorensen / Simpson (turnover) /
#     nestedness-resultant dissimilarity between two sites, weighting
#     within- and between-group pairs as they occur in the design, from the
#     same independent pairs. turnover_share = E[beta_sim] / E[beta_sor].
# These describe the OBSERVED community at this sampling effort (detection is
# not modelled): an undersampled survey has more apparent turnover than the
# latent community, and the truth reflects that by construction.

`%||%` <- function(a, b) if (is.null(a)) b else a

com_species_pool <- function(params, S) {
  list(o1 = stats::runif(S, -0.1, 1.1), o2 = stats::runif(S, -0.1, 1.1),
       t1 = params$tol1 * exp(stats::rnorm(S, 0, 0.3)),
       t2 = params$tol2 * exp(stats::rnorm(S, 0, 0.3)),
       h = exp(stats::rnorm(S, log(params$abund_max), params$abund_sd)))
}

# Draw gradient positions u (n x 2) on [0, 1].
com_draw_u <- function(n, rho) {
  z1 <- stats::rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(n)
  cbind(stats::pnorm(z1), stats::pnorm(z2))
}

# Simulate counts for sites with given group indices g (1..G); empty sites are
# redrawn (env and counts). Returns list(u, y).
com_draw_sites <- function(g, pool, params, design, delta, n_grad) {
  n <- length(g); S <- length(pool$h)
  G <- nrow(delta)
  sg <- rep(params$site_noise_sd, G); sg[G] <- sg[G] * (params$dispersion_mult %||% 1)
  u <- matrix(NA_real_, n, 2); y <- matrix(0L, n, S)
  todo <- seq_len(n); it <- 0
  while (length(todo) && it < 200) {
    it <- it + 1
    m <- length(todo)
    uu <- com_draw_u(m, design$env_cor %||% 0)
    eta <- matrix(log((design$effort %||% 1) * pool$h), m, S, byrow = TRUE) -
      outer(uu[, 1], pool$o1, "-")^2 / matrix(2 * pool$t1^2, m, S, byrow = TRUE)
    if (n_grad >= 2) eta <- eta - outer(uu[, 2], pool$o2, "-")^2 / matrix(2 * pool$t2^2, m, S, byrow = TRUE)
    gg <- g[todo]
    eta <- eta + delta[gg, , drop = FALSE]
    s <- sg[gg]
    eta <- eta + matrix(stats::rnorm(m * S, -s^2 / 2, s), m, S)  # recycles s by row (column-major)
    yy <- matrix(stats::rpois(m * S, exp(eta)), m, S)
    ok <- rowSums(yy) > 0
    u[todo[ok], ] <- uu[ok, , drop = FALSE]; y[todo[ok], ] <- yy[ok, , drop = FALSE]
    todo <- todo[!ok]
  }
  if (length(todo)) stop("simulator could not produce non-empty sites: raise abund_max or effort")
  list(u = u, y = y)
}

com_group_sizes <- function(design) {
  G <- as.integer(design$n_groups %||% 0)
  if (G < 1) return(integer())
  if (!is.null(design$group_sizes)) {
    gs <- as.integer(unlist(design$group_sizes)); stopifnot(length(gs) == G); return(gs)
  }
  n <- as.integer(design$n_sites)
  base <- rep(n %/% G, G); base[seq_len(n %% G)] <- base[seq_len(n %% G)] + 1L
  base
}

com_simulate <- function(params, design, seed) {
  set.seed(seed)
  S <- as.integer(params$n_species)
  n_grad <- as.integer(design$n_gradients %||% 2)
  gs <- com_group_sizes(design)
  G <- max(1L, length(gs))
  n <- if (length(gs)) sum(gs) else as.integer(design$n_sites)
  g <- if (length(gs)) rep(seq_len(G), gs) else rep(1L, n)
  pool <- com_species_pool(params, S)
  delta <- matrix(0, G, S)
  if (G > 1) delta[-1, ] <- stats::rnorm((G - 1) * S, 0, params$group_effect_sd)
  sim <- com_draw_sites(g, pool, params, design, delta, n_grad)
  # names and units
  env_names <- unlist(design$env_names %||% c("elevation", "soil_moisture"))[seq_len(n_grad)]
  rng <- design$env_ranges %||% list(c(200, 1800), c(5, 45))
  sp <- unlist(design$species_names %||% sprintf("sp%02d", seq_len(S)))[seq_len(S)]
  site <- sprintf("S%03d", seq_len(n))
  env <- data.frame(site = site, stringsAsFactors = FALSE)
  for (k in seq_len(n_grad)) env[[env_names[k]]] <- round(rng[[k]][1] + sim$u[, k] * diff(rng[[k]]), 2)
  gname <- design$group_name %||% "habitat"
  glev <- unlist(design$group_levels %||% LETTERS[seq_len(G)])[seq_len(G)]
  if (length(gs)) env[[gname]] <- factor(glev[g], levels = glev)
  # shuffle row order so groups are not blocked in the file
  ord <- sample.int(n)
  env <- env[ord, , drop = FALSE]; env$site <- site; rownames(env) <- NULL
  Y <- sim$y[ord, , drop = FALSE]
  dtype <- design$data_type %||% "abundance"
  if (dtype == "presence") Y[] <- as.integer(Y > 0)
  dimnames(Y) <- list(site, sp)
  data <- list(comm = Y, env = env, env_vars = env_names,
               group = if (length(gs)) gname else NULL,
               partition = stats::setNames(as.list(env_names), env_names),
               data_type = dtype)
  truth <- com_truth(pool, params, design, delta, n_grad, g, env_names, gname, dtype,
                     seed = seed)
  list(data = data, truth = truth)
}

# ---- truth by large independent Monte Carlo ---------------------------------

hellinger <- function(Y) { rs <- rowSums(Y); sqrt(Y / rs) }

r2_lm <- function(Yc, X) {
  # Yc column-centred; X without intercept (centred inside)
  X <- scale(as.matrix(X), scale = FALSE)
  q <- qr(X)
  fitv <- qr.fitted(q, Yc)
  sum(fitv^2) / sum(Yc^2)
}

com_truth <- function(pool, params, design, delta, n_grad, g, env_names, gname, dtype, seed) {
  set.seed(seed + 7777777L)
  G <- nrow(delta)
  pi_g <- tabulate(g, G) / length(g)
  bin <- function(y) if (dtype == "presence") (y > 0) * 1 else y
  # (1) RDA R2 truths
  Nt <- as.integer(design$truth_n %||% 40000)
  gt <- if (G > 1) sample.int(G, Nt, replace = TRUE, prob = pi_g) else rep(1L, Nt)
  big <- com_draw_sites(gt, pool, params, design, delta, n_grad)
  Yh <- hellinger(bin(big$y)); Yc <- scale(Yh, scale = FALSE)
  U <- big$u[, seq_len(n_grad), drop = FALSE]
  tr <- c()
  if (n_grad == 1) {
    tr[paste0("R2adj_", env_names[1])] <- r2_lm(Yc, U[, 1])
    tr["R2adj_env_total"] <- tr[[1]]
  } else {
    r1 <- r2_lm(Yc, U[, 1]); r2 <- r2_lm(Yc, U[, 2]); r12 <- r2_lm(Yc, U)
    tr[paste0("R2adj_", env_names[1])] <- r1
    tr[paste0("R2adj_", env_names[2])] <- r2
    tr[paste0("R2adj_unique_", env_names[1])] <- r12 - r2
    tr[paste0("R2adj_unique_", env_names[2])] <- r12 - r1
    tr["R2adj_shared"] <- r1 + r2 - r12
    tr["R2adj_env_total"] <- r12
  }
  # (2) pair-based truths: independent pairs for every group combination
  np <- as.integer(design$truth_pairs %||% 60000)
  E <- list()
  for (a in seq_len(G)) for (b in a:G) {
    s1 <- com_draw_sites(rep(a, np), pool, params, design, delta, n_grad)$y
    s2 <- com_draw_sites(rep(b, np), pool, params, design, delta, n_grad)$y
    E[[paste(a, b)]] <- pair_stats(bin(s1), bin(s2))
  }
  # weights of within/between pairs as they occur in the sample design
  ng <- tabulate(g, G); n <- length(g)
  w <- matrix(0, G, G)
  for (a in seq_len(G)) for (b in a:G)
    w[a, b] <- if (a == b) ng[a] * (ng[a] - 1) / 2 else ng[a] * ng[b]
  w <- w / sum(w)
  avg <- function(k) sum(vapply(names(E), function(nm) {
    ab <- as.integer(strsplit(nm, " ")[[1]]); w[ab[1], ab[2]] * E[[nm]][[k]] }, 0))
  tr["beta_sor_pair_mean"] <- avg("sor")
  tr["beta_sim_pair_mean"] <- avg("sim")
  tr["beta_sne_pair_mean"] <- avg("sne")
  tr["turnover_share"] <- tr[["beta_sim_pair_mean"]] / tr[["beta_sor_pair_mean"]]
  if (G > 1) {
    # Bray-Curtis PERMANOVA R2 limit, groups at proportions pi_g
    Ed2 <- function(a, b) E[[paste(min(a, b), max(a, b))]][["bray2"]]
    tot <- 0; for (a in seq_len(G)) for (b in seq_len(G)) tot <- tot + pi_g[a] * pi_g[b] * Ed2(a, b)
    wit <- sum(vapply(seq_len(G), function(a) pi_g[a] * Ed2(a, a), 0))
    tr[paste0("R2adj_", gname, "_permanova")] <- 1 - wit / tot
  }
  tr
}

# Per-pair dissimilarities for row-aligned matrices A, B (same species).
pair_stats <- function(A, B) {
  bray <- rowSums(abs(A - B)) / rowSums(A + B)
  pa <- A > 0; pb <- B > 0
  a <- rowSums(pa & pb); b <- rowSums(pa & !pb); c <- rowSums(!pa & pb)
  mn <- pmin(b, c); mx <- pmax(b, c)
  sor <- (b + c) / (2 * a + b + c)
  sim <- ifelse(a + mn > 0, mn / (a + mn), 0)
  list(bray2 = mean(bray^2), sor = mean(sor), sim = mean(sim), sne = mean(sor - sim))
}
