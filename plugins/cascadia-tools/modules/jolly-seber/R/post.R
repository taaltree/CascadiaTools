# jolly-seber: from parameter draws to reported quantities.
#
# Latent counts (N_t, entrants B_t, Nsuper) are drawn exactly from their
# conditional posterior given each parameter draw. Every individual's latent
# history is one of T(T+1)/2 + 1 "paths": never entered, or alive from entry
# occasion e to last-alive occasion d (e <= d). For each unique capture
# pattern u (including the M - n all-zero augmented rows) the posterior over
# (group, path) is prior(path | group) x pi_group x P(y_u | path); the w_u
# animals with that pattern are allocated to paths with one multinomial draw.
# Summing the path probabilities gives the forward-algorithm likelihood used
# in the Nimble model (checked in dev/jolly-seber/check_paths.R).

js_paths <- function(T) {
  pe <- unlist(lapply(seq_len(T), function(e) rep(e, T - e + 1)))
  pd <- unlist(lapply(seq_len(T), function(e) e:T))
  alive <- outer(seq_along(pe), seq_len(T), function(k, t) pe[k] <= t & t <= pd[k])
  list(e = pe, d = pd, P = length(pe), alive = alive * 1)
}

# Parameter draws (matrix rows) -> list of natural-scale quantities.
js_param_draws <- function(prep, S) {
  T <- prep$T
  col <- function(v, k) S[, sprintf("%s[%d]", v, k)]
  b <- vapply(seq_len(T), function(t) col("b", t), numeric(nrow(S)))
  if (nrow(S) == 1) b <- matrix(b, 1)
  psi <- S[, "psi"]
  bp <- vapply(seq_len(ncol(prep$Xp)), function(k) col("beta_p", k), numeric(nrow(S)))
  bphi <- vapply(seq_len(ncol(prep$Xphi)), function(k) col("beta_phi", k), numeric(nrow(S)))
  if (nrow(S) == 1) { bp <- matrix(bp, 1); bphi <- matrix(bphi, 1) }
  lp <- bp %*% t(prep$Xp)            # draws x T
  lphi <- bphi %*% t(prep$Xphi)      # draws x (T-1)
  bg <- if (prep$G == 2) S[, "beta_grp"] else rep(0, nrow(S))
  pig2 <- if (prep$G == 2) S[, "pi_grp"] else rep(0, nrow(S))
  cb <- t(apply(cbind(0, b[, -T, drop = FALSE]), 1, cumsum))
  if (T == 1) cb <- matrix(0, nrow(S), 1)
  gamma <- psi * b / (1 - psi * cb)
  list(b = b, psi = psi, p = stats::plogis(lp), lphi = lphi, bg = bg, pig2 = pig2,
       gamma = gamma, beta_p = bp, beta_phi = bphi)
}

# One draw: latent path counts per pattern. Returns list(N, B, Nsuper, cells)
js_sample_paths <- function(prep, pr, k, paths, need_cells = FALSE) {
  T <- prep$T; P <- paths$P
  p <- pr$p[k, ]; gam <- pr$gamma[k, ]
  lq <- log1p(-p); lpp <- log(p)
  # cumulative log-lik per pattern; C[, t+1] = sum_{s<=t}
  Y <- prep$pat
  ll <- Y * rep(lpp, each = nrow(Y)) + (1 - Y) * rep(lq, each = nrow(Y))
  C <- cbind(0, t(apply(ll, 1, cumsum)))
  if (T == 1) C <- cbind(0, ll)
  LL <- C[, paths$d + 1, drop = FALSE] - C[, paths$e, drop = FALSE]
  LL[!prep$valid] <- -Inf
  lne <- c(0, cumsum(log1p(-gam)))                 # log P(not entered by t), lne[t+1]
  lent <- lne[paths$e] + log(gam[paths$e])
  out_lp <- matrix(0, 2, P)
  for (g in 1:2) {
    lph <- if (g == 1) pr$lphi[k, ] else pr$lphi[k, ] + pr$bg[k]
    lsv <- log(stats::plogis(lph)) * prep$dt       # log survival per interval
    csv <- c(0, cumsum(lsv))
    ldie <- c(log1p(-exp(lsv)), 0)                 # death after d (none after T)
    out_lp[g, ] <- lent + (csv[paths$d] - csv[paths$e]) + ldie[paths$d]
  }
  lpi <- c(log1p(-pr$pig2[k]), log(pr$pig2[k]))
  never_lp <- lne[T + 1]
  U <- prep$U
  cnt <- matrix(0, 2, P); nev <- 0
  cells <- if (need_cells) vector("list", U) else NULL
  for (u in seq_len(U)) {
    gs <- if (prep$pat_g[u] == 0) (if (prep$G == 2) 1:2 else 1) else prep$pat_g[u]
    lw <- unlist(lapply(gs, function(g) lpi[g] * (prep$G == 2) + out_lp[g, ] + LL[u, ]))
    isz <- u == U
    if (isz) lw <- c(lw, never_lp)
    m <- max(lw); pw <- exp(lw - m)
    dr <- as.vector(stats::rmultinom(1, prep$pat_w[u], pw))
    if (isz) { nev <- dr[length(dr)]; dr <- dr[-length(dr)] }
    dm <- matrix(dr, length(gs), P, byrow = TRUE)
    cnt[gs, ] <- cnt[gs, ] + dm
    if (need_cells) cells[[u]] <- colSums(dm)
  }
  tot <- colSums(cnt)
  N <- as.vector(tot %*% paths$alive)
  B <- vapply(2:T, function(t) sum(tot[paths$e == t]), 1)
  list(N = N, B = B, Nsuper = sum(tot), Nsuper2 = sum(cnt[2, ]), cells = cells, by_group = cnt)
}

js_valid_mask <- function(prep, paths) {
  f <- apply(prep$pat, 1, function(r) if (any(r > 0)) min(which(r > 0)) else NA)
  l <- apply(prep$pat, 1, function(r) if (any(r > 0)) max(which(r > 0)) else NA)
  v <- outer(seq_len(prep$U), seq_len(paths$P), function(u, k) is.na(f[u]) | (paths$e[k] <= f[u] & paths$d[k] >= l[u]))
  v
}

js_thin_idx <- function(n, m) if (n <= m) seq_len(n) else unique(round(seq(1, n, length.out = m)))

# Reported-quantity matrix for one chain (rows = draws).
js_report_chain <- function(prep, S, latent = NULL) {
  T <- prep$T
  pr <- js_param_draws(prep, S)
  nd <- nrow(S)
  if (is.null(latent)) {
    paths <- js_paths(T)
    prep$valid <- js_valid_mask(prep, paths)
    N <- matrix(0, nd, T); B <- matrix(0, nd, T - 1); Ns <- numeric(nd); Ns2 <- numeric(nd)
    for (k in seq_len(nd)) {
      sp <- js_sample_paths(prep, pr, k, paths)
      N[k, ] <- sp$N; B[k, ] <- sp$B; Ns[k] <- sp$Nsuper; Ns2[k] <- sp$Nsuper2
    }
  } else { N <- latent$N; B <- latent$B; Ns <- latent$Nsuper; Ns2 <- latent$Nsuper2 }
  colnames(N) <- paste0("abundance_", seq_len(T))
  colnames(B) <- paste0("entrants_", 2:T)
  bb <- pr$b; colnames(bb) <- paste0("entry_prob_", seq_len(T))
  f <- B / N[, -T, drop = FALSE]; colnames(f) <- paste0("recruitment_rate_", seq_len(T - 1))
  lam <- N[, -1, drop = FALSE] / N[, -T, drop = FALSE]; colnames(lam) <- paste0("growth_rate_", seq_len(T - 1))
  f[!is.finite(f)] <- NA; lam[!is.finite(lam)] <- NA
  gm <- (N[, T] / N[, 1])^(1 / sum(prep$dt)); gm[!is.finite(gm)] <- NA
  out <- cbind(superpopulation_size = Ns, N, B, bb, f, lam, growth_rate_mean = gm)
  # survival (per time unit) and detection on natural scales
  sv <- stats::plogis(pr$beta_phi); colnames(sv) <- colnames(prep$Xphi)
  out <- cbind(out, sv)
  if (prep$G == 2) {
    l2 <- prep$group_levels[2]
    out <- cbind(out, stats::setNames(data.frame(pr$bg), paste0("survival_beta_", l2)) |> as.matrix(),
                 # realised share of level 2 among the animals in the superpopulation
                 stats::setNames(data.frame(Ns2 / Ns), paste0("prop_", l2)) |> as.matrix())
    if (ncol(prep$Xphi) == 1) out <- cbind(out, stats::setNames(data.frame(stats::plogis(pr$beta_phi[, 1] + pr$bg)),
                                                                paste0("survival_", l2)) |> as.matrix())
  }
  dp <- pr$beta_p
  if (prep$Kp_lev == T) { dp <- stats::plogis(dp) }
  else { dp[, 1] <- stats::plogis(dp[, 1]) }
  colnames(dp) <- colnames(prep$Xp)
  cbind(out, dp)
}

js_pool <- function(chains, m) {
  per <- max(1, ceiling(m / length(chains)))
  do.call(rbind, lapply(chains, function(x) x[js_thin_idx(nrow(x), per), , drop = FALSE]))
}

# Posterior predictive checks (Freeman-Tukey discrepancies):
#  individual: animals' total captures vs expected sum_t z_it p_t (all M rows)
#  occasion:   animals caught per occasion vs N_t p_t
js_ppc <- function(prep, S) {
  T <- prep$T
  pr <- js_param_draws(prep, S)
  paths <- js_paths(T)
  prep$valid <- js_valid_mask(prep, paths)
  ytot <- rowSums(prep$pat)
  nt_obs <- colSums(prep$ch)
  res <- matrix(NA_real_, nrow(S), 4, dimnames = list(NULL, c("ind_obs", "ind_rep", "occ_obs", "occ_rep")))
  for (k in seq_len(nrow(S))) {
    sp <- js_sample_paths(prep, pr, k, paths, need_cells = TRUE)
    p <- pr$p[k, ]
    Ep <- as.vector(paths$alive %*% p)               # expected captures per path
    fi_obs <- 0; fi_rep <- 0; nt_rep <- numeric(T)
    for (u in seq_len(prep$U)) {
      cc <- sp$cells[[u]]
      for (j in which(cc > 0)) {
        c0 <- cc[j]
        fi_obs <- fi_obs + c0 * (sqrt(ytot[u]) - sqrt(Ep[j]))^2
        tt <- which(paths$alive[j, ] > 0)
        yr <- matrix(stats::rbinom(c0 * length(tt), 1, rep(p[tt], each = c0)), c0, length(tt))
        fi_rep <- fi_rep + sum((sqrt(rowSums(yr)) - sqrt(Ep[j]))^2)
        nt_rep[tt] <- nt_rep[tt] + colSums(yr)
      }
    }
    E_t <- sp$N * p
    res[k, ] <- c(fi_obs, fi_rep, sum((sqrt(nt_obs) - sqrt(E_t))^2), sum((sqrt(nt_rep) - sqrt(E_t))^2))
  }
  list(draws = res, p = c(individual = mean(res[, "ind_rep"] >= res[, "ind_obs"]),
                          occasion = mean(res[, "occ_rep"] >= res[, "occ_obs"])))
}
