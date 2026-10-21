# scr-open: from parameter draws to the reported quantities.
#
# For engines that integrate out activity centres and alive states (nimble
# marginal, nimble hmm, numpyro), realised quantities are drawn per posterior
# draw from their exact conditional posterior given the parameters and the
# data: each animal's activity-centre pixel is drawn from P(g | history), then
# its alive-state path by forward-filtering backward-sampling (FFBS) given
# that pixel. Augmented (never detected) rows are handled the same way with an
# all-zero history. This yields joint draws of N_t, recruits B_t, the
# superpopulation, per-capita recruitment B_t / N_{t-1} and growth N_t / N_{t-1}
# (between-period correlation preserved), and realised density surfaces.

so_lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
so_log1mexp <- function(lp) ifelse(lp > -0.693, log(-expm1(lp)), log1p(-exp(lp)))

so_lpi <- function(prep, th) {
  PD <- ncol(prep$XD)
  if (!PD) return(rep(-log(prep$G), prep$G))
  lw <- as.vector(prep$XD %*% th[paste0("bD", seq_len(PD))])
  lw - so_lse(lw)
}

# log P(detect at j | pixel g) pieces: lg = logit-ish (G x J), base (G x T)
so_emis <- function(prep, th) {
  lp <- log(th[["p0"]]) - prep$d2 / (2 * th[["sigma"]]^2)
  l1 <- so_log1mexp(lp)
  list(lg = lp - l1, base = l1 %*% prep$usage)
}

# log alive-emission for detected individuals: array n x G x T
so_e_det <- function(prep, em) {
  n <- prep$n; G <- prep$G; TT <- prep$TT
  e <- array(rep(as.vector(t(em$base)), each = n), c(n, TT, G))   # n x T x G
  if (length(prep$nzc)) {
    contrib <- prep$nzc * t(em$lg[, prep$nzj, drop = FALSE])       # nnz x G
    key <- (prep$nzt - 1) * n + prep$nzi                             # index into n x T
    add <- rowsum(contrib, key, reorder = TRUE)
    kk <- as.integer(rownames(add))
    em2 <- matrix(e, n * TT, G)
    em2[kk, ] <- em2[kk, ] + add
    e <- array(em2, c(n, TT, G))
  }
  aperm(e, c(1, 3, 2))                                              # n x G x T
}

# Forward pass in probability space. E: R x G x T alive emission (already
# shifted for seen periods), U: R x T (1 = period unseen). Returns alpha
# array R x G x T x 3.
so_forward <- function(E, U, gam, phi) {
  d <- dim(E); R <- d[1]; G <- d[2]; TT <- d[3]
  A <- array(0, c(R, G, TT, 3))
  a1 <- (1 - gam[1]) * matrix(U[, 1], R, G); a2 <- gam[1] * E[, , 1]; a3 <- matrix(0, R, G)
  A[, , 1, 1] <- a1; A[, , 1, 2] <- a2; A[, , 1, 3] <- a3
  for (t in 2:TT) {
    u <- matrix(U[, t], R, G)
    n1 <- a1 * (1 - gam[t]) * u
    n2 <- (a1 * gam[t] + a2 * phi) * E[, , t]
    n3 <- (a2 * (1 - phi) + a3) * u
    a1 <- n1; a2 <- n2; a3 <- n3
    A[, , t, 1] <- a1; A[, , t, 2] <- a2; A[, , t, 3] <- a3
  }
  A
}

# Backward sampling of state paths. al: R x T x 3 forward probabilities at
# each row's chosen pixel. Returns R x T matrix of states 1/2/3.
so_backsample <- function(al, gam, phi) {
  R <- dim(al)[1]; TT <- dim(al)[2]
  st <- matrix(0L, R, TT)
  samp <- function(P) {             # R x 3 unnormalised -> state
    P <- P / rowSums(P); u <- stats::runif(nrow(P))
    1L + (u > P[, 1]) + (u > P[, 1] + P[, 2])
  }
  st[, TT] <- samp(matrix(al[, TT, ], R, 3))
  for (t in (TT - 1):1) {
    nx <- st[, t + 1]
    # Tr(s -> nx) for transition into period t+1
    tr <- cbind(ifelse(nx == 1, 1 - gam[t + 1], ifelse(nx == 2, gam[t + 1], 0)),
                ifelse(nx == 2, phi, ifelse(nx == 3, 1 - phi, 0)),
                ifelse(nx == 3, 1, 0))
    P <- matrix(al[, t, ], R, 3) * tr
    st[, t] <- samp(P)
  }
  st
}

so_gumbel_pick <- function(LW) {     # rows of log-weights -> sampled column
  LW <- LW - log(-log(matrix(stats::runif(length(LW)), nrow(LW))))
  max.col(LW, ties.method = "first")
}

# One parameter draw -> realised states for all M rows (+ pixel draws).
so_realise <- function(prep, th, M, keep_pc = FALSE) {
  gam <- th[paste0("gam", seq_len(prep$TT))]; phi <- th[["phi"]]
  lpi <- so_lpi(prep, th); em <- so_emis(prep, th)
  n <- prep$n; G <- prep$G; TT <- prep$TT
  st_det <- NULL; g_det <- integer(); Pc <- NULL
  if (n) {
    U <- 1 - prep$seen * 1
    Et <- vector("list", TT)
    a1 <- NULL; a2 <- NULL; a3 <- NULL
    for (t in seq_len(TT)) {
      et <- matrix(em$base[, t], n, G, byrow = TRUE)
      k <- which(prep$nzt == t)
      if (length(k)) {
        add <- rowsum(prep$nzc[k] * t(em$lg[, prep$nzj[k], drop = FALSE]), prep$nzi[k])
        ii <- as.integer(rownames(add))
        sub <- et[ii, , drop = FALSE] + add
        mx <- sub[cbind(seq_along(ii), max.col(sub, "first"))]
        et[ii, ] <- sub - mx                 # shift seen periods (per-row max)
      }
      E <- exp(et); Et[[t]] <- E
      u <- U[, t]
      if (t == 1) {
        a1 <- (1 - gam[1]) * u; a2 <- gam[1] * E; a3 <- matrix(0, n, G)
      } else {
        n2 <- (a1 * gam[t] + a2 * phi) * E
        a3 <- (a2 * (1 - phi) + a3) * u
        a1 <- a1 * (1 - gam[t]) * u
        a2 <- n2
      }
    }
    L <- a1 + a2 + a3
    LW <- log(L) + matrix(lpi, n, G, byrow = TRUE)
    g_det <- so_gumbel_pick(LW)
    if (keep_pc) { Pc <- exp(LW - apply(LW, 1, so_lse)) }
    Es <- vapply(seq_len(TT), function(t) Et[[t]][cbind(seq_len(n), g_det)], numeric(n))
    if (n == 1) Es <- matrix(Es, 1)
    A <- so_forward(array(Es, c(n, 1, TT)), U, gam, phi)
    st_det <- so_backsample(array(A[, 1, , ], c(n, TT, 3)), gam, phi)
  }
  # undetected rows: one all-zero history per pixel
  nu <- M - n
  E0 <- array(exp(em$base), c(1, G, TT))
  A0 <- so_forward(E0, matrix(1, 1, TT), gam, phi)
  L0 <- A0[1, , TT, 1] + A0[1, , TT, 2] + A0[1, , TT, 3]
  w0 <- exp(lpi) * L0
  g_u <- sample.int(G, nu, replace = TRUE, prob = w0)
  al0 <- array(0, c(nu, TT, 3))
  for (s in 1:3) al0[, , s] <- matrix(A0[1, , , s][g_u, ], nu, TT)
  st_u <- so_backsample(al0, gam, phi)
  list(states = rbind(st_det, st_u), g = c(g_det, g_u), Pc = Pc, lf0 = log(sum(w0)))
}

so_names <- function(prep) {
  lab <- prep$period_labels; TT <- prep$TT
  list(N = paste0("abundance_period_", lab), D = paste0(so_density_stem(prep$units), "_period_", lab),
       B = paste0("recruits_period_", lab[-1]), f = paste0("per_capita_recruitment_period_", lab[-1]),
       lam = paste0("growth_rate_period_", lab[-1]), gam = paste0("entry_prob_period_", lab))
}

# Reported draws for one chain. m: engine draws (internal names); idx: rows.
so_report_chain <- function(prep, m, idx, M, surf = NULL) {
  m <- m[idx, , drop = FALSE]
  S <- nrow(m); TT <- prep$TT
  nm <- so_names(prep)
  if (all(paste0("N", seq_len(TT)) %in% colnames(m))) {
    Nt <- m[, paste0("N", seq_len(TT)), drop = FALSE]
    Bt <- m[, paste0("B", seq_len(TT)), drop = FALSE]
    Ns <- m[, "Nsuper"]
  } else {
    Nt <- matrix(0, S, TT); Bt <- matrix(0, S, TT); Ns <- numeric(S)
    for (r in seq_len(S)) {
      re <- so_realise(prep, m[r, ], M, keep_pc = !is.null(surf))
      st <- re$states
      al <- st == 2
      Nt[r, ] <- colSums(al)
      Bt[r, ] <- c(Nt[r, 1], colSums(st[, -1, drop = FALSE] == 2 & st[, -TT, drop = FALSE] == 1))
      Ns[r] <- sum(rowSums(al) > 0)
      if (!is.null(surf)) {
        for (t in seq_len(TT)) surf$real[, t] <- surf$real[, t] + tabulate(re$g[al[, t]], prep$G)
        if (!is.null(re$Pc)) surf$Pc <- surf$Pc + re$Pc
        surf$R <- surf$R + 1
      }
    }
  }
  Ak <- prep$area * so_area_conv(prep)
  prevN <- pmax(Nt[, -TT, drop = FALSE], 1)
  out <- cbind(Nt, Nt / Ak, Ns, Bt[, -1, drop = FALSE], Bt[, -1, drop = FALSE] / prevN,
               Nt[, -1, drop = FALSE] / prevN, m[, "phi"], m[, "p0"], m[, "sigma"] * prep$cscale,
               pi * (sqrt(stats::qchisq(0.95, 2)) * m[, "sigma"])^2 * so_area_conv(prep),
               m[, paste0("gam", seq_len(TT)), drop = FALSE])
  colnames(out) <- c(nm$N, nm$D, "superpopulation_size", nm$B, nm$f, nm$lam, "survival", "baseline_detection",
                     "sigma", "home_range_area_95", nm$gam)
  for (k in seq_along(prep$dcov)) { out <- cbind(out, m[, paste0("bD", k)]); colnames(out)[ncol(out)] <- paste0("density_beta_", prep$dcov[k]) }
  attr(out, "surf") <- surf
  out
}

so_thin_idx <- function(nr, max_draws) if (nr <= max_draws) seq_len(nr) else unique(round(seq(1, nr, length.out = max_draws)))
so_pool <- function(chains, ndraw) {
  all <- do.call(rbind, chains)
  all[so_thin_idx(nrow(all), ndraw), , drop = FALSE]
}

# Realised density surfaces per period (expected animals per pixel), from a
# pooled set of parameter draws (works for every engine).
so_surfaces <- function(prep, draws, M) {
  surf <- list(real = matrix(0, prep$G, prep$TT), Pc = if (prep$n) matrix(0, prep$n, prep$G) else NULL, R = 0)
  dd <- draws[, setdiff(colnames(draws), c(paste0("N", seq_len(prep$TT)), paste0("B", seq_len(prep$TT)), "Nsuper")), drop = FALSE]
  r <- so_report_chain(prep, dd, seq_len(nrow(dd)), M, surf = surf)
  s <- attr(r, "surf")
  list(realised = s$real / s$R, Pc = if (!is.null(s$Pc)) s$Pc / s$R, R = s$R)
}

# Expected density (reported units) at standardised covariate xs, per draw,
# for period 1: E[N_1] = M gamma_1 spread over pixels by pi.
so_density_at <- function(prep, pd, xs, M) {
  PD <- ncol(prep$XD)
  vapply(seq_len(nrow(pd)), function(r) {
    b <- if (PD) pd[r, paste0("bD", seq_len(PD))] else numeric(0)
    lw <- if (PD) as.vector(prep$XD %*% b) else rep(0, prep$G)
    pd[r, "gam1"] * M * exp(sum(b * xs) - so_lse(lw)) / prep$pix^2 / so_area_conv(prep)
  }, 1)
}
