# popgen-ne engine: effective population size (LD and temporal methods),
# pairwise relatedness, heterozygosity and an optional bottleneck test.
#
# Genotypes everywhere are an integer array `geno[n, L, 2]` of allele codes
# (any integers; NA = missing). Rows are individuals, columns loci.
#
# Formulas follow the NeEstimator v2 documentation ("Ne Calculations",
# "On Burrows coefficients", "Implementation of Jackknife Process on Samples")
# and source code (Do et al. 2014), which implement:
#   LD method   : Waples (2006) bias correction, Waples & Do (2008) LDNe,
#                 Peel et al. (2013) missing-data weighting, jackknife over
#                 individuals of Jones, Ovenden & Wang (2016).
#   Temporal    : Nei & Tajima (1981) Fc, Pollak (1983) Fk, Waples (1989)
#                 sampling-plan corrections, Jorde & Ryman (2007) Fs'.
#   Relatedness : Queller & Goodnight (1989); Lynch & Ritland (1999).

# ---------------------------------------------------------------------------
# Allele dosage matrix
# ---------------------------------------------------------------------------

# Dosage (0/1/2) of every allele at every locus for individuals `idx`.
# Returns X (n x A, NA where the locus is untyped), loc (locus index of each
# column), allele code, freq (over typed individuals), typed (n x L logical).
dosage_matrix <- function(geno, idx = seq_len(dim(geno)[1]), loci = seq_len(dim(geno)[2])) {
  g1 <- geno[idx, loci, 1, drop = FALSE]; g2 <- geno[idx, loci, 2, drop = FALSE]
  dim(g1) <- dim(g2) <- c(length(idx), length(loci))
  typed <- !is.na(g1) & !is.na(g2)
  cols <- list(); loc <- integer(); al <- integer()
  for (j in seq_along(loci)) {
    a <- sort(unique(c(g1[typed[, j], j], g2[typed[, j], j])))
    if (!length(a)) next
    for (k in a) {
      v <- (g1[, j] == k) + (g2[, j] == k)
      v[!typed[, j]] <- NA
      cols[[length(cols) + 1]] <- v
    }
    loc <- c(loc, rep(j, length(a))); al <- c(al, a)
  }
  X <- if (length(cols)) do.call(cbind, cols) else matrix(numeric(), length(idx), 0)
  storage.mode(X) <- "double"
  nt <- colSums(typed)
  freq <- colSums(X, na.rm = TRUE) / (2 * nt[loc])
  list(X = X, loc = loc, allele = al, freq = freq, typed = typed, n_typed = nt, loci = loci)
}

# ---------------------------------------------------------------------------
# LD method (LDNe)
# ---------------------------------------------------------------------------

# Expected r^2 from sampling S individuals (Waples 2006; NeEstimator ExpR2Samp)
exp_r2_sample <- function(S) ifelse(S >= 30, 1 / S + 3.19 / S^2, 0.0018 + 0.907 / S + 4.44 / S^2)

# Bias-corrected Ne from r^2-drift (NeEstimator eq. 1.7 / LD_Ne()).
# Returns the INVERSE 1/Ne, which is continuous through r' = 0 (infinite Ne)
# and negative when r' < 0; Ne = 1 / inv.
ld_inv_ne <- function(rprime, S, mating = "random") {
  if (mating == "monogamy") {
    a <- ifelse(S >= 30, 2 / 3, 0.618); b <- ifelse(S >= 30, 7.2, 5.24)
  } else {
    a <- ifelse(S >= 30, 1 / 3, 0.308); b <- ifelse(S >= 30, 2.76, 2.08)
  }
  2 * rprime / (a + sqrt(pmax(a^2 - b * rprime, 0)))
}

inv_to_ne <- function(inv) ifelse(inv > 0, 1 / inv, Inf)

# Pcrit rule of Waples & Do (2010): 0.02 when S > 25, otherwise a value that
# removes singletons only (1/(2S) < Pcrit <= 1/S).
auto_pcrit <- function(S) if (S > 25) 0.02 else 1 / (2 * S - 1)

# Choose the allele columns used for LD at a given Pcrit.
ld_columns <- function(d, pcrit) {
  L <- length(d$loci)
  keep <- logical(length(d$loc)); n_ind <- integer(L); n_eff <- integer(L)
  for (j in seq_len(L)) {
    k <- which(d$loc == j)
    if (length(k) < 2) next
    f <- d$freq[k]
    if (pcrit > 0) {
      if (any(f > 1 - pcrit)) next          # locus nearly fixed
      el <- f >= pcrit - 1e-9
    } else el <- f > 0 & f < 1
    if (!any(el)) next
    dropped <- any(!el)
    ne <- sum(el)
    n_eff[j] <- ne
    n_ind[j] <- if (dropped) ne else ne - 1
    if (n_ind[j] < 1) next
    kk <- k[el]
    if (!dropped && ne == 2) kk <- kk[1]    # biallelic: one column carries all r^2
    keep[kk] <- TRUE
  }
  list(cols = which(keep), n_ind = n_ind, n_eff = n_eff)
}

# Pairwise r^2 (Burrows composite, S/(S-1) adjusted) between allele columns,
# aggregated to locus pairs, for the full sample and (optionally) every
# leave-one-individual-out sample.
# Returns locus-level matrices r2 (mean over allele pairs) and S (sample size
# at both loci), plus jackknife locus-pair arrays summarised to overall R2.
ld_core <- function(X, loc, n_ind, jackknife = TRUE, mating = "random", pcrit = 0, dropped = NULL) {
  n <- nrow(X); A <- ncol(X)
  ul <- sort(unique(loc)); L <- length(ul)
  li <- match(loc, ul)
  Gm <- matrix(0, A, L); Gm[cbind(seq_len(A), li)] <- 1
  ncol_l <- colSums(Gm)
  M <- !is.na(X); X0 <- X; X0[!M] <- 0
  Mi <- M * 1
  # locus-level typing mask (all columns of a locus share it)
  ML <- Mi[, match(seq_len(L), li), drop = FALSE]
  complete <- all(M)
  same_locus <- outer(li, li, "==")
  pair_r2 <- function(N, Sa, Sb, Saa, Sbb, Sab) {
    cv <- Sab / N - (Sa / N) * (Sb / N)
    va <- Saa / N - (Sa / N)^2; vb <- Sbb / N - (Sb / N)^2
    r2 <- cv^2 / (va * vb) * (N / (N - 1))^2
    r2[!(va > 1e-10 & vb > 1e-10) | N < 2] <- 0
    r2[r2 > 1] <- 1
    r2[same_locus] <- 0
    r2
  }
  agg <- if (A == L && all(li == seq_len(L))) function(r2) r2 else
    function(r2) (crossprod(Gm, r2 %*% Gm)) / outer(ncol_l, ncol_l)
  if (complete) {
    N <- n
    s <- colSums(X0); q <- colSums(X0^2); P <- crossprod(X0)
    r2 <- pair_r2(N, outer(s, rep(1, A)), outer(rep(1, A), s), outer(q, rep(1, A)), outer(rep(1, A), q), P)
    SL <- matrix(n, L, L)
  } else {
    N <- crossprod(Mi); Sa <- crossprod(X0, Mi); Saa <- crossprod(X0^2, Mi); P <- crossprod(X0)
    r2 <- pair_r2(N, Sa, t(Sa), Saa, t(Saa), P)
    SL <- crossprod(ML)
  }
  r2L <- agg(r2)
  out <- list(r2 = r2L, S = SL, n_ind = n_ind[ul], missing = !complete, loci = ul)
  if (jackknife && n > 3) {
    # Leave-one-individual-out r^2 (NeEstimator "jackknife on samples"):
    # allele columns stay those eligible in the full sample, but a column
    # that becomes rare (< Pcrit), fixed, or absent without individual k is
    # dropped for that replicate and the number of independent alleles of
    # its locus is recounted (NeEstimator SampJack notes, eq. 10).
    ut <- upper.tri(r2L)
    cnt <- colSums(X0); ntp <- colSums(Mi)
    if (is.null(dropped)) dropped <- rep(FALSE, L)
    single <- ncol_l == 1
    rk <- numeric(n); ok <- logical(n)
    for (k in seq_len(n)) {
      xk <- X0[k, ]; mk <- Mi[k, ]
      fk <- (cnt - xk) / (2 * (ntp - mk))
      el <- if (pcrit > 0) fk >= pcrit - 1e-9 & fk <= 1 - pcrit + 1e-9 else fk > 0 & fk < 1
      el[!is.finite(fk)] <- FALSE
      ne_l <- as.vector(crossprod(Gm, el * 1))
      fs_l <- as.vector(crossprod(Gm, fk * el))
      nind_k <- ifelse(single, ne_l, ifelse(!dropped & fs_l > 1 - 1e-6, ne_l - 1, ne_l))
      if (complete) {
        Nk <- n - 1
        sk <- s - xk; qk <- q - xk^2
        r2k <- pair_r2(Nk, outer(sk, rep(1, A)), outer(rep(1, A), sk),
                       outer(qk, rep(1, A)), outer(rep(1, A), qk), P - outer(xk, xk))
        Sk <- matrix(n - 1, L, L)
      } else {
        Nk <- N - outer(mk, mk); Sak <- Sa - outer(xk, mk); Saak <- Saa - outer(xk^2, mk)
        r2k <- pair_r2(Nk, Sak, t(Sak), Saak, t(Saak), P - outer(xk, xk))
        mlk <- ML[k, ]; Sk <- SL - outer(mlk, mlk)
      }
      if (A == L && all(li == seq_len(L))) {
        rL <- r2k[ut]
      } else {
        sums <- crossprod(Gm, (r2k * outer(el, el)) %*% Gm)
        rL <- (sums / outer(ne_l, ne_l))[ut]
      }
      w <- outer(nind_k, nind_k)[ut]
      if (!complete) w <- w * Sk[ut]^2
      keep <- w > 0 & is.finite(rL)
      ok[k] <- any(keep)
      rk[k] <- sum(w[keep] * rL[keep]) / sum(w[keep])
    }
    out$r2_jack <- rk[ok]
  }
  out
}

# Full LDNe estimate for one sample.
ld_ne <- function(geno, idx, pcrit = "auto", jackknife = TRUE, mating = "random",
                  max_loci = 3000, seed = 1, ci = 0.95) {
  d <- dosage_matrix(geno, idx)
  S0 <- length(idx)
  if (identical(pcrit, "auto")) pcrit <- auto_pcrit(S0)
  pcrit <- as.numeric(pcrit)
  sel <- ld_columns(d, pcrit)
  usable <- which(sel$n_ind > 0)
  note <- character()
  if (length(usable) > max_loci) {
    set.seed(seed)
    sub <- sort(sample(usable, max_loci))
    note <- sprintf("LD method used a random subset of %d of %d usable loci (max_loci).", max_loci, length(usable))
    sel$cols <- sel$cols[d$loc[sel$cols] %in% sub]
  }
  na_row <- function(msg) list(ne = NA_real_, inv = NA_real_, inv_lo = NA_real_, inv_hi = NA_real_,
                               ne_lo = NA_real_, ne_hi = NA_real_, S = S0, n_loci = 0, n_comp = 0,
                               r2 = NA_real_, er2 = NA_real_, pcrit = pcrit, note = msg)
  if (length(unique(d$loc[sel$cols])) < 2) return(na_row("fewer than 2 polymorphic loci pass Pcrit"))
  X <- d$X[, sel$cols, drop = FALSE]
  ul <- sort(unique(d$loc[sel$cols]))
  core <- ld_core(X, d$loc[sel$cols], sel$n_ind, jackknife = jackknife, mating = mating, pcrit = pcrit,
                  dropped = sel$n_ind[ul] == sel$n_eff[ul])
  ut <- upper.tri(core$r2)
  r2 <- core$r2[ut]; S <- core$S[ut]; nij <- outer(core$n_ind, core$n_ind)[ut]
  good <- S >= 2
  r2 <- r2[good]; S <- S[good]; nij <- nij[good]
  Ncomp <- sum(nij)
  Sh <- Ncomp / sum(nij / S)                     # weighted harmonic S (eq. 1.2)
  w <- if (core$missing) nij * S^2 else nij       # eq. 1.4
  R2 <- sum(w * r2) / sum(w)                      # eq. 1.5
  Eij <- exp_r2_sample(S)
  ER2 <- if (core$missing) sum(w * Eij) / sum(w) else exp_r2_sample(Sh)
  rprime <- R2 - ER2
  if (core$missing) {                             # eqs 1.10-1.12
    ne0 <- inv_to_ne(ld_inv_ne(rprime, Sh, mating))
    if (is.finite(ne0) && ne0 > 0) {
      w2 <- nij * S^2 / (S + 3 * ne0)^2
      rprime <- sum(w2 * (r2 - Eij)) / sum(w2)
      ER2 <- sum(w2 * Eij) / sum(w2)
      R2 <- rprime + ER2
    }
  }
  inv <- ld_inv_ne(rprime, Sh, mating)
  a <- (1 - ci) / 2
  # Parametric CI: chi-square on R2 with df = number of independent comparisons
  df_par <- Ncomp
  par_lo <- df_par * R2 / stats::qchisq(1 - a, df_par); par_hi <- df_par * R2 / stats::qchisq(a, df_par)
  out <- list(R2 = R2, er2 = ER2, rprime = rprime, inv = inv, ne = inv_to_ne(inv), S = Sh,
              n_loci = length(core$loci), n_comp = Ncomp, pcrit = pcrit,
              inv_par = sort(ld_inv_ne(c(par_lo, par_hi) - ER2, Sh, mating)))
  if (jackknife && length(core$r2_jack) > 3) {
    rk <- core$r2_jack; nj <- length(rk)
    vj <- (nj - 1) / nj * sum((rk - mean(rk))^2) * 0.84^2   # Jones et al. 2016 correction
    phi <- vj / mean(rk)^2
    df_j <- max(1, round(2 / phi))
    lo <- df_j * R2 / stats::qchisq(1 - a, df_j); hi <- df_j * R2 / stats::qchisq(a, df_j)
    out$inv_lo <- ld_inv_ne(lo - ER2, Sh, mating); out$inv_hi <- ld_inv_ne(hi - ER2, Sh, mating)
    out$df_jack <- df_j
    out$ci_type <- "jackknife over individuals"
  } else {
    out$inv_lo <- out$inv_par[1]; out$inv_hi <- out$inv_par[2]
    out$ci_type <- "parametric chi-square"
  }
  out$ne_lo <- inv_to_ne(out$inv_hi); out$ne_hi <- inv_to_ne(out$inv_lo)
  out$note <- note
  out
}

# ---------------------------------------------------------------------------
# Temporal method
# ---------------------------------------------------------------------------

# r2_drift: drift LD (R^2' of the LD method, ~1/(3Ne)) used to inflate the
# variance for correlation of frequency change among unlinked loci in LD:
# DEFF = 1 + deff_c * (D - D/L) * r2_drift (calibrated by simulation, VALIDATION.md).
temporal_ne <- function(geno, idx1, idx2, t, census_N = NA, pcrit = "auto", ci = 0.95,
                        r2_drift = 0, deff_c = 0.5) {
  L <- dim(geno)[2]
  if (identical(pcrit, "auto")) pcrit <- auto_pcrit(min(length(idx1), length(idx2)))
  q <- as.numeric(pcrit)
  plan1 <- is.finite(census_N) && census_N > 0
  ext <- if (plan1) 1 / census_N else 0
  rows <- vector("list", L)
  for (j in seq_len(L)) {
    a1 <- geno[idx1, j, ]; a2 <- geno[idx2, j, ]
    t1 <- stats::complete.cases(a1); t2 <- stats::complete.cases(a2)
    m <- sum(t1); nn <- sum(t2)
    if (m < 1 || nn < 1) next
    al <- sort(unique(c(a1[t1, ], a2[t2, ])))
    x <- vapply(al, function(k) sum(a1[t1, ] == k), 0) / (2 * m)
    y <- vapply(al, function(k) sum(a2[t2, ] == k), 0) / (2 * nn)
    zb <- (m * x + nn * y) / (m + nn)
    if (q > 0) {
      if (any(zb >= 1 - q)) next                        # nearly monomorphic: dropped
      rare <- zb < q
      if (any(rare)) { x <- c(x[!rare], sum(x[rare])); y <- c(y[!rare], sum(y[rare])) }
    }
    z <- (x + y) / 2
    if (length(x) < 2 || all(z %in% c(0, 1))) next
    K <- length(x)
    rows[[j]] <- c(K = K, m = m, n = nn, s = 2 / (1 / m + 1 / nn),
                   fc = sum((x - y)^2 / (z - x * y)), fk = sum((x - y)^2 / z),
                   top = sum((x - y)^2), bot = sum(z * (1 - z)))
  }
  r <- do.call(rbind, rows)
  if (is.null(r) || nrow(r) < 2) return(NULL)
  r <- as.data.frame(r)
  a <- (1 - ci) / 2
  df_ind <- sum(r$K - 1)
  deff <- 1 + deff_c * (df_ind - df_ind / nrow(r)) * max(0, r2_drift, na.rm = TRUE)
  df_eff <- max(1, df_ind / deff)
  chi <- function(Fv) c(df_eff * Fv / stats::qchisq(1 - a, df_eff), df_eff * Fv / stats::qchisq(a, df_eff))
  invS_K <- sum(r$K / r$s) / sum(r$K); invS_K1 <- sum((r$K - 1) / r$s) / sum(r$K - 1)
  # Nei & Tajima Fc
  Fc <- sum(r$fc) / sum(r$K); Fcp <- Fc - invS_K + ext
  # Pollak Fk
  Fk <- sum(r$fk) / sum(r$K - 1); Fkp <- Fk - invS_K1 + ext
  # Jorde & Ryman Fs'
  fs_prime <- function(top, bot, invh, inv2) {
    Fs <- top / bot
    (Fs * (1 - invh / 4 + ext / 4) - invh + ext) / ((1 + Fs / 4) * (1 - inv2))
  }
  invh_all <- mean(1 / r$s); inv2_all <- mean(1 / (2 * r$n))
  Fsp <- fs_prime(sum(r$top), sum(r$bot), invh_all, inv2_all)
  nl <- nrow(r)
  jk <- vapply(seq_len(nl), function(l)
    fs_prime(sum(r$top[-l]), sum(r$bot[-l]), mean(1 / r$s[-l]), mean(1 / (2 * r$n[-l]))), 0)
  fbar <- mean(jk); se <- sqrt(deff * (nl - 1) / nl * sum((jk - fbar)^2))
  z <- stats::qnorm(1 - a)
  # chi-square interval on the raw Fs with df = 2/CV^2 from the jackknife
  # over loci (x DEFF), bounds mapped to Fs' (monotone): keeps the interval
  # right-skewed like the sampling distribution of Fs.
  Fs_all <- sum(r$top) / sum(r$bot)
  jkF <- (sum(r$top) - r$top) / (sum(r$bot) - r$bot)
  vF <- deff * (nl - 1) / nl * sum((jkF - mean(jkF))^2)
  df_F <- if (vF > 0) max(1, 2 * Fs_all^2 / vF) else 1e6
  map_F <- function(Fs) (Fs * (1 - invh_all / 4 + ext / 4) - invh_all + ext) / ((1 + Fs / 4) * (1 - inv2_all))
  Fs_ci <- c(df_F * Fs_all / stats::qchisq(1 - a, df_F), df_F * Fs_all / stats::qchisq(a, df_F))
  list(
    JR = list(F = sum(r$top) / sum(r$bot), Fp = Fsp, Fp_lo = map_F(Fs_ci[1]), Fp_hi = map_F(Fs_ci[2]),
              ci_type = "jackknife over loci (chi-square)", df_chisq = df_F,
              Fp_lo_normal = fbar - z * se, Fp_hi_normal = fbar + z * se),   # NeEstimator's interval
    NT = list(F = Fc, Fp = Fcp, Fp_lo = chi(Fc)[1] - invS_K + ext, Fp_hi = chi(Fc)[2] - invS_K + ext,
              ci_type = "parametric chi-square"),
    Pollak = list(F = Fk, Fp = Fkp, Fp_lo = chi(Fk)[1] - invS_K1 + ext, Fp_hi = chi(Fk)[2] - invS_K1 + ext,
                  ci_type = "parametric chi-square"),
    t = t, S = 1 / invS_K, n_loci = nl, n_alleles_indep = df_ind, deff = deff, plan = if (plan1) "I" else "II",
    census_N = census_N, pcrit = q)
}

# F' -> 1/(2Ne) per generation = F'/t ; Ne = t/(2F')
temporal_rows <- function(tm, method, t) {
  m <- tm[[method]]
  c(het = m$Fp / t, het_lo = m$Fp_lo / t, het_hi = m$Fp_hi / t)
}

# ---------------------------------------------------------------------------
# Relatedness (moment estimators)
# ---------------------------------------------------------------------------

# Allele-frequency reference for a group of individuals.
rel_prep <- function(geno, idx) {
  d <- dosage_matrix(geno, idx)
  p <- d$freq
  p[p <= 0] <- NA
  list(X = d$X, loc = d$loc, p = p, L = length(d$loci))
}

# All pairwise estimates for the group (matrix algebra; see engine notes in
# README). Returns symmetric matrices for "QG" and "LR".
rel_matrix <- function(rp) {
  X <- rp$X; M <- !is.na(X); X0 <- X; X0[!M] <- 0
  p <- rp$p; p[is.na(p)] <- 1
  A <- ncol(X); Lm <- matrix(0, A, rp$L); Lm[cbind(seq_len(A), rp$loc)] <- 1
  ML <- (M * 1) %*% Lm; ML <- (ML > 0) * 1         # n x L typed loci
  Mc <- ML %*% t(Lm)                                # n x A: allele's locus typed
  Lxy <- ML %*% t(ML)                               # loci typed in both
  # Queller & Goodnight 1989 (ratio of sums over loci), x as reference
  numQ <- 0.5 * tcrossprod(X0) - (X0 %*% (p * t(Mc)))
  denQ <- (0.5 * X0^2 - sweep(X0, 2, p, "*")) %*% t(Mc)
  rQ <- numQ / denQ
  QG <- (rQ + t(rQ)) / 2
  # Lynch & Ritland 1999 (weighted per-locus, eq. 5a with weights 7a)
  WR <- 0.5 * (X0 %*% (t(X0) / p)) - 2 * Lxy
  W <- (sweep(0.5 * X0^2, 2, p, "/")) %*% t(Mc) - 2 * Lxy
  rL <- WR / W
  LR <- (rL + t(rL)) / 2
  diag(QG) <- diag(LR) <- NA
  list(QG = QG, LR = LR)
}

# Per-locus contributions for a set of pairs (rows i,j), used for the locus
# jackknife of class means and for simulated pairs. Ga, Gb: pairs x A dosage.
rel_pair_terms <- function(Ga, Gb, loc, p, L) {
  Ma <- !is.na(Ga); Mb <- !is.na(Gb)
  A <- ncol(Ga); Lm <- matrix(0, A, L); Lm[cbind(seq_len(A), loc)] <- 1
  both <- ((Ma * 1) %*% Lm > 0) & ((Mb * 1) %*% Lm > 0)   # pairs x L
  Ga[!Ma] <- 0; Gb[!Mb] <- 0
  P <- matrix(p, nrow(Ga), A, byrow = TRUE)
  sh <- (Ga * Gb) %*% Lm
  term <- function(Gx) list(
    qn = 0.5 * sh - (Gx * P) %*% Lm, qd = (0.5 * Gx^2 - Gx * P) %*% Lm,
    lw = (0.5 * Gx^2 / P) %*% Lm - 2)
  lr_num <- 0.5 * ((Ga * Gb / P) %*% Lm) - 2
  ta <- term(Ga); tb <- term(Gb)
  z <- function(m) { m[!both] <- 0; m }
  list(qn_a = z(ta$qn), qn_b = z(tb$qn), qd_a = z(ta$qd), qd_b = z(tb$qd), lr = z(lr_num),
       lw_a = z(ta$lw), lw_b = z(tb$lw))
}

rel_from_terms <- function(tt) {
  s <- rowSums
  cbind(QG = (s(tt$qn_a) / s(tt$qd_a) + s(tt$qn_b) / s(tt$qd_b)) / 2,
        LR = (s(tt$lr) / s(tt$lw_a) + s(tt$lr) / s(tt$lw_b)) / 2)
}

# Leave-one-locus-out class means (pairs x L sums minus each locus' column).
rel_jack_terms <- function(tt) {
  loo <- function(m) rowSums(m) - m          # pairs x L: sum without locus l
  QG <- (loo(tt$qn_a) / loo(tt$qd_a) + loo(tt$qn_b) / loo(tt$qd_b)) / 2
  lr <- loo(tt$lr)
  LR <- (lr / loo(tt$lw_a) + lr / loo(tt$lw_b)) / 2
  cbind(QG = colMeans(QG, na.rm = TRUE), LR = colMeans(LR, na.rm = TRUE))
}

# Mean relatedness of a class of pairs with locus-jackknife 95% CI.
rel_class_mean <- function(rp, i, j, ci = 0.95, se_method = "max") {
  Ga <- rp$X[i, , drop = FALSE]; Gb <- rp$X[j, , drop = FALSE]
  p <- rp$p; p[is.na(p)] <- 1
  tt <- rel_pair_terms(Ga, Gb, rp$loc, p, rp$L)
  full <- colMeans(rel_from_terms(tt), na.rm = TRUE)
  L <- rp$L
  jk <- rel_jack_terms(tt)
  se_loc <- sqrt((L - 1) / L * colSums(sweep(jk, 2, colMeans(jk))^2))
  # jackknife over individuals: drop every pair containing individual k
  pr <- rel_from_terms(tt)
  u <- unique(c(i, j))
  if (length(u) > 2) {
    jki <- t(vapply(u, function(k) { keep <- i != k & j != k
      if (!any(keep)) return(c(NA_real_, NA_real_)); colMeans(pr[keep, , drop = FALSE], na.rm = TRUE) }, numeric(2)))
    jki <- jki[stats::complete.cases(jki), , drop = FALSE]
    nu <- nrow(jki)
    se_ind <- sqrt((nu - 1) / nu * colSums(sweep(jki, 2, colMeans(jki))^2))
  } else se_ind <- c(NA, NA)
  se <- if (identical(se_method, "loci")) se_loc else if (identical(se_method, "individuals")) se_ind else pmax(se_loc, se_ind, na.rm = TRUE)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  data.frame(estimator = names(full), estimate = full, sd = se, lower = full - z * se, upper = full + z * se,
             n_pairs = length(i), se_loci = se_loc, se_ind = se_ind, row.names = NULL)
}

# Simulate genotype pairs of known relationship from allele frequencies and
# estimate their relatedness: misclassification rates at the user's markers.
rel_sim_pairs <- function(rp, n_pairs = 200, seed = 1) {
  set.seed(seed)
  p <- rp$p; p[is.na(p)] <- 0
  L <- rp$L; A <- length(p)
  draw <- function(np) {   # np x L allele indices (column ids into p)
    sapply(seq_len(L), function(l) { k <- which(rp$loc == l); if (length(k) == 1) rep(k, np) else sample(k, np, TRUE, prob = p[k]) })
  }
  to_dos <- function(a1, a2) {
    G <- matrix(0, nrow(a1), A)
    for (l in seq_len(L)) { G[cbind(seq_len(nrow(a1)), a1[, l])] <- G[cbind(seq_len(nrow(a1)), a1[, l])] + 1
                            G[cbind(seq_len(nrow(a1)), a2[, l])] <- G[cbind(seq_len(nrow(a1)), a2[, l])] + 1 }
    G
  }
  ibd <- list(unrelated = c(1, 0, 0), halfsib = c(0.5, 0.5, 0), fullsib = c(0.25, 0.5, 0.25), parent_offspring = c(0, 1, 0))
  pp <- p; pp[pp <= 0] <- 1
  out <- lapply(names(ibd), function(cl) {
    np <- n_pairs
    a1 <- draw(np); a2 <- draw(np)
    b1 <- draw(np); b2 <- draw(np)
    k <- matrix(sample(0:2, np * L, TRUE, prob = ibd[[cl]]), np, L)
    b1[k >= 1] <- a1[k >= 1]; b2[k == 2] <- a2[k == 2]
    tt <- rel_pair_terms(to_dos(a1, a2), to_dos(b1, b2), rp$loc, pp, L)
    data.frame(class = cl, rel_from_terms(tt))
  })
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Heterozygosity and bottleneck (heterozygosity-excess) test
# ---------------------------------------------------------------------------

het_summary <- function(geno, idx, ci = 0.95) {
  L <- dim(geno)[2]
  st <- t(vapply(seq_len(L), function(j) {
    a <- geno[idx, j, , drop = FALSE]; dim(a) <- c(length(idx), 2)
    ok <- stats::complete.cases(a); n <- sum(ok)
    if (n < 2) return(c(NA, NA, NA, NA))
    al <- c(a[ok, 1], a[ok, 2]); f <- table(al) / (2 * n)
    c(Ho = mean(a[ok, 1] != a[ok, 2]), He = (2 * n / (2 * n - 1)) * (1 - sum(f^2)), k = length(f), n = n)
  }, numeric(4)))
  colnames(st) <- c("Ho", "He", "k", "n")
  st <- st[stats::complete.cases(st), , drop = FALSE]
  st <- st[st[, "k"] > 1, , drop = FALSE]
  Lu <- nrow(st)
  f <- function(s) c(Ho = mean(s[, "Ho"]), He = mean(s[, "He"]), FIS = 1 - mean(s[, "Ho"]) / mean(s[, "He"]))
  full <- f(st)
  jk <- t(vapply(seq_len(Lu), function(l) f(st[-l, , drop = FALSE]), numeric(3)))
  se <- sqrt((Lu - 1) / Lu * colSums(sweep(jk, 2, colMeans(jk))^2))
  z <- stats::qnorm(1 - (1 - ci) / 2)
  list(table = data.frame(stat = names(full), estimate = full, sd = se, lower = full - z * se, upper = full + z * se,
                          row.names = NULL), per_locus = st)
}

# Ewens sampling formula: distribution of He given k alleles in a sample of n
# genes under the infinite-alleles model at mutation-drift equilibrium
# (Cornuet & Luikart 1996). Returns the simulated He values (Hoppe urn,
# conditioned on k by rejection).
.heq_cache <- new.env(parent = emptyenv())
heq_iam <- function(n, k, nrep = 2000) {
  key <- paste(n, k, nrep)
  if (!is.null(.heq_cache[[key]])) return(.heq_cache[[key]])
  if (k <= 1) return(rep(0, nrep))
  th <- stats::uniroot(function(th) sum(th / (th + 0:(n - 1))) - k, c(1e-6, 1e4))$root
  he <- numeric(0); tries <- 0
  while (length(he) < nrep && tries < 200) {
    tries <- tries + 1
    for (r in seq_len(nrep)) {   # Hoppe urn
      cnt <- integer(0)
      for (i in 0:(n - 1)) {
        if (stats::runif(1) < th / (th + i)) cnt <- c(cnt, 1L)
        else { j <- sample.int(length(cnt), 1, prob = cnt); cnt[j] <- cnt[j] + 1L }
      }
      if (length(cnt) == k) he <- c(he, (n / (n - 1)) * (1 - sum((cnt / n)^2)))
      if (length(he) >= nrep) break
    }
  }
  .heq_cache[[key]] <- he
  he
}

# Heterozygosity-excess test. The P-value is a Monte Carlo test: the observed
# mean standardised difference (He - E[He|k]) / SD[He|k] over loci is compared
# with its distribution when every locus draws He from its own equilibrium
# distribution given k. This is calibrated by construction; the Wilcoxon
# signed-rank test used by BOTTLENECK (Piry et al. 1999) assumes symmetric
# differences, which the skewed He|k distribution violates (false-positive
# rate 0.18 at nominal 0.05 in our check), so it is reported only for reference.
bottleneck_test <- function(per_locus, nrep = 1000, seed = 1, nnull = 4999) {
  set.seed(seed)
  sims <- lapply(seq_len(nrow(per_locus)), function(l) heq_iam(2 * per_locus[l, "n"], per_locus[l, "k"], nrep))
  m <- vapply(sims, mean, 0); s <- vapply(sims, stats::sd, 0)
  ok <- is.finite(s) & s > 0
  d <- per_locus[, "He"] - m
  z <- (d / s)[ok]
  Tobs <- mean(z)
  Tnull <- rowMeans(vapply(which(ok), function(l) (sample(sims[[l]], nnull, TRUE) - m[l]) / s[l], numeric(nnull)))
  wt <- suppressWarnings(stats::wilcox.test(d, alternative = "greater"))
  list(n_loci = length(d), n_excess = sum(d > 0), expected_excess = length(d) * 0.5,
       p_value = (1 + sum(Tnull >= Tobs)) / (nnull + 1), p_wilcoxon = wt$p.value,
       mean_diff = mean(d), mean_std_diff = Tobs)
}
