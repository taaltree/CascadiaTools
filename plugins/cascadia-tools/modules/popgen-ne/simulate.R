# Wright-Fisher simulation with a known pedigree.
#
# An ideal population of N diploid individuals (random union of gametes,
# selfing allowed with probability 1/N, discrete generations, Poisson-like
# family sizes) has Ne = N. Founders are in Hardy-Weinberg and linkage
# equilibrium; loci are unlinked (each gamete takes either parental allele
# with probability 1/2). After `burn_in` generations (drift-LD equilibrium for
# unlinked loci is reached within a few generations) we take
#   - an "early" sample of S individuals (generation 0),
#   - a "late" sample of S individuals t generations later,
#   - a kin panel of known full-sib and half-sib pairs born in the late
#     generation (extra offspring of that generation's parents; not part of
#     the census and not used for Ne).
# Sampling is non-destructive (sampled individuals still reproduce), so the
# temporal method uses sampling plan I with census N (Waples 1989).
#
# Pedigree kinship is tracked exactly (recursion on the kinship matrix), so
# the truth for relatedness is the pedigree relatedness of each pair
# expressed relative to the sample's own allele frequencies (the reference
# that moment estimators use):
#   r*_ij = 2 (th_ij - th_i. - th_.j + th_..) / (1 - th_..),
# where th is kinship relative to the founders and dots are means over the
# sampled group (self-kinship included).

wf_simulate <- function(params, design, seed) {
  set.seed(seed)
  N <- as.integer(params$Ne)
  L <- as.integer(design$n_loci)
  S <- as.integer(design$n_ind)
  tgap <- as.integer(design$generations_between)
  burn <- as.integer(design$burn_in)
  marker <- design$marker
  # founders
  if (marker == "snp") {
    p <- stats::runif(L, design$snp_maf_min, 1 - design$snp_maf_min)
    H1 <- matrix(stats::rbinom(N * L, 1, rep(p, each = N)) + 1L, N, L)
    H2 <- matrix(stats::rbinom(N * L, 1, rep(p, each = N)) + 1L, N, L)
  } else {
    k <- as.integer(design$n_alleles)
    H1 <- H2 <- matrix(0L, N, L)
    for (l in seq_len(L)) {
      f <- stats::rgamma(k, design$allele_dirichlet); f <- f / sum(f)
      H1[, l] <- sample.int(k, N, TRUE, f) + 100L; H2[, l] <- sample.int(k, N, TRUE, f) + 100L
    }
  }
  storage.mode(H1) <- storage.mode(H2) <- "integer"
  need_kin <- isTRUE(design$n_fullsib_pairs > 0) || isTRUE(design$n_halfsib_pairs > 0) || isTRUE(design$n_unrelated_pairs > 0)
  phi <- if (need_kin) diag(0.5, N) else matrix(0, 1, 1)   # kinship among founders
  par_prev <- NULL
  gamete <- function(A1, A2, par) {
    u <- matrix(stats::runif(length(par) * L) < 0.5, length(par), L)
    out <- A1[par, , drop = FALSE]; out[u] <- A2[par, , drop = FALSE][u]; out
  }
  kin_full <- function(phi, a, b) {
    m <- 0.25 * (phi[a, a, drop = FALSE] + phi[a, b, drop = FALSE] + phi[b, a, drop = FALSE] + phi[b, b, drop = FALSE])
    diag(m) <- 0.5 * (1 + phi[cbind(a, b)])
    m
  }
  # track kinship over the last `kin_gens` generations only when N is large
  ngen <- burn + tgap
  kin_from <- if (N > 600) max(1, ngen - 12) else 1
  early <- NULL
  for (g in seq_len(ngen)) {
    a <- sample.int(N, N, TRUE); b <- sample.int(N, N, TRUE)
    n1 <- gamete(H1, H2, a); n2 <- gamete(H1, H2, b)
    if (need_kin && g >= kin_from) phi <- kin_full(phi, a, b) else if (need_kin) phi <- diag(0.5, N)
    H1_old <- H1; H2_old <- H2
    H1 <- n1; H2 <- n2
    if (g == burn) {
      ie <- sort(sample.int(N, S))
      early <- list(H1 = H1[ie, , drop = FALSE], H2 = H2[ie, , drop = FALSE])
    }
    if (g == ngen) {
      # the parents' generation is (H1_old, H2_old) with kinship phi_par
      last <- list(a = a, b = b, Hp1 = H1_old, Hp2 = H2_old, gp_prev = par_prev)
    }
    par_prev <- list(a = a, b = b)
    if (g == ngen - 1) phi_par <- phi
  }
  if (!need_kin) phi_par <- NULL else if (ngen == 1) phi_par <- diag(0.5, N)
  il <- sort(sample.int(N, S))
  # kin panel from the parents of the late generation
  nfs <- as.integer(design$n_fullsib_pairs); nhs <- as.integer(design$n_halfsib_pairs)
  pa <- integer(); pb <- integer(); kin <- NULL
  if (nfs > 0) for (f in seq_len(nfs)) {
    pr <- sample.int(N, 2); pa <- c(pa, pr[1], pr[1]); pb <- c(pb, pr[2], pr[2])
  }
  if (nhs > 0) for (f in seq_len(nhs)) {
    pr <- sample.int(N, 3); pa <- c(pa, pr[1], pr[1]); pb <- c(pb, pr[2], pr[3])
  }
  npan <- length(pa)
  if (npan) {
    k1 <- gamete(last$Hp1, last$Hp2, pa); k2 <- gamete(last$Hp1, last$Hp2, pb)
  }
  # kinship of the relatedness group (late sample + panel) from parental kinship
  ga <- c(last$a[il], pa); gb <- c(last$b[il], pb)
  K <- if (need_kin) kin_full(phi_par, ga, gb) else diag(0.5, length(ga))
  m_i <- rowMeans(K); m_all <- mean(K)
  Rtrue <- 2 * (K - outer(m_i, m_i, "+") + m_all) / (1 - m_all)
  th_s <- mean(K[upper.tri(K)])
  Rwg <- 2 * (K - th_s) / (1 - th_s)          # Weir & Goudet (2017) beta-style reference
  Rraw <- 2 * K                               # relative to the founders
  alt <- list()
  nl <- length(il)
  ids_late <- sprintf("L%03d", seq_len(nl)); ids_pan <- if (npan) sprintf("K%03d", seq_len(npan)) else character()
  ids_early <- sprintf("E%03d", seq_len(S))
  pairs <- NULL; tr <- c()
  if (nfs > 0) {
    i <- nl + seq(1, 2 * nfs, by = 2); j <- i + 1
    pairs <- rbind(pairs, data.frame(id1 = c(ids_late, ids_pan)[i], id2 = c(ids_late, ids_pan)[j], relationship = "fullsib"))
    tr["fullsib"] <- mean(Rtrue[cbind(i, j)])
    alt$fullsib <- c(wg = mean(Rwg[cbind(i, j)]), raw = mean(Rraw[cbind(i, j)]))
  }
  if (nhs > 0) {
    i <- nl + 2 * nfs + seq(1, 2 * nhs, by = 2); j <- i + 1
    pairs <- rbind(pairs, data.frame(id1 = c(ids_late, ids_pan)[i], id2 = c(ids_late, ids_pan)[j], relationship = "halfsib"))
    tr["halfsib"] <- mean(Rtrue[cbind(i, j)])
    alt$halfsib <- c(wg = mean(Rwg[cbind(i, j)]), raw = mean(Rraw[cbind(i, j)]))
  }
  # unrelated: late-sample pairs sharing no parent or grandparent
  if (isTRUE(design$n_unrelated_pairs > 0)) {
    pl <- cbind(last$a[il], last$b[il])
    gpar <- function(x) if (is.null(last$gp_prev)) c(x, x) else c(last$gp_prev$a[x], last$gp_prev$b[x])
    gpl <- t(apply(pl, 1, function(r) c(gpar(r[1]), gpar(r[2]))))
    cand <- which(upper.tri(diag(nl)), arr.ind = TRUE)
    okp <- apply(cand, 1, function(r) !any(pl[r[1], ] %in% pl[r[2], ]) && !any(gpl[r[1], ] %in% gpl[r[2], ]))
    cand <- cand[okp, , drop = FALSE]
    take <- cand[sample.int(nrow(cand), min(nrow(cand), design$n_unrelated_pairs)), , drop = FALSE]
    pairs <- rbind(pairs, data.frame(id1 = ids_late[take[, 1]], id2 = ids_late[take[, 2]], relationship = "unrelated"))
    tr["unrelated"] <- mean(Rtrue[take])
    alt$unrelated <- c(wg = mean(Rwg[take]), raw = mean(Rraw[take]))
  }
  # assemble genotype array
  late1 <- H1[il, , drop = FALSE]; late2 <- H2[il, , drop = FALSE]
  G1 <- rbind(early$H1, late1, if (npan) k1); G2 <- rbind(early$H2, late2, if (npan) k2)
  n <- nrow(G1)
  geno <- array(NA_integer_, c(n, L, 2)); geno[, , 1] <- G1; geno[, , 2] <- G2
  if (isTRUE(design$missing_rate > 0)) {
    miss <- matrix(stats::runif(n * L) < design$missing_rate, n, L)
    g1 <- geno[, , 1]; g2 <- geno[, , 2]; g1[miss] <- NA; g2[miss] <- NA
    geno[, , 1] <- g1; geno[, , 2] <- g2
  }
  loci <- sprintf("%s%04d", if (marker == "snp") "snp" else "msat", seq_len(L))
  dimnames(geno) <- list(c(ids_early, ids_late, ids_pan), loci, NULL)
  ind <- data.frame(id = c(ids_early, ids_late, ids_pan), pop = "sim",
                    sample = c(rep("early", S), rep("late", nl), rep("late", npan)),
                    time = c(rep(0, S), rep(tgap, nl + npan)),
                    use_for_ne = c(rep(TRUE, S + nl), rep(isTRUE(design$panel_in_ne), npan)), stringsAsFactors = FALSE)
  list(geno = geno, ind = ind, known_kin = pairs,
       census_N = if (isFALSE(design$census_known)) NULL else c(sim = N), generation_time = 1,
       marker = marker, rel_truth = tr, rel_truth_alt = alt)
}
