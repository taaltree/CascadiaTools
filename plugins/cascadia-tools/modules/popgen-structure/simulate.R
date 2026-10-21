# Wright-Fisher simulator for subdivided populations (island, 1-D and 2-D
# stepping-stone, panmixia) and the EXACT expected F-statistics for the same
# model from an identity-by-descent recursion. The recursion gives the truth
# that simulation recovery compares the estimators with.
#
# Life cycle per generation (all demes, non-overlapping generations):
#   each offspring in deme i is selfed with probability `selfing`;
#   each parent slot independently comes from deme k with probability M[i,k]
#   (backward migration matrix) and is a uniformly chosen adult of deme k
#   (random union of gametes; a selfed offspring uses one parent twice);
#   each gene is a random one of the parent's two genes (unlinked loci);
#   each gene mutates with probability mu to an allele drawn uniformly from
#   the K possible alleles (K-allele, parent-independent mutation).
# Founders: per locus an ancestral frequency vector p ~ Dirichlet(alpha,..);
# every founder gene is an independent draw from p (no initial structure).
#
# Identity recursion (Malecot / Rousset 1996 style, exact for this model):
#   F_i  = P(two genes within an individual of deme i are identical in state)
#   Q_ij = P(two genes in distinct individuals of demes i, j are identical)
#   G_k  = (1/N_k)(1+F_k)/2 + (1-1/N_k) Q_kk      (two independent parent draws in deme k)
#   H    = Q with diagonal replaced by G
#   Q'   = M H M^T
#   F'_i = s * sum_k M_ik (1+F_k)/2 + (1-s) Q'_ii
#   X''  = c X' + (1-c)/K,  c = (1-mu)^2          (mutation, for F and Q)
# Start: F = Q = E[sum p^2] = (alpha+1)/(K alpha+1).
# Estimand definitions (what WC84 / Nei estimators target):
#   Fst (pair i,j) = ((Q_ii+Q_jj)/2 - Q_ij) / (1 - Q_ij)
#   Fst (global)   = (Qw - Qb) / (1 - Qb), Qw = mean_i Q_ii, Qb = mean_{i<j} Q_ij
#   Fis (global)   = (Fbar - Qw) / (1 - Qw);  He_i = 1 - Q_ii;  Ho_i = 1 - F_i

ps_lattice <- function(model, n_pops, spacing) {
  if (model %in% c("stepping_stone_2d", "island", "panmixia")) {
    nc <- ceiling(sqrt(n_pops))
    xy <- data.frame(x = ((seq_len(n_pops) - 1) %% nc) * spacing,
                     y = ((seq_len(n_pops) - 1) %/% nc) * spacing)
  } else xy <- data.frame(x = (seq_len(n_pops) - 1) * spacing, y = 0)
  xy
}

ps_migration_matrix <- function(model, n_pops, m, xy = NULL, spacing = 1) {
  d <- n_pops
  M <- matrix(0, d, d)
  if (model == "panmixia") return(matrix(1 / d, d, d))
  if (model == "island") {
    M[] <- m / (d - 1); diag(M) <- 1 - m
    return(M)
  }
  if (is.null(xy)) xy <- ps_lattice(model, d, spacing)
  D <- as.matrix(stats::dist(xy))
  nb <- abs(D - spacing) < 1e-8 * spacing
  nmax <- if (model == "stepping_stone_2d") 4 else 2
  for (i in seq_len(d)) {
    k <- which(nb[i, ])
    M[i, k] <- m / nmax
    M[i, i] <- 1 - length(k) * m / nmax   # edge demes keep the unused share (reflecting boundary)
  }
  M
}

# Exact expected F, Q after T generations.
ps_identity_recursion <- function(M, N, selfing, mu, K, alpha, T_gen) {
  d <- nrow(M)
  N <- rep_len(N, d)
  q0 <- (alpha + 1) / (K * alpha + 1)
  Fv <- rep(q0, d); Q <- matrix(q0, d, d)
  cmu <- (1 - mu)^2
  for (t in seq_len(T_gen)) {
    G <- (1 / N) * (1 + Fv) / 2 + (1 - 1 / N) * diag(Q)
    H <- Q; diag(H) <- G
    Qn <- M %*% H %*% t(M)
    Fn <- selfing * as.vector(M %*% ((1 + Fv) / 2)) + (1 - selfing) * diag(Qn)
    Q <- cmu * Qn + (1 - cmu) / K
    Fv <- cmu * Fn + (1 - cmu) / K
  }
  list(F = Fv, Q = Q)
}

ps_truth_from_identity <- function(id, pop_names, xy = NULL, ibd_distance = "log") {
  Q <- id$Q; Fv <- id$F; d <- nrow(Q)
  Qw <- mean(diag(Q)); Qb <- mean(Q[upper.tri(Q)])
  out <- c(Fst_global = (Qw - Qb) / (1 - Qb),
           Fis_global = (mean(Fv) - Qw) / (1 - Qw),
           He_mean = 1 - Qw, Ho_mean = 1 - mean(Fv))
  pw <- outer(diag(Q), diag(Q), "+") / 2
  Fp <- (pw - Q) / (1 - Q)
  pr <- which(upper.tri(Q), arr.ind = TRUE)
  fst_pairs <- stats::setNames(Fp[pr], paste0("Fst_", pop_names[pr[, 1]], "_vs_", pop_names[pr[, 2]]))
  slope <- NA_real_
  if (!is.null(xy)) {
    D <- as.matrix(stats::dist(xy))[pr]
    x <- if (ibd_distance == "log") log(D) else D
    y <- Fp[pr] / (1 - Fp[pr])
    slope <- if (stats::var(x) > 0) stats::cov(x, y) / stats::var(x) else 0
  }
  list(global = out, pairs = fst_pairs, ibd_slope = slope,
       He_pop = stats::setNames(1 - diag(Q), pop_names),
       Fis_pop = stats::setNames((Fv - diag(Q)) / (1 - diag(Q)), pop_names), Fpair = Fp)
}

# Forward Wright-Fisher simulation. Returns allele matrices (individuals x
# loci) for the whole population in the final generation, with deme labels.
ps_wf_forward <- function(M, N, selfing, mu, K, alpha, T_gen, n_loci) {
  d <- nrow(M); N <- rep_len(N, d)
  deme <- rep(seq_len(d), N); ntot <- sum(N)
  start <- c(0, cumsum(N))[seq_len(d)]
  # founders
  A1 <- A2 <- matrix(0L, ntot, n_loci)
  for (l in seq_len(n_loci)) {
    g <- stats::rgamma(K, alpha); p <- g / sum(g)
    A1[, l] <- sample.int(K, ntot, replace = TRUE, prob = p)
    A2[, l] <- sample.int(K, ntot, replace = TRUE, prob = p)
  }
  cumM <- t(apply(M, 1, cumsum))
  pick_parent <- function(dm) {
    u <- stats::runif(length(dm))
    k <- rowSums(u > cumM[dm, , drop = FALSE]) + 1L
    k <- pmin(k, d)
    start[k] + ceiling(stats::runif(length(dm)) * N[k])
  }
  for (t in seq_len(T_gen)) {
    p1 <- pick_parent(deme)
    p2 <- pick_parent(deme)
    self <- stats::runif(ntot) < selfing
    p2[self] <- p1[self]
    s1 <- matrix(stats::runif(ntot * n_loci) < 0.5, ntot, n_loci)
    s2 <- matrix(stats::runif(ntot * n_loci) < 0.5, ntot, n_loci)
    g1 <- ifelse(s1, A1[p1, , drop = FALSE], A2[p1, , drop = FALSE])
    g2 <- ifelse(s2, A1[p2, , drop = FALSE], A2[p2, , drop = FALSE])
    if (mu > 0) {
      m1 <- stats::runif(ntot * n_loci) < mu
      if (any(m1)) g1[m1] <- sample.int(K, sum(m1), replace = TRUE)
      m2 <- stats::runif(ntot * n_loci) < mu
      if (any(m2)) g2[m2] <- sample.int(K, sum(m2), replace = TRUE)
    }
    A1 <- matrix(as.integer(g1), ntot, n_loci); A2 <- matrix(as.integer(g2), ntot, n_loci)
  }
  list(A1 = A1, A2 = A2, deme = deme)
}

ps_simulate <- function(params, design, seed) {
  set.seed(seed)
  model <- design$model
  d <- as.integer(design$n_pops)
  marker <- if (is.null(design$marker)) "microsat" else design$marker
  K <- if (marker == "snp") 2L else as.integer(params$n_alleles)
  mu <- if (marker == "snp") 0 else params$mu
  alpha <- if (marker == "snp") params$snp_alpha else params$alpha
  spacing <- design$spacing_km
  xy <- ps_lattice(model, d, spacing)
  M <- ps_migration_matrix(model, d, params$m, xy, spacing)
  N <- rep_len(params$Ne, d)
  pops <- paste0("P", seq_len(d))
  n_shadow <- if (is.null(design$n_shadow_loci)) 300L else as.integer(design$n_shadow_loci)
  truth_type <- if (is.null(design$truth)) "realised" else design$truth
  if (truth_type != "realised") n_shadow <- 0L
  L <- as.integer(design$n_loci)
  sim <- ps_wf_forward(M, N, params$selfing, mu, K, alpha, params$T_gen, L + n_shadow)
  # sample individuals per deme
  n_ind <- rep_len(design$n_ind, d)
  idx <- unlist(lapply(seq_len(d), function(i) {
    w <- which(sim$deme == i); w[sample.int(length(w), min(n_ind[i], length(w)))]
  }))
  a1 <- sim$A1[idx, seq_len(L), drop = FALSE]; a2 <- sim$A2[idx, seq_len(L), drop = FALSE]
  if (marker == "snp") { a1 <- a1 - 1L; a2 <- a2 - 1L }   # 0/1 allele codes
  else { a1 <- 100L + 2L * a1; a2 <- 100L + 2L * a2 }     # allele sizes in bp (dinucleotide)
  miss <- if (is.null(design$missing_rate)) 0 else design$missing_rate
  if (miss > 0) {
    z <- matrix(stats::runif(length(a1)) < miss, nrow(a1)); a1[z] <- NA; a2[z] <- NA
  }
  loci <- if (marker == "snp") sprintf("snp%04d", seq_len(L)) else sprintf("Loc%02d", seq_len(L))
  colnames(a1) <- colnames(a2) <- loci
  pop <- factor(pops[sim$deme[idx]], levels = pops)
  ind <- sprintf("%s_%03d", as.character(pop), stats::ave(seq_along(idx), pop, FUN = seq_along))
  coords <- data.frame(pop = pops, x = xy$x, y = xy$y)
  ibd_x <- if (is.null(design$ibd_distance)) "log" else design$ibd_distance
  data <- ps_make_data(a1, a2, ind, pop, marker = marker, coords = coords, coord_type = "km",
                       ibd_distance = ibd_x, source = sprintf("simulated %s, %d demes", model, d))
  keep_pairs <- c(paste0("Fst_P1_vs_P2"), paste0("Fst_P1_vs_P", d))
  if (truth_type == "realised") {
    # Realised (genome-wide) values for THESE populations: census of every individual alive in the
    # final generation at the sampled + shadow loci (same pedigree), computed with the same estimators.
    cen <- ps_make_data(sim$A1, sim$A2, seq_len(nrow(sim$A1)), factor(pops[sim$deme], levels = pops),
                        coords = coords, ibd_distance = ibd_x)
    colnames(cen$a1) <- colnames(cen$a2) <- cen$loci <- paste0("L", seq_len(ncol(sim$A1)))
    cc <- ps_components(cen)
    st <- ps_stats_W(cc, matrix(1, ncol(sim$A1), 1), ps_pair_distances(cen, cc), ibd_x)[, 1]
    truth <- st[c("Fst_global", "Fis_global", "He_mean", "Ho_mean", "ibd_slope", keep_pairs)]
  } else {
    # Expected values over replicate histories (exact identity recursion)
    id <- ps_identity_recursion(M, N, params$selfing, mu, K, alpha, params$T_gen)
    tr <- ps_truth_from_identity(id, pops, xy, ibd_x)
    truth <- c(tr$global, ibd_slope = tr$ibd_slope)
    truth <- c(truth, tr$pairs[intersect(keep_pairs, names(tr$pairs))])
  }
  list(data = data, truth = truth)
}

# Equilibrium approximation for the finite island model (Wright 1931; with
# the finite-island correction, e.g. Rousset 1996 / Whitlock & McCauley 1999).
ps_wright_island_fst <- function(Ne, m, n_pops) 1 / (1 + 4 * Ne * m * (n_pops / (n_pops - 1))^2)
