# Estimation engine for popgen-structure.
#
# Internal data structure (returned by both simulate() and read_data()):
#   list(a1, a2      integer matrices individuals x loci: the two allele codes
#                    (microsat: allele size; SNP: 0/1), NA = missing
#        ind         individual ids;  pop  factor of populations
#        loci        locus names;     marker "microsat" | "snp"
#        coords      data.frame(pop, x, y) or NULL (x/y in km, or lon/lat)
#        coord_type  "km" | "lonlat";  ibd_distance "log" | "linear"
#        warnings    character: non-fatal issues found while reading
#        problems    character: fatal issues found while reading)
#
# Estimators: Weir & Cockerham (1984) theta (FST) and f (FIS) as ratios of
# sums of per-allele variance components over loci; Nei (1987) unbiased gene
# diversity Hs (He) and observed heterozygosity Ho; rarefied allelic
# richness (El Mousadik & Petit 1996, Hurlbert 1971); private alleles.
# Uncertainty: delete-one-locus jackknife over loci (loci are the replicates of
# the evolutionary process; Weir 1996) plus the cross-locus covariance caused by
# loci sharing the same individuals (delete-one-individual jackknife, see ps_uncertainty);
# intervals estimate +/- t(L-1) SE. Validation showed the plain percentile
# bootstrap over loci undercovers (VALIDATION.md).

ps_make_data <- function(a1, a2, ind, pop, marker = "microsat", coords = NULL, coord_type = "km",
                         ibd_distance = "log", source = "", warnings = character(),
                         problems = character()) {
  storage.mode(a1) <- "integer"; storage.mode(a2) <- "integer"
  bad <- is.na(a1) | is.na(a2); a1[bad] <- NA; a2[bad] <- NA
  pop <- droplevels(as.factor(pop))
  list(a1 = a1, a2 = a2, ind = as.character(ind), pop = pop, loci = colnames(a1),
       marker = marker, coords = coords, coord_type = coord_type, ibd_distance = ibd_distance,
       source = source, warnings = warnings, problems = problems)
}

# ---- per-locus tables ------------------------------------------------------

ps_tables <- function(data) {
  pops <- levels(data$pop); r <- length(pops); pi <- as.integer(data$pop)
  L <- ncol(data$a1)
  lapply(seq_len(L), function(l) {
    g1 <- data$a1[, l]; g2 <- data$a2[, l]
    ok <- !is.na(g1)
    al <- sort(unique(c(g1[ok], g2[ok])))
    k <- length(al)
    n <- tabulate(pi[ok], r)
    if (k == 0) return(list(alleles = al, n = n, x = matrix(0, r, 0), h = matrix(0, r, 0), nhet = rep(0, r)))
    i1 <- match(g1[ok], al); i2 <- match(g2[ok], al); pp <- pi[ok]
    x <- matrix(tabulate((i1 - 1) * r + pp, r * k) + tabulate((i2 - 1) * r + pp, r * k), r, k)
    het <- i1 != i2
    hc <- matrix(tabulate((i1[het] - 1) * r + pp[het], r * k) + tabulate((i2[het] - 1) * r + pp[het], r * k), r, k)
    list(alleles = al, n = n, x = x, h = hc, nhet = tabulate(pp[het], r))
  })
}

# Weir & Cockerham (1984) variance components summed over alleles, for the
# populations in `rows` (populations with no typed individuals are skipped).
ps_wc_abc <- function(tb, rows) {
  n <- tb$n[rows]; keep <- n > 0
  if (sum(keep) < 2 || ncol(tb$x) < 2) return(c(0, 0, 0))
  n <- n[keep]; x <- tb$x[rows, , drop = FALSE][keep, , drop = FALSE]
  hc <- tb$h[rows, , drop = FALSE][keep, , drop = FALSE]
  r <- length(n); nbar <- sum(n) / r
  if (nbar <= 1) return(c(0, 0, 0))
  nc <- (r * nbar - sum(n^2) / (r * nbar)) / (r - 1)
  p <- x / (2 * n); h <- hc / n
  pbar <- colSums(n * p) / (r * nbar)
  s2 <- colSums(n * sweep(p, 2, pbar)^2) / ((r - 1) * nbar)
  hbar <- colSums(n * h) / (r * nbar)
  pq <- pbar * (1 - pbar)
  a <- nbar / nc * (s2 - 1 / (nbar - 1) * (pq - (r - 1) / r * s2 - hbar / 4))
  b <- nbar / (nbar - 1) * (pq - (r - 1) / r * s2 - (2 * nbar - 1) / (4 * nbar) * hbar)
  c <- hbar / 2
  c(sum(a), sum(b), sum(c))
}

# Same components for every pair of populations at once (vectorised).
# pr: 2 x P matrix of population indices. Returns P x 3 (a, b, c).
ps_wc_pairs <- function(tb, pr) {
  P <- ncol(pr); out <- matrix(0, P, 3)
  if (ncol(tb$x) < 2) return(out)
  n1 <- tb$n[pr[1, ]]; n2 <- tb$n[pr[2, ]]
  ok <- n1 > 0 & n2 > 0 & (n1 + n2) / 2 > 1
  if (!any(ok)) return(out)
  n1 <- n1[ok]; n2 <- n2[ok]
  x1 <- tb$x[pr[1, ok], , drop = FALSE]; x2 <- tb$x[pr[2, ok], , drop = FALSE]
  h1 <- tb$h[pr[1, ok], , drop = FALSE] / n1; h2 <- tb$h[pr[2, ok], , drop = FALSE] / n2
  p1 <- x1 / (2 * n1); p2 <- x2 / (2 * n2)
  nbar <- (n1 + n2) / 2
  nc <- 2 * nbar - (n1^2 + n2^2) / (2 * nbar)
  pbar <- (n1 * p1 + n2 * p2) / (2 * nbar)
  s2 <- (n1 * (p1 - pbar)^2 + n2 * (p2 - pbar)^2) / nbar
  hbar <- (n1 * h1 + n2 * h2) / (2 * nbar)
  pq <- pbar * (1 - pbar)
  a <- nbar / nc * (s2 - 1 / (nbar - 1) * (pq - s2 / 2 - hbar / 4))
  b <- nbar / (nbar - 1) * (pq - s2 / 2 - (2 * nbar - 1) / (4 * nbar) * hbar)
  out[ok, ] <- cbind(rowSums(a), rowSums(b), rowSums(hbar / 2))
  out
}

# Naive (deliberately WRONG, for the negative control): Nei's GST from
# sample frequencies with no sample-size corrections; Hs = 1 - sum p^2.
ps_naive_parts <- function(tb, rows) {
  n <- tb$n[rows]; keep <- n > 0
  if (sum(keep) < 2 || ncol(tb$x) < 1) return(c(0, 0))
  x <- tb$x[rows, , drop = FALSE][keep, , drop = FALSE]; n <- n[keep]
  p <- x / (2 * n)
  hs <- mean(1 - rowSums(p^2)); ht <- 1 - sum(colMeans(p)^2)
  c(ht - hs, ht)
}

# All per-locus components the bootstrap needs, as matrices with loci in
# columns, so any statistic is a ratio of weighted sums over loci.
ps_components <- function(data, tabs = ps_tables(data), estimator = "WC", g = NULL) {
  pops <- levels(data$pop); r <- length(pops); L <- length(tabs)
  pr <- utils::combn(r, 2)
  pair_names <- paste0(pops[pr[1, ]], "_vs_", pops[pr[2, ]])
  P <- ncol(pr)
  glob <- matrix(0, 3, L); pa <- pb <- pc <- matrix(0, P, L)
  Hs <- Ho <- HsN <- AR <- mask <- matrix(0, r, L)
  nmin <- min(vapply(tabs, function(t) { n <- t$n[t$n > 0]; if (length(n)) min(n) else Inf }, 1))
  if (is.null(g)) g <- 2 * nmin  # rarefaction size (genes)
  for (l in seq_len(L)) {
    tb <- tabs[[l]]
    n <- tb$n; ok <- n >= 2
    mask[, l] <- ok
    if (ncol(tb$x)) {
      p <- tb$x / pmax(2 * n, 1)
      ho <- tb$nhet / pmax(n, 1)
      sp2 <- rowSums(p^2)
      Ho[, l] <- ifelse(ok, ho, 0)
      Hs[, l] <- ifelse(ok, n / (n - 1) * (1 - sp2 - ho / (2 * n)), 0)
      HsN[, l] <- ifelse(ok, 1 - sp2, 0)
      if (is.finite(g)) {
        N2 <- 2 * n
        AR[, l] <- ifelse(ok, rowSums(1 - exp(lchoose(N2 - tb$x, g) - lchoose(N2, g))), 0)
      }
    }
    if (estimator == "WC") {
      glob[, l] <- ps_wc_abc(tb, seq_len(r))
      v <- ps_wc_pairs(tb, pr); pa[, l] <- v[, 1]; pb[, l] <- v[, 2]; pc[, l] <- v[, 3]
    } else {
      v <- ps_naive_parts(tb, seq_len(r)); glob[, l] <- c(v[1], v[2] - v[1], 0)
      for (q in seq_len(P)) { v <- ps_naive_parts(tb, pr[, q]); pa[q, l] <- v[1]; pb[q, l] <- v[2] - v[1]; pc[q, l] <- 0 }
    }
  }
  # private alleles: alleles seen in exactly one population
  pa_count <- rowSums(vapply(tabs, function(tb) {
    if (!ncol(tb$x)) return(rep(0, r))
    pres <- tb$x > 0; one <- colSums(pres) == 1
    rowSums(pres[, one, drop = FALSE])
  }, numeric(r)))
  list(pops = pops, pairs = pr, pair_names = pair_names, glob = glob, pa = pa, pb = pb, pc = pc,
       Hs = Hs, HsN = HsN, Ho = Ho, AR = AR, mask = mask, rarefaction_genes = g,
       private = stats::setNames(pa_count, pops), n_alleles = vapply(tabs, function(t) length(t$alleles), 1),
       n_typed = vapply(tabs, function(t) t$n, numeric(r)), estimator = estimator)
}

ps_pair_distances <- function(data, comps) {
  if (is.null(data$coords)) return(NULL)
  cc <- data$coords[match(comps$pops, data$coords$pop), ]
  if (anyNA(cc$x)) return(NULL)
  if (identical(data$coord_type, "lonlat")) {
    eco_require("geosphere")
    D <- geosphere::distm(cbind(cc$x, cc$y), fun = geosphere::distHaversine) / 1000
  } else D <- as.matrix(stats::dist(cbind(cc$x, cc$y)))
  D[t(comps$pairs)]
}

# Statistics for a matrix of locus weights W (L x B): column b gives the
# multiplicity of each locus in bootstrap sample b (all ones = point estimate).
ps_stats_W <- function(comps, W, dist = NULL, ibd_distance = "log") {
  W <- as.matrix(W)
  sg <- comps$glob %*% W                     # 3 x B
  if (comps$estimator == "WC") {
    fst <- sg[1, ] / colSums(sg)
    fis <- 1 - sg[3, ] / (sg[2, ] + sg[3, ])
  } else {
    fst <- sg[1, ] / (sg[1, ] + sg[2, ])
    fis <- 1 - colSums((comps$Ho %*% W)) / colSums((comps$HsN %*% W))
  }
  mw <- comps$mask %*% W                     # r x B: loci used per pop
  Hs_src <- if (comps$estimator == "WC") comps$Hs else comps$HsN
  he_pop <- (Hs_src %*% W) / mw
  ho_pop <- (comps$Ho %*% W) / mw
  ar_pop <- (comps$AR %*% W) / mw
  fis_pop <- 1 - (comps$Ho %*% W) / (Hs_src %*% W)
  A <- comps$pa %*% W; Bm <- comps$pb %*% W; C <- comps$pc %*% W
  fst_pair <- A / (A + Bm + C)
  out <- rbind(Fst_global = fst, Fis_global = fis,
               He_mean = colSums(Hs_src %*% W) / colSums(mw),
               Ho_mean = colSums(comps$Ho %*% W) / colSums(mw))
  rownames(he_pop) <- paste0("He_", comps$pops); rownames(ho_pop) <- paste0("Ho_", comps$pops)
  rownames(ar_pop) <- paste0("AR_", comps$pops); rownames(fis_pop) <- paste0("Fis_", comps$pops)
  rownames(fst_pair) <- paste0("Fst_", comps$pair_names)
  out <- rbind(out, he_pop, ho_pop, ar_pop, fis_pop, fst_pair)
  if (!is.null(dist)) {
    ok <- dist > 0
    if (sum(ok) >= 3) {
      x <- if (ibd_distance == "log") log(dist[ok]) else dist[ok]
      y <- fst_pair[ok, , drop = FALSE] / (1 - fst_pair[ok, , drop = FALSE])
      xc <- x - mean(x)
      slope <- colSums(xc * y) / sum(xc^2)
      out <- rbind(out, ibd_slope = slope)
    }
  }
  out
}

# Multiallelic Hardy-Weinberg test: Monte Carlo version of the conditional
# exact test (Guo & Thompson 1992), permuting alleles among individuals;
# exact enumeration for biallelic loci (Wigginton et al. 2005).
ps_hwe_one <- function(g1, g2, B = 1000) {
  ok <- !is.na(g1); g1 <- g1[ok]; g2 <- g2[ok]; n <- length(g1)
  if (n < 5) return(c(p = NA_real_, Fhat = NA_real_))
  al <- unique(c(g1, g2)); k <- length(al)
  if (k < 2) return(c(p = NA_real_, Fhat = NA_real_))
  ho <- mean(g1 != g2); pf <- table(c(g1, g2)) / (2 * n)
  he <- (1 - sum(pf^2)) * 2 * n / (2 * n - 1)
  Fhat <- 1 - ho / he
  if (k == 2) {
    na <- sum(c(g1, g2) == al[1]); nb <- 2 * n - na; nmin <- min(na, nb)
    hs <- seq(nmin %% 2, nmin, by = 2)
    lp <- hs * log(2) + lfactorial(n) - lfactorial((nmin - hs) / 2) - lfactorial(hs) -
      lfactorial(n - hs - (nmin - hs) / 2) + lfactorial(na) + lfactorial(nb) - lfactorial(2 * n)
    pr <- exp(lp); hobs <- sum(g1 != g2)
    pobs <- pr[hs == hobs]
    return(c(p = min(1, sum(pr[pr <= pobs * (1 + 1e-7)])), Fhat = Fhat))
  }
  stat <- function(a, b) {
    lo <- pmin(a, b); hi <- pmax(a, b)
    sum(a != b) * log(2) - sum(lfactorial(table(paste(lo, hi))))
  }
  s0 <- stat(g1, g2)
  pool <- c(g1, g2)
  sb <- vapply(seq_len(B), function(i) { z <- sample(pool); stat(z[1:n], z[(n + 1):(2 * n)]) }, 1)
  c(p = (sum(sb <= s0 + 1e-9) + 1) / (B + 1), Fhat = Fhat)
}

ps_hwe <- function(data, B = 1000) {
  pops <- levels(data$pop)
  res <- do.call(rbind, lapply(pops, function(pp) {
    w <- data$pop == pp
    do.call(rbind, lapply(seq_along(data$loci), function(l) {
      v <- ps_hwe_one(data$a1[w, l], data$a2[w, l], B)
      data.frame(pop = pp, locus = data$loci[l], p = v[["p"]], Fis = v[["Fhat"]])
    }))
  }))
  res$q <- stats::ave(res$p, res$pop, FUN = function(p) { q <- p; ok <- !is.na(p); q[ok] <- stats::p.adjust(p[ok], "BH"); q })
  res
}

# ---- adegenet objects, PCA, DAPC, AMOVA -----------------------------------

ps_genind <- function(data) {
  eco_require("adegenet")
  g <- matrix(paste(data$a1, data$a2, sep = "/"), nrow(data$a1))
  g[is.na(data$a1)] <- NA
  colnames(g) <- make.names(data$loci, unique = TRUE)
  suppressWarnings(adegenet::df2genind(as.data.frame(g, stringsAsFactors = FALSE), sep = "/",
                                        ind.names = data$ind, pop = data$pop, ploidy = 2, NA.char = NA))
}

ps_pca <- function(gi, nf = 4) {
  X <- adegenet::tab(gi, freq = TRUE, NA.method = "mean")
  X <- X[, apply(X, 2, stats::var) > 0, drop = FALSE]
  pc <- ade4::dudi.pca(X, center = TRUE, scale = FALSE, scannf = FALSE, nf = nf)
  list(scores = pc$li, eig = pc$eig, pct = 100 * pc$eig / sum(pc$eig))
}

ps_dapc <- function(gi, settings, n_pops) {
  set.seed(if (is.null(settings$seed)) 1 else settings$seed)
  ni <- adegenet::nInd(gi)
  maxk <- min(if (is.null(settings$max_k)) max(2 * n_pops, 10) else settings$max_k, ni - 1, 40)
  X <- adegenet::tab(gi, freq = TRUE, NA.method = "mean")
  npca <- min(ni - 1, ncol(X))
  fc <- adegenet::find.clusters(gi, n.pca = npca, max.n.clust = maxk, choose.n.clust = FALSE,
                                criterion = if (is.null(settings$k_criterion)) "diffNgroup" else settings$k_criterion,
                                n.iter = 1e5, n.start = 20)
  K <- length(levels(fc$grp))
  bic <- data.frame(K = as.integer(sub("K=", "", names(fc$Kstat))), BIC = as.numeric(fc$Kstat))
  out <- list(K = K, bic = bic, grp = fc$grp, posterior = NULL, n_pca = NA)
  if (K >= 2) {
    # Thia (2023): retain only K - 1 PCs (the axes carrying among-group structure)
    np <- max(1, K - 1)
    dp <- adegenet::dapc(gi, fc$grp, n.pca = np, n.da = K - 1)
    out$posterior <- dp$posterior; out$n_pca <- np
  }
  out
}

ps_amova <- function(gi, nperm = 999) {
  eco_require(c("poppr", "ade4"))
  suppressPackageStartupMessages(require("ade4", character.only = TRUE))  # poppr.amova calls quasieuclid() unqualified
  adegenet::strata(gi) <- data.frame(pop = adegenet::pop(gi))
  am <- poppr::poppr.amova(gi, ~pop, quiet = TRUE, within = FALSE)
  tst <- ade4::randtest(am, nrepet = nperm)
  list(table = am$componentsofcovariance, phi = am$statphi, p = tst$pvalue)
}

# ---- main fit ---------------------------------------------------------------

ps_fit <- function(data, settings) {
  st <- settings
  est_name <- if (is.null(st$estimator)) "WC" else st$estimator
  set.seed(if (is.null(st$seed)) 1 else st$seed)
  tabs <- ps_tables(data)
  comps <- ps_components(data, tabs, est_name)
  L <- length(tabs)
  dist <- ps_pair_distances(data, comps)
  ibd_x <- if (is.null(data$ibd_distance)) "log" else data$ibd_distance
  point <- ps_stats_W(comps, matrix(1, L, 1), dist, ibd_x)[, 1]
  unc <- ps_uncertainty(data, comps, dist, ibd_x, st)
  prob <- if (is.null(st$ci)) 0.95 else st$ci
  a <- (1 - prob) / 2
  tq <- stats::qt(1 - a, max(L - 1, 1))
  se <- unc$se[names(point)]
  est <- data.frame(param = names(point), estimate = unname(point), sd = unname(se),
                    lower = unname(point - tq * se), upper = unname(point + tq * se), row.names = NULL)
  # Private alleles (counts; no interval) and island-model Nm (derived)
  est <- rbind(est, data.frame(param = paste0("PrivateAlleles_", comps$pops), estimate = unname(comps$private),
                               sd = NA_real_, lower = NA_real_, upper = NA_real_))
  nm <- function(f) ifelse(f > 0, (1 / f - 1) / 4, Inf)
  fr <- est[est$param == "Fst_global", ]
  est <- rbind(est, data.frame(param = "Nm_island_model", estimate = nm(point[["Fst_global"]]), sd = NA_real_,
                               lower = nm(fr$upper), upper = nm(fr$lower)))
  notes <- character()
  mono <- sum(comps$n_alleles < 2)
  notes <- c(notes, sprintf("%d loci, %d populations, %d individuals; %d monomorphic loci (they add nothing to FST but count as He = 0 in diversity means).",
                            L, length(comps$pops), length(data$ind), mono))
  notes <- c(notes, sprintf("Intervals: estimate +/- t(%d df) x SE, SE from a delete-one-locus jackknife plus the cross-locus covariance from a delete-one-individual jackknife over %d individuals (median share of variance from shared individuals: %.0f%%). Allelic richness rarefied to %d genes (the smallest population sample).",
                            L - 1, unc$B, 100 * stats::median(unc$cov_share[c("Fst_global", "Fis_global", "He_mean")], na.rm = TRUE), comps$rarefaction_genes))
  if (est_name != "WC") notes <- c(notes, "WARNING: naive (uncorrected) estimator requested; for validation only.")
  if (length(data$warnings)) notes <- c(notes, data$warnings)
  obj <- list(comps = comps, jack = unc$jack, se = unc$se, se_jack = unc$se_jack, tq = tq,
              dist = dist, ibd_distance = ibd_x, settings = st)
  # optional heavier analyses (real-data fits)
  if (!is.null(dist) && isTRUE((st$mantel_perm %||% 999) > 0)) {
    eco_require("vegan")
    pr <- comps$pairs; r <- length(comps$pops)
    Fm <- matrix(0, r, r); f <- point[paste0("Fst_", comps$pair_names)]
    Fm[t(pr)] <- f / (1 - f); Fm <- Fm + t(Fm)
    Dm <- matrix(0, r, r); Dm[t(pr)] <- if (ibd_x == "log") log(pmax(dist, 1e-9)) else dist; Dm <- Dm + t(Dm)
    mt <- vegan::mantel(stats::as.dist(Fm), stats::as.dist(Dm), permutations = st$mantel_perm %||% 999)
    obj$mantel <- list(r = unname(mt$statistic), p = mt$signif, perm = mt$permutations)
  }
  if (isTRUE(st$hwe %||% TRUE)) obj$hwe <- ps_hwe(data, st$hwe_B %||% 1000)
  if (isTRUE(st$pca %||% TRUE) || isTRUE(st$dapc %||% TRUE) || isTRUE(st$amova)) gi <- ps_genind(data)
  if (isTRUE(st$pca %||% TRUE)) obj$pca <- tryCatch(ps_pca(gi), error = function(e) { notes <<- c(notes, paste("PCA failed:", conditionMessage(e))); NULL })
  if (isTRUE(st$dapc %||% TRUE)) obj$dapc <- tryCatch(ps_dapc(gi, st, length(comps$pops)),
                                                     error = function(e) { notes <<- c(notes, paste("DAPC failed:", conditionMessage(e))); NULL })
  if (isTRUE(st$amova)) obj$amova <- tryCatch(ps_amova(gi, st$amova_perm %||% 999),
                                              error = function(e) { notes <<- c(notes, paste("AMOVA failed:", conditionMessage(e))); NULL })
  if (!is.null(obj$hwe)) {
    h <- obj$hwe; ns <- sum(h$q < 0.05, na.rm = TRUE)
    byloc <- tapply(h$q < 0.05, h$locus, function(z) mean(z, na.rm = TRUE))
    badloc <- names(byloc)[byloc > 0.5 & !is.na(byloc)]
    notes <- c(notes, sprintf("Hardy-Weinberg: %d of %d population x locus tests deviate after Benjamini-Hochberg correction within population (q < 0.05).%s",
                              ns, sum(!is.na(h$q)),
                              if (length(badloc)) paste0(" Loci deviating in more than half the populations (check for null alleles or scoring errors): ", paste(badloc, collapse = ", "), ".") else ""))
  }
  ok <- all(is.finite(point[c("Fst_global", "Fis_global", "He_mean")])) && sum(comps$n_alleles >= 2) >= 2
  if (!ok) notes <- c(notes, "Fewer than two polymorphic loci or non-finite estimates: results are not usable.")
  list(estimates = est,
       diagnostics = list(converged = ok, notes = notes),   # likelihood-free estimators: no R-hat/ESS
       object = obj)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Standard errors.
#  (1) delete-one-locus jackknife: genetic (evolutionary) + per-locus sampling variance,
#      treating loci as independent;
#  (2) plus the covariance BETWEEN loci created because all loci are measured on the same
#      individuals (inbred individuals are homozygous at many loci; relatives share alleles at
#      many loci). Estimated with a delete-one-individual jackknife stratified by population:
#      V_joint (individual deleted at all loci) - V_sep (deleted at one locus at a time, summed
#      over loci). Pseudo-values are centred within populations so the common effect of n_i - 1
#      cancels. Bootstrapping individuals with replacement was rejected: the number of duplicated
#      individuals is shared by all loci and inflated FST intervals near 0 (VALIDATION.md).
ps_uncertainty <- function(data, comps, dist, ibd_x, st) {
  L <- ncol(data$a1)
  J <- matrix(1, L, L); diag(J) <- 0
  jk <- ps_stats_W(comps, J, dist, ibd_x)
  se_jack <- sqrt((L - 1) / L * rowSums((jk - rowMeans(jk, na.rm = TRUE))^2, na.rm = TRUE))
  extra <- rep(0, nrow(jk))
  if (!isFALSE(st$individual_covariance)) {
    mats <- c("glob", "pa", "pb", "pc", "Hs", "HsN", "Ho", "AR", "mask")
    # mixing weights: column l uses the reduced data at locus l only
    Wmix <- rbind(1 - diag(L), diag(L))
    pops <- as.integer(data$pop); n_i <- tabulate(pops, nlevels(data$pop))
    S <- nrow(jk)
    th_joint <- matrix(NA_real_, S, length(pops)); th_sep <- array(NA_real_, c(S, L, length(pops)))
    for (j in seq_along(pops)) {
      dd <- data; dd$a1 <- data$a1[-j, , drop = FALSE]; dd$a2 <- data$a2[-j, , drop = FALSE]; dd$pop <- data$pop[-j]
      cj <- ps_components(dd, estimator = comps$estimator, g = comps$rarefaction_genes)
      th_joint[, j] <- ps_stats_W(cj, matrix(1, L, 1), dist, ibd_x)[, 1]
      aug <- comps
      for (k in mats) aug[[k]] <- cbind(comps[[k]], cj[[k]])
      th_sep[, , j] <- ps_stats_W(aug, Wmix, dist, ibd_x)
    }
    vj <- vs <- rep(0, S)
    for (i in seq_along(n_i)) {
      w <- which(pops == i); ni <- length(w); if (ni < 2) next
      dj <- th_joint[, w, drop = FALSE]; dj <- dj - rowMeans(dj)
      vj <- vj + (ni - 1) / ni * rowSums(dj^2)
      ds <- th_sep[, , w, drop = FALSE]
      ds <- sweep(ds, c(1, 2), apply(ds, c(1, 2), mean))
      vs <- vs + (ni - 1) / ni * apply(ds^2, 1, sum)
    }
    extra <- pmax(vj - vs, 0)
    extra[!is.finite(extra) | grepl("^AR_", rownames(jk))] <- 0   # AR: rarefaction size fixed by the smallest sample; loci jackknife only
  }
  se <- sqrt(se_jack^2 + extra)
  names(se) <- names(se_jack) <- rownames(jk)
  list(se = se, se_jack = se_jack, jack = jk, B = length(data$ind),
       cov_share = stats::setNames(extra / se^2, rownames(jk)))
}
