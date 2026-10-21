# Engine for the community module.
#
# Estimands with intervals are computed by hand (closed-form least squares on
# the Hellinger-transformed matrix; sums of squared distances for PERMANOVA;
# pairwise Baselga components) so that a delete-one-site jackknife is cheap and
# exact. Point estimates agree with vegan::varpart / RsquareAdj, adonis2 and
# betapart::beta.pair to numerical precision (VALIDATION.md, cross-check).
# Permutation tests, ordinations, betadisper and indicator species use vegan
# and indicspecies.

com_transform <- function(Y, method = "hellinger") {
  switch(method,
    hellinger = sqrt(Y / rowSums(Y)),
    chord = Y / sqrt(rowSums(Y^2)),
    log = log1p(Y),
    none = Y,
    stop("transform must be hellinger, chord, log or none"))
}

adj_r2 <- function(r2, n, p) 1 - (1 - r2) * (n - 1) / (n - p - 1)

# Model matrix (no intercept column) for a set of env variables.
com_design <- function(env, vars, poly = 1) {
  rhs <- vapply(vars, function(v) {
    if (is.numeric(env[[v]]) && poly > 1) sprintf("poly(%s, %d, raw = TRUE)", v, poly) else v
  }, "")
  X <- stats::model.matrix(stats::as.formula(paste("~", paste(rhs, collapse = " + "))), env)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

# R2 (unadjusted) of the multivariate regression of Y on X, and its rank.
r2_rows <- function(Y, X) {
  Yc <- scale(Y, scale = FALSE)
  if (!ncol(X)) return(c(r2 = 0, p = 0))
  Xc <- scale(X, scale = FALSE)
  q <- qr(Xc)
  c(r2 = sum(qr.fitted(q, Yc)^2) / sum(Yc^2), p = q$rank)
}

# Variation partitioning (Borcard et al. 1992; Peres-Neto et al. 2006) for a
# named list of model matrices. Returns a named vector of fractions:
#   R2adj_<set> (marginal), R2adj_unique_<set> (conditional on all others),
#   R2adj_shared (total - sum of unique; for 2 sets this is fraction [b]),
#   R2adj_env_total.
com_varpart_vec <- function(Y, Xs, adjusted = TRUE) {
  n <- nrow(Y)
  fr <- function(X) { r <- r2_rows(Y, X); if (adjusted) adj_r2(r[["r2"]], n, r[["p"]]) else r[["r2"]] }
  k <- length(Xs); nm <- names(Xs)
  out <- c()
  marg <- vapply(Xs, fr, 0)
  tot <- fr(do.call(cbind, Xs))
  out[paste0("R2adj_", nm)] <- marg
  if (k > 1) {
    uni <- vapply(seq_len(k), function(i) tot - fr(do.call(cbind, Xs[-i])), 0)
    out[paste0("R2adj_unique_", nm)] <- uni
    out["R2adj_shared"] <- tot - sum(uni)
  }
  out["R2adj_env_total"] <- tot
  out
}

# Sum of squares decomposition for one-way PERMANOVA from squared distances.
permanova_r2 <- function(D2, g) {
  n <- length(g)
  sst <- sum(D2[upper.tri(D2)]) / n
  ssw <- 0
  for (lv in unique(g)) { i <- which(g == lv); if (length(i) > 1) ssw <- ssw + sum(D2[i, i][upper.tri(D2[i, i])]) / length(i) }
  1 - ssw / sst
}

# Jackknife of permanova adjusted R2 (vectorised delete-one-site).
permanova_jack <- function(D2, g, adjusted = TRUE) {
  g <- as.integer(factor(g)); n <- length(g); G <- max(g)
  U <- D2; U[lower.tri(U, diag = TRUE)] <- 0
  tot <- sum(U); rs <- rowSums(D2)            # D2 symmetric, zero diagonal
  ng <- tabulate(g, G)
  wsum <- vapply(seq_len(G), function(a) { i <- which(g == a); sum(D2[i, i]) / 2 }, 0)
  rsw <- vapply(seq_len(n), function(i) sum(D2[i, g == g[i]]), 0)
  f <- function(r2, nn) if (adjusted) adj_r2(r2, nn, G - 1) else r2
  full <- f(1 - sum(wsum / ng) / (tot / n), n)
  jk <- vapply(seq_len(n), function(i) {
    ws <- wsum; nn <- ng; ws[g[i]] <- ws[g[i]] - rsw[i]; nn[g[i]] <- nn[g[i]] - 1
    ssw <- sum(ifelse(nn > 0, ws / pmax(nn, 1), 0))
    f(1 - ssw / ((tot - rs[i]) / (n - 1)), n - 1)
  }, 0)
  list(estimate = full, jack = jk)
}

# Pairwise Baselga (2010) Sorensen components for all site pairs (matrices).
beta_pair_mats <- function(P) {
  P <- (P > 0) * 1
  a <- tcrossprod(P); sr <- rowSums(P)
  b <- sr - a; c <- t(b)           # b[i,j] = species only in i; c = only in j
  mn <- pmin(b, c)
  sor <- (b + c) / (2 * a + b + c)
  sim <- ifelse(a + mn > 0, mn / (a + mn), 0)
  diag(sor) <- 0; diag(sim) <- 0
  list(sor = sor, sim = sim, sne = sor - sim)
}

# Mean over pairs and delete-one-site jackknife values of that mean.
pair_mean_jack <- function(M) {
  n <- nrow(M); tot <- sum(M[upper.tri(M)]); rs <- rowSums(M) - diag(M)
  list(estimate = tot / (n * (n - 1) / 2),
       jack = (tot - rs) / ((n - 1) * (n - 2) / 2))
}

jack_rows <- function(name, est, jk, level = 0.95) {
  n <- length(jk)
  se <- sqrt((n - 1) / n * sum((jk - mean(jk))^2))
  tq <- stats::qt(1 - (1 - level) / 2, n - 1)
  data.frame(param = name, estimate = est, sd = se, lower = est - tq * se, upper = est + tq * se,
             row.names = NULL)
}

# ---- main fit ---------------------------------------------------------------

com_fit <- function(data, settings) {
  st <- merge_list(com_default_settings(), settings)
  # reproducible permutations without resetting the caller's random stream
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) rm(list = ".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(st$seed %||% 1)
  Y <- as.matrix(data$comm); env <- data$env
  notes <- as.character(data$notes %||% character())   # read_data notes (ignored sites, dropped species...) reach the report
  zero <- colSums(Y > 0, na.rm = TRUE) == 0
  if (any(zero)) {   # never-recorded species carry no information and break IndVal tables
    notes <- c(notes, sprintf("%d species with no records were ignored.", sum(zero)))
    Y <- Y[, !zero, drop = FALSE]
  }
  presence <- isTRUE(st$binary) || identical(data$data_type, "presence")
  Yb <- if (presence) (Y > 0) * 1 else Y
  n <- nrow(Y); level <- st$ci
  adjusted <- !identical(st$r2, "unadjusted")
  rows <- list(); pv <- c()
  # --- (1) RDA + variation partitioning, jackknife intervals
  Yt <- com_transform(Yb, st$transform)
  part <- data$partition
  if (!is.null(st$partition)) part <- st$partition
  if (is.null(part) || !length(part)) part <- stats::setNames(as.list(data$env_vars), data$env_vars)
  part <- part[vapply(part, length, 1L) > 0]
  Xs <- lapply(part, function(v) com_design(env, unlist(v), st$poly %||% 1))
  vp <- if (length(Xs)) com_varpart_vec(Yt, Xs, adjusted) else c()
  if (length(vp)) {
    jk <- t(vapply(seq_len(n), function(i)
      com_varpart_vec(Yt[-i, , drop = FALSE], lapply(Xs, function(X) X[-i, , drop = FALSE]), adjusted),
      vp))
    if (length(vp) == 1) jk <- matrix(jk, ncol = 1, dimnames = list(NULL, names(vp)))
    for (k in names(vp)) rows[[k]] <- jack_rows(k, vp[[k]], jk[, k], level)
  }
  # --- (2) PERMANOVA effect size on the chosen distance, jackknife interval
  gname <- data$group
  D <- vegan::vegdist(Yb, method = st$distance, binary = presence)
  if (!is.null(gname)) {
    gr <- droplevels(factor(env[[gname]]))
    pj <- permanova_jack(as.matrix(D)^2, gr, adjusted)
    nm <- paste0("R2adj_", gname, "_permanova")
    rows[[nm]] <- jack_rows(nm, pj$estimate, pj$jack, level)
  }
  # --- (3) Baselga beta components (presence/absence), jackknife intervals
  bp <- beta_pair_mats(Y)
  js <- lapply(bp, pair_mean_jack)
  for (k in c("sor", "sim", "sne")) {
    nm <- paste0("beta_", k, "_pair_mean"); rows[[nm]] <- jack_rows(nm, js[[k]]$estimate, js[[k]]$jack, level)
  }
  rows$turnover_share <- jack_rows("turnover_share", js$sim$estimate / js$sor$estimate,
                                   js$sim$jack / js$sor$jack, level)
  est <- do.call(rbind, rows)
  est$p_value <- NA_real_
  obj <- list(settings = st, presence = presence, beta_pairs = bp, D = D, partition = part)
  # --- permutation tests and extras (skipped in fast settings)
  perm <- as.integer(st$n_perm %||% 0)
  if (perm > 0 && length(Xs)) {
    ctrl <- permute::how(nperm = perm)
    Xall <- as.data.frame(do.call(cbind, Xs))
    names(Xall) <- make.names(names(Xall), unique = TRUE)
    setcols <- split(names(Xall), rep(names(Xs), vapply(Xs, ncol, 1L)))
    r_all <- vegan::rda(Yt ~ ., data = Xall)
    est$p_value[est$param == "R2adj_env_total"] <- stats::anova(r_all, permutations = ctrl)[1, "Pr(>F)"]
    for (k in names(Xs)) {
      f1 <- stats::as.formula(paste("Yt ~", paste(setcols[[k]], collapse = " + ")))
      est$p_value[est$param == paste0("R2adj_", k)] <- stats::anova(vegan::rda(f1, data = Xall), permutations = ctrl)[1, "Pr(>F)"]
      if (length(Xs) > 1) {
        others <- unlist(setcols[setdiff(names(Xs), k)])
        f2 <- stats::as.formula(paste("Yt ~", paste(setcols[[k]], collapse = " + "),
                                      "+ Condition(", paste(others, collapse = " + "), ")"))
        est$p_value[est$param == paste0("R2adj_unique_", k)] <- stats::anova(vegan::rda(f2, data = Xall), permutations = ctrl)[1, "Pr(>F)"]
      }
    }
    obj$rda <- r_all
    obj$rda_axes <- tryCatch(stats::anova(r_all, by = "axis", permutations = permute::how(nperm = min(perm, 499))),
                             error = function(e) NULL)
    obj$varpart_vegan <- tryCatch(if (length(Xs) %in% 2:4) do.call(vegan::varpart, c(list(Yt), lapply(Xs, as.matrix))) else NULL,
                                  error = function(e) NULL)
  }
  if (perm > 0 && !is.null(gname)) {
    ctrl <- permute::how(nperm = perm)
    df <- data.frame(g = gr)
    a2 <- vegan::adonis2(D ~ g, data = df, permutations = ctrl)
    obj$permanova <- a2
    est$p_value[est$param == paste0("R2adj_", gname, "_permanova")] <- a2[1, "Pr(>F)"]
    bd <- vegan::betadisper(D, gr, type = st$dispersion_type %||% "centroid")
    bt <- vegan::permutest(bd, permutations = ctrl)
    obj$betadisper <- bd; obj$betadisper_test <- bt
    md <- tapply(bd$distances, gr, mean)
    jd <- jack_rows(paste0("dispersion_ratio_", gname), max(md) / min(md), rep(NA_real_, 2), level)
    jd$sd <- NA_real_; jd$lower <- NA_real_; jd$upper <- NA_real_
    jd$p_value <- bt$tab[1, "Pr(>F)"]
    est <- rbind(est, jd)
    if (nlevels(gr) > 2 && isTRUE(st$pairwise)) {
      prs <- utils::combn(levels(gr), 2)
      pw <- do.call(rbind, lapply(seq_len(ncol(prs)), function(k) {
        keep <- gr %in% prs[, k]
        Dk <- stats::as.dist(as.matrix(D)[keep, keep])
        ak <- vegan::adonis2(Dk ~ g, data = data.frame(g = droplevels(gr[keep])), permutations = ctrl)
        data.frame(pair = paste(prs[, k], collapse = " vs "), R2 = ak[1, "R2"], F = ak[1, "F"], p = ak[1, "Pr(>F)"])
      }))
      pw$p_holm <- stats::p.adjust(pw$p, "holm")
      obj$permanova_pairwise <- pw
    }
    if (length(Xs)) {
      Xall2 <- cbind(df, env[, unique(unlist(part)), drop = FALSE])
      fm <- stats::as.formula(paste("D ~ g +", paste(unique(unlist(part)), collapse = " + ")))
      obj$permanova_env <- tryCatch(vegan::adonis2(fm, data = Xall2, permutations = ctrl, by = "margin"),
                                    error = function(e) NULL)
    }
    if (isTRUE(st$indicator)) {
      ip <- tryCatch(indicspecies::multipatt(Yb, gr, func = "IndVal.g", duleg = TRUE,
                                             control = ctrl), error = function(e) NULL)
      if (!is.null(ip)) {
        sg <- ip$sign
        grp_names <- levels(gr)
        tab <- data.frame(species = rownames(sg),
                          group = grp_names[vapply(seq_len(nrow(sg)), function(i) { w <- which.max(unlist(sg[i, seq_along(grp_names)])); if (length(w)) w else NA_integer_ }, 1L)],
                          indval = sg$stat, p = sg$p.value, stringsAsFactors = FALSE)
        A <- ip$A; B <- ip$B
        tab <- tab[!is.na(tab$group), , drop = FALSE]
        tab$A <- vapply(seq_len(nrow(tab)), function(i) A[tab$species[i], tab$group[i]], 0)
        tab$B <- vapply(seq_len(nrow(tab)), function(i) B[tab$species[i], tab$group[i]], 0)
        tab$p_fdr <- stats::p.adjust(tab$p, "BH")
        obj$indicators <- tab[order(tab$p, -tab$indval), ]
      }
    }
  }
  # Ordinations
  conv <- TRUE
  if (isTRUE(st$nmds)) {
    nm <- tryCatch(suppressWarnings(vegan::metaMDS(Yb, distance = st$distance, k = st$nmds_k,
                                                   trymax = st$nmds_trymax, autotransform = FALSE,
                                                   binary = presence, trace = 0)),
                   error = function(e) { notes <<- c(notes, paste("NMDS failed:", conditionMessage(e))); NULL })
    if (!is.null(nm)) {
      obj$nmds <- nm
      if (!isTRUE(nm$converged)) notes <- c(notes, sprintf(
        "NMDS did not repeat its best solution in %d random starts (stress %.3f); treat the configuration as approximate.", st$nmds_trymax, nm$stress))
      est <- rbind(est, data.frame(param = "nmds_stress", estimate = nm$stress, sd = NA, lower = NA, upper = NA, p_value = NA))
      if (length(data$env_vars))
        obj$envfit <- tryCatch(vegan::envfit(nm, env[, data$env_vars, drop = FALSE],
                                             permutations = max(perm, 99)), error = function(e) NULL)
    }
  }
  if (isTRUE(st$dca)) {
    dc <- tryCatch(vegan::decorana(Yb), error = function(e) NULL)
    if (!is.null(dc)) {
      len <- diff(range(vegan::scores(dc, display = "sites")[, 1]))
      obj$dca <- dc
      est <- rbind(est, data.frame(param = "dca_axis1_length_sd", estimate = len, sd = NA, lower = NA, upper = NA, p_value = NA))
    }
    obj$pca <- vegan::rda(Yt)
  }
  if (perm > 0 && length(Xs)) {
    Xdf <- env[, unique(unlist(part)), drop = FALSE]
    fm <- stats::as.formula(paste("Yb ~", paste(names(Xdf), collapse = " + ")))
    if (isTRUE(st$cca)) {
      cc <- tryCatch(vegan::cca(fm, data = Xdf), error = function(e) NULL)
      if (!is.null(cc)) {
        ra <- vegan::RsquareAdj(cc, permutations = min(perm, 999))
        obj$cca <- cc
        est <- rbind(est, data.frame(param = "cca_R2adj_env_total", estimate = ra$adj.r.squared, sd = NA, lower = NA, upper = NA,
                                     p_value = stats::anova(cc, permutations = permute::how(nperm = perm))[1, "Pr(>F)"]))
      }
    }
    if (isTRUE(st$dbrda)) {
      fmd <- stats::as.formula(paste("D ~", paste(names(Xdf), collapse = " + ")))
      db <- tryCatch(vegan::dbrda(fmd, data = Xdf), error = function(e) NULL)
      if (!is.null(db)) {
        obj$dbrda <- db
        est <- rbind(est, data.frame(param = "dbrda_R2adj_env_total", estimate = vegan::RsquareAdj(db)$adj.r.squared,
                                     sd = NA, lower = NA, upper = NA,
                                     p_value = stats::anova(db, permutations = permute::how(nperm = perm))[1, "Pr(>F)"]))
      }
    }
  }
  if (isTRUE(st$multisite)) {
    bm <- unlist(betapart::beta.multi((Y > 0) * 1, index.family = "sorensen"))
    obj$beta_multi <- bm
    est <- rbind(est, data.frame(param = c("beta_SOR_multisite", "beta_SIM_multisite", "beta_SNE_multisite"),
                                 estimate = unname(bm[c("beta.SOR", "beta.SIM", "beta.SNE")]),
                                 sd = NA, lower = NA, upper = NA, p_value = NA))
  }
  if (isTRUE(st$species_curves) && length(data$env_vars)) obj$curves <- com_species_curves(Yb, env, data$env_vars, presence, level)
  rownames(est) <- NULL
  list(estimates = est,
       diagnostics = list(converged = conv, rhat_max = NA_real_, ess_min = NA_real_, notes = notes),
       object = obj)
}

com_default_settings <- function() list(
  transform = "hellinger", distance = "bray", binary = FALSE, poly = 1, r2 = "adjusted",
  ci = 0.95, n_perm = 999, nmds = TRUE, nmds_k = 2, nmds_trymax = 50, dca = TRUE, cca = TRUE,
  dbrda = TRUE, indicator = TRUE, pairwise = TRUE, multisite = TRUE, species_curves = TRUE,
  dispersion_type = "centroid", seed = 1)

# Per-species unimodal (quadratic) GLMs along each gradient, other gradients
# held at their means: the classical Gaussian-logit/log response curve
# (ter Braak & Looman 1986). Quasi-likelihood for overdispersed counts.
com_species_curves <- function(Yb, env, vars, presence, level = 0.95, n_top = 6) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  fam <- if (presence) stats::quasibinomial() else stats::quasipoisson()
  out <- list()
  for (v in vars) {
    others <- setdiff(vars[vapply(vars, function(w) is.numeric(env[[w]]), TRUE)], v)
    rhs <- paste(c(sprintf("poly(%s, 2, raw = TRUE)", v), others), collapse = " + ")
    res <- lapply(colnames(Yb), function(sp) {
      y <- Yb[, sp]
      if (sum(y > 0) < 5) return(NULL)
      d <- env; d$y <- y
      m <- tryCatch(suppressWarnings(stats::glm(stats::as.formula(paste("y ~", rhs)), family = fam, data = d)),
                    error = function(e) NULL)
      if (is.null(m)) return(NULL)
      m0 <- suppressWarnings(stats::glm(stats::as.formula(paste("y ~", if (length(others)) paste(others, collapse = "+") else "1")), family = fam, data = d))
      xs <- seq(min(env[[v]]), max(env[[v]]), length.out = 80)
      nd <- data.frame(xs); names(nd) <- v
      for (w in others) nd[[w]] <- mean(env[[w]])
      pr <- stats::predict(m, newdata = nd, se.fit = TRUE, type = "link")
      inv <- fam$linkinv
      list(dev_gain = (m0$deviance - m$deviance) / m0$null.deviance,
           pred = data.frame(x = xs, estimate = inv(pr$fit), lower = inv(pr$fit - z * pr$se.fit),
                             upper = inv(pr$fit + z * pr$se.fit), species = sp))
    })
    names(res) <- colnames(Yb)
    res <- Filter(Negate(is.null), res)
    if (!length(res)) next
    gain <- vapply(res, function(r) r$dev_gain, 0)
    top <- names(sort(gain, decreasing = TRUE))[seq_len(min(n_top, length(gain)))]
    out[[v]] <- list(pred = do.call(rbind, lapply(res[top], function(r) r$pred)), gain = sort(gain, decreasing = TRUE))
  }
  out
}
