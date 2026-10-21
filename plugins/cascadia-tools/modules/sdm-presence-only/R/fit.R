# fit() for sdm-presence-only: engine dispatch, estimates, derived
# quantities, spatial block cross-validation, quadrature check, GOF.

sdm_settings <- list(
  engine = "auto",            # auto (integrated if survey data, else ppm) | ppm | integrated | survey_only | maxnet
  bias_correction = TRUE,     # FALSE ignores the bias covariates (only for demonstration / negative controls)
  perfect_detection = FALSE,  # TRUE: treat survey detections as true presence/absence (p = 1)
  draws = 2000,               # parametric draws (ML engines) for derived quantities and bands
  evaluate = TRUE, cv_folds = 5, cv_block_size = NULL,   # spatial block CV (block side in coordinate units)
  quadrature_check = TRUE, gof = TRUE, gof_reps = 1000,
  maxnet_classes = "l", maxnet_regmult = 1, maxnet_boot = 100,
  compare_engines = character(), seed = 1)

sdm_fast <- merge_list(sdm_settings, list(draws = 500, evaluate = FALSE, quadrature_check = FALSE,
                                          gof = FALSE, maxnet_boot = 50))

sdm_mvn <- function(mu, V, n) {
  L <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(L)) { e <- eigen(V, symmetric = TRUE); L <- t(e$vectors %*% diag(sqrt(pmax(e$values, 0)))) }
  z <- matrix(stats::rnorm(n * length(mu)), n)
  sweep(z %*% L, 2, mu, "+")
}

sdm_run_engine <- function(engine, prep, settings) {
  switch(engine,
    ppm = sdm_engine_ppm(prep, settings),
    integrated = sdm_engine_ml(prep, settings, use_po = TRUE, use_survey = TRUE),
    survey_only = sdm_engine_ml(prep, settings, use_po = FALSE, use_survey = TRUE),
    maxnet = sdm_engine_maxnet(prep, settings),
    stop("unknown engine '", engine, "' (ppm, integrated, survey_only, maxnet)"))
}

# Draws of the parameter vector (ML: multivariate normal; maxnet: bootstrap)
sdm_draws <- function(eng, n) {
  if (!is.null(eng$boot)) return(eng$boot)
  if (is.null(eng$vcov)) return(matrix(eng$coef, 1, dimnames = list(NULL, names(eng$coef))))
  d <- sdm_mvn(eng$coef, eng$vcov, n); colnames(d) <- names(eng$coef); d
}

# habitat linear predictor (log intensity relative to mean covariates) for a
# matrix of standardised habitat terms; returns draws x rows (row 1 = point estimate
# when est = TRUE)
sdm_hab_link <- function(obj, Xnew, Znew = NULL, draws = TRUE, bias = FALSE) {
  eng <- obj$eng
  if (!is.null(eng$betas)) {   # maxnet
    nd <- as.data.frame(Xnew[, obj$prep$habitat, drop = FALSE])
    for (b in obj$prep$bias) nd[[b]] <- if (bias && !is.null(Znew)) Znew[, b] else 0
    B <- if (draws) rbind(eng$betas, eng$boot_betas) else matrix(eng$betas, 1, dimnames = list(NULL, names(eng$betas)))
    base <- nd[1, , drop = FALSE]; base[] <- 0
    out <- t(apply(B, 1, function(b) { b <- b[b != 0]; if (!length(b)) return(rep(0, nrow(nd)))
      sdm_maxnet_link(b, nd) - sdm_maxnet_link(b, base) }))
    return(out)
  }
  bn <- paste0("intensity_beta_", colnames(obj$prep$X))
  D <- if (draws) rbind(eng$coef, obj$draws) else matrix(eng$coef, 1, dimnames = list(NULL, names(eng$coef)))
  out <- D[, bn, drop = FALSE] %*% t(Xnew)
  if (bias && ncol(obj$prep$Z)) out <- out + D[, paste0("bias_beta_", colnames(obj$prep$Z)), drop = FALSE] %*% t(Znew)
  if ("log_density" %in% colnames(D) && !bias) out <- out + D[, "log_density"]
  out
}

sdm_fit <- function(data, settings) {
  settings <- merge_list(sdm_settings, settings)
  set.seed(settings$seed)
  prep <- sdm_prepare(data, settings)
  engine <- settings$engine
  if (engine == "auto") engine <- if (!is.null(prep$survey)) "integrated" else "ppm"
  if (prep$n_pres == 0 && engine != "survey_only") stop("no presence records inside the study region")
  eng <- sdm_run_engine(engine, prep, settings)
  draws <- sdm_draws(eng, settings$draws)
  obj <- list(engine = engine, eng = eng, prep = prep, draws = draws, settings = settings)
  # ---- estimates table
  shown <- grep("^(intensity_beta_|bias_beta_)", names(eng$coef), value = TRUE)
  if (!is.null(eng$boot)) {
    est <- data.frame(param = shown, estimate = unname(eng$coef[shown]),
                      sd = apply(draws[, shown, drop = FALSE], 2, stats::sd),
                      lower = apply(draws[, shown, drop = FALSE], 2, stats::quantile, 0.025),
                      upper = apply(draws[, shown, drop = FALSE], 2, stats::quantile, 0.975), row.names = NULL)
  } else {
    se <- if (is.null(eng$vcov)) rep(NA_real_, length(shown)) else sqrt(pmax(diag(eng$vcov)[shown], 0))
    est <- wald_rows(eng$coef[shown], se, shown)
  }
  if ("log_density" %in% names(eng$coef)) {
    A <- prep$A; Amed <- stats::median(prep$survey$A)
    der <- function(th) {
      th <- stats::setNames(th, names(eng$coef))
      lam <- exp(th[["log_density"]] + drop(prep$X %*% th[paste0("intensity_beta_", colnames(prep$X))]))
      c(density_at_mean_covariates = exp(th[["log_density"]]),
        expected_total_abundance = sum(A * lam),
        mean_occupancy_prob = mean(-expm1(-lam * Amed)),
        detection_prob = if ("logit_detection" %in% names(th)) stats::plogis(th[["logit_detection"]]) else NA)
    }
    pt <- der(eng$coef); dd <- t(apply(draws, 1, der))
    keep <- names(pt)[is.finite(pt)]
    est <- rbind(est, data.frame(param = keep, estimate = unname(pt[keep]),
                                 sd = apply(dd[, keep, drop = FALSE], 2, stats::sd),
                                 lower = apply(dd[, keep, drop = FALSE], 2, stats::quantile, 0.025),
                                 upper = apply(dd[, keep, drop = FALSE], 2, stats::quantile, 0.975), row.names = NULL))
    obj$survey_area_km2 <- Amed
  }
  # ---- diagnostics
  notes <- prep$notes
  conv <- isTRUE(eng$converged)
  # boundary / non-identification check on every link-scale parameter except the
  # PO intercept (whose level is arbitrary): |estimate| > 10 or SE > 5 on the log/logit
  # scale means a 150-fold multiplicative uncertainty, i.e. not estimable from these data
  chk <- setdiff(names(eng$coef), "po_intercept")
  se_chk <- if (!is.null(eng$vcov)) sqrt(pmax(diag(eng$vcov)[chk], 0)) else apply(draws[, intersect(chk, colnames(draws)), drop = FALSE], 2, stats::sd)[chk]
  big <- stats::setNames(abs(eng$coef[chk]) > 10 | (!is.na(se_chk) & se_chk > 5), chk)
  if (any(big)) { conv <- FALSE; notes <- c(notes, paste("Estimates on a boundary or not identified (|estimate| > 10 or SE > 5 on the link scale):", paste(chk[big], collapse = ", "))) }
  if (!conv && !any(big)) notes <- c(notes, "The optimiser did not converge or the information matrix is singular: estimates are not reliable.")
  if (engine == "maxnet") notes <- c(notes, sprintf(
    "maxnet: intervals are percentile bootstrap intervals over %d resamples of the presence records; maxnet shrinks coefficients towards zero (regularisation multiplier %g), so effects are conservative.",
    nrow(eng$boot), settings$maxnet_regmult),
    if (grepl("[hpt]", eng$classes)) sprintf("maxnet features '%s' include hinge/product/threshold terms: the table reports only the linear (and quadratic) coefficients, which are NOT the full effect of a covariate; read the response curves instead.", eng$classes))
  # ---- naive comparison: what if the bias covariate were ignored?
  if (length(prep$bias) && engine %in% c("ppm", "integrated", "maxnet")) {
    s2 <- settings; s2$bias_correction <- FALSE
    p2 <- sdm_prepare(data, s2)
    e2 <- tryCatch(if (engine == "integrated") sdm_engine_ml(p2, s2, TRUE, TRUE) else sdm_engine_ppm(p2, s2), error = function(e) NULL)
    if (!is.null(e2)) obj$naive <- e2$coef
  }
  # ---- extras for real-data reports
  if (isTRUE(settings$gof) && engine != "survey_only") {
    obj$gof <- tryCatch(sdm_gof(obj, settings), error = function(e) { message("gof failed: ", conditionMessage(e)); NULL })
    if (!is.null(obj$gof)) notes <- c(notes, sprintf(
      "Goodness of fit: records counted in %d spatial blocks vs model expectation, Pearson dispersion %.2f (1 = Poisson; P(replicate >= observed) = %.2f from %d parametric replicates).%s",
      obj$gof$nblocks, obj$gof$chat, obj$gof$p, length(obj$gof$rep),
      if (obj$gof$p < 0.05) " The records are more clustered than the model allows (unmodelled habitat, residual sampling bias or duplicate records), so the intervals are too narrow; treat them as a lower bound on uncertainty." else ""))
  }
  if (isTRUE(settings$quadrature_check) && engine != "survey_only") {
    obj$quad <- tryCatch(sdm_quadrature_check(obj, settings), error = function(e) { message("quadrature check failed: ", conditionMessage(e)); NULL })
    if (!is.null(obj$quad)) {
      qt <- obj$quad$table; qt <- qt[qt$spacing > 1, ]
      mx <- tapply(abs(qt$change_se), qt$spacing, max)
      mx <- mx[order(as.numeric(names(mx)))]
      notes <- c(notes, sprintf(
        "Quadrature check (Renner et al. 2015): largest change in any coefficient when the integral uses only %s, relative to the full grid of %d cells used for the estimates.%s",
        paste(sprintf("every %s cell: %.2f SE", c("1" = "", "2" = "2nd", "4" = "4th", "8" = "8th")[names(mx)], mx), collapse = "; "), nrow(prep$X),
        if (length(mx) && mx[1] > 0.5) " Even halving the quadrature changes the estimates by more than 0.5 SE: use a finer covariate grid." else " Changes shrink towards 0 as the grid is refined, so the full grid approximates the integral well."))
    }
  }
  if (isTRUE(settings$evaluate) && engine != "survey_only") {
    obj$cv <- tryCatch(sdm_block_cv(data, obj, settings), error = function(e) { message("block CV failed: ", conditionMessage(e)); NULL })
    if (!is.null(obj$cv)) notes <- c(notes, sdm_cv_note(obj$cv))
  }
  if (length(settings$compare_engines)) {
    cmp <- list(data.frame(engine = engine, est))
    for (e2 in setdiff(unlist(settings$compare_engines), engine)) {
      s2 <- settings; s2$engine <- e2; s2$compare_engines <- character(); s2$evaluate <- FALSE
      s2$quadrature_check <- FALSE; s2$gof <- FALSE
      r2 <- tryCatch(sdm_fit(data, s2), error = function(e) { message("compare engine ", e2, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r2)) cmp[[length(cmp) + 1]] <- data.frame(engine = e2, r2$estimates)
    }
    obj$comparison <- do.call(rbind, cmp)
  }
  dg <- list(converged = conv, notes = notes, loglik = eng$loglik)
  list(estimates = est, diagnostics = dg, object = obj)
}

# ---- spatial blocks ----------------------------------------------------------
sdm_blocks <- function(grid, size = NULL, k = 5) {
  ext <- max(diff(range(grid$x)), diff(range(grid$y)))
  if (is.null(size)) size <- ext / 5
  bx <- floor((grid$x - min(grid$x)) / size); by <- floor((grid$y - min(grid$y)) / size)
  as.integer(factor(paste(bx, by)))
}

sdm_subset_prep <- function(prep, cells_keep) {
  p2 <- prep
  p2[["quad_idx"]] <- which(cells_keep); p2[["quad_w"]] <- prep$A[which(cells_keep)]
  p2$pres_cell <- prep$pres_cell[cells_keep[prep$pres_cell]]
  p2$n_pres <- length(p2$pres_cell)
  if (!is.null(prep$survey)) {
    k <- cells_keep[prep$survey$cell]
    s <- prep$survey
    p2$survey <- list(cell = s$cell[k], A = s$A[k], ydet = s$ydet[k], J = s$J[k], site = s$site[k],
                      y = s$y[k, , drop = FALSE], X = s$X[k, , drop = FALSE], estimate_p = s$estimate_p)
    if (!length(p2$survey$cell)) p2$survey <- NULL
  }
  p2
}

sdm_auc <- function(pres, bg) {
  r <- rank(c(pres, bg)); np <- length(pres); nb <- length(bg)
  (sum(r[seq_len(np)]) - np * (np + 1) / 2) / (np * nb)
}

# Continuous Boyce index (Hirzel et al. 2006): moving window of width 10% of
# the prediction range, predicted-to-expected ratio vs window, Spearman rho.
sdm_boyce <- function(pres, bg, nwin = 100, width = 0.1) {
  r <- range(c(pres, bg)); w <- width * diff(r)
  mids <- seq(r[1] + w / 2, r[2] - w / 2, length.out = nwin)
  P <- vapply(mids, function(m) mean(pres >= m - w / 2 & pres <= m + w / 2), 0)
  E <- vapply(mids, function(m) mean(bg >= m - w / 2 & bg <= m + w / 2), 0)
  ok <- E > 0
  f <- P[ok] / E[ok]
  rho <- if (sum(ok) > 2 && stats::sd(f) > 0) suppressWarnings(stats::cor(mids[ok], f, method = "spearman")) else NA_real_
  list(index = rho, curve = data.frame(mid = mids[ok], pe = f))
}

sdm_block_cv <- function(data, obj, settings) {
  prep <- obj$prep
  k <- settings$cv_folds
  blk <- sdm_blocks(prep$grid, settings$cv_block_size, k)
  set.seed(settings$seed + 7)
  ub <- unique(blk); fold_of_block <- stats::setNames(sample(rep(seq_len(k), length.out = length(ub))), ub)
  fold <- fold_of_block[as.character(blk)]
  s2 <- settings; s2$maxnet_boot <- 0
  rows <- list(); curves <- list()
  for (f in seq_len(k)) {
    test <- fold == f
    tp <- prep$pres_cell[test[prep$pres_cell]]
    if (length(tp) < 3 || sum(!test) < 10) next
    ptr <- sdm_subset_prep(prep, !test)
    eng_name <- if (obj$engine == "integrated" && is.null(ptr$survey)) "ppm" else obj$engine
    e <- tryCatch(sdm_run_engine(eng_name, ptr, s2), error = function(e) NULL)
    if (is.null(e)) next
    o2 <- list(eng = e, prep = prep)
    tc <- which(test)
    hab <- drop(sdm_hab_link(o2, prep$X[tc, , drop = FALSE], prep$Z[tc, , drop = FALSE], draws = FALSE))
    rec <- if (ncol(prep$Z)) drop(sdm_hab_link(o2, prep$X[tc, , drop = FALSE], prep$Z[tc, , drop = FALSE], draws = FALSE, bias = TRUE)) else hab
    pos <- match(tp, tc)
    bh <- sdm_boyce(hab[pos], hab); br <- sdm_boyce(rec[pos], rec)
    sv_auc <- NA_real_
    if (!is.null(prep$survey)) {
      ts <- which(test[prep$survey$cell])
      yy <- prep$survey$ydet[ts] > 0
      if (length(ts) >= 5 && any(yy) && any(!yy)) {
        sh <- drop(sdm_hab_link(o2, prep$survey$X[ts, , drop = FALSE], draws = FALSE))
        sv_auc <- sdm_auc(sh[yy], sh[!yy])
      }
    }
    rows[[length(rows) + 1]] <- data.frame(fold = f, n_test_records = length(tp), n_test_cells = length(tc),
                                           auc_habitat = sdm_auc(hab[pos], hab), boyce_habitat = bh$index,
                                           auc_records = sdm_auc(rec[pos], rec), boyce_records = br$index,
                                           auc_survey = sv_auc)
    curves[[length(curves) + 1]] <- data.frame(fold = f, bh$curve)
  }
  if (!length(rows)) return(NULL)
  list(folds = do.call(rbind, rows), curves = do.call(rbind, curves), nblocks = length(ub), block = blk, fold = fold)
}

sdm_cv_note <- function(cv) {
  f <- cv$folds
  s <- sprintf("Spatial block cross-validation (%d folds of %d blocks): mean test AUC %.2f and Boyce index %.2f for the bias-corrected habitat map; AUC %.2f and Boyce %.2f for predicting where records occur (habitat x sampling effort).",
               nrow(f), cv$nblocks, mean(f$auc_habitat), mean(f$boyce_habitat, na.rm = TRUE),
               mean(f$auc_records), mean(f$boyce_records, na.rm = TRUE))
  if (any(is.finite(f$auc_survey))) s <- paste(s, sprintf("Against held-out survey sites (detected vs not), AUC %.2f.", mean(f$auc_survey, na.rm = TRUE)))
  s
}

# ---- quadrature convergence (Renner et al. 2015) ---------------------------
sdm_quadrature_check <- function(obj, settings) {
  prep <- obj$prep
  sp <- sdm_grid_spacing(prep$grid)
  ix <- round((prep$grid$x - min(prep$grid$x)) / sp["dx"]); iy <- round((prep$grid$y - min(prep$grid$y)) / sp["dy"])
  s2 <- settings; s2$maxnet_boot <- 0
  out <- list()
  for (k in c(8, 4, 2, 1)) {
    keep <- (ix %% k == 0) & (iy %% k == 0)
    if (sum(keep) < 50 && k > 1) next
    p2 <- prep; p2[["quad_idx"]] <- which(keep); p2[["quad_w"]] <- rep(sum(prep$A) / sum(keep), sum(keep))
    e <- tryCatch(sdm_run_engine(obj$engine, p2, s2), error = function(e) NULL)
    if (is.null(e)) next
    nm <- grep("^(intensity_beta_|bias_beta_|log_density)", names(e$coef), value = TRUE)
    out[[length(out) + 1]] <- data.frame(spacing = k, n_quadrature = sum(keep), param = nm, estimate = unname(e$coef[nm]))
  }
  d <- do.call(rbind, out)
  full <- obj$eng$coef
  se <- if (!is.null(obj$eng$vcov)) sqrt(diag(obj$eng$vcov)) else apply(obj$draws, 2, stats::sd)
  d$change_se <- (d$estimate - full[d$param]) / se[d$param]
  sp2 <- sort(unique(d$spacing))
  last <- if (length(sp2) > 1) max(abs(d$change_se[d$spacing == sp2[2]])) else NA_real_
  list(table = d, spacing = sort(unique(d$spacing), decreasing = TRUE), max_change_se_last = last)
}

# ---- goodness of fit: record counts per spatial block ----------------------
sdm_gof <- function(obj, settings) {
  prep <- obj$prep
  blk <- sdm_blocks(prep$grid, settings$cv_block_size)
  lin <- drop(sdm_hab_link(obj, prep$X, prep$Z, draws = FALSE, bias = TRUE))
  mu <- prep$A * exp(lin - max(lin))
  share <- tapply(mu, blk, sum); share <- share / sum(share)
  n <- prep$n_pres
  O <- tabulate(blk[prep$pres_cell], nbins = length(share))
  E <- n * share
  X2 <- function(o) sum((o - E)^2 / E)
  obs <- X2(O)
  set.seed(settings$seed + 11)
  rep <- vapply(seq_len(settings$gof_reps), function(i) X2(stats::rmultinom(1, n, share)[, 1]), 0)
  npar <- length(grep("^(intensity_beta_|bias_beta_)", names(obj$eng$coef))) + 1
  list(obs = obs, rep = rep, p = mean(rep >= obs), chat = obs / max(length(share) - npar, 1),
       nblocks = length(share), blocks = data.frame(block = seq_along(share), observed = O, expected = E))
}
