# Hierarchical distance sampling with unmarked (maximum likelihood).
#   hds: distsamp (Royle, Dawson & Bates 2004). Counts per transect and
#        distance bin; N_i ~ Poisson(D_i a_i), log D_i = X_i beta;
#        y_ij ~ Poisson(D_i a_i pi_ij(sigma_i)), log sigma_i = Z_i alpha.
#   gds: gdistsamp (Chandler, Royle & King 2011) for repeated visits:
#        M_i ~ Poisson(D_i a_i), N_it ~ Bin(M_i, phi), y_ijt ~ Multinom(N_it, pi_ij).
# Derived quantities (density at an average site, mean density over the
# surveyed transects or a prediction grid, abundance, average detection
# probability) are computed from multivariate-normal draws of the MLEs
# (parametric bootstrap of the estimates; Lande et al. 2003 style).

ds_ukey <- c(hn = "halfnorm", hr = "hazard")

# TRUE when a point survey has effort other than one visit per point (HDS needs an offset).
ds_point_offset <- function(data) data$survey == "point" && any(abs(data$transects$effort * data$n_visits - 1) > 1e-9)

ds_umf <- function(data, cut, gds = FALSE) {
  tr <- data$tz
  sc <- tr[, unique(c(data$st_terms, data$dt_terms)), drop = FALSE]
  if (!ncol(sc)) sc <- data.frame(.dummy = rep(1, nrow(tr)))
  nV <- data$n_visits
  # unmarked takes the area of a point as pi w^2 whatever tlength says, so point
  # effort other than 1 (several visits pooled per point) enters as an offset
  # log(effort) on density (ds_point_offset); found in independent verification.
  if (!gds && ds_point_offset(data)) sc$.log_eff <- log(data$transects$effort * nV)
  if (!gds) {
    y <- ds_bin_counts(data, cut)
    unmarked::unmarkedFrameDS(y = y, siteCovs = sc, dist.breaks = cut,
                              tlength = data$transects$effort * 1000 * nV, survey = data$survey, unitsIn = "m")
  } else {
    y <- ds_bin_counts(data, cut, by_visit = TRUE)
    unmarked::unmarkedFrameGDS(y = y, siteCovs = sc, dist.breaks = cut, numPrimary = nV,
                               tlength = data$transects$effort * 1000, survey = data$survey, unitsIn = "m")
  }
}

ds_engine_unmarked <- function(data, settings, gds = FALSE) {
  eco_require("unmarked")
  cut <- ds_hds_cut(data, settings)
  w <- max(cut)
  umf <- ds_umf(data, cut, gds)
  if (!data$dt_site_level) stop("hierarchical distance sampling needs detection covariates measured per transect (columns of the transect file)")
  fst <- stats::as.formula(ds_rhs(data$formulas$density))
  fdt <- stats::as.formula(ds_rhs(data$formulas$detection))
  y <- if (gds) ds_bin_counts(data, cut, by_visit = TRUE) else NULL
  a_i <- ds_cov_area(data, w) / (if (gds) data$n_visits else 1)  # km2 covered per visit (gds) / in total (hds)
  fits <- list()
  st_txt <- ds_rhs(data$formulas$density)
  if (!gds && ds_point_offset(data)) st_txt <- paste(st_txt, "+ offset(.log_eff)")
  nbs <- ncol(stats::model.matrix(fst, data$tz)); nas <- ncol(stats::model.matrix(fdt, data$tz))
  for (k in settings$keys) {
    if (!k %in% names(ds_ukey)) next
    fm <- tryCatch(suppressWarnings(if (!gds) {
      # unmarked's default start (log sigma = 0, i.e. 1 m) sometimes ends on a flat
      # detection function (sigma -> 10 km, NaN SEs; 8/60 exact-distance point data sets in
      # independent verification). Also fit from data-based starts and keep the better fit.
      ff <- stats::as.formula(paste(ds_rhs(data$formulas$detection), st_txt))
      s0 <- c(log(max(sum(umf@y) / sum(a_i) / 0.5, 0.01)), rep(0, nbs - 1), log(w / 2), rep(0, nas - 1), if (k == "hr") log(2.5))
      cand <- list(tryCatch(unmarked::distsamp(ff, umf, keyfun = ds_ukey[[k]], output = "density", unitsOut = "kmsq", starts = s0), error = function(e) NULL),
                   tryCatch(unmarked::distsamp(ff, umf, keyfun = ds_ukey[[k]], output = "density", unitsOut = "kmsq"), error = function(e) NULL))
      cand <- Filter(Negate(is.null), cand)
      if (!length(cand)) stop("distsamp failed")
      cand[[which.min(vapply(cand, function(f) f@negLogLike, 1))]]
    } else {
      maxc <- max(vapply(seq_len(data$n_visits), function(t) max(rowSums(y[, (t - 1) * (length(cut) - 1) + seq_len(length(cut) - 1), drop = FALSE])), 1))
      K <- settings$K %||% NA
      if (is.na(K)) {
        # start generous, then check the Poisson tail at the fitted densities
        K <- max(50, ceiling(maxc * 6 + 20))
      }
      for (it in 1:3) {
        fm0 <- unmarked::gdistsamp(fst, ~ 1, fdt, umf, keyfun = ds_ukey[[k]], output = "density", unitsOut = "kmsq",
                                   mixture = "P", K = K)
        lam <- exp(stats::model.matrix(fst, data$tz) %*% fm0@estimates@estimates$lambda@estimates)
        if (is.null(settings$K) && any(stats::ppois(K, max(lam * a_i), lower.tail = FALSE) > 1e-6)) K <- 2 * K else break
      }
      attr(fm0, "K") <- K
      fm0
    }), error = function(e) e)
    if (inherits(fm, "error")) next
    fits[[k]] <- fm
  }
  if (!length(fits)) stop("unmarked ", if (gds) "gdistsamp" else "distsamp", " failed for every key")
  # Stability screen (as for Distance::ds): drop a key whose fit has a link-scale SE > 3 or a
  # non-invertible Hessian when another key fitted stably.
  # Also (same rule as for Distance::ds) drop a key whose average detection probability
  # has CV > 50%: a hazard-rate fitted as a near-spike to coarse bins can pass the SE
  # screen yet leave P_a, and hence density, essentially unidentified (seen in point
  # transects with ~40 detections: density estimates 100-1000x the truth).
  stable <- vapply(names(fits), function(k) { f <- fits[[k]]
    v <- tryCatch(diag(unmarked::vcov(f)), error = function(e) NA)
    all(is.finite(v)) && all(v > 0) && all(sqrt(v) < 3) &&
      ds_hds_pa_cv(f, k, w, data$survey) <= (settings$max_pa_cv %||% 0.5) }, TRUE)
  hds_dropped <- character()
  if (any(stable) && !all(stable)) { hds_dropped <- names(fits)[!stable]; fits <- fits[stable] }
  aic <- vapply(fits, function(f) f@AIC, 1)
  best <- names(which.min(aic))
  fm <- fits[[best]]
  E <- fm@estimates@estimates
  sname <- if (gds) "lambda" else "state"
  beta <- E[[sname]]@estimates; alpha <- E$det@estimates
  shape <- if (best == "hr") E$scale@estimates else NULL
  phi <- if (gds) E$phi@estimates else NULL
  th <- c(beta, alpha, shape, phi)
  nb <- length(beta); na <- length(alpha)
  V <- tryCatch(unmarked::vcov(fm), error = function(e) matrix(NA_real_, length(th), length(th)))
  # unmarked's vcov order: state/lambda, [phi], det, [scale] for gds; state, det, [scale] for distsamp
  vn <- colnames(V)
  ord <- c(grep("^lam", vn), grep("^p\\(", vn)[!grepl("scale", vn[grep("^p\\(", vn)])],
           grep("scale", vn), grep("^phi\\(", vn))
  if (length(ord) != length(th)) stop("internal: could not match unmarked vcov names: ", paste(vn, collapse = ", "))
  V <- V[ord, ord, drop = FALSE]
  dimnames(V) <- list(NULL, NULL)
  se <- sqrt(pmax(diag(V), 0))
  conv <- fm@opt$convergence == 0 && all(is.finite(diag(V))) && all(diag(V) > 0) && all(se < 5)
  notes <- character()
  if (!conv) notes <- c(notes, "unmarked fit: optimiser did not converge, Hessian not invertible, or a link-scale SE > 5 (parameters poorly identified).")
  if (length(hds_dropped)) notes <- c(notes, sprintf("Hierarchical model: %s key dropped (parameters not identified).", paste(vapply(hds_dropped, function(k) c(hn = "half-normal", hr = "hazard-rate")[[k]], ""), collapse = ", ")))

  Xs <- stats::model.matrix(fst, data$tz)
  Xd <- stats::model.matrix(fdt, data$tz)
  gs <- ds_group_size(data)
  nd <- settings$n_draws
  dr <- ds_mvn(th, V, nd)
  ib <- seq_len(nb); ia <- nb + seq_len(na)
  ish <- if (!is.null(shape)) nb + na + 1 else NULL
  iph <- if (!is.null(phi)) length(th) else NULL
  sz <- if (isTRUE(data$clusters)) pmax(stats::rnorm(nd, gs$mean, gs$se), 1) else rep(1, nd)
  # group density per transect, per draw (nd x n_transects)
  Dg <- exp(dr[, ib, drop = FALSE] %*% t(Xs))
  sig <- exp(dr[, ia, drop = FALSE] %*% t(Xd))
  bsh <- if (!is.null(ish)) exp(dr[, ish]) else NULL
  Pa <- matrix(0, nd, nrow(Xd))
  # P_a only varies with detection covariates; compute for unique rows
  ukey <- unique(Xd)
  for (u in seq_len(nrow(ukey))) {
    cols <- which(apply(Xd, 1, function(r) all(r == ukey[u, ])))
    pa_u <- ds_pa(best, sig[, cols[1]], w, data$survey, bsh)
    Pa[, cols] <- pa_u
  }
  wts <- Dg * matrix(a_i, nd, length(a_i), byrow = TRUE)
  pa_avg <- rowSums(Pa * wts) / rowSums(wts)
  # average site: numeric covariates at their mean (0 after standardising), factors at reference level
  D0 <- exp(dr[, 1]) * sz
  out <- list(density_at_average_site = D0)
  if (isTRUE(data$clusters)) out$group_density_at_average_site <- exp(dr[, 1])
  # density over the study area: prediction grid if supplied, else per-stratum mean of transects
  if (!is.null(data$grid)) {
    G <- data$grid
    Xg <- stats::model.matrix(stats::delete.response(stats::terms(fst)), G$z)
    Dgrid <- exp(dr[, ib, drop = FALSE] %*% t(Xg)) * sz
    Ntot <- Dgrid %*% G$area_km2
    out$abundance <- as.numeric(Ntot)
    out$density_per_km2 <- as.numeric(Ntot) / sum(G$area_km2)
  } else {
    str <- data$transects$stratum
    if (!is.null(data$strata)) {
      Ntot <- 0
      for (s in data$strata$stratum) Ntot <- Ntot + rowMeans(Dg[, str == s, drop = FALSE]) * sz * data$strata$area_km2[data$strata$stratum == s]
      out$abundance <- Ntot
      out$density_per_km2 <- Ntot / sum(data$strata$area_km2)
    } else out$density_per_km2 <- rowMeans(Dg) * sz
  }
  if (isTRUE(data$clusters)) out$group_density_per_km2 <- out$density_per_km2 / sz
  out$detection_prob <- pa_avg
  out[[if (data$survey == "line") "esw_m" else "edr_m"]] <- ds_esw(pa_avg, w, data$survey)
  if (best == "hn") out$sigma_m <- exp(dr[, ia[1]]) else { out$hr_scale_m <- exp(dr[, ia[1]]); out$hr_shape <- bsh }
  if (gds) out$availability_prob <- stats::plogis(dr[, iph])
  est <- ds_draw_rows(do.call(cbind, out), names(out))
  # covariate effects: Wald on the link scale (+ rate ratios)
  cn <- colnames(Xs)
  for (j in seq_len(nb)[-1]) {
    nmj <- ds_clean(cn[j])
    est <- rbind(est, wald_rows(beta[j], se[j], paste0("density_beta_", nmj)),
                 wald_rows(beta[j], se[j], paste0("rate_ratio_", nmj), transform = exp))
  }
  dn <- colnames(Xd)
  for (j in seq_len(na)[-1]) est <- rbind(est, wald_rows(alpha[j], se[nb + j], paste0("detection_beta_", ds_clean(dn[j]))))
  # Freeman-Tukey parametric bootstrap GOF (optional)
  gof <- NULL
  if (settings$gof_nboot > 0) gof <- tryCatch({
    ft <- function(f) { o <- unmarked::getY(f@data); e <- unmarked::fitted(f); sum((sqrt(o) - sqrt(e))^2, na.rm = TRUE) }
    pb <- suppressWarnings(unmarked::parboot(fm, statistic = ft, nsim = settings$gof_nboot))
    list(type = "freeman_tukey", stat = pb@t0[[1]], rep = as.numeric(pb@t.star[, 1]), p = mean(pb@t.star[, 1] >= pb@t0[[1]]))
  }, error = function(e) NULL)
  list(engine = if (gds) "unmarked::gdistsamp" else "unmarked::distsamp", estimates = est, fit = fm, key = best,
       selection = data.frame(model = paste(c(hn = "half-normal", hr = "hazard-rate")[names(aic)], "(binned HDS)"),
                              AIC = aic, dAIC = aic - min(aic), row.names = NULL),
       theta = th, vcov = V, nb = nb, na = na, draws = dr, size_draws = sz, cut = cut,
       gof = gof, converged = conv, notes = notes, K = attr(fm, "K"))
}

# CV of the average detection probability (at the reference level of detection
# covariates) of an unmarked distsamp/gdistsamp fit, from 400 multivariate-normal
# draws of the detection parameters. Inf when the Hessian is unusable.
ds_hds_pa_cv <- function(f, key, w, survey) {
  V <- tryCatch(unmarked::vcov(f), error = function(e) NULL)
  if (is.null(V)) return(Inf)
  vn <- colnames(V)
  ip <- grep("^p\\(", vn)
  isc <- grep("scale", vn)
  ia <- setdiff(ip, isc)[1]
  E <- f@estimates@estimates
  mu <- E$det@estimates[1]
  idx <- ia
  if (key == "hr") { mu <- c(mu, E$scale@estimates[1]); idx <- c(ia, isc[1]) }
  Vs <- V[idx, idx, drop = FALSE]
  if (anyNA(idx) || !all(is.finite(Vs)) || any(diag(Vs) <= 0)) return(Inf)
  dr <- ds_mvn(stats::setNames(mu, NULL), Vs, 400)
  pa <- ds_pa(key, exp(dr[, 1]), w, survey, if (key == "hr") exp(dr[, 2]) else NULL)
  stats::sd(pa) / mean(pa)
}
