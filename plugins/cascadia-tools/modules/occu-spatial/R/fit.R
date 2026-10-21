# occu-spatial: fitting through spOccupancy (Doser et al. 2022) and turning
# its output into ecotools estimands with R-hat / bulk ESS per estimand.

clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

resolve_model <- function(model, type) {
  if (is.null(model) || identical(model, "auto"))
    model <- switch(type, single = "spPGOcc", multi = "sfMsPGOcc", integrated = "intPGOcc")
  ok <- list(single = c("spPGOcc", "PGOcc"), multi = c("msPGOcc", "sfMsPGOcc"),
             integrated = c("intPGOcc", "spIntPGOcc"))[[type]]
  if (!model %in% ok)
    stop(sprintf("settings$model = '%s' does not fit %s data; use one of: %s", model, type,
                 paste(ok, collapse = ", ")))
  model
}

is_spatial_model <- function(model) model %in% c("spPGOcc", "sfMsPGOcc", "spIntPGOcc")

# stacked-chain draws matrix -> iterations x chains x variables array
chains_array <- function(M, n_chains) {
  M <- as.matrix(M); n <- nrow(M) %/% n_chains
  arr <- array(NA_real_, c(n, n_chains, ncol(M)), dimnames = list(NULL, NULL, colnames(M)))
  for (k in seq_len(n_chains)) arr[, k, ] <- M[(k - 1) * n + seq_len(n), , drop = FALSE]
  arr
}

summarise_estimands <- function(D, n_chains, prob = 0.95) {
  eco_require("posterior")
  arr <- chains_array(D, n_chains)
  a <- (1 - prob) / 2
  rows <- lapply(seq_len(dim(arr)[3]), function(v) {
    x <- arr[, , v, drop = TRUE]; x <- matrix(x, ncol = n_chains)
    q <- stats::quantile(x, c(a, 1 - a), names = FALSE)
    data.frame(param = dimnames(arr)[[3]][v], estimate = mean(x), sd = stats::sd(as.vector(x)),
               lower = q[1], upper = q[2],
               rhat = if (n_chains > 1) posterior::rhat(x) else NA_real_,
               ess_bulk = posterior::ess_bulk(x))
  })
  do.call(rbind, rows)
}

# max/min pairwise distance (km) for the default phi prior
# Default bounds on the effective range: the 5% quantile and the maximum of
# the inter-site distances. (spOccupancy's own default uses the MINIMUM
# distance; with irregular sites that puts ~95% of the uniform-on-phi prior
# mass on ranges shorter than the typical site spacing, where a spatial
# effect is indistinguishable from independent site noise.)
dist_range <- function(coords) {
  d <- as.numeric(stats::dist(coords)); d <- d[d > 0]
  c(min = unname(stats::quantile(d, 0.05)), max = max(d))
}

sp_priors <- function(data, settings, spatial) {
  pr <- list(beta.normal = list(mean = 0, var = 2.72))
  if (data$type != "integrated") pr$alpha.normal <- list(mean = 0, var = 2.72)  # integrated: package default N(0, 2.72) per source
  if (spatial) {
    dr <- dist_range(data$coords)
    rk <- if (!is.null(settings$range_prior_km)) unlist(settings$range_prior_km) else c(dr[["min"]], dr[["max"]])
    pr$phi.unif <- c(3 / rk[2], 3 / rk[1])
    pr$sigma.sq.ig <- unlist(d_or(settings$sigma_sq_ig, c(2, 1)))
  }
  pr
}

ms_priors <- function(data, settings, spatial) {
  pr <- list(beta.comm.normal = list(mean = 0, var = 2.72), alpha.comm.normal = list(mean = 0, var = 2.72),
             tau.sq.beta.ig = list(a = 0.1, b = 0.1), tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
  if (spatial) {
    dr <- dist_range(data$coords)
    rk <- if (!is.null(settings$range_prior_km)) unlist(settings$range_prior_km) else c(dr[["min"]], dr[["max"]])
    q <- d_or(settings$n_factors, 1)
    pr$phi.unif <- list(a = rep(3 / rk[2], q), b = rep(3 / rk[1], q))
  }
  pr
}

covs_df <- function(occ_X, J) {
  if (ncol(occ_X)) occ_X else data.frame(.dummy = rep(0, J))
}

run_spoccupancy <- function(data, settings, model) {
  eco_require("spOccupancy")
  st <- settings
  bl <- d_or(st$batch_length, 25)
  n_iter <- d_or(st$n_iter, 10000)
  n_batch <- ceiling(n_iter / bl)
  n_burn <- d_or(st$n_burn, round(n_iter / 2))
  n_thin <- d_or(st$n_thin, 1)
  n_ch <- d_or(st$n_chains, 3)
  nn <- d_or(st$n_neighbors, 15)
  thr <- d_or(st$n_omp_threads, 1)
  J <- nrow(data$coords)
  psi_f <- stats::as.formula(data$psi_formula)
  set.seed(d_or(st$seed, 1))
  co <- unname(data$coords)
  if (data$type == "single") {
    d <- list(y = unname(data$y), occ.covs = covs_df(data$occ_X, J), det.covs = data$det_X, coords = co)
    if (!length(d$det.covs)) d$det.covs <- list(.dummy = matrix(0, nrow(data$y), ncol(data$y)))
    p_f <- stats::as.formula(data$p_formula)
    if (model == "spPGOcc")
      spOccupancy::spPGOcc(psi_f, p_f, data = d, priors = sp_priors(data, st, TRUE), cov.model = "exponential",
                           NNGP = TRUE, n.neighbors = nn, n.batch = n_batch, batch.length = bl,
                           n.burn = n_burn, n.thin = n_thin, n.chains = n_ch, n.omp.threads = thr,
                           tuning = list(phi = 1), verbose = FALSE)
    else spOccupancy::PGOcc(psi_f, p_f, data = d, priors = sp_priors(data, st, FALSE), n.samples = n_batch * bl,
                            n.burn = n_burn, n.thin = n_thin, n.chains = n_ch, n.omp.threads = thr, verbose = FALSE)
  } else if (data$type == "multi") {
    d <- list(y = unname(data$y), occ.covs = covs_df(data$occ_X, J), det.covs = data$det_X, coords = co)
    dimnames(d$y)[[1]] <- data$species
    if (!length(d$det.covs)) d$det.covs <- list(.dummy = matrix(0, dim(data$y)[2], dim(data$y)[3]))
    p_f <- stats::as.formula(data$p_formula)
    if (model == "sfMsPGOcc")
      spOccupancy::sfMsPGOcc(psi_f, p_f, data = d, priors = ms_priors(data, st, TRUE), cov.model = "exponential",
                             NNGP = TRUE, n.neighbors = nn, n.factors = d_or(st$n_factors, 1),
                             n.batch = n_batch, batch.length = bl, n.burn = n_burn, n.thin = n_thin,
                             n.chains = n_ch, n.omp.threads = thr, tuning = list(phi = 1), verbose = FALSE)
    else spOccupancy::msPGOcc(psi_f, p_f, data = d[c("y", "occ.covs", "det.covs")], priors = ms_priors(data, st, FALSE),
                              n.samples = n_batch * bl, n.burn = n_burn, n.thin = n_thin, n.chains = n_ch,
                              n.omp.threads = thr, verbose = FALSE)
  } else {
    src <- data$sources
    dcov <- lapply(src, function(s) { x <- data$det_X[[s]]
      if (!length(x)) list(.dummy = matrix(0, nrow(data$y[[s]]), ncol(data$y[[s]]))) else x })
    d <- list(y = unname(lapply(src, function(s) unname(data$y[[s]]))), occ.covs = covs_df(data$occ_X, J),
              det.covs = unname(dcov), sites = unname(lapply(src, function(s) data$sites[[s]])))
    p_f <- unname(lapply(src, function(s) stats::as.formula(data$p_formula[[s]])))
    if (model == "spIntPGOcc") {
      d$coords <- co
      spOccupancy::spIntPGOcc(psi_f, p_f, data = d, priors = sp_priors(data, st, TRUE), cov.model = "exponential",
                              NNGP = TRUE, n.neighbors = nn, n.batch = n_batch, batch.length = bl,
                              n.burn = n_burn, n.thin = n_thin, n.chains = n_ch, n.omp.threads = thr,
                              tuning = list(phi = 1), verbose = FALSE)
    } else spOccupancy::intPGOcc(psi_f, p_f, data = d, priors = sp_priors(data, st, FALSE), n.samples = n_batch * bl,
                                 n.burn = n_burn, n.thin = n_thin, n.chains = n_ch, n.omp.threads = thr,
                                 verbose = FALSE)
  }
}

# name the columns of a coefficient-draw matrix: intercept -> probability
coef_block <- function(M, prob_name, prefix, suffix = "") {
  M <- as.matrix(M); cn <- colnames(M)
  out <- list()
  for (k in seq_along(cn)) {
    if (cn[k] == "(Intercept)") out[[paste0(prob_name, suffix)]] <- stats::plogis(M[, k])
    else if (cn[k] != ".dummy") out[[paste0(prefix, clean_name(cn[k]), suffix)]] <- M[, k]
  }
  out
}

estimand_draws <- function(out, data, model) {
  L <- list()
  if (data$type == "single") {
    L <- c(L, coef_block(out$beta.samples, "occupancy_prob", "psi_beta_"),
           coef_block(out$alpha.samples, "detection_prob", "p_beta_"))
    L$mean_occupancy <- rowMeans(out$psi.samples)
    L$prop_sites_occupied <- rowMeans(out$z.samples)
    if (model == "spPGOcc") {
      th <- as.matrix(out$theta.samples)
      L$spatial_variance <- th[, "sigma.sq"]
      L$spatial_decay_per_km <- th[, "phi"]
      L$effective_range_km <- 3 / th[, "phi"]
    }
  } else if (data$type == "multi") {
    bc <- as.matrix(out$beta.comm.samples); ac <- as.matrix(out$alpha.comm.samples)
    tb <- sqrt(as.matrix(out$tau.sq.beta.samples)); ta <- sqrt(as.matrix(out$tau.sq.alpha.samples))
    L <- c(L, coef_block(bc, "comm_occupancy_prob", "comm_psi_beta_"))
    colnames(tb) <- ifelse(colnames(tb) == "(Intercept)", "occupancy_logit", paste0("psi_beta_", clean_name(colnames(tb))))
    for (k in seq_len(ncol(tb))) if (!grepl("dummy", colnames(tb)[k])) L[[paste0("comm_sd_", colnames(tb)[k])]] <- tb[, k]
    L <- c(L, coef_block(ac, "comm_detection_prob", "comm_p_beta_"))
    colnames(ta) <- ifelse(colnames(ta) == "(Intercept)", "detection_logit", paste0("p_beta_", clean_name(colnames(ta))))
    for (k in seq_len(ncol(ta))) if (!grepl("dummy", colnames(ta)[k])) L[[paste0("comm_sd_", colnames(ta)[k])]] <- ta[, k]
    rich <- apply(out$z.samples, c(1, 3), sum)       # draws x sites
    L$mean_site_richness <- rowMeans(rich)
    if (model == "sfMsPGOcc") {
      th <- as.matrix(out$theta.samples)
      for (r in seq_len(ncol(th))) L[[paste0("effective_range_km_f", r)]] <- 3 / th[, r]
    }
    # species-level (not in truth; reported in the table)
    bs <- as.matrix(out$beta.samples); as_ <- as.matrix(out$alpha.samples)
    for (sp in data$species) {
      cols <- grep(paste0("-", sp, "$"), colnames(bs)); m <- bs[, cols, drop = FALSE]
      colnames(m) <- sub(paste0("-", sp, "$"), "", colnames(m))
      L <- c(L, coef_block(m, "occupancy_prob", "psi_beta_", paste0("_", sp)))
      cols <- grep(paste0("-", sp, "$"), colnames(as_)); m <- as_[, cols, drop = FALSE]
      colnames(m) <- sub(paste0("-", sp, "$"), "", colnames(m))
      L <- c(L, coef_block(m, "detection_prob", "p_beta_", paste0("_", sp)))
    }
  } else {
    L <- c(L, coef_block(out$beta.samples, "occupancy_prob", "psi_beta_"))
    A <- as.matrix(out$alpha.samples); pos <- 0
    for (s in data$sources) {
      nt <- 1 + length(data$p_terms[[s]])
      m <- A[, pos + seq_len(nt), drop = FALSE]; pos <- pos + nt
      colnames(m) <- c("(Intercept)", data$p_terms[[s]])
      L <- c(L, coef_block(m, "detection_prob", "p_beta_", paste0("_", s)))
    }
    L$mean_occupancy <- rowMeans(out$psi.samples)
    L$prop_sites_occupied <- rowMeans(out$z.samples)
    if (model == "spIntPGOcc") {
      th <- as.matrix(out$theta.samples)
      L$spatial_variance <- th[, "sigma.sq"]; L$effective_range_km <- 3 / th[, "phi"]
    }
  }
  do.call(cbind, L)
}

# quantities that are weakly identified and NOT used for the convergence flag
# when settings$converge_on = "core" (see VALIDATION.md)
spatial_hyper <- function(p) grepl("^(spatial_variance|spatial_decay|effective_range)", p)
species_level <- function(p, species) if (length(species)) grepl(paste0("_(", paste(species, collapse = "|"), ")$"), p) else rep(FALSE, length(p))

site_summaries <- function(out, data, model) {
  qs <- function(M, f = identity) { M <- f(M)
    data.frame(mean = colMeans(M), sd = apply(M, 2, stats::sd),
               lower = apply(M, 2, stats::quantile, 0.025), upper = apply(M, 2, stats::quantile, 0.975)) }
  st <- data.frame(site = data$site_ids, x = data$coords[, 1], y = data$coords[, 2])
  if (data$type == "multi") {
    rich <- apply(out$z.samples, c(1, 3), sum)
    r <- qs(rich); st$richness <- r$mean; st$richness_sd <- r$sd; st$richness_lower <- r$lower; st$richness_upper <- r$upper
    st$naive_richness <- apply(data$y, 2, function(m) sum(apply(m, 1, function(v) any(v %in% 1))))
    if (model == "sfMsPGOcc") {
      # w.samples: draws x factors x sites; species-level effect = lambda %*% w
      st$factor1 <- colMeans(out$w.samples[, 1, , drop = TRUE])
    }
  } else {
    p <- qs(out$psi.samples); st$psi <- p$mean; st$psi_sd <- p$sd; st$psi_lower <- p$lower; st$psi_upper <- p$upper
    st$z <- colMeans(out$z.samples)
    if (!is.null(out$w.samples)) { w <- qs(out$w.samples); st$w <- w$mean; st$w_sd <- w$sd }
    st$detected <- if (data$type == "single") apply(data$y, 1, function(v) any(v %in% 1))
                   else vapply(seq_len(nrow(data$coords)), function(j) any(vapply(data$sources, function(s) {
                     i <- which(data$sites[[s]] == j); length(i) && any(data$y[[s]][i, ] %in% 1) }, logical(1))), logical(1))
  }
  st
}

gof_check <- function(out, data) {
  g <- tryCatch(spOccupancy::ppcOcc(out, fit.stat = "freeman-tukey", group = 1),
                error = function(e) NULL)
  if (is.null(g)) return(NULL)
  if (is.list(g$fit.y)) {   # integrated: one per source
    return(lapply(seq_along(g$fit.y), function(r) list(obs = g$fit.y[[r]], rep = g$fit.y.rep[[r]],
                                                       p = mean(g$fit.y.rep[[r]] > g$fit.y[[r]]),
                                                       label = data$sources[r])))
  }
  if (is.matrix(g$fit.y)) {  # multi: draws x species
    return(lapply(seq_len(ncol(g$fit.y)), function(i) list(obs = g$fit.y[, i], rep = g$fit.y.rep[, i],
                                                            p = mean(g$fit.y.rep[, i] > g$fit.y[, i]),
                                                            label = data$species[i])))
  }
  list(list(obs = g$fit.y, rep = g$fit.y.rep, p = mean(g$fit.y.rep > g$fit.y), label = "all sites"))
}

predict_grid <- function(out, data, model, settings) {
  pg <- data$pred
  if (is.null(pg) || isFALSE(settings$predict)) return(NULL)
  X0 <- cbind(1, as.matrix(pg[, data$psi_terms, drop = FALSE]))
  colnames(X0)[1] <- "(Intercept)"
  co0 <- as.matrix(pg[, c("x", "y")])
  args <- list(object = out, X.0 = X0, verbose = FALSE)
  if (is_spatial_model(model)) { args$coords.0 <- co0; args$n.omp.threads <- d_or(settings$n_omp_threads, 1) }
  pr <- tryCatch(do.call(stats::predict, args), error = function(e) { message("prediction failed: ", conditionMessage(e)); NULL })
  if (is.null(pr)) return(NULL)
  g <- data.frame(x = pg$x, y = pg$y)
  if (data$type == "multi") {
    rich <- apply(pr$z.0.samples, c(1, 3), sum)
    g$richness <- colMeans(rich); g$richness_sd <- apply(rich, 2, stats::sd)
  } else {
    g$psi <- colMeans(pr$psi.0.samples); g$psi_sd <- apply(pr$psi.0.samples, 2, stats::sd)
    if (!is.null(pr$w.0.samples)) g$w <- colMeans(pr$w.0.samples)
  }
  g
}

fit_occu_spatial <- function(data, settings) {
  model <- resolve_model(settings$model, data$type)
  t0 <- proc.time()[["elapsed"]]
  out <- run_spoccupancy(data, settings, model)
  D <- estimand_draws(out, data, model)
  n_ch <- out$n.chains
  est <- summarise_estimands(D, n_ch)
  conv_on <- d_or(settings$converge_on, "all")
  core <- if (conv_on == "core") !spatial_hyper(est$param) & !species_level(est$param, data$species) else rep(TRUE, nrow(est))
  rmax <- max(est$rhat[core], na.rm = TRUE)
  notes <- sprintf("spOccupancy::%s; %d chains x %d post-burn-in draws kept (thin %d); %.0f s.",
                   model, n_ch, out$n.post, out$n.thin, proc.time()[["elapsed"]] - t0)
  sh <- spatial_hyper(est$param)
  if (any(sh)) notes <- c(notes, sprintf("Spatial range/variance: max R-hat %.3f, min bulk ESS %.0f. These are weakly identified; read them as descriptive (see Known limits).",
                                         max(est$rhat[sh], na.rm = TRUE), min(est$ess_bulk[sh], na.rm = TRUE)))
  if (conv_on == "core") notes <- c(notes, "Convergence flag uses all estimands except the spatial range/variance and species-level coefficients.")
  gof <- if (isTRUE(settings$gof)) gof_check(out, data) else NULL
  if (!is.null(gof)) {
    ps <- vapply(gof, `[[`, numeric(1), "p")
    notes <- c(notes, sprintf("Posterior predictive check (Freeman-Tukey, sites): Bayesian p = %s%s.",
                              paste(sprintf("%.2f", ps[seq_len(min(8, length(ps)))]), collapse = ", "),
                              if (length(ps) > 8) sprintf(" ... (%d species; %d outside 0.1-0.9)", length(ps), sum(ps < 0.1 | ps > 0.9)) else ""))
  }
  waic <- tryCatch(spOccupancy::waicOcc(out), error = function(e) NULL)
  sites <- site_summaries(out, data, model)
  pred <- predict_grid(out, data, model, settings)
  ess_min <- min(est$ess_bulk[core], na.rm = TRUE)
  if (conv_on == "all" && is.finite(rmax) && (rmax >= 1.01 || ess_min < 400))
    notes <- c(notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): rerun with more iterations (settings: n_iter, n_burn, n_thin).",
                              rmax, ess_min))
  if (length(data$notes)) notes <- c(paste("Data:", data$notes), notes)
  keep <- if (isTRUE(settings$keep_fit)) out else NULL
  list(estimates = est,
       diagnostics = list(converged = is.finite(rmax) && rmax < 1.05, rhat_max = rmax,
                          ess_min = ess_min, notes = notes),
       object = list(model = model, type = data$type, draws = D, n_chains = n_ch,
                     beta = as.matrix(if (data$type == "multi") out$beta.comm.samples else out$beta.samples),
                     beta_species = if (data$type == "multi") as.matrix(out$beta.samples),
                     alpha = as.matrix(if (data$type == "multi") out$alpha.comm.samples else out$alpha.samples),
                     sites = sites, pred = pred, gof = gof, waic = waic,
                     lambda = if (!is.null(out$lambda.samples)) colMeans(out$lambda.samples),
                     spfit = keep, settings = settings))
}
