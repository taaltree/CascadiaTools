# nmixture: from engine output to the contract estimates table, plus derived
# abundance, goodness of fit, K sensitivity and mixture comparison.

# Public (biological-scale) quantities from a matrix of theta draws.
# Returns list(values = matrix, link = named vector of "log"/"logit"/"identity").
nm_public <- function(TH, des, mixture, dyn = "constant") {
  TH <- as.matrix(TH)
  if (is.null(dim(TH)) || nrow(TH) == 0) stop("no draws")
  col <- function(n) TH[, n]
  zi <- if (mixture == "ZIP") stats::plogis(col("logit_zi")) else 0
  out <- list(); link <- character()
  add <- function(n, v, l) { out[[n]] <<- v; link[[n]] <<- l }
  ab_name <- if (des$model == "static") "abundance_at_average_site" else "initial_abundance_per_site"
  add(ab_name, exp(col("lam[1]")) * (1 - zi), "log")
  add("detection_prob", stats::plogis(col("p[1]")), "logit")
  if (des$model == "dynamic") {
    if (dyn != "notrend") add("recruitment_per_site", exp(col("gam[1]")), "log")
    add("apparent_survival", stats::plogis(col("om[1]")), "logit")
  }
  if (mixture == "NB") add("nb_size", exp(col("log_size")), "log")
  if (mixture == "ZIP") add("zero_inflation_prob", zi, "logit")
  slopes <- function(pre, cn, lab, rr = NULL) {
    if (length(cn) < 2) return()
    for (k in 2:length(cn)) {
      v <- col(sprintf("%s[%d]", pre, k)); nm <- nm_clean(cn[k])
      add(paste0(lab, "_beta_", nm), v, "identity")
      if (!is.null(rr)) add(paste0(rr, "_", nm), exp(v), "log")
    }
  }
  slopes("lam", des$colnames$lam, "lambda", "rate_ratio")
  slopes("p", des$colnames$p, "p", "odds_ratio_detection")
  if (des$model == "dynamic") {
    if (dyn != "notrend") slopes("gam", des$colnames$gam, "gamma", "rate_ratio_recruitment")
    slopes("om", des$colnames$om, "omega", "odds_ratio_survival")
  }
  list(values = do.call(cbind, out), link = unlist(link))
}

nm_link <- function(l) switch(l, log = list(f = log, inv = exp), logit = list(f = stats::qlogis, inv = stats::plogis),
                              list(f = identity, inv = identity))

# Wald / delta-method rows for the public quantities of a likelihood fit.
nm_wald_public <- function(eng, des, mixture, dyn) {
  th <- eng$theta; V <- eng$vcov
  p0 <- nm_public(matrix(th, 1, dimnames = list(NULL, names(th))), des, mixture, dyn)
  est <- p0$values[1, ]
  rows <- lapply(names(est), function(n) {
    lk <- nm_link(p0$link[[n]])
    gfun <- function(t) { names(t) <- names(th); lk$f(nm_public(matrix(t, 1, dimnames = list(NULL, names(th))), des, mixture, dyn)$values[1, n]) }
    g0 <- gfun(th)
    grad <- vapply(seq_along(th), function(k) { h <- 1e-5 * max(1, abs(th[k])); e <- th; e[k] <- e[k] + h
      (gfun(e) - g0) / h }, 1)
    se_g <- sqrt(max(drop(t(grad) %*% V %*% grad), 0))
    if (!is.finite(se_g)) se_g <- NA_real_
    lo <- lk$inv(g0 - 1.96 * se_g); hi <- lk$inv(g0 + 1.96 * se_g)
    dinv <- (lk$inv(g0 + 1e-6) - lk$inv(g0 - 1e-6)) / 2e-6
    data.frame(param = n, estimate = est[[n]], sd = abs(dinv) * se_g, lower = min(lo, hi), upper = max(lo, hi))
  })
  do.call(rbind, rows)
}

# Derived abundance quantities from N draws (nd x M x T)
nm_N_public <- function(Nd, des) {
  if (des$model == "static") {
    tot <- apply(Nd[, , 1, drop = FALSE], 1, sum)
    cbind(total_abundance = tot, mean_abundance_per_site = tot / des$M)
  } else {
    tot <- apply(Nd, c(1, 3), sum)
    colnames(tot) <- paste0("total_abundance_t", seq_len(des$T))
    cbind(tot, abundance_ratio_last_first = tot[, des$T] / pmax(tot[, 1], 1))
  }
}

nm_site_summary <- function(Nd, des) {
  rows <- list()
  for (t in seq_len(des$T)) {
    m <- Nd[, , t, drop = FALSE]; dim(m) <- dim(Nd)[1:2]
    yy <- matrix(des$y[, , t], des$M, des$J)
    rows[[t]] <- data.frame(site = des$site_ids, season = des$season_labels[t],
                            max_count = suppressWarnings(apply(yy, 1, max, na.rm = TRUE)),
                            estimate = colMeans(m), lower = apply(m, 2, stats::quantile, 0.025),
                            upper = apply(m, 2, stats::quantile, 0.975))
  }
  out <- do.call(rbind, rows); out$max_count[!is.finite(out$max_count)] <- NA
  out
}

nm_choose_K_post <- function(eng, des, mixture, dyn, tol) {
  K <- eng$K
  u <- nm_unpack(eng$theta, des, mixture, dyn)
  while (nm_tail_prob(u, K, mixture, des) > tol && K < 2000) K <- ceiling(K * 1.5)
  K
}

nm_fit <- function(data, settings) {
  settings <- merge_list(nm_settings, settings)
  set.seed(settings$seed)
  mixture <- match.arg(settings$mixture, c("P", "NB", "ZIP"))
  settings$mixture <- mixture
  dyn <- match.arg(settings$dynamics, c("constant", "notrend"))
  des <- nm_prepare(data, settings)
  if (des$model == "dynamic" && des$T < 2) stop("dynamic model needs at least 2 seasons")
  engine <- match.arg(settings$engine, c("unmarked", "nimble", "rlik"))
  eng <- switch(engine, unmarked = nm_engine_unmarked(des, settings),
                rlik = nm_engine_rlik(des, settings),
                nimble = nm_engine_nimble(des, settings))
  if (engine == "nimble") eng$K <- nm_choose_K_post(eng, des, mixture, dyn, settings$K_tail_tol)
  notes <- eng$notes
  # ---- derived abundance with parameter uncertainty propagated
  if (engine == "nimble") {
    pub_chains <- lapply(eng$chains, function(TH) {
      Nd <- nm_N_draws(TH, des, eng$K, mixture, dyn)
      cbind(nm_public(TH, des, mixture, dyn)$values, nm_N_public(Nd, des))
    })
    s <- summarise_draws_eco(pub_chains)
    est <- s$estimates
    Nd_all <- nm_N_draws(eng$draws[seq(1, nrow(eng$draws), length.out = min(nrow(eng$draws), 1000)), , drop = FALSE],
                         des, eng$K, mixture, dyn)
    dg <- list(converged = s$diagnostics$converged, rhat_max = s$diagnostics$rhat_max, ess_min = s$diagnostics$ess_min)
    if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
      notes <- c(notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): run longer (settings niter/nburnin).",
                                dg$rhat_max, dg$ess_min))
  } else {
    est <- nm_wald_public(eng, des, mixture, dyn)
    Nd_all <- nm_N_draws(eng$draws, des, eng$K, mixture, dyn)
    est <- rbind(est, derived_rows(nm_N_public(Nd_all, des)))
    dg <- list(converged = isTRUE(eng$converged), rhat_max = NA_real_, ess_min = NA_real_)
  }
  rownames(est) <- NULL
  site_N <- nm_site_summary(Nd_all, des)
  obj <- list(engine = engine, des = des, theta = eng$theta, vcov = eng$vcov, draws = eng$draws,
              chains = eng$chains, K = eng$K, loglik = eng$loglik, aic = eng$aic, fm = eng$fm,
              settings = settings, site_N = site_N, N_draws_total = nm_N_public(Nd_all, des))
  # ---- K sensitivity
  if (isTRUE(settings$k_sensitivity)) {
    obj$k_sens <- tryCatch(nm_k_sensitivity(eng, des, settings), error = function(e) { notes <<- c(notes, paste("K sensitivity failed:", conditionMessage(e))); NULL })
    if (!is.null(obj$k_sens)) {
      ks <- obj$k_sens; rel <- max(abs(ks$total_abundance / ks$total_abundance[1] - 1))
      notes <- c(notes, sprintf("K sensitivity: K = %s; total abundance changes by at most %.2f%% when K is increased%s.",
                                paste(ks$K, collapse = ", "), 100 * rel,
                                if (rel > 0.01) " - MORE THAN 1%: abundance is sensitive to the upper bound and may be poorly identified" else ""))
    }
  }
  # ---- mixture comparison (closed model, likelihood engines)
  if (isTRUE(settings$compare_mixtures) && engine != "nimble" && des$model == "static") {
    obj$mixtures <- tryCatch(nm_compare_mixtures(des, settings), error = function(e) NULL)
    if (!is.null(obj$mixtures)) {
      mx <- obj$mixtures
      best <- mx$mixture[which.min(mx$AIC)]
      notes <- c(notes, sprintf("Abundance distribution comparison (AIC): %s. Lowest AIC: %s%s.",
                                paste(sprintf("%s %.1f", mx$mixture, mx$AIC), collapse = ", "), best,
                                if (best != mixture) sprintf(" - the fitted model used %s; consider refitting with mixture: %s and compare abundance estimates", mixture, best) else ""))
    }
  }
  # ---- goodness of fit
  if (isTRUE(settings$gof_nboot > 0)) {
    obj$gof <- tryCatch(nm_gof(eng, des, settings), error = function(e) { notes <<- c(notes, paste("GOF failed:", conditionMessage(e))); NULL })
    if (!is.null(obj$gof)) {
      g <- obj$gof
      notes <- c(notes, sprintf("Goodness of fit (%s, %d replicates): chi-square %.1f, c-hat %.2f, p = %.3f; Freeman-Tukey p = %.3f.%s",
                                g$method, g$nboot, g$chisq_obs, g$c_hat, g$p_chisq, g$p_ft,
                                if (min(g$p_chisq, g$p_ft) < 0.05) " LACK OF FIT: abundance estimates may be biased (overdispersion, unmodelled heterogeneity or closure violation)." else ""))
    }
  }
  dg$notes <- notes
  list(estimates = est, diagnostics = dg, object = obj)
}

nm_k_sensitivity <- function(eng, des, settings) {
  mixture <- settings$mixture; dyn <- settings$dynamics
  Ks <- unique(c(eng$K, ceiling(eng$K * settings$k_mult)))
  plug <- function(th, K) {
    u <- nm_unpack(th, des, mixture, dyn)
    if (des$model == "static") {
      pp <- nm_static_post(u, des, K, mixture)$post
      sum(pp %*% (0:K))
    } else NA_real_
  }
  rows <- lapply(Ks, function(K) {
    if (eng$engine == "nimble" || K == eng$K) th <- eng$theta
    else {
      e2 <- if (eng$engine == "unmarked") nm_engine_unmarked(des, merge_list(settings, list(n_draws = 2)), K = K)
            else nm_engine_rlik(des, merge_list(settings, list(n_draws = 2)), K = K)
      th <- e2$theta
    }
    pub <- nm_public(matrix(th, 1, dimnames = list(NULL, names(th))), des, mixture, dyn)$values[1, 1]
    tot <- if (des$model == "static") plug(th, K) else {
      names(th) <- names(eng$theta)
      sum(nm_dyn_ffbs(nm_unpack(th, des, mixture, dyn), des, K, mixture)$FT %*% (0:K))
    }
    data.frame(K = K, abundance = pub, total_abundance = tot)
  })
  do.call(rbind, rows)
}

nm_compare_mixtures <- function(des, settings) {
  rows <- lapply(c("P", "NB", "ZIP"), function(mx) {
    st <- merge_list(settings, list(mixture = mx, n_draws = 2))
    e <- tryCatch(if (settings$engine == "rlik") nm_engine_rlik(des, st) else nm_engine_unmarked(des, st), error = function(e) NULL)
    if (is.null(e)) return(data.frame(mixture = mx, AIC = NA, logLik = NA, abundance_at_average_site = NA, converged = FALSE))
    ab <- nm_public(matrix(e$theta, 1, dimnames = list(NULL, names(e$theta))), des, mx)$values[1, 1]
    data.frame(mixture = mx, AIC = e$aic, logLik = e$loglik, abundance_at_average_site = ab, converged = e$converged)
  })
  out <- do.call(rbind, rows)
  out$delta_AIC <- out$AIC - min(out$AIC, na.rm = TRUE)
  out
}
