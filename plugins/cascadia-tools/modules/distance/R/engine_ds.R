# Conventional / multiple-covariate distance sampling with Distance::ds
# (Miller et al. 2019). Candidate detection functions: each key in
# settings$keys; with settings$adjustment ("cos" or "herm") ds() runs
# sequential forward AIC selection of adjustment terms (0, 1, 2, ...), so the
# no-adjustment model is always a candidate. With detection covariates no
# adjustments are fitted (Marques et al. 2007). The model with the lowest AIC
# is used. Density and abundance come from dht with the Fewster et al. (2009)
# R2 (lines) / P2 (points) encounter-rate variance, the detection-function
# variance (delta method) and the cluster-size variance, combined with
# Satterthwaite df and log-normal intervals (Buckland et al. 2001).

ds_engine_distance <- function(data, settings) {
  eco_require("Distance")
  w <- data$truncation_m
  o <- data$obs_w
  tr <- data$transects
  nV <- data$n_visits
  flat <- data.frame(object = seq_len(nrow(o)), Sample.Label = o$transect,
                     Region.Label = tr$stratum[match(o$transect, tr$transect)],
                     stringsAsFactors = FALSE)
  if (data$distance_type == "binned") { flat$distbegin <- o$distbegin; flat$distend <- o$distend; flat$distance <- (o$distbegin + o$distend) / 2 }
  else flat$distance <- o$distance
  if (isTRUE(data$clusters)) flat$size <- o$size
  # covariates get a prefix: Distance/mrds reserve names such as 'observer', 'detected', 'size'
  for (v in data$dt_terms) flat[[paste0("cv_", v)]] <- o[[v]]
  has_area <- !is.null(data$strata)
  rt <- if (has_area) data.frame(Region.Label = data$strata$stratum, Area = data$strata$area_km2)
        else data.frame(Region.Label = unique(tr$stratum), Area = 0)
  st <- data.frame(Sample.Label = tr$transect, Region.Label = tr$stratum, Effort = tr$effort * nV)
  ot <- flat[, c("object", "Region.Label", "Sample.Label")]
  dd <- flat[, setdiff(names(flat), c("Region.Label", "Sample.Label"))]
  ftxt <- ds_rhs(data$formulas$detection)
  for (v in data$dt_terms) ftxt <- gsub(paste0("\\b", v, "\\b"), paste0("cv_", v), ftxt, perl = TRUE)
  form <- stats::as.formula(ftxt)
  cov_model <- length(data$dt_terms) > 0
  cands <- list()
  for (k in settings$keys) {
    adj <- if (cov_model || identical(settings$adjustment, "none")) NULL else settings$adjustment
    args <- list(data = dd, truncation = w, transect = data$survey, formula = form, key = k,
                 adjustment = adj, region_table = rt, sample_table = st, obs_table = ot,
                 convert_units = 0.001, quiet = TRUE, max_adjustments = settings$max_adjustments)
    if (data$distance_type == "binned") args$cutpoints <- data$cutpoints[data$cutpoints <= w + 1e-9]
    m <- tryCatch(suppressWarnings(suppressMessages(do.call(Distance::ds, args))), error = function(e) e)
    if (inherits(m, "error") || is.null(m$ddf) || !is.finite(m$ddf$criterion)) next
    cands[[k]] <- m
  }
  if (!length(cands)) stop("Distance::ds failed for every candidate detection function")
  # Stability screen: a candidate whose average detection probability has CV > 50% (e.g. a
  # hazard-rate fitted as a spike at zero to binned data, P_a -> 0, density -> infinity) is not
  # identified by the data. It is dropped when at least one stable candidate remains; if none
  # remains the fit is flagged as not converged.
  pa_cv <- vapply(cands, function(m) { sm <- tryCatch(summary(m$ddf), error = function(e) NULL)
    if (is.null(sm) || !is.finite(sm$average.p.se)) Inf else sm$average.p.se / sm$average.p }, 1)
  unstable <- names(cands)[pa_cv > settings$max_pa_cv %||% 0.5]
  dropped <- character()
  if (length(unstable) && length(unstable) < length(cands)) { dropped <- unstable; cands <- cands[setdiff(names(cands), unstable)] }
  all_unstable <- length(unstable) == length(pa_cv)
  aic <- vapply(cands, function(m) m$ddf$criterion, 1)
  best <- names(which.min(aic))
  m <- cands[[best]]
  nadj <- if (is.null(m$ddf$ds$aux$ddfobj$adjustment)) 0L else length(m$ddf$ds$aux$ddfobj$adjustment$order)
  sel <- data.frame(model = vapply(names(cands), function(k) {
      a <- cands[[k]]$ddf$ds$aux$ddfobj$adjustment
      paste0(c(hn = "half-normal", hr = "hazard-rate", unif = "uniform")[[k]],
             if (!is.null(a)) paste0(" + ", length(a$order), " ", a$series, " adj.") else "",
             if (cov_model) paste0(" (", ds_rhs(data$formulas$detection), ")") else "")
    }, ""), AIC = aic, dAIC = aic - min(aic), row.names = NULL)
  sel <- sel[order(sel$AIC), ]

  ests <- lapply(names(cands), function(k) ds_model_rows(cands[[k]], k, data, has_area, w))
  est <- if (isTRUE(settings$model_average) && length(cands) > 1) ds_model_average(ests, aic, data, w) else ests[[which.min(aic)]]

  if (isTRUE(settings$.boot)) return(list(estimates = est))
  gof <- tryCatch({
    g <- suppressWarnings(Distance::gof_ds(m, plot = FALSE, chisq = data$distance_type == "binned"))
    if (data$distance_type == "binned" || is.null(g$dsgof$CvM)) list(type = "chisq", stat = g$chisquare$chi1$chisq, p = g$chisquare$chi1$p, df = g$chisquare$chi1$df)
    else list(type = "cvm", stat = g$dsgof$CvM$W, p = g$dsgof$CvM$p, edf = g$dsgof$edf, cdf = g$dsgof$cdf)
  }, error = function(e) NULL)

  conv <- all(is.finite(est$sd[est$param %in% c("density_per_km2", "detection_prob")])) &&
          isTRUE(m$ddf$ds$converge == 0 || is.null(m$ddf$ds$converge))
  notes <- character()
  if (!conv) notes <- c(notes, "Distance::ds did not converge or gave non-finite standard errors.")
  if (length(dropped)) notes <- c(notes, sprintf("Candidate detection function(s) %s dropped: detection probability not identified (CV > 50%%, e.g. a spike at zero distance).",
                                                 paste(vapply(dropped, dist_keyname, ""), collapse = ", ")))
  if (all_unstable) { conv <- FALSE; notes <- c(notes, "Every candidate detection function has CV(detection probability) > 50%: the distances do not determine the detection function; density is not reliable.") }
  if (sum(o$visit >= 0) < 60) notes <- c(notes, sprintf("Only %d detections within truncation; Buckland et al. (2001, 2015) recommend at least 60-80 for a reliable detection function.", nrow(o)))
  boot <- NULL
  if (identical(settings$variance, "bootstrap") && (settings$nboot %||% 0) > 0) {
    boot <- ds_bootstrap(data, settings, est)
    if (boot$n_ok >= max(20, 0.8 * boot$n)) {
      for (p in colnames(boot$draws)) {
        i <- which(est$param == p)
        est$sd[i] <- stats::sd(boot$draws[, p])
        est[i, c("lower", "upper")] <- stats::quantile(boot$draws[, p], c(0.025, 0.975), names = FALSE)
      }
    } else {
      conv <- FALSE
      notes <- c(notes, sprintf("Bootstrap failed in %d of %d replicates; analytic intervals kept.", boot$n - boot$n_ok, boot$n))
    }
  }
  list(engine = "Distance::ds", estimates = est, fit = m, boot = boot, key = best, nadj = nadj,
       selection = sel, gof = gof, converged = conv, notes = notes)
}

# Estimates table from one fitted ds model.
ds_model_rows <- function(m, key, data, has_area, w) {
  nadj0 <- if (is.null(m$ddf$ds$aux$ddfobj$adjustment)) 0L else length(m$ddf$ds$aux$ddfobj$adjustment$order)
  sm <- summary(m$ddf)
  ind <- m$dht$individuals
  Dt <- ind$D; Nt <- ind$N
  tot <- if (nrow(Dt) > 1) Dt[Dt$Label == "Total", ] else Dt[1, ]
  rows <- list(ds_row("density_per_km2", tot$Estimate, tot$se, tot$lcl, tot$ucl))
  if (has_area) {
    Ntot <- if (nrow(Nt) > 1) Nt[Nt$Label == "Total", ] else Nt[1, ]
    rows[[length(rows) + 1]] <- ds_row("abundance", Ntot$Estimate, Ntot$se, Ntot$lcl, Ntot$ucl)
  }
  strat <- Dt[Dt$Label != "Total", ]
  if (nrow(strat) > 1) for (i in seq_len(nrow(strat))) {
    rows[[length(rows) + 1]] <- ds_row(paste0("density_per_km2_", ds_clean(strat$Label[i])), strat$Estimate[i], strat$se[i], strat$lcl[i], strat$ucl[i])
    if (has_area) { nn <- Nt[Nt$Label == strat$Label[i], ]
      rows[[length(rows) + 1]] <- ds_row(paste0("abundance_", ds_clean(strat$Label[i])), nn$Estimate, nn$se, nn$lcl, nn$ucl) }
  }
  if (isTRUE(data$clusters)) {
    cl <- m$dht$clusters$D; ct <- if (nrow(cl) > 1) cl[cl$Label == "Total", ] else cl[1, ]
    rows[[length(rows) + 1]] <- ds_row("group_density_per_km2", ct$Estimate, ct$se, ct$lcl, ct$ucl)
    es <- m$dht$Expected.S; et <- if (nrow(es) > 1) es[es$Region == "Total", ] else es[1, ]
    rows[[length(rows) + 1]] <- ds_lognorm_row("mean_group_size", et$Expected.S, et$se.Expected.S)
  }
  pa <- sm$average.p; pa_se <- sm$average.p.se
  rows[[length(rows) + 1]] <- ds_logit_row("detection_prob", pa, pa_se)
  if (data$survey == "line") {
    rows[[length(rows) + 1]] <- ds_row("esw_m", pa * w, pa_se * w, NA, NA)
  } else {
    rows[[length(rows) + 1]] <- ds_row("edr_m", w * sqrt(pa), w * pa_se / (2 * sqrt(pa)), NA, NA)
  }
  r <- rows[[length(rows)]]; pr <- rows[[length(rows) - 1]]
  rows[[length(rows)]]$lower <- if (data$survey == "line") pr$lower * w else w * sqrt(pr$lower)
  rows[[length(rows)]]$upper <- if (data$survey == "line") pr$upper * w else w * sqrt(pr$upper)
  ks <- sm$coeff$key.scale
  if (!is.null(ks) && nadj0 == 0) {
    nm0 <- if (key == "hn") "sigma_m" else "hr_scale_m"
    rows[[length(rows) + 1]] <- ds_row(nm0, exp(ks$estimate[1]), exp(ks$estimate[1]) * ks$se[1],
                                       exp(ks$estimate[1] - 1.96 * ks$se[1]), exp(ks$estimate[1] + 1.96 * ks$se[1]))
  }
  if (!is.null(ks) && nrow(ks) > 1) for (i in 2:nrow(ks))
    rows[[length(rows) + 1]] <- ds_row(paste0("detection_beta_", ds_clean(sub("^cv_", "", rownames(ks)[i]))), ks$estimate[i], ks$se[i],
                                       ks$estimate[i] - 1.96 * ks$se[i], ks$estimate[i] + 1.96 * ks$se[i])
  if (key == "hr" && nadj0 == 0 && !is.null(sm$coeff$key.shape)) {
    sh <- sm$coeff$key.shape
    rows[[length(rows) + 1]] <- ds_row("hr_shape", exp(sh$estimate[1]), exp(sh$estimate[1]) * sh$se[1],
                                       exp(sh$estimate[1] - 1.96 * sh$se[1]), exp(sh$estimate[1] + 1.96 * sh$se[1]))
  }
  isum <- ind$summary; itot <- if (nrow(isum) > 1) isum[isum$Region == "Total", ] else isum[1, ]
  er_name <- if (data$survey == "line") "encounter_rate_per_km" else "encounter_rate_per_point"
  rows[[length(rows) + 1]] <- ds_lognorm_row(er_name, itot$ER, itot$se.ER)
  est <- do.call(rbind, rows)
  est
}

# Model averaging over the candidate key-function families (Buckland, Burnham &
# Augustin 1997): theta_MA = sum w_k theta_k with Akaike weights w_k, and
# se_MA = sum w_k sqrt(se_k^2 + (theta_k - theta_MA)^2), which adds
# model-selection uncertainty. Applied to the quantities every candidate
# estimates (density, abundance, P_a, ESW/EDR, encounter rate, group size);
# key-specific parameters (sigma, shape, covariate effects) come from the
# AIC-best model.
ds_model_average <- function(ests, aic, data, w) {
  wt <- exp(-(aic - min(aic)) / 2); wt <- wt / sum(wt)
  best <- ests[[which.min(aic)]]
  common <- Reduce(intersect, lapply(ests, function(e) e$param))
  common <- intersect(common, c("density_per_km2", "abundance", "group_density_per_km2", "mean_group_size", "detection_prob",
                                grep("^(density_per_km2_|abundance_|encounter_rate_)", best$param, value = TRUE)))
  out <- best
  for (p in common) {
    th <- vapply(ests, function(e) e$estimate[e$param == p], 1)
    se <- vapply(ests, function(e) e$sd[e$param == p], 1)
    ma <- sum(wt * th); sma <- sum(wt * sqrt(se^2 + (th - ma)^2))
    r <- if (p == "detection_prob") ds_logit_row(p, ma, sma) else ds_lognorm_row(p, ma, sma)
    out[out$param == p, c("estimate", "sd", "lower", "upper")] <- r[, c("estimate", "sd", "lower", "upper")]
  }
  pr <- out[out$param == "detection_prob", ]
  ie <- which(out$param %in% c("esw_m", "edr_m"))
  if (length(ie)) {
    if (data$survey == "line") out[ie, c("estimate", "sd", "lower", "upper")] <- c(pr$estimate, pr$sd, pr$lower, pr$upper) * w
    else out[ie, c("estimate", "sd", "lower", "upper")] <- c(w * sqrt(pr$estimate), w * pr$sd / (2 * sqrt(pr$estimate)), w * sqrt(pr$lower), w * sqrt(pr$upper))
  }
  out
}

# Nonparametric bootstrap over transects (Buckland et al. 2001 sec. 3.6.4;
# the approach of Distance::bootdht): resample transects (lines/points) with
# replacement within each stratum, refit the whole detection-function
# selection (keys, adjustment terms, AIC / model averaging) on every
# replicate and recompute density. The bootstrap SD replaces the analytic SE
# and the 2.5/97.5% percentiles give the interval, so the uncertainty from
# choosing the detection function is included. Key-specific parameters
# (sigma, shape, covariate effects) keep their analytic (Wald) intervals.
ds_boot_params <- function(p) grepl("^(density_per_km2|abundance|group_density_per_km2|mean_group_size|detection_prob|esw_m|edr_m|encounter_rate_)", p)

ds_bootstrap <- function(data, settings, est) {
  B <- settings$nboot
  tr <- data$transects
  o <- data$obs_w
  st_idx <- split(seq_len(nrow(tr)), tr$stratum)
  obs_by_tr <- split(seq_len(nrow(o)), factor(o$transect, levels = tr$transect))
  bs <- merge_list(settings, list(.boot = TRUE))
  params <- est$param[ds_boot_params(est$param)]
  out <- matrix(NA_real_, B, length(params), dimnames = list(NULL, params))
  for (b in seq_len(B)) {
    pick <- unlist(lapply(st_idx, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]), use.names = FALSE)
    newid <- sprintf("b%04d", seq_along(pick))
    t2 <- tr[pick, , drop = FALSE]; t2$transect <- newid
    rows <- lapply(seq_along(pick), function(k) { r <- obs_by_tr[[pick[k]]]; if (!length(r)) return(NULL)
      x <- o[r, , drop = FALSE]; x$transect <- newid[k]; x })
    o2 <- do.call(rbind, rows)
    if (is.null(o2) || nrow(o2) < 10) next
    d2 <- data; d2$transects <- t2; d2$obs_w <- o2
    e <- tryCatch(ds_engine_distance(d2, bs)$estimates, error = function(e) NULL)
    if (is.null(e)) next
    out[b, ] <- e$estimate[match(params, e$param)]
  }
  ok <- stats::complete.cases(out)
  list(draws = out[ok, , drop = FALSE], n_ok = sum(ok), n = B)
}
