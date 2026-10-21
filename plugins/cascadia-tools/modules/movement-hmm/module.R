# movement-hmm: hidden Markov models for animal movement tracks
# (Langrock et al. 2012; Michelot et al. 2016; McClintock & Michelot 2018).
#
# Data structure (what simulate()$data and read_data() both return):
#   list(track = data.frame(ID, time [POSIXct], x, y [metres, projected],
#                           <covariate columns>),
#        covariates = character()   transition-probability covariates,
#        tod = TRUE/FALSE           model time of day on transitions,
#        tz = "UTC"                 time zone used for the clock hour,
#        nstates = 2,               default number of states,
#        interval_h = NULL,         nominal fix interval (NULL = median gap),
#        crs = character            CRS of x/y (documentation only))

source(file.path(MODULE_DIR, "hmm_engine.R"), local = TRUE)
source(file.path(MODULE_DIR, "output.R"), local = TRUE)

simulate <- function(params, design, seed) {
  sim <- simulate_tracks(params, design, seed)
  data <- list(track = sim$track, covariates = "forest", tod = TRUE, tz = "UTC",
               nstates = length(params$step_mean), interval_h = design$fix_interval_h %||% 1,
               crs = "local metres")
  # truth must use the same covariate centre as fit(): mean over the regular grid
  g <- regularise_grid(data$track[order(data$track$ID, data$track$time), ], data$interval_h)
  g$forest <- stats::ave(g$forest, g$ID, FUN = fill_locf)
  list(data = data, truth = truth_from_params(params, mean(g$forest, na.rm = TRUE), data$interval_h))
}

fit <- function(data, settings) {
  eco_require("momentuHMM")
  K <- as.integer(settings$nstates %||% data$nstates %||% 2)
  prep <- prep_track(data, settings)
  step <- prep$pd$step
  zm <- any(step == 0, na.rm = TRUE)
  nstarts <- settings$nstarts %||% 5
  set.seed(settings$seed %||% 1)
  t0 <- proc.time()[["elapsed"]]
  fits <- lapply(seq_len(nstarts), function(i) {
    sv <- start_values(step, K, jitter = i > 1, zm = zm)
    tryCatch(hmm_fit_once(prep, K, sv, retryFits = settings$retryFits %||% 0), error = function(e) NULL)
  })
  ok <- !vapply(fits, is.null, logical(1))
  if (!any(ok)) stop("every starting value failed; check the data (very few fixes, all steps zero?)")
  nll <- vapply(fits, function(f) if (is.null(f)) Inf else f$mod$minimum, 0)
  best <- fits[[which.min(nll)]]
  n_best <- sum(nll - min(nll) < 0.5)
  m <- order_states(best, prep, K)
  if (m$mod$minimum > min(nll) + 0.01) m <- best   # relabel refit drifted: keep original
  ord_ok <- !is.unsorted(m$mle$step["mean", ])
  nd <- settings$ndraws %||% 2000
  dr <- hmm_draws(m, K, nd, zm)
  bdim <- dim(m$mle$beta)
  terms <- prep$terms
  rn <- rownames(m$mle$beta)
  if (!identical(rn[-1], terms)) {
    # momentuHMM may rename formula terms; map by position (same order as formula)
    rn <- c("(Intercept)", terms)
  }
  nat_hat <- w2nat(dr$est, K, zm, bdim)
  # sanity: our working->natural map must reproduce momentuHMM's MLE
  map_err <- max(abs(nat_hat$mu - m$mle$step["mean", ]) / m$mle$step["mean", ],
                 abs(nat_hat$beta - m$mle$beta))
  nmT <- offdiag_names(K)
  qfun <- function(f) {
    hat <- f(nat_hat)
    if (!dr$ok) return(data.frame(estimate = hat, sd = NA_real_, lower = NA_real_, upper = NA_real_))
    D <- t(apply(dr$W, 1, function(w) f(w2nat(w, K, zm, bdim))))
    if (length(hat) == 1) D <- t(D)
    data.frame(estimate = hat, sd = apply(D, 2, stats::sd),
               lower = apply(D, 2, stats::quantile, 0.025), upper = apply(D, 2, stats::quantile, 0.975))
  }
  # slopes per ORIGINAL covariate unit; intercept at covariate means
  slope_scale <- c(1, vapply(terms, function(tm) {
    v <- sub("_z$", "", tm)
    if (grepl("_z$", tm)) prep$scaling[[v]]$scale else 1
  }, 0))
  beta_names <- as.vector(outer(c("intercept", sub("_z$", "", terms)), nmT,
                                function(a, b) paste0("tpm_", b, "_", a)))
  est <- rbind(
    cbind(param = paste0("step_mean_m_state", 1:K), qfun(function(n) n$mu)),
    cbind(param = paste0("step_sd_m_state", 1:K), qfun(function(n) n$sd)),
    if (zm) cbind(param = paste0("zero_step_prob_state", 1:K), qfun(function(n) n$zm)),
    cbind(param = paste0("turn_concentration_state", 1:K), qfun(function(n) n$kappa)),
    cbind(param = beta_names, qfun(function(n) as.vector(n$beta / slope_scale))),
    cbind(param = paste0("activity_budget_state", 1:K),
          qfun(function(n) activity_budget(n$beta, terms, K, prep$interval_h, prep$phase0)))
  )
  est$param <- as.character(est$param); rownames(est) <- NULL
  notes <- character()
  code <- m$mod$code %||% NA
  conv <- isTRUE(code %in% 1:3) && dr$ok
  if (!dr$ok) notes <- c(notes, "Hessian not invertible: standard errors unavailable (parameter on a boundary or not identifiable). Do not interpret.")
  if (isTRUE(dr$conditional)) notes <- c(notes, "The initial-state distribution is on a boundary (singular full Hessian); intervals condition on its estimate. This nuisance parameter does not affect the reported quantities.")
  if (length(dr$zm_boundary)) notes <- c(notes, sprintf("The probability of a zero-length step in state(s) %s is estimated at ~0 (on a boundary), so it is held fixed when computing intervals; the data contain %d zero-length step(s) in total.",
                                                        paste(dr$zm_boundary, collapse = ", "), sum(step == 0, na.rm = TRUE)))
  if (!isTRUE(code %in% 1:3)) notes <- c(notes, sprintf("Optimiser stopped with nlm code %s (not converged).", code))
  notes <- c(notes, sprintf("%d of %d starting values reached the best log-likelihood (%.2f) within 0.5 units%s.",
                            n_best, sum(ok), -min(nll),
                            if (n_best < 2 && nstarts > 1) "; only one start found it, so the likelihood surface is multimodal: increase nstarts and check" else ""))
  if (!ord_ok) notes <- c(notes, "States could not be relabelled in order of mean step length.")
  if (map_err > 1e-6) notes <- c(notes, sprintf("Internal parameter map mismatch %.2g (report this).", map_err))
  if (is.null(settings$interval_h %||% data$interval_h))
    notes <- c(notes, sprintf("Fix interval taken as %s (most common gap between fixes); if the collar schedule differs, set data: interval_h.",
                              if (prep$interval_h < 1) sprintf("%g min", round(prep$interval_h * 60, 2)) else sprintf("%g h", round(prep$interval_h, 3))))
  if (prep$n_missing > 0) notes <- c(notes, sprintf("%d of %d regular time steps had no fix (treated as missing; covariates carried forward).",
                                                    prep$n_missing, nrow(prep$grid)))
  for (v in prep$covs) {
    nmiss <- sum(is.na(data$track[[v]]) & !is.na(data$track$x))
    if (nmiss > 0) notes <- c(notes, sprintf("%d fixes with coordinates had no value of %s; the previous fix's value was carried forward.", nmiss, v))
  }
  if (prep$n_dropped > 0) notes <- c(notes, sprintf("%d fixes were off the %.2g-h grid by more than the tolerance and were dropped%s.",
                                                    prep$n_dropped, prep$interval_h,
                                                    if (prep$n_dropped > 0.05 * nrow(data$track)) "; this is a lot: check the fix schedule and set data: interval_h, or use settings: regularise: crawl" else ""))
  obj <- list(m = m, prep = prep, K = K, zm = zm, draws = dr, terms = terms, bdim = bdim,
              llk = -m$mod$minimum, n_best = n_best, nstarts = sum(ok), aic = 2 * m$mod$minimum + 2 * length(m$mod$estimate),
              time = proc.time()[["elapsed"]] - t0, state_names = settings$state_names)
  if (isTRUE(settings$checks %||% TRUE)) {
    obj$viterbi <- tryCatch(momentuHMM::viterbi(m), error = function(e) NULL)
    obj$pres <- tryCatch(suppressWarnings(momentuHMM::pseudoRes(m)), error = function(e) NULL)
    if (!is.null(obj$pres)) {
      sr <- obj$pres$stepRes[is.finite(obj$pres$stepRes)]
      ar <- obj$pres$angleRes[is.finite(obj$pres$angleRes)]
      a1 <- stats::acf(obj$pres$stepRes, na.action = stats::na.pass, plot = FALSE, lag.max = 1)$acf[2]
      obj$gof <- list(step_ks_p = suppressWarnings(stats::ks.test(sr, "pnorm")$p.value),
                      angle_ks_p = suppressWarnings(stats::ks.test(ar, "pnorm")$p.value),
                      step_acf1 = a1)
      notes <- c(notes, sprintf("Pseudo-residuals: step KS p = %.3f, angle KS p = %.3f, lag-1 autocorrelation of step residuals = %.2f.",
                                obj$gof$step_ks_p, obj$gof$angle_ks_p, a1))
    }
  }
  list(estimates = est,
       diagnostics = list(converged = conv, rhat_max = NA_real_, ess_min = NA_real_, notes = notes),
       object = obj)
}

read_data <- function(config) {
  d <- config$data
  f <- file.path(config$.dir %||% ".", d$file)
  raw <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  col <- function(k, def) d[[k]] %||% def
  need <- c(col("id", "id"), col("time", "timestamp"),
            if (!is.null(d$lon) || !is.null(d$lat)) c(d$lon %||% "lon", d$lat %||% "lat") else c(col("x", "x"), col("y", "y")))
  miss <- setdiff(need, names(raw))
  if (length(miss)) stop("columns not found in ", d$file, ": ", paste(miss, collapse = ", "),
                         ". Available: ", paste(names(raw), collapse = ", "),
                         ". Fix the column name in the config (data: id / time / lon / lat / x / y) or in the file.")
  tz <- d$tz %||% "UTC"
  traw <- raw[[col("time", "timestamp")]]
  fmt <- d$time_format %||% "%Y-%m-%d %H:%M:%S"
  tm <- as.POSIXct(traw, format = fmt, tz = tz)
  if (all(is.na(tm))) {
    # Common exports the default format misses. ISO 8601 ("2024-06-01T12:00:33",
    # optional fractional seconds; a trailing "Z" means UTC and is converted to
    # the local tz). Tried before as.POSIXct's own guesses, which silently fall
    # back to date-only parsing and drop the clock time (verification 2026-10-19).
    t2 <- sub("Z$", "", trimws(traw)); isz <- grepl("Z$", trimws(traw))
    for (f2 in c("%Y-%m-%dT%H:%M:%OS", "%Y-%m-%d %H:%M:%OS", "%Y/%m/%d %H:%M:%OS", "%Y-%m-%dT%H:%M", "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M")) {
      tu <- as.POSIXct(t2, format = f2, tz = "UTC"); tl <- as.POSIXct(t2, format = f2, tz = tz)
      cand <- tl; cand[isz] <- tu[isz]; attr(cand, "tzone") <- tz
      if (mean(!is.na(cand[nzchar(t2)])) > 0.5) { tm <- cand; break }
    }
  }
  if (all(is.na(tm))) stop(sprintf("could not read any timestamp in column '%s' with time_format '%s' (first value: '%s'). Set data: time_format to match the file, e.g. '%%m/%%d/%%Y %%H:%%M' for 06/01/2024 13:00 or '%%d.%%m.%%Y %%H:%%M:%%S' for 01.06.2024 13:00:00.",
                               col("time", "timestamp"), fmt, traw[which(nzchar(traw))[1]]))
  okt <- !is.na(tm)
  if (any(okt) && all(format(tm[okt], "%H:%M:%S") == "00:00:00") && any(grepl("[0-9]{1,2}:[0-9]{2}", traw[okt])))
    stop(sprintf("timestamps in column '%s' were read as dates only (the clock time was lost; first value: '%s'). Set data: time_format to match the file exactly, e.g. '%%Y-%%m-%%dT%%H:%%M:%%S' for 2024-06-01T13:00:00.",
                 col("time", "timestamp"), traw[which(okt)[1]]))
  if (!is.null(d$lon)) {
    eco_require("sf")
    ll <- data.frame(lon = suppressWarnings(as.numeric(raw[[d$lon]])), lat = suppressWarnings(as.numeric(raw[[d$lat]])))
    okll <- stats::complete.cases(ll)
    if (any(abs(ll$lat[okll]) > 90) || any(abs(ll$lon[okll]) > 180))
      stop(sprintf("column '%s' has values outside -90..90 (e.g. %g) or '%s' outside -180..180, so these are not WGS84 latitude/longitude. Are longitude and latitude swapped in the config? If the coordinates are already in metres (UTM), use data: x / y instead of lon / lat.",
                   d$lat, ll$lat[okll][which.max(abs(ll$lat[okll]))], d$lon))
    crs <- d$crs_out
    if (is.null(crs)) {
      zone <- floor((mean(ll$lon[okll]) + 180) / 6) + 1
      crs <- (if (mean(ll$lat[okll]) >= 0) 32600 else 32700) + zone
    }
    xy <- matrix(NA_real_, nrow(ll), 2)
    pts <- sf::st_as_sf(ll[okll, ], coords = c("lon", "lat"), crs = 4326)
    xy[okll, ] <- sf::st_coordinates(sf::st_transform(pts, crs))
    x <- xy[, 1]; y <- xy[, 2]; crs_txt <- paste0("EPSG:", crs)
  } else {
    x <- raw[[col("x", "x")]]; y <- raw[[col("y", "y")]]; crs_txt <- d$crs %||% "projected metres (as supplied)"
  }
  covs <- unlist(d$covariates %||% list())
  mc <- setdiff(covs, names(raw))
  if (length(mc)) stop(sprintf("covariate column(s) %s not found in %s. Available: %s. Check the spelling under data: covariates.",
                               paste0("'", mc, "'", collapse = ", "), d$file, paste(names(raw), collapse = ", ")))
  tr <- data.frame(ID = as.character(raw[[col("id", "id")]]), time = tm, x = x, y = y)
  tr$row_in_file <- seq_len(nrow(raw)) + 1L   # +1 for the header line
  for (v in covs) tr[[v]] <- raw[[v]]
  list(track = tr, covariates = covs, tod = !isFALSE(d$time_of_day), tz = tz,
       nstates = config$settings$nstates %||% 2, interval_h = d$interval_h,
       crs = crs_txt, raw_n = nrow(raw), file = d$file)
}

check_data <- function(data) {
  pr <- character()
  tr <- data$track
  if (anyNA(tr$time)) pr <- c(pr, sprintf("%d timestamps could not be read%s; fix or delete those rows, or set data: time_format: to match the file (e.g. '%%d/%%m/%%Y %%H:%%M').",
                                          sum(is.na(tr$time)),
                                          if (!is.null(tr$row_in_file)) paste0(" (file lines ", paste(utils::head(tr$row_in_file[is.na(tr$time)], 8), collapse = ", "), ")") else ""))
  for (v in data$covariates) {
    if (!v %in% names(tr)) { pr <- c(pr, sprintf("covariate '%s' is not a column of the track file (columns: %s); check the spelling under data: covariates.",
                                                 v, paste(setdiff(names(tr), c("ID", "time", "x", "y")), collapse = ", "))); next }
    if (!is.numeric(tr[[v]])) {
      bad <- unique(tr[[v]][is.na(suppressWarnings(as.numeric(tr[[v]]))) & !is.na(tr[[v]]) & nzchar(trimws(tr[[v]]))])
      pr <- c(pr, if (length(bad)) sprintf("covariate '%s' contains text that is not a number (%s); leave missing values empty or NA, and code categories as separate 0/1 columns.",
                                           v, paste0("'", utils::head(bad, 5), "'", collapse = ", "))
                  else sprintf("covariate '%s' is not numeric; code categories as 0/1 columns.", v))
    }
  }
  n <- table(tr$ID)
  if (any(n < 30)) {
    few <- names(n)[n < 30]; big <- names(n)[n >= 30]
    hint <- vapply(few, function(a) {
      if (!length(big)) return("")
      dd <- utils::adist(a, big)[1, ]
      if (min(dd) <= 2) sprintf(" ('%s' looks like a typo for '%s')", a, big[which.min(dd)]) else ""
    }, "")
    pr <- c(pr, sprintf("animal(s) %s have fewer than 30 fixes%s; HMMs need long tracks (correct misspelled ids, or remove these animals).",
                        paste(few, collapse = ", "), paste(hint, collapse = "")))
  }
  okt <- !is.na(tr$time)
  dup <- duplicated(tr[okt, c("ID", "time")])
  if (any(dup)) pr <- c(pr, sprintf("%d duplicated timestamps within an animal (same ID and time); remove duplicates first.", sum(dup)))
  if (all(okt)) {
    gaps <- lapply(split(as.numeric(tr$time), tr$ID), function(t) diff(sort(t)) / 3600)
    med <- vapply(gaps, function(g) if (length(g)) stats::median(g) else NA_real_, 0)
    if (any(med <= 0, na.rm = TRUE))
      pr <- c(pr, sprintf("animal(s) %s have most fixes at identical times (median gap 0); check that the time column includes the clock time and that ids are not merged.",
                          paste(names(med)[which(med <= 0)], collapse = ", ")))
    else if (sum(is.finite(med)) > 1 && (max(med, na.rm = TRUE) / min(med, na.rm = TRUE)) > 1.5)
      pr <- c(pr, sprintf("animals have different fix schedules (median intervals %s h); fit them separately or subsample to a common interval.",
                          paste(sprintf("%s=%.2g", names(med), med), collapse = ", ")))
  }
  if (sum(stats::complete.cases(tr[, c("x", "y")])) < 0.5 * nrow(tr))
    pr <- c(pr, "more than half of the fixes have no coordinates.")
  pr
}

module <- ecomodule(
  id = "movement-hmm",
  title = "Hidden Markov models for animal movement (behavioural states)",
  description = paste(
    "Which behavioural modes (e.g. resting/foraging versus travelling) does an animal switch between,",
    "how do they differ in step length and turning, how much time is spent in each, and how do time of day",
    "and habitat change the chance of switching between them? Fitted to regular-interval GPS tracks with a",
    "gamma step-length / von Mises turning-angle hidden Markov model whose transition probabilities depend on covariates."),
  defaults = list(
    params = list(step_mean = c(60, 450), step_sd = c(50, 300), angle_kappa = c(0.5, 2.5),
                  tpm_intercept = c(-1.6, -1.2), tpm_tod_cos = c(-0.9, 0.7), tpm_tod_sin = c(0.4, -0.3),
                  tpm_forest = c(-0.025, 0.025), forest_ref = 50),
    design = list(n_animals = 6, n_fixes = 600, fix_interval_h = 1, missing_prop = 0.05,
                  jitter_min = 2, landscape_km = 20)),
  simulate = simulate, fit = fit, read_data = read_data, check_data = check_data,
  interpret = interpret, figures = figures, plots = plots, methods_text = methods_text,
  settings = list(nstarts = 10, retryFits = 0, ndraws = 2000, checks = TRUE),
  fast_settings = list(nstarts = 4, retryFits = 0, ndraws = 1000, checks = FALSE),
  power_params = c("tpm_1to2_forest", "tpm_2to1_forest"),
  engine = "momentuHMM (maximum likelihood, nlm; fallback moveHMM for cross-checks)",
  references = c(
    "Langrock R, King R, Matthiopoulos J, Thomas L, Fortin D, Morales JM (2012) Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions. Ecology 93(11):2336-2342. doi:10.1890/11-2241.1",
    "McClintock BT, Michelot T (2018) momentuHMM: R package for generalized hidden Markov models of animal movement. Methods in Ecology and Evolution 9(6):1518-1530. doi:10.1111/2041-210X.12995",
    "Michelot T, Langrock R, Patterson TA (2016) moveHMM: an R package for the statistical modelling of animal movement data using hidden Markov models. Methods in Ecology and Evolution 7:1308-1315. doi:10.1111/2041-210X.12578",
    "McClintock BT, Langrock R, Gimenez O, Cam E, Borchers DL, Glennie R, Patterson TA (2020) Uncovering ecological state dynamics with hidden Markov models. Ecology Letters 23(12):1878-1903. doi:10.1111/ele.13610",
    "Pohle J, Langrock R, van Beest FM, Schmidt NM (2017) Selecting the number of states in hidden Markov models: pragmatic solutions illustrated using animal movement. Journal of Agricultural, Biological and Environmental Statistics 22(3):270-293. doi:10.1007/s13253-017-0283-8",
    "Johnson DS, London JM, Lea MA, Durban JW (2008) Continuous-time correlated random walk model for animal telemetry data. Ecology 89(5):1208-1215. doi:10.1890/07-1032.1",
    "Koslik J-O, Feldmann CC, Mews S, Michels R, Langrock R (2025) Inference on the state process of periodically inhomogeneous hidden Markov models for animal behavior. Annals of Applied Statistics 19(4):2724-2737. arXiv:2312.14583",
    "Zucchini W, MacDonald IL, Langrock R (2016) Hidden Markov Models for Time Series: An Introduction Using R, 2nd edn. Chapman & Hall/CRC."
  )
)
