# Engine for movement-ssf.
#
#   ssf_build_steps()  amt: regularise tracks, steps by burst, pooled tentative
#                      gamma / von Mises distributions, random (available)
#                      steps, covariates at step ends.
#   ssf_fit_issa()     glmmTMB Poisson formulation of the conditional logistic
#                      likelihood with step-specific intercepts whose variance
#                      is fixed at 10^6 and independent random slopes per
#                      animal (Muff et al. 2020); movement-kernel updates
#                      (Avgar et al. 2016; Fieberg et al. 2021) with delta-
#                      method intervals.
#   ssf_fit_rsf()      use-availability RSF: weighted logistic glmmTMB with
#                      available points in each animal's MCP (weight 5000,
#                      Fieberg et al. 2021), random intercept variance fixed
#                      large and random slopes (Muff et al. 2020).
#   ssf_home_ranges()  AKDE (ctmm, Fleming et al. 2015) + KDE and MCP (amt).
#   ssf_uhc()          out-of-sample used-habitat calibration check
#                      (Fieberg et al. 2018).

ssf_z <- function(land, v, x, y) {
  s <- land$scaling[[v]]
  (ssf_lookup(land, x, y, v) - s$center) / s$scale
}

ssf_interval_hours <- function(fixes) {
  d <- unlist(lapply(split(fixes$t, fixes$id), function(t) diff(as.numeric(sort(t)))))
  round(stats::median(d) / 60) / 60   # nearest minute, in hours
}

ssf_build_steps <- function(data, st) {
  eco_require(c("amt", "lubridate"))
  suppressPackageStartupMessages(requireNamespace("amt"))
  fx <- data$fixes
  h <- st$interval_hours %||% data$interval_hours %||% ssf_interval_hours(fx)
  tol <- st$tolerance_minutes %||% max(1, round(60 * h * 0.1))
  notes <- character()
  rate <- lubridate::minutes(round(60 * h)); tolp <- lubridate::minutes(tol)
  ids <- unique(fx$id)
  stp <- list(); trks <- list()
  for (i in ids) {
    d <- fx[fx$id == i, ]
    d <- d[order(d$t), ]
    tr <- amt::make_track(d, x, y, t)
    tr <- amt::track_resample(tr, rate = rate, tolerance = tolp)
    tr <- amt::filter_min_n_burst(tr, min_n = 3)
    trks[[i]] <- tr
    if (nrow(tr) < 3) { notes <- c(notes, sprintf("%s: fewer than 3 regular fixes; dropped.", i)); next }
    s <- amt::steps_by_burst(tr)
    s <- s[!is.na(s$ta_) & s$sl_ > 0, ]
    if (nrow(s) < (st$min_steps %||% 20)) {
      notes <- c(notes, sprintf("%s: only %d usable steps (< %d); dropped from the step-selection fit.", i, nrow(s), st$min_steps %||% 20))
      next
    }
    s$id <- i
    stp[[i]] <- s
  }
  if (!length(stp)) stop("no animal has enough regular steps; check interval_hours / tolerance_minutes")
  all_s <- do.call(rbind, lapply(stp, function(s) data.frame(sl = s$sl_, ta = s$ta_)))
  sl_d <- suppressMessages(amt::fit_distr(all_s$sl, "gamma"))
  ta_d <- amt::fit_distr(all_s$ta, "vonmises")
  set.seed(st$seed %||% 1)
  rs <- lapply(names(stp), function(i) {
    s <- stp[[i]]
    r <- suppressWarnings(amt::random_steps(s, n_control = st$n_control, sl_distr = sl_d, ta_distr = ta_d))
    r <- as.data.frame(r)
    r$id <- i
    r$stratum <- paste(i, r$step_id_, sep = "_")
    r
  })
  rs <- do.call(rbind, rs)
  land <- data$land
  for (v in data$covariates) rs[[paste0("z_", v)]] <- ssf_z(land, v, rs$x2_, rs$y2_)
  for (v in data$covariates) rs[[paste0("start_", v)]] <- ssf_lookup(land, rs$x1_, rs$y1_, v)
  cen <- do.call(rbind, lapply(split(fx, fx$id), function(d) data.frame(id = d$id[1], cx = mean(d$x), cy = mean(d$y))))
  m <- match(rs$id, cen$id)
  # Range residency: -a * |end - c_i|^2 with an unknown centre c_i expands to
  # -a * (hx^2 + hy^2) + 2a c_ix hx + 2a c_iy hy + const (const cancels within a
  # stratum), so a common coefficient on home_d2 plus animal-specific linear
  # terms in hx, hy estimates the centre jointly with selection. hx, hy are
  # km from the animal's mean location (only for numerical conditioning).
  rs$hx <- (rs$x2_ - cen$cx[m]) / 1000; rs$hy <- (rs$y2_ - cen$cy[m]) / 1000
  rs$home_d2 <- rs$hx^2 + rs$hy^2
  rs$sl_km <- rs$sl_ / 1000
  rs$log_sl <- log(rs$sl_)
  rs$cos_ta <- cos(rs$ta_)
  rs$case <- as.integer(rs$case_)
  zc <- paste0("z_", data$covariates)
  bad_used <- unique(rs$stratum[rs$case == 1 & !stats::complete.cases(rs[, zc, drop = FALSE])])
  n_av_na <- sum(rs$case == 0 & !stats::complete.cases(rs[, zc, drop = FALSE]))
  if (length(bad_used)) notes <- c(notes, sprintf("%d used steps ended outside the covariate rasters (or on NA cells) and were dropped.", length(bad_used)))
  if (n_av_na) notes <- c(notes, sprintf("%d available steps (%.1f%%) fell outside the rasters and were dropped; extend the rasters if this is large.", n_av_na, 100 * n_av_na / sum(rs$case == 0)))
  rs <- rs[!(rs$stratum %in% bad_used) & stats::complete.cases(rs[, zc, drop = FALSE]), ]
  list(steps = rs, sl_distr = sl_d, ta_distr = ta_d, interval_hours = h, tolerance_minutes = tol,
       tracks = trks, centres = cen, notes = notes,
       n_used = tapply(rs$case, rs$id, sum))
}

# Fit a glmmTMB model with the first random-effect SD fixed (Muff et al. 2020).
ssf_tmb <- function(fo, d, family, n_re, weights = NULL, fix_first_sd = TRUE, REML = FALSE) {
  eco_require("glmmTMB")
  args <- list(formula = fo, data = d, family = family, REML = isTRUE(REML))
  if (!is.null(weights)) args$weights <- weights
  if (fix_first_sd) {
    args$map <- list(theta = factor(c(NA, seq_len(n_re - 1))))
    args$start <- list(theta = c(log(1e3), rep(0, n_re - 1)))
  }
  w <- character()
  m <- withCallingHandlers(do.call(glmmTMB::glmmTMB, args),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  attr(m, "eco_warnings") <- unique(w)
  m
}

ssf_theta <- function(m) {
  pf <- m$sdr$par.fixed; V <- m$sdr$cov.fixed
  k <- which(names(pf) == "theta")
  list(est = unname(pf[k]), se = sqrt(pmax(diag(V)[k], 0)))
}

ssf_home_terms <- function(st, n_id) {
  if (!isTRUE(st$home_attraction)) return(character())
  if (identical(st$home_centre, "mean") ) return("home_d2")
  if (n_id > 1) c("home_d2", "id:hx", "id:hy") else c("home_d2", "hx", "hy")
}

ssf_fit_issa <- function(S, data, st) {
  covs <- data$covariates; zc <- paste0("z_", covs)
  d <- S$steps
  d$id <- factor(d$id); d$stratum <- factor(d$stratum)
  mv <- c("sl_km", "log_sl", "cos_ta")
  fx <- c(zc, mv, ssf_home_terms(st, nlevels(d$id)))
  rs <- isTRUE(st$random_slopes) && nlevels(d$id) >= 2
  re <- if (rs) paste0(" + (0 + ", zc, " | id)", collapse = "") else ""
  fo <- stats::as.formula(paste("case ~ -1 +", paste(fx, collapse = " + "), "+ (1 | stratum)", re))
  m <- ssf_tmb(fo, d, stats::poisson(), n_re = 1 + if (rs) length(zc) else 0, REML = isTRUE(st$reml) && rs)
  b <- glmmTMB::fixef(m)$cond
  V <- as.matrix(stats::vcov(m)$cond)
  se <- sqrt(pmax(diag(V), 0))
  z <- stats::qnorm(0.975)
  rows <- list()
  add <- function(p, e, s, lo, hi) rows[[length(rows) + 1]] <<- data.frame(param = p, estimate = e, sd = s, lower = lo, upper = hi)
  for (v in covs) {
    k <- paste0("z_", v)
    add(paste0("sel_", v, "_mean"), b[[k]], se[[k]], b[[k]] - z * se[[k]], b[[k]] + z * se[[k]])
  }
  if (rs) {
    th <- ssf_theta(m)
    for (j in seq_along(covs)) {
      e <- exp(th$est[j]); s <- th$se[j]
      add(paste0("sel_", covs[j], "_sd"), e, e * s, exp(th$est[j] - z * s), exp(th$est[j] + z * s))
    }
  }
  for (v in covs) {
    k <- paste0("z_", v)
    add(paste0("rss_", v, "_per_sd"), exp(b[[k]]), exp(b[[k]]) * se[[k]], exp(b[[k]] - z * se[[k]]), exp(b[[k]] + z * se[[k]]))
    u <- data$units[[v]]
    if (!is.null(u)) {
      f <- u / data$land$scaling[[v]]$scale
      add(paste0("rss_", v, "_per_unit"), exp(f * b[[k]]), exp(f * b[[k]]) * f * se[[k]],
          exp(f * (b[[k]] - z * se[[k]])), exp(f * (b[[k]] + z * se[[k]])))
    }
  }
  # movement kernel updates (Avgar et al. 2016; amt::update_gamma / update_vonmises)
  kt <- S$sl_distr$params$shape; stt <- S$sl_distr$params$scale; kap_t <- S$ta_distr$params$kappa
  bl <- b[["log_sl"]]; bs <- b[["sl_km"]] / 1000; bc <- b[["cos_ta"]]
  Vm <- V[mv, mv]
  shape <- kt + bl
  rate <- 1 / stt - bs
  notes <- character()
  if (rate <= 0 || shape <= 0) {
    notes <- c(notes, "The updated step-length distribution is not a valid gamma (shape or rate <= 0); movement parameters are not reported. Use more available steps or check the step-length data.")
  } else {
    scale <- 1 / rate; mean_sl <- shape * scale
    g_shape <- c(0, 1, 0); g_scale <- c(scale^2 / 1000, 0, 0); g_mean <- c(shape * scale^2 / 1000, scale, 0)
    dm <- function(g) sqrt(max(0, drop(t(g) %*% Vm %*% g)))
    pos <- function(p, e, g) { s <- dm(g); sl <- s / e; add(p, e, s, e * exp(-z * sl), e * exp(z * sl)) }
    pos("step_shape", shape, g_shape)
    pos("step_scale_m", scale, g_scale)
    pos("step_mean_m", mean_sl, g_mean)
  }
  kap <- kap_t + bc; sk <- se[["cos_ta"]]
  add("turn_kappa", kap, sk, kap - z * sk, kap + z * sk)
  if (kap < 0) notes <- c(notes, "Updated von Mises concentration is negative: animals tend to reverse direction between fixes.")
  if (isTRUE(st$home_attraction)) {
    add("home_attraction_per_km2", -b[["home_d2"]], se[["home_d2"]], -b[["home_d2"]] - z * se[["home_d2"]], -b[["home_d2"]] + z * se[["home_d2"]])
  }
  est <- do.call(rbind, rows)
  # individual-level coefficients (population mean + conditional mode)
  indiv <- NULL
  if (rs) {
    rr <- glmmTMB::ranef(m, condVar = TRUE)$cond$id
    cv <- attr(rr, "condVar")
    indiv <- do.call(rbind, lapply(seq_along(covs), function(j) {
      k <- zc[j]
      csd <- if (is.list(cv)) sqrt(as.vector(cv[[j]])) else if (!is.null(cv)) sqrt(cv[j, j, ]) else NA_real_
      e <- b[[k]] + rr[[k]]; s <- sqrt(se[[k]]^2 + csd^2)
      data.frame(id = rownames(rr), covariate = covs[j], estimate = e, lower = e - z * s, upper = e + z * s)
    }))
  }
  conv <- isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess) && all(is.finite(se))
  list(estimates = est, model = m, coef = b, vcov = V, indiv = indiv, converged = conv,
       notes = c(notes, if (!conv) sprintf("glmmTMB did not converge cleanly (convergence code %s, Hessian positive definite: %s).", m$fit$convergence, m$sdr$pdHess),
                 if (length(attr(m, "eco_warnings"))) paste("glmmTMB warnings:", paste(attr(m, "eco_warnings"), collapse = " | "))),
       random_slopes = rs, formula = fo)
}

# Use-availability RSF within each animal's 100% MCP.
ssf_fit_rsf <- function(data, st, prefix = "rsf_") {
  eco_require(c("amt", "sf", "glmmTMB"))
  covs <- data$covariates; zc <- paste0("z_", covs)
  set.seed((st$seed %||% 1) + 7)
  fx <- data$fixes
  pts <- do.call(rbind, lapply(split(fx, fx$id), function(d) {
    tr <- amt::make_track(d, x, y, t)
    hr <- amt::hr_mcp(tr, levels = 1)
    rp <- amt::random_points(hr, n = st$rsf_n_avail * nrow(d), presence = tr)
    data.frame(id = d$id[1], case = as.integer(rp$case_), x = rp$x_, y = rp$y_)
  }))
  for (v in covs) pts[[paste0("z_", v)]] <- ssf_z(data$land, v, pts$x, pts$y)
  pts <- pts[stats::complete.cases(pts[, zc, drop = FALSE]), ]
  pts$w <- ifelse(pts$case == 1, 1, st$rsf_weight)
  pts$id <- factor(pts$id)
  rs <- isTRUE(st$random_slopes) && nlevels(pts$id) >= 2
  re <- if (rs) paste0(" + (0 + ", zc, " | id)", collapse = "") else ""
  fo <- stats::as.formula(paste("case ~", paste(zc, collapse = " + "), "+ (1 | id)", re))
  m <- ssf_tmb(fo, pts, stats::binomial(), n_re = 1 + if (rs) length(zc) else 0, weights = pts$w)
  b <- glmmTMB::fixef(m)$cond; se <- sqrt(pmax(diag(as.matrix(stats::vcov(m)$cond)), 0))
  z <- stats::qnorm(0.975)
  est <- data.frame(param = paste0(prefix, covs, "_mean"), estimate = b[zc], sd = se[zc],
                    lower = b[zc] - z * se[zc], upper = b[zc] + z * se[zc], row.names = NULL)
  if (rs) {
    th <- ssf_theta(m)
    est <- rbind(est, data.frame(param = paste0(prefix, covs, "_sd"), estimate = exp(th$est), sd = exp(th$est) * th$se,
                                 lower = exp(th$est - z * th$se), upper = exp(th$est + z * th$se)))
  }
  conv <- isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess)
  list(estimates = est, model = m, converged = conv, n_used = sum(pts$case), n_avail = sum(pts$case == 0))
}

# Home ranges: AKDE (ctmm) with CI, plus KDE (reference bandwidth) and MCP (amt).
ssf_home_ranges <- function(data, st) {
  eco_require(c("ctmm", "amt", "sf"))
  lev <- st$hr_level %||% 0.95
  crs <- data$land$crs
  fx <- data$fixes
  one <- function(i) {
    d <- fx[fx$id == i, ]; d <- d[order(d$t), ]
    row <- data.frame(id = i, n_fixes = nrow(d), akde_km2 = NA, akde_lower = NA, akde_upper = NA,
                      kde_km2 = NA, mcp_km2 = NA, model = NA, n_eff = NA, note = "")
    tr <- amt::make_track(d, x, y, t, crs = crs)
    row$kde_km2 <- tryCatch(amt::hr_area(amt::hr_kde(tr, levels = lev), units = FALSE)$area / 1e6, error = function(e) NA)
    row$mcp_km2 <- tryCatch(amt::hr_area(amt::hr_mcp(tr, levels = lev), units = FALSE)$area / 1e6, error = function(e) NA)
    ak <- tryCatch({
      ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(d, coords = c("x", "y"), crs = crs), 4326))
      df <- data.frame(individual.local.identifier = i, timestamp = format(d$t, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                       location.long = ll[, 1], location.lat = ll[, 2])
      tel <- suppressMessages(suppressWarnings(ctmm::as.telemetry(df, timezone = "UTC", projection = sf::st_crs(crs)$proj4string)))
      gu <- ctmm::ctmm.guess(tel, interactive = FALSE)
      fm <- suppressWarnings(ctmm::ctmm.select(tel, gu, trace = 0))
      ud <- ctmm::akde(tel, fm)
      sm <- summary(ud, level.UD = lev, units = FALSE)
      ci <- sm$CI[grep("^area", rownames(sm$CI)), ] / 1e6
      pol <- ctmm::as.sf(ud, level.UD = lev)
      pol <- sf::st_transform(pol, crs)
      list(ci = ci, model = summary(fm)$name %||% "", neff = unname(sm$DOF["area"]), pol = pol)
    }, error = function(e) e)
    if (inherits(ak, "error")) row$note <- paste("AKDE failed:", conditionMessage(ak))
    else {
      row$akde_km2 <- ak$ci[2]; row$akde_lower <- ak$ci[1]; row$akde_upper <- ak$ci[3]
      row$model <- ak$model; row$n_eff <- ak$neff
      cc <- as.data.frame(sf::st_coordinates(ak$pol))
      fcol <- grep("^L[0-9]$", names(cc), value = TRUE)
      feat <- cc[[fcol[length(fcol)]]]
      nm <- as.character(ak$pol$name)[feat]
      cc$ci <- ifelse(grepl("low", nm), "lower", ifelse(grepl("high", nm), "upper", "estimate"))
      cc$grp <- interaction(cc$ci, do.call(paste, cc[fcol]), drop = TRUE)
      cc$id <- i
      pol <- cc[, c("X", "Y", "ci", "grp", "id")]
    }
    list(row = row, pol = if (inherits(ak, "error")) NULL else pol)
  }
  ids <- unique(fx$id)
  nc <- max(1, as.integer(st$cores %||% 1))
  res <- if (nc > 1) parallel::mclapply(ids, one, mc.cores = nc) else lapply(ids, one)
  list(table = do.call(rbind, lapply(res, `[[`, "row")), polygons = do.call(rbind, lapply(res, `[[`, "pol")), level = lev)
}

# Out-of-sample used-habitat calibration (Fieberg et al. 2018): fit a
# population-level conditional logistic model to ~80% of strata, then for the
# held-out strata choose one step per stratum with probability proportional
# to exp(x'beta*) for M draws beta* ~ N(beta_hat, V) and compare the
# distribution of each covariate at chosen steps with the observed used steps.
ssf_uhc <- function(S, data, st, fit) {
  eco_require("survival")
  d <- S$steps
  covs <- data$covariates; zc <- paste0("z_", covs)
  set.seed((st$seed %||% 1) + 11)
  str <- unique(d$stratum)
  test <- sample(str, max(10, round(0.2 * length(str))))
  tr <- d[!d$stratum %in% test, ]; te <- d[d$stratum %in% test, ]
  terms <- c(zc, "sl_km", "log_sl", "cos_ta", ssf_home_terms(st, length(unique(d$id))))
  d$id <- factor(d$id)
  # conditional logit = stratified Cox model with one event per stratum (exact with Breslow ties here)
  tr$.time <- 1
  fo <- stats::as.formula(paste("survival::Surv(.time, case) ~", paste(terms, collapse = " + "), "+ survival::strata(stratum)"))
  cm <- survival::coxph(fo, data = tr, ties = "breslow")
  b <- stats::coef(cm); V <- stats::vcov(cm)
  keep <- !is.na(b); V <- V[keep, keep]; b <- b[keep]
  M <- st$uhc_draws %||% 200
  L <- chol(V)
  X <- stats::model.matrix(stats::as.formula(paste("~ -1 +", paste(terms, collapse = " + "))), te)
  X <- X[, names(b), drop = FALSE]
  b[is.na(b)] <- 0
  grp <- as.integer(factor(te$stratum))
  orig <- function(v, zz) zz * data$land$scaling[[v]]$scale + data$land$scaling[[v]]$center
  sims <- matrix(NA, M, length(covs), dimnames = list(NULL, covs))
  dens <- list()
  for (k in seq_len(M)) {
    bb <- b + drop(stats::rnorm(length(b)) %*% L)
    eta <- drop(X %*% bb)
    ch <- vapply(split(seq_along(grp), grp), function(ix) { p <- exp(eta[ix] - max(eta[ix])); ix[sample.int(length(ix), 1, prob = p)] }, 1L)
    for (v in covs) {
      val <- orig(v, te[[paste0("z_", v)]][ch])
      sims[k, v] <- mean(val)
      if (k <= 50) dens[[length(dens) + 1]] <- data.frame(covariate = v, draw = k, value = val)
    }
  }
  obs <- vapply(covs, function(v) mean(orig(v, te[[paste0("z_", v)]][te$case == 1])), 0)
  avail <- vapply(covs, function(v) mean(orig(v, te[[paste0("z_", v)]][te$case == 0])), 0)
  tab <- data.frame(covariate = covs, observed_used_mean = obs, available_mean = avail,
                    predicted_mean = colMeans(sims),
                    pred_lower = apply(sims, 2, stats::quantile, 0.025), pred_upper = apply(sims, 2, stats::quantile, 0.975),
                    row.names = NULL)
  tab$p_two_sided <- vapply(covs, function(v) { q <- mean(sims[, v] >= obs[[v]]); 2 * min(q, 1 - q) }, 0)
  tab$ok <- tab$observed_used_mean >= tab$pred_lower & tab$observed_used_mean <= tab$pred_upper
  used <- do.call(rbind, lapply(covs, function(v) data.frame(covariate = v, value = orig(v, te[[paste0("z_", v)]][te$case == 1]))))
  availd <- do.call(rbind, lapply(covs, function(v) data.frame(covariate = v, value = orig(v, te[[paste0("z_", v)]][te$case == 0]))))
  list(table = tab, sim = do.call(rbind, dens), used = used, avail = availd, n_test = length(test))
}
