# scr-search: maximum-likelihood engines in secr (Efford 2011; secr >= 5).
#
#  secr_transect : INDEPENDENT ROUTE. Each search track is a secr "transect"
#                  detector (Efford 2011; secr Appendix D), used only on its
#                  occasion; samples enter with their exact locations (snapped
#                  onto the track); Poisson counts (binomN = 0), hazard
#                  half-normal (HHN). The expected number of samples is the
#                  hazard integrated along the track: secr normalises the
#                  kernel along the line by sqrt(2 pi) sigma, so
#                  lambda0_secr = lambda0_per_km * sqrt(2 pi) * sigma_km
#                  (checked empirically with sim.capthist; see VALIDATION.md).
#                  No cell discretisation, no local evaluation.
#  secr_cells    : the SAME discretised model as the Nimble engine (one
#                  "count" detector per (cell, track) unit at the detector location used by Nimble (centroid or centre),
#                  usage = km searched), fitted by ML: checks the Nimble
#                  likelihood and its local-evaluation truncation.

ss_secr_mask <- function(prep, f) {
  msk <- secr::read.mask(data = data.frame(x = prep$mxy[, 1] * f, y = prep$mxy[, 2] * f), spacing = prep$pix * f)
  if (length(prep$dcov)) secr::covariates(msk) <- as.data.frame(prep$XD)
  msk
}

ss_secr_rows <- function(prep, fit, lam_fun, f) {
  cf <- stats::coef(fit); z <- stats::qnorm(0.975)
  rn <- secr::region.N(fit, se.N = TRUE)
  aconv <- ss_area_conv(prep)
  EN <- rn["E.N", ]; RN <- rn["R.N", ]
  rows <- list(
    data.frame(param = ss_density_name(prep$units), estimate = EN$estimate / prep$area * aconv, sd = EN$SE.estimate / prep$area * aconv,
               lower = EN$lcl / prep$area * aconv, upper = EN$ucl / prep$area * aconv),
    data.frame(param = "abundance", estimate = RN$estimate, sd = RN$SE.estimate, lower = RN$lcl, upper = RN$ucl))
  # lambda0 and derived quantities by the delta method on the link scale (log lambda0, log sigma)
  V <- fit$beta.vcv; nm <- rownames(cf)
  il <- match("lambda0", nm); is <- match("sigma", nm)
  bl <- cf$beta[il]; bs <- cf$beta[is]
  vl <- V[il, il]; vs <- V[is, is]; cls <- V[il, is]
  lam <- lam_fun(bl, bs)   # list(est, grad wrt (bl, bs))
  se_l <- sqrt(max(0, c(lam$g %*% matrix(c(vl, cls, cls, vs), 2) %*% lam$g)))
  rows[[3]] <- data.frame(param = ss_lam_name(prep$units), estimate = lam$est, sd = se_l,
                          lower = lam$est * exp(-z * se_l / lam$est), upper = lam$est * exp(z * se_l / lam$est))
  sg <- function(v) exp(v) / f * prep$cscale
  rows[[4]] <- data.frame(param = "sigma", estimate = sg(bs), sd = sg(bs) * sqrt(vs), lower = sg(bs - z * sqrt(vs)), upper = sg(bs + z * sqrt(vs)))
  hr <- function(v) pi * (sqrt(stats::qchisq(0.95, 2)) * exp(v) / f)^2 * ss_hr_conv(prep)
  rows[[5]] <- data.frame(param = "home_range_area_95", estimate = hr(bs), sd = 2 * hr(bs) * sqrt(vs),
                          lower = hr(bs - z * sqrt(vs)), upper = hr(bs + z * sqrt(vs)))
  # samples from a straight track through the centre = lambda0 sqrt(2 pi) sigma
  lt <- log(lam$est / ss_len_conv(prep)) + log(sqrt(2 * pi)) + (bs - log(f))
  gl <- c(lam$g[1] / lam$est, lam$g[2] / lam$est + 1)
  se_t <- sqrt(max(0, c(gl %*% matrix(c(vl, cls, cls, vs), 2) %*% gl)))
  rows[[6]] <- data.frame(param = "samples_per_track_through_centre", estimate = exp(lt),
                          sd = exp(lt) * se_t, lower = exp(lt - z * se_t), upper = exp(lt + z * se_t))
  for (v in prep$dcov) { r <- cf[paste0("D.", v), ]
    rows[[length(rows) + 1]] <- data.frame(param = paste0("density_beta_", v), estimate = r$beta, sd = r$SE.beta,
                                           lower = r$beta - z * r$SE.beta, upper = r$beta + z * r$SE.beta) }
  for (v in prep$wcov) { r <- cf[paste0("lambda0.", v), ]
    rows[[length(rows) + 1]] <- data.frame(param = paste0("detection_beta_", v), estimate = r$beta, sd = r$SE.beta,
                                           lower = r$beta - z * r$SE.beta, upper = r$beta + z * r$SE.beta) }
  est <- do.call(rbind, rows); rownames(est) <- NULL
  est
}

ss_secr_fit <- function(ch, msk, prep, start, method, timecov = NULL) {
  fD <- stats::as.formula(paste("D ~", if (length(prep$dcov)) paste(prep$dcov, collapse = " + ") else "1"))
  fl <- stats::as.formula(paste("lambda0 ~", if (length(prep$wcov)) paste(prep$wcov, collapse = " + ") else "1"))
  suppressWarnings(secr::secr.fit(ch, mask = msk, model = list(fD, fl), detectfn = "HHN", binomN = 0,
                                  trace = FALSE, start = start, method = method, ncores = 1, timecov = timecov))
}

ss_engine_secr_cells <- function(prep, settings) {
  eco_require("secr"); suppressPackageStartupMessages(requireNamespace("secr"))
  if (prep$observation != "poisson") stop("the secr_cells engine supports observation = poisson only")
  f <- 1000
  ud <- prep$units_df
  # one "count" detector per cell; every search unit's track is an "occasion"
  # (usage = length searched in the cell on that track); track-level
  # detection covariates become time covariates. Cell-level detection
  # covariates (from the mask) become detector covariates.
  tkey <- if (all(is.na(ud$track_id))) rep("all", nrow(ud)) else ud$track_id
  tl <- unique(tkey); nt <- length(tl)
  tdf <- data.frame(x = prep$cxy[, 1] * f, y = prep$cxy[, 2] * f)
  rownames(tdf) <- paste0("C", seq_len(prep$C))
  tr <- secr::read.traps(data = tdf, detector = "count", binary.usage = FALSE)
  us <- matrix(0, prep$C, nt)
  us[cbind(prep$ucell, match(tkey, tl))] <- exp(prep$logE)
  secr::usage(tr) <- us
  tcov <- NULL; dvars <- character()
  if (length(prep$wcov)) {
    first <- match(tl, tkey)
    constant_in_track <- vapply(seq_along(prep$wcov), function(k) all(tapply(prep$W[, k], tkey, function(v) length(unique(v)) == 1)), TRUE)
    if (all(constant_in_track)) {
      tcov <- as.data.frame(prep$W[first, , drop = FALSE])
    } else {
      constant_in_cell <- vapply(seq_along(prep$wcov), function(k) all(tapply(prep$W[, k], prep$ucell, function(v) length(unique(v)) == 1)), TRUE)
      if (!all(constant_in_cell)) stop("secr_cells: each detection covariate must be constant within a track or within a cell")
      dc <- as.data.frame(prep$W[match(seq_len(prep$C), prep$ucell), , drop = FALSE])
      secr::covariates(tr) <- dc
    }
  }
  sm <- prep$samples
  socc <- match(if (all(is.na(ud$track_id))) rep("all", nrow(sm)) else ud$track_id[match(sm$unit, ud$unit)], tl)
  ch <- secr::make.capthist(data.frame(Session = 1, ID = sm$individual, Occasion = socc,
                                       TrapID = paste0("C", ud$cellkey[match(sm$unit, ud$unit)] |> match(prep$ckeys))),
                            tr, fmt = "trapID", noccasions = nt)
  msk <- ss_secr_mask(prep, f)
  ini <- ss_inits(prep, prep$M)
  start <- list(D = max(prep$n, 1) * 2 / (prep$area * f^2 / 1e4), lambda0 = ini$lam0, sigma = ini$sigma * f)
  t0 <- proc.time()[["elapsed"]]
  fit <- ss_secr_fit(ch, msk, prep, start, if (is.null(settings$secr_method)) "Newton-Raphson" else settings$secr_method,
                     timecov = tcov)
  lam_fun <- function(bl, bs) list(est = exp(bl) * ss_len_conv(prep), g = c(exp(bl) * ss_len_conv(prep), 0))
  est <- ss_secr_rows(prep, fit, lam_fun, f)
  ok <- all(is.finite(est$sd)) && all(est$sd > 0)
  list(estimates = est, fit = fit, t_total = proc.time()[["elapsed"]] - t0, converged = ok, engine = "secr_cells")
}

# project points onto a polyline: returns the nearest point on the track
ss_project <- function(px, py, x, y) {
  ax <- x[-length(x)]; ay <- y[-length(y)]; bx <- x[-1]; by <- y[-1]
  dx <- bx - ax; dy <- by - ay; L2 <- dx^2 + dy^2; L2[L2 == 0] <- 1e-12
  out <- matrix(NA_real_, length(px), 2)
  for (k in seq_along(px)) {
    t <- pmin(1, pmax(0, ((px[k] - ax) * dx + (py[k] - ay) * dy) / L2))
    qx <- ax + t * dx; qy <- ay + t * dy
    j <- which.min((qx - px[k])^2 + (qy - py[k])^2)
    out[k, ] <- c(qx[j], qy[j])
  }
  out
}

ss_engine_secr_transect <- function(prep, settings) {
  eco_require("secr"); suppressPackageStartupMessages(requireNamespace("secr"))
  if (prep$effort_model != "track_length") stop("secr_transect uses the tracks themselves; effort_model must be track_length")
  if (prep$observation != "poisson") stop("secr_transect supports observation = poisson only")
  f <- 1000
  tr <- prep$tracks; tinfo <- prep$tinfo
  tid <- unique(tr$track_id)
  tl <- lapply(tid, function(t) { d <- tr[tr$track_id == t, ]
    xy <- data.frame(x = (d$x - prep$origin[1]) / prep$cscale * f, y = (d$y - prep$origin[2]) / prep$cscale * f)
    xy[c(TRUE, rowSums(abs(diff(as.matrix(xy)))) > 0), ] })
  # thin very dense GPS tracks to <= 400 vertices (integration cost scales with vertices)
  tl <- lapply(tl, function(d) if (nrow(d) > 400) d[unique(c(round(seq(1, nrow(d), length.out = 400)), nrow(d))), ] else d)
  names(tl) <- tid
  trs <- secr::make.transect(transectlist = tl, exclusive = FALSE)
  K <- prep$K
  occ <- tinfo$occasion[match(tid, tinfo$track_id)]
  us <- matrix(0, length(tid), K); us[cbind(seq_along(tid), occ)] <- 1
  secr::usage(trs) <- us
  if (length(prep$wcov)) {
    W <- prep$W[match(tid, prep$units_df$track_id), , drop = FALSE]
    secr::covariates(trs) <- as.data.frame(W)
  }
  sm <- prep$samples
  sx <- (sm$x - prep$origin[1]) / prep$cscale * f; sy <- (sm$y - prep$origin[2]) / prep$cscale * f
  pxy <- matrix(NA_real_, nrow(sm), 2)
  for (t in unique(sm$track_assigned)) { ii <- which(sm$track_assigned == t); d <- tl[[t]]
    pxy[ii, ] <- ss_project(sx[ii], sy[ii], d$x, d$y) }
  socc <- tinfo$occasion[match(sm$track_assigned, tinfo$track_id)]
  ch <- secr::make.capthist(data.frame(Session = 1, ID = sm$individual, Occasion = socc, x = pxy[, 1], y = pxy[, 2]),
                            trs, fmt = "XY", noccasions = K, snapXY = TRUE)
  msk <- ss_secr_mask(prep, f)
  ini <- ss_inits(prep, prep$M)
  sgm <- ini$sigma
  start <- list(D = max(prep$n, 1) * 2 / (prep$area * f^2 / 1e4), lambda0 = ini$lam0 * sqrt(2 * pi) * sgm, sigma = sgm * f)
  t0 <- proc.time()[["elapsed"]]
  fit <- ss_secr_fit(ch, msk, prep, start, if (is.null(settings$secr_method)) "Newton-Raphson" else settings$secr_method)
  # lambda0 per internal length = lambda0_secr / (sqrt(2 pi) sigma_internal); sigma_internal = exp(bs)/f
  lam_fun <- function(bl, bs) { e <- exp(bl) / (sqrt(2 * pi) * exp(bs) / f) * ss_len_conv(prep); list(est = e, g = c(e, -e)) }
  est <- ss_secr_rows(prep, fit, lam_fun, f)
  ok <- all(is.finite(est$sd)) && all(est$sd > 0)
  list(estimates = est, fit = fit, t_total = proc.time()[["elapsed"]] - t0, converged = ok, engine = "secr_transect")
}
