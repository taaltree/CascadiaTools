# Simulation for movement-ssf: telemetry from a known integrated step-selection
# process on a simulated landscape.
#
# Landscape: smooth Gaussian random fields (FFT convolution, periodic) on a
# grid of `cell_m` metres -> "forest" (% canopy cover, 0-100) and "elevation"
# (m). Further covariates can be requested with design$covariates; every layer
# is standardised with its landscape-wide mean and SD (exactly as fit() does
# with a researcher's raster), so selection coefficients are per 1 SD of the
# raster.
#
# Movement (Avgar et al. 2016; Signer et al. 2017): at every step an animal
# draws n_candidates candidate end points from the selection-free movement
# kernel  step length ~ Gamma(shape, scale),  turn angle ~ von Mises(0, kappa)
# and picks one with probability proportional to
#   w = exp( sum_k beta_ik * z_k(end) - attraction * d2(end) )
# where z_k is the standardised covariate, beta_ik ~ Normal(mu_k, sd_k) is
# animal i's selection coefficient, and d2 is the squared distance (km^2) to
# the animal's home-range centre (range residency). With many candidates this
# is a Monte Carlo version of the iSSF redistribution kernel. Fixes are
# recorded every `interval_hours`, then a share (1 - fix_rate) of fixes is
# lost at random, as with real GPS collars.

ssf_grf <- function(nr, nc, range_cells) {
  # Gaussian random field with Gaussian covariance, sd of kernel = range_cells
  z <- matrix(stats::rnorm(nr * nc), nr, nc)
  kx <- pmin(0:(nr - 1), nr - 0:(nr - 1)); ky <- pmin(0:(nc - 1), nc - 0:(nc - 1))
  k <- exp(-outer(kx^2, ky^2, "+") / (2 * range_cells^2))
  f <- Re(stats::fft(stats::fft(z) * stats::fft(k), inverse = TRUE))
  (f - mean(f)) / stats::sd(as.vector(f))
}

ssf_sim_landscape <- function(design) {
  nr <- nc <- as.integer(design$land_cells)
  res <- design$cell_m
  rc <- design$field_range_m / res
  z_el <- ssf_grf(nr, nc, 2 * rc)
  z_fo <- 0.9 * ssf_grf(nr, nc, rc) - 0.3 * z_el
  z_fo <- (z_fo - mean(z_fo)) / stats::sd(as.vector(z_fo))
  layers <- list(forest = 100 * stats::plogis(1.3 * z_fo + 0.3),
                 elevation = 1250 + 220 * z_el)
  extra <- setdiff(design$covariates, names(layers))
  for (v in extra) layers[[v]] <- ssf_grf(nr, nc, rc)
  layers <- layers[design$covariates]
  list(xmin = design$x0, ymax = design$y0 + nr * res, resx = res, resy = res,
       nrow = nr, ncol = nc, layers = layers, crs = design$crs,
       scaling = lapply(layers, function(m) list(center = mean(m), scale = stats::sd(as.vector(m)))))
}

# Fast covariate lookup on the stored landscape (NA outside).
ssf_cell_index <- function(land, x, y) {
  col <- floor((x - land$xmin) / land$resx) + 1
  row <- floor((land$ymax - y) / land$resy) + 1
  ok <- !is.na(col) & !is.na(row) & col >= 1 & col <= land$ncol & row >= 1 & row <= land$nrow
  idx <- rep(NA_integer_, length(x))
  idx[ok] <- as.integer((col[ok] - 1) * land$nrow + row[ok])
  idx
}
ssf_lookup <- function(land, x, y, layer) land$layers[[layer]][ssf_cell_index(land, x, y)]

ssf_simulate <- function(params, design, seed) {
  set.seed(seed)
  covs <- design$covariates
  na <- as.integer(design$n_animals); ns <- as.integer(design$n_steps)
  land <- ssf_sim_landscape(design)
  Z <- lapply(covs, function(v) (land$layers[[v]] - land$scaling[[v]]$center) / land$scaling[[v]]$scale)
  names(Z) <- covs
  mu <- unlist(params$sel_mean)[covs]; sdv <- unlist(params$sel_sd)[covs]
  B <- sapply(covs, function(v) stats::rnorm(na, mu[[v]], sdv[[v]]))
  B <- matrix(B, na, length(covs), dimnames = list(NULL, covs))
  # home-range centres inside the central part of the landscape
  W <- land$ncol * land$resx; H <- land$nrow * land$resy
  mar <- design$centre_margin_m
  cx <- land$xmin + stats::runif(na, mar, W - mar); cy <- land$ymax - stats::runif(na, mar, H - mar)
  nc <- as.integer(design$n_candidates)
  nb <- as.integer(design$burnin_steps)
  x <- cx; y <- cy; head <- stats::runif(na, -pi, pi)
  X <- matrix(NA_real_, na, ns + 1); Y <- X
  att <- params$attraction %||% 0
  for (s in seq_len(nb + ns + 1)) {
    sl <- matrix(stats::rgamma(na * nc, shape = params$step_shape, scale = params$step_scale), na, nc)
    ta <- matrix(ssf_rvonmises(na * nc, params$turn_kappa), na, nc)
    dir <- head + ta
    ex <- x + sl * cos(dir); ey <- y + sl * sin(dir)
    idx <- ssf_cell_index(land, as.vector(ex), as.vector(ey))
    lw <- matrix(0, na, nc)
    for (v in covs) lw <- lw + B[, v] * matrix(Z[[v]][idx], na, nc)
    lw <- lw - att * ((ex - cx)^2 + (ey - cy)^2) / 1e6
    lw[is.na(lw)] <- -Inf                  # off the landscape: never chosen
    pick <- apply(lw, 1, function(l) { p <- exp(l - max(l)); sample.int(nc, 1, prob = p) })
    ii <- cbind(seq_len(na), pick)
    x <- ex[ii]; y <- ey[ii]; head <- dir[ii]
    if (s > nb) { X[, s - nb] <- x; Y[, s - nb] <- y }
  }
  ids <- sprintf("%s%02d", design$id_prefix, seq_len(na))
  t0 <- as.POSIXct(design$start_time, tz = "UTC")
  dt <- design$interval_hours * 3600
  fixes <- do.call(rbind, lapply(seq_len(na), function(i) {
    tt <- t0 + (0:ns) * dt + stats::runif(ns + 1, 0, design$time_jitter_s)
    data.frame(id = ids[i], t = tt, x = X[i, ], y = Y[i, ], stringsAsFactors = FALSE)
  }))
  if (design$fix_rate < 1) fixes <- fixes[stats::runif(nrow(fixes)) < design$fix_rate, ]
  rownames(fixes) <- NULL
  truth <- c(stats::setNames(mu, paste0("sel_", covs, "_mean")),
             stats::setNames(sdv, paste0("sel_", covs, "_sd")),
             stats::setNames(exp(mu), paste0("rss_", covs, "_per_sd")),
             step_shape = params$step_shape, step_scale_m = params$step_scale,
             step_mean_m = params$step_shape * params$step_scale,
             turn_kappa = params$turn_kappa)
  if (att > 0) truth <- c(truth, home_attraction_per_km2 = att)
  data <- list(fixes = fixes, land = land, covariates = covs,
               labels = list(forest = "Forest cover (%)", elevation = "Elevation (m)"),
               names = list(forest = "forest cover", elevation = "elevation"),
               units = list(forest = 10, elevation = 100),
               unit_text = list(forest = "10 percentage points of forest cover", elevation = "100 m of elevation"),
               interval_hours = design$interval_hours, animal_label = design$animal_label,
               truth_individual = data.frame(id = ids, B, centre_x = cx, centre_y = cy))
  list(data = data, truth = truth)
}

# von Mises(0, kappa) draws (Best & Fisher 1979), vectorised.
ssf_rvonmises <- function(n, kappa) {
  if (kappa < 1e-8) return(stats::runif(n, -pi, pi))
  out <- numeric(n); todo <- seq_len(n)
  tau <- 1 + sqrt(1 + 4 * kappa^2); rho <- (tau - sqrt(2 * tau)) / (2 * kappa); r <- (1 + rho^2) / (2 * rho)
  while (length(todo)) {
    m <- length(todo)
    u1 <- stats::runif(m); u2 <- stats::runif(m); u3 <- stats::runif(m)
    z <- cos(pi * u1); f <- (1 + r * z) / (r + z); c <- kappa * (r - f)
    acc <- (c * (2 - c) - u2 > 0) | (log(c / u2) + 1 - c >= 0)
    th <- sign(u3 - 0.5) * acos(f)
    out[todo[acc]] <- th[acc]
    todo <- todo[!acc]
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
