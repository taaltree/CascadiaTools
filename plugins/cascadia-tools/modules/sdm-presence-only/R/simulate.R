# Simulation for sdm-presence-only.
#
# A gridded landscape (nx x ny cells of cell_km x cell_km) with two habitat
# covariates and one sampling-bias covariate (distance to the nearest road):
#   forest_cover_pct : smooth Gaussian random field mapped to 0-100 %
#   elevation_m      : smooth field, positively correlated with distance to
#                      road (roads follow valleys), so bias and habitat are
#                      confounded, the situation where ignoring bias fails
#   dist_road_km     : distance from the cell centre to the nearest of
#                      n_roads straight roads
# Species intensity (individuals per km2):
#   log lambda(s) = log(density_at_mean_covariates) + sum_k beta_k z_k(s)
# Presence-only records: an inhomogeneous Poisson process with intensity
#   lambda(s) b(s),  log b(s) = gamma + delta * z_road(s)
# with gamma set so that the expected number of records is n_presence
# (Fithian et al. 2015 eq. 1-2). z_k are covariates standardised over the
# landscape grid, exactly as fit() standardises them, so the slopes are
# "per 1 SD" on both sides.
# Optional survey sites (Koshkina et al. 2017; Dorazio 2014): a random sample
# of cells, each with a plot of survey_area_km2; N ~ Poisson(lambda A),
# occupied if N > 0, detected on each of survey_visits visits with
# probability detection_prob.

sdm_grf <- function(nx, ny, range_cells) {
  # Smooth Gaussian random field by FFT convolution of white noise with a
  # Gaussian kernel (padded so edges don't wrap). Returns an nx x ny matrix.
  pad <- ceiling(3 * range_cells)
  NX <- nx + 2 * pad; NY <- ny + 2 * pad
  w <- matrix(stats::rnorm(NX * NY), NX, NY)
  dx <- pmin(0:(NX - 1), NX - 0:(NX - 1)); dy <- pmin(0:(NY - 1), NY - 0:(NY - 1))
  k <- exp(-outer(dx^2, dy^2, "+") / (2 * range_cells^2))
  f <- Re(stats::fft(stats::fft(w) * stats::fft(k), inverse = TRUE))
  f[(pad + 1):(pad + nx), (pad + 1):(pad + ny)]
}

sdm_zs <- function(v) (v - mean(v)) / stats::sd(v)

sdm_landscape <- function(design) {
  nx <- design$nx; ny <- design$ny; ck <- design$cell_km
  xs <- (seq_len(nx) - 0.5) * ck; ys <- (seq_len(ny) - 0.5) * ck
  g <- expand.grid(x = xs, y = ys)               # x fastest: cell = ix + (iy-1)*nx
  rc <- design$range_km / ck
  # roads: straight lines through random points at random angles
  W <- nx * ck; H <- ny * ck
  d <- rep(Inf, nrow(g))
  for (r in seq_len(design$n_roads)) {
    px <- stats::runif(1, 0.15, 0.85) * W; py <- stats::runif(1, 0.15, 0.85) * H
    th <- stats::runif(1, 0, pi)
    d <- pmin(d, abs(-(g$x - px) * sin(th) + (g$y - py) * cos(th)))
  }
  dz <- sdm_zs(sqrt(d))                          # sqrt: less skewed driver of valleys
  rho <- design$bias_habitat_cor
  e1 <- sdm_zs(as.vector(sdm_grf(nx, ny, rc)))
  e2 <- sdm_zs(as.vector(sdm_grf(nx, ny, rc)))
  elev <- 900 + 280 * sdm_zs(rho * dz + sqrt(1 - rho^2) * e1)
  forest <- 100 * stats::plogis(1.3 * sdm_zs(0.5 * rho * dz + sqrt(1 - 0.25 * rho^2) * e2) + 0.4)
  g$forest_cover_pct <- round(forest, 1)
  g$elevation_m <- round(elev)
  g$dist_road_km <- round(d, 3)
  g$area_km2 <- ck^2
  g
}

sdm_simulate <- function(params, design, seed) {
  set.seed(seed)
  g <- sdm_landscape(design)
  hab <- c("forest_cover_pct", "elevation_m"); bias <- "dist_road_km"
  Z <- sapply(c(hab, bias), function(v) sdm_zs(g[[v]]))
  # habitat linear predictor (quadratic terms if params name *_sq)
  eta <- rep(log(params$density_at_mean_covariates), nrow(g))
  quad <- character()
  for (v in hab) {
    b <- params[[paste0("intensity_beta_", v)]]; if (!is.null(b)) eta <- eta + b * Z[, v]
    b2 <- params[[paste0("intensity_beta_", v, "_sq")]]
    if (!is.null(b2)) { eta <- eta + b2 * Z[, v]^2; quad <- c(quad, v) }
  }
  lam <- exp(eta)
  zeta <- params$bias_beta_dist_road_km * Z[, bias]
  A <- g$area_km2
  c0 <- log(design$n_presence / sum(A * lam * exp(zeta)))
  mu <- A * lam * exp(c0 + zeta)
  n <- stats::rpois(nrow(g), mu)
  ck <- design$cell_km
  idx <- rep(seq_len(nrow(g)), n)
  pres <- data.frame(x = g$x[idx] + stats::runif(length(idx), -ck / 2, ck / 2),
                     y = g$y[idx] + stats::runif(length(idx), -ck / 2, ck / 2))
  survey <- NULL
  if (isTRUE(design$survey_sites > 0)) {
    ns <- design$survey_sites; J <- design$survey_visits
    cells <- sample.int(nrow(g), ns)
    As <- rep(design$survey_area_km2, ns)
    N <- stats::rpois(ns, lam[cells] * As)
    p <- if (J > 1) params$detection_prob else 1
    y <- matrix(stats::rbinom(ns * J, 1, p * (N > 0)), ns, J)
    survey <- list(sites = data.frame(site = sprintf("S%03d", seq_len(ns)), x = g$x[cells], y = g$y[cells],
                                      area_km2 = As),
                   y = y)
  }
  truth <- c()
  for (v in hab) {
    for (sfx in c("", "_sq")) {
      nm <- paste0("intensity_beta_", v, sfx)
      if (!is.null(params[[nm]])) truth[nm] <- params[[nm]]
    }
  }
  truth["bias_beta_dist_road_km"] <- params$bias_beta_dist_road_km
  if (!is.null(survey)) {
    truth["density_at_mean_covariates"] <- params$density_at_mean_covariates
    truth["expected_total_abundance"] <- sum(A * lam)
    truth["mean_occupancy_prob"] <- mean(-expm1(-lam * design$survey_area_km2))
    if (design$survey_visits > 1) truth["detection_prob"] <- params$detection_prob
  }
  data <- list(grid = g, presence = pres, survey = survey,
               habitat = hab, bias = bias, quadratic = quad,
               coord_unit = "km",
               labels = list(forest_cover_pct = "Forest cover (%)", elevation_m = "Elevation (m)",
                             dist_road_km = "Distance to road (km)"))
  list(data = data, truth = truth)
}
