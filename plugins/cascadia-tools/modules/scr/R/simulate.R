# scr: simulate a camera-trap / hair-snare style SCR study (Efford 2004;
# Royle & Young 2008; Borchers & Efford 2008; Royle et al. 2014 ch. 9).
#
# Activity centres: inhomogeneous Poisson process on a rectangular state space
# (trap extent + buffer), optionally with a non-habitat lake removed,
# log D(s) = beta0 + beta_canopy * canopy_std(s) [+ beta_elev * elevation_std(s)];
# covariates are smooth Gaussian random fields on the habitat mask (canopy
# optionally with a centre-high trend so the trap array sits in good habitat).
# Detection on each occasion k at trap j (only while trap j operates, with
# effort u_jk):
#   binomial / Bernoulli (proximity): p = 1 - (1 - p0_ijk exp(-d^2 / 2 sigma_ij^2))^u_jk
#   Poisson (counts of samples):      y ~ Poisson(u_jk lambda0_ijk exp(-d^2 / 2 sigma_ij^2))
# link(p0_ijk) = link(p0) + beta_trail on_trail_j + beta_rain rain_k
#               + beta_male male_i + beta_b (caught at j before k)
# log sigma_ij = log sigma + beta_lure lure_j + beta_sigma_male male_i.
# The default design (no extras) reproduces the round-1 simulator exactly.

scr_smooth_field <- function(nx, ny, sd_px) {
  pad <- ceiling(3 * sd_px) + 1
  z <- matrix(stats::rnorm((nx + 2 * pad) * (ny + 2 * pad)), nx + 2 * pad, ny + 2 * pad)
  if (sd_px > 0) {
    k <- stats::dnorm(seq(-pad, pad), 0, sd_px); k <- k / sum(k)
    z <- apply(z, 2, function(c) stats::filter(c, k, circular = TRUE))
    z <- t(apply(z, 1, function(r) stats::filter(r, k, circular = TRUE)))
  }
  z[pad + seq_len(nx), pad + seq_len(ny)]
}

scr_zstd <- function(v) (v - mean(v)) / stats::sd(v)

scr_simulate <- function(params, design, seed) {
  set.seed(seed)
  isT <- function(k) isTRUE(design[[k]])
  obs <- if (is.null(design$obs_model)) "binomial" else design$obs_model
  pois <- obs == "poisson"
  side <- design$trap_grid_side; sp <- design$trap_spacing; K <- design$n_occasions
  buf <- design$buffer; pix <- design$mask_spacing
  g1 <- buf + sp * (seq_len(side) - 1)
  traps <- expand.grid(x = g1, y = g1)
  J <- nrow(traps)
  traps <- data.frame(id = sprintf("%s%02d", if (pois) "HS" else "CAM", seq_len(J)), x = traps$x, y = traps$y)
  has_trail <- isT("detection_covariate")
  if (has_trail) traps$on_trail <- as.integer(seq_len(J) %in% sample.int(J, round(design$prop_on_trail * J)))
  # trap operation: each trap fails from a random occasion onward with prob trap_failure
  usage <- matrix(1L, J, K, dimnames = list(traps$id, NULL))
  fail <- which(stats::runif(J) < design$trap_failure)
  for (j in fail) { k0 <- sample.int(K, 1); if (k0 > 1) usage[j, k0:K] <- 0L }
  # habitat mask (rectangle) and canopy-cover field
  xlim <- (side - 1) * sp + 2 * buf
  nxp <- round(xlim / pix)
  gx <- (seq_len(nxp) - 0.5) * pix
  mask <- expand.grid(x = gx, y = gx)
  G <- nrow(mask)
  f <- as.vector(scr_smooth_field(nxp, nxp, design$covariate_range / pix))
  f <- (f - mean(f)) / stats::sd(f)
  if (!is.null(design$covariate_trend) && design$covariate_trend != 0) {
    cen <- mean(g1)
    r <- sqrt((mask$x - cen)^2 + (mask$y - cen)^2)
    tr <- -(r - mean(r)) / stats::sd(r)
    f <- f + design$covariate_trend * tr
    f <- (f - mean(f)) / stats::sd(f)
  }
  mask$canopy <- 50 + 12 * f
  has_cov <- isTRUE(design$density_covariate)
  beta <- if (has_cov) params$density_beta_canopy else 0
  # ---- extras (each draws random numbers only when switched on)
  has_elev <- isT("density_covariate2")
  hab <- rep(TRUE, G)
  if (!is.null(design$nonhabitat) && design$nonhabitat > 0) {
    # a lake in the buffer on the west side: disk covering `nonhabitat` of the state space
    rad <- sqrt(design$nonhabitat * xlim^2 / pi)
    hab <- (mask$x - 0)^2 + (mask$y - xlim / 2)^2 > rad^2
  }
  lw <- beta * f
  if (has_elev) {
    f2 <- as.vector(scr_smooth_field(nxp, nxp, 2 * design$covariate_range / pix))
    f2 <- scr_zstd(f2[hab]); fe <- rep(NA_real_, G); fe[hab] <- f2
    mask$elevation <- 900 + 150 * fe
    lw <- lw + params$density_beta_elevation * ifelse(hab, fe, 0)
  }
  if (!all(hab)) {
    # re-standardise canopy over the habitat pixels so the per-SD effect is defined on the fitted mask
    fh <- scr_zstd(f[hab]); f[hab] <- fh; mask$canopy[hab] <- 50 + 12 * fh
    lw <- (if (has_cov) beta * f else 0) + (if (has_elev) params$density_beta_elevation * ifelse(hab, fe, 0) else 0)
  }
  a_km2 <- sum(hab) * pix^2 / 1e6
  EN <- params$density * a_km2
  w <- exp(lw) * hab; w <- w / sum(w)
  N <- stats::rpois(1, EN)
  cell <- sample.int(G, N, replace = TRUE, prob = w)
  s <- cbind(mask$x[cell] + stats::runif(N, -pix / 2, pix / 2), mask$y[cell] + stats::runif(N, -pix / 2, pix / 2))
  has_sex <- isT("sex"); has_rain <- isT("occasion_covariate"); has_beh <- isT("behaviour")
  has_lure <- isT("sigma_covariate"); cont_eff <- isT("continuous_effort")
  extras <- has_sex || has_rain || has_beh || has_lure || cont_eff || pois
  det0 <- if (pois) params$baseline_encounter_rate else params$baseline_detection
  link <- if (pois) log else stats::qlogis
  lp0 <- link(det0) + if (has_trail) params$detection_beta_on_trail * traps$on_trail else 0
  caps <- list()
  if (!extras) {
    p0j <- stats::plogis(lp0)
    d2 <- outer(s[, 1], traps$x, "-")^2 + outer(s[, 2], traps$y, "-")^2
    p <- sweep(exp(-d2 / (2 * params$sigma^2)), 2, p0j, "*")
    for (k in seq_len(K)) {
      yk <- matrix(stats::rbinom(N * J, 1, p * rep(usage[, k], each = N)), N, J)
      w1 <- which(yk == 1, arr.ind = TRUE)
      if (nrow(w1)) caps[[k]] <- data.frame(ind = w1[, 1], occasion = k, trap = traps$id[w1[, 2]])
    }
    sex <- NULL
  } else {
    if (has_lure) traps$lure <- as.integer(seq_len(J) %in% sample.int(J, round(J / 2)))
    if (cont_eff) {
      # effort = days a detector was active / nominal days (hair snares, scat stations), plus random inactive checks
      eff <- matrix(round(stats::runif(J * K, 0.5, 1.5), 1), J, K)
      off <- matrix(stats::runif(J * K) < (if (is.null(design$prop_inactive)) 0.1 else design$prop_inactive), J, K)
      usage <- usage * eff * (!off)
      dimnames(usage) <- list(traps$id, NULL)
    }
    rain <- NULL
    if (has_rain) { z <- scr_zstd(stats::rnorm(K)); rain <- 12 + 6 * z }
    sex <- if (has_sex) stats::rbinom(N, 1, 1 - params$proportion_female) else rep(0L, N)
    sig_i <- if (has_sex) ifelse(sex == 1, params$sigma_male, params$sigma) else rep(params$sigma, N)
    dm <- if (has_sex) (if (pois) params$baseline_encounter_rate_male else params$baseline_detection_male) else det0
    lsex <- if (has_sex) ifelse(sex == 1, link(dm) - link(det0), 0) else rep(0, N)
    sigij <- matrix(sig_i, N, J)
    if (has_lure) sigij <- sigij * exp(params$sigma_beta_lure * matrix(traps$lure, N, J, byrow = TRUE))
    d2 <- outer(s[, 1], traps$x, "-")^2 + outer(s[, 2], traps$y, "-")^2
    kern <- exp(-d2 / (2 * sigij^2))
    caught <- matrix(FALSE, N, J)
    for (k in seq_len(K)) {
      eta <- matrix(lp0, N, J, byrow = TRUE) + lsex
      if (has_rain) eta <- eta + params$detection_beta_rain * z[k]
      if (has_beh) eta <- eta + params$detection_beta_trap_response * caught
      U <- matrix(usage[, k], N, J, byrow = TRUE)
      if (pois) {
        yk <- matrix(stats::rpois(N * J, U * exp(eta) * kern), N, J)
      } else {
        pk <- stats::plogis(eta) * kern
        yk <- matrix(stats::rbinom(N * J, 1, 1 - (1 - pk)^U), N, J) * (U > 0)
      }
      caught <- caught | yk > 0
      w1 <- which(yk > 0, arr.ind = TRUE)
      if (nrow(w1)) {
        rp <- yk[w1]   # one row per genetic sample (Poisson) or detection
        caps[[k]] <- data.frame(ind = rep(w1[, 1], rp), occasion = k, trap = traps$id[rep(w1[, 2], rp)])
      }
    }
  }
  cap <- do.call(rbind, caps)
  if (is.null(cap)) cap <- data.frame(ind = integer(), occasion = integer(), trap = character())
  det <- sort(unique(cap$ind))
  cap$individual <- sprintf("ID%03d", match(cap$ind, det))
  ind_tab <- NULL
  if (has_sex && length(det)) {
    sx <- ifelse(sex[det] == 1, "M", "F")
    pu <- if (is.null(design$prop_sex_unknown)) 0 else design$prop_sex_unknown
    sx[stats::runif(length(det)) < pu] <- NA
    ind_tab <- data.frame(individual = sprintf("ID%03d", seq_along(det)), sex = sx)
  }
  cap <- cap[order(cap$individual, cap$occasion), c("individual", "occasion", "trap")]
  rownames(cap) <- NULL
  mask <- mask[hab, , drop = FALSE]; rownames(mask) <- NULL
  dcovs <- c(if (has_cov) "canopy", if (has_elev) "elevation")
  data <- list(traps = traps, usage = usage, captures = cap, n_occasions = K,
               mask = mask, mask_spacing = pix, buffer = buf, units = "m", mask_kind = "csv",
               density_covariates = dcovs,
               detection_covariates = if (has_trail) "on_trail" else character(),
               covariate_labels = list(canopy = "Canopy cover (%)", on_trail = "Camera on trail", elevation = "Elevation (m)",
                                       rain = "Rainfall (mm)", lure = "Scent lure at detector"))
  if (extras) {
    data$obs_model <- obs
    if (has_rain) { data$occasion_covariates <- data.frame(occasion = seq_len(K), rain = rain); data$occasion_detection_covariates <- "rain" }
    if (has_lure) data$sigma_covariates <- "lure"
    if (has_sex) { data$individuals <- ind_tab; data$sex_effects <- c("p0", "sigma") }
    if (has_beh) data$behaviour <- "trap"
  } else if (obs != "binomial") data$obs_model <- obs
  dn <- if (pois) "baseline_encounter_rate" else "baseline_detection"
  hr <- function(sg) pi * (sqrt(stats::qchisq(0.95, 2)) * sg / 1000)^2
  truth <- c(density_per_km2 = params$density, abundance = N)
  if (!has_sex) {
    truth[dn] <- det0; truth["sigma"] <- params$sigma; truth["home_range_area_95"] <- hr(params$sigma)
  } else {
    truth[paste0(dn, "_female")] <- det0; truth[paste0(dn, "_male")] <- dm
    truth["sigma_female"] <- params$sigma; truth["sigma_male"] <- params$sigma_male
    truth["home_range_area_95_female"] <- hr(params$sigma); truth["home_range_area_95_male"] <- hr(params$sigma_male)
    truth["proportion_female"] <- params$proportion_female
  }
  if (has_cov) truth["density_beta_canopy"] <- params$density_beta_canopy
  if (has_elev) truth["density_beta_elevation"] <- params$density_beta_elevation
  if (has_trail) truth["detection_beta_on_trail"] <- params$detection_beta_on_trail
  if (has_rain) truth["detection_beta_rain"] <- params$detection_beta_rain
  if (has_lure) truth["sigma_beta_lure"] <- params$sigma_beta_lure
  if (has_sex) { truth["detection_beta_male"] <- link(dm) - link(det0); truth["sigma_beta_male"] <- log(params$sigma_male / params$sigma) }
  if (has_beh) truth["detection_beta_trap_response"] <- params$detection_beta_trap_response
  list(data = data, truth = truth,
       latent = list(s = s, N = N, n = length(det), EN = EN, sex = sex))
}
