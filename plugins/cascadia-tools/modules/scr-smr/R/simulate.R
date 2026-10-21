# scr-smr: simulate a camera-trap spatial mark-resight study.
#
# Population: N ~ Poisson(D * A) activity centres on a rectangular state
# space (camera grid + buffer), log-linear in a smooth "canopy" field
# (as modules/scr). Half-normal space use with scale sigma.
#
# Marking:
#   smr  : n_marked animals chosen at random from the population (e.g. collared
#          earlier, marked by an unmodelled process unrelated to space).
#   gsmr : a marking grid (live traps) of marking_grid_side^2 traps centred on
#          (camera-grid centre + marking_offset), marking_occasions occasions,
#          capture prob pm0 exp(-d^2 / 2 sigma^2) per trap and occasion; every
#          animal caught at least once is marked. Marked animals are therefore
#          concentrated near the marking grid (Whittington et al. 2018).
# Resighting (cameras, K occasions): sightings of animal i at camera j on
# occasion k ~ Poisson(lambda0_j exp(-d^2 / 2 sigma^2)) (resight_model
# "poisson") or Bernoulli with that probability ("bernoulli"); log(lambda0_j)
# [logit for Bernoulli] = log(lambda0) + beta_trail * on_trail_j. A sighting of
# a marked animal is identified with probability id_probability, otherwise it
# is recorded as "marked_unidentified". Optional telemetry: telemetry_fixes
# bivariate-normal fixes N(s_i, sigma^2 I) for a share of the marked animals.

smr_smooth_field <- function(nx, ny, sd_px) {
  pad <- ceiling(3 * sd_px) + 1
  z <- matrix(stats::rnorm((nx + 2 * pad) * (ny + 2 * pad)), nx + 2 * pad, ny + 2 * pad)
  if (sd_px > 0) {
    k <- stats::dnorm(seq(-pad, pad), 0, sd_px); k <- k / sum(k)
    z <- apply(z, 2, function(c) stats::filter(c, k, circular = TRUE))
    z <- t(apply(z, 1, function(r) stats::filter(r, k, circular = TRUE)))
  }
  z[pad + seq_len(nx), pad + seq_len(ny)]
}

smr_simulate <- function(params, design, seed) {
  set.seed(seed)
  side <- design$trap_grid_side; sp <- design$trap_spacing; K <- design$n_occasions
  buf <- design$buffer; pix <- design$mask_spacing
  model <- if (is.null(design$model)) "smr" else design$model
  pois <- !identical(design$resight_model, "bernoulli")
  g1 <- buf + sp * (seq_len(side) - 1)
  cg <- expand.grid(x = g1, y = g1)
  J <- nrow(cg)
  cams <- data.frame(id = sprintf("CAM%02d", seq_len(J)), x = cg$x, y = cg$y)
  has_trail <- isTRUE(design$detection_covariate)
  if (has_trail) cams$on_trail <- as.integer(seq_len(J) %in% sample.int(J, round(design$prop_on_trail * J)))
  usage <- matrix(1L, J, K, dimnames = list(cams$id, NULL))
  fail <- which(stats::runif(J) < design$trap_failure)
  for (j in fail) { k0 <- sample.int(K, 1); if (k0 > 1) usage[j, k0:K] <- 0L }
  # state space and habitat covariate
  xlim <- (side - 1) * sp + 2 * buf
  nxp <- round(xlim / pix)
  gx <- (seq_len(nxp) - 0.5) * pix
  mask <- expand.grid(x = gx, y = gx)
  G <- nrow(mask)
  f <- as.vector(smr_smooth_field(nxp, nxp, design$covariate_range / pix))
  f <- (f - mean(f)) / stats::sd(f)
  mask$canopy <- 50 + 12 * f
  has_cov <- isTRUE(design$density_covariate)
  beta <- if (has_cov) params$density_beta_canopy else 0
  A_km2 <- G * pix^2 / 1e6
  w <- exp(beta * f); w <- w / sum(w)
  N <- stats::rpois(1, params$density * A_km2)
  cell <- sample.int(G, N, replace = TRUE, prob = w)
  s <- cbind(mask$x[cell] + stats::runif(N, -pix / 2, pix / 2), mask$y[cell] + stats::runif(N, -pix / 2, pix / 2))
  sig <- params$sigma
  # ---- marking
  marking <- NULL; marked <- logical(N)
  if (model == "gsmr") {
    ms <- design$marking_grid_side; msp <- design$marking_spacing
    cen <- mean(g1) + if (is.null(design$marking_offset)) c(0, 0) else unlist(design$marking_offset)
    off <- msp * (seq_len(ms) - 1) - msp * (ms - 1) / 2
    mg <- expand.grid(x = cen[1] + off, y = cen[2] + off)
    L <- nrow(mg); Km <- design$marking_occasions
    mtr <- data.frame(id = sprintf("LT%02d", seq_len(L)), x = mg$x, y = mg$y)
    d2m <- outer(s[, 1], mtr$x, "-")^2 + outer(s[, 2], mtr$y, "-")^2
    pm <- params$marking_detection * exp(-d2m / (2 * sig^2))
    caps <- list()
    for (k in seq_len(Km)) {
      yk <- matrix(stats::rbinom(N * L, 1, pm), N, L)
      w1 <- which(yk == 1, arr.ind = TRUE)
      if (nrow(w1)) caps[[k]] <- data.frame(ind = w1[, 1], occasion = k, trap = mtr$id[w1[, 2]])
    }
    cap <- do.call(rbind, caps)
    if (is.null(cap)) cap <- data.frame(ind = integer(), occasion = integer(), trap = character())
    marked[unique(cap$ind)] <- TRUE
  } else {
    nm <- min(design$n_marked, N)
    marked[sample.int(N, nm)] <- TRUE
  }
  mid <- which(marked)
  lab <- rep(NA_character_, N); lab[mid] <- sprintf("M%03d", seq_along(mid))
  if (model == "gsmr") {
    cap$individual <- lab[cap$ind]
    cap <- cap[order(cap$individual, cap$occasion), c("individual", "occasion", "trap")]; rownames(cap) <- NULL
    marking <- list(traps = mtr, usage = matrix(1L, nrow(mtr), Km, dimnames = list(mtr$id, NULL)),
                    captures = cap, n_occasions = Km)
  }
  # ---- resighting
  lr0 <- if (pois) log(params$resight_rate) else stats::qlogis(params$resight_rate)
  lr0 <- lr0 + if (has_trail) params$detection_beta_on_trail * cams$on_trail else 0
  r0 <- if (pois) exp(lr0) else stats::plogis(lr0)
  d2 <- outer(s[, 1], cams$x, "-")^2 + outer(s[, 2], cams$y, "-")^2
  lam <- sweep(exp(-d2 / (2 * sig^2)), 2, r0, "*")
  delta <- if (is.null(params$id_probability)) 1 else params$id_probability
  rows <- list()
  for (k in seq_len(K)) {
    lk <- lam * rep(usage[, k], each = N)
    yk <- if (pois) matrix(stats::rpois(N * J, lk), N, J) else matrix(stats::rbinom(N * J, 1, lk), N, J)
    w1 <- which(yk > 0, arr.ind = TRUE)
    if (!nrow(w1)) next
    cnts <- yk[w1]
    ind <- rep(w1[, 1], cnts); jj <- rep(w1[, 2], cnts)
    st <- ifelse(marked[ind], ifelse(stats::runif(length(ind)) < delta, "marked", "marked_unidentified"), "unmarked")
    rows[[length(rows) + 1]] <- data.frame(camera = cams$id[jj], occasion = k, status = st,
                                           animal_id = ifelse(st == "marked", lab[ind], NA_character_))
  }
  si <- do.call(rbind, rows)
  if (is.null(si)) si <- data.frame(camera = character(), occasion = integer(), status = character(), animal_id = character())
  si <- si[order(si$occasion, si$camera), ]; rownames(si) <- NULL
  # ---- telemetry
  tel <- NULL
  nfix <- if (is.null(design$telemetry_fixes)) 0 else design$telemetry_fixes
  if (nfix > 0 && length(mid)) {
    fr <- if (is.null(design$telemetry_fraction)) 1 else design$telemetry_fraction
    tid <- mid[seq_len(max(1, round(fr * length(mid))))]
    tel <- do.call(rbind, lapply(tid, function(i) data.frame(animal_id = lab[i], x = stats::rnorm(nfix, s[i, 1], sig),
                                                             y = stats::rnorm(nfix, s[i, 2], sig))))
  }
  data <- list(cameras = cams, usage = usage, n_occasions = K, sightings = si,
               marked = if (model == "smr") data.frame(id = lab[mid]) else NULL,
               unidentified_recorded = delta < 1, telemetry = tel, marking = marking,
               mask = mask, mask_spacing = pix, buffer = buf, units = "m",
               density_covariates = if (has_cov) "canopy" else character(),
               detection_covariates = if (has_trail) "on_trail" else character(),
               covariate_labels = list(canopy = "Canopy cover (%)", on_trail = "Camera on trail"))
  m <- length(mid)
  truth <- c(density_per_km2 = N / A_km2, abundance = N, marked_fraction = if (N > 0) m / N else NA,
             sigma = sig, home_range_area_95 = pi * (sqrt(stats::qchisq(0.95, 2)) * sig / 1000)^2)
  truth[if (pois) "resight_rate" else "resight_detection"] <- params$resight_rate
  if (delta < 1) truth["id_probability"] <- delta
  if (model == "gsmr") truth["marking_detection"] <- params$marking_detection
  if (has_cov) truth["density_beta_canopy"] <- params$density_beta_canopy
  if (has_trail) truth["detection_beta_on_trail"] <- params$detection_beta_on_trail
  list(data = data, truth = truth, latent = list(s = s, N = N, marked = marked, m = m))
}
