# scr-search: simulate a search-encounter (detection-dog / snow-tracking)
# noninvasive genetic survey.
#
# Landscape: square searched region (search_side) + buffer = state space,
# a few winding roads crossing it, and a smooth habitat covariate ("forest
# cover") that can be correlated with proximity to roads.
# Search tracks: per occasion, tracks_per_occasion search bouts of 4-10 km.
# A share road_bias of them follow a road (searchers walk / drive roads);
# the rest are correlated random walks off-road. GPS fixes every ~50 m.
# Each track has a team covariate (experienced dog team 0/1).
# Animals: Poisson number of activity centres, log-linear in forest cover.
# Samples: deposited and found CONTINUOUSLY along the walked path:
# the number of samples of animal i on a short piece of track of length l
# at distance d from its centre is Poisson(lambda0 * l * exp(-d^2/2sigma^2) *
# exp(beta_team * team)); each sample's GPS location is the piece centre
# (pieces are <= 50 m) plus 10 m jitter. The fitted model only sees tracks and samples,
# so the cell discretisation is tested, not assumed.

ss_smooth_field <- function(nx, ny, sd_px) {
  pad <- ceiling(3 * sd_px) + 1
  z <- matrix(stats::rnorm((nx + 2 * pad) * (ny + 2 * pad)), nx + 2 * pad, ny + 2 * pad)
  if (sd_px > 0) {
    k <- stats::dnorm(seq(-pad, pad), 0, sd_px); k <- k / sum(k)
    z <- apply(z, 2, function(c) stats::filter(c, k, circular = TRUE))
    z <- t(apply(z, 1, function(r) stats::filter(r, k, circular = TRUE)))
  }
  z[pad + seq_len(nx), pad + seq_len(ny)]
}

# A winding road crossing the searched region [lo, hi]^2, points every `step`
ss_road <- function(lo, hi, step) {
  cen <- c(stats::runif(1, lo + 0.3 * (hi - lo), hi - 0.3 * (hi - lo)), stats::runif(1, lo + 0.3 * (hi - lo), hi - 0.3 * (hi - lo)))
  th <- stats::runif(1, 0, pi)
  L <- 1.5 * (hi - lo)
  t <- seq(-L, L, by = step)
  amp <- 0.06 * (hi - lo); wl <- stats::runif(1, 0.3, 0.6) * (hi - lo); ph <- stats::runif(1, 0, 2 * pi)
  off <- amp * sin(2 * pi * t / wl + ph)
  x <- cen[1] + t * cos(th) - off * sin(th); y <- cen[2] + t * sin(th) + off * cos(th)
  ok <- x >= lo & x <= hi & y >= lo & y <= hi
  data.frame(x = x[ok], y = y[ok])
}

ss_simulate <- function(params, design, seed) {
  set.seed(seed)
  S <- design$search_side; buf <- design$buffer; pix <- design$mask_spacing
  lo <- buf; hi <- buf + S
  side <- S + 2 * buf
  nxp <- round(side / pix)
  gx <- (seq_len(nxp) - 0.5) * pix
  mask <- expand.grid(x = gx, y = gx)
  G <- nrow(mask)
  step <- 50
  # ---- roads
  roads <- lapply(seq_len(design$n_roads), function(k) { r <- ss_road(lo, hi, step); while (nrow(r) < 40) r <- ss_road(lo, hi, step); r })
  rp <- do.call(rbind, roads)
  droad <- ss_min_dist(mask$x, mask$y, rp$x, rp$y)
  # ---- habitat covariate, optionally correlated with closeness to roads
  f <- as.vector(ss_smooth_field(nxp, nxp, design$covariate_range / pix)); f <- (f - mean(f)) / stats::sd(f)
  rho <- if (is.null(design$cov_road_corr)) 0 else design$cov_road_corr
  if (rho != 0) {
    cl <- -log(droad + 500); cl <- (cl - mean(cl)) / stats::sd(cl)
    f <- sqrt(1 - rho^2) * f + rho * cl; f <- (f - mean(f)) / stats::sd(f)
  }
  mask$forest <- round(55 + 12 * f, 1)
  fz <- (mask$forest - mean(mask$forest)) / stats::sd(mask$forest)
  has_cov <- isTRUE(design$density_covariate)
  beta <- if (has_cov) params$density_beta_forest else 0
  # ---- animals
  a_km2 <- G * pix^2 / 1e6
  N <- stats::rpois(1, params$density * a_km2)
  w <- exp(beta * fz); w <- w / sum(w)
  cell <- sample.int(G, N, replace = TRUE, prob = w)
  s <- cbind(mask$x[cell] + stats::runif(N, -pix / 2, pix / 2), mask$y[cell] + stats::runif(N, -pix / 2, pix / 2))
  # ---- search tracks
  K <- design$n_occasions; nt <- design$tracks_per_occasion
  trk <- list(); tinfo <- list(); id <- 0
  for (k in seq_len(K)) for (t in seq_len(nt)) {
    id <- id + 1
    L <- stats::runif(1, design$track_km[1], design$track_km[2]) * 1000
    nstep <- max(2, round(L / step))
    if (stats::runif(1) < design$road_bias) {
      r <- roads[[sample.int(length(roads), 1)]]
      i0 <- sample.int(nrow(r), 1); dir <- sample(c(-1, 1), 1)
      ii <- i0 + dir * (0:(nstep - 1)); ii <- ii[ii >= 1 & ii <= nrow(r)]
      if (length(ii) < nstep / 2) { ii <- i0 - dir * (0:(nstep - 1)); ii <- ii[ii >= 1 & ii <= nrow(r)] }
      xy <- cbind(r$x[ii], r$y[ii]) + matrix(stats::rnorm(2 * length(ii), 0, 5), ncol = 2)
    } else {
      xy <- matrix(NA_real_, nstep, 2)
      xy[1, ] <- c(stats::runif(1, lo, hi), stats::runif(1, lo, hi)); h <- stats::runif(1, 0, 2 * pi)
      for (q in 2:nstep) {
        h <- h + stats::rnorm(1, 0, 0.35)
        nx <- xy[q - 1, 1] + step * cos(h); ny <- xy[q - 1, 2] + step * sin(h)
        if (nx < lo || nx > hi) { h <- pi - h; nx <- xy[q - 1, 1] + step * cos(h) }
        if (ny < lo || ny > hi) { h <- -h; ny <- xy[q - 1, 2] + step * sin(h) }
        xy[q, ] <- c(min(max(nx, lo), hi), min(max(ny, lo), hi))
      }
    }
    tid <- sprintf("T%03d", id)
    trk[[id]] <- data.frame(track_id = tid, occasion = k, x = round(xy[, 1], 1), y = round(xy[, 2], 1))
    tinfo[[id]] <- data.frame(track_id = tid, occasion = k, experienced_team = stats::rbinom(1, 1, 0.5))
  }
  tracks <- do.call(rbind, trk); tinfo <- do.call(rbind, tinfo)
  has_team <- isTRUE(design$detection_covariate)
  tracks$experienced_team <- tinfo$experienced_team[match(tracks$track_id, tinfo$track_id)]
  # ---- samples deposited along the walked path
  pcs <- ss_pieces(tracks[, c("track_id", "x", "y")], step)
  pcs$team <- tinfo$experienced_team[match(pcs$track_id, tinfo$track_id)]
  pcs$occasion <- tinfo$occasion[match(pcs$track_id, tinfo$track_id)]
  mult <- params$lambda0_per_km * pcs$len / 1000 * if (has_team) exp(params$detection_beta_experienced_team * pcs$team) else 1
  sg2 <- 2 * params$sigma^2
  out <- list()
  for (i in seq_len(N)) {
    lam <- mult * exp(-((pcs$x - s[i, 1])^2 + (pcs$y - s[i, 2])^2) / sg2)
    k <- stats::rpois(length(lam), lam)
    if (any(k > 0)) {
      ix <- rep(which(k > 0), k[k > 0])
      out[[length(out) + 1]] <- data.frame(ind = i, piece = ix)
    }
  }
  smp <- do.call(rbind, out)
  if (is.null(smp)) smp <- data.frame(ind = integer(), piece = integer())
  # position: centre of the (<= 50 m) piece + 10 m GPS / off-trail jitter
  if (nrow(smp)) {
    pp <- pcs[smp$piece, ]
    sx <- pp$x + stats::rnorm(nrow(smp), 0, 10)
    sy <- pp$y + stats::rnorm(nrow(smp), 0, 10)
    det <- sort(unique(smp$ind))
    samples <- data.frame(sample_id = sprintf("S%04d", seq_len(nrow(smp))), individual = sprintf("ID%03d", match(smp$ind, det)),
                          x = round(sx, 1), y = round(sy, 1), track_id = pp$track_id, occasion = pp$occasion)
    samples <- samples[order(samples$occasion, samples$track_id), ]
    samples$sample_id <- sprintf("S%04d", seq_len(nrow(samples))); rownames(samples) <- NULL
  } else samples <- data.frame(sample_id = character(), individual = character(), x = numeric(), y = numeric(),
                               track_id = character(), occasion = integer())
  data <- list(tracks = if (has_team) tracks else tracks[, c("track_id", "occasion", "x", "y")],
               samples = samples, n_occasions = K, mask = mask[, c("x", "y", "forest")], mask_spacing = pix, buffer = buf,
               cell_size = design$cell_size, units = "m",
               density_covariates = if (has_cov) "forest" else character(),
               detection_covariates = if (has_team) "experienced_team" else character(),
               covariate_labels = list(forest = "Forest cover (%)", experienced_team = "Experienced dog team"),
               roads = rp)
  truth <- c(density_per_km2 = params$density, abundance = N, lambda0_per_km = params$lambda0_per_km,
             sigma = params$sigma, home_range_area_95 = pi * (sqrt(stats::qchisq(0.95, 2)) * params$sigma / 1000)^2,
             samples_per_track_through_centre = params$lambda0_per_km * sqrt(2 * pi) * params$sigma / 1000)
  if (has_cov) truth["density_beta_forest"] <- params$density_beta_forest
  if (has_team) truth["detection_beta_experienced_team"] <- params$detection_beta_experienced_team
  list(data = data, truth = truth, latent = list(s = s, N = N, n = length(unique(smp$ind))))
}
