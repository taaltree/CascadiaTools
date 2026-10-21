# scr-open: simulate a multi-year spatial capture-recapture study of an open
# population (R port of ../../python/simulate_open.py, in map units of metres,
# with trap outages and named covariates added).
#
# Population: N_1 ~ Poisson(D * A) founders alive in period 1; between periods
# each animal alive at t-1 survives with probability phi; B_t ~ Poisson(f *
# N_{t-1}) recruits enter at t (per-capita recruitment f). Activity centres are
# static and drawn once per animal from pi_g ∝ exp(beta * X_g) over the mask
# pixels (uniform within pixel); X is a smooth Gaussian random field
# ("canopy cover"), standardised over the mask.
# Detection: proximity detectors; on each secondary occasion an animal alive
# in period t is detected at trap j with p = p0 exp(-d^2 / (2 sigma^2)) while
# the trap operates.

so_smooth_field <- function(nx, ny, sd_px) {
  pad <- ceiling(3 * sd_px) + 1
  z <- matrix(stats::rnorm((nx + 2 * pad) * (ny + 2 * pad)), nx + 2 * pad, ny + 2 * pad)
  if (sd_px > 0) {
    k <- stats::dnorm(seq(-pad, pad), 0, sd_px); k <- k / sum(k)
    z <- apply(z, 2, function(c) stats::filter(c, k, circular = TRUE))
    z <- t(apply(z, 1, function(r) stats::filter(r, k, circular = TRUE)))
  }
  z[pad + seq_len(nx), pad + seq_len(ny)]
}

so_simulate <- function(params, design, seed) {
  set.seed(seed)
  side <- design$trap_grid_side; sp <- design$trap_spacing
  K <- design$n_occasions; TT <- design$n_periods
  buf <- design$buffer; pix <- design$mask_spacing
  g1 <- buf + sp * (seq_len(side) - 1)
  tg <- expand.grid(x = g1, y = g1)
  J <- nrow(tg)
  traps <- data.frame(id = sprintf("T%02d", seq_len(J)), x = tg$x, y = tg$y)
  # trap outages: in each period a trap fails with prob trap_failure from a
  # random occasion onward
  usage <- matrix(as.integer(K), J, TT, dimnames = list(traps$id, NULL))
  occ_on <- array(1L, c(J, K, TT))
  tf <- if (is.null(design$trap_failure)) 0 else design$trap_failure
  for (t in seq_len(TT)) for (j in which(stats::runif(J) < tf)) {
    k0 <- sample.int(K, 1); occ_on[j, k0:K, t] <- 0L; usage[j, t] <- as.integer(k0 - 1)
  }
  xlim <- (side - 1) * sp + 2 * buf
  nxp <- round(xlim / pix)
  gx <- (seq_len(nxp) - 0.5) * pix
  mask <- expand.grid(x = gx, y = gx)
  G <- nrow(mask)
  f <- as.vector(so_smooth_field(nxp, nxp, design$covariate_range / pix))
  f <- (f - mean(f)) / stats::sd(f)
  mask$canopy <- 50 + 12 * f
  has_cov <- isTRUE(design$density_covariate)
  beta <- if (has_cov) params$density_beta_canopy else 0
  a_km2 <- G * pix^2 / 1e6
  EN1 <- params$density * a_km2
  # population process
  N1 <- stats::rpois(1, EN1)
  alive <- matrix(0L, N1, TT); alive[, 1] <- 1L
  for (t in 2:TT) {
    nprev <- sum(alive[, t - 1])
    surv <- alive[, t - 1] == 1
    alive[surv, t] <- stats::rbinom(sum(surv), 1, params$survival)
    B <- stats::rpois(1, params$per_capita_recruitment * nprev)
    if (B > 0) { nw <- matrix(0L, B, TT); nw[, t] <- 1L; alive <- rbind(alive, nw) }
  }
  Ns <- nrow(alive)
  entry <- apply(alive, 1, function(a) which(a == 1)[1])
  w <- exp(beta * f); w <- w / sum(w)
  cell <- sample.int(G, Ns, replace = TRUE, prob = w)
  s <- cbind(mask$x[cell] + stats::runif(Ns, -pix / 2, pix / 2), mask$y[cell] + stats::runif(Ns, -pix / 2, pix / 2))
  d2 <- outer(s[, 1], traps$x, "-")^2 + outer(s[, 2], traps$y, "-")^2
  p <- params$baseline_detection * exp(-d2 / (2 * params$sigma^2))
  caps <- list()
  for (t in seq_len(TT)) {
    al <- which(alive[, t] == 1)
    if (!length(al)) next
    for (k in seq_len(K)) {
      pk <- p[al, , drop = FALSE] * rep(occ_on[, k, t], each = length(al))
      yk <- matrix(stats::rbinom(length(pk), 1, pk), length(al), J)
      w1 <- which(yk == 1, arr.ind = TRUE)
      if (nrow(w1)) caps[[length(caps) + 1]] <- data.frame(ind = al[w1[, 1]], period = t, occasion = k, trap = traps$id[w1[, 2]])
    }
  }
  cap <- do.call(rbind, caps)
  if (is.null(cap)) cap <- data.frame(ind = integer(), period = integer(), occasion = integer(), trap = character())
  det <- sort(unique(cap$ind))
  cap$individual <- sprintf("ID%03d", match(cap$ind, det))
  cap <- cap[order(cap$individual, cap$period, cap$occasion), c("individual", "period", "occasion", "trap")]
  rownames(cap) <- NULL
  data <- list(traps = traps, usage = usage, captures = cap, n_occasions = K, n_periods = TT,
               period_labels = as.character(seq_len(TT)),
               mask = mask, mask_spacing = pix, buffer = buf, units = "m",
               density_covariates = if (has_cov) "canopy" else character(),
               covariate_labels = list(canopy = "Canopy cover (%)"))
  Nt <- colSums(alive)
  Bt <- vapply(seq_len(TT), function(t) sum(entry == t), 1)
  lab <- data$period_labels
  truth <- c(stats::setNames(Nt, paste0("abundance_period_", lab)),
             stats::setNames(Nt / a_km2, paste0("density_per_km2_period_", lab)),
             superpopulation_size = Ns,
             stats::setNames(Bt[-1], paste0("recruits_period_", lab[-1])),
             stats::setNames(Bt[-1] / pmax(Nt[-TT], 1), paste0("per_capita_recruitment_period_", lab[-1])),
             stats::setNames(Nt[-1] / pmax(Nt[-TT], 1), paste0("growth_rate_period_", lab[-1])),
             survival = params$survival,
             baseline_detection = params$baseline_detection, sigma = params$sigma,
             home_range_area_95 = pi * (sqrt(stats::qchisq(0.95, 2)) * params$sigma / 1000)^2)
  if (has_cov) truth["density_beta_canopy"] <- params$density_beta_canopy
  list(data = data, truth = truth,
       latent = list(s = s, alive = alive, entry = entry, cell = cell, n = length(det), det = det, EN1 = EN1, X = f))
}
