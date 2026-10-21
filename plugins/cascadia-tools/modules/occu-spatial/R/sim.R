# occu-spatial: data simulation written from the model equations (not from
# spOccupancy::simOcc, which is used only as an independent cross-check).
#
# Spatial random effect: w ~ MVN(0, sigma2 * exp(-phi * d)), exponential
# covariance, d = distance in km, phi = 3 / effective range (correlation 0.05
# at the effective range). Occupancy logit(psi_j) = b0 + x_j'b + w_j,
# z_j ~ Bern(psi_j), y_jk ~ Bern(z_j p_jk), logit(p_jk) = a0 + v_jk'a.
#
# All three data types return the SAME structure as read_data() (see io.R):
#   type, y, occ_X (data.frame of standardised model terms), det_X, coords (km),
#   site_ids, species/sources, terms + scaling, site_table, formulas.

# exponential-covariance Gaussian process at coords (km)
gp_draw <- function(coords, sigma_sq, range_km) {
  if (sigma_sq <= 0) return(rep(0, nrow(coords)))
  D <- as.matrix(stats::dist(coords))
  C <- sigma_sq * exp(-(3 / range_km) * D)
  L <- chol(C + diag(1e-8 * sigma_sq, nrow(C)))
  drop(crossprod(L, stats::rnorm(nrow(C))))
}

std_exact <- function(v, use = rep(TRUE, length(v))) (v - mean(v[use])) / stats::sd(v[use])

# site-level covariate: a mix of spatial structure (GP) and white noise,
# standardised exactly so the true slope is "per 1 SD" as in fit().
sim_site_cov <- function(coords, frac_spatial, range_km) {
  g <- if (frac_spatial > 0) gp_draw(coords, 1, range_km) else 0
  std_exact(sqrt(frac_spatial) * g + sqrt(1 - frac_spatial) * stats::rnorm(nrow(coords)))
}

sim_coords <- function(n, extent) {
  cbind(x = stats::runif(n, 0, extent), y = stats::runif(n, 0, extent))
}

# survey availability: first survey always done, later ones missing w.p. prop_missing
sim_avail <- function(J, K, prop_missing) {
  a <- matrix(stats::runif(J * K) > prop_missing, J, K)
  a[, 1] <- TRUE
  a
}

param_terms <- function(params, prefix) {
  nm <- grep(paste0("^", prefix), names(params), value = TRUE)
  stats::setNames(unlist(params[nm]), sub(paste0("^", prefix), "", nm))
}

d_or <- function(x, default) if (is.null(x)) default else x

base_data <- function(type, coords, occ_X, site_ids) {
  st <- data.frame(site = site_ids, x = coords[, 1], y = coords[, 2], occ_X, check.names = FALSE)
  terms <- names(occ_X)
  list(type = type, coords = coords, occ_X = occ_X, site_ids = site_ids,
       psi_terms = terms, psi_formula = if (length(terms)) paste("~", paste(terms, collapse = " + ")) else "~ 1",
       scaling = id_scaling(terms),
       term_source = stats::setNames(terms, terms), site_table = st, coord_units = "km",
       simulated = TRUE)
}

id_scaling <- function(terms) stats::setNames(lapply(terms, function(t) list(center = 0, scale = 1, label = t, binary = FALSE)), terms)

sim_det_covs <- function(names_, avail) {
  out <- lapply(names_, function(n) {
    m <- matrix(stats::rnorm(length(avail)), nrow(avail))
    m[avail] <- std_exact(m[avail]); m[!avail] <- NA; m
  })
  stats::setNames(out, names_)
}

lin_det <- function(a0, a, det_X, J, K) {
  eta <- matrix(a0, J, K)
  for (n in names(a)) eta <- eta + a[[n]] * det_X[[n]]
  stats::plogis(eta)
}

# ---- single species (spPGOcc / PGOcc) -------------------------------------
sim_single <- function(params, design) {
  J <- design$n_sites; K <- design$n_surveys
  coords <- sim_coords(J, design$extent_km)
  b <- param_terms(params, "psi_beta_"); a <- param_terms(params, "p_beta_")
  occ_X <- as.data.frame(lapply(stats::setNames(names(b), names(b)), function(n)
    sim_site_cov(coords, design$cov_spatial_frac, design$cov_range_km)))
  s2 <- d_or(params$spatial_variance, 0)
  w <- gp_draw(coords, s2, params$effective_range_km)
  eta <- stats::qlogis(params$occupancy_prob) + w
  for (n in names(b)) eta <- eta + b[[n]] * occ_X[[n]]
  psi <- stats::plogis(eta)
  z <- stats::rbinom(J, 1, psi)
  avail <- sim_avail(J, K, design$prop_missing)
  det_X <- sim_det_covs(names(a), avail)
  p <- lin_det(stats::qlogis(params$detection_prob), a, det_X, J, K)
  y <- matrix(stats::rbinom(J * K, 1, ifelse(avail, z * p, 0)), J, K); y[!avail] <- NA
  dat <- base_data("single", coords, occ_X, sprintf("S%03d", seq_len(J)))
  dat$y <- y; dat$det_X <- det_X; dat$p_terms <- names(a); dat$p_scaling <- id_scaling(names(a))
  dat$p_formula <- if (length(a)) paste("~", paste(names(a), collapse = " + ")) else "~ 1"
  truth <- c(occupancy_prob = params$occupancy_prob,
             stats::setNames(b, paste0("psi_beta_", names(b))),
             detection_prob = params$detection_prob,
             stats::setNames(a, paste0("p_beta_", names(a))),
             mean_occupancy = mean(psi), prop_sites_occupied = mean(z))
  if (s2 > 0) truth <- c(truth, spatial_variance = s2, effective_range_km = params$effective_range_km)
  dat$truth_latent <- list(psi = psi, z = z, w = w)
  list(data = dat, truth = truth)
}

# ---- multispecies (msPGOcc: n_factors = 0; sfMsPGOcc: n_factors >= 1) -----
sim_multi <- function(params, design) {
  J <- design$n_sites; K <- design$n_surveys; N <- design$n_species
  q <- d_or(design$n_factors, 0)
  coords <- sim_coords(J, design$extent_km)
  mb <- param_terms(params, "comm_psi_beta_"); sb <- param_terms(params, "comm_sd_psi_beta_")
  ma <- param_terms(params, "comm_p_beta_"); sa <- param_terms(params, "comm_sd_p_beta_")
  occ_X <- as.data.frame(lapply(stats::setNames(names(mb), names(mb)), function(n)
    sim_site_cov(coords, design$cov_spatial_frac, design$cov_range_km)))
  avail <- sim_avail(J, K, design$prop_missing)
  det_X <- sim_det_covs(names(ma), avail)
  sp <- sprintf("SP%02d", seq_len(N))
  b0 <- stats::rnorm(N, stats::qlogis(params$comm_occupancy_prob), params$comm_sd_occupancy_logit)
  a0 <- stats::rnorm(N, stats::qlogis(params$comm_detection_prob), params$comm_sd_detection_logit)
  B <- vapply(names(mb), function(n) stats::rnorm(N, mb[[n]], sb[[n]]), numeric(N))
  A <- vapply(names(ma), function(n) stats::rnorm(N, ma[[n]], sa[[n]]), numeric(N))
  B <- matrix(B, N); A <- matrix(A, N)
  W <- matrix(0, N, J)
  if (q > 0) {
    # factor model: w_i(s) = sum_r lambda_ir f_r(s), f_r ~ GP(1, phi_r)
    Lam <- matrix(stats::rnorm(N * q, 0, d_or(design$loading_sd, 1)), N, q)
    Lam[upper.tri(Lam)] <- 0; diag(Lam) <- 1
    rng <- rep_len(unlist(params$factor_range_km), q)
    Fm <- t(vapply(seq_len(q), function(r) gp_draw(coords, 1, rng[r]), numeric(J)))
    W <- Lam %*% matrix(Fm, q)
  }
  eta <- b0 + W
  for (k in seq_along(mb)) eta <- eta + outer(B[, k], occ_X[[k]])
  psi <- stats::plogis(eta)
  z <- matrix(stats::rbinom(N * J, 1, psi), N, J)
  y <- array(NA_integer_, c(N, J, K))
  for (i in seq_len(N)) {
    a_i <- stats::setNames(as.list(A[i, ]), names(ma))
    p <- lin_det(a0[i], a_i, det_X, J, K)
    yi <- matrix(stats::rbinom(J * K, 1, ifelse(avail, z[i, ] * p, 0)), J, K); yi[!avail] <- NA
    y[i, , ] <- yi
  }
  dimnames(y) <- list(sp, NULL, NULL)
  dat <- base_data("multi", coords, occ_X, sprintf("S%03d", seq_len(J)))
  dat$y <- y; dat$det_X <- det_X; dat$p_terms <- names(ma); dat$species <- sp; dat$p_scaling <- id_scaling(names(ma))
  dat$p_formula <- if (length(ma)) paste("~", paste(names(ma), collapse = " + ")) else "~ 1"
  truth <- c(comm_occupancy_prob = params$comm_occupancy_prob,
             comm_sd_occupancy_logit = params$comm_sd_occupancy_logit,
             stats::setNames(mb, paste0("comm_psi_beta_", names(mb))),
             stats::setNames(sb, paste0("comm_sd_psi_beta_", names(sb))),
             comm_detection_prob = params$comm_detection_prob,
             comm_sd_detection_logit = params$comm_sd_detection_logit,
             stats::setNames(ma, paste0("comm_p_beta_", names(ma))),
             stats::setNames(sa, paste0("comm_sd_p_beta_", names(sa))),
             mean_site_richness = mean(colSums(z)))
  if (q > 0) truth <- c(truth, stats::setNames(rep_len(unlist(params$factor_range_km), q),
                                               paste0("effective_range_km_f", seq_len(q))))
  dat$truth_latent <- list(psi = psi, z = z, w = W, richness = colSums(z))
  list(data = dat, truth = truth)
}

# ---- integrated: several data sources, shared occupancy (intPGOcc) --------
sim_integrated <- function(params, design) {
  src <- unlist(design$sources); R <- length(src)
  nj <- rep_len(unlist(design$n_sites_source), R); Kr <- rep_len(unlist(design$n_surveys_source), R)
  ov <- d_or(design$n_overlap, 0)
  # first source sites 1..nj[1]; source 2 overlaps the LAST `ov` sites of source 1, etc.
  sites <- list(); nxt <- 0
  for (r in seq_len(R)) {
    if (r == 1) { sites[[r]] <- seq_len(nj[1]); nxt <- nj[1] }
    else { o <- min(ov, nj[r]); sh <- if (o > 0) (nj[1] - o + 1):nj[1] else integer()
           new <- if (nj[r] - o > 0) (nxt + 1):(nxt + nj[r] - o) else integer()
           sites[[r]] <- c(sh, new); nxt <- nxt + length(new) }
  }
  J <- nxt
  coords <- sim_coords(J, design$extent_km)
  b <- param_terms(params, "psi_beta_")
  occ_X <- as.data.frame(lapply(stats::setNames(names(b), names(b)), function(n)
    sim_site_cov(coords, design$cov_spatial_frac, design$cov_range_km)))
  s2 <- d_or(params$spatial_variance, 0)
  w <- if (s2 > 0 && isTRUE(design$integrated_spatial)) gp_draw(coords, s2, params$effective_range_km) else rep(0, J)
  eta <- stats::qlogis(params$occupancy_prob) + w
  for (n in names(b)) eta <- eta + b[[n]] * occ_X[[n]]
  psi <- stats::plogis(eta); z <- stats::rbinom(J, 1, psi)
  y <- det_X <- p_terms <- list(); truth_det <- c()
  for (r in seq_len(R)) {
    s <- src[r]; Jr <- length(sites[[r]])
    allp <- param_terms(params, "p_beta_")
    a <- allp[grepl(paste0("_", s, "$"), names(allp))]
    names(a) <- sub(paste0("_", s, "$"), "", names(a))
    avail <- sim_avail(Jr, Kr[r], design$prop_missing)
    dx <- sim_det_covs(names(a), avail)
    p0 <- params[[paste0("detection_prob_", s)]]
    p <- lin_det(stats::qlogis(p0), a, dx, Jr, Kr[r])
    yr <- matrix(stats::rbinom(Jr * Kr[r], 1, ifelse(avail, z[sites[[r]]] * p, 0)), Jr, Kr[r]); yr[!avail] <- NA
    y[[s]] <- yr; det_X[[s]] <- dx; p_terms[[s]] <- names(a)
    truth_det <- c(truth_det, stats::setNames(p0, paste0("detection_prob_", s)),
                   stats::setNames(a, paste0("p_beta_", names(a), "_", s)))
  }
  dat <- base_data("integrated", coords, occ_X, sprintf("S%03d", seq_len(J)))
  dat$y <- y; dat$det_X <- det_X; dat$p_terms <- p_terms; dat$sources <- src
  dat$p_scaling <- lapply(p_terms, id_scaling)
  dat$sites <- stats::setNames(sites, src)
  dat$p_formula <- lapply(p_terms, function(t) if (length(t)) paste("~", paste(t, collapse = " + ")) else "~ 1")
  truth <- c(occupancy_prob = params$occupancy_prob,
             stats::setNames(b, paste0("psi_beta_", names(b))), truth_det,
             mean_occupancy = mean(psi), prop_sites_occupied = mean(z))
  dat$truth_latent <- list(psi = psi, z = z, w = w)
  list(data = dat, truth = truth)
}

simulate_occu_spatial <- function(params, design, seed) {
  set.seed(seed)
  type <- d_or(design$type, "single")
  switch(type,
         single = sim_single(params, design),
         multi = sim_multi(params, design),
         integrated = sim_integrated(params, design),
         stop("design$type must be single, multi or integrated"))
}
