# Simulation: animals (or groups) placed by an inhomogeneous Poisson process
# whose intensity varies between survey cells with a habitat covariate, then
# detected from line or point transects with a known detection function.
#
# Study area = n_transects cells; each cell holds one transect (a line of
# transect_length_km down its middle in a cell transect_length_km x
# cell_width_km, or a point in a cell of point_cell_km2). Groups in cell i:
#   N_i ~ Poisson(D_i * a_cell),  log D_i = log D0 + beta * habitat_i
# Groups inside the covered strip/disc of half-width W_sim = 1.5 x truncation
# get uniform positions (perpendicular distance x ~ U(0, W_sim) for lines,
# radial r with density 2r / W_sim^2 for points), so distances beyond the
# truncation distance are recorded as they would be in the field.
# Optional: group sizes (1 + Poisson(mean - 1)), observer covariate on log
# sigma, repeated visits with availability (temporary emigration; Chandler et
# al. 2011), and violations used as negative controls: responsive movement
# away from the observer before detection, heaping of small distances at 0,
# rounding. With repeated visits, positions are redrawn on every visit unless
# design$fixed_positions = TRUE.

dist_simulate <- function(params, design, seed) {
  set.seed(seed)
  p <- params; d <- design
  survey <- d$survey
  w <- d$truncation_m
  Wsim <- 1.5 * w
  nT <- d$n_transects
  nV <- max(1, d$n_visits)
  key <- p$detection_key
  mu_s <- p$mean_group_size
  clusters <- mu_s > 1
  hab <- as.numeric(scale(stats::rnorm(nT)))          # standardised: mean 0, SD 1
  D <- p$group_density_at_average_site * exp(p$density_beta_habitat * hab)  # groups / km2
  a_cell <- if (survey == "line") d$transect_length_km * d$cell_width_km else d$point_cell_km2
  a_cov_sim <- if (survey == "line") d$transect_length_km * 2 * Wsim / 1000 else pi * (Wsim / 1000)^2
  if (a_cov_sim > a_cell) stop("covered strip wider than the survey cell; increase cell_width_km / point_cell_km2")
  nstrat <- max(1, d$n_strata)
  stratum <- if (nstrat > 1) paste0("S", (seq_len(nT) - 1) %% nstrat + 1) else rep("all", nT)
  observer <- if (isTRUE(d$observers)) rep(c("A", "B"), length.out = nT) else rep("A", nT)
  sig_t <- p$sigma_m * exp(p$detection_beta_observer * (observer == "B"))
  phi <- if (nV > 1) p$availability_prob else 1

  Ncell <- stats::rpois(nT, D * a_cell)
  Ncov <- stats::rbinom(nT, Ncell, a_cov_sim / a_cell)
  ind_total <- if (clusters) Ncell + stats::rpois(nT, Ncell * (mu_s - 1)) else Ncell
  obs <- list()
  pa_cov <- numeric()   # realised P_a of every group inside truncation (for the truth)
  for (i in seq_len(nT)) {
    if (!Ncov[i]) next
    x <- if (survey == "line") stats::runif(Ncov[i], 0, Wsim) else Wsim * sqrt(stats::runif(Ncov[i]))
    sz <- if (clusters) 1 + stats::rpois(Ncov[i], mu_s - 1) else rep(1, Ncov[i])
    inw <- x <= w
    pa_cov <- c(pa_cov, rep(ds_pa(key, sig_t[i], w, survey, p$hr_shape), sum(inw)))
    for (v in seq_len(nV)) {
      avail <- stats::runif(Ncov[i]) < phi
      # repeated visits: animals move around their home range between visits, so
      # positions are redrawn each visit (the gdistsamp assumption), each animal
      # within the band where it lives (inside or outside the truncation distance),
      # so the superpopulation of the covered area stays the same animals.
      # fixed_positions = TRUE keeps them (sedentary animals; a known violation).
      if (v > 1 && !isTRUE(d$fixed_positions)) {
        u <- stats::runif(Ncov[i])
        x <- if (survey == "line") ifelse(inw, u * w, w + u * (Wsim - w))
             else ifelse(inw, w * sqrt(u), sqrt(w^2 + u * (Wsim^2 - w^2)))
      }
      xo <- x
      if (d$responsive_move_m > 0) xo <- x + d$responsive_move_m * exp(-x / d$responsive_move_m)
      det <- avail & stats::runif(Ncov[i]) < ds_g(xo, key, sig_t[i], p$hr_shape)
      if (!any(det)) next
      xr <- xo[det]
      if (d$heap_zero_m > 0) xr[xr < d$heap_zero_m] <- 0
      if (d$round_m > 0) xr <- round(xr / d$round_m) * d$round_m
      obs[[length(obs) + 1]] <- data.frame(transect = sprintf("T%02d", i), visit = v,
                                           distance = xr, size = sz[det])
    }
  }
  obs <- if (length(obs)) do.call(rbind, obs) else data.frame(transect = character(), visit = integer(), distance = numeric(), size = numeric())
  tr <- data.frame(transect = sprintf("T%02d", seq_len(nT)), stratum = stratum,
                   effort = if (survey == "line") d$transect_length_km else 1,
                   habitat = hab, observer = observer)
  A_str <- tapply(rep(a_cell, nT), stratum, sum)
  strata <- data.frame(stratum = names(A_str), area_km2 = as.numeric(A_str))
  cut <- NULL
  if (d$distance_type == "binned") {
    cut <- seq(0, w, length.out = d$n_bins + 1)
    obs <- obs[obs$distance <= w, ]
    j <- findInterval(obs$distance, cut, rightmost.closed = TRUE, all.inside = TRUE)
    obs$distbegin <- cut[j]; obs$distend <- cut[j + 1]
    obs$distance <- (obs$distbegin + obs$distend) / 2
  }
  if (!clusters) obs$size <- NULL
  data <- list(
    survey = survey, distance_type = d$distance_type, truncation_m = w, cutpoints = cut,
    obs = obs, transects = tr, strata = strata, n_visits = nV, clusters = clusters,
    formulas = list(density = "~ habitat",
                    detection = if (isTRUE(d$observers)) "~ observer" else "~ 1"),
    state_covs = "habitat", det_covs = if (isTRUE(d$observers)) "observer" else character(),
    scaling = list(habitat = list(center = 0, scale = 1)),
    species = "simulated species", group_word = if (clusters) "group" else "animal",
    grid = NULL)

  A <- nT * a_cell
  pa_true <- if (length(pa_cov)) mean(pa_cov) else ds_pa(key, p$sigma_m, w, survey, p$hr_shape)
  truth <- c(density_per_km2 = sum(ind_total) / A,
             abundance = sum(ind_total),
             density_at_average_site = p$group_density_at_average_site * mu_s,
             density_beta_habitat = p$density_beta_habitat,
             rate_ratio_habitat = exp(p$density_beta_habitat),
             detection_prob = pa_true)
  truth[[if (survey == "line") "esw_m" else "edr_m"]] <- ds_esw(pa_true, w, survey)
  if (key == "hn") truth[["sigma_m"]] <- p$sigma_m
  if (key == "hr") { truth[["hr_scale_m"]] <- p$sigma_m; truth[["hr_shape"]] <- p$hr_shape }
  if (clusters) {
    truth[["group_density_per_km2"]] <- sum(Ncell) / A
    truth[["group_density_at_average_site"]] <- p$group_density_at_average_site
    truth[["mean_group_size"]] <- mu_s
  }
  if (isTRUE(d$observers)) truth[["detection_beta_observerB"]] <- p$detection_beta_observer
  if (nV > 1) truth[["availability_prob"]] <- phi
  if (nstrat > 1) for (s in unique(stratum))
    truth[[paste0("density_per_km2_", s)]] <- sum(ind_total[stratum == s]) / A_str[[s]]
  list(data = data, truth = truth)
}
