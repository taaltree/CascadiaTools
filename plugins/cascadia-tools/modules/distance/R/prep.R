# Shared data preparation for all engines.
#
# data (from dist_simulate() or dist_read_data()):
#   survey "line"/"point"; distance_type "exact"/"binned"; truncation_m;
#   cutpoints (m, binned only); obs data.frame(transect, distance [m],
#   [distbegin, distend], [size], [visit], [obs-level detection covariates]);
#   transects data.frame(transect, effort [km of line per visit, or 1 per
#   point visit], stratum, covariates...); strata data.frame(stratum, area_km2)
#   or NULL; n_visits; clusters; formulas list(density, detection);
#   state_covs, det_covs; scaling; grid (optional prediction grid).

`%||%` <- function(a, b) if (is.null(a)) b else a

ds_rhs <- function(f) { f <- trimws(f %||% "~ 1"); if (!startsWith(f, "~")) f <- paste("~", f); f }
ds_terms <- function(f) all.vars(stats::as.formula(ds_rhs(f)))

# Standardise numeric site covariates (mean 0, SD 1), keep factors; record
# the scaling so figures can draw original units. Returns the prepared list.
dist_prep <- function(data) {
  w <- data$truncation_m
  o <- data$obs
  if (data$distance_type == "binned" && !all(c("distbegin", "distend") %in% names(o)))
    stop("binned data need distbegin and distend columns")
  o <- if (data$distance_type == "binned") o[o$distend <= w + 1e-9, , drop = FALSE] else o[o$distance <= w, , drop = FALSE]
  if (!"visit" %in% names(o)) o$visit <- 1
  if (isTRUE(data$clusters) && !"size" %in% names(o)) stop("clusters = TRUE but no size column")
  tr <- data$transects
  if (!"stratum" %in% names(tr)) tr$stratum <- "all"
  tr$stratum <- as.character(tr$stratum)
  st_terms <- ds_terms(data$formulas$density)
  dt_terms <- ds_terms(data$formulas$detection)
  scaling <- data$scaling %||% list()
  zs <- tr
  for (v in unique(c(st_terms, dt_terms))) {
    if (!v %in% names(tr)) next
    x <- tr[[v]]
    if (is.character(x) || is.logical(x)) { zs[[v]] <- factor(x); next }
    if (is.factor(x)) next
    sc <- scaling[[v]]
    if (is.null(sc)) { sc <- list(center = mean(x), scale = if (stats::sd(x) > 0) stats::sd(x) else 1); scaling[[v]] <- sc }
    zs[[v]] <- (x - sc$center) / sc$scale
  }
  # observation-level detection covariates (MCDS): copy transect covariates
  for (v in dt_terms) {
    if (!v %in% names(o)) {
      o[[v]] <- zs[[v]][match(o$transect, tr$transect)]
    } else if (is.numeric(o[[v]])) {
      if (is.null(scaling[[v]])) scaling[[v]] <- list(center = mean(o[[v]]), scale = max(stats::sd(o[[v]]), 1e-8))
      o[[v]] <- (o[[v]] - scaling[[v]]$center) / scaling[[v]]$scale
    } else o[[v]] <- factor(o[[v]])
  }
  if (!is.null(data$grid)) {
    gz <- data$grid
    for (v in st_terms) {
      if (!v %in% names(gz)) stop("prediction grid lacks the density covariate '", v, "'")
      if (is.factor(zs[[v]])) gz[[v]] <- factor(gz[[v]], levels = levels(zs[[v]]))
      else gz[[v]] <- (gz[[v]] - scaling[[v]]$center) / scaling[[v]]$scale
    }
    data$grid$z <- gz
  }
  data$obs_w <- o
  data$tz <- zs
  data$scaling <- scaling
  data$st_terms <- st_terms; data$dt_terms <- dt_terms
  data$dt_site_level <- all(dt_terms %in% names(tr))
  data
}

# Covered area per transect (km2) for all visits.
ds_cov_area <- function(data, w = data$truncation_m) {
  if (data$survey == "line") 2 * (w / 1000) * data$transects$effort * data$n_visits
  else pi * (w / 1000)^2 * data$transects$effort * data$n_visits
}

ds_total_area <- function(data) if (is.null(data$strata)) NA_real_ else sum(data$strata$area_km2)

# Wald / log-normal helpers -------------------------------------------------
ds_row <- function(name, est, se, lo, hi) data.frame(param = name, estimate = est, sd = se, lower = lo, upper = hi)
ds_lognorm_row <- function(name, est, se) {
  cv <- se / est; C <- exp(1.96 * sqrt(log(1 + cv^2)))
  ds_row(name, est, se, est / C, est * C)
}
ds_logit_row <- function(name, est, se) {
  l <- stats::qlogis(est); sl <- se / (est * (1 - est))
  ds_row(name, est, se, stats::plogis(l - 1.96 * sl), stats::plogis(l + 1.96 * sl))
}
ds_draw_rows <- function(dr, names) {
  dr <- as.matrix(dr); colnames(dr) <- names
  derived_rows(dr)
}
ds_clean <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

ds_mvn <- function(theta, V, n) {
  ok <- all(is.finite(V)) && all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > -1e-10)
  if (!ok) return(matrix(theta, n, length(theta), byrow = TRUE, dimnames = list(NULL, names(theta))))
  ev <- eigen(V, symmetric = TRUE)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(theta))
  out <- sweep(matrix(stats::rnorm(n * length(theta)), n) %*% t(A), 2, theta, "+")
  colnames(out) <- names(theta)
  out
}

# Counts per transect x distance bin (x visit) within truncation.
ds_bin_counts <- function(data, cut, by_visit = FALSE) {
  o <- data$obs_w
  x <- if (data$distance_type == "binned") (o$distbegin + o$distend) / 2 else o$distance
  j <- findInterval(x, cut, rightmost.closed = TRUE, all.inside = TRUE)
  J <- length(cut) - 1
  ids <- data$transects$transect
  i <- match(o$transect, ids)
  if (!by_visit) {
    y <- matrix(0L, length(ids), J)
    for (k in seq_along(j)) y[i[k], j[k]] <- y[i[k], j[k]] + 1L
  } else {
    T <- data$n_visits
    y <- matrix(0L, length(ids), J * T)
    for (k in seq_along(j)) { col <- (o$visit[k] - 1) * J + j[k]; y[i[k], col] <- y[i[k], col] + 1L }
  }
  rownames(y) <- ids
  y
}

ds_hds_cut <- function(data, settings) {
  if (data$distance_type == "binned") return(data$cutpoints[data$cutpoints <= data$truncation_m + 1e-9])
  seq(0, data$truncation_m, length.out = settings$hds_bins + 1)
}

# Mean group size (observed groups within truncation) with SE.
ds_group_size <- function(data) {
  if (!isTRUE(data$clusters)) return(list(mean = 1, se = 0))
  s <- data$obs_w$size
  list(mean = mean(s), se = stats::sd(s) / sqrt(length(s)))
}

# Fewster et al. (2009) R2 (lines) / P2 (points) encounter-rate variance:
# var(n/L) = K / (L^2 (K - 1)) * sum_k l_k^2 (n_k / l_k - n / L)^2
ds_er_var <- function(n_k, l_k) {
  K <- length(l_k); L <- sum(l_k); n <- sum(n_k)
  K / (L^2 * (K - 1)) * sum(l_k^2 * (n_k / l_k - n / L)^2)
}
