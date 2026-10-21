# scr-open: data structure, habitat mask and the internal representation
# shared by every engine. Conventions follow modules/scr (copied, not sourced,
# so the two modules can evolve independently).
#
# Data structure (returned by simulate()$data and read_data()):
#   traps        data.frame(id, x, y)
#   usage        J x T integer matrix: number of secondary occasions (0..K)
#                trap j operated in primary period t (rows = traps$id)
#   captures     data.frame(individual, period, occasion, trap): one row per
#                detection; period is the index 1..T, occasion 1..K. Duplicate
#                rows (same animal, period, occasion, trap) are collapsed:
#                proximity detectors record at most one detection per occasion.
#   n_occasions  K (secondary occasions per primary period)
#   n_periods    T (primary periods; population open between, closed within)
#   period_labels character(T), e.g. years
#   mask         data.frame(x, y, <density covariates>) of habitat pixel
#                centres on a regular grid, or NULL (built from traps + buffer)
#   mask_spacing, buffer   numeric (map units) or NULL
#   units        "m", "km" or NA
#   density_covariates     character
#   covariate_labels       named list of display labels
#
# Internally distances are in km when units are m/km (cscale = 1000 or 1),
# otherwise divided by the median nearest-neighbour trap spacing.

so_cscale <- function(units, traps) {
  if (identical(units, "m")) return(1000)
  if (identical(units, "km")) return(1)
  so_trap_spacing(traps)
}

so_trap_spacing <- function(traps) {
  xy <- as.matrix(traps[, c("x", "y")])
  if (nrow(xy) < 2) return(1)
  d <- as.matrix(stats::dist(xy)); diag(d) <- Inf
  stats::median(apply(d, 1, min))
}

# Pooled RPSV-type sigma guess from repeat detections within the same period
# (Calhoun & Casby 1958; secr::RPSV(CC = TRUE)); used to size a default mask.
so_sigma_guess <- function(captures, traps) {
  xy <- traps[match(captures$trap, traps$id), c("x", "y")]
  sp <- split(seq_len(nrow(captures)), captures$individual)
  ss <- 0; df <- 0
  for (ix in sp) if (length(ix) > 1) {
    m <- colMeans(xy[ix, , drop = FALSE])
    ss <- ss + sum((xy$x[ix] - m[1])^2 + (xy$y[ix] - m[2])^2)
    df <- df + length(ix) - 1
  }
  if (df < 1 || ss <= 0) return(NA_real_)
  sqrt(ss / (2 * df))
}

so_make_mask <- function(traps, buffer, spacing) {
  xr <- range(traps$x) + c(-buffer, buffer); yr <- range(traps$y) + c(-buffer, buffer)
  gx <- seq(xr[1] + spacing / 2, xr[2], by = spacing)
  gy <- seq(yr[1] + spacing / 2, yr[2], by = spacing)
  m <- expand.grid(x = gx, y = gy)
  keep <- so_dist_to_traps(m, traps) <= buffer
  m[keep, , drop = FALSE]
}

so_dist_to_traps <- function(pts, traps) {
  d2 <- outer(pts$x, traps$x, "-")^2 + outer(pts$y, traps$y, "-")^2
  sqrt(apply(d2, 1, min))
}

so_infer_spacing <- function(mask) {
  dx <- diff(sort(unique(round(mask$x, 6)))); dy <- diff(sort(unique(round(mask$y, 6))))
  s <- c(dx[dx > 0], dy[dy > 0])
  if (!length(s)) stop("habitat mask needs at least two distinct x and y values")
  min(s)
}

so_std <- function(v) {
  u <- unique(v[!is.na(v)])
  if (all(u %in% c(0, 1))) return(list(x = v, center = 0, scale = 1, binary = TRUE))
  ce <- mean(v, na.rm = TRUE); sc <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  list(x = (v - ce) / sc, center = ce, scale = sc, binary = FALSE)
}

# Collapse capture rows to an n x J x T array of counts (occasions detected).
so_counts <- function(data) {
  tr <- as.character(data$traps$id)
  TT <- data$n_periods
  cap <- unique(data.frame(individual = as.character(data$captures$individual),
                           period = as.integer(data$captures$period),
                           occasion = as.integer(data$captures$occasion),
                           trap = as.character(data$captures$trap)))
  ids <- unique(cap$individual)
  y <- array(0L, c(length(ids), length(tr), TT), dimnames = list(ids, tr, NULL))
  ijt <- cbind(match(cap$individual, ids), match(cap$trap, tr), cap$period)
  for (r in seq_len(nrow(ijt))) y[ijt[r, , drop = FALSE]] <- y[ijt[r, , drop = FALSE]] + 1L
  list(y = y, captures = cap, ids = ids)
}

so_period_labels <- function(data) {
  if (!is.null(data$period_labels)) as.character(data$period_labels) else as.character(seq_len(data$n_periods))
}

so_prepare <- function(data, settings = list()) {
  traps <- data$traps; traps$id <- as.character(traps$id)
  J <- nrow(traps); K <- data$n_occasions; TT <- data$n_periods
  if (TT < 2) stop("an open-population model needs at least 2 primary periods (use module 'scr' for one period)")
  usage <- if (is.null(data$usage)) matrix(as.integer(K), J, TT) else as.matrix(data$usage)
  storage.mode(usage) <- "integer"
  cn <- so_counts(data)
  y <- cn$y; n <- dim(y)[1]
  dcov <- if (!is.null(settings$density_covariates)) unlist(settings$density_covariates) else unlist(data$density_covariates)
  dcov <- as.character(dcov)
  notes <- character()
  sig_guess <- so_sigma_guess(cn$captures[, c("individual", "trap")], traps)
  mask <- data$mask
  if (is.null(mask)) {
    sp0 <- so_trap_spacing(traps)
    sg <- if (is.finite(sig_guess) && sig_guess > 0) sig_guess else sp0
    buffer <- if (!is.null(data$buffer)) data$buffer else 4 * sg
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else sg / 2
    mask <- so_make_mask(traps, buffer, spacing)
    mx <- if (is.null(settings$max_mask_points)) 2500 else settings$max_mask_points
    while (nrow(mask) > mx) { spacing <- spacing * 1.25; mask <- so_make_mask(traps, buffer, spacing) }
    notes <- c(notes, sprintf("Habitat mask built from the traps: buffer %s, pixel spacing %s (%s units), %d pixels.",
                              format(signif(buffer, 3)), format(signif(spacing, 3)), so_unit_word(data$units), nrow(mask)))
    if (length(dcov) && !all(dcov %in% names(mask)))
      stop("density covariates need a habitat mask file with those columns: ", paste(dcov, collapse = ", "))
  } else {
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else so_infer_spacing(mask)
    buffer <- min(min(traps$x) - min(mask$x), max(mask$x) - max(traps$x),
                  min(traps$y) - min(mask$y), max(mask$y) - max(traps$y)) + spacing / 2
  }
  if (!is.null(settings$buffer)) {
    keep <- so_dist_to_traps(mask, traps) <= settings$buffer
    mask <- mask[keep, , drop = FALSE]; buffer <- settings$buffer
    notes <- c(notes, sprintf("State space trimmed to %s units from the nearest trap (settings$buffer).", format(settings$buffer)))
  }
  if (length(dcov)) {
    bad <- rowSums(is.na(mask[, dcov, drop = FALSE])) > 0
    if (any(bad)) {
      mask <- mask[!bad, , drop = FALSE]
      notes <- c(notes, sprintf("%d habitat pixels with a missing density covariate were treated as non-habitat.", sum(bad)))
    }
  }
  rownames(mask) <- NULL
  G <- nrow(mask)
  cscale <- so_cscale(data$units, traps)
  ox <- min(mask$x) - spacing / 2; oy <- min(mask$y) - spacing / 2
  mxy <- cbind((mask$x - ox) / cscale, (mask$y - oy) / cscale)
  txy <- cbind((traps$x - ox) / cscale, (traps$y - oy) / cscale)
  d2 <- outer(mxy[, 1], txy[, 1], "-")^2 + outer(mxy[, 2], txy[, 2], "-")^2
  pix <- spacing / cscale
  dsc <- list(); XD <- matrix(0, G, 0)
  for (v in dcov) { s <- so_std(mask[[v]]); dsc[[v]] <- s[c("center", "scale", "binary")]; XD <- cbind(XD, s$x) }
  colnames(XD) <- dcov
  # sparse detections sorted by individual then period
  nz <- which(y > 0, arr.ind = TRUE)
  nz <- nz[order(nz[, 1], nz[, 3], nz[, 2]), , drop = FALSE]
  nzc <- y[nz]
  nzstart <- match(seq_len(n), nz[, 1])
  nzend <- vapply(seq_len(n), function(i) max(which(nz[, 1] == i)), 1L)
  seen <- apply(y, c(1, 3), sum) > 0
  if (n == 1) seen <- matrix(seen, 1)
  M <- settings$M
  if (is.null(M)) M <- max(6 * n, n + 300)
  M <- as.integer(M)
  width <- max(diff(range(mxy[, 1])), diff(range(mxy[, 2]))) + pix
  sigma_max <- if (!is.null(settings$sigma_max)) settings$sigma_max / cscale else width / 2
  list(n = n, J = J, K = K, TT = TT, G = G, M = M, y = y, usage = usage, ids = cn$ids,
       captures = cn$captures, traps = traps, mask = mask, mxy = mxy, txy = txy, d2 = d2,
       pix = pix, area = G * pix^2, spacing = spacing, buffer = buffer, cscale = cscale, origin = c(ox, oy),
       units = if (is.null(data$units)) NA_character_ else data$units,
       XD = XD, dcov = dcov, dscale = dsc, period_labels = so_period_labels(data),
       nzi = nz[, 1], nzj = nz[, 2], nzt = nz[, 3], nzc = nzc, nzstart = nzstart, nzend = nzend,
       seen = seen, sigma_max = sigma_max, sigma_guess = sig_guess, trap_spacing = so_trap_spacing(traps),
       n_spatial_recap = sum(apply(y, 1, function(a) sum(rowSums(a) > 0)) > 1),
       n_by_period = colSums(seen), n_recaptured_between = sum(rowSums(seen) > 1),
       notes = notes)
}

so_unit_word <- function(units) if (is.null(units) || is.na(units)) "coordinate" else units
so_unit_len <- function(prep) if (is.na(prep$units)) "map units" else prep$units
so_known_units <- function(units) !is.null(units) && !is.na(units) && units %in% c("m", "km")
so_area_label <- function(units) if (so_known_units(units)) "km²" else "map unit²"
so_density_stem <- function(units) if (so_known_units(units)) "density_per_km2" else "density_per_unit2"
# internal area -> reported area (km2 when units known, else map unit2)
so_area_conv <- function(prep) if (so_known_units(prep$units)) 1 else prep$cscale^2
