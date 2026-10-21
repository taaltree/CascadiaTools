# scr-search: data structure and the internal representation shared by all
# engines (search tracks -> detector cells with effort; samples -> records).
#
# Data structure (returned by simulate()$data and read_data()):
#   tracks    data.frame(track_id, occasion, x, y [, track-level detection
#             covariates]): GPS fixes of each search route in the order walked
#             (one row per fix). One track = one search bout (a dog team on one
#             day, one snow-tracking transect...).
#   samples   data.frame(sample_id, individual, x, y [, track_id] [, occasion]):
#             one row per genotyped sample (scat, hair) at its GPS location.
#             track_id is optional: if missing, a sample is attributed to the
#             nearest search track (on its occasion when occasion is given).
#   n_occasions K
#   mask      data.frame(x, y, <density covariates> [, <cell detection covariates>])
#             habitat pixel centres on a regular grid, or NULL (built from the
#             tracks + buffer)
#   mask_spacing, buffer, cell_size    numeric (map units) or NULL
#   units     "m", "km" or NA
#   density_covariates, detection_covariates (track columns or mask columns)
#   covariate_labels  optional list(column = "Nice label")
#
# Internally coordinates are km (units m -> /1000) relative to the mask
# origin; with unknown units they are divided by the detector cell size.

ss_bucket <- function(x, b = 50) {
  # geometric buckets (factor ~1.25) so replicates of one design share a compiled model
  x <- max(x, 1)
  v <- b; while (v < x) v <- ceiling(v * 1.25)
  as.integer(v)
}

ss_cscale <- function(units, cell_size) {
  if (identical(units, "m")) return(1000)
  if (identical(units, "km")) return(1)
  cell_size
}

ss_std <- function(v) {
  u <- unique(v[!is.na(v)])
  if (all(u %in% c(0, 1))) return(list(x = v, center = 0, scale = 1, binary = TRUE))
  ce <- mean(v, na.rm = TRUE); sc <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  list(x = (v - ce) / sc, center = ce, scale = sc, binary = FALSE)
}

ss_infer_spacing <- function(mask) {
  dx <- diff(sort(unique(round(mask$x, 6)))); dy <- diff(sort(unique(round(mask$y, 6))))
  s <- c(dx[dx > 0], dy[dy > 0])
  if (!length(s)) stop("habitat mask needs at least two distinct x and y values")
  min(s)
}

# Split each track into short pieces (<= maxlen) and return their midpoints
# and lengths. Consecutive fixes of the same track are joined; tracks never
# join each other.
ss_pieces <- function(tracks, maxlen) {
  out <- list()
  sp <- split(seq_len(nrow(tracks)), factor(tracks$track_id, levels = unique(tracks$track_id)))
  for (tid in names(sp)) {
    ix <- sp[[tid]]
    if (length(ix) < 2) next
    x <- tracks$x[ix]; y <- tracks$y[ix]
    dx <- diff(x); dy <- diff(y); L <- sqrt(dx^2 + dy^2)
    k <- pmax(1L, ceiling(L / maxlen))
    seg <- rep(seq_along(L), k)
    pos <- (sequence(k) - 0.5) / rep(k, k)
    out[[tid]] <- data.frame(track_id = tid, x = x[seg] + pos * dx[seg], y = y[seg] + pos * dy[seg],
                             len = L[seg] / rep(k, k))
  }
  p <- do.call(rbind, out); rownames(p) <- NULL
  p[p$len > 0, , drop = FALSE]
}

# Quick sigma guess from the spread of each individual's samples (pooled
# RPSV). Only used for default mask buffer, local-evaluation radius and inits.
ss_sigma_guess <- function(samples) {
  sp <- split(seq_len(nrow(samples)), samples$individual)
  ss <- 0; df <- 0
  for (ix in sp) if (length(ix) > 1) {
    mx <- mean(samples$x[ix]); my <- mean(samples$y[ix])
    ss <- ss + sum((samples$x[ix] - mx)^2 + (samples$y[ix] - my)^2); df <- df + length(ix) - 1
  }
  if (df < 1 || ss <= 0) return(NA_real_)
  sqrt(ss / (2 * df))
}

ss_make_mask <- function(pts, buffer, spacing) {
  xr <- range(pts$x) + c(-buffer, buffer); yr <- range(pts$y) + c(-buffer, buffer)
  m <- expand.grid(x = seq(xr[1] + spacing / 2, xr[2], by = spacing), y = seq(yr[1] + spacing / 2, yr[2], by = spacing))
  sub <- pts[unique(round(seq(1, nrow(pts), length.out = min(nrow(pts), 3000)))), ]
  d <- ss_min_dist(m$x, m$y, sub$x, sub$y)
  m[d <= buffer, , drop = FALSE]
}

ss_min_dist <- function(x1, y1, x2, y2, chunk = 2000) {
  out <- numeric(length(x1))
  for (s in seq(1, length(x1), by = chunk)) {
    ii <- s:min(length(x1), s + chunk - 1)
    d2 <- outer(x1[ii], x2, "-")^2 + outer(y1[ii], y2, "-")^2
    out[ii] <- sqrt(apply(d2, 1, min))
  }
  out
}

# Everything an engine needs, on internal scale.
ss_prepare <- function(data, settings = list()) {
  notes <- character()
  units <- if (is.null(data$units)) NA_character_ else data$units
  tracks <- data$tracks; tracks$track_id <- as.character(tracks$track_id)
  samples <- data$samples; samples$individual <- as.character(samples$individual)
  if (!"track_id" %in% names(samples)) samples$track_id <- NA_character_
  samples$track_id <- as.character(samples$track_id)
  K <- if (!is.null(data$n_occasions)) data$n_occasions else max(tracks$occasion)
  dcov <- as.character(unlist(if (!is.null(settings$density_covariates)) settings$density_covariates else data$density_covariates))
  wcov <- as.character(unlist(if (!is.null(settings$detection_covariates)) settings$detection_covariates else data$detection_covariates))
  effort_model <- if (is.null(settings$effort_model)) "track_length" else settings$effort_model
  observation <- if (is.null(settings$observation)) "poisson" else settings$observation
  if (!observation %in% c("poisson", "bernoulli")) stop("settings$observation must be 'poisson' or 'bernoulli'")
  # ---- track-level table (one row per track)
  tinfo <- tracks[!duplicated(tracks$track_id), c("track_id", "occasion", intersect(wcov, names(tracks))), drop = FALSE]
  # ---- cell size, mask
  cs <- if (!is.null(settings$cell_size)) settings$cell_size else data$cell_size
  sig_guess <- ss_sigma_guess(samples)
  if (is.null(cs)) {
    cs <- if (is.finite(sig_guess)) signif(sig_guess / 1.5, 2) else stop("give data: cell_size (detector cell side, map units)")
    notes <- c(notes, sprintf("Detector cell size not given: using %s %s (about 2/3 of a rough sigma guess).", format(cs), ss_unit_word(units)))
  }
  mask <- data$mask
  if (is.null(mask)) {
    sg <- if (is.finite(sig_guess)) sig_guess else 2 * cs
    buffer <- if (!is.null(data$buffer)) data$buffer else 4 * sg
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else min(cs, sg / 2)
    mask <- ss_make_mask(tracks, buffer, spacing)
    while (nrow(mask) > (if (is.null(settings$max_mask_points)) 3000 else settings$max_mask_points)) {
      spacing <- spacing * 1.2; mask <- ss_make_mask(tracks, buffer, spacing)
    }
    notes <- c(notes, sprintf("Habitat mask built from the search tracks: buffer %s, pixel spacing %s %s, %d pixels.",
                              format(signif(buffer, 3)), format(signif(spacing, 3)), ss_unit_word(units), nrow(mask)))
    if (length(dcov)) stop("density covariates need a habitat mask CSV or raster with columns: ", paste(dcov, collapse = ", "))
  } else {
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else ss_infer_spacing(mask)
    buffer <- min(min(tracks$x) - min(mask$x), max(mask$x) - max(tracks$x),
                  min(tracks$y) - min(mask$y), max(mask$y) - max(tracks$y)) + spacing / 2
  }
  if (length(dcov)) {
    miss <- setdiff(dcov, names(mask)); if (length(miss)) stop("habitat mask has no column(s): ", paste(miss, collapse = ", "))
    bad <- rowSums(is.na(mask[, dcov, drop = FALSE])) > 0
    if (any(bad)) { mask <- mask[!bad, , drop = FALSE]
      notes <- c(notes, sprintf("%d habitat pixels with a missing density covariate were treated as non-habitat.", sum(bad))) }
  }
  G <- nrow(mask)
  cscale <- ss_cscale(units, cs)
  ox <- min(mask$x) - spacing / 2; oy <- min(mask$y) - spacing / 2
  mxy <- cbind((mask$x - ox) / cscale, (mask$y - oy) / cscale)
  pix <- spacing / cscale; csz <- cs / cscale
  # ---- pieces -> detector cells and units (cell x track)
  tr_int <- data.frame(track_id = tracks$track_id, x = (tracks$x - ox) / cscale, y = (tracks$y - oy) / cscale)
  pcs <- ss_pieces(tr_int, csz / 20)
  pcs$cx <- floor(pcs$x / csz); pcs$cy <- floor(pcs$y / csz)
  pcs$cellkey <- paste(pcs$cx, pcs$cy)
  pcs$unitkey <- paste(pcs$cellkey, pcs$track_id, sep = "|")
  ukeys <- unique(pcs$unitkey)
  ueff <- tapply(pcs$len, factor(pcs$unitkey, levels = ukeys), sum)
  first <- pcs[match(ukeys, pcs$unitkey), ]
  units_df <- data.frame(unit = ukeys, cellkey = first$cellkey, cx = first$cx, cy = first$cy,
                         track_id = first$track_id, effort = as.numeric(ueff), stringsAsFactors = FALSE)
  units_df$occasion <- tinfo$occasion[match(units_df$track_id, tinfo$track_id)]
  # ---- samples -> units (nearest piece of the sample's track, or of any track)
  s_int <- data.frame(x = (samples$x - ox) / cscale, y = (samples$y - oy) / cscale)
  pcs_occ <- tinfo$occasion[match(pcs$track_id, tinfo$track_id)]
  su <- character(nrow(samples)); sdist <- numeric(nrow(samples))
  for (r in seq_len(nrow(samples))) {
    cand <- if (!is.na(samples$track_id[r]) && samples$track_id[r] %in% pcs$track_id) which(pcs$track_id == samples$track_id[r])
            else if ("occasion" %in% names(samples) && !is.na(samples$occasion[r]) && any(pcs_occ == samples$occasion[r])) which(pcs_occ == samples$occasion[r])
            else seq_len(nrow(pcs))
    d2 <- (pcs$x[cand] - s_int$x[r])^2 + (pcs$y[cand] - s_int$y[r])^2
    k <- cand[which.min(d2)]
    su[r] <- pcs$unitkey[k]; sdist[r] <- sqrt(min(d2))
  }
  far <- sdist > csz / 2
  if (any(far)) notes <- c(notes, sprintf("%d samples were more than half a detector cell (%s %s) from the track they were attributed to (max %s %s); check sample coordinates and track IDs.",
                                          sum(far), format(signif(cs / 2, 3)), ss_unit_word(units), format(signif(max(sdist) * cscale, 3)), ss_unit_word(units)))
  samples$unit <- su; samples$track_assigned <- units_df$track_id[match(su, units_df$unit)]
  samples$dist_to_track <- sdist * cscale
  # ---- effort model: how effort enters (track_length = correct; others are
  # the shortcuts used as negative controls / sensitivity analyses)
  if (effort_model == "visits") {
    units_df$effort <- csz
    notes <- c(notes, "Effort model 'visits': every track passing through a cell counts as one standard visit, regardless of the length searched.")
  } else if (effort_model == "uniform") {
    # every cell in the bounding box of the searches gets the same effort, track identity ignored
    cxr <- range(units_df$cx); cyr <- range(units_df$cy)
    grid <- expand.grid(cx = cxr[1]:cxr[2], cy = cyr[1]:cyr[2])
    grid$cellkey <- paste(grid$cx, grid$cy)
    units_df <- data.frame(unit = grid$cellkey, cellkey = grid$cellkey, cx = grid$cx, cy = grid$cy, track_id = NA_character_,
                           effort = csz, occasion = NA, stringsAsFactors = FALSE)
    samples$unit <- sub("\\|.*$", "", samples$unit)
    if (length(wcov)) notes <- c(notes, "Effort model 'uniform': detection covariates dropped (search effort ignored).")
    wcov <- character()
    notes <- c(notes, "Effort model 'uniform': search effort IGNORED (every cell of the searched rectangle treated as searched equally). Use only as a sensitivity analysis.")
  } else if (effort_model != "track_length") stop("settings$effort_model must be track_length, visits or uniform")
  U <- nrow(units_df)
  ckeys <- unique(units_df$cellkey)
  C <- length(ckeys)
  ucell <- match(units_df$cellkey, ckeys)
  cxy <- cbind((units_df$cx[match(ckeys, units_df$cellkey)] + 0.5) * csz, (units_df$cy[match(ckeys, units_df$cellkey)] + 0.5) * csz)
  ccentre <- cxy   # geometric cell centres (maps, cell-level covariates)
  # detector location of each cell: the geometric centre, or the length-weighted
  # centroid of the track searched inside the cell (a closer approximation of the
  # hazard integrated along the track; Efford 2011)
  detloc <- if (is.null(settings$detector_location)) "centroid" else settings$detector_location
  if (!detloc %in% c("centre", "centroid")) stop("settings$detector_location must be 'centroid' or 'centre'")
  if (detloc == "centroid" && effort_model != "uniform") {
    fk <- factor(pcs$cellkey, levels = ckeys)
    wl <- tapply(pcs$len, fk, sum)
    cxy <- cbind(as.numeric(tapply(pcs$x * pcs$len, fk, sum) / wl), as.numeric(tapply(pcs$y * pcs$len, fk, sum) / wl))
  } else detloc <- "centre"
  # ---- detection covariates per unit: track columns, else mask columns at the cell centre
  wsc <- list(); W <- matrix(0, U, 0)
  for (v in wcov) {
    if (v %in% names(tinfo)) val <- tinfo[[v]][match(units_df$track_id, tinfo$track_id)]
    else if (v %in% names(mask)) {
      ccx <- ccentre[ucell, 1]; ccy <- ccentre[ucell, 2]
      nn <- vapply(seq_len(U), function(u) which.min((mxy[, 1] - ccx[u])^2 + (mxy[, 2] - ccy[u])^2), 1L)
      val <- mask[[v]][nn]
    } else stop("detection covariate '", v, "' is neither a column of the track file nor of the habitat mask")
    if (is.character(val) || is.factor(val)) {
      lv <- sort(unique(as.character(val)))
      if (length(lv) != 2) stop("categorical detection covariate '", v, "' must have exactly 2 levels (code others as 0/1 columns)")
      notes <- c(notes, sprintf("Detection covariate '%s': '%s' = 1, '%s' = 0 (reference).", v, lv[2], lv[1]))
      val <- as.integer(as.character(val) == lv[2])
    }
    if (anyNA(val)) stop("detection covariate '", v, "' has missing values on some searched cells/tracks")
    s <- ss_std(as.numeric(val)); wsc[[v]] <- s[c("center", "scale", "binary")]; W <- cbind(W, s$x)
  }
  colnames(W) <- wcov
  # ---- density covariates
  dsc <- list(); XD <- matrix(0, G, 0)
  for (v in dcov) { s <- ss_std(mask[[v]]); dsc[[v]] <- s[c("center", "scale", "binary")]; XD <- cbind(XD, s$x) }
  colnames(XD) <- dcov
  # ---- records: individual x unit counts
  ids <- unique(samples$individual)
  n <- length(ids)
  rec <- stats::aggregate(list(count = rep(1L, nrow(samples))),
                          list(ind = match(samples$individual, ids), unit = match(samples$unit, units_df$unit)), sum)
  rec <- rec[order(rec$ind, rec$unit), ]
  if (observation == "bernoulli") rec$count <- 1L
  rec$cell <- ucell[rec$unit]
  # ---- local evaluation radius and pixel-cell pairs
  sg_int <- if (is.finite(sig_guess)) sig_guess / cscale else 2 * csz
  d2gc <- outer(mxy[, 1], cxy[, 1], "-")^2 + outer(mxy[, 2], cxy[, 2], "-")^2
  span <- 0
  for (i in seq_len(n)) { cc <- unique(rec$cell[rec$ind == i]); if (length(cc) > 1) span <- max(span, sqrt(max(stats::dist(cxy[cc, , drop = FALSE])))) }
  R <- if (!is.null(settings$local_radius)) settings$local_radius / cscale else 6 * sg_int
  R <- max(R, span / 2 + 3 * csz)
  if (isTRUE(settings$no_local)) R <- Inf
  pr <- which(d2gc <= R^2, arr.ind = TRUE)
  pr <- pr[order(pr[, 1], pr[, 2]), , drop = FALSE]
  # per-individual candidate pixels: within R of every cell where it was detected
  sets <- vector("list", n)
  for (i in seq_len(n)) {
    cc <- unique(rec$cell[rec$ind == i])
    ok <- rep(TRUE, G); for (c1 in cc) ok <- ok & d2gc[, c1] <= R^2
    sets[[i]] <- which(ok)
    if (!length(sets[[i]])) stop("internal: empty local window for an individual; increase settings$local_radius")
  }
  # ---- data augmentation size, priors
  M <- settings$M
  if (is.null(M)) M <- max(4 * n, n + 150)
  width <- max(diff(range(mxy[, 1])), diff(range(mxy[, 2]))) + pix
  sigma_max <- if (!is.null(settings$sigma_max)) settings$sigma_max / cscale else width / 2
  tl <- tapply(units_df$effort, units_df$occasion, sum)
  list(n = n, ids = ids, K = K, G = G, C = C, U = U, M = M, mask = mask, mxy = mxy, pix = pix, area = G * pix^2,
       spacing = spacing, buffer = buffer, cscale = cscale, origin = c(ox, oy), units = units, cell_size = cs, csz = csz,
       cxy = cxy, ccentre = ccentre, detector_location = detloc, ckeys = ckeys, units_df = units_df, ucell = ucell, logE = log(units_df$effort), W = W, wcov = wcov, wscale = wsc,
       XD = XD, dcov = dcov, dscale = dsc, rec = rec, sets = sets, pairs = pr, pd2 = d2gc[pr], d2gc = d2gc, R = R,
       samples = samples, tracks = tracks, tinfo = tinfo, effort_model = effort_model, observation = observation,
       total_effort = sum(units_df$effort) * cscale, effort_by_occasion = tl * cscale,
       sigma_guess = sig_guess, sigma_max = sigma_max, n_samples = nrow(samples),
       n_spatial_recap = sum(tapply(rec$cell, rec$ind, function(v) length(unique(v))) > 1), notes = notes)
}

ss_unit_word <- function(units) if (is.null(units) || is.na(units)) "map units" else units
ss_density_name <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "density_per_km2" else "density_per_unit2"
ss_area_label <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "km²" else "map unit²"
ss_len_label <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "km" else "map unit"
# lambda0 is per internal length unit: per km when units known
ss_lam_name <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "lambda0_per_km" else "lambda0_per_unit"
