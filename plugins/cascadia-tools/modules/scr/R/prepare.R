# scr: data structure, habitat mask, and the internal representation shared
# by every engine.
#
# Data structure (returned by simulate()$data and read_data()):
#   traps     data.frame(id, x, y, <detector covariates>)
#   usage     J x K matrix (rows = traps$id): 0/1 operating, or non-negative
#             effort (e.g. days a hair snare was active / nominal days).
#             Effort acts on the hazard scale (Efford, Borchers & Mowat 2013):
#             Bernoulli p = 1 - (1 - p)^effort, Poisson mean = effort x lambda.
#   captures  data.frame(individual, occasion, trap [, count]): one row per
#             detection (camera/live trap) or per genetic sample (NGS). For
#             binomial/Bernoulli models repeat rows of the same individual x
#             trap x occasion collapse to one detection; for the Poisson
#             model they are counted (number of samples).
#   n_occasions K
#   individuals optional data.frame(individual, sex) with sex "F", "M" or NA
#   occasion_covariates optional data.frame with K rows (one per occasion)
#   mask      data.frame(x, y, <density covariates>) of habitat pixel centres
#             on a regular grid, or NULL (built from traps [+ polygon/raster])
#   mask_kind "csv", "traps", "polygon" or "raster" (how the mask was made)
#   habitat_polygon (sf), raster_file, habitat_layer: sources to rebuild the
#             mask with another buffer (buffer check)
#   mask_spacing, buffer   numeric (coordinate units) or NULL
#   units     "m", "km" or NA (unknown)
#   Model choices (each can be overridden by the same-named settings key):
#   obs_model ("binomial", "bernoulli", "poisson"), density_covariates,
#   detection_covariates (detector-level, on p0/lambda0),
#   occasion_detection_covariates (occasion-level, on p0/lambda0),
#   sigma_covariates (detector-level, on sigma), sex_effects (subset of
#   c("p0", "sigma")), behaviour ("none" or "trap": trap-specific response)
#
# Internally all distances are rescaled so numbers are O(1): coordinates in
# metres become km (cscale = 1000), km stay km, unknown units are divided by
# the median nearest-neighbour trap spacing.

scr_cscale <- function(units, traps) {
  if (identical(units, "m")) return(1000)
  if (identical(units, "km")) return(1)
  scr_trap_spacing(traps)
}

scr_trap_spacing <- function(traps) {
  xy <- as.matrix(traps[, c("x", "y")])
  if (nrow(xy) < 2) return(1)
  d <- as.matrix(stats::dist(xy)); diag(d) <- Inf
  stats::median(apply(d, 1, min))
}

# Quick sigma guess from the spread of repeat detections (pooled RPSV,
# Calhoun & Casby 1958 as used by secr::RPSV(CC = TRUE)). Used only to size
# a default mask when the researcher does not provide one.
scr_sigma_guess <- function(captures, traps) {
  xy <- traps[match(captures$trap, traps$id), c("x", "y")]
  sp <- split(seq_len(nrow(captures)), captures$individual)
  ss <- 0; df <- 0
  for (ix in sp) if (length(ix) > 1) {
    m <- colMeans(xy[ix, , drop = FALSE])
    ss <- ss + sum((xy$x[ix] - m[1])^2 + (xy$y[ix] - m[2])^2)
    df <- df + length(ix) - 1
  }
  # all repeat detections at the same trap -> no spatial information (NA, not 0)
  if (df < 1 || ss <= 0) return(NA_real_)
  sqrt(ss / (2 * df))
}

# Habitat mask from the trap layout: grid of pixel centres at `spacing`
# within `buffer` of the nearest trap (secr "trapbuffer" type).
scr_make_mask <- function(traps, buffer, spacing) {
  xr <- range(traps$x) + c(-buffer, buffer); yr <- range(traps$y) + c(-buffer, buffer)
  gx <- seq(xr[1] + spacing / 2, xr[2], by = spacing)
  gy <- seq(yr[1] + spacing / 2, yr[2], by = spacing)
  m <- expand.grid(x = gx, y = gy)
  keep <- scr_dist_to_traps(m, traps) <= buffer
  m[keep, , drop = FALSE]
}

scr_dist_to_traps <- function(pts, traps) {
  d2 <- outer(pts$x, traps$x, "-")^2 + outer(pts$y, traps$y, "-")^2
  sqrt(apply(d2, 1, min))
}

# Build the habitat mask from its sources (traps + buffer, optionally
# intersected with a habitat polygon and/or a habitat raster layer, with
# density covariates extracted from a raster). Used by read_data() and by the
# buffer check (rebuilds the mask with a larger buffer).
scr_build_mask <- function(data, buffer, spacing) {
  traps <- data$traps
  m <- scr_make_mask(traps, buffer, spacing)
  notes <- character()
  if (!is.null(data$habitat_polygon)) {
    eco_require("sf")
    poly <- sf::st_union(sf::st_geometry(data$habitat_polygon))
    pts <- sf::st_as_sf(m, coords = c("x", "y"), crs = sf::st_crs(poly))
    inside <- lengths(sf::st_intersects(pts, poly)) > 0
    notes <- c(notes, sprintf("%d of %d candidate mask pixels lie outside the habitat polygon and were dropped.", sum(!inside), length(inside)))
    m <- m[inside, , drop = FALSE]
  }
  if (!is.null(data$raster_file)) {
    eco_require("terra")
    r <- terra::rast(data$raster_file)
    v <- terra::extract(r, as.matrix(m[, c("x", "y")]))
    hl <- data$habitat_layer
    if (!is.null(hl)) {
      if (!hl %in% names(v)) stop("raster has no habitat layer named '", hl, "' (layers: ", paste(names(r), collapse = ", "), ")", call. = FALSE)
      h <- v[[hl]]
      keep <- !is.na(h) & h > 0
      notes <- c(notes, sprintf("%d mask pixels are non-habitat in raster layer '%s' and were dropped.", sum(!keep), hl))
      m <- m[keep, , drop = FALSE]; v <- v[keep, , drop = FALSE]
    }
    for (k in names(v)) if (k %in% data$density_covariates) m[[k]] <- v[[k]]
    miss <- setdiff(data$density_covariates, names(m))
    if (length(miss)) stop("raster has no layer named: ", paste(miss, collapse = ", "), " (layers: ", paste(names(r), collapse = ", "), ")", call. = FALSE)
  }
  rownames(m) <- NULL
  attr(m, "notes") <- notes
  m
}

scr_infer_spacing <- function(mask) {
  dx <- diff(sort(unique(round(mask$x, 6)))); dy <- diff(sort(unique(round(mask$y, 6))))
  s <- c(dx[dx > 0], dy[dy > 0])
  if (!length(s)) stop("habitat mask needs at least two distinct x and y values")
  min(s)
}

# Model choice: settings override the data-level default
scr_choice <- function(settings, data, key, default = NULL) {
  if (!is.null(settings[[key]])) return(settings[[key]])
  if (!is.null(data[[key]])) return(data[[key]])
  default
}

# Detection records -> per-individual x trap x occasion counts.
# obs "poisson": rows are counted (samples); otherwise collapsed to 0/1.
scr_counts <- function(data, obs = "binomial") {
  tr <- as.character(data$traps$id)
  K <- data$n_occasions
  cap <- data.frame(individual = as.character(data$captures$individual),
                    occasion = as.integer(data$captures$occasion),
                    trap = as.character(data$captures$trap))
  cnt <- if (!is.null(data$captures$count)) as.numeric(data$captures$count) else rep(1, nrow(cap))
  if (obs != "poisson") { keep <- !duplicated(cap); cap <- cap[keep, , drop = FALSE]; cnt <- rep(1, nrow(cap)) }
  ids <- unique(cap$individual)
  y3 <- array(0, c(length(ids), length(tr), K), dimnames = list(ids, tr, NULL))
  ix <- cbind(match(cap$individual, ids), match(cap$trap, tr), cap$occasion)
  for (r in seq_len(nrow(ix))) y3[ix[r, 1], ix[r, 2], ix[r, 3]] <- y3[ix[r, 1], ix[r, 2], ix[r, 3]] + cnt[r]
  y <- apply(y3, c(1, 2), sum)
  if (!is.matrix(y)) y <- matrix(y, length(ids), length(tr), dimnames = list(ids, tr))
  list(y = y, y3 = y3, captures = unique(cap), ids = ids)
}

scr_std <- function(v) {
  # 0/1 indicators stay 0/1 (reference = 0); continuous -> mean 0, SD 1
  u <- unique(v[!is.na(v)])
  if (all(u %in% c(0, 1))) return(list(x = v, center = 0, scale = 1, binary = TRUE))
  ce <- mean(v, na.rm = TRUE); sc <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  list(x = (v - ce) / sc, center = ce, scale = sc, binary = FALSE)
}

scr_sex_code <- function(s) {
  s <- toupper(trimws(as.character(s)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("F", "FEMALE", "0", "FEM")] <- "F"
  out[s %in% c("M", "MALE", "1", "MAL")] <- "M"
  out
}

# Everything an engine needs, on internal scale.
scr_prepare <- function(data, settings = list()) {
  traps <- data$traps
  traps$id <- as.character(traps$id)
  J <- nrow(traps)
  K <- data$n_occasions
  usage <- if (is.null(data$usage)) matrix(1, J, K) else { u <- as.matrix(data$usage); storage.mode(u) <- "double"; u }
  notes <- character()
  if (isTRUE(settings$ignore_effort)) {
    # negative control only: pretend every detector worked with full effort on every occasion
    usage[] <- 1
    notes <- c(notes, "WARNING: settings ignore_effort = true: detector effort / operation was ignored (negative control; do not use for inference).")
  }
  obs <- scr_choice(settings, data, "obs_model", "binomial")
  if (!obs %in% c("binomial", "bernoulli", "poisson")) stop("obs_model must be binomial, bernoulli or poisson (got '", obs, "')")
  binary_usage <- all(usage %in% c(0, 1))
  if (obs == "binomial" && !binary_usage) {
    obs <- "bernoulli"
    notes <- c(notes, "Usage contains partial effort (values other than 0/1), so the binomial model was fitted occasion by occasion as Bernoulli with p = 1 - (1 - p)^effort.")
  }
  cn <- scr_counts(data, obs)
  y <- cn$y; y3 <- cn$y3; n <- nrow(y)
  Kj <- rowSums(usage > 0)
  dcov <- as.character(unlist(scr_choice(settings, data, "density_covariates", character())))
  pcov <- as.character(unlist(scr_choice(settings, data, "detection_covariates", character())))
  ocov <- as.character(unlist(scr_choice(settings, data, "occasion_detection_covariates", character())))
  scov <- as.character(unlist(scr_choice(settings, data, "sigma_covariates", character())))
  sexeff <- as.character(unlist(scr_choice(settings, data, "sex_effects", character())))
  sexeff <- intersect(sexeff, c("p0", "sigma"))
  behav <- scr_choice(settings, data, "behaviour", "none")
  if (!behav %in% c("none", "trap")) stop("behaviour must be 'none' or 'trap' (trap-specific learned response, secr 'bk')")
  # ---- mask
  sig_guess <- scr_sigma_guess(cn$captures, traps)
  mask <- data[["mask"]]
  if (is.null(mask)) {
    sp0 <- scr_trap_spacing(traps)
    sg <- if (is.finite(sig_guess) && sig_guess > 0) sig_guess else sp0
    buffer <- if (!is.null(data[["buffer"]])) data[["buffer"]] else 4 * sg
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else sg / 2
    mask <- scr_build_mask(data, buffer, spacing)
    while (nrow(mask) > (if (is.null(settings$max_mask_points)) 4000 else settings$max_mask_points)) {
      spacing <- spacing * 1.25; mask <- scr_build_mask(data, buffer, spacing)
    }
    notes <- c(notes, attr(mask, "notes"), sprintf("Habitat mask built from the traps: buffer %s, pixel spacing %s (%s units), %d pixels.",
                              format(signif(buffer, 3)), format(signif(spacing, 3)), scr_unit_word(data$units), nrow(mask)))
    if (length(dcov) && !all(dcov %in% names(mask)))
      stop("density covariates need a habitat mask file or raster with those columns: ", paste(dcov, collapse = ", "))
  } else {
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else scr_infer_spacing(mask)
    # effective buffer: smallest distance from the trap array's bounding box to the mask's outer edge
    buffer <- if (!is.null(data[["buffer"]]) && !identical(data$mask_kind, "csv")) data[["buffer"]] else
      min(min(traps$x) - min(mask$x), max(mask$x) - max(traps$x),
          min(traps$y) - min(mask$y), max(mask$y) - max(traps$y)) + spacing / 2
  }
  # Optional trimming of the state space (settings buffer): used to show what a
  # too-small buffer does (negative control), by the buffer check, or to trim a large mask.
  if (!is.null(settings[["buffer"]])) {
    keep <- scr_dist_to_traps(mask, traps) <= settings[["buffer"]]
    mask <- mask[keep, , drop = FALSE]
    buffer <- settings[["buffer"]]
    notes <- c(notes, sprintf("State space trimmed to %s units from the nearest trap (settings: buffer).", format(settings[["buffer"]])))
  }
  if (length(dcov)) {
    miss <- setdiff(dcov, names(mask))
    if (length(miss)) stop("habitat mask has no column for density covariate(s): ", paste(miss, collapse = ", "))
    bad <- rowSums(is.na(mask[, dcov, drop = FALSE])) > 0
    if (any(bad)) {
      mask <- mask[!bad, , drop = FALSE]
      notes <- c(notes, sprintf("%d habitat pixels with a missing density covariate were treated as non-habitat.", sum(bad)))
    }
  }
  G <- nrow(mask)
  cscale <- scr_cscale(data$units, traps)
  ox <- min(mask$x) - spacing / 2; oy <- min(mask$y) - spacing / 2
  mxy <- cbind((mask$x - ox) / cscale, (mask$y - oy) / cscale)
  txy <- cbind((traps$x - ox) / cscale, (traps$y - oy) / cscale)
  d2 <- outer(mxy[, 1], txy[, 1], "-")^2 + outer(mxy[, 2], txy[, 2], "-")^2
  pix <- spacing / cscale
  # ---- covariates
  dsc <- list(); XD <- matrix(0, G, 0)
  for (v in dcov) { s <- scr_std(mask[[v]]); dsc[[v]] <- s[c("center", "scale", "binary")]; XD <- cbind(XD, s$x) }
  colnames(XD) <- dcov
  tcov <- function(v, what) {
    if (!v %in% names(traps)) stop(what, " covariate '", v, "' is not a column of the trap file")
    if (anyNA(traps[[v]])) stop(what, " covariate '", v, "' has missing values in the trap file")
    scr_std(as.numeric(traps[[v]]))
  }
  psc <- list(); XP <- matrix(0, J, 0)
  for (v in pcov) { s <- tcov(v, "detection"); psc[[v]] <- s[c("center", "scale", "binary")]; XP <- cbind(XP, s$x) }
  colnames(XP) <- pcov
  osc <- list(); XO <- matrix(0, K, 0)
  if (length(ocov)) {
    oc <- data$occasion_covariates
    if (is.null(oc) || nrow(oc) != K) stop("occasion-level detection covariates need an occasion table with one row per occasion (", K, ")")
    for (v in ocov) {
      if (!v %in% names(oc)) stop("occasion covariate '", v, "' is not a column of the occasion table")
      if (anyNA(oc[[v]])) stop("occasion covariate '", v, "' has missing values")
      s <- scr_std(as.numeric(oc[[v]])); osc[[v]] <- s[c("center", "scale", "binary")]; XO <- cbind(XO, s$x)
    }
  }
  colnames(XO) <- ocov
  ssc <- list(); XS <- matrix(0, J, 0)
  for (v in scov) { s <- tcov(v, "sigma"); ssc[[v]] <- s[c("center", "scale", "binary")]; XS <- cbind(XS, s$x) }
  colnames(XS) <- scov
  # ---- occasion structure: aggregate over occasions when nothing varies by occasion
  needK <- length(ocov) > 0 || behav == "trap" || (obs == "bernoulli" && !binary_usage)
  if (needK) {
    Kq <- K; use <- usage; obs_code <- if (obs == "poisson") 3L else 2L
    if (obs == "binomial") notes <- c(notes, "Detection varies by occasion, so the binomial model was fitted occasion by occasion (Bernoulli; identical likelihood).")
  } else {
    Kq <- 1L; use <- matrix(rowSums(usage), J, 1); obs_code <- if (obs == "poisson") 3L else 1L
  }
  XPq <- matrix(0, J * Kq, length(pcov) + length(ocov))
  for (k in seq_len(Kq)) {
    rows <- (k - 1) * J + seq_len(J)
    if (length(pcov)) XPq[rows, seq_along(pcov)] <- XP
    if (length(ocov)) XPq[rows, length(pcov) + seq_along(ocov)] <- matrix(XO[k, ], J, length(ocov), byrow = TRUE)
  }
  colnames(XPq) <- c(pcov, ocov)
  # ---- sex
  S <- if (length(sexeff)) 2L else 1L
  sexi <- rep(-1, n)
  if (S == 2) {
    ind <- data$individuals
    if (!is.null(ind)) {
      sx <- scr_sex_code(ind$sex[match(cn$ids, as.character(ind$individual))])
      sexi[sx %in% "F"] <- 0; sexi[sx %in% "M"] <- 1
    }
    if (all(sexi < 0)) notes <- c(notes, "WARNING: sex-specific detection requested but no detected individual has a known sex; the sexes are then a latent mixture and weakly identified.")
  }
  # ---- detected individual x trap pairs (sorted by individual), per-occasion counts and behaviour state
  pr <- which(y > 0, arr.ind = TRUE)
  pr <- pr[order(pr[, 1], pr[, 2]), , drop = FALSE]
  np <- nrow(pr)
  yp <- matrix(0, np, Kq); bp <- matrix(0, np, Kq)
  for (p in seq_len(np)) {
    yk <- y3[pr[p, 1], pr[p, 2], ]
    if (Kq == 1) yp[p, 1] <- sum(yk) else {
      yp[p, ] <- yk
      if (behav == "trap") { f <- min(which(yk > 0)); if (f < K) bp[p, (f + 1):K] <- 1 * (usage[pr[p, 2], (f + 1):K] > 0) }
    }
  }
  pstart <- match(seq_len(n), pr[, 1]); pend <- vapply(seq_len(n), function(i) max(which(pr[, 1] == i)), 1L)
  # ---- parameter layout of the coefficient vector beta
  PD <- length(dcov); PP <- ncol(XPq); PS <- length(scov)
  nb <- PD + PP + PS
  iSP <- if ("p0" %in% sexeff) (nb <- nb + 1) else 0
  iSS <- if ("sigma" %in% sexeff) (nb <- nb + 1) else 0
  iB <- if (behav == "trap") (nb <- nb + 1) else 0
  bidx <- c(0, PD, PD, PP, PD + PP, PS, iSP, iSS, iB)
  bnames <- c(if (PD) paste0("density_beta_", dcov), if (PP) paste0("detection_beta_", c(pcov, ocov)),
              if (PS) paste0("sigma_beta_", scov), if (iSP) "detection_beta_male", if (iSS) "sigma_beta_male",
              if (iB) "detection_beta_trap_response")
  # ---- data augmentation size
  M <- settings$M
  if (is.null(M)) M <- max(5 * n, n + 150)
  width <- max(diff(range(mxy[, 1])), diff(range(mxy[, 2]))) + pix
  sigma_max <- if (!is.null(settings$sigma_max)) settings$sigma_max / cscale else width / 2
  sp_traps <- scr_trap_spacing(traps)
  nz <- which(y > 0, arr.ind = TRUE); nz <- nz[order(nz[, 1], nz[, 2]), , drop = FALSE]
  list(n = n, J = J, K = K, G = G, M = M, y = y, y3 = y3, Kj = Kj, usage = usage, ids = cn$ids,
       captures = cn$captures, traps = traps, mask = mask, mxy = mxy, txy = txy, d2 = d2,
       pix = pix, area = G * pix^2, spacing = spacing, buffer = buffer, cscale = cscale, origin = c(ox, oy),
       units = if (is.null(data$units)) NA_character_ else data$units,
       obs = obs, obs_code = obs_code, Kq = Kq, use = use, binary_usage = binary_usage,
       XD = XD, XP = XP, XO = XO, XPq = XPq, XS = XS, dcov = dcov, pcov = pcov, ocov = ocov, scov = scov,
       dscale = dsc, pscale = psc, oscale = osc, sscale = ssc,
       S = S, sexeff = sexeff, sexi = sexi, behav = behav,
       pi = pr[, 1], pj = pr[, 2], yp = yp, bp = bp, pstart = pstart, pend = pend,
       nb = nb, bidx = bidx, bnames = bnames,
       nzi = nz[, 1], nzj = nz[, 2], nzc = y[nz],
       lam_max = if (is.null(settings$lambda0_max)) 10 else settings$lambda0_max,
       sigma_max = sigma_max, sigma_guess = sig_guess, trap_spacing = sp_traps,
       n_spatial_recap = sum(rowSums(y > 0) > 1), notes = notes)
}

# Is this a "basic" model (binomial/Poisson aggregated, detector-level p0
# covariates only)? The latent Nimble engine and NumPyro support only these.
scr_is_basic <- function(prep) prep$Kq == 1 && prep$S == 1 && !length(prep$scov) && prep$behav == "none"

scr_unit_word <- function(units) if (is.null(units) || is.na(units)) "coordinate" else units

# Names of the reported quantities (depend on whether map units are known)
scr_density_name <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "density_per_km2" else "density_per_unit2"
scr_area_label <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "km²" else "map unit²"
scr_det_name <- function(prep) if (prep$obs == "poisson") "baseline_encounter_rate" else "baseline_detection"
