# scr-smr: data structure, habitat mask and the internal representation
# shared by the engines. Geometry helpers are copied from modules/scr
# (prefix smr_ so both modules can be loaded in one session).
#
# Data structure (returned by simulate()$data and read_data()):
#   cameras     data.frame(id, x, y, <detection covariates>): resighting detectors
#   usage       J x K 0/1 matrix: camera j operating on resight occasion k
#   n_occasions K (resighting occasions)
#   sightings   data.frame(camera, occasion, status, animal_id): one row per
#               sighting event. status is "marked" (marked animal, identity
#               read; animal_id required), "unmarked" (animal without marks),
#               or "marked_unidentified" (marked but the mark could not be
#               read; secr's Tm). animal_id is NA unless status = "marked".
#   marked      data.frame(id): every animal carrying a mark during
#               resighting (SMR), including marked animals never resighted.
#               NULL for gSMR (marked = animals caught in the marking phase).
#   unidentified_recorded  TRUE if marked-but-unidentified sightings were
#               recorded (then the identification probability is estimated);
#               FALSE: every sighting of a marked animal was identified.
#   telemetry   optional data.frame(animal_id, x, y): telemetry fixes of
#               marked animals (same coordinate system).
#   marking     optional (gSMR): list(traps = data.frame(id, x, y),
#               usage = L x Km 0/1, captures = data.frame(individual,
#               occasion, trap), n_occasions = Km). Animals captured here are
#               the marked animals; the capture process is modelled jointly.
#   mask, mask_spacing, buffer, units, density_covariates,
#   detection_covariates, covariate_labels   as in modules/scr
#
# Internally coordinates are in km (units m or km) or in trap-spacing units.

smr_cscale <- function(units, traps) {
  if (identical(units, "m")) return(1000)
  if (identical(units, "km")) return(1)
  smr_trap_spacing(traps)
}

smr_trap_spacing <- function(traps) {
  xy <- as.matrix(traps[, c("x", "y")])
  if (nrow(xy) < 2) return(1)
  d <- as.matrix(stats::dist(xy)); diag(d) <- Inf
  stats::median(apply(d, 1, min))
}

smr_make_mask <- function(traps, buffer, spacing) {
  xr <- range(traps$x) + c(-buffer, buffer); yr <- range(traps$y) + c(-buffer, buffer)
  gx <- seq(xr[1] + spacing / 2, xr[2], by = spacing)
  gy <- seq(yr[1] + spacing / 2, yr[2], by = spacing)
  m <- expand.grid(x = gx, y = gy)
  keep <- smr_dist_to_traps(m, traps) <= buffer
  m[keep, , drop = FALSE]
}

smr_dist_to_traps <- function(pts, traps) {
  d2 <- outer(pts$x, traps$x, "-")^2 + outer(pts$y, traps$y, "-")^2
  sqrt(apply(d2, 1, min))
}

smr_infer_spacing <- function(mask) {
  dx <- diff(sort(unique(round(mask$x, 6)))); dy <- diff(sort(unique(round(mask$y, 6))))
  s <- c(dx[dx > 0], dy[dy > 0])
  if (!length(s)) stop("habitat mask needs at least two distinct x and y values")
  min(s)
}

smr_std <- function(v) {
  u <- unique(v[!is.na(v)])
  if (all(u %in% c(0, 1))) return(list(x = v, center = 0, scale = 1, binary = TRUE))
  ce <- mean(v, na.rm = TRUE); sc <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  list(x = (v - ce) / sc, center = ce, scale = sc, binary = FALSE)
}

smr_unit_word <- function(units) if (is.null(units) || is.na(units)) "coordinate" else units
smr_density_name <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "density_per_km2" else "density_per_unit2"
smr_area_label <- function(units) if (!is.null(units) && !is.na(units) && units %in% c("m", "km")) "km²" else "map unit²"
# internal area -> reported area (km2 when units known, else map unit2)
smr_area_conv <- function(prep) if (!is.na(prep$units) && prep$units %in% c("m", "km")) 1 else prep$cscale^2

# Pooled RPSV-type sigma guess from repeat locations of the same animal
# (identified sightings, marking captures, telemetry). Only sizes the mask.
smr_sigma_guess <- function(pts) {
  if (is.null(pts) || !nrow(pts)) return(NA_real_)
  sp <- split(seq_len(nrow(pts)), pts$id)
  ss <- 0; df <- 0
  for (ix in sp) if (length(ix) > 1) {
    mx <- mean(pts$x[ix]); my <- mean(pts$y[ix])
    ss <- ss + sum((pts$x[ix] - mx)^2 + (pts$y[ix] - my)^2); df <- df + length(ix) - 1
  }
  if (df < 1 || ss <= 0) return(NA_real_)
  sqrt(ss / (2 * df))
}

smr_model_type <- function(data, settings) {
  m <- if (is.null(settings$model)) "auto" else settings$model
  if (m == "auto") m <- if (!is.null(data$marking)) "gsmr" else "smr"
  if (!m %in% c("smr", "gsmr")) stop("settings$model must be auto, smr or gsmr")
  if (m == "gsmr" && is.null(data$marking)) stop("model gsmr needs marking-phase data (data: marking_traps and marking_captures)")
  m
}

smr_prepare <- function(data, settings = list()) {
  model <- smr_model_type(data, settings)
  pois <- !identical(settings$resight_model, "bernoulli")
  notes <- character()
  cams <- data$cameras; cams$id <- as.character(cams$id)
  J <- nrow(cams); K <- data$n_occasions
  usage <- if (is.null(data$usage)) matrix(1L, J, K) else as.matrix(data$usage)
  Kj <- rowSums(usage)
  si <- data$sightings
  si$camera <- as.character(si$camera); si$status <- as.character(si$status)
  si$animal_id <- as.character(si$animal_id)
  # ---- who is marked
  mk_cap <- NULL
  if (!is.null(data$marking)) {
    mk <- data$marking
    mk$traps$id <- as.character(mk$traps$id)
    mcap <- unique(data.frame(individual = as.character(mk$captures$individual),
                              occasion = as.integer(mk$captures$occasion), trap = as.character(mk$captures$trap)))
    mk_cap <- mcap
  }
  ids <- if (!is.null(mk_cap)) unique(mk_cap$individual) else as.character(data$marked$id)
  if (!is.null(mk_cap) && !is.null(data$marked)) ids <- unique(c(ids, as.character(data$marked$id)))
  m <- length(ids)
  # ---- identified sightings: m x J counts (Bernoulli: distinct occasions)
  sid <- si[si$status == "marked", , drop = FALSE]
  if (!pois) {
    dup <- duplicated(sid[, c("animal_id", "camera", "occasion")])
    if (any(dup)) notes <- c(notes, sprintf("Bernoulli resight model: %d repeat sightings of the same marked animal at the same camera on the same occasion were collapsed to one.", sum(dup)))
    sid <- sid[!dup, , drop = FALSE]
  }
  yID <- matrix(0, m, J, dimnames = list(ids, cams$id))
  if (nrow(sid)) {
    ij <- cbind(match(sid$animal_id, ids), match(sid$camera, cams$id))
    for (r in seq_len(nrow(ij))) yID[ij[r, 1], ij[r, 2]] <- yID[ij[r, 1], ij[r, 2]] + 1
  }
  cnt <- function(st) { v <- si$camera[si$status == st]; as.numeric(table(factor(v, levels = cams$id))) }
  nU <- cnt("unmarked"); nTm <- cnt("marked_unidentified")
  has_tm <- isTRUE(data$unidentified_recorded) || sum(nTm) > 0
  if (isFALSE(settings$use_unidentified)) has_tm <- FALSE
  if (!has_tm && sum(nTm) > 0) notes <- c(notes, "Marked-but-unidentified sightings were ignored (settings$use_unidentified = FALSE).")
  # ---- telemetry
  tel <- data$telemetry
  if (!is.null(tel) && nrow(tel)) {
    tel$animal_id <- as.character(tel$animal_id)
    tel <- tel[tel$animal_id %in% ids, , drop = FALSE]
  }
  has_tel <- !is.null(tel) && nrow(tel) > 0 && !isFALSE(settings$use_telemetry)
  # ---- mask
  pts <- rbind(
    if (nrow(sid)) data.frame(id = sid$animal_id, x = cams$x[match(sid$camera, cams$id)], y = cams$y[match(sid$camera, cams$id)]),
    if (!is.null(mk_cap)) data.frame(id = mk_cap$individual, x = data$marking$traps$x[match(mk_cap$trap, data$marking$traps$id)],
                                     y = data$marking$traps$y[match(mk_cap$trap, data$marking$traps$id)]))
  sig_guess <- smr_sigma_guess(pts)
  if (has_tel) { sg_t <- smr_sigma_guess(data.frame(id = tel$animal_id, x = tel$x, y = tel$y)); if (is.finite(sg_t)) sig_guess <- sg_t }
  dcov <- as.character(if (!is.null(settings$density_covariates)) unlist(settings$density_covariates) else unlist(data$density_covariates))
  pcov <- as.character(if (!is.null(settings$detection_covariates)) unlist(settings$detection_covariates) else unlist(data$detection_covariates))
  alltr <- rbind(cams[, c("x", "y")], if (!is.null(data$marking)) data$marking$traps[, c("x", "y")])
  mask <- data$mask
  if (is.null(mask)) {
    sp0 <- smr_trap_spacing(cams)
    sg <- if (is.finite(sig_guess) && sig_guess > 0) sig_guess else sp0
    buffer <- if (!is.null(data$buffer)) data$buffer else 4 * sg
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else sg / 2
    mask <- smr_make_mask(alltr, buffer, spacing)
    while (nrow(mask) > (if (is.null(settings$max_mask_points)) 2500 else settings$max_mask_points)) {
      spacing <- spacing * 1.25; mask <- smr_make_mask(alltr, buffer, spacing)
    }
    notes <- c(notes, sprintf("Habitat mask built from the detectors: buffer %s, pixel spacing %s (%s units), %d pixels.",
                              format(signif(buffer, 3)), format(signif(spacing, 3)), smr_unit_word(data$units), nrow(mask)))
    if (length(dcov)) stop("density covariates need a habitat mask file or raster with those columns: ", paste(dcov, collapse = ", "))
  } else {
    spacing <- if (!is.null(data$mask_spacing)) data$mask_spacing else smr_infer_spacing(mask)
    buffer <- min(min(alltr$x) - min(mask$x), max(mask$x) - max(alltr$x),
                  min(alltr$y) - min(mask$y), max(mask$y) - max(alltr$y)) + spacing / 2
  }
  if (length(dcov)) {
    bad <- rowSums(is.na(mask[, dcov, drop = FALSE])) > 0
    if (any(bad)) {
      mask <- mask[!bad, , drop = FALSE]
      notes <- c(notes, sprintf("%d habitat pixels with a missing density covariate were treated as non-habitat.", sum(bad)))
    }
  }
  G <- nrow(mask)
  cscale <- smr_cscale(data$units, cams)
  ox <- min(mask$x) - spacing / 2; oy <- min(mask$y) - spacing / 2
  tox <- function(x) (x - ox) / cscale; toy <- function(y) (y - oy) / cscale
  mxy <- cbind(tox(mask$x), toy(mask$y))
  txy <- cbind(tox(cams$x), toy(cams$y))
  pix <- spacing / cscale
  # ---- covariates
  dsc <- list(); XD <- matrix(0, G, 0)
  for (v in dcov) { s <- smr_std(mask[[v]]); dsc[[v]] <- s[c("center", "scale", "binary")]; XD <- cbind(XD, s$x) }
  colnames(XD) <- dcov
  psc <- list(); XP <- matrix(0, J, 0)
  for (v in pcov) {
    if (!v %in% names(cams)) stop("detection covariate '", v, "' is not a column of the camera file")
    s <- smr_std(cams[[v]]); psc[[v]] <- s[c("center", "scale", "binary")]; XP <- cbind(XP, s$x)
  }
  colnames(XP) <- pcov
  # ---- telemetry sufficient statistics (n fixes, mean x, mean y, sum of squares)
  telS <- matrix(0, m, 4)
  if (has_tel) {
    for (k in seq_len(m)) {
      w <- tel$animal_id == ids[k]
      if (any(w)) {
        tx <- tox(tel$x[w]); ty <- toy(tel$y[w])
        telS[k, ] <- c(sum(w), mean(tx), mean(ty), sum((tx - mean(tx))^2 + (ty - mean(ty))^2))
      }
    }
  }
  # ---- marking phase (gSMR)
  mkp <- NULL
  if (!is.null(data$marking)) {
    mk <- data$marking
    L <- nrow(mk$traps); Km <- mk$n_occasions
    musage <- if (is.null(mk$usage)) matrix(1L, L, Km) else as.matrix(mk$usage)
    xm <- matrix(0, m, L, dimnames = list(ids, mk$traps$id))
    ij <- cbind(match(mk_cap$individual, ids), match(mk_cap$trap, mk$traps$id))
    for (r in seq_len(nrow(ij))) xm[ij[r, 1], ij[r, 2]] <- xm[ij[r, 1], ij[r, 2]] + 1
    mkp <- list(L = L, Km = Km, Kmj = rowSums(musage), usage = musage, mtxy = cbind(tox(mk$traps$x), toy(mk$traps$y)),
                xm = xm, traps = mk$traps, captures = mk_cap, n_marked_captured = length(unique(mk_cap$individual)))
  }
  width <- max(diff(range(mxy[, 1])), diff(range(mxy[, 2]))) + pix
  sigma_max <- if (!is.null(settings$sigma_max)) settings$sigma_max / cscale else width / 2
  list(model = model, pois = pois, has_tm = has_tm, has_tel = has_tel, has_marking = model == "gsmr",
       J = J, K = K, Kj = Kj, usage = usage, cameras = cams, txy = txy, m = m, ids = ids,
       yID = yID, nU = nU, nTm = nTm, telS = telS, telemetry = if (has_tel) tel else NULL, mk = mkp,
       sightings = si, G = G, mask = mask, mxy = mxy, pix = pix, area = G * pix^2, spacing = spacing,
       buffer = buffer, cscale = cscale, origin = c(ox, oy),
       units = if (is.null(data$units)) NA_character_ else data$units,
       XD = XD, XP = XP, dcov = dcov, pcov = pcov, dscale = dsc, pscale = psc,
       sigma_max = sigma_max, sigma_guess = sig_guess, trap_spacing = smr_trap_spacing(cams),
       n_id_sightings = sum(yID), n_unmarked_sightings = sum(nU), n_unid_sightings = sum(nTm),
       n_marked_resighted = sum(rowSums(yID) > 0), n_spatial_resight = sum(rowSums(yID > 0) > 1),
       notes = notes)
}
