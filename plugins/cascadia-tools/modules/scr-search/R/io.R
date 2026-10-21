# scr-search: reading the researcher's files, data checks, writing simulated
# data as the CSVs a field project would have.

ss_col <- function(df, name, what, file) {
  if (!name %in% names(df)) stop(sprintf("%s: column '%s' (%s) not found; columns are: %s", file, name, what,
                                         paste(names(df), collapse = ", ")), call. = FALSE)
  df[[name]]
}

# GPX track(s) -> data.frame(track_id, x, y) in the projected CRS `crs`
ss_read_gpx <- function(path, crs) {
  eco_require("sf")
  if (is.null(crs)) stop("GPX tracks are in longitude/latitude: give data: crs (e.g. 'EPSG:32610' for UTM zone 10N) so they can be projected to metres")
  pts <- sf::st_read(path, layer = "track_points", quiet = TRUE)
  pts <- sf::st_transform(pts, crs)
  xy <- sf::st_coordinates(pts)
  stem <- tools::file_path_sans_ext(basename(path))
  tid <- if (length(unique(pts$track_fid)) > 1) paste0(stem, "_", pts$track_fid) else stem
  data.frame(track_id = tid, x = xy[, 1], y = xy[, 2],
             time = if ("time" %in% names(pts)) as.character(pts$time) else NA_character_,
             order = seq_len(nrow(pts)), stringsAsFactors = FALSE)
}

ss_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc)) stop("config needs a data: block (tracks, samples, ...)")
  g <- function(k, d) if (is.null(dc[[k]])) d else dc[[k]]
  rd <- function(f) utils::read.csv(file.path(config$.dir, f), check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
  notes <- character()
  # ---- tracks: one CSV of GPS fixes, or one/several GPX files
  tf <- unlist(dc$tracks)
  if (all(grepl("\\.gpx$", tf, ignore.case = TRUE))) {
    tracks <- do.call(rbind, lapply(tf, function(f) ss_read_gpx(file.path(config$.dir, f), dc$crs)))
    tracks$occasion <- 1L
  } else {
    if (length(tf) != 1) stop("data: tracks must be one CSV file or a list of GPX files")
    tdf <- rd(tf)
    tracks <- data.frame(track_id = as.character(ss_col(tdf, g("track_id_col", "track_id"), "track ID", tf)),
                         x = as.numeric(ss_col(tdf, g("x_col", "x"), "x coordinate", tf)),
                         y = as.numeric(ss_col(tdf, g("y_col", "y"), "y coordinate", tf)), stringsAsFactors = FALSE)
    oc <- g("occasion_col", "occasion")
    tracks$occasion <- if (oc %in% names(tdf)) suppressWarnings(as.integer(tdf[[oc]])) else 1L
    for (v in names(tdf)) if (!v %in% names(tracks)) tracks[[v]] <- tdf[[v]]
    if (!is.null(dc$order_col)) {   # e.g. a GPS time stamp: sort fixes within each track
      o <- order(factor(tracks$track_id, levels = unique(tracks$track_id)), ss_col(tdf, dc$order_col, "fix order / time", tf))
      tracks <- tracks[o, ]
    }
  }
  # optional per-track table (occasion, team, snow...) joined by track ID
  if (!is.null(dc$track_info)) {
    ti <- rd(dc$track_info)
    tic <- ss_col(ti, g("track_id_col", "track_id"), "track ID", dc$track_info)
    m <- match(tracks$track_id, as.character(tic))
    if (anyNA(m)) stop(sprintf("%s: %d tracks have no row in the track_info file (e.g. %s)", dc$track_info,
                               length(unique(tracks$track_id[is.na(m)])), tracks$track_id[is.na(m)][1]), call. = FALSE)
    for (v in setdiff(names(ti), g("track_id_col", "track_id"))) tracks[[v]] <- ti[[v]][m]
    oc <- g("occasion_col", "occasion"); if (oc %in% names(ti)) tracks$occasion <- as.integer(ti[[oc]][m])
  }
  # ---- samples (genotyped; failed genotypes have an empty individual)
  sf_ <- dc$samples
  sdf <- rd(sf_)
  ind <- as.character(ss_col(sdf, g("individual_col", "individual"), "individual (genotype) ID", sf_))
  ind <- trimws(ind)
  failed <- is.na(ind) | ind %in% c("", "NA", "fail", "failed", "FAILED", "unknown")
  samples <- data.frame(sample_id = if (!is.null(dc$sample_id_col)) as.character(ss_col(sdf, dc$sample_id_col, "sample ID", sf_)) else paste0("S", seq_len(nrow(sdf))),
                        individual = ind, x = as.numeric(ss_col(sdf, g("sample_x_col", g("x_col", "x")), "sample x", sf_)),
                        y = as.numeric(ss_col(sdf, g("sample_y_col", g("y_col", "y")), "sample y", sf_)), stringsAsFactors = FALSE)
  stc <- g("sample_track_col", g("track_id_col", "track_id"))
  samples$track_id <- if (stc %in% names(sdf)) as.character(sdf[[stc]]) else NA_character_
  soc <- g("sample_occasion_col", g("occasion_col", "occasion"))
  samples$occasion <- if (soc %in% names(sdf)) suppressWarnings(as.integer(sdf[[soc]])) else NA_integer_
  if (any(failed)) notes <- c(notes, sprintf("%d of %d samples had no individual ID (failed or non-target genotypes) and were dropped; they carry no information on individuals but WERE part of the search effort, which is already counted from the tracks.",
                                             sum(failed), nrow(samples)))
  samples <- samples[!failed, , drop = FALSE]
  # ---- habitat mask
  mask <- NULL
  if (!is.null(dc$mask)) {
    mdf <- rd(dc$mask)
    mask <- data.frame(x = as.numeric(ss_col(mdf, g("mask_x_col", "x"), "x", dc$mask)),
                       y = as.numeric(ss_col(mdf, g("mask_y_col", "y"), "y", dc$mask)))
    for (v in unique(c(unlist(dc$density_covariates), intersect(unlist(dc$detection_covariates), names(mdf)))))
      mask[[v]] <- as.numeric(ss_col(mdf, v, "covariate", dc$mask))
    if (!is.null(dc$habitat_col)) { h <- ss_col(mdf, dc$habitat_col, "habitat flag", dc$mask); mask <- mask[h %in% c(1, TRUE), , drop = FALSE] }
  }
  data <- list(tracks = tracks, samples = samples, n_occasions = if (!is.null(dc$n_occasions)) dc$n_occasions else max(tracks$occasion, na.rm = TRUE),
               mask = mask, mask_spacing = dc$mask_spacing, buffer = dc$buffer, cell_size = dc$cell_size,
               units = if (is.null(dc$units)) NA_character_ else dc$units,
               density_covariates = as.character(unlist(dc$density_covariates)),
               detection_covariates = as.character(unlist(dc$detection_covariates)),
               covariate_labels = dc$labels, read_notes = notes)
  if (!is.null(dc$raster)) {
    eco_require("terra")
    r <- terra::rast(file.path(config$.dir, dc$raster))
    if (is.null(data$mask)) {
      sg <- ss_sigma_guess(samples); cs <- if (!is.null(dc$cell_size)) dc$cell_size else 1
      buf <- if (!is.null(dc$buffer)) dc$buffer else 4 * (if (is.finite(sg)) sg else 2 * cs)
      sp <- if (!is.null(dc$mask_spacing)) dc$mask_spacing else cs
      data$mask <- ss_make_mask(tracks, buf, sp); data$mask_spacing <- sp; data$buffer <- buf
    }
    v <- terra::extract(r, as.matrix(data$mask[, c("x", "y")]))
    for (k in names(v)) if (k %in% c(data$density_covariates, data$detection_covariates)) data$mask[[k]] <- v[[k]]
  }
  data
}

ss_check_data <- function(data) {
  pr <- character()
  tr <- data$tracks; sm <- data$samples
  if (!nrow(tr)) return("track file has no rows")
  if (anyNA(tr$x) || anyNA(tr$y)) {
    bad <- which(is.na(tr$x) | is.na(tr$y))
    pr <- c(pr, sprintf("track file: %d GPS fixes have missing or non-numeric coordinates (e.g. track %s, fix %d of that track); delete those rows or correct the coordinates",
                        length(bad), tr$track_id[bad[1]], sum(tr$track_id[seq_len(bad[1])] == tr$track_id[bad[1]])))
  }
  if (anyNA(tr$occasion)) pr <- c(pr, sprintf("track file / field sheet: %d tracks have no occasion (session) number, e.g. track %s; give whole numbers 1, 2, ...",
                                              length(unique(tr$track_id[is.na(tr$occasion)])), tr$track_id[is.na(tr$occasion)][1]))
  np <- table(tr$track_id); if (any(np < 2)) pr <- c(pr, sprintf("%d tracks have a single GPS fix (no length searched), e.g. %s; remove them or add fixes", sum(np < 2), names(np)[np < 2][1]))
  occ_per <- tapply(tr$occasion, tr$track_id, function(v) length(unique(v)))
  if (any(occ_per > 1, na.rm = TRUE)) pr <- c(pr, sprintf("%d tracks span several occasions (e.g. %s); give each search bout its own track ID", sum(occ_per > 1, na.rm = TRUE), names(occ_per)[which(occ_per > 1)[1]]))
  # big jumps between fixes = two tracks merged under one ID, or a coordinate typo
  jumps <- unlist(lapply(split(tr, tr$track_id), function(d) if (nrow(d) > 1) sqrt(diff(d$x)^2 + diff(d$y)^2) else numeric()))
  jumps <- jumps[is.finite(jumps)]
  cs <- data$cell_size
  if (!is.null(cs) && length(jumps) && any(jumps > 5 * cs)) pr <- c(pr, sprintf("track file: %d jumps between consecutive GPS fixes exceed 5 detector cells (%s); two tracks share an ID or a coordinate is wrong",
                                                                               sum(jumps > 5 * cs), format(5 * cs)))
  if (!nrow(sm)) pr <- c(pr, "sample file has no genotyped samples")
  if (anyNA(sm$x) || anyNA(sm$y)) pr <- c(pr, sprintf("sample file: %d samples have missing coordinates (e.g. sample %s); add the GPS location or remove the sample",
                                                       sum(is.na(sm$x) | is.na(sm$y)), sm$sample_id[which(is.na(sm$x) | is.na(sm$y))[1]]))
  if (nrow(sm)) {   # genotype IDs that differ only in case / spacing are almost always typos (they create ghost individuals)
    key <- toupper(gsub("[[:space:]_-]", "", sm$individual)); ids <- unique(sm$individual)
    clash <- tapply(sm$individual, key, function(v) length(unique(v)) > 1)
    if (any(clash)) { k <- names(clash)[clash][1]
      pr <- c(pr, sprintf("sample file: %d individual IDs differ only in upper/lower case, spaces or dashes (e.g. %s); each spelling counts as a different animal, so make them identical",
                          sum(clash), paste(unique(sm$individual[key == k]), collapse = " vs "))) }
  }
  if (any(!is.na(sm$track_id))) {
    unk <- setdiff(unique(sm$track_id[!is.na(sm$track_id) & sm$track_id != ""]), unique(tr$track_id))
    if (length(unk)) pr <- c(pr, sprintf("sample file: %d track IDs are not in the track file (e.g. %s); check spelling", length(unk), paste(utils::head(unk, 3), collapse = ", ")))
  }
  okt <- !is.na(tr$x) & !is.na(tr$y); oks <- !is.na(sm$x) & !is.na(sm$y)
  if (any(oks) && any(okt)) {
    sm <- sm[oks, , drop = FALSE]
    dd <- ss_min_dist(sm$x, sm$y, tr$x[okt], tr$y[okt])
    lim <- if (!is.null(cs)) cs else stats::quantile(dd, 0.99) * 2
    own <- !is.na(sm$track_id) & sm$track_id %in% tr$track_id
    if (any(own) && !is.null(cs)) {
      down <- vapply(which(own), function(r) { k <- okt & tr$track_id == sm$track_id[r]; min(sqrt((tr$x[k] - sm$x[r])^2 + (tr$y[k] - sm$y[r])^2)) }, 1)
      if (any(down > cs)) pr <- c(pr, sprintf("sample file: %d samples lie more than one detector cell (%s) from the track recorded for them (e.g. sample %s on track %s, %s away); check the sample coordinates or the track ID",
                                             sum(down > cs), format(cs), sm$sample_id[which(own)[which.max(down)]], sm$track_id[which(own)[which.max(down)]], format(signif(max(down), 3))))
    }
    if (any(dd > lim)) pr <- c(pr, sprintf("%d samples lie more than %s map units from any search track (e.g. sample %s, %s away); samples must come from the recorded searches (opportunistic samples need their own effort or must be excluded)",
                                            sum(dd > lim), format(signif(lim, 3)), sm$sample_id[which.max(dd)], format(signif(max(dd), 3))))
  }
  for (v in data$detection_covariates) {
    inT <- v %in% names(tr); inM <- !is.null(data$mask) && v %in% names(data$mask)
    if (!inT && !inM) pr <- c(pr, sprintf("detection covariate '%s' is not a column of the track file (or track_info) or the habitat mask", v))
    if (inT && anyNA(tr[[v]])) pr <- c(pr, sprintf("track covariate '%s' has missing values on %d tracks (e.g. track %s); fill them in on the field sheet or drop the covariate",
                                                   v, length(unique(tr$track_id[is.na(tr[[v]])])), tr$track_id[is.na(tr[[v]])][1]))
  }
  if (!is.null(data$mask)) {
    m <- data$mask
    if (nrow(m) < 100) pr <- c(pr, sprintf("habitat mask has only %d pixels; use a finer grid (aim for 1000-3000)", nrow(m)))
    for (v in data$density_covariates) if (!v %in% names(m)) pr <- c(pr, sprintf("habitat mask has no column '%s'", v))
    sp <- tryCatch(ss_infer_spacing(m), error = function(e) NA)
    if (is.finite(sp)) {
      dt <- ss_min_dist(tr$x[okt], tr$y[okt], m$x, m$y)
      if (any(dt > sp)) pr <- c(pr, sprintf("%d GPS fixes of the tracks lie outside the habitat mask; the mask must cover all searched ground plus a buffer", sum(dt > sp)))
    }
  } else if (length(data$density_covariates)) pr <- c(pr, "density covariates need a habitat mask CSV (data: mask) or a raster (data: raster)")
  if (!is.na(data$units) && !data$units %in% c("m", "km")) pr <- c(pr, "data: units must be m or km (or omitted if unknown)")
  if (is.null(cs)) pr <- c(pr, "data: cell_size (side of the square detector cells that discretise the search tracks, map units) is required; about sigma/2 to sigma (e.g. 1000 m for a mid-sized carnivore)")
  pr
}

ss_write_data <- function(data, out) {
  utils::write.csv(data$tracks, file.path(out, "search_tracks.csv"), row.names = FALSE)
  utils::write.csv(data$samples, file.path(out, "genetic_samples.csv"), row.names = FALSE)
  if (!is.null(data$mask)) utils::write.csv(data$mask, file.path(out, "habitat_grid.csv"), row.names = FALSE)
  invisible(out)
}
