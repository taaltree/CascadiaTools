# scr-smr: reading the researcher's files, data checks, writing simulated data.

smr_col <- function(df, name, what, file) {
  if (!name %in% names(df)) stop(sprintf("%s: column '%s' (%s) not found; columns are: %s", file, name, what,
                                         paste(names(df), collapse = ", ")), call. = FALSE)
  df[[name]]
}

# usage strings ("1101") -> 0/1 matrix
smr_usage <- function(u, K, ids, file) {
  u <- gsub("[^01]", "", as.character(u))
  if (is.null(K)) K <- max(nchar(u))
  if (any(nchar(u) != K)) stop(sprintf("%s: usage strings must have one 0/1 character per occasion (%d); row(s) %s differ",
                                       file, K, paste(utils::head(ids[nchar(u) != K], 5), collapse = ", ")), call. = FALSE)
  um <- t(vapply(u, function(s) as.integer(strsplit(s, "")[[1]])[seq_len(K)], integer(K)))
  rownames(um) <- ids
  um
}

smr_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc)) stop("config needs a data: block (cameras, sightings, n_occasions, marked or marking_captures, ...)")
  g <- function(k, d) if (is.null(dc[[k]])) d else dc[[k]]
  rd <- function(f, cc = NA) {
    fp <- file.path(config$.dir, f)
    if (!file.exists(fp)) stop(sprintf("data file '%s' not found in %s", f, config$.dir), call. = FALSE)
    # only force column classes for columns that exist (a misnamed column then gets the clear smr_col() message)
    if (!all(is.na(cc))) { hdr <- names(utils::read.csv(fp, nrows = 1, check.names = FALSE)); cc <- cc[names(cc) %in% hdr]; if (!length(cc)) cc <- NA }
    utils::read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE, colClasses = cc, na.strings = c("", "NA"))
  }
  # ---- cameras
  uc <- dc[["usage_col"]]
  cdf <- rd(dc[["cameras"]], if (!is.null(uc)) stats::setNames("character", uc) else NA)
  cams <- data.frame(id = as.character(smr_col(cdf, g("camera_id_col", "camera_id"), "camera ID", dc[["cameras"]])),
                     x = as.numeric(smr_col(cdf, g("x_col", "x"), "x coordinate", dc[["cameras"]])),
                     y = as.numeric(smr_col(cdf, g("y_col", "y"), "y coordinate", dc[["cameras"]])))
  pcov <- as.character(unlist(dc[["detection_covariates"]]))
  for (v in pcov) cams[[v]] <- smr_col(cdf, v, "detection covariate", dc[["cameras"]])
  K <- dc[["n_occasions"]]
  if (!is.null(uc)) {
    usage <- smr_usage(smr_col(cdf, uc, "usage (0/1 per occasion)", dc[["cameras"]]), K, cams$id, dc[["cameras"]])
    K <- ncol(usage)
  } else {
    if (is.null(K)) stop("data: n_occasions (number of resighting occasions) is required unless a usage column is given")
    usage <- matrix(1L, nrow(cams), K, dimnames = list(cams$id, NULL))
  }
  # ---- sightings (one row per sighting event, or with a count column)
  sdf <- rd(dc[["sightings"]])
  st_raw <- as.character(smr_col(sdf, g("status_col", "status"), "marked / unmarked / marked_unidentified", dc[["sightings"]]))
  mv <- unlist(g("marked_values", c("marked", "M", "marked_identified")))
  uv <- unlist(g("unmarked_values", c("unmarked", "U")))
  nv <- unlist(g("unidentified_values", c("marked_unidentified", "marked_unknown_id", "MU")))
  status <- ifelse(st_raw %in% mv, "marked", ifelse(st_raw %in% uv, "unmarked", ifelse(st_raw %in% nv, "marked_unidentified", NA)))
  si <- data.frame(camera = as.character(smr_col(sdf, g("camera_col", "camera_id"), "camera ID", dc[["sightings"]])),
                   occasion = suppressWarnings(as.integer(smr_col(sdf, g("occasion_col", "occasion"), "occasion number", dc[["sightings"]]))),
                   status = status, raw_status = st_raw,
                   animal_id = as.character(if (g("animal_col", "animal_id") %in% names(sdf)) sdf[[g("animal_col", "animal_id")]] else NA))
  si$animal_id[si$status != "marked" & !is.na(si$status)] <- NA
  if (!is.null(dc[["count_col"]])) {
    cn <- suppressWarnings(as.integer(smr_col(sdf, dc[["count_col"]], "number of animals in the sighting", dc[["sightings"]])))
    if (anyNA(cn) || any(cn < 1)) stop(dc[["sightings"]], ": count column must be whole numbers >= 1", call. = FALSE)
    si <- si[rep(seq_len(nrow(si)), cn), ]
  }
  rownames(si) <- NULL
  # ---- marked animals (SMR) / marking phase (gSMR)
  marked <- NULL
  if (!is.null(dc[["marked"]])) {
    mdf <- rd(dc[["marked"]])
    marked <- data.frame(id = as.character(smr_col(mdf, g("marked_id_col", "animal_id"), "marked animal ID", dc[["marked"]])))
    if (!is.null(dc[["marked_in_area_col"]])) {   # optional 0/1: collar functioning / animal present during resighting
      keep <- smr_col(mdf, dc[["marked_in_area_col"]], "marked animal available (0/1)", dc[["marked"]]) %in% c(1, TRUE)
      marked <- marked[keep, , drop = FALSE]
    }
  }
  marking <- NULL
  if (!is.null(dc[["marking_captures"]])) {
    if (is.null(dc[["marking_traps"]])) stop("data: marking_captures needs marking_traps (CSV of marking trap locations)")
    muc <- dc[["marking_usage_col"]]
    tdf <- rd(dc[["marking_traps"]], if (!is.null(muc)) stats::setNames("character", muc) else NA)
    mtr <- data.frame(id = as.character(smr_col(tdf, g("marking_trap_id_col", "trap_id"), "marking trap ID", dc[["marking_traps"]])),
                      x = as.numeric(smr_col(tdf, g("marking_x_col", g("x_col", "x")), "x coordinate", dc[["marking_traps"]])),
                      y = as.numeric(smr_col(tdf, g("marking_y_col", g("y_col", "y")), "y coordinate", dc[["marking_traps"]])))
    Km <- dc[["marking_occasions"]]
    if (!is.null(muc)) { musage <- smr_usage(smr_col(tdf, muc, "usage", dc[["marking_traps"]]), Km, mtr$id, dc[["marking_traps"]]); Km <- ncol(musage) }
    else {
      if (is.null(Km)) stop("data: marking_occasions is required unless marking_usage_col is given")
      musage <- matrix(1L, nrow(mtr), Km, dimnames = list(mtr$id, NULL))
    }
    cdf2 <- rd(dc[["marking_captures"]])
    cap <- data.frame(individual = as.character(smr_col(cdf2, g("marking_animal_col", "animal_id"), "animal ID", dc[["marking_captures"]])),
                      occasion = suppressWarnings(as.integer(smr_col(cdf2, g("marking_occasion_col", "occasion"), "marking occasion", dc[["marking_captures"]]))),
                      trap = as.character(smr_col(cdf2, g("marking_trap_col", "trap_id"), "marking trap ID", dc[["marking_captures"]])))
    marking <- list(traps = mtr, usage = musage, captures = cap, n_occasions = Km)
  }
  tel <- NULL
  if (!is.null(dc[["telemetry"]])) {
    tdf <- rd(dc[["telemetry"]])
    tel <- data.frame(animal_id = as.character(smr_col(tdf, g("telemetry_animal_col", "animal_id"), "animal ID", dc[["telemetry"]])),
                      x = as.numeric(smr_col(tdf, g("telemetry_x_col", "x"), "x", dc[["telemetry"]])),
                      y = as.numeric(smr_col(tdf, g("telemetry_y_col", "y"), "y", dc[["telemetry"]])))
  }
  mask <- NULL
  if (!is.null(dc[["mask"]])) {
    mdf <- rd(dc[["mask"]])
    mask <- data.frame(x = as.numeric(smr_col(mdf, g("mask_x_col", "x"), "x", dc[["mask"]])),
                       y = as.numeric(smr_col(mdf, g("mask_y_col", "y"), "y", dc[["mask"]])))
    for (v in unlist(dc[["density_covariates"]])) mask[[v]] <- as.numeric(smr_col(mdf, v, "density covariate", dc[["mask"]]))
    if (!is.null(dc[["habitat_col"]])) {
      h <- smr_col(mdf, dc[["habitat_col"]], "habitat flag", dc[["mask"]])
      mask <- mask[h %in% c(1, TRUE), , drop = FALSE]
    }
  }
  uir <- dc[["unidentified_recorded"]]
  if (is.null(uir)) uir <- any(si$status == "marked_unidentified", na.rm = TRUE)
  data <- list(cameras = cams, usage = usage, n_occasions = K, sightings = si, marked = marked,
               unidentified_recorded = isTRUE(uir), telemetry = tel, marking = marking,
               mask = mask, mask_spacing = dc[["mask_spacing"]], buffer = dc[["buffer"]],
               units = if (is.null(dc[["units"]])) NA_character_ else dc[["units"]],
               density_covariates = as.character(unlist(dc[["density_covariates"]])),
               detection_covariates = pcov, covariate_labels = dc[["labels"]],
               species = dc[["species"]])
  if (!is.null(dc[["raster"]])) {
    eco_require("terra")
    r <- terra::rast(file.path(config$.dir, dc[["raster"]]))
    if (is.null(data$mask)) {
      sp <- if (!is.null(dc[["mask_spacing"]])) dc[["mask_spacing"]] else smr_trap_spacing(cams) / 3
      buf <- if (!is.null(dc[["buffer"]])) dc[["buffer"]] else 3 * smr_trap_spacing(cams)
      data$mask <- smr_make_mask(rbind(cams[, c("x", "y")], if (!is.null(marking)) marking$traps[, c("x", "y")]), buf, sp)
      data$mask_spacing <- sp; data$buffer <- buf
    }
    v <- terra::extract(r, as.matrix(data$mask[, c("x", "y")]))
    for (k in names(v)) if (k %in% data$density_covariates) data$mask[[k]] <- v[[k]]
    miss <- setdiff(data$density_covariates, names(data$mask))
    if (length(miss)) stop("raster has no layer named: ", paste(miss, collapse = ", "))
  }
  data
}

smr_check_data <- function(data) {
  pr <- character()
  cams <- data$cameras
  if (anyNA(cams$x) || anyNA(cams$y)) pr <- c(pr, sprintf("camera file: %d cameras have missing or non-numeric coordinates", sum(is.na(cams$x) | is.na(cams$y))))
  if (any(duplicated(cams$id))) pr <- c(pr, paste("camera file: duplicated camera IDs:", paste(unique(cams$id[duplicated(cams$id)]), collapse = ", ")))
  if (any(duplicated(cams[, c("x", "y")]))) pr <- c(pr, "camera file: two cameras share the same coordinates; merge them or check the coordinates")
  si <- data$sightings
  if (any(is.na(si$status))) {
    bad <- unique(si$raw_status[is.na(si$status)])
    pr <- c(pr, sprintf("sightings file: %d rows have an unrecognised status (%s); use marked / unmarked / marked_unidentified or set marked_values/unmarked_values/unidentified_values in the config",
                        sum(is.na(si$status)), paste(utils::head(bad, 5), collapse = ", ")))
  }
  unk <- setdiff(unique(si$camera), cams$id)
  if (length(unk)) pr <- c(pr, sprintf("sightings file: %d camera IDs are not in the camera file (e.g. %s); check spelling/case",
                                        length(unk), paste(utils::head(unk, 5), collapse = ", ")))
  K <- data$n_occasions
  if (anyNA(si$occasion)) pr <- c(pr, "sightings file: some occasions are missing or not whole numbers")
  badk <- si$occasion[!is.na(si$occasion) & (si$occasion < 1 | si$occasion > K)]
  if (length(badk)) pr <- c(pr, sprintf("sightings file: occasions must be 1..%d but found %s", K, paste(utils::head(unique(badk), 5), collapse = ", ")))
  if (!length(unk) && !anyNA(si$occasion) && !length(badk)) {
    off <- data$usage[cbind(match(si$camera, cams$id), si$occasion)] == 0
    if (any(off)) pr <- c(pr, sprintf("%d sightings are at cameras marked as not operating on that occasion (e.g. %s occasion %d); fix the usage column",
                                      sum(off), si$camera[off][1], si$occasion[off][1]))
  }
  idm <- si$animal_id[which(si$status == "marked")]
  if (any(is.na(idm) | idm == "")) pr <- c(pr, sprintf("sightings file: %d sightings have status 'marked' but no animal ID; if the mark could not be read use status 'marked_unidentified'",
                                                     sum(is.na(idm) | idm == "")))
  marked_ids <- if (!is.null(data$marking)) unique(as.character(data$marking$captures$individual)) else as.character(data$marked$id)
  if (!is.null(data$marking) && !is.null(data$marked)) marked_ids <- unique(c(marked_ids, as.character(data$marked$id)))
  if (is.null(data$marking) && is.null(data$marked))
    pr <- c(pr, "no marked-animal list: give data: marked (CSV of every animal carrying a mark during resighting, including those never photographed) or the marking-phase captures (marking_traps + marking_captures) for generalized SMR")
  if (length(marked_ids) && any(duplicated(marked_ids)) && is.null(data$marking)) pr <- c(pr, "marked-animal file: duplicated animal IDs")
  stray <- setdiff(unique(idm[!is.na(idm)]), marked_ids)
  if (length(stray)) pr <- c(pr, sprintf("sightings file: %d identified marked animals are not in the marked-animal list (e.g. %s). Every marked animal must be listed (SMR) or appear in the marking captures (gSMR).",
                                          length(stray), paste(utils::head(stray, 5), collapse = ", ")))
  if (!length(marked_ids)) pr <- c(pr, "there are no marked animals: spatial mark-resight needs at least a few marked animals (with none, density is only weakly identified; see Chandler & Royle 2013)")
  if (!is.null(data$marking)) {
    mk <- data$marking; mc <- mk$captures
    unk2 <- setdiff(unique(mc$trap), mk$traps$id)
    if (length(unk2)) pr <- c(pr, sprintf("marking captures: %d trap IDs not in the marking-trap file (e.g. %s)", length(unk2), paste(utils::head(unk2, 5), collapse = ", ")))
    if (anyNA(mc$occasion) || any(mc$occasion < 1 | mc$occasion > mk$n_occasions, na.rm = TRUE))
      pr <- c(pr, sprintf("marking captures: occasions must be whole numbers 1..%d", mk$n_occasions))
  }
  if (!is.null(data$telemetry)) {
    tel <- data$telemetry
    if (anyNA(tel$x) || anyNA(tel$y)) pr <- c(pr, "telemetry file: some fixes have missing coordinates")
    st <- setdiff(unique(tel$animal_id), marked_ids)
    if (length(st)) pr <- c(pr, sprintf("telemetry file: %d animals are not marked animals (e.g. %s); telemetry is only used for marked animals", length(st), paste(utils::head(st, 3), collapse = ", ")))
  }
  for (v in data$detection_covariates) if (anyNA(cams[[v]])) pr <- c(pr, sprintf("camera file: detection covariate '%s' has missing values", v))
  if (!is.null(data$mask)) {
    m <- data$mask
    if (nrow(m) < 50) pr <- c(pr, sprintf("habitat mask has only %d pixels; use a finer grid (aim for 800-2500)", nrow(m)))
    if (nrow(m) > 4000) pr <- c(pr, sprintf("habitat mask has %d pixels; the latent-centre model needs <= ~2500 (use a coarser grid, spacing about sigma/2)", nrow(m)))
    for (v in data$density_covariates) if (!v %in% names(m)) pr <- c(pr, sprintf("habitat mask has no column '%s'", v))
    sp <- tryCatch(smr_infer_spacing(m), error = function(e) NA)
    if (is.finite(sp)) {
      off <- abs((m$x - min(m$x)) / sp - round((m$x - min(m$x)) / sp)) > 0.01 | abs((m$y - min(m$y)) / sp - round((m$y - min(m$y)) / sp)) > 0.01
      if (any(off)) pr <- c(pr, sprintf("habitat mask: %d points are not on a regular grid of spacing %g", sum(off), sp))
      okc <- !is.na(cams$x) & !is.na(cams$y)
      dt <- smr_dist_to_traps(cams[okc, , drop = FALSE], data.frame(x = m$x, y = m$y))
      if (any(dt > sp)) pr <- c(pr, sprintf("%d cameras lie outside the habitat mask; the mask must cover the cameras plus a buffer", sum(dt > sp)))
    }
  } else if (length(data$density_covariates)) pr <- c(pr, "density covariates need a habitat mask CSV (data: mask) or a raster (data: raster)")
  if (!is.na(data$units) && !data$units %in% c("m", "km")) pr <- c(pr, "data: units must be m or km (or omitted if unknown)")
  pr
}

smr_write_data <- function(data, out) {
  cams <- data$cameras
  cams$usage <- apply(data$usage, 1, paste, collapse = "")
  utils::write.csv(cams, file.path(out, "smr_cameras.csv"), row.names = FALSE)
  utils::write.csv(data$sightings, file.path(out, "smr_sightings.csv"), row.names = FALSE)
  if (!is.null(data$marked)) utils::write.csv(data.frame(animal_id = data$marked$id), file.path(out, "smr_marked.csv"), row.names = FALSE)
  if (!is.null(data$telemetry)) utils::write.csv(data$telemetry, file.path(out, "smr_telemetry.csv"), row.names = FALSE)
  if (!is.null(data$marking)) {
    utils::write.csv(data.frame(trap_id = data$marking$traps$id, x = data$marking$traps$x, y = data$marking$traps$y),
                     file.path(out, "smr_marking_traps.csv"), row.names = FALSE)
    utils::write.csv(data.frame(animal_id = data$marking$captures$individual, occasion = data$marking$captures$occasion,
                                trap_id = data$marking$captures$trap), file.path(out, "smr_marking_captures.csv"), row.names = FALSE)
  }
  if (!is.null(data$mask)) utils::write.csv(data$mask, file.path(out, "smr_habitat.csv"), row.names = FALSE)
  invisible(out)
}
