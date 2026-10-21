# scr: reading the researcher's files, data checks, writing simulated data.

scr_col <- function(df, name, what, file, key = NULL) {
  if (!name %in% names(df)) stop(sprintf("%s: column '%s' (%s) not found; columns are: %s.%s", file, name, what,
                                         paste(names(df), collapse = ", "),
                                         if (!is.null(key)) sprintf(" Rename the column, or set `%s:` in the data: block of config.yaml to the right column name.", key) else
                                           " Check the spelling in config.yaml."), call. = FALSE)
  df[[name]]
}

scr_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc)) stop("config needs a data: block (traps, captures, n_occasions, ...)")
  path <- function(f) file.path(config$.dir, f)
  rd <- function(f, cc = NA) utils::read.csv(path(f), check.names = FALSE, stringsAsFactors = FALSE,
                                             strip.white = TRUE, colClasses = cc, na.strings = c("NA", ""))
  g <- function(k, d) if (is.null(dc[[k]])) d else dc[[k]]
  # usage strings such as "01111111" must stay text (leading zeros)
  tdf <- rd(dc$traps, if (!is.null(dc$usage_col)) stats::setNames("character", dc$usage_col) else NA)
  traps <- data.frame(id = as.character(scr_col(tdf, g("trap_id_col", "id"), "trap ID", dc$traps, "trap_id_col")),
                      x = suppressWarnings(as.numeric(scr_col(tdf, g("x_col", "x"), "x coordinate", dc$traps, "x_col"))),
                      y = suppressWarnings(as.numeric(scr_col(tdf, g("y_col", "y"), "y coordinate", dc$traps, "y_col"))))
  pcov <- as.character(unlist(dc$detection_covariates)); scov <- as.character(unlist(dc$sigma_covariates))
  for (v in unique(c(pcov, scov))) traps[[v]] <- scr_col(tdf, v, "detector covariate", dc$traps)
  # trap operation / effort: usage string per trap ("110111"), an effort table, or all operating
  K <- dc$n_occasions
  uc <- dc$usage_col
  if (!is.null(dc[["effort"]])) {
    edf <- rd(dc[["effort"]])
    eid <- g("effort_trap_col", g("trap_id_col", "id"))
    if (all(c(eid, "occasion", "effort") %in% names(edf)) && ncol(edf) <= 4) {
      # long format: trap, occasion, effort
      if (is.null(K)) K <- max(edf$occasion, na.rm = TRUE)
      usage <- matrix(0, nrow(traps), K)
      ii <- match(as.character(edf[[eid]]), traps$id)
      if (anyNA(ii)) stop(sprintf("%s: detector(s) not in the trap file: %s", dc[["effort"]], paste(utils::head(unique(edf[[eid]][is.na(ii)]), 5), collapse = ", ")), call. = FALSE)
      usage[cbind(ii, as.integer(edf$occasion))] <- as.numeric(edf$effort)
    } else {
      # wide format: detector id + one column per occasion (in order)
      ids <- as.character(scr_col(edf, eid, "detector ID", dc[["effort"]], "effort_trap_col"))
      em0 <- as.matrix(edf[, setdiff(names(edf), eid), drop = FALSE])
      em <- em0; suppressWarnings(storage.mode(em) <- "double")
      txt <- vapply(seq_len(ncol(em)), function(k) all(is.na(em[, k])) && any(!is.na(em0[, k]) & nzchar(trimws(em0[, k]))), TRUE)
      if (any(txt)) stop(sprintf("%s: column(s) %s are not numbers. The effort table must hold the detector ID column ('%s') and one numeric column per occasion (days active or 0/1), nothing else; remove or rename the other columns.",
                                 dc[["effort"]], paste(sprintf("'%s'", colnames(em)[txt]), collapse = ", "), eid), call. = FALSE)
      if (is.null(K)) K <- ncol(em)
      if (ncol(em) != K) stop(sprintf("%s: %d occasion columns but n_occasions is %d", dc[["effort"]], ncol(em), K), call. = FALSE)
      ii <- match(traps$id, ids)
      if (anyNA(ii)) stop(sprintf("%s: no effort row for detector(s) %s", dc[["effort"]], paste(utils::head(traps$id[is.na(ii)], 5), collapse = ", ")), call. = FALSE)
      usage <- em[ii, , drop = FALSE]
    }
    if (!is.null(dc$effort_scale)) usage <- usage / dc$effort_scale   # e.g. days active / 7 nominal days
  } else if (!is.null(uc)) {
    u <- as.character(scr_col(tdf, uc, "usage (0/1 per occasion)", dc$traps))
    u <- gsub("[^01]", "", u)
    if (is.null(K)) K <- max(nchar(u))
    if (any(nchar(u) != K)) stop(sprintf("%s: usage strings must have one 0/1 character per occasion (%d); trap(s) %s differ",
                                         dc$traps, K, paste(utils::head(traps$id[nchar(u) != K], 5), collapse = ", ")), call. = FALSE)
    usage <- t(vapply(u, function(s) as.integer(strsplit(s, "")[[1]])[seq_len(K)], integer(K)))
  } else {
    if (is.null(K)) stop("data: n_occasions is required (number of sampling occasions) unless a usage column or effort table is given")
    usage <- matrix(1L, nrow(traps), K)
  }
  rownames(usage) <- traps$id
  # detections (camera/live-trap records) or genetic samples
  cfile <- if (!is.null(dc$samples)) dc$samples else dc$captures
  if (is.null(cfile)) stop("data: needs captures: (detections) or samples: (genetic samples) file")
  cdf <- rd(cfile)
  qc <- NULL
  if (!is.null(dc$qc_col)) {
    q <- scr_col(cdf, dc$qc_col, "sample QC result", cfile, "qc_col")
    pass <- as.character(unlist(g("qc_pass", c("pass", "PASS", "Pass", "ok", "OK", "1", "TRUE"))))
    keep <- as.character(q) %in% pass
    if (!any(keep)) {
      tq <- table(as.character(q), useNA = "ifany")
      stop(sprintf("%s: no sample passed QC. Column '%s' holds %s, but qc_pass is [%s] (matching is case-sensitive). Set `qc_pass:` in the data: block of config.yaml to the value(s) that mean a sample passed.",
                   cfile, dc$qc_col, paste(sprintf("'%s' (%d)", names(tq), as.integer(tq)), collapse = ", "), paste(pass, collapse = ", ")), call. = FALSE)
    }
    qc <- list(n_samples = nrow(cdf), n_pass = sum(keep), n_fail = sum(!keep), column = dc$qc_col)
    cdf <- cdf[keep, , drop = FALSE]
  }
  cap <- data.frame(individual = as.character(scr_col(cdf, g("individual_col", "individual"), "animal ID", cfile, "individual_col")),
                    occasion = suppressWarnings(as.integer(scr_col(cdf, g("occasion_col", "occasion"), "occasion number", cfile, "occasion_col"))),
                    trap = as.character(scr_col(cdf, g("trap_col", "trap"), "trap ID", cfile, "trap_col")))
  cap$row <- suppressWarnings(as.integer(rownames(cdf)))   # data row in the file (kept through QC filtering) for messages
  if (!is.null(dc$count_col)) cap$count <- suppressWarnings(as.numeric(scr_col(cdf, dc$count_col, "number of samples", cfile, "count_col")))
  # sex: individuals file or a column of the detection file
  ind <- NULL; sex_bad_calls <- NULL; read_notes <- character()
  if (!is.null(dc$individuals)) {
    idf <- rd(dc$individuals)
    ind <- data.frame(individual = as.character(scr_col(idf, g("individual_col", "individual"), "animal ID", dc$individuals, "individual_col")),
                      sex_raw = as.character(scr_col(idf, g("sex_col", "sex"), "sex", dc$individuals, "sex_col")))
  } else if (!is.null(dc$sex_col)) {
    sx <- data.frame(individual = cap$individual, sex_raw = as.character(scr_col(cdf, dc$sex_col, "sex", cfile, "sex_col")))
    sx$sex_raw[is.na(sx$sex_raw)] <- ""
    sx <- sx[!is.na(sx$individual) & nzchar(sx$individual), , drop = FALSE]
    # sample-level calls that are neither F/M nor an "unknown" code are reported by check_data()
    code <- scr_sex_code(sx$sex_raw)
    unrec <- nzchar(sx$sex_raw) & is.na(code) & !toupper(sx$sex_raw) %in% c("U", "UNK", "UNKNOWN", "?")
    if (any(unrec)) sex_bad_calls <- unique(sx$sex_raw[unrec])
    # one sex per animal: the most frequent recognised call (F/M); conflicting calls are noted
    ind <- do.call(rbind, lapply(split(code, sx$individual), function(v) {
      v <- v[!is.na(v)]
      data.frame(sex_raw = if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else NA_character_,
                 conflict = length(unique(v)) > 1)
    }))
    ind$individual <- rownames(ind)
    if (any(ind$conflict)) read_notes <- c(read_notes, sprintf(
      "Sex calls disagree among the samples of %d individual(s) (%s); the most frequent call was used. Conflicting sex calls can indicate a genotyping error or a mixed sample: check them in the lab data.",
      sum(ind$conflict), paste(utils::head(ind$individual[ind$conflict], 5), collapse = ", ")))
  }
  if (!is.null(ind)) { ind$sex <- scr_sex_code(ind$sex_raw); ind <- ind[!is.na(ind$individual), c("individual", "sex", "sex_raw")] }
  # occasion-level covariates
  occ <- NULL
  if (!is.null(dc$occasions)) {
    odf <- rd(dc$occasions)
    oc <- g("occasion_id_col", "occasion")
    o <- as.integer(scr_col(odf, oc, "occasion number", dc$occasions, "occasion_id_col"))
    occ <- odf[order(o), , drop = FALSE]; occ$occasion <- sort(o)
  }
  data <- list(traps = traps, usage = usage, captures = cap, n_occasions = K, individuals = ind,
               occasion_covariates = occ, mask = NULL, mask_kind = "traps",
               mask_spacing = dc$mask_spacing, buffer = dc[["buffer"]],
               units = if (is.null(dc$units)) NA_character_ else dc$units,
               density_covariates = as.character(unlist(dc$density_covariates)),
               detection_covariates = pcov, sigma_covariates = scov,
               occasion_detection_covariates = as.character(unlist(dc$occasion_covariates)),
               obs_model = if (is.null(dc$obs_model)) NULL else tolower(trimws(dc$obs_model)),
               behaviour = if (is.null(dc$behaviour)) NULL else tolower(trimws(dc$behaviour)),
               sex_effects = if (!is.null(dc$sex_effects)) as.character(unlist(dc$sex_effects)) else if (!is.null(ind)) c("p0", "sigma") else NULL,
               covariate_labels = dc$labels, qc = qc)
  if (!is.null(dc[["mask"]])) {
    mdf <- rd(dc[["mask"]])
    mask <- data.frame(x = as.numeric(scr_col(mdf, g("mask_x_col", "x"), "x", dc[["mask"]])),
                       y = as.numeric(scr_col(mdf, g("mask_y_col", "y"), "y", dc[["mask"]])))
    for (v in data$density_covariates) mask[[v]] <- as.numeric(scr_col(mdf, v, "density covariate", dc[["mask"]]))
    if (!is.null(dc$habitat_col)) {   # optional 0/1 habitat flag: drop non-habitat
      h <- scr_col(mdf, dc$habitat_col, "habitat flag", dc[["mask"]])
      mask <- mask[h %in% c(1, TRUE), , drop = FALSE]
    }
    data[["mask"]] <- mask; data$mask_kind <- "csv"
  } else if (!is.null(dc$habitat_polygon) || !is.null(dc$raster)) {
    # polygon (shapefile / GeoPackage) of suitable habitat and/or raster (GeoTIFF) of covariates [+ habitat layer]
    if (!is.null(dc$habitat_polygon)) {
      eco_require("sf")
      data$habitat_polygon <- sf::st_read(path(dc$habitat_polygon), quiet = TRUE)
      data$mask_kind <- "polygon"
    }
    if (!is.null(dc$raster)) {
      data$raster_file <- normalizePath(path(dc$raster))
      data$habitat_layer <- dc$habitat_layer
      if (is.null(dc$habitat_polygon)) data$mask_kind <- "raster"
    }
    # Build the mask from the valid rows only, so that a missing coordinate, a misspelled
    # detector ID or an out-of-range occasion is reported by check_data() with a clear
    # message instead of breaking the mask construction here.
    okt <- is.finite(traps$x) & is.finite(traps$y)
    if (!any(okt)) stop(sprintf("%s: no detector has numeric coordinates (columns '%s', '%s')", dc$traps, g("x_col", "x"), g("y_col", "y")), call. = FALSE)
    if (!is.na(data$units) && data$units %in% c("m", "km") && all(abs(traps$x[okt]) <= 180) && all(abs(traps$y[okt]) <= 90))
      stop(sprintf("%s: the detector coordinates look like longitude/latitude in degrees, but units is '%s'. Project the detectors (and the habitat polygon/raster) to a metric CRS such as UTM first.",
                   dc$traps, data$units), call. = FALSE)
    dsub <- data; dsub$traps <- traps[okt, , drop = FALSE]
    okc <- cap$trap %in% dsub$traps$id & !is.na(cap$occasion) & cap$occasion >= 1 & cap$occasion <= K &
      !is.na(cap$individual) & nzchar(cap$individual)
    dsub$captures <- cap[okc, , drop = FALSE]
    sg <- if (nrow(dsub$captures)) scr_sigma_guess(scr_counts(dsub)$captures, dsub$traps) else NA_real_
    if (!is.finite(sg)) sg <- scr_trap_spacing(dsub$traps)
    buf <- if (!is.null(dc[["buffer"]])) dc[["buffer"]] else 4 * sg
    sp <- if (!is.null(dc$mask_spacing)) dc$mask_spacing else sg / 2
    data[["mask"]] <- scr_build_mask(dsub, buf, sp); data$mask_spacing <- sp; data[["buffer"]] <- buf
    data$mask_notes <- attr(data[["mask"]], "notes")
  }
  data$sex_bad_calls <- sex_bad_calls
  if (length(read_notes)) data$mask_notes <- c(data$mask_notes, read_notes)
  data
}

scr_check_data <- function(data) {
  pr <- character()
  tr <- data$traps
  badxy <- is.na(tr$x) | is.na(tr$y)
  if (any(badxy)) pr <- c(pr, sprintf("trap file: %d trap(s) have missing or non-numeric coordinates (%s); coordinates must be plain numbers without commas, units or letters",
                                      sum(badxy), paste(utils::head(tr$id[badxy], 5), collapse = ", ")))
  if (any(duplicated(tr$id))) pr <- c(pr, paste("trap file: duplicated trap IDs:", paste(unique(tr$id[duplicated(tr$id)]), collapse = ", ")))
  if (any(duplicated(tr[, c("x", "y")]))) pr <- c(pr, "trap file: two traps share the same coordinates; merge them or check the coordinates")
  cap <- data$captures
  crow <- if (!is.null(cap$row)) cap$row else seq_len(nrow(cap))
  if (!nrow(cap)) pr <- c(pr, "capture file has no detections")
  unk <- setdiff(unique(cap$trap), tr$id)
  if (length(unk)) pr <- c(pr, sprintf("capture file: %d trap IDs are not in the trap file (e.g. %s); check spelling/case",
                                        length(unk), paste(utils::head(unk, 5), collapse = ", ")))
  if (anyNA(cap$occasion)) pr <- c(pr, sprintf("capture file: %d row(s) have a missing or non-whole-number occasion (data row(s) %s)",
                                                sum(is.na(cap$occasion)), paste(utils::head(crow[is.na(cap$occasion)], 5), collapse = ", ")))
  K <- data$n_occasions
  badk <- cap$occasion[!is.na(cap$occasion) & (cap$occasion < 1 | cap$occasion > K)]
  if (length(badk)) pr <- c(pr, sprintf("capture file: occasions must be 1..%d (n_occasions) but found %s", K, paste(utils::head(unique(badk), 5), collapse = ", ")))
  noid <- is.na(cap$individual) | cap$individual %in% c("", "NA")
  if (any(noid)) pr <- c(pr, sprintf("capture file: %d row(s) have no animal ID (data row(s) %s); remove them or fill in the ID%s",
                                     sum(noid), paste(utils::head(crow[noid], 5), collapse = ", "),
                                     if (is.null(data$qc)) " (for genetic samples, set qc_col: to drop samples that failed genotyping)" else ""))
  if (!is.null(cap$count) && (anyNA(cap$count) || any(cap$count < 0) || any(cap$count != round(cap$count))))
    pr <- c(pr, "capture file: the sample-count column must hold whole numbers >= 0")
  u <- data$usage
  if (anyNA(u) || any(u < 0, na.rm = TRUE)) {
    bad <- which(is.na(u) | u < 0, arr.ind = TRUE)
    pr <- c(pr, sprintf("trap operation / effort: %d value(s) are missing, not numbers or negative (e.g. %s); every detector x occasion needs a number >= 0 (0 = detector not operating)",
                        nrow(bad), paste(utils::head(sprintf("detector %s occasion %d", tr$id[bad[, 1]], bad[, 2]), 3), collapse = "; ")))
  }
  if (!length(unk) && !anyNA(cap$occasion) && !length(badk) && !anyNA(u)) {
    off <- u[cbind(match(cap$trap, tr$id), cap$occasion)] == 0
    if (any(off)) pr <- c(pr, sprintf("%d detections are at traps marked as not operating on that occasion (e.g. %s occasion %d); fix the usage/effort table",
                                      sum(off), cap$trap[off][1], cap$occasion[off][1]))
  }
  obs <- if (is.null(data$obs_model)) "binomial" else data$obs_model
  if (!obs %in% c("binomial", "bernoulli", "poisson")) pr <- c(pr, sprintf("obs_model must be binomial, bernoulli or poisson (lower case; got '%s')", obs))
  if (!is.null(data$behaviour) && !data$behaviour %in% c("none", "trap"))
    pr <- c(pr, sprintf("behaviour must be none or trap (trap-specific learned response at the detector where the animal was first detected; got '%s')", data$behaviour))
  if (length(data$sex_bad_calls)) pr <- c(pr, sprintf("sex column: value(s) not recognised (use F/M, female/male, or leave empty/U for unknown): %s",
                                                      paste(utils::head(data$sex_bad_calls, 5), collapse = ", ")))
  if (!is.na(data$units) && data$units %in% c("m", "km") && !any(badxy) && all(abs(tr$x) <= 180) && all(abs(tr$y) <= 90))
    pr <- c(pr, sprintf("trap file: the coordinates look like longitude/latitude in degrees, but units is '%s'; project them (and any mask, polygon or raster) to a metric CRS such as UTM", data$units))
  for (v in unique(c(data$detection_covariates, data$sigma_covariates))) if (anyNA(tr[[v]])) pr <- c(pr, sprintf("trap file: detector covariate '%s' has missing values", v))
  if (length(data$occasion_detection_covariates)) {
    oc <- data$occasion_covariates
    if (is.null(oc)) pr <- c(pr, "occasion covariates listed but no occasions: file given")
    else {
      if (nrow(oc) != K || !all(sort(oc$occasion) == seq_len(K))) pr <- c(pr, sprintf("occasions file must have exactly one row per occasion 1..%d", K))
      for (v in data$occasion_detection_covariates) {
        if (!v %in% names(oc)) pr <- c(pr, sprintf("occasions file has no column '%s'", v))
        else if (anyNA(suppressWarnings(as.numeric(oc[[v]])))) pr <- c(pr, sprintf("occasions file: covariate '%s' has missing or non-numeric values", v))
      }
    }
  }
  if (!is.null(data$individuals) && !is.null(data$individuals$sex_raw)) {
    ind <- data$individuals   # animals missing from this table, or with empty/U sex, are treated as unknown sex
    bad <- !is.na(ind$sex_raw) & nzchar(ind$sex_raw) & is.na(ind$sex) & !toupper(ind$sex_raw) %in% c("U", "UNK", "UNKNOWN", "?")
    if (any(bad)) pr <- c(pr, sprintf("sex values not recognised (use F/M, female/male, or leave empty/U for unknown): %s", paste(utils::head(unique(ind$sex_raw[bad]), 5), collapse = ", ")))
  }
  if (!is.null(data[["mask"]])) {
    m <- data[["mask"]]
    if (nrow(m) < 50) pr <- c(pr, sprintf("habitat mask has only %d pixels; use a finer grid (aim for 1000-4000)", nrow(m)))
    for (v in data$density_covariates) if (!v %in% names(m)) pr <- c(pr, sprintf("habitat mask has no column '%s'", v))
    sp <- tryCatch(scr_infer_spacing(m), error = function(e) NA)
    if (is.finite(sp) && !any(badxy)) {
      off <- abs((m$x - min(m$x)) / sp - round((m$x - min(m$x)) / sp)) > 0.01 | abs((m$y - min(m$y)) / sp - round((m$y - min(m$y)) / sp)) > 0.01
      if (any(off)) pr <- c(pr, sprintf("habitat mask: %d points are not on a regular grid of spacing %g", sum(off), sp))
      dt <- scr_dist_to_traps(tr, data.frame(id = seq_len(nrow(m)), x = m$x, y = m$y))
      if (any(dt > sp)) pr <- c(pr, sprintf("%d traps lie outside the habitat mask (more than one pixel from any habitat pixel: %s); the mask must cover the traps plus a buffer, and traps should be in habitat",
                                           sum(dt > sp), paste(utils::head(tr$id[dt > sp], 5), collapse = ", ")))
    }
  } else if (length(data$density_covariates)) pr <- c(pr, "density covariates need a habitat mask CSV (data: mask) or a raster (data: raster)")
  if (!is.na(data$units) && !data$units %in% c("m", "km")) pr <- c(pr, "data: units must be m or km (or omitted if unknown)")
  pr
}

# Write simulated data as the CSVs a researcher would have (used by `core/run.R simulate`)
scr_write_data <- function(data, out) {
  tr <- data$traps
  if (all(data$usage %in% c(0, 1))) tr$usage <- apply(data$usage, 1, paste, collapse = "")
  else utils::write.csv(data.frame(id = tr$id, data$usage, check.names = FALSE), file.path(out, "scr_effort.csv"), row.names = FALSE)
  utils::write.csv(tr, file.path(out, "scr_traps.csv"), row.names = FALSE)
  utils::write.csv(data$captures, file.path(out, "scr_captures.csv"), row.names = FALSE)
  if (!is.null(data$individuals)) utils::write.csv(data$individuals, file.path(out, "scr_individuals.csv"), row.names = FALSE)
  if (!is.null(data$occasion_covariates)) utils::write.csv(data$occasion_covariates, file.path(out, "scr_occasions.csv"), row.names = FALSE)
  if (!is.null(data[["mask"]])) utils::write.csv(data[["mask"]], file.path(out, "scr_habitat.csv"), row.names = FALSE)
  invisible(out)
}
