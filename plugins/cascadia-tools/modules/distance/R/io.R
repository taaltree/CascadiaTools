# Reading researchers' files.
#
# config$data:
#   survey:          line | point
#   observations:    CSV, one row per detected animal or group (rows with an empty
#                    distance are treated as effort-only rows and dropped)
#   transects:       CSV, one row per transect line / point, INCLUDING those with no detections
#   transect_col:    id column present in both files (default "transect")
#   distance_col:    perpendicular (line) or radial (point) distance (default "distance_m")
#   distance_units:  m (default) | km | ft | yd
#   distance_type:   exact (default) | binned
#   cutpoints_m:     bin edges for binned data (e.g. [0, 25, 50, 100, 200]); with
#                    distbegin_col / distend_col, or distances are assigned to bins
#   size_col:        group (cluster) size, optional
#   visit_col:       visit number for repeated surveys (optional; enables gdistsamp)
#   length_col:      line length (lines; default "length_km"); length_units km | m
#   visits_col:      number of times each transect/point was surveyed (default 1)
#   stratum_col:     stratum of each transect (optional)
#   study_area_km2:  number, or a map {stratum: area} (optional; gives abundance)
#   truncation_m:    number, or "10%" to drop the farthest 10% of distances
#   prediction_grid: CSV of the study area: x, y, area_km2 and the density covariates (optional)
#   species, group_word ("animal"/"group"/"covey"/"pod" ...)
# config$formulas: density ("~ shrub_cover_pct"), detection ("~ 1" or "~ observer")

dist_unit_m <- c(m = 1, km = 1000, ft = 0.3048, yd = 0.9144)

# Numeric column with a clear error naming the rows and values that are not numbers
# (e.g. "85m", "1,20"), instead of silently turning them into NA.
dist_num <- function(x, col, file) {
  if (is.numeric(x)) return(x)
  xs <- trimws(as.character(x))
  v <- suppressWarnings(as.numeric(xs))
  bad <- is.na(v) & !is.na(xs) & xs != ""
  if (any(bad)) stop(sprintf("column '%s' of the %s file must contain plain numbers, but %d value(s) are not numbers: %s. Remove units or text (write 85, not 85m) and use a dot as the decimal mark.",
                             col, file, sum(bad), paste(sprintf("row %d = '%s'", utils::head(which(bad), 5) + 1, utils::head(xs[bad], 5)), collapse = ", ")), call. = FALSE)
  v
}

dist_read_data <- function(config) {
  dc <- config$data
  rd <- function(f) utils::read.csv(file.path(config$.dir, f), stringsAsFactors = FALSE, check.names = FALSE)
  ob <- rd(dc$observations)
  tr <- rd(dc$transects)
  tcol <- dc$transect_col %||% "transect"
  dcol <- dc$distance_col %||% "distance_m"
  cols <- function(d) paste0(" Columns found: ", paste(names(d), collapse = ", "), ".")
  if (!tcol %in% names(ob)) stop("observations file has no column '", tcol, "' (set data$transect_col).", cols(ob))
  if (!tcol %in% names(tr)) stop("transects file has no column '", tcol, "' (set data$transect_col).", cols(tr))
  survey <- dc$survey %||% "line"
  dtype <- dc$distance_type %||% "exact"
  um <- dist_unit_m[[dc$distance_units %||% "m"]]
  keycol <- if (dtype == "binned" && !is.null(dc$distbegin_col)) dc$distbegin_col else dcol
  if (!keycol %in% names(ob)) stop("observations file has no column '", keycol, "' (set data$distance_col).", cols(ob))
  blank <- is.na(ob[[keycol]]) | ob[[keycol]] == ""
  if (any(blank)) message(sprintf("Note: %d observation row(s) with an empty distance (file rows %s) were treated as effort-only rows and dropped; if these were real detections, fill in the distance.",
                                  sum(blank), paste(utils::head(which(blank) + 1, 10), collapse = ", ")))
  ob <- ob[!blank, , drop = FALSE]
  o <- data.frame(transect = as.character(ob[[tcol]]),
                  distance = if (dcol %in% names(ob)) dist_num(ob[[dcol]], dcol, "observations") * um else NA_real_, stringsAsFactors = FALSE)
  cut <- NULL
  if (dtype == "binned") {
    cut <- as.numeric(unlist(dc$cutpoints_m))
    if (!length(cut)) stop("binned distances need data$cutpoints_m")
    if (!is.null(dc$distbegin_col)) {
      o$distbegin <- dist_num(ob[[dc$distbegin_col]], dc$distbegin_col, "observations") * um; o$distend <- dist_num(ob[[dc$distend_col]], dc$distend_col, "observations") * um
    } else {
      j <- findInterval(o$distance, cut, rightmost.closed = TRUE, all.inside = TRUE)
      o$distbegin <- cut[j]; o$distend <- cut[j + 1]
    }
    o$distance <- (o$distbegin + o$distend) / 2
  }
  clusters <- !is.null(dc$size_col)
  if (clusters) {
    if (!dc$size_col %in% names(ob)) stop("observations file has no group-size column '", dc$size_col, "' (set data$size_col, or remove it if animals were counted singly).", cols(ob))
    o$size <- dist_num(ob[[dc$size_col]], dc$size_col, "observations")
  }
  if (!is.null(dc$visit_col)) o$visit <- as.integer(ob[[dc$visit_col]])
  # observation-level detection covariates
  fdt <- ds_terms(config$formulas$detection)
  for (v in fdt) if (v %in% names(ob) && !v %in% names(tr)) o[[v]] <- ob[[v]]

  t2 <- data.frame(transect = as.character(tr[[tcol]]), stringsAsFactors = FALSE)
  vis <- if (!is.null(dc$visits_col)) as.numeric(tr[[dc$visits_col]]) else rep(1, nrow(tr))
  if (survey == "line") {
    lc <- dc$length_col %||% "length_km"
    if (!lc %in% names(tr)) stop("transects file has no line-length column '", lc, "' (set data$length_col).", cols(tr))
    len <- dist_num(tr[[lc]], lc, "transects") * (if ((dc$length_units %||% "km") == "m") 0.001 else 1)
  } else len <- rep(1, nrow(tr))
  nV <- 1
  if (!is.null(dc$visit_col) && length(unique(vis)) == 1 && vis[1] > 1) { nV <- vis[1]; t2$effort <- len }
  else { t2$effort <- len * vis; o$visit <- 1 }
  t2$stratum <- if (!is.null(dc$stratum_col)) as.character(tr[[dc$stratum_col]]) else "all"
  for (v in setdiff(names(tr), c(tcol, dc$stratum_col))) t2[[v]] <- tr[[v]]
  strata <- NULL
  if (!is.null(dc$study_area_km2)) {
    a <- dc$study_area_km2
    strata <- if (is.list(a)) data.frame(stratum = names(a), area_km2 = as.numeric(unlist(a)))
              else data.frame(stratum = unique(t2$stratum)[1], area_km2 = as.numeric(a))
    if (!is.list(a) && length(unique(t2$stratum)) > 1) stop("study_area_km2 must give an area for every stratum")
  }
  w <- dc$truncation_m
  if (is.null(w)) w <- if (dtype == "binned") max(cut) else max(o$distance)
  if (is.character(w) && grepl("%$", w)) w <- stats::quantile(o$distance, 1 - as.numeric(sub("%", "", w)) / 100, names = FALSE)
  w <- as.numeric(w)
  grid <- if (!is.null(dc$prediction_grid)) rd(dc$prediction_grid) else NULL
  list(survey = survey, distance_type = dtype, truncation_m = w, cutpoints = cut, obs = o, transects = t2,
       strata = strata, n_visits = nV, clusters = clusters,
       formulas = list(density = config$formulas$density %||% "~ 1", detection = config$formulas$detection %||% "~ 1"),
       state_covs = ds_terms(config$formulas$density), det_covs = fdt, scaling = list(),
       species = dc$species %||% "animals", group_word = dc$group_word %||% (if (clusters) "group" else "animal"),
       grid = grid, labels = config$figures$labels, raw_distance_units = dc$distance_units %||% "m")
}

# Fatal problems only (fit stops); softer warnings are added to the fit notes.
dist_check_data <- function(data) {
  pr <- character()
  o <- data$obs; tr <- data$transects
  if (!data$survey %in% c("line", "point")) pr <- c(pr, "data$survey must be 'line' or 'point'")
  miss <- setdiff(unique(o$transect), tr$transect)
  if (length(miss)) {
    hint <- vapply(utils::head(miss, 3), function(x) { a <- agrep(x, tr$transect, max.distance = 1, value = TRUE)
      if (length(a)) sprintf("%s (did you mean %s?)", x, a[1]) else x }, "")
    pr <- c(pr, sprintf("%d transect id(s) in the observations are missing from the transects file (e.g. %s). Every surveyed transect, including those with no detections, must be listed there, with the same spelling in both files.",
                        length(miss), paste(hint, collapse = ", ")))
  }
  if (any(is.na(o$distance) & !(data$distance_type == "binned" & !is.null(o$distbegin)))) pr <- c(pr, sprintf("%d observation(s) have no usable distance", sum(is.na(o$distance))))
  if (anyDuplicated(tr$transect)) pr <- c(pr, "transect ids are duplicated in the transects file")
  if (any(o$distance < 0, na.rm = TRUE)) pr <- c(pr, sprintf("%d negative distance(s)", sum(o$distance < 0, na.rm = TRUE)))
  if (any(!is.finite(tr$effort)) || any(tr$effort <= 0)) pr <- c(pr, "every transect needs a positive length (lines) / number of visits (points)")
  if (isTRUE(data$clusters) && any(bad <- (!is.finite(o$size) | o$size < 1))) pr <- c(pr, sprintf("group sizes must be whole numbers >= 1: %d missing or invalid (transect(s) %s; fill in missing sizes or drop size_col)",
                                                                                       sum(bad), paste(utils::head(unique(o$transect[bad]), 5), collapse = ", ")))
  if (!is.null(data$strata)) {
    ms <- setdiff(unique(tr$stratum), data$strata$stratum)
    if (length(ms)) pr <- c(pr, paste("study_area_km2 has no area for stratum:", paste(ms, collapse = ", ")))
  }
  for (v in c(ds_terms(data$formulas$density))) if (!v %in% names(tr)) pr <- c(pr, sprintf("density covariate '%s' is not a column of the transects file", v))
  for (v in ds_terms(data$formulas$detection)) if (!v %in% c(names(tr), names(o))) pr <- c(pr, sprintf("detection covariate '%s' not found in either file", v))
  for (v in intersect(ds_terms(data$formulas$density), names(tr))) if (anyNA(tr[[v]])) pr <- c(pr, sprintf("covariate '%s' has missing values for transect(s) %s; fill them in or remove the covariate from the density formula", v, paste(utils::head(tr$transect[is.na(tr[[v]])], 5), collapse = ", ")))
  nin <- if (data$distance_type == "binned") sum(o$distend <= data$truncation_m + 1e-9, na.rm = TRUE) else sum(o$distance <= data$truncation_m, na.rm = TRUE)
  if (nin < 20) pr <- c(pr, sprintf("only %d detections within the truncation distance: too few to fit a detection function (aim for 60-80+)", nin))
  pr
}

# Non-fatal warnings about data quality (heaping, few detections, spikes).
dist_data_warnings <- function(data) {
  o <- data$obs; w <- data$truncation_m
  out <- character()
  if (data$distance_type == "exact") {
    x <- o$distance[o$distance <= w]
    z <- mean(x == 0)
    if (z > 0.05) out <- c(out, sprintf("%.0f%% of distances are exactly 0: possible rounding/heaping at zero, which biases density upward. Consider binning the data (distance_type: binned) with a first bin that absorbs the heaping.", 100 * z))
    r <- x[x > 0]
    if (length(r) > 20 && mean(abs(r / 10 - round(r / 10)) < 1e-8) > 0.6) out <- c(out, "Most distances are multiples of 10: heaping to round numbers. If the Q-Q plot shows steps, analyse the data in bins.")
  }
  out
}
