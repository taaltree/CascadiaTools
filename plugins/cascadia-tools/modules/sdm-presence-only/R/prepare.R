# Data input, checks and model matrices for sdm-presence-only.
#
# data (from read_data() or simulate()):
#   grid      data.frame(x, y, area_km2, <habitat covariates>, <bias covariates>)
#             one row per cell of a REGULAR grid covering the study region
#             (the region the presence records could have come from)
#   presence  data.frame(x, y) presence-only records (same coordinate system)
#   survey    NULL or list(sites = data.frame(site, x, y, area_km2),
#                          y = sites x visits matrix of 0/1/NA)
#   habitat   names of habitat covariates (drive species intensity)
#   bias      names of sampling-bias covariates (drive where records are made)
#   quadratic habitat covariates that also get a squared term
#   coord_unit "m" or "km"; labels: named list of axis labels

sdm_cell_index <- function(grid, x, y) {
  ux <- sort(unique(grid$x)); uy <- sort(unique(grid$y))
  dx <- if (length(ux) > 1) min(diff(ux)) else 1; dy <- if (length(uy) > 1) min(diff(uy)) else 1
  gx <- round((grid$x - ux[1]) / dx); gy <- round((grid$y - uy[1]) / dy)
  nxx <- max(gx) + 1
  key <- gx + gy * nxx
  lut <- rep(NA_integer_, (max(gy) + 1) * nxx); lut[key + 1] <- seq_len(nrow(grid))
  ix <- round((x - ux[1]) / dx); iy <- round((y - uy[1]) / dy)
  ok <- ix >= 0 & ix < nxx & iy >= 0 & iy <= max(gy)
  out <- rep(NA_integer_, length(x))
  out[ok] <- lut[ix[ok] + iy[ok] * nxx + 1]
  out
}

sdm_grid_spacing <- function(grid) {
  ux <- sort(unique(grid$x)); uy <- sort(unique(grid$y))
  c(dx = if (length(ux) > 1) min(diff(ux)) else NA, dy = if (length(uy) > 1) min(diff(uy)) else NA)
}

# Standardise covariates over the grid, build habitat (X) and bias (Z)
# matrices for cells, presences and survey sites.
sdm_prepare <- function(data, settings = list()) {
  g <- data$grid
  hab <- data$habitat; bias <- if (isFALSE(settings$bias_correction)) character() else data$bias
  quad <- intersect(data$quadratic, hab)
  vars <- unique(c(hab, data$bias))
  keep <- stats::complete.cases(g[, c("x", "y", "area_km2", vars), drop = FALSE])
  notes <- if (is.null(data$read_notes)) character() else data$read_notes
  for (b in data$bias) for (v in hab) {
    r <- suppressWarnings(stats::cor(g[[b]], g[[v]], use = "complete.obs"))
    if (is.finite(r) && abs(r) > 0.85) notes <- c(notes, sprintf(
      "Habitat covariate '%s' is strongly correlated with the sampling-bias covariate '%s' (r = %.2f): where people looked and this habitat can barely be told apart, so its effect has a wide interval. Do not drop the bias covariate to 'fix' this: that would credit sampling effort to habitat. Survey data with non-detections (engine: integrated) help most.", v, b, r))
  }
  if (any(!keep)) notes <- c(notes, sprintf("%d grid cells with missing covariate values were dropped from the study region.", sum(!keep)))
  g <- g[keep, , drop = FALSE]
  scaling <- lapply(stats::setNames(vars, vars), function(v) list(center = mean(g[[v]]), scale = stats::sd(g[[v]])))
  zf <- function(df, v) (df[[v]] - scaling[[v]]$center) / scaling[[v]]$scale
  mkX <- function(df) {
    X <- matrix(sapply(hab, function(v) zf(df, v)), nrow = nrow(df))
    for (v in quad) X <- cbind(X, zf(df, v)^2)
    colnames(X) <- c(hab, if (length(quad)) paste0(quad, "_sq"))
    X
  }
  mkZ <- function(df) {
    if (!length(bias)) return(matrix(0, nrow(df), 0))
    Z <- sapply(bias, function(v) zf(df, v)); matrix(Z, nrow = nrow(df), dimnames = list(NULL, bias))
  }
  pc <- sdm_cell_index(g, data$presence$x, data$presence$y)
  if (any(is.na(pc))) notes <- c(notes, sprintf("%d presence records fell outside the covariate grid (or in cells with missing covariates) and were dropped.", sum(is.na(pc))))
  pc <- pc[!is.na(pc)]
  prep <- list(grid = g, X = mkX(g), Z = mkZ(g), A = g$area_km2, pres_cell = pc,
               habitat = hab, bias = bias, bias_all = data$bias, quadratic = quad, scaling = scaling,
               notes = notes, n_pres = length(pc), mkX = mkX, mkZ = mkZ, zf = zf)
  if (!is.null(data$survey)) {
    s <- data$survey
    sc <- sdm_cell_index(g, s$sites$x, s$sites$y)
    ok <- !is.na(sc) & rowSums(!is.na(s$y)) > 0
    if (any(!ok)) prep$notes <- c(prep$notes, sprintf("%d survey sites outside the grid or never visited were dropped.", sum(!ok)))
    Y <- s$y[ok, , drop = FALSE]
    prep$survey <- list(cell = sc[ok], A = s$sites$area_km2[ok], ydet = rowSums(Y, na.rm = TRUE),
                        J = rowSums(!is.na(Y)), site = s$sites$site[ok], y = Y,
                        X = mkX(g[sc[ok], , drop = FALSE]))
    prep$survey$estimate_p <- max(prep$survey$J) > 1 && !isTRUE(settings$perfect_detection)
  }
  prep
}

# ---- reading researchers' files --------------------------------------------

sdm_read_data <- function(config) {
  d <- config$data; dir <- config$.dir
  path <- function(f) if (grepl("^/", f)) f else file.path(dir, f)
  hab <- unlist(d$habitat_covariates); bias <- unlist(d$bias_covariates)
  vars <- c(hab, bias)
  # covariate grid: CSV (x, y, covariates) or rasters
  if (!is.null(d$covariate_grid)) {
    g <- utils::read.csv(path(d$covariate_grid), check.names = FALSE)
    xc <- if (is.null(d$grid_x)) "x" else d$grid_x; yc <- if (is.null(d$grid_y)) "y" else d$grid_y
    miss <- setdiff(c(xc, yc, vars), names(g))
    if (length(miss)) stop("covariate grid ", d$covariate_grid, " lacks column(s): ", paste(miss, collapse = ", "),
                           ". Columns present: ", paste(names(g), collapse = ", "))
    g <- data.frame(x = g[[xc]], y = g[[yc]], g[, vars, drop = FALSE], check.names = FALSE)
  } else if (!is.null(d$rasters)) {
    eco_require("terra")
    r <- terra::rast(lapply(unlist(d$rasters), function(f) terra::rast(path(f))))
    names(r) <- names(unlist(d$rasters))
    miss <- setdiff(vars, names(r))
    if (length(miss)) stop("rasters: no layer named ", paste(miss, collapse = ", "), " (name each raster in config data$rasters)")
    g <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
    g <- g[, c("x", "y", vars)]
    if (is.null(d$coord_unit) && isTRUE(terra::is.lonlat(r))) stop("rasters are in longitude/latitude: project them to a metric CRS (e.g. UTM) first")
  } else stop("config data needs covariate_grid (CSV) or rasters")
  unit <- if (is.null(d$coord_unit)) "m" else d$coord_unit
  sp <- sdm_grid_spacing(g)
  f <- switch(unit, m = 1e-6, km = 1, stop("coord_unit must be m or km"))
  g$area_km2 <- unname(sp["dx"] * sp["dy"] * f)
  # presences
  p <- utils::read.csv(path(d$presence), check.names = FALSE)
  px <- if (is.null(d$presence_x)) "x" else d$presence_x; py <- if (is.null(d$presence_y)) "y" else d$presence_y
  if (!all(c(px, py) %in% names(p))) stop("presence file ", d$presence, " needs columns ", px, " and ", py,
                                          " (set presence_x / presence_y in config). Columns present: ", paste(names(p), collapse = ", "))
  if (!is.null(d$species_col) && !is.null(d$species)) p <- p[p[[d$species_col]] %in% d$species, , drop = FALSE]
  pres <- data.frame(x = p[[px]], y = p[[py]])
  read_notes <- character()
  if (any(!stats::complete.cases(pres))) {
    bad <- which(!stats::complete.cases(pres))
    read_notes <- c(read_notes, sprintf("%d presence record(s) with a missing coordinate were dropped (file rows %s).",
                                        length(bad), paste(utils::head(bad + 1, 10), collapse = ", ")))
  }
  pres <- pres[stats::complete.cases(pres), , drop = FALSE]
  # survey (optional): long (site, x, y, visit, detected) or wide
  survey <- NULL; survey_issues <- character()
  if (!is.null(d$survey)) {
    s <- utils::read.csv(path(d$survey), check.names = FALSE)
    sc <- if (is.null(d$survey_site_col)) "site" else d$survey_site_col
    sx <- if (is.null(d$survey_x)) px else d$survey_x; sy <- if (is.null(d$survey_y)) py else d$survey_y
    if (identical(d$survey_format, "wide")) {
      vc <- unlist(d$survey_visit_cols)
      Y <- as.matrix(s[, vc, drop = FALSE]); sites <- s[, c(sc, sx, sy)]
    } else {
      vcol <- if (is.null(d$survey_visit_col)) "visit" else d$survey_visit_col
      dcol <- if (is.null(d$survey_detection_col)) "detected" else d$survey_detection_col
      miss <- setdiff(c(sc, sx, sy, vcol, dcol), names(s))
      if (length(miss)) stop("survey file lacks column(s): ", paste(miss, collapse = ", "))
      ids <- unique(s[[sc]])
      vis <- sort(unique(s[[vcol]]))
      # catch typos in site ids: one id at several locations, or one location under several ids
      xy <- paste(s[[sx]], s[[sy]])
      nxy <- tapply(xy, s[[sc]], function(v) length(unique(v)))
      for (k in names(nxy)[nxy > 1]) survey_issues <- c(survey_issues, sprintf(
        "survey site '%s' appears with %d different coordinates; each site needs one location (check for a typo in the site id or coordinates)", k, nxy[[k]]))
      nid <- tapply(as.character(s[[sc]]), xy, function(v) unique(v))
      for (k in nid[lengths(nid) > 1]) survey_issues <- c(survey_issues, sprintf(
        "survey sites %s have exactly the same coordinates: is one of them a misspelled site id? Fix the id (or, if they are distinct sites, give each its own coordinates)",
        paste0("'", k, "'", collapse = " and ")))
      dup <- duplicated(paste(s[[sc]], s[[vcol]]))
      if (any(dup)) survey_issues <- c(survey_issues, sprintf("survey file has %d duplicated site x visit rows (e.g. site '%s' visit %s)",
                                                              sum(dup), s[[sc]][dup][1], s[[vcol]][dup][1]))
      Y <- matrix(NA_real_, length(ids), length(vis))
      Y[cbind(match(s[[sc]], ids), match(s[[vcol]], vis))] <- s[[dcol]]
      sites <- s[match(ids, s[[sc]]), c(sc, sx, sy)]
    }
    a <- if (!is.null(d$survey_area_col) && d$survey_area_col %in% names(s)) s[match(sites[[sc]], s[[sc]]), d$survey_area_col]
         else rep(if (is.null(d$survey_area_km2)) stop("give survey_area_km2 (area each survey covers, km2) or survey_area_col") else d$survey_area_km2, nrow(sites))
    rownames(Y) <- as.character(sites[[sc]])
    survey <- list(sites = data.frame(site = sites[[sc]], x = sites[[sx]], y = sites[[sy]], area_km2 = a), y = Y,
                   issues = survey_issues)
  }
  labs <- d$labels; if (is.null(labs)) labs <- list()
  list(grid = g, presence = pres, survey = survey, habitat = hab, bias = bias,
       quadratic = unlist(d$quadratic), coord_unit = unit, labels = labs, read_notes = read_notes)
}

sdm_check_data <- function(data) {
  pr <- character()
  g <- data$grid
  if (!length(data$habitat)) pr <- c(pr, "no habitat covariates given (data: habitat_covariates)")
  vars <- c(data$habitat, data$bias)
  for (v in vars) {
    if (!is.numeric(g[[v]])) pr <- c(pr, sprintf("covariate '%s' is not numeric", v))
    else if (stats::sd(g[[v]], na.rm = TRUE) == 0) pr <- c(pr, sprintf("covariate '%s' does not vary across the grid", v))
  }
  if (all(abs(range(g$x, na.rm = TRUE)) <= 180) && all(abs(range(g$y, na.rm = TRUE)) <= 90) && identical(data$coord_unit, "m"))
    pr <- c(pr, "grid coordinates look like longitude/latitude (degrees). Project grid, presences and survey sites to a metric CRS (e.g. UTM, metres) first; cell areas are computed from the coordinates.")
  sp <- sdm_grid_spacing(g)
  if (any(!is.finite(sp))) pr <- c(pr, "the covariate grid needs at least 2 distinct x and y values")
  else {
    offx <- abs(((g$x - min(g$x)) / sp["dx"]) - round((g$x - min(g$x)) / sp["dx"]))
    offy <- abs(((g$y - min(g$y)) / sp["dy"]) - round((g$y - min(g$y)) / sp["dy"]))
    if (max(offx, offy) > 0.01) pr <- c(pr, "covariate grid is not a regular grid (cell centres are not evenly spaced); resample it to a regular grid")
    if (anyDuplicated(paste(g$x, g$y))) pr <- c(pr, "covariate grid has duplicated x,y cell centres")
  }
  if (nrow(data$presence) < 10) pr <- c(pr, sprintf("only %d presence records: at least ~20-30 are needed for even one covariate effect", nrow(data$presence)))
  if (length(vars) > 1 && nrow(g) > 2) {
    cm <- stats::cor(g[, vars], use = "complete.obs"); cm[upper.tri(cm, TRUE)] <- 0
    hi <- which(abs(cm) > 0.85, arr.ind = TRUE)
    for (k in seq_len(nrow(hi))) {
      a <- rownames(cm)[hi[k, 1]]; b <- colnames(cm)[hi[k, 2]]
      # a habitat covariate collinear with a BIAS covariate is not a data error: the model stays
      # unbiased with wide intervals (verified by simulation), and dropping the bias covariate
      # would wrongly attribute effort to habitat. fit() reports it as a note instead.
      if (xor(a %in% data$bias, b %in% data$bias)) next
      pr <- c(pr, sprintf("covariates '%s' and '%s' are almost collinear (r = %.2f): their effects cannot be separated; drop one",
                          a, b, cm[hi[k, 1], hi[k, 2]]))
    }
  }
  if (!is.null(data$survey)) {
    y <- data$survey$y
    if (!all(y %in% c(0, 1, NA))) {
      bad <- which(!(y %in% c(0, 1, NA)) , arr.ind = FALSE)
      bi <- arrayInd(bad, dim(y))
      sid <- if (!is.null(rownames(y))) rownames(y)[bi[, 1]] else as.character(data$survey$sites$site[bi[, 1]])
      pr <- c(pr, sprintf("survey detections must be 0/1 (NA or absent row for visits not done); found %s",
                          paste(utils::head(sprintf("'%s' at site %s visit %d", y[bad], sid, bi[, 2]), 5), collapse = ", ")))
    }
    pr <- c(pr, data$survey$issues)
    if (any(rowSums(!is.na(y)) == 0)) pr <- c(pr, sprintf("%d survey sites have no visits recorded", sum(rowSums(!is.na(y)) == 0)))
    if (any(!is.finite(data$survey$sites$area_km2) | data$survey$sites$area_km2 <= 0)) pr <- c(pr, "survey area must be a positive number of km2")
  }
  pr
}

# write simulated data in the researcher-facing CSV format (core/run.R simulate)
sdm_write_data <- function(data, out) {
  utils::write.csv(data$presence, file.path(out, "presence.csv"), row.names = FALSE)
  utils::write.csv(data$grid, file.path(out, "covariate_grid.csv"), row.names = FALSE)
  if (!is.null(data$survey)) {
    s <- data$survey
    long <- do.call(rbind, lapply(seq_len(ncol(s$y)), function(j)
      data.frame(s$sites[, c("site", "x", "y")], visit = j, detected = s$y[, j])))
    utils::write.csv(long[!is.na(long$detected), ], file.path(out, "survey.csv"), row.names = FALSE)
  }
}
