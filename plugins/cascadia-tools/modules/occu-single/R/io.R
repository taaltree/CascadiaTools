# occu-single: reading researchers' files and checking them.
#
# Long format (default; one row per survey actually done):
#   site, [season], survey, detected, [survey covariates...]
# Wide format (one row per site; detection columns in survey order, NA = not surveyed):
#   site, <detection columns...>; survey covariates as extra wide files.
# Site covariates: one row per site. Season covariates (dynamic): site, season, covariates.

occ_read_csv <- function(dir, f) {
  if (is.null(f)) return(NULL)
  p <- if (file.exists(f)) f else file.path(dir, f)
  if (!file.exists(p)) stop("file not found: ", f, " (looked in ", dir, ")")
  # Blank cells are missing values in every column. Without na.strings = "" a blank
  # cell in a TEXT column (observer, habitat) is read as "", which silently becomes
  # an extra factor level (and, sorting first, the reference level), so the model
  # gets a fake category instead of a "missing value" message from check_data().
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE,
                  na.strings = c("NA", ""))
}

# Detection values that are text ("yes", "Y", "present") would otherwise become
# NA silently, i.e. be treated as surveys not done. Report them instead.
occ_text_detections <- function(raw, num, what, where) {
  raw <- as.vector(raw); num <- as.vector(num)
  bad <- !is.na(raw) & trimws(as.character(raw)) != "" & is.na(num)
  if (!any(bad)) return(character())
  sprintf("%s has %d value(s) that are not numbers, e.g. '%s' at %s (values found: %s). Use 1 = detected, 0 = not detected, and leave the cell blank (or drop the row) for a survey that was not done.",
          what, sum(bad), raw[bad][1], where[bad][1], paste(utils::head(unique(raw[bad]), 6), collapse = ", "))
}

# Time order of seasons. Numbers (2021, 2022) sort numerically. Text labels
# that share one pattern with a single number ("Y1".."Y12", "season 3") sort by
# that number. Anything else ("spring", "autumn", "2021-22") needs an explicit
# data$season_order in the config, because alphabetical order would silently
# scramble the time order of a dynamic model.
occ_season_order <- function(x, season_order, col) {
  u <- unique(x[!is.na(x)])
  if (!is.null(season_order)) {
    so <- as.character(unlist(season_order))
    miss <- setdiff(as.character(u), so)
    if (length(miss)) stop(sprintf("season(s) %s in column '%s' are not listed in data$season_order (%s)",
                                   paste(miss, collapse = ", "), col, paste(so, collapse = ", ")))
    return(so[so %in% as.character(u)])
  }
  if (is.numeric(u)) return(sort(u))
  u <- as.character(u)
  num <- suppressWarnings(as.numeric(u))
  if (!anyNA(num)) return(u[order(num)])
  digits <- regmatches(u, gregexpr("[0-9]+", u))
  stem <- gsub("[0-9]+", "#", u)
  if (all(lengths(digits) == 1) && length(unique(stem)) == 1) return(u[order(as.numeric(unlist(digits)))])
  stop(sprintf("cannot tell the time order of the seasons in column '%s' (%s). Add the order to the config, e.g.\n  data:\n    season_order: [%s]\n(earliest first), or use numbers such as years.",
               col, paste(utils::head(u, 8), collapse = ", "), paste(utils::head(u, 8), collapse = ", ")))
}

occ_read_data <- function(config) {
  dc <- config$data; dir <- if (is.null(config$.dir)) "." else config$.dir
  if (is.null(dc)) stop("config needs a 'data:' section")
  model <- if (!is.null(dc$model)) dc$model else if (!is.null(config$settings$model)) config$settings$model else "static"
  model <- match.arg(model, c("static", "dynamic"))
  fmt <- if (is.null(dc$format)) "long" else dc$format
  site_col <- if (is.null(dc$site_col)) "site" else dc$site_col
  issues <- character()
  det <- occ_read_csv(dir, dc$detections)
  if (is.null(det)) stop("config data$detections (the detection file) is required")
  if (!site_col %in% names(det)) stop(sprintf("detection file has no '%s' column (columns: %s); set data$site_col",
                                              site_col, paste(names(det), collapse = ", ")))
  det[[site_col]] <- as.character(det[[site_col]])
  sc_raw <- occ_read_csv(dir, dc$site_covs)
  sites <- unique(c(det[[site_col]], if (!is.null(sc_raw)) as.character(sc_raw[[site_col]])))
  sites <- sort(unique(det[[site_col]]))
  obs_names <- unlist(dc$obs_covs)
  obs <- list(); season_labels <- NULL

  if (fmt == "long") {
    sv <- if (is.null(dc$survey_col)) "survey" else dc$survey_col
    dv <- if (is.null(dc$detection_col)) "detected" else dc$detection_col
    se <- if (is.null(dc$season_col)) "season" else dc$season_col
    need <- c(sv, dv, if (model == "dynamic") se, obs_names)
    miss <- setdiff(need, names(det))
    if (length(miss)) stop("detection file is missing column(s): ", paste(miss, collapse = ", "),
                           " (columns present: ", paste(names(det), collapse = ", "), ")")
    sn <- suppressWarnings(as.numeric(det[[sv]]))
    if (anyNA(sn) || any(sn != round(sn)) || any(sn < 1))
      stop(sprintf("survey column '%s' must hold visit numbers 1, 2, 3, ... within each site%s (found e.g. '%s')",
                   sv, if (model == "dynamic") " and season" else "",
                   det[[sv]][which(is.na(sn) | sn != round(sn) | sn < 1)[1]]))
    J <- max(sn)
    if (model == "dynamic") {
      season_labels <- occ_season_order(det[[se]], dc$season_order, se)
      tt <- match(as.character(det[[se]]), as.character(season_labels))
      if (anyNA(tt)) stop(sprintf("season column '%s' is blank in %d row(s), e.g. row %d of the detection file; every survey needs a season",
                                  se, sum(is.na(tt)), which(is.na(tt))[1] + 1))
    } else tt <- rep(1L, nrow(det))
    T <- max(tt)
    ii <- match(det[[site_col]], sites)
    key <- paste(ii, sn, tt)
    if (any(duplicated(key))) {
      d <- which(duplicated(key))[1]
      issues <- c(issues, sprintf("%d duplicated survey rows (same site, %ssurvey), e.g. site %s survey %s. Each survey must appear once.%s",
                                  sum(duplicated(key)), if (model == "dynamic") "season, " else "",
                                  det[[site_col]][d], det[[sv]][d],
                                  if (model == "static") " If the repeated rows are different seasons or years at the same site, this is multi-season data: set data: model: dynamic and data: season_col: <column>." else ""))
    }
    yv <- suppressWarnings(as.numeric(det[[dv]]))
    issues <- c(issues, occ_text_detections(det[[dv]], yv, sprintf("detection column '%s'", dv),
                                            sprintf("site %s survey %s", det[[site_col]], det[[sv]])))
    y <- array(NA_real_, c(length(sites), J, T))
    y[cbind(ii, sn, tt)] <- yv
    for (k in obs_names) {
      v <- det[[k]]
      a <- array(if (is.numeric(v)) NA_real_ else NA_character_, c(length(sites), J, T))
      a[cbind(ii, sn, tt)] <- v
      obs[[k]] <- a
    }
  } else if (fmt == "wide") {
    dcols <- if (!is.null(dc$detection_cols)) unlist(dc$detection_cols) else setdiff(names(det), site_col)
    miss <- setdiff(dcols, names(det)); if (length(miss)) stop("detection columns not found: ", paste(miss, collapse = ", "))
    det <- det[match(sites, det[[site_col]]), , drop = FALSE]
    ym0 <- as.matrix(det[, dcols, drop = FALSE])
    ym <- suppressWarnings(matrix(as.numeric(ym0), nrow(ym0), dimnames = dimnames(ym0)))
    issues <- c(issues, occ_text_detections(ym0, ym, "detection columns",
                                            sprintf("site %s column %s", rep(sites, ncol(ym0)), rep(dcols, each = nrow(ym0)))))
    if (model == "dynamic") {
      J <- dc$surveys_per_season
      if (is.null(J)) stop("wide format for a dynamic model needs data$surveys_per_season")
      T <- length(dcols) / J
      if (T != round(T)) stop("number of detection columns is not a multiple of surveys_per_season")
      season_labels <- if (!is.null(dc$season_labels)) unlist(dc$season_labels) else seq_len(T)
    } else { J <- length(dcols); T <- 1 }
    y <- array(ym, c(length(sites), J, T))
    for (k in names(dc$obs_cov_files)) {
      oc <- occ_read_csv(dir, dc$obs_cov_files[[k]])
      oc <- oc[match(sites, as.character(oc[[site_col]])), setdiff(names(oc), site_col), drop = FALSE]
      if (ncol(oc) != J * T) stop(sprintf("survey covariate file for '%s' has %d columns but there are %d detection columns",
                                          k, ncol(oc), J * T))
      obs[[k]] <- array(as.matrix(oc), c(length(sites), J, T))
    }
  } else stop("data$format must be 'long' or 'wide'")

  if (model == "static") { y <- matrix(y, length(sites), J); obs <- lapply(obs, function(a) matrix(a, length(sites), J)) }
  # site covariates
  sc <- data.frame(row.names = seq_along(sites))
  if (!is.null(sc_raw)) {
    if (!site_col %in% names(sc_raw)) stop("site covariate file has no '", site_col, "' column")
    sc_raw[[site_col]] <- as.character(sc_raw[[site_col]])
    nos <- setdiff(sites, sc_raw[[site_col]])
    if (length(nos)) issues <- c(issues, sprintf("%d site(s) with detection data have no row in the site covariate file: %s",
                                                length(nos), paste(utils::head(nos, 8), collapse = ", ")))
    extra <- setdiff(sc_raw[[site_col]], sites)
    if (length(extra)) issues <- c(issues, sprintf("note: %d site(s) in the site covariate file have no surveys and are ignored: %s",
                                                  length(extra), paste(utils::head(extra, 8), collapse = ", ")))
    use <- if (!is.null(dc$site_cov_cols)) unlist(dc$site_cov_cols) else setdiff(names(sc_raw), site_col)
    sc <- sc_raw[match(sites, sc_raw[[site_col]]), use, drop = FALSE]
    rownames(sc) <- NULL
  }
  # season (site x season) covariates
  yc <- list()
  if (model == "dynamic" && !is.null(dc$season_covs)) {
    se <- if (is.null(dc$season_col)) "season" else dc$season_col
    sr <- occ_read_csv(dir, dc$season_covs)
    for (cc in c(site_col, se)) if (!cc %in% names(sr))
      stop(sprintf("season covariate file has no '%s' column (columns: %s); it needs %s, %s and the covariates",
                   cc, paste(names(sr), collapse = ", "), site_col, se))
    sr[[site_col]] <- as.character(sr[[site_col]])
    use <- if (!is.null(dc$season_cov_cols)) unlist(dc$season_cov_cols) else setdiff(names(sr), c(site_col, se))
    ii <- match(sr[[site_col]], sites); tt <- match(as.character(sr[[se]]), as.character(season_labels))
    ok <- !is.na(ii) & !is.na(tt)
    unk <- setdiff(unique(as.character(sr[[se]][!is.na(sr[[se]])])), as.character(season_labels))
    if (length(unk)) issues <- c(issues, sprintf(
      "season covariate file has season value(s) %s that do not occur in the detection file's '%s' column (%s). Use the same season labels in both files.",
      paste(utils::head(unk, 6), collapse = ", "), se, paste(utils::head(season_labels, 6), collapse = ", ")))
    for (k in use) {
      a <- matrix(if (is.numeric(sr[[k]])) NA_real_ else NA_character_, length(sites), length(season_labels))
      a[cbind(ii[ok], tt[ok])] <- sr[[k]][ok]
      yc[[k]] <- a
    }
  }
  forms <- if (!is.null(config$formulas)) config$formulas else dc$formulas
  list(model = model, y = y, site_covs = sc, obs_covs = obs, yearly_covs = yc,
       formulas = merge_list(occ_default_formulas, forms), site_ids = sites,
       season_labels = if (is.null(season_labels)) NULL else as.character(season_labels),
       read_issues = issues)
}

occ_check_data <- function(data) {
  pr <- grep("^note:", data$read_issues, value = TRUE, invert = TRUE)
  y <- data$y; dm <- occ_dims(data); M <- dm[["M"]]; J <- dm[["J"]]; T <- dm[["T"]]
  ids <- if (is.null(data$site_ids)) seq_len(M) else data$site_ids
  y3 <- array(y, c(M, J, T))
  bad <- !is.na(y3) & !(y3 %in% c(0, 1))
  if (any(bad)) {
    k <- which(bad, arr.ind = TRUE)[1, ]
    pr <- c(pr, sprintf("detections must be 0/1 (or blank for no survey); %d value(s) are not, e.g. %s at site %s survey %d. If these are counts, convert to 1 = detected.",
                        sum(bad), format(y3[bad][1]), ids[k[1]], k[2]))
  }
  if (J < 2) pr <- c(pr, "only one survey per site: occupancy and detection cannot be separated. At least 2 (preferably 3+) repeat surveys are needed.")
  ns <- apply(!is.na(y3), c(1, 3), sum)
  if (all(ns <= 1)) pr <- c(pr, "no site has more than one survey in a season; repeat surveys are required to estimate detection.")
  none <- which(rowSums(ns) == 0)
  if (length(none)) pr <- c(pr, sprintf("%d site(s) have no surveys at all: %s. Remove them.", length(none), paste(utils::head(ids[none], 8), collapse = ", ")))
  if (sum(y3, na.rm = TRUE) == 0) pr <- c(pr, "the species was never detected: occupancy cannot be estimated.")
  for (k in names(data$site_covs)) {
    v <- data$site_covs[[k]]
    if (anyNA(v)) pr <- c(pr, sprintf("site covariate '%s' is missing at %d site(s), e.g. %s. Fill in or drop the covariate.",
                                      k, sum(is.na(v)), paste(utils::head(ids[is.na(v)], 5), collapse = ", ")))
    else if (length(unique(v)) < 2) pr <- c(pr, sprintf("site covariate '%s' has the same value at every site.", k))
  }
  for (k in names(data$obs_covs)) {
    a <- array(data$obs_covs[[k]], c(M, J, T))
    m <- !is.na(y3) & is.na(a)
    if (any(m)) { w <- which(m, arr.ind = TRUE)[1, ]
      pr <- c(pr, sprintf("survey covariate '%s' is missing for %d survey(s) that have detection data, e.g. site %s survey %d%s. Fill it in or blank the detection for that survey.",
                          k, sum(m), ids[w[1]], w[2], if (T > 1) sprintf(" season %d", w[3]) else "")) }
  }
  for (k in names(data$yearly_covs)) {
    a <- data$yearly_covs[[k]][, seq_len(max(T - 1, 1)), drop = FALSE]
    if (anyNA(a)) pr <- c(pr, sprintf("season covariate '%s' is missing for %d site-season(s) (needed for seasons 1..%d).", k, sum(is.na(a)), T - 1))
  }
  if (data$model == "dynamic" && T < 2) pr <- c(pr, "dynamic model needs at least 2 seasons.")
  # Text covariates that are really numbers ("39.8%", "1,200", "12 m") or that have
  # so many distinct values that they cannot be categories. Either would silently
  # become a factor with one coefficient per value.
  f <- merge_list(occ_default_formulas, data$formulas)
  used <- unique(unlist(lapply(f, function(x) all.vars(stats::as.formula(x)))))
  covs <- c(as.list(data$site_covs), data$obs_covs, data$yearly_covs)
  for (k in intersect(used, names(covs))) {
    v <- covs[[k]]; if (is.numeric(v)) next
    v <- trimws(as.character(v)); v <- v[!is.na(v)]; if (!length(v)) next
    numlike <- grepl("^[-+]?[0-9][0-9.,]*[ ]*(%|[A-Za-z]{1,3}[.]?)?$", v)  # 39.8%, 1,200, 12 m, 5 km
    num <- suppressWarnings(as.numeric(gsub("[^0-9.+-]", "", v)))
    nlev <- length(unique(v))
    if (mean(numlike) >= 0.8 && any(numlike & is.na(suppressWarnings(as.numeric(v)))))
      pr <- c(pr, sprintf("covariate '%s' is stored as text, not numbers (e.g. '%s'), probably because of units, symbols or commas. Remove them so the column holds plain numbers (e.g. '%s'); otherwise every value becomes its own category.",
                          k, v[numlike & is.na(suppressWarnings(as.numeric(v)))][1], format(num[numlike & is.na(suppressWarnings(as.numeric(v)))][1])))
    else if (nlev > max(12, M / 4))
      pr <- c(pr, sprintf("text covariate '%s' has %d different values; a category covariate needs few levels with several sites each. Check that it is not an ID, date or number column.", k, nlev))
  }
  # formula variables exist
  av <- c(names(data$site_covs), names(data$obs_covs), names(data$yearly_covs))
  for (s in names(f)) {
    v <- all.vars(stats::as.formula(f[[s]])); miss <- setdiff(v, av)
    if (length(miss)) pr <- c(pr, sprintf("formula for %s uses %s, which is not a column in the covariate files (available: %s)",
                                          s, paste(miss, collapse = ", "), paste(av, collapse = ", ")))
  }
  pr
}
