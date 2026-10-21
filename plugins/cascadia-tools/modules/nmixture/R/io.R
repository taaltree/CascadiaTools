# nmixture: reading researchers' count files and checking them.
#
# Long format (default; one row per visit actually done):
#   site, [season], visit, count, [visit covariates...]
# Wide format (one row per site; count columns in visit order, blank = not visited):
#   site, <count columns...>; visit covariates as extra wide files (same layout).
# Site covariates: one row per site. Season covariates (dynamic): site, season, covariates.

nm_read_csv <- function(dir, f) {
  if (is.null(f)) return(NULL)
  p <- if (file.exists(f)) f else file.path(dir, f)
  if (!file.exists(p)) stop("file not found: ", f, " (looked in ", dir, ")")
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
}

# " (did you mean 'Count'?)" when a missing column differs only in case/spacing from a present one
nm_col_hint <- function(miss, present) {
  key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  h <- unlist(lapply(miss, function(m) { k <- present[key(present) == key(m)]
    if (length(k)) sprintf("'%s' (did you mean '%s'?)", m, k[1]) }))
  if (length(h)) paste0("; ", paste(h, collapse = ", ")) else ""
}

# " (closest in the covariate file: 'PT-007')" for a probably misspelled site name
nm_near_hint <- function(s, pool) {
  pool <- unique(pool[!is.na(pool)]); if (!length(pool)) return("")
  d <- utils::adist(tolower(s), tolower(pool))[1, ]
  if (min(d) <= 2) sprintf(" (closest name in the site covariate file: '%s')", pool[which.min(d)]) else ""
}

# Text in a covariate the formulas use: numbers with symbols ("55%"), a few words among
# numbers ("calm"), or very many categories. Returns a message or NULL.
nm_text_cov_issue <- function(v, name, kind, n_units) {
  if (!is.character(v)) return(NULL)
  vv <- v[!is.na(v) & nzchar(trimws(v))]; if (!length(vv)) return(NULL)
  num <- suppressWarnings(as.numeric(vv))
  strp <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", vv)))
  nl <- length(unique(vv))
  if (mean(!is.na(num)) >= 0.5) {
    bad <- unique(vv[is.na(num)])
    return(sprintf("%s covariate '%s' is mostly numbers but %d value(s) are text (e.g. '%s'), so the whole column would be treated as a category. Replace the text with numbers, or leave the cell blank if the value is unknown.",
                   kind, name, sum(is.na(num)), bad[1]))
  }
  if (mean(!is.na(strp)) >= 0.9 && nl > 5)
    return(sprintf("%s covariate '%s' holds numbers with symbols or units (e.g. '%s'), so it would be treated as a category with %d levels. Remove the symbols so the column holds plain numbers.",
                   kind, name, vv[1], nl))
  if (nl > max(10, n_units / 2))
    return(sprintf("%s covariate '%s' is text with %d different values; as a category it needs one parameter per value. If it is a measurement, store it as plain numbers; otherwise group the values into fewer categories.",
                   kind, name, nl))
  NULL
}

nm_read_data <- function(config) {
  dc <- config$data; dir <- config$.dir %||% "."
  if (is.null(dc)) stop("config needs a 'data:' section")
  model <- dc$model %||% config$settings$model %||% "static"
  model <- match.arg(model, c("static", "dynamic"))
  fmt <- dc$format %||% "long"
  site_col <- dc$site_col %||% "site"
  issues <- character()
  cnt <- nm_read_csv(dir, dc$counts)
  if (is.null(cnt)) stop("config data$counts (the count file) is required")
  if (!site_col %in% names(cnt)) stop(sprintf("count file has no '%s' column (columns: %s); set data$site_col",
                                              site_col, paste(names(cnt), collapse = ", ")))
  cnt[[site_col]] <- as.character(cnt[[site_col]])
  sites <- sort(unique(cnt[[site_col]]))
  obs_names <- unlist(dc$obs_covs)
  obs <- list(); season_labels <- NULL
  if (fmt == "long") {
    vv <- dc$visit_col %||% "visit"; cv <- dc$count_col %||% "count"; se <- dc$season_col %||% "season"
    need <- c(vv, cv, if (model == "dynamic") se, obs_names)
    miss <- setdiff(need, names(cnt))
    if (length(miss)) stop("count file is missing column(s): ", paste(miss, collapse = ", "),
                           nm_col_hint(miss, names(cnt)),
                           " (columns present: ", paste(names(cnt), collapse = ", "), ")")
    vn <- suppressWarnings(as.numeric(cnt[[vv]]))
    if (anyNA(vn) || any(vn != round(vn)) || any(vn < 1))
      stop(sprintf("visit column '%s' must hold visit numbers 1, 2, 3, ... within each site%s (found e.g. '%s')",
                   vv, if (model == "dynamic") " and season" else "", cnt[[vv]][which(is.na(vn) | vn != round(vn) | vn < 1)[1]]))
    J <- max(vn)
    if (model == "dynamic") { season_labels <- sort(unique(cnt[[se]])); tt <- match(cnt[[se]], season_labels) }
    else tt <- rep(1L, nrow(cnt))
    T <- max(tt); ii <- match(cnt[[site_col]], sites)
    key <- paste(ii, vn, tt)
    if (any(duplicated(key))) {
      d <- which(duplicated(key))[1]
      issues <- c(issues, sprintf("%d duplicated visit rows (same site, %svisit), e.g. site %s visit %s. Each visit must appear once.",
                                  sum(duplicated(key)), if (model == "dynamic") "season, " else "", cnt[[site_col]][d], cnt[[vv]][d]))
    }
    yv <- suppressWarnings(as.numeric(cnt[[cv]]))
    bad <- which(is.na(yv) & !is.na(cnt[[cv]]) & nzchar(as.character(cnt[[cv]])))
    if (length(bad)) issues <- c(issues, sprintf("count column '%s' has %d non-numeric value(s), e.g. '%s' (row %d).", cv, length(bad), cnt[[cv]][bad[1]], bad[1] + 1))
    y <- array(NA_real_, c(length(sites), J, T)); y[cbind(ii, vn, tt)] <- yv
    for (k in obs_names) {
      v <- cnt[[k]]
      a <- array(if (is.numeric(v)) NA_real_ else NA_character_, c(length(sites), J, T))
      a[cbind(ii, vn, tt)] <- v
      obs[[k]] <- a
    }
  } else if (fmt == "wide") {
    ccols <- if (!is.null(dc$count_cols)) unlist(dc$count_cols) else setdiff(names(cnt), site_col)
    miss <- setdiff(ccols, names(cnt)); if (length(miss)) stop("count columns not found: ", paste(miss, collapse = ", "))
    cnt <- cnt[match(sites, cnt[[site_col]]), , drop = FALSE]
    ym <- as.matrix(cnt[, ccols, drop = FALSE]); storage.mode(ym) <- "numeric"
    if (model == "dynamic") {
      J <- dc$visits_per_season
      if (is.null(J)) stop("wide format for a dynamic model needs data$visits_per_season")
      T <- length(ccols) / J
      if (T != round(T)) stop("number of count columns is not a multiple of visits_per_season")
      season_labels <- if (!is.null(dc$season_labels)) unlist(dc$season_labels) else seq_len(T)
    } else { J <- length(ccols); T <- 1 }
    y <- array(ym, c(length(sites), J, T))
    for (k in names(dc$obs_cov_files)) {
      oc <- nm_read_csv(dir, dc$obs_cov_files[[k]])
      oc <- oc[match(sites, as.character(oc[[site_col]])), setdiff(names(oc), site_col), drop = FALSE]
      if (ncol(oc) != J * T) stop(sprintf("visit covariate file for '%s' has %d columns but there are %d count columns", k, ncol(oc), J * T))
      obs[[k]] <- array(as.matrix(oc), c(length(sites), J, T))
    }
  } else stop("data$format must be 'long' or 'wide'")
  if (model == "static") { y <- matrix(y, length(sites), J); obs <- lapply(obs, function(a) matrix(a, length(sites), J)) }
  sc <- data.frame(row.names = seq_along(sites))
  sc_raw <- nm_read_csv(dir, dc$site_covs)
  if (!is.null(sc_raw)) {
    if (!site_col %in% names(sc_raw)) stop(sprintf("site covariate file has no '%s' column%s (columns present: %s). The site column must have the same name in the count file and the site covariate file.",
                                                   site_col, nm_col_hint(site_col, names(sc_raw)), paste(names(sc_raw), collapse = ", ")))
    sc_raw[[site_col]] <- as.character(sc_raw[[site_col]])
    nos <- setdiff(sites, sc_raw[[site_col]])
    extra <- setdiff(sc_raw[[site_col]], sites)
    if (length(nos)) issues <- c(issues, sprintf("%d site(s) with counts have no row in the site covariate file: %s. Check the spelling of these site names in the two files.",
                                                length(nos), paste(vapply(utils::head(nos, 8), function(s) paste0(s, nm_near_hint(s, sc_raw[[site_col]])), ""), collapse = ", ")))
    if (length(extra)) issues <- c(issues, sprintf("note: %d site(s) in the site covariate file have no counts and are ignored: %s",
                                                  length(extra), paste(utils::head(extra, 8), collapse = ", ")))
    use <- if (!is.null(dc$site_cov_cols)) unlist(dc$site_cov_cols) else setdiff(names(sc_raw), site_col)
    sc <- sc_raw[match(sites, sc_raw[[site_col]]), use, drop = FALSE]; rownames(sc) <- NULL
  }
  yc <- list()
  if (model == "dynamic" && !is.null(dc$season_covs)) {
    se <- dc$season_col %||% "season"
    sr <- nm_read_csv(dir, dc$season_covs)
    miss <- setdiff(c(site_col, se), names(sr))
    if (length(miss)) stop(sprintf("season covariate file is missing column(s) %s%s (columns present: %s). It needs the same site and season columns as the count file.",
                                   paste(miss, collapse = ", "), nm_col_hint(miss, names(sr)), paste(names(sr), collapse = ", ")))
    sr[[site_col]] <- as.character(sr[[site_col]])
    use <- if (!is.null(dc$season_cov_cols)) unlist(dc$season_cov_cols) else setdiff(names(sr), c(site_col, se))
    ii <- match(sr[[site_col]], sites); tt <- match(sr[[se]], season_labels); ok <- !is.na(ii) & !is.na(tt)
    for (k in use) {
      a <- matrix(if (is.numeric(sr[[k]])) NA_real_ else NA_character_, length(sites), length(season_labels))
      a[cbind(ii[ok], tt[ok])] <- sr[[k]][ok]; yc[[k]] <- a
    }
  }
  forms <- config$formulas %||% dc$formulas
  list(model = model, y = y, site_covs = sc, obs_covs = obs, yearly_covs = yc,
       formulas = merge_list(nm_default_formulas, forms), site_ids = sites,
       season_labels = if (is.null(season_labels)) NULL else as.character(season_labels),
       read_issues = issues, units = dc$units, species = dc$species,
       labels = config$figures$labels %||% dc$labels)
}

nm_check_data <- function(data) {
  pr <- grep("^note:", data$read_issues, value = TRUE, invert = TRUE)
  dm <- nm_dims(data); M <- dm[["M"]]; J <- dm[["J"]]; T <- dm[["T"]]
  ids <- data$site_ids %||% seq_len(M)
  y3 <- array(data$y, c(M, J, T))
  bad <- !is.na(y3) & (y3 < 0 | y3 != round(y3))
  if (any(bad)) { k <- which(bad, arr.ind = TRUE)[1, ]
    pr <- c(pr, sprintf("counts must be whole numbers >= 0 (blank for a visit not done); %d value(s) are not, e.g. %s at site %s visit %d.",
                        sum(bad), format(y3[bad][1]), ids[k[1]], k[2])) }
  if (J < 2) pr <- c(pr, "only one visit per site: abundance and detection cannot be separated. N-mixture models need at least 2 (preferably 3+) repeat counts during a period when the population is closed.")
  nv <- apply(!is.na(y3), c(1, 3), sum)
  if (all(nv <= 1)) pr <- c(pr, "no site has more than one visit in a season; repeat counts are required to estimate detection.")
  none <- which(rowSums(nv) == 0)
  if (length(none)) pr <- c(pr, sprintf("%d site(s) have no visits at all: %s. Remove them.", length(none), paste(utils::head(ids[none], 8), collapse = ", ")))
  if (sum(y3, na.rm = TRUE) == 0) pr <- c(pr, "the species was never counted: abundance cannot be estimated.")
  for (k in names(data$site_covs)) {
    v <- data$site_covs[[k]]
    if (anyNA(v)) pr <- c(pr, sprintf("site covariate '%s' is missing at %d site(s), e.g. %s. Fill in or drop the covariate.",
                                      k, sum(is.na(v)), paste(utils::head(ids[is.na(v)], 5), collapse = ", ")))
    else if (length(unique(v)) < 2) pr <- c(pr, sprintf("site covariate '%s' has the same value at every site.", k))
  }
  f <- merge_list(nm_default_formulas, data$formulas)
  pvars <- all.vars(stats::as.formula(f$p))
  for (k in intersect(names(data$obs_covs), pvars)) {
    a <- array(data$obs_covs[[k]], c(M, J, T)); m <- !is.na(y3) & is.na(a)
    if (any(m)) { w <- which(m, arr.ind = TRUE)[1, ]
      pr <- c(pr, sprintf("visit covariate '%s' is missing for %d visit(s) that have counts, e.g. site %s visit %d%s. Fill it in or blank the count for that visit.",
                          k, sum(m), ids[w[1]], w[2], if (T > 1) sprintf(" season %d", w[3]) else "")) }
  }
  for (k in names(data$yearly_covs)) {
    a <- data$yearly_covs[[k]][, seq_len(max(T - 1, 1)), drop = FALSE]
    if (anyNA(a)) pr <- c(pr, sprintf("season covariate '%s' is missing for %d site-season(s) (needed for seasons 1..%d).", k, sum(is.na(a)), T - 1))
  }
  if ((data$model %||% "static") == "dynamic" && T < 2) pr <- c(pr, "dynamic model needs at least 2 seasons.")
  av <- c(names(data$site_covs), names(data$obs_covs), names(data$yearly_covs))
  use <- if ((data$model %||% "static") == "dynamic") names(f) else c("lambda", "p")
  for (s in use) {
    v <- all.vars(stats::as.formula(f[[s]])); miss <- setdiff(v, av)
    if (length(miss)) pr <- c(pr, sprintf("formula for %s uses %s, which is not a column in the covariate files (available: %s)",
                                          s, paste(miss, collapse = ", "), paste(av, collapse = ", ")))
    if (s %in% c("lambda", "gamma", "omega")) {
      bad <- intersect(v, names(data$obs_covs))
      if (length(bad)) pr <- c(pr, sprintf("formula for %s uses visit-level covariate(s) %s; abundance can only depend on site (or season) covariates.", s, paste(bad, collapse = ", ")))
    }
  }
  fv <- unique(unlist(lapply(f[use], function(x) all.vars(stats::as.formula(x)))))
  for (k in intersect(fv, names(data$site_covs))) pr <- c(pr, nm_text_cov_issue(data$site_covs[[k]], k, "site", M))
  for (k in intersect(fv, names(data$obs_covs))) pr <- c(pr, nm_text_cov_issue(as.vector(data$obs_covs[[k]]), k, "visit", sum(!is.na(y3))))
  for (k in intersect(fv, names(data$yearly_covs))) pr <- c(pr, nm_text_cov_issue(as.vector(data$yearly_covs[[k]]), k, "season", M * T))
  ymax <- max(y3, na.rm = TRUE)
  if (is.finite(ymax) && ymax > 300) pr <- c(pr, sprintf("note: largest count is %d; computation scales with the upper bound K, so the fit may be slow.", ymax))
  grep("^note:", pr, value = TRUE, invert = TRUE)
}
