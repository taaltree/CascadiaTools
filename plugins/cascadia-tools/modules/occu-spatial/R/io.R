# occu-spatial: reading researchers' files and checking them.
#
# Files (CSV):
#   detections: one row per survey actually done (long format)
#     site, survey, detected (0/1) [, species] [, source] [, survey covariates]
#   sites: one row per site: site, x, y (projected coordinates, m or km), site covariates
#   prediction (optional): x, y and the same site covariates on a grid -> maps
# The data type follows the columns: species_col given -> multi-species,
# source_col given -> integrated (several survey methods), else single species.

`%||%` <- function(a, b) if (is.null(a)) b else a

# model terms from a formula, numeric terms standardised (mean 0, SD 1)
build_terms <- function(formula, df, scaling = NULL, xlev = NULL, what = "data") {
  f <- stats::as.formula(formula)
  miss <- setdiff(all.vars(f), names(df))
  if (length(miss))
    stop(sprintf("the formula '%s' uses %s, but the %s file has no such column. Its columns are: %s. Fix the spelling in the config (formulas:) or rename the column.",
                 paste(deparse(f), collapse = ""), paste0("'", miss, "'", collapse = ", "), what,
                 paste(names(df), collapse = ", ")), call. = FALSE)
  mf <- stats::model.frame(f, df, na.action = stats::na.pass, xlev = xlev)
  mm <- stats::model.matrix(f, mf)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  raw <- colnames(mm); nm <- clean_name(raw)
  colnames(mm) <- nm
  new_scaling <- is.null(scaling)
  if (new_scaling) scaling <- list()
  for (k in seq_along(nm)) {
    v <- mm[, k]
    binary <- all(v[!is.na(v)] %in% c(0, 1))
    if (new_scaling) scaling[[nm[k]]] <- if (binary) list(center = 0, scale = 1, label = raw[k], binary = TRUE)
                                         else list(center = mean(v, na.rm = TRUE), scale = stats::sd(v, na.rm = TRUE), label = raw[k], binary = FALSE)
    sc <- scaling[[nm[k]]]
    mm[, k] <- (v - sc$center) / sc$scale
  }
  xl <- lapply(mf, function(c) if (is.factor(c) || is.character(c)) levels(factor(c)))
  list(X = as.data.frame(mm), terms = nm, scaling = scaling,
       formula = if (length(nm)) paste("~", paste(nm, collapse = " + ")) else "~ 1",
       xlev = Filter(Negate(is.null), xl))
}

survey_index <- function(site, survey) {
  # within-site order of survey labels -> 1..K
  o <- stats::ave(seq_along(site), site, FUN = function(i) rank(survey[i], ties.method = "first"))
  as.integer(o)
}

det_matrices <- function(rows, site_idx, k_idx, J, K, p_formula) {
  if (is.null(p_formula) || identical(trimws(p_formula), "~ 1") || identical(trimws(p_formula), "~1"))
    return(list(X = list(), terms = character(), scaling = list(), formula = "~ 1"))
  bt <- build_terms(p_formula, rows, what = "detections")
  X <- lapply(bt$X, function(v) { m <- matrix(NA_real_, J, K); m[cbind(site_idx, k_idx)] <- v; m })
  list(X = X, terms = bt$terms, scaling = bt$scaling, formula = bt$formula)
}

read_occu_spatial <- function(config) {
  dd <- config$data; dir <- config$.dir %||% "."
  rd <- function(f) {
    fp <- if (grepl("^(/|~|[A-Za-z]:)", f)) path.expand(f) else file.path(dir, f)
    if (!file.exists(fp)) stop(sprintf("data file '%s' not found (looked for %s). Paths in the config are relative to the config file's folder unless absolute.", f, fp), call. = FALSE)
    utils::read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  }
  det <- rd(dd$detections); st <- rd(dd$sites)
  sc <- dd$site_col %||% "site"; vc <- dd$survey_col %||% "survey"; yc <- dd$detection_col %||% "detected"
  xc <- dd$x_col %||% "x"; ycol <- dd$y_col %||% "y"
  spc <- dd$species_col; src <- dd$source_col
  keys <- c(site_col = sc, survey_col = vc, detection_col = yc, species_col = spc, source_col = src)
  for (k in names(keys)) if (!keys[[k]] %in% names(det))
    stop(sprintf("the detections file (%s) has no column '%s'. Its columns are: %s. Set 'data: %s:' in the config to the right column name.",
                 dd$detections, keys[[k]], paste(names(det), collapse = ", "), k), call. = FALSE)
  keys <- c(site_col = sc, x_col = xc, y_col = ycol)
  for (k in names(keys)) if (!keys[[k]] %in% names(st))
    stop(sprintf("the sites file (%s) has no column '%s'. Its columns are: %s. Set 'data: %s:' in the config to the right column name.",
                 dd$sites, keys[[k]], paste(names(st), collapse = ", "), k), call. = FALSE)
  fac <- if (identical(tolower(dd$coord_units %||% "m"), "m")) 1000 else 1
  det[[sc]] <- as.character(det[[sc]]); st[[sc]] <- as.character(st[[sc]])
  type <- if (!is.null(spc)) "multi" else if (!is.null(src)) "integrated" else "single"
  # rows with an empty detection value are dropped (treated as surveys not done), and the user is told
  na_det <- is.na(det[[yc]]) | trimws(as.character(det[[yc]])) == ""
  read_notes <- character()
  if (any(na_det)) {
    ex <- utils::head(paste(det[[sc]][na_det], det[[vc]][na_det]), 5)
    read_notes <- c(read_notes, sprintf("%d detection row(s) had an empty/NA detection value (e.g. %s) and were treated as surveys NOT done. If the survey was done, enter 0 or 1.",
                                        sum(na_det), paste(ex, collapse = "; ")))
    det <- det[!na_det, , drop = FALSE]
  }
  dkey <- do.call(paste, c(det[c(sc, vc, spc, src)], sep = " / "))
  dup_rows <- unique(dkey[duplicated(dkey)])
  dup_sites <- st[[sc]][duplicated(st[[sc]])]
  raw_xy <- cbind(st[[xc]], st[[ycol]])
  missing_sites <- setdiff(unique(det[[sc]]), st[[sc]])
  if (length(missing_sites)) {
    sugg <- vapply(utils::head(missing_sites, 10), function(s) {
      dd_ <- utils::adist(s, st[[sc]])[1, ]
      if (min(dd_) <= 2) sprintf("'%s' (did you mean '%s'?)", s, st[[sc]][which.min(dd_)]) else sprintf("'%s'", s) }, "")
    stop(sprintf("%d site id(s) appear in the detections file but not in the sites file: %s. Site ids must match exactly in both files (check spelling, leading zeros and spaces).",
                 length(missing_sites), paste(sugg, collapse = ", ")), call. = FALSE)
  }
  st_keep <- st[[sc]] %in% det[[sc]]
  if (any(!st_keep)) read_notes <- c(read_notes, sprintf("%d site(s) in the sites file have no detection rows and were left out (e.g. %s).",
                                                       sum(!st_keep), paste(utils::head(st[[sc]][!st_keep], 5), collapse = ", ")))
  st <- st[st_keep, , drop = FALSE]
  J <- nrow(st); site_ids <- st[[sc]]
  coords <- cbind(x = st[[xc]] / fac, y = st[[ycol]] / fac)
  psi_f <- config$formulas$psi %||% "~ 1"
  occ <- build_terms(psi_f, st, what = "sites")
  out <- list(type = type, coords = coords, occ_X = occ$X, site_ids = site_ids,
              psi_terms = occ$terms, psi_formula = occ$formula, scaling = occ$scaling,
              site_table = data.frame(site = site_ids, x = coords[, 1], y = coords[, 2], st, check.names = FALSE),
              coord_units = "km", simulated = FALSE, raw_psi_formula = psi_f, notes = read_notes,
              dup_rows = dup_rows, dup_sites = dup_sites, raw_xy = raw_xy[st_keep, , drop = FALSE],
              coord_units_raw = tolower(dd$coord_units %||% "m"))
  p_f <- config$formulas$p %||% "~ 1"
  if (type %in% c("single", "multi")) {
    sv <- unique(det[, c(sc, vc), drop = FALSE])
    sv$j <- match(sv[[sc]], site_ids); sv$k <- survey_index(sv[[sc]], sv[[vc]])
    K <- max(sv$k)
    key <- paste(det[[sc]], det[[vc]], sep = "\r")
    svkey <- paste(sv[[sc]], sv[[vc]], sep = "\r")
    # survey covariates: one row per site-survey (first row if repeated per species)
    rows <- det[!duplicated(key), , drop = FALSE]
    rows <- rows[match(svkey, key[!duplicated(key)]), , drop = FALSE]
    dm <- det_matrices(rows, sv$j, sv$k, J, K, if (is.list(p_f)) p_f[[1]] else p_f)
    out$det_X <- dm$X; out$p_terms <- dm$terms; out$p_formula <- dm$formula; out$p_scaling <- dm$scaling
    out$raw_p_formula <- p_f
    if (type == "single") {
      y <- matrix(NA_real_, J, K)
      y[cbind(match(det[[sc]], site_ids), sv$k[match(key, svkey)])] <- det[[yc]]
      out$y <- y
    } else {
      sp <- sort(unique(as.character(det[[spc]])))
      y <- array(NA_real_, c(length(sp), J, K), dimnames = list(sp, NULL, NULL))
      # a surveyed site-survey with no row for a species = not detected
      for (i in seq_along(sp)) { m <- matrix(NA_real_, J, K); m[cbind(sv$j, sv$k)] <- 0; y[i, , ] <- m }
      y[cbind(match(det[[spc]], sp), match(det[[sc]], site_ids), sv$k[match(key, svkey)])] <- det[[yc]]
      n_imp <- length(sp) * nrow(sv) - nrow(unique(det[, c(spc, sc, vc)]))
      if (n_imp > 0) out$notes <- c(out$notes, sprintf("%d species x survey combinations had no row and were treated as 'not detected'.", n_imp))
      out$y <- y; out$species <- sp
    }
  } else {
    sources <- unique(as.character(det[[src]]))
    out$sources <- sources; out$y <- out$det_X <- out$p_terms <- out$p_formula <- out$sites <- out$p_scaling <- list()
    for (s in sources) {
      d <- det[det[[src]] == s, , drop = FALSE]
      ss <- unique(d[[sc]]); ss <- site_ids[site_ids %in% ss]
      j <- match(d[[sc]], ss); k <- survey_index(d[[sc]], d[[vc]]); K <- max(k)
      y <- matrix(NA_real_, length(ss), K); y[cbind(j, k)] <- d[[yc]]
      pf <- if (is.list(p_f)) (p_f[[s]] %||% "~ 1") else p_f
      dm <- det_matrices(d, j, k, length(ss), K, pf)
      out$y[[s]] <- y; out$det_X[[s]] <- dm$X; out$p_terms[[s]] <- dm$terms
      out$p_formula[[s]] <- dm$formula; out$p_scaling[[s]] <- dm$scaling
      out$sites[[s]] <- match(ss, site_ids)
    }
    out$raw_p_formula <- p_f
  }
  if (!is.null(dd$prediction)) {
    pg <- rd(dd$prediction)
    pt <- build_terms(psi_f, pg, scaling = occ$scaling, xlev = occ$xlev, what = "prediction")
    out$pred <- data.frame(x = pg[[xc]] / fac, y = pg[[ycol]] / fac, pt$X, check.names = FALSE)
    if (anyNA(out$pred)) { out$notes <- c(out$notes, sprintf("%d prediction cells with missing covariates dropped.", sum(!stats::complete.cases(out$pred))))
                           out$pred <- out$pred[stats::complete.cases(out$pred), ] }
  }
  out$labels <- config$figures$labels
  out
}

check_occu_spatial <- function(data) {
  pr <- character()
  ys <- if (data$type == "integrated") data$y else list(all = data$y)
  for (nm in names(ys)) {
    v <- as.vector(ys[[nm]])
    bad <- v[!is.na(v) & !v %in% c(0, 1)]
    if (length(bad)) pr <- c(pr, sprintf("detections%s must be 0 or 1; found %s", if (nm == "all") "" else paste0(" (", nm, ")"),
                                         paste(utils::head(unique(bad), 5), collapse = ", ")))
  }
  if (anyNA(data$occ_X)) {
    bad <- data$site_ids[!stats::complete.cases(data$occ_X)]
    pr <- c(pr, sprintf("site covariates missing at %d site(s) (e.g. %s); fill them in or drop those sites.", length(bad),
                        paste(utils::head(bad, 5), collapse = ", ")))
  }
  if (length(data$dup_rows))
    pr <- c(pr, sprintf("%d site/survey combination(s) appear more than once in the detections file (e.g. %s). Each survey must be one row%s; remove or merge the duplicates.",
                        length(data$dup_rows), paste(utils::head(data$dup_rows, 3), collapse = "; "),
                        if (data$type == "multi") " per species" else if (data$type == "integrated") " per method" else ""))
  if (length(data$dup_sites))
    pr <- c(pr, sprintf("site id(s) listed more than once in the sites file: %s.", paste(utils::head(unique(data$dup_sites), 5), collapse = ", ")))
  if (!is.null(data$raw_xy) && !anyNA(data$raw_xy) && all(abs(data$raw_xy[, 1]) <= 180) && all(abs(data$raw_xy[, 2]) <= 90) &&
      !(identical(data$coord_units_raw, "km") && diff(range(data$raw_xy[, 1])) > 5))
    pr <- c(pr, sprintf("the site coordinates look like longitude/latitude in degrees (x from %.3f to %.3f, y from %.3f to %.3f). Distances must be in projected units: convert to UTM or another projected system (metres) and set 'coord_units: m'.",
                        min(data$raw_xy[, 1]), max(data$raw_xy[, 1]), min(data$raw_xy[, 2]), max(data$raw_xy[, 2])))
  if (anyNA(data$coords)) pr <- c(pr, "some sites have no x/y coordinates.")
  else if (any(duplicated(round(data$coords, 6))))
    pr <- c(pr, sprintf("%d sites share identical coordinates with another site; the spatial models need distinct locations (merge them or jitter by a few metres).",
                        sum(duplicated(round(data$coords, 6)))))
  chk_det <- function(X, y, lab) {
    out <- character()
    for (n in names(X)) {
      m <- X[[n]]; y2 <- if (length(dim(y)) == 3) y[1, , ] else y
      bad <- is.na(m) & !is.na(y2)
      if (any(bad)) out <- c(out, sprintf("survey covariate '%s'%s is missing for %d survey(s) that were done (e.g. site %s); fill it in, or blank the detection value if the survey was not done.",
                                          n, lab, sum(bad), paste(utils::head(unique(sites_lab[which(bad, arr.ind = TRUE)[, 1]]), 3), collapse = ", ")))
    }
    out
  }
  if (data$type == "integrated") for (s in data$sources) { sites_lab <- data$site_ids[data$sites[[s]]]; pr <- c(pr, chk_det(data$det_X[[s]], data$y[[s]], paste0(" (", s, ")"))) }
  else { sites_lab <- data$site_ids; pr <- c(pr, chk_det(data$det_X, data$y, "")) }
  J <- nrow(data$coords)
  if (J < 20) pr <- c(pr, sprintf("only %d sites: too few to estimate occupancy covariate effects or spatial structure.", J))
  K <- if (data$type == "single") ncol(data$y) else if (data$type == "multi") dim(data$y)[3] else max(vapply(data$y, ncol, 1L))
  if (K < 2) pr <- c(pr, "every site has a single survey: detection cannot be separated from occupancy (need >= 2 surveys at some sites).")
  pr
}
