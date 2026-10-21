# abund-interactions: reading researchers' files and checking them.
#
# Two input formats (config data$format):
#
# records (default; what camera-trap software exports)
#   deployments: one row per camera deployment
#       site, start_date, end_date   (dates as YYYY-MM-DD; a site may have several deployments)
#   records: one row per detection (photo / sequence / independent event)
#       site, date (or date-time), species
#   Deployment days are cut into occasions of `occasion_days` days starting at
#   each site's first day; K = number of days the camera ran in the occasion,
#   y = number of those days with >= 1 record of the species.
#
# occasions (already summarised): one row per site x occasion
#   site, occasion (1, 2, ...), effort_days, <dominant detection days>, <subordinate detection days>,
#   [occasion covariates...]
#
# Site covariates: one row per site (site, covariates..., optional x / y coordinates).

abi_csv <- function(dir, f) {
  if (is.null(f)) return(NULL)
  p <- if (file.exists(f)) f else file.path(dir, f)
  if (!file.exists(p)) stop("file not found: ", f, " (looked in ", dir, ")")
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
}

abi_need_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(sprintf("%s is missing column(s): %s (columns present: %s). Set the column names in the config data: section.",
                                 what, paste(miss, collapse = ", "), paste(names(df), collapse = ", ")))
}

abi_read_data <- function(config) {
  dc <- config$data; dir <- if (is.null(config$.dir)) "." else config$.dir
  if (is.null(dc)) stop("config needs a 'data:' section")
  g <- function(k, d) if (is.null(dc[[k]])) d else dc[[k]]
  site_col <- g("site_col", "site")
  dom <- dc$dominant; sub <- dc$subordinate
  if (is.null(dom) || is.null(sub)) stop("config data$dominant and data$subordinate (species names as written in your data) are required")
  fmt <- g("format", "records")
  warn <- character()
  obs <- list()
  if (fmt == "records") {
    dep <- abi_csv(dir, dc$deployments); rec <- abi_csv(dir, dc$records)
    if (is.null(dep) || is.null(rec)) stop("format 'records' needs data$deployments and data$records files")
    start_col <- g("start_col", "start_date"); end_col <- g("end_col", "end_date")
    date_col <- g("date_col", "date"); sp_col <- g("species_col", "species")
    abi_need_cols(dep, c(site_col, start_col, end_col), "deployment file")
    abi_need_cols(rec, c(site_col, date_col, sp_col), "records file")
    dfmt <- g("date_format", "%Y-%m-%d")
    pd <- function(x, what) {
      d <- as.Date(substr(as.character(x), 1, 10), format = dfmt)
      if (anyNA(d)) stop(sprintf("could not read %d %s value(s) as dates with format '%s' (e.g. '%s'); set data$date_format",
                                 sum(is.na(d)), what, dfmt, x[which(is.na(d))[1]]))
      d
    }
    # stray spaces in site ids / species names (common in hand-edited sheets) silently split a
    # site or species in two: trim them and say so
    for (cl in list(list("dep", site_col), list("rec", site_col), list("rec", sp_col))) {
      v <- get(cl[[1]])[[cl[[2]]]]
      if (is.character(v)) {
        tv <- trimws(v); nch <- sum(tv != v, na.rm = TRUE)
        if (nch) warn <- c(warn, sprintf("%d value(s) in column '%s' of the %s file had leading/trailing spaces; they were removed.",
                                         nch, cl[[2]], if (cl[[1]] == "dep") "deployment" else "records"))
        if (cl[[1]] == "dep") dep[[cl[[2]]]] <- tv else rec[[cl[[2]]]] <- tv
      }
    }
    nod <- is.na(rec[[date_col]]) | trimws(as.character(rec[[date_col]])) == ""
    if (any(nod)) stop(sprintf("%d record(s) have an empty date (first on line %d of the records file): fill or delete them",
                               sum(nod), which(nod)[1] + 1))
    blank <- is.na(rec[[site_col]]) | rec[[site_col]] == ""
    if (any(blank)) {
      warn <- c(warn, sprintf("%d record(s) have an empty site id (first on line %d of the records file) and were ignored.",
                              sum(blank), which(blank)[1] + 1))
      rec <- rec[!blank, , drop = FALSE]
    }
    dep$.s <- pd(dep[[start_col]], "deployment start"); dep$.e <- pd(dep[[end_col]], "deployment end")
    if (any(dep$.e < dep$.s)) stop(sprintf("deployment(s) end before they start, e.g. site %s", dep[[site_col]][which(dep$.e < dep$.s)[1]]))
    rec$.d <- pd(rec[[date_col]], "record date")
    dep[[site_col]] <- as.character(dep[[site_col]]); rec[[site_col]] <- as.character(rec[[site_col]])
    sites <- sort(unique(dep[[site_col]]))
    L <- as.integer(g("occasion_days", 7))
    # days operating per site
    days <- lapply(sites, function(s) {
      d <- dep[dep[[site_col]] == s, ]
      sort(unique(do.call(c, lapply(seq_len(nrow(d)), function(k) seq(d$.s[k], d$.e[k], by = "day")))))
    })
    nocc <- vapply(days, function(d) as.integer(ceiling((as.numeric(max(d) - min(d)) + 1) / L)), 1L)
    J <- if (!is.null(dc$n_occasions)) as.integer(dc$n_occasions) else max(nocc)
    M <- length(sites)
    K <- matrix(0, M, J); yA <- matrix(0, M, J); yS <- matrix(0, M, J)
    outside <- 0
    for (i in seq_len(M)) {
      d0 <- min(days[[i]])
      occ <- as.integer(as.numeric(days[[i]] - d0) %/% L) + 1L
      ok <- occ <= J
      tab <- tabulate(occ[ok], J); K[i, ] <- tab
      r <- rec[rec[[site_col]] == sites[i], ]
      inop <- !(r$.d %in% days[[i]])
      outside <- outside + sum(inop)
      r <- r[!inop, ]
      for (sp in c("A", "S")) {
        nm <- if (sp == "A") dom else sub
        dd <- unique(r$.d[r[[sp_col]] == nm])
        o <- as.integer(as.numeric(dd - d0) %/% L) + 1L
        o <- o[o <= J]
        cnt <- tabulate(o, J)
        if (sp == "A") yA[i, ] <- cnt else yS[i, ] <- cnt
      }
    }
    if (outside > 0) warn <- c(warn, sprintf("%d record(s) fall on days when no camera was deployed at that site and were ignored: check deployment dates.", outside))
    unk <- setdiff(unique(rec[[site_col]]), sites)
    # An unknown id that is one edit away from a deployed site is almost certainly a typo: its
    # detections would otherwise be dropped and that site read as "never detected". Stop instead.
    near <- vapply(unk, function(u) { h <- sites[utils::adist(tolower(u), tolower(sites))[1, ] <= 1]; if (length(h)) h[1] else NA_character_ }, "")
    if (any(!is.na(near)))
      stop(sprintf("site id(s) in the records file do not match any deployment but look like typos: %s. Make the site ids identical in the records and deployments files (otherwise those detections would be dropped and the site treated as never detected).",
                   paste(sprintf("'%s' (did you mean '%s'?)", unk[!is.na(near)], near[!is.na(near)]), collapse = ", ")))
    if (length(unk)) warn <- c(warn, sprintf("%d site(s) in the records file have no deployment and were ignored: %s", length(unk), paste(utils::head(unk, 6), collapse = ", ")))
    spp <- unique(rec[[sp_col]])
    for (nm in c(dom, sub)) if (!nm %in% spp)
      stop(sprintf("species '%s' never appears in column '%s' of the records file (species present: %s)", nm, sp_col, paste(utils::head(spp, 15), collapse = ", ")))
    yA[K == 0] <- NA; yS[K == 0] <- NA
  } else if (fmt == "occasions") {
    det <- abi_csv(dir, dc$occasions)
    if (is.null(det)) stop("format 'occasions' needs data$occasions (one row per site x occasion)")
    oc_col <- g("occasion_col", "occasion"); ef_col <- g("effort_col", "effort_days")
    a_col <- g("dominant_col", dom); s_col <- g("subordinate_col", sub)
    obs_names <- unlist(dc$obs_covs)
    abi_need_cols(det, c(site_col, oc_col, ef_col, a_col, s_col, obs_names), "occasion file")
    det[[site_col]] <- trimws(as.character(det[[site_col]]))
    sites <- sort(unique(det[[site_col]])); M <- length(sites)
    jj <- suppressWarnings(as.integer(det[[oc_col]]))
    if (anyNA(jj) || any(jj < 1)) stop(sprintf("occasion column '%s' must hold occasion numbers 1, 2, 3, ...", oc_col))
    J <- max(jj); ii <- match(det[[site_col]], sites)
    if (any(duplicated(paste(ii, jj)))) stop("duplicated site x occasion rows in the occasion file; each occasion must appear once per site")
    K <- matrix(0, M, J); yA <- matrix(NA_real_, M, J); yS <- matrix(NA_real_, M, J)
    K[cbind(ii, jj)] <- as.numeric(det[[ef_col]])
    yA[cbind(ii, jj)] <- as.numeric(det[[a_col]]); yS[cbind(ii, jj)] <- as.numeric(det[[s_col]])
    for (k in obs_names) { a <- matrix(if (is.numeric(det[[k]])) NA_real_ else NA_character_, M, J); a[cbind(ii, jj)] <- det[[k]]; obs[[k]] <- a }
  } else stop("data$format must be 'records' or 'occasions'")

  # site covariates
  scr <- abi_csv(dir, dc$site_covs)
  sc <- data.frame(row.names = seq_len(M)); coords <- NULL
  if (!is.null(scr)) {
    abi_need_cols(scr, site_col, "site covariate file")
    scr[[site_col]] <- trimws(as.character(scr[[site_col]]))
    if (any(duplicated(scr[[site_col]]))) stop(sprintf("site covariate file has more than one row for site(s): %s",
      paste(utils::head(unique(scr[[site_col]][duplicated(scr[[site_col]])]), 8), collapse = ", ")))
    nos <- setdiff(sites, scr[[site_col]])
    if (length(nos)) stop(sprintf("%d site(s) have survey data but no row in the site covariate file: %s", length(nos), paste(utils::head(nos, 8), collapse = ", ")))
    scr <- scr[match(sites, scr[[site_col]]), , drop = FALSE]
    cc <- unlist(dc$coords)
    if (length(cc) == 2) { abi_need_cols(scr, cc, "site covariate file"); coords <- data.frame(x = scr[[cc[1]]], y = scr[[cc[2]]]) }
    use <- if (!is.null(dc$site_cov_cols)) unlist(dc$site_cov_cols) else setdiff(names(scr), c(site_col, cc))
    abi_need_cols(scr, use, "site covariate file (site_cov_cols in the config)")
    sc <- scr[, use, drop = FALSE]; rownames(sc) <- NULL
  }
  fm <- config[["formulas"]]
  if (is.null(fm))
    warn <- c(warn, "No 'formulas:' section in the config: every submodel (dominant abundance and detection, subordinate occupancy and detection) is intercept-only. If you meant to include covariates, add them under 'formulas:' (keys lambda, r, psi, p).")
  bad <- setdiff(names(fm), c("lambda", "r", "psi", "p"))
  if (length(bad)) stop(sprintf("unknown key(s) under 'formulas:' in the config: %s (allowed: lambda, r, psi, p)", paste(bad, collapse = ", ")))
  sp_lab <- c(dominant = g("dominant_label", dom), subordinate = g("subordinate_label", sub))
  list(yA = yA, yS = yS, K = K, site_covs = sc, obs_covs = obs,
       # [[ ]] not $: fm$p would partial-match fm$psi when no p formula is given
       formulas = list(lambda = fm[["lambda"]] %||% "~ 1", r = fm[["r"]] %||% "~ 1", psi = fm[["psi"]] %||% "~ 1", p = fm[["p"]] %||% "~ 1"),
       modifier = dc$interaction_modifier, species = sp_lab, site_ids = sites, coords = coords,
       occasion_days = if (fmt == "records") L else NULL, warnings = warn)
}

abi_check_data <- function(data) {
  pr <- character()
  K <- data$K; K[is.na(K)] <- 0
  M <- nrow(K)
  if (any(K < 0) || any(K != round(K))) pr <- c(pr, "effort (days per occasion) must be whole numbers >= 0")
  for (s in c("yA", "yS")) {
    y <- data[[s]]; nm <- data$species[[if (s == "yA") "dominant" else "subordinate"]]
    w <- which(!is.na(y) & y > K, arr.ind = TRUE)
    if (length(w)) pr <- c(pr, sprintf("%s: %d occasion(s) have more detection days than days surveyed, e.g. site %s occasion %d",
                                       nm, nrow(w), data$site_ids[w[1, 1]], w[1, 2]))
    if (any(!is.na(y) & (y < 0 | y != round(y)))) pr <- c(pr, sprintf("%s: detections must be whole numbers >= 0", nm))
    if (sum(y, na.rm = TRUE) == 0) pr <- c(pr, sprintf("%s was never detected: nothing to estimate", nm))
  }
  nos <- which(rowSums(K) == 0)
  if (length(nos)) pr <- c(pr, sprintf("%d site(s) have no survey effort at all: %s", length(nos), paste(utils::head(data$site_ids[nos], 8), collapse = ", ")))
  vars <- unique(unlist(lapply(data$formulas, function(f) all.vars(stats::as.formula(f)))))
  miss <- setdiff(vars, c(names(data$site_covs), names(data$obs_covs)))
  if (length(miss)) pr <- c(pr, sprintf("formula variable(s) not found in the site or occasion covariates: %s (available: %s)",
                                        paste(miss, collapse = ", "), paste(c(names(data$site_covs), names(data$obs_covs)), collapse = ", ")))
  for (v in intersect(vars, names(data$site_covs))) if (anyNA(data$site_covs[[v]]))
    pr <- c(pr, sprintf("site covariate '%s' has %d missing value(s): fill them or drop those sites", v, sum(is.na(data$site_covs[[v]]))))
  if (!is.null(data$modifier) && !data$modifier %in% names(data$site_covs))
    pr <- c(pr, sprintf("interaction_modifier '%s' is not a column of the site covariate file", data$modifier))
  if (M < 50) pr <- c(pr, sprintf("only %d sites: abundance-mediated interactions need many sites (Twining et al. 2025 used 300-1000); results would not be reliable", M))
  pr
}
