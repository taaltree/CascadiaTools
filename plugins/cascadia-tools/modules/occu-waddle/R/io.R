# occu-waddle: reading researchers' files and checking them.
#
# config$data keys (paths relative to the config file):
#   detections      CSV with one row per site x survey.
#   format          "wide" (default): one 0/1 column per species;
#                   "long": one row per site x survey x species with a 0/1 column.
#   site_col, survey_col           site id and survey number columns
#   dominant_col, subordinate_col  (wide) detection columns of the two species
#   species_col, detection_col     (long) species name column and 0/1 column
#   dominant, subordinate          (long) species names as written in species_col;
#                                  (wide, optional) display names
#   surveyed_col    optional: rows where this is 0 or NA count as "survey not done"
#   obs_covs        survey-level covariate columns of the detections file
#   site_covs       CSV with one row per site (site_col + covariates)
# config$formulas: psiA, psiBA, psiBa, pA, pB, rA, rBA, rBa (see README).

wd_read_data <- function(config) {
  d <- config$data
  if (is.null(d)) stop("config has no 'data' section")
  rd <- function(f) utils::read.csv(file.path(config$.dir, f), stringsAsFactors = FALSE, check.names = FALSE)
  det <- rd(d$detections)
  sc_ <- if (is.null(d$site_col)) "site" else d$site_col
  sv_ <- if (is.null(d$survey_col)) "survey" else d$survey_col
  miss <- setdiff(c(sc_, sv_), names(det))
  if (length(miss)) stop("detections file lacks column(s): ", paste(miss, collapse = ", "),
                         " (set data$site_col / data$survey_col in the config). Columns in the file are: ",
                         paste(names(det), collapse = ", "))
  fmt <- if (is.null(d$format)) "wide" else d$format
  if (fmt == "long") {
    spc <- if (is.null(d$species_col)) "species" else d$species_col
    dc <- if (is.null(d$detection_col)) "detected" else d$detection_col
    if (!all(c(spc, dc) %in% names(det))) stop("long format needs columns '", spc, "' and '", dc, "'")
    if (is.null(d$dominant) || is.null(d$subordinate)) stop("long format needs data$dominant and data$subordinate (species names)")
    sp <- c(d$dominant, d$subordinate)
    nf <- setdiff(sp, det[[spc]])
    if (length(nf)) stop("species not found in column '", spc, "': ", paste(nf, collapse = ", "))
    key <- paste(det[[sc_]], det[[sv_]], sep = "\r")
    first <- !duplicated(key)
    base <- det[first, setdiff(names(det), c(spc, dc)), drop = FALSE]
    for (k in 1:2) {
      sub <- det[det[[spc]] == sp[k], ]
      base[[paste0(".y", k)]] <- sub[[dc]][match(key[first], paste(sub[[sc_]], sub[[sv_]], sep = "\r"))]
    }
    det <- base; cols <- c(".y1", ".y2"); labels <- sp
  } else {
    cols <- c(d$dominant_col, d$subordinate_col)
    if (length(cols) != 2) stop("data$dominant_col and data$subordinate_col are required (the two detection columns)")
    miss <- setdiff(cols, names(det))
    if (length(miss)) stop("detections file lacks species column(s): ", paste(miss, collapse = ", "),
                           ". Columns in the file are: ", paste(names(det), collapse = ", "),
                           ". Set data$dominant_col / data$subordinate_col in the config to the right names.")
    labels <- c(if (is.null(d$dominant)) cols[1] else d$dominant, if (is.null(d$subordinate)) cols[2] else d$subordinate)
  }
  sc <- NULL
  if (!is.null(d$site_covs)) {
    sc <- rd(d$site_covs)
    if (!sc_ %in% names(sc)) stop("site covariate file lacks the site column '", sc_, "'")
    if (anyDuplicated(sc[[sc_]])) stop("site covariate file has duplicated site ids: ",
                                       paste(unique(sc[[sc_]][duplicated(sc[[sc_]])]), collapse = ", "))
    site_ids <- as.character(sc[[sc_]])
    extra <- setdiff(unique(as.character(det[[sc_]])), site_ids)
    if (length(extra)) stop("site(s) in the detections file but not in the site covariate file: ",
                            paste(utils::head(extra, 10), collapse = ", "))
  } else site_ids <- sort(unique(as.character(det[[sc_]])))
  surveys <- sort(unique(det[[sv_]]))
  M <- length(site_ids); J <- length(surveys)
  ii <- match(as.character(det[[sc_]]), site_ids); jj <- match(det[[sv_]], surveys)
  dup <- duplicated(cbind(ii, jj))
  if (any(dup)) stop("detections file has more than one row for site ", det[[sc_]][dup][1], " survey ",
                     det[[sv_]][dup][1], " (one row per site x survey is expected)")
  done <- rep(TRUE, nrow(det))
  if (!is.null(d$surveyed_col)) {
    v <- det[[d$surveyed_col]]
    if (is.null(v)) stop("surveyed_col '", d$surveyed_col, "' not in the detections file")
    done <- !is.na(v) & v != 0
  }
  unreadable <- list()
  mk <- function(col) {
    raw <- det[[col]]
    if (is.numeric(raw) || is.logical(raw)) v <- as.numeric(raw)
    else {                                  # text: accept yes/no style entries, flag anything else
      s <- tolower(trimws(as.character(raw)))
      v <- suppressWarnings(as.numeric(s))
      v[s %in% c("yes", "y", "true", "t", "detected", "present")] <- 1
      v[s %in% c("no", "n", "false", "f", "not detected", "absent")] <- 0
      bad <- is.na(v) & !is.na(raw) & !(s %in% c("", "na", "nan", "."))
      if (any(bad & done)) unreadable[[col]] <<- unique(as.character(raw[bad & done]))
    }
    v[!done] <- NA
    m <- matrix(NA_real_, M, J); m[cbind(ii, jj)] <- v; m
  }
  yA <- mk(cols[1]); yB <- mk(cols[2])
  oc <- list()
  for (k in unlist(d$obs_covs)) {
    if (!k %in% names(det)) stop("survey covariate '", k, "' not in the detections file")
    m <- matrix(if (is.numeric(det[[k]])) NA_real_ else NA_character_, M, J)
    m[cbind(ii, jj)] <- det[[k]]
    m[is.na(yA) | is.na(yB)] <- NA
    oc[[k]] <- m
  }
  site_covs <- if (is.null(sc)) data.frame(row.names = seq_len(M)) else sc[, setdiff(names(sc), sc_), drop = FALSE]
  forms <- if (!is.null(config$formulas)) config$formulas else d$formulas
  list(species = c(dominant = labels[1], subordinate = labels[2]), yA = yA, yB = yB,
       site_covs = site_covs, obs_covs = oc, formulas = forms, site_ids = site_ids,
       survey_labels = as.character(surveys), unreadable = unreadable)
}

wd_check_data <- function(data) {
  pr <- character()
  yA <- data$yA; yB <- data$yB
  if (!identical(dim(yA), dim(yB))) return("the two species' detection tables have different sizes")
  ids <- if (is.null(data$site_ids)) seq_len(nrow(yA)) else data$site_ids
  for (col in names(data$unreadable))
    pr <- c(pr, sprintf("detection column '%s' contains entries that are not 0/1 (or yes/no): %s. Correct them; leave the cell blank only if that survey was not done.",
                        col, paste(sprintf("'%s'", utils::head(data$unreadable[[col]], 5)), collapse = ", ")))
  for (k in 1:2) {
    y <- list(yA, yB)[[k]]; nm <- data$species[[k]]
    v <- y[!is.na(y)]
    if (any(!v %in% c(0, 1))) pr <- c(pr, sprintf("detections of '%s' must be 0 or 1; found %s. Convert counts to detected (1) / not detected (0).",
                                                   nm, paste(utils::head(unique(v[!v %in% c(0, 1)]), 5), collapse = ", ")))
    if (sum(v) == 0) pr <- c(pr, sprintf("'%s' was never detected; its occupancy cannot be estimated.", nm))
  }
  if (!identical(is.na(yA), is.na(yB))) {
    bad <- which(is.na(yA) != is.na(yB), arr.ind = TRUE)[1, ]
    pr <- c(pr, sprintf("site %s survey %d records one species but not the other; every survey must record both species (0 or 1) or neither.",
                        ids[bad[1]], bad[2]))
  }
  ns <- rowSums(!is.na(yA))
  if (any(ns == 0)) pr <- c(pr, sprintf("%d site(s) have no surveys at all (e.g. %s); remove them from the site file.",
                                         sum(ns == 0), paste(utils::head(ids[ns == 0], 5), collapse = ", ")))
  if (stats::median(ns) < 2) pr <- c(pr, "most sites have fewer than 2 surveys; detection cannot be separated from occupancy.")
  rs <- tryCatch(wd_resolve(data$formulas), error = function(e) { pr <<- c(pr, conditionMessage(e)); NULL })
  if (is.null(rs)) return(pr)
  sc <- data$site_covs
  for (p in rs$free) {
    v <- all.vars(stats::as.formula(rs$forms[[p]]))
    if (p %in% wd_occ_pars) {
      miss <- setdiff(v, names(sc))
      if (length(miss)) pr <- c(pr, sprintf("the %s formula uses %s, which is not a column of the site file (occupancy covariates must be site-level).",
                                            p, paste(miss, collapse = ", ")))
      for (k in intersect(v, names(sc))) if (anyNA(sc[[k]]))
        pr <- c(pr, sprintf("site covariate '%s' is missing at %d site(s) (e.g. %s); fill or drop those sites.",
                            k, sum(is.na(sc[[k]])), paste(utils::head(ids[is.na(sc[[k]])], 5), collapse = ", ")))
    } else {
      miss <- setdiff(v, c(names(sc), names(data$obs_covs)))
      if (length(miss)) pr <- c(pr, sprintf("the %s formula uses unknown covariate(s): %s", p, paste(miss, collapse = ", ")))
      for (k in intersect(v, names(data$obs_covs))) {
        bad <- is.na(data$obs_covs[[k]]) & !is.na(yA)
        if (any(bad)) pr <- c(pr, sprintf("survey covariate '%s' is missing for %d survey(s) that were done (e.g. site %s); fill them or mark the survey as not done.",
                                          k, sum(bad), ids[which(bad, arr.ind = TRUE)[1, 1]]))
      }
      for (k in intersect(v, names(sc))) if (anyNA(sc[[k]]))
        pr <- c(pr, sprintf("site covariate '%s' (used for detection) has missing values", k))
    }
  }
  unique(pr)
}

# Non-blocking warnings about thin data, added to the fit notes.
wd_data_warnings <- function(des) {
  eA <- rowSums(des$yA, na.rm = TRUE) > 0; eB <- rowSums(des$yB, na.rm = TRUE) > 0
  out <- character()
  dn <- des$species[["dominant"]]; sn <- des$species[["subordinate"]]
  if (sum(eA & eB) < 10)
    out <- c(out, sprintf("Only %d site(s) had both %s and %s detected; occupancy of %s where %s is present (psiBA) and the detection parameters rA, rBA, rBa rest on very few sites, so their intervals are wide.",
                          sum(eA & eB), dn, sn, sn, dn))
  if (sum(eA) < 10) out <- c(out, sprintf("%s was detected at only %d site(s).", dn, sum(eA)))
  if (sum(eB) < 10) out <- c(out, sprintf("%s was detected at only %d site(s).", sn, sum(eB)))
  out
}
