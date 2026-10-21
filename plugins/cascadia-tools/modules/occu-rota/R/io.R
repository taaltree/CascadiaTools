# occu-rota: reading researchers' files and checking them.
#
# config$data keys (paths relative to the config file):
#   detections     CSV of surveys.
#   format         "wide_species" (default): one row per site x survey, one 0/1
#                  column per species (typical camera-trap / point-count sheet);
#                  "long": one row per site x survey x species with a 0/1 column.
#   site_col, survey_col           column names (site id; survey/occasion number)
#   species_cols   (wide_species) the species detection columns, 2 or 3
#   species_col, detection_col     (long) species name column and 0/1 column
#   species        (long, optional) which species to keep, in this order
#   surveyed_col   optional: survey effort column; rows where it is 0 or NA
#                  count as "not surveyed" (e.g. camera not working that week)
#   obs_covs       survey-level covariate columns of the detections file
#   site_covs      CSV with one row per site (site_col + covariates)
#   max_order      2 (pairwise, default) or 3 (adds the three-way term; 3 species only)
# config$formulas: occupancy, interaction, detection (see README).

rota_read_data <- function(config) {
  d <- config$data
  if (is.null(d)) stop("config has no 'data' section")
  rd <- function(f) utils::read.csv(file.path(config$.dir, f), stringsAsFactors = FALSE, check.names = FALSE)
  det <- rd(d$detections)
  fmt <- if (is.null(d$format)) "wide_species" else d$format
  sc_ <- if (is.null(d$site_col)) "site" else d$site_col
  sv_ <- if (is.null(d$survey_col)) "survey" else d$survey_col
  need <- c(sc_, sv_)
  miss <- setdiff(need, names(det))
  if (length(miss)) stop("detections file lacks column(s): ", paste(miss, collapse = ", "),
                         " (set data$site_col / data$survey_col in the config). Columns in the file: ",
                         paste(names(det), collapse = ", "))
  if (fmt == "long") {
    spc <- if (is.null(d$species_col)) "species" else d$species_col
    dc <- if (is.null(d$detection_col)) "detected" else d$detection_col
    if (!all(c(spc, dc) %in% names(det))) stop("long format needs columns '", spc, "' and '", dc, "'")
    species_raw <- if (!is.null(d$species)) as.character(unlist(d$species)) else sort(unique(det[[spc]]))
    det <- det[det[[spc]] %in% species_raw, ]
    key <- paste(det[[sc_]], det[[sv_]], sep = "\r")
    first <- !duplicated(key)
    base <- det[first, setdiff(names(det), c(spc, dc)), drop = FALSE]
    for (s in species_raw) {
      sub <- det[det[[spc]] == s, ]
      base[[s]] <- sub[[dc]][match(key[first], paste(sub[[sc_]], sub[[sv_]], sep = "\r"))]
    }
    det <- base; species_cols <- species_raw
  } else {
    species_cols <- as.character(unlist(d$species_cols))
    if (length(species_cols) < 2) stop("data$species_cols must list the 2 or 3 species detection columns")
    miss <- setdiff(species_cols, names(det))
    if (length(miss)) stop("detections file lacks species column(s): ", paste(miss, collapse = ", "),
                           ". Columns in the file: ", paste(names(det), collapse = ", "),
                           " (check spelling and capitalisation in data$species_cols)")
  }
  species <- rota_clean(species_cols)
  site_ids <- if (!is.null(d$site_covs)) NULL else sort(unique(as.character(det[[sc_]])))
  sc <- NULL
  if (!is.null(d$site_covs)) {
    sc <- rd(d$site_covs)
    if (!sc_ %in% names(sc)) stop("site covariate file lacks the site column '", sc_, "'")
    if (anyDuplicated(sc[[sc_]])) stop("site covariate file has duplicated site ids: ",
                                       paste(unique(sc[[sc_]][duplicated(sc[[sc_]])]), collapse = ", "))
    site_ids <- as.character(sc[[sc_]])
    extra <- setdiff(unique(as.character(det[[sc_]])), site_ids)
    if (length(extra)) {
      # suggest the closest site id in the site file (typos, case, O vs 0)
      hint <- vapply(utils::head(extra, 10), function(x) {
        dd <- utils::adist(tolower(x), tolower(site_ids))[1, ]
        if (min(dd) <= 2) sprintf("%s (did you mean %s?)", x, site_ids[which.min(dd)]) else x
      }, "")
      stop("site(s) in the detections file but not in the site covariate file: ",
           paste(hint, collapse = ", "), ". Site ids must match exactly (including capitalisation) in both files.")
    }
  }
  surveys <- sort(unique(det[[sv_]]))
  M <- length(site_ids); J <- length(surveys)
  ii <- match(as.character(det[[sc_]]), site_ids); jj <- match(det[[sv_]], surveys)
  dup <- duplicated(cbind(ii, jj))
  if (any(dup)) stop("detections file has more than one row for site ", det[[sc_]][dup][1], " survey ", det[[sv_]][dup][1],
                     " (one row per site x survey is expected)")
  done <- rep(TRUE, nrow(det))
  if (!is.null(d$surveyed_col)) {
    v <- det[[d$surveyed_col]]
    if (is.null(v)) stop("surveyed_col '", d$surveyed_col, "' not in the detections file")
    done <- !is.na(v) & v != 0
    has_det <- Reduce(`|`, lapply(species_cols, function(k) !is.na(det[[k]])))
    nodone <- which(is.na(v) & has_det)
    if (length(nodone)) attr(done, "na_effort") <- sprintf("site %s survey %s", det[[sc_]][nodone], det[[sv_]][nodone])
  }
  y <- array(NA_real_, c(M, J, length(species)))
  for (s in seq_along(species_cols)) {
    v <- suppressWarnings(as.numeric(det[[species_cols[s]]]))
    v[!done] <- NA
    y[cbind(ii, jj, s)] <- v
  }
  # a survey with any species missing is treated as not done for all; a survey
  # recorded for some species but blank for others is almost always a data-entry
  # slip, so it is remembered and reported by check_data()
  anyna <- apply(is.na(y), c(1, 2), any); somena <- anyna & !apply(is.na(y), c(1, 2), all)
  partial <- which(somena, arr.ind = TRUE)
  partial_na <- if (nrow(partial)) sprintf("site %s survey %s", site_ids[partial[, 1]], surveys[partial[, 2]]) else character()
  allna <- anyna
  for (s in seq_along(species)) { m <- y[, , s]; m[allna] <- NA; y[, , s] <- m }
  oc <- list()
  for (k in unlist(d$obs_covs)) {
    if (!k %in% names(det)) stop("survey covariate '", k, "' not in the detections file")
    m <- matrix(if (is.numeric(det[[k]])) NA_real_ else NA_character_, M, J)
    m[cbind(ii, jj)] <- det[[k]]
    m[allna] <- NA
    oc[[k]] <- m
  }
  site_covs <- if (is.null(sc)) data.frame(row.names = seq_len(M)) else sc[, setdiff(names(sc), sc_), drop = FALSE]
  forms <- if (!is.null(config$formulas)) config$formulas else d$formulas
  list(species = species, y = y, site_covs = site_covs, obs_covs = oc, formulas = forms,
       max_order = if (is.null(d$max_order)) 2L else as.integer(d$max_order),
       site_ids = site_ids, survey_labels = as.character(surveys), species_labels = species_cols,
       partial_na = partial_na, na_effort = attr(done, "na_effort"), surveyed_col = d$surveyed_col)
}

rota_check_data <- function(data) {
  pr <- character()
  S <- length(data$species)
  if (!S %in% 2:3) return(sprintf("occu-rota needs 2 or 3 species; got %d (%s). For more species use a joint species distribution model.",
                                  S, paste(data$species, collapse = ", ")))
  y <- data$y
  if (length(dim(y)) != 3 || dim(y)[3] != S) return("y must be a sites x surveys x species array")
  v <- y[!is.na(y)]
  if (any(!v %in% c(0, 1))) pr <- c(pr, sprintf("detection values must be 0 or 1; found %s. Convert counts to detected (1) / not detected (0).",
                                                 paste(utils::head(unique(v[!v %in% c(0, 1)]), 5), collapse = ", ")))
  na1 <- is.na(y[, , 1])
  for (s in seq_len(S)[-1]) if (!identical(na1, is.na(y[, , s])))
    pr <- c(pr, sprintf("species '%s' has missing surveys where '%s' does not; every survey must record all species (or none).",
                        data$species[s], data$species[1]))
  if (length(data$partial_na))
    pr <- c(pr, sprintf("%d survey(s) have a 0/1 value for some species but are blank for others (e.g. %s). Every survey must record all species: enter 0 (not detected) or 1, or, if the survey was not done, leave all species blank or set the effort column to 0.",
                        length(data$partial_na), paste(utils::head(data$partial_na, 3), collapse = "; ")))
  if (length(data$na_effort))
    pr <- c(pr, sprintf("%d survey(s) have detection values but a blank '%s' (e.g. %s); they were treated as NOT surveyed. If the camera was working, fill in '%s'; if not, this is fine (blank the detections to silence this message).",
                        length(data$na_effort), data$surveyed_col, paste(utils::head(data$na_effort, 3), collapse = "; "), data$surveyed_col))
  ns <- rowSums(!na1)
  ids <- if (is.null(data$site_ids)) seq_len(dim(y)[1]) else data$site_ids
  if (any(ns == 0)) pr <- c(pr, sprintf("%d site(s) have no surveys at all (e.g. %s); remove them from the site file.",
                                         sum(ns == 0), paste(utils::head(ids[ns == 0], 5), collapse = ", ")))
  if (stats::median(ns) < 2) pr <- c(pr, "most sites have fewer than 2 surveys; detection cannot be separated from occupancy.")
  for (s in seq_len(S)) if (sum(y[, , s], na.rm = TRUE) == 0)
    pr <- c(pr, sprintf("species '%s' was never detected; its occupancy and interactions cannot be estimated.", data$species[s]))
  if (anyDuplicated(data$species)) pr <- c(pr, "species names must be unique")
  sc <- data$site_covs
  fm <- data$formulas
  vars_state <- unique(unlist(lapply(c(unlist(fm$occupancy), unlist(fm$interaction)), function(f) {
    f <- rota_form_str(f); if (rota_is_zero(f)) character() else all.vars(stats::as.formula(f)) })))
  vars_det <- unique(unlist(lapply(unlist(fm$detection), function(f) all.vars(stats::as.formula(rota_form_str(f))))))
  miss <- setdiff(vars_state, names(sc))
  if (length(miss)) pr <- c(pr, sprintf("occupancy/interaction formulas use covariate(s) not in the site file: %s (occupancy covariates must be site-level).",
                                         paste(miss, collapse = ", ")))
  for (k in intersect(vars_state, names(sc))) if (anyNA(sc[[k]]))
    pr <- c(pr, sprintf("site covariate '%s' is missing at %d site(s) (e.g. %s); fill or drop those sites.",
                        k, sum(is.na(sc[[k]])), paste(utils::head(ids[is.na(sc[[k]])], 5), collapse = ", ")))
  # log() of a covariate needs strictly positive values
  all_f <- vapply(c(unlist(fm$occupancy), unlist(fm$interaction), unlist(fm$detection)), rota_form_str, "")
  for (k in unique(c(vars_state, vars_det))) {
    if (!any(grepl(paste0("log\\(\\s*", k, "\\s*\\)"), all_f))) next
    v <- if (k %in% names(sc)) sc[[k]] else if (k %in% names(data$obs_covs)) as.vector(data$obs_covs[[k]][!na1]) else NULL
    if (is.numeric(v) && any(v <= 0, na.rm = TRUE))
      pr <- c(pr, sprintf("covariate '%s' is used as log(%s) but has %d value(s) <= 0; use log(%s + 1) (or another offset) in the formula, or check the data.",
                          k, k, sum(v <= 0, na.rm = TRUE), k))
  }
  # factor levels that differ only by capitalisation or spaces ("yes" vs "Yes ")
  for (k in intersect(unique(c(vars_state, vars_det)), names(sc))) {
    v <- sc[[k]]
    if (!(is.character(v) || is.factor(v))) next
    lv <- unique(as.character(stats::na.omit(v))); key <- tolower(trimws(lv))
    if (anyDuplicated(key)) pr <- c(pr, sprintf("site covariate '%s' has levels that differ only by capitalisation or spaces: %s. Make them identical.",
                                                k, paste(sprintf("'%s'", lv[key %in% key[duplicated(key)]]), collapse = ", ")))
  }
  missd <- setdiff(vars_det, c(names(sc), names(data$obs_covs)))
  if (length(missd)) pr <- c(pr, sprintf("detection formulas use unknown covariate(s): %s", paste(missd, collapse = ", ")))
  for (k in intersect(vars_det, names(data$obs_covs))) {
    m <- data$obs_covs[[k]]
    bad <- is.na(m) & !na1
    if (any(bad)) pr <- c(pr, sprintf("survey covariate '%s' is missing for %d survey(s) that were done (e.g. site %s); fill them or mark the survey as not done.",
                                      k, sum(bad), ids[which(bad, arr.ind = TRUE)[1, 1]]))
  }
  for (k in intersect(vars_det, names(sc))) if (anyNA(sc[[k]]))
    pr <- c(pr, sprintf("site covariate '%s' (used for detection) has missing values", k))
  pr
}

# Non-blocking warnings about sparse data, added to the fit notes.
rota_data_warnings <- function(des) {
  y0 <- des$y; y0[is.na(y0)] <- 0
  det <- apply(y0, c(1, 3), max)
  out <- character()
  for (a in seq_len(des$S - 1)) for (b in (a + 1):des$S) {
    both <- sum(det[, a] * det[, b])
    if (both < 5)
      out <- c(out, sprintf("Only %d site(s) had both %s and %s detected; the %s-%s interaction is weakly identified and its interval will be wide or unstable.",
                            both, des$species[a], des$species[b], des$species[a], des$species[b]))
  }
  nd <- colSums(det)
  for (s in seq_len(des$S)) if (nd[s] < 10)
    out <- c(out, sprintf("%s was detected at only %d site(s); its occupancy and interactions are poorly informed.", des$species[s], nd[s]))
  out
}
