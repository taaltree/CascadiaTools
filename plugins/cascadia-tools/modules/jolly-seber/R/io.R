# jolly-seber: reading researchers' files, data checks, writing simulated data.
#
# Two input formats (config data: block):
#  1. captures: long CSV, one row per capture (animal ID + occasion, optional
#     group column such as sex, any other columns ignored). Several captures
#     of one animal within an occasion count once.
#  2. histories: one row per animal with a capture-history string ("0110100",
#     MARK style), optional freq and group columns.
# Optional occasions: CSV with one row per occasion (also occasions with no
# captures) giving its key, an interval to the next occasion or a numeric
# time / date, and occasion covariates for detection (e.g. effort).

js_col <- function(df, name, what, file) {
  if (!name %in% names(df)) {
    near <- names(df)[tolower(gsub("[^A-Za-z0-9]", "", names(df))) == tolower(gsub("[^A-Za-z0-9]", "", name))]
    stop(sprintf("%s: column '%s' (%s) not found; columns are: %s%s", file, name, what,
                 paste(names(df), collapse = ", "),
                 if (length(near)) sprintf(". Did you mean '%s'? Rename the column or change the config entry to that name", near[1]) else ""),
         call. = FALSE)
  }
  df[[name]]
}

js_unknown_values <- c("", "U", "u", "unk", "unknown", "Unknown", "UNK", "?", "NA", "na")

js_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc)) stop("config needs a data: block (captures or histories, occasions, ...)")
  g <- function(k, d) if (is.null(dc[[k]])) d else dc[[k]]
  rd <- function(f) utils::read.csv(file.path(config$.dir, f), check.names = FALSE, stringsAsFactors = FALSE,
                                    strip.white = TRUE, colClasses = "character")
  notes <- character()
  # occasions table
  occ <- NULL
  if (!is.null(dc$occasions)) {
    odf <- rd(dc$occasions)
    key <- trimws(as.character(js_col(odf, g("occasion_key", g("occasion_col", "occasion")), "occasion key", dc$occasions)))
    occ <- data.frame(key = key, stringsAsFactors = FALSE)
    if (!is.null(dc$time_col)) {
      tv <- js_col(odf, dc$time_col, "occasion time", dc$occasions)
      tn <- suppressWarnings(as.numeric(tv))
      if (anyNA(tn)) {   # dates
        d <- as.Date(tv, tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y"))
        div <- switch(g("time_unit", "year"), year = 365.25, month = 365.25 / 12, week = 7, day = 1, 365.25)
        tn <- as.numeric(d) / div
      }
      occ$time <- tn
      occ <- occ[order(occ$time), , drop = FALSE]
      occ$interval <- c(diff(occ$time), NA)
    } else if (!is.null(dc$interval_col)) {
      occ$interval <- suppressWarnings(as.numeric(js_col(odf, dc$interval_col, "interval to next occasion", dc$occasions)))
    } else occ$interval <- c(rep(1, nrow(occ) - 1), NA)
    for (v in unlist(dc$detection_covariates)) {
      raw <- trimws(as.character(js_col(odf, v, "detection covariate", dc$occasions)))
      num <- suppressWarnings(as.numeric(raw))
      nonnum <- !is.na(raw) & nzchar(raw) & is.na(num)
      if (any(nonnum)) stop(sprintf("%s: detection covariate '%s' must be a number; found %s at occasion(s) %s", dc$occasions, v,
                                    paste0("'", raw[nonnum], "'", collapse = ", "), paste(key[nonnum], collapse = ", ")), call. = FALSE)
      occ[[v]] <- num[match(occ$key, key)]
    }
    lab <- if (!is.null(dc$label_col)) as.character(js_col(odf, dc$label_col, "occasion label", dc$occasions)[match(occ$key, key)]) else occ$key
    occ$label <- lab
  }
  grp_col <- dc$group_col
  if (!is.null(dc$captures)) {
    cdf <- rd(dc$captures)
    id <- trimws(as.character(js_col(cdf, g("id_col", "animal_id"), "animal ID", dc$captures)))
    oc <- trimws(as.character(js_col(cdf, g("occasion_col", "occasion"), "occasion", dc$captures)))
    # rows without an animal ID cannot enter a capture history: drop them, say which
    noid <- is.na(id) | !nzchar(id)
    if (any(noid)) {
      notes <- c(notes, sprintf("%d capture row(s) in %s have no %s (file rows %s) and were left out; add the IDs if they are known.",
                                sum(noid), dc$captures, g("id_col", "animal_id"), paste(utils::head(which(noid) + 1, 10), collapse = ", ")))
      cdf <- cdf[!noid, , drop = FALSE]; id <- id[!noid]; oc <- oc[!noid]
    }
    if (anyNA(oc) || any(!nzchar(oc)))
      stop(sprintf("%s: %d capture row(s) have no %s value (file rows %s); every capture needs its occasion", dc$captures,
                   sum(is.na(oc) | !nzchar(oc)), g("occasion_col", "occasion"),
                   paste(utils::head(which(is.na(oc) | !nzchar(oc)) + 1, 10), collapse = ", ")), call. = FALSE)
    if (is.null(occ)) {
      u <- unique(oc); un <- suppressWarnings(as.numeric(u))
      u <- if (!anyNA(un)) u[order(un)] else sort(u)
      occ <- data.frame(key = u, interval = if (!anyNA(un)) c(diff(sort(un)), NA) else c(rep(1, length(u) - 1), NA),
                        label = u, stringsAsFactors = FALSE)
      notes <- c(notes, "No occasions file: occasions are the distinct values in the capture file; an occasion with no captures at all would be missed.")
    }
    bad <- setdiff(unique(oc), occ$key)
    if (length(bad)) stop(sprintf("%s: occasion value(s) %s are not in the occasions table (%s)", dc$captures,
                                  paste(utils::head(bad, 5), collapse = ", "), paste(occ$key, collapse = ", ")), call. = FALSE)
    ids <- unique(id[nzchar(id)])
    ch <- matrix(0L, length(ids), nrow(occ), dimnames = list(ids, occ$key))
    ok <- nzchar(id)
    ch[cbind(match(id[ok], ids), match(oc[ok], occ$key))] <- 1L
    dup <- sum(duplicated(paste(id, oc)))
    if (dup) notes <- c(notes, sprintf("%d repeat captures of the same animal within an occasion were counted once.", dup))
    grp <- NULL
    if (!is.null(grp_col)) {
      gv <- as.character(js_col(cdf, grp_col, "group (e.g. sex)", dc$captures))
      gv[gv %in% js_unknown_values] <- NA
      # one value per animal: the most frequent non-missing record
      grp <- vapply(ids, function(a) { x <- gv[id == a & !is.na(gv)]; if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1] }, "")
      conf <- vapply(ids, function(a) length(unique(stats::na.omit(gv[id == a]))) > 1, TRUE)
      if (any(conf)) notes <- c(notes, sprintf("%d animals were recorded with conflicting %s; the most frequent value was used (%s).",
                                               sum(conf), grp_col, paste(utils::head(ids[conf], 5), collapse = ", ")))
    }
  } else if (!is.null(dc$histories)) {
    hdf <- rd(dc$histories)
    s <- gsub("[^01]", "", as.character(js_col(hdf, g("history_col", "ch"), "capture history", dc$histories)))
    T <- max(nchar(s))
    if (any(nchar(s) != T)) stop(sprintf("%s: capture histories must all have %d characters", dc$histories, T), call. = FALSE)
    fr <- if (!is.null(dc$freq_col)) as.integer(js_col(hdf, dc$freq_col, "frequency", dc$histories)) else rep(1L, length(s))
    idx <- rep(seq_along(s), fr)
    ch <- t(vapply(s[idx], function(x) as.integer(strsplit(x, "")[[1]]), integer(T)))
    rownames(ch) <- if (!is.null(dc$id_col)) make.unique(as.character(hdf[[dc$id_col]][idx])) else sprintf("H%04d", seq_along(idx))
    grp <- if (!is.null(grp_col)) { gv <- as.character(js_col(hdf, grp_col, "group", dc$histories))[idx]; gv[gv %in% js_unknown_values] <- NA; gv } else NULL
    if (is.null(occ)) occ <- data.frame(key = as.character(seq_len(T)), interval = c(rep(1, T - 1), NA), label = as.character(seq_len(T)))
    if (nrow(occ) != T) stop(sprintf("occasions table has %d rows but capture histories have %d occasions", nrow(occ), T), call. = FALSE)
  } else stop("data: needs either captures: (one row per capture) or histories: (one capture-history string per animal)")
  ch <- ch[rowSums(ch) > 0, , drop = FALSE]
  lev <- NULL
  if (!is.null(grp) && !is.null(dc$group_recode)) {   # e.g. {F: female, M: male}
    rc <- unlist(dc$group_recode)
    hit <- !is.na(grp) & grp %in% names(rc)
    grp[hit] <- rc[grp[hit]]
  }
  if (!is.null(grp)) {
    grp <- grp[rownames(ch)]
    lev <- if (!is.null(dc$group_levels)) as.character(unlist(dc$group_levels)) else sort(unique(stats::na.omit(grp)))
  }
  occ$occasion <- seq_len(nrow(occ))
  colnames(ch) <- occ$label
  list(ch = ch, group = if (!is.null(grp)) unname(grp) else NULL, group_levels = lev,
       group_name = if (!is.null(grp_col)) g("group_name", grp_col) else NULL,
       occasions = occ, detection_covariates = as.character(unlist(dc$detection_covariates)),
       time_unit = g("time_unit", "year"), species = dc$species, read_notes = notes,
       covariate_labels = dc$labels)
}

js_check_data <- function(data) {
  pr <- character()
  ch <- data$ch; T <- ncol(ch)
  if (T < 3) pr <- c(pr, sprintf("only %d sampling occasions: Jolly-Seber needs at least 3 (4+ to estimate survival and recruitment usefully)", T))
  if (!all(ch %in% c(0, 1))) pr <- c(pr, "capture histories must contain only 0 and 1")
  if (nrow(ch) < 10) pr <- c(pr, sprintf("only %d animals were captured; too few for an open-population model", nrow(ch)))
  dt <- data$occasions$interval[seq_len(T - 1)]
  if (any(!is.finite(dt) | dt <= 0)) pr <- c(pr, sprintf("intervals between occasions must be positive numbers; check occasion(s) %s",
                                                         paste(data$occasions$label[which(!is.finite(dt) | dt <= 0)], collapse = ", ")))
  for (v in data$detection_covariates) {
    x <- data$occasions[[v]]
    if (is.null(x) || anyNA(x)) pr <- c(pr, sprintf("detection covariate '%s' is missing for occasion(s) %s (every occasion needs a value)", v,
                                                    paste(data$occasions$label[is.na(x)], collapse = ", ")))
  }
  if (!is.null(data$group)) {
    lev <- data$group_levels
    if (length(lev) != 2) pr <- c(pr, sprintf("the group column must have exactly 2 levels (found: %s); set data: group_levels to the two to compare", paste(lev, collapse = ", ")))
    bad <- setdiff(unique(stats::na.omit(data$group)), lev)
    if (length(bad)) pr <- c(pr, sprintf("group values %s are not among group_levels (%s); recode them or mark them unknown", paste(bad, collapse = ", "), paste(lev, collapse = ", ")))
  }
  pr
}

# Design warnings that don't block the fit.
js_design_notes <- function(prep, est) {
  ch <- prep$ch; T <- prep$T
  nt <- colSums(ch)
  out <- character()
  if (any(nt == 0)) out <- c(out, sprintf("No animals were caught at occasion(s) %s: abundance there rests entirely on the model.", paste(prep$labels[nt == 0], collapse = ", ")))
  rec <- sum(rowSums(ch) > 1)
  if (rec < 20) out <- c(out, sprintf("Only %d animals were recaptured; survival and capture probability are weakly informed.", rec))
  if (prep$Kp_lev == T) out <- c(out, "Capture probability varies freely by occasion: abundance and entry at the FIRST occasion and survival into the LAST occasion are only weakly identified (classic Jolly-Seber confounding); interpret N_1, entry_prob_1 and the last survival with care.")
  if (prep$Kphi_lev > 1 && prep$Kp_lev == T) out <- c(out, "With both survival and capture probability varying by occasion, the last survival and last capture probability are confounded (only their product is estimable).")
  out
}

js_write_data <- function(data, out) {
  ch <- data$ch
  lab <- data$occasions$label
  w <- which(ch > 0, arr.ind = TRUE)
  cap <- data.frame(animal_id = rownames(ch)[w[, 1]], occasion = lab[w[, 2]], stringsAsFactors = FALSE)
  if (!is.null(data$group)) cap[[if (is.null(data$group_name)) "group" else data$group_name]] <- data$group[w[, 1]]
  cap <- cap[order(cap$occasion, cap$animal_id), ]
  utils::write.csv(cap, file.path(out, "captures.csv"), row.names = FALSE, na = "")
  utils::write.csv(data$occasions, file.path(out, "occasions.csv"), row.names = FALSE, na = "")
  invisible(TRUE)
}

# R2ucare goodness-of-fit tests of the CJS part (Pradel et al. 2005;
# Gimenez et al. 2018) on the captured animals (per group if grouped).
js_ucare <- function(prep) {
  if (!requireNamespace("R2ucare", quietly = TRUE)) return(NULL)
  ch <- prep$ch
  if (prep$T < 3 || nrow(ch) < 10) return(NULL)
  gs <- if (prep$G == 2) split(seq_len(nrow(ch)), factor(prep$gcode, levels = 0:2, labels = c("unknown", prep$group_levels))) else list(all = seq_len(nrow(ch)))
  gs <- gs[vapply(gs, length, 1) >= 10]
  res <- lapply(names(gs), function(k) {
    X <- ch[gs[[k]], , drop = FALSE]; fr <- rep(1, nrow(X))
    t3 <- tryCatch(suppressWarnings(R2ucare::test3sr(X, fr)), error = function(e) NULL)
    t2 <- tryCatch(suppressWarnings(R2ucare::test2ct(X, fr)), error = function(e) NULL)
    ov <- tryCatch(suppressWarnings(R2ucare::overall_CJS(X, fr)), error = function(e) NULL)
    rows <- list()
    if (!is.null(t3)) rows[[1]] <- data.frame(group = k, test = "3.SR (transients)", stat = t3$test3sr[["stat"]], df = t3$test3sr[["df"]], p = t3$test3sr[["p_val"]])
    if (!is.null(t2)) rows[[2]] <- data.frame(group = k, test = "2.CT (trap dependence)", stat = t2$test2ct[["stat"]], df = t2$test2ct[["df"]], p = t2$test2ct[["p_val"]])
    if (!is.null(ov)) rows[[3]] <- data.frame(group = k, test = "overall CJS", stat = ov$chi2, df = ov$degree_of_freedom, p = ov$p_value)
    do.call(rbind, rows)
  })
  out <- do.call(rbind, res)
  if (is.null(out) || !nrow(out)) return(NULL)
  rownames(out) <- NULL
  out
}

js_ucare_text <- function(g) {
  s <- sprintf("%s [%s]: chi2 = %.1f, df = %s, p = %.3f", g$test, g$group, g$stat, format(g$df), g$p)
  flag <- any(g$p < 0.05, na.rm = TRUE)
  c(paste0("R2ucare goodness-of-fit tests of the survival/recapture part (captured animals): ", paste(s, collapse = "; "), "."),
    if (flag) "At least one test has p < 0.05: transients (3.SR), trap-happiness/shyness (2.CT) or other heterogeneity may bias survival and abundance; consider a transient or trap-dependence model and treat estimates with caution."
    else "No test flagged lack of fit (p >= 0.05 for all), although these tests have low power with few occasions.")
}

`%||%` <- function(a, b) if (is.null(a)) b else a
