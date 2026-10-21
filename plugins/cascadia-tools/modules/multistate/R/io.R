# multistate: reading researchers' files and checking them.
#
# Capture file, one of three formats (data$format):
#   long (default): one row per capture: id, occasion, state
#                   state = a state name (e.g. "infected"), or the unknown code
#                   (data$unknown_code, default "U") when the animal was caught but
#                   its state could not be determined.
#   wide:           one row per animal: id, then one column per occasion holding
#                   0 (not seen), a state name or number 1..S, or the unknown code.
#   ch:             one row per animal: id and a capture-history string such as
#                   "0A0BU0" (one character per occasion; 0 = not seen).
# Individual covariates (optional): id + one column per covariate.
# Occasion covariates (optional): occasion + one column per covariate. A value
#   at occasion t is used for detection at occasion t and for survival and
#   transitions over the interval t -> t+1.

ms_read_csv <- function(dir, f) {
  if (is.null(f)) return(NULL)
  p <- if (file.exists(f)) f else file.path(dir, f)
  if (!file.exists(p)) stop("file not found: ", f, " (looked in ", dir, ")")
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE,
                  na.strings = c("NA", ""), colClasses = "character")
}

ms_occ_order <- function(x, given) {
  if (!is.null(given)) return(as.character(unlist(given)))
  u <- unique(as.character(x))
  num <- suppressWarnings(as.numeric(u))
  if (!anyNA(num)) return(u[order(num)])
  if (mean(!is.na(num)) >= 0.5)
    stop(sprintf(paste("occasion value(s) %s are not numbers, unlike the other occasions (%s). Check for typos",
                       "(e.g. the letter O instead of zero) in the occasion column, or list the occasions in data$occasion_order, earliest first."),
                 paste0("'", u[is.na(num)], "'", collapse = ", "), paste(utils::head(u[!is.na(num)], 5), collapse = ", ")))
  d <- regmatches(u, gregexpr("[0-9]+", u)); stem <- gsub("[0-9]+", "#", u)
  if (all(lengths(d) == 1) && length(unique(stem)) == 1) return(u[order(as.numeric(unlist(d)))])
  stop(sprintf("cannot tell the time order of the occasions (%s); list them in data$occasion_order, earliest first",
               paste(utils::head(u, 8), collapse = ", ")))
}

ms_code_state <- function(v, states, unk, where, issues, flag_empty = FALSE) {
  v <- trimws(as.character(v))
  out <- rep(NA_integer_, length(v))
  out[v %in% c("0", "")] <- 0L
  out[v == unk] <- length(states) + 1L
  m <- match(v, states); out[!is.na(m)] <- m[!is.na(m)]
  num <- suppressWarnings(as.integer(v))
  okn <- is.na(out) & !is.na(num) & num >= 1 & num <= length(states) + 1
  out[okn] <- num[okn]
  if (flag_empty && any(is.na(v) | v == "")) issues <- c(issues, sprintf(
    "%d capture row(s) have an empty state, e.g. %s. If the animal was caught but its state is unknown, enter the unknown code '%s'; if it was not caught, delete the row.",
    sum(is.na(v) | v == ""), where[is.na(v) | v == ""][1], unk))
  bad <- is.na(out) & !is.na(v)
  if (any(bad)) issues <- c(issues, sprintf(
    "%d state value(s) are not one of the states (%s), the unknown code '%s' or 0, e.g. '%s' at %s. Check spelling/capitalisation or add the state to data$states.",
    sum(bad), paste(states, collapse = ", "), unk, v[bad][1], where[bad][1]))
  out[is.na(out)] <- 0L
  list(code = out, issues = issues)
}

ms_read_data <- function(config) {
  dc <- config$data; dir <- if (is.null(config$.dir)) "." else config$.dir
  if (is.null(dc)) stop("config needs a 'data:' section")
  states <- as.character(unlist(dc$states))
  if (length(states) < 2) stop("data$states must list at least two states in order, e.g. states: [uninfected, infected]")
  unk <- if (is.null(dc$unknown_code)) "U" else as.character(dc$unknown_code)
  idc <- if (is.null(dc$id_col)) "id" else dc$id_col
  fmt <- if (is.null(dc$format)) "long" else dc$format
  cap <- ms_read_csv(dir, dc$captures)
  if (is.null(cap)) stop("config data$captures (the capture file) is required")
  if (!idc %in% names(cap)) stop(sprintf("capture file has no '%s' column (columns: %s); set data$id_col", idc, paste(names(cap), collapse = ", ")))
  issues <- character()
  if (fmt == "long") {
    oc <- if (is.null(dc$occasion_col)) "occasion" else dc$occasion_col
    sc <- if (is.null(dc$state_col)) "state" else dc$state_col
    miss <- setdiff(c(oc, sc), names(cap))
    if (length(miss)) stop(sprintf("capture file lacks column(s) %s (columns: %s); set data$occasion_col / data$state_col",
                                   paste(miss, collapse = ", "), paste(names(cap), collapse = ", ")))
    occ_tab <- ms_read_csv(dir, dc$occasions)
    occs <- ms_occ_order(c(cap[[oc]], if (!is.null(occ_tab)) occ_tab[[if (is.null(dc$occasion_col)) "occasion" else dc$occasion_col]]), dc$occasion_order)
    miss_o <- setdiff(unique(cap[[oc]]), occs)
    if (length(miss_o)) stop("occasion(s) not in data$occasion_order: ", paste(miss_o, collapse = ", "))
    ids <- unique(cap[[idc]])
    cs <- ms_code_state(cap[[sc]], states, unk, sprintf("row %d", seq_len(nrow(cap)) + 1), issues, flag_empty = TRUE)
    issues <- cs$issues
    y <- matrix(0L, length(ids), length(occs))
    ii <- match(cap[[idc]], ids); tt <- match(cap[[oc]], occs)
    for (r in seq_len(nrow(cap))) {
      v <- cs$code[r]; old <- y[ii[r], tt[r]]
      if (old == 0L || old == length(states) + 1L) y[ii[r], tt[r]] <- v
      else if (v <= length(states) && v != old)
        issues <- c(issues, sprintf("animal %s was recorded in two different states (%s and %s) at occasion %s; one state per occasion is needed (use the unknown code '%s' if it cannot be resolved).",
                                    ids[ii[r]], states[old], states[v], occs[tt[r]], unk))
    }
  } else if (fmt == "wide") {
    occs <- setdiff(names(cap), c(idc, unlist(dc$ignore_cols)))
    if (!is.null(dc$occasion_order)) occs <- as.character(unlist(dc$occasion_order))
    ids <- cap[[idc]]
    y <- matrix(0L, nrow(cap), length(occs))
    for (j in seq_along(occs)) {
      cs <- ms_code_state(cap[[occs[j]]], states, unk, sprintf("animal %s, column %s", ids, occs[j]), issues)
      issues <- cs$issues; y[, j] <- cs$code
    }
  } else if (fmt == "ch") {
    chc <- if (is.null(dc$ch_col)) "ch" else dc$ch_col
    ch <- cap[[chc]]; ids <- cap[[idc]]
    n <- nchar(ch)
    if (length(unique(n)) > 1) issues <- c(issues, sprintf("capture histories have different lengths (%s); every history needs one character per occasion.",
                                                           paste(sort(unique(n)), collapse = ", ")))
    Tn <- max(n)
    chars <- t(vapply(ch, function(s) { v <- strsplit(s, "")[[1]]; c(v, rep("0", Tn - length(v))) }, character(Tn)))
    y <- matrix(0L, length(ch), Tn)
    for (j in seq_len(Tn)) {
      cs <- ms_code_state(chars[, j], states, unk, sprintf("animal %s, occasion %d", ids, j), issues)
      issues <- cs$issues; y[, j] <- cs$code
    }
    occs <- if (!is.null(dc$occasion_order)) as.character(unlist(dc$occasion_order)) else as.character(seq_len(Tn))
  } else stop("data$format must be long, wide or ch")
  if (anyDuplicated(ids)) issues <- c(issues, sprintf("animal id(s) appear on more than one row of the wide/ch file: %s",
                                                      paste(utils::head(ids[duplicated(ids)], 5), collapse = ", ")))
  # covariates
  ind <- ms_read_csv(dir, dc$individuals)
  ind_covs <- data.frame(id = ids)
  wanted <- unique(c(unlist(dc$survival), unlist(dc$transition), unlist(dc$detection)))
  if (!is.null(ind)) {
    if (!idc %in% names(ind)) stop(sprintf("individual covariate file has no '%s' column", idc))
    m <- match(ids, ind[[idc]])
    if (anyNA(m)) issues <- c(issues, sprintf("%d animal(s) in the capture file are missing from the individual covariate file, e.g. %s.",
                                              sum(is.na(m)), ids[is.na(m)][1]))
    for (k in setdiff(names(ind), idc)) if (k %in% wanted) ind_covs[[k]] <- suppressWarnings(as.numeric(ind[[k]][m]))
  }
  occ <- ms_read_csv(dir, dc$occasions)
  occ_covs <- data.frame(occasion = occs)
  if (!is.null(occ)) {
    ocn <- if (is.null(dc$occasion_col)) "occasion" else dc$occasion_col
    if (!ocn %in% names(occ)) stop(sprintf("occasion covariate file has no '%s' column", ocn))
    m <- match(occs, as.character(occ[[ocn]]))
    if (anyNA(m)) issues <- c(issues, sprintf("occasion(s) %s have no row in the occasion covariate file.",
                                              paste(occs[is.na(m)], collapse = ", ")))
    for (k in setdiff(names(occ), ocn)) if (k %in% wanted) occ_covs[[k]] <- suppressWarnings(as.numeric(occ[[k]][m]))
  }
  for (k in wanted) if (!k %in% c(names(ind_covs), names(occ_covs)))
    issues <- c(issues, sprintf("covariate '%s' (listed under data$survival/transition/detection) is not a column of the individual or occasion covariate files.", k))
  model <- if (is.null(dc$model)) "auto" else dc$model
  if (model == "auto") model <- if (any(y == length(states) + 1L)) "multievent" else "AS"
  list(y = y, states = states, ids = ids, occasion_labels = occs, ind_covs = ind_covs, occ_covs = occ_covs,
       model = model, unknown_code = unk, interval_unit = if (is.null(dc$interval_unit)) NULL else as.character(dc$interval_unit),
       formulas = list(survival = as.character(unlist(dc$survival)), transition = as.character(unlist(dc$transition)),
                       detection = as.character(unlist(dc$detection))),
       issues = issues)
}

ms_check_data <- function(data) {
  pr <- data$issues
  y <- data$y; S <- length(data$states); T <- ncol(y)
  if (T < 3) pr <- c(pr, sprintf("only %d capture occasions: multistate models need at least 3 (survival and detection are confounded at the last interval).", T))
  if (any(y < 0 | y > S + 1)) pr <- c(pr, "state codes must be 0 (not seen), 1..S or S+1 (state unknown)")
  never <- rowSums(y > 0) == 0
  if (any(never)) pr <- c(pr, sprintf("%d animal(s) have no captures at all (all zeros); remove them or check ids.", sum(never)))
  f <- apply(y, 1, function(r) which(r > 0)[1])
  if (sum(f < T, na.rm = TRUE) < 10) pr <- c(pr, "fewer than 10 animals were marked before the last occasion: nothing to estimate.")
  recap <- sum(vapply(seq_len(nrow(y)), function(i) if (is.na(f[i])) 0 else sum(y[i, -seq_len(f[i])] > 0), 0))
  if (recap < 20) pr <- c(pr, sprintf("only %d recaptures in total: survival, transitions and detection cannot be separated with so few.", recap))
  seen_state <- tabulate(y[y >= 1 & y <= S], S)
  if (any(seen_state == 0)) pr <- c(pr, sprintf("state(s) %s never recorded: drop them from data$states or check the codes.",
                                                paste(data$states[seen_state == 0], collapse = ", ")))
  if (identical(data$model, "AS") && any(y == S + 1)) pr <- c(pr, "the data contain 'state unknown' records but data$model is AS; use model: multievent (or auto).")
  for (k in setdiff(names(data$ind_covs), "id")) if (anyNA(data$ind_covs[[k]]))
    pr <- c(pr, sprintf("individual covariate '%s' is missing (or not a number) for %d animal(s), e.g. %s; fill it in or drop the covariate.",
                        k, sum(is.na(data$ind_covs[[k]])), data$ids[is.na(data$ind_covs[[k]])][1]))
  for (k in setdiff(names(data$occ_covs), "occasion")) {
    v <- data$occ_covs[[k]]
    need <- c(if (k %in% data$formulas$detection) seq_len(T)[-1],
              if (k %in% c(data$formulas$survival, data$formulas$transition)) seq_len(T - 1))
    need <- sort(unique(need))
    if (anyNA(v[need])) pr <- c(pr, sprintf("occasion covariate '%s' is missing at occasion(s) %s.", k, paste(data$occasion_labels[need][is.na(v[need])], collapse = ", ")))
  }
  unique(pr)
}

# Write simulated data in the long researcher format (core/run.R simulate).
ms_write_data <- function(data, out) {
  S <- length(data$states); unk <- if (is.null(data$unknown_code)) "U" else data$unknown_code
  w <- which(data$y > 0, arr.ind = TRUE)
  st <- ifelse(data$y[w] == S + 1, unk, data$states[pmin(data$y[w], S)])
  cap <- data.frame(id = data$ids[w[, 1]], occasion = data$occasion_labels[w[, 2]], state = st)
  cap <- cap[order(cap$id, as.numeric(cap$occasion)), ]
  utils::write.csv(cap, file.path(out, "multistate_captures.csv"), row.names = FALSE)
  utils::write.csv(data$ind_covs, file.path(out, "multistate_individuals.csv"), row.names = FALSE)
  utils::write.csv(data$occ_covs, file.path(out, "multistate_occasions.csv"), row.names = FALSE)
  invisible(out)
}
