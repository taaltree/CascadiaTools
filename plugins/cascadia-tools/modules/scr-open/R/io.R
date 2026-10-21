# scr-open: reading the researcher's files, data checks, writing simulated data.
#
# config$data keys:
#   traps: CSV, one row per trap (trap_id_col, x_col, y_col)
#   usage_cols: optional list, one column of the trap file per period, each
#               either the number of occasions the trap operated (0..K) or a
#               0/1 string with one character per occasion ("11011")
#   captures: CSV, one row per detection (individual_col, period_col,
#             occasion_col, trap_col; extra columns ignored)
#   periods: optional list of period values in time order (e.g. [2019, 2020,
#            ...]); default = sorted unique values of the period column.
#            Give it whenever a period had no detections.
#   n_occasions: K, secondary occasions per period (required)
#   units: m | km (omit if unknown)
#   mask: optional habitat CSV (mask_x_col, mask_y_col, habitat_col,
#         density_covariates); else built from traps + buffer
#   buffer, mask_spacing: optional (map units)
#   labels: display labels for covariates

so_col <- function(df, name, what, file) {
  if (!name %in% names(df)) stop(sprintf("%s: column '%s' (%s) not found; columns are: %s", file, name, what,
                                         paste(names(df), collapse = ", ")), call. = FALSE)
  df[[name]]
}

so_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc)) stop("config needs a data: block (traps, captures, n_occasions, ...)")
  g <- function(k, d) if (is.null(dc[[k]])) d else dc[[k]]
  rd <- function(f) utils::read.csv(file.path(config$.dir, f), check.names = FALSE, stringsAsFactors = FALSE,
                                    strip.white = TRUE, colClasses = "character")
  tdf <- rd(dc$traps)
  traps <- data.frame(id = as.character(so_col(tdf, g("trap_id_col", "id"), "trap ID", dc$traps)),
                      x = as.numeric(so_col(tdf, g("x_col", "x"), "x coordinate", dc$traps)),
                      y = as.numeric(so_col(tdf, g("y_col", "y"), "y coordinate", dc$traps)))
  K <- dc$n_occasions
  if (is.null(K)) stop("data: n_occasions (secondary occasions per primary period) is required")
  cdf <- rd(dc$captures)
  pcol <- so_col(cdf, g("period_col", "period"), "primary period (e.g. year)", dc$captures)
  periods <- if (!is.null(dc$periods)) as.character(unlist(dc$periods)) else {
    u <- unique(pcol); nu <- suppressWarnings(as.numeric(u))
    if (all(!is.na(nu))) u[order(nu)] else sort(u)
  }
  TT <- length(periods)
  pidx <- match(as.character(pcol), periods)
  cap <- data.frame(individual = as.character(so_col(cdf, g("individual_col", "individual"), "animal ID", dc$captures)),
                    period = pidx,
                    occasion = suppressWarnings(as.integer(so_col(cdf, g("occasion_col", "occasion"), "occasion within period", dc$captures))),
                    trap = as.character(so_col(cdf, g("trap_col", "trap"), "trap ID", dc$captures)))
  attr(cap, "raw_period") <- as.character(pcol)
  uc <- unlist(dc$usage_cols)
  if (!is.null(uc)) {
    if (length(uc) != TT) stop(sprintf("data: usage_cols has %d columns but there are %d periods (%s)", length(uc), TT, paste(periods, collapse = ", ")))
    occ_on <- array(NA_integer_, c(nrow(traps), K, TT))
    usage <- sapply(seq_along(uc), function(t) {
      u <- trimws(as.character(so_col(tdf, uc[[t]], "usage (occasions operated)", dc$traps)))
      isstr <- grepl("^[01]+$", u) & nchar(u) == K & K > 1
      for (j in which(isstr)) occ_on[j, , t] <<- as.integer(strsplit(u[j], "")[[1]])
      ifelse(isstr, nchar(gsub("0", "", u)), suppressWarnings(as.integer(u)))
    })
    usage <- matrix(as.integer(usage), nrow(traps), TT)
  } else { usage <- matrix(as.integer(K), nrow(traps), TT); occ_on <- NULL }
  rownames(usage) <- traps$id
  mask <- NULL
  dcov <- as.character(unlist(dc$density_covariates))
  if (!is.null(dc$mask)) {
    mdf <- rd(dc$mask)
    mask <- data.frame(x = as.numeric(so_col(mdf, g("mask_x_col", "x"), "x", dc$mask)),
                       y = as.numeric(so_col(mdf, g("mask_y_col", "y"), "y", dc$mask)))
    for (v in dcov) mask[[v]] <- as.numeric(so_col(mdf, v, "density covariate", dc$mask))
    if (!is.null(dc$habitat_col)) {
      h <- so_col(mdf, dc$habitat_col, "habitat flag", dc$mask)
      mask <- mask[h %in% c("1", "TRUE", "true"), , drop = FALSE]
    }
  }
  list(traps = traps, usage = usage, occ_on = occ_on, captures = cap, n_occasions = as.integer(K), n_periods = TT,
       period_labels = periods, mask = mask, mask_spacing = dc$mask_spacing, buffer = dc$buffer,
       units = if (is.null(dc$units)) NA_character_ else dc$units,
       density_covariates = dcov, covariate_labels = dc$labels)
}

so_check_data <- function(data) {
  pr <- character()
  tr <- data$traps; K <- data$n_occasions; TT <- data$n_periods
  if (anyNA(tr$x) || anyNA(tr$y)) pr <- c(pr, sprintf("trap file: %d traps have missing or non-numeric coordinates (e.g. %s); fix or remove them",
                                                      sum(is.na(tr$x) | is.na(tr$y)), paste(utils::head(tr$id[is.na(tr$x) | is.na(tr$y)], 3), collapse = ", ")))
  if (any(duplicated(tr$id))) pr <- c(pr, paste("trap file: duplicated trap IDs:", paste(unique(tr$id[duplicated(tr$id)]), collapse = ", ")))
  if (any(duplicated(tr[, c("x", "y")]))) pr <- c(pr, "trap file: two traps share the same coordinates; merge them or check the coordinates")
  if (TT < 2) pr <- c(pr, "only one primary period: use module 'scr' (closed population) instead")
  u <- data$usage
  badu <- which(is.na(u) | u < 0 | u > K, arr.ind = TRUE)
  if (length(badu)) {
    if (!is.matrix(badu)) badu <- matrix(badu, ncol = 2)
    ex <- utils::head(sprintf("trap %s, period %s", tr$id[badu[, 1]], (if (is.null(data$period_labels)) seq_len(TT) else data$period_labels)[badu[, 2]]), 3)
    pr <- c(pr, sprintf("trap usage must be the number of occasions operated (0 to %d) or a 0/1 string of length %d for every trap and period; %d entries are not (e.g. %s)",
                        K, K, nrow(badu), paste(ex, collapse = "; ")))
  }
  cap <- data$captures
  if (!nrow(cap)) pr <- c(pr, "capture file has no detections")
  if (anyNA(cap$period)) {
    rp <- attr(cap, "raw_period")
    if (is.null(rp)) rp <- rep(NA_character_, nrow(cap))
    miss <- is.na(rp) | trimws(rp) == ""
    if (any(miss)) pr <- c(pr, sprintf("capture file: %d rows have no period (year) value (e.g. row %d); fill them in or remove those rows", sum(miss), which(miss)[1]))
    unl <- is.na(cap$period) & !miss
    if (any(unl)) pr <- c(pr, sprintf("capture file: %d rows have period values that are not in the period list (%s), e.g. %s; check for typos or add them to data: periods",
                                      sum(unl), paste(data$period_labels, collapse = ", "), paste(utils::head(unique(rp[unl]), 3), collapse = ", ")))
  }
  unk <- setdiff(unique(cap$trap), tr$id)
  if (length(unk)) pr <- c(pr, sprintf("capture file: %d trap IDs are not in the trap file (e.g. %s); check spelling/case",
                                        length(unk), paste(utils::head(unk, 5), collapse = ", ")))
  if (anyNA(cap$occasion)) pr <- c(pr, "capture file: some occasions are missing or not whole numbers")
  badk <- cap$occasion[!is.na(cap$occasion) & (cap$occasion < 1 | cap$occasion > K)]
  if (length(badk)) pr <- c(pr, sprintf("capture file: occasions must be 1..%d (n_occasions, numbered within each period) but found %s", K, paste(utils::head(unique(badk), 5), collapse = ", ")))
  if (anyNA(cap$individual) || any(cap$individual == "")) pr <- c(pr, "capture file: some rows have no animal ID")
  if (!length(unk) && !anyNA(cap$occasion) && !length(badk) && !anyNA(cap$period) && !anyNA(u)) {
    jj <- match(cap$trap, tr$id)
    oo <- data$occ_on
    on_str <- if (is.null(oo)) rep(NA_integer_, nrow(cap)) else oo[cbind(jj, cap$occasion, cap$period)]
    # 0/1 strings: the detection must fall on an operating occasion; counts: occasion <= count
    off <- ifelse(is.na(on_str), cap$occasion > u[cbind(jj, cap$period)], on_str == 0L)
    if (any(off)) pr <- c(pr, sprintf("%d detections are at traps recorded as not operating on that occasion (e.g. trap %s, period %s, occasion %d); fix the usage columns (a count k means occasions 1..k; a 0/1 string gives each occasion)",
                                      sum(off), cap$trap[off][1], data$period_labels[cap$period[off][1]], cap$occasion[off][1]))
  }
  if (!anyNA(cap$period) && nrow(cap)) {
    np <- tabulate(cap$period, TT)
    if (any(np == 0)) pr <- c(pr, sprintf("no detections at all in period(s) %s: check the data (if real, the model can still run but those periods are estimated from the priors and neighbours only; remove this check by fixing the period list)",
                                         paste(data$period_labels[np == 0], collapse = ", ")))
  }
  if (!is.null(data$mask)) {
    m <- data$mask
    if (nrow(m) < 50) pr <- c(pr, sprintf("habitat mask has only %d pixels; use a finer grid (aim for 500-2500)", nrow(m)))
    for (v in data$density_covariates) if (!v %in% names(m)) pr <- c(pr, sprintf("habitat mask has no column '%s'", v))
    sp <- tryCatch(so_infer_spacing(m), error = function(e) NA)
    if (is.finite(sp) && !anyNA(tr$x) && !anyNA(tr$y)) {
      off <- abs((m$x - min(m$x)) / sp - round((m$x - min(m$x)) / sp)) > 0.01 | abs((m$y - min(m$y)) / sp - round((m$y - min(m$y)) / sp)) > 0.01
      if (any(off)) pr <- c(pr, sprintf("habitat mask: %d points are not on a regular grid of spacing %g", sum(off), sp))
      dt <- so_dist_to_traps(tr, data.frame(id = seq_len(nrow(m)), x = m$x, y = m$y))
      if (any(dt > sp)) pr <- c(pr, sprintf("%d traps lie outside the habitat mask; the mask must cover the traps plus a buffer", sum(dt > sp)))
    }
  } else if (length(data$density_covariates)) pr <- c(pr, "density covariates need a habitat mask CSV (data: mask)")
  if (!is.na(data$units) && !data$units %in% c("m", "km")) pr <- c(pr, "data: units must be m or km (or omitted if unknown)")
  pr
}

so_write_data <- function(data, out) {
  tr <- data$traps
  for (t in seq_len(data$n_periods)) tr[[paste0("occasions_", data$period_labels[t])]] <- data$usage[, t]
  utils::write.csv(tr, file.path(out, "scr-open_traps.csv"), row.names = FALSE)
  cp <- data$captures; cp$period <- data$period_labels[cp$period]
  utils::write.csv(cp, file.path(out, "scr-open_captures.csv"), row.names = FALSE)
  if (!is.null(data$mask)) utils::write.csv(data$mask, file.path(out, "scr-open_habitat.csv"), row.names = FALSE)
  invisible(out)
}
