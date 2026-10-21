# cjs: reading researchers' capture data and checking it.
#
# Three capture formats (config data$format):
#   wide  : one row per individual; a capture-history string column (ch_col, e.g. "0110")
#           or one 0/1 column per occasion (occasion_cols)
#   long  : one row per capture: id_col, occasion_col (year or occasion number / label)
#   dates : one row per capture: id_col, date_col; dates are binned into occasions by
#           occasion_unit (year | month | week) or explicit occasion_breaks (start dates)
# Optional tables: individuals (one row per individual: group, age at marking,
# individual covariates), occasion_covs (one row per occasion).

cjs_path <- function(config, f) if (is.null(f)) NULL else if (file.exists(f)) f else file.path(config$.dir, f)

cjs_read_data <- function(config) {
  d <- config$data
  if (is.null(d)) stop("config needs a `data:` block")
  fmt <- d$format %||% "long"
  idc <- d$id_col %||% "id"
  cap <- utils::read.csv(cjs_path(config, d$captures), stringsAsFactors = FALSE, check.names = FALSE)
  if (!idc %in% names(cap)) stop(sprintf("captures file has no id column '%s' (columns: %s)", idc, paste(names(cap), collapse = ", ")))
  cap[[idc]] <- trimws(as.character(cap[[idc]]))
  bad_id <- which(is.na(cap[[idc]]) | cap[[idc]] == "")
  if (length(bad_id)) stop(sprintf("%d rows of the captures file have an empty id in column '%s' (rows %s): fill in or delete them",
                                   length(bad_id), idc, paste(utils::head(bad_id + 1L, 5), collapse = ", ")))
  cols_txt <- paste(names(cap), collapse = ", ")
  notes <- character()
  ind_extra <- NULL
  if (fmt == "wide") {
    if (!is.null(d$ch_col)) {
      chs <- gsub("[^01]", "", as.character(cap[[d$ch_col]]))
      if (length(unique(nchar(chs))) != 1) stop("capture histories in column '", d$ch_col, "' have different lengths")
      ch <- do.call(rbind, lapply(strsplit(chs, ""), as.integer))
      T <- ncol(ch)
      occ_lab <- if (!is.null(d$occasions)) as.character(unlist(d$occasions)) else as.character(seq_len(T))
      used <- d$ch_col
    } else {
      oc <- unlist(d$occasion_cols)
      if (is.null(oc)) stop("wide format needs ch_col (history string) or occasion_cols (one 0/1 column per occasion)")
      ch <- as.matrix(cap[, oc]); storage.mode(ch) <- "integer"
      occ_lab <- if (!is.null(d$occasions)) as.character(unlist(d$occasions)) else oc
      used <- oc
    }
    ids <- cap[[idc]]
    ind_extra <- cap[, setdiff(names(cap), c(idc, used)), drop = FALSE]
  } else {
    if (fmt == "long") {
      oc <- d$occasion_col %||% "occasion"
      if (!oc %in% names(cap)) stop("long format needs occasion_col (column '", oc, "' not found; columns in the captures file: ", cols_txt, ")")
      lab <- as.character(cap[[oc]])
      num <- suppressWarnings(as.numeric(lab))
      levs <- if (!is.null(d$occasions)) as.character(unlist(d$occasions))
              else if (all(is.finite(num))) as.character(sort(unique(num))) else sort(unique(lab))
    } else if (fmt == "dates") {
      dc <- d$date_col %||% "date"
      if (!dc %in% names(cap)) stop("dates format needs date_col (column '", dc, "' not found; columns in the captures file: ", cols_txt, ")")
      dt <- as.Date(as.character(cap[[dc]]), format = d$date_format %||% "%Y-%m-%d")
      if (anyNA(dt)) stop(sprintf("%d capture dates could not be read with format '%s' (e.g. '%s')", sum(is.na(dt)),
                                  d$date_format %||% "%Y-%m-%d", as.character(cap[[dc]])[which(is.na(dt))[1]]))
      if (!is.null(d$season_months)) {
        inm <- as.integer(format(dt, "%m")) %in% unlist(d$season_months)
        if (any(!inm)) notes <- c(notes, sprintf("%d captures outside the sampling months (%s) were dropped.", sum(!inm), paste(unlist(d$season_months), collapse = ",")))
        cap <- cap[inm, , drop = FALSE]; dt <- dt[inm]
      }
      if (!is.null(d$occasion_breaks)) {
        br <- as.Date(as.character(unlist(d$occasion_breaks)))
        k <- findInterval(as.numeric(dt), as.numeric(br))
        if (any(k == 0)) stop("some captures are earlier than the first occasion_breaks date")
        lab <- format(br[k], "%Y-%m-%d"); levs <- format(br, "%Y-%m-%d")
      } else {
        u <- d$occasion_unit %||% "year"
        lab <- switch(u, year = format(dt, "%Y"), month = format(dt, "%Y-%m"), week = format(dt, "%Y-W%V"),
                      stop("occasion_unit must be year, month or week"))
        levs <- if (!is.null(d$occasions)) as.character(unlist(d$occasions)) else sort(unique(lab))
      }
    } else stop("data format must be wide, long or dates")
    if (any(!lab %in% levs)) stop("captures at occasions not listed in `occasions`: ", paste(utils::head(unique(lab[!lab %in% levs])), collapse = ", "))
    ids <- unique(cap[[idc]])
    ch <- matrix(0L, length(ids), length(levs))
    ch[cbind(match(cap[[idc]], ids), match(lab, levs))] <- 1L
    nd <- sum(duplicated(paste(cap[[idc]], lab)))
    if (nd) notes <- c(notes, sprintf("%d repeat captures of the same individual within an occasion were collapsed (one detection per occasion).", nd))
    occ_lab <- levs
  }
  T <- ncol(ch)
  # time intervals
  iv <- if (!is.null(d$time_intervals)) as.numeric(unlist(d$time_intervals)) else {
    num <- suppressWarnings(as.numeric(occ_lab))
    if (all(is.finite(num)) && (fmt != "wide" || !is.null(d$occasions)) && all(diff(num) > 0) && (d$occasion_unit %||% "year") == "year") diff(num) else rep(1, T - 1)
  }
  if (length(iv) != T - 1) stop(sprintf("time_intervals needs %d values (one per interval between %d occasions)", T - 1, T))
  if (any(iv != 1)) notes <- c(notes, sprintf("Unequal time intervals (%s): survival is reported per unit time (per %s).", paste(iv, collapse = ", "), d$occasion_unit %||% "year"))
  # individual table
  ind <- data.frame(row.names = seq_along(ids))
  if (!is.null(d$individuals)) {
    it <- utils::read.csv(cjs_path(config, d$individuals), stringsAsFactors = FALSE, check.names = FALSE)
    if (!idc %in% names(it)) stop(sprintf("individuals file has no id column '%s'", idc))
    it[[idc]] <- trimws(as.character(it[[idc]]))
    if (anyDuplicated(it[[idc]])) stop("individuals file has duplicated ids: ", paste(utils::head(it[[idc]][duplicated(it[[idc]])]), collapse = ", "))
    miss <- setdiff(ids, it[[idc]])
    if (length(miss)) stop(sprintf("%d captured individuals are not in the individuals file (e.g. %s)", length(miss), paste(utils::head(miss, 5), collapse = ", ")))
    nc <- setdiff(it[[idc]], ids)
    if (length(nc)) notes <- c(notes, sprintf("%d individuals in the individuals file were never captured and were ignored.", length(nc)))
    ind <- it[match(ids, it[[idc]]), setdiff(names(it), idc), drop = FALSE]
  } else if (!is.null(ind_extra)) {
    ind <- ind_extra
  } else if (fmt != "wide") {
    # individual-level columns repeated on capture rows (constant per id)
    cols <- unique(c(d$group_col, d$age_col, unlist(d$individual_covs)))
    if (length(cols)) {
      mc <- setdiff(cols, names(cap))
      if (length(mc)) stop(sprintf("column(s) %s named in the config not found in the captures file (columns: %s); give an individuals file or fix the names",
                                   paste0("'", mc, "'", collapse = ", "), cols_txt))
      ind <- cap[match(ids, cap[[idc]]), cols, drop = FALSE]
      for (k in cols) {
        nv <- tapply(cap[[k]], cap[[idc]], function(v) length(unique(v)))
        if (any(nv > 1)) stop(sprintf("column '%s' changes between captures of the same individual (e.g. %s); individual covariates must be constant", k, names(nv)[nv > 1][1]))
      }
    }
  }
  keepc <- c(d$group_col, d$age_col, unlist(d$individual_covs))
  miss_c <- setdiff(keepc, names(ind))
  if (length(miss_c)) stop(sprintf("column(s) %s named in the config (group_col / age_col / individual_covs) not found in the %s; available columns: %s",
                                   paste0("'", miss_c, "'", collapse = ", "), if (!is.null(d$individuals)) "individuals file" else "captures file",
                                   paste(names(ind), collapse = ", ")))
  if (!is.null(d$individual_covs) || !is.null(d$group_col) || !is.null(d$age_col)) ind <- ind[, intersect(keepc, names(ind)), drop = FALSE]
  rownames(ind) <- NULL
  if (!is.null(d$group_col)) {
    if (!d$group_col %in% names(ind)) stop("group_col '", d$group_col, "' not found")
    g <- as.character(ind[[d$group_col]])
    lv <- if (!is.null(d$group_levels)) as.character(unlist(d$group_levels)) else sort(unique(g[!is.na(g)]))
    odd <- setdiff(unique(g[!is.na(g) & g != ""]), lv)
    if (length(odd)) stop(sprintf("group column '%s' has value(s) %s that are not in group_levels (%s) for %d individuals (e.g. %s): correct the typo or add the level to group_levels",
                                  d$group_col, paste0("'", odd, "'", collapse = ", "), paste(lv, collapse = ", "),
                                  sum(g %in% odd), paste(utils::head(ids[g %in% odd], 3), collapse = ", ")))
    ind[[d$group_col]] <- factor(g, levels = lv)
  }
  # occasion covariates
  occ <- data.frame(row.names = seq_len(T))
  if (!is.null(d$occasion_covs)) {
    ot <- utils::read.csv(cjs_path(config, d$occasion_covs), stringsAsFactors = FALSE, check.names = FALSE)
    oc <- d$occasion_cov_col %||% names(ot)[1]
    if (!oc %in% names(ot)) stop(sprintf("occasion_cov_col '%s' not found in the occasion covariate file (columns: %s)", oc, paste(names(ot), collapse = ", ")))
    m <- match(occ_lab, as.character(ot[[oc]]))
    if (anyNA(m[-T])) stop("occasion covariate file lacks rows for occasions: ", paste(occ_lab[is.na(m)], collapse = ", "))
    occ <- ot[m, setdiff(names(ot), oc), drop = FALSE]; rownames(occ) <- NULL
  }
  fo <- config$formulas %||% d$formulas %||% list(phi = "~ 1", p = "~ 1")
  list(ch = ch, ids = ids, occasions = occ_lab, intervals = iv, ind = ind, occ = occ,
       group = d$group_col, age_col = d$age_col, formulas = fo, species = d$species %||% "individuals",
       units = d$covariate_units %||% list(), occasion_unit = d$occasion_unit %||% (if (fmt == "dates" || cjs_year_labels(occ_lab)) "year" else "occasion"),
       notes = notes, labels = config$figures$labels %||% list())
}

cjs_check_data <- function(data) {
  pr <- character()
  ch <- data$ch; T <- ncol(ch)
  if (!all(ch %in% c(0, 1))) pr <- c(pr, "capture histories must contain only 0 and 1")
  if (T < 3) pr <- c(pr, sprintf("only %d occasions: a CJS model needs at least 3 (survival and detection are not separable with 2)", T))
  nz <- rowSums(ch) == 0
  if (any(nz)) pr <- c(pr, sprintf("%d individuals have no captures (e.g. %s): remove them", sum(nz), paste(utils::head(data$ids[nz], 3), collapse = ", ")))
  if (anyDuplicated(data$ids)) pr <- c(pr, "duplicated individual ids")
  first <- cjs_first(ch)
  for (t in 2:T) if (sum(ch[, t]) == 0 && any(first < t, na.rm = TRUE))
    pr <- c(pr, sprintf("occasion '%s' has no captures at all: if it was not sampled, drop it from `occasions` and give time_intervals (the gap becomes a longer interval)", data$occasions[t]))
  if (length(data$intervals) != T - 1 || any(!is.finite(data$intervals)) || any(data$intervals <= 0))
    pr <- c(pr, "time intervals must be T-1 positive numbers")
  bad <- intersect(setdiff(c(names(data$ind), names(data$occ)), data$age_col), c(cjs_marked_reserved, "tsm", "td", "age_class"))
  if (length(bad)) pr <- c(pr, sprintf("column name(s) %s are reserved (they are built-in model variables); rename them, e.g. '%s_class'", paste(bad, collapse = ", "), bad[1]))
  fo <- data$formulas
  vars <- tryCatch(unique(c(all.vars(stats::as.formula(cjs_parse_formula(fo$phi %||% "~1")$formula)), all.vars(stats::as.formula(cjs_rhs(fo$p %||% "~1"))))),
                   error = function(e) { pr <<- c(pr, paste("formula problem:", conditionMessage(e))); character() })
  for (proc in c("phi", "p")) {
    f <- if (proc == "phi") tryCatch(cjs_parse_formula(fo$phi %||% "~1")$formula, error = function(e) NULL) else cjs_rhs(fo$p %||% "~1")
    if (!is.null(f)) tryCatch(cjs_check_vars(f, proc, data), error = function(e) pr <<- c(pr, conditionMessage(e)))
  }
  for (v in intersect(vars, names(data$ind))) if (anyNA(data$ind[[v]]))
    pr <- c(pr, sprintf("individual covariate '%s' is missing for %d individuals (e.g. %s): fill in, or drop those individuals", v, sum(is.na(data$ind[[v]])), paste(utils::head(data$ids[is.na(data$ind[[v]])], 3), collapse = ", ")))
  for (v in intersect(vars, names(data$occ))) {
    if (anyNA(data$occ[[v]][-T])) pr <- c(pr, sprintf("occasion covariate '%s' is missing for occasions %s", v, paste(data$occasions[-T][is.na(data$occ[[v]][-T])], collapse = ", ")))
  }
  if (!is.null(data$group)) {
    tab <- table(data$ind[[data$group]])
    if (anyNA(data$ind[[data$group]])) pr <- c(pr, sprintf("%d individuals have no value for group '%s'", sum(is.na(data$ind[[data$group]])), data$group))
    if (any(tab < 3)) pr <- c(pr, sprintf("group level(s) %s have fewer than 3 individuals: merge or drop them", paste(names(tab)[tab < 3], collapse = ", ")))
  }
  if (!is.null(data$age_col)) {
    a <- tolower(as.character(data$ind[[data$age_col]]))
    ok <- a %in% c("adult", "ad", "a", "juvenile", "juv", "j", "young", "hy", "ahy", "asy", "subadult")
    if (any(!ok)) pr <- c(pr, sprintf("age column '%s' has values other than adult/juvenile: %s", data$age_col, paste(utils::head(unique(a[!ok])), collapse = ", ")))
  }
  pr
}

# Write simulated data in the researcher format (long captures + tables).
cjs_write_data <- function(data, out) {
  w <- which(data$ch > 0, arr.ind = TRUE)
  utils::write.csv(data.frame(id = data$ids[w[, 1]], occasion = data$occasions[w[, 2]])[order(w[, 1], w[, 2]), ],
                   file.path(out, "cjs_sim_captures.csv"), row.names = FALSE)
  utils::write.csv(cbind(id = data$ids, data$ind), file.path(out, "cjs_sim_individuals.csv"), row.names = FALSE)
  utils::write.csv(cbind(occasion = data$occasions, data$occ), file.path(out, "cjs_sim_occasions.csv"), row.names = FALSE)
}

# Occasion labels that are calendar years (e.g. 1981, 1982, ...): survival is then annual.
cjs_year_labels <- function(lab) {
  num <- suppressWarnings(as.numeric(lab))
  length(num) > 1 && all(is.finite(num)) && all(num == round(num)) && all(num >= 1800 & num <= 2200) && all(diff(num) > 0)
}
