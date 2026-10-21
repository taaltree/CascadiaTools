# ipm: reading researchers' files, data checks, m-array construction.
#
# Three data sets (Schaub & Kery 2022, ch. 4-5):
#   counts       one row per year: year, count (breeding females / pairs / territories)
#   captures     capture-recapture of individually marked animals, either
#                long (one row per capture or resighting: id, year, age at that
#                capture) or wide (one row per animal: id, age at marking,
#                one 0/1 column per year)
#   productivity one row per monitored brood (year, fledglings[, mother_age])
#                or one row per year (year, broods, fledglings[, mother_age])
#   covariates   optional, one row per year
# Years present in any file define the study period (contiguous).

ipm_age_code <- function(x, what = c("capture", "mother")) {
  what <- match.arg(what)
  v <- gsub("[.]+$", "", tolower(trimws(as.character(x))))   # "juv." -> "juv"
  out <- rep(NA_integer_, length(v))
  if (what == "capture") {
    # 1 = marked as young of the year (first-year survival applies until the
    # next breeding season); 2 = marked as a breeder (yearling or older).
    # A bird first caught as a breeding yearling ("SY") has already survived
    # its first year, so it is an adult for the capture-recapture model.
    out[v %in% c("1", "0", "j", "juv", "juvenile", "juveniles", "nestling", "nestlings", "chick", "chicks", "pullus",
                 "pulli", "young", "hatch-year", "hy", "fledgling", "fledglings", "hatchling", "calf", "pup", "fawn", "yoy")] <- 1L
    out[v %in% c("2", "2+", "a", "ad", "adult", "adults", "breeder", "yearling", "sy", "second-year", "asy", "ahy",
                 "after-hatch-year", "after-second-year", "older", "old")] <- 2L
  } else {
    # mother age class for fecundity: 1 = first-time (1-year-old) breeder, 2 = older
    out[v %in% c("1", "yearling", "sy", "second-year", "1y", "first-time", "first", "young", "1-year-old")] <- 1L
    out[v %in% c("2", "2+", "a", "ad", "adult", "adults", "older", "asy", "after-second-year", "experienced", "old")] <- 2L
  }
  out
}

# Age-dependent m-arrays (Burnham et al. 1987; Schaub & Kery 2022 4.5.1.2).
# age = 1: marked as juvenile (young of the year); 2: marked as adult. After
# its first recapture a juvenile is re-released as an adult. Rows: release
# occasion 1..T-1; columns 1..T-1: first recapture at occasion 2..T; column
# T: never recaptured.
ipm_marray <- function(ch, age) {
  T <- ncol(ch)
  mj <- matrix(0, T - 1, T); ma <- matrix(0, T - 1, T)
  for (i in seq_len(nrow(ch))) {
    occ <- which(ch[i, ] > 0)
    if (!length(occ)) next
    for (k in seq_along(occ)) {
      r <- occ[k]
      if (r >= T) next
      juv <- (k == 1 && age[i] == 1)
      col <- if (k < length(occ)) occ[k + 1] - 1 else T
      if (juv) mj[r, col] <- mj[r, col] + 1 else ma[r, col] <- ma[r, col] + 1
    }
  }
  list(marr_j = mj, marr_a = ma)
}

ipm_read_csv <- function(config, key) {
  f <- config$data[[key]]
  if (is.null(f)) return(NULL)
  p <- if (file.exists(f)) f else file.path(config$.dir, f)
  if (!file.exists(p)) stop("cannot find the ", key, " file '", f, "' (looked in ", config$.dir, ")")
  utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
}

ipm_col <- function(df, name, file, default = NULL, key = NULL) {
  nm <- if (is.null(name)) default else name
  if (is.null(nm) || !nm %in% names(df))
    stop("the ", file, " file has no column '", nm, "'. Columns found: ", paste(names(df), collapse = ", "),
         ". Rename the column in the file, or name the column to use in the config under data: ",
         if (is.null(key)) "(README, section Config, lists every *_col key)" else
           paste0("e.g.  ", key, ": \"<column name in your file>\""), ".")
  df[[nm]]
}

ipm_read_data <- function(config) {
  d <- config$data
  cnt <- ipm_read_csv(config, "counts")
  cap <- ipm_read_csv(config, "captures")
  pro <- ipm_read_csv(config, "productivity")
  cov <- ipm_read_csv(config, "covariates")
  if (is.null(cnt) || is.null(cap) || is.null(pro))
    stop("an IPM needs three files: data: counts, captures and productivity")
  cy <- as.integer(ipm_col(cnt, d$year_col, "counts", "year", "year_col"))
  cc <- as.numeric(ipm_col(cnt, d$count_col, "counts", "count", "count_col"))
  # captures
  fmt <- if (is.null(d$capture_format)) "long" else d$capture_format
  idc <- if (is.null(d$id_col)) "id" else d$id_col
  agc <- if (is.null(d$age_col)) "age" else d$age_col
  if (fmt == "marray") {
    # m-array table (as written by write_data / IPMbook::marrayAge): columns
    # age (juvenile/adult), release_year, one column per recapture year, never
    ag <- ipm_age_code(ipm_col(cap, agc, "captures", key = "age_col"))
    ry <- as.integer(ipm_col(cap, if (is.null(d$release_year_col)) "release_year" else d$release_year_col, "captures", key = "release_year_col"))
    rc <- grep("[0-9]{4}$", names(cap), value = TRUE)
    nv <- if ("never" %in% names(cap)) "never" else stop("m-array file needs a 'never' column (released but never recaptured)")
    myrs <- sort(unique(c(ry, as.integer(sub("^[^0-9]*", "", rc)))))
  }
  if (fmt == "long") {
    cid <- as.character(ipm_col(cap, idc, "captures", key = "id_col"))
    cyr <- as.integer(ipm_col(cap, if (is.null(d$capture_year_col)) "year" else d$capture_year_col, "captures", key = "capture_year_col"))
    craw <- as.character(ipm_col(cap, agc, "captures", key = "age_col"))
    cag <- ipm_age_code(craw)
  } else if (fmt == "wide") {
    ycols <- grep("^[^0-9]*[0-9]{4}$", names(cap), value = TRUE)
    if (!length(ycols)) stop("wide capture file: expected one column per year named like 2005 or y2005")
    yrs <- as.integer(sub("^[^0-9]*", "", ycols))
    long <- do.call(rbind, lapply(seq_along(ycols), function(k) {
      w <- which(as.numeric(cap[[ycols[k]]]) > 0)
      data.frame(id = as.character(cap[[idc]][w]), year = rep(yrs[k], length(w)),
                 first_age = ipm_age_code(cap[[agc]][w]), raw_age = as.character(cap[[agc]][w]))
    }))
    cid <- long$id; cyr <- long$year
    cag <- long$first_age; craw <- long$raw_age
    # wide: age column is age at first marking -> only first capture counts
  } else if (fmt != "marray") stop("data: capture_format must be long, wide or marray")
  # productivity
  pfmt <- if (!is.null(d$productivity_format)) d$productivity_format else
    if (any(c(d$broods_col, "broods", "n_broods") %in% names(pro))) "aggregated" else "broods"
  py <- as.integer(ipm_col(pro, d$productivity_year_col, "productivity", "year", "productivity_year_col"))
  pf <- as.numeric(ipm_col(pro, d$fledglings_col, "productivity", "fledglings", "fledglings_col"))
  pb <- if (pfmt == "aggregated") as.numeric(ipm_col(pro, d$broods_col, "productivity",
                                                    if ("broods" %in% names(pro)) "broods" else "n_broods", "broods_col")) else rep(1, length(pf))
  mac <- if (is.null(d$mother_age_col)) "mother_age" else d$mother_age_col
  page <- if (mac %in% names(pro)) ipm_age_code(pro[[mac]], "mother") else rep(NA_integer_, length(pf))
  # study period
  # duplicated count years: ambiguous, the researcher must decide
  dupc <- unique(cy[duplicated(cy) & !is.na(cy)])
  if (length(dupc))
    stop("the counts file has more than one row for year(s) ", paste(dupc, collapse = ", "), " (values: ",
         paste(vapply(dupc, function(y) paste(cc[cy %in% y], collapse = " and "), ""), collapse = "; "),
         "). Keep one row per year (if the rows are separate plots or observers, combine them into one annual count first).")
  yr <- sort(unique(c(cy, if (fmt == "marray") myrs else cyr, py)))
  if (!is.null(d$first_year)) yr <- yr[yr >= d$first_year]
  if (!is.null(d$last_year)) yr <- yr[yr <= d$last_year]
  # a mistyped year (2106 for 2016) silently stretches the study period:
  # stop when the data years fall into blocks separated by > 3 empty years
  if (length(yr) > 1 && any(diff(yr) > 4)) {
    blk <- cumsum(c(1, diff(yr) > 4)); main <- as.integer(names(which.max(table(blk))))
    odd <- yr[blk != main]
    where <- function(y) paste(c(if (any(cy %in% y)) "counts", if (fmt != "marray" && any(cyr %in% y)) "captures",
                                 if (fmt == "marray" && any(myrs %in% y)) "captures", if (any(py %in% y)) "productivity"), collapse = "/")
    stop("year(s) ", paste(sprintf("%d (%s file)", odd, vapply(odd, where, "")), collapse = ", "),
         " lie far outside the other years (", min(yr[blk == main]), "-", max(yr[blk == main]),
         "). Probably a typing error: correct the year, or if the gap is real, restrict the study period with data: first_year / last_year.")
  }
  years <- seq(min(yr), max(yr))
  T <- length(years)
  counts <- rep(NA_real_, T); counts[match(cy, years)[!is.na(match(cy, years))]] <- cc[!is.na(match(cy, years))]
  if (fmt == "marray") {
    mj <- matrix(0, T - 1, T); mA <- matrix(0, T - 1, T)
    for (r in seq_len(nrow(cap))) {
      i <- match(ry[r], years)
      if (is.na(i) || i >= T || is.na(ag[r])) next
      v <- numeric(T)
      for (k in rc) { j <- match(as.integer(sub("^[^0-9]*", "", k)), years); if (!is.na(j) && j > i) v[j - 1] <- v[j - 1] + as.numeric(cap[[k]][r]) }
      v[T] <- as.numeric(cap[[nv]][r])
      if (ag[r] == 1) mj[i, ] <- mj[i, ] + v else mA[i, ] <- mA[i, ] + v
    }
    ma <- list(marr_j = mj, marr_a = mA)
    # juvenile releases are first markings; adult releases mix new marks and re-releases
    first_age <- integer(); nm_j <- sum(mj); nm_a <- NA_real_
    bad_age <- 0; ch <- NULL; typo <- character(); bad_vals <- character(); bad_ids <- character()
  } else {
  # capture histories
  keep <- cyr %in% years & !is.na(cid)
  ids <- unique(cid[keep])
  ch <- matrix(0L, length(ids), T)
  ch[cbind(match(cid[keep], ids), match(cyr[keep], years))] <- 1L
  # age at first capture
  first_age <- rep(NA_integer_, length(ids))
  o <- order(match(cid[keep], ids), cyr[keep])
  ki <- match(cid[keep], ids)[o]; ka <- cag[keep][o]
  first_age[ki[!duplicated(ki)]] <- ka[!duplicated(ki)]
  bad_age <- sum(is.na(first_age))
  bad_vals <- unique(craw[keep][o][!duplicated(ki)][is.na(ka[!duplicated(ki)])])
  bad_ids <- ids[is.na(first_age)]
  first_age[is.na(first_age)] <- 2L
  # possible ring-reading / typing errors: ids whose letter/digit format is rare
  # in this file (e.g. 'AF7339l' among 'AF99999'-style rings: l/1, O/0 mix-ups,
  # a dropped or extra character), reported with the closest common id. A
  # misread ring splits one bird into a 'death' plus a new adult.
  typo <- character()
  if (length(ids) > 20) {
    pat <- gsub("[A-Za-z]", "A", gsub("[0-9]", "9", ids))
    tb <- table(pat)
    rare <- names(tb)[tb <= max(2, 0.01 * length(ids))]
    odd <- ids[pat %in% rare]
    common <- ids[!pat %in% rare]
    for (x in utils::head(odd, 20)) {
      dd <- if (length(common)) utils::adist(x, common)[1, ] else numeric()
      typo <- c(typo, if (length(dd) && min(dd) <= 1) sprintf("'%s' (close to '%s')", x, common[which.min(dd)]) else sprintf("'%s'", x))
    }
  }
  ma <- ipm_marray(ch, first_age)
  nm_j <- sum(first_age == 1); nm_a <- sum(first_age == 2)
  }
  # productivity table
  agef_req <- d$age_specific_fecundity
  have_age <- any(!is.na(page))
  agef <- if (is.null(agef_req) || identical(agef_req, "auto")) have_age && mean(!is.na(page)) > 0.8 else isTRUE(agef_req)
  pk <- py %in% years & !is.na(pf)
  n_na_fl <- sum(py %in% years & is.na(pf))
  bad_mage <- if (agef) unique(as.character(pro[[mac]])[pk & is.na(page)]) else character()
  n_bad_mage <- if (agef) sum(pk & is.na(page)) else 0
  tab <- function(sel) {
    J <- tapply(pf[sel], factor(py[sel], levels = years), sum); B <- tapply(pb[sel], factor(py[sel], levels = years), sum)
    J[is.na(J)] <- 0; B[is.na(B)] <- 0; list(J = as.numeric(J), B = as.numeric(B))
  }
  if (agef) {
    a1 <- tab(pk & page %in% 1L); a2 <- tab(pk & page %in% 2L)
    prod <- data.frame(year = years, J1 = a1$J, B1 = a1$B, J2 = a2$J, B2 = a2$B)
  } else {
    a <- tab(pk)
    prod <- data.frame(year = years, J1 = a$J, B1 = a$B, J2 = 0, B2 = 0)
  }
  # covariates
  covs <- NULL
  if (!is.null(cov)) {
    cyc <- as.integer(ipm_col(cov, d$covariate_year_col, "covariates", "year", "covariate_year_col"))
    covs <- data.frame(year = years)
    for (nm in setdiff(names(cov), c("year", d$covariate_year_col)))
      covs[[nm]] <- as.numeric(cov[[nm]])[match(years, cyc)]
  }
  eff <- if (is.null(d$effects)) list() else lapply(d$effects, unlist)
  list(species = if (is.null(d$species)) "the population" else d$species,
       years = years, counts = counts, marr_j = ma$marr_j, marr_a = ma$marr_a,
       prod = prod, age_fecundity = agef, covariates = covs, effects = eff,
       n_marked = c(juveniles = nm_j, adults = nm_a),
       count_unit = if (is.null(d$count_unit)) "breeding females" else d$count_unit,
       figure_labels = config$figures$labels,   # used for covariate names in the interpretation text
       notes = c(if (bad_age) sprintf("%d marked animals had an age at first capture that was not recognised (value(s): %s; e.g. ids %s) and were treated as adults (marked as breeders). If they were young of the year, recode the age column (juvenile / nestling / pullus vs adult) and refit.",
                                      bad_age, paste0("'", utils::head(bad_vals, 6), "'", collapse = ", "), paste(utils::head(bad_ids, 5), collapse = ", ")),
                 if (length(typo)) sprintf("%d animal id(s) have an unusual format compared with the other ids: %s%s. If these are typing or reading errors (l/1, O/0, missing character), correct them: a misread ring turns one bird into a 'death' plus a new adult and biases adult survival down.",
                                           length(typo), paste(utils::head(typo, 5), collapse = "; "), if (length(typo) > 5) "; ..." else ""),
                 if (n_bad_mage) sprintf("%d brood(s) had a mother age that was not recognised (value(s): %s) and were left out of the productivity data. Use yearling / first-time vs adult / older in the mother-age column and refit.",
                                         n_bad_mage, paste0("'", utils::head(bad_mage, 6), "'", collapse = ", ")),
                 if (n_na_fl) sprintf("%d brood(s) with no fledgling number (empty/NA) were left out of the productivity data (year(s) %s). If they failed, enter 0: failed broods must be included.",
                                      n_na_fl, paste(unique(py[is.na(pf) & py %in% years]), collapse = ", ")),
                 if (have_age && !agef) "Mother age was given for too few broods (< 80%): fecundity is modelled as the same for yearlings and adults."),
       raw = list(ch = ch, age = first_age, brood_year = py, fledglings = pf, broods = pb, mother_age = page))
}

ipm_check_data <- function(data) {
  pr <- character()
  T <- length(data$years)
  if (T < 5) pr <- c(pr, sprintf("only %d years: an IPM needs at least 5 (better 10+) years of counts to estimate population change and immigration.", T))
  if (sum(!is.na(data$counts)) < 3) pr <- c(pr, "fewer than 3 years have a population count")
  if (any(data$counts < 0, na.rm = TRUE)) pr <- c(pr, "negative counts found")
  if (sum(data$marr_j[, -ncol(data$marr_j)]) == 0)
    pr <- c(pr, "no animal marked as a juvenile was ever recaptured: juvenile (first-year) survival is not identifiable from capture data. Check the age column (juvenile/nestling vs adult).")
  if (sum(data$marr_a[, -ncol(data$marr_a)]) == 0)
    pr <- c(pr, "no adult was ever recaptured: adult survival cannot be estimated.")
  if (sum(data$prod$B1) + sum(data$prod$B2) == 0) pr <- c(pr, "the productivity file contains no broods within the study years")
  if (any(data$prod$J1 < 0 | data$prod$J2 < 0, na.rm = TRUE)) pr <- c(pr, "negative fledgling numbers in the productivity file")
  if (isTRUE(data$age_fecundity) && (sum(data$prod$B1) == 0 || sum(data$prod$B2) == 0))
    pr <- c(pr, "age-specific fecundity requested but there are no broods of yearling (or of older) mothers; set data: age_specific_fecundity: false")
  for (rate in names(data$effects)) {
    if (!rate %in% c("juvenile_survival", "adult_survival", "fecundity", "immigration"))
      pr <- c(pr, sprintf("effects: '%s' is not a rate; use juvenile_survival, adult_survival, fecundity or immigration", rate))
    for (cv in data$effects[[rate]]) {
      if (is.null(data$covariates) || !cv %in% names(data$covariates)) {
        pr <- c(pr, sprintf("covariate '%s' (effect on %s) is not a column of the covariates file%s. Check the spelling under data: effects: in the config.", cv, rate,
                            if (is.null(data$covariates)) " (no covariates file was given: data: covariates)" else
                              paste0(" (columns found: ", paste(setdiff(names(data$covariates), "year"), collapse = ", "), ")"))); next
      }
      x <- data$covariates[[cv]]
      need <- if (rate == "fecundity") seq_len(T) else seq_len(T - 1)
      if (anyNA(x[need])) pr <- c(pr, sprintf("covariate '%s' is missing for year(s) %s", cv,
                                             paste(data$years[need][is.na(x[need])], collapse = ", ")))
    }
  }
  pr
}

# write simulated data in the researcher-facing CSV formats
ipm_write_data <- function(data, out) {
  T <- length(data$years)
  utils::write.csv(data.frame(year = data$years, count = data$counts), file.path(out, "ipm_counts.csv"), row.names = FALSE)
  pr <- if (isTRUE(data$age_fecundity))
    rbind(data.frame(year = data$years, mother_age = "yearling", broods = data$prod$B1, fledglings = data$prod$J1),
          data.frame(year = data$years, mother_age = "adult", broods = data$prod$B2, fledglings = data$prod$J2))
  else data.frame(year = data$years, broods = data$prod$B1, fledglings = data$prod$J1)
  utils::write.csv(pr[order(pr$year), ], file.path(out, "ipm_productivity_by_year.csv"), row.names = FALSE)
  m <- rbind(data.frame(age = "juvenile", release_year = data$years[-T], data$marr_j),
             data.frame(age = "adult", release_year = data$years[-T], data$marr_a))
  names(m)[-(1:2)] <- c(paste0("recap_", data$years[-1]), "never")
  utils::write.csv(m, file.path(out, "ipm_marray.csv"), row.names = FALSE)
  if (!is.null(data$covariates)) utils::write.csv(data$covariates, file.path(out, "ipm_covariates.csv"), row.names = FALSE)
  invisible(out)
}
