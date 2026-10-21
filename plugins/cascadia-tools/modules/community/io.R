# Reading researchers' files for the community module.
#
# Species data (config data$species), one of:
#   format: wide  - one row per site, first column (site_col) = site id, one
#                   column per species (counts, cover, biomass or 0/1)
#   format: long  - one row per site x species record: site_col, species_col,
#                   count_col (missing combinations are zeros). Surveys that
#                   were repeated at a site are summed.
# Site data (config data$env): one row per site with env_site_col (default:
# site_col), the environmental variables (env_vars) and optionally a grouping
# column (group: habitat type, treatment, region...).

com_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc$species)) stop("config data$species (site x species CSV) is required")
  path <- function(f) if (grepl("^/", f)) f else file.path(config$.dir %||% ".", f)
  sp <- utils::read.csv(path(dc$species), check.names = FALSE, stringsAsFactors = FALSE)
  site_col <- dc$site_col %||% names(sp)[1]
  fmt <- dc$format %||% "wide"
  notes <- character(); problems <- character()
  if (!site_col %in% names(sp)) stop("site column '", site_col, "' not found in ", dc$species,
                                     ". Its first columns are: ", paste(utils::head(names(sp), 6), collapse = ", "),
                                     ". Set data$site_col to the right name (names are case-sensitive).")
  if (fmt == "long") {
    spc <- dc$species_col %||% "species"; cc <- dc$count_col %||% "count"
    miss <- setdiff(c(spc, cc), names(sp))
    if (length(miss)) stop("long-format species file lacks column(s): ", paste(miss, collapse = ", "))
    sp[[cc]] <- suppressWarnings(as.numeric(sp[[cc]]))
    Y <- stats::xtabs(stats::as.formula(paste0("`", cc, "` ~ `", site_col, "` + `", spc, "`")), data = sp)
    Y <- matrix(as.numeric(Y), nrow(Y), ncol(Y), dimnames = dimnames(Y))
  } else {
    ids <- as.character(sp[[site_col]])
    if (anyDuplicated(ids)) problems <- c(problems, sprintf(
      "Site id(s) appear more than once in the species file: %s. Combine repeated surveys (sum or average) or use format: long.",
      paste(unique(ids[duplicated(ids)]), collapse = ", ")))
    M <- sp[, setdiff(names(sp), c(site_col, unlist(dc$drop_cols))), drop = FALSE]
    nonnum <- names(M)[!vapply(M, function(v) is.numeric(v) || all(is.na(suppressWarnings(as.numeric(v))) == is.na(v)), TRUE)]
    if (length(nonnum)) problems <- c(problems, sprintf(
      "Column(s) %s in the species file contain text; every column except the site id must be a number (use drop_cols to ignore notes columns).",
      paste(nonnum, collapse = ", ")))
    M[] <- lapply(M, function(v) suppressWarnings(as.numeric(v)))
    Y <- as.matrix(M); rownames(Y) <- ids
  }
  env <- NULL
  env_vars <- unlist(dc$env_vars)
  group <- dc$group
  if (!is.null(dc$env)) {
    env <- utils::read.csv(path(dc$env), stringsAsFactors = FALSE, check.names = FALSE)
    esc <- dc$env_site_col %||% site_col
    if (!esc %in% names(env)) stop("site column '", esc, "' not found in ", dc$env, ". Its columns are: ",
                                   paste(names(env), collapse = ", "), ". Set data$env_site_col (names are case-sensitive).")
    env[[esc]] <- as.character(env[[esc]])
    if (is.null(env_vars)) env_vars <- setdiff(names(env)[vapply(env, is.numeric, TRUE)], esc)
    miss <- setdiff(c(env_vars, group), names(env))
    if (length(miss)) stop("environment file ", dc$env, " lacks column(s): ", paste(miss, collapse = ", "),
                           ". Its columns are: ", paste(names(env), collapse = ", "),
                           ". Fix env_vars / group in the config (names are case-sensitive).")
    only_sp <- setdiff(rownames(Y), env[[esc]]); only_env <- setdiff(env[[esc]], rownames(Y))
    if (length(only_sp)) problems <- c(problems, sprintf(
      "Site(s) %s are in the species file but not in the environment file. Check spelling of site ids (e.g. 'S01' vs 'S1').",
      paste(utils::head(only_sp, 10), collapse = ", ")))
    if (length(only_env)) notes <- c(notes, sprintf(
      "Site(s) %s are in the environment file but have no species data; they were ignored.", paste(utils::head(only_env, 10), collapse = ", ")))
    keep <- intersect(rownames(Y), env[[esc]])
    Y <- Y[keep, , drop = FALSE]
    env <- env[match(keep, env[[esc]]), , drop = FALSE]
    names(env)[names(env) == esc] <- "site"
    if (!is.null(group)) env[[group]] <- factor(env[[group]])
    for (v in env_vars) if (is.character(env[[v]])) env[[v]] <- factor(env[[v]])
  } else {
    env <- data.frame(site = rownames(Y)); env_vars <- character(); group <- NULL
  }
  rownames(env) <- NULL
  # data type
  dtype <- dc$data_type %||% "auto"
  vals <- Y[!is.na(Y)]
  if (dtype == "auto") dtype <- if (length(vals) && all(vals %in% c(0, 1))) "presence" else "abundance"
  # drop species never recorded (they carry no information)
  zero <- colnames(Y)[colSums(Y > 0, na.rm = TRUE) == 0]
  if (length(zero)) { notes <- c(notes, sprintf("%d species with no records were dropped: %s.", length(zero), paste(utils::head(zero, 10), collapse = ", ")))
                      Y <- Y[, !colnames(Y) %in% zero, drop = FALSE] }
  minocc <- as.integer(dc$min_sites %||% 0)
  if (minocc > 1) {
    rare <- colSums(Y > 0, na.rm = TRUE) < minocc
    if (any(rare)) { notes <- c(notes, sprintf("%d species found at fewer than %d sites were dropped (min_sites).", sum(rare), minocc))
                     Y <- Y[, !rare, drop = FALSE] }
  }
  # species names that look like typos of each other (case/punctuation, or one edit apart)
  spn <- colnames(Y); key <- toupper(gsub("[^[:alnum:]]", "", spn))
  dup <- unique(key[duplicated(key)])
  if (length(dup)) problems <- c(problems, sprintf(
    "Species names that differ only in case, spaces or punctuation: %s. Merge them (same species) or rename them.",
    paste(vapply(dup, function(k) paste(spn[key == k], collapse = " / "), ""), collapse = "; ")))
  if (length(spn) > 1 && length(spn) <= 2000) {
    dm <- utils::adist(key); dm[lower.tri(dm, diag = TRUE)] <- NA
    nr <- nchar(key)[row(dm)]; nc <- nchar(key)[col(dm)]
    # short codes (e.g. 4-letter bird codes) routinely differ by one letter: only flag an extra/missing
    # character there, or a one-letter difference between longer names
    near <- which(dm == 1 & pmin(nr, nc) >= 3 & (nr != nc | pmin(nr, nc) >= 7), arr.ind = TRUE)
    if (nrow(near)) notes <- c(notes, sprintf(
      "Possible species-name typos (one character apart): %s. If these are the same species, fix the names and rerun.",
      paste(utils::head(sprintf("%s / %s", spn[near[, 1]], spn[near[, 2]]), 10), collapse = "; ")))
  }
  part <- dc$partition
  if (is.null(part)) part <- stats::setNames(as.list(env_vars), make.names(env_vars))
  names(part) <- make.names(names(part))
  list(comm = Y, env = env, env_vars = env_vars, group = group, partition = part,
       data_type = dtype, notes = notes, problems = problems,
       # report text uses data$labels; fall back on the figures: labels block so the
       # interpretation names species/levels the same way the figures do
       labels = utils::modifyList(as.list(config$figures$labels %||% list()), as.list(dc$labels %||% list())))
}

com_check_data <- function(data) {
  pr <- data$problems %||% character()
  Y <- data$comm
  if (anyNA(Y)) {
    w <- which(is.na(Y), arr.ind = TRUE)
    pr <- c(pr, sprintf("Missing species values (NA) at %d cells, e.g. site %s, species %s. Use 0 for 'not recorded'; drop sites that were not surveyed.",
                        nrow(w), rownames(Y)[w[1, 1]], colnames(Y)[w[1, 2]]))
  }
  if (any(Y < 0, na.rm = TRUE)) pr <- c(pr, "Negative species values found; counts, cover and biomass must be >= 0.")
  empty <- rownames(Y)[rowSums(Y, na.rm = TRUE) == 0]
  if (length(empty)) pr <- c(pr, sprintf(
    "Site(s) %s recorded no species. Hellinger, Bray-Curtis and Sorensen are undefined for empty sites: remove them (and report how many) before analysis.",
    paste(utils::head(empty, 10), collapse = ", ")))
  if (nrow(Y) < 8) pr <- c(pr, sprintf("Only %d sites: multivariate tests and ordinations need at least about 10 sites.", nrow(Y)))
  if (ncol(Y) < 3) pr <- c(pr, "Fewer than 3 species: use a univariate model (glm module) instead.")
  if (identical(data$data_type, "presence") && !all(Y[!is.na(Y)] %in% c(0, 1)))
    pr <- c(pr, "data_type is 'presence' but values other than 0/1 occur.")
  env <- data$env
  for (v in data$env_vars) {
    if (anyNA(env[[v]])) pr <- c(pr, sprintf("Environmental variable '%s' has missing values at site(s) %s.", v,
                                            paste(utils::head(env$site[is.na(env[[v]])], 10), collapse = ", ")))
    else if (is.numeric(env[[v]]) && stats::sd(env[[v]]) == 0) pr <- c(pr, sprintf("Environmental variable '%s' is constant.", v))
  }
  if (length(data$env_vars) > 1) {
    num <- data$env_vars[vapply(data$env_vars, function(v) is.numeric(env[[v]]), TRUE)]
    if (length(num) > 1) {
      cm <- suppressWarnings(stats::cor(env[, num], use = "pairwise"))
      hi <- which(abs(cm) > 0.95 & upper.tri(cm), arr.ind = TRUE)
      if (nrow(hi)) pr <- c(pr, sprintf("Environmental variables %s and %s are almost identical (r = %.2f); keep one.",
                                       num[hi[1, 1]], num[hi[1, 2]], cm[hi[1, 1], hi[1, 2]]))
    }
  }
  for (k in names(data$partition)) {
    miss <- setdiff(unlist(data$partition[[k]]), names(env))
    if (length(miss)) pr <- c(pr, sprintf("partition set '%s' names unknown variable(s): %s.", k, paste(miss, collapse = ", ")))
  }
  if (!is.null(data$group)) {
    g <- env[[data$group]]
    if (anyNA(g)) pr <- c(pr, sprintf("Group column '%s' has missing values.", data$group))
    tb <- table(g)
    if (length(tb) < 2) pr <- c(pr, sprintf("Group column '%s' has only one level; PERMANOVA needs at least two groups.", data$group))
    small <- names(tb)[tb < 3]
    if (length(small)) pr <- c(pr, sprintf("Group(s) %s have fewer than 3 sites; merge or drop them.", paste(small, collapse = ", ")))
  }
  pr
}

com_write_data <- function(data, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  Y <- data.frame(site = rownames(data$comm), data$comm, check.names = FALSE)
  utils::write.csv(Y, file.path(out, "community_species.csv"), row.names = FALSE)
  utils::write.csv(data$env, file.path(out, "community_sites.csv"), row.names = FALSE)
  invisible(out)
}
