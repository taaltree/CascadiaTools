# Shared helpers for the timeseries module: data reading/checking, the
# year x site matrix every engine uses, formatting.
#
# Data structure (what simulate()$data and read_data() both return):
#   list(model      = "exponential" | "gompertz" | "marss" | "bayes" | "loglinear" (or NULL),
#        df         = data.frame(year, site, count, <covariates...>) one row per
#                     site-year; count NA = not surveyed that year,
#        covariates = character names of covariate columns (raw units),
#        response   = "count" | "index" (index: any positive abundance index),
#        meta       = list(source, count_label, units, ...))

`%||%` <- function(a, b) if (is.null(a)) b else a

ts_fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
ts_sig <- function(x, d = 3) format(signif(x, d), big.mark = ",", scientific = FALSE, trim = TRUE)
ts_pct <- function(x, d = 1) sprintf(paste0("%+.", d, "f%%"), x)

TS_MODELS <- c("exponential", "gompertz", "marss", "bayes", "loglinear")

ts_model <- function(data, settings = list()) {
  m <- settings$model %||% data$model %||% "exponential"
  if (!m %in% TS_MODELS)
    stop("settings: model must be one of ", paste(TS_MODELS, collapse = ", "), call. = FALSE)
  m
}

ts_sites <- function(df) {
  s <- unique(as.character(df$site))
  s[order(s)]
}

# Full year sequence (min..max) x sites matrix of the response. Years absent
# from the file are inserted as missing (NA) so engines see the true spacing.
ts_matrix <- function(df) {
  yrs <- seq(min(df$year), max(df$year))
  sites <- ts_sites(df)
  Y <- matrix(NA_real_, length(sites), length(yrs), dimnames = list(sites, yrs))
  for (i in seq_len(nrow(df))) Y[as.character(df$site[i]), as.character(df$year[i])] <- df$count[i]
  list(Y = Y, years = yrs, sites = sites)
}

# Covariates per year (same for all sites; single-site models). Standardised;
# centre/scale stored so figures can return to original units.
ts_cov_matrix <- function(df, covs, years) {
  if (!length(covs)) return(list(Z = matrix(0, length(years), 0), scaling = list(), raw = NULL))
  raw <- matrix(NA_real_, length(years), length(covs), dimnames = list(years, covs))
  for (k in covs) {
    v <- tapply(df[[k]], df$year, function(z) mean(z, na.rm = TRUE))
    raw[names(v), k] <- v
  }
  sc <- lapply(covs, function(k) list(center = mean(raw[, k], na.rm = TRUE),
                                      scale = stats::sd(raw[, k], na.rm = TRUE)))
  names(sc) <- covs
  Z <- raw
  for (k in covs) Z[, k] <- (raw[, k] - sc[[k]]$center) / sc[[k]]$scale
  list(Z = Z, scaling = sc, raw = raw)
}

ts_col <- function(df, col, what) {
  if (is.null(col) || !col %in% names(df))
    stop(sprintf("Column '%s' (%s) not found in the data file. Columns present: %s. Set data: %s in the config.",
                 col %||% "<unset>", what, paste(names(df), collapse = ", "), what), call. = FALSE)
  df[[col]]
}

ts_num <- function(x, col) {
  if (is.numeric(x)) return(as.numeric(x))
  s <- trimws(as.character(x)); s[s %in% c("", "NA", "na", "-", ".")] <- NA
  v <- suppressWarnings(as.numeric(s))
  bad <- which(!is.na(s) & is.na(v))
  if (length(bad))
    stop(sprintf("Column '%s' has non-numeric values in data row(s) %s (spreadsheet line %s, counting the header) (e.g. '%s'). Fix them, or leave the cell blank for a year that was not surveyed.",
                 col, paste(utils::head(bad, 5), collapse = ", "), paste(utils::head(bad, 5) + 1, collapse = ", "), s[bad[1]]), call. = FALSE)
  v
}

ts_read_data <- function(config) {
  dc <- config$data %||% list()
  f <- dc$file %||% stop("config: data: file: is required", call. = FALSE)
  path <- if (file.exists(f)) f else file.path(config$.dir %||% ".", f)
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  year <- ts_num(ts_col(raw, dc$year %||% "year", "year"), dc$year %||% "year")
  count <- ts_num(ts_col(raw, dc$count %||% "count", "count"), dc$count %||% "count")
  site <- if (!is.null(dc$site)) as.character(ts_col(raw, dc$site, "site")) else rep("all", nrow(raw))
  covs <- unlist(dc$covariates %||% list())
  df <- data.frame(year = year, site = site, count = count)
  for (k in covs) df[[k]] <- ts_num(ts_col(raw, k, "covariates"), k)
  list(model = config$settings$model %||% dc$model, df = df, covariates = covs,
       response = dc$response %||% "count",
       meta = list(source = basename(path), count_label = dc$count_label %||% "Count",
                   covariate_labels = dc$covariate_labels %||% list(),
                   year_label = "Year"))
}

ts_check_data <- function(data, settings = list()) {
  df <- data$df; pr <- character()
  m <- tryCatch(ts_model(data, settings), error = function(e) { pr <<- c(pr, conditionMessage(e)); NULL })
  if (anyNA(df$year)) pr <- c(pr, "Some rows have no year. Every row needs a year.")
  if (any(df$year != round(df$year), na.rm = TRUE)) pr <- c(pr, "Years must be whole numbers (one survey per year). For seasonal data aggregate to one value per year first.")
  dup <- duplicated(df[, c("year", "site")])
  if (any(dup)) pr <- c(pr, sprintf("Year %s appears more than once for site '%s'. Give one count per site per year (average or sum replicate counts first).",
                                    df$year[dup][1], df$site[dup][1]))
  if (any(df$count < 0, na.rm = TRUE)) pr <- c(pr, "Negative counts found.")
  # Site names: catch typos ("Nrth Fork" vs "North Fork", "north fork ") that
  # would silently create an extra, very short site.
  sites <- unique(as.character(df$site))
  if (length(sites) > 1) {
    nyr_site <- vapply(sites, function(s) length(unique(df$year[df$site == s & !is.na(df$count)])), 1L)
    norm <- tolower(gsub("[^[:alnum:]]", "", sites))
    for (i in seq_along(sites)) for (j in seq_along(sites)) if (i < j) {
      same <- norm[i] == norm[j]
      d <- utils::adist(norm[i], norm[j])[1, 1]
      close <- d <= 2 && min(nchar(norm[c(i, j)])) >= 4
      few <- min(nyr_site[c(i, j)]) < 0.5 * max(nyr_site[c(i, j)])
      if (same || (close && few)) {
        k <- if (nyr_site[i] <= nyr_site[j]) i else j; o <- setdiff(c(i, j), k)
        pr <- c(pr, sprintf("Site '%s' (%d year(s) with counts) looks like a misspelling of site '%s' (%d years). If they are the same site, correct the name in the data file; otherwise rename one so the names are clearly different.",
                            sites[k], nyr_site[k], sites[o], nyr_site[o]))
      }
    }
    short <- sites[nyr_site < 3]
    if (length(short)) pr <- c(pr, sprintf("Site(s) %s have counts in fewer than 3 years. Each site needs at least 3 counted years (check the site names for typos, or drop the site).",
                                          paste0("'", short, "'", collapse = ", ")))
  }
  nyr <- length(unique(df$year[!is.na(df$count)]))
  if (nyr < 8) pr <- c(pr, sprintf("Only %d years with counts. State-space trend models need at least ~10 years (15+ to separate process and observation error).", nyr))
  if (!is.null(m) && m %in% c("exponential", "gompertz", "marss", "bayes") && any(df$count == 0, na.rm = TRUE)) {
    z <- df[!is.na(df$count) & df$count == 0, ]
    pr <- c(pr, sprintf("Count of 0 in year %s (site '%s'). The log-scale state-space models need positive counts; use model: loglinear (negative binomial counts) or set settings: zero_offset: 0.5 (adds 0.5 to all counts; biases small counts).",
                        z$year[1], z$site[1]))
  }
  if (!is.null(m) && m %in% c("exponential", "gompertz", "bayes") && length(unique(df$site)) > 1)
    pr <- c(pr, sprintf("Model '%s' fits one population; the file has %d sites. Use model: marss (several sub-populations) or loglinear, or sum the sites first.",
                        m, length(unique(df$site))))
  for (k in data$covariates) {
    if (!is.null(m) && m == "bayes") {
      yrs <- seq(min(df$year), max(df$year))
      have <- unique(df$year[!is.na(df[[k]])])
      miss <- setdiff(yrs[-1], have)
      if (length(miss)) pr <- c(pr, sprintf("Covariate '%s' is missing in year(s) %s. Growth into every year needs a covariate value, even in years without a count: add rows with the covariate and a blank count.",
                                            k, paste(utils::head(miss, 6), collapse = ", ")))
    } else if (anyNA(df[[k]][!is.na(df$count)])) pr <- c(pr, sprintf("Covariate '%s' has missing values in surveyed years.", k))
    if (stats::sd(df[[k]], na.rm = TRUE) == 0) pr <- c(pr, sprintf("Covariate '%s' never changes.", k))
  }
  if (!is.null(m) && m %in% c("exponential", "gompertz", "marss") && length(data$covariates))
    pr <- c(pr, sprintf("Model '%s' does not use covariates; use model: bayes (covariates on growth rate) or loglinear (covariates on abundance), or remove data: covariates.", m))
  pr
}

# log response matrix with optional zero offset
ts_logY <- function(Y, settings) {
  off <- settings$zero_offset %||% 0
  log(Y + off)
}

# Quantile summary of a draws matrix into contract rows.
ts_draw_rows <- function(draws, est = NULL, prob = 0.95) {
  draws <- as.matrix(draws)
  a <- (1 - prob) / 2
  qs <- t(apply(draws, 2, function(v) stats::quantile(v, c(a, 1 - a), na.rm = TRUE, names = FALSE)))
  data.frame(param = colnames(draws),
             estimate = if (is.null(est)) colMeans(draws, na.rm = TRUE) else unname(est[colnames(draws)]),
             sd = apply(draws, 2, stats::sd, na.rm = TRUE),
             lower = qs[, 1], upper = qs[, 2], row.names = NULL)
}

ts_rows <- function(param, estimate, sd, lower, upper)
  data.frame(param = param, estimate = estimate, sd = sd, lower = lower, upper = upper, row.names = NULL)

ts_est <- function(fit, p) { e <- fit$estimates; e[e$param == p, , drop = FALSE] }
