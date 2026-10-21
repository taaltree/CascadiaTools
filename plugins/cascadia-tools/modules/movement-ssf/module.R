# ecotools module `movement-ssf`: habitat selection and movement from GPS
# telemetry. Integrated step-selection analysis (iSSA; Avgar et al. 2016) with
# individual random slopes fitted by the Poisson trick in glmmTMB (Muff et al.
# 2020), relative selection strength (Avgar et al. 2017), movement-kernel
# updates (Fieberg et al. 2021), a use-availability RSF for comparison and home
# ranges (AKDE via ctmm; KDE and MCP via amt).
#
# Data structure (what simulate()$data and read_data() both return):
#   list(fixes      = data.frame(id, t = POSIXct UTC, x, y) in a projected CRS (m),
#        land       = list(xmin, ymax, resx, resy, nrow, ncol, crs,
#                          layers  = named list of nrow x ncol matrices (row 1 = north),
#                          scaling = list(<cov> = list(center, scale))),
#        covariates = character, names of the layers used in the model,
#        labels, units, unit_text = per-covariate axis label, meaningful
#                     change and its wording, interval_hours (NULL = auto),
#        animal_label = e.g. "elk")

for (f in c("simulate.R", "engine.R", "output.R")) sys.source(file.path(MODULE_DIR, f), envir = environment())

SSF_SETTINGS <- list(method = "issa", n_control = 20, interval_hours = NULL, tolerance_minutes = NULL,
                     random_slopes = TRUE, reml = TRUE, home_attraction = TRUE, home_centre = "estimated", min_steps = 20,
                     rsf = TRUE, rsf_n_avail = 10, rsf_weight = 5000,
                     home_range = TRUE, hr_level = 0.95, uhc = TRUE, uhc_draws = 200, seed = 1,
                     cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

ssf_fit <- function(data, settings) {
  st <- merge_list(SSF_SETTINGS, settings)
  notes <- character()
  if (identical(st$method, "rsf")) {
    # Negative-control route: RSF coefficients reported under the SSF names.
    R <- ssf_fit_rsf(data, st, prefix = "sel_")
    return(list(estimates = R$estimates,
                diagnostics = list(converged = R$converged, notes = "RSF fitted in place of the iSSA (method: rsf)."),
                object = R$model, settings = st))
  }
  S <- ssf_build_steps(data, st)
  F <- ssf_fit_issa(S, data, st)
  est <- F$estimates
  notes <- c(notes, data$intake_notes, S$notes, F$notes)
  n_an <- length(unique(S$steps$id))
  notes <- c(notes, sprintf("%d animals, %d used steps at a %s-h interval (tolerance %d min), %d available steps per used step; tentative step-length gamma(shape %.2f, scale %.0f m), turn-angle von Mises(kappa %.2f).",
                            n_an, sum(S$steps$case), format(S$interval_hours), S$tolerance_minutes, st$n_control,
                            S$sl_distr$params$shape, S$sl_distr$params$scale, S$ta_distr$params$kappa))
  if (F$random_slopes && n_an < 6) notes <- c(notes, sprintf("Only %d animals: the among-animal SD of selection is poorly estimated (validation: biased low and under-covered with 5 animals).", n_an))
  R <- NULL; H <- NULL; U <- NULL
  if (isTRUE(st$rsf)) {
    R <- tryCatch(ssf_fit_rsf(data, st), error = function(e) { notes <<- c(notes, paste("RSF failed:", conditionMessage(e))); NULL })
    if (!is.null(R)) est <- rbind(est, R$estimates)
  }
  if (isTRUE(st$home_range)) {
    H <- tryCatch(ssf_home_ranges(data, st), error = function(e) { notes <<- c(notes, paste("Home ranges failed:", conditionMessage(e))); NULL })
    if (!is.null(H)) {
      ht <- H$table
      ok <- is.finite(ht$akde_km2)
      if (any(ok)) est <- rbind(est, data.frame(param = paste0("hr_akde_km2_", ht$id[ok]), estimate = ht$akde_km2[ok], sd = NA,
                                                lower = ht$akde_lower[ok], upper = ht$akde_upper[ok]))
      if (any(!ok)) notes <- c(notes, ht$note[!ok])
    }
  }
  if (isTRUE(st$uhc)) {
    U <- tryCatch(ssf_uhc(S, data, st, F), error = function(e) { notes <<- c(notes, paste("UHC check failed:", conditionMessage(e))); NULL })
    if (!is.null(U)) {
      bad <- U$table$covariate[!U$table$ok]
      notes <- c(notes, if (length(bad)) sprintf("Used-habitat calibration (held-out 20%% of steps): observed use of %s lies OUTSIDE the 95%% predictive interval; the model misrepresents selection for it (missing covariate, non-linearity, or individual differences).", paste(bad, collapse = ", "))
                        else "Used-habitat calibration (held-out 20% of steps): observed mean of every covariate at used steps lies within the model's 95% predictive interval.")
    }
  }
  rownames(est) <- NULL
  keep <- S$steps[, c("id", "stratum", "case", "sl_", "ta_", "x1_", "y1_", "x2_", "y2_", paste0("z_", data$covariates))]
  list(estimates = est,
       diagnostics = list(converged = F$converged, notes = notes, model_checks_ok = if (!is.null(U)) all(U$table$ok) else NA),
       object = F$model, issa = F[c("coef", "vcov", "indiv", "random_slopes", "formula")],
       steps = keep, sl_distr = S$sl_distr$params, ta_distr = S$ta_distr$params,
       interval_hours = S$interval_hours, tolerance_minutes = S$tolerance_minutes, rsf = if (!is.null(R)) R[c("estimates", "n_used", "n_avail")],
       home_ranges = H, uhc = U, settings = st)
}

ssf_read_data <- function(config) {
  eco_require(c("terra", "sf"))
  dc <- config$data
  # YAML 1.1 reads an unquoted key `y:` (and `n:`, `yes:`...) as the boolean TRUE/FALSE
  if (is.null(dc$y) && !is.null(dc[["TRUE"]])) dc$y <- dc[["TRUE"]]
  rp <- function(f) if (grepl("^/", f)) f else file.path(config$.dir, f)
  raw <- utils::read.csv(rp(dc$fixes), stringsAsFactors = FALSE, check.names = FALSE)
  col <- function(k, def) { n <- dc[[k]] %||% def; if (!n %in% names(raw)) stop(sprintf("column '%s' (data: %s) not found in %s; columns are: %s", n, k, dc$fixes, paste(names(raw), collapse = ", "))); raw[[n]] }
  tt <- col("time", "timestamp")
  fmt <- dc$time_format %||% c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%d %H:%M", "%m/%d/%Y %H:%M")
  t <- as.POSIXct(tt, tz = dc$timezone %||% "UTC", tryFormats = fmt, optional = TRUE)
  fixes <- data.frame(id = as.character(col("id", "animal_id")), t = t,
                      x = suppressWarnings(as.numeric(col("x", "x"))), y = suppressWarnings(as.numeric(col("y", "y"))),
                      stringsAsFactors = FALSE)
  # coordinates that are present but not numbers (e.g. "511,014.9") are not failed fixes
  bad_xy <- character(); bad_row <- rep(FALSE, nrow(fixes))
  for (k in c("x", "y")) { r <- as.character(col(k, k)); b <- !is.na(r) & nzchar(trimws(r)) & is.na(fixes[[k]]); bad_xy <- c(bad_xy, r[b]); bad_row <- bad_row | b }
  n_bad_time <- sum(is.na(fixes$t) & !is.na(tt) & nzchar(as.character(tt)))
  no_id <- is.na(fixes$id) | !nzchar(trimws(fixes$id))
  n_no_id <- sum(no_id)
  fixes <- fixes[!no_id, ]; bad_row <- bad_row[!no_id]
  n_raw <- nrow(fixes) + n_no_id; n_na <- sum((is.na(fixes$x) | is.na(fixes$y)) & !bad_row)
  intake <- sprintf("%d GPS rows read; %d without coordinates (failed fixes) were dropped.", n_raw, n_na)
  if (n_bad_time > 0.5 * n_raw) stop(sprintf("Most timestamps (%d of %d) could not be read, e.g. '%s'. Set data: time_format (R strptime format, e.g. '%%d/%%m/%%Y %%H:%%M') and timezone.", n_bad_time, n_raw, tt[1]))
  if (any(bad_row)) intake <- paste(intake, sprintf("%d rows had coordinates that are not numbers (e.g. '%s') and were dropped.", sum(bad_row), bad_xy[1]))
  if (n_no_id) intake <- paste(intake, sprintf("%d rows had no animal id and were dropped.", n_no_id))
  fixes <- fixes[stats::complete.cases(fixes), ]
  attr(fixes, "n_bad_time") <- n_bad_time
  attr(fixes, "n_no_id") <- n_no_id
  attr(fixes, "bad_xy") <- bad_xy
  attr(fixes, "n_failed") <- n_na
  rasters <- dc$rasters
  if (is.null(rasters) || !length(rasters)) stop("data: rasters: give one GeoTIFF per covariate, e.g. rasters: {forest: forest.tif}")
  rl <- lapply(rasters, function(f) terra::rast(rp(f)))
  r1 <- rl[[1]]
  crs <- if (!is.null(dc$crs)) dc$crs else paste0("EPSG:", terra::crs(r1, describe = TRUE)$code)
  if (isTRUE(dc$lonlat)) {
    ok <- stats::complete.cases(fixes[, c("x", "y")])
    p <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(fixes[ok, ], coords = c("x", "y"), crs = 4326), terra::crs(r1)))
    fixes$x[ok] <- p[, 1]; fixes$y[ok] <- p[, 2]
  }
  covs <- unlist(dc$covariates %||% names(rasters))
  # Standardise by the landscape-wide mean and SD of each raster (before cropping)
  scaling <- lapply(rl[covs], function(r) {
    g <- terra::global(r, c("mean", "sd"), na.rm = TRUE)
    list(center = g[1, 1], scale = g[1, 2])
  })
  # crop to the data extent + buffer so lookups are cheap
  ok <- stats::complete.cases(fixes[, c("x", "y")])
  buf <- dc$buffer_m %||% 5000
  e <- terra::ext(min(fixes$x[ok]) - buf, max(fixes$x[ok]) + buf, min(fixes$y[ok]) - buf, max(fixes$y[ok]) + buf)
  layers <- list(); geo <- NULL
  for (v in covs) {
    r <- rl[[v]]
    if (!terra::compareGeom(r, r1, stopOnError = FALSE)) r <- terra::resample(r, r1)
    rc <- terra::crop(r, e)
    layers[[v]] <- terra::as.matrix(rc, wide = TRUE)
    if (is.null(geo)) geo <- rc
  }
  land <- list(xmin = terra::xmin(geo), ymax = terra::ymax(geo), resx = terra::res(geo)[1], resy = terra::res(geo)[2],
               nrow = terra::nrow(geo), ncol = terra::ncol(geo), layers = layers, crs = crs, scaling = scaling,
               lonlat_raster = terra::is.lonlat(r1))
  labs <- dc$labels %||% list(); un <- dc$units %||% list(); ut <- dc$unit_text %||% list()
  for (v in covs) {
    if (is.null(labs[[v]])) labs[[v]] <- v
    if (is.null(ut[[v]]) && !is.null(un[[v]])) ut[[v]] <- paste(un[[v]], "units of", v)
  }
  nms <- dc$names %||% list(); for (v in covs) if (is.null(nms[[v]])) nms[[v]] <- v
  list(fixes = fixes, land = land, covariates = covs, labels = labs, units = un, unit_text = ut, names = nms,
       animal_plural = dc$animal_plural, intake_notes = intake, n_failed_fixes = n_na, n_fix_attempts = n_raw - n_no_id,
       interval_hours = dc$interval_hours, animal_label = dc$animal_label %||% "animal")
}

ssf_check_data <- function(data) {
  pr <- character()
  fx <- data$fixes
  nb <- attr(fx, "n_bad_time") %||% sum(is.na(fx$t))
  if (nb) pr <- c(pr, sprintf("%d timestamps could not be read; set data: time_format (e.g. '%%d/%%m/%%Y %%H:%%M') and timezone.", nb))
  if (isTRUE(data$land$lonlat_raster)) pr <- c(pr, "The covariate rasters are in longitude/latitude. Reproject them to a projected CRS in metres (e.g. UTM) first; step lengths and home-range areas need metres.")
  ok <- stats::complete.cases(fx)
  f <- fx[ok, ]
  if (nrow(f) && all(abs(f$x) <= 180) && all(abs(f$y) <= 90)) pr <- c(pr, "Coordinates look like longitude/latitude; set data: lonlat: true so they are projected to the raster CRS.")
  if (isTRUE(attr(fx, "n_no_id") > 0)) pr <- c(pr, sprintf("%d rows have an empty animal id; they were dropped. Fill in the id if those fixes belong to a collared animal.", attr(fx, "n_no_id")))
  bx <- attr(fx, "bad_xy")
  if (length(bx)) pr <- c(pr, sprintf("%d coordinate values are not plain numbers (e.g. '%s'); those rows were dropped. Remove thousands separators/units or fix the decimal mark in the CSV.", length(bx), bx[1]))
  # ids that differ only by case, spaces, punctuation or leading zeros ("E05" vs "e05 ", "E3" vs "E03")
  uid <- unique(f$id)
  key <- gsub("(?<![0-9])0+(?=[0-9])", "", toupper(gsub("[^A-Za-z0-9]", "", uid)), perl = TRUE)
  for (k in unique(key[duplicated(key)])) {
    same <- uid[key == k]
    pr <- c(pr, sprintf("Animal ids %s look like the same animal written differently (case, spaces or leading zeros: %s fixes). Correct the id column so each animal has one spelling, or rename them if they really are different animals.",
                        paste0("'", same, "'", collapse = " and "), paste(table(f$id)[same], collapse = " and ")))
  }
  dup <- duplicated(f[, c("id", "t")])
  if (any(dup)) pr <- c(pr, sprintf("%d duplicated (animal, timestamp) fixes, e.g. %s at %s; remove duplicates (keep the better fix).", sum(dup), f$id[dup][1], format(f$t[dup][1])))
  n <- table(f$id)
  if (any(n < 30)) pr <- c(pr, sprintf("Animals with fewer than 30 fixes: %s. Drop them or merge deployments.", paste(names(n)[n < 30], collapse = ", ")))
  if (length(n) < 2) pr <- c(pr, "Only one animal: population-level selection and among-animal variation need several animals (set settings: random_slopes: false to describe this one animal).")
  for (v in data$covariates) {
    z <- ssf_lookup(data$land, f$x, f$y, v)
    out <- mean(is.na(z))
    if (out > 0.05) {
      by_id <- tapply(is.na(z), f$id, mean); worst <- names(by_id)[by_id > 0.5]
      pr <- c(pr, sprintf("%.0f%% of fixes fall outside (or on NA cells of) the '%s' raster%s. Check that the GPS coordinates and the raster use the same CRS (%s).", 100 * out, v,
                          if (length(worst)) sprintf(" (mostly animal(s) %s, whose fixes are largely off the raster: wrong coordinates or CRS for these animals?)", paste(worst, collapse = ", ")) else "",
                          data$land$crs))
    }
    if (isTRUE(data$land$scaling[[v]]$scale == 0)) pr <- c(pr, sprintf("Raster '%s' is constant.", v))
  }
  pr
}

module <- ecomodule(
  id = "movement-ssf",
  title = "Habitat selection from telemetry: integrated step selection (iSSA), RSF and home ranges",
  description = "Where do GPS-collared animals choose to go, given where they can reach in one step? Estimates habitat-selection coefficients and relative selection strength (RSS) for landscape covariates, how much individuals differ, the selection-free movement kernel (step length and turning), a resource selection function for comparison, and home-range areas (AKDE, KDE, MCP).",
  defaults = list(
    params = list(sel_mean = list(forest = 0.8, elevation = -0.5), sel_sd = list(forest = 0.4, elevation = 0.3),
                  step_shape = 2, step_scale = 120, turn_kappa = 1, attraction = 0.25),
    design = list(n_animals = 12, n_steps = 400, interval_hours = 2, fix_rate = 0.95, time_jitter_s = 120,
                  covariates = c("forest", "elevation"), land_cells = 320, cell_m = 50, field_range_m = 400,
                  centre_margin_m = 5000, n_candidates = 200, burnin_steps = 100,
                  x0 = 500000, y0 = 5000000, crs = "EPSG:32611", start_time = "2024-05-01 00:00:00",
                  id_prefix = "E", animal_label = "elk")),
  simulate = ssf_simulate, fit = ssf_fit, read_data = ssf_read_data, check_data = ssf_check_data,
  interpret = ssf_interpret, plots = ssf_plots, figures = ssf_figures, methods_text = ssf_methods_text,
  write_data = ssf_write_data,
  settings = SSF_SETTINGS,
  fast_settings = merge_list(SSF_SETTINGS, list(rsf = FALSE, home_range = FALSE, uhc = FALSE)),
  power_params = c("sel_forest_mean", "sel_elevation_mean"),
  engine = "amt (steps, available steps) + glmmTMB (Poisson-trick mixed conditional logit); ctmm (AKDE)",
  references = c(
    "Avgar T, Potts JR, Lewis MA, Boyce MS (2016) Integrated step selection analysis: bridging the gap between resource selection and animal movement. Methods in Ecology and Evolution 7:619-630. doi:10.1111/2041-210X.12528",
    "Avgar T, Lele SR, Keim JL, Boyce MS (2017) Relative selection strength: quantifying effect size in habitat- and step-selection inference. Ecology and Evolution 7:5322-5330. doi:10.1002/ece3.3122",
    "Fieberg J, Signer J, Smith B, Avgar T (2021) A 'How to' guide for interpreting parameters in habitat-selection analyses. Journal of Animal Ecology 90:1027-1043. doi:10.1111/1365-2656.13441",
    "Muff S, Signer J, Fieberg J (2020) Accounting for individual-specific variation in habitat-selection studies: efficient estimation of mixed-effects models using Bayesian or frequentist computation. Journal of Animal Ecology 89:80-92. doi:10.1111/1365-2656.13087",
    "Signer J, Fieberg J, Avgar T (2019) Animal movement tools (amt): R package for managing tracking data and conducting habitat selection analyses. Ecology and Evolution 9:880-890. doi:10.1002/ece3.4823",
    "Signer J, Fieberg J, Avgar T (2017) Estimating utilization distributions from fitted step-selection functions. Ecosphere 8:e01771. doi:10.1002/ecs2.1771",
    "Fieberg JR, Forester JD, Street GM, Johnson DH, ArchMiller AA, Matthiopoulos J (2018) Used-habitat calibration plots: a new procedure for validating species distribution, resource selection, and step-selection models. Ecography 41:737-752. doi:10.1111/ecog.03123",
    "Fleming CH, Fagan WF, Mueller T, Olson KA, Leimgruber P, Calabrese JM (2015) Rigorous home range estimation with movement data: a new autocorrelated kernel density estimator. Ecology 96:1182-1188. doi:10.1890/14-2010.1",
    "Calabrese JM, Fleming CH, Gurarie E (2016) ctmm: an R package for analyzing animal relocation data as a continuous-time stochastic process. Methods in Ecology and Evolution 7:1124-1132. doi:10.1111/2041-210X.12559",
    "Brooks ME, Kristensen K, van Benthem KJ, et al. (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9:378-400.",
    "Fortin D, Beyer HL, Boyce MS, Smith DW, Duchesne T, Mao JS (2005) Wolves influence elk movements: behavior shapes a trophic cascade in Yellowstone National Park. Ecology 86:1320-1330."))
