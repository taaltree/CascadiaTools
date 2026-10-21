# fisheries-growth: von Bertalanffy growth (nls + Kimura LRT, or Bayesian
# fixed/hierarchical groups in Nimble with group covariates on Linf and K),
# length-weight (log-log with bias correction) and condition.
#
# Data structure (simulate()$data and read_data() both return this):
#   list(df = data.frame(fish_id, age, length, group (factor), weight (NA ok)),
#        group_cov = data.frame(group, <cov>) or NULL, cov_name = NULL or name,
#        ref_group, reference_length, units = list(length, weight, age),
#        ws = NULL or list(log10_a, b, min_length, source), species)

for (f in c("simulate.R", "engine_nls.R", "engine_nimble.R", "lw.R", "output.R"))
  source(file.path(MODULE_DIR, f), local = TRUE)

fg_fit <- function(data, settings) {
  eco_require(c("FSA"))
  settings <- merge_list(fg_default_settings, settings)
  df <- data$df
  df$group <- droplevels(factor(df$group))
  if (!is.null(data$ref_group) && data$ref_group %in% levels(df$group))
    df$group <- stats::relevel(df$group, ref = as.character(data$ref_group))
  data$df <- df
  gf <- if (settings$engine == "nimble") fit_nimble_fg(data, settings) else fit_nls_fg(data, settings)
  est <- gf$estimates
  for (k in c("rhat", "ess_bulk")) if (!k %in% names(est)) est[[k]] <- NA_real_
  lw <- if (any(is.finite(df$weight))) fit_lw(data) else NULL
  if (!is.null(lw)) est <- rbind(est, cbind(lw$rows, rhat = NA_real_, ess_bulk = NA_real_))
  notes <- gf$diagnostics$notes
  # Is Linf supported by old fish? Compare the oldest fish to the age at 90% Linf.
  gn <- levels(df$group); sg <- fg_safe(gn)
  for (g in seq_along(gn)) {
    a90 <- est$estimate[est$param == paste0("age_at_90pct_Linf_", sg[g])]
    amax <- max(df$age[df$group == gn[g]]); nold <- sum(df$age[df$group == gn[g]] >= a90)
    if (length(a90) && is.finite(a90) && (amax < a90 || nold < 5))
      notes <- c(notes, sprintf("Group %s: only %d fish are at or beyond the estimated age at 90%% of Linf (%.1f yr; oldest fish %.0f yr). Linf is extrapolated beyond most of the data; treat Linf (and K, which is correlated with it) with caution.",
                                gn[g], nold, a90, amax))
  }
  if (settings$engine == "nimble") {
    if (gf$diagnostics$rhat_max > 1.01 || gf$diagnostics$ess_min < 400)
      notes <- c(notes, sprintf("Real-data standard not met (R-hat %.3f, min bulk ESS %.0f; need < 1.01 and >= 400): rerun with more iterations (settings: niter).",
                                gf$diagnostics$rhat_max, gf$diagnostics$ess_min))
  }
  if (!is.null(lw) && length(lw$outliers))
    notes <- c(notes, sprintf("%d fish have length-weight residuals > 4 SD (possible data-entry errors): %s",
                              length(lw$outliers), paste(utils::head(lw$outliers, 10), collapse = ", ")))
  dg <- gf$diagnostics; dg$notes <- notes
  list(estimates = est, diagnostics = dg,
       object = list(engine = settings$engine, settings = settings, growth = gf, lw = lw))
}

fg_default_settings <- list(engine = "nls", error = "multiplicative", group_effects = "random",
                            kimura = TRUE, niter = 20000, nburnin = 5000, nchains = 3, thin = 1,
                            seed = 1)

fg_read_data <- function(config) {
  dc <- config$data
  f <- dc$file; if (!grepl("^/", f)) f <- file.path(config$.dir, f)
  raw <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  col <- function(key, default = NULL, required = TRUE) {
    nm <- if (!is.null(dc[[key]])) dc[[key]] else default
    if (is.null(nm) || !nm %in% names(raw)) {
      if (required) stop(sprintf("column '%s' (config data:%s) not found in %s. Columns are: %s",
                                 nm, key, basename(f), paste(names(raw), collapse = ", ")))
      return(NULL)
    }
    raw[[nm]]
  }
  grp <- col("group", required = !is.null(dc$group))
  if (!is.null(dc$filter)) {  # e.g. filter: {year: 2014}
    keep <- rep(TRUE, nrow(raw))
    for (k in names(dc$filter)) {
      if (!k %in% names(raw))
        stop(sprintf("filter column '%s' (config data:filter) not found in %s. Columns are: %s",
                     k, basename(f), paste(names(raw), collapse = ", ")))
      kk <- raw[[k]] %in% unlist(dc$filter[[k]])
      if (!any(kk))
        stop(sprintf("filter %s: %s matches no rows of %s. Values present: %s",
                     k, paste(unlist(dc$filter[[k]]), collapse = ", "), basename(f),
                     paste(utils::head(sort(unique(raw[[k]])), 20), collapse = ", ")))
      keep <- keep & kk
    }
  } else keep <- rep(TRUE, nrow(raw))
  df <- data.frame(
    fish_id = { v <- col("fish_id", required = FALSE); if (is.null(v)) seq_len(nrow(raw)) else v },
    age = suppressWarnings(as.numeric(col("age", "age"))),
    length = suppressWarnings(as.numeric(col("length", "length"))),
    group = if (is.null(grp)) "all" else as.character(grp),
    weight = { v <- col("weight", required = FALSE); if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v)) },
    stringsAsFactors = FALSE)[keep, ]
  if (!is.null(dc$age_offset)) df$age <- df$age + as.numeric(dc$age_offset)
  lv <- if (!is.null(dc$group_levels)) unlist(dc$group_levels) else sort(unique(df$group))
  unlisted <- setdiff(unique(df$group), lv)  # values not in data:group_levels (typos?)
  df$group <- factor(df$group, levels = lv)
  gcov <- NULL; covname <- NULL
  if (!is.null(dc$group_covariates)) {
    gf <- dc$group_covariates$file; if (!grepl("^/", gf)) gf <- file.path(config$.dir, gf)
    gcov <- utils::read.csv(gf, stringsAsFactors = FALSE, check.names = FALSE)
    gcol <- if (is.null(dc$group_covariates$group)) dc$group else dc$group_covariates$group
    names(gcov)[names(gcov) == gcol] <- "group"
    gcov$group <- as.character(gcov$group)
    covname <- dc$group_covariates$covariate
  }
  ws <- NULL; sw <- dc$standard_weight
  if (!is.null(sw)) {
    if (!is.null(sw$species)) {
      w <- FSA::wsVal(sw$species, units = if (is.null(sw$units)) "metric" else sw$units,
                      group = sw$group)
      if (nrow(w) != 1) stop(sprintf("standard_weight: %d equations match species '%s'; add group: (one of %s)",
                                     nrow(w), sw$species, paste(unique(w$group), collapse = ", ")))
      ml <- w[[grep("^min\\.", names(w))[1]]]
      ws <- list(log10_a = w$int, b = w$slope, min_length = if (is.finite(ml)) ml else NULL,
                 source = sprintf("%s standard weight equation (%s, %s method; via FSA::wsVal)", sw$species, w$source, w$method))
    } else ws <- list(log10_a = sw$log10_a, b = sw$b, min_length = sw$min_length, source = sw$source)
  }
  lref <- if (!is.null(dc$reference_length)) dc$reference_length else
    signif(stats::median(df$length, na.rm = TRUE), 2)
  list(df = df, group_cov = gcov, cov_name = covname,
       ref_group = if (!is.null(dc$reference_group)) dc$reference_group else levels(df$group)[1],
       reference_length = lref,
       units = list(length = if (is.null(dc$length_units)) "mm" else dc$length_units,
                    weight = if (is.null(dc$weight_units)) "g" else dc$weight_units, age = "years"),
       ws = ws, species = if (!is.null(dc$species)) dc$species else sw$species,
       unlisted_groups = unlisted)
}

# Suggest the group a rare level was probably meant to be (case / small typos).
fg_typo_hint <- function(g, others) {
  others <- setdiff(others, g)
  if (!length(others)) return("")
  hit <- others[tolower(others) == tolower(g)]
  if (!length(hit)) hit <- others[utils::adist(tolower(g), tolower(others)) <= max(1, floor(nchar(g) / 4))]
  if (length(hit)) sprintf(" Is it a misspelling of '%s'? Correct it in the data file.", hit[1]) else ""
}

fg_check_data <- function(data) {
  df <- data$df; pr <- character()
  if (!nrow(df)) return("no rows left after reading/filtering the file")
  ids <- function(i) paste(utils::head(df$fish_id[i], 5), collapse = ", ")
  bad <- which(!is.finite(df$age)); if (length(bad)) pr <- c(pr, sprintf("%d fish have a missing or non-numeric age (fish %s%s). Remove them or fill them in.", length(bad), ids(bad), if (length(bad) > 5) ", ..." else ""))
  bad <- which(!is.finite(df$length)); if (length(bad)) pr <- c(pr, sprintf("%d fish have a missing or non-numeric length (fish %s%s). Remove them or fill them in (lengths must be plain numbers, e.g. 400 not '40O' or '400mm').", length(bad), ids(bad), if (length(bad) > 5) ", ..." else ""))
  if (any(df$length <= 0, na.rm = TRUE)) pr <- c(pr, sprintf("some lengths are zero or negative (fish %s)", ids(which(df$length <= 0))))
  if (any(df$age < 0, na.rm = TRUE)) pr <- c(pr, sprintf("some ages are negative (fish %s)", ids(which(df$age < 0))))
  if (any(df$age > 150, na.rm = TRUE)) pr <- c(pr, sprintf("ages up to %g: ages must be in years (are these months or days?).", max(df$age, na.rm = TRUE)))
  if (anyNA(df$group)) {
    ul <- data$unlisted_groups; ul <- ul[!is.na(ul) & nzchar(ul)]
    hints <- vapply(ul, function(u) { h <- fg_typo_hint(u, levels(df$group)); if (nzchar(h)) sub(".* of ('[^']*').*", " (did you mean \\1?)", h) else "" }, "")
    pr <- c(pr, sprintf("%d fish have a group value that is blank or not listed in data:group_levels (%s). Fix the spelling in the data file, fill in blank groups, or add the value to group_levels.",
                        sum(is.na(df$group)),
                        if (length(ul)) paste0("'", ul, "'", hints, collapse = ", ") else "blank group"))
  }
  if (any(duplicated(df$fish_id))) pr <- c(pr, sprintf("fish_id is duplicated (e.g. %s): each row should be one fish", df$fish_id[duplicated(df$fish_id)][1]))
  for (g in levels(df$group)) {
    d <- df[df$group == g & is.finite(df$age), ]
    if (nrow(d) < 15) pr <- c(pr, sprintf("group '%s' has only %d aged fish; growth curves need at least 15 (preferably 50+) per group.%s Otherwise merge groups or drop it.", g, nrow(d),
                                          fg_typo_hint(g, levels(df$group))))
    else if (length(unique(d$age)) < 4) pr <- c(pr, sprintf("group '%s' has only %d distinct ages; a 3-parameter growth curve needs at least 4.", g, length(unique(d$age))))
  }
  w <- df$weight[is.finite(df$weight)]
  if (length(w) && any(w <= 0)) pr <- c(pr, "some weights are zero or negative (leave missing weights blank)")
  if (length(w) > 10) {
    b <- stats::coef(stats::lm(log(df$weight) ~ log(df$length), subset = is.finite(df$weight) & df$weight > 0))[2]
    if (b < 1.5 || b > 4.5) pr <- c(pr, sprintf("length-weight slope b = %.2f is outside the biologically plausible 2-4 range: check that length and weight columns and units are right (Froese 2006).", b))
    # Gross weight errors (e.g. some weights entered in kg, or a misplaced decimal):
    # robust fit, then flag fish whose weight is > 2.5x off the length-weight line.
    ok <- which(is.finite(df$weight) & df$weight > 0 & is.finite(df$length) & df$length > 0)
    lL <- log(df$length[ok]); lW <- log(df$weight[ok])
    rb <- stats::median(lW - 3 * lL); r <- lW - 3 * lL - rb  # b = 3 start, then one refit on the inner 90%
    inn <- abs(r) < stats::quantile(abs(r), 0.9)
    cf <- stats::coef(stats::lm(lW[inn] ~ lL[inn])); r <- lW - cf[1] - cf[2] * lL
    far <- ok[abs(r) > log(2.5)]
    if (length(far)) pr <- c(pr, sprintf("%d fish weigh less than 1/2.5 or more than 2.5 times the weight expected for their length (fish %s%s): check for weights in the wrong units (kg vs g) or typing errors, then fix or blank them.",
                                         length(far), ids(far), if (length(far) > 5) ", ..." else ""))
  }
  if (!is.null(data$cov_name)) {
    if (is.null(data$group_cov) || !data$cov_name %in% names(data$group_cov))
      pr <- c(pr, sprintf("group covariate '%s' not found in the group covariate file", data$cov_name))
    else {
      miss <- setdiff(levels(df$group), data$group_cov$group)
      if (length(miss)) pr <- c(pr, sprintf("groups missing from the group covariate file: %s", paste(miss, collapse = ", ")))
    }
  }
  if (!is.null(data$ref_group) && !data$ref_group %in% levels(df$group))
    pr <- c(pr, sprintf("reference_group '%s' is not one of the groups (%s)", data$ref_group, paste(levels(df$group), collapse = ", ")))
  pr
}

module <- ecomodule(
  id = "fisheries-growth",
  title = "Fish growth (von Bertalanffy) and length-weight condition",
  description = paste(
    "How fast do fish grow and how large do they get, and do growth and condition differ among",
    "populations, sexes or years (or with a group covariate such as temperature)? Fits the von Bertalanffy",
    "growth curve to age-length data (nls with Kimura likelihood-ratio tests, or a Bayesian fixed- or",
    "random-effects model in Nimble), plus the length-weight relationship W = aL^b and condition."),
  defaults = list(
    params = list(Linf = c(600, 650, 600), K = c(0.30, 0.25, 0.30), t0 = -0.5,
                  sigma = 0.08, error = "multiplicative",
                  Linf_typical = 600, K_typical = 0.30, sd_log_Linf = 0.08, sd_log_K = 0.15,
                  Linf_beta = -0.05, K_beta = 0.10,
                  lw_a = 3.5e-6, lw_b = 3.18, lw_sigma = 0.08, condition = c(1, 1.05, 1)),
    design = list(n_groups = 3, n_per_group = 150, Z = 0.35, age_min = 1, max_age = 15,
                  age_offset = 0, prop_weighed = 1, reference_length = 500, min_length = 0,
                  hierarchical = FALSE, cov_name = "temperature", cov_mean = 20, cov_sd = 2)),
  simulate = simulate_fg, fit = fg_fit,
  settings = fg_default_settings,
  fast_settings = merge_list(fg_default_settings, list(niter = 3000, nburnin = 1000)),
  read_data = fg_read_data, check_data = fg_check_data,
  interpret = fg_interpret, plots = function(fit, data) fg_figures(fit, data, fig_opts()),
  figures = fg_figures, methods_text = fg_methods,
  power_params = c("Linf_diff_B_vs_A", "K_diff_B_vs_A", "t0_diff_B_vs_A",
                   "Linf_diff_C_vs_A", "K_diff_C_vs_A", "t0_diff_C_vs_A",
                   "condition_ratio_B_vs_A", "condition_ratio_C_vs_A",
                   "Linf_beta_temperature", "K_beta_temperature"),
  null_values = c(Linf_diff_B_vs_A = 0, K_diff_B_vs_A = 0, t0_diff_B_vs_A = 0,
                  Linf_diff_C_vs_A = 0, K_diff_C_vs_A = 0, t0_diff_C_vs_A = 0,
                  condition_ratio_B_vs_A = 1, condition_ratio_C_vs_A = 1,
                  Linf_beta_temperature = 0, K_beta_temperature = 0),
  engine = "stats::nls + FSA starts (Kimura LRT) or nimble (Bayesian fixed/hierarchical); lm for length-weight",
  references = c(
    "von Bertalanffy L (1938) A quantitative theory of organic growth. Human Biology 10:181-213.",
    "Beverton RJH, Holt SJ (1957) On the Dynamics of Exploited Fish Populations. Fishery Investigations Series II 19. HMSO, London.",
    "Kimura DK (1980) Likelihood methods for the von Bertalanffy growth curve. Fishery Bulletin 77:765-776.",
    "Helser TE, Lai HL (2004) A Bayesian hierarchical meta-analysis of fish growth: with an example for North American largemouth bass, Micropterus salmoides. Ecological Modelling 178:399-416.",
    "Ogle DH (2016) Introductory Fisheries Analyses with R. Chapman & Hall/CRC, Boca Raton.",
    "Ogle DH, Doll JC, Wheeler AP, Dinno A (2025) FSA: Simple Fisheries Stock Assessment Methods. R package version 0.10.x.",
    "Froese R (2006) Cube law, condition factor and weight-length relationships: history, meta-analysis and recommendations. Journal of Applied Ichthyology 22:241-253.",
    "Sprugel DG (1983) Correcting for bias in log-transformed allometric equations. Ecology 64:209-210.",
    "Le Cren ED (1951) The length-weight relationship and seasonal cycle in gonad weight and condition in the perch (Perca fluviatilis). Journal of Animal Ecology 20:201-219.",
    "Wege GJ, Anderson RO (1978) Relative weight (Wr): a new index of condition for largemouth bass. Pages 79-91 in Novinger GD, Dillard JG (eds) New Approaches to the Management of Small Impoundments. American Fisheries Society, North Central Division, Special Publication 5.",
    "Murphy BR, Brown ML, Springer TA (1990) Evaluation of the relative weight (Wr) index, with new applications to walleye. North American Journal of Fisheries Management 10:85-97.",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.")
)
