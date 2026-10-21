# ecotools module `glm`: generalised linear (mixed) models made easy.
#
# Data structure (what simulate()$data and read_data() both return):
#   list(df   = data.frame with the response, predictors, grouping columns,
#               optional effort (offset) and trials columns,
#        spec = list(formula, family, link, offset, offset_per, trials, scale,
#                    units, unit_text, labels, response_label, zi, reference))
# See README.md for every key. Engine: glmmTMB (Laplace ML); stats::glm /
# MASS::glm.nb as a fallback for fixed-effects-only models.

for (f in c("engine.R", "simulate.R", "output.R", "figures.R")) sys.source(file.path(MODULE_DIR, f), envir = environment())

GLM_DEFAULT_SETTINGS <- list(engine = "glmmTMB", REML = TRUE, diagnostics = TRUE, n_sim_resid = 250, n_draws = 2000,
                             seed = 1, force_family = NULL, drop_random = FALSE)

glm_fit <- function(data, settings) {
  eco_require(c("glmmTMB", "reformulas", "MASS"))
  st <- merge_list(GLM_DEFAULT_SETTINGS, settings)
  prep <- glm_prepare(data, st)
  can_glm <- !length(prep$groups) && prep$fam %in% c("gaussian", "binomial", "poisson", "nbinom2", "Gamma")
  notes <- character()
  if (identical(st$engine, "glm")) {
    if (!can_glm) stop("engine 'glm' only handles models without random effects and families gaussian, binomial, poisson, nbinom2, Gamma")
    fo <- glm_fit_glm(prep)
  } else {
    fo <- tryCatch(glm_fit_tmb(prep, REML = st$REML), error = function(e) e)
    if (inherits(fo, "error")) {
      if (!can_glm) stop("glmmTMB failed: ", conditionMessage(fo))
      notes <- c(notes, sprintf("glmmTMB failed (%s); the model was refitted with stats::glm / MASS::glm.nb.", conditionMessage(fo)))
      fo <- glm_fit_glm(prep)
    }
  }
  est <- glm_estimates(fo, prep, st)
  notes <- c(notes, glm_fit_notes(fo, prep, est))
  dh <- NULL; r2 <- NULL
  if (isTRUE(st$diagnostics)) {
    dh <- glm_dharma(fo, prep, st)
    r2 <- glm_r2(fo, prep)
    notes <- c(notes, dharma_notes(dh, prep))
    if (length(r2)) notes <- c(notes, paste0("R2: ", paste(sprintf("%s = %s", names(r2), fmt(r2, 3)), collapse = "; "), "."))
  }
  checks_ok <- if (is.null(dh) || !is.null(dh$error)) NA else all(dh$p >= 0.05)
  list(estimates = est,
       diagnostics = list(converged = fo$converged, notes = notes, engine = fo$engine, r2 = r2,
                          dharma_p = if (!is.null(dh)) dh$p, model_checks_ok = checks_ok,
                          logLik = fo$logLik, AIC = tryCatch(stats::AIC(fo$model), error = function(e) NA_real_)),
       object = fo$model, prep = prep, fo = fo, dharma = dh, settings = st)
}

glm_fit_notes <- function(fo, prep, est) {
  n <- character()
  n <- c(n, sprintf("Engine %s; %d observations%s.", fo$engine, fo$nobs,
                    if (length(prep$groups)) paste0(", ", paste(sprintf("%d levels of %s", vapply(prep$groups, function(g) nlevels(prep$df[[g]]), 1L), prep$groups), collapse = ", ")) else ""))
  if (!fo$converged) n <- c(n, sprintf("NOT CONVERGED: optimiser code/message '%s'; Hessian positive definite: %s. Simplify the model (fewer terms, drop a random effect) or check for predictors that perfectly separate the outcome.",
                                       fo$message %||% "", fo$pdHess))
  if (length(fo$warnings)) n <- c(n, paste("Engine warnings:", paste(fo$warnings, collapse = " | ")))
  for (g in prep$groups) {
    k <- nlevels(prep$df[[g]])
    if (k < 5) n <- c(n, sprintf("Only %d levels of %s: a random-intercept SD is poorly estimated with fewer than ~5 groups; consider %s as a fixed factor.", k, g, g))
    tb <- table(prep$df[[g]])
    if (stats::median(tb) >= 3 && any(tb == 1))
      n <- c(n, sprintf("%d level(s) of %s have a single observation while the typical %s has %d (%s). Check that these ids are not misspelled.",
                        sum(tb == 1), g, g, as.integer(stats::median(tb)), paste(utils::head(names(tb)[tb == 1], 6), collapse = ", ")))
    r <- est[est$param == paste0("re_sd_", clean_name(g)), ]
    if (nrow(r) && (r$estimate < 1e-3 || !is.finite(r$upper))) n <- c(n, sprintf("Random-intercept SD for %s is at the boundary (~0): singular fit. Fixed-effect estimates are still usable; the SD is not.", g))
  }
  if (prep$fam == "binomial" && prep$link == "logit") {
    big <- est$param[grepl("^beta_", est$param) & abs(est$estimate) > 8 & est$param != "beta_intercept"]
    if (length(big)) n <- c(n, paste("Very large logit coefficients (", paste(big, collapse = ", "), "): possible complete separation (a predictor perfectly predicts 0/1); intervals are unreliable."))
  }
  if (length(prep$num_vars) > 1) {
    cm <- suppressWarnings(stats::cor(prep$raw[prep$num_vars]))
    hi <- which(abs(cm) > 0.7 & upper.tri(cm), arr.ind = TRUE)
    if (nrow(hi)) n <- c(n, paste0("Correlated predictors (|r| > 0.7): ", paste(sprintf("%s & %s (r = %s)", rownames(cm)[hi[, 1]], colnames(cm)[hi[, 2]], fmt(cm[hi], 2)), collapse = "; "),
                                   ". Their separate effects are hard to tell apart; intervals will be wide and each effect is 'holding the other constant'."))
  }
  npar <- length(fo$beta)
  if (fo$nobs / npar < 10) n <- c(n, sprintf("Only %.1f observations per fixed-effect coefficient: estimates may be unstable; consider fewer predictors.", fo$nobs / npar))
  n
}

# ---- real data ------------------------------------------------------------------
glm_read_data <- function(config) {
  dc <- config$data
  if (is.null(dc$file)) stop("config needs data: file: <your csv>")
  path <- dc$file
  if (!grepl("^/", path)) path <- file.path(config$.dir %||% ".", path)
  if (!file.exists(path)) stop("data file not found: ", path)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE, strip.white = TRUE,
                        na.strings = c("", "NA", "na", "N/A", "n/a", "."))
  spec <- list(formula = dc$formula, family = dc$family, link = dc$link, offset = dc[["offset"]],
               offset_per = dc$offset_per %||% 1, trials = dc$trials, scale = dc$scale %||% "all",
               units = dc$units, unit_text = dc$unit_text, labels = dc$labels,
               response_label = dc$response_label, zi = dc$zi, reference = dc$reference, source = path)
  if (isTRUE(dc$drop_missing)) {
    used <- tryCatch(glm_used_columns(spec), error = function(e) character())
    used <- intersect(used, names(df))
    keep <- stats::complete.cases(df[used])
    spec$dropped_rows <- which(!keep)
    df <- df[keep, , drop = FALSE]
  }
  list(df = df, spec = spec)
}

glm_used_columns <- function(spec) {
  pf <- parse_glm_formula(spec$formula)
  unique(c(pf$response, pf$vars, pf$groups, spec[["offset"]], spec$trials))
}

glm_check_data <- function(data) {
  pr <- character()
  sp <- data$spec; df <- data$df
  fam <- tryCatch(canon_family(sp$family), error = function(e) { pr <<- c(pr, conditionMessage(e)); NULL })
  pf <- tryCatch(parse_glm_formula(sp$formula), error = function(e) { pr <<- c(pr, conditionMessage(e)); NULL })
  if (is.null(fam) || is.null(pf)) return(pr)
  if (!is.null(sp$link) && !sp$link %in% c("identity", "log", "logit", "probit", "cloglog", "inverse", "sqrt"))
    pr <- c(pr, sprintf("link '%s' not recognised", sp$link))
  cols <- unique(c(pf$response, pf$vars, pf$groups, sp[["offset"]], sp$trials))
  miss <- setdiff(cols, names(df))
  if (length(miss)) return(c(pr, sprintf("column(s) %s not found in the data file. Columns present: %s. (Column names are made syntactic by R: spaces and symbols become dots.)",
                                         paste(miss, collapse = ", "), paste(names(df), collapse = ", "))))
  for (v in cols) {
    k <- which(is.na(df[[v]]))
    if (length(k)) pr <- c(pr, sprintf("column %s has %d missing value(s) (data rows %s%s, not counting the header line). Fill them, or set data: drop_missing: true to drop those rows.",
                                       v, length(k), paste(utils::head(k, 8), collapse = ", "), if (length(k) > 8) ", ..." else ""))
  }
  if (length(pr)) return(pr)
  y <- df[[pf$response]]
  if (!is.numeric(y)) return(c(pr, sprintf("response %s must be numeric; it contains text such as '%s'", pf$response, utils::head(y[is.na(suppressWarnings(as.numeric(y)))], 1))))
  isint <- function(x) all(abs(x - round(x)) < 1e-8)
  if (fam == "binomial") {
    if (is.null(sp$trials)) {
      if (!all(y %in% c(0, 1))) pr <- c(pr, sprintf("binomial response %s must be 0/1. If it is a number of successes, give the number of trials under data: trials: <column>.", pf$response))
    } else {
      tr <- df[[sp$trials]]
      if (any(y < 0) || !isint(y) || !isint(tr)) pr <- c(pr, "successes and trials must be whole numbers >= 0")
      bad <- which(y > tr)
      if (length(bad)) pr <- c(pr, sprintf("successes exceed trials in row(s) %s", paste(utils::head(bad, 8), collapse = ", ")))
      if (any(tr == 0)) pr <- c(pr, sprintf("%d row(s) have 0 trials; remove them", sum(tr == 0)))
    }
  }
  if (fam %in% COUNT_FAMILIES) {
    if (any(y < 0)) pr <- c(pr, sprintf("count response %s has negative values", pf$response))
    if (!isint(y)) pr <- c(pr, sprintf("count response %s has non-integer values (e.g. %s). Counts must be whole numbers; for densities use the count with an offset for area/effort, or a Gamma/gaussian family for continuous data.",
                                       pf$response, utils::head(y[abs(y - round(y)) > 1e-8], 1)))
  }
  if (fam == "beta" && any(y <= 0 | y >= 1))
    pr <- c(pr, sprintf("beta response %s must be strictly between 0 and 1 (found %d values at 0 or 1, or outside). Percentages must be divided by 100. If exact 0s/1s occur, the common fix is y' = (y*(n-1) + 0.5)/n (Smithson & Verkuilen 2006), or use binomial with successes/trials if the proportion comes from counts.",
                                    pf$response, sum(y <= 0 | y >= 1)))
  if (fam == "Gamma" && any(y <= 0)) pr <- c(pr, sprintf("Gamma response %s must be > 0 (found %d values <= 0)", pf$response, sum(y <= 0)))
  if (!is.null(sp[["offset"]])) {
    if (any(df[[sp[["offset"]]]] <= 0)) pr <- c(pr, sprintf("effort/offset column %s must be > 0 (log is taken)", sp[["offset"]]))
    lk <- sp$link %||% default_link(fam)
    if (lk != "log") pr <- c(pr, sprintf("an effort offset needs a log link; this model uses '%s'. For binomial data with unequal effort, include effort as a predictor (or use cloglog link) instead.", lk))
  } else if (!is.null(sp$offset_per) && !isTRUE(all.equal(as.numeric(sp$offset_per), 1)))
    pr <- c(pr, sprintf("offset_per: %s is set but no effort column is given under offset:. Add offset: <effort column>, or remove offset_per.", sp$offset_per))
  for (g in pf$groups) if (length(unique(df[[g]])) < 2) pr <- c(pr, sprintf("grouping column %s has a single value; remove (1|%s)", g, g))
  # likely typos in group ids: labels that are identical once case, punctuation and leading zeros are ignored
  for (g in pf$groups) {
    lv <- unique(as.character(df[[g]]))
    key <- gsub("(?<![0-9])0+(?=[0-9])", "", tolower(gsub("[^A-Za-z0-9]", "", lv)), perl = TRUE)
    dup <- split(lv, key); dup <- dup[lengths(dup) > 1]
    if (length(dup)) pr <- c(pr, sprintf("grouping column %s has labels that look like the same %s spelled differently: %s. If they are the same %s, correct the spelling in the data file; each spelling is otherwise treated as a separate %s.",
                                         g, g, paste(vapply(dup, function(d) paste(sprintf("'%s' (%d rows)", d, as.integer(table(df[[g]])[d])), collapse = " vs "), ""), collapse = "; "), g, g))
  }
  for (v in pf$vars) {
    x <- df[[v]]
    if (is.character(x) || is.logical(x) || is.factor(x)) {
      tb <- table(x)
      if (length(tb) < 2) pr <- c(pr, sprintf("factor %s has only one level (%s)", v, names(tb)))
      if (length(tb) > 30) pr <- c(pr, sprintf("factor %s has %d levels: if it identifies sites/individuals, use it as a random intercept (1|%s) instead", v, length(tb), v))
      if (!is.null(sp$reference[[v]]) && !sp$reference[[v]] %in% names(tb))
        pr <- c(pr, sprintf("reference level '%s' for %s not found (levels: %s)", sp$reference[[v]], v, paste(names(tb), collapse = ", ")))
    } else if (is.numeric(x) && stats::sd(x) == 0) pr <- c(pr, sprintf("predictor %s is constant", v))
  }
  for (v in names(sp$units)) if (!v %in% pf$vars) pr <- c(pr, sprintf("units given for %s, which is not a predictor in the formula", v))
  if (nrow(df) < 10) pr <- c(pr, sprintf("only %d rows of data", nrow(df)))
  pr
}

module <- ecomodule(
  id = "glm",
  title = "GLM / GLMM with biological interpretation",
  description = paste(
    "How does a response (a count, presence/absence, a proportion, a size or mass, a continuous measurement) change with",
    "habitat, environmental or treatment predictors, allowing for repeated measures within sites or individuals?",
    "Fits generalised linear (mixed) models for Gaussian, binomial, Poisson, negative-binomial, beta, Gamma and",
    "zero-inflated count data, and reports effects as odds ratios, rate ratios and changes on the response scale per SD",
    "and per biologically meaningful unit, with conditional-effects plots and residual checks."),
  defaults = list(
    params = list(beta_intercept = 0.7, beta_forest = 0.3, beta_temp = -0.2,
                  beta_habitat_riparian = 0.4, beta_habitat_wetland = -0.3,
                  re_sd = 0.5, dispersion = 2, zi_prob = 0.2),
    design = list(family = "nbinom2", n_groups = 30, n_per_group = 8, group_name = "site",
                  covariates = list(
                    forest = list(type = "normal", mean = 55, sd = 20, level = "group"),
                    temp = list(type = "normal", mean = 15, sd = 5),
                    habitat = list(type = "factor", levels = c("upland", "riparian", "wetland"), level = "group")),
                  interactions = NULL,
                  offset = list(name = "effort", min = 0.5, max = 2),
                  units = list(forest = 10, temp = 5), response = "count")),
  simulate = glm_simulate, fit = glm_fit,
  settings = GLM_DEFAULT_SETTINGS,
  fast_settings = merge_list(GLM_DEFAULT_SETTINGS, list(diagnostics = FALSE)),
  read_data = glm_read_data, check_data = glm_check_data,
  interpret = glm_interpret, plots = glm_plots, figures = glm_figures, methods_text = glm_methods,
  power_params = c("beta_forest", "beta_temp", "beta_habitat_riparian", "beta_habitat_wetland"),
  engine = "glmmTMB (REML for mixed models, ML otherwise; Laplace approximation); stats::glm / MASS::glm.nb fallback",
  references = c(
    "Brooks ME, Kristensen K, van Benthem KJ, Magnusson A, Berg CW, Nielsen A, Skaug HJ, Maechler M, Bolker BM (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9(2):378-400. doi:10.32614/RJ-2017-066",
    "Bolker BM, Brooks ME, Clark CJ, Geange SW, Poulsen JR, Stevens MHH, White J-SS (2009). Generalized linear mixed models: a practical guide for ecology and evolution. Trends in Ecology and Evolution 24(3):127-135.",
    "McCullagh P, Nelder JA (1989). Generalized Linear Models, 2nd ed. Chapman & Hall.",
    "Hartig F (2026). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package 0.5.0. doi:10.32614/CRAN.package.DHARMa",
    "Dunn PK, Smyth GK (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5(3):236-244.",
    "Nakagawa S, Johnson PCD, Schielzeth H (2017). The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. Journal of the Royal Society Interface 14:20170213. doi:10.1098/rsif.2017.0213",
    "Luedecke D, Ben-Shachar MS, Patil I, Waggoner P, Makowski D (2021). performance: an R package for assessment, comparison and testing of statistical models. Journal of Open Source Software 6(60):3139.",
    "Schielzeth H (2010). Simple means to improve the interpretability of regression coefficients. Methods in Ecology and Evolution 1(2):103-113. doi:10.1111/j.2041-210X.2010.00012.x",
    "Ferrari SLP, Cribari-Neto F (2004). Beta regression for modelling rates and proportions. Journal of Applied Statistics 31(7):799-815.",
    "Hardin JW, Hilbe JM (2007). Generalized Linear Models and Extensions, 2nd ed. Stata Press. (NB1/NB2 parameterisations)")
)
