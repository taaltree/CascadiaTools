# Estimation for the causal module: GLM (any family) with the DAG's
# adjustment set, optional unit / unit+time fixed effects for panels,
# g-computation of the average causal effect with marginaleffects, the
# "what if" comparison under other adjustment strategies, and the DAG's
# testable implications (dagitty::localTests).

CAUSAL_DEFAULT_SETTINGS <- list(
  whatif = TRUE,             # estimate the effect under alternative adjustment strategies
  ci_method = "unconditional", # unconditional | delta | bootstrap
  n_boot = 500,
  vcov = "auto",             # auto (cluster-robust by unit for panels, sandwich HC3 for poisson, model-based otherwise) | model | robust
  local_tests = TRUE, ci_test_type = "auto",  # auto -> "cis" (numeric) or "cis.pillai" (factors)
  diagnostics = TRUE, n_sim_resid = 250,
  force_adjust = NULL,       # NEGATIVE CONTROLS ONLY: adjust for this set even if invalid
  force_panel = NULL,        # NEGATIVE CONTROLS ONLY: override the panel estimator (e.g. "none")
  seed = 1)

# ---- data preparation ----------------------------------------------------------
causal_prepare <- function(data, st) {
  sp <- data$spec; df <- data$df
  fam <- causal_family_info(sp[["family"]], sp[["link"]])
  id <- causal_identify(sp, names(df))
  x <- sp[["exposure"]]; y <- sp[["outcome"]]
  if (!is.null(st$force_panel)) id$panel_used <- st$force_panel else id$panel_used <- id$panel
  # exposure type
  xv <- df[[x]]
  if (is.logical(xv)) xv <- as.numeric(xv)
  if (is.character(xv) || is.factor(xv)) {
    lv <- if (is.factor(xv)) levels(xv) else sort(unique(xv))
    if (length(lv) != 2) stop("a categorical exposure must have exactly 2 levels (got ", length(lv), "); code it as several 0/1 exposures or pick one contrast")
    ref <- sp[["exposure_reference"]] %||% lv[1]
    if (!ref %in% lv) stop("exposure_reference '", ref, "' is not a level of ", x)
    sp[["exposure_levels"]] <- c(ref, setdiff(lv, ref))
    xv <- as.numeric(as.character(xv) != ref)
  }
  df[[x]] <- xv
  ux <- sort(unique(stats::na.omit(xv)))
  ex_type <- sp[["exposure_type"]] %||% if (length(ux) == 2 && all(ux %in% c(0, 1))) "binary" else "continuous"
  # character adjusters -> factors
  for (v in setdiff(intersect(names(df), id$observed), c(x, y))) if (is.character(df[[v]])) df[[v]] <- factor(df[[v]])
  # complete cases on the variables any strategy may use
  use <- unique(c(id$observed, sp[["unit"]], sp[["time"]], sp[["trials"]], sp[["offset"]]))
  use <- intersect(use, names(df))
  cc <- stats::complete.cases(df[, use, drop = FALSE])
  notes <- character()
  if (any(!cc)) notes <- c(notes, sprintf("%d of %d rows had missing values in DAG variables and were dropped (complete-case analysis assumes missingness does not depend on the outcome).", sum(!cc), length(cc)))
  df <- df[cc, , drop = FALSE]
  if (!is.null(sp[["unit"]])) df$.unit <- factor(df[[sp[["unit"]]]])
  if (!is.null(sp[["time"]])) df$.time <- factor(df[[sp[["time"]]]])
  if (!is.null(sp[["offset"]])) {
    if (any(df[[sp[["offset"]]]] <= 0)) stop("offset (effort) column '", sp[["offset"]], "' must be > 0")
    df$.log_offset <- log(df[[sp[["offset"]]]])
  }
  contrast <- if (ex_type == "binary") 1 else (sp[["contrast"]] %||% signif(stats::sd(df[[x]]), 2))
  list(df = df, sp = sp, fam = fam, id = id, x = x, y = y, ex_type = ex_type, contrast = contrast, notes = notes)
}

causal_formula <- function(pr, set, panel) {
  x <- pr$x; y <- pr$y
  set <- setdiff(set, pr$id$absorbed)
  lhs <- if (pr$fam$family == "binomial" && !is.null(pr$sp[["trials"]]))
    sprintf("cbind(%s, %s - %s)", y, pr$sp[["trials"]], y) else y
  rhs <- c(x, set)
  if (panel %in% c("fixed", "twoway")) rhs <- c(rhs, ".unit")
  if (panel %in% c("twoway", "time")) rhs <- c(rhs, ".time")
  if (!is.null(pr$sp[["offset"]])) rhs <- c(rhs, "offset(.log_offset)")
  stats::as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")))
}

causal_glm <- function(pr, fo, df = pr$df) {
  w <- character()
  m <- withCallingHandlers(
    if (pr$fam$family == "negbin") MASS::glm.nb(fo, data = df, link = log)
    else stats::glm(fo, family = switch(pr$fam$family,
                                        gaussian = stats::gaussian(link = pr$fam$link),
                                        poisson = stats::poisson(link = pr$fam$link),
                                        binomial = stats::binomial(link = pr$fam$link),
                                        Gamma = stats::Gamma(link = pr$fam$link)), data = df),
    warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  m <- tryCatch(marginaleffects::set_modeldata(m, df), error = function(e) m)
  attr(m, "eco_warnings") <- unique(w)
  m
}

causal_vcov <- function(m, pr, panel, st) {
  how <- st$vcov
  # Poisson: the total-effect model marginalises over mediators and other
  # causes of the outcome, which makes counts overdispersed given (X, Z) even
  # when the outcome is Poisson given ALL its parents -> sandwich SEs.
  if (how == "auto") how <- if (panel != "none" && !is.null(pr$sp[["unit"]])) "cluster"
                            else if (pr$fam$family == "poisson") "robust" else "model"
  V <- switch(how,
    model = stats::vcov(m),
    robust = { eco_require("sandwich"); sandwich::vcovHC(m, type = "HC3") },
    cluster = { eco_require("sandwich"); sandwich::vcovCL(m, cluster = pr$df$.unit, type = "HC1") })
  attr(V, "how") <- how
  V
}

# Counterfactual datasets: lo/hi exposure (and offset fixed at offset_per).
causal_cf <- function(pr, df = pr$df) {
  x <- pr$x
  lo <- df; hi <- df
  if (pr$ex_type == "binary") { lo[[x]] <- 0; hi[[x]] <- 1 } else hi[[x]] <- df[[x]] + pr$contrast
  if (!is.null(pr$sp[["offset"]])) {
    per <- pr$sp[["offset_per"]] %||% 1
    lo$.log_offset <- hi$.log_offset <- log(per)
  }
  list(lo = lo, hi = hi)
}

# Point estimates by hand (used for the bootstrap and influence terms).
causal_point <- function(m, pr, df = pr$df) {
  cf <- causal_cf(pr, df)
  a <- stats::predict(m, newdata = cf$hi, type = "response")
  b <- stats::predict(m, newdata = cf$lo, type = "response")
  list(a = a, b = b, beta = unname(stats::coef(m)[pr$x]), ate = mean(a) - mean(b), ratio = mean(a) / mean(b),
       m1 = mean(a), m0 = mean(b))
}

# Effect estimates from one fitted model: beta (link scale, per unit of x),
# ate (response scale, per contrast), ratio and potential-outcome means.
causal_effects <- function(m, pr, panel, st, names) {
  eco_require("marginaleffects")
  V <- causal_vcov(m, pr, panel, st)
  z <- stats::qnorm(0.975)
  x <- pr$x
  b <- unname(stats::coef(m)[x]); seb <- sqrt(V[x, x])
  cf <- causal_cf(pr)
  nd <- pr$df
  if (!is.null(pr$sp[["offset"]])) nd$.log_offset <- log(pr$sp[["offset_per"]] %||% 1)
  vars <- if (pr$ex_type == "binary") stats::setNames(list(c(0, 1)), x) else stats::setNames(list(pr$contrast), x)
  ac <- as.data.frame(marginaleffects::avg_comparisons(m, variables = vars, newdata = nd, vcov = V))
  rc <- if (!is.na(pr$fam$ratio))
    as.data.frame(marginaleffects::avg_comparisons(m, variables = vars, newdata = nd, vcov = V, comparison = "ratioavg"))
  pt <- causal_point(m, pr)
  # influence-function term for sampling of the covariate distribution
  # (makes the interval target the population, not just these sites)
  n <- nrow(pr$df)
  agg <- function(infl) {
    if (!is.null(pr$sp[["unit"]]) && panel != "none") infl <- tapply(infl, pr$df$.unit, sum)
    sum(infl^2) / n^2
  }
  extra <- if (st$ci_method == "unconditional") list(
    ate = agg((pt$a - pt$b) - pt$ate),
    ratio = agg((pt$a - pt$ratio * pt$b) / pt$m0),
    m1 = agg(pt$a - pt$m1), m0 = agg(pt$b - pt$m0)) else list(ate = 0, ratio = 0, m1 = 0, m0 = 0)
  rows <- list()
  rows$beta <- data.frame(param = names[["beta"]], estimate = b, sd = seb, lower = b - z * seb, upper = b + z * seb)
  se_ate <- sqrt(ac$std.error[1]^2 + extra$ate)
  rows$ate <- data.frame(param = names[["ate"]], estimate = ac$estimate[1], sd = se_ate,
                         lower = ac$estimate[1] - z * se_ate, upper = ac$estimate[1] + z * se_ate)
  if (!is.null(rc)) {
    r <- rc$estimate[1]; se_r <- sqrt(rc$std.error[1]^2 + extra$ratio)
    rows$ratio <- data.frame(param = names[["ratio"]], estimate = r, sd = se_r,
                             lower = exp(log(r) - z * se_r / r), upper = exp(log(r) + z * se_r / r))
  }
  if ("mean0" %in% names(names)) {
    ap <- as.data.frame(marginaleffects::avg_predictions(m, variables = stats::setNames(list(c(0, 1)), x), newdata = nd, vcov = V, by = x))
    for (k in c(0, 1)) {
      e <- ap$estimate[ap[[x]] == k]; s <- sqrt(ap$std.error[ap[[x]] == k]^2 + extra[[paste0("m", k)]])
      pos <- e > 0 && pr$fam$family != "gaussian"
      rows[[paste0("mean", k)]] <- data.frame(param = names[[paste0("mean", k)]], estimate = e, sd = s,
        lower = if (pos) exp(log(e) - z * s / e) else e - z * s, upper = if (pos) exp(log(e) + z * s / e) else e + z * s)
    }
  }
  est <- do.call(rbind, rows)
  if (st$ci_method == "bootstrap") est <- causal_bootstrap(m, pr, panel, est, names, st)
  attr(est, "vcov_how") <- attr(V, "how")
  attr(est, "point") <- pt
  est
}

# Nonparametric bootstrap (rows, or whole units for panels): percentile CIs.
causal_bootstrap <- function(m, pr, panel, est, names, st) {
  set.seed(st$seed)
  fo <- stats::formula(m)
  df <- pr$df
  clus <- !is.null(pr$sp[["unit"]]) && panel != "none"
  B <- st$n_boot
  res <- matrix(NA_real_, B, 5, dimnames = list(NULL, c("beta", "ate", "ratio", "m0", "m1")))
  for (bb in seq_len(B)) {
    if (clus) {
      u <- levels(df$.unit); s <- sample(u, length(u), replace = TRUE)
      d <- do.call(rbind, lapply(seq_along(s), function(i) { di <- df[df$.unit == s[i], , drop = FALSE]; di$.unit <- paste0("b", i); di }))
      d$.unit <- factor(d$.unit)
    } else d <- df[sample(nrow(df), replace = TRUE), , drop = FALSE]
    mb <- tryCatch(suppressWarnings(causal_glm(pr, fo, d)), error = function(e) NULL)
    if (is.null(mb)) next
    p <- tryCatch(causal_point(mb, pr, d), error = function(e) NULL)
    if (!is.null(p)) res[bb, ] <- c(p$beta, p$ate, p$ratio, p$m0, p$m1)
  }
  map <- c(beta = "beta", ate = "ate", ratio = "ratio", mean0 = "m0", mean1 = "m1")
  for (k in names(map)) if (k %in% names(names)) {
    i <- est$param == names[[k]]; v <- res[, map[[k]]]; v <- v[is.finite(v)]
    est$sd[i] <- stats::sd(v); est$lower[i] <- stats::quantile(v, 0.025); est$upper[i] <- stats::quantile(v, 0.975)
  }
  est
}

# ---- local tests of the DAG ----------------------------------------------------
causal_local_tests <- function(pr, st) {
  g <- pr$id$g_latent
  t <- tryCatch(dagitty::impliedConditionalIndependencies(g, type = "missing.edge"), error = function(e) NULL)
  if (is.null(t) || !length(t)) return(NULL)
  key <- vapply(t, function(ci) paste(sort(c(ci$X, ci$Y)), collapse = "|"), "")
  t <- t[!duplicated(key)]
  obs <- pr$id$observed
  t <- t[vapply(t, function(ci) all(c(ci$X, ci$Y, unlist(ci$Z)) %in% obs), logical(1))]
  if (!length(t)) return(NULL)
  class(t) <- "dagitty.cis"
  d <- pr$df[, obs, drop = FALSE]
  has_f <- any(vapply(d, function(v) is.factor(v) || is.character(v), logical(1)))
  type <- if (st$ci_test_type == "auto") (if (has_f) "cis.pillai" else "cis") else st$ci_test_type
  if (type == "cis") d[] <- lapply(d, as.numeric)
  r <- tryCatch(dagitty::localTests(g, d, type = type, tests = t, abbreviate.names = FALSE),
                error = function(e) { attr(e, "msg") <- conditionMessage(e); NULL })
  if (is.null(r)) return(NULL)
  out <- data.frame(test = rownames(r), estimate = r$estimate, p_value = r$p.value,
                    lower = r[["2.5%"]], upper = r[["97.5%"]], stringsAsFactors = FALSE)
  out$p_holm <- stats::p.adjust(out$p_value, "holm")
  pv <- pmax(out$p_value, 1e-300)
  C <- -2 * sum(log(pv)); k <- nrow(out)
  attr(out, "fisher_c") <- c(C = C, df = 2 * k, p = stats::pchisq(C, 2 * k, lower.tail = FALSE))
  attr(out, "type") <- type
  out
}

causal_dharma <- function(m, st) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) return(NULL)
  tryCatch({
    set.seed(st$seed)
    r <- DHARMa::simulateResiduals(m, n = st$n_sim_resid, plot = FALSE)
    list(res = r, p = c(uniformity = DHARMa::testUniformity(r, plot = FALSE)$p.value,
                        dispersion = DHARMa::testDispersion(r, plot = FALSE)$p.value,
                        outliers = DHARMa::testOutliers(r, plot = FALSE, type = "bootstrap", nBoot = 100)$p.value))
  }, error = function(e) list(error = conditionMessage(e)))
}

# E-value (VanderWeele & Ding 2017) for a ratio estimate and its CI limit
# closest to 1.
evalue <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
causal_evalue <- function(est_row, pr, ratio_row = NULL) {
  if (!is.null(ratio_row) && nrow(ratio_row)) {
    rr <- ratio_row$estimate; lo <- ratio_row$lower; hi <- ratio_row$upper; kind <- "ratio"
  } else {
    # continuous outcome: RR ~ exp(0.91 d), d = effect / SD(outcome) (VanderWeele & Ding 2017)
    s <- stats::sd(pr$df[[pr$y]])
    f <- function(v) exp(0.91 * v / s)
    rr <- f(est_row$estimate); lo <- f(est_row$lower); hi <- f(est_row$upper); kind <- "approx (standardised difference)"
  }
  e_ci <- if (lo <= 1 && hi >= 1) 1 else evalue(if (rr > 1) lo else hi)
  c(point = evalue(rr), ci = e_ci, kind = kind)
}

# ---- main fit -------------------------------------------------------------------
causal_fit <- function(data, settings) {
  eco_require(c("dagitty", "marginaleffects", "MASS"))
  st <- merge_list(CAUSAL_DEFAULT_SETTINGS, settings)
  pr <- causal_prepare(data, st)
  id <- pr$id
  notes <- c(pr$notes, id$total$notes)
  panel <- id$panel_used
  set <- if (!is.null(st$force_adjust)) unlist(st$force_adjust) else id$total$set
  if (is.null(set)) stop(causal_not_identified_msg(id))
  if (!is.null(st$force_adjust)) notes <- c(notes, sprintf("NEGATIVE CONTROL: adjustment set forced to {%s}%s.", paste(set, collapse = ", "),
                                                         if (!isTRUE(dagitty::isAdjustmentSet(id$g, set, pr$x, pr$y))) " (the DAG says this set does NOT identify the effect)" else ""))
  if (!is.null(st$force_panel)) notes <- c(notes, sprintf("NEGATIVE CONTROL: panel estimator forced to '%s'.", panel))
  if (panel != "none" && is.null(pr$sp[["unit"]])) stop("panel estimator '", panel, "' needs `unit:` (the column identifying sites/individuals)")
  if (panel == "twoway" && is.null(pr$sp[["time"]])) stop("panel 'twoway' needs `time:` (the column identifying years/periods)")
  unabsorbed <- intersect(set, id$latent)
  if (panel == "none" && length(intersect(set, id$absorbed))) {
    # absorbed latents only make sense with fixed effects
    set <- setdiff(set, id$absorbed)
  }
  if (length(setdiff(unabsorbed, id$absorbed)) && is.null(st$force_adjust)) stop("adjustment set includes unmeasured variables: ", paste(unabsorbed, collapse = ", "))
  # time-invariant adjusters are collinear with unit fixed effects
  if (panel %in% c("fixed", "twoway")) {
    inv <- Filter(function(v) v %in% names(pr$df) && all(tapply(pr$df[[v]], pr$df$.unit, function(z) length(unique(z)) == 1)), setdiff(set, id$absorbed))
    if (length(inv)) notes <- c(notes, sprintf("%s constant within every %s: already controlled by the %s fixed effects, so not entered separately.", paste(inv, collapse = ", "), pr$sp[["unit"]], pr$sp[["unit"]]))
    set <- setdiff(set, inv)
  }
  nm <- causal_estimand_names(pr$x, pr$y, pr$fam, pr$ex_type)
  fo <- causal_formula(pr, set, panel)
  m <- causal_glm(pr, fo)
  est <- causal_effects(m, pr, panel, st, nm)
  pt_main <- attr(est, "point"); vc_how <- attr(est, "vcov_how")
  est$adjusted_for <- paste(c(setdiff(set, id$absorbed), if (panel %in% c("fixed", "twoway")) paste0(pr$sp[["unit"]], " fixed effects"),
                              if (panel == "twoway") paste0(pr$sp[["time"]], " fixed effects")), collapse = ", ")
  if (est$adjusted_for[1] == "") est$adjusted_for <- "(nothing)"
  est$role <- "total effect"
  est$strategy <- "DAG set (main estimate)"
  est$valid_by_dag <- if (!is.null(st$force_adjust)) isTRUE(dagitty::isAdjustmentSet(id$g, unlist(st$force_adjust), pr$x, pr$y)) else TRUE
  if (pr$fam$link == "logit") est$role[est$param == nm[["beta"]]] <- "total effect (conditional log odds ratio; non-collapsible)"
  conv <- isTRUE(m$converged) && all(is.finite(est$estimate)) && all(is.finite(est$sd))
  warn <- attr(m, "eco_warnings")
  models <- list(total = m)
  # direct effect
  if (id$effect %in% c("direct", "both")) {
    dset <- id$direct$set
    if (is.null(dset)) notes <- c(notes, "The direct effect is not identifiable from this DAG by regression adjustment (dagitty found no valid set); only the total effect is reported.")
    else {
      nmd <- causal_estimand_names(pr$x, pr$y, pr$fam, "continuous", "direct_")
      if (pr$ex_type == "binary") nmd <- nmd[c("beta", "ate", if (!is.na(pr$fam$ratio)) "ratio")]
      md <- causal_glm(pr, causal_formula(pr, dset, panel))
      ed <- causal_effects(md, pr, panel, st, nmd)
      ed$adjusted_for <- paste(setdiff(dset, id$absorbed), collapse = ", "); ed$role <- "controlled direct effect"
      ed$strategy <- "Direct-effect set"
      ed$valid_by_dag <- TRUE
      est <- rbind(est, ed); models$direct <- md
      conv <- conv && isTRUE(md$converged)
    }
  }
  # what if
  wi <- NULL
  if (isTRUE(st$whatif) && is.null(st$force_adjust)) {
    strat <- causal_strategies(id$g, pr$x, pr$y, id$total$set, id$observed, panel)
    if (panel != "none") {
      strat <- rbind(strat, data.frame(tag = "no_fixed_effects", label = sprintf("Ignore %s (pooled GLM, no fixed effects)", pr$sp[["unit"]]), short = "No fixed effects (pooled)",
                                       kind = "panel", valid = FALSE, set = I(list(setdiff(id$total$set, id$absorbed)))))
    }
    rows <- list()
    for (i in seq_len(nrow(strat))) {
      if (strat$tag[i] == "dag") next
      s <- strat$set[[i]]
      pnl <- if (strat$tag[i] == "no_fixed_effects") "none" else panel
      if (pnl != "none") {
        # fixed effects follow the absorbed confounders in this strategy's set
        has_u <- any(unlist(pr$sp[["unit_constant"]]) %in% s); has_t <- any(unlist(pr$sp[["time_constant"]]) %in% s)
        pnl <- if (has_u && has_t) "twoway" else if (has_u) "fixed" else if (has_t) "time" else "none"
      }
      if (length(intersect(s, setdiff(id$latent, if (pnl != "none") id$absorbed)))) next
      if (pnl %in% c("fixed", "twoway")) s <- setdiff(s, Filter(function(v) v %in% names(pr$df) && all(tapply(pr$df[[v]], pr$df$.unit, function(z) length(unique(z)) == 1)), s))
      mi <- tryCatch(causal_glm(pr, causal_formula(pr, s, pnl)), error = function(e) NULL)
      if (is.null(mi)) next
      ei <- tryCatch(causal_effects(mi, pr, pnl, merge_list(st, list(ci_method = if (st$ci_method == "bootstrap") "unconditional" else st$ci_method)),
                                    nm[intersect(c("beta", "ate", "ratio"), names(nm))]), error = function(e) NULL)
      if (is.null(ei)) next
      ei <- ei[ei$param == nm[["ate"]], ]
      ei$param <- paste0("whatif_", strat$tag[i], "__", nm[["ate"]])
      ei$adjusted_for <- if (length(setdiff(s, id$absorbed))) paste(setdiff(s, id$absorbed), collapse = ", ") else "(nothing)"
      if (pnl %in% c("fixed", "twoway")) ei$adjusted_for <- paste0(ei$adjusted_for, " + ", pr$sp[["unit"]], " fixed effects")
      if (pnl %in% c("time", "twoway")) ei$adjusted_for <- paste0(ei$adjusted_for, " + ", pr$sp[["time"]], " fixed effects")
      ei$role <- paste("what if:", strat$label[i])
      ei$strategy <- strat$short[i]
      ei$valid_by_dag <- strat$valid[i]
      rows[[length(rows) + 1]] <- ei
    }
    wi <- if (length(rows)) do.call(rbind, rows) else NULL
    if (!is.null(wi)) est <- rbind(est, wi)
  }
  rownames(est) <- NULL
  # checks
  lt <- if (isTRUE(st$local_tests)) tryCatch(causal_local_tests(pr, st), error = function(e) NULL) else NULL
  dh <- if (isTRUE(st$diagnostics)) causal_dharma(m, st) else NULL
  ev <- tryCatch(causal_evalue(est[est$param == nm[["ate"]], ], pr,
                               if (!is.na(pr$fam$ratio)) est[est$param == nm[["ratio"]], ]), error = function(e) NULL)
  st$vcov_used <- vc_how
  notes <- c(notes, causal_fit_notes(m, pr, id, set, panel, est, nm, lt, dh, ev, warn, st))
  list(estimates = est,
       diagnostics = list(converged = conv, notes = notes, local_tests = lt,
                          dharma_p = if (!is.null(dh)) dh$p, evalue = ev,
                          model_checks_ok = if (is.null(dh$p)) NA else all(dh$p >= 0.05)),
       object = m, models = models, prep = pr[c("sp", "fam", "x", "y", "ex_type", "contrast")],
       id = id[c("x", "y", "absorbed", "panel", "panel_used", "latent", "roles", "observed", "effect")],
       sets = list(chosen = set, how = id$total$how, minimal = id$total$sets, direct = id$direct$set,
                   optimal = tryCatch(causal_optimal_set(id$g, pr$x, pr$y), error = function(e) NULL)),
       dag_string = causal_dag_string(pr$sp[["dag"]], pr$sp[["edges"]]), dharma = dh, settings = st,
       data_used = pr$df, names = nm, point = list(m0 = pt_main$m0, m1 = pt_main$m1))
}

causal_not_identified_msg <- function(id) {
  lat <- id$latent
  sprintf(paste0("The causal effect of %s on %s is NOT identifiable by adjustment with the variables you measured: ",
                 "every valid adjustment set includes an unmeasured variable (%s). Options: measure it; if it is constant within ",
                 "sites and you have repeated measures, declare it in `unit_constant:` and use `panel: fixed`; or revise the DAG."),
          id$x, id$y, if (length(lat)) paste(lat, collapse = ", ") else "none declared")
}

causal_fit_notes <- function(m, pr, id, set, panel, est, nm, lt, dh, ev, warn, st) {
  n <- character()
  n <- c(n, sprintf("Adjustment set for the total effect: {%s} (%s; %d minimal set%s found by dagitty%s).",
                    paste(setdiff(set, id$absorbed), collapse = ", "), switch(id$total$how, minimal = "smallest minimal set", optimal = "optimal set", user = "set given in config", none = "none"),
                    length(id$total$sets), if (length(id$total$sets) == 1) "" else "s",
                    if (length(id$absorbed) && panel != "none") paste0("; ", paste(c(
                      if (length(intersect(id$absorbed, unlist(pr$sp[["unit_constant"]])))) paste(intersect(id$absorbed, unlist(pr$sp[["unit_constant"]])), collapse = ", ") |> paste("absorbed by", pr$sp[["unit"]], "fixed effects"),
                      if (panel == "twoway" && length(intersect(id$absorbed, unlist(pr$sp[["time_constant"]])))) paste(intersect(id$absorbed, unlist(pr$sp[["time_constant"]])), collapse = ", ") |> paste("absorbed by", pr$sp[["time"]], "fixed effects")), collapse = "; ")) else ""))
  if (identical(id$total$how, "user") && !isTRUE(id$total$valid))
    n <- c(n, "WARNING: the adjustment set given in config is NOT valid for this DAG (dagitty::isAdjustmentSet is FALSE). The estimate is not the causal effect.")
  n <- c(n, sprintf("%s fit (%s family, %s link) to %d rows; standard errors: %s; intervals: %s.",
                    if (pr$fam$family == "negbin") "MASS::glm.nb" else "stats::glm", pr$fam$family, pr$fam$link, nrow(pr$df),
                    st$vcov_used %||% st$vcov, st$ci_method))
  if (!isTRUE(m$converged)) n <- c(n, "NOT CONVERGED: the GLM optimiser did not converge; do not interpret.")
  if (length(warn)) n <- c(n, paste("Engine warnings:", paste(warn, collapse = " | ")))
  if (pr$fam$family == "binomial" && any(abs(stats::coef(m)[!grepl("^\\.unit|^\\.time", names(stats::coef(m)))]) > 10, na.rm = TRUE))
    n <- c(n, "Very large logit coefficients: possible separation; intervals unreliable.")
  # positivity / overlap
  x <- pr$df[[pr$x]]
  if (pr$ex_type == "continuous") {
    out <- mean(x + pr$contrast > max(x))
    if (out > 0.1) n <- c(n, sprintf("Extrapolation: for %.0f%% of rows, exposure + %s is above the largest observed value; the effect for them relies on the model's shape.", 100 * out, format(pr$contrast)))
  } else {
    adj <- setdiff(set, c(id$absorbed, id$latent))
    if (length(adj)) {
      ps <- tryCatch(stats::fitted(stats::glm(stats::reformulate(adj, pr$x), stats::binomial, pr$df)), error = function(e) NULL)
      if (!is.null(ps)) {
        ext <- mean(ps < 0.05 | ps > 0.95)
        n <- c(n, sprintf("Overlap (positivity): propensity of exposure given the adjustment set ranges %.2f-%.2f; %.0f%% of rows have propensity < 0.05 or > 0.95%s.",
                          min(ps), max(ps), 100 * ext, if (ext > 0.1) " - poor overlap: the effect partly extrapolates between groups that never overlap" else ""))
      }
    }
  }
  if (!is.null(lt)) {
    fc <- attr(lt, "fisher_c")
    bad <- lt$test[lt$p_holm < 0.05]
    n <- c(n, sprintf("DAG check (%d testable conditional independencies, %s tests): Fisher's C = %.1f, df = %d, p = %.3f.%s",
                      nrow(lt), attr(lt, "type"), fc["C"], fc["df"], fc["p"],
                      if (length(bad)) paste0(" Contradicted by the data (Holm-adjusted p < 0.05): ", paste(bad, collapse = "; "),
                                              ". The DAG is missing at least one arrow or common cause; revise it before trusting the estimate.") else " No implied independence is contradicted by the data (this cannot prove the DAG is right)."))
    if (panel != "none") n <- c(n, "Local tests treat rows as independent; with repeated measures their p-values are too small, so a rejection is weaker evidence than it looks.")
  } else n <- c(n, "The DAG implies no testable conditional independencies among measured variables (it is saturated or the tests involve unmeasured nodes), so the data cannot check it.")
  if (!is.null(dh$p)) {
    bad <- names(dh$p)[dh$p < 0.05]
    n <- c(n, sprintf("Outcome-model residual checks (DHARMa): uniformity p = %.3f, dispersion p = %.3f, outliers p = %.3f.%s",
                      dh$p[["uniformity"]], dh$p[["dispersion"]], dh$p[["outliers"]],
                      if (length(bad)) paste0(" Flagged: ", paste(bad, collapse = ", "), if (pr$fam$family == "poisson" && "dispersion" %in% bad) ". Overdispersed counts: use family negbin." else ". Inspect the residual plot; with several tests one p < 0.05 also occurs by chance.") else ""))
  }
  if (!is.null(ev)) n <- c(n, sprintf("E-value: an unmeasured confounder would need to be associated with both %s and %s by a risk ratio of %.2f each (%.2f to move the 95%% CI to include no effect) to explain away the estimate (%s).",
                                      pr$x, pr$y, as.numeric(ev[["point"]]), as.numeric(ev[["ci"]]), ev[["kind"]]))
  n
}
