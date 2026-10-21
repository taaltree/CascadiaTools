# Reading data, checks, interpretation, figures and methods text for `sem`.

# ---- real data ----------------------------------------------------------------
sem_read_data <- function(config) {
  dc <- config$data %||% list()
  if (is.null(dc$file)) stop("config needs data: file: <your csv>")
  path <- dc$file
  if (!grepl("^/", path)) path <- file.path(config$.dir %||% ".", path)
  if (!file.exists(path)) stop("data file not found: ", path)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE, strip.white = TRUE,
                        na.strings = c("", "NA", "na", "N/A", "n/a", "."))
  for (nm in names(dc$transform)) {
    ex <- dc$transform[[nm]]
    v <- tryCatch(eval(parse(text = ex), envir = df, enclos = baseenv()),
                  error = function(e) stop(sprintf("transform %s = '%s' failed: %s", nm, ex, conditionMessage(e))))
    df[[nm]] <- v
  }
  engine <- config$engine %||% (if (!is.null(config$lavaan) && is.null(config$model)) "lavaan" else "piecewise")
  spec <- list(engine = engine, models = config$model, lavaan = config$lavaan, group = config$group,
               estimator = config$estimator, effects = config$effects, labels = config$labels,
               units = config$units, drop_missing = isTRUE(dc$drop_missing), source = path)
  list(df = df, spec = spec)
}

sem_model_vars <- function(spec, df = NULL) {
  if (identical(spec$engine, "lavaan")) {
    pt <- lavaan::lavaanify(spec$lavaan, auto = FALSE)
    lv <- unique(pt$lhs[pt$op == "=~"])
    v <- unique(c(pt$lhs[pt$op %in% c("~", "=~", "~~")], pt$rhs[pt$op %in% c("~", "=~", "~~")]))
    return(list(vars = setdiff(v[nzchar(v)], lv), latent = lv))
  }
  comps <- sem_parse_models(spec$models)
  list(vars = unique(unlist(lapply(comps, function(c) c(c$response, c$preds, c$groups)))), latent = character(), comps = comps)
}

sem_check_data <- function(data) {
  pr <- character(); sp <- data$spec; df <- data$df
  if (!sp$engine %in% c("piecewise", "lavaan")) return(sprintf("engine '%s' unknown: use piecewise or lavaan", sp$engine))
  mv <- tryCatch(sem_model_vars(sp), error = function(e) { pr <<- c(pr, paste("model specification problem:", conditionMessage(e))); NULL })
  if (is.null(mv)) return(pr)
  if (sp$engine == "piecewise") {
    r <- tryCatch({ sem_topo_order(mv$comps); NULL }, error = function(e) conditionMessage(e))
    if (!is.null(r)) pr <- c(pr, r)
  }
  miss <- setdiff(c(mv$vars, sp$group), names(df))
  if (length(miss)) {
    hint <- unlist(lapply(miss, function(v) {
      cand <- names(df)[agrepl(v, names(df), max.distance = 0.25, ignore.case = TRUE)]
      cand <- unique(c(cand, names(df)[vapply(names(df), function(n) agrepl(n, v, max.distance = 0.25, ignore.case = TRUE), logical(1))]))
      if (length(cand)) sprintf(" Did you mean %s for %s? Rename the column in the CSV or change the config to match.", paste(cand, collapse = " or "), v)
    }))
    return(c(pr, sprintf("column(s) %s not found in the data file. Columns present: %s. (R turns spaces and symbols in column names into dots.)%s",
                         paste(miss, collapse = ", "), paste(names(df), collapse = ", "), paste(hint, collapse = ""))))
  }
  # grouping columns (random intercepts, lavaan group:): IDs that differ only by
  # case, spaces, punctuation or leading zeros are almost always typos
  gcols <- unique(c(if (sp$engine == "piecewise") unlist(lapply(mv$comps, `[[`, "groups")), sp$group))
  for (gc in gcols) {
    ids <- unique(as.character(stats::na.omit(df[[gc]])))
    key <- gsub("(?<![0-9])0+(?=[0-9])", "", gsub("[^a-z0-9]", "", tolower(ids)), perl = TRUE)
    dup <- key[duplicated(key)]
    for (k in unique(dup)) {
      v <- ids[key == k]
      pr <- c(pr, sprintf("grouping column %s has IDs that look like the same site spelled differently: %s (rows %s). Correct the spelling so each site has one ID; otherwise they are treated as different sites.",
                          gc, paste0("'", v, "'", collapse = ", "), paste(utils::head(which(df[[gc]] %in% v[-which.max(table(df[[gc]])[v])]), 8), collapse = ", ")))
    }
  }
  if (!isTRUE(sp$drop_missing)) for (v in mv$vars) {
    k <- which(is.na(df[[v]]))
    if (length(k)) pr <- c(pr, sprintf("column %s has %d missing value(s) (rows %s%s). Every component of a piecewise SEM must use the same rows: fill them, or set data: drop_missing: true to drop incomplete rows.",
                                       v, length(k), paste(utils::head(k, 8), collapse = ", "), if (length(k) > 8) ", ..." else ""))
  }
  if (length(pr)) return(pr)
  isint <- function(x) all(abs(x - round(x)) < 1e-8, na.rm = TRUE)
  if (sp$engine == "piecewise") {
    for (c in mv$comps) {
      y <- df[[c$response]]
      if (!is.numeric(y)) {
        vals <- utils::head(unique(as.character(stats::na.omit(y))), 6)
        pr <- c(pr, sprintf("response %s must be numeric but contains text (values: %s).%s", c$response, paste0("'", vals, "'", collapse = ", "),
                            if (c$family == "binomial") " For a presence/absence or yes/no variable, recode it as 1 and 0 (e.g. 1 = present, 0 = absent) in the CSV, or add it under data: transform: (e.g. shrub01: as.integer(shrub == 'yes'))." else " Remove units or text from the numbers in the CSV."))
        next
      }
      if (c$family == "binomial" && !all(y %in% c(0, 1, NA))) pr <- c(pr, sprintf("binomial response %s must be 0/1; rows %s have other values", c$response, paste(utils::head(which(!y %in% c(0, 1, NA)), 8), collapse = ", ")))
      if (c$family %in% c("poisson", "nbinom") && (any(y < 0, na.rm = TRUE) || !isint(y)))
        pr <- c(pr, sprintf("count response %s must be whole numbers >= 0 (family %s); rows %s have %s",
                            c$response, c$family, paste(utils::head(which(y < 0 | abs(y - round(y)) > 1e-8), 8), collapse = ", "),
                            paste(utils::head(y[which(y < 0 | abs(y - round(y)) > 1e-8)], 8), collapse = ", ")))
      for (x in c$preds) {
        v <- df[[x]]
        if (is.numeric(v) && isTRUE(stats::sd(v, na.rm = TRUE) == 0)) pr <- c(pr, sprintf("predictor %s is constant", x))
        if (!is.numeric(v) && length(unique(stats::na.omit(v))) > 20) pr <- c(pr, sprintf("predictor %s is text with %d distinct values: if it identifies sites, use (1|%s) instead", x, length(unique(v)), x))
        if (!is.numeric(v) && x %in% vapply(mv$comps, `[[`, "", "response")) pr <- c(pr, sprintf("%s is a response and a predictor but is not numeric", x))
      }
      for (g in c$groups) if (length(unique(df[[g]])) < 3) pr <- c(pr, sprintf("grouping column %s has fewer than 3 groups: use it as a fixed factor instead of (1|%s)", g, g))
    }
  } else {
    num <- intersect(mv$vars, names(df))
    nn <- num[!vapply(df[num], is.numeric, logical(1))]
    if (length(nn)) pr <- c(pr, sprintf("lavaan needs numeric variables; %s are not numeric (use 0/1 dummies for categories, or group: for multigroup)", paste(nn, collapse = ", ")))
    vv <- vapply(df[setdiff(num, nn)], stats::var, 1, na.rm = TRUE)
    if (length(vv) > 1 && max(vv) / min(vv) > 1000)
      pr <- c(pr, sprintf("variances differ by a factor of %.0f (%s largest, %s smallest); lavaan's optimiser can fail. Rescale with data: transform: (e.g. dnbr_100: dnbr / 100).",
                          max(vv) / min(vv), names(vv)[which.max(vv)], names(vv)[which.min(vv)]))
  }
  n <- sum(stats::complete.cases(df[, intersect(mv$vars, names(df)), drop = FALSE]))
  if (n < 20) pr <- c(pr, sprintf("only %d complete rows: too few for a structural equation model", n))
  pr
}

# ---- labels -------------------------------------------------------------------
sem_lab <- function(v, data, opts = NULL) {
  vapply(as.character(v), function(k) {
    d <- data$spec$labels[[k]] %||% gsub("_", " ", k)
    if (is.null(opts)) d else fig_label(k, d, opts)
  }, "", USE.NAMES = FALSE)
}

sem_numeric_paths <- function(fit) {
  t <- fit$table
  if (fit$engine == "lavaan") t <- t[t$op == "~", ] else t <- t[!is.na(t$std_factor), ]
  t
}

# ---- interpretation -----------------------------------------------------------
sem_effect_phrase <- function(b, lo, hi, link, ylab, verb = "is associated with") {
  out <- sem_effect_phrase0(b, lo, hi, link, ylab)
  sub("^is associated with", verb, out)
}

sem_effect_phrase0 <- function(b, lo, hi, link, ylab) {
  if (link == "log") {
    sprintf("is associated with the expected %s multiplied by %s (95%% CI %s to %s; %s%%)", ylab, fmt(exp(b), 2), fmt(exp(lo), 2), fmt(exp(hi), 2),
            sprintf("%+.0f", 100 * (exp(b) - 1)))
  } else if (link == "logit") {
    sprintf("is associated with the odds of %s multiplied by %s (95%% CI %s to %s)", ylab, fmt(exp(b), 2), fmt(exp(lo), 2), fmt(exp(hi), 2))
  } else if (link == "probit") {
    sprintf("is associated with a shift of %s (95%% CI %s to %s) in the probit of %s", fmt_sig(b), fmt_sig(lo), fmt_sig(hi), ylab)
  } else sprintf("is associated with a difference of %s (95%% CI %s to %s) in %s", fmt_sig(b), fmt_sig(lo), fmt_sig(hi), ylab)
}

sem_certainty <- function(lo, hi) if (lo > 0 || hi < 0) "the 95% interval excludes zero" else "the 95% interval includes zero, so the direction is uncertain"

sem_interpret <- function(fit, data) {
  out <- character(); est <- fit$estimates; df <- fit$df
  g <- function(p) est[est$param == p, , drop = FALSE]
  if (fit$engine == "piecewise") {
    ds <- fit$dsep
    out <- c(out, if (isTRUE(ds$saturated)) "The model is saturated (every pair of variables is linked), so the data cannot test it; the path estimates are still valid regressions."
             else if (!is.null(ds$table)) sprintf("Global fit. Shipley's d-separation test checks the %d pairs of variables the model says are unlinked once their causes are accounted for. Fisher's C = %.2f (df = %d), P = %.3f: %s",
                                                 ds$k, ds$C, as.integer(ds$df), ds$p,
                                                 if (ds$p < 0.05) paste0("the data contradict the model. The strongest missing links are ",
                                                                         paste(sprintf("%s -> %s (P = %s)", sem_lab(ds$table$missing_from, data)[order(ds$table$p_value)][1:min(3, sum(ds$table$p_value < 0.05))],
                                                                                       sem_lab(ds$table$missing_to, data)[order(ds$table$p_value)][1:min(3, sum(ds$table$p_value < 0.05))],
                                                                                       fmt_sig(sort(ds$table$p_value)[1:min(3, sum(ds$table$p_value < 0.05))], 2)), collapse = "; "),
                                                                         ". Path and indirect-effect estimates from a rejected model can be biased; revise the model (add the link if it is biologically plausible) before interpreting.")
                                                 else "no evidence that the model leaves out a needed link. This supports, but does not prove, the causal structure drawn.")
             else "The d-separation test could not be computed.")
  } else {
    fm <- fit$fit_measures
    if (!is.null(fm)) out <- c(out, if (fm[["df"]] == 0) "The model is saturated (df = 0), so its global fit cannot be tested."
      else sprintf("Global fit. Chi-square = %.2f (df = %d, P = %.3f); CFI = %.3f, RMSEA = %.3f (90%% CI %.3f-%.3f), SRMR = %.3f. %s Conventional guides (Hu & Bentler 1999) for good fit are CFI >= 0.95, RMSEA <= 0.06 and SRMR <= 0.08; the chi-square test is the formal one.",
                   fm[["chisq"]], as.integer(fm[["df"]]), fm[["pvalue"]], fm[["cfi"]], fm[["rmsea"]], fm[["rmsea.ci.lower"]], fm[["rmsea.ci.upper"]], fm[["srmr"]],
                   if (fm[["pvalue"]] < 0.05) "P < 0.05: the model's implied covariances differ from the data; see the modification indices for what may be missing, and interpret paths cautiously." else "P >= 0.05: no evidence of misfit."))
  }
  # paths
  tp <- sem_numeric_paths(fit)
  for (i in seq_len(nrow(tp))) {
    r <- tp[i, ]
    x <- r$from; y <- r$to
    lk <- if (fit$engine == "piecewise") r$link else "identity"
    sdx <- if (x %in% names(df) && is.numeric(df[[x]])) stats::sd(df[[x]]) else NA
    if (is.na(sdx) && fit$engine == "lavaan" && x %in% fit$latent && y %in% names(df) && is.finite(r$std_factor))
      sdx <- r$std_factor * stats::sd(df[[y]], na.rm = TRUE)   # model-implied SD of the latent variable
    u <- data$spec$units[[x]]
    ylab <- sem_lab(y, data); xlab <- sem_lab(x, data)
    std <- g(paste0("std_", x, "_to_", y))
    grp <- if (fit$engine == "lavaan" && length(fit$group_labels) > 1) sprintf(" [group %s]", fit$group_labels[r$group]) else ""
    if (is.na(sdx)) {
      out <- c(out, sprintf("%s -> %s%s: each 1-unit increase in %s %s; %s.", xlab, ylab, grp, xlab,
                            sem_effect_phrase(r$estimate, r$lower, r$upper, lk, ylab), sem_certainty(r$lower, r$upper)))
      next
    }
    s1 <- sprintf("%s -> %s%s: a 1-SD increase in %s (%s) %s", xlab, ylab, grp, xlab,
                  if (fit$engine == "lavaan" && x %in% fit$latent) sprintf("%s units of its first indicator, the latent variable's scale", fmt_sig(sdx)) else paste(fmt_sig(sdx), "units"),
                  sem_effect_phrase(r$estimate * sdx, r$lower * sdx, r$upper * sdx, lk, ylab))
    if (!is.null(u)) s1 <- paste0(s1, sprintf("; an increase of %s units %s", fmt_sig(u), sem_effect_phrase(r$estimate * u, r$lower * u, r$upper * u, lk, ylab)))
    if (nrow(std)) s1 <- paste0(s1, sprintf(". Standardised coefficient %s (95%% CI %s to %s)", fmt(std$estimate, 2), fmt(std$lower, 2), fmt(std$upper, 2)))
    out <- c(out, paste0(s1, "; ", sem_certainty(r$lower, r$upper), "."))
  }
  # indirect / total
  ei <- fit$effects_info
  for (i in seq_len(if (is.null(ei)) 0 else nrow(ei))) {
    a <- ei$from[i]; b <- ei$to[i]
    gs <- if (!is.null(ei$group) && fit$engine == "lavaan" && length(fit$group_labels) > 1) paste0("_", make.names(fit$group_labels[ei$group[i]])) else ""
    ind <- g(paste0("indirect_", a, "_to_", b, gs)); tot <- g(paste0("total_", a, "_to_", b, gs)); dir <- g(paste0("direct_", a, "_to_", b, gs))
    sind <- g(paste0("std_indirect_", a, "_to_", b, gs)); stot <- g(paste0("std_total_", a, "_to_", b, gs))
    if (!nrow(ind)) next
    lk <- if (fit$engine == "piecewise") ei$sink_link[i] else "identity"
    sda <- if (a %in% names(df) && is.numeric(df[[a]])) stats::sd(df[[a]]) else NA
    al <- sem_lab(a, data); bl <- sem_lab(b, data)
    routes <- paste(vapply(strsplit(ei$routes[i], "; ")[[1]], function(r) paste(sem_lab(strsplit(r, " -> ")[[1]], data), collapse = " -> "), ""), collapse = "; ")
    nongm <- isTRUE(nzchar(ei$nongaussian_mediators[i] %||% ""))
    if (!nongm) {
      s <- sprintf("Indirect and total effects of %s on %s (routes: %s). The indirect (mediated) association is %s (95%% CI %s to %s) on the %s scale per unit of %s; standardised %s (95%% CI %s to %s), %s.",
                   al, bl, routes, fmt_sig(ind$estimate), fmt_sig(ind$lower), fmt_sig(ind$upper),
                   if (lk == "identity") paste0("original (units of ", bl, ")") else if (lk == "log") "log expected-count" else "log-odds",
                   al, fmt(sind$estimate, 2), fmt(sind$lower, 2), fmt(sind$upper, 2), sem_certainty(ind$lower, ind$upper))
    } else {
      # Non-Gaussian mediator: the raw product a*b multiplies a coefficient on the
      # mediator's link scale by one per observed unit of the mediator, so it is
      # not an effect on any one scale (and can even differ in sign from the
      # standardised total). Describe and compare standardised values only.
      s <- sprintf("Indirect and total effects of %s on %s (routes: %s). The route passes through non-Gaussian mediator(s) (%s), so the raw product of path coefficients (%s, 95%% CI %s to %s) combines the mediator's link scale with its observed units and is not an effect on any single scale; it is reported only as a test of the mediated pathway (%s). Use the standardised values (Grace et al. 2018): standardised indirect association %s (95%% CI %s to %s), standardised total %s (95%% CI %s to %s).",
                   al, bl, routes, ei$nongaussian_mediators[i], fmt_sig(ind$estimate), fmt_sig(ind$lower), fmt_sig(ind$upper),
                   sem_certainty(ind$lower, ind$upper), fmt(sind$estimate, 2), fmt(sind$lower, 2), fmt(sind$upper, 2),
                   fmt(stot$estimate, 2), fmt(stot$lower, 2), fmt(stot$upper, 2))
      dir <- g(paste0("std_direct_", a, "_to_", b, gs)); ind <- sind; tot <- stot
    }
    if (!is.na(sda) && !nongm)
      s <- paste(s, sprintf("Per 1-SD increase in %s (%s units), the indirect route alone %s; the total association (direct + indirect) %s.",
                            al, fmt_sig(sda), sem_effect_phrase(ind$estimate * sda, ind$lower * sda, ind$upper * sda, lk, bl),
                            sem_effect_phrase(tot$estimate * sda, tot$lower * sda, tot$upper * sda, lk, bl, "corresponds to")))
    sw <- if (nongm) "standardised " else ""
    if (nrow(dir) && dir$estimate != 0 && sign(dir$estimate) != sign(ind$estimate) && abs(ind$estimate) > 0)
      s <- paste(s, sprintf("The %sdirect (%s) and indirect (%s) associations have opposite signs and partly cancel, so the total (%s) is smaller in magnitude than the stronger of the two; a weak total does not mean the driver has little influence.",
                            sw, if (nongm) fmt(dir$estimate, 2) else fmt_sig(dir$estimate), if (nongm) fmt(ind$estimate, 2) else fmt_sig(ind$estimate), if (nongm) fmt(tot$estimate, 2) else fmt_sig(tot$estimate)))
    else if (!nrow(dir)) s <- paste(s, sprintf("The model has no direct path from %s to %s, so the whole association is indirect.", al, bl))
    else if (abs(tot$estimate) > 0 && sign(tot$estimate) == sign(ind$estimate))
      s <- paste(s, sprintf("About %.0f%% of the %stotal association runs through the mediators (point estimate only).", 100 * min(1, ind$estimate / tot$estimate), sw))
    out <- c(out, s)
  }
  if (fit$engine == "lavaan" && length(fit$latent)) {
    ld <- est[grepl("^std_loading_", est$param), ]
    if (nrow(ld)) out <- c(out, sprintf("Measurement model: standardised loadings of indicators on their latent variable(s) range from %s to %s. Loadings above ~0.7 mean the indicator mostly reflects the latent variable; low loadings mean it is a weak measure.",
                                        fmt(min(ld$estimate), 2), fmt(max(ld$estimate), 2)))
  }
  if (fit$engine == "piecewise") {
    nong <- Filter(function(c) c$family != "gaussian", fit$comps)
    if (length(nong)) out <- c(out, sprintf("Detection versus ecology. %s %s analysed as recorded (%s). If not every individual or presence is detected, these responses mix abundance or occurrence with detectability, so a path into them can reflect an effect on detection (e.g. dense vegetation hiding animals) as well as on the animals; the SEM cannot separate the two. %s",
      paste(sem_lab(vapply(nong, `[[`, "", "response"), data), collapse = " and "), if (length(nong) > 1) "are" else "is",
      paste(vapply(nong, function(c) switch(c$family, poisson = "Poisson counts", nbinom = "negative-binomial counts", "0/1 outcomes"), ""), collapse = ", "),
      if (any(vapply(nong, function(c) length(c$groups) > 0, logical(1)))) "Effects in components with a random intercept are for a typical group (random effect at 0)." else ""))
  }
  if (length(fit$r2)) out <- c(out, paste0("Variance explained (R2): ", paste(sprintf("%s %.2f", sem_lab(names(fit$r2), data), fit$r2), collapse = "; "), "."))
  out <- c(out, "Causal reading. A structural equation model encodes a causal hypothesis, and the global test can reject it, but passing the test does not prove it: the coefficients are causal effects only if the diagram is right and no unmeasured variable drives two of the measured ones. Without that justification (an experiment, or strong prior knowledge stated in the methods), read each path as 'is associated with, holding the other parents constant'.")
  out
}

# ---- figures ------------------------------------------------------------------
# Publication suite (docs/FIGURES.md, row `sem`): path_diagram, indirect_effects,
# dsep (piecewise) / missing_paths + loadings (lavaan), partial_<x>_to_<y> for
# every numeric path, effects_panel. Every label goes through fig_label() (keys
# are the variable names, plus the keys listed in README "Figures"); covariates
# are drawn in the units of the data columns (the module never standardises;
# figures: units: {<var>: {center, scale}} back-transforms columns that the user
# supplied already standardised).

# Interval at level `ci` for rows of fit$estimates: the stored intervals are
# 95% (Monte Carlo / delta / bootstrap); any other level uses estimate +- z*SE.
sem_fig_ci <- function(d, ci) {
  if (abs(ci - 0.95) > 1e-9) {
    z <- stats::qnorm(1 - (1 - ci) / 2)
    d$lower <- d$estimate - z * d$sd; d$upper <- d$estimate + z * d$sd
  }
  d
}
sem_ci_note <- function(fit, ci) {
  meth <- switch(fit$settings$ci_method %||% "montecarlo", montecarlo = "Monte Carlo", delta = "delta-method",
                 bootstrap = "bootstrap percentile", "")
  if (abs(ci - 0.95) > 1e-9) sprintf("%d%% normal-approximation confidence", round(100 * ci))
  else sprintf("95%% %s confidence", meth)
}

# Semantic colours drawn from the active palette (so greyscale / custom
# palettes apply): positive path, negative path, flagged missing link.
sem_colours <- function(opts) {
  pal <- opts$palette %||% "okabe_ito"
  p <- if (length(pal) > 1) pal else eco_pal[[pal]]
  if (is.null(p)) p <- eco_pal$okabe_ito
  i <- if (length(p) >= 7) c(5, 6, 7) else c(1, min(3, length(p)), min(2, length(p)))
  stats::setNames(p[i], c("positive", "negative", "missing"))
}

# Text sizes in mm for geoms, from the figure base size (pt).
sem_txt <- function(opts, rel = 0.8) (opts$base_size %||% 9) * rel / ggplot2::.pt

# greyscale: the core constructors draw single-series lines/bands in the
# (green) accent colour whatever the palette; recolour them to ink/grey so a
# print-safe export carries no colour at all.
sem_grey <- function(p, opts) {
  if (!isTRUE(opts$greyscale) && !identical(opts$palette, "grey")) return(p)
  for (i in seq_along(p$layers)) {
    g <- class(p$layers[[i]]$geom)[1]
    if (g == "GeomRibbon" && !is.null(p$layers[[i]]$aes_params$fill)) p$layers[[i]]$aes_params$fill <- "grey45"
    if (g == "GeomLine" && !is.null(p$layers[[i]]$aes_params$colour)) p$layers[[i]]$aes_params$colour <- eco_pal$ink
  }
  p
}

sem_wrap <- function(x, width) vapply(x, function(z) paste(strwrap(z, width), collapse = "\n"), "", USE.NAMES = FALSE)

sem_layout <- function(edges, nodes) {
  lvl <- stats::setNames(rep(0, length(nodes)), nodes)
  for (it in seq_along(nodes)) for (k in seq_len(nrow(edges))) lvl[edges$to[k]] <- max(lvl[edges$to[k]], lvl[edges$from[k]] + 1)
  # exogenous variables sit one level before their earliest child, so their
  # arrows do not have to cross intermediate levels
  for (v in setdiff(nodes, edges$to)) {
    ch <- edges$to[edges$from == v]
    if (length(ch)) lvl[v] <- max(0, min(lvl[ch]) - 1)
  }
  d <- data.frame(name = nodes, x = 3 * lvl[nodes], stringsAsFactors = FALSE)
  d$y <- stats::ave(d$x, d$x, FUN = function(z) rev(seq_along(z) - (length(z) + 1) / 2))
  d$y <- d$y * 1.1
  d
}

sem_path_plot <- function(edges, nodes_df, missing = NULL, latent = character(), data = NULL, opts = NULL, ci = 0.95) {
  eco_require(c("ggraph", "igraph"))
  cols <- sem_colours(opts)
  nd <- nodes_df
  nd$label <- sem_wrap(sem_lab(nd$name, data, opts), 13)
  nd$type <- ifelse(nd$name %in% latent, "latent", "observed")
  lev_in <- fig_label("ci_excludes_0", sprintf("%d%% CI excludes 0", round(100 * ci)), opts)
  lev_out <- fig_label("ci_includes_0", sprintf("%d%% CI includes 0", round(100 * ci)), opts)
  lab_pos <- fig_label("sign_positive", "Positive", opts); lab_neg <- fig_label("sign_negative", "Negative", opts)
  e <- edges[, c("from", "to", "estimate", "lower", "upper")]
  e$kind <- "path"
  e$sign <- ifelse(e$estimate >= 0, "positive", "negative")
  e$support <- ifelse(e$lower > 0 | e$upper < 0, lev_in, lev_out)
  e$lab <- fmt(e$estimate, 2)
  # bend omitted-link arcs away from the middle of the diagram so they do not
  # cross the fitted paths
  ymid <- mean(range(nd$y))
  if (!is.null(missing) && nrow(missing)) {
    yy <- (nd$y[match(missing$from, nd$name)] + nd$y[match(missing$to, nd$name)]) / 2
    xs <- nd$x[match(missing$to, nd$name)] - nd$x[match(missing$from, nd$name)]
    # ggraph: positive strength bends to the right of the travel direction
    bend <- ifelse(yy <= ymid, -1, 1) * ifelse(xs >= 0, 1, -1)
    e <- rbind(cbind(e, bend = 0), data.frame(from = missing$from, to = missing$to, estimate = 0, lower = NA, upper = NA, kind = "missing",
                             sign = "missing", support = lev_out, bend = bend,
                             lab = paste0(fig_label("missing_link", "missing?", opts), " P = ", fmt_sig(missing$p_value, 2))))
  } else e$bend <- 0
  e$support <- factor(e$support, levels = c(lev_in, lev_out))
  nx <- stats::setNames(nd$x, nd$name); ny <- stats::setNames(nd$y, nd$name)
  e$skip <- e$kind == "path" & abs(nx[e$to] - nx[e$from]) > 3.5
  e$side <- "above"
  yb <- (ny[e$from] + ny[e$to]) / 2
  e$pbend <- ifelse(yb >= ymid, 1, -1) * ifelse(nx[e$to] >= nx[e$from], 1, -1)
  for (k in which(e$kind == "path" & !e$skip)) {
    oth <- e$from[e$kind == "path" & e$to == e$to[k] & e$from != e$from[k]]
    up_arc <- any(e$skip & e$from == e$from[k] & e$pbend * sign(nx[e$to] - nx[e$from]) > 0)
    if ((length(oth) && any(ny[oth] > ny[e$from[k]])) || up_arc) e$side <- replace(e$side, k, "below")
  }
  g <- igraph::graph_from_data_frame(e, vertices = nd[, c("name", "x", "y", "label", "type")], directed = TRUE)
  arw <- grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
  # caps sized from the node-label font (label_rect defaults to 12 pt, which left
  # arrows stopping well short of the smaller node boxes)
  node_pt <- sem_txt(opts, 0.85) * ggplot2::.pt
  cap <- function(lbl) ggraph::label_rect(lbl, padding = ggplot2::margin(1.6, 1.6, 1.6, 1.6, "mm"), fontsize = node_pt, lineheight = 0.9)
  p <- ggraph::ggraph(g, layout = "manual", x = nd$x, y = nd$y) +
    ggraph::geom_edge_link(ggplot2::aes(edge_width = abs(estimate), edge_colour = sign, edge_linetype = support,
                                        start_cap = cap(node1.label), end_cap = cap(node2.label),
                                        filter = kind == "path" & !skip),
                           arrow = arw, show.legend = c(edge_width = FALSE, edge_colour = TRUE, edge_linetype = TRUE))
  # coefficient labels. Arrows that skip a level: white box on the arrow a
  # quarter of the way along (their middle passes close to intermediate
  # nodes). Arrows between adjacent levels: beside the arrow, on the side away
  # from the other arrows entering the same node (below when another parent
  # sits higher), so labels never sit on an arrowhead or another arrow.
  for (sd in c("above", "below"))
    p <- p + ggraph::geom_edge_link(ggplot2::aes(label = lab, start_cap = cap(node1.label), end_cap = cap(node2.label),
                                                 filter = kind == "path" & !skip & side == !!sd),
                                    edge_alpha = 0, angle_calc = "along", label_size = sem_txt(opts, 0.8), label_colour = eco_pal$ink,
                                    label_dodge = grid::unit(if (sd == "above") 3.8 else -3.8, "mm"), show.legend = FALSE)
  # arrows that skip a level are drawn as arcs bent away from the middle of
  # the diagram, so they do not run through intermediate nodes
  for (sg in c(1, -1)) if (any(e$skip & e$pbend == sg))
    p <- p + ggraph::geom_edge_arc(ggplot2::aes(edge_width = abs(estimate), edge_colour = sign, edge_linetype = support,
                                                label = lab, start_cap = cap(node1.label), end_cap = cap(node2.label),
                                                filter = skip & pbend == !!sg),
                                   strength = 0.3 * sg, arrow = arw, angle_calc = "along", label_size = sem_txt(opts, 0.8),
                                   label_colour = eco_pal$ink, label_dodge = grid::unit(3.4 * (if (sg > 0) -1 else 1), "mm"),
                                   show.legend = FALSE)
  for (b in setdiff(unique(e$bend[e$kind == "missing"]), NA)) {
    p <- p + ggraph::geom_edge_arc(ggplot2::aes(label = lab, start_cap = cap(node1.label), end_cap = cap(node2.label),
                                                filter = kind == "missing" & bend == !!b), strength = 0.18 * b, edge_colour = cols[["missing"]],
                                   edge_linetype = "dotted", edge_width = 0.6, arrow = grid::arrow(length = grid::unit(1.8, "mm"), type = "open"),
                                   angle_calc = "along", label_dodge = grid::unit(3.4 * b, "mm"), label_size = sem_txt(opts, 0.75),
                                   label_colour = cols[["missing"]])
  }
  p <- p + ggraph::geom_node_label(ggplot2::aes(label = label, fill = type), size = sem_txt(opts, 0.85), lineheight = 0.9,
                                   colour = eco_pal$ink, label.padding = grid::unit(0.3, "lines"), show.legend = FALSE) +  # latent boxes are shaded; explained in the caption
    ggraph::scale_edge_width(range = c(0.4, 2.2), guide = "none") +
    ggraph::scale_edge_colour_manual(values = c(positive = cols[["positive"]], negative = cols[["negative"]], missing = cols[["missing"]]),
                                     breaks = c("positive", "negative"), labels = c(lab_pos, lab_neg), name = NULL) +
    ggraph::scale_edge_linetype_manual(values = stats::setNames(c("solid", "dashed"), c(lev_in, lev_out)), name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = c(observed = "white", latent = "#F0E9D2"), name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 1.1)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = if (any(e$kind == "missing") || any(e$skip)) 1 else 0.6)) +
    ggplot2::guides(x = "none", y = "none", edge_colour = ggplot2::guide_legend(order = 1, override.aes = list(edge_width = 1)),
                    edge_linetype = ggplot2::guide_legend(order = 2, override.aes = list(edge_width = 0.6, edge_colour = eco_pal$ink),
                                                          theme = ggplot2::theme(legend.key.width = grid::unit(10, "mm")))) +
    ggplot2::labs(x = NULL, y = NULL)
  nlev <- length(unique(nd$x))
  attr(p, "fig_size") <- list(width = if (nlev >= 4) "double" else "onehalf",
                              height = 40 + 22 * max(table(nd$x)))
  p
}

sem_path_diagram <- function(fit, data, opts, ci = 0.95) {
  est <- fit$estimates
  z <- stats::qnorm(1 - (1 - ci) / 2)
  if (fit$engine == "piecewise") {
    t <- sem_numeric_paths(fit)
    s <- sem_fig_ci(est[match(paste0("std_", t$from, "_to_", t$to), est$param), ], ci)
    edges <- data.frame(from = t$from, to = t$to, estimate = s$estimate, lower = s$lower, upper = s$upper)
    fac <- fit$table[is.na(fit$table$std_factor), ]
    if (nrow(fac)) {
      ff <- unique(fac[, c("from", "to")])
      for (k in seq_len(nrow(ff))) {
        zz <- fac[fac$from == ff$from[k] & fac$to == ff$to[k], ]
        edges <- rbind(edges, data.frame(from = ff$from[k], to = ff$to[k], estimate = NA, lower = min(zz$lower), upper = max(zz$upper)))
      }
    }
    nodes <- unique(c(edges$from, edges$to))
    miss <- if (!is.null(fit$dsep$table)) { m <- fit$dsep$table[fit$dsep$table$p_value < 0.05, ]; data.frame(from = m$missing_from, to = m$missing_to, p_value = m$p_value) }
    latent <- character()
  } else {
    t <- fit$table[fit$table$op == "~" & fit$table$group == 1, ]
    edges <- data.frame(from = t$from, to = t$to, estimate = t$std, lower = t$std - z * t$std_se, upper = t$std + z * t$std_se)
    nodes <- unique(c(edges$from, edges$to)); latent <- fit$latent
    mi <- fit$modindices
    miss <- if (!is.null(mi) && nrow(mi)) { m <- mi[mi$op == "~", , drop = FALSE]; m <- m[m$lhs %in% nodes & m$rhs %in% nodes, ]
      if (nrow(m)) data.frame(from = m$rhs, to = m$lhs, p_value = stats::pchisq(m$mi, 1, lower.tail = FALSE)) }
  }
  edges$estimate[is.na(edges$estimate)] <- 0
  lay <- sem_layout(edges, nodes)
  sem_path_plot(edges, lay, miss, latent, data, opts, ci)
}

# Conditional effect of x on comp$response (other parents at their means /
# modal level) with partial residuals. Gaussian: y - fitted + fitted(x only).
# Count (log link): working residuals (y - mu)/mu added to the conditional
# linear predictor and back-transformed (visreg-style partial residuals).
# Binary: no residual points (a rug of observed x).
sem_effect_pred <- function(fit, comp, x, n = 100, ci = 0.95) {
  m <- fit$object[[comp$response]]; df <- fit$df
  ff <- lme4::nobars(comp$formula)
  tt <- stats::delete.response(stats::terms(ff))
  nd <- df[rep(1, n), , drop = FALSE]
  for (v in comp$preds) nd[[v]] <- if (is.numeric(df[[v]])) mean(df[[v]]) else factor(names(which.max(table(df[[v]]))), levels = levels(factor(df[[v]])))
  nd[[x]] <- seq(min(df[[x]]), max(df[[x]]), length.out = n)
  cv <- sem_coef_vcov(fit$fos[[comp$response]])
  xl <- lapply(df[vapply(df, is.factor, logical(1))], levels)
  X <- stats::model.matrix(tt, nd, xlev = xl)
  X <- X[, names(cv$b), drop = FALSE]
  eta <- as.numeric(X %*% cv$b); se <- sqrt(pmax(0, rowSums((X %*% cv$V) * X)))
  z <- stats::qnorm(1 - (1 - ci) / 2)
  inv <- switch(comp$link, identity = identity, log = exp, logit = stats::plogis, probit = stats::pnorm)
  Xo <- stats::model.matrix(tt, df, xlev = xl)[, names(cv$b), drop = FALSE]
  Xm <- Xo; for (v in comp$preds) if (v != x) { cols <- grep(paste0("^", v), colnames(Xm)); Xm[, cols] <- matrix(X[1, cols], nrow(Xm), length(cols), byrow = TRUE) }
  y <- df[[comp$response]]
  obs <- NULL
  if (comp$family == "gaussian") {
    obs <- data.frame(x = df[[x]], y = y - as.numeric(Xo %*% cv$b) + as.numeric(Xm %*% cv$b))
  } else if (comp$family != "binomial" && comp$link == "log") {
    mu <- tryCatch(as.numeric(stats::fitted(m)), error = function(e) NULL)
    if (is.null(mu) || length(mu) != length(y)) mu <- exp(as.numeric(Xo %*% cv$b))
    obs <- data.frame(x = df[[x]], y = exp(as.numeric(Xm %*% cv$b) + (y - mu) / mu))
  } else obs <- data.frame(x = df[[x]])
  list(pred = data.frame(x = nd[[x]], estimate = inv(eta), lower = inv(eta - z * se), upper = inv(eta + z * se)), obs = obs)
}

# lavaan: partial effect of each regression path (single-group models). Latent
# variables enter as factor scores (lavPredict, units of the first indicator).
# Line through the means with slope b; band from SE(b) and the residual SE of
# the mean; points are partial residuals y - sum_{k != x} b_k (v_k - mean v_k).
sem_lavaan_partials <- function(fit, ci = 0.95) {
  m <- fit$object
  if (length(fit$group_labels %||% "") > 1) return(list())
  pe <- lavaan::parameterEstimates(m)
  V <- fit$df
  if (length(fit$latent)) {
    sc <- tryCatch(lavaan::lavPredict(m, type = "lv"), error = function(e) NULL)
    if (is.null(sc) || nrow(sc) != nrow(V)) return(list())
    V <- cbind(V, as.data.frame(sc))
  }
  z <- stats::qnorm(1 - (1 - ci) / 2)
  out <- list()
  reg <- pe[pe$op == "~", ]
  for (yv in unique(reg$lhs)) {
    r <- reg[reg$lhs == yv, ]
    if (!all(c(yv, r$rhs) %in% names(V))) next
    rv <- pe$est[pe$op == "~~" & pe$lhs == yv & pe$rhs == yv][1]
    for (k in seq_len(nrow(r))) {
      xv <- r$rhs[k]; b <- r$est[k]; sb <- r$se[k]
      xs <- seq(min(V[[xv]]), max(V[[xv]]), length.out = 100)
      my <- mean(V[[yv]]); mx <- mean(V[[xv]])
      se <- sqrt(rv / nrow(V) + (xs - mx)^2 * sb^2)
      fitl <- my + b * (xs - mx)
      others <- setdiff(r$rhs, xv)
      adj <- V[[yv]]
      for (o in others) adj <- adj - r$est[r$rhs == o] * (V[[o]] - mean(V[[o]]))
      out[[paste0(xv, "_to_", yv)]] <- list(x = xv, y = yv, latent_x = xv %in% fit$latent, latent_y = yv %in% fit$latent,
        pred = data.frame(x = xs, estimate = fitl, lower = fitl - z * se, upper = fitl + z * se),
        obs = data.frame(x = V[[xv]], y = adj))
    }
  }
  out
}

sem_figures <- function(fit, data, opts = fig_opts()) {
  figs <- list()
  ci <- opts$ci %||% 0.95
  pct <- round(100 * ci)
  scal <- data$scaling
  cols <- sem_colours(opts)
  figs$path_diagram <- with_caption(sem_path_diagram(fit, data, opts, ci),
    paste0("Path diagram with standardised coefficients", if (fit$engine == "lavaan") " (lavaan std.all)" else " (piecewiseSEM conventions: SD scaling; latent-theoretic for binary, observation-empirical for count responses)",
           sprintf(". Arrow width is proportional to |coefficient|; colour gives the sign; dashed arrows have a %d%% confidence interval that includes zero.", pct),
           if (fit$engine == "piecewise") " Dotted arrows (labelled with P) are links the model omits but the d-separation test flags (P < 0.05)." else " Dotted arrows (labelled with P) are omitted regressions with modification index > 3.84.",
           if (fit$engine == "lavaan" && length(fit$latent %||% character()) > 0) " Shaded boxes are latent variables (measured through their indicators)." else "",
           if (fit$engine == "lavaan" && length(fit$group_labels) > 1) sprintf(" Shown for group %s.", fit$group_labels[1]) else ""))
  est <- fit$estimates
  eff <- est[grepl("^std_(direct|indirect|total)_", est$param), ]
  if (nrow(eff)) {
    eff <- sem_fig_ci(eff, ci)
    eff$type <- sub("^std_([a-z]+)_.*", "\\1", eff$param)
    pr <- sub("^std_[a-z]+_", "", eff$param)
    eff$param <- vapply(pr, function(p) { i <- regexpr("_to_", p)
      sem_wrap(paste(sem_lab(substr(p, 1, i - 1), data, opts), "→", sem_lab(substr(p, i + 4, nchar(p)), data, opts)), 24) }, "")
    eff$param <- factor(eff$param, levels = rev(unique(eff$param)))
    tl <- c(direct = fig_label("effect_direct", "Direct", opts), indirect = fig_label("effect_indirect", "Indirect", opts),
            total = fig_label("effect_total", "Total", opts))
    eff$group <- factor(tl[eff$type], levels = tl)
    p <- fig_coef(eff[, c("param", "estimate", "lower", "upper", "group")], null_value = 0, order = FALSE,
                  xlab = fig_label("axis_std_effect", sprintf("Standardised effect (%d%% CI)", pct), opts)) +
      ggplot2::aes(shape = group) + ggplot2::scale_shape_manual(values = c(16, 17, 15), name = NULL) +
      ggplot2::scale_x_continuous(n.breaks = 5) +
      ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE), shape = ggplot2::guide_legend(reverse = TRUE))
    # legend sits to the right (core default), so give the panel room
    attr(p, "fig_size") <- list(width = "onehalf", height = 34 + 9 * length(unique(eff$param)) * 1.4)
    figs$indirect_effects <- with_caption(p, sprintf("Direct, indirect (summed over all mediated routes) and total standardised effects (points) with %s intervals (lines); dashed line = no effect. Pairs with no direct path show only indirect and total (equal) effects.",
                                                     sem_ci_note(fit, ci)))
  }
  if (fit$engine == "piecewise" && !is.null(fit$dsep$table) && nrow(fit$dsep$table)) {
    d <- fit$dsep$table
    d$claim <- sem_wrap(paste(sem_lab(d$missing_from, data, opts), "→", sem_lab(d$missing_to, data, opts)), 24)
    d$claim <- stats::reorder(factor(d$claim), -d$p_value)
    l_flag <- fig_label("dsep_flagged", "P < 0.05 (link needed?)", opts); l_ok <- fig_label("dsep_consistent", "Consistent with model", opts)
    d$flag <- factor(ifelse(d$p_value < 0.05, l_flag, l_ok), levels = c(l_flag, l_ok))
    p <- ggplot2::ggplot(d, ggplot2::aes(p_value, claim, colour = flag, shape = flag)) +
      ggplot2::geom_vline(xintercept = 0.05, linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_point(size = 2) + ggplot2::scale_x_log10(limits = c(min(1e-3, min(d$p_value)), 1), breaks = c(1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.2, 1),
                                                                   labels = function(b) format(b, scientific = FALSE, drop0trailing = TRUE)) +
      ggplot2::scale_colour_manual(values = stats::setNames(c(cols[["missing"]], eco_pal$null), c(l_flag, l_ok)), name = NULL, drop = FALSE) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(17, 16), c(l_flag, l_ok)), name = NULL, drop = FALSE) +
      ggplot2::guides(colour = "none", shape = "none") +   # meaning given in the caption; a side legend squeezes the panel
      ggplot2::labs(x = fig_label("axis_dsep_p", "P value (log scale)", opts),
                    y = fig_label("axis_dsep_claim", "Omitted link tested", opts))
    attr(p, "fig_size") <- list(width = "single", height = 38 + 11 * nrow(d))
    figs$dsep <- with_caption(p, sprintf("d-separation tests (Shipley 2009): each row tests whether a link the model omits is needed, conditioning on both variables' parents (points = P values; highlighted triangles flag P < 0.05, circles are consistent with the model). Fisher's C = %.2f, df = %d, P = %.3f. Dashed line P = 0.05.",
                                         fit$dsep$C, as.integer(fit$dsep$df), fit$dsep$p))
  }
  if (fit$engine == "lavaan") {
    mi <- fit$modindices
    if (!is.null(mi) && nrow(mi)) {
      d <- utils::head(mi, 10)
      lw <- function(v) sem_lab(v, data, opts)
      d$what <- sem_wrap(ifelse(d$op == "~", paste(lw(d$rhs), "→", lw(d$lhs)),
                         ifelse(d$op == "~~", paste(lw(d$lhs), "↔", lw(d$rhs)), paste(lw(d$lhs), "→", lw(d$rhs), "(loading)"))), 28)
      d$what <- stats::reorder(factor(d$what), d$mi)
      p <- ggplot2::ggplot(d, ggplot2::aes(mi, what)) + ggplot2::geom_col(fill = "#bfcfc7", colour = eco_pal$muted, linewidth = 0.2) +
        ggplot2::geom_vline(xintercept = 3.84, linetype = 2, colour = eco_pal$null) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(x = fig_label("axis_mi", "Modification index (expected chi-square drop)", opts), y = NULL)
      attr(p, "fig_size") <- list(height = 25 + 8 * nrow(d))
      figs$missing_paths <- with_caption(p, "Largest modification indices: parameters the model fixes at zero whose addition would most improve fit (dashed line 3.84 = P 0.05). Add only with a biological reason.")
    }
    lt <- fit$table[fit$table$op == "=~" & fit$table$group == 1, ]
    if (nrow(lt)) {
      z <- stats::qnorm(1 - (1 - ci) / 2)
      ld <- data.frame(param = sem_wrap(paste(sem_lab(lt$from, data, opts), "→", sem_lab(lt$to, data, opts)), 30),
                       estimate = lt$std, lower = lt$std - z * lt$std_se, upper = lt$std + z * lt$std_se)
      p <- fig_coef(ld, null_value = 0, xlab = fig_label("axis_std_loading", sprintf("Standardised loading (%d%% CI)", pct), opts)) +
        ggplot2::geom_vline(xintercept = 0.7, linetype = 3, colour = eco_pal$null, linewidth = 0.4) + ggplot2::expand_limits(x = 1)
      figs$loadings <- with_caption(p, sprintf("Standardised factor loadings (latent variable → indicator; points) with %d%% confidence intervals (lines). Dotted line 0.7: values above it indicate a strong indicator; dashed line = 0.", pct))
    }
  }
  eps <- list()
  if (fit$engine == "piecewise") {
    for (c in fit$comps) for (x in c$preds) {
      if (!is.numeric(fit$df[[x]])) next
      pp <- tryCatch(sem_effect_pred(fit, c, x, ci = ci), error = function(e) NULL)
      if (is.null(pp)) next
      pp$pred$x <- fig_x(pp$pred$x, x, opts, scal); pp$obs$x <- fig_x(pp$obs$x, x, opts, scal)
      ylab <- sem_lab(c$response, data, opts)
      # "Expected <label>" rather than "<label> (expected)": labels often end in a
      # unit in parentheses already ("Bird detections (per point count)")
      exp_lab <- if (grepl("^[A-Z][a-z]", ylab)) paste0(tolower(substr(ylab, 1, 1)), substring(ylab, 2)) else ylab
      ylab2 <- switch(c$family, gaussian = ylab, binomial = paste0("P(", ylab, ")"),
                      paste(fig_label("expected_prefix", "Expected", opts), exp_lab))
      p <- fig_cond_effect(pp$pred, sem_wrap(sem_lab(x, data, opts), 34), sem_wrap(ylab2, 24), obs = if (isTRUE(opts$show_data)) pp$obs else NULL,
                           ylim01 = c$family == "binomial")
      if (c$family %in% c("poisson", "nbinom")) p <- p + ggplot2::expand_limits(y = 0)
      p <- sem_grey(p, opts)
      nm <- paste0("partial_", x, "_to_", c$response)
      eps[[nm]] <- p
      figs[[nm]] <- with_caption(p, sprintf("Conditional effect of %s on %s (response scale; line) with %d%% confidence band; other predictors of %s held at their means%s. %s",
                                            sem_lab(x, data, opts), sem_lab(c$response, data, opts), pct, sem_lab(c$response, data, opts),
                                            if (length(c$groups)) ", random effects at zero" else "",
                                            switch(c$family, gaussian = "Points are partial residuals (observed value minus the fitted contribution of the other predictors).",
                                                   binomial = "Rug shows observed values of the predictor.",
                                                   "Points are partial residuals: working residuals added to the conditional linear predictor and back-transformed to counts.")))
    }
  } else {
    lp <- tryCatch(sem_lavaan_partials(fit, ci), error = function(e) list())
    for (k in names(lp)) {
      q <- lp[[k]]
      q$pred$x <- fig_x(q$pred$x, q$x, opts, scal); q$obs$x <- fig_x(q$obs$x, q$x, opts, scal)
      xl <- sem_lab(q$x, data, opts); yl <- sem_lab(q$y, data, opts)
      if (q$latent_x) xl <- paste(xl, fig_label("factor_score", "(factor score)", opts))
      if (q$latent_y) yl <- paste(yl, fig_label("factor_score", "(factor score)", opts))
      p <- fig_cond_effect(q$pred, sem_wrap(xl, 34), sem_wrap(yl, 24), obs = if (isTRUE(opts$show_data)) q$obs else NULL)
      p <- sem_grey(p, opts)
      nm <- paste0("partial_", k)
      eps[[nm]] <- p
      figs[[nm]] <- with_caption(p, sprintf("Partial effect of %s on %s (line, unstandardised lavaan coefficient through the means) with %d%% confidence band; points are partial residuals (other predictors' contributions removed).%s",
                                            sem_lab(q$x, data, opts), sem_lab(q$y, data, opts), pct,
                                            if (q$latent_x || q$latent_y) " Latent variables are shown as factor scores (lavPredict), centred, in the units of their first indicator; factor scores are shrunk towards the mean, so the points are less spread than the true latent values." else ""))
    }
  }
  if (length(eps) > 1) {
    # theme each panel: a theme added to the patchwork later reaches only the last panel
    eps <- lapply(eps, function(p) p + theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic", opts$legend %||% "right"))
    figs$effects_panel <- with_caption(compose_figures(eps, ncol = 2), "All path-level conditional effects on the response scale (panels as in the individual partial_ figures; same bands and points).")
    attr(figs$effects_panel, "fig_size") <- list(width = "double", height = 55 * ceiling(length(eps) / 2))
  }
  figs
}

# ---- methods ------------------------------------------------------------------
sem_methods <- function(fit, data, settings) {
  ci <- switch(settings$ci_method %||% "montecarlo",
               montecarlo = sprintf("Monte Carlo confidence intervals (Preacher & Selig 2012; %d draws from the estimated sampling distribution of the path coefficients)", settings$n_mc %||% 4000),
               delta = "delta-method (first-order Taylor) standard errors",
               bootstrap = sprintf("percentile intervals from %d nonparametric bootstrap resamples", settings$n_boot %||% 500))
  if (fit$engine == "piecewise") {
    comps <- vapply(fit$comps, function(c) sprintf("%s (%s, %s link%s)", c$text, c$family, c$link, if (length(c$groups)) ", random intercept" else ""), "")
    sprintf(paste(
      "We fitted a piecewise structural equation model (Shipley 2009; Lefcheck 2016) consisting of %d component models: %s.",
      "Gaussian components were fitted by least squares (lm) or REML (lme4::lmer), Poisson and binomial components by maximum likelihood (glm / lme4::glmer), and negative-binomial components with MASS::glm.nb or lme4::glmer.nb.",
      "The overall model was evaluated with Shipley's test of directed separation: each independence claim of the basis set was tested with the partial coefficient of the omitted link in the model of the causally later variable, conditioning on the parents of both variables (%s), and the P values were combined as Fisher's C = -2 sum ln(P), compared with a chi-square distribution on 2k degrees of freedom (piecewiseSEM %s).",
      "Standardised coefficients scaled each coefficient by SD(x)/SD(y), using the latent-theoretic SD of the linear predictor plus the link variance for binary responses and the observation-empirical SD for count responses (Grace et al. 2018).",
      "Indirect effects were computed as the sum over all directed routes of the products of the path coefficients and total effects as direct plus indirect; their 95%% intervals are %s. Analyses were run in R %s with ecotools module sem."),
      length(comps), paste(comps, collapse = "; "),
      if (identical(fit$dsep$direction, "conserve")) "for claims involving non-Gaussian variables both directions were tested and the smaller P value kept, piecewiseSEM conserve = TRUE" else "Shipley 2009",
      as.character(utils::packageVersion("piecewiseSEM")), ci, getRversion())
  } else {
    sprintf(paste(
      "We fitted a structural equation model in lavaan %s (Rosseel 2012) by %s%s.",
      "Global fit was assessed with the model chi-square test, CFI, RMSEA (with 90%% CI) and SRMR (Hu & Bentler 1999), and modification indices were inspected for omitted links.",
      "Standardised coefficients are lavaan's std.all solution. Indirect effects were computed as the sum over directed routes of products of regression coefficients, and total effects as direct plus indirect, with 95%% intervals from %s. Latent variables were scaled by fixing the first loading to 1.",
      "Analyses were run in R %s with ecotools module sem."),
      as.character(utils::packageVersion("lavaan")), lavaan::lavInspect(fit$object, "options")$estimator,
      if (length(fit$group_labels) > 1) paste0(" as a multigroup model across ", paste(fit$group_labels, collapse = ", ")) else "",
      ci, getRversion())
  }
}
