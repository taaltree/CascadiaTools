# lavaan engine (Rosseel 2012): global-covariance SEM for latent variables,
# multigroup models and models with correlated errors. Gaussian (ML / robust
# MLR) estimation. Indirect and total effects are computed from the
# regression paths with the same graph code as the piecewise engine, using
# the joint sampling covariance of all free parameters (Monte Carlo or delta
# method) or lavaan's own bootstrap draws.

sem_lavaan_edges_syntax <- function(syntax) {
  pt <- lavaan::lavaanify(syntax, auto = FALSE)
  r <- pt[pt$op == "~", ]; l <- pt[pt$op == "=~", ]
  e <- rbind(data.frame(from = r$rhs, to = r$lhs, op = rep("~", nrow(r)), stringsAsFactors = FALSE),
             data.frame(from = l$lhs, to = l$rhs, op = rep("=~", nrow(l)), stringsAsFactors = FALSE))
  e <- unique(e)
  e$key <- paste0(e$from, "_to_", e$to)
  e
}

sem_fit_lavaan <- function(data, st) {
  sp <- data$spec
  syntax <- sp$lavaan
  if (is.null(syntax) || !nzchar(syntax)) stop("engine lavaan needs `lavaan:` model syntax in the config")
  df <- data$df
  for (v in names(df)) if (is.character(df[[v]])) df[[v]] <- factor(df[[v]])
  boot <- identical(st$ci_method, "bootstrap")
  args <- list(model = syntax, data = df, estimator = sp$estimator %||% st$estimator %||% "ML",
               missing = st$missing %||% "listwise")
  if (!is.null(sp$group)) args$group <- sp$group
  if (boot) { args$se <- "bootstrap"; args$bootstrap <- st$n_boot %||% 500 }
  w <- character()
  set.seed(st$seed %||% 1)
  fit <- withCallingHandlers(do.call(lavaan::sem, args),
                             warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  conv <- isTRUE(lavaan::lavInspect(fit, "converged"))
  postok <- tryCatch(suppressWarnings(lavaan::lavInspect(fit, "post.check")), error = function(e) FALSE)
  pe <- lavaan::parameterEstimates(fit, ci = TRUE, level = 0.95, standardized = FALSE,
                                   boot.ci.type = if (boot) "perc" else "norm")
  ss <- lavaan::standardizedSolution(fit, type = "std.all", level = 0.95)
  pt <- lavaan::parTable(fit)
  ngroups <- lavaan::lavInspect(fit, "ngroups")
  glabs <- if (ngroups > 1) lavaan::lavInspect(fit, "group.label") else ""
  gsuf <- function(g) if (ngroups > 1) paste0("_", make.names(glabs[g])) else ""
  if (is.null(pe$group)) pe$group <- 1L
  if (is.null(ss$group)) ss$group <- 1L
  if (is.null(pt$group)) pt$group <- 1L
  lv <- lavaan::lavNames(fit, "lv")
  ind <- unique(pt$rhs[pt$op == "=~"])
  rows <- list()
  for (i in seq_len(nrow(pe))) {
    r <- pe[i, ]
    nm <- switch(r$op, "~" = paste0("path_", r$rhs, "_to_", r$lhs), "=~" = paste0("loading_", r$lhs, "_", r$rhs),
                 "~~" = if (r$lhs == r$rhs) paste0(if (r$lhs %in% lv) "var_" else "resvar_", r$lhs) else paste0("cov_", r$lhs, "_", r$rhs),
                 ":=" = r$lhs, "~1" = paste0("intercept_", r$lhs), NULL)
    if (is.null(nm)) next
    if (r$op == "~~" && r$lhs == r$rhs && !isTRUE(st$report_variances)) next
    if (r$op == "~1") next
    s <- ss[ss$lhs == r$lhs & ss$op == r$op & ss$rhs == r$rhs & ss$group == r$group, ]
    kind <- switch(r$op, "~" = "path (unstandardised)", "=~" = "loading (unstandardised)", "~~" = "covariance", ":=" = "defined parameter")
    rows[[length(rows) + 1]] <- data.frame(param = paste0(nm, gsuf(r$group)), estimate = r$est, sd = r$se,
                                           lower = r$ci.lower, upper = r$ci.upper, kind = kind,
                                           p_value = r$pvalue, from = if (r$op == "~") r$rhs else r$lhs,
                                           to = if (r$op == "~") r$lhs else r$rhs, op = r$op, group = r$group,
                                           std = if (nrow(s)) s$est.std[1] else NA_real_,
                                           std_se = if (nrow(s)) s$se[1] else NA_real_, stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  tab$std_factor <- ifelse(abs(tab$estimate) > 1e-12, tab$std / tab$estimate, NA_real_)
  # standardised rows (std.all with lavaan delta-method SE)
  sr <- tab[tab$op %in% c("~", "=~") & !is.na(tab$std), ]
  z <- stats::qnorm(0.975)
  stdrows <- data.frame(param = sub("^(path|loading)_", "std_\\1_", sr$param), estimate = sr$std, sd = sr$std_se,
                        lower = sr$std - z * sr$std_se, upper = sr$std + z * sr$std_se,
                        kind = sub("unstandardised", "standardised, std.all", sr$kind), p_value = NA_real_,
                        stringsAsFactors = FALSE)
  stdrows$param <- sub("^std_path_", "std_", stdrows$param)
  # effects
  theta_hat <- lavaan::coef(fit)
  V <- tryCatch(lavaan::lavInspect(fit, "vcov"), error = function(e) NULL)
  method <- st$ci_method %||% "montecarlo"
  draws_all <- NULL
  if (method == "montecarlo" && !is.null(V)) { set.seed(st$seed %||% 1); draws_all <- sem_rmvn(st$n_mc %||% 4000, theta_hat, V) }
  if (method == "bootstrap") draws_all <- lavaan::lavInspect(fit, "boot")
  eff <- list(); effinfo <- list()
  for (g in seq_len(ngroups)) {
    eg <- pt[pt$op %in% c("~", "=~") & pt$group == g, ]
    edges <- data.frame(from = ifelse(eg$op == "~", eg$rhs, eg$lhs), to = ifelse(eg$op == "~", eg$lhs, eg$rhs),
                        free = eg$free, est = eg$est, op = eg$op, stringsAsFactors = FALSE)
    edges$key <- paste0(edges$from, "_to_", edges$to)
    struct <- edges[!(edges$op == "=~"), ]   # effects only through structural paths
    # treat latent->indicator edges as not part of effects
    pairs <- sem_pairs_from_spec(sp$effects, struct)
    if (!nrow(pairs) || is.null(V) && method != "bootstrap") next
    keep <- pairs$from %in% c(struct$from, struct$to) & pairs$to %in% c(struct$from, struct$to)
    pairs <- pairs[keep, , drop = FALSE]
    if (!nrow(pairs)) next
    trow <- tab[tab$op == "~" & tab$group == g, ]
    fac <- stats::setNames(trow$std_factor, paste0(trow$from, "_to_", trow$to))[struct$key]
    fac[is.na(fac)] <- 0
    fn <- sem_effect_fun(struct, pairs); fns <- sem_effect_fun(struct, pairs, "std_")
    map <- function(th_all) { v <- struct$est; f <- struct$free > 0; v[f] <- th_all[struct$free[f]]; stats::setNames(v, struct$key) }
    gu <- function(th_all) fn(map(th_all)); gsd <- function(th_all) fns(map(th_all) * fac)
    e1 <- sem_ci_rows(gu, theta_hat, if (method == "delta") "delta" else "montecarlo", V, draws_all)
    e2 <- sem_ci_rows(gsd, theta_hat, if (method == "delta") "delta" else "montecarlo", V, draws_all)
    e1$kind <- paste0(sub("_.*", "", e1$param), " effect (unstandardised)")
    e2$kind <- paste0(sub("^std_([a-z]+)_.*", "\\1", e2$param), " effect (standardised)")
    ee <- rbind(e1, e2)
    for (s in attr(fn, "specs")) if (!any(s$direct)) ee <- ee[!ee$param %in% paste0(c("direct_", "std_direct_"), s$from, "_to_", s$to), ]
    ee$param <- paste0(ee$param, gsuf(g)); ee$p_value <- NA_real_
    eff[[g]] <- ee
    effinfo[[g]] <- do.call(rbind, lapply(attr(fn, "specs"), function(s)
      data.frame(from = s$from, to = s$to, n_paths = length(s$paths), n_indirect = sum(!s$direct),
                 routes = paste(vapply(s$paths, paste, "", collapse = " -> "), collapse = "; "),
                 mediators = paste(sem_mediators(s), collapse = ", "), nongaussian_mediators = "",
                 sink_family = "gaussian", sink_link = "identity", group = g, stringsAsFactors = FALSE)))
  }
  keep <- c("param", "estimate", "sd", "lower", "upper", "kind", "p_value")
  est <- rbind(tab[, keep], stdrows[, keep], do.call(rbind, eff))
  rownames(est) <- NULL
  fm <- tryCatch(lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower",
                                            "rmsea.ci.upper", "srmr", "aic", "bic")), error = function(e) NULL)
  if (!is.null(fm) && grepl("MLR|MLM", args$estimator)) {
    fr <- tryCatch(lavaan::fitMeasures(fit, c("chisq.scaled", "pvalue.scaled", "cfi.robust", "rmsea.robust")), error = function(e) NULL)
    fm <- c(fm, fr)
  }
  mi <- tryCatch(lavaan::modindices(fit, sort. = TRUE, minimum.value = 3.84), error = function(e) NULL)
  if (!is.null(mi) && nrow(mi)) {
    mi <- mi[mi$op %in% c("~", "~~", "=~"), , drop = FALSE]
    mi <- mi[!(mi$op == "~~" & mi$lhs == mi$rhs), , drop = FALSE]
  }
  notes <- sem_lavaan_notes(fit, fm, mi, conv, postok, w, df, ngroups, glabs)
  list(estimates = est,
       diagnostics = list(converged = conv && isTRUE(postok), notes = notes, engine = "lavaan",
                          fit_measures = fm,
                          model_checks_ok = if (is.null(fm) || !is.finite(fm[["pvalue"]]) || fm[["df"]] == 0) NA else fm[["pvalue"]] >= 0.05),
       object = fit, df = df, table = tab, fit_measures = fm, modindices = mi,
       effects_info = if (length(effinfo)) do.call(rbind, effinfo) else NULL,
       latent = lv, indicators = ind, group_labels = glabs, settings = st, engine = "lavaan")
}

sem_lavaan_notes <- function(fit, fm, mi, conv, postok, w, df, ngroups, glabs) {
  n <- character()
  n <- c(n, sprintf("Engine lavaan %s, estimator %s; %d observations%s.", as.character(utils::packageVersion("lavaan")),
                    lavaan::lavInspect(fit, "options")$estimator, sum(lavaan::lavInspect(fit, "nobs")),
                    if (ngroups > 1) paste0(" in ", ngroups, " groups (", paste(glabs, collapse = ", "), ")") else ""))
  if (!conv) n <- c(n, "NOT CONVERGED: the optimiser did not converge. Check the model is identified (each latent variable needs >= 3 indicators, or 2 with constraints) and variables are on similar scales.")
  if (!isTRUE(postok)) n <- c(n, "Inadmissible solution (negative variance or non-positive-definite covariance matrix, a 'Heywood case'): do not interpret; the model is misspecified or the sample is too small.")
  if (length(w)) n <- c(n, paste("lavaan warnings:", paste(utils::head(unique(w), 3), collapse = " | ")))
  if (!is.null(fm)) {
    if (fm[["df"]] == 0) n <- c(n, "The model is saturated (df = 0): global fit cannot be tested.")
    else n <- c(n, sprintf("Global fit: chi-square = %.2f, df = %d, P = %.3f; CFI = %.3f; TLI = %.3f; RMSEA = %.3f (90%% CI %.3f-%.3f); SRMR = %.3f. %s",
                           fm[["chisq"]], as.integer(fm[["df"]]), fm[["pvalue"]], fm[["cfi"]], fm[["tli"]], fm[["rmsea"]],
                           fm[["rmsea.ci.lower"]], fm[["rmsea.ci.upper"]], fm[["srmr"]],
                           if (fm[["pvalue"]] < 0.05) "P < 0.05: the model-implied covariances differ from the observed ones (the model omits something)."
                           else "P >= 0.05: no evidence of misfit (this does not prove the model is right)."))
  }
  if (!is.null(mi) && nrow(mi)) {
    top <- utils::head(mi, 5)
    n <- c(n, paste0("Largest modification indices (possible missing links; add only with a biological reason): ",
                     paste(sprintf("%s %s %s (MI = %.1f)", top$lhs, top$op, top$rhs, top$mi), collapse = "; "), "."))
  }
  n
}
