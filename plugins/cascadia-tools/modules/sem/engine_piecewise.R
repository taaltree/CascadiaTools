# Piecewise SEM engine (Shipley 2000, 2009; Lefcheck 2016): each endogenous
# variable gets its own (G)LM(M); the global model is tested with Shipley's
# d-separation test (Fisher's C) via piecewiseSEM::dSep. Path coefficients,
# standardised coefficients (piecewiseSEM conventions, recomputed here so the
# Monte Carlo / bootstrap intervals use the same definition) and indirect /
# total effects (sums over directed paths of products of coefficients).

sem_fit_component <- function(c, df) {
  f <- c$formula
  environment(f) <- globalenv()
  re <- length(c$groups) > 0
  # Mixed models: fit with centred/scaled numeric predictors (lme4 convergence
  # checks are scale-sensitive) and back-transform coefficients afterwards.
  scaling <- NULL
  if (re) {
    num <- c$preds[vapply(c$preds, function(x) is.numeric(df[[x]]), logical(1))]
    if (length(num)) {
      scaling <- list(center = vapply(num, function(x) mean(df[[x]]), 1), scale = vapply(num, function(x) stats::sd(df[[x]]), 1))
      for (x in num) df[[x]] <- (df[[x]] - scaling$center[[x]]) / scaling$scale[[x]]
    }
  }
  fam <- switch(c$family, gaussian = stats::gaussian(), poisson = stats::poisson(),
                binomial = stats::binomial(link = c$link), nbinom = NULL)
  w <- character()
  m <- withCallingHandlers(
    if (c$family == "gaussian" && !re) do.call(stats::lm, list(formula = f, data = df))
    else if (c$family == "gaussian") do.call(lme4::lmer, list(formula = f, data = df, REML = TRUE))
    else if (c$family == "nbinom" && !re) do.call(MASS::glm.nb, list(formula = f, data = df))
    else if (c$family == "nbinom") do.call(lme4::glmer.nb, list(formula = f, data = df))
    else if (!re) do.call(stats::glm, list(formula = f, family = fam, data = df))
    else do.call(lme4::glmer, list(formula = f, family = fam, data = df)),
    warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") },
    message = function(e) invokeRestart("muffleMessage"))
  conv <- TRUE
  if (inherits(m, "glm")) conv <- isTRUE(m$converged)
  if (inherits(m, "merMod")) {
    msg <- m@optinfo$conv$lme4$messages
    if (length(msg) && any(!grepl("singular", msg))) conv <- FALSE
  }
  list(model = m, warnings = unique(w), converged = conv, scaling = scaling,
       singular = inherits(m, "merMod") && lme4::isSingular(m))
}

sem_coef_vcov <- function(fo) {
  m <- fo$model
  out <- if (inherits(m, "merMod")) list(b = lme4::fixef(m), V = as.matrix(stats::vcov(m)), df = Inf)
         else list(b = stats::coef(m), V = stats::vcov(m), df = if (inherits(m, "lm") && !inherits(m, "glm")) m$df.residual else Inf)
  sc <- fo$scaling
  if (!is.null(sc)) {
    nm <- names(out$b); A <- diag(length(nm)); dimnames(A) <- list(nm, nm)
    for (x in names(sc$center)) if (x %in% nm) {
      A[x, x] <- 1 / sc$scale[[x]]
      if ("(Intercept)" %in% nm) A["(Intercept)", x] <- -sc$center[[x]] / sc$scale[[x]]
    }
    out$b <- stats::setNames(as.numeric(A %*% out$b), nm)
    out$V <- A %*% out$V %*% t(A); dimnames(out$V) <- list(nm, nm)
  }
  out
}

sem_link_pred <- function(m) {
  if (inherits(m, "glm") || inherits(m, "merMod")) stats::predict(m, type = "link") else stats::fitted(m)
}

# Per-component results: coefficient rows, numeric-edge theta/V, std factors
sem_component_table <- function(c, fo, df) {
  cv <- sem_coef_vcov(fo)
  y <- c$response
  eta <- as.numeric(sem_link_pred(fo$model))
  yobs <- df[[y]]
  if (length(yobs) != length(eta)) yobs <- stats::model.frame(fo$model)[[y]]
  sdy <- sem_sd_y(yobs, eta, c$family, c$link)
  rows <- list(); theta <- c(); fac <- c(); keys <- c()
  qz <- if (is.finite(cv$df)) stats::qt(0.975, cv$df) else stats::qnorm(0.975)
  for (x in c$preds) {
    cn <- if (x %in% names(cv$b)) x else grep(paste0("^", x), names(cv$b), value = TRUE)
    cn <- setdiff(cn, "(Intercept)")
    numeric_x <- is.numeric(df[[x]]) && identical(cn, x)
    for (k in cn) {
      b <- cv$b[[k]]; se <- sqrt(cv$V[k, k])
      lev <- if (numeric_x) "" else sub(paste0("^", x), "", k)
      nm <- if (numeric_x) paste0(x, "_to_", y) else paste0(x, "_", make.names(lev), "_to_", y)
      st <- if (numeric_x) stats::sd(df[[x]]) / sdy else NA_real_
      p <- 2 * (if (is.finite(cv$df)) stats::pt(-abs(b / se), cv$df) else stats::pnorm(-abs(b / se)))
      rows[[length(rows) + 1]] <- data.frame(
        param = paste0("path_", nm), estimate = b, sd = se, lower = b - qz * se, upper = b + qz * se,
        kind = "path (unstandardised, link scale)", from = x, to = y, level = lev, p_value = p,
        std_factor = st, family = c$family, link = c$link, stringsAsFactors = FALSE)
      if (numeric_x) { theta[nm] <- b; fac[nm] <- st; keys[nm] <- k }
    }
  }
  tab <- if (length(rows)) do.call(rbind, rows) else NULL
  Vn <- if (length(keys)) cv$V[keys, keys, drop = FALSE] else matrix(0, 0, 0)
  if (length(keys)) dimnames(Vn) <- list(names(keys), names(keys))
  list(tab = tab, theta = theta, V = Vn, fac = fac, sd_y = sdy, df = cv$df)
}

sem_blockdiag <- function(mats) {
  k <- sum(vapply(mats, nrow, 1L)); out <- matrix(0, k, k); i <- 0
  nm <- unlist(lapply(mats, rownames))
  for (m in mats) if (nrow(m)) { j <- i + seq_len(nrow(m)); out[j, j] <- m; i <- i + nrow(m) }
  dimnames(out) <- list(nm, nm); out
}

sem_prepare_df <- function(df, comps) {
  for (v in names(df)) if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  used <- unique(unlist(lapply(comps, function(c) c(c$response, c$preds, c$groups))))
  keep <- stats::complete.cases(df[, intersect(used, names(df)), drop = FALSE])
  df[keep, , drop = FALSE]
}

sem_fit_all_components <- function(comps, df) {
  fos <- lapply(comps, sem_fit_component, df = df)
  names(fos) <- vapply(comps, `[[`, "", "response")
  fos
}

sem_collect <- function(comps, fos, df) {
  parts <- lapply(seq_along(comps), function(i) sem_component_table(comps[[i]], fos[[i]], df))
  list(tab = do.call(rbind, lapply(parts, `[[`, "tab")),
       theta = unlist(lapply(parts, function(p) p$theta)),
       fac = unlist(lapply(parts, function(p) p$fac)),
       V = sem_blockdiag(lapply(parts, `[[`, "V")),
       sd_y = vapply(parts, `[[`, 1, "sd_y"))
}

sem_dsep <- function(models, df, comps, st) {
  out <- list(table = NULL, C = NA_real_, df = NA_real_, p = NA_real_, k = 0L, note = NULL, saturated = FALSE)
  # d-sep refits use psem's data: give it a copy with Gaussian / predictor
  # columns standardised (t and z tests are invariant to linear rescaling;
  # count and binary responses keep their original values).
  keepraw <- unlist(lapply(comps, function(c) if (c$family != "gaussian") c$response))
  dsd <- df
  for (v in names(dsd)) if (is.numeric(dsd[[v]]) && !v %in% keepraw && isTRUE(stats::sd(dsd[[v]]) > 0)) dsd[[v]] <- as.numeric(scale(dsd[[v]]))
  ps <- tryCatch(suppressWarnings(suppressMessages(do.call(piecewiseSEM::psem, c(unname(models), list(data = dsd))))),
                 error = function(e) e)
  if (inherits(ps, "error")) { out$note <- paste("d-separation test could not be run:", conditionMessage(ps)); return(out) }
  bs <- tryCatch(suppressWarnings(piecewiseSEM::basisSet(ps)), error = function(e) list())
  if (!length(bs)) { out$saturated <- TRUE; out$note <- "The model is saturated (every pair of variables is linked), so it has no testable independence claims: d-separation cannot check it."; out$psem <- ps; return(out) }
  # Direction of each independence claim. piecewiseSEM duplicates claims that
  # involve non-Gaussian variables (testing y1 ~ y2 and y2 ~ y1). "causal"
  # (default) keeps only the claim whose response is the causally later
  # variable, as in Shipley's (2000, 2009) basis set; "conserve" keeps both
  # and uses the smaller P value (piecewiseSEM conserve = TRUE; anti-
  # conservative, see VALIDATION.md).
  dirn <- st$dsep_direction %||% "causal"
  if (dirn == "causal") {
    pos <- stats::setNames(seq_along(comps), vapply(comps, `[[`, "", "response"))
    tp <- function(v) if (v %in% names(pos)) pos[[v]] else 0
    key <- vapply(bs, function(b) paste(sort(b[1:2]), collapse = "|"), "")
    keep <- vapply(seq_along(bs), function(i) {
      dup <- sum(key == key[i]) > 1
      !dup || tp(bs[[i]][2]) > tp(bs[[i]][1])
    }, logical(1))
    bs2 <- bs[keep]; class(bs2) <- "basisSet"
    ds <- tryCatch(suppressWarnings(suppressMessages(piecewiseSEM::dSep(ps, basis.set = bs2, conserve = FALSE, .progressBar = FALSE))),
                   error = function(e) e)
  } else {
    ds <- tryCatch(suppressWarnings(suppressMessages(piecewiseSEM::dSep(ps, conserve = TRUE, .progressBar = FALSE))),
                   error = function(e) e)
  }
  if (inherits(ds, "error")) { out$note <- paste("d-separation test failed:", conditionMessage(ds)); return(out) }
  ds <- as.data.frame(ds)
  names(ds) <- gsub("\\.$", "", names(ds))
  pcol <- grep("^P\\.Value", names(ds), value = TRUE)[1]
  p <- as.numeric(ds[[pcol]])
  claim <- trimws(as.character(ds[[1]]))
  resp <- trimws(sub("~.*", "", claim)); pred <- trimws(sub("\\+.*", "", sub(".*~", "", claim)))
  tab <- data.frame(claim = paste0(pred, " -> ", resp, " | parents"), missing_from = pred, missing_to = resp,
                    test = ds$Test.Type, df = as.numeric(ds$DF), statistic = as.numeric(ds$Crit.Value),
                    p_value = p, stringsAsFactors = FALSE)
  pc <- pmax(p, 1e-300)
  out$table <- tab; out$k <- nrow(tab)
  out$C <- -2 * sum(log(pc)); out$df <- 2 * nrow(tab)
  out$p <- stats::pchisq(out$C, out$df, lower.tail = FALSE)
  out$psem <- ps; out$direction <- dirn
  out
}

sem_r2 <- function(ps) {
  r <- tryCatch(suppressWarnings(suppressMessages(piecewiseSEM::rsquared(ps))), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) return(NULL)
  rc <- intersect(c("R.squared", "Marginal"), names(r))[1]
  if (is.na(rc)) return(NULL)
  out <- stats::setNames(as.numeric(r[[rc]]), r$Response)
  if ("Conditional" %in% names(r)) attr(out, "conditional") <- stats::setNames(as.numeric(r$Conditional), r$Response)
  out
}

sem_fit_piecewise <- function(data, st) {
  comps <- sem_topo_order(sem_parse_models(data$spec$models))
  df <- sem_prepare_df(data$df, comps)
  fos <- sem_fit_all_components(comps, df)
  models <- lapply(fos, `[[`, "model")
  col <- sem_collect(comps, fos, df)
  tab <- col$tab
  # standardised path rows
  stdrows <- tab[!is.na(tab$std_factor), ]
  if (nrow(stdrows)) {
    s <- stdrows$std_factor
    stdrows <- transform(stdrows, param = sub("^path_", "std_", param), estimate = estimate * s, sd = sd * s,
                         lower = lower * s, upper = upper * s, kind = "path (standardised)")
  }
  # odds / rate ratios for non-identity links
  ratrows <- tab[tab$link %in% c("log", "logit"), ]
  if (nrow(ratrows))
    ratrows <- transform(ratrows, param = paste0(ifelse(link == "log", "rate_ratio_", "odds_ratio_"), sub("^path_", "", param)),
                         estimate = exp(estimate), sd = exp(estimate) * sd, lower = exp(lower), upper = exp(upper),
                         kind = ifelse(link == "log", "rate ratio per unit", "odds ratio per unit"))
  num <- tab[!is.na(tab$std_factor), ]
  edges <- data.frame(from = num$from, to = num$to, key = paste0(num$from, "_to_", num$to), stringsAsFactors = FALSE)
  pairs <- sem_pairs_from_spec(data$spec$effects, edges)
  eff <- NULL; effinfo <- NULL
  if (nrow(pairs)) {
    for (i in seq_len(nrow(pairs))) for (v in c(pairs$from[i], pairs$to[i]))
      if (!v %in% c(edges$from, edges$to)) stop(sprintf("effects: variable '%s' is not a numeric variable in the model", v))
    g <- sem_effect_fun(edges, pairs)
    gs <- sem_effect_fun(edges, pairs, "std_")
    fac <- col$fac[edges$key]
    gstd <- function(th) gs(th * fac)
    method <- st$ci_method %||% "montecarlo"
    draws <- NULL; draws_std <- NULL
    if (method == "montecarlo") {
      set.seed(st$seed %||% 1)
      draws <- sem_rmvn(st$n_mc %||% 4000, col$theta[edges$key], col$V[edges$key, edges$key, drop = FALSE])
      colnames(draws) <- edges$key
    } else if (method == "bootstrap") {
      bt <- sem_bootstrap(comps, df, edges, st)
      draws <- bt$theta; draws_std <- bt$theta_std
      if (nrow(draws) < 0.8 * (st$n_boot %||% 500)) warning("many bootstrap refits failed")
    }
    th <- col$theta[edges$key]
    e1 <- sem_ci_rows(g, th, method, col$V[edges$key, edges$key, drop = FALSE], draws)
    e2 <- if (method == "bootstrap") sem_ci_rows(gs, th * fac, "montecarlo", draws = draws_std) else sem_ci_rows(gstd, th, method, col$V[edges$key, edges$key, drop = FALSE], draws)
    e1$kind <- paste0(sub("_.*", "", e1$param), " effect (unstandardised, link scale of response)")
    e2$kind <- paste0(sub("^std_([a-z]+)_.*", "\\1", e2$param), " effect (standardised)")
    eff <- rbind(e1, e2)
    specs <- attr(g, "specs")
    nodir <- vapply(specs, function(s) !any(s$direct), logical(1))
    for (s in specs[nodir]) eff <- eff[!eff$param %in% paste0(c("direct_", "std_direct_"), s$from, "_to_", s$to), ]
    fam <- stats::setNames(vapply(comps, `[[`, "", "family"), vapply(comps, `[[`, "", "response"))
    lnk <- stats::setNames(vapply(comps, `[[`, "", "link"), vapply(comps, `[[`, "", "response"))
    effinfo <- do.call(rbind, lapply(specs, function(s) {
      med <- sem_mediators(s)
      data.frame(from = s$from, to = s$to, n_paths = length(s$paths), n_indirect = sum(!s$direct),
                 routes = paste(vapply(s$paths, paste, "", collapse = " -> "), collapse = "; "),
                 mediators = paste(med, collapse = ", "),
                 nongaussian_mediators = paste(med[fam[med] != "gaussian"], collapse = ", "),
                 sink_family = fam[[s$to]], sink_link = lnk[[s$to]], stringsAsFactors = FALSE)
    }))
    # ratio for total effects on log/logit sinks with gaussian mediators
    rr <- list()
    for (i in seq_len(nrow(effinfo))) {
      if (!effinfo$sink_link[i] %in% c("log", "logit") || nzchar(effinfo$nongaussian_mediators[i])) next
      b <- paste0(effinfo$from[i], "_to_", effinfo$to[i])
      r <- eff[eff$param == paste0("total_", b), ]
      pre <- if (effinfo$sink_link[i] == "log") "rate_ratio_total_" else "odds_ratio_total_"
      rr[[length(rr) + 1]] <- data.frame(param = paste0(pre, b), estimate = exp(r$estimate), sd = exp(r$estimate) * r$sd,
                                         lower = exp(r$lower), upper = exp(r$upper),
                                         kind = if (pre == "rate_ratio_total_") "rate ratio per unit (total effect)" else "odds ratio per unit (total effect)",
                                         stringsAsFactors = FALSE)
    }
    if (length(rr)) eff <- rbind(eff, do.call(rbind, rr))
  }
  ds <- if (isTRUE(st$dsep %||% TRUE)) sem_dsep(models, df, comps, st) else list(note = "d-separation test skipped (settings: dsep: false)")
  r2 <- if (!is.null(ds$psem) && isTRUE(st$r2 %||% TRUE)) sem_r2(ds$psem) else NULL
  keep <- c("param", "estimate", "sd", "lower", "upper", "kind")
  est <- rbind(tab[, c(keep, "p_value")], cbind(stdrows[, keep], p_value = stdrows$p_value),
               if (nrow(ratrows)) cbind(ratrows[, keep], p_value = ratrows$p_value),
               if (!is.null(eff)) cbind(eff[, keep], p_value = NA_real_))
  rownames(est) <- NULL
  conv <- all(vapply(fos, `[[`, TRUE, "converged"))
  notes <- sem_piecewise_notes(comps, fos, ds, r2, effinfo, df, data, st)
  list(estimates = est,
       diagnostics = list(converged = conv, notes = notes, engine = "piecewiseSEM",
                          fisher_C = ds$C, fisher_df = ds$df, fisher_p = ds$p, n_claims = ds$k,
                          model_checks_ok = if (is.na(ds$p %||% NA)) NA else ds$p >= 0.05),
       object = models, fos = fos, comps = comps, df = df, table = tab, dsep = ds, r2 = r2, effects_info = effinfo,
       sd_y = stats::setNames(col$sd_y, names(fos)), settings = st, engine = "piecewise")
}

sem_bootstrap <- function(comps, df, edges, st) {
  set.seed(st$seed %||% 1)
  nb <- st$n_boot %||% 500
  grp <- unique(unlist(lapply(comps, `[[`, "groups")))
  res <- lapply(seq_len(nb), function(b) {
    if (length(grp)) {
      lv <- unique(df[[grp[1]]]); pick <- sample(lv, length(lv), replace = TRUE)
      d <- do.call(rbind, lapply(seq_along(pick), function(j) { x <- df[df[[grp[1]]] == pick[j], , drop = FALSE]; x[[grp[1]]] <- paste0(pick[j], "_", j); x }))
      d[[grp[1]]] <- factor(d[[grp[1]]])
    } else d <- df[sample(nrow(df), replace = TRUE), , drop = FALSE]
    out <- tryCatch({
      fos <- sem_fit_all_components(comps, d)
      col <- sem_collect(comps, fos, d)
      list(th = col$theta[edges$key], ths = col$theta[edges$key] * col$fac[edges$key])
    }, error = function(e) NULL)
    out
  })
  res <- Filter(Negate(is.null), res)
  list(theta = do.call(rbind, lapply(res, `[[`, "th")), theta_std = do.call(rbind, lapply(res, `[[`, "ths")))
}

sem_piecewise_notes <- function(comps, fos, ds, r2, effinfo, df, data, st) {
  n <- character()
  n <- c(n, sprintf("Engine piecewiseSEM: %d component models fitted to %d rows%s.", length(comps), nrow(df),
                    if (nrow(df) < nrow(data$df)) sprintf(" (%d rows with missing values in model variables were dropped)", nrow(data$df) - nrow(df)) else ""))
  for (i in seq_along(comps)) {
    fo <- fos[[i]]; c <- comps[[i]]
    if (!fo$converged) n <- c(n, sprintf("NOT CONVERGED: component %s. Simplify it (fewer predictors, drop the random intercept).", c$text))
    if (isTRUE(fo$singular)) n <- c(n, sprintf("Component %s: random-intercept variance estimated at ~0 (singular fit); fixed effects are still usable.", c$text))
    if (length(fo$warnings)) n <- c(n, sprintf("Component %s warnings: %s", c$text, paste(utils::head(fo$warnings, 3), collapse = " | ")))
    if (c$family == "poisson") {
      m <- fo$model
      disp <- sum(stats::residuals(m, type = "pearson")^2) / max(1, stats::df.residual(m))
      if (disp > 1.5) n <- c(n, sprintf("Component %s: Pearson dispersion %.2f > 1.5, counts are overdispersed; use family: nbinom (intervals are too narrow otherwise).", c$text, disp))
    }
  }
  if (!is.null(ds$note)) n <- c(n, ds$note)
  if (!is.null(ds$table)) {
    n <- c(n, sprintf("d-separation test (Shipley 2009): Fisher's C = %.2f, df = %d, P = %.3f over %d independence claims. %s",
                      ds$C, as.integer(ds$df), ds$p, ds$k,
                      if (ds$p < 0.05) "P < 0.05: the data contradict the model; at least one missing link is needed (see the missing-path table)."
                      else "P >= 0.05: no evidence that the model omits a needed link (this does not prove the model is right)."))
    bad <- ds$table[ds$table$p_value < 0.05, ]
    if (nrow(bad)) n <- c(n, paste0("Missing links with P < 0.05 in the d-sep test: ",
                                    paste(sprintf("%s -> %s (P = %s)", bad$missing_from, bad$missing_to, fmt_sig(bad$p_value, 2)), collapse = "; "),
                                    ". (Each claim is a separate test, so with several claims one P < 0.05 can occur by chance; judge the model by Fisher's C first.) Consider adding them if they make biological sense, or a correlated error (%~~%) if both share an unmeasured cause."))
  }
  if (length(r2)) n <- c(n, paste0("R2 by response: ", paste(sprintf("%s = %.2f", names(r2), r2), collapse = "; "), "."))
  if (!is.null(effinfo)) for (i in seq_len(nrow(effinfo)))
    if (nzchar(effinfo$nongaussian_mediators[i]))
      n <- c(n, sprintf("Indirect effect %s -> %s passes through non-Gaussian mediator(s) %s: the unstandardised product mixes link and response scales and is not a well-defined effect; compare the standardised values only as relative importance (Grace et al. 2018).",
                        effinfo$from[i], effinfo$to[i], effinfo$nongaussian_mediators[i]))
  nr <- nrow(df); np <- max(vapply(comps, function(c) length(c$preds) + 1, 1))
  if (nr < 5 * np * 2) n <- c(n, sprintf("Only %d rows for components with up to %d coefficients: path estimates are imprecise and the d-sep test has little power to detect missing links.", nr, np))
  n
}
