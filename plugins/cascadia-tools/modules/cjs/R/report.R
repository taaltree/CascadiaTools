# cjs: predictions, interpretation, figures and methods text.

cjs_fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
cjs_ci <- function(r, bayes, d = 2) sprintf("%s (95%% %s %s-%s)", cjs_fmt(r$estimate, d), if (bayes) "CrI" else "CI",
                                           cjs_fmt(r$lower, d), cjs_fmt(r$upper, d))
cjs_unit <- function(data, v) { u <- data$units[[v]]; if (is.null(u) || !nzchar(u)) "" else paste0(" ", u) }
cjs_nice <- function(v) { s <- gsub("_", " ", v); paste0(toupper(substr(s, 1, 1)), substring(s, 2)) }
# Lower-case the first letter for use mid-sentence, unless it starts an acronym ("PIT tag").
cjs_lc1 <- function(s) ifelse(grepl("^[A-Z][a-z]", s), paste0(tolower(substr(s, 1, 1)), substring(s, 2)), s)
# User labels (config figures: labels:, carried in data$labels and opts) for variables and levels.
cjs_L <- function(key, default, data = NULL, opts = NULL) {
  if (!is.null(opts$labels[[key]])) return(opts$labels[[key]])
  if (!is.null(data$labels[[key]])) return(data$labels[[key]])
  default
}
cjs_vlab <- function(data, v, opts = NULL) cjs_L(v, cjs_nice(v), data, opts)
# Map each "_"-separated token of an estimand suffix ("F", "male_adult") through the labels.
cjs_glab <- function(x, data, opts = NULL) vapply(x, function(z) {
  if (!nzchar(z)) return("all")
  if (!is.null(cjs_L(z, NULL, data, opts))) return(cjs_L(z, z, data, opts))
  paste(vapply(strsplit(z, "_")[[1]], function(t) cjs_L(t, gsub("_", " ", t), data, opts), ""), collapse = " ")
}, "", USE.NAMES = FALSE)

# Readable name of a design-matrix column (coefficient), e.g. "sex_M" -> "Sex: Male vs Female".
cjs_term_label <- function(term, des, data, process = "phi", opts = NULL) {
  if (term == "tsm_first") return("First interval after marking")
  if (term == "td") return("Caught at previous occasion")
  if (term == "age_class_juvenile") return("Juvenile (first interval)")
  if (term == "Time") return("Linear trend (per occasion)")
  if (grepl("^tsm_first_", term))
    return(paste0("First interval x ", cjs_lc1(cjs_term_label(sub("^tsm_first_", "", term), des, data, process, opts))))
  if (grepl("^time_", term)) {
    lab <- sub("^time_", "", term)
    return(if (process == "phi") paste("Interval starting", lab) else paste("Occasion", lab))
  }
  mmo <- if (process == "phi") des$mphi else des$mp
  for (v in setdiff(names(mmo$xlev), c("time", "tsm", "age_class"))) {
    lv <- sub("^_", "", mmo$xlev[[v]])
    for (l in lv[-1]) if (cjs_clean(paste0(v, "_", l)) == term)
      return(sprintf("%s: %s vs %s", cjs_vlab(data, v, opts), cjs_L(l, l, data, opts), cjs_L(lv[1], lv[1], data, opts)))
  }
  if (term %in% names(des$scaling))
    return(if (isTRUE(des$scaling[[term]]$binary)) cjs_vlab(data, term, opts) else paste0(cjs_vlab(data, term, opts), ", per SD"))
  term
}

# Response-scale predictions for design rows X of one process. For Nimble,
# `eps` (draws x rows) adds random time effects.
cjs_pred_X <- function(obj, process, X, ci = 0.95, eps = NULL) {
  nm <- names(obj$theta)
  idx <- setdiff(grep(paste0("^", process, ":"), nm), which(nm == "phi:sd_time"))
  a <- (1 - ci) / 2
  if (obj$engine == "marked") {
    V <- obj$vcov[idx, idx, drop = FALSE]
    lp <- drop(X %*% obj$theta[idx]); se <- sqrt(pmax(rowSums((X %*% V) * X), 0))
    z <- stats::qnorm(1 - a); z5 <- stats::qnorm(0.75)
    return(data.frame(estimate = stats::plogis(lp), lower = stats::plogis(lp - z * se), upper = stats::plogis(lp + z * se),
                      lower50 = stats::plogis(lp - z5 * se), upper50 = stats::plogis(lp + z5 * se)))
  }
  D <- obj$draws
  sel <- if (nrow(D) > 2000) round(seq(1, nrow(D), length.out = 2000)) else seq_len(nrow(D))
  eta <- X %*% t(D[sel, idx, drop = FALSE])
  if (!is.null(eps)) eta <- eta + t(eps[sel, , drop = FALSE])
  val <- stats::plogis(eta)
  data.frame(estimate = rowMeans(val), lower = apply(val, 1, stats::quantile, a), upper = apply(val, 1, stats::quantile, 1 - a),
             lower50 = apply(val, 1, stats::quantile, 0.25), upper50 = apply(val, 1, stats::quantile, 0.75))
}

cjs_proc_df <- function(des, process) if (process == "phi") des$dd$phi else des$dd$p[des$dd$p$.t > 1, , drop = FALSE]
cjs_proc_mm <- function(des, process) if (process == "phi") des$mphi else des$mp

# Conditional effect of a numeric covariate (original units), one curve per group.
cjs_cond_effect <- function(obj, process, var, ci = 0.95, n = 100) {
  des <- obj$des; df <- cjs_proc_df(des, process); mmo <- cjs_proc_mm(des, process)
  tmpl <- cjs_template(df, mmo)
  grp <- des$data$group
  glev <- if (!is.null(grp) && grp %in% mmo$vars) levels(df[[grp]]) else NA
  sc <- des$scaling[[var]]
  raw <- if (sc$type == "individual") des$data$ind[[var]] else des$data$occ[[var]][if (process == "phi") seq_len(des$T - 1) else 2:des$T]
  xs <- seq(min(raw, na.rm = TRUE), max(raw, na.rm = TRUE), length.out = n)
  out <- lapply(glev, function(g) {
    rows <- tmpl[rep(1, n), , drop = FALSE]
    rows[[var]] <- (xs - sc$center) / sc$scale
    if (!is.na(g)) rows[[grp]] <- factor(g, levels = levels(df[[grp]]))
    pr <- cjs_pred_X(obj, process, cjs_newX(mmo, rows), ci)
    pr$x <- xs
    if (!is.na(g)) pr$group <- cjs_lab(g)
    pr
  })
  list(pred = do.call(rbind, out), raw = raw[is.finite(raw)])
}

# Survival (or detection) per interval/occasion, one row per group.
cjs_by_time <- function(obj, process, ci = 0.95) {
  des <- obj$des; df <- cjs_proc_df(des, process); mmo <- cjs_proc_mm(des, process)
  T <- des$T; tt <- if (process == "phi") seq_len(T - 1) else 2:T
  tmpl <- cjs_template(df, mmo)
  grp <- des$data$group
  glev <- if (!is.null(grp) && grp %in% mmo$vars) levels(df[[grp]]) else NA
  occs <- names(des$data$occ)
  out <- lapply(glev, function(g) {
    rows <- tmpl[rep(1, length(tt)), , drop = FALSE]
    if ("time" %in% mmo$vars) rows$time <- factor(levels(df$time)[seq_along(tt)], levels = levels(df$time))
    if ("Time" %in% mmo$vars) rows$Time <- seq_along(tt) - 1
    for (v in intersect(occs, mmo$vars)) {
      sc <- des$scaling[[v]]
      rows[[v]] <- if (is.null(sc)) factor(df[[v]][match(tt, df$.t)], levels = levels(df[[v]])) else (des$data$occ[[v]][tt] - sc$center) / sc$scale
    }
    if (!is.na(g)) rows[[grp]] <- factor(g, levels = levels(df[[grp]]))
    eps <- if (process == "phi" && des$random_time && !is.null(obj$eps_chains)) do.call(rbind, obj$eps_chains) else NULL
    pr <- cjs_pred_X(obj, process, cjs_newX(mmo, rows), ci, eps = eps)
    pr$time <- if (process == "phi") paste0(des$data$occasions[tt], "-", des$data$occasions[tt + 1]) else des$data$occasions[tt]
    pr$t <- tt
    if (!is.na(g)) pr$group <- cjs_lab(g)
    pr
  })
  d <- do.call(rbind, out)
  d$time <- factor(d$time, levels = unique(d$time))
  if (des$confounded_last) d <- d[d$t != (if (process == "phi") T - 1 else T), , drop = FALSE]
  d
}

cjs_time_varies <- function(des, process) {
  mmo <- cjs_proc_mm(des, process)
  any(c("time", "Time", names(des$data$occ)) %in% mmo$vars) || (process == "phi" && des$random_time)
}

# When variation over time comes ONLY from occasion covariates (no time / Time / (1 | time)),
# the per-interval values are model predictions from those covariates, not year-specific
# estimates. Returns the covariate names in that case, otherwise NULL.
cjs_time_from_covs <- function(des, process) {
  mmo <- cjs_proc_mm(des, process)
  if (any(c("time", "Time") %in% mmo$vars) || (process == "phi" && des$random_time)) return(NULL)
  v <- intersect(names(des$data$occ), mmo$vars)
  if (length(v)) v else NULL
}

cjs_num_covs <- function(des, process) {
  mmo <- cjs_proc_mm(des, process)
  v <- intersect(mmo$vars, names(des$scaling))
  v <- v[!vapply(v, function(k) isTRUE(des$scaling[[k]]$binary), TRUE)]
  if (process == "phi") v else v[vapply(v, function(k) des$scaling[[k]]$type == "individual", TRUE)]
}

cjs_interpret <- function(fit, data) {
  obj <- fit$object; des <- obj$des; est <- fit$estimates
  bayes <- obj$engine != "marked"
  row <- function(n) { r <- est[est$param == n, , drop = FALSE]; if (nrow(r)) r else NULL }
  unit <- data$occasion_unit %||% "occasion"
  per <- if (unit == "year") "annual" else paste0("per-", unit)
  out <- character()
  ch <- data$ch
  out <- c(out, sprintf("Data: %d marked %s over %d %s occasions (%s to %s); %d captures in total, %d individuals (%.0f%%) were recaptured at least once after marking.",
                        nrow(ch), data$species %||% "individuals", ncol(ch), unit, data$occasions[1], data$occasions[ncol(ch)], sum(ch),
                        sum(rowSums(ch) > 1), 100 * mean(rowSums(ch) > 1)))
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge (see Diagnostics). The numbers below should not be interpreted.")
  # survival rows (not time-specific)
  sv <- est[grepl("^survival", est$param) & !grepl("temporal_sd|_x_detection", est$param), , drop = FALSE]
  covtxt <- c(cjs_num_covs(des, "phi"))
  bin <- intersect(des$mphi$vars, names(des$scaling)[vapply(des$scaling, function(z) isTRUE(z$binary), TRUE)])
  avg <- if (length(covtxt)) paste0(" at average ", paste(cjs_lc1(vapply(covtxt, function(v) cjs_vlab(data, v), "")), collapse = " and ")) else ""
  if (length(bin)) avg <- paste0(avg, " with ", paste(bin, collapse = ", "), " = 0")
  if (nrow(sv) && nrow(sv) <= 8) {
    out <- c(out, sprintf("Apparent %s survival%s%s: %s. Apparent survival is the probability of surviving AND staying in the study area; animals that emigrate permanently are indistinguishable from deaths, so true survival is at least this high (apparent survival = true survival x site fidelity).",
                          per, avg, if ("tsm" %in% des$mphi$vars) " of residents (intervals after the first one following marking)" else "",
                          paste(sprintf("%s %s", cjs_glab(sub("^survival_?", "", sv$param), data), cjs_ci(sv, bayes)), collapse = "; ")))
  } else if (nrow(sv)) {
    out <- c(out, sprintf("Apparent %s survival ranged from %.2f to %.2f across the %d groups/intervals estimated (table). Apparent survival = true survival x site fidelity.",
                          per, min(sv$estimate), max(sv$estimate), nrow(sv)))
  }
  if (!identical(des$formulas$phi, des$formulas$phi_user))
    out <- c(out, "Because transients are modelled (tsm), the survival effects below refer to residents (intervals after the first one following marking); the first interval after marking has its own effects, since transients dilute it.")
  # group / factor contrasts (survival and recapture)
  for (proc in c("phi", "p")) {
    pre <- paste0("^", proc, "_beta_")
    co <- est[grepl(pre, est$param) & !grepl(paste0(pre, "(time_|tsm|td$)"), est$param), , drop = FALSE]
    co <- co[!sub(pre, "", co$param) %in% c(cjs_num_covs(des, proc), "Time"), , drop = FALSE]
    what <- if (proc == "phi") "apparent survival" else "recapture probability"
    for (k in seq_len(nrow(co))) {
      r <- co[k, ]
      term <- sub(pre, "", r$param)
      # 0/1 covariates (e.g. a flood-year indicator) are not groups
      where <- if (isTRUE(des$scaling[[term]]$binary)) sprintf("when %s = 1 than when it is 0", cjs_lc1(cjs_vlab(data, term))) else "in this group"
      out <- c(out, sprintf("%s (%s): odds ratio %s (logit coefficient %.2f). %s",
                            cjs_term_label(sub(pre, "", r$param), des, data, proc), what,
                            cjs_ci(data.frame(estimate = exp(r$estimate), lower = exp(r$lower), upper = exp(r$upper)), bayes),
                            r$estimate, if (r$lower > 0) sprintf("The interval is above 1: %s is higher %s.", what, where)
                                        else if (r$upper < 0) sprintf("The interval is below 1: %s is lower %s.", what, where)
                                        else "The interval includes 1: the data do not show a clear difference."))
    }
  }
  # numeric covariates on survival / detection
  for (proc in c("phi", "p")) for (v in cjs_num_covs(des, proc)) {
    r <- row(paste0(proc, "_beta_", v)); if (is.null(r)) next
    sc <- des$scaling[[v]]
    ce <- cjs_cond_effect(obj, proc, v)
    raw <- ce$raw; q <- stats::quantile(raw, c(0.1, 0.9), names = FALSE)
    g1 <- ce$pred; if ("group" %in% names(g1)) g1 <- g1[g1$group == g1$group[1], ]
    lo <- g1[which.min(abs(g1$x - q[1])), ]; hi <- g1[which.min(abs(g1$x - q[2])), ]
    what <- if (proc == "phi") "apparent survival" else "detection probability"
    out <- c(out, sprintf("%s: each 1 SD (%s%s) increase multiplies the odds of %s by %s. Holding other covariates at their means%s, %s is %.2f at the 10th percentile (%s%s) and %.2f at the 90th percentile (%s%s). %s This is an association, not proof of a causal effect.",
                          cjs_vlab(data, v), cjs_fmt(sc$scale, 2), cjs_unit(data, v), what,
                          cjs_ci(data.frame(estimate = exp(r$estimate), lower = exp(r$lower), upper = exp(r$upper)), bayes),
                          if ("group" %in% names(ce$pred)) paste0(" (", cjs_glab(ce$pred$group[1], data), ")") else "",
                          what, lo$estimate, cjs_fmt(q[1], 1), cjs_unit(data, v), hi$estimate, cjs_fmt(q[2], 1), cjs_unit(data, v),
                          if (r$lower > 0) "The interval excludes no effect (positive association)." else if (r$upper < 0) "The interval excludes no effect (negative association)." else "The interval includes no effect."))
  }
  if (!is.null(row("phi_beta_Time"))) {
    r <- row("phi_beta_Time")
    out <- c(out, sprintf("Linear trend in survival: logit change per %s %s (odds ratio per %s %s).", unit, cjs_ci(r, bayes), unit,
                          cjs_ci(data.frame(estimate = exp(r$estimate), lower = exp(r$lower), upper = exp(r$upper)), bayes)))
  }
  if (cjs_time_varies(des, "phi")) {
    bt <- cjs_by_time(obj, "phi")
    tc <- cjs_time_from_covs(des, "phi")
    out <- c(out, sprintf("%s from %.2f (%s) to %.2f (%s)%s%s.",
                          if (is.null(tc)) "Survival varied between intervals" else "Predicted survival per interval ranged",
                          min(bt$estimate), bt$time[which.min(bt$estimate)],
                          max(bt$estimate), bt$time[which.max(bt$estimate)], if ("group" %in% names(bt)) paste0(" (across ", paste(cjs_glab(unique(bt$group), data), collapse = " and "), ")") else "",
                          if (is.null(tc)) "" else sprintf("; these values follow only from the interval's %s (the model has no separate year effect, so other year-to-year variation is not captured)",
                                                           paste(cjs_lc1(vapply(tc, function(v) cjs_vlab(data, v), "")), collapse = " and "))))
  }
  if (!is.null(row("survival_temporal_sd"))) {
    r <- row("survival_temporal_sd")
    out <- c(out, sprintf("Temporal variability: the SD of survival among intervals on the logit scale is %s. With few intervals this variance is imprecise and partly prior-driven (see Known limits).", cjs_ci(r, bayes)))
  }
  tr <- est[grepl("^proportion_transient", est$param), , drop = FALSE]
  # 1 - phi_first/phi_later is unbounded when estimated; a proportion is shown within [0, 1]
  for (k in c("estimate", "lower", "upper")) tr[[k]] <- pmin(pmax(tr[[k]], 0), 1)
  if (nrow(tr)) out <- c(out, sprintf("Transients: an estimated %s of newly marked animals were transients (passing through and never returning), computed as 1 - survival(first interval)/survival(later intervals) (Pradel et al. 1997; values truncated to 0-1 for display).",
                                      paste(sprintf("%s%s", cjs_ci(tr, bayes), ifelse(grepl("_", sub("^proportion_transient", "", tr$param)), paste0(" (", cjs_glab(sub("^proportion_transient_", "", tr$param), data), ")"), "")), collapse = "; ")))
  dv <- est[grepl("^detection", est$param), , drop = FALSE]
  if (nrow(dv) && nrow(dv) <= 6) out <- c(out, sprintf("Detection (recapture) probability per occasion: %s.",
                                                       paste(sprintf("%s %s", ifelse(dv$param == "detection", "all", cjs_glab(sub("^detection_", "", dv$param), data)), cjs_ci(dv, bayes)), collapse = "; ")))
  else if (nrow(dv)) out <- c(out, sprintf("Detection probability ranged from %.2f to %.2f across occasions/groups.", min(dv$estimate), max(dv$estimate)))
  if (!is.null(row("p_beta_td"))) {
    r <- row("p_beta_td")
    out <- c(out, sprintf("Trap dependence: animals caught at the previous occasion had odds of recapture multiplied by %s (%s).",
                          cjs_ci(data.frame(estimate = exp(r$estimate), lower = exp(r$lower), upper = exp(r$upper)), bayes),
                          if (r$lower > 0) "trap-happy" else if (r$upper < 0) "trap-shy" else "no clear trap response"))
  }
  if (des$confounded_last) out <- c(out, "In this fully time-dependent model the last survival and last detection probabilities are confounded; only their product is reported.")
  g <- obj$gof
  if (!is.null(g)) {
    gm <- if (is.null(g$model)) cjs_gof_model(g, des) else g$model
    tb <- g$table
    comp <- function(lbl) tb[tb$test == lbl, ]
    txt <- sprintf("Goodness of fit (R2ucare tests of the time-dependent CJS model): overall chi-square %.1f on %d df, p = %.3f, c-hat = %.2f.",
                   tb$chi2[nrow(tb)], as.integer(tb$df[nrow(tb)]), g$p, g$chat)
    if ("3.SR (transience)" %in% gm$flagged) {
      r <- comp("3.SR (transience)")
      txt <- c(txt, sprintf("Test 3.SR detected transients (p = %.3f): %s", r$p_value,
                            if ("tsm" %in% des$mphi$vars) "the fitted model accounts for them (tsm)." else "add 'tsm' to the survival formula; otherwise survival is biased low."))
    }
    if ("2.CT (trap dependence)" %in% gm$flagged) {
      r <- comp("2.CT (trap dependence)")
      txt <- c(txt, sprintf("Test 2.CT detected trap dependence (p = %.3f): %s", r$p_value,
                            if ("td" %in% des$mp$vars) "the fitted model accounts for it (td)." else "add 'td' to the detection formula."))
    }
    if (length(gm$excluded))
      txt <- c(txt, sprintf("Excluding the component(s) the model accounts for, c-hat = %.2f (p = %.3f).", gm$chat, gm$p))
    unaddressed <- setdiff(gm$flagged, gm$excluded)
    txt <- c(txt, if (is.finite(gm$chat) && gm$chat <= 1.5 && isTRUE(gm$p >= 0.05 || is.na(gm$p)) && !length(unaddressed))
                    "No evidence of remaining lack of fit."
                  else if (is.finite(gm$chat) && gm$chat > 3) "c-hat > 3: a structural problem (unmodelled transients, trap dependence or groups); revise the model before relying on these estimates."
                  else "Some lack of fit remains (see the component tests in Diagnostics and the gof figure); consider settings chat_adjust: true, which widens intervals by sqrt(c-hat).")
    out <- c(out, paste(txt, collapse = " "))
  }
  if (nrow(ch) < 150 || ncol(ch) < 5)
    out <- c(out, sprintf("Caution: this study (%d individuals, %d occasions) is at or below the smallest design validated for this module (see VALIDATION.md, harder design): intervals may be too narrow or estimates biased; run `recover` at your own design before relying on them.", nrow(ch), ncol(ch)))
  out
}

cjs_methods_text <- function(fit, data, settings) {
  obj <- fit$object; des <- obj$des; s <- obj$settings
  eng <- switch(obj$engine,
    marked = "by maximum likelihood with the R package marked (Laake et al. 2013); 95% confidence intervals for probabilities are Wald intervals on the logit scale",
    nimble = sprintf("in a Bayesian framework with NIMBLE (de Valpine et al. 2017) using the marginal CJS likelihood (dCJS distribution, nimbleEcology; Goldstein et al. 2024; Turek et al. 2016), %d chains of %d iterations after %d burn-in (automated-factor slice samplers); priors: Logistic(0,1) on logit-scale intercepts (uniform on the probability scale), Normal(0, SD %s) on other logit-scale coefficients%s",
                     s$nchains, s$niter - s$nburnin, s$nburnin, s$prior_sd, if (des$random_time) ", half-Normal(0, 1) on the SD of random time effects" else ""),
    nimble_latent = sprintf("in a Bayesian framework with NIMBLE using the state-space (latent alive state) formulation of Gimenez et al. (2007) and Kery & Schaub (2012), %d chains of %d iterations after %d burn-in; priors as in the marginal version", s$nchains, s$niter - s$nburnin, s$nburnin))
  covs <- intersect(c(des$mphi$vars, des$mp$vars), names(des$scaling))
  paste0(sprintf("Apparent survival (phi) and recapture probability (p) were estimated with Cormack-Jolly-Seber models (Lebreton et al. 1992) from the capture histories of %d individually marked %s over %d occasions. ",
                 nrow(data$ch), data$species %||% "individuals", ncol(data$ch)),
         sprintf("Survival was modelled as logit(phi) %s%s and detection as logit(p) %s. ", des$formulas$phi,
                 if (des$random_time) " plus a normally distributed random effect of time interval" else "", des$formulas$p),
         if (length(covs)) sprintf("Continuous covariates (%s) were standardised to mean 0 and SD 1 before fitting. ", paste(covs, collapse = ", ")) else "",
         if ("tsm" %in% des$mphi$vars) "Transients were accounted for by letting survival in the first interval after marking differ from later intervals (Pradel et al. 1997). " else "",
         if ("td" %in% des$mp$vars) "Trap dependence was modelled by letting recapture probability depend on capture at the previous occasion (Pradel 1993). " else "",
         sprintf("Models were fitted %s. ", eng),
         if (!is.null(obj$gof)) "Goodness of fit of the fully time-dependent CJS model was assessed with the U-CARE component tests (3.SR transience, 3.SM, 2.CT trap dependence, 2.CL; Pradel et al. 2005) using R2ucare (Gimenez et al. 2018), with c-hat = overall chi-square / df. " else "",
         if (isTRUE(s$chat_adjust)) "Standard errors were inflated by sqrt(c-hat) when c-hat > 1. " else "",
         "Survival is apparent survival: the product of true survival and the probability of remaining in the study area.")
}

# ---------------------------------------------------------------- figures
cjs_figures <- function(fit, data, opts = fig_opts()) {
  obj <- fit$object; des <- obj$des; est <- fit$estimates
  ci <- opts$ci %||% 0.95
  bayes <- obj$engine != "marked"
  itype <- if (bayes) "credible" else "confidence"
  pct <- round(100 * ci)
  unit <- data$occasion_unit %||% "occasion"
  ylab_phi <- fig_label("survival", if (unit == "year") "Apparent annual survival" else sprintf("Apparent survival\nper %s", unit), opts)
  ylab_p <- fig_label("detection", "Recapture probability", opts)
  figs <- list()
  covlab <- function(v) fig_label(v, paste0(cjs_nice(v), if (nzchar(cjs_unit(data, v))) paste0(" (", trimws(cjs_unit(data, v)), ")") else ""), opts)
  glab <- function(g) vapply(g, function(z) fig_label(z, z, opts), "")
  grp_cap <- function(pred) if ("group" %in% names(pred)) sprintf(", one curve per %s", fig_label(des$data$group, des$data$group, opts)) else ""
  effect_figs <- list()
  for (proc in c("phi", "p")) for (v in cjs_num_covs(des, proc)) {
    ce <- cjs_cond_effect(obj, proc, v, ci)
    pred <- ce$pred
    if ("group" %in% names(pred)) pred$group <- glab(pred$group)
    p <- fig_cond_effect(pred, covlab(v), if (proc == "phi") ylab_phi else ylab_p,
                         obs = if (isTRUE(opts$show_data)) data.frame(x = ce$raw) else NULL, ylim01 = TRUE)
    nm <- paste0(if (proc == "phi") "survival_" else "detection_", v)
    figs[[nm]] <- with_caption(p, sprintf("%s versus %s (line: %s; dark band: 50%%, light band: %d%% %s interval)%s; other covariates at their means. Rug: observed %s values%s.",
                                          if (proc == "phi") "Apparent survival" else "Recapture probability", cjs_lc1(covlab(v)),
                                          if (bayes) "posterior mean" else "maximum-likelihood estimate", pct, itype, grp_cap(pred),
                                          cjs_lc1(covlab(v)), if (des$scaling[[v]]$type == "occasion") " (one per interval)" else " (one per individual)"))
    if (proc == "phi") effect_figs[[nm]] <- p
  }
  if (length(effect_figs) > 1) {
    # the core styles only top-level ggplots: style each panel before composing
    sub <- lapply(names(effect_figs), function(n) apply_fig_opts(effect_figs[[n]], n, opts))
    pn <- compose_figures(sub, ncol = min(2, length(sub)))
    attr(pn, "fig_size") <- list(width = "double", aspect = if (length(sub) > 2) 0.75 else 0.45)
    figs$survival_effects_panel <- with_caption(pn,
      sprintf("Conditional effects of all continuous covariates on apparent survival (line: %s; bands: 50%% and %d%% %s intervals); other covariates at their means.",
              if (bayes) "posterior mean" else "estimate", pct, itype))
  }
  for (proc in c("phi", "p")) if (cjs_time_varies(des, proc)) {
    bt <- cjs_by_time(obj, proc, ci)
    if (proc == "phi") {
      # "2016-2017" -> "2016-17" for year labels (keeps the axis legible)
      yy <- des$data$occasions
      if (all(grepl("^[0-9]{4}$", yy)))
        levels(bt$time) <- sub("^([0-9]{4})-[0-9]{2}([0-9]{2})$", "\\1-\\2", levels(bt$time))
    }
    has_g <- "group" %in% names(bt)
    if (has_g) bt$group <- glab(bt$group)
    pd <- ggplot2::position_dodge(width = if (has_g) 0.45 else 0)
    m <- if (has_g) ggplot2::aes(time, estimate, colour = group, shape = group) else ggplot2::aes(time, estimate)
    p <- ggplot2::ggplot(bt, m) +
      do.call(ggplot2::geom_linerange, c(list(ggplot2::aes(ymin = lower, ymax = upper), linewidth = 0.45, position = pd), eco_fixed(has_g, colour = eco_pal$ink))) +
      do.call(ggplot2::geom_linerange, c(list(ggplot2::aes(ymin = lower50, ymax = upper50), linewidth = 1.2, position = pd), eco_fixed(has_g, colour = eco_pal$ink))) +
      do.call(ggplot2::geom_point, c(list(size = 1.8, position = pd), eco_fixed(has_g, colour = eco_pal$ink))) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = fig_label(if (proc == "phi") "interval" else "occasion", if (proc == "phi") "Interval" else "Occasion", opts),
                    y = if (proc == "phi") ylab_phi else ylab_p, colour = NULL, shape = NULL) +
      # a guide (not theme()) so the core theme applied at export does not undo the rotation
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = if (length(levels(bt$time)) > 4) 45 else 0))
    nm <- if (proc == "phi") "survival_time" else "detection_time"
    figs[[nm]] <- with_caption(p, sprintf("%s per %s (points: %s; thick bar: 50%%, thin bar: %d%% %s interval)%s%s%s.",
                                          if (proc == "phi") "Apparent survival" else "Recapture probability",
                                          if (proc == "phi") "interval between occasions" else "occasion",
                                          if (bayes) "posterior means" else "estimates", pct, itype,
                                          if (has_g) sprintf("; colours and symbols = %s", fig_label(des$data$group, des$data$group, opts)) else "",
                                          if (des$confounded_last) "; the last interval/occasion is omitted because it is confounded with the other process" else "",
                                          if (!is.null(tcv <- cjs_time_from_covs(des, proc))) sprintf("; values are predicted from %s only (no separate year effect in the model)",
                                                                                                     paste(cjs_lc1(vapply(tcv, function(v) fig_label(v, cjs_vlab(des$data, v), opts), "")), collapse = " and ")) else ""))
  }
  # survival / detection by group (categorical)
  for (proc in c("phi", "p")) {
    stem <- if (proc == "phi") "^survival" else "^detection"
    r <- est[grepl(stem, est$param) & !grepl("temporal_sd|_x_detection", est$param), , drop = FALSE]
    if (cjs_time_varies(des, proc) || nrow(r) < 2) next
    r$level <- glab(gsub("_", " ", sub(paste0(stem, "_?"), "", r$param)))
    p <- fig_cat_effect(r, if (!is.null(data$group) && is.null(data$age_col)) cjs_vlab(data, data$group, opts) else fig_label("group", "Group", opts), if (proc == "phi") ylab_phi else ylab_p) +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    figs[[if (proc == "phi") "survival_groups" else "detection_groups"]] <- with_caption(p,
      sprintf("%s by group (point: %s; bar: %d%% %s interval)%s.", if (proc == "phi") "Apparent survival" else "Recapture probability",
              if (bayes) "posterior mean" else "estimate", pct, itype,
              if (length(cjs_num_covs(des, proc))) "; continuous covariates at their means" else ""))
  }
  if (!cjs_time_varies(des, "p") && !"detection_groups" %in% names(figs)) {
    r <- est[grepl("^detection", est$param), , drop = FALSE]
    if (nrow(r) == 1) {
      r$level <- fig_label("all", "All occasions", opts)
      figs$detection <- with_caption(fig_cat_effect(r, NULL, ylab_p) + ggplot2::scale_y_continuous(limits = c(0, 1)),
                                     sprintf("Recapture probability per occasion, constant over occasions (point: estimate; bar: %d%% %s interval).", pct, itype))
    }
  }
  # coefficients (odds ratios)
  co <- est[grepl("^(phi|p)_beta_", est$param) & !grepl("^phi_beta_tsm_first_", est$param), , drop = FALSE]
  if (nrow(co)) {
    proc <- ifelse(startsWith(co$param, "phi_"), "phi", "p")
    co$process <- factor(ifelse(proc == "phi", fig_label("Survival", "Survival", opts), fig_label("Detection", "Recapture", opts)),
                         levels = c(fig_label("Survival", "Survival", opts), fig_label("Detection", "Recapture", opts)))
    co$param <- mapply(function(term, pr) fig_label(paste0(pr, "_beta_", term), cjs_term_label(term, des, data, pr, opts), opts),
                       sub("^(phi|p)_beta_", "", co$param), proc, USE.NAMES = FALSE)
    # one forest plot (no legend/facets: they squeeze the panel at single-column width);
    # recapture rows are prefixed, survival rows first
    co$param <- ifelse(proc == "p", paste0(fig_label("Detection", "Recapture", opts), ": ", co$param), co$param)
    co <- co[order(proc != "phi", -co$estimate), ]
    # interval for the plot at opts$ci (table rows are 95%): logit-scale normal approximation from sd
    if (abs(ci - 0.95) > 1e-9 && all(is.finite(co$sd))) { z <- stats::qnorm(1 - (1 - ci) / 2); co$lower <- co$estimate - z * co$sd; co$upper <- co$estimate + z * co$sd }
    wrap <- function(x) vapply(x, function(z) paste(strwrap(z, 30), collapse = "\n"), "")
    co$param <- factor(co$param, levels = rev(unique(co$param)))
    pc <- fig_coef(co[, c("param", "estimate", "lower", "upper")], null_value = 1, transform = exp, log_x = TRUE,
                   xlab = sprintf("Odds ratio (%d%% %s interval)", pct, itype), order = FALSE, opts = NULL) +
      ggplot2::scale_y_discrete(labels = wrap)
    attr(pc, "fig_size") <- list(width = "onehalf", height = 30 + 11 * nrow(co))
    figs$coefficients <- with_caption(pc, sprintf("Effects on the odds of apparent survival and of recapture: odds ratio per 1 SD for continuous covariates, versus the reference level for factors (point: %s; bar: %d%% %s interval; log scale). Dashed line: no effect.%s",
                                                  if (bayes) "posterior mean" else "estimate", pct, itype,
                                                  if (any(grepl("^phi_beta_tsm_first_", est$param))) " Survival effects are for residents; first-interval (transient) deviations are in the estimates table." else ""))
  }
  # goodness of fit
  g <- obj$gof
  if (!is.null(g)) {
    t <- g$table
    t$ratio <- ifelse(t$df > 0, t$chi2 / t$df, NA_real_)
    t$lab <- ifelse(t$df > 0, sprintf("p = %.3f, df %d", t$p_value, as.integer(t$df)), "not testable (df 0)")
    t$test <- factor(t$test, levels = rev(t$test))
    t$flag <- ifelse(is.finite(t$p_value) & t$p_value < 0.05, "p < 0.05", "p >= 0.05")
    p <- ggplot2::ggplot(t, ggplot2::aes(ratio, test)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = ratio, yend = test), colour = eco_pal$muted, linewidth = 0.5, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(shape = flag), size = 2.2, colour = eco_pal$ink, na.rm = TRUE) +
      ggplot2::geom_text(ggplot2::aes(x = pmax(ifelse(is.finite(ratio), ratio, 0), 1.02), label = lab), hjust = -0.15, size = 2.3, lineheight = 0.85, colour = eco_pal$ink) +
      ggplot2::scale_shape_manual(values = c("p < 0.05" = 17, "p >= 0.05" = 16), name = NULL, drop = FALSE) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0)), limits = c(0, 1.9 * max(c(1.2, t$ratio), na.rm = TRUE))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = "Chi-square / df (1 = expected under good fit)", y = NULL)
    attr(p, "fig_size") <- list(width = "onehalf", aspect = 0.5)
    figs$gof <- with_caption(p, sprintf("Goodness-of-fit tests of the time-dependent CJS model (U-CARE components via R2ucare%s). Values well above 1 indicate lack of fit: 3.SR = transients, 2.CT = trap dependence. Overall c-hat = %.2f.",
                                        if (length(g$groups) > 1) ", summed over groups" else "", g$chat))
  }
  # captures per occasion (data)
  ch <- data$ch; first <- cjs_first(ch)
  T <- ncol(ch)
  cd <- rbind(data.frame(occasion = data$occasions, n = tabulate(first, T), type = "Newly marked"),
              data.frame(occasion = data$occasions, n = colSums(ch) - tabulate(first, T), type = "Recaptured"))
  cd$occasion <- factor(cd$occasion, levels = data$occasions)
  cd$type <- factor(glab(cd$type), levels = glab(c("Newly marked", "Recaptured")))
  figs$captures <- with_caption(
    ggplot2::ggplot(cd, ggplot2::aes(occasion, n, fill = type)) + ggplot2::geom_col(width = 0.75) +
      ggplot2::labs(x = fig_label("occasion", "Occasion", opts), y = fig_label("n_captured", "Individuals captured", opts), fill = NULL) +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = if (T > 6) 45 else 0)),
    "Individuals captured per occasion: newly marked (first capture) and recaptures of previously marked animals.")
  if (bayes && !is.null(obj$chains)) {
    keep <- colnames(obj$chains[[1]])
    # readable facet titles: "Survival: Sex: Male vs Female" instead of "phi:sex_M"
    plab <- vapply(keep, function(k) {
      pr <- sub(":.*$", "", k); term <- sub("^[^:]*:", "", k)
      head <- if (pr == "phi") "Survival" else "Recapture"
      if (term == "(Intercept)") return(paste(head, "intercept (logit)"))
      if (term == "sd_time") return(paste(head, "SD of year effects (logit)"))
      paste0(head, ": ", cjs_term_label(term, des, data, pr, opts))
    }, "", USE.NAMES = FALSE)
    plab <- make.unique(plab)
    rn <- function(m) { m <- as.matrix(m)[, keep, drop = FALSE]; colnames(m) <- plab; m }
    np <- length(keep); nc <- if (np > 4) 3 else 2
    wrapf <- ggplot2::label_wrap_gen(32)
    pt <- fig_trace(lapply(obj$chains, rn)) +
      ggplot2::facet_wrap(~ind, scales = "free_y", ncol = nc, labeller = ggplot2::label_wrap_gen(24)) +
      ggplot2::scale_x_continuous(n.breaks = 3) + ggplot2::scale_y_continuous(n.breaks = 3)
    attr(pt, "fig_size") <- list(width = "double", height = 20 + 30 * ceiling(np / nc))
    figs$trace <- with_caption(pt, sprintf("MCMC traces of the logit-scale coefficients after burn-in (%d chains, one colour per chain); continuous covariates per 1 SD. Chains should overlap with no trend.", length(obj$chains)))
    pp <- fig_posterior(rn(do.call(rbind, obj$chains))) +
      ggplot2::facet_wrap(~param, scales = "free", ncol = nc, labeller = wrapf) +
      ggplot2::scale_x_continuous(n.breaks = 4) + ggplot2::scale_y_continuous(n.breaks = 3)
    attr(pp, "fig_size") <- list(width = "double", height = 20 + 30 * ceiling(np / nc))
    figs$posterior <- with_caption(pp, "Posterior densities of the logit-scale coefficients (continuous covariates per 1 SD; odds ratio = exp(value)).")
  }
  figs
}

cjs_plots <- function(fit, data) {
  f <- cjs_figures(fit, data, fig_opts())
  f[setdiff(names(f), c("trace", "posterior", "survival_effects_panel"))]
}
