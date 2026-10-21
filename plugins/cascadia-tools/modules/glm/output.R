# Interpretation sentences, figures and methods text for the glm module.

plural <- function(x) ifelse(grepl("[^aeiou]y$", x), sub("y$", "ies", x), ifelse(grepl("(s|x|ch|sh)$", x), paste0(x, "es"), paste0(x, "s")))

vlab <- function(prep, v) vapply(v, function(x) as.character(prep$spec$labels[[x]] %||% x), "", USE.NAMES = FALSE)

resp_phrase <- function(prep) {
  lab <- prep$spec$response_label %||% prep$response
  if (prep$fam %in% COUNT_FAMILIES) paste("expected", lab)
  else if (prep$fam == "binomial") paste("probability of", lab)
  else paste("mean", lab)
}

response_axis_label <- function(prep) {
  lab <- prep$spec$response_label %||% prep$response
  base <- if (prep$fam %in% COUNT_FAMILIES) paste("Expected", lab)
          else if (prep$fam == "binomial") paste("Probability of", lab)
          else paste("Mean", lab)
  if (!is.null(prep$offset)) base <- sprintf("%s\nper %s %s", base, fmt(prep$offset_per), vlab(prep, prep$offset))
  base
}

# variables that share an interaction term with v
interaction_partners <- function(prep, v) {
  fac <- attr(prep$tt, "factors")
  if (is.null(fac) || !length(fac)) return(character())
  rv <- lapply(rownames(fac), function(r) all.vars(str2lang(r)))
  rows_v <- which(vapply(rv, function(a) v %in% a, logical(1)))
  if (!length(rows_v)) return(character())
  cols <- which(colSums(fac[rows_v, , drop = FALSE] > 0) > 0 & colSums(fac > 0) >= 2)
  out <- character()
  for (cc in cols) for (r in which(fac[, cc] > 0)) out <- c(out, rv[[r]])
  setdiff(unique(out), v)
}

# point + 95% band of expected response for model-scale newdata
pred_band <- function(prep, fo, nd, draws) {
  est <- predict_mu(prep, fo, nd)
  if (is.null(draws)) return(data.frame(estimate = est, lower = est, upper = est))
  mu <- predict_mu(prep, fo, nd, draws)
  data.frame(estimate = est, lower = apply(mu, 1, stats::quantile, 0.025, na.rm = TRUE),
             upper = apply(mu, 1, stats::quantile, 0.975, na.rm = TRUE))
}

held_text <- function(prep, except) {
  others <- setdiff(prep$vars, except)
  parts <- character()
  if (any(others %in% prep$num_vars)) parts <- c(parts, "other continuous predictors at their mean")
  for (v in intersect(others, prep$fac_vars)) parts <- c(parts, sprintf("%s = %s", vlab(prep, v), prep$xlev[[v]][1]))
  if (length(prep$groups)) parts <- c(parts, paste("typical", paste(prep$groups, collapse = "/")))
  paste(parts, collapse = "; ")
}

# Figure names double as file names in core export_figures(); '%' and '/' break the
# saved file path (ggsave treats '%' as a format string), so spell them out.
fig_name <- function(x) gsub("\\( ", "(", gsub("\\s+", " ", gsub("/", " per ", gsub("%", " percent", x, fixed = TRUE), fixed = TRUE)))

# Several curves (interactions). Drawn here rather than with core cond_effect_plot(),
# whose grouped branch passes colour/fill = NULL and loses the grouping in ggplot2 4.0.
grouped_effect_plot <- function(pred, xlab, ylab, observed_x) {
  ggplot2::ggplot(pred, ggplot2::aes(x, estimate, colour = group, fill = group, group = group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_rug(data = data.frame(x = observed_x), ggplot2::aes(x = x), inherit.aes = FALSE, alpha = 0.4, sides = "b") +
    ggplot2::scale_colour_manual(values = c("#2f6f4f", "#c47f1a", "#3b6ea8", "#9c3d54", "#6b5b95", "#7a7a7a")[(seq_along(levels(pred$group)) - 1) %% 6 + 1]) +
    ggplot2::scale_fill_manual(values = c("#2f6f4f", "#c47f1a", "#3b6ea8", "#9c3d54", "#6b5b95", "#7a7a7a")[(seq_along(levels(pred$group)) - 1) %% 6 + 1]) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")
}

glm_plots <- function(fit, data) {
  eco_require(c("ggplot2", "patchwork"))
  prep <- fit$prep; fo <- fit$fo; st <- fit$settings
  draws <- coef_draws(fo, st$n_draws, st$seed)
  ylab <- response_axis_label(prep)
  out <- list()
  for (v in prep$vars) {
    part <- interaction_partners(prep, v)[1]
    if (is.na(part)) part <- NULL
    base <- prep$ref_raw
    if (v %in% prep$num_vars) {
      xr <- prep$raw[[v]]
      grid <- seq(min(xr), max(xr), length.out = 100)
      nd0 <- base[rep(1, 100), , drop = FALSE]
      nd0[[v]] <- grid
      gvals <- list(NULL)
      if (!is.null(part)) {
        if (part %in% prep$num_vars) {
          m <- mean(prep$raw[[part]]); s <- stats::sd(prep$raw[[part]])
          gvals <- stats::setNames(as.list(m + c(-1, 0, 1) * s),
                                   sprintf("%s = %s (%s)", vlab(prep, part), fmt(m + c(-1, 0, 1) * s), c("mean - 1 SD", "mean", "mean + 1 SD")))
        } else gvals <- stats::setNames(as.list(prep$xlev[[part]]), sprintf("%s = %s", vlab(prep, part), prep$xlev[[part]]))
      }
      pred <- do.call(rbind, lapply(seq_along(gvals), function(i) {
        nd <- nd0
        if (!is.null(gvals[[i]])) nd[[part]] <- if (part %in% prep$num_vars) gvals[[i]] else factor(gvals[[i]], levels = prep$xlev[[part]])
        nd <- to_model_scale(prep, nd)
        pb <- pred_band(prep, fo, nd, draws)
        d <- data.frame(x = grid, pb)
        if (!is.null(gvals[[i]])) d$group <- names(gvals)[i]
        d
      }))
      if ("group" %in% names(pred)) pred$group <- factor(pred$group, levels = names(gvals))
      p <- if ("group" %in% names(pred)) grouped_effect_plot(pred, vlab(prep, v), ylab, xr) else cond_effect_plot(pred, vlab(prep, v), ylab, observed_x = xr)
      p <- p +
        ggplot2::labs(subtitle = paste("Held at:", held_text(prep, c(v, part)))) +
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"))
      if (prep$fam == "binomial") p <- p + ggplot2::coord_cartesian(ylim = c(0, NA))
      out[[fig_name(paste("Effect of", vlab(prep, v)))]] <- p
    } else {
      lv <- prep$xlev[[v]]
      gvals <- list(NULL)
      if (!is.null(part) && part %in% prep$fac_vars)
        gvals <- stats::setNames(as.list(prep$xlev[[part]]), sprintf("%s = %s", vlab(prep, part), prep$xlev[[part]]))
      pred <- do.call(rbind, lapply(seq_along(gvals), function(i) {
        nd <- base[rep(1, length(lv)), , drop = FALSE]
        nd[[v]] <- factor(lv, levels = lv)
        if (!is.null(gvals[[i]])) nd[[part]] <- factor(gvals[[i]], levels = prep$xlev[[part]])
        d <- data.frame(level = factor(lv, levels = lv), pred_band(prep, fo, to_model_scale(prep, nd), draws))
        d$group <- if (!is.null(gvals[[i]])) names(gvals)[i] else "all"
        d
      }))
      pd <- ggplot2::position_dodge(width = if (length(gvals) > 1) 0.4 else 0)
      p <- ggplot2::ggplot(pred, ggplot2::aes(level, estimate, colour = group)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.8, position = pd) +
        ggplot2::geom_point(size = 3, position = pd) +
        ggplot2::labs(x = vlab(prep, v), y = ylab, colour = NULL,
                      subtitle = paste("Held at:", held_text(prep, c(v, if (length(gvals) > 1) part)))) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"))
      if (length(gvals) == 1) p <- p + ggplot2::scale_colour_manual(values = "#2f6f4f", guide = "none")
      out[[fig_name(paste("Effect of", vlab(prep, v)))]] <- p
    }
  }
  # effect-size (coefficient) plot
  est <- fit$estimates
  rp <- ratio_prefix(prep$fam, prep$link)
  pref <- if (is.na(rp)) "beta_" else paste0(rp, "_")
  labs_nonint <- prep$labels_coef[names(prep$labels_coef) != "(Intercept)"]
  es <- est[est$param %in% paste0(pref, labs_nonint), ]
  if (nrow(es)) {
    es$term <- factor(es$per, levels = rev(unique(es$per)))
    nullv <- if (is.na(rp)) 0 else 1
    p <- ggplot2::ggplot(es, ggplot2::aes(estimate, term)) +
      ggplot2::geom_vline(xintercept = nullv, linetype = 2, colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), width = 0.2, colour = "#2f6f4f", orientation = "y") +
      ggplot2::geom_point(size = 2.6, colour = "#2f6f4f") +
      ggplot2::labs(x = if (is.na(rp)) "Difference in mean (95% CI)" else paste0(tools::toTitleCase(ratio_words(prep$fam, prep$link)), " (95% CI, log axis)"),
                    y = NULL, subtitle = sprintf("Dashed line = no effect (%s)", nullv)) +
      ggplot2::theme_minimal(base_size = 12)
    if (!is.na(rp)) p <- p + ggplot2::scale_x_log10()
    out[["Effect sizes"]] <- p
  }
  dh <- fit$dharma
  if (!is.null(dh) && is.null(dh$error)) {
    r <- sort(dh$resid); n <- length(r)
    q <- ggplot2::ggplot(data.frame(e = (seq_len(n) - 0.5) / n, o = r), ggplot2::aes(e, o)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#c0392b") +
      ggplot2::geom_point(size = 1, alpha = 0.6) +
      ggplot2::labs(x = "Expected (uniform)", y = "Observed scaled residual",
                    title = "QQ plot of scaled residuals",
                    subtitle = sprintf("Uniformity p = %s; points should follow the line", fmt(dh$p[["uniformity"]], 2))) +
      ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10.5), plot.subtitle = ggplot2::element_text(size = 8.5))
    rk <- rank(dh$pred, ties.method = "average") / length(dh$pred)
    rv <- ggplot2::ggplot(data.frame(x = rk, y = dh$resid), ggplot2::aes(x, y)) +
      ggplot2::geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 2, colour = "grey55") +
      ggplot2::geom_point(size = 1, alpha = 0.5) +
      ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), colour = "#c0392b", se = FALSE) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = "Model prediction (rank-transformed)", y = "Scaled residual",
                    title = "Residuals vs predicted",
                    subtitle = sprintf("Dispersion p = %s; red line should be flat", fmt(dh$p[["dispersion"]], 2))) +
      ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10.5), plot.subtitle = ggplot2::element_text(size = 8.5))
    out[["Residual checks (DHARMa)"]] <- q + rv
  }
  out
}

# ---- interpretation ------------------------------------------------------------
effect_sentence <- function(prep, b, lo, hi, what, resp) {
  rp <- ratio_prefix(prep$fam, prep$link)
  if (identical(rp, "odds_ratio"))
    return(sprintf("%s multiplies the odds %s by %s (95%% CI %s to %s)", what,
                   if (prep$fam == "beta") sprintf("of the %s, mean/(1 - mean),", resp) else paste("of", sub("^probability of ", "", resp)),
                   fmt(exp(b)), fmt(exp(lo)), fmt(exp(hi))))
  if (!is.na(rp)) {
    pct <- if (b >= 0) sprintf("%s%% higher", fmt(100 * (exp(b) - 1), 2)) else sprintf("%s%% lower", fmt(100 * (1 - exp(b)), 2))
    art <- if (grepl("^(8|11|18)([^0-9]|$)|^8[0-9]", pct)) "an" else "a"   # "an 11%", "an 8.5%", "an 80%"
    return(sprintf("%s is associated with %s %s %s (%s %s, 95%% CI %s to %s)", what, art, pct, resp,
                   ratio_words(prep$fam, prep$link), fmt(exp(b)), fmt(exp(lo)), fmt(exp(hi))))
  }
  if (prep$link == "identity")
    return(sprintf("%s is associated with a change of %s (95%% CI %s to %s) in the %s", what, fmt(b), fmt(lo), fmt(hi), resp))
  sprintf("%s has coefficient %s (95%% CI %s to %s) on the %s-link scale", what, fmt(b), fmt(lo), fmt(hi), prep$link)
}

certainty_sentence <- function(prep, lo, hi) {
  if (lo > 0 || hi < 0) {
    return(sprintf("The 95%% interval excludes no effect, so the data clearly indicate a %s association.",
                   if (lo > 0) "positive" else "negative"))
  }
  rp <- ratio_prefix(prep$fam, prep$link)
  rng <- if (!is.na(rp)) sprintf("from %s%% lower to %s%% higher", fmt(100 * (1 - exp(lo)), 2), fmt(100 * (exp(hi) - 1), 2))
         else sprintf("from %s to %s", fmt(lo), fmt(hi))
  sprintf("The 95%% interval includes no effect: the data are compatible with effects %s%s, including none.",
          rng, if (identical(rp, "odds_ratio")) " in the odds" else "")
}

# Plain-language reading of a two-way interaction coefficient
describe_term <- function(prep, p) {
  if (p %in% prep$num_vars) return(list(type = "num", var = p, txt = if (p %in% prep$scaled)
    sprintf("1 SD of %s (%s)", vlab(prep, p), fmt(prep$scaling[[p]][["sd"]])) else sprintf("1 unit of %s", vlab(prep, p))))
  for (v in prep$fac_vars) if (startsWith(p, v) && substring(p, nchar(v) + 1) %in% prep$xlev[[v]])
    return(list(type = "fac", var = v, lvl = substring(p, nchar(v) + 1), ref = prep$xlev[[v]][1]))
  list(type = "other", txt = p)
}

interaction_sentence <- function(prep, cn, b, lo, hi) {
  parts <- strsplit(cn, ":", fixed = TRUE)[[1]]
  tail_txt <- paste0(" (see the conditional-effects figure with several curves). ",
                     sub("association", "interaction", certainty_sentence(prep, lo, hi)))
  if (length(parts) != 2) return(sprintf("Higher-order interaction %s: coefficient %s (95%% CI %s to %s) on the %s scale; read it from the figures%s",
                                         cn, fmt(b), fmt(lo), fmt(hi), prep$link, tail_txt))
  d <- lapply(parts, describe_term, prep = prep)
  if (d[[1]]$type != "num" && d[[2]]$type == "num") d <- rev(d)
  rp <- ratio_prefix(prep$fam, prep$link)
  resp <- resp_phrase(prep)
  mult <- !is.na(rp)
  eff_of <- function(x) if (x$type == "num") sprintf("the %s per %s", if (mult) ratio_words(prep$fam, prep$link) else paste0("change in the ", resp), x$txt)
                        else if (x$type == "fac") sprintf("the %s for %s = %s vs %s", if (mult) ratio_words(prep$fam, prep$link) else paste0("difference in the ", resp), vlab(prep, x$var), x$lvl, x$ref)
                        else x$txt
  cond <- if (d[[2]]$type == "num") sprintf("for each %s increase", d[[2]]$txt)
          else if (d[[2]]$type == "fac") sprintf("when %s = %s rather than %s", vlab(prep, d[[2]]$var), d[[2]]$lvl, d[[2]]$ref)
          else paste("with", d[[2]]$txt)
  core <- if (mult) sprintf("%s is multiplied by %s (95%% CI %s to %s)", eff_of(d[[1]]), fmt(exp(b)), fmt(exp(lo)), fmt(exp(hi)))
          else sprintf("%s changes by %s (95%% CI %s to %s)", eff_of(d[[1]]), fmt(b), fmt(lo), fmt(hi))
  first <- d[[1]]; nm <- if (first$type == "fac") vlab(prep, first$var) else if (first$type == "num") vlab(prep, first$var) else first$txt
  other <- if (d[[2]]$type %in% c("num", "fac")) vlab(prep, d[[2]]$var) else d[[2]]$txt
  paste0(sprintf("Interaction of %s with %s: %s, %s", nm, other, cond, core), tail_txt)
}

glm_interpret <- function(fit, data) {
  prep <- fit$prep; fo <- fit$fo; est <- fit$estimates
  z <- stats::qnorm(0.975)
  resp <- resp_phrase(prep)
  out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge cleanly (see Diagnostics). The numbers below are not reliable until this is fixed.")
  g <- prep$groups
  out <- c(out, sprintf("Model: %s GLM%s (%s link) for %s, fitted to %d observations%s%s.",
    switch(prep$fam, gaussian = "Gaussian (normal)", binomial = if (is.null(prep$trials)) "Binomial (0/1)" else "Binomial (successes out of trials)",
           poisson = "Poisson", nbinom2 = "Negative-binomial", beta = "Beta", Gamma = "Gamma",
           zipoisson = "Zero-inflated Poisson", zinb = "Zero-inflated negative-binomial"),
    if (length(g)) "M" else "", prep$link, prep$spec$response_label %||% prep$response, fo$nobs,
    if (length(g)) paste0(" from ", paste(sprintf("%d %s groups", vapply(g, function(x) nlevels(prep$df[[x]]), 1L), g), collapse = " and ")) else "",
    if (length(prep$scaled)) paste0("; continuous predictors were standardised (", paste(vlab(prep, prep$scaled), collapse = ", "), ") so effects below are per 1 standard deviation (SD) and per the stated units") else ""))
  r0 <- est[est$param == "mean_response_at_reference", ]
  if (nrow(r0)) out <- c(out, sprintf("Baseline: at the reference point (%s), the %s is %s (95%% CI %s to %s).",
                                      r0$per, resp, fmt(r0$estimate), fmt(r0$lower), fmt(r0$upper)))
  ref <- prep$ref_raw
  pred_at <- function(v, val) { nd <- ref; nd[[v]] <- val; predict_mu(prep, fo, to_model_scale(prep, nd)) }
  se_all <- sqrt(pmax(diag(fo$V), 0))
  for (j in seq_along(fo$beta)) {
    cn <- names(fo$beta)[j]; if (cn == "(Intercept)") next
    b <- fo$beta[[j]]; s <- se_all[[j]]; lo <- b - z * s; hi <- b + z * s
    if (cn %in% prep$num_vars) {
      v <- cn; part <- interaction_partners(prep, v)
      what <- if (v %in% prep$scaled) sprintf("Each 1-SD increase in %s (1 SD = %s)", vlab(prep, v), fmt(prep$scaling[[v]][["sd"]]))
              else sprintf("Each 1-unit increase in %s", vlab(prep, v))
      txt <- paste0(effect_sentence(prep, b, lo, hi, what, resp), ".")
      u <- unit_value(prep, v)
      if (!is.null(u)) {
        k <- if (v %in% prep$scaled) u / prep$scaling[[v]][["sd"]] else u
        txt <- paste(txt, paste0(effect_sentence(prep, b * k, lo * k, hi * k, paste("Per", unit_text(prep, v), "this"), resp), "."))
      }
      xr <- range(prep$raw[[v]])
      m1 <- pred_at(v, xr[1]); m2 <- pred_at(v, xr[2])
      txt <- paste(txt, sprintf("Across the observed range of %s (%s to %s), with %s, the %s goes from %s to %s.",
                                vlab(prep, v), fmt(xr[1]), fmt(xr[2]), held_text(prep, v), resp, fmt(m1), fmt(m2)))
      if (length(part)) txt <- paste(txt, sprintf("Because %s is in an interaction with %s, this is its effect when %s is at its mean or reference level.",
                                                  vlab(prep, v), paste(vlab(prep, part), collapse = ", "), paste(vlab(prep, part), collapse = "/")))
      out <- c(out, paste(txt, certainty_sentence(prep, lo, hi)))
    } else if (!grepl(":", cn) && any(vapply(prep$fac_vars, function(v) startsWith(cn, v) && substring(cn, nchar(v) + 1) %in% prep$xlev[[v]], logical(1)))) {
      v <- prep$fac_vars[vapply(prep$fac_vars, function(v) startsWith(cn, v) && substring(cn, nchar(v) + 1) %in% prep$xlev[[v]], logical(1))][1]
      lvl <- substring(cn, nchar(v) + 1); r <- prep$xlev[[v]][1]
      what <- sprintf("Compared with %s = %s (reference), %s = %s", vlab(prep, v), r, vlab(prep, v), lvl)
      m0 <- pred_at(v, factor(r, levels = prep$xlev[[v]])); m1 <- pred_at(v, factor(lvl, levels = prep$xlev[[v]]))
      txt <- paste0(effect_sentence(prep, b, lo, hi, what, resp), sprintf(": %s vs %s.", fmt(m1), fmt(m0)))
      if (length(interaction_partners(prep, v))) txt <- paste(txt, sprintf("(This contrast is for the mean/reference level of %s, because of an interaction.)",
                                                                       paste(vlab(prep, interaction_partners(prep, v)), collapse = ", ")))
      out <- c(out, paste(txt, certainty_sentence(prep, lo, hi)))
    } else if (grepl(":", cn)) {
      out <- c(out, interaction_sentence(prep, cn, b, lo, hi))
    } else {
      out <- c(out, sprintf("Term %s (coefficient %s, 95%% CI %s to %s on the link scale) is a transformed or polynomial term; read its shape from the conditional-effects figure rather than from the coefficient.",
                            cn, fmt(b), fmt(lo), fmt(hi)))
    }
  }
  if (!is.null(fo$re)) for (i in seq_len(nrow(fo$re))) {
    r <- est[est$param == paste0("re_sd_", clean_name(fo$re$group[i])), ]
    gname <- fo$re$group[i]
    if (r$estimate < 1e-3 || !is.finite(r$upper)) {
      out <- c(out, sprintf("Between-%s variation was estimated at essentially zero (SD %s): %s do not differ beyond what the predictors explain, so the random intercept adds little (the fit is 'singular').", gname, fmt(r$estimate), plural(gname)))
    } else if (prep$link %in% c("log", "logit")) {
      f <- exp(z * r$estimate)
      out <- c(out, sprintf("%s differ beyond the predictors: SD of %s intercepts = %s (95%% CI %s to %s) on the %s scale. The middle 95%% of %s have %s between %s and %s times that of a typical %s.",
                            tools::toTitleCase(plural(gname)), gname, fmt(r$estimate), fmt(r$lower), fmt(r$upper), prep$link, plural(gname),
                            if (prep$link == "logit") "odds" else "expected values", fmt(1 / f), fmt(f), gname))
    } else {
      out <- c(out, sprintf("%s differ beyond the predictors: SD of %s intercepts = %s (95%% CI %s to %s) in response units; 95%% of %s lie within +/- %s of a typical %s.",
                            tools::toTitleCase(plural(gname)), gname, fmt(r$estimate), fmt(r$lower), fmt(r$upper), plural(gname), fmt(z * r$estimate), gname))
    }
    out[length(out)] <- paste(out[length(out)], "Predictions in the figures are for a typical group (random effect = 0).")
  }
  dn <- disp_name(prep$fam)
  d <- est[est$param %in% dn, ]
  if (nrow(d)) {
    mresp <- mean(prep$raw[[prep$response]])
    out <- c(out, switch(prep$fam,
      gaussian = sprintf("Residual SD = %s (95%% CI %s to %s): individual observations scatter around the predicted mean by about this much.", fmt(d$estimate), fmt(d$lower), fmt(d$upper)),
      nbinom2 = , zinb = sprintf("Counts are more variable than Poisson: NB size theta = %s (95%% CI %s to %s); variance = mean + mean^2/theta, so at the average observed count (%s) the variance is about %s times the mean.",
                                 fmt(d$estimate), fmt(d$lower), fmt(d$upper), fmt(mresp), fmt(1 + mresp / d$estimate)),
      beta = sprintf("Beta precision phi = %s (95%% CI %s to %s): variance = mean(1 - mean)/(1 + phi); at a mean proportion of 0.5 the SD is about %s.",
                     fmt(d$estimate), fmt(d$lower), fmt(d$upper), fmt(sqrt(0.25 / (1 + d$estimate)))),
      Gamma = sprintf("Coefficient of variation = %s (95%% CI %s to %s): at any predicted mean, observations typically scatter by about %s%% of that mean.",
                      fmt(d$estimate), fmt(d$lower), fmt(d$upper), fmt(100 * d$estimate, 2))))
  }
  zr <- est[est$param == "zero_inflation_prob", ]
  if (nrow(zr)) out <- c(out, sprintf("An estimated %s%% of observations (95%% CI %s%% to %s%%) are excess ('structural') zeros beyond those the count distribution produces, e.g. visits when the species could not be present or detected. The expected values above include these zeros.",
                                      fmt(100 * zr$estimate, 2), fmt(100 * zr$lower, 2), fmt(100 * zr$upper, 2)))
  r2 <- fit$diagnostics$r2
  if (length(r2)) {
    if (all(c("R2_marginal", "R2_conditional") %in% names(r2)))
      out <- c(out, sprintf("Variance explained (Nakagawa R2): predictors alone %s%% (marginal); predictors plus %s differences %s%% (conditional).",
                            fmt(100 * r2[["R2_marginal"]], 2), paste(g, collapse = "/"), fmt(100 * r2[["R2_conditional"]], 2)))
    else out <- c(out, sprintf("Variance explained (%s): %s%%.", names(r2)[1], fmt(100 * r2[[1]], 2)))
  }
  dh <- fit$dharma
  if (!is.null(dh) && is.null(dh$error)) {
    bad <- names(dh$p)[dh$p < 0.05]
    out <- c(out, if (length(bad)) sprintf("Model checks: the DHARMa residual tests flagged %s (p < 0.05); see Diagnostics before relying on the intervals.", paste(bad, collapse = ", "))
                  else "Model checks: DHARMa simulated residuals showed no detectable problems with distribution, dispersion, outliers or zeros (all p >= 0.05).")
  }
  if (prep$fam %in% COUNT_FAMILIES || (prep$fam == "binomial" && is.null(prep$trials)))
    out <- c(out, sprintf("These are effects on the %s as recorded, which combines how many/whether animals are present with how likely they are to be detected or recorded. If detectability changes with a predictor (e.g. dense canopy hides birds, observer or survey conditions), its effect on %s is not purely an effect on abundance or occurrence; repeated-survey designs can separate the two (occupancy or N-mixture modules).",
                          if (prep$fam == "binomial") "probability of a recorded 1" else "expected count",
                          if (prep$fam == "binomial") "the recorded outcome" else "counts"))
  out <- c(out, "These are associations estimated from the data; unless the predictors were experimentally manipulated, they do not show that a predictor causes the change in the response (other correlated factors could be responsible). Predictions outside the observed range of a predictor are extrapolation.")
  out
}

glm_methods <- function(fit, data, settings) {
  prep <- fit$prep; fo <- fit$fo
  famtxt <- switch(prep$fam, gaussian = "Gaussian", binomial = "binomial", poisson = "Poisson",
                   nbinom2 = "negative binomial (quadratic variance, NB2: Var = mu + mu^2/theta)",
                   beta = "beta (mean-precision parameterisation; Ferrari & Cribari-Neto 2004)",
                   Gamma = "Gamma", zipoisson = "zero-inflated Poisson (constant zero-inflation probability)",
                   zinb = "zero-inflated negative binomial (NB2, constant zero-inflation probability)")
  paste0(
    sprintf("We modelled %s with a %s generalised linear %smodel with a %s link", prep$spec$response_label %||% prep$response, famtxt,
            if (length(prep$groups)) "mixed " else "", prep$link),
    if (!is.null(prep$trials)) sprintf(" (successes out of %s trials)", prep$trials) else "",
    sprintf(". Fixed effects: %s.", if (prep$rhs_fixed == "1") "intercept only" else prep$rhs_fixed),
    if (length(prep$groups)) sprintf(" Random intercepts for %s accounted for repeated observations within groups.", paste(prep$groups, collapse = " and ")) else "",
    if (!is.null(prep$offset)) sprintf(" log(%s) was included as an offset so that the response is modelled per unit of effort.", prep$offset) else "",
    if (length(prep$scaled)) sprintf(" Continuous predictors (%s) were centred and divided by their standard deviation before fitting (Schielzeth 2010).", paste(prep$scaled, collapse = ", ")) else "",
    if (fo$engine == "glmmTMB") sprintf(" The model was fitted by %s in glmmTMB %s (Brooks et al. 2017) in R %s.",
                                      if (isTRUE(fo$REML) && isTRUE(fo$reml_ml_fixef)) "restricted maximum likelihood for the random-effect SD(s) and dispersion (REML: fixed effects integrated out together with the random effects by the Laplace approximation, which reduces the downward bias of random-effect SDs), with the fixed effects then estimated by maximum likelihood (Laplace approximation) with those variance parameters held at their REML values (the analogue of REML in linear mixed models; the REML mode of the fixed effects is biased for non-Gaussian data)"
                                      else if (isTRUE(fo$REML)) "restricted maximum likelihood (REML: fixed effects integrated out together with the random effects by the Laplace approximation, which reduces the downward bias of random-effect SDs)"
                                      else if (length(prep$groups)) "maximum likelihood (Laplace approximation for the random effects)" else "maximum likelihood",
                                      as.character(utils::packageVersion("glmmTMB")), paste(R.version$major, R.version$minor, sep = "."))
    else sprintf(" The model was fitted by maximum likelihood with %s in R %s.", if (prep$fam == "nbinom2") "MASS::glm.nb" else "stats::glm", paste(R.version$major, R.version$minor, sep = ".")),
    " 95% confidence intervals are Wald intervals on the link scale, back-transformed to odds/rate ratios or to the response scale; intervals for SDs and dispersion are Wald intervals on the log scale.",
    if (isTRUE(settings$diagnostics)) sprintf(" Model fit was checked with simulated quantile residuals (DHARMa; Hartig 2026; %d simulations): tests of uniformity, dispersion, outliers%s. Variance explained is reported as R2 (Nakagawa et al. 2017, via performance; Luedecke et al. 2021).",
                                              settings$n_sim_resid, if (prep$fam %in% COUNT_FAMILIES) " and zero inflation" else "") else "",
    " Conditional-effects curves show the expected response across the observed range of each predictor with other predictors at their mean or reference level, with 95% intervals from 2000 draws from the sampling distribution of the estimates.",
    " The module was validated by simulation (planted-truth recovery of all reported estimands; see VALIDATION.md).")
}
