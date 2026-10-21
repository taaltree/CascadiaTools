# Interpretation, figures and methods text for fisheries-growth.

fg_p <- function(p) ifelse(p < 1e-4, "< 0.0001", paste("=", sprintf("%.2g", p)))
fg_row <- function(est, p) est[est$param == p, , drop = FALSE]
fg_fmt <- function(r, d = 0, unit = "") {
  if (!nrow(r)) return("NA")
  sprintf(paste0("%.", d, "f%s (95%% %s %.", d, "f to %.", d, "f)"), r$estimate, unit,
          "CI", r$lower, r$upper)
}

fg_curves <- function(fit, data, ci = 0.95, inner = 0.5) {
  ob <- fit$object; df <- fg_prep_df(data, ob)
  if (ob$engine == "nimble") fg_curve_nimble(ob$growth, df, ob$settings$error, ci, inner)
  else fg_curve_nls(ob$growth, df, ob$settings$error, ci)
}

fg_prep_df <- function(data, ob) {
  df <- data$df; df$group <- droplevels(factor(df$group))
  if (!is.null(data$ref_group) && data$ref_group %in% levels(df$group))
    df$group <- stats::relevel(df$group, ref = as.character(data$ref_group))
  df
}

fg_interpret <- function(fit, data) {
  est <- fit$estimates; ob <- fit$object; df <- fg_prep_df(data, ob)
  gn <- levels(df$group); sg <- fg_safe(gn); G <- length(gn)
  lu <- data$units$length; wu <- data$units$weight
  out <- sprintf("Data: %d aged fish in %d group%s (%s), ages %g-%g years, lengths %g-%g %s. Error model: %s.",
                 nrow(df), G, if (G > 1) "s" else "",
                 paste(sprintf("%s n = %d", gn, as.integer(table(df$group))), collapse = "; "),
                 min(df$age), max(df$age), min(df$length), max(df$length), lu,
                 if (ob$settings$error == "multiplicative")
                   "multiplicative (lognormal), so curves describe the median length at age"
                 else "additive (normal), so curves describe the mean length at age")
  for (g in seq_len(G)) {
    L <- fg_row(est, paste0("Linf_", sg[g])); K <- fg_row(est, paste0("K_", sg[g]))
    A <- fg_row(est, paste0("age_at_90pct_Linf_", sg[g]))
    out <- c(out, sprintf(
      "%s: fish approach an asymptotic length (Linf) of %s and grow toward it at K = %s per year. They reach 90%% of Linf (about %.0f %s) at age %s years.",
      gn[g], fg_fmt(L, 0, paste0(" ", lu)), fg_fmt(K, 3), 0.9 * L$estimate, lu, fg_fmt(A, 1)))
  }
  if (G > 1) {
    dtxt <- c()
    for (g in 2:G) for (p in c("Linf", "K")) {
      r <- fg_row(est, paste0(p, "_diff_", sg[g], "_vs_", sg[1])); if (!nrow(r)) next
      u <- if (p == "Linf") paste0(" ", lu) else " per year"; d <- if (p == "Linf") 0 else 3
      dir <- if (r$lower > 0) "larger" else if (r$upper < 0) "smaller" else NA
      dtxt <- c(dtxt, if (is.na(dir))
        sprintf("%s in %s differs from %s by %s; the interval includes zero, so the data cannot tell these groups apart on %s.",
                p, gn[g], gn[1], fg_fmt(r, d, u), p)
        else sprintf("%s in %s is %s than in %s by %s.", p, gn[g], dir, gn[1],
                     sub("^-", "", fg_fmt(transform(r, estimate = abs(estimate), lower = min(abs(c(lower, upper))),
                                                    upper = max(abs(c(lower, upper)))), d, u))))
    }
    out <- c(out, paste0("Group differences (reference group: ", gn[1], "). ", paste(dtxt, collapse = " ")))
    kim <- ob$growth$kimura
    if (!is.null(kim)) {
      k <- kim[!is.na(kim$p_value), ]
      out <- c(out, paste0("Kimura (1980) likelihood-ratio tests against the model with all parameters group-specific: ",
        paste(sprintf("%s: chi-square = %.1f, df = %d, p %s", tolower(k$hypothesis), k$chisq, k$df,
                      fg_p(k$p_value)), collapse = "; "),
        ". A small p-value means the constrained model fits worse than separate curves, i.e. evidence that the groups differ in that parameter; a large p-value means no detectable difference, not proof that the groups are identical. Linf, K and t0 are strongly correlated, so a difference often shows up in more than one test."))
    }
  }
  if (isTRUE(ob$growth$random)) {
    Lt <- fg_row(est, "Linf_typical"); Kt <- fg_row(est, "K_typical")
    sL <- fg_row(est, "sd_log_Linf"); sK <- fg_row(est, "sd_log_K")
    out <- c(out, sprintf(
      "Hierarchical model: a typical group%s has Linf %s and K %s. Among groups, Linf varies with a standard deviation of %s on the log scale (about %.0f%% from group to group) and K with %s (about %.0f%%).",
      if (isTRUE(ob$growth$hascov)) " (at the average covariate value)" else "",
      fg_fmt(Lt, 0, paste0(" ", lu)), fg_fmt(Kt, 3), fg_fmt(sL, 3), 100 * sL$estimate, fg_fmt(sK, 3), 100 * sK$estimate))
    if (G < 10) out <- c(out, sprintf(
      "With only %d groups the among-group SDs are imprecise and tend to be overestimated (simulation: SD of log Linf +15%% with 8 groups, +25%% with 6); read them through their 95%% intervals rather than as point values.", G))
    if (isTRUE(ob$growth$hascov)) {
      cv <- data$cov_name; sc <- ob$growth$cov_scaling
      for (p in c("Linf", "K")) {
        b <- fg_row(est, paste0(p, "_beta_", cv))
        pc <- 100 * (exp(c(b$estimate, b$lower, b$upper)) - 1)
        dr <- fit$object$growth$derived; pr <- mean(unlist(lapply(dr, function(m) m[, paste0(p, "_beta_", cv)])) > 0)
        out <- c(out, sprintf(
          "Each 1 SD increase in %s (%.2f units) is associated with a %+.1f%% change in %s (95%% CI %+.1f%% to %+.1f%%; posterior probability that the effect is positive %.2f). This is an among-group association across %d groups, not proof that %s causes the difference.",
          cv, sc$scale, pc[1], p, pc[2], pc[3], pr, G, cv))
      }
    }
  }
  lw <- ob$lw
  if (!is.null(lw)) {
    b <- fg_row(est, "lw_b")
    allo <- if (b$lower > 3) "positive allometry: fish become relatively heavier (deeper-bodied) as they grow"
            else if (b$upper < 3) "negative allometry: fish become relatively more slender as they grow"
            else "consistent with isometric growth (b = 3)"
    mw <- est[grepl("^mean_weight_at_", est$param), ]
    out <- c(out, sprintf(
      "Length-weight (%d weighed fish): W = %.3g x L^%.3f, with b %s; this is %s. Typical fish b values are 2.5-3.5 (Froese 2006). A %s %s fish of group %s weighs on average %s (back-transformed with the lognormal bias correction).",
      nrow(lw$df), fg_row(est, "lw_a")$estimate, b$estimate, sprintf("95%% CI %.3f to %.3f", b$lower, b$upper),
      allo, format(data$reference_length), lu, gn[1], fg_fmt(mw, 0, paste0(" ", wu))))
    if (!is.null(lw$slope_test))
      out <- c(out, sprintf("Do groups share one slope b? F = %.2f (df %d, %d), p %s; group-specific b: %s. %s",
        lw$slope_test$F, lw$slope_test$df1, lw$slope_test$df2, fg_p(lw$slope_test$p_value),
        paste(sprintf("%s %.3f", names(lw$slope_test$b_by_group), lw$slope_test$b_by_group), collapse = ", "),
        if (lw$slope_test$p_value < 0.05) "Slopes differ, so the condition ratios below compare weight at the average log length only." else ""))
    if (G > 1) for (g in 2:G) {
      r <- fg_row(est, paste0("condition_ratio_", sg[g], "_vs_", sg[1]))
      pc <- 100 * (c(r$estimate, r$lower, r$upper) - 1)
      out <- c(out, sprintf("At the same length, fish in %s weigh %.1f%% %s than fish in %s (95%% CI %+.1f%% to %+.1f%%)%s.",
                            gn[g], abs(pc[1]), if (pc[1] >= 0) "more" else "less", gn[1], pc[2], pc[3],
                            if (r$lower <= 1 && r$upper >= 1) "; the interval includes no difference" else ""))
    }
    ct <- lw$cond_tab
    out <- c(out, sprintf("%s by group: %s.%s", lw$cond_name,
      paste(sprintf(if (!is.null(data$ws)) "%s %.1f (95%% CI %.1f-%.1f, n = %d)" else "%s %.3f (95%% CI %.3f-%.3f, n = %d)", ct$group, ct$mean, ct$lower, ct$upper, ct$n), collapse = "; "),
      if (!is.null(data$ws)) sprintf(" Wr = 100 means a fish weighs exactly the standard weight for its length (%s); standards built with the regression-line-percentile method represent the 75th percentile of population mean weights, so many healthy populations sit below 100.", data$ws$source)
      else " Kn = 1 is the average weight for that length in these data; Kn compares fish within this data set only."))
  }
  out <- c(out, "Caveats: ages are treated as exact (ageing error is not modelled), and size-selective gear catches the faster-growing fish of young age classes, which inflates length at young ages and biases K down and t0 low. Linf is only as reliable as the number of old fish: see the notes above.")
  if (length(fit$diagnostics$notes)) out <- c(out, paste("Notes:", paste(fit$diagnostics$notes, collapse = " ")))
  out
}

fg_figures <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2"))
  ob <- fit$object; est <- fit$estimates; df <- fg_prep_df(data, ob)
  gn <- levels(df$group); sg <- fg_safe(gn); G <- length(gn)
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  inner <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci
  show_data <- is.null(opts$show_data) || isTRUE(opts$show_data)
  lu <- data$units$length; wu <- data$units$weight
  bayes <- ob$engine == "nimble"
  itype <- if (bayes) "credible" else "confidence"
  pct <- function(p) round(100 * p)
  # Group (entity) labels and encodings: every figure uses the same factor
  # levels with drop = FALSE, so a group keeps its colour, shape and linetype
  # in every figure even when it is missing from one data subset (e.g. no
  # weighed fish). Shape/linetype are the second encoding for greyscale.
  glab <- function(g) vapply(as.character(g), fig_label, "", opts = opts)
  gl <- unname(glab(gn)); gf <- function(g) factor(unname(glab(g)), levels = gl)
  facet <- G > 8
  pal <- if (is.null(opts$palette)) "okabe_ito" else opts$palette
  gshapes <- stats::setNames(rep(c(16, 17, 15, 18, 1, 2, 0, 5), length.out = G), gl)
  gscales <- list(scale_colour_eco(pal, drop = FALSE, limits = gl), scale_fill_eco(pal, drop = FALSE, limits = gl),
                  ggplot2::scale_shape_manual(values = gshapes, drop = FALSE, limits = gl),
                  ggplot2::scale_linetype_discrete(drop = FALSE, limits = gl))
  xl <- fig_label("age", "Age (years)", opts)
  yl <- fig_label("length", sprintf("Length (%s)", lu), opts)
  figs <- list()
  # Size hint for figures that want a non-default aspect. The core only reads
  # width/height from a module hint (it ignores a hinted aspect), so turn the
  # aspect into a height here, honouring per_figure width/aspect overrides.
  fsz <- function(name, width, aspect) {
    u <- opts$per_figure[[name]]
    w <- if (!is.null(u$width)) u$width else width; a <- if (!is.null(u$aspect)) u$aspect else aspect
    wmm <- tryCatch(fig_width_mm(w, opts), error = function(e) NULL)
    if (is.null(wmm) || !is.numeric(wmm)) list(width = width, aspect = aspect) else list(width = width, height = wmm * a)
  }

  # ---- growth_curves --------------------------------------------------------
  cv <- fg_curves(fit, data, ci, inner)
  # Draw each group's curve only over that group's observed ages (no
  # extrapolation past its oldest fish).
  amax <- tapply(df$age, df$group, max); amin <- tapply(df$age, df$group, min)
  cv <- cv[cv$x <= amax[as.character(cv$group)] + 0.25 & cv$x >= amin[as.character(cv$group)] - 0.25, ]
  cv$group <- gf(cv$group)
  dd <- df; dd$group <- gf(dd$group)
  aesg <- if (facet) ggplot2::aes(x, estimate) else ggplot2::aes(x, estimate, colour = group, fill = group, linetype = group)
  p1 <- ggplot2::ggplot(cv, aesg)
  if (show_data)
    p1 <- p1 + (if (facet) ggplot2::geom_point(data = dd, ggplot2::aes(age, length), inherit.aes = FALSE, alpha = 0.35, size = 0.7,
                                               colour = eco_pal$muted, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1))
                else ggplot2::geom_point(data = dd, ggplot2::aes(age, length, colour = group, shape = group), inherit.aes = FALSE,
                                         alpha = 0.35, size = 0.8, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1)))
  p1 <- p1 + do.call(ggplot2::geom_ribbon, c(list(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.22, colour = NA),
                                             if (facet) list(fill = eco_pal$band)))
  if (all(c("lower50", "upper50") %in% names(cv)))
    p1 <- p1 + do.call(ggplot2::geom_ribbon, c(list(ggplot2::aes(ymin = lower50, ymax = upper50), alpha = 0.3, colour = NA),
                                               if (facet) list(fill = eco_pal$band)))
  p1 <- p1 +
    do.call(ggplot2::geom_line, c(list(ggplot2::aes(y = pi_lower), linewidth = 0.3, linetype = 3), if (facet) list(colour = eco_pal$band))) +
    do.call(ggplot2::geom_line, c(list(ggplot2::aes(y = pi_upper), linewidth = 0.3, linetype = 3), if (facet) list(colour = eco_pal$band))) +
    do.call(ggplot2::geom_line, c(list(linewidth = 0.8), if (facet) list(colour = eco_pal$band))) +
    ggplot2::labs(x = xl, y = yl, colour = NULL, fill = NULL, linetype = NULL, shape = NULL)
  p1 <- if (facet) p1 + ggplot2::facet_wrap(~group) + ggplot2::scale_x_continuous(n.breaks = 3) else p1 + gscales
  figs$growth_curves <- with_caption(p1, sprintf(
    "Length at age with the fitted von Bertalanffy curve per group (line: %s length at age, drawn over each group's observed ages), %d%% %s band for the curve%s and %d%% prediction interval for individual fish (dotted). %s",
    if (ob$settings$error == "multiplicative") "median" else "mean", pct(ci), itype,
    if (bayes) sprintf(" (outer shading; inner shading %d%%)", pct(inner)) else " (shaded)", pct(ci),
    if (show_data) "Points are individual fish (ages jittered)." else ""))

  # ---- parameter_comparison -------------------------------------------------
  pars <- c("Linf", "K", if (!isTRUE(ob$growth$random)) "t0", "age_at_90pct_Linf")
  plab <- c(Linf = sprintf("Linf (%s)", lu), K = "K (per year)", t0 = "t0 (years)", age_at_90pct_Linf = "Age at 90% Linf (years)")
  dr <- if (bayes) stack_draws(ob$growth$derived) else NULL
  pc <- do.call(rbind, lapply(pars, function(p) {
    nm <- paste0(p, "_", sg)
    r <- est[match(nm, est$param), c("estimate", "sd", "lower", "upper")]
    if (bayes && all(nm %in% colnames(dr))) {
      q <- function(pp) apply(dr[, nm, drop = FALSE], 2, stats::quantile, pp)
      r$lower <- q((1 - ci) / 2); r$upper <- q(1 - (1 - ci) / 2)
      r$lower50 <- q((1 - inner) / 2); r$upper50 <- q(1 - (1 - inner) / 2)
    } else if (!bayes && abs(ci - 0.95) > 1e-9) {
      # Same Wald construction as the estimates table, at the requested level:
      # Linf and K on the log scale, t0 and age at 90% Linf on the natural scale.
      z <- stats::qnorm(1 - (1 - ci) / 2)
      if (p %in% c("Linf", "K")) { r$lower <- r$estimate * exp(-z * r$sd / r$estimate); r$upper <- r$estimate * exp(z * r$sd / r$estimate) }
      else { r$lower <- r$estimate - z * r$sd; r$upper <- r$estimate + z * r$sd }
    }
    data.frame(par = fig_label(p, plab[[p]], opts), group = gf(gn), r)
  }))
  pc$par <- factor(pc$par, levels = unique(pc$par))
  pc$gy <- factor(pc$group, levels = rev(gl))
  ink <- if (facet) list(colour = eco_pal$ink) else list()
  p2 <- ggplot2::ggplot(pc, if (facet) ggplot2::aes(estimate, gy) else ggplot2::aes(estimate, gy, colour = group, shape = group)) +
    do.call(ggplot2::geom_linerange, c(list(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.6), ink))
  if ("lower50" %in% names(pc)) p2 <- p2 + do.call(ggplot2::geom_linerange, c(list(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.4), ink))
  p2 <- p2 + do.call(ggplot2::geom_point, c(list(size = 1.8), ink)) +
    ggplot2::facet_wrap(~par, scales = "free_x", nrow = 1) + ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_x_continuous(n.breaks = 4) + ggplot2::theme(panel.spacing.x = ggplot2::unit(6, "mm")) +
    (if (!facet) gscales) + ggplot2::guides(colour = "none", shape = "none", fill = "none", linetype = "none")
  figs$parameter_comparison <- with_caption(p2, sprintf(
    "Growth parameters by group (panel titles give units): %s with %d%% %s interval%s. Linf, K and t0 are strongly correlated; compare whole curves and age at 90%% Linf rather than one parameter alone.",
    if (bayes) "posterior median" else "estimate", pct(ci), if (bayes) "credible" else "Wald confidence",
    if (bayes) sprintf(" (thin line) and %d%% interval (thick line)", pct(inner)) else ""))
  attr(figs$parameter_comparison, "fig_size") <- list(width = "double", height = 28 + 7 * G)

  # ---- residuals ------------------------------------------------------------
  fitted_v <- fg_fitted(fit, df)
  rd <- data.frame(age = df$age, fitted = fitted_v, group = gf(df$group),
                   resid = if (ob$settings$error == "multiplicative") log(df$length) - log(fitted_v) else df$length - fitted_v)
  rlab <- fig_label("residual", if (ob$settings$error == "multiplicative") "Residual (log length)" else sprintf("Residual (%s)", lu), opts)
  th <- theme_eco(if (is.null(opts$base_size)) 9 else opts$base_size, opts$font, opts$style, opts$legend)
  r1 <- ggplot2::ggplot(rd, ggplot2::aes(fitted, resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = eco_pal$null) +
    ggplot2::geom_point(if (facet) NULL else ggplot2::aes(colour = group, shape = group), alpha = 0.45, size = 0.8) +
    ggplot2::geom_smooth(ggplot2::aes(group = 1), method = "loess", formula = y ~ x, se = FALSE, colour = eco_pal$ink, linewidth = 0.6) +
    ggplot2::labs(x = fig_label("fitted_length", sprintf("Fitted length (%s)", lu), opts), y = rlab, colour = NULL, shape = NULL) + th +
    (if (!facet) gscales[c(1, 3)])
  # Continuous age axis (one box per age), so gaps between sampled ages stay visible
  r2 <- ggplot2::ggplot(rd, ggplot2::aes(age, resid, group = age)) + ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = eco_pal$null) +
    ggplot2::geom_boxplot(outlier.size = 0.5, fill = NA, colour = eco_pal$ink, linewidth = 0.3, width = 0.6) + ggplot2::labs(x = xl, y = rlab) + th
  figs$residuals <- with_caption(compose_figures(list(r1, r2), ncol = 2),
    sprintf("Growth-model residuals (%s). A: against fitted length; points are individual fish, line is a loess smooth (a trend means the curve shape is wrong). B: by age; boxes show median and quartiles, whiskers 1.5 x IQR (boxes should centre on zero with similar spread; a funnel means the error model is wrong). Ages with one fish appear as a single bar.",
            if (ob$settings$error == "multiplicative") "log observed minus log fitted length" else sprintf("observed minus fitted length, %s", lu)))
  attr(figs$residuals, "fig_size") <- fsz("residuals", "double", 0.42)

  # ---- length_weight and condition -----------------------------------------
  lw <- ob$lw
  if (!is.null(lw)) {
    m <- lw$lm; z <- stats::qt(1 - (1 - ci) / 2, m$df.residual)
    lgn <- levels(lw$df$group)
    nd <- do.call(rbind, lapply(lgn, function(g) {
      r <- range(lw$df$length[lw$df$group == g]); data.frame(logL = log(seq(r[1], r[2], length.out = 60)), group = factor(g, levels = lgn))
    }))
    pr <- stats::predict(m, newdata = nd, se.fit = TRUE); bc <- exp(lw$sigma^2 / 2)
    lp <- data.frame(x = exp(nd$logL), estimate = exp(pr$fit) * bc, lower = exp(pr$fit - z * pr$se.fit) * bc,
                     upper = exp(pr$fit + z * pr$se.fit) * bc, group = gf(nd$group))
    ld <- lw$df; ld$group <- gf(ld$group)
    p4 <- if (facet) {
      ggplot2::ggplot(lp, ggplot2::aes(x, estimate)) +
        (if (show_data) ggplot2::geom_point(data = ld, ggplot2::aes(length, weight), inherit.aes = FALSE, alpha = 0.3, size = 0.5, colour = eco_pal$muted)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.3, fill = eco_pal$band) +
        ggplot2::geom_line(linewidth = 0.7, colour = eco_pal$band) + ggplot2::facet_wrap(~group)
    } else {
      ggplot2::ggplot(lp, ggplot2::aes(x, estimate, colour = group, fill = group, linetype = group)) +
        (if (show_data) ggplot2::geom_point(data = ld, ggplot2::aes(length, weight, colour = group, shape = group), inherit.aes = FALSE, alpha = 0.3, size = 0.7)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) + ggplot2::geom_line(linewidth = 0.7) + gscales
    }
    p4 <- p4 + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::labs(x = yl, y = fig_label("weight", sprintf("Weight (%s)", wu), opts), colour = NULL, fill = NULL, linetype = NULL, shape = NULL)
    figs$length_weight <- with_caption(p4, sprintf(
      "Length-weight relationship on log-log axes: mean weight at length (bias-corrected back-transform; line, drawn over each group's observed lengths) with %d%% confidence band, common slope b = %.3f with group-specific intercepts.%s",
      pct(ci), stats::coef(m)[2], if (show_data) " Points are weighed fish." else ""))
    # Condition: group means with t intervals at opts$ci, from the per-fish index
    cd <- lw$df[is.finite(lw$df$cond_index), ]
    ct <- do.call(rbind, lapply(lgn, function(g) {
      v <- cd$cond_index[cd$group == g]; n <- length(v); mn <- mean(v)
      h <- if (n > 1) stats::qt(1 - (1 - ci) / 2, n - 1) * stats::sd(v) / sqrt(n) else NA
      data.frame(group = g, n = n, mean = mn, lower = mn - h, upper = mn + h)
    }))
    ct <- ct[ct$n > 0, ]; ct$group <- gf(ct$group); cd$group <- gf(cd$group)
    ref <- if (!is.null(data$ws)) 100 else 1
    p5 <- ggplot2::ggplot(ct, ggplot2::aes(group, mean)) +
      ggplot2::geom_hline(yintercept = ref, linetype = 2, colour = eco_pal$null)
    if (show_data) p5 <- p5 + (if (facet) ggplot2::geom_jitter(data = cd, ggplot2::aes(group, cond_index), width = 0.15, height = 0, alpha = 0.25, size = 0.6, colour = eco_pal$muted)
                               else ggplot2::geom_jitter(data = cd, ggplot2::aes(group, cond_index, colour = group, shape = group), width = 0.15, height = 0, alpha = 0.25, size = 0.7))
    p5 <- p5 + ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$ink, linewidth = 0.7) +
      (if (facet) ggplot2::geom_point(size = 2, colour = eco_pal$ink)
       else ggplot2::geom_point(ggplot2::aes(fill = group, shape = group), size = 2.2, colour = eco_pal$ink, stroke = 0.5)) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      (if (!facet) c(gscales[1:2], list(ggplot2::scale_shape_manual(values = stats::setNames(c(21, 24, 22, 23, 21, 24, 22, 23)[(seq_len(G) - 1) %% 8 + 1], gl), drop = FALSE, limits = gl)))) +
      ggplot2::guides(colour = "none", fill = "none", shape = "none") +
      ggplot2::labs(x = NULL, y = fig_label("condition", lw$cond_name, opts))
    # Many groups: stagger the labels (a guide survives the core theme, theme() tweaks do not)
    if (facet) p5 <- p5 + ggplot2::guides(x = ggplot2::guide_axis(n.dodge = 2))
    figs$condition <- with_caption(p5, sprintf(
      "%s by group: mean (large symbol) with %d%% confidence interval (vertical line; may be shorter than the symbol when n is large; n = %s fish).%s Dashed line: %s.", lw$cond_name, pct(ci),
      paste(ct$n, collapse = ", "), if (show_data) " Small points are individual fish." else "",
      if (!is.null(data$ws)) sprintf("Wr = 100, the standard weight for the fish's length (%s)", data$ws$source)
      else "Kn = 1 (average weight for length in these data)"))
  }

  # ---- covariate effects (hierarchical model with a group covariate) ---------
  if (isTRUE(ob$growth$hascov)) {
    gfit <- ob$growth; M <- stack_draws(gfit$samples); sc <- gfit$cov_scaling; cvn <- data$cov_name
    scl <- stats::setNames(list(sc), cvn)
    zs <- seq(min(gfit$x), max(gfit$x), length.out = 60)
    a2 <- (1 - ci) / 2; i2 <- (1 - inner) / 2
    for (p in c("Linf", "K")) {
      mu <- if (p == "Linf") M[, "mu_lL"] else M[, "mu_lK"]; b <- if (p == "Linf") M[, "bL"] else M[, "bK"]
      P <- exp(outer(mu, rep(1, 60)) + outer(b, zs))
      q <- function(pp) apply(P, 2, stats::quantile, pp)
      pred <- data.frame(x = fig_x(zs, cvn, opts, scl), estimate = apply(P, 2, stats::median), lower = q(a2), upper = q(1 - a2),
                         lower50 = q(i2), upper50 = q(1 - i2))
      gp <- est[match(paste0(p, "_", sg), est$param), ]
      obs <- data.frame(x = fig_x(gfit$x, cvn, opts, scl), y = gp$estimate, lo = gp$lower, hi = gp$upper)
      pp <- fig_cond_effect(pred, fig_label(cvn, cvn, opts), fig_label(p, if (p == "Linf") sprintf("Linf (%s)", lu) else "K (per year)", opts), obs = NULL) +
        (if (all(abs(obs$x - round(obs$x)) < 1e-8)) ggplot2::scale_x_continuous(breaks = function(l) { b <- pretty(l); b[b == round(b)] })) +
        ggplot2::geom_pointrange(data = obs, ggplot2::aes(x, y, ymin = lo, ymax = hi),
                                 inherit.aes = FALSE, colour = eco_pal$muted, size = 0.25)
      figs[[paste0(p, "_", cvn)]] <- with_caption(pp, sprintf(
        "%s of a typical group against %s (line: posterior median; bands: %d%% and %d%% credible intervals). Points: group-level posterior means with 95%% credible intervals (shrunk toward the line by the hierarchical model).",
        p, fig_label(cvn, cvn, opts), pct(ci), pct(inner)))
    }
  }

  # ---- MCMC diagnostics (Bayesian) ------------------------------------------
  if (bayes) {
    keep <- intersect(c(paste0("Linf_", sg[1:min(G, 3)]), paste0("K_", sg[1:min(G, 3)]), "Linf_typical", "K_typical", "t0",
                        paste0("t0_", sg[1]), "sd_log_Linf", "sd_log_K", "residual_sd"), colnames(ob$growth$derived[[1]]))
    rn <- function(mm) { mm <- as.matrix(mm)[, keep, drop = FALSE]; colnames(mm) <- vapply(keep, fig_label, "", opts = opts); mm }
    dl <- lapply(ob$growth$derived, rn)
    figs$trace <- with_caption(fig_trace(dl) + ggplot2::scale_x_continuous(n.breaks = 3, labels = function(b) ifelse(b >= 1000, paste0(b / 1000, "k"), b)), "MCMC traces (one colour per chain) for key parameters; chains should overlap with no trends.")
    attr(figs$trace, "fig_size") <- fsz("trace", "double", 0.6)
    figs$posterior <- with_caption(fig_posterior(do.call(rbind, dl)), "Posterior densities of key parameters (all chains combined).")
    attr(figs$posterior, "fig_size") <- fsz("posterior", "double", 0.6)
  }
  figs
}

fg_fitted <- function(fit, df) {
  ob <- fit$object; gi <- as.integer(df$group)
  if (ob$engine == "nimble") {
    M <- stack_draws(ob$growth$samples); G <- nlevels(df$group)
    lL <- colMeans(M[, sprintf("lL[%d]", seq_len(G)), drop = FALSE]); lK <- colMeans(M[, sprintf("lK[%d]", seq_len(G)), drop = FALSE])
    t0 <- colMeans(M[, sprintf("t0[%d]", seq_len(G)), drop = FALSE])
    fg_vb(df$age, exp(lL[gi]), exp(lK[gi]), t0[gi])
  } else {
    th <- t(vapply(ob$growth$gp, function(x) x$theta, numeric(3)))
    fg_vb(df$age, exp(th[gi, 1]), exp(th[gi, 2]), th[gi, 3])
  }
}

fg_methods <- function(fit, data, settings) {
  ob <- fit$object; st <- ob$settings; G <- nlevels(factor(data$df$group))
  err <- if (st$error == "multiplicative") "multiplicative lognormal error (log length ~ Normal(log VB(age), sigma))" else "additive normal error"
  g <- if (st$engine == "nimble") {
    if (st$group_effects == "random") sprintf(
      "We fitted a Bayesian hierarchical von Bertalanffy model (Helser & Lai 2004) in NIMBLE (de Valpine et al. 2017) with %s. Group-level log Linf and log K were drawn from normal distributions%s, with t0 shared among groups. Priors: log Linf mean ~ N(log max length, 1); log K mean ~ N(log 0.3, 1.5); among-group SDs ~ half-normal(0, 0.5); covariate effects ~ N(0, 1); t0 ~ N(0, 2) truncated below the youngest age; sigma ~ U(0, %s). We ran %d chains of %d iterations (%d burn-in) with automated-factor slice samplers on each group's (log Linf, log K) pair and judged convergence by R-hat < 1.01 and bulk ESS >= 400. We report posterior means with 95%% equal-tailed credible intervals, except for the among-group SDs, whose right-skewed posteriors are summarised by the median.",
      err, if (isTRUE(ob$growth$hascov)) sprintf(" whose means depended linearly on the standardised group covariate %s", data$cov_name) else "",
      if (st$error == "multiplicative") "1" else "max length", st$nchains, st$niter, st$nburnin)
    else sprintf(
      "We fitted a Bayesian von Bertalanffy model in NIMBLE (de Valpine et al. 2017) with group-specific Linf, K and t0 and %s. Priors: log Linf ~ N(log max length, 1), log K ~ N(log 0.3, 1.5), t0 ~ N(0, 2) truncated below each group's youngest age, sigma ~ U(0, %s). Groups were compared through posterior differences. We ran %d chains of %d iterations (%d burn-in); convergence was judged by R-hat < 1.01 and bulk ESS >= 400.",
      err, if (st$error == "multiplicative") "1" else "max length", st$nchains, st$niter, st$nburnin)
  } else sprintf(
    "We fitted the typical von Bertalanffy growth function L = Linf(1 - exp(-K(age - t0))) (von Bertalanffy 1938; Beverton & Holt 1957) by nonlinear least squares (R stats::nls, port algorithm) with %s and group-specific parameters, using starting values from FSA::findGrowthStarts (Ogle 2016). Parameters were estimated as log Linf, log K and t0; 95%% Wald intervals were back-transformed, and differences between groups and the age at 90%% of Linf (t0 + ln10/K) used the delta method.%s",
    err, if (G > 1) " Growth was compared among groups with Kimura's (1980) likelihood-ratio tests of the general model against models with a common Linf, K, t0, or all three (-2 ln Lambda = N ln(RSS_H/RSS_general), chi-square with (groups - 1) df per constrained parameter)." else "")
  l <- if (!is.null(ob$lw)) sprintf(
    " The length-weight relationship W = aL^b was fitted by least squares on log-transformed data%s; predicted mean weights were back-transformed with the lognormal bias correction exp(s^2/2) (Sprugel 1983). Condition was summarised as %s.",
    if (G > 1) " with a common slope and group-specific intercepts (the interaction model with group-specific slopes was compared by an F test)" else "",
    if (!is.null(data$ws)) sprintf("relative weight Wr = 100 W/Ws (Wege & Anderson 1978) using the %s", data$ws$source)
    else "relative condition Kn = W / predicted W (Le Cren 1951)") else ""
  paste0(g, l, " Analyses used the ecotools fisheries-growth module (R ", getRversion(), ").")
}
