# Publication figure suite for the glm module (docs/FIGURES.md, skills/eco-figures).
#
# glm_figures(fit, data, opts) returns a named list of ggplots, each with a
# caption, built from the core constructors in core/R/figures.R:
#   effect_<var>     conditional effect of every predictor on the response scale
#                    (factors via fig_cat_effect; interactions: one curve per
#                    level of the moderator, or at mean - 1 SD / mean / mean + 1 SD)
#   effects_panel    all effect_* figures composed with A/B/C tags
#   coefficients     forest plot of odds / rate ratios (or differences for an
#                    identity link) per 1 SD and per factor contrast
#   residuals        DHARMa scaled residuals: QQ plot + residuals vs prediction
#   random_effects   caterpillar plot of the random intercepts (mixed models only)
# The core applies theme, palette and journal sizes; nothing is hard-coded here.
# Read-only with respect to the fit: nothing here changes estimation.

# centre/scale of the standardised predictors in the form fig_x() expects
glm_fig_scaling <- function(prep)
  lapply(prep$scaling, function(s) list(center = s[["mean"]], scale = s[["sd"]]))

# Readable breaks for a log axis spanning [lo, hi] (fig_coef's default log10 scale
# can leave a whole side of the null line without a tick).
# panel_mm (optional): approximate panel width; the densest tick set whose labels
# do not collide at that width is used (at single-column width with long row
# labels 0.5 and 0.7 would otherwise overprint).
glm_log_breaks <- function(lo, hi, panel_mm = NULL, char_mm = 1.6) {
  sets <- list(c(1, 1.5, 2, 3, 5, 7), c(1, 2, 3, 5), c(1, 2, 5), c(1, 3), 1)
  pick <- function(m) {
    cand <- as.vector(outer(m, 10^(-4:4)))
    b <- sort(cand[cand >= lo & cand <= hi])
    while (length(b) > 7) b <- b[c(TRUE, FALSE)]
    if (length(b) < 3) b <- signif(exp(seq(log(lo), log(hi), length.out = 4)), 2)
    unique(sort(c(b, 1)[c(b, 1) >= lo & c(b, 1) <= hi]))
  }
  if (is.null(panel_mm)) return(pick(sets[[1]]))
  mm_per <- panel_mm / max(log10(hi) - log10(lo), 1e-6)
  for (m in sets) {
    b <- pick(m)
    if (length(b) < 2) return(b)
    w <- (nchar(format(b, drop0trailing = TRUE, trim = TRUE)) + 1.5) * char_mm
    gap <- diff(log10(b)) * mm_per
    if (all(gap >= (head(w, -1) + tail(w, -1)) / 2)) return(b)
  }
  pick(sets[[length(sets)]])
}

# Approximate width (mm) of one character of axis text at the figure's base size
glm_fig_char_mm <- function(opts) (opts$base_size %||% 9) * 0.9 * 0.3528 * 0.55

# Final width (mm) of figure `id`: per_figure width, else the module's hint, else the global width
glm_fig_width_mm <- function(id, hint, opts) {
  w <- opts$per_figure[[id]]$width %||% hint %||% opts$width %||% "single"
  tryCatch(if (is.null(opts$.journal)) fig_widths[[w]] %||% as.numeric(w) else fig_width_mm(w, opts), error = function(e) 84)
}

# Wrap an axis title that is wider than the panel
glm_fig_wrap <- function(txt, panel_mm, opts) {
  if (grepl("\n", txt)) return(txt)
  n <- max(12, floor(panel_mm / (glm_fig_char_mm(opts) / 0.9)))
  if (nchar(txt) <= n) return(txt)
  # keep short bracketed groups such as "(90% CI)" on one line
  prot <- gsub("\\(([^() ]{1,6}) ([^() ]{1,6})\\)", "(\\1\u00a0\\2)", txt)
  gsub("\u00a0", " ", paste(strwrap(prot, n), collapse = "\n"), fixed = TRUE)
}

# User axis renames given through per_figure: {<id>: {axis: {x: ...}}} also reach the caption
glm_fig_axis <- function(opts, id, k, default) opts$per_figure[[id]]$axis[[k]] %||% default

# greyscale: true -> the single-series colour (eco_pal$band, green) becomes grey;
# the core greyscale option only swaps the discrete palette.
glm_fig_grey <- function(p, opts) {
  if (!isTRUE(opts$greyscale) || !inherits(p, "ggplot")) return(p)
  greens <- toupper(c(eco_pal$band, eco_pal$accent))
  for (i in seq_along(p$layers)) {
    ap <- p$layers[[i]]$aes_params
    for (k in intersect(names(ap), c("colour", "fill")))
      if (is.character(ap[[k]]) && toupper(ap[[k]]) %in% greens)
        p$layers[[i]]$aes_params[[k]] <- if (k == "fill") "#4d4d4d" else "#333333"
  }
  p
}
glm_ucfirst <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))

# Labels: the user's figures: labels: map wins, then data: labels:, then the column name
glm_fig_var <- function(prep, v, opts) fig_label(v, glm_ucfirst(vlab(prep, v)), opts)
glm_fig_level <- function(v, lvl, opts) fig_label(paste0(v, "_", lvl), fig_label(lvl, lvl, opts), opts)
glm_fig_group <- function(prep, g, opts) fig_label(g, glm_ucfirst(prep$spec$labels[[g]] %||% g), opts)

glm_fig_ylab <- function(prep, opts) fig_label("response", response_axis_label(prep), opts)

# Expected response (incl. zero inflation, per offset_per effort) with a ci band
glm_fig_band <- function(prep, fo, nd, draws, ci) {
  est <- predict_mu(prep, fo, nd)
  if (is.null(draws)) return(data.frame(estimate = est, lower = est, upper = est))
  mu <- predict_mu(prep, fo, nd, draws)
  a <- (1 - ci) / 2
  data.frame(estimate = est, lower = apply(mu, 1, stats::quantile, a, na.rm = TRUE),
             upper = apply(mu, 1, stats::quantile, 1 - a, na.rm = TRUE))
}

# Observed response on the plotted scale, or NULL when raw points are not meaningful
# (0/1 data: a rug of the covariate is used instead).
# Counts are not drawn as raw points either: individual counts span a far wider
# range than the expected value, so they flatten the fitted curve (seen on the
# Wood Thrush example: raw counts to 12, effects 0.5-3); counts get a rug.
glm_fig_obs_y <- function(prep) {
  y <- prep$raw[[prep$response]]
  if (prep$fam == "binomial") return(if (is.null(prep$trials)) NULL else y / prep$raw[[prep$trials]])
  if (prep$fam %in% COUNT_FAMILIES) return(NULL)
  if (!is.null(prep$offset)) y <- y / prep$raw[[prep$offset]] * prep$offset_per
  y
}
glm_fig_obs_text <- function(prep) {
  if (prep$fam == "binomial") return("observed proportions (successes / trials)")
  if (!is.null(prep$offset)) return(sprintf("observed values rescaled to %s %s", fmt(prep$offset_per), vlab(prep, prep$offset)))
  "observed values"
}

# x in original units: model-scale z -> fig_x (user units: entry, else recorded scaling)
glm_fig_x <- function(prep, v, raw, opts, scaling) {
  if (!v %in% prep$scaled) return(raw)
  s <- prep$scaling[[v]]
  fig_x((raw - s[["mean"]]) / s[["sd"]], v, opts, scaling)
}

glm_fig_held <- function(prep, except, opts) {
  others <- setdiff(prep$vars, except)
  parts <- character()
  num <- intersect(others, prep$num_vars)
  if (length(num)) parts <- c(parts, sprintf("%s at %s", paste(vapply(num, glm_fig_var, "", prep = prep, opts = opts), collapse = ", "),
                                             if (length(num) > 1) "their means" else "its mean"))
  for (v in intersect(others, prep$fac_vars))
    parts <- c(parts, sprintf("%s = %s", glm_fig_var(prep, v, opts), glm_fig_level(v, prep$xlev[[v]][1], opts)))
  if (length(prep$groups)) parts <- c(parts, sprintf("random effects at zero (a typical %s)", paste(prep$groups, collapse = "/")))
  if (!length(parts)) "" else paste(parts, collapse = "; ")
}

# moderator values for curves: factor levels, or the observed 10th/50th/90th
# percentiles of a numeric partner (mean +/- 1 SD can fall outside the data for
# skewed predictors, e.g. a negative distance)
glm_fig_moderator <- function(prep, part, opts, scaling) {
  if (is.null(part)) return(NULL)
  if (part %in% prep$num_vars) {
    vals <- unname(stats::quantile(prep$raw[[part]], c(0.1, 0.5, 0.9)))
    shown <- glm_fig_x(prep, part, vals, opts, scaling)
    list(values = as.list(vals), labels = sprintf("%s (%s)", fmt(shown), c("10th percentile", "median", "90th percentile")),
         title = glm_fig_var(prep, part, opts), numeric = TRUE)
  } else {
    lv <- prep$xlev[[part]]
    list(values = as.list(lv), labels = vapply(lv, glm_fig_level, "", v = part, opts = opts, USE.NAMES = FALSE),
         title = glm_fig_var(prep, part, opts), numeric = FALSE)
  }
}

glm_fig_effect_num <- function(prep, fo, draws, v, part, opts, ci, scaling) {
  xr <- prep$raw[[v]]
  grid <- seq(min(xr), max(xr), length.out = 100)
  nd0 <- prep$ref_raw[rep(1, length(grid)), , drop = FALSE]
  nd0[[v]] <- grid
  mod <- glm_fig_moderator(prep, part, opts, scaling)
  vals <- if (is.null(mod)) list(NULL) else mod$values
  pred <- do.call(rbind, lapply(seq_along(vals), function(i) {
    nd <- nd0
    if (!is.null(vals[[i]])) nd[[part]] <- if (mod$numeric) vals[[i]] else factor(vals[[i]], levels = prep$xlev[[part]])
    d <- data.frame(x = glm_fig_x(prep, v, grid, opts, scaling), glm_fig_band(prep, fo, to_model_scale(prep, nd), draws, ci))
    if (!is.null(mod)) d$group <- mod$labels[i]
    d
  }))
  if (!is.null(mod)) pred$group <- factor(pred$group, levels = mod$labels)
  xlab <- glm_fig_var(prep, v, opts); ylab <- glm_fig_ylab(prep, opts)
  show <- isTRUE(opts$show_data)
  oy <- glm_fig_obs_y(prep)
  obs <- if (show) data.frame(x = glm_fig_x(prep, v, xr, opts, scaling)) else NULL
  if (show && !is.null(oy)) obs$y <- oy
  p <- fig_cond_effect(pred, xlab, ylab, obs = obs, ylim01 = prep$fam %in% c("binomial", "beta"))
  if (show && !is.null(oy))   # rug as well as points: shows where the data are
    p <- p + ggplot2::geom_rug(data = obs, ggplot2::aes(x = x), inherit.aes = FALSE, alpha = 0.35,
                               sides = "b", length = ggplot2::unit(0.02, "npc"))
  if (!is.null(mod)) {  # second encoding so groups survive greyscale printing
    p <- p + ggplot2::aes(linetype = group) + ggplot2::labs(colour = mod$title, fill = mod$title, linetype = mod$title)
    # a side legend squeezes the panel to ~1/3 of a single column; default to 1.5 columns
    # (override with per_figure: {<id>: {width: single}} plus legend: bottom)
    attr(p, "fig_size") <- list(width = "onehalf", aspect = 0.55)
  }
  if (glm_fig_zero_base(prep)) p <- p + ggplot2::expand_limits(y = 0)
  attr(p, "glm_obs") <- show
  p
}

# Ratio-scale responses (counts, Gamma means): the y axis starts at 0 so the slope
# is not exaggerated by a truncated axis.
glm_fig_zero_base <- function(prep) prep$fam %in% c(COUNT_FAMILIES, "Gamma") || (prep$fam == "gaussian" && prep$link == "log")

glm_fig_effect_fac <- function(prep, fo, draws, v, part, opts, ci) {
  lv <- prep$xlev[[v]]
  lvl_lab <- vapply(lv, glm_fig_level, "", v = v, opts = opts, USE.NAMES = FALSE)
  mod <- if (!is.null(part) && part %in% prep$fac_vars) glm_fig_moderator(prep, part, opts, NULL) else NULL
  vals <- if (is.null(mod)) list(NULL) else mod$values
  pred <- do.call(rbind, lapply(seq_along(vals), function(i) {
    nd <- prep$ref_raw[rep(1, length(lv)), , drop = FALSE]
    nd[[v]] <- factor(lv, levels = lv)
    if (!is.null(vals[[i]])) nd[[part]] <- factor(vals[[i]], levels = prep$xlev[[part]])
    d <- data.frame(level = factor(lvl_lab, levels = lvl_lab), glm_fig_band(prep, fo, to_model_scale(prep, nd), draws, ci))
    if (!is.null(mod)) d$group <- mod$labels[i]
    d
  }))
  xlab <- glm_fig_var(prep, v, opts); ylab <- glm_fig_ylab(prep, opts)
  oy <- glm_fig_obs_y(prep)
  show <- isTRUE(opts$show_data) && !is.null(oy)
  obs <- if (show) data.frame(level = factor(lvl_lab[match(as.character(prep$raw[[v]]), lv)], levels = lvl_lab), y = oy) else NULL
  if (is.null(mod)) {
    p <- fig_cat_effect(pred, xlab, ylab, obs = obs)
  } else {
    pred$group <- factor(pred$group, levels = mod$labels)
    pd <- ggplot2::position_dodge(width = 0.45)
    p <- ggplot2::ggplot(pred, ggplot2::aes(level, estimate, colour = group, shape = group))
    if (show) p <- p + ggplot2::geom_point(data = obs, ggplot2::aes(level, y), inherit.aes = FALSE,
                                           position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
                                           alpha = 0.3, size = 0.9, colour = eco_pal$muted)
    p <- p + ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), position = pd, size = 0.35, linewidth = 0.6) +
      ggplot2::labs(x = xlab, y = ylab, colour = mod$title, shape = mod$title)
    attr(p, "fig_size") <- list(width = "onehalf", aspect = 0.55)
  }
  if (prep$fam %in% c("binomial", "beta")) p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02)))
  if (glm_fig_zero_base(prep)) p <- p + ggplot2::expand_limits(y = 0)
  attr(p, "glm_obs") <- show
  p
}

# Descriptive default label for a fixed-effect coefficient
glm_fig_coef_text <- function(prep, cn, opts) {
  one <- function(p) {
    if (p %in% prep$num_vars) return(list(name = glm_fig_var(prep, p, opts), var = p))
    for (v in prep$fac_vars) if (startsWith(p, v) && substring(p, nchar(v) + 1) %in% prep$xlev[[v]])
      return(list(name = sprintf("%s: %s", glm_fig_var(prep, v, opts), glm_fig_level(v, substring(p, nchar(v) + 1), opts)),
                  ref = sprintf(" vs %s", glm_fig_level(v, prep$xlev[[v]][1], opts))))
    list(name = p)
  }
  parts <- lapply(strsplit(cn, ":", fixed = TRUE)[[1]], one)
  if (length(parts) == 1) {
    q <- parts[[1]]
    if (!is.null(q$var)) return(if (q$var %in% prep$scaled) sprintf("%s, per SD (%s)", q$name, fmt(prep$scaling[[q$var]][["sd"]]))
                                else sprintf("%s, per unit", q$name))
    return(paste0(q$name, q$ref %||% ""))
  }
  paste(vapply(parts, `[[`, "", "name"), collapse = " × ")
}

glm_fig_coefficients <- function(prep, fo, opts, ci) {
  b <- fo$beta; se <- sqrt(pmax(diag(fo$V), 0))
  keep <- names(b) != "(Intercept)"
  if (!any(keep)) return(NULL)
  z <- stats::qnorm(1 - (1 - ci) / 2)
  rp <- ratio_prefix(prep$fam, prep$link)
  tr <- if (is.na(rp)) identity else exp
  cn <- names(b)[keep]
  L <- prep$labels_coef[cn]
  lab <- vapply(seq_along(cn), function(j) {
    d <- glm_fig_coef_text(prep, cn[j], opts)
    d <- fig_label(paste0("beta_", L[[j]]), d, opts)
    # not keyed by the bare column name: that key renames the variable itself (and the
    # row would lose its "per SD" scale)
    if (!is.na(rp)) d <- fig_label(paste0(rp, "_", L[[j]]), d, opts)
    d
  }, "")
  est <- data.frame(param = lab, estimate = tr(b[keep]), lower = tr(b[keep] - z * se[keep]), upper = tr(b[keep] + z * se[keep]))
  words <- if (is.na(rp)) sprintf("Difference in %s", if (prep$link == "identity") "mean" else paste(prep$link, "-link value")) else glm_ucfirst(ratio_words(prep$fam, prep$link))
  xlab <- fig_label("coefficients_axis", sprintf("%s (%d%% CI)", words, round(100 * ci)), opts)
  # wide row labels need more than a single column
  hint_w <- if (max(nchar(lab)) > 30) "onehalf" else NULL
  panel_mm <- max(20, glm_fig_width_mm("coefficients", hint_w, opts) - max(nchar(lab)) * glm_fig_char_mm(opts) - 8)
  p <- fig_coef(est, null_value = if (is.na(rp)) 0 else 1, log_x = FALSE, xlab = glm_fig_wrap(xlab, panel_mm, opts), opts = NULL)
  if (!is.na(rp)) {
    rg <- range(c(est$lower, est$upper, 1), finite = TRUE)
    p <- p + ggplot2::scale_x_log10(breaks = glm_log_breaks(rg[1], rg[2], panel_mm, glm_fig_char_mm(opts)),
                                    labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE))
  }
  if (!is.null(hint_w)) attr(p, "fig_size")$width <- hint_w
  list(plot = p, rp = rp, words = words)
}

glm_fig_residuals <- function(dh, prep, opts) {
  r <- sort(dh$resid); n <- length(r)
  ann <- function(txt) ggplot2::annotate("text", x = 0, y = 1.13, label = txt, hjust = 0, vjust = 1,
                                         size = 0.8 * (opts$base_size %||% 9) / ggplot2::.pt, colour = eco_pal$muted)
  q <- ggplot2::ggplot(data.frame(e = (seq_len(n) - 0.5) / n, o = r), ggplot2::aes(e, o)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = eco_pal$null, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_point(size = 0.7, alpha = 0.6, colour = eco_pal$ink) +
    ann(sprintf("KS test p = %s", fmt(dh$p[["uniformity"]], 2))) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.13), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = fig_label("qq_expected", "Expected quantile (uniform)", opts),
                  y = fig_label("qq_observed", "Scaled residual", opts))
  rk <- rank(dh$pred, ties.method = "average") / length(dh$pred)
  d <- data.frame(x = rk, y = dh$resid)
  nb <- max(3, min(8, floor(n / 40)))   # >= ~40 residuals per bin keeps the quantiles stable
  d$bin <- cut(d$x, breaks = seq(0, 1, length.out = nb + 1), include.lowest = TRUE)
  qs <- do.call(rbind, lapply(split(d, d$bin), function(s) if (nrow(s) >= 5)
    data.frame(x = mean(s$x), q = c("0.25", "0.50", "0.75"), y = stats::quantile(s$y, c(0.25, 0.5, 0.75), names = FALSE))))
  rv <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    ggplot2::geom_point(size = 0.7, alpha = 0.45, colour = eco_pal$ink)
  if (!is.null(qs)) rv <- rv + ggplot2::geom_line(data = qs, ggplot2::aes(x, y, group = q), colour = eco_pal$band, linewidth = 0.8)
  rv <- rv + ann(sprintf("Dispersion test p = %s", fmt(dh$p[["dispersion"]], 2))) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.13), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = fig_label("resid_predicted", "Prediction (rank-transformed)", opts),
                  y = fig_label("qq_observed", "Scaled residual", opts))
  list(q, rv)
}

glm_fig_random_effects <- function(prep, fo, opts, ci) {
  if (is.null(fo$re) || fo$engine != "glmmTMB") return(NULL)
  re <- tryCatch(as.data.frame(glmmTMB::ranef(fo$model, condVar = TRUE)), error = function(e) NULL)
  if (is.null(re) || !nrow(re)) return(NULL)
  re <- re[re$component == "cond", ]
  z <- stats::qnorm(1 - (1 - ci) / 2)
  mult <- prep$link %in% c("log", "logit")
  tr <- if (mult) exp else identity
  re$estimate <- tr(re$condval); re$lower <- tr(re$condval - z * re$condsd); re$upper <- tr(re$condval + z * re$condsd)
  re$grpvar <- as.character(re$grpvar)
  re$facet <- vapply(re$grpvar, glm_fig_group, "", prep = prep, opts = opts)
  re$lab <- mapply(function(g, l) fig_label(paste0(g, "_", l), l, opts), re$grpvar, as.character(re$grp))
  re <- re[order(re$facet, re$estimate), ]
  re$key <- factor(paste(re$facet, re$lab, sep = "\r"), levels = unique(paste(re$facet, re$lab, sep = "\r")))
  nmax <- max(table(re$facet))
  gwords <- paste(unique(re$facet), collapse = " / ")
  # interval level is stated in the caption; a longer title is clipped at single-column width
  xlab <- if (mult) sprintf("%s vs a typical %s", if (prep$link == "logit") "Odds multiplier" else "Multiplier", tolower(gwords))
          else sprintf("Difference from a typical %s", tolower(gwords))
  p <- ggplot2::ggplot(re, ggplot2::aes(estimate, key)) +
    ggplot2::geom_vline(xintercept = if (mult) 1 else 0, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), colour = eco_pal$band, linewidth = 0.45) +
    ggplot2::geom_point(size = 1, colour = eco_pal$ink) +
    ggplot2::scale_y_discrete(labels = if (nmax > 50) NULL else function(k) sub("^.*\r", "", k)) +
    ggplot2::labs(x = fig_label("random_effects_axis", xlab, opts), y = if (length(unique(re$facet)) == 1) re$facet[1] else NULL)
  if (mult) { rg <- range(c(re$lower, re$upper, 1), finite = TRUE)
              pm <- max(20, glm_fig_width_mm("random_effects", NULL, opts) - max(nchar(as.character(re$lab))) * glm_fig_char_mm(opts) - 14)
              p <- p + ggplot2::scale_x_log10(breaks = glm_log_breaks(rg[1], rg[2], pm, glm_fig_char_mm(opts)),
                                              labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE)) }
  if (length(unique(re$facet)) > 1) p <- p + ggplot2::facet_wrap(~facet, scales = "free_y", ncol = 1)
  attr(p, "fig_size") <- list(height = 24 + min(nrow(re), 200) * (if (nmax > 50) 1.2 else 2.6))
  list(plot = p, mult = mult, n = nrow(re))
}

# Core styling for panels composed into one patchwork figure (the core's theme
# reaches only the last panel of a patchwork, so each panel is styled with the
# user's options first).
glm_fig_panel <- function(p, name, opts) apply_fig_opts(p, name, opts)

glm_figures <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2", "patchwork"))
  if (is.null(opts$ci)) opts$ci <- 0.95
  prep <- fit$prep; fo <- fit$fo; st <- fit$settings %||% GLM_DEFAULT_SETTINGS
  ci <- opts$ci; pc <- round(100 * ci)
  scaling <- data$scaling %||% glm_fig_scaling(prep)
  draws <- coef_draws(fo, st$n_draws %||% 2000, st$seed %||% 1)
  band_txt <- if (is.null(draws)) "no band: the covariance of the estimates could not be computed"
              else sprintf("%d%% confidence band from %d draws from the multivariate-normal sampling distribution of the fixed effects",
                           pc, st$n_draws %||% 2000)
  ylab <- glm_fig_ylab(prep, opts)
  figs <- list(); eff <- list()
  for (v in prep$vars) {
    part <- interaction_partners(prep, v)[1]
    if (is.na(part)) part <- NULL
    id <- paste0("effect_", clean_name(v))
    if (v %in% prep$num_vars) {
      p <- glm_fig_effect_num(prep, fo, draws, v, part, opts, ci, scaling)
      cap <- sprintf("%s (line) versus %s, with %s.%s Other predictors held at: %s.%s%s",
                     gsub("\n", " ", glm_fig_axis(opts, id, "y", ylab)), glm_fig_axis(opts, id, "x", glm_fig_var(prep, v, opts)), band_txt,
                     if (!is.null(part)) sprintf(" One curve per %s (%s), because the model includes their interaction; the band is shown for each curve.",
                                                 glm_fig_var(prep, part, opts), if (part %in% prep$num_vars) "at its observed 10th, 50th and 90th percentiles" else "level") else "",
                     glm_fig_held(prep, c(v, part), opts),
                     if (isTRUE(attr(p, "glm_obs"))) {
                       if (!is.null(glm_fig_obs_y(prep))) sprintf(" Grey points: %s, not adjusted for the other predictors or random effects; rug: observed values of %s.",
                                                                  glm_fig_obs_text(prep), glm_fig_axis(opts, id, "x", glm_fig_var(prep, v, opts)))
                       else sprintf(" Rug: observed values of %s.", glm_fig_axis(opts, id, "x", glm_fig_var(prep, v, opts)))
                     } else "",
                     " The curve spans the observed range only.")
    } else {
      p <- glm_fig_effect_fac(prep, fo, draws, v, part, opts, ci)
      fpart <- !is.null(part) && part %in% prep$fac_vars
      cap <- sprintf("%s for each level of %s (points) with %s intervals from %d draws from the sampling distribution of the fixed effects.%s Other predictors held at: %s.%s",
                     gsub("\n", " ", glm_fig_axis(opts, id, "y", ylab)), glm_fig_axis(opts, id, "x", glm_fig_var(prep, v, opts)), sprintf("%d%% confidence", pc), st$n_draws %||% 2000,
                     if (fpart) sprintf(" Shown separately for each level of %s (colour and symbol), because the model includes their interaction.", glm_fig_var(prep, part, opts))
                     else if (!is.null(part)) sprintf(" %s interacts with %s; see effect_%s for the curves per level.", glm_fig_var(prep, v, opts), glm_fig_var(prep, part, opts), clean_name(part)) else "",
                     glm_fig_held(prep, c(v, if (fpart) part), opts),
                     if (isTRUE(attr(p, "glm_obs"))) sprintf(" Grey points: %s, not adjusted for the other predictors or random effects.", glm_fig_obs_text(prep)) else "")
    }
    attr(p, "glm_obs") <- NULL
    p <- glm_fig_grey(p, opts)
    figs[[id]] <- with_caption(p, cap)
    eff[[id]] <- p
  }
  if (length(eff) > 1) {
    pan <- lapply(eff, glm_fig_panel, name = "effects_panel", opts = opts)
    nc <- if (length(pan) <= 4) 2 else 3
    pp <- compose_figures(pan, ncol = nc)
    attr(pp, "fig_size") <- list(width = "double", height = 8 + 62 * ceiling(length(pan) / nc))
    figs$effects_panel <- with_caption(pp, sprintf(
      "Conditional effects of every predictor on the %s (%s). Panels: %s. Lines/points: estimates; bands/bars: %d%% confidence intervals from the sampling distribution of the fixed effects; in each panel the other predictors are held at their mean (continuous) or reference level (factors) and random effects at zero.%s",
      sub("^(.)", "\\L\\1", gsub("\n", " ", ylab), perl = TRUE), if (length(prep$groups)) sprintf("typical %s", paste(tolower(vapply(prep$groups, glm_fig_group, "", prep = prep, opts = opts)), collapse = "/")) else "population level",
      paste(sprintf("%s, %s", LETTERS[seq_along(eff)], vapply(prep$vars, glm_fig_var, "", prep = prep, opts = opts)), collapse = "; "), pc,
      if (!isTRUE(opts$show_data)) "" else if (is.null(glm_fig_obs_y(prep))) " Rug: observed predictor values." else sprintf(" Grey points: %s (not adjusted for other predictors); rug: observed predictor values.", glm_fig_obs_text(prep))))
  } else if (length(eff) == 1) {
    figs$effects_panel <- figs[[1]]
  }
  cf <- glm_fig_coefficients(prep, fo, opts, ci)
  if (!is.null(cf)) figs$coefficients <- with_caption(cf$plot, sprintf(
    "%s for each predictor (points) with %d%% Wald confidence intervals (bars)%s: continuous predictors per 1 SD (SD in brackets, original units), factor levels against the reference level. Dashed line: no effect (%s).%s%s",
    cf$words, pc, if (!is.na(cf$rp)) " on a log axis" else "", if (is.na(cf$rp)) 0 else 1,
    if (!is.na(cf$rp)) sprintf(" A value of 1.5 means a 50%% higher %s.", if (cf$rp == "odds_ratio") "odds" else if (cf$rp == "rate_ratio") "expected count" else "mean") else "",
    if (any(grepl(":", names(fo$beta), fixed = TRUE))) sprintf(" Rows with %s are interaction terms: how much the effect of the first predictor (per SD, or vs reference) %s at the stated level (or per SD) of the second; main-effect rows of interacting predictors are their effects at the mean or reference level of the partner.",
      "\u00d7", if (is.na(cf$rp)) "changes" else "is multiplied") else ""))
  dh <- fit$dharma
  if (is.null(dh)) dh <- tryCatch(glm_dharma(fo, prep, st), error = function(e) NULL)
  if (!is.null(dh) && is.null(dh$error)) {
    pr <- lapply(lapply(glm_fig_residuals(dh, prep, opts), glm_fig_grey, opts = opts), glm_fig_panel, name = "residuals", opts = opts)
    pres <- compose_figures(pr, ncol = 2)
    attr(pres, "fig_size") <- list(width = "onehalf", aspect = 0.5)
    figs$residuals <- with_caption(pres, sprintf(
      "DHARMa simulated quantile residuals (%d simulations). (A) QQ plot against the uniform distribution expected under a correct model; points should follow the dashed line (Kolmogorov-Smirnov uniformity p = %s). (B) Residuals against the rank-transformed model prediction; solid lines are the 0.25, 0.5 and 0.75 quantiles of the residuals in bins of the prediction and should lie flat on the dashed lines (dispersion test p = %s; outlier test p = %s%s).",
      dh$nsim %||% st$n_sim_resid, fmt(dh$p[["uniformity"]], 2), fmt(dh$p[["dispersion"]], 2), fmt(dh$p[["outliers"]], 2),
      if (!is.na(dh$p["zero_inflation"])) sprintf("; zero-inflation p = %s", fmt(dh$p[["zero_inflation"]], 2)) else ""))
  }
  rf <- glm_fig_random_effects(prep, fo, opts, ci)
  if (!is.null(rf)) rf$plot <- glm_fig_grey(rf$plot, opts)
  if (!is.null(rf)) figs$random_effects <- with_caption(rf$plot, sprintf(
    "Estimated random intercepts (conditional modes, points) of %s with %d%% intervals from their conditional SDs (bars), sorted by estimate%s. Dashed line: a typical %s (random effect 0). These intervals do not include uncertainty in the fixed effects.",
    paste(unique(vapply(prep$groups, glm_fig_group, "", prep = prep, opts = opts)), collapse = " and "), pc,
    if (rf$mult) sprintf(", shown as %s on a log axis", if (prep$link == "logit") "multipliers of the odds" else "multipliers of the expected value") else "",
    paste(prep$groups, collapse = "/")))
  # a per_figure aspect: must win over the row-count height hint (the core keeps the hint)
  for (id in names(figs)) {
    u <- opts$per_figure[[id]]
    if (!is.null(u$aspect) && is.null(u$height) && !is.null(attr(figs[[id]], "fig_size")$height))
      attr(figs[[id]], "fig_size")$height <- NULL
  }
  figs
}
