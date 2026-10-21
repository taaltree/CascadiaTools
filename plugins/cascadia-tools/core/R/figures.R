# Publication figures shared by every module.
#
# Modules build their figure suite with the fig_* constructors below and
# return it from figures(fit, data, opts) as a named list of ggplots. The
# core then applies the user's customisation (fig_opts, from the `figures:`
# block of a config) and exports each figure at journal sizes in vector and
# raster formats (save_figure / export_figures).
#
# Design rules (see skills/eco-figures):
#   - categorical colour: Okabe-Ito order, fixed per entity, never cycled;
#     > 8 groups -> facets, not more hues
#   - sequential: one hue light->dark (viridis family) for magnitudes/maps
#   - diverging: two hues with a neutral midpoint, centred on the null value
#   - no dual y-axes; uncertainty always drawn (bands/intervals)
#   - covariates plotted in ORIGINAL units (back-transform scaled ones)
#   - text in ink colours, not series colours; recessive grid

# ---- palettes -------------------------------------------------------------

eco_pal <- list(
  # Okabe & Ito (2008) colour-universal design; black moved last so it can be
  # reserved for "observed data"/reference elements.
  okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7", "#999999", "#000000"),
  # Print-safe greys for grayscale journals
  grey = c("#1a1a1a", "#666666", "#a6a6a6", "#d9d9d9"),
  ink = "#1d2321", muted = "#5d6762", accent = "#2f6f4f",
  band = "#2f6f4f", truth = "#c0392b", null = "#7a7a7a",
  diverging = c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")
)

eco_categorical <- function(n, palette = "okabe_ito") {
  p <- if (length(palette) > 1) palette else eco_pal[[palette]]
  if (n > length(p))
    warning("More groups (", n, ") than palette colours (", length(p),
            "): use facets or group small categories into 'Other'.")
  p[seq_len(min(n, length(p)))]
}

scale_colour_eco <- function(palette = "okabe_ito", ...)
  ggplot2::discrete_scale("colour", palette = function(n) eco_categorical(n, palette), ...)
scale_fill_eco <- function(palette = "okabe_ito", ...)
  ggplot2::discrete_scale("fill", palette = function(n) eco_categorical(n, palette), ...)
scale_fill_eco_seq <- function(option = "viridis", direction = 1, ...)
  ggplot2::scale_fill_viridis_c(option = option, direction = direction, ...)
scale_fill_eco_div <- function(midpoint = 0, ...)
  ggplot2::scale_fill_gradient2(low = eco_pal$diverging[["low"]], mid = eco_pal$diverging[["mid"]],
                                high = eco_pal$diverging[["high"]], midpoint = midpoint, ...)

# ---- theme ----------------------------------------------------------------

# style: "classic" (axis lines, no grid: most journals), "minimal" (light
# grid), "grey" (classic + greyscale-friendly).
theme_eco <- function(base_size = 9, base_family = "", style = "classic", legend = "right") {
  t <- ggplot2::theme_classic(base_size = base_size, base_family = base_family)
  if (style == "minimal")
    t <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(colour = "#e6e6e6", linewidth = 0.3))
  t + ggplot2::theme(
    text = ggplot2::element_text(colour = eco_pal$ink),
    axis.text = ggplot2::element_text(colour = eco_pal$ink, size = ggplot2::rel(0.9)),
    axis.line = if (style == "minimal") ggplot2::element_blank() else ggplot2::element_line(linewidth = 0.35),
    axis.ticks = ggplot2::element_line(linewidth = 0.35),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", hjust = 0),
    legend.position = legend, legend.key.size = ggplot2::unit(0.9, "lines"),
    legend.title = ggplot2::element_text(size = ggplot2::rel(0.95)),
    plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
    plot.tag = ggplot2::element_text(face = "bold"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA))
}

# ---- sizes and export -----------------------------------------------------

# Widths in mm. Approximate; always check the target journal's author guide.
fig_widths <- c(single = 84, onehalf = 129, double = 174, slide = 254)
journal_presets <- list(
  generic = list(single = 84, onehalf = 129, double = 174, max_h = 234, font = 9),
  esa     = list(single = 84, onehalf = 129, double = 174, max_h = 234, font = 9),   # Ecology, Ecosphere, Ecol Appl
  besjournals = list(single = 80, onehalf = 125, double = 166, max_h = 230, font = 9), # J Anim Ecol, MEE, J Appl Ecol
  wiley_jwm = list(single = 85, onehalf = 129, double = 178, max_h = 230, font = 9),
  cjfas   = list(single = 86, onehalf = 129, double = 178, max_h = 230, font = 9),
  plos    = list(single = 132, onehalf = 132, double = 190, max_h = 222, font = 10),
  nature  = list(single = 89, onehalf = 120, double = 183, max_h = 247, font = 7),
  slide   = list(single = 254, onehalf = 254, double = 254, max_h = 143, font = 16)
)

# Default options; every key can be overridden from config `figures:`.
fig_defaults <- function() list(
  journal = "generic", width = "single", height = NULL, aspect = 0.75,
  formats = c("pdf", "png"), dpi = 600, base_size = NULL, font = "",
  style = "classic", palette = "okabe_ito", legend = "right",
  titles = FALSE, panel_tags = TRUE, ci = 0.95, inner_ci = 0.5,
  show_data = TRUE, labels = list(), units = list(), include = NULL,
  exclude = NULL, per_figure = list(), greyscale = FALSE)

fig_opts <- function(cfg_fig = NULL) {
  o <- merge_list(fig_defaults(), cfg_fig)
  j <- journal_presets[[o$journal]]
  if (is.null(j)) { warning("unknown journal preset '", o$journal, "'; using generic"); j <- journal_presets$generic }
  o$.journal <- j
  if (is.null(o$base_size)) o$base_size <- j$font
  if (isTRUE(o$greyscale)) o$palette <- "grey"
  o
}

# Width in mm from "single"/"onehalf"/"double" or a number (mm).
fig_width_mm <- function(width, opts) {
  if (is.numeric(width)) return(width)
  w <- opts$.journal[[width]]
  if (is.null(w)) w <- fig_widths[[width]]
  w
}

# Label lookup: user `labels:` map wins, then the module default.
fig_label <- function(key, default = key, opts = NULL) {
  if (!is.null(opts$labels[[key]])) opts$labels[[key]] else default
}

# Apply global + per-figure customisation to one ggplot.
apply_fig_opts <- function(p, name, opts) {
  pf <- merge_list(opts, opts$per_figure[[name]])
  hint0 <- attr(p, "fig_size")
  p <- p + theme_eco(pf$base_size, pf$font, pf$style, pf$legend)
  # Only add discrete scales if the plot maps colour/fill to a discrete var
  if (fig_maps_discrete(p, "colour")) p <- p + scale_colour_eco(pf$palette)
  if (fig_maps_discrete(p, "fill")) p <- p + scale_fill_eco(pf$palette)
  if (!isTRUE(pf$titles)) p <- p + ggplot2::labs(title = NULL, subtitle = NULL)
  # YAML 1.1 reads an unquoted `y:` key as TRUE; treat it as the y axis
  if (length(pf$axis)) names(pf$axis)[names(pf$axis) %in% c("TRUE", "true")] <- "y"
  for (k in intersect(names(pf$axis), c("x", "y")))
    p <- p + do.call(ggplot2::labs, stats::setNames(list(pf$axis[[k]]), k))
  if (!is.null(pf$xlim)) p <- p + ggplot2::coord_cartesian(xlim = unlist(pf$xlim))
  if (!is.null(pf$ylim)) p <- p + ggplot2::coord_cartesian(ylim = unlist(pf$ylim))
  hint <- hint0
  user <- opts$per_figure[[name]]
  attr(p, "fig_size") <- list(
    width = if (!is.null(user$width)) user$width else if (!is.null(hint$width)) hint$width else pf$width,
    height = if (!is.null(user$height)) user$height else hint$height,
    aspect = if (!is.null(user$aspect)) user$aspect else if (!is.null(hint$aspect)) hint$aspect else pf$aspect)
  p
}

fig_maps_discrete <- function(p, aes_name) {
  m <- c(p$mapping, unlist(lapply(p$layers, function(l) l$mapping), recursive = FALSE))
  m <- m[names(m) %in% c(aes_name, if (aes_name == "colour") "color")]
  if (!length(m)) return(FALSE)
  # Has a scale already been set by the module? Then respect it.
  if (any(vapply(p$scales$scales, function(s) aes_name %in% s$aesthetics, logical(1)))) return(FALSE)
  v <- tryCatch(rlang::eval_tidy(m[[1]], p$data), error = function(e) NULL)
  is.null(v) || is.character(v) || is.factor(v) || is.logical(v)
}

save_figure <- function(p, file_stem, opts = fig_opts(), width = NULL, height = NULL) {
  sz <- attr(p, "fig_size")
  wmm <- fig_width_mm(if (!is.null(width)) width else if (!is.null(sz$width)) sz$width else opts$width, opts)
  hmm <- if (!is.null(height)) height else if (!is.null(sz$height)) sz$height else wmm * (if (!is.null(sz$aspect)) sz$aspect else opts$aspect)
  hmm <- min(hmm, opts$.journal$max_h)
  dir.create(dirname(file_stem), recursive = TRUE, showWarnings = FALSE)
  out <- character()
  for (fmt in opts$formats) {
    f <- paste0(file_stem, ".", fmt)
    dev <- switch(fmt,
      pdf = grDevices::cairo_pdf,
      svg = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg,
      png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png",
      tiff = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff",
      eps = grDevices::cairo_ps, NULL)
    args <- list(filename = f, plot = p, width = wmm, height = hmm, units = "mm", dpi = opts$dpi, bg = "white")
    if (!is.null(dev)) args$device <- dev
    if (fmt == "tiff") args$compression <- "lzw"
    ok <- tryCatch({ do.call(ggplot2::ggsave, args); TRUE },
                   error = function(e) { message("save ", f, " failed: ", conditionMessage(e)); FALSE })
    if (ok) out <- c(out, f)
  }
  invisible(out)
}

# Export a module's whole figure suite + a gallery page + the option sheet.
export_figures <- function(module, fit, data, opts, outdir) {
  figs <- if (is.function(module$figures)) module$figures(fit, data, opts)
          else if (is.function(module$plots)) module$plots(fit, data) else list()
  figs <- Filter(Negate(is.null), figs)
  if (!is.null(opts$include)) figs <- figs[intersect(names(figs), opts$include)]
  if (!is.null(opts$exclude)) figs <- figs[setdiff(names(figs), opts$exclude)]
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (n in names(figs)) {
    p <- figs[[n]]
    cap <- attr(p, "caption")
    if (inherits(p, "ggplot")) p <- apply_fig_opts(p, n, opts)
    files <- save_figure(p, file.path(outdir, n), opts)
    rows[[n]] <- data.frame(figure = n, files = paste(basename(files), collapse = ", "),
                            caption = if (is.null(cap)) "" else cap)
  }
  idx <- do.call(rbind, rows)
  if (!is.null(idx)) {
    utils::write.csv(idx, file.path(outdir, "figures_index.csv"), row.names = FALSE)
    figure_gallery(module, idx, outdir)
  }
  invisible(idx)
}

figure_gallery <- function(module, idx, outdir) {
  items <- vapply(seq_len(nrow(idx)), function(i) {
    png <- file.path(outdir, paste0(idx$figure[i], ".png"))
    img <- if (file.exists(png)) sprintf("<img src='%s' alt='%s'/>", basename(png), idx$figure[i]) else ""
    sprintf("<section><h2>%s</h2>%s<p class='meta'>%s</p><p>%s</p></section>",
            html_escape(idx$figure[i]), img, html_escape(idx$files[i]), html_escape(idx$caption[i]))
  }, "")
  writeLines(html_page(paste(module$title, "- figures"), items,
                       "Edit the `figures:` block of your config and re-run `core/run.R figures` to customise."),
             file.path(outdir, "gallery.html"))
}

# Compose several figures into one multi-panel figure with A/B/C tags.
compose_figures <- function(plots, ncol = NULL, tags = TRUE, guides = "collect") {
  eco_require("patchwork")
  p <- patchwork::wrap_plots(plots, ncol = ncol, guides = guides)
  if (tags) p <- p + patchwork::plot_annotation(tag_levels = "A")
  p
}

with_caption <- function(p, caption) { attr(p, "caption") <- caption; p }

# Back-transform a scaled covariate to original units for x axes:
# data should record centre/scale when modules standardise covariates.
unscale <- function(x, center = 0, scale = 1) x * scale + center

# ---- figure constructors (generic, reused by modules) ---------------------

# Coefficient / effect-size dot-and-whisker ("forest") plot.
# est: data.frame(param, estimate, lower, upper[, group][, lower50, upper50]).
fig_coef <- function(est, null_value = 0, xlab = "Effect (95% interval)", transform = NULL,
                     log_x = FALSE, order = TRUE, opts = NULL) {
  eco_require("ggplot2")
  d <- est
  # User relabelling of parameter names (figures: labels: {beta1: "Forest cover"})
  if (!is.null(opts)) d$param <- vapply(as.character(d$param), fig_label, "", opts = opts)
  if (is.function(transform)) for (k in intersect(c("estimate", "lower", "upper", "lower50", "upper50"), names(d))) d[[k]] <- transform(d[[k]])
  if (order) d$param <- stats::reorder(factor(d$param), d$estimate)
  has_g <- "group" %in% names(d)
  pd <- ggplot2::position_dodge(width = if (has_g) 0.5 else 0)
  base <- if (has_g) ggplot2::aes(estimate, param, colour = group) else ggplot2::aes(estimate, param)
  p <- ggplot2::ggplot(d, base) +
    ggplot2::geom_vline(xintercept = null_value, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    do.call(ggplot2::geom_linerange, c(list(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.5, position = pd),
                                       eco_fixed(has_g, colour = eco_pal$ink)))
  if (all(c("lower50", "upper50") %in% names(d)))
    p <- p + do.call(ggplot2::geom_linerange, c(list(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.3, position = pd),
                                                eco_fixed(has_g, colour = eco_pal$ink)))
  p <- p + do.call(ggplot2::geom_point, c(list(size = 1.8, position = pd), eco_fixed(has_g, colour = eco_pal$ink))) +
    ggplot2::labs(x = xlab, y = NULL, colour = NULL)
  if (log_x) p <- p + ggplot2::scale_x_log10()
  # Height scales with the number of rows (mm) so short tables aren't stretched
  attr(p, "fig_size") <- list(height = 22 + 6 * length(unique(d$param)) * (if (has_g) 1.4 else 1))
  p
}

# Conditional effect: response vs covariate, others fixed; band(s), optional
# observed data (points/rug) and groups. pred: x, estimate, lower, upper
# [, lower50, upper50][, group]. obs: data.frame(x[, y][, group]).
fig_cond_effect <- function(pred, xlab, ylab, obs = NULL, rug = TRUE, ylim01 = FALSE) {
  eco_require("ggplot2")
  has_g <- "group" %in% names(pred)
  m <- if (has_g) ggplot2::aes(x, estimate, colour = group, fill = group) else ggplot2::aes(x, estimate)
  p <- ggplot2::ggplot(pred, m)
  if (!is.null(obs) && "y" %in% names(obs))
    p <- p + ggplot2::geom_point(data = obs, ggplot2::aes(x, y), inherit.aes = FALSE,
                                 alpha = 0.35, size = 0.9, colour = eco_pal$muted)
  p <- p + do.call(ggplot2::geom_ribbon, c(list(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA),
                                           eco_fixed(has_g, fill = eco_pal$band)))
  if (all(c("lower50", "upper50") %in% names(pred)))
    p <- p + do.call(ggplot2::geom_ribbon, c(list(ggplot2::aes(ymin = lower50, ymax = upper50), alpha = 0.25, colour = NA),
                                             eco_fixed(has_g, fill = eco_pal$band)))
  p <- p + do.call(ggplot2::geom_line, c(list(linewidth = 0.8), eco_fixed(has_g, colour = eco_pal$band))) +
    ggplot2::labs(x = xlab, y = ylab, colour = NULL, fill = NULL)
  if (!is.null(obs) && rug && !"y" %in% names(obs))
    p <- p + ggplot2::geom_rug(data = obs, ggplot2::aes(x = x), inherit.aes = FALSE, alpha = 0.35,
                               sides = "b", length = ggplot2::unit(0.02, "npc"))
  if (ylim01) p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02)))
  p
}

# Categorical-predictor effect: point + interval per level (+ raw data jitter).
fig_cat_effect <- function(pred, xlab, ylab, obs = NULL) {
  p <- ggplot2::ggplot(pred, ggplot2::aes(level, estimate))
  if (!is.null(obs)) p <- p + ggplot2::geom_jitter(data = obs, ggplot2::aes(level, y), width = 0.15, height = 0,
                                                   alpha = 0.3, size = 0.9, colour = eco_pal$muted)
  p + ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$ink, size = 0.35, linewidth = 0.6) +
    ggplot2::labs(x = xlab, y = ylab)
}

# Posterior densities (one panel per parameter), optional prior and truth.
# draws: matrix/data.frame with one column per parameter.
fig_posterior <- function(draws, params = colnames(draws), truth = NULL, prior = NULL, ncol = NULL) {
  d <- utils::stack(as.data.frame(draws)[, params, drop = FALSE])
  names(d) <- c("value", "param")
  p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
    ggplot2::geom_density(fill = eco_pal$band, alpha = 0.3, colour = eco_pal$band, linewidth = 0.5) +
    ggplot2::facet_wrap(~param, scales = "free", ncol = ncol) +
    ggplot2::labs(x = "Value", y = "Posterior density")
  if (!is.null(prior)) {
    pd <- utils::stack(as.data.frame(prior)[, intersect(params, colnames(prior)), drop = FALSE]); names(pd) <- c("value", "param")
    p <- p + ggplot2::geom_density(data = pd, linetype = 2, colour = eco_pal$muted, linewidth = 0.4)
  }
  if (!is.null(truth)) p <- p + ggplot2::geom_vline(data = data.frame(param = names(truth), t = truth),
                                                    ggplot2::aes(xintercept = t), colour = eco_pal$truth, linewidth = 0.6)
  p
}

# MCMC traces: samples = list of per-chain matrices.
fig_trace <- function(samples, params = NULL) {
  d <- do.call(rbind, lapply(seq_along(samples), function(k) {
    m <- as.matrix(samples[[k]]); if (!is.null(params)) m <- m[, intersect(params, colnames(m)), drop = FALSE]
    data.frame(iter = seq_len(nrow(m)), chain = factor(k), utils::stack(as.data.frame(m)))
  }))
  ggplot2::ggplot(d, ggplot2::aes(iter, values, colour = chain)) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.8) + ggplot2::facet_wrap(~ind, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = NULL, colour = "Chain")
}

# Gridded surface (density, occupancy, intensity) + optional points.
# grid: x, y, value. points: x, y[, group].
fig_map <- function(grid, fill_lab, points = NULL, diverging = FALSE, midpoint = 0, option = "viridis") {
  p <- ggplot2::ggplot(grid, ggplot2::aes(x, y)) + ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::coord_equal(expand = FALSE) + ggplot2::labs(x = NULL, y = NULL, fill = fill_lab)
  p <- p + if (diverging) scale_fill_eco_div(midpoint) else scale_fill_eco_seq(option)
  if (!is.null(points)) p <- p + ggplot2::geom_point(data = points, ggplot2::aes(x, y), inherit.aes = FALSE,
                                                     shape = 21, fill = "white", colour = "black", size = 1.2, stroke = 0.3)
  p
}

# Time series: estimated latent trajectory with band and observed values.
# est: time, estimate, lower, upper[, group]; obs: time, y[, group]
fig_timeseries <- function(est, ylab, obs = NULL, xlab = "Year") {
  has_g <- "group" %in% names(est)
  m <- if (has_g) ggplot2::aes(time, estimate, colour = group, fill = group) else ggplot2::aes(time, estimate)
  p <- ggplot2::ggplot(est, m) +
    do.call(ggplot2::geom_ribbon, c(list(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA),
                                    eco_fixed(has_g, fill = eco_pal$band))) +
    do.call(ggplot2::geom_line, c(list(linewidth = 0.8), eco_fixed(has_g, colour = eco_pal$band)))
  if (!is.null(obs)) p <- p + ggplot2::geom_point(data = obs, ggplot2::aes(time, y), inherit.aes = FALSE,
                                                  size = 1.2, colour = eco_pal$ink)
  p + ggplot2::labs(x = xlab, y = ylab, colour = NULL, fill = NULL)
}

# Matrix heatmap (correlations, transition probs, pairwise FST). Diverging
# around `midpoint` when values are signed, sequential otherwise.
fig_matrix <- function(mat, fill_lab, diverging = NULL, midpoint = 0, digits = 2, label = TRUE) {
  d <- as.data.frame(as.table(as.matrix(mat))); names(d) <- c("row", "col", "value")
  if (is.null(diverging)) diverging <- any(d$value < 0, na.rm = TRUE)
  p <- ggplot2::ggplot(d, ggplot2::aes(col, row, fill = value)) + ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::coord_equal() + ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(x = NULL, y = NULL, fill = fill_lab) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  p <- p + if (diverging) scale_fill_eco_div(midpoint) else scale_fill_eco_seq()
  if (label) p <- p + ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(value), "", formatC(value, digits = digits, format = "f"))),
                                         size = 2.3, colour = eco_pal$ink)
  p
}

# Posterior / parametric-bootstrap predictive check: observed statistic vs
# replicated distribution.
fig_ppc <- function(rep_stat, obs_stat, stat_lab) {
  pval <- mean(rep_stat >= obs_stat)
  ggplot2::ggplot(data.frame(v = rep_stat), ggplot2::aes(v)) +
    ggplot2::geom_histogram(bins = 30, fill = "#bfcfc7", colour = "white") +
    ggplot2::geom_vline(xintercept = obs_stat, colour = eco_pal$ink, linewidth = 0.8) +
    ggplot2::labs(x = stat_lab, y = "Replicates",
                  caption = sprintf("Observed (line); P(replicate >= observed) = %.2f", pval))
}

# Simulation recovery: estimate vs truth per replicate with coverage colour.
fig_recovery <- function(reps) {
  d <- reps[!is.na(reps$param), ]
  d$covered <- ifelse(d$lower <= d$truth & d$upper >= d$truth, "Interval covers truth", "Misses truth")
  d <- d[order(d$param, d$estimate), ]; d$rank <- stats::ave(d$estimate, d$param, FUN = seq_along)
  ggplot2::ggplot(d, ggplot2::aes(rank, estimate, colour = covered)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), linewidth = 0.3) +
    ggplot2::geom_point(size = 0.6) + ggplot2::geom_hline(ggplot2::aes(yintercept = truth), colour = eco_pal$truth) +
    ggplot2::facet_wrap(~param, scales = "free") +
    ggplot2::scale_colour_manual(values = c("Interval covers truth" = "#7a9a8c", "Misses truth" = "#D55E00")) +
    ggplot2::labs(x = "Replicate (sorted by estimate)", y = "Estimate with 95% interval", colour = NULL) +
    ggplot2::theme(legend.position = "bottom")
}

# Power curve with Monte Carlo error bars (binomial SE).
fig_power <- function(pw, xlab = NULL, target = 0.8) {
  d <- pw[!is.na(pw$power), ]
  d$se <- sqrt(d$power * (1 - d$power) / pmax(d$n_converged, 1))
  ggplot2::ggplot(d, ggplot2::aes(value, power, colour = param)) +
    ggplot2::geom_hline(yintercept = target, linetype = 2, colour = eco_pal$null) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = pmax(0, power - 2 * se), ymax = pmin(1, power + 2 * se)), size = 0.25) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = if (is.null(xlab)) d$vary[1] else xlab, y = "Power", colour = NULL)
}

# Path diagram for DAGs / SEMs from an edge list: from, to, estimate
# [, lower, upper]; node positions optional (x, y) else a layered layout.
fig_path <- function(edges, nodes = NULL, digits = 2, signed = TRUE) {
  if (is.null(nodes)) {
    v <- unique(c(edges$from, edges$to)); lvl <- stats::setNames(rep(0, length(v)), v)
    for (i in seq_along(v)) for (k in seq_len(nrow(edges))) lvl[edges$to[k]] <- max(lvl[edges$to[k]], lvl[edges$from[k]] + 1)
    nodes <- data.frame(name = v, x = lvl[v])
    nodes$y <- stats::ave(nodes$x, nodes$x, FUN = function(z) seq_along(z) - (length(z) + 1) / 2)
  }
  e <- merge(merge(edges, stats::setNames(nodes[, c("name", "x", "y")], c("from", "x0", "y0"))),
             stats::setNames(nodes[, c("name", "x", "y")], c("to", "x1", "y1")))
  has_est <- "estimate" %in% names(e)
  if (has_est) { e$sign <- ifelse(e$estimate >= 0, "positive", "negative"); e$lab <- formatC(e$estimate, digits = digits, format = "f") }
  shrink <- 0.18
  e$xa <- e$x0 + shrink * (e$x1 - e$x0); e$ya <- e$y0 + shrink * (e$y1 - e$y0)
  e$xb <- e$x1 - shrink * (e$x1 - e$x0); e$yb <- e$y1 - shrink * (e$y1 - e$y0)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = e, ggplot2::aes(x = xa, y = ya, xend = xb, yend = yb,
                                                 linewidth = if (has_est) abs(estimate) else 0.6,
                                                 colour = if (has_est && signed) sign else "path"),
                          arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"), type = "closed")) +
    ggplot2::geom_label(data = nodes, ggplot2::aes(x, y, label = name), size = 2.8, label.size = 0.3) +
    ggplot2::scale_linewidth(range = c(0.3, 1.8), guide = "none") +
    ggplot2::scale_colour_manual(values = c(positive = "#0072B2", negative = "#D55E00", path = eco_pal$ink), name = NULL) +
    ggplot2::theme_void()
  if (has_est) p <- p + ggplot2::geom_label(data = e, ggplot2::aes((x0 + x1) / 2, (y0 + y1) / 2, label = lab),
                                            size = 2.4, label.size = 0, fill = "white")
  p
}

# Covariate values for an x axis in ORIGINAL units. Priority: the user's
# `units:` entry in config, then the scaling the module recorded when it
# standardised the covariate (scaling = list(<key> = list(center, scale))).
fig_x <- function(x, key, opts = NULL, scaling = NULL) {
  u <- if (!is.null(opts$units[[key]])) opts$units[[key]] else scaling[[key]]
  if (is.null(u)) return(x)
  x * u$scale + u$center
}
