# nmixture: interpretation, figures and methods text.

nm_fmt <- function(x, d = 2) formatC(x, digits = d, format = "f", big.mark = ",")
nm_ci <- function(r, d = 2) sprintf("%s (95%% %s %s to %s)", nm_fmt(r$estimate, d),
                                    if (isTRUE(r$.bayes)) "CrI" else "CI", nm_fmt(r$lower, d), nm_fmt(r$upper, d))
nm_row <- function(est, n, bayes = FALSE) { r <- est[est$param == n, , drop = FALSE]; if (!nrow(r)) return(NULL); r$.bayes <- bayes; r }

# Variables (raw names) used by a formula, classified numeric / factor in a df
nm_formula_vars <- function(form, df) {
  v <- all.vars(stats::as.formula(form)); v <- intersect(v, names(df))
  list(num = v[vapply(v, function(k) is.numeric(df[[k]]), TRUE)], fac = v[vapply(v, function(k) is.factor(df[[k]]), TRUE)])
}

# Conditional effect of one variable on lambda (expected N per site) or p.
# Other numeric variables at their mean, factors at their reference level.
nm_cond_effect <- function(obj, process, var, n = 100, ci = 0.95, inner = 0.5) {
  des <- obj$des; mix <- obj$settings$mixture
  src <- switch(process, lam = des$sc, p = des$odf, gam = des$ydf, om = des$ydf)
  form <- switch(process, lam = des$formulas$lambda, p = des$formulas$p, gam = des$formulas$gamma, om = des$formulas$omega)
  base <- src[rep(1, n), , drop = FALSE]
  for (k in names(src)) {
    if (is.numeric(src[[k]])) base[[k]] <- if (isTRUE(obj$settings$standardize) && !is.null(des$scaling[[k]])) 0 else mean(src[[k]], na.rm = TRUE)
    else if (is.factor(src[[k]])) base[[k]] <- factor(levels(src[[k]])[1], levels = levels(src[[k]]))
  }
  is_fac <- is.factor(src[[var]])
  if (is_fac) { lv <- levels(src[[var]]); base <- base[rep(1, length(lv)), , drop = FALSE]; base[[var]] <- factor(lv, levels = lv) }
  else { r <- range(src[[var]], na.rm = TRUE); base[[var]] <- seq(r[1], r[2], length.out = n) }
  X <- stats::model.matrix(stats::as.formula(form), stats::model.frame(stats::as.formula(form), base, na.action = stats::na.pass,
                                                                       xlev = lapply(Filter(is.factor, src[intersect(names(src), all.vars(stats::as.formula(form)))]), levels)))
  TH <- obj$draws; if (nrow(TH) > 1000) TH <- TH[seq(1, nrow(TH), length.out = 1000), , drop = FALSE]
  B <- TH[, grep(paste0("^", process, "\\["), colnames(TH)), drop = FALSE]
  eta <- X %*% t(B)
  val <- switch(process,
    lam = exp(eta) * matrix(if (mix == "ZIP") 1 - stats::plogis(TH[, "logit_zi"]) else 1, nrow(eta), ncol(eta), byrow = TRUE),
    p = stats::plogis(eta), gam = exp(eta), om = stats::plogis(eta))
  pt_th <- obj$theta[grep(paste0("^", process, "\\["), names(obj$theta))]
  pt_eta <- drop(X %*% pt_th)
  pt <- switch(process, lam = exp(pt_eta) * (if (mix == "ZIP") 1 - stats::plogis(obj$theta[["logit_zi"]]) else 1),
               p = stats::plogis(pt_eta), gam = exp(pt_eta), om = stats::plogis(pt_eta))
  if (obj$engine == "nimble") pt <- rowMeans(val)
  a <- (1 - ci) / 2; b <- (1 - inner) / 2
  out <- data.frame(estimate = pt, lower = apply(val, 1, stats::quantile, a), upper = apply(val, 1, stats::quantile, 1 - a),
                    lower50 = apply(val, 1, stats::quantile, b), upper50 = apply(val, 1, stats::quantile, 1 - b))
  if (is_fac) out$level <- levels(src[[var]])
  else out$x <- fig_x(base[[var]], var, NULL, des$scaling)
  out
}

nm_proc_label <- c(lam = "Expected abundance per site", p = "Detection probability per visit",
                   gam = "Recruits per site per season", om = "Apparent survival between seasons")
nm_proc_prefix <- c(lam = "abundance", p = "detection", gam = "recruitment", om = "survival")

# Readable default for a raw column name: canopy_cover_pct -> "Canopy cover (%)".
nm_pretty <- function(v) {
  units <- c(pct = "%", m = "m", km = "km", c = "\u00b0C", mm = "mm", ha = "ha", beaufort = "Beaufort")
  w <- strsplit(v, "_")[[1]]; u <- NULL
  if (length(w) > 1 && tolower(w[length(w)]) %in% names(units)) { u <- units[[tolower(w[length(w)])]]; w <- w[-length(w)] }
  if (length(w) > 1 && tolower(w[1]) == "julian" && tolower(w[2]) == "day") w <- c("day", "of", "year")
  out <- paste(w, collapse = " "); out <- paste0(toupper(substr(out, 1, 1)), substr(out, 2, nchar(out)))
  if (!is.null(u)) paste0(out, " (", u, ")") else out
}

# Wrap a long y-axis label onto two lines when it would not fit the panel height
# of figure `nm` (after per_figure width / aspect / height overrides).
nm_fit_ylab <- function(lab, nm, opts) {
  pf <- merge_list(opts, opts$per_figure[[nm]])
  w <- tryCatch(fig_width_mm(pf$width %||% "single", pf), error = function(e) 80)
  h <- if (!is.null(pf$height)) pf$height else w * (pf$aspect %||% 0.75)
  if (grepl("\n", lab)) return(lab)
  # measured text width (mm) at the figure's font size; fallback ~0.5 em per character
  lw <- tryCatch({
    # measure on the device used for export (ragg fonts are wider than the null pdf device's metrics)
    if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(tempfile(fileext = ".png"), 50, 50, units = "mm", res = 72)
    else grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    g <- grid::textGrob(lab, gp = grid::gpar(fontsize = pf$base_size %||% 9, fontfamily = pf$font %||% ""))
    grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE)
  }, error = function(e) nchar(lab) * (pf$base_size %||% 9) * 0.3528 * 0.5)
  if (lw <= h - 12) return(lab)                        # ~12 mm go to the x axis and margins
  paste(strwrap(lab, width = ceiling(nchar(lab) / 2) + 2), collapse = "\n")
}

# Greyscale: the core constructors draw bands and lines in the fixed accent green
# (eco_pal$band), which the core does not re-map; swap it for print greys here.
nm_greyify <- function(p) {
  for (i in seq_along(p$layers)) {
    ap <- p$layers[[i]]$aes_params
    for (k in intersect(names(ap), c("colour", "fill")))
      if (identical(toupper(ap[[k]]), toupper(eco_pal$band))) p$layers[[i]]$aes_params[[k]] <- if (k == "fill") eco_pal$grey[2] else eco_pal$grey[1]
  }
  p
}

# Publication figure suite (docs/FIGURES.md, row "nmixture"). Every label goes
# through fig_label(); covariates are drawn in original units (des$scaling,
# overridable with the figures: units: block); intervals at opts$ci and
# opts$inner_ci; theme, palette and size are applied by the core.
nm_figures <- function(fit, data, opts = fig_opts()) {
  obj <- fit$object; des <- obj$des; est <- fit$estimates
  # labels given in the data config act as defaults; the figures: block wins
  dl <- data$labels %||% obj$labels
  if (length(dl)) for (k in names(dl)) if (is.null(opts$labels[[k]])) opts$labels[[k]] <- dl[[k]]
  ci <- opts$ci %||% 0.95; inner <- opts$inner_ci %||% 0.5
  pc <- function(x) round(100 * x)
  bayes <- obj$engine == "nimble"
  itype <- if (bayes) "credible" else "confidence"
  iabbr <- if (bayes) "CrI" else "CI"
  dyn <- obj$settings$dynamics %||% "constant"
  scaling <- data$scaling %||% des$scaling
  sp <- data$species %||% "individuals"
  # module-drawn accents: greys when the user asks for greyscale (the core only re-maps discrete scales)
  grey <- isTRUE(opts$greyscale) || identical(opts$palette, "grey")
  col_count <- if (grey) eco_pal$grey[2] else eco_pal$okabe_ito[6]   # largest / minimum counts
  col_band <- if (grey) eco_pal$grey[3] else eco_pal$band
  figs <- list(); panels <- list(); plab <- list()
  procs <- c(lam = "lambda", p = "p")
  if (des$model == "dynamic") procs <- c(procs, if (dyn != "notrend") c(gam = "gamma"), om = "omega")
  # ---- conditional effects: abundance_<cov>, detection_<cov> (+ dynamic recruitment_/survival_)
  for (pr in names(procs)) {
    src <- switch(pr, lam = des$sc, p = des$odf, gam = des$ydf, om = des$ydf)
    vv <- nm_formula_vars(des$formulas[[procs[[pr]]]], src)
    for (v in c(vv$num, vv$fac)) {
      ce <- tryCatch(nm_cond_effect(obj, pr, v, ci = ci, inner = inner), error = function(e) NULL)
      if (is.null(ce)) next
      nm <- paste0(nm_proc_prefix[[pr]], "_", nm_clean(v))
      ylab <- fig_label(nm_proc_prefix[[pr]], nm_proc_label[[pr]], opts)
      plab[[nm_proc_prefix[[pr]]]] <- ylab
      ylab_f <- nm_fit_ylab(ylab, nm, opts)
      xlab <- fig_label(v, nm_pretty(v), opts)
      jit <- FALSE
      if ("x" %in% names(ce)) {
        if (!is.null(opts$units[[v]])) ce$x <- fig_x((ce$x - (scaling[[v]]$center %||% 0)) / (scaling[[v]]$scale %||% 1), v, opts, NULL)
        obsx <- fig_x(src[[v]][!is.na(src[[v]])], v, opts, scaling)
        # few distinct values (wind force, visit number): jitter the rug so ticks don't hide it
        # integer-valued covariates (day of year) get the same treatment, kept inside +/- 0.4 units
        intval <- all(abs(obsx - round(obsx)) < 1e-6)
        if (length(unique(obsx)) <= 12 || (intval && length(unique(obsx)) <= 60)) {
          jit <- TRUE
          amp <- 0.03 * diff(range(obsx)); if (intval) amp <- min(amp, 0.4)
          obsx <- obsx + ((seq_along(obsx) %% 9) - 4) / 4 * amp
        }
        obs <- if (isTRUE(opts$show_data)) data.frame(x = obsx) else NULL
        p <- fig_cond_effect(ce, xlab, ylab_f, obs = obs, ylim01 = pr %in% c("p", "om"))
        if (pr %in% c("lam", "gam")) p <- p + ggplot2::expand_limits(y = 0)
      } else {
        ce$level <- vapply(as.character(ce$level), function(l) fig_label(l, l, opts), "")   # level renames
        ce$level <- factor(ce$level, levels = unique(ce$level))
        p <- fig_cat_effect(ce, xlab, ylab_f) +
          ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = eco_pal$ink, linewidth = 1.4) +
          # long (renamed) level labels are angled so they don't collide in narrow composed panels
          (if (max(nchar(levels(ce$level))) > 6) ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 30)) else NULL)
        if (pr %in% c("lam", "gam")) p <- p + ggplot2::expand_limits(y = 0)
        else p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02)))
      }
      if (grey) p <- nm_greyify(p)
      figs[[nm]] <- with_caption(p, sprintf("%s versus %s (%s: estimate; %s: %d%% and %d%% %s intervals; other covariates at their mean or reference level%s).%s",
        ylab, sub(" \\(.*\\)$", "", xlab), if ("x" %in% names(ce)) "line" else "point",
        if ("x" %in% names(ce)) "outer and inner bands" else "thin and thick bars", pc(ci), pc(inner), itype,
        if (isTRUE(opts$show_data) && "x" %in% names(ce)) paste0("; rug: observed values", if (jit) ", jittered slightly" else "") else "",
        if (pr == "p") " This is a detection effect: it changes how many animals are counted, not how many are there." else ""))
      panels[[nm_proc_prefix[[pr]]]][[nm]] <- p
    }
  }
  # multi-panel: all effects on one process (only when there are two or more)
  for (pp in names(panels)) {
    if (length(panels[[pp]]) < 2) next
    nc <- if (length(panels[[pp]]) <= 2) 2 else min(3, length(panels[[pp]]))
    th <- theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic", opts$legend %||% "right")
    pg <- compose_figures(unname(panels[[pp]]), ncol = nc) & th
    attr(pg, "fig_size") <- list(width = "double", aspect = 0.42 * ceiling(length(panels[[pp]]) / nc))
    pid <- paste0(pp, "_effects")
    figs[[pid]] <- with_caption(pg, sprintf("All covariate effects on %s (panels %s; same conventions as the single-effect figures: %d%% and %d%% %s intervals, other covariates at their mean or reference level).",
                                            sub("^(\\w)", "\\L\\1", plab[[pp]], perl = TRUE),
                                            paste(LETTERS[seq_along(panels[[pp]])], collapse = ", "), pc(ci), pc(inner), itype))
  }
  # ---- coefficients: rate ratios and odds ratios, "Process: covariate" at opts$ci (+ inner)
  rr <- est[grepl("^(rate_ratio|odds_ratio)", est$param), c("param", "estimate", "lower", "upper")]
  if (nrow(rr)) {
    if (bayes) {
      pv <- tryCatch(nm_public(obj$draws, des, obj$settings$mixture, dyn)$values, error = function(e) NULL)
      q <- function(n, a) if (!is.null(pv) && n %in% colnames(pv)) stats::quantile(pv[, n], a, names = FALSE) else NA_real_
      rr$lower <- vapply(rr$param, q, 1, a = (1 - ci) / 2); rr$upper <- vapply(rr$param, q, 1, a = 1 - (1 - ci) / 2)
      rr$lower50 <- vapply(rr$param, q, 1, a = (1 - inner) / 2); rr$upper50 <- vapply(rr$param, q, 1, a = 1 - (1 - inner) / 2)
    } else {   # Wald intervals are symmetric on the log scale: rescale the stored 95% interval to opts$ci / inner_ci
      se <- (log(rr$upper) - log(rr$lower)) / (2 * stats::qnorm(0.975))
      zz <- function(l) stats::qnorm(1 - (1 - l) / 2)
      rr$lower <- exp(log(rr$estimate) - zz(ci) * se); rr$upper <- exp(log(rr$estimate) + zz(ci) * se)
      rr$lower50 <- exp(log(rr$estimate) - zz(inner) * se); rr$upper50 <- exp(log(rr$estimate) + zz(inner) * se)
    }
    pre <- c(rate_ratio_recruitment_ = "Recruitment", odds_ratio_survival_ = "Survival",
             odds_ratio_detection_ = "Detection", rate_ratio_ = "Abundance")
    srcs <- list(des$sc, des$odf, des$ydf)
    allv <- unique(unlist(lapply(srcs, names)))
    rr$param <- vapply(rr$param, function(pn) {
      k <- names(pre)[vapply(names(pre), function(x) startsWith(pn, x), TRUE)][1]
      cv <- substring(pn, nchar(k) + 1)
      hit <- allv[vapply(allv, function(v) nm_clean(v) == cv || startsWith(cv, nm_clean(v)), TRUE)]
      nice <- if (length(hit)) {
        v <- hit[which.max(nchar(hit))]; lv <- substring(cv, nchar(nm_clean(v)) + 1)
        lab <- sub(" \\(.*\\)$", "", fig_label(v, nm_pretty(v), opts))
        if (nzchar(lv)) {
          s <- Filter(function(d) !is.null(d) && is.factor(d[[v]]), srcs)
          ref <- if (length(s)) levels(s[[1]][[v]])[1] else NULL
          ll <- fig_label(lv, lv, opts)
          # a renamed level that already names the variable ("Observer 2") is not prefixed again
          paste0(if (startsWith(tolower(ll), tolower(lab))) "" else paste0(lab, " "), ll,
                 if (!is.null(ref)) paste0(" vs ", fig_label(ref, ref, opts)) else "")
        } else lab
      } else cv
      fig_label(pn, paste0(pre[[k]], ": ", nice), opts)
    }, "", USE.NAMES = FALSE)
    rr$param <- factor(rr$param, levels = rev(rr$param))
    figs$coefficients <- with_caption(
      fig_coef(rr, null_value = 1, log_x = TRUE, order = FALSE,
               xlab = fig_label("coefficients", sprintf("Ratio (%d%% %s)", pc(ci), iabbr), opts), opts = NULL),
      sprintf("Rate ratios (abundance, recruitment) and odds ratios (detection, survival) per 1 SD of each numeric covariate, or relative to the reference level of a categorical covariate (points: estimates; thin and thick bars: %d%% and %d%% %s intervals%s); log scale; dashed line = no effect (1).",
              pc(ci), pc(inner), itype, if (bayes) "" else ", Wald intervals on the log scale"))
  }
  # ---- abundance by site (caterpillar)
  sn <- obj$site_N
  if (des$model == "dynamic") sn <- sn[sn$season == des$season_labels[des$T], ]
  sn <- sn[order(sn$estimate), ]; sn$rank <- seq_len(nrow(sn))
  ytop <- max(sn$upper, sn$max_count, na.rm = TRUE)
  pcat <- ggplot2::ggplot(sn, ggplot2::aes(rank, estimate)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = col_band, alpha = 0.55, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(y = max_count), shape = 4, size = 0.9, colour = col_count, na.rm = TRUE) +
    ggplot2::geom_point(size = 0.8, colour = eco_pal$ink) +
    ggplot2::annotate("point", x = 0.03 * nrow(sn), y = c(0.97, 0.89) * ytop, shape = c(16, 4), size = c(0.9, 1.1),
                      colour = c(eco_pal$ink, col_count)) +
    ggplot2::annotate("text", x = 0.06 * nrow(sn), y = c(0.97, 0.89) * ytop, hjust = 0, size = 2.5, colour = eco_pal$ink,
                      label = c(fig_label("site_estimate", "Estimated abundance (95% interval)", opts),
                                fig_label("site_max_count", "Largest count", opts))) +
    ggplot2::scale_x_continuous(breaks = function(l) { b <- pretty(c(1, nrow(sn))); b[b >= 1 & b <= nrow(sn)] }) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::labs(x = fig_label("site_rank", "Site (sorted by estimated abundance)", opts),
                  y = nm_fit_ylab(fig_label("abundance_by_site", "Abundance per site", opts), "abundance_by_site", opts))
  figs$abundance_by_site <- with_caption(pcat, sprintf("Estimated abundance of %s at each surveyed site%s (points; bars: 95%% %s intervals including parameter uncertainty; the per-site interval level is fixed at fit time) and the largest count made there (crosses). The gap between them is what the model attributes to imperfect detection.",
    sp, if (des$model == "dynamic") paste0(" in season ", des$season_labels[des$T]) else "", itype))
  # ---- total abundance distribution
  Nt <- obj$N_draws_total
  tcol <- if (des$model == "static") "total_abundance" else paste0("total_abundance_t", des$T)
  td <- data.frame(v = Nt[, tcol])
  r <- est[est$param == tcol, ]
  qo <- stats::quantile(td$v, c((1 - ci) / 2, 1 - (1 - ci) / 2), names = FALSE)
  qi <- stats::quantile(td$v, c((1 - inner) / 2, 1 - (1 - inner) / 2), names = FALSE)
  sumMax <- sum(sn$max_count, na.rm = TRUE)
  ptot <- ggplot2::ggplot(td, ggplot2::aes(v)) +
    ggplot2::geom_histogram(bins = 40, fill = col_band, alpha = 0.35, colour = "white", linewidth = 0.2) +
    ggplot2::annotate("rect", xmin = qi[1], xmax = qi[2], ymin = -Inf, ymax = Inf, fill = col_band, alpha = 0.15) +
    ggplot2::geom_vline(xintercept = r$estimate, colour = eco_pal$ink, linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = qo, colour = eco_pal$ink, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = sumMax, colour = col_count, linetype = 4, linewidth = 0.6) +
    ggplot2::annotate("text", x = sumMax, y = Inf, label = fig_label("min_count", "Sum of largest counts", opts),
                      angle = 90, hjust = 1.05, vjust = -0.5, size = 2.4, colour = eco_pal$ink) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = fig_label("total_abundance", sprintf("Total abundance over the %d surveyed sites", des$M), opts),
                  y = fig_label("draws", if (bayes) "Posterior draws" else "Draws", opts))
  figs$total_abundance <- with_caption(ptot,
    sprintf("Distribution of total abundance of %s across the surveyed sites%s (%s). Solid line: estimate; dashed lines: %d%% interval; shaded: %d%% interval. Dot-dash line: sum of the largest count at each site (%d), a minimum count.",
            sp, if (des$model == "dynamic") paste0(" in season ", des$season_labels[des$T]) else "",
            if (bayes) "posterior draws" else "draws from the sampling distribution of the parameters combined with the conditional distribution of each site's N",
            pc(ci), pc(inner), sumMax))
  # ---- trajectory (dynamic)
  if (des$model == "dynamic") {
    cols <- paste0("total_abundance_t", seq_len(des$T))
    tr <- data.frame(time = seq_len(des$T), estimate = est$estimate[match(cols, est$param)],
                     lower = apply(Nt[, cols, drop = FALSE], 2, stats::quantile, (1 - ci) / 2),
                     upper = apply(Nt[, cols, drop = FALSE], 2, stats::quantile, 1 - (1 - ci) / 2),
                     lower50 = apply(Nt[, cols, drop = FALSE], 2, stats::quantile, (1 - inner) / 2),
                     upper50 = apply(Nt[, cols, drop = FALSE], 2, stats::quantile, 1 - (1 - inner) / 2))
    obsd <- data.frame(time = seq_len(des$T), y = vapply(seq_len(des$T), function(t) sum(apply(matrix(des$y[, , t], des$M), 1, function(v) if (all(is.na(v))) 0 else max(v, na.rm = TRUE))), 1))
    p <- fig_timeseries(tr, fig_label("trajectory", "Total abundance (all sites)", opts), obs = obsd, xlab = fig_label("season", "Season", opts)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = col_band, alpha = 0.25) +
      ggplot2::scale_x_continuous(breaks = seq_len(des$T), labels = des$season_labels) +
      ggplot2::expand_limits(y = 0)
    if (grey) p <- nm_greyify(p)
    figs$trajectory <- with_caption(p, sprintf("Estimated total abundance of %s per season over the %d surveyed sites (line; bands: %d%% and %d%% %s intervals); points: sum over sites of the largest count (a minimum count index).",
                                               sp, des$M, pc(ci), pc(inner), itype))
  }
  # ---- GOF
  if (!is.null(obj$gof)) {
    g <- obj$gof
    figs$gof <- with_caption(fig_ppc(g$chisq_rep, g$chisq_obs, fig_label("gof", "Chi-square discrepancy", opts)) +
                               ggplot2::labs(y = fig_label("replicates", "Simulated data sets", opts)),
                             sprintf("Goodness of fit, %s: histogram of the chi-square statistic for data simulated from the fitted model; vertical line: observed data (c-hat %.2f, p = %.2f). Small p (< 0.05) means the counts vary more than the model allows.",
                                     g$method, g$c_hat, g$p_chisq))
  }
  # ---- K sensitivity
  if (!is.null(obj$k_sens)) {
    ks <- obj$k_sens
    figs$K_sensitivity <- with_caption(
      ggplot2::ggplot(ks, ggplot2::aes(K, total_abundance)) + ggplot2::geom_line(colour = col_band) +
        ggplot2::geom_point(colour = eco_pal$ink, size = 1.5) +
        ggplot2::scale_y_continuous(limits = range(c(ks$total_abundance * 0.95, ks$total_abundance * 1.05))) +
        ggplot2::labs(x = fig_label("K", "Upper bound K on abundance per site", opts), y = fig_label("K_sensitivity", "Total abundance (plug-in)", opts)),
      "Sensitivity of estimated total abundance to the upper bound K used to sum over possible abundances (Dennis et al. 2015); the model is refitted at each K and total abundance is the sum of each site's expected N given its counts at the point estimates (a plug-in value, so it can differ slightly from the reported total, which also averages over parameter uncertainty). The y axis spans +/- 5% of the estimate; a flat line means K is large enough.")
  }
  figs
}

nm_plots <- function(fit, data) nm_figures(fit, data, fig_opts())

nm_interpret <- function(fit, data) {
  obj <- fit$object; des <- obj$des; est <- fit$estimates; st <- obj$settings
  bayes <- obj$engine == "nimble"
  R <- function(n) nm_row(est, n, bayes)
  sp <- data$species %||% "individuals"
  labs <- data$labels %||% obj$labels
  lab <- function(v) { l <- labs[[v]]; if (is.null(l)) v else sub(" \\(.*\\)$", "", l) }   # readable name, unit suffix dropped
  out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge or some parameters are not identifiable (see Diagnostics). Do not interpret the numbers below until this is resolved.")
  if (des$model == "static") {
    a <- R("abundance_at_average_site"); tt <- R("total_abundance"); mm <- R("mean_abundance_per_site")
    smax <- sum(apply(matrix(des$y, des$M), 1, function(v) if (all(is.na(v))) 0 else max(v, na.rm = TRUE)))
    lf <- nm_formula_vars(des$formulas$lambda, des$sc)$fac
    lref <- if (length(lf)) sprintf(" (%s)", paste(vapply(lf, function(v) sprintf("%s '%s'", lab(v), levels(des$sc[[v]])[1]), ""), collapse = ", ")) else ""
    out <- c(out, sprintf("At a site with average covariate values%s the model expects %s %s. Summed over the %d surveyed sites, estimated abundance is %s (an average of %s per site). The sum of the largest count at each site was %d, so the model infers that about %.0f%% of the animals present were missed even on the visit with the highest count.",
                          lref, nm_ci(a), sp, des$M, nm_ci(tt, 0), nm_ci(mm), smax, 100 * (1 - smax / tt$estimate)))
  } else {
    a <- R("initial_abundance_per_site"); g <- R("recruitment_per_site"); om <- R("apparent_survival")
    out <- c(out, sprintf("In the first season (%s) the model expects %s %s at a site with average covariates.", des$season_labels[1], nm_ci(a), sp))
    out <- c(out, sprintf("Between seasons, %s of animals survive and stay at a site (apparent survival: deaths and permanent emigration are not separated)%s.",
                          paste0(nm_fmt(100 * om$estimate, 0), "% (95% ", if (bayes) "CrI" else "CI", " ", nm_fmt(100 * om$lower, 0), "-", nm_fmt(100 * om$upper, 0), "%)"),
                          if (!is.null(g)) sprintf(", and each site gains %s new animals per season (births plus immigration)", nm_ci(g)) else ""))
    t1 <- R("total_abundance_t1"); tT <- R(paste0("total_abundance_t", des$T)); rr <- R("abundance_ratio_last_first")
    out <- c(out, sprintf("Total abundance over the %d surveyed sites was %s in %s and %s in %s, a ratio of %s (%s).",
                          des$M, nm_ci(t1, 0), des$season_labels[1], nm_ci(tT, 0), des$season_labels[des$T], nm_ci(rr),
                          if (rr$lower > 1) "an increase" else if (rr$upper < 1) "a decline" else "the data cannot distinguish an increase from a decline"))
  }
  # covariate effects on abundance (and dynamics)
  eff <- function(pre, rrpre, proc, what, src) {
    o <- character()
    vv <- nm_formula_vars(des$formulas[[proc]], src)
    for (v in vv$num) {
      r <- R(paste0(rrpre, "_", nm_clean(v))); if (is.null(r)) next
      s <- des$scaling[[v]]
      per <- if (!is.null(s)) sprintf("1 SD (%s units)", nm_fmt(s$scale, 2)) else "1 unit"
      pct <- function(x) sprintf("%+.0f%%", 100 * (x - 1))
      txt <- sprintf("Each increase of %s in %s is associated with %s being multiplied by %s, i.e. %s (95%% %s %s to %s).",
                     per, lab(v), what, nm_ci(r), pct(r$estimate), if (bayes) "CrI" else "CI", pct(r$lower), pct(r$upper))
      if (!is.null(s)) {
        nice <- signif(s$scale, 1)
        b <- R(paste0(pre, "_beta_", nm_clean(v)))
        f <- exp(c(b$estimate, b$lower, b$upper) * nice / s$scale)
        txt <- paste(txt, sprintf("Per %s units of %s the factor is %s (%s to %s).", format(nice), lab(v), nm_fmt(f[1]), nm_fmt(f[2]), nm_fmt(f[3])))
      }
      ce <- tryCatch(nm_cond_effect(obj, c(lambda = "lam", gamma = "gam")[[proc]], v), error = function(e) NULL)
      if (!is.null(ce) && proc == "lambda")
        txt <- paste(txt, sprintf("Across the observed range (%s to %s) expected abundance per site goes from %s to %s, other covariates held at their mean.",
                                  nm_fmt(min(ce$x)), nm_fmt(max(ce$x)), nm_fmt(ce$estimate[1]), nm_fmt(ce$estimate[nrow(ce)])))
      if ((r$lower - 1) * (r$upper - 1) < 0) txt <- paste(txt, "The interval includes no effect (1), so the direction is uncertain.")
      o <- c(o, txt)
    }
    for (v in vv$fac) {
      rows <- est[grepl(paste0("^", rrpre, "_", nm_clean(v)), est$param), ]
      for (k in seq_len(nrow(rows))) {
        lv <- sub(paste0("^", rrpre, "_", nm_clean(v)), "", rows$param[k])
        o <- c(o, sprintf("Sites with %s = '%s' have %s multiplied by %s relative to the reference level '%s' (%s).",
                          lab(v), lv, what, nm_ci(cbind(rows[k, ], .bayes = bayes)), levels(src[[v]])[1],
                          if ((rows$lower[k] - 1) * (rows$upper[k] - 1) < 0) "the interval includes no difference (1)" else sprintf("%+.0f%%", 100 * (rows$estimate[k] - 1))))
      }
    }
    o
  }
  out <- c(out, eff("lambda", "rate_ratio", "lambda", "expected abundance", des$sc))
  if (des$model == "dynamic") {
    if (st$dynamics != "notrend") out <- c(out, eff("gamma", "rate_ratio_recruitment", "gamma", "recruitment", des$ydf))
    vo <- nm_formula_vars(des$formulas$omega, des$ydf)
    for (v in vo$num) {
      r <- R(paste0("odds_ratio_survival_", nm_clean(v))); if (is.null(r)) next
      out <- c(out, sprintf("Each increase of %s in %s (in the season before the transition) multiplies the odds that an animal survives and stays at the site by %s%s. The survival figure above refers to average %s.",
                            if (!is.null(des$scaling[[v]])) sprintf("1 SD (%s units)", nm_fmt(des$scaling[[v]]$scale, 2)) else "1 unit", lab(v), nm_ci(r),
                            if ((r$lower - 1) * (r$upper - 1) < 0) "; the interval includes no effect (1)" else "", lab(v)))
    }
    for (v in vo$fac) {
      rows <- est[grepl(paste0("^odds_ratio_survival_", nm_clean(v)), est$param), ]
      for (k in seq_len(nrow(rows))) out <- c(out, sprintf("Survival odds for %s = '%s' relative to '%s': %s.", lab(v),
                                                           sub(paste0("^odds_ratio_survival_", nm_clean(v)), "", rows$param[k]), levels(des$ydf[[v]])[1], nm_ci(cbind(rows[k, ], .bayes = bayes))))
    }
  }
  # detection
  p <- R("detection_prob")
  cum <- 1 - (1 - c(p$estimate, p$lower, p$upper))^des$J
  vv <- nm_formula_vars(des$formulas$p, des$odf)
  # detection_prob is at the mean of numeric covariates and the reference level of factors: say which
  refl <- if (length(vv$fac)) sprintf(" (%s)", paste(vapply(vv$fac, function(v) sprintf("%s '%s'", lab(v), levels(des$odf[[v]])[1]), ""), collapse = ", ")) else ""
  out <- c(out, sprintf("Detection: on a visit with average covariate values%s each animal present had a %s probability of being counted; over %d such visits the probability of counting an animal at least once is %s (%s to %s). Detection effects describe the counting process, not where animals are.",
                        refl, nm_ci(p), des$J, nm_fmt(cum[1]), nm_fmt(cum[2]), nm_fmt(cum[3])))
  for (v in vv$num) {
    r <- R(paste0("odds_ratio_detection_", nm_clean(v))); if (is.null(r)) next
    out <- c(out, sprintf("Each increase of %s in %s multiplies the odds of detecting an animal by %s.",
                          if (!is.null(des$scaling[[v]])) sprintf("1 SD (%s units)", nm_fmt(des$scaling[[v]]$scale, 2)) else "1 unit", lab(v), nm_ci(r)))
  }
  for (v in vv$fac) {
    rows <- est[grepl(paste0("^odds_ratio_detection_", nm_clean(v)), est$param), ]
    for (k in seq_len(nrow(rows))) out <- c(out, sprintf("Detection odds for %s '%s' relative to '%s': %s.", lab(v), sub(paste0("^odds_ratio_detection_", nm_clean(v)), "", rows$param[k]), levels(des$odf[[v]])[1], nm_ci(cbind(rows[k, ], .bayes = bayes))))
  }
  if (st$mixture == "NB") { r <- R("nb_size")
    out <- c(out, sprintf("Abundance varies among sites more than a Poisson distribution allows (negative binomial size %s; variance = mean + mean^2/size). Negative-binomial N-mixture models are the most prone to unstable, inflated abundance (Kery 2018).", nm_ci(r))) }
  if (st$mixture == "ZIP") { r <- R("zero_inflation_prob")
    out <- c(out, sprintf("An estimated %s of sites are structurally unoccupied (zero inflation probability %s).", paste0(nm_fmt(100 * r$estimate, 0), "%"), nm_ci(r))) }
  if (!is.null(obj$gof)) {
    g <- obj$gof
    out <- c(out, if (min(g$p_chisq, g$p_ft) < 0.05)
      sprintf("Model fit: the goodness-of-fit test indicates lack of fit (c-hat %.2f, p = %.3f). Counts vary more than the model allows; absolute abundance is likely biased. Try mixture: NB, add detection covariates, or treat abundance as a relative index.", g$c_hat, min(g$p_chisq, g$p_ft))
      else sprintf("Model fit: %s (c-hat %.2f, smaller of the chi-square and Freeman-Tukey p-values = %.2f). In our simulations this test caught strong site-level detection heterogeneity and overdispersion, but it can miss double counting, movement between visits and individual-level detection heterogeneity (Duarte et al. 2018; Link et al. 2018), so it does not prove the assumptions hold.",
                   if (min(g$p_chisq, g$p_ft) < 0.10) "no significant lack of fit at the 5% level, but the p-value is borderline" else "no evidence of lack of fit",
                   g$c_hat, min(g$p_chisq, g$p_ft)))
  }
  if (!is.null(obj$mixtures)) {
    mx <- obj$mixtures; best <- mx$mixture[which.min(mx$AIC)]
    out <- c(out, sprintf("Choice of abundance distribution (AIC): %s. Expected abundance at an average site under each: %s.",
                          paste(sprintf("%s dAIC %.1f", mx$mixture, mx$delta_AIC), collapse = "; "),
                          paste(sprintf("%s %s", mx$mixture, nm_fmt(mx$abundance_at_average_site)), collapse = ", ")))
  }
  out <- c(out, "Caution on absolute abundance: N-mixture abundance relies on closure within a season, no double counting, independent detection of individuals and correctly modelled detection. Small violations can bias absolute N substantially without failing goodness-of-fit tests (Barker et al. 2018; Link et al. 2018). Relative differences among sites and covariate effects (rate ratios) are more robust than the absolute totals. All covariate relationships are associations, not demonstrated causes.")
  out
}

nm_methods_text <- function(fit, data, settings) {
  obj <- fit$object; des <- obj$des; st <- obj$settings
  mixname <- c(P = "Poisson", NB = "negative binomial", ZIP = "zero-inflated Poisson")[[st$mixture]]
  f <- des$formulas
  eng <- switch(obj$engine,
    unmarked = sprintf("by maximum likelihood with the R package unmarked (%s; Fiske & Chandler 2011), summing over abundance up to K = %d per site", if (des$model == "static") "pcount" else "pcountOpen", obj$K),
    rlik = sprintf("by maximum likelihood with a hand-written likelihood (K = %d)", obj$K),
    nimble = sprintf("in a Bayesian framework with NIMBLE (de Valpine et al. 2017) using the marginalised N-mixture distributions of nimbleEcology (Goldstein & de Valpine 2022), %d chains of %d iterations (%d burn-in, thinning %d) with automated factor slice sampling. Priors: log expected abundance intercept ~ Normal(0, %g), logit detection intercept ~ Normal(0, %g), covariate coefficients ~ Normal(0, %g)%s",
                     st$nchains, st$niter, st$nburnin, st$thin, st$prior_sd_lambda_int, st$prior_sd_p_int, st$prior_sd_beta,
                     switch(st$mixture, NB = sprintf(", log negative-binomial size ~ Normal(0, %g)", st$prior_sd_log_size), ZIP = ", zero-inflation probability ~ Uniform(0, 1)", "")))
  base <- if (des$model == "static")
    sprintf("Repeated counts at %d sites (up to %d visits each) were analysed with a binomial N-mixture model (Royle 2004). Site abundance N_i followed a %s distribution with log(expected abundance) %s, and each count was Binomial(N_i, p_ij) with logit(p) %s. The population at each site was assumed closed across visits.",
            des$M, des$J, mixname, f$lambda, f$p)
  else sprintf("Repeated counts at %d sites over %d seasons (up to %d visits per season) were analysed with the open-population N-mixture model of Dail & Madsen (2011). Initial abundance followed a %s distribution with log(expected abundance) %s; between seasons, survivors were Binomial(N_t-1, omega) with logit(omega) %s and recruits Poisson(gamma)%s; counts were Binomial(N_t, p) with logit(p) %s. Populations were assumed closed among visits within a season.",
               des$M, des$T, des$J, mixname, f$lambda, f$omega,
               if (st$dynamics == "notrend") " with gamma = lambda(1 - omega) (no-trend dynamics)" else sprintf(" with log(gamma) %s", f$gamma), f$p)
  std <- if (isTRUE(st$standardize) && length(des$scaling)) " Numeric covariates were standardised (mean 0, SD 1), so effects are per 1 SD; rate ratios exp(beta) are reported for abundance and odds ratios for detection." else ""
  der <- if (obj$engine == "nimble") " Site abundances and totals were drawn from their conditional distribution given the counts at each posterior draw."
         else sprintf(" Site and total abundances were obtained by drawing %d parameter vectors from the asymptotic normal distribution of the estimates and, for each, drawing every site's abundance from its conditional distribution given the counts, so intervals include parameter uncertainty.", st$n_draws)
  gof <- if (!is.null(obj$gof)) sprintf(" Goodness of fit was assessed with a %s of the chi-square and Freeman-Tukey statistics (%d replicates).", obj$gof$method, obj$gof$nboot) else ""
  ks <- if (!is.null(obj$k_sens)) sprintf(" Sensitivity to the truncation bound K was checked by refitting with K = %s (Dennis et al. 2015).", paste(obj$k_sens$K, collapse = ", ")) else ""
  paste0(base, " The model was fitted ", eng, ".", std, der, gof, ks)
}
