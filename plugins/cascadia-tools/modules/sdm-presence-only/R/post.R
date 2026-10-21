# Interpretation, figures and methods text for sdm-presence-only.

sdm_lab <- function(v, data, opts = NULL) {
  d <- if (!is.null(data$labels[[v]])) data$labels[[v]] else v
  fig_label(v, d, opts)
}

sdm_th <- function(opts) theme_eco(if (is.null(opts$base_size)) 9 else opts$base_size, opts$font, opts$style, opts$legend)

# unit from a label such as "Elevation (m)" -> "m"
sdm_unit <- function(lab) { m <- regmatches(lab, regexpr("\\(([^()]*)\\)\\s*$", lab)); if (length(m)) gsub("^\\(|\\)\\s*$", "", m) else "units" }
# step = TRUE: an increment (for % covariates that is percentage points, not a relative %)
sdm_amount <- function(x, lab, step = FALSE) { u <- sdm_unit(lab); if (u == "%") { if (step) paste(x, "percentage points") else paste0(x, "%") } else paste(x, u) }

sdm_nice_unit <- function(s) { u <- 10^floor(log10(s)); if (s / u >= 5) 5 * u else u }

sdm_fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")

sdm_ratio_text <- function(r, lo, hi) {
  pct <- function(z) sprintf("%+.0f%%", 100 * (z - 1))
  sprintf("x%s (95%% CI %s to %s; %s)", sdm_fmt(r), sdm_fmt(lo), sdm_fmt(hi), pct(r))
}

sdm_interpret <- function(fit, data) {
  o <- fit$object; prep <- o$prep; est <- fit$estimates
  integ <- "log_density" %in% names(o$eng$coef)
  cu <- if (identical(data$coord_unit, "km")) "km" else "m"
  area <- sum(prep$A)
  out <- sprintf("The analysis used %d presence records inside a study region of %s km2 (%d grid cells)%s, fitted with the %s engine.",
                 prep$n_pres, format(round(area), big.mark = ","), nrow(prep$X),
                 if (!is.null(prep$survey)) sprintf(" and %d survey sites (%d with the species detected)", length(prep$survey$cell), sum(prep$survey$ydet > 0)) else "",
                 o$eng$engine)
  resp <- if (integ) "expected density (individuals per km2)" else "relative intensity of occurrence (how densely the species' locations are packed in the landscape)"
  for (v in prep$habitat) {
    sc <- prep$scaling[[v]]; lab <- sdm_lab(v, data)
    r <- est[est$param == paste0("intensity_beta_", v), ]
    rq <- est[est$param == paste0("intensity_beta_", v, "_sq"), ]
    if (!nrow(r)) next
    if (!nrow(rq)) {
      u <- sdm_nice_unit(sc$scale)
      txt <- sprintf("%s: holding the other covariates at their landscape means, each 1 SD increase (%s) is associated with a change in %s of %s; per %s, x%s (95%% CI %s to %s).",
                     lab, sdm_amount(signif(sc$scale, 3), lab, TRUE), resp, sdm_ratio_text(exp(r$estimate), exp(r$lower), exp(r$upper)),
                     sdm_amount(format(u), lab, TRUE), sdm_fmt(exp(r$estimate * u / sc$scale)), sdm_fmt(exp(r$lower * u / sc$scale)), sdm_fmt(exp(r$upper * u / sc$scale)))
      rg <- range(prep$grid[[v]])
      zr <- (rg - sc$center) / sc$scale
      txt <- paste(txt, sprintf("Across the observed range (%s to %s) this is a %s-fold difference.", sdm_amount(signif(rg[1], 3), lab), sdm_amount(signif(rg[2], 3), lab),
                                sdm_fmt(exp(abs(r$estimate * diff(zr))), 1)))
      if (r$lower < 0 && r$upper > 0) txt <- paste(txt, "The interval includes no change, so the data are consistent with no association.")
    } else {
      b1 <- r$estimate; b2 <- rq$estimate
      zopt <- -b1 / (2 * b2); xopt <- sc$center + zopt * sc$scale
      rg <- range(prep$grid[[v]])
      shape <- if (b2 < 0) "peaks (hump-shaped response)" else "is lowest (U-shaped response)"
      txt <- sprintf("%s has a curved association (quadratic term %s, 95%% CI %s to %s on the log scale): %s %s at about %s%s. See the conditional-effect plot for the full curve with its band.",
                     lab, sdm_fmt(b2), sdm_fmt(rq$lower), sdm_fmt(rq$upper), resp, shape, signif(xopt, 3),
                     if (xopt < rg[1] || xopt > rg[2]) " (outside the observed range, so within the data the response is monotone)" else "")
    }
    out <- c(out, txt)
  }
  for (b in prep$bias) {
    r <- est[est$param == paste0("bias_beta_", b), ]; sc <- prep$scaling[[b]]; u <- sdm_nice_unit(sc$scale)
    blab <- sdm_lab(b, data)
    txt <- sprintf("Sampling bias, %s: the rate at which records were made changed by %s per 1 SD (%s) and x%s (95%% CI %s to %s) per %s. This describes where people looked, not where the species lives. The habitat effects above and the map are corrected for it by holding %s at a common value for prediction (Warton et al. 2013).",
                   blab, sdm_ratio_text(exp(r$estimate), exp(r$lower), exp(r$upper)), sdm_amount(signif(sc$scale, 3), blab, TRUE),
                   sdm_fmt(exp(r$estimate * u / sc$scale)), sdm_fmt(exp(r$lower * u / sc$scale)), sdm_fmt(exp(r$upper * u / sc$scale)),
                   sdm_amount(format(u), blab, TRUE), blab)
    if (!is.null(o$naive)) {
      hb <- paste0("intensity_beta_", prep$habitat)
      cm <- stats::cor(prep$grid[, c(b, prep$habitat)])[1, -1]
      sh <- vapply(hb, function(k) sprintf("%s %s instead of %s", sdm_lab(sub("intensity_beta_", "", k), data),
                                           sdm_fmt(o$naive[[k]]), sdm_fmt(o$eng$coef[[k]])), "")
      txt <- paste(txt, sprintf("For comparison, ignoring the bias covariate would give log-intensity slopes per SD of %s (correlation of %s with the habitat covariates: %s). Differences show how much the sampling pattern would otherwise be mistaken for habitat preference.",
                                paste(sh, collapse = "; "), sdm_lab(b, data), paste(sprintf("%.2f", cm), collapse = ", ")))
    }
    out <- c(out, txt)
  }
  if (!length(prep$bias) && !isFALSE(o$settings$bias_correction))
    out <- c(out, "No sampling-bias covariate was given. Presence-only records almost always follow observer effort (roads, towns, trails, protected areas). If effort is correlated with any habitat covariate, the habitat effects above are biased and cannot be checked from these data. Supply a bias covariate (e.g. distance to road) or survey data.")
  if (integ) {
    g <- function(p) est[est$param == p, ]
    d <- g("density_at_mean_covariates"); N <- g("expected_total_abundance"); oc <- g("mean_occupancy_prob"); p <- g("detection_prob")
    txt <- sprintf("Because survey data with non-detections were included, the absolute level is estimable: density at average habitat is %s individuals per km2 (95%% CI %s to %s), giving an expected %s individuals (95%% CI %s to %s) across the %s km2 region, and an average probability of %s (95%% CI %s to %s) that a %s km2 survey plot is occupied.",
                   signif(d$estimate, 3), signif(d$lower, 3), signif(d$upper, 3), format(round(N$estimate), big.mark = ","),
                   format(round(N$lower), big.mark = ","), format(round(N$upper), big.mark = ","), format(round(area), big.mark = ","),
                   sdm_fmt(oc$estimate), sdm_fmt(oc$lower), sdm_fmt(oc$upper), signif(o$survey_area_km2, 3))
    if (nrow(p)) txt <- paste(txt, sprintf("Detection probability per survey visit at an occupied site was %s (95%% CI %s to %s).", sdm_fmt(p$estimate), sdm_fmt(p$lower), sdm_fmt(p$upper)))
    txt <- paste(txt, "These absolute numbers assume that an occupied plot holds at least one individual of a Poisson process, that the survey sites and presence records respond to habitat in the same way, and that detection does not vary among sites; abundance is the most assumption-dependent of these and is better read as an index than a census.")
    out <- c(out, txt)
  } else {
    out <- c(out, "What presence-only data cannot tell you: how common the species is. The overall rate of records mixes abundance with sampling effort, so prevalence, occupancy probability and abundance are not identifiable (Ward et al. 2009; Fithian et al. 2015). The map is RELATIVE: a value of 2 means twice the landscape-average intensity of occurrence, not a 2% or 20% chance of presence. To estimate absolute occurrence or abundance, add survey data with non-detections (engine: integrated).")
  }
  if (!is.null(o$cv)) {
    f <- o$cv$folds
    out <- c(out, sprintf("Out-of-sample evaluation by spatial block cross-validation: the bias-corrected habitat map ranked held-out records above random background with mean AUC %s (0.5 = no better than random) and continuous Boyce index %s (-1 to 1; > 0 means records are more frequent than expected where predicted intensity is high). Predicting where records occur (habitat x effort) gave AUC %s. Held-out records are themselves biased towards high-effort areas, so the habitat-map scores understate how well it describes the species' distribution.%s",
                          sdm_fmt(mean(f$auc_habitat)), sdm_fmt(mean(f$boyce_habitat, na.rm = TRUE)), sdm_fmt(mean(f$auc_records)),
                          if (any(is.finite(f$auc_survey))) sprintf(" Against held-out survey sites, which are not biased by effort, the AUC was %s.", sdm_fmt(mean(f$auc_survey, na.rm = TRUE))) else ""))
  }
  if (!is.null(o$gof) && o$gof$p < 0.05)
    out <- c(out, sprintf("Model check: records are more clustered among spatial blocks than the model expects (dispersion %.1f, P = %.3f). Intervals are therefore too narrow; possible causes are a missing habitat covariate, sampling effort that the bias covariate does not capture, or repeated records of the same animals/sites.", o$gof$chat, o$gof$p))
  out <- c(out, "All effects are associations between covariates and where the species was recorded, after allowing for the stated sampling bias. They are not evidence that the covariate causes the pattern.")
  out
}

# ---- prediction helpers for figures -----------------------------------------
sdm_band <- function(L, est_row = 1, ci = 0.95, inner = 0.5, trans = exp) {
  D <- L[-est_row, , drop = FALSE]; if (!nrow(D)) D <- L
  q <- function(p) apply(D, 2, stats::quantile, p, na.rm = TRUE)
  a <- (1 - ci) / 2; b <- (1 - inner) / 2
  data.frame(estimate = trans(L[est_row, ]), lower = trans(q(a)), upper = trans(q(1 - a)),
             lower50 = trans(q(b)), upper50 = trans(q(1 - b)))
}

sdm_curve <- function(fit, data, v, opts, kind = c("habitat", "bias"), n = 100) {
  kind <- match.arg(kind)
  o <- fit$object; prep <- o$prep
  rg <- range(prep$grid[[v]])
  xs <- seq(rg[1], rg[2], length.out = n)
  zx <- (xs - prep$scaling[[v]]$center) / prep$scaling[[v]]$scale
  Xn <- matrix(0, n, ncol(prep$X), dimnames = list(NULL, colnames(prep$X)))
  Zn <- matrix(0, n, ncol(prep$Z), dimnames = list(NULL, colnames(prep$Z)))
  if (kind == "habitat") {
    Xn[, v] <- zx; if (paste0(v, "_sq") %in% colnames(Xn)) Xn[, paste0(v, "_sq")] <- zx^2
    L <- sdm_hab_link(o, Xn, Zn)
  } else {
    Zn[, v] <- zx
    L <- sdm_hab_link(o, Xn, Zn, bias = TRUE)   # habitat at means (0) -> relative recording rate
  }
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; inner <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci
  cbind(x = xs, sdm_band(L, ci = ci, inner = inner))
}

sdm_figures <- function(fit, data, opts = fig_opts()) {
  o <- fit$object; prep <- o$prep; est <- fit$estimates
  integ <- "log_density" %in% names(o$eng$coef)
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  cil <- round(100 * ci)
  itype <- if (o$engine == "maxnet") "percentile bootstrap" else "confidence (parametric draws from the ML fit)"
  # maps in km (metre coordinates give unreadable 7-digit tick labels)
  cf_ <- if (identical(data$coord_unit, "km")) 1 else 1e-3
  cu <- "km"
  figs <- list()
  # ---- maps
  link <- drop(sdm_hab_link(o, prep$X, prep$Z, draws = FALSE))
  val <- exp(link)
  if (!integ) val <- val / stats::weighted.mean(val, prep$A)
  grid <- data.frame(x = prep$grid$x * cf_, y = prep$grid$y * cf_, value = val)
  pts <- if (isTRUE(opts$show_data)) data.frame(x = prep$grid$x[prep$pres_cell], y = prep$grid$y[prep$pres_cell]) else NULL
  if (isTRUE(opts$show_data)) pts <- data.frame(x = data$presence$x * cf_, y = data$presence$y * cf_)[!is.na(sdm_cell_index(prep$grid, data$presence$x, data$presence$y)), ]
  flab <- if (integ) fig_label("intensity_map_fill", "Density\n(per km²)", opts) else fig_label("intensity_map_fill", "Relative\nintensity", opts)
  m <- suppressMessages(fig_map(grid, flab, points = pts) +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = flab) +
    ggplot2::labs(x = fig_label("x", paste0("Easting (", cu, ")"), opts), y = fig_label("y", paste0("Northing (", cu, ")"), opts)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)))
  figs$intensity_map <- with_caption(m, paste0(
    if (integ) "Estimated density of the species (individuals per km², log colour scale) " else "Relative intensity of occurrence (1 = landscape average; log colour scale) ",
    "predicted from the habitat covariates with sampling bias removed", if (!is.null(pts)) "; white points are the presence records, which also reflect where people looked." else "."))
  if (ncol(prep$Z)) {
    D <- o$eng$coef[paste0("bias_beta_", colnames(prep$Z))]
    bv <- exp(drop(prep$Z %*% D)); bv <- bv / stats::weighted.mean(bv, prep$A)
    bm <- suppressMessages(fig_map(data.frame(x = prep$grid$x * cf_, y = prep$grid$y * cf_, value = bv),
                                   fig_label("bias_map_fill", "Relative\nrecording\neffort", opts), points = pts) +
      ggplot2::scale_fill_viridis_c(trans = "log10", option = "magma", name = fig_label("bias_map_fill", "Relative\nrecording\neffort", opts)) +
      ggplot2::labs(x = fig_label("x", paste0("Easting (", cu, ")"), opts), y = fig_label("y", paste0("Northing (", cu, ")"), opts)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)))
    figs$bias_map <- with_caption(bm, "Estimated relative recording effort (sampling bias; 1 = landscape average, log colour scale) from the bias covariates, with the presence records.")
  }
  # ---- conditional effects
  ylab <- if (integ) fig_label("density_axis", "Density (individuals per km²)", opts) else fig_label("intensity_axis", "Relative intensity (1 = average)", opts)
  eff <- list()
  for (v in prep$habitat) {
    pr <- sdm_curve(fit, data, v, opts, "habitat")
    obs <- if (isTRUE(opts$show_data)) data.frame(x = prep$grid[[v]][prep$pres_cell]) else NULL
    p <- fig_cond_effect(pr, sdm_lab(v, data, opts), ylab, obs = obs)
    if (!integ) p <- p + ggplot2::geom_hline(yintercept = 1, linetype = 3, colour = eco_pal$null, linewidth = 0.3)
    eff[[v]] <- with_caption(p, sprintf("%s vs %s, other covariates at their landscape means and sampling bias removed. Line: estimate; bands: %d%% and %d%% %s intervals; rug: covariate value at each presence record.",
                                        if (integ) "Expected density" else "Relative intensity of occurrence (1 = intensity at average conditions)",
                                        tolower(sdm_lab(v, data, opts)), cil, round(100 * opts$inner_ci), itype))
    figs[[paste0("effect_", v)]] <- eff[[v]]
  }
  bl <- list()
  for (b in prep$bias) {
    pr <- sdm_curve(fit, data, b, opts, "bias")
    obs <- if (isTRUE(opts$show_data)) data.frame(x = prep$grid[[b]][prep$pres_cell]) else NULL
    bl[[b]] <- fig_cond_effect(pr, sdm_lab(b, data, opts), fig_label("bias_axis", "Relative recording rate", opts), obs = obs) +
      ggplot2::geom_hline(yintercept = 1, linetype = 3, colour = eco_pal$null, linewidth = 0.3)
  }
  if (length(bl)) figs$bias_effect <- with_caption(if (length(bl) == 1) bl[[1]] else compose_figures(lapply(bl, function(p) p + sdm_th(opts))),
    sprintf("Sampling-bias (observer effort) effect: rate at which records are made relative to the mean of the bias covariate, habitat held constant. Bands: %d%% and %d%% %s intervals; rug: values at presence records. This is a property of the observers, not the species.",
            cil, round(100 * opts$inner_ci), itype))
  # ---- evaluation
  if (!is.null(o$cv)) {
    f <- o$cv$folds
    long <- rbind(data.frame(fold = f$fold, metric = "AUC", map = "Habitat map", value = f$auc_habitat),
                  data.frame(fold = f$fold, metric = "AUC", map = "Records map", value = f$auc_records),
                  data.frame(fold = f$fold, metric = "Boyce index", map = "Habitat map", value = f$boyce_habitat),
                  data.frame(fold = f$fold, metric = "Boyce index", map = "Records map", value = f$boyce_records))
    if (any(is.finite(f$auc_survey))) long <- rbind(long, data.frame(fold = f$fold, metric = "AUC", map = "Survey sites", value = f$auc_survey))
    long <- long[is.finite(long$value), ]
    ref <- data.frame(metric = c("AUC", "Boyce index"), r = c(0.5, 0))
    pa <- ggplot2::ggplot(long, ggplot2::aes(value, map)) +
      ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = r), linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_point(colour = eco_pal$muted, size = 1.2, alpha = 0.7) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = eco_pal$accent, colour = eco_pal$ink) +
      ggplot2::facet_wrap(~metric, ncol = 1, scales = "free") + ggplot2::labs(x = "Score on held-out spatial blocks", y = NULL)
    cv <- o$cv$curves; cv$fold <- factor(cv$fold)
    pb <- ggplot2::ggplot(cv, ggplot2::aes(mid, pe, group = fold)) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_line(colour = eco_pal$band, alpha = 0.7, linewidth = 0.5) +
      ggplot2::labs(x = "Predicted log relative intensity (habitat map)", y = "Predicted /\nexpected ratio")
    ev <- compose_figures(list(pa + sdm_th(opts), pb + sdm_th(opts)), ncol = 1, guides = "keep") + patchwork::plot_layout(heights = c(1.3, 1))
    attr(ev, "fig_size") <- list(width = "onehalf", height = 130)
    figs$evaluation <- with_caption(ev, sprintf("Spatial block cross-validation (%d folds). (A) AUC and continuous Boyce index on held-out blocks for the habitat map (sampling bias removed), the records map (habitat x recording effort) and, if present, held-out survey sites (detected vs not); grey points = folds, diamonds = mean, dashed = no skill. (B) Predicted-to-expected ratio of held-out records along the habitat-map gradient, one line per fold (Hirzel et al. 2006); a rising curve means more records than expected where the map predicts high intensity.", nrow(f)))
  }
  if (length(eff)) {
    pl <- lapply(seq_along(eff), function(i) eff[[i]] + sdm_th(opts) + if (i > 1) ggplot2::labs(y = NULL) else NULL)
    rp <- compose_figures(pl, ncol = min(3, length(eff)))
    attr(rp, "fig_size") <- list(width = "double", height = if (length(eff) <= 3) 75 else 140)
    figs$response_curves_panel <- with_caption(rp, sprintf("Response curves for all habitat covariates (y axis: %s), others at their means, sampling bias removed; bands %d%% and %d%% %s intervals.",
                                                           tolower(ylab), cil, round(100 * opts$inner_ci), itype))
  }
  # ---- coefficients (rate ratios per SD)
  cf <- est[grepl("^(intensity_beta_|bias_beta_)", est$param), ]
  cf$group <- ifelse(grepl("^bias", cf$param), "Sampling bias", "Habitat")
  cf$param <- vapply(cf$param, function(p) {
    v <- sub("^(intensity_beta_|bias_beta_)", "", p); sq <- grepl("_sq$", v); v0 <- sub("_sq$", "", v)
    fig_label(p, paste0(sdm_lab(v0, data), if (sq) " (squared)" else ""), opts) }, "")
  figs$coefficients <- with_caption(fig_coef(cf, null_value = 1, transform = exp, log_x = TRUE,
                                             xlab = "Multiplicative change per 1 SD (95% interval)"),
    sprintf("Habitat effects on the intensity of occurrence and sampling-bias effects on the recording rate, as multiplicative changes per 1 SD of each covariate (log scale) with %s intervals; dashed line = no effect.", if (o$engine == "maxnet") "95% percentile bootstrap" else "95% Wald confidence"))
  # ---- checks
  if (!is.null(o$quad)) {
    q <- o$quad$table
    se_all <- if (!is.null(o$eng$vcov)) sqrt(diag(o$eng$vcov)) else apply(o$draws, 2, stats::sd)
    q$full <- o$eng$coef[q$param]; q$se <- se_all[q$param]
    q$param <- vapply(q$param, function(p) {
      v <- sub("^(intensity_beta_|bias_beta_)", "", p)
      if (p == "log_density") "log density" else paste0(sdm_lab(sub("_sq$", "", v), data), if (grepl("_sq$", v)) " (sq.)" else "") }, "")
    qp <- 
      ggplot2::ggplot(q, ggplot2::aes(n_quadrature, estimate)) + ggplot2::geom_line(colour = eco_pal$band) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = full - 0.5 * se, ymax = full + 0.5 * se), fill = eco_pal$band, alpha = 0.12) +
        ggplot2::geom_point(colour = eco_pal$ink, size = 1.2) +
        ggplot2::scale_x_log10(breaks = c(10, 100, 1000, 10000, 1e5), labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)) +
        ggplot2::facet_wrap(~param, scales = "free_y", ncol = 2, labeller = ggplot2::label_wrap_gen(18)) +
        ggplot2::labs(x = "Number of quadrature points (log scale)", y = "Estimate") +
        ggplot2::theme(strip.clip = "off")
    attr(qp, "fig_size") <- list(width = "onehalf", height = 45 * ceiling(length(unique(q$param)) / 2) + 15)
    figs$quadrature_convergence <- with_caption(qp,
      "Quadrature convergence (Renner et al. 2015): coefficient estimates as the quadrature grid is refined from every 8th cell to every cell. Shaded band: full-grid estimate +/- 0.5 standard errors; points inside it at the right mean the integral is well approximated.")
  }
  if (!is.null(o$gof))
    figs$gof <- with_caption(fig_ppc(o$gof$rep, o$gof$obs, "Pearson X² of record counts per spatial block"),
      sprintf("Goodness of fit: observed Pearson statistic (line) for record counts in %d spatial blocks vs %d datasets simulated from the fitted model. An observed value in the far right tail means records are more clustered than the model allows.", o$gof$nblocks, length(o$gof$rep)))
  if (!is.null(o$comparison)) {
    cm <- o$comparison[grepl("^(intensity_beta_|bias_beta_)", o$comparison$param), ]
    cm$group <- cm$engine
    cm$param <- vapply(cm$param, function(p) {
      v <- sub("^(intensity_beta_|bias_beta_)", "", p); sq <- grepl("_sq$", v); v0 <- sub("_sq$", "", v)
      fig_label(p, paste0(sdm_lab(v0, data), if (sq) " (squared)" else "", if (grepl("^bias", p)) " [bias]" else ""), opts) }, "")
    ec <- fig_coef(cm, null_value = 1, transform = exp, log_x = TRUE, xlab = "Multiplicative change per 1 SD (95% interval)")
    attr(ec, "fig_size")$width <- "onehalf"
    figs$engine_comparison <- with_caption(ec,
      "The same effects estimated by different engines on the same data, as multiplicative changes per 1 SD (log scale); dashed line = no effect. ppm and maxnet use the presence-only records only; integrated adds the survey sites. maxnet intervals are percentile bootstrap, the others Wald.")
  }
  figs
}

sdm_methods_text <- function(fit, data, settings) {
  o <- fit$object; prep <- o$prep; s <- o$settings
  hab <- paste(vapply(prep$habitat, sdm_lab, "", data = data), collapse = ", ")
  quad <- if (length(prep$quadratic)) paste0(" with quadratic terms for ", paste(vapply(prep$quadratic, sdm_lab, "", data = data), collapse = ", ")) else ""
  bias <- if (length(prep$bias)) paste(vapply(prep$bias, sdm_lab, "", data = data), collapse = ", ") else NULL
  core <- sprintf("We modelled the %d presence records as an inhomogeneous Poisson point process (Warton & Shepherd 2010; Renner et al. 2015) over a %s km2 study region represented by %d grid cells. The log intensity of the species was a linear function of %s%s, each standardised to mean 0 and SD 1 over the region.",
                  prep$n_pres, format(round(sum(prep$A)), big.mark = ","), nrow(prep$X), hab, quad)
  b <- if (!is.null(bias)) sprintf(" To account for uneven sampling effort, the presence records were modelled as a thinned version of this process whose log thinning rate was a linear function of %s (Warton et al. 2013; Fithian et al. 2015); predictions and habitat effects hold the bias covariate constant.", bias)
       else " No sampling-bias covariate was included, so effort is assumed uniform across the region."
  e <- switch(o$engine,
    ppm = " The model was fitted by down-weighted Poisson regression using every grid cell centre as a quadrature point weighted by cell area (Renner et al. 2015), with Wald 95% confidence intervals.",
    maxnet = sprintf(" The model was fitted with maxnet (Phillips et al. 2017), the open-source implementation of Maxent, which is equivalent to an inhomogeneous Poisson process (Fithian & Hastie 2013; Renner & Warton 2013), using '%s' features, regularisation multiplier %g and all grid cells as background; 95%% intervals are percentile intervals from %d bootstrap resamples of the presence records.",
                     o$eng$classes, s$maxnet_regmult, nrow(o$eng$boot)),
    integrated = sprintf(" Presence-only records were combined with %d survey sites (%s) in a joint likelihood (Fithian et al. 2015; Koshkina et al. 2017; Dorazio 2014): a survey plot of area A was occupied with probability 1 - exp(-lambda A) (complementary log-log link)%s. The integrated likelihood was maximised numerically and 95%% confidence intervals for derived quantities (density at mean covariates, expected abundance, mean plot occupancy) were obtained from %d draws from the asymptotic normal distribution of the estimates.",
                         length(prep$survey$cell), if (isTRUE(prep$survey$estimate_p)) "repeated visits" else "single visits treated as presence/absence",
                         if (isTRUE(prep$survey$estimate_p)) " and the species was detected on each visit to an occupied plot with constant probability p" else "", s$draws),
    survey_only = " Only the survey data were used: occupancy with a complementary log-log link to the Poisson intensity.")
  ev <- if (isTRUE(s$evaluate)) sprintf(" Predictive performance was assessed by spatial block cross-validation (%d folds of square blocks; Roberts et al. 2017) using AUC and the continuous Boyce index (Hirzel et al. 2006).", s$cv_folds) else ""
  ck <- paste0(if (isTRUE(s$quadrature_check)) " Convergence of the quadrature approximation was checked by refitting with coarser quadrature grids." else "",
               if (isTRUE(s$gof)) " Goodness of fit was assessed by comparing record counts in spatial blocks with parametric replicates from the fitted model." else "")
  paste0(core, b, e, ev, ck, " The method was validated by simulation (ecotools sdm-presence-only VALIDATION.md). Analyses used R ", getRversion(), ".")
}

sdm_references <- c(
  "Warton DI, Shepherd LC (2010) Poisson point process models solve the 'pseudo-absence problem' for presence-only data in ecology. Annals of Applied Statistics 4:1383-1402.",
  "Renner IW, Elith J, Baddeley A, Fithian W, Hastie T, Phillips SJ, Popovic G, Warton DI (2015) Point process models for presence-only analysis. Methods in Ecology and Evolution 6:366-379. doi:10.1111/2041-210X.12352",
  "Fithian W, Hastie T (2013) Finite-sample equivalence in statistical models for presence-only data. Annals of Applied Statistics 7:1917-1939.",
  "Renner IW, Warton DI (2013) Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics 69:274-281.",
  "Phillips SJ, Anderson RP, Dudik M, Schapire RE, Blair ME (2017) Opening the black box: an open-source release of Maxent. Ecography 40:887-893. doi:10.1111/ecog.03049",
  "Warton DI, Renner IW, Ramp D (2013) Model-based control of observer bias for the analysis of presence-only data in ecology. PLoS ONE 8:e79168.",
  "Fithian W, Elith J, Hastie T, Keith DA (2015) Bias correction in species distribution models: pooling survey and collection data for multiple species. Methods in Ecology and Evolution 6:424-438. doi:10.1111/2041-210X.12242",
  "Koshkina V, Wang Y, Gordon A, Dorazio RM, White M, Stone L (2017) Integrated species distribution models: combining presence-background data and site-occupancy data with imperfect detection. Methods in Ecology and Evolution 8:420-430. doi:10.1111/2041-210X.12738",
  "Dorazio RM (2014) Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. Global Ecology and Biogeography 23:1472-1484. doi:10.1111/geb.12216",
  "Ward G, Hastie T, Barry S, Elith J, Leathwick JR (2009) Presence-only data and the EM algorithm. Biometrics 65:554-563.",
  "Berman M, Turner TR (1992) Approximating point process likelihoods with GLIM. Applied Statistics 41:31-38.",
  "Hirzel AH, Le Lay G, Helfer V, Randin C, Guisan A (2006) Evaluating the ability of habitat suitability models to predict species presences. Ecological Modelling 199:142-152.",
  "Roberts DR, Bahn V, Ciuti S, et al. (2017) Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography 40:913-929.",
  "Baddeley A, Rubak E, Turner R (2015) Spatial Point Patterns: Methodology and Applications with R. Chapman & Hall/CRC (spatstat; independent cross-check).")
