# Output for movement-ssf: interpretation, methods text, figures, data writer.

ssf_fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
ssf_row <- function(est, p) est[est$param == p, , drop = FALSE]
ssf_lab <- function(data, v, opts = NULL) fig_label(v, data$labels[[v]] %||% v, opts)

ssf_interpret <- function(fit, data) {
  est <- fit$estimates
  an <- data$animal_label %||% "animal"
  ans <- data$animal_plural %||% (if (an %in% c("elk", "deer", "moose", "caribou", "bison", "sheep", "fish")) an else paste0(an, "s"))
  out <- character()
  if (is.null(fit$issa)) {
    return("This fit used the resource selection function only (method: rsf); see the estimates table. RSF coefficients describe selection at the home-range scale and their intervals ignore autocorrelation between successive fixes.")
  }
  h <- fit$interval_hours
  for (v in data$covariates) {
    m <- ssf_row(est, paste0("sel_", v, "_mean")); s <- ssf_row(est, paste0("sel_", v, "_sd"))
    r <- ssf_row(est, paste0("rss_", v, "_per_sd")); ru <- ssf_row(est, paste0("rss_", v, "_per_unit"))
    sc <- data$land$scaling[[v]]$scale
    nm <- data$names[[v]] %||% v
    dir <- if (m$lower > 0) "selected locations with higher" else if (m$upper < 0) "avoided locations with higher" else "showed no clear selection or avoidance of"
    txt <- sprintf("%s: %s %s %s at the scale of one %s-h step (selection coefficient %s per SD, 95%% CI %s to %s). Of two locations an animal could equally reach that differ by 1 SD of %s (SD = %s in the raster's units), the one with higher %s is %s times as likely to be chosen (relative selection strength, RSS; 95%% CI %s to %s).",
                   ssf_lab(data, v), tools::toTitleCase(ans), dir, nm, format(h), ssf_fmt(m$estimate), ssf_fmt(m$lower), ssf_fmt(m$upper),
                   nm, ssf_fmt(sc, 1), nm, ssf_fmt(r$estimate), ssf_fmt(r$lower), ssf_fmt(r$upper))
    if (nrow(ru)) txt <- paste(txt, sprintf("Per %s, RSS = %s (%s to %s).", data$unit_text[[v]] %||% paste(data$units[[v]], "units"),
                                            ssf_fmt(ru$estimate), ssf_fmt(ru$lower), ssf_fmt(ru$upper)))
    txt <- paste(txt, sprintf("In the RSS figure for %s, each value is compared with an equally reachable location at the landscape-average value (%s); the curve is for a typical animal (population mean), shown on a log scale.",
                              nm, ssf_fmt(data$land$scaling[[v]]$center, 1)))
    if (nrow(s)) {
      pr <- stats::pnorm(m$estimate / max(s$estimate, 1e-8))
      txt <- paste(txt, sprintf("Individuals differ: the among-animal SD of the selection coefficient is %s (95%% CI %s to %s) on the log-RSS scale per SD, so an estimated %s of %s in this population would show positive selection (log-RSS > 0) for %s.",
                                ssf_fmt(s$estimate), ssf_fmt(s$lower), ssf_fmt(s$upper),
                                if (pr > 0.99) "more than 99%" else if (pr < 0.01) "less than 1%" else sprintf("%.0f%%", 100 * pr), ans, nm))
    }
    out <- c(out, txt)
  }
  ha <- ssf_row(est, "home_attraction_per_km2")
  if (nrow(ha)) out <- c(out, sprintf("Range residency: step ends farther from each animal's centre of activity were chosen less often (attraction coefficient %s per km^2, 95%% CI %s to %s). This term is included so that the habitat coefficients do not absorb the animals' tendency to stay in their ranges; treat its size only as an indication. With the centre estimated from the same track it is overestimated (in this module's validation by about 25%% with 400 steps per animal, 6%% with 1500 steps) and its interval is too narrow, so it is not interpreted quantitatively.",
                                      ssf_fmt(ha$estimate, 3), ssf_fmt(ha$lower, 3), ssf_fmt(ha$upper, 3)))
  mm <- ssf_row(est, "step_mean_m"); k <- ssf_row(est, "turn_kappa")
  if (nrow(mm)) {
    obs <- mean(fit$steps$sl_[fit$steps$case == 1])
    out <- c(out, sprintf("Movement: after removing the effect of habitat selection, the selection-free step length (the distance an animal would move in %s h in homogeneous habitat) has mean %s m (95%% CI %s to %s; gamma shape %s, scale %s m), compared with an observed mean step of %s m. The turn-angle concentration is %s (95%% CI %s to %s); values above 0 mean animals tend to keep their heading (directional persistence), 0 means random turning.",
                            format(h), ssf_fmt(mm$estimate, 0), ssf_fmt(mm$lower, 0), ssf_fmt(mm$upper, 0),
                            ssf_fmt(ssf_row(est, "step_shape")$estimate), ssf_fmt(ssf_row(est, "step_scale_m")$estimate, 0), ssf_fmt(obs, 0),
                            ssf_fmt(k$estimate), ssf_fmt(k$lower), ssf_fmt(k$upper)))
  }
  if (!is.null(fit$rsf)) {
    rr <- fit$rsf$estimates
    txt <- vapply(data$covariates, function(v) { r <- rr[rr$param == paste0("rsf_", v, "_mean"), ]
      sprintf("%s %s (95%% CI %s to %s)", v, ssf_fmt(r$estimate), ssf_fmt(r$lower), ssf_fmt(r$upper)) }, "")
    out <- c(out, sprintf("Resource selection function (use vs. availability inside each animal's minimum convex polygon, ignoring movement): coefficients per SD %s. These describe selection at the home-range scale and are not interchangeable with the step-scale coefficients above. Their intervals assume independent fixes; with autocorrelated GPS data they are too narrow (in this module's validation, with 12 animals and 2-h fixes, the RSF intervals covered the true step-selection coefficient in only 52%% of replicates for elevation and 91%% for forest cover; RSF and step-scale coefficients estimate different things), so use the iSSA results for inference.",
                          paste(txt, collapse = "; ")))
  }
  if (!is.null(fit$home_ranges)) {
    ht <- fit$home_ranges$table; ok <- is.finite(ht$akde_km2)
    if (any(ok)) out <- c(out, sprintf("Home ranges (%.0f%% AKDE, which accounts for autocorrelation): median %s km^2 across %d %s (range %s to %s km^2). Conventional estimators that assume independent fixes give %s: median KDE %s km^2 and MCP %s km^2. AKDE is the recommended estimate (Fleming et al. 2015); its effective sample size per animal (median %.0f) shows how many independent range crossings the track contains; below ~10 the area is very uncertain. In this module's validation (simulated range-resident tracks with about 13 effective range crossings) AKDE areas were on average 9%% too small and the 95%% intervals contained the true area in 85%% of simulations, so treat the AKDE intervals as somewhat too narrow for short deployments.",
                                       100 * fit$home_ranges$level, ssf_fmt(stats::median(ht$akde_km2[ok]), 1), sum(ok), ans,
                                       ssf_fmt(min(ht$akde_km2[ok]), 1), ssf_fmt(max(ht$akde_km2[ok]), 1),
                                       { mk <- stats::median(ht$kde_km2, na.rm = TRUE); mm <- stats::median(ht$mcp_km2, na.rm = TRUE); ma <- stats::median(ht$akde_km2[ok])
                                         if (isTRUE(mk < ma && mm < ma)) "smaller areas, as expected with autocorrelated fixes" else if (isTRUE(mk > ma && mm > ma)) "larger areas here (unusual; check the AKDE model and effective sample size)" else "different areas" },
                                       ssf_fmt(stats::median(ht$kde_km2, na.rm = TRUE), 1), ssf_fmt(stats::median(ht$mcp_km2, na.rm = TRUE), 1),
                                       stats::median(ht$n_eff[ok], na.rm = TRUE)))
  }
  if (!is.null(fit$uhc)) {
    t <- fit$uhc$table
    out <- c(out, paste0("Model check (used-habitat calibration on held-out steps): ", paste(sprintf("%s observed mean at used steps %s vs model-predicted %s (95%% interval %s to %s)%s",
      t$covariate, ssf_fmt(t$observed_used_mean, 1), ssf_fmt(t$predicted_mean, 1), ssf_fmt(t$pred_lower, 1), ssf_fmt(t$pred_upper, 1),
      ifelse(t$ok, "", " -- OUTSIDE")), collapse = "; "), "."))
  }
  nf <- data$n_failed_fixes
  if (!is.null(nf) && nf > 0) out <- c(out, sprintf("Fix success: %d of %d fix attempts (%.0f%%) had no location and were dropped, which breaks the track into separate bursts. The analysis assumes fixes were lost independently of habitat. If collars fail more often under dense canopy or in steep terrain, steps into those places are under-recorded and selection for them is underestimated (a detection problem, not animal behaviour); compare fix success across habitats before interpreting a weak or negative effect of such a covariate.",
                                           nf, data$n_fix_attempts %||% NA_integer_, 100 * nf / (data$n_fix_attempts %||% NA_integer_)))
  out <- c(out, "These are associations between where animals stepped and the mapped covariates, conditional on what was reachable; they are not proof that the covariate causes the choice (unmapped features correlated with the covariate produce the same signal).")
  out
}

ssf_methods_text <- function(fit, data, settings) {
  st <- fit$settings %||% settings
  rs <- !is.null(fit$issa) && isTRUE(fit$issa$random_slopes)
  paste0(
    "GPS fixes were regularised to a ", format(fit$interval_hours %||% NA), "-h interval (tolerance ", fit$tolerance_minutes %||% st$tolerance_minutes %||% "10% of the interval",
    " min) and split into bursts of consecutive fixes with amt (Signer et al. 2019). For each observed step we drew ", st$n_control,
    " available steps from a gamma step-length and von Mises turn-angle distribution fitted to all observed steps (the tentative movement kernel). ",
    "Covariates (", paste(data$covariates, collapse = ", "), ") were extracted at step end points and standardised by the mean and SD of each raster across the study landscape. ",
    "We fitted an integrated step-selection function (Avgar et al. 2016) with step length, log step length and the cosine of the turn angle as movement terms",
    if (isTRUE(st$home_attraction)) ", plus the squared distance (km^2) of the step end from each animal's centre of activity to represent range residency (the centre was estimated jointly through animal-specific linear terms in the coordinates)" else "",
    ". The conditional logistic likelihood was fitted as a Poisson model with step-specific intercepts whose variance was fixed at 10^6",
    if (rs) " and independent animal-level random slopes for each habitat covariate (Muff et al. 2020)" else "",
    " in glmmTMB (Brooks et al. 2017)", if (rs && isTRUE(st$reml)) ", with variance components estimated by restricted maximum likelihood (REML), which in validation removed most of the downward bias of the maximum-likelihood among-animal SD" else "", ". Selection is reported as relative selection strength (RSS = exp(beta); Avgar et al. 2017) per SD and per a biologically meaningful unit, with Wald 95% confidence intervals; among-animal SDs have intervals from the log-SD scale. ",
    "Selection-free movement parameters were obtained by updating the tentative gamma shape and scale and von Mises concentration with the movement coefficients (Avgar et al. 2016; Fieberg et al. 2021), with delta-method intervals. ",
    if (!is.null(fit$rsf)) paste0("For comparison we fitted a use-availability resource selection function with ", st$rsf_n_avail, " random available points per fix inside each animal's minimum convex polygon, weighted logistic regression (available weight ", st$rsf_weight, ") with a random intercept of fixed large variance and random slopes (Muff et al. 2020; Fieberg et al. 2021). ") else "",
    if (!is.null(fit$home_ranges)) paste0("Home ranges (", 100 * (st$hr_level %||% 0.95), "%) were estimated per animal by autocorrelated kernel density estimation (AKDE) with the continuous-time movement model chosen by AICc in ctmm (Fleming et al. 2015; Calabrese et al. 2016), and for comparison by Gaussian-reference KDE and MCP in amt. ") else "",
    if (!is.null(fit$uhc)) "Model fit was checked with out-of-sample used-habitat calibration (Fieberg et al. 2018) on 20% of steps held out at random. " else "",
    "The analysis was validated by simulation (ecotools movement-ssf; see VALIDATION.md).")
}

# ---- figures -----------------------------------------------------------------

# Back-transform a standardised covariate to original units. A user `figures: units:`
# entry ({center, scale}) wins over the landscape scaling recorded by read_data/simulate.
ssf_orig <- function(data, v, z, opts = NULL) {
  u <- opts$units[[v]]
  if (!is.null(u) && !is.list(u)) opts$units[[v]] <- NULL   # data `units:` (RSS change) is not a scaling
  fig_x(z, v, opts, data$land$scaling)
}

# greyscale: fig_cond_effect draws the single-series line/band in the green accent whatever
# the palette; recolour to ink/grey so a print-safe export carries no colour.
ssf_grey <- function(p, opts) {
  if (!isTRUE(opts$greyscale) && !identical(opts$palette, "grey")) return(p)
  for (i in seq_along(p$layers)) {
    g <- class(p$layers[[i]]$geom)[1]
    if (g == "GeomRibbon" && !is.null(p$layers[[i]]$aes_params$fill)) p$layers[[i]]$aes_params$fill <- "grey45"
    if (g == "GeomLine" && !is.null(p$layers[[i]]$aes_params$colour)) p$layers[[i]]$aes_params$colour <- eco_pal$ink
  }
  p
}

# Size hint with a default aspect. The core's apply_fig_opts ignores an `aspect` hint (only a
# per_figure aspect or the global one is used), so pass a height computed from the width;
# per_figure width/aspect/height overrides still win.
ssf_size <- function(name, width, aspect, opts) {
  pf <- opts$per_figure[[name]]
  w <- fig_width_mm(pf$width %||% width, opts)
  list(width = width, height = w * (pf$aspect %||% aspect))
}

# Figure suite (docs/FIGURES.md, movement-ssf): rss_<cov>, selection_coefficients,
# movement_kernel, used_available, home_ranges, hr_area (+ uhc, the model check).
# Fixed colour roles shared by every figure: Used = vermillion, Available = blue,
# population mean = ink, individuals = muted grey; AKDE/KDE/MCP = Okabe-Ito 1-3.
ssf_figures <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2", "patchwork"))
  if (is.null(fit$issa)) return(list())
  ci <- opts$ci %||% 0.95; z <- stats::qnorm(1 - (1 - ci) / 2); cip <- round(100 * ci)
  L <- function(key, def) fig_label(key, def, opts)
  covlab <- function(v) vapply(v, function(x) L(x, data$labels[[x]] %||% x), "", USE.NAMES = FALSE)
  est <- fit$estimates; b <- fit$issa$coef; V <- fit$issa$vcov
  st <- fit$steps
  covs <- data$covariates
  out <- list()
  an <- L("animal", data$animal_label %||% "animal")
  an_a <- paste(if (grepl("^[aeiouAEIOU]", an)) "an" else "a", an)
  an_pl <- L("animals", data$animal_plural %||% if (tolower(an) %in% c("elk", "deer", "moose", "caribou", "sheep", "bison")) an else paste0(an, "s"))
  lab_used <- L("used", "Used"); lab_avail <- L("available", "Available")
  # Used/available colours honour greyscale and user palettes (linetype is the second encoding)
  pal <- opts$palette %||% "okabe_ito"
  ua_pair <- if (identical(pal, "grey")) eco_pal$grey[c(1, 2)] else if (length(pal) > 1) rep_len(pal, 2) else eco_pal$okabe_ito[c(6, 5)]
  ua_cols <- stats::setNames(ua_pair, c(lab_used, lab_avail))
  cov_levels <- covlab(covs)   # panels follow the covariate order of the config (as the rss_<cov> figures)
  ylab_dens <- L("density", "Probability density")
  show <- isTRUE(opts$show_data %||% TRUE)
  # RSS curves (Avgar et al. 2017): location x vs a reference location at the landscape mean
  for (v in covs) {
    k <- paste0("z_", v)
    zr <- stats::quantile(st[[k]], c(0.01, 0.99), na.rm = TRUE)
    zz <- seq(zr[1], zr[2], length.out = 120)
    se <- sqrt(V[k, k])
    pred <- data.frame(x = ssf_orig(data, v, zz, opts), estimate = exp(b[[k]] * zz),
                       lower = exp(zz * b[[k]] - z * abs(zz) * se), upper = exp(zz * b[[k]] + z * abs(zz) * se))
    used <- st[st$case == 1, ]
    uo <- data.frame(x = ssf_orig(data, v, used[[k]], opts))
    uo <- uo[uo$x >= min(pred$x) & uo$x <= max(pred$x), , drop = FALSE]
    if (nrow(uo) > 400) { set.seed(1); uo <- uo[sort(sample.int(nrow(uo), 400)), , drop = FALSE] }   # readable rug
    p <- ssf_grey(fig_cond_effect(pred, covlab(v), L("rss", "Relative selection strength (RSS)"),
                         obs = if (show) uo else NULL), opts) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
      ggplot2::scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.2, 0.5, 1, 2, 5, 10, 30, 100), labels = function(x) format(x, drop0trailing = TRUE))
    out[[paste0("rss_", v)]] <- with_caption(p, sprintf(
      "Relative selection strength (Avgar et al. 2017) for %s: how many times more likely %s is to end a step at a location with the value on the x axis than at an equally reachable location with the landscape-average value (%s), other covariates equal. Line: population-level (typical %s) estimate; band: %d%% confidence interval; log scale, so the line is the log-RSS; dashed line = no selection (RSS = 1). Curve spans the 1st-99th percentile of the covariate at used and available steps. Rug: covariate at a random sample of up to 400 used step ends.",
      tolower(data$names[[v]] %||% covlab(v)), an_a, ssf_fmt(data$land$scaling[[v]]$center, 1), an, cip))
  }
  # Selection coefficients: population mean + individual
  lab_pop <- L("population_mean", "Population mean"); lab_ind <- L("individual", "Individual")
  pop <- do.call(rbind, lapply(covs, function(v) {
    r <- ssf_row(est, paste0("sel_", v, "_mean"))
    data.frame(covariate = covlab(v), id = lab_pop, estimate = r$estimate,
               lower = r$estimate - z * r$sd, upper = r$estimate + z * r$sd, level = lab_pop)
  }))
  ind <- fit$issa$indiv
  if (!is.null(ind)) {
    zz <- z / stats::qnorm(0.975)
    ind$lower <- ind$estimate - (ind$estimate - ind$lower) * zz; ind$upper <- ind$estimate + (ind$upper - ind$estimate) * zz
    ind$covariate <- covlab(ind$covariate)
    ind$id <- covlab(as.character(ind$id))
    ind$level <- lab_ind
  }
  dd <- rbind(pop, if (!is.null(ind)) ind[, names(pop)])
  dd$id <- factor(dd$id, levels = c(rev(sort(unique(setdiff(dd$id, lab_pop)))), lab_pop))
  dd$level <- factor(dd$level, levels = c(lab_pop, lab_ind))
  dd$covariate <- factor(dd$covariate, levels = cov_levels)
  p <- ggplot2::ggplot(dd, ggplot2::aes(estimate, id, shape = level)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper, colour = level), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = level), size = 1.8) +
    ggplot2::facet_wrap(~covariate, scales = "free_x") +
    ggplot2::scale_colour_manual(values = stats::setNames(c(eco_pal$ink, eco_pal$muted), c(lab_pop, lab_ind))) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(18, 16), c(lab_pop, lab_ind))) +
    ggplot2::labs(x = L("sel_axis", "Selection coefficient (log-RSS per 1 SD of covariate)"), y = NULL, colour = NULL, shape = NULL)
  attr(p, "fig_size") <- list(width = "double", height = 30 + 4.5 * length(unique(dd$id)))
  sd_txt <- paste(vapply(covs, function(v) {
    lb <- covlab(v); u <- if (grepl("\\(.*\\)\\s*$", lb)) paste0(" ", sub(".*\\((.*)\\)\\s*$", "\\1", lb)) else ""
    sprintf("1 SD of %s = %s%s", data$names[[v]] %||% tolower(lb), ssf_fmt(data$land$scaling[[v]]$scale, 1), u) }, ""), collapse = "; ")
  out$selection_coefficients <- with_caption(p, sprintf(
    "Step-selection coefficients (log-RSS per 1 SD of each covariate; %s). Diamond: population mean with %d%% confidence interval; dots: individual %s (population mean plus the animal's conditional random slope) with approximate %d%% intervals. Dashed line: no selection.",
    sd_txt, cip, an_pl, cip))
  # Movement kernel: tentative vs updated
  slmax <- stats::quantile(st$sl_[st$case == 1], 0.99)
  xs <- seq(1, slmax, length.out = 300)
  sh <- ssf_row(est, "step_shape"); sc <- ssf_row(est, "step_scale_m"); kp <- ssf_row(est, "turn_kappa")
  k_tent <- L("kernel_tentative", "Tentative (fitted to observed steps)"); k_upd <- L("kernel_updated", "Selection-free (updated)")
  kl <- c(k_upd, k_tent)
  sld <- data.frame(x = rep(xs, 2), density = c(stats::dgamma(xs, fit$sl_distr$shape, scale = fit$sl_distr$scale),
                                                  if (nrow(sh)) stats::dgamma(xs, sh$estimate, scale = sc$estimate) else NA),
                    kernel = factor(rep(c(k_tent, k_upd), each = length(xs)), levels = kl))
  ang <- seq(-pi, pi, length.out = 300)
  dvm <- function(a, k) exp(k * cos(a)) / (2 * pi * besselI(k, 0))
  tad <- data.frame(x = rep(ang, 2), density = c(dvm(ang, fit$ta_distr$kappa), dvm(ang, max(kp$estimate, 0))),
                    kernel = factor(rep(c(k_tent, k_upd), each = length(ang)), levels = kl))
  obs_sl <- data.frame(x = st$sl_[st$case == 1 & st$sl_ <= slmax]); obs_ta <- data.frame(x = st$ta_[st$case == 1])
  kcols <- stats::setNames(eco_categorical(2, opts$palette %||% "okabe_ito"), kl)
  mk <- function(d, o, xl, bins, yl) ggplot2::ggplot(d, ggplot2::aes(x, density)) +
    ggplot2::geom_histogram(data = o, ggplot2::aes(x = x, y = ggplot2::after_stat(density)), inherit.aes = FALSE,
                            bins = bins, fill = "grey85", colour = "white", linewidth = 0.1) +
    ggplot2::geom_line(ggplot2::aes(colour = kernel, linetype = kernel), linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = kcols) + ggplot2::scale_linetype_manual(values = stats::setNames(c(1, 2), kl)) +
    ggplot2::labs(x = xl, y = yl, colour = NULL, linetype = NULL)
  th <- function(p) p + theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic", "bottom")
  out$movement_kernel <- with_caption(compose_figures(list(
    th(mk(sld, obs_sl, L("step_length", sprintf("Step length (m per %s h)", format(fit$interval_hours))), 40, L("density_sl", "Probability density (per m)"))),
    th(mk(tad, obs_ta, L("turn_angle", "Turn angle (radians)"), 36, L("density_ta", "Probability density (per radian)")) +
         ggplot2::guides(colour = "none", linetype = "none") +
         ggplot2::scale_x_continuous(breaks = c(-pi, -pi / 2, 0, pi / 2, pi), labels = c("-π", "-π/2", "0", "π/2", "π")))),
    ncol = 2) & ggplot2::theme(legend.position = "bottom"),
    "Movement kernel. A: step length; B: turn angle (0 = straight ahead). Grey histograms: observed steps. Dashed: tentative gamma (step) and von Mises (turn) distributions fitted to the observed steps and used to draw available steps; solid: selection-free (updated) distributions, which remove the effect of habitat selection on observed movement (Avgar et al. 2016). Point estimates; step lengths above the 99th percentile are not shown.")
  attr(out$movement_kernel, "fig_size") <- ssf_size("movement_kernel", "double", 0.45, opts)
  # Used vs available
  ua <- do.call(rbind, lapply(covs, function(v) data.frame(covariate = covlab(v),
                                                           value = ssf_orig(data, v, st[[paste0("z_", v)]], opts),
                                                           type = factor(ifelse(st$case == 1, lab_used, lab_avail), levels = c(lab_used, lab_avail)))))
  ua$covariate <- factor(ua$covariate, levels = cov_levels)
  out$used_available <- with_caption(
    ggplot2::ggplot(ua, ggplot2::aes(value, colour = type, linetype = type)) +
      ggplot2::geom_density(linewidth = 0.7, key_glyph = "path") + ggplot2::facet_wrap(~covariate, scales = "free") +
      ggplot2::scale_colour_manual(values = ua_cols) +
      ggplot2::scale_linetype_manual(values = stats::setNames(c(1, 2), c(lab_used, lab_avail))) +
      ggplot2::labs(x = L("covariate_value", "Covariate value (units in panel title)"), y = ylab_dens, colour = NULL, linetype = NULL),
    sprintf("Distribution (kernel density) of each covariate at used step ends (solid: observed GPS steps) and at available step ends (dashed: %d per used step, drawn from the tentative movement kernel). A shift of the used curve relative to the available curve is the raw signal of selection; the y scale differs between panels because covariate units differ.",
            fit$settings$n_control %||% 20L))
  attr(out$used_available, "fig_size") <- ssf_size("used_available", "double", 0.45, opts)
  # UHC
  if (!is.null(fit$uhc)) {
    u <- fit$uhc
    sim <- u$sim; sim$covariate <- covlab(sim$covariate)
    us <- u$used; us$covariate <- covlab(us$covariate); av <- u$avail; av$covariate <- covlab(av$covariate)
    sim$covariate <- factor(sim$covariate, levels = cov_levels); us$covariate <- factor(us$covariate, levels = cov_levels)
    av$covariate <- factor(av$covariate, levels = cov_levels)
    lab_pred <- L("uhc_predicted", "Model-predicted use"); lab_uho <- L("uhc_used", "Used (held out)")
    lt <- stats::setNames(c(1, 1, 2), c(lab_pred, lab_uho, lab_avail))
    cl <- stats::setNames(c(if (identical(pal, "grey")) eco_pal$grey[4] else "grey70", ua_pair), c(lab_pred, lab_uho, lab_avail))
    p <- ggplot2::ggplot() +
      ggplot2::geom_density(data = sim, ggplot2::aes(value, group = draw, colour = lab_pred, linetype = lab_pred), linewidth = 0.2, key_glyph = "path") +
      ggplot2::geom_density(data = av, ggplot2::aes(value, colour = lab_avail, linetype = lab_avail), linewidth = 0.6, key_glyph = "path") +
      ggplot2::geom_density(data = us, ggplot2::aes(value, colour = lab_uho, linetype = lab_uho), linewidth = 0.8, key_glyph = "path") +
      ggplot2::facet_wrap(~covariate, scales = "free") +
      ggplot2::scale_linetype_manual(values = lt, breaks = names(lt)) + ggplot2::scale_colour_manual(values = cl, breaks = names(cl)) +
      ggplot2::labs(x = L("covariate_value", "Covariate value (units in panel title)"), y = ylab_dens, linetype = NULL, colour = NULL)
    out$uhc <- with_caption(p, sprintf("Used-habitat calibration (Fieberg et al. 2018) on %d strata (about 20%% of steps) held out from fitting. Grey: covariate distributions at steps chosen by the fitted model (one line per parameter draw, 50 of %d draws shown); solid colour: observed used steps; dashed: available steps. A well-calibrated model has the observed used curve inside the grey envelope.",
                                       u$n_test %||% NA, fit$settings$uhc_draws %||% 200))
    attr(out$uhc, "fig_size") <- ssf_size("uhc", "double", 0.45, opts)
  }
  # Home ranges: one panel per animal (facets, as there can be more animals than palette colours)
  H <- fit$home_ranges
  lab_e <- L("easting", "Easting (km)"); lab_n <- L("northing", "Northing (km)")
  if (!is.null(H) && !is.null(H$polygons) && nrow(H$polygons)) {
    pg <- H$polygons; pg$grp <- interaction(pg$id, pg$grp, drop = TRUE)
    lab_est <- L("akde_estimate", sprintf("%.0f%% AKDE estimate", 100 * H$level)); lab_ci <- L("akde_ci", "95% CI of the contour")
    pg$ci <- factor(ifelse(pg$ci == "estimate", lab_est, lab_ci), levels = c(lab_est, lab_ci))
    fx <- data$fixes[stats::complete.cases(data$fixes), ]
    fx <- fx[fx$id %in% pg$id, ]
    pg$id <- covlab(as.character(pg$id)); fx$id <- covlab(as.character(fx$id))
    nid <- length(unique(pg$id)); nc <- if (nid <= 4) nid else if (nid <= 9) 3 else if (nid <= 16) 4 else 5
    nr <- ceiling(nid / nc)
    p <- ggplot2::ggplot() +
      (if (show) ggplot2::geom_point(data = fx, ggplot2::aes(x / 1000, y / 1000), size = 0.08, colour = "grey55", alpha = 0.5)) +
      ggplot2::geom_polygon(data = pg, ggplot2::aes(X / 1000, Y / 1000, group = grp, linetype = ci), fill = NA, colour = eco_pal$ink, linewidth = 0.35) +
      ggplot2::facet_wrap(~id, ncol = nc) +
      ggplot2::scale_linetype_manual(values = stats::setNames(c(1, 3), c(lab_est, lab_ci))) +
      ggplot2::scale_x_continuous(n.breaks = 4) + ggplot2::scale_y_continuous(n.breaks = 4) +
      ggplot2::coord_equal() +
      ggplot2::labs(x = lab_e, y = lab_n, linetype = NULL)
    asp <- diff(range(c(pg$Y, fx$y))) / diff(range(c(pg$X, fx$x)))
    attr(p, "fig_size") <- list(width = "double", height = min(230, 22 + nr * (174 - 45) / nc * asp))
    out$home_ranges <- with_caption(p, sprintf("%.0f%% autocorrelated kernel density (AKDE) home range of each %s (solid) with the 95%% confidence contours of that home range (dotted); grey points are the %s's GPS fixes. All panels share the same scale. Coordinates in km (%s).",
                                               100 * H$level, an, an, data$land$crs))
  }
  if (!is.null(H)) {
    ht <- H$table
    mth <- c(L("akde", "AKDE"), L("kde", "KDE (href)"), L("mcp", "MCP"))
    hd <- rbind(data.frame(id = ht$id, method = mth[1], estimate = ht$akde_km2, lower = ht$akde_lower, upper = ht$akde_upper),
                data.frame(id = ht$id, method = mth[2], estimate = ht$kde_km2, lower = NA, upper = NA),
                data.frame(id = ht$id, method = mth[3], estimate = ht$mcp_km2, lower = NA, upper = NA))
    hd$method <- factor(hd$method, levels = mth)
    hd$id <- covlab(as.character(hd$id))
    hd$id <- factor(hd$id, levels = rev(sort(unique(hd$id))))
    p <- ggplot2::ggplot(hd, ggplot2::aes(estimate, id, colour = method, shape = method)) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.5, na.rm = TRUE, show.legend = FALSE,
                              position = ggplot2::position_dodge(width = 0.6)) +
      ggplot2::geom_point(size = 1.6, position = ggplot2::position_dodge(width = 0.6), na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = stats::setNames(eco_categorical(3, opts$palette %||% "okabe_ito"), mth)) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17, 15), mth)) +
      ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::labs(x = L("hr_area", sprintf("%.0f%% home-range area (km²)", 100 * H$level)), y = NULL, colour = NULL, shape = NULL)
    attr(p, "fig_size") <- list(height = 30 + 6 * length(unique(hd$id)))
    out$hr_area <- with_caption(p, sprintf("%.0f%% home-range area per %s: AKDE (circles) with 95%% confidence interval (ctmm), and KDE with reference bandwidth (triangles) and MCP (squares), which assume independent fixes and typically underestimate area for autocorrelated GPS data. KDE and MCP have no interval.", 100 * H$level, an))
  }
  out
}

ssf_plots <- function(fit, data) {
  f <- ssf_figures(fit, data, fig_opts())
  f[intersect(c(paste0("rss_", data$covariates), "selection_coefficients", "movement_kernel", "used_available", "uhc", "home_ranges", "hr_area"), names(f))]
}

# Write simulated data in the researcher format (GPS CSV + one GeoTIFF per covariate).
ssf_write_data <- function(data, outdir, stem = "gps_fixes") {
  eco_require("terra")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fx <- data$fixes
  utils::write.csv(data.frame(animal_id = fx$id, timestamp = format(fx$t, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                              utm_e = round(fx$x, 1), utm_n = round(fx$y, 1)),
                   file.path(outdir, paste0(stem, ".csv")), row.names = FALSE)
  L <- data$land
  for (v in names(L$layers)) {
    r <- terra::rast(L$layers[[v]], extent = terra::ext(L$xmin, L$xmin + L$ncol * L$resx, L$ymax - L$nrow * L$resy, L$ymax), crs = L$crs)
    terra::writeRaster(r, file.path(outdir, paste0(v, ".tif")), overwrite = TRUE)
  }
  invisible(outdir)
}
