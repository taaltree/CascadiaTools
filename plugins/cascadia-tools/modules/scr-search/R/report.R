# scr-search: design notes, interpretation, figures, methods text, references.

ss_references <- c(
  "Royle JA, Kery M, Guelat J (2011) Spatial capture-recapture models for search-encounter data. Methods in Ecology and Evolution 2:602-611. doi:10.1111/j.2041-210X.2011.00116.x",
  "Thompson CM, Royle JA, Garner JD (2012) A framework for inference about carnivore density from unstructured spatial sampling of scat using detector dogs. Journal of Wildlife Management 76:863-871. doi:10.1002/jwmg.317",
  "Russell RE, Royle JA, Desimone R, Schwartz MK, Edwards VL, Pilgrim KP, McKelvey KS (2012) Estimating abundance of mountain lions from unstructured spatial sampling. Journal of Wildlife Management 76:1551-1561. doi:10.1002/jwmg.412",
  "Efford MG (2011) Estimation of population density by spatially explicit capture-recapture analysis of data from area searches. Ecology 92:2202-2207. doi:10.1890/11-0332.1",
  "Efford MG (2025) Spatially explicit capture-recapture in R (SECRbook), Appendix D: Area and transect search; and secr: Spatially explicit capture-recapture models, R package version 5.x. https://murrayefford.github.io/SECRbook/",
  "Milleret C, Dupont P, Bonenfant C, Broseth H, Flagstad O, Sutherland C, Bischof R (2019) A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363. doi:10.1002/ece3.4751",
  "Moqanaki EM, Milleret C, Tourani M, Dupont P, Bischof R (2021) Consequences of ignoring variable and spatially autocorrelated detection probability in spatial capture-recapture. Landscape Ecology 36:2879-2895. doi:10.1007/s10980-021-01283-x",
  "Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press, Waltham MA. (Ch. 15: search-encounter models)",
  "Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(Suppl 2):S521-S537. doi:10.1007/s10336-010-0619-4",
  "Borchers DL, Efford MG (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64:377-385. doi:10.1111/j.1541-0420.2007.00927.x",
  "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487")

# ---------------------------------------------------------------- draws
# Common parameter draws for every engine (internal units) with attribute M.
# secr: pseudo-draws from the multivariate normal of the link-scale MLEs,
# expressed as psi = E[N] / M with large M.
ss_pdraws <- function(fit, ndraw = 1000) {
  ob <- fit$object
  if (!is.null(ob$chains)) { d <- ss_pool(ob$chains, ndraw); attr(d, "M") <- ob$prep$M_used; return(d) }
  prep <- ob$prep; sf <- ob$secr
  b <- sf$fit$par; V <- sf$beta.vcv
  rn <- rownames(stats::coef(sf))
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  B <- sweep(matrix(stats::rnorm(ndraw * length(b)), ndraw) %*% t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(b))), 2, b, "+")
  colnames(B) <- rn
  PD <- length(prep$dcov); PW <- length(prep$wcov)
  Mbig <- 1000 * max(prep$n, 1)
  cell_ha <- (prep$pix * 1000)^2 / 1e4
  out <- matrix(NA_real_, ndraw, 3 + PD + PW)
  colnames(out) <- c("psi", "lam0", "sigma", if (PD) paste0("bD", seq_len(PD)), if (PW) paste0("bW", seq_len(PW)))
  for (r in seq_len(ndraw)) {
    bD <- if (PD) B[r, paste0("D.", prep$dcov)] else numeric(0)
    lw <- if (PD) as.vector(prep$XD %*% bD) else rep(0, prep$G)
    out[r, "psi"] <- min(sum(exp(B[r, "D"] + lw)) * cell_ha / Mbig, 0.999)
    sg <- exp(B[r, "sigma"]) / 1000
    out[r, "sigma"] <- sg
    out[r, "lam0"] <- if (ob$engine == "secr_transect") exp(B[r, "lambda0"]) / (sqrt(2 * pi) * sg) else exp(B[r, "lambda0"])
    if (PD) out[r, paste0("bD", seq_len(PD))] <- bD
    if (PW) out[r, paste0("bW", seq_len(PW))] <- B[r, paste0("lambda0.", prep$wcov)]
  }
  attr(out, "M") <- Mbig
  out
}

ss_density_at <- function(prep, pd, xs) {
  M <- attr(pd, "M"); PD <- ncol(prep$XD)
  vapply(seq_len(nrow(pd)), function(r) {
    b <- if (PD) pd[r, paste0("bD", seq_len(PD))] else numeric(0)
    lw <- if (PD) as.vector(prep$XD %*% b) else rep(0, prep$G)
    pd[r, "psi"] * M * exp(sum(b * xs) - ss_lse(lw)) / prep$pix^2 * ss_area_conv(prep)
  }, 1)
}

ss_q <- function(x, p = c(0.025, 0.975)) stats::quantile(x, p, na.rm = TRUE, names = FALSE)
ss_lab <- function(v, data, opts = NULL) {
  d <- if (!is.null(data$covariate_labels[[v]])) data$covariate_labels[[v]] else v
  fig_label(v, d, opts)
}
ss_fmt <- function(x, d = 2) sub("\\.$", "", trimws(formatC(x, digits = d, format = "fg", flag = "#")))

# Descriptive, relabellable names for reported parameters (figures)
ss_param_label <- function(p, prep, data = NULL, opts = NULL) {
  al <- ss_area_label(prep$units); ll <- ss_len_label(prep$units); ul <- ss_unit_word(prep$units)
  vapply(p, function(k) {
    d <- switch(k,
      density_per_km2 = , density_per_unit2 = sprintf("Density (per %s)", al),
      abundance = "Abundance (N)",
      lambda0_per_km = , lambda0_per_unit = sprintf("lambda0 (samples per %s)", ll),
      sigma = sprintf("sigma (%s)", ul),
      home_range_area_95 = sprintf("95%% home range (%s)", al),
      samples_per_track_through_centre = "Samples, track via centre",
      if (startsWith(k, "density_beta_")) sprintf("Density: %s", ss_lab(sub("^density_beta_", "", k), data))
      else if (startsWith(k, "detection_beta_")) sprintf("Detection: %s", ss_lab(sub("^detection_beta_", "", k), data))
      else k)
    fig_label(k, d, opts)
  }, "")
}
# map axes in km when coordinates are metres
ss_map_axes <- function(prep, opts) {
  km <- identical(prep$units, "m")
  f <- if (km) function(b) format(b / 1000, trim = TRUE, scientific = FALSE) else ggplot2::waiver()
  list(ggplot2::scale_x_continuous(labels = f, breaks = scales::breaks_pretty(4), expand = c(0, 0)),
       ggplot2::scale_y_continuous(labels = f, breaks = scales::breaks_pretty(4), expand = c(0, 0)),
       ggplot2::labs(x = fig_label("easting", if (km) "Easting (km)" else sprintf("Easting (%s)", ss_unit_word(prep$units)), opts),
                     y = fig_label("northing", if (km) "Northing (km)" else sprintf("Northing (%s)", ss_unit_word(prep$units)), opts)))
}
ss_bin_levels <- function(v, data, opts) c(fig_label(paste0(v, "_0"), "No", opts), fig_label(paste0(v, "_1"), "Yes", opts))

# ---------------------------------------------------------------- design notes
ss_design_notes <- function(prep, est) {
  g <- function(p) est$estimate[est$param == p]
  sig <- g("sigma"); out <- character()
  if (!length(sig) || !is.finite(sig)) return(out)
  ul <- ss_unit_word(prep$units)
  out <- c(out, sprintf("Search effort: %s %s of track over %d occasions in %d detector cells of %s %s (%d cell x track units); %d genotyped samples from %d individuals.",
                        ss_fmt(prep$total_effort / (if (identical(prep$units, "m")) 1000 else 1), 3), ss_len_label(prep$units),
                        prep$K, prep$C, format(prep$cell_size), ul, prep$U, prep$n_samples, prep$n))
  r <- prep$cell_size / sig
  out <- c(out, sprintf("Detector cell side is %.2f x sigma (recommended <= 1; coarser cells blur where samples were found and bias sigma upward).", r))
  if (r > 1) out <- c(out, "WARNING: detector cells are larger than sigma; use smaller cells (settings/data cell_size).")
  if (is.finite(prep$buffer)) {
    br <- prep$buffer / sig
    out <- c(out, sprintf("State-space buffer is %.1f x sigma beyond the searched area.", br))
    if (br < 3) out <- c(out, "WARNING: the buffer is less than 3 sigma; animals living outside the state space can be sampled, inflating density. Use a buffer of at least 3-4 sigma.")
  }
  out <- c(out, sprintf("%d of %d identified individuals were sampled in more than one detector cell (spatial recaptures).", prep$n_spatial_recap, prep$n))
  if (prep$n_spatial_recap < 5) out <- c(out, "WARNING: fewer than 5 spatial recaptures: sigma (and therefore density) is weakly identified.")
  if (prep$n < 15) out <- c(out, sprintf("WARNING: only %d individuals identified; estimates will be imprecise and may be biased.", prep$n))
  if (prep$effort_model != "track_length") out <- c(out, sprintf("NOTE: effort model '%s' is a simplification; the default 'track_length' (km searched per cell) is the validated model.", prep$effort_model))
  if (is.finite(prep$R)) out <- c(out, sprintf("Local evaluation: pixel-cell pairs farther apart than %s %s (%.1f x sigma) were ignored (contribution < exp(-%.0f)).",
                                               format(signif(prep$R * prep$cscale, 3)), ul, prep$R * prep$cscale / sig, (prep$R * prep$cscale / sig)^2 / 2))
  out
}

# ---------------------------------------------------------------- interpret
ss_interpret <- function(fit, data) {
  est <- fit$estimates; prep <- fit$object$prep
  row <- function(p) est[est$param == p, ]
  dn <- ss_density_name(prep$units); al <- ss_area_label(prep$units); ll <- ss_len_label(prep$units); ul <- ss_unit_word(prep$units)
  D <- row(dn); N <- row("abundance"); lam <- row(ss_lam_name(prep$units)); sg <- row("sigma"); hr <- row("home_range_area_95")
  tc <- row("samples_per_track_through_centre")
  isml <- grepl("^secr", fit$object$engine)
  ci <- function(r, d = 2) sprintf("%s (95%% %s %s to %s)", ss_fmt(r$estimate, d), if (isml) "CI" else "CrI", ss_fmt(r$lower, d), ss_fmt(r$upper, d))
  A_rep <- prep$area * ss_hr_conv(prep)
  eff <- prep$total_effort / (if (identical(prep$units, "m")) 1000 else 1)
  out <- c(
    sprintf("Data: %s %s of search tracks over %d occasions yielded %d genotyped samples from %d individuals; %d individuals were sampled in more than one detector cell, which is what lets the model estimate how far animals range.",
            ss_fmt(eff, 3), ll, prep$K, prep$n_samples, prep$n, prep$n_spatial_recap),
    sprintf("Density was estimated at %s animals per %s%s. Over the %s %s state space (the searched area plus a buffer), that is %s animals, of which %d (%.0f%%) were identified from samples.",
            ci(D, 2), al, if (dn == "density_per_km2") sprintf(", i.e. %s per 100 km2", ss_fmt(100 * D$estimate, 3)) else "",
            ss_fmt(A_rep, 3), al, ci(N, 3), prep$n, 100 * prep$n / N$estimate),
    sprintf("The spatial scale of space use, sigma, was %s %s. If space use is a circular bivariate normal around each activity centre, 95%% of an animal's activity falls within %s %s of its centre, a home-range area of about %s %s; treat this as an index of ranging scale.",
            ci(sg, 3), ul, ss_fmt(sqrt(stats::qchisq(0.95, 2)) * sg$estimate, 3), ul, ci(hr, 2), al),
    sprintf("Encounter rate: a searcher covering 1 %s right at an animal's activity centre finds on average %s of its samples%s; a long straight search route passing through the activity centre would yield %s samples of that animal. The rate falls off with distance as a half-normal curve (see the detection-function figure), so the expected number of samples from an animal is proportional to the search effort near its activity centre.",
            ll, ci(lam, 2), if (length(prep$wcov)) " (at the reference level of the detection covariates)" else "", ci(tc, 2))
  )
  pd <- tryCatch(ss_pdraws(fit, 1000), error = function(e) NULL)
  for (k in seq_along(prep$wcov)) {
    v <- prep$wcov[k]; r <- row(paste0("detection_beta_", v)); sc <- prep$wscale[[v]]; lab <- ss_lab(v, data)
    rr <- exp(c(r$estimate, r$lower, r$upper))
    out <- c(out, sprintf("Detection: %s multiplied the number of samples found per %s searched by %s (95%% interval %s to %s)%s. This is a detection (search efficiency) effect, not a difference in where animals live.",
                          if (isTRUE(sc$binary)) sprintf("'%s' (1 vs 0)", lab) else sprintf("each 1 SD (%s) increase in '%s'", ss_fmt(sc$scale), lab),
                          ll, ss_fmt(rr[1], 3), ss_fmt(rr[2], 3), ss_fmt(rr[3], 3),
                          if (!is.null(pd) && isTRUE(sc$binary)) sprintf("; per-km encounter rate at the activity centre %s with %s = 1 versus %s with %s = 0",
                                                                         ss_fmt(mean(pd[, "lam0"] * exp(pd[, paste0("bW", k)])) * ss_len_conv(prep)), lab, ss_fmt(mean(pd[, "lam0"]) * ss_len_conv(prep)), lab) else ""))
  }
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; r <- row(paste0("density_beta_", v)); sc <- prep$dscale[[v]]; lab <- ss_lab(v, data)
    rr <- exp(c(r$estimate, r$lower, r$upper))
    s1 <- sprintf("Density and %s: each 1 SD (%s units) higher %s is associated with %s%% %s density (rate ratio %s, 95%% interval %s to %s).",
                  lab, ss_fmt(sc$scale), lab, ss_fmt(abs(100 * (rr[1] - 1)), 2), if (rr[1] >= 1) "higher" else "lower",
                  ss_fmt(rr[1], 3), ss_fmt(rr[2], 3), ss_fmt(rr[3], 3))
    if (!is.null(pd)) {
      xr <- stats::quantile(prep$mask[[v]], c(0.1, 0.9))
      xs <- function(x) { z <- rep(0, length(prep$dcov)); z[k] <- (x - sc$center) / sc$scale; z }
      lo <- ss_density_at(prep, pd, xs(xr[1])); hi <- ss_density_at(prep, pd, xs(xr[2]))
      s1 <- paste(s1, sprintf("Expected density is %s per %s (%s to %s) where %s = %s (10th percentile of the state space) and %s (%s to %s) where it = %s (90th percentile).",
                              ss_fmt(mean(lo)), al, ss_fmt(ss_q(lo)[1]), ss_fmt(ss_q(lo)[2]), lab, ss_fmt(xr[1], 3),
                              ss_fmt(mean(hi)), ss_fmt(ss_q(hi)[1]), ss_fmt(ss_q(hi)[2]), ss_fmt(xr[2], 3)))
      if (!isml) s1 <- paste(s1, sprintf("Posterior probability that the association is positive: %.2f.", mean(pd[, paste0("bD", k)] > 0)))
    }
    s1 <- paste(s1, "This is a correlation between density and habitat across the state space, not proof that the habitat causes higher density. Because search effort was modelled explicitly (km searched per cell), more samples in heavily searched areas are not mistaken for higher density.")
    out <- c(out, s1)
  }
  out <- c(out,
           sprintf("Design: detector cells were %s %s (%.2f sigma) and the state space extends %.1f sigma beyond the searched area (recommended >= 3-4 sigma).",
                   format(prep$cell_size), ul, prep$cell_size / sg$estimate, prep$buffer / sg$estimate),
           "Assumptions: the population is closed during the survey, activity centres are fixed, samples of different individuals are found independently, the number of samples found per km of search declines with distance from the activity centre as a half-normal curve shared by all individuals, and every genotype is correctly assigned. Samples must come from the recorded searches: opportunistic samples found without recorded effort bias density upward where they cluster.",
           "Validation: at the module's default design (16 x 16 km searched area, ~250 km of dog-team tracks over 3 occasions, ~80 animals) the Nimble engine recovers density, abundance, lambda0, sigma and both covariate effects; when effort is concentrated along roads in good habitat, ignoring effort biases the habitat effect while modelling km searched per cell does not. With sparser data (posterior CV of density above about 0.2, fewer than about 35 animals identified) the posterior mean of density and N ran about 10% high in simulations (intervals kept nominal coverage); see the diagnostics note and modules/scr-search/VALIDATION.md.")
  out
}

# ---------------------------------------------------------------- figures
ss_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  ob <- fit$object; prep <- ob$prep
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; a <- (1 - ci) / 2
  ici <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci; ai <- (1 - ici) / 2
  al <- ss_area_label(prep$units); ul <- ss_unit_word(prep$units); ll <- ss_len_label(prep$units)
  isml <- grepl("^secr", ob$engine); band <- if (isml) "confidence (parametric bootstrap)" else "credible"
  pd <- ss_pdraws(fit, 1000)
  M <- attr(pd, "M")
  figs <- list()
  to_x <- function(v) v * prep$cscale + prep$origin[1]; to_y <- function(v) v * prep$cscale + prep$origin[2]
  xl <- fig_label("easting", sprintf("Easting (%s)", ul), opts); yl <- fig_label("northing", sprintf("Northing (%s)", ul), opts)
  trk <- prep$tracks
  smp <- prep$samples
  # ---- effort map
  ce <- tapply(exp(prep$logE) * prep$cscale / (if (identical(prep$units, "m")) 1000 else 1), factor(prep$ucell, levels = seq_len(prep$C)), sum)
  cc0 <- if (!is.null(prep$ccentre)) prep$ccentre else prep$cxy
  cg <- data.frame(x = to_x(cc0[, 1]), y = to_y(cc0[, 2]), value = as.numeric(ce))
  pe <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = cg, ggplot2::aes(x, y, fill = value), width = prep$cell_size, height = prep$cell_size, colour = "white", linewidth = 0.1) +
    scale_fill_eco_seq(option = "mako", direction = -1) +
    ggplot2::geom_path(data = trk, ggplot2::aes(x, y, group = track_id), colour = "grey20", linewidth = 0.15, alpha = 0.6) +
    ggplot2::coord_equal() + ss_map_axes(prep, opts) + ggplot2::labs(fill = sprintf("%s searched\nper cell", ll))
  if (isTRUE(opts$show_data) && nrow(smp))
    pe <- pe + ggplot2::geom_point(data = smp, ggplot2::aes(x, y), shape = 21, fill = "#E69F00", colour = "black", size = 0.9, stroke = 0.2)
  figs$effort_map <- with_caption(pe, sprintf(
    "Search effort: length of search track (%s) inside each %s %s detector cell (fill), GPS tracks (grey lines) and genotyped samples (orange points; %d samples of %d individuals). Unsearched ground has no detector.",
    ll, format(prep$cell_size), ul, prep$n_samples, prep$n))
  # ---- density surfaces
  surf <- ob$surfaces
  if (is.null(surf)) surf <- tryCatch(ss_surfaces(prep, pd[ss_thin_idx(nrow(pd), 60), , drop = FALSE], M), error = function(e) NULL)
  if (!is.null(surf)) {
    cellA <- prep$pix^2 / ss_area_conv(prep)
    grid <- data.frame(x = prep$mask$x, y = prep$mask$y, value = surf$realised / cellA)
    p <- fig_map(grid, sprintf("Animals per %s", al)) + ss_map_axes(prep, opts) +
      ggplot2::geom_path(data = trk, ggplot2::aes(x, y, group = track_id), inherit.aes = FALSE, colour = "white", linewidth = 0.12, alpha = 0.5)
    if (isTRUE(opts$show_data) && prep$n) {
      cen <- surf$Pc %*% prep$mxy
      p <- p + ggplot2::geom_point(data = data.frame(x = to_x(cen[, 1]), y = to_y(cen[, 2])), ggplot2::aes(x, y),
                                   inherit.aes = FALSE, shape = 4, colour = "#D55E00", size = 1.1, stroke = 0.5)
    }
    figs$density_surface <- with_caption(p, sprintf(
      "Posterior mean realised density (animals per %s): identified animals placed by their estimated activity centres plus the expected density of animals never sampled. White lines: search tracks; orange crosses: posterior mean activity centres of the %d identified individuals.", al, prep$n))
    if (length(prep$dcov)) {
      g2 <- data.frame(x = prep$mask$x, y = prep$mask$y, value = surf$expected / cellA)
      figs$density_expected_surface <- with_caption(
        fig_map(g2, sprintf("Expected animals per %s", al)) + ss_map_axes(prep, opts),
        sprintf("Expected density surface implied by the habitat model (%s).", paste(vapply(prep$dcov, ss_lab, "", data = data, opts = opts), collapse = ", ")))
    }
    g3 <- data.frame(x = prep$mask$x, y = prep$mask$y, value = surf$pdetect)
    figs$detection_map <- with_caption(
      fig_map(g3, "P(sampled at\nleast once)") + ss_map_axes(prep, opts) +
        ggplot2::geom_path(data = trk, ggplot2::aes(x, y, group = track_id), inherit.aes = FALSE, colour = "white", linewidth = 0.12, alpha = 0.5),
      "Posterior mean probability that an animal whose activity centre is in each pixel left at least one sample that was found, given the realised search effort. Low values: animals there are effectively unsampled, so density there rests on the habitat model and the buffer.")
  }
  # ---- conditional effects of density covariates
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; sc <- prep$dscale[[v]]
    xs <- seq(min(prep$mask[[v]]), max(prep$mask[[v]]), length.out = 60)
    Dm <- sapply(xs, function(x) { z <- rep(0, length(prep$dcov)); z[k] <- (x - sc$center) / sc$scale; ss_density_at(prep, pd, z) })
    pred <- data.frame(x = xs, estimate = colMeans(Dm), lower = apply(Dm, 2, stats::quantile, a), upper = apply(Dm, 2, stats::quantile, 1 - a),
                       lower50 = apply(Dm, 2, stats::quantile, ai), upper50 = apply(Dm, 2, stats::quantile, 1 - ai))
    obs <- if (isTRUE(opts$show_data)) data.frame(x = prep$mask[[v]][ss_thin_idx(prep$G, 300)]) else NULL
    figs[[paste0("density_", v)]] <- with_caption(
      fig_cond_effect(pred, ss_lab(v, data, opts), fig_label("density", sprintf("Density (animals per %s)", al), opts), obs = obs) +
        ggplot2::scale_y_continuous(limits = c(0, NA)),
      sprintf("Expected density versus %s (other density covariates at their state-space mean), with %d%% and %d%% %s bands. Rug: covariate values across the state space.",
              ss_lab(v, data, opts), round(100 * ci), round(100 * ici), band))
  }
  # ---- detection function: samples per km searched vs distance (by binary detection covariate)
  dmax <- 4 * stats::median(pd[, "sigma"]) * prep$cscale
  dd <- seq(0, dmax, length.out = 80)
  kmx <- identical(prep$units, "m"); dsc <- if (kmx) 1000 else 1
  bin <- prep$wcov[vapply(prep$wcov, function(v) isTRUE(prep$wscale[[v]]$binary), TRUE)]
  levs <- if (length(bin)) c(0, 1) else 0
  rows <- list()
  for (lv in levs) {
    l0 <- pd[, "lam0"] * ss_len_conv(prep) * if (length(bin)) exp(lv * pd[, paste0("bW", match(bin[1], prep$wcov))]) else 1
    P <- sapply(dd, function(x) l0 * exp(-(x / prep$cscale)^2 / (2 * pd[, "sigma"]^2)))
    r <- data.frame(x = dd / dsc, estimate = colMeans(P), lower = apply(P, 2, stats::quantile, a), upper = apply(P, 2, stats::quantile, 1 - a),
                    lower50 = apply(P, 2, stats::quantile, ai), upper50 = apply(P, 2, stats::quantile, 1 - ai))
    if (length(bin)) r$group <- sprintf("%s: %s", ss_lab(bin[1], data, opts), ss_bin_levels(bin[1], data, opts)[lv + 1])
    rows[[length(rows) + 1]] <- r
  }
  pdet <- fig_cond_effect(do.call(rbind, rows), fig_label("distance", sprintf("Distance from activity centre (%s)", if (kmx) "km" else ul), opts),
                          fig_label("encounter_rate", sprintf("Samples found per %s searched", ll), opts)) +
    ggplot2::geom_vline(xintercept = prep$cell_size / dsc, linetype = 2, colour = "grey45") +
    ggplot2::annotate("text", x = prep$cell_size / dsc, y = Inf, label = " cell size", hjust = 0, vjust = 1.3, size = 2.6, colour = "grey30") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
  if (length(bin)) {
    pdet <- pdet + ggplot2::guides(colour = ggplot2::guide_legend(position = "bottom", nrow = 2, title = NULL),
                                   fill = ggplot2::guide_legend(position = "bottom", nrow = 2, title = NULL))
    attr(pdet, "fig_size") <- list(width = "onehalf", aspect = 0.7)
  }
  figs$detection_function <- with_caption(pdet, sprintf(
    "Expected number of an animal's samples found per %s of search, versus the distance of the search from its activity centre: lambda0 exp(-d^2 / 2 sigma^2), with %d%% and %d%% %s bands%s. Dashed line: detector cell size.",
    ll, round(100 * ci), round(100 * ici), band, if (length(bin)) sprintf(", by %s", ss_lab(bin[1], data, opts)) else ""))
  # ---- detection covariates (response scale: lambda0 per km)
  for (k in seq_along(prep$wcov)) {
    v <- prep$wcov[k]; sc <- prep$wscale[[v]]
    if (isTRUE(sc$binary)) {
      L <- cbind(pd[, "lam0"], pd[, "lam0"] * exp(pd[, paste0("bW", k)])) * ss_len_conv(prep)
      pr <- data.frame(level = factor(ss_bin_levels(v, data, opts), levels = ss_bin_levels(v, data, opts)), estimate = colMeans(L), lower = apply(L, 2, stats::quantile, a), upper = apply(L, 2, stats::quantile, 1 - a))
      figs[[paste0("detection_", v)]] <- with_caption(
        fig_cat_effect(pr, ss_lab(v, data, opts), fig_label("lambda0", sprintf("Samples per %s at centre", ll), opts)),
        sprintf("Encounter rate lambda0 (samples found per %s searched at an animal's activity centre) without (No) and with (Yes) %s, with %d%% %s intervals.",
                ll, tolower(ss_lab(v, data, opts)), round(100 * ci), band))
    } else {
      wv <- unscale(prep$W[, k], sc$center, sc$scale)
      xs <- seq(min(wv), max(wv), length.out = 60)
      L <- sapply(xs, function(x) pd[, "lam0"] * exp(pd[, paste0("bW", k)] * (x - sc$center) / sc$scale) * ss_len_conv(prep))
      pred <- data.frame(x = xs, estimate = colMeans(L), lower = apply(L, 2, stats::quantile, a), upper = apply(L, 2, stats::quantile, 1 - a),
                         lower50 = apply(L, 2, stats::quantile, ai), upper50 = apply(L, 2, stats::quantile, 1 - ai))
      figs[[paste0("detection_", v)]] <- with_caption(
        fig_cond_effect(pred, ss_lab(v, data, opts), fig_label("lambda0", sprintf("Samples per %s at activity centre", ll), opts),
                        obs = if (isTRUE(opts$show_data)) data.frame(x = wv) else NULL) + ggplot2::scale_y_continuous(limits = c(0, NA)),
        sprintf("Encounter rate lambda0 versus %s (other detection covariates at reference/mean), %d%% and %d%% %s bands; rug: values over searched cells.",
                ss_lab(v, data, opts), round(100 * ci), round(100 * ici), band))
    }
  }
  # ---- abundance
  if (!is.null(ob$rep_chains)) {
    Nd <- unlist(lapply(ob$rep_chains, function(m) m[, "abundance"]))
    qq <- ss_q(Nd, c(a, 1 - a))
    figs$abundance_posterior <- with_caption(
      ggplot2::ggplot(data.frame(N = Nd), ggplot2::aes(N)) +
        ggplot2::geom_histogram(binwidth = max(1, round(diff(range(Nd)) / 40)), fill = eco_pal$band, colour = "white", alpha = 0.7) +
        ggplot2::geom_vline(xintercept = prep$n, linetype = 3, colour = "grey30") +
        ggplot2::geom_vline(xintercept = qq, linetype = 2, colour = eco_pal$ink) +
        ggplot2::labs(x = fig_label("abundance", "Abundance in the state space (N)", opts), y = "Posterior draws"),
      sprintf("Posterior distribution of the number of activity centres in the %s %s state space; dashed lines: %d%% credible interval; dotted line: number of individuals identified (%d).",
              format(round(prep$area * ss_hr_conv(prep))), al, round(100 * ci), prep$n))
  } else {
    r <- fit$estimates[fit$estimates$param == "abundance", ]
    figs$abundance_posterior <- with_caption(
      ggplot2::ggplot(r, ggplot2::aes(estimate, "N")) + ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper)) +
        ggplot2::geom_vline(xintercept = prep$n, linetype = 3) + ggplot2::labs(x = "Abundance in the state space (N)", y = NULL),
      "Realised abundance in the state space (secr region.N) with 95% confidence interval; dotted line: individuals identified.")
  }
  # ---- samples per individual and GOF
  ppc <- ob$ppc
  if (!is.null(ppc) && prep$n) {
    obsn <- table(factor(pmin(ppc$obs_caps, 10), levels = 1:10))
    repn <- t(vapply(ppc$caps_rep, function(v) as.numeric(table(factor(pmin(v, 10), levels = 1:10))), numeric(10)))
    dfc <- data.frame(k = factor(c(1:9, "10+"), levels = c(1:9, "10+")), obs = as.numeric(obsn),
                      est = colMeans(repn), lo = apply(repn, 2, stats::quantile, a), hi = apply(repn, 2, stats::quantile, 1 - a))
    figs$captures_per_individual <- with_caption(
      ggplot2::ggplot(dfc, ggplot2::aes(k)) +
        ggplot2::geom_col(ggplot2::aes(y = obs), fill = "grey75", width = 0.7) +
        ggplot2::geom_pointrange(ggplot2::aes(y = est, ymin = lo, ymax = hi), colour = eco_pal$ink, size = 0.3) +
        ggplot2::labs(x = if (prep$observation == "bernoulli") "Search units with samples per individual" else "Samples per identified individual", y = "Number of individuals"),
      sprintf("Bars: observed number of individuals with 1, 2, ... samples. Points and lines: posterior predictive mean and %d%% interval from data simulated under the fitted model with the realised search effort. Bars far outside the intervals indicate unmodelled heterogeneity (e.g. sex differences, latrines producing clusters of scats).", round(100 * ci)))
    nm <- c(n_detected = "Individuals identified", total_samples = "Total samples", cells_per_individual = "Cells per individual",
            disp_samples = "Samples/individual: var/mean")
    pl <- lapply(names(nm), function(s) fig_ppc(ppc$rep[, s], stats::median(ppc$obs[, s], na.rm = TRUE), nm[[s]]))
    ftd <- data.frame(o = ppc$obs[, "ft_cell_totals"], r = ppc$rep[, "ft_cell_totals"])
    pl[[5]] <- ggplot2::ggplot(ftd, ggplot2::aes(o, r)) + ggplot2::geom_point(alpha = 0.4, size = 0.8, colour = eco_pal$ink) +
      ggplot2::geom_abline(linetype = 2, colour = "grey45") +
      ggplot2::labs(x = "FT statistic, observed", y = "FT, replicated",
                    caption = sprintf("P(replicate >= observed) = %.2f", ppc$p[["ft_cell_totals"]]))
    pl <- lapply(pl, function(p) p + theme_eco(opts$base_size, opts$font, opts$style, opts$legend))
    figs$gof <- with_caption(compose_figures(pl, ncol = 2),
      "Posterior predictive checks: histograms of statistics from data simulated under the fitted model, line = observed value; last panel: Freeman-Tukey discrepancy of samples per detector cell, observed vs replicated (points above the 1:1 line = data consistent with the model).")
    attr(figs$gof, "fig_size") <- list(width = "double", aspect = 1)
  }
  # ---- coefficients
  ec <- fit$estimates[grepl("^(density|detection)_beta_", fit$estimates$param), ]
  if (nrow(ec)) {
    ec$param <- ss_param_label(ec$param, prep, data, opts)
    pc <- fig_coef(ec, null_value = 0, xlab = "Log rate ratio (95% interval)")
    attr(pc, "fig_size") <- c(attr(pc, "fig_size"), list(width = "onehalf"))
    figs$coefficients <- with_caption(pc,
      "Log-scale coefficients: density covariates per 1 SD of the covariate; detection (encounter-rate) covariates per 1 SD, or level 1 vs 0 for 0/1 covariates. exp(coefficient) is the rate ratio. Dashed line: no effect.")
  }
  # ---- engine comparison
  if (!is.null(ob$comparison)) {
    cm <- ob$comparison
    cm$engine <- factor(cm$engine, levels = unique(cm$engine))
    cm$param <- factor(ss_param_label(as.character(cm$param), prep, data, opts), levels = unique(ss_param_label(as.character(cm$param), prep, data, opts)))
    figs$engine_comparison <- with_caption(
      ggplot2::ggplot(cm, ggplot2::aes(estimate, engine, colour = engine, shape = engine)) +
        ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.25) +
        ggplot2::facet_wrap(~param, scales = "free_x", ncol = 3, labeller = ggplot2::label_wrap_gen(22)) +
        ggplot2::scale_x_continuous(breaks = function(l) { r <- diff(l); b <- pretty(l, 3); b[b >= l[1] + 0.1 * r & b <= l[2] - 0.1 * r] },
                                    expand = ggplot2::expansion(mult = 0.12)) +
        ggplot2::guides(colour = "none", shape = "none") +
        ggplot2::labs(x = "Estimate (95% interval)", y = NULL) +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(size = 6),
                       panel.spacing.x = ggplot2::unit(4, "mm")),
      "The same data fitted with each engine: Bayesian marginal likelihood on detector cells (Nimble, local evaluation), maximum likelihood of the same cell model (secr_cells) and maximum likelihood with the tracks as continuous transect detectors and exact sample locations (secr_transect). Agreement indicates the result depends neither on the software nor on discretising the tracks into cells.")
    attr(figs$engine_comparison, "fig_size") <- list(width = "double", aspect = 0.8)
  }
  if (!is.null(ob$rep_chains)) {
    rc <- lapply(ob$rep_chains, function(m) { colnames(m) <- vapply(ss_param_label(colnames(m), prep, data, opts), function(z) if (nchar(z) > 26) sub(": ", ":\n", z) else z, ""); m })
    figs$trace <- with_caption(fig_trace(rc), "MCMC traces of the reported quantities, one colour per chain (thinned).")
    attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.8)
    figs$posterior <- with_caption(fig_posterior(do.call(rbind, rc)), "Posterior densities of the reported quantities.")
    attr(figs$posterior, "fig_size") <- list(width = "double", aspect = 0.8)
  }
  figs
}

# ---------------------------------------------------------------- methods
ss_methods_text <- function(fit, data, settings) {
  ob <- fit$object; prep <- ob$prep; st <- ob$settings
  al <- ss_area_label(prep$units); ul <- ss_unit_word(prep$units); ll <- ss_len_label(prep$units)
  dens <- if (length(prep$dcov)) sprintf("log-linear in %s (standardised to mean 0, SD 1 over the state space)", paste(prep$dcov, collapse = ", ")) else "homogeneous"
  detc <- if (length(prep$wcov)) sprintf(" multiplied by exp(alpha' w_u) for search covariate(s) %s", paste(prep$wcov, collapse = ", ")) else ""
  obs <- if (prep$observation == "poisson") "the number of genotyped samples of individual i in unit u was Poisson with mean lambda_iu"
         else "whether individual i was sampled at least once in unit u was Bernoulli with probability 1 - exp(-lambda_iu)"
  base <- sprintf(paste(
    "We estimated population density with a search-encounter spatial capture-recapture model (Royle et al. 2011; Thompson et al. 2012; Royle et al. 2014).",
    "GPS tracks of the searches were discretised into %s %s square detector cells (Milleret et al. 2019; Moqanaki et al. 2021); the effort E_u of each search unit u (a track within a cell) was the length searched (%s km in total over %d occasions, %d cells).",
    "Each animal has a fixed activity centre in a state space represented by a habitat mask of %d pixels (%s %s spacing, %s %s, extending %s %s beyond the searched area); activity centres follow a Poisson point process whose intensity is %s.",
    "The expected number of samples of an animal with activity centre s found in unit u was lambda_iu = lambda0 E_u exp(-d(s, c_u)^2 / (2 sigma^2))%s, where c_u is %s, and %s.",
    "The 95%% home-range area was derived as pi (2.448 sigma)^2 and the expected samples from a long straight track through the activity centre as lambda0 sqrt(2 pi) sigma."),
    format(prep$cell_size), ul, ss_fmt(prep$total_effort / (if (identical(prep$units, "m")) 1000 else 1), 3), prep$K, prep$C,
    prep$G, format(signif(prep$spacing, 3)), ul, ss_fmt(prep$area * ss_hr_conv(prep), 3), al,
    format(signif(prep$buffer, 3)), ul, dens, detc,
    if (identical(prep$detector_location, "centroid")) "the length-weighted centroid of the track searched inside the cell" else "the cell centre", obs)
  eng <- if (ob$engine == "secr_transect") sprintf("The model was fitted by maximum likelihood in the R package secr (version %s; Efford 2011) with each search track as a transect detector used on its own occasion, exact sample locations, Poisson counts (binomN = 0) and a hazard half-normal detection function; realised abundance came from region.N.",
                                                   as.character(utils::packageVersion("secr")))
  else if (ob$engine == "secr_cells") sprintf("The model was fitted by maximum likelihood in the R package secr (version %s) with one count detector per cell (at the same location as in the Bayesian model), each search track treated as an occasion with usage equal to the length searched in the cell, track-level search covariates as time covariates, and a hazard half-normal detection function (binomN = 0).", as.character(utils::packageVersion("secr")))
  else sprintf(paste(
    "The model was fitted in a Bayesian framework with data augmentation (Royle & Dorazio 2012; M = %d, psi ~ Uniform(0,1)) and priors lambda0 ~ Uniform(0, %g per %s), sigma ~ Uniform(0, %s %s) and Normal(0, SD 2) for covariate coefficients.",
    "Activity centres and inclusion indicators were integrated out over the habitat mask (marginal likelihood), with local evaluation: pixel-cell pairs more than %s %s apart (%.1f sigma) were ignored and each identified animal's centre was evaluated only on pixels within that radius of all its sample cells (Milleret et al. 2019).",
    "The remaining parameters were sampled jointly with an adaptive block random-walk sampler in NIMBLE %s (de Valpine et al. 2017), %d chains of %d iterations after %d burn-in.",
    "Abundance was derived per draw as n + Binomial(M - n, P(z = 1 | never sampled)) and density as psi M / area.",
    "Convergence was assessed with rank-normalised R-hat and bulk effective sample size; model fit with posterior predictive checks of individuals identified, total samples, cells per individual, dispersion of samples per individual and a Freeman-Tukey statistic of samples per cell."),
    prep$M_used, ss_lam_max(prep) * ss_len_conv(prep), ll, format(signif(prep$sigma_max * prep$cscale, 3)), ul,
    format(signif(prep$R * prep$cscale, 3)), ul, prep$R * prep$cscale / fit$estimates$estimate[fit$estimates$param == "sigma"],
    as.character(utils::packageVersion("nimble")), st$nchains, st$niter - st$nburnin, st$nburnin)
  paste(base, eng)
}
