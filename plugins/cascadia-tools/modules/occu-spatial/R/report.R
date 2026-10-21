# occu-spatial: interpretation, methods text and the publication figure suite
# (docs/FIGURES.md row "occu-spatial").

est_row <- function(fit, p) fit$estimates[fit$estimates$param == p, , drop = FALSE]
fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
ci_txt <- function(r, d = 2, f = identity) sprintf("%s (95%% CrI %s to %s)", fmt(f(r$estimate), d), fmt(f(r$lower), d), fmt(f(r$upper), d))

term_label <- function(t, data, opts = NULL, sc = data$scaling[[t]]) {
  def <- if (!is.null(data$labels[[t]])) data$labels[[t]] else if (!is.null(sc$label)) sc$label else t
  fig_label(t, def, opts)
}
term_sd_txt <- function(t, sc) {
  if (is.null(sc) || isTRUE(sc$binary)) return("")
  if (isTRUE(all.equal(sc$scale, 1)) && isTRUE(all.equal(sc$center, 0))) return(" (1 SD)")
  sprintf(" (1 SD = %s in original units)", trimws(formatC(sc$scale, digits = 3, format = "fg")))
}
# a round step near 1 SD in original units, and the odds ratio for it
nice_step <- function(sdv) { b <- 10^floor(log10(sdv)); s <- b * c(1, 2, 5, 10); s[which.min(abs(log(s / sdv)))] }
step_txt <- function(draws, sc, what = "occupancy") {
  if (is.null(sc) || isTRUE(sc$binary) || (isTRUE(all.equal(sc$scale, 1)) && isTRUE(all.equal(sc$center, 0)))) return("")
  st <- nice_step(sc$scale); or <- exp(draws * st / sc$scale)
  # exp(posterior mean of the slope), the same convention as the per-SD odds ratio
  sprintf(" Per %s-unit increase the odds of %s are multiplied by %s (95%% CrI %s to %s).", format(st), what,
          fmt(exp(mean(draws) * st / sc$scale)), fmt(stats::quantile(or, 0.025)), fmt(stats::quantile(or, 0.975)))
}

# posterior of a linear predictor along one standardised term, others at 0
curve_draws <- function(B, term, xs) {
  b0 <- B[, "(Intercept)"]; b1 <- B[, term]
  stats::plogis(outer(b0, rep(1, length(xs))) + outer(b1, xs))
}
band <- function(P, xs, ci = 0.95, inner = 0.5) {
  a <- (1 - ci) / 2; ai <- (1 - inner) / 2
  data.frame(x = xs, estimate = colMeans(P),
             lower = apply(P, 2, stats::quantile, a), upper = apply(P, 2, stats::quantile, 1 - a),
             lower50 = apply(P, 2, stats::quantile, ai), upper50 = apply(P, 2, stats::quantile, 1 - ai))
}

interpret_occu_spatial <- function(fit, data) {
  o <- fit$object; out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the MCMC did not converge (R-hat >= 1.05 for at least one estimand). Do not interpret the numbers below; rerun with more iterations (settings: n_iter, n_burn).")
  model_txt <- switch(o$model,
    spPGOcc = "a single-species spatial occupancy model (spatially structured site effect, exponential covariance, nearest-neighbour Gaussian process)",
    PGOcc = "a single-species occupancy model without a spatial effect",
    msPGOcc = "a multi-species (community) occupancy model without spatial effects",
    sfMsPGOcc = sprintf("a multi-species spatial factor occupancy model (%d spatial factor%s shared among species)", o$settings$n_factors %||% 1, if ((o$settings$n_factors %||% 1) > 1) "s" else ""),
    intPGOcc = "an integrated occupancy model that combines several survey methods with a shared occupancy process",
    spIntPGOcc = "a spatial integrated occupancy model")
  out <- c(out, paste0("Model: ", model_txt, ". All effects are associations from observational data, not causal effects. Covariate slopes are per 1 standard deviation (SD) of the covariate across sites (or surveys)."))
  if (o$type %in% c("single", "integrated")) {
    r0 <- est_row(fit, "occupancy_prob"); mo <- est_row(fit, "mean_occupancy"); pso <- est_row(fit, "prop_sites_occupied")
    naive <- mean(o$sites$detected)
    out <- c(out, sprintf("Occupancy probability at average covariate values%s is %s. Averaged over the %d surveyed sites, mean occupancy probability is %s, and the estimated proportion of surveyed sites that are occupied is %s, compared with %.2f of sites where the species was detected at least once (the difference is the correction for missed detections).",
                          if (is_spatial_model(o$model)) " and zero spatial effect" else "", ci_txt(r0), nrow(data$coords), ci_txt(mo), ci_txt(pso), naive))
    for (t in data$psi_terms) {
      r <- est_row(fit, paste0("psi_beta_", t)); if (!nrow(r)) next
      P <- curve_draws(o$beta, t, c(-1, 1)); lo <- colMeans(P)
      pp <- mean(o$draws[, paste0("psi_beta_", t)] > 0)
      out <- c(out, sprintf("Occupancy and %s: each 1-SD increase%s multiplies the odds of occupancy by %s; predicted occupancy goes from %.2f at 1 SD below the mean to %.2f at 1 SD above it (other covariates at their means).%s Posterior probability that the association is positive: %.2f.%s",
                            term_label(t, data), term_sd_txt(t, data$scaling[[t]]), ci_txt(r, 2, exp), lo[1], lo[2],
                            step_txt(o$draws[, paste0("psi_beta_", t)], data$scaling[[t]]), pp,
                            if (r$lower > 0 || r$upper < 0) "" else " The 95% interval includes no association."))
    }
    srcs <- if (o$type == "integrated") data$sources else ""
    for (s in srcs) {
      suf <- if (nzchar(s)) paste0("_", s) else ""
      rp <- est_row(fit, paste0("detection_prob", suf))
      out <- c(out, sprintf("Detection%s: when the species is present, a single survey%s detects it with probability %s (at average survey conditions).",
                            if (nzchar(s)) paste0(" by ", s) else "", if (nzchar(s)) paste0(" by ", s) else "", ci_txt(rp)))
      pt <- if (o$type == "integrated") data$p_terms[[s]] else data$p_terms
      psc <- if (o$type == "integrated") data$p_scaling[[s]] else data$p_scaling
      for (t in pt) {
        r <- est_row(fit, paste0("p_beta_", t, suf)); if (!nrow(r)) next
        out <- c(out, sprintf("Detection and %s%s: each 1-SD increase%s multiplies the odds of detection by %s.%s This is a detection effect (how easily the species is found), not a habitat effect.",
                              term_label(t, data, sc = psc[[t]]), if (nzchar(s)) paste0(" (", s, ")") else "",
                              term_sd_txt(t, psc[[t]]), ci_txt(r, 2, exp), step_txt(o$draws[, paste0("p_beta_", t, suf)], psc[[t]], "detection")))
      }
    }
    if (o$type == "integrated")
      out <- c(out, "Integrating the sources assumes they all observe the SAME occupancy state at a shared site (same species, same season, closure across methods); detection is estimated separately for each source.")
  } else {
    r0 <- est_row(fit, "comm_occupancy_prob"); sd0 <- est_row(fit, "comm_sd_occupancy_logit")
    N <- length(data$species)
    out <- c(out, sprintf("Across the %d species, the community-average occupancy probability at average covariates is %s; species vary around it with SD %s on the logit scale.", N, ci_txt(r0), ci_txt(sd0)))
    for (t in data$psi_terms) {
      r <- est_row(fit, paste0("comm_psi_beta_", t)); rs <- est_row(fit, paste0("comm_sd_psi_beta_", t))
      sp <- fit$estimates[fit$estimates$param %in% paste0("psi_beta_", t, "_", data$species), ]
      npos <- sum(sp$lower > 0); nneg <- sum(sp$upper < 0)
      out <- c(out, sprintf("Community response to %s: the average species' odds of occupancy are multiplied by %s per 1-SD increase%s (among-species SD of the slope %s on the logit scale). %d of %d species have a 95%% interval entirely above zero and %d entirely below zero; the rest are uncertain.",
                            term_label(t, data), ci_txt(r, 2, exp), term_sd_txt(t, data$scaling[[t]]), ci_txt(rs), npos, N, nneg))
    }
    rd <- est_row(fit, "comm_detection_prob")
    out <- c(out, sprintf("Detection: the community-average per-survey detection probability is %s.", ci_txt(rd)))
    rr <- est_row(fit, "mean_site_richness")
    out <- c(out, sprintf("Species richness: on average %s of the %d modelled species are estimated to occur at a surveyed site, versus %.1f detected per site (naive count). Richness here counts only species in the data set; species never detected anywhere are not included.",
                          ci_txt(rr, 1), N, mean(o$sites$naive_richness)))
  }
  if (is_spatial_model(o$model)) {
    rg <- fit$estimates[grepl("^effective_range_km", fit$estimates$param), ]
    sv <- est_row(fit, "spatial_variance")
    rng_txt <- paste(vapply(seq_len(nrow(rg)), function(i) sprintf("%s km (95%% CrI %s to %s km)", fmt(rg$estimate[i], 1), fmt(rg$lower[i], 1), fmt(rg$upper[i], 1)), ""), collapse = "; ")
    sv_txt <- ""
    if (nrow(sv)) {
      ig <- unlist(o$settings$sigma_sq_ig %||% c(2, 1))
      pq <- 1 / stats::qgamma(c(0.975, 0.025), shape = ig[1], rate = ig[2])
      near_prior <- (sv$upper - sv$lower) > 0.5 * diff(pq)
      sv_txt <- sprintf(" The residual spatial variance is %s on the logit scale; its prior alone gives a 95%% interval of %s to %s%s.",
                        ci_txt(sv), fmt(pq[1]), fmt(pq[2]),
                        if (near_prior) ", so the data add little information: these numbers mostly restate the prior and do not show whether residual clustering is strong or weak. In simulations without real residual clustering, the spatial model's occupancy slopes were about 10% larger than the non-spatial model's, so also fit the non-spatial model (settings: model: PGOcc) and report both if they differ" else "")
    }
    out <- c(out, sprintf("Residual spatial structure: the spatial random effect has an estimated effective range (distance at which the residual correlation drops to 0.05) of %s.%s The spatial variance and range are weakly identified by detection/non-detection data (simulation shows wide intervals and prior influence), so read them descriptively, not as a measured clustering or dispersal scale. Including the spatial term mainly guards the covariate intervals against being too narrow if residual clustering is present.",
                          rng_txt, sv_txt))
  }
  if (!is.null(o$gof)) {
    ps <- vapply(o$gof, `[[`, numeric(1), "p")
    bad <- sum(ps < 0.1 | ps > 0.9)
    out <- c(out, if (bad == 0) sprintf("Model check: posterior predictive p-values (Freeman-Tukey, by site) are %s; no sign of lack of fit (values near 0 or 1 would indicate it).", paste(fmt(ps), collapse = ", "))
                  else sprintf("Model check: %d of %d posterior predictive p-values are outside 0.1-0.9, suggesting lack of fit (e.g. unmodelled detection heterogeneity or missing covariates). Interpret with caution.", bad, length(ps)))
  }
  out <- c(out, paste(
    "Validated limits (VALIDATION.md): in simulation, covariate effects, baseline occupancy and detection were recovered within 10% with nominal 95% coverage at 200 sites x 4 surveys (single species), at 150 sites x 15 species (community means) and in the integrated design.",
    "With about 100-120 sites and per-survey detection of 0.25-0.3, posterior means of occupancy slopes and baseline occupancy were pulled away from zero by about 10-28%, while the 95% intervals still covered the truth about as often as they should. Treat effect sizes from small surveys as likely overstated.",
    "Between-species SDs in community models were overestimated by 15-30% with 10-15 species. The spatial variance and range are not reliably recovered and are reported descriptively only.",
    if (identical(o$model, "sfMsPGOcc")) "For the spatial factor model (sfMsPGOcc), 24 of 100 simulated data sets did not converge at short run lengths and species-level intercepts can mix slowly; check that every reported species-level estimate has R-hat < 1.01 and bulk ESS >= 400, and use msPGOcc if they do not." else NULL))
  out
}

methods_occu_spatial <- function(fit, data, settings) {
  o <- fit$object; st <- o$settings
  base <- "We modelled detection/non-detection data with hierarchical occupancy models that separate the occupancy state of each site from imperfect detection (MacKenzie et al. 2002), fitted by Polya-Gamma data-augmentation Gibbs samplers in the R package spOccupancy (Doser et al. 2022; Polson et al. 2013)."
  mod <- switch(o$model,
    spPGOcc = "Occupancy probability was modelled on the logit scale as a linear function of standardised site covariates plus a spatially structured random effect, a zero-mean Gaussian process with exponential covariance sigma^2 exp(-phi d), approximated by a nearest-neighbour Gaussian process (Datta et al. 2016) with %d neighbours. The effective range was computed as 3/phi.",
    PGOcc = "Occupancy probability was modelled on the logit scale as a linear function of standardised site covariates.",
    msPGOcc = "Species-specific occupancy and detection coefficients were treated as random effects drawn from community-level normal distributions (Dorazio & Royle 2005), giving community means and among-species standard deviations.",
    sfMsPGOcc = "Species-specific coefficients were drawn from community-level normal distributions, and residual spatial structure was modelled with a spatial factor model (Doser et al. 2023): %d latent spatial factors, each a nearest-neighbour Gaussian process with exponential correlation (%d neighbours), with species-specific loadings (upper triangle fixed at 0 and diagonal at 1 for identifiability).",
    intPGOcc = "Data from several survey methods were combined in an integrated occupancy model (Doser, Leuenberger et al. 2022), in which all methods observe a shared latent occupancy state and each method has its own detection model.",
    spIntPGOcc = "Data from several survey methods were combined in an integrated occupancy model with a spatial random effect (NNGP, exponential covariance).")
  mod <- if (o$model == "sfMsPGOcc") sprintf(mod, st$n_factors %||% 1, st$n_neighbors %||% 15) else if (grepl("%d", mod)) sprintf(mod, st$n_neighbors %||% 15) else mod
  pri <- paste0("Priors: Normal(0, variance 2.72) on occupancy and detection coefficients (and on community means)",
                if (o$type == "multi") "; inverse-gamma(0.1, 0.1) on among-species variances" else "",
                if (is_spatial_model(o$model)) "; uniform on phi, with the effective range bounded between the 5% quantile and the maximum of the inter-site distances" else "",
                if (o$model %in% c("spPGOcc", "spIntPGOcc")) "; inverse-gamma(2, 1) on the spatial variance" else "", ".")
  mc <- sprintf("We ran %d chains of %d iterations (burn-in %d, thinning %d) and checked convergence with the rank-normalised R-hat and bulk effective sample size (posterior package).",
                st$n_chains, st$n_iter, st$n_burn, st$n_thin)
  gof <- if (isTRUE(st$gof)) "Fit was assessed with a posterior predictive check using the Freeman-Tukey statistic on site-level detection counts." else ""
  paste(base, mod, "Numeric covariates were standardised to mean 0 and SD 1.", pri, mc, gof,
        "Analyses used the ecotools occu-spatial module, validated by simulation recovery (see VALIDATION.md).")
}

# ---- figures ---------------------------------------------------------------

site_map <- function(df, col, fill_lab, diverging = FALSE, opts = NULL, pred = NULL, pcol = col) {
  if (!is.null(pred) && pcol %in% names(pred)) {
    g <- data.frame(x = pred$x, y = pred$y, value = pred[[pcol]])
    p <- fig_map(g, fill_lab, points = df[, c("x", "y")], diverging = diverging)
  } else {
    d <- data.frame(x = df$x, y = df$y, value = df[[col]])
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y, colour = value)) + ggplot2::geom_point(size = 0.9) +
      ggplot2::coord_equal() + ggplot2::labs(colour = fill_lab)
    p <- p + if (diverging) ggplot2::scale_colour_gradient2(low = eco_pal$diverging[["low"]], mid = eco_pal$diverging[["mid"]], high = eco_pal$diverging[["high"]], midpoint = 0)
             else ggplot2::scale_colour_viridis_c()
  }
  # greyscale: sequential surfaces in one grey ramp (signed surfaces keep the diverging scale)
  # greyscale: sequential surfaces in one grey ramp; signed surfaces run dark (below 0) to
  # light (above 0) through mid grey at 0, so the sign stays readable without hue
  isg <- !is.null(pred) && pcol %in% names(pred)
  if (isTRUE(opts$greyscale) && !diverging)
    p <- p + if (isg) ggplot2::scale_fill_gradient(low = "grey92", high = "grey10")
             else ggplot2::scale_colour_gradient(low = "grey85", high = "grey10")
  if (isTRUE(opts$greyscale) && diverging)
    p <- p + if (isg) ggplot2::scale_fill_gradient2(low = "grey5", mid = "grey55", high = "grey97", midpoint = 0)
             else ggplot2::scale_colour_gradient2(low = "grey5", mid = "grey55", high = "grey97", midpoint = 0)
  p + ggplot2::scale_x_continuous(n.breaks = 3) + ggplot2::scale_y_continuous(n.breaks = 3) +
    ggplot2::labs(x = fig_label("easting", "Easting (km)", opts), y = fig_label("northing", "Northing (km)", opts))
}
# effective export width: per_figure width wins, then a non-default global width, else the module default
fig_eff_width <- function(name, default, opts) {
  g <- opts$width %||% "single"
  opts$per_figure[[name]]$width %||% (if (!identical(g, "single")) g else default)
}
map_pair <- function(a, b, opts = fig_opts(), name = "occupancy_map") {
  # The core themes only the last patchwork panel (with opts$legend), so give
  # panel A the same theme and legend position; both legends then sit alike.
  leg <- opts$legend %||% "right"
  w <- fig_eff_width(name, "double", opts); wmm <- fig_width_mm(w, opts)
  stack <- wmm < 120   # side by side needs ~1.5 columns; stack the panels at single-column width
  th <- theme_eco(opts$base_size, opts$font, opts$style, leg)
  # horizontal colourbars need room for their tick labels: title above, bar ~45% of the
  # panel width. Set on the guide itself, because the core re-themes the last panel.
  if (leg %in% c("bottom", "top")) {
    gb <- ggplot2::guide_colourbar(theme = ggplot2::theme(
      legend.key.width = ggplot2::unit(0.45 * wmm / if (stack) 1 else 2, "mm"),
      legend.key.height = ggplot2::unit(2.5, "mm"), legend.title.position = "top"))
    th <- list(th, ggplot2::guides(fill = gb, colour = gb))
  }
  p <- compose_figures(list(a + th, b + th), ncol = if (stack) 1 else 2, guides = "keep")
  # the core honours fig_size width/height hints (not aspect), so give a height
  # that matches the width actually used
  ratio <- if (stack) (if (leg %in% c("bottom", "top")) 2.1 else 1.6) else (if (leg %in% c("bottom", "top")) 0.6 else 0.44)
  attr(p, "fig_size") <- list(width = w, height = wmm * ratio)
  p
}
single_map <- function(p, opts = fig_opts(), name = "spatial_effect_map") {
  w <- fig_eff_width(name, "onehalf", opts)
  leg <- opts$legend %||% "right"
  attr(p, "fig_size") <- list(width = w, height = fig_width_mm(w, opts) * if (leg %in% c("bottom", "top")) 1.1 else 0.78); p }

# readable default species name ("long_tailed_weasel" -> "Long tailed weasel"); figures: labels: wins
sp_label <- function(sp, opts = NULL) {
  z <- gsub("_", " ", sp); fig_label(sp, paste0(toupper(substr(z, 1, 1)), substring(z, 2)), opts) }

# readable default label for an estimand (users can override via figures: labels:)
param_label <- function(p, data, opts = NULL) {
  lab <- function(t, sc) fig_label(t, if (!is.null(data$labels[[t]])) data$labels[[t]] else if (!is.null(sc[[t]]$label)) sc[[t]]$label else t, opts)
  if (grepl("^(comm_)?psi_beta_", p)) { t <- sub("^(comm_)?psi_beta_", "", p); return(paste0("Occupancy: ", lab(t, data$scaling))) }
  if (grepl("^(comm_)?p_beta_", p)) {
    t <- sub("^(comm_)?p_beta_", "", p)
    if (data$type == "integrated") for (s in data$sources) if (grepl(paste0("_", s, "$"), t)) {
      t2 <- sub(paste0("_", s, "$"), "", t); return(paste0("Detection (", s, "): ", lab(t2, data$p_scaling[[s]]))) }
    return(paste0("Detection: ", lab(t, data$p_scaling)))
  }
  known <- c(occupancy_prob = "Occupancy probability (mean covariates)",
             detection_prob = "Detection probability (mean covariates)",
             mean_occupancy = "Mean occupancy over sites",
             prop_sites_occupied = "Proportion of sites occupied",
             spatial_variance = "Spatial variance (logit scale)",
             spatial_decay_per_km = "Spatial decay phi (per km)",
             effective_range_km = "Effective range (km)")
  if (p %in% names(known)) return(known[[p]])
  z <- gsub("_", " ", p); paste0(toupper(substr(z, 1, 1)), substring(z, 2))
}

figures_occu_spatial <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2", "patchwork"))
  o <- fit$object; ci <- opts$ci %||% 0.95; inner <- opts$inner_ci %||% 0.5
  figs <- list()
  pct <- round(100 * ci)
  sites <- o$sites
  # occupancy / detection conditional effects
  if (o$type %in% c("single", "integrated")) {
    for (t in data$psi_terms) {
      xs <- seq(min(data$occ_X[[t]]), max(data$occ_X[[t]]), length.out = 100)
      b <- band(curve_draws(o$beta, t, xs), xs, ci, inner)
      sc <- data$scaling[[t]]; b$x <- unscale(b$x, sc$center, sc$scale)
      obs <- data.frame(x = unscale(data$occ_X[[t]], sc$center, sc$scale))
      figs[[paste0("effect_", t)]] <- with_caption(
        fig_cond_effect(b, term_label(t, data, opts), fig_label("occupancy", "Occupancy probability", opts),
                        obs = if (isTRUE(opts$show_data)) obs else NULL, ylim01 = TRUE),
        sprintf("%s vs %s (posterior mean, %d%% and %d%% credible bands), other covariates at their means%s; rug = surveyed sites.",
                fig_label("occupancy", "Occupancy probability", opts), term_label(t, data, opts), pct, round(100 * inner), if (is_spatial_model(o$model)) " and spatial effect at zero" else ""))
    }
    srcs <- if (o$type == "integrated") data$sources else ""
    A <- o$alpha; pos <- 0
    for (s in srcs) {
      pt <- if (o$type == "integrated") data$p_terms[[s]] else data$p_terms
      psc <- if (o$type == "integrated") data$p_scaling[[s]] else data$p_scaling
      dx <- if (o$type == "integrated") data$det_X[[s]] else data$det_X
      nt <- 1 + length(pt); As <- A[, pos + seq_len(nt), drop = FALSE]; pos <- pos + nt
      colnames(As) <- c("(Intercept)", pt)
      for (t in pt) {
        v <- dx[[t]][!is.na(dx[[t]])]; xs <- seq(min(v), max(v), length.out = 100)
        b <- band(curve_draws(As, t, xs), xs, ci, inner)
        sc <- psc[[t]] %||% list(center = 0, scale = 1, label = t)
        b$x <- unscale(b$x, sc$center, sc$scale)
        nm <- paste0("detection_", t, if (nzchar(s)) paste0("_", s) else "")
        figs[[nm]] <- with_caption(
          fig_cond_effect(b, fig_label(t, sc$label %||% t, opts), fig_label("detection", paste0("Detection probability per survey", if (nzchar(s)) paste0(" (", s, ")")), opts),
                          obs = if (isTRUE(opts$show_data)) data.frame(x = unscale(v, sc$center, sc$scale)) else NULL, ylim01 = TRUE),
          sprintf("Per-survey detection probability vs %s (posterior mean, %d%% and %d%% credible bands), other detection covariates at their means; rug = surveys.",
                  fig_label(t, sc$label %||% t, opts), pct, round(100 * inner)))
      }
    }
  } else {
    # community: community-mean curve with band + species posterior-mean curves
    for (t in data$psi_terms) {
      xs <- seq(min(data$occ_X[[t]]), max(data$occ_X[[t]]), length.out = 100)
      b <- band(curve_draws(o$beta, t, xs), xs, ci, inner)
      sc <- data$scaling[[t]]; xo <- unscale(xs, sc$center, sc$scale); b$x <- xo
      Bs <- o$beta_species
      spc <- do.call(rbind, lapply(data$species, function(sp) {
        b0 <- mean(Bs[, paste0("(Intercept)-", sp)]); b1 <- mean(Bs[, paste0(t, "-", sp)])
        data.frame(x = xo, estimate = stats::plogis(b0 + b1 * xs), species = sp_label(sp, opts))
      }))
      p <- ggplot2::ggplot() +
        ggplot2::geom_line(data = spc, ggplot2::aes(x, estimate, group = species), colour = eco_pal$muted, alpha = 0.45, linewidth = 0.35) +
        ggplot2::geom_ribbon(data = b, ggplot2::aes(x, ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.2) +
        ggplot2::geom_ribbon(data = b, ggplot2::aes(x, ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25) +
        ggplot2::geom_line(data = b, ggplot2::aes(x, estimate), colour = eco_pal$band, linewidth = 1) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = term_label(t, data, opts), y = fig_label("occupancy", "Occupancy probability", opts))
      if (isTRUE(opts$show_data)) p <- p + ggplot2::geom_rug(data = data.frame(x = unscale(data$occ_X[[t]], sc$center, sc$scale)),
                                                            ggplot2::aes(x = x), inherit.aes = FALSE, alpha = 0.3, sides = "b")
      figs[[paste0("effect_", t)]] <- with_caption(p, sprintf("Occupancy vs %s for the average species of the community (thick line, %d%% and %d%% credible bands) and each species (thin grey lines, posterior means); other covariates at their means.",
                                                              term_label(t, data, opts), pct, round(100 * inner)))
    }
    # community_effects forest plot
    for (t in data$psi_terms) {
      e <- fit$estimates
      sp <- e[e$param %in% paste0("psi_beta_", t, "_", data$species), ]
      sp$param <- vapply(sub(paste0("^psi_beta_", t, "_"), "", sp$param), function(s) sp_label(s, opts), "")
      cm <- e[e$param == paste0("comm_psi_beta_", t), ]; cm$param <- "Community mean"
      d <- rbind(cm, sp); d$group <- ifelse(d$param == "Community mean", "Community mean", "Species")
      d$param <- factor(d$param, levels = c(sp$param[order(sp$estimate)], "Community mean"))
      cml <- fig_label("community_mean", "Community mean", opts); spl <- fig_label("species", "Species", opts)
      d$param <- as.character(d$param); d$param[d$param == "Community mean"] <- cml
      d$param <- factor(d$param, levels = c(sp$param[order(sp$estimate)], cml))
      d$group <- factor(ifelse(d$group == "Community mean", cml, spl), levels = c(cml, spl))
      p <- ggplot2::ggplot(d, ggplot2::aes(estimate, param, colour = group, shape = group)) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = eco_pal$null) +
        ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper)) + ggplot2::geom_point(size = 1.6) +
        ggplot2::scale_colour_manual(values = stats::setNames(c(if (isTRUE(opts$greyscale)) "grey55" else eco_pal$truth, eco_pal$ink), c(cml, spl)), name = NULL) +
        ggplot2::scale_shape_manual(values = stats::setNames(c(18, 16), c(cml, spl)), name = NULL) +
        ggplot2::labs(x = fig_label(paste0("community_effects_", t, "_x"),
                                    sprintf("Effect of %s on occupancy\n(logit scale, per 1 SD; %d%% CrI)", term_label(t, data, opts), pct), opts), y = NULL)
      attr(p, "fig_size") <- list(height = 25 + 5 * nrow(d))
      figs[[if (length(data$psi_terms) > 1) paste0("community_effects_", t) else "community_effects"]] <-
        with_caption(p, sprintf("Species-specific effects of %s on occupancy (circles) and the community mean (diamond), posterior means with %d%% credible intervals; dashed line = no effect.", term_label(t, data, opts), pct))
    }
    figs$richness_map <- with_caption(
      map_pair(site_map(sites, "richness", fig_label("richness", "Species\nrichness", opts), opts = opts, pred = o$pred),
               site_map(sites, "richness_sd", fig_label("richness_sd", "Posterior SD\n(species)", opts), opts = opts, pred = o$pred, pcol = "richness_sd"), opts, "richness_map"),
      paste0("Estimated species richness (number of modelled species present; posterior mean, A) and its posterior SD (B)", if (!is.null(o$pred)) " on the prediction grid; points = surveyed sites." else " at surveyed sites."))
  }
  # coefficients
  e <- fit$estimates
  co <- e[grepl("^(psi_beta_|p_beta_|comm_psi_beta_|comm_p_beta_)", e$param) & !species_level(e$param, data$species), ]
  if (nrow(co)) {
    co$group <- ifelse(grepl("psi_beta", co$param), fig_label("occupancy_process", "Occupancy", opts), fig_label("detection_process", "Detection", opts))
    ai <- (1 - inner) / 2
    has_d <- co$param %in% colnames(o$draws)
    co$lower50 <- co$upper50 <- NA_real_
    co$lower50[has_d] <- vapply(co$param[has_d], function(p) stats::quantile(o$draws[, p], ai), 0)
    co$upper50[has_d] <- vapply(co$param[has_d], function(p) stats::quantile(o$draws[, p], 1 - ai), 0)
    if (!all(has_d)) co$lower50 <- co$upper50 <- NULL
    co$param <- vapply(co$param, function(p) fig_label(p, param_label(p, data, opts), opts), "")
    # rows grouped by process (occupancy on top), ordered by estimate within a process
    occ <- grepl("psi_beta", e$param[match(rownames(co), rownames(e))])
    lv <- c(co$param[!occ][order(co$estimate[!occ])], co$param[occ][order(co$estimate[occ])])
    co$param <- factor(co$param, levels = unique(lv))
    gl <- c(fig_label("occupancy_process", "Occupancy", opts), fig_label("detection_process", "Detection", opts))
    co$group <- factor(co$group, levels = gl)
    pc <- fig_coef(co, null_value = 0, order = FALSE,
                   xlab = fig_label("coefficients_x", sprintf("Coefficient (logit scale, per 1 SD)"), opts)) +
      ggplot2::aes(shape = group) + ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17), gl), name = NULL)
    attr(pc, "fig_size") <- list(width = fig_eff_width("coefficients", "onehalf", opts),
                                 height = 30 + 9 * nrow(co))
    figs$coefficients <- with_caption(pc,
                                      sprintf("Occupancy (circles) and detection (triangles) coefficients on the logit scale (per 1 SD of the covariate): posterior means with %d%% (thin) and %d%% (thick) credible intervals; dashed line = no association.", pct, round(100 * inner)))
  }
  # maps
  if (o$type != "multi") {
    figs$occupancy_map <- with_caption(
      map_pair(site_map(sites, "psi", fig_label("occupancy_map_fill", "Occupancy\nprobability", opts), opts = opts, pred = o$pred),
               site_map(sites, "psi_sd", fig_label("occupancy_sd_fill", "Posterior SD\n(probability)", opts), opts = opts, pred = o$pred, pcol = "psi_sd"), opts, "occupancy_map"),
      paste0("Predicted occupancy probability (posterior mean, A) and its uncertainty (posterior SD, B)",
             if (!is.null(o$pred)) " on the prediction grid; points = surveyed sites." else " at the surveyed sites."))
    if ("w" %in% names(sites))
      figs$spatial_effect_map <- with_caption(single_map(site_map(sites, "w", fig_label("spatial_effect", "Spatial effect\n(logit scale)", opts), diverging = TRUE, opts = opts, pred = o$pred), opts),
                                              paste0("Posterior mean of the spatial random effect (logit scale)", if (!is.null(o$pred)) " interpolated to the prediction grid" else "", ": ", if (isTRUE(opts$greyscale)) "light = more occupied than the covariates predict, dark = less, mid grey = as predicted" else "red = more occupied than the covariates predict, blue = less, white = as predicted", "; points = surveyed sites."))
  } else if (o$model == "sfMsPGOcc") {
    figs$spatial_effect_map <- with_caption(single_map(site_map(sites, "factor1", fig_label("spatial_factor1", "Spatial\nfactor 1", opts), diverging = TRUE, opts = opts), opts),
                                            "Posterior mean of the first latent spatial factor at surveyed sites; species load on it with their own sign and strength.")
  }
  # range posterior
  rc <- grep("^effective_range_km", colnames(o$draws), value = TRUE)
  if (length(rc)) {
    dr <- dist_range(data$coords)
    rk <- if (!is.null(o$settings$range_prior_km)) unlist(o$settings$range_prior_km) else c(dr[["min"]], dr[["max"]])
    # phi ~ U(3/rk2, 3/rk1) implies range r = 3/phi has density 3 / (r^2 (3/rk1 - 3/rk2)) on [rk1, rk2];
    # drawn exactly (a KDE of prior draws invents a mode at the lower bound)
    rr <- seq(rk[1], rk[2], length.out = 400)
    pd <- do.call(rbind, lapply(rc, function(r) data.frame(param = r, value = rr, dens = 3 / (rr^2 * (3 / rk[1] - 3 / rk[2])))))
    dd <- do.call(rbind, lapply(rc, function(r) data.frame(param = r, value = o$draws[, r])))
    rlab <- vapply(rc, function(r) fig_label(r, if (r == "effective_range_km") "Spatial random effect" else sub("^effective_range_km_f", "Spatial factor ", r), opts), "")
    pr <- ggplot2::ggplot(dd, ggplot2::aes(value)) +
      ggplot2::geom_density(fill = eco_pal$band, alpha = 0.3, colour = eco_pal$band, linewidth = 0.5, bounds = rk) +
      ggplot2::geom_line(data = pd, ggplot2::aes(value, dens), linetype = 2, colour = eco_pal$muted, linewidth = 0.4) +
      (if (length(rc) > 1) ggplot2::facet_wrap(~param, scales = "free", labeller = ggplot2::as_labeller(rlab)) else ggplot2::facet_null()) +
      ggplot2::labs(x = fig_label("range_x", "Effective spatial range (km)", opts), y = fig_label("range_y", "Density", opts))
    figs$range <- with_caption(pr,
                               "Posterior of the effective spatial range (distance at which residual correlation falls to 0.05; filled) and the exact prior density implied by the uniform prior on phi (dashed), both on the prior's support. Overlap with the prior shows how little the data say about the range.")
  }
  if (!is.null(o$gof)) {
    gd <- do.call(rbind, lapply(o$gof, function(g) data.frame(obs = g$obs, rep = g$rep, group = sprintf("%s (p = %.2f)", g$label, g$p))))
    p <- ggplot2::ggplot(gd, ggplot2::aes(obs, rep)) + ggplot2::geom_point(alpha = 0.15, size = 0.4, colour = eco_pal$ink) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = if (isTRUE(opts$greyscale)) "grey45" else eco_pal$truth) +
      ggplot2::labs(x = fig_label("gof_x", "Freeman-Tukey discrepancy, observed data", opts),
                    y = fig_label("gof_y", "Discrepancy, replicated data", opts))
    if (length(o$gof) > 1) p <- p + ggplot2::facet_wrap(~group, scales = "free")
    pv <- paste(vapply(o$gof, function(g) sprintf("%sBayesian p = %.2f", if (length(o$gof) > 1) paste0(g$label, ": ") else "", g$p), ""), collapse = "; ")
    figs$gof <- with_caption(p, sprintf("Posterior predictive check: Freeman-Tukey discrepancy of replicated vs observed site-level detection counts, one point per posterior draw; diagonal line = 1:1. Points spread evenly around the line (Bayesian p near 0.5) indicate adequate fit (%s).", pv))
  }
  core <- colnames(o$draws)[!species_level(colnames(o$draws), data$species)][1:min(9, sum(!species_level(colnames(o$draws), data$species)))]
  pt <- fig_trace(lapply(seq_len(o$n_chains), function(k) {
    n <- nrow(o$draws) / o$n_chains; o$draws[(k - 1) * n + seq_len(n), core, drop = FALSE] }), core) +
    ggplot2::facet_wrap(~ind, scales = "free_y", ncol = 3,
                        labeller = ggplot2::as_labeller(function(x) vapply(x, function(k)
                          paste(strwrap(fig_label(k, param_label(k, data, opts), opts), 26), collapse = "\n"), ""))) +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = ggplot2::rel(0.75)))
  attr(pt, "fig_size") <- list(width = "double", aspect = 0.75)
  pt <- pt + ggplot2::labs(x = fig_label("trace_x", "Saved iteration (after burn-in and thinning)", opts))
  figs$trace <- with_caption(pt, "MCMC traces per chain (post burn-in, thinned) for the main estimands; well-mixed chains overlap with no trend.")
  figs
}

plots_occu_spatial <- function(fit, data) {
  f <- figures_occu_spatial(fit, data, fig_opts())
  f[intersect(names(f), c(grep("^(effect|detection)_", names(f), value = TRUE), "occupancy_map", "richness_map", "community_effects", "range"))]
}
