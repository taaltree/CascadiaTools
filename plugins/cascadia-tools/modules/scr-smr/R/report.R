# scr-smr: design notes, interpretation, figures, methods text, references.

smr_references <- c(
  "Chandler RB, Royle JA (2013) Spatially explicit models for inference about density in unmarked or partially marked populations. Annals of Applied Statistics 7:936-954. doi:10.1214/12-AOAS610",
  "Sollmann R, Gardner B, Parsons AW, Stocking JJ, McClintock BT, Simons TR, Pollock KH, O'Connell AF (2013) A spatial mark-resight model augmented with telemetry data. Ecology 94:553-559. doi:10.1890/12-1256.1",
  "Whittington J, Hebblewhite M, Chandler RB (2018) Generalized spatial mark-resight models with an application to grizzly bears. Journal of Applied Ecology 55:157-168. doi:10.1111/1365-2664.12954",
  "Jimenez J, Chandler R, Tobajas J, Descalzo E, Mateo R, Ferreras P (2019) Generalized spatial mark-resight models with incomplete identification: an application to red fox density estimates. Ecology and Evolution 9:4739-4748. doi:10.1002/ece3.5077",
  "Efford MG, Hunter CM (2018) Spatial capture-mark-resight estimation of animal population density. Biometrics 74:411-420. doi:10.1111/biom.12766",
  "Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press, Waltham MA (ch. 18-19: unmarked and partially marked populations).",
  "Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(Suppl 2):S521-S537. doi:10.1007/s10336-010-0619-4",
  "Efford MG (2025) secr: Spatially explicit capture-recapture models. R package version 5.x (mark-resight vignette secr-markresight.pdf). doi:10.32614/CRAN.package.secr",
  "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487",
  "Sun CC, Fuller AK, Royle JA (2014) Trap configuration and spacing influences parameter estimates in spatial capture-recapture models. PLoS ONE 9:e88025. doi:10.1371/journal.pone.0088025")

smr_fmt <- function(x, d = 2) sub("\\.$", "", trimws(formatC(x, digits = d, format = "fg", flag = "#")))
smr_unit_len <- function(prep) if (is.na(prep$units)) "map units" else prep$units
smr_lab <- function(v, data, opts = NULL) {
  d <- if (!is.null(data$covariate_labels[[v]])) data$covariate_labels[[v]] else v
  fig_label(v, d, opts)
}
smr_ci_word <- function(fit) if (identical(fit$object$engine, "secr")) "CI" else "CrI"

# ---------------------------------------------------------------- design notes
smr_design_notes <- function(prep, est) {
  g <- function(p) est$estimate[est$param == p]
  out <- character()
  out <- c(out, sprintf("Marked animals: %d (%s); %d of them were resighted and identified at least once, %d at more than one camera. Identified sightings: %d; unmarked sightings: %d%s.",
                        prep$m, if (prep$has_marking) "caught in the modelled marking phase" else if (identical(prep$model, "smr") && !is.null(prep$mk)) "marking phase present but IGNORED (marked treated as a random sample)" else "marked before resighting, assumed a random sample of the population",
                        prep$n_marked_resighted, prep$n_spatial_resight, prep$n_id_sightings, prep$n_unmarked_sightings,
                        if (prep$has_tm) sprintf("; marked but unidentified sightings: %d", prep$n_unid_sightings) else ""))
  if (prep$n_marked_resighted < 5) out <- c(out, sprintf("WARNING: only %d marked animals were resighted; sigma (and therefore density) rests on very few animals. Expect wide intervals; see VALIDATION.md (harder design).", prep$n_marked_resighted))
  if (prep$has_tel) out <- c(out, sprintf("Telemetry: %d fixes from %d marked animals inform their activity centres and sigma.", sum(prep$telS[, 1]), sum(prep$telS[, 1] > 0)))
  sig <- g("sigma")
  if (length(sig) && is.finite(sig)) {
    ratio <- prep$trap_spacing / sig
    out <- c(out, sprintf("Camera spacing (%s %s) is %.1f x sigma (recommended <= 2; Sun et al. 2014). State-space buffer is %.1f x sigma (recommended >= 3).",
                          format(signif(prep$trap_spacing, 3)), smr_unit_len(prep), ratio, prep$buffer / sig))
    if (ratio > 2) out <- c(out, "WARNING: cameras are more than 2 sigma apart; few marked animals are photographed at several cameras.")
    if (prep$buffer / sig < 3) out <- c(out, "WARNING: the buffer is less than 3 sigma; density is inflated when animals living outside the state space are photographed.")
  }
  if (identical(prep$model, "smr") && !is.null(prep$mk))
    out <- c(out, "WARNING: marking-phase captures were supplied but the conventional SMR model was requested: if marking was spatially concentrated the density estimate is biased (Whittington et al. 2018). Use model: gsmr.")
  out
}

# ---------------------------------------------------------------- interpret
smr_interpret <- function(fit, data) {
  est <- fit$estimates; prep <- fit$object$prep
  row <- function(p) est[est$param == p, ]
  cw <- smr_ci_word(fit)
  ci <- function(r, d = 2) sprintf("%s (95%% %s %s to %s)", smr_fmt(r$estimate, d), cw, smr_fmt(r$lower, d), smr_fmt(r$upper, d))
  dn <- smr_density_name(prep$units); al <- smr_area_label(prep$units)
  A <- prep$area * smr_area_conv(prep)
  D <- row(dn); N <- row("abundance"); mf <- row("marked_fraction"); sg <- row("sigma"); hr <- row("home_range_area_95")
  rr <- row(smr_rate_name(prep))
  sp <- if (!is.null(data$species)) data$species else "animals"
  out <- c(
    sprintf("Data: %d marked %s (%s); cameras at %d sites over %d resighting occasions recorded %d sightings of identified marked animals, %d of unmarked animals%s. %d of the %d marked animals were photographed.",
            prep$m, sp, if (prep$model == "gsmr") sprintf("caught at %d marking traps over %d occasions", prep$mk$L, prep$mk$Km) else "marked before the camera survey",
            prep$J, prep$K, prep$n_id_sightings, prep$n_unmarked_sightings,
            if (prep$has_tm) sprintf(" and %d of marked animals whose mark could not be read", prep$n_unid_sightings) else "",
            prep$n_marked_resighted, prep$m),
    sprintf("Density was estimated at %s %s per %s%s. Over the %s %s state space (camera array plus buffer) that is %s %s, of which the %d marked animals are a fraction %s.",
            ci(D, 2), sp, al, if (dn == "density_per_km2") sprintf(" (%s per 100 km2)", smr_fmt(100 * D$estimate, 3)) else "",
            smr_fmt(A, 3), al, ci(N, 3), sp, prep$m, ci(mf, 2)),
    sprintf("Sigma, the spatial scale of space use, was %s %s: 95%% of an animal's activity falls within about %s %s of its activity centre (circular home range of %s %s; treat as an index of ranging scale).",
            ci(sg, 3), smr_unit_len(prep), smr_fmt(sqrt(stats::qchisq(0.95, 2)) * sg$estimate, 3), smr_unit_len(prep), ci(hr, 2), al),
    if (prep$pois) sprintf("A camera at an animal's activity centre records it on average %s times per occasion%s; the rate falls off with distance as a half-normal curve (detection-function figure).",
                           ci(rr, 2), if (length(prep$pcov)) " at the reference level of the camera covariates" else "")
    else sprintf("A camera at an animal's activity centre records it on a given occasion with probability %s%s.", ci(rr, 2),
                 if (length(prep$pcov)) " at the reference level of the camera covariates" else "")
  )
  if (prep$has_tm) { r <- row("id_probability")
    out <- c(out, sprintf("Of the photographs of marked animals, a share %s could be identified to individual; the rest were still used as 'marked, identity unknown' counts.", ci(r, 2))) }
  if (prep$model == "gsmr") { r <- row("marking_detection")
    out <- c(out, sprintf("Marking phase: a marking trap at an animal's activity centre caught it on a given occasion with probability %s. Because marking was modelled, marked animals are not assumed to be a random sample: animals living near the marking traps are more likely to carry marks, and the model accounts for that (generalized SMR; Whittington et al. 2018).", ci(r, 2))) }
  else out <- c(out, "Conventional SMR assumes the marked animals are a random sample of the population with respect to where they live. If animals were marked at traps concentrated in part of the area, this assumption fails and density is biased (in our simulations about -30% when marking was centred on the camera grid); supply the marking captures and use model gsmr instead.")
  for (k in seq_along(prep$pcov)) {
    v <- prep$pcov[k]; r <- row(paste0("detection_beta_", v)); sc <- prep$pscale[[v]]; lab <- smr_lab(v, data)
    eff <- exp(c(r$estimate, r$lower, r$upper))
    what <- if (prep$pois) "the sighting rate" else "the odds of a sighting"
    out <- c(out, sprintf("Detection: %s %s multiplied %s by %s (95%% interval %s to %s). This is a detection effect (how easily animals are photographed), not a difference in where animals live.",
                          if (isTRUE(sc$binary)) sprintf("'%s' (1 vs 0)", lab) else sprintf("each 1 SD (%s) increase in '%s'", smr_fmt(sc$scale), lab),
                          "", what, smr_fmt(eff[1]), smr_fmt(eff[2]), smr_fmt(eff[3])))
    out[length(out)] <- gsub("  ", " ", out[length(out)])
  }
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; r <- row(paste0("density_beta_", v)); sc <- prep$dscale[[v]]; lab <- smr_lab(v, data)
    rr2 <- exp(c(r$estimate, r$lower, r$upper))
    s1 <- sprintf("Density and %s: each 1 SD (%s units) higher %s is associated with %s%% %s density (rate ratio %s, 95%% interval %s to %s).",
                  lab, smr_fmt(sc$scale), lab, smr_fmt(abs(100 * (rr2[1] - 1)), 2), if (rr2[1] >= 1) "higher" else "lower",
                  smr_fmt(rr2[1]), smr_fmt(rr2[2]), smr_fmt(rr2[3]))
    pd <- smr_pdraws(fit, 1000)
    if (!is.null(pd)) {
      xr <- stats::quantile(prep$mask[[v]], c(0.1, 0.9))
      lo <- smr_density_at(prep, pd, k, xr[1]); hi <- smr_density_at(prep, pd, k, xr[2])
      s1 <- paste(s1, sprintf("Expected density is %s per %s (%s to %s) where %s = %s (10th percentile of the state space) and %s (%s to %s) where it = %s (90th percentile).",
                              smr_fmt(mean(lo)), al, smr_fmt(smr_q(lo)[1]), smr_fmt(smr_q(lo)[2]), lab, smr_fmt(xr[1], 3),
                              smr_fmt(mean(hi)), smr_fmt(smr_q(hi)[1]), smr_fmt(smr_q(hi)[2]), smr_fmt(xr[2], 3)))
      if (!identical(fit$object$engine, "secr")) s1 <- paste(s1, sprintf("Posterior probability that the association is positive: %.2f.", mean(pd[, paste0("bD", k)] > 0)))
    }
    out <- c(out, paste(s1, "This is an association between density and habitat across the state space, not proof that the habitat causes higher density."))
  }
  out <- c(out,
    "Assumptions: closed population during resighting; marked and unmarked animals share the same detection function and sigma; marks are not lost; every marked animal present is listed (including those never photographed); unmarked photos are counted independently (no double counting of a single event); for conventional SMR, marked animals are a spatially random sample.",
    paste("Validation (modules/scr-smr/VALIDATION.md, 100 simulated datasets per design): with 20 randomly marked animals of ~73 and 49 cameras, SMR estimated density with +5% bias and 98% interval coverage; gSMR with marking at central box traps gave +4.5% bias and 96% coverage.",
          "Fitting conventional SMR to those gSMR data (marked animals treated as a random sample) underestimated density by 31% (coverage 10%), which is why the marking captures should be modelled when marking was done at traps.",
          "With only 5 marked animals of ~41 the posterior mean of density was 15% too high (coverage still 95%); with so few marks, rely on the interval rather than the point estimate.",
          if (identical(prep$model, "gsmr"))
            "Cross-check: secr's joint marking + sighting fit was biased low (-26% abundance, coverage 60%) on these gSMR simulations while the Bayesian gSMR was unbiased, so secr's abundance is not a reference for gSMR; abundance agreed with secr SCR fitted to the marking captures alone (median ratio 0.93), and sigma and detection rates agreed with secr's joint fit (ratios 0.97-1.08)."
          else "Cross-check: for conventional SMR the secr pseudo-likelihood fit agreed with the Bayesian SMR on density (mean ratio 1.02 over 48 datasets)."))
  out
}

# ---------------------------------------------------------------- draws for figures
# Top-level draws in a common form: N, sigma, lam0, delta, pm0, bD1.., bP1..
smr_pdraws <- function(fit, ndraw = 1000) {
  ob <- fit$object; prep <- ob$prep
  if (!is.null(ob$chains)) {
    d <- smr_pool(ob$chains, ndraw)
    cn <- colnames(d); cn <- sub("^bD\\[(\\d+)\\]$", "bD\\1", cn); cn <- sub("^bP\\[(\\d+)\\]$", "bP\\1", cn); colnames(d) <- cn
    return(d)
  }
  # secr: parametric draws from the MLE (multivariate normal on the link scale)
  sf <- ob$secr; b <- sf$fit$par; V <- sf$beta.vcv
  if (is.null(V) || any(!is.finite(V))) return(NULL)
  rn <- rownames(stats::coef(sf))
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  B <- sweep(matrix(stats::rnorm(ndraw * length(b)), ndraw) %*% t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(b))), 2, b, "+")
  colnames(B) <- rn
  f <- 1000; A_ha <- prep$area * f^2 / 1e4
  PD <- length(prep$dcov); PP <- length(prep$pcov)
  out <- matrix(NA_real_, ndraw, 3 + PD + PP)
  colnames(out) <- c("N", "sigma", "lam0", if (PD) paste0("bD", seq_len(PD)), if (PP) paste0("bP", seq_len(PP)))
  for (r in seq_len(ndraw)) {
    bD <- if (PD) B[r, paste0("D.", prep$dcov)] else numeric(0)
    lw <- if (PD) as.vector(prep$XD %*% bD) else rep(0, prep$G)
    out[r, "N"] <- sum(exp(B[r, "D"] + lw)) * (prep$pix * f)^2 / 1e4
    out[r, "sigma"] <- exp(B[r, "sigma"]) / f
    g0 <- B[r, "g0"] + if (prep$has_marking) B[r, "g0.tssighting"] else 0
    out[r, "lam0"] <- stats::plogis(g0)
    if (PD) out[r, paste0("bD", seq_len(PD))] <- bD
    if (PP) out[r, paste0("bP", seq_len(PP))] <- B[r, paste0("g0.", prep$pcov)]
  }
  out
}

# Expected density (reported units) at covariate k = x (others at their mean, i.e. 0)
smr_density_at <- function(prep, pd, k, x) {
  PD <- length(prep$dcov); sc <- prep$dscale[[prep$dcov[k]]]
  xs <- (x - sc$center) / sc$scale
  vapply(seq_len(nrow(pd)), function(r) {
    b <- pd[r, paste0("bD", seq_len(PD))]
    lw <- as.vector(prep$XD %*% b)
    pd[r, "N"] * exp(b[k] * xs - smr_lse(lw)) / (prep$pix^2 * smr_area_conv(prep))
  }, 1)
}

# ---------------------------------------------------------------- figures
smr_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  ob <- fit$object; prep <- ob$prep
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; a <- (1 - ci) / 2
  ici <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci; ai <- (1 - ici) / 2
  al <- smr_area_label(prep$units); ul <- smr_unit_len(prep)
  bayes <- !identical(ob$engine, "secr")
  bw <- if (bayes) "credible" else "confidence (parametric bootstrap)"
  pd <- smr_pdraws(fit, 1000)
  figs <- list()
  to_x <- function(v) v * prep$cscale + prep$origin[1]; to_y <- function(v) v * prep$cscale + prep$origin[2]
  xl <- fig_label("easting", sprintf("Easting (%s)", ul), opts); yl <- fig_label("northing", sprintf("Northing (%s)", ul), opts)
  cellA <- prep$pix^2 * smr_area_conv(prep)
  dev_pts <- data.frame(x = prep$cameras$x, y = prep$cameras$y, type = fig_label("cameras", "Camera", opts))
  if (prep$model == "gsmr") dev_pts <- rbind(dev_pts, data.frame(x = prep$mk$traps$x, y = prep$mk$traps$y, type = fig_label("marking_traps", "Marking trap", opts)))
  add_dev <- function(p) p + ggplot2::geom_point(data = dev_pts, ggplot2::aes(x, y, shape = type), inherit.aes = FALSE,
                                                 fill = "white", colour = "black", size = 1.1, stroke = 0.3) +
    ggplot2::scale_shape_manual(values = c(21, 24)[seq_along(unique(dev_pts$type))], name = NULL)
  surf <- ob$surfaces
  if (!is.null(surf)) {
    grid <- data.frame(x = prep$mask$x, y = prep$mask$y, value = surf$realised / cellA)
    p <- add_dev(fig_map(grid, sprintf("Animals per %s", al)) + ggplot2::labs(x = xl, y = yl))
    if (isTRUE(opts$show_data) && prep$m) {
      cen <- surf$centres
      p <- p + ggplot2::geom_point(data = data.frame(x = to_x(cen[, 1]), y = to_y(cen[, 2])), ggplot2::aes(x, y),
                                   inherit.aes = FALSE, shape = 4, colour = "#D55E00", size = 1.3, stroke = 0.6)
    }
    figs$density_surface <- with_caption(p, sprintf(
      "Posterior mean realised density (animals per %s): marked and unmarked activity centres combined, averaged over %d posterior draws. Circles: cameras%s; orange crosses: posterior mean activity centres of the %d marked animals.",
      al, ob$settings$surface_draws, if (prep$model == "gsmr") ", triangles: marking traps" else "", prep$m))
    attr(figs$density_surface, "fig_size") <- list(width = "onehalf", aspect = 0.8)
    if (length(prep$dcov)) {
      g2 <- data.frame(x = prep$mask$x, y = prep$mask$y, value = surf$expected / cellA)
      figs$density_expected_surface <- with_caption(add_dev(fig_map(g2, sprintf("Expected animals per %s", al)) + ggplot2::labs(x = xl, y = yl)),
        sprintf("Expected density surface implied by the habitat model (%s), posterior mean.", paste(vapply(prep$dcov, smr_lab, "", data = data, opts = opts), collapse = ", ")))
      attr(figs$density_expected_surface, "fig_size") <- list(width = "onehalf", aspect = 0.8)
    }
  }
  # conditional effects of density covariates
  if (!is.null(pd)) for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]
    xs <- seq(min(prep$mask[[v]]), max(prep$mask[[v]]), length.out = 60)
    Dm <- sapply(xs, function(x) smr_density_at(prep, pd, k, x))
    pred <- data.frame(x = xs, estimate = colMeans(Dm), lower = apply(Dm, 2, stats::quantile, a), upper = apply(Dm, 2, stats::quantile, 1 - a),
                       lower50 = apply(Dm, 2, stats::quantile, ai), upper50 = apply(Dm, 2, stats::quantile, 1 - ai))
    obs <- if (isTRUE(opts$show_data)) data.frame(x = prep$mask[[v]][smr_thin_idx(prep$G, 300)]) else NULL
    figs[[paste0("density_", v)]] <- with_caption(
      fig_cond_effect(pred, smr_lab(v, data, opts), fig_label("density", sprintf("Density (animals per %s)", al), opts), obs = obs) +
        ggplot2::scale_y_continuous(limits = c(0, NA)),
      sprintf("Expected density versus %s (other density covariates at their state-space mean), with %d%% and %d%% %s bands. Rug: covariate values across the state space.",
              smr_lab(v, data, opts), round(100 * ci), round(100 * ici), bw))
  }
  # detection function(s): resighting (and marking for gSMR)
  if (!is.null(pd)) {
    dmax <- 4 * stats::median(pd[, "sigma"]) * prep$cscale
    dd <- seq(0, dmax, length.out = 80)
    bin <- prep$pcov[vapply(prep$pcov, function(v) isTRUE(prep$pscale[[v]]$binary), TRUE)]
    levs <- if (length(bin)) c(0, 1) else 0
    rows <- list()
    for (lv in levs) {
      eta <- if (length(bin)) lv * pd[, paste0("bP", match(bin[1], prep$pcov))] else 0
      r0 <- if (prep$pois) pd[, "lam0"] * exp(eta) else stats::plogis(stats::qlogis(pd[, "lam0"]) + eta)
      P <- sapply(dd, function(x) r0 * exp(-(x / prep$cscale)^2 / (2 * pd[, "sigma"]^2)))
      r <- data.frame(x = dd, estimate = colMeans(P), lower = apply(P, 2, stats::quantile, a), upper = apply(P, 2, stats::quantile, 1 - a),
                      lower50 = apply(P, 2, stats::quantile, ai), upper50 = apply(P, 2, stats::quantile, 1 - ai))
      r$group <- if (length(bin)) sprintf("Camera, %s = %d", smr_lab(bin[1], data, opts), lv) else "Camera resighting"
      rows[[length(rows) + 1]] <- r
    }
    if (prep$model == "gsmr" && "pm0" %in% colnames(pd)) {
      P <- sapply(dd, function(x) pd[, "pm0"] * exp(-(x / prep$cscale)^2 / (2 * pd[, "sigma"]^2)))
      rows[[length(rows) + 1]] <- data.frame(x = dd, estimate = colMeans(P), lower = apply(P, 2, stats::quantile, a), upper = apply(P, 2, stats::quantile, 1 - a),
                                             lower50 = apply(P, 2, stats::quantile, ai), upper50 = apply(P, 2, stats::quantile, 1 - ai),
                                             group = "Marking trap (capture probability)")
    }
    dfd <- do.call(rbind, rows)
    if (length(unique(dfd$group)) == 1) dfd$group <- NULL
    pdet <- fig_cond_effect(dfd, fig_label("distance", sprintf("Distance from activity centre to detector (%s)", ul), opts),
                            fig_label("detection_rate", paste0(if (prep$pois) "Expected sightings per occasion" else "Detection probability per occasion",
                                   if (prep$model == "gsmr" && "pm0" %in% colnames(pd)) "\n(trap: capture probability)" else ""), opts)) +
      ggplot2::geom_vline(xintercept = prep$trap_spacing, linetype = 2, colour = "grey45") +
      ggplot2::annotate("text", x = prep$trap_spacing, y = Inf, label = " camera spacing", hjust = 0, vjust = 1.3, size = 2.6, colour = "grey30") +
      ggplot2::scale_y_continuous(limits = c(0, NA))
    # one panel per curve (direct strip labels; no legend squeezing a single-column figure)
    if (!is.null(dfd$group)) {
      dfd$group <- factor(dfd$group, levels = unique(dfd$group))
      pdet <- pdet + ggplot2::facet_wrap(~group, ncol = 1, scales = "free_y") +
        ggplot2::guides(colour = "none", fill = "none")
    }
    pdet <- pdet + ggplot2::scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
                                               breaks = scales::breaks_pretty(n = 4))
    figs$detection_function <- with_caption(pdet, sprintf(
      "Half-normal detection functions with %d%% and %d%% %s bands: %s%s. Dashed line: camera spacing (recommended <= 2 sigma).",
      round(100 * ci), round(100 * ici), bw,
      if (prep$pois) "expected number of camera sightings per occasion" else "camera detection probability per occasion",
      if (prep$model == "gsmr") "; and per-occasion capture probability at a marking trap" else ""))
    if (!is.null(dfd$group)) attr(figs$detection_function, "fig_size") <- list(aspect = 0.45 * length(levels(dfd$group)) + 0.1)
  }
  # abundance and marked fraction posteriors
  if (!is.null(ob$rep_chains)) {
    R <- do.call(rbind, ob$rep_chains)
    Nd <- R[, "abundance"]; qq <- smr_q(Nd, c(a, 1 - a))
    figs$abundance_posterior <- with_caption(
      ggplot2::ggplot(data.frame(N = Nd), ggplot2::aes(N)) +
        ggplot2::geom_histogram(binwidth = 1, fill = eco_pal$band, colour = NA, alpha = 0.75) +
        ggplot2::geom_vline(xintercept = prep$m, linetype = 3, colour = "grey30") +
        ggplot2::geom_vline(xintercept = qq, linetype = 2, colour = eco_pal$ink) +
        ggplot2::labs(x = fig_label("abundance", "Abundance in the state space (N)", opts), y = "Posterior draws"),
      sprintf("Posterior distribution of the number of animals in the %s %s state space; dashed lines: %d%% credible interval; dotted line: number of marked animals (%d), the lower bound.",
              smr_fmt(prep$area * smr_area_conv(prep), 3), al, round(100 * ci), prep$m))
    mfd <- R[, "marked_fraction"]
    figs$marked_fraction <- with_caption(
      ggplot2::ggplot(data.frame(v = mfd), ggplot2::aes(v)) +
        ggplot2::geom_density(fill = eco_pal$band, alpha = 0.3, colour = eco_pal$band) +
        ggplot2::geom_vline(xintercept = smr_q(mfd, c(a, 1 - a)), linetype = 2, colour = eco_pal$ink) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = fig_label("marked_fraction", "Fraction of the population carrying marks", opts), y = "Posterior density"),
      sprintf("Posterior of the marked fraction m/N (m = %d marked animals); dashed lines: %d%% credible interval.", prep$m, round(100 * ci)))
  } else {
    r <- fit$estimates[fit$estimates$param == "abundance", ]
    figs$abundance_posterior <- with_caption(
      ggplot2::ggplot(r, ggplot2::aes(estimate, "N")) + ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper)) +
        ggplot2::geom_vline(xintercept = prep$m, linetype = 3) + ggplot2::labs(x = "Abundance in the state space (N)", y = NULL),
      "Expected abundance in the state space (secr) with 95% confidence interval; dotted line: number of marked animals.")
  }
  # sightings per marked individual (observed vs posterior predictive) + PPC panel
  ppc <- ob$ppc
  if (!is.null(ppc) && prep$m) {
    cap10 <- function(v) table(factor(pmin(v, 10), levels = 0:10))
    obsn <- cap10(rowSums(prep$yID))
    repn <- t(vapply(ppc$sight_rep, function(z) as.numeric(cap10(z$per_ind)), numeric(11)))
    dfc <- data.frame(k = factor(c(0:9, "10+"), levels = c(0:9, "10+")), obs = as.numeric(obsn),
                      est = colMeans(repn), lo = apply(repn, 2, stats::quantile, a), hi = apply(repn, 2, stats::quantile, 1 - a))
    figs$captures_per_individual <- with_caption(
      ggplot2::ggplot(dfc, ggplot2::aes(k)) +
        ggplot2::geom_col(ggplot2::aes(y = obs), fill = "grey75", width = 0.7) +
        ggplot2::geom_pointrange(ggplot2::aes(y = est, ymin = lo, ymax = hi), colour = eco_pal$ink, size = 0.3) +
        ggplot2::labs(x = "Identified sightings per marked animal", y = "Number of marked animals"),
      sprintf("Bars: observed number of marked animals identified 0, 1, 2, ... times on camera. Points and lines: posterior predictive mean and %d%% interval. Bars outside the intervals indicate unmodelled heterogeneity among marked animals.", round(100 * ci)))
    # unmarked sightings by camera: observed vs expected
    EU <- sapply(ppc$sight_rep, function(z) z$EU)
    RU <- sapply(ppc$sight_rep, function(z) z$per_cam_unmarked)
    dcam <- data.frame(obs = prep$nU, est = rowMeans(EU), lo = apply(RU, 1, stats::quantile, a), hi = apply(RU, 1, stats::quantile, 1 - a))
    figs$unmarked_by_camera <- with_caption(
      ggplot2::ggplot(dcam, ggplot2::aes(est, obs)) +
        ggplot2::geom_abline(linetype = 2, colour = "grey45") +
        ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi), colour = "grey65", linewidth = 0.4) +
        ggplot2::geom_point(size = 1.2, colour = eco_pal$ink) +
        ggplot2::labs(x = "Expected unmarked sightings\n(posterior mean per camera)", y = "Observed unmarked sightings"),
      sprintf("Unmarked sightings at each camera: observed versus expected under the fitted model; horizontal lines: %d%% posterior predictive interval for that camera; dashed: 1:1.", round(100 * ci)))
    nm <- c(unmarked_total = "Unmarked sightings (total)", identified_total = "Identified sightings (total)",
            marked_resighted = "Marked animals resighted", unidentified_total = "Marked, unidentified sightings")
    nm <- nm[names(nm) %in% colnames(ppc$rep)]
    int_breaks <- function(l) { b <- unique(round(pretty(l))); b[b >= l[1] & b <= l[2]] }
    pl <- lapply(names(nm), function(s) fig_ppc(ppc$rep[, s], stats::median(ppc$obs[, s], na.rm = TRUE), nm[[s]]) +
                   ggplot2::scale_x_continuous(breaks = int_breaks))
    for (s in c("ft_unmarked_cameras", "ft_marked_individuals")) {
      ftd <- data.frame(o = ppc$obs[, s], r = ppc$rep[, s])
      pl[[length(pl) + 1]] <- ggplot2::ggplot(ftd, ggplot2::aes(o, r)) + ggplot2::geom_point(alpha = 0.4, size = 0.8, colour = eco_pal$ink) +
        ggplot2::geom_abline(linetype = 2, colour = "grey45") +
        ggplot2::labs(x = sprintf("FT discrepancy, observed\n(%s)", if (s == "ft_unmarked_cameras") "unmarked by camera" else "marked by animal"),
                      y = "FT, replicated", caption = sprintf("P(replicate >= observed) = %.2f", ppc$p[[s]]))
    }
    figs$gof <- with_caption(compose_figures(pl, ncol = 2) & theme_eco(opts$base_size, opts$font, opts$style, opts$legend),
      "Posterior predictive checks: histograms of statistics from data replicated under the fitted model (line = observed); scatter panels: Freeman-Tukey discrepancy of unmarked sightings by camera and of identified sightings by marked animal, observed vs replicated (points spread on both sides of the 1:1 line = consistent with the model).")
    attr(figs$gof, "fig_size") <- list(width = "double", aspect = 1.15)
  }
  # coefficients
  ec <- fit$estimates[grepl("^(density|detection)_beta_", fit$estimates$param), ]
  if (nrow(ec)) {
    cv_of <- sub("^(density|detection)_beta_", "", ec$param)
    ec$param <- sprintf("%s: %s", ifelse(grepl("^density", ec$param), "Density (per SD)",
                                         if (prep$pois) "Sighting rate" else "Detection"),
                        vapply(cv_of, function(v) smr_lab(v, data, opts), ""))
    pco <- fig_coef(ec, null_value = 0, xlab = "Coefficient (link scale)", opts = opts) +
      ggplot2::scale_y_discrete(labels = function(x) vapply(strwrap(x, 28, simplify = FALSE), paste, "", collapse = "\n"))
    attr(pco, "fig_size") <- list(height = 45 + 12 * nrow(ec))
    figs$coefficients <- with_caption(pco,
      "Link-scale coefficients with 95% intervals: density covariates on the log scale per 1 SD (exp(beta) = multiplicative change in density); camera covariates on the log (Poisson) or logit (Bernoulli) scale (0/1 covariates: level 1 vs 0). Dashed line: no effect.")
  }
  if (!is.null(ob$comparison)) {
    cm <- ob$comparison
    cm <- cm[cm$param %in% c(smr_density_name(prep$units), "abundance", "sigma", smr_rate_name(prep), "id_probability",
                             "marking_detection", "marked_fraction", grep("^density_beta_", cm$param, value = TRUE)), ]
    plab <- c(density_per_km2 = "Density (per km\u00b2)", density_per_unit2 = "Density (per map unit\u00b2)", abundance = "Abundance (N)",
              sigma = sprintf("sigma (%s)", ul), resight_rate = "Resight rate lambda0", resight_detection = "Resight probability p0",
              id_probability = "P(mark identified)", marking_detection = "Marking capture prob. p0", marked_fraction = "Marked fraction")
    cm$plab <- ifelse(cm$param %in% names(plab), plab[cm$param],
                      paste("Density effect:", vapply(sub("^density_beta_", "", cm$param), function(v) smr_lab(v, data, opts), "")))
    figs$engine_comparison <- with_caption(
      ggplot2::ggplot(cm, ggplot2::aes(estimate, engine, colour = engine)) +
        ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.25) +
        ggplot2::facet_wrap(~plab, scales = "free_x", ncol = 3, labeller = ggplot2::label_wrap_gen(26)) + ggplot2::labs(x = "Estimate (95% interval)", y = NULL) +
        ggplot2::scale_x_continuous(breaks = function(l) { b <- pretty(l, 3); r <- diff(l); b[b >= l[1] + 0.08 * r & b <= l[2] - 0.08 * r] }) +
        ggplot2::theme(legend.position = "none", panel.spacing.x = grid::unit(1.5, "lines")),
      paste("The same data fitted by the Bayesian data-augmentation model (Nimble) and by maximum pseudo-likelihood (secr mark-resight; Efford & Hunter 2018). secr estimates expected abundance and treats unmarked counts as independent Poisson counts with an overdispersion adjustment.",
            if (identical(prep$model, "gsmr")) "For gSMR, secr's joint fit underestimated abundance by about 26% in validation simulations (VALIDATION.md section 8), so a lower secr abundance here is expected; compare sigma and detection rates." else "Small differences are expected."))
    attr(figs$engine_comparison, "fig_size") <- list(width = "double", aspect = 0.6)
  }
  if (!is.null(ob$rep_chains)) {
    keep <- setdiff(colnames(ob$rep_chains[[1]]), c("abundance_unmarked", "home_range_area_95"))
    wrapnm <- function(mm) { colnames(mm) <- sub("^(density|detection)_beta_", "\\1_beta\n", colnames(mm)); mm }
    figs$trace <- with_caption(fig_trace(lapply(ob$rep_chains, function(mm) wrapnm(mm[smr_thin_idx(nrow(mm), 500), keep, drop = FALSE]))),
                               "MCMC traces of the reported quantities, one colour per chain (thinned for display).")
    attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.75)
    figs$posterior <- with_caption(fig_posterior(wrapnm(do.call(rbind, ob$rep_chains)[, keep, drop = FALSE])), "Posterior densities of the reported quantities.")
    attr(figs$posterior, "fig_size") <- list(width = "double", aspect = 0.75)
  }
  figs
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------- methods
smr_methods_text <- function(fit, data, settings) {
  ob <- fit$object; prep <- ob$prep; st <- ob$settings
  al <- smr_area_label(prep$units)
  dens <- if (length(prep$dcov)) sprintf("log-linear in %s (standardised over the state space)", paste(prep$dcov, collapse = ", ")) else "homogeneous"
  detc <- if (length(prep$pcov)) sprintf(", with camera-level covariate(s) %s on the %s scale", paste(prep$pcov, collapse = ", "), if (prep$pois) "log" else "logit") else ""
  mod <- if (prep$model == "gsmr") sprintf(paste(
    "We used a generalized spatial mark-resight model (Whittington et al. 2018; Jimenez et al. 2019) in which the marking process is modelled jointly with resighting:",
    "captures at %d marking traps over %d occasions were Binomial with probability pm0 exp(-d^2 / (2 sigma^2)), for every individual in the population (zero histories for animals never caught), so that the spatial distribution of marked animals is not assumed random."),
    prep$mk$L, prep$mk$Km)
  else "We used a spatial mark-resight model (Chandler & Royle 2013; Sollmann et al. 2013) in which the marked animals are a known set assumed to be a spatially random sample of the population."
  base <- sprintf(paste(
    "%s Activity centres follow a Poisson point process over a state space of %d habitat pixels (%s %s spacing, %s %s, extending %s %s beyond the detectors) with intensity %s.",
    "At %d cameras operated for up to %d occasions, sightings of an individual at camera j were %s with %s lambda0_j exp(-d^2 / (2 sigma^2))%s.",
    "Sightings of marked animals were identified to individual%s; sightings of unmarked animals were modelled as camera-level counts with expectation the sum of the rates of all unmarked individuals%s.%s"),
    mod, prep$G, format(signif(prep$spacing, 3)), smr_unit_len(prep), smr_fmt(prep$area * smr_area_conv(prep), 3), al,
    format(signif(prep$buffer, 3)), smr_unit_len(prep), dens, prep$J, prep$K,
    if (prep$pois) "Poisson" else "Bernoulli per occasion", if (prep$pois) "expected count" else "probability", detc,
    if (prep$has_tm) " with probability delta; sightings of marked animals whose identity could not be read were modelled as camera-level Poisson counts with expectation (1 - delta) times the summed rates of marked animals (Sollmann et al. 2013)" else "",
    if (!prep$pois) " (a Poisson approximation to the sum of Bernoulli detections; Chandler & Royle 2013)" else "",
    if (prep$has_tel) " Telemetry fixes of marked animals were bivariate normal around their activity centres with the same sigma (Sollmann et al. 2013)." else "")
  eng <- if (ob$engine == "secr") sprintf("The model was fitted by maximum pseudo-likelihood in secr %s (Efford & Hunter 2018), with variances of unmarked and unidentified counts adjusted for overdispersion estimated from %d simulations.",
                                          as.character(utils::packageVersion("secr")), st$secr_nsim)
  else sprintf(paste(
    "Inference was Bayesian with data augmentation (Royle & Dorazio 2012; M = %d rows%s), latent activity centres and inclusion indicators, in NIMBLE %s (de Valpine et al. 2017): %d chains of %d iterations after %d burn-in, thinned by %d,",
    "with priors psi ~ U(0,1), sigma ~ U(0, %s %s), %s ~ U(0, %s)%s%s and Normal(0, SD 2) for covariate coefficients.",
    "Abundance N is the number of activity centres in the state space, density N / area and the marked fraction m / N.",
    "Convergence was assessed with rank-normalised R-hat and bulk effective sample size; model fit with posterior predictive Freeman-Tukey checks of unmarked sightings by camera and identified sightings by marked animal."),
    ob$M, if (prep$model == "smr") sprintf("; the %d marked animals are included with certainty and psi applies to the remaining rows", prep$m) else "",
    as.character(utils::packageVersion("nimble")), st$nchains, st$niter - st$nburnin, st$nburnin, st$thin,
    format(signif(prep$sigma_max * prep$cscale, 3)), smr_unit_len(prep), if (prep$pois) "lambda0" else "p0", if (prep$pois) st$lam0_max else 1,
    if (prep$has_tm) ", delta ~ U(0,1)" else "", if (prep$model == "gsmr") ", pm0 ~ U(0,1)" else "")
  paste(base, eng)
}
