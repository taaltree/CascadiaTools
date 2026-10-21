# axis breaks that keep tick labels away from the panel edges (labels at the
# edge collide with the neighbouring facet or are clipped at the figure border)
so_inner_breaks <- function(l) {
  b <- pretty(l, 4); r <- diff(l)
  b[b >= l[1] + 0.05 * r & b <= l[2] - 0.05 * r]
}

# scr-open: design notes, interpretation, figures, methods text, references.

so_references <- c(
  "Gardner B, Reppucci J, Lucherini M, Royle JA (2010) Spatially explicit inference for open populations: estimating demographic parameters from camera-trap studies. Ecology 91:3376-3383. doi:10.1890/09-0804.1",
  "Royle JA, Chandler RB, Sollmann R, Gardner B (2014) Spatial Capture-Recapture. Academic Press, Waltham MA (chapter 16: open population models).",
  "Efford MG, Schofield MR (2020) A spatial open-population capture-recapture model. Biometrics 76:392-402. doi:10.1111/biom.13150",
  "Glennie R, Borchers DL, Murchie M, Harmsen BJ, Foster RJ (2019) Open population maximum likelihood spatial capture-recapture. Biometrics 75:1345-1355. doi:10.1111/biom.13078",
  "Schwarz CJ, Arnason AN (1996) A general methodology for the analysis of capture-recapture experiments in open populations. Biometrics 52:860-873. doi:10.2307/2533048",
  "Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(Suppl 2):S521-S537. doi:10.1007/s10336-010-0619-4",
  "Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS: a hierarchical perspective. Academic Press (Jolly-Seber model via data augmentation, restricted dynamic occupancy formulation).",
  "Borchers DL, Efford MG (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64:377-385. doi:10.1111/j.1541-0420.2007.00927.x",
  "Turek D, de Valpine P, Paciorek CJ (2016) Efficient Markov chain Monte Carlo sampling for hierarchical hidden Markov models. Environmental and Ecological Statistics 23:549-564. doi:10.1007/s10651-016-0353-z",
  "Milleret C, Dupont P, Bonenfant C, Broseth H, Flagstad O, Sutherland C, Bischof R (2019) A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363. doi:10.1002/ece3.4751",
  "Sun CC, Fuller AK, Royle JA (2014) Trap configuration and spacing influences parameter estimates in spatial capture-recapture models. PLoS ONE 9:e88025. doi:10.1371/journal.pone.0088025",
  "Efford MG (2024) openCR: Open population capture-recapture. R package version 2.2. doi:10.32614/CRAN.package.openCR",
  "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487",
  "Phan D, Pradhan N, Jankowiak M (2019) Composable effects for flexible and accelerated probabilistic programming in NumPyro. arXiv:1912.11554")

so_q <- function(x, p = c(0.025, 0.975)) stats::quantile(x, p, na.rm = TRUE, names = FALSE)
so_fmt <- function(x, d = 2) sub("\\.$", "", formatC(x, digits = d, format = "fg", flag = "#"))
so_lab <- function(v, data, opts = NULL) {
  d <- if (!is.null(data$covariate_labels[[v]])) data$covariate_labels[[v]] else v
  fig_label(v, d, opts)
}
so_draws <- function(fit) do.call(rbind, fit$object$rep_chains)

so_design_notes <- function(prep, est) {
  g <- function(p) est$estimate[est$param == p]
  sig <- g("sigma"); out <- character()
  if (!length(sig) || !is.finite(sig)) return(out)
  ratio <- prep$trap_spacing / sig
  out <- c(out, sprintf("Trap spacing (%s %s) is %.1f x sigma; spacing of about 1-2 sigma or less is recommended (Sun et al. 2014).",
                        format(signif(prep$trap_spacing, 3)), so_unit_len(prep), ratio))
  if (ratio > 2) out <- c(out, "WARNING: traps are more than 2 sigma apart; few spatial recaptures, sigma and density may be imprecise or biased.")
  if (!is.null(prep$buffer) && is.finite(prep$buffer)) {
    br <- prep$buffer / sig
    out <- c(out, sprintf("State-space buffer is %.1f x sigma beyond the outermost traps.", br))
    if (br < 3) out <- c(out, "WARNING: the buffer is less than 3 sigma; animals living outside the state space can be detected, inflating density. Use a buffer of at least 3-4 sigma.")
  }
  out <- c(out, sprintf("%d of %d detected individuals were detected in more than one period (between-period recaptures carry the information on survival).",
                        prep$n_recaptured_between, prep$n))
  if (prep$n_recaptured_between < 10) out <- c(out, "WARNING: fewer than 10 individuals were detected in more than one period: survival and recruitment are weakly identified.")
  if (any(prep$n_by_period < 10)) out <- c(out, sprintf("WARNING: fewer than 10 individuals detected in period(s) %s; abundance there leans on the survival/recruitment model.",
                                                        paste(prep$period_labels[prep$n_by_period < 10], collapse = ", ")))
  # Does the trapped area represent the habitat range? (density covariate
  # effects and total abundance extrapolate to habitat the traps never sampled)
  if (length(prep$dcov) && !is.null(prep$mxy) && !is.null(prep$txy)) {
    s_int <- sig / prep$cscale
    dmin <- sqrt(apply(prep$mxy, 1, function(p) min((prep$txy[, 1] - p[1])^2 + (prep$txy[, 2] - p[2])^2)))
    near <- dmin <= 2 * s_int
    if (sum(near) >= 5) for (k in seq_along(prep$dcov)) {
      x <- prep$XD[, k]
      qa <- stats::quantile(x, c(0.1, 0.9)); qn <- stats::quantile(x[near], c(0.1, 0.9))
      cover <- max(0, min(qa[2], qn[2]) - max(qa[1], qn[1])) / (qa[2] - qa[1])
      dm <- mean(x[near]) - mean(x)
      out <- c(out, sprintf("Habitat sampled by the traps (within 2 sigma): mean %s is %+.2f SD from the state-space mean; the traps cover %.0f%% of its 10th-90th percentile range.",
                            prep$dcov[k], dm, 100 * cover))
    }
    b <- est$estimate[match(paste0("density_beta_", prep$dcov), est$param)]
    if (sum(near) >= 5 && all(is.finite(b))) {
      w <- exp(as.vector(prep$XD %*% b))
      ratio <- mean(w) / mean(w[near]); far <- sum(w[!near]) / sum(w)
      out <- c(out, sprintf("Under the fitted density model, %.0f%% of animals have activity centres more than 2 sigma from any trap, and average density over the state space is %.2f x the density around the traps.", 100 * far, ratio))
      if (ratio > 1.5 || ratio < 1 / 1.5)
        out <- c(out, sprintf("WARNING: the traps sample habitat of unrepresentative density (state-space average %.2f x the trapped area), so total abundance rests on extrapolating the density covariate effect to habitat the traps did not cover; a homogeneous-density model would likely give abundance roughly %.1f times %s (a rough guide, not a prediction). Check the covariate map and effect carefully and, in future, place traps across the full habitat range.", ratio, max(ratio, 1 / ratio), if (ratio > 1) "lower" else "higher"))
    }
  }
  out
}

# ---------------------------------------------------------------- interpret
so_interpret <- function(fit, data) {
  est <- fit$estimates; prep <- fit$object$prep; TT <- prep$TT; lab <- prep$period_labels
  row <- function(p) est[est$param == p, ]
  nm <- so_names(prep); al <- so_area_label(prep$units)
  ci <- function(r, d = 2) sprintf("%s (95%% CrI %s to %s)", so_fmt(r$estimate, d), so_fmt(r$lower, d), so_fmt(r$upper, d))
  dr <- so_draws(fit)
  Nt <- dr[, nm$N, drop = FALSE]
  chg <- Nt[, TT] / pmax(Nt[, 1], 1)
  pdec <- mean(chg < 1)
  Ak <- prep$area * so_area_conv(prep)
  out <- c(
    sprintf("Data: %d individuals were detected (%s per period: %s) at %d traps over %d primary periods of %d occasions; %d individuals were detected in more than one period.",
            prep$n, "individuals", paste(sprintf("%s: %d", lab, prep$n_by_period), collapse = ", "), prep$J, TT, prep$K, prep$n_recaptured_between),
    sprintf("Abundance in the %s %s state space (traps plus buffer) was %s in period %s and %s in period %s (density %s and %s animals per %s). Over the study the population changed by a factor of %s (95%% CrI %s to %s); posterior probability of a decline from the first to the last period: %s.",
            so_fmt(Ak, 3), al, ci(row(nm$N[1]), 3), lab[1], ci(row(nm$N[TT]), 3), lab[TT],
            so_fmt(row(nm$D[1])$estimate, 2), so_fmt(row(nm$D[TT])$estimate, 2), al,
            so_fmt(mean(chg), 3), so_fmt(so_q(chg)[1], 3), so_fmt(so_q(chg)[2], 3),
            if (pdec < 0.005) "< 0.01" else if (pdec > 0.995) "> 0.99" else sprintf("%.2f", pdec)),
    sprintf("Apparent survival between consecutive periods was %s: an animal alive in one period had an estimated %.0f%% chance of still being alive and resident in the next. 'Apparent' because death and permanent emigration from the study area cannot be told apart. Expected residence time is about %s periods (1 / (1 - survival)).",
            ci(row("survival")), 100 * row("survival")$estimate, so_fmt(1 / (1 - row("survival")$estimate), 2)),
    sprintf("Recruitment (births plus immigration of new animals into the state space): %s new animals per period (%s); per-capita recruitment rates (recruits per animal alive in the previous period) were %s.",
            paste(sprintf("%s in %s", vapply(nm$B, function(p) so_fmt(row(p)$estimate, 3), ""), lab[-1]), collapse = ", "),
            "posterior means",
            paste(sprintf("%s: %s", lab[-1], vapply(nm$f, function(p) ci(row(p)), "")), collapse = "; ")),
    sprintf("Period-to-period growth rates (abundance this period / abundance in the previous period): %s. A growth rate of 1 means a stable population; the change reflects survival (%s) plus per-capita recruitment.",
            paste(sprintf("%s: %s", lab[-1], vapply(nm$lam, function(p) ci(row(p), 3), "")), collapse = "; "), so_fmt(row("survival")$estimate)),
    sprintf("In total %s different animals were alive in the state space at some point during the study (superpopulation), of whom %d (%.0f%%) were ever detected.",
            ci(row("superpopulation_size"), 3), prep$n, 100 * prep$n / row("superpopulation_size")$estimate),
    sprintf("Detection: a trap at an animal's activity centre detects it on a given occasion with probability %s; the spatial scale sigma was %s %s, so 95%% of activity falls within %s %s of the centre (home-range area about %s %s, assuming circular bivariate-normal space use).",
            ci(row("baseline_detection")), ci(row("sigma"), 3), so_unit_len(prep), so_fmt(2.448 * row("sigma")$estimate, 3), so_unit_len(prep),
            ci(row("home_range_area_95")), al))
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; r <- row(paste0("density_beta_", v)); sc <- prep$dscale[[v]]; lb <- so_lab(v, data)
    rr <- exp(c(r$estimate, r$lower, r$upper))
    b <- dr[, paste0("density_beta_", v)]
    xr <- stats::quantile(prep$mask[[v]], c(0.1, 0.9))
    ratio <- exp(b * (xr[2] - xr[1]) / sc$scale)
    if (isTRUE(sc$binary)) {
      out <- c(out, sprintf("Density and %s (0/1): density where %s = 1 is %s times that where it is 0 (95%% CrI %s to %s); posterior probability the association is positive: %.2f. This describes where animals' activity centres are, not proof that %s causes the difference.",
                            lb, lb, so_fmt(rr[1]), so_fmt(rr[2]), so_fmt(rr[3]), mean(b > 0), lb))
      next
    }
    out <- c(out, sprintf("Density and %s: each 1 SD (%s) higher %s is associated with %s%% %s density (rate ratio %s, 95%% CrI %s to %s). Between the 10th (%s) and 90th (%s) percentile of %s across the state space, expected density differs by a factor of %s (%s to %s). Posterior probability the association is positive: %.2f. This describes where animals' activity centres are, not proof that %s causes higher density.",
                          lb, so_fmt(sc$scale), lb, so_fmt(abs(100 * (rr[1] - 1)), 2), if (rr[1] >= 1) "higher" else "lower",
                          so_fmt(rr[1]), so_fmt(rr[2]), so_fmt(rr[3]), so_fmt(xr[1], 3), so_fmt(xr[2], 3), lb,
                          so_fmt(stats::median(ratio)), so_fmt(so_q(ratio)[1]), so_fmt(so_q(ratio)[2]), mean(b > 0), lb))
  }
  out <- c(out,
    "Assumptions: closure within each primary period; activity centres fixed across periods (animals that shift their range between years bias survival low and can bias density; Glennie et al. 2019); survival and detection the same for all animals and periods (sex or age differences, trap-happiness and transients violate this; check the posterior predictive tests); detection a half-normal function of distance; recruits enter the state space between periods.",
    "Validation: see modules/scr-open/VALIDATION.md for simulation-recovery results (bias and interval coverage of every quantity above) at the default design and at harder designs, and the designs where the model fails.")
  out
}

# ---------------------------------------------------------------- figures
# core apply_fig_opts() adds the theme with `+`, which on a patchwork reaches
# only the last panel; theme every panel here before composing.
so_compose <- function(plots, opts, ncol = 2) {
  th <- theme_eco(opts$base_size, opts$font, opts$style, opts$legend)
  compose_figures(lapply(plots, function(p) p + th), ncol = ncol)
}

so_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  ob <- fit$object; prep <- ob$prep; TT <- prep$TT; lab <- prep$period_labels
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; a <- (1 - ci) / 2
  ici <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci; ai <- (1 - ici) / 2
  al <- so_area_label(prep$units); ul <- so_unit_len(prep)
  if (identical(prep$units, "m")) { kmu <- "km"; kmlab <- function(v) format(v / 1000, trim = TRUE, scientific = FALSE) } else { kmu <- ul; kmlab <- ggplot2::waiver() }
  pretty <- function(x) so_pretty_param(x, prep, data, opts)
  nm <- so_names(prep)
  dr <- so_draws(fit)
  pl <- fig_label("period", "Primary period", opts)
  tnum <- suppressWarnings(as.numeric(lab)); if (anyNA(tnum)) tnum <- seq_len(TT)
  band <- function(D, x) data.frame(time = x, estimate = colMeans(D), lower = apply(D, 2, stats::quantile, a),
                                    upper = apply(D, 2, stats::quantile, 1 - a),
                                    lower50 = apply(D, 2, stats::quantile, ai), upper50 = apply(D, 2, stats::quantile, 1 - ai))
  addin <- function(p, d) p + ggplot2::geom_ribbon(data = d, ggplot2::aes(time, ymin = lower50, ymax = upper50), inherit.aes = FALSE,
                                                    fill = eco_pal$band, alpha = 0.25)
  figs <- list()
  # abundance trajectory
  bN <- band(dr[, nm$N, drop = FALSE], tnum)
  obs <- if (isTRUE(opts$show_data)) data.frame(time = tnum, y = prep$n_by_period) else NULL
  p <- addin(fig_timeseries(bN, fig_label("abundance", "Abundance in the state space (N)", opts), obs = obs, xlab = pl), bN) +
    ggplot2::geom_point(data = bN, ggplot2::aes(time, estimate), colour = eco_pal$band, size = 1.4) +
    ggplot2::scale_x_continuous(breaks = tnum, labels = lab) + ggplot2::scale_y_continuous(limits = c(0, NA))
  figs$abundance_trajectory <- with_caption(p, sprintf(
    "Realised abundance per period (posterior mean, line) with %d%% and %d%% credible bands. Black points: number of individuals detected in each period (a minimum count). Density per %s = abundance / %s %s.",
    round(100 * ci), round(100 * ici), al, format(signif(prep$area * so_area_conv(prep), 3), big.mark = ","), al))
  # recruitment: recruits and per-capita rate
  bB <- band(dr[, nm$B, drop = FALSE], tnum[-1]); bf <- band(dr[, nm$f, drop = FALSE], tnum[-1])
  pB <- ggplot2::ggplot(bB, ggplot2::aes(time, estimate)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$ink, linewidth = 0.5) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = eco_pal$ink, linewidth = 1.3) +
    ggplot2::geom_point(size = 1.8, colour = eco_pal$ink) + ggplot2::scale_x_continuous(breaks = tnum[-1], labels = lab[-1]) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(x = pl, y = fig_label("recruits", "New recruits", opts))
  pf <- ggplot2::`%+%`(pB, bf) + ggplot2::labs(y = fig_label("per_capita_recruitment", "Per-capita recruitment", opts))
  figs$recruitment <- with_caption(so_compose(list(pB, pf), opts, ncol = 2),
    sprintf("(A) Number of new animals entering the state space in each period and (B) per-capita recruitment (recruits per animal alive in the previous period); points are posterior means, thick and thin bars %d%% and %d%% credible intervals.",
            round(100 * ici), round(100 * ci)))
  attr(figs$recruitment, "fig_size") <- list(width = "double", aspect = 0.4)
  # growth rate
  bl <- band(dr[, nm$lam, drop = FALSE], tnum[-1])
  figs$growth_rate <- with_caption(
    ggplot2::`%+%`(pB, bl) + ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null) +
      ggplot2::scale_y_continuous() + ggplot2::labs(y = fig_label("growth_rate", "Growth rate\n(N this year / N last year)", opts)),
    sprintf("Realised population growth rate between consecutive periods (posterior mean, %d%% and %d%% credible intervals); dashed line: stable population.", round(100 * ici), round(100 * ci)))
  # survival
  phd <- dr[, "survival"]
  figs$survival <- with_caption(
    ggplot2::ggplot(data.frame(v = phd), ggplot2::aes(v)) +
      ggplot2::geom_density(fill = eco_pal$band, alpha = 0.3, colour = eco_pal$band) +
      ggplot2::geom_vline(xintercept = so_q(phd, c(a, 1 - a)), linetype = 2, colour = eco_pal$ink) +
      ggplot2::geom_hline(yintercept = 0, colour = eco_pal$muted, linewidth = 0.3) +
      ggplot2::scale_x_continuous(breaks = so_inner_breaks, expand = ggplot2::expansion(mult = 0.06)) +
      ggplot2::labs(x = fig_label("survival", "Apparent survival between periods", opts), y = "Posterior density"),
    sprintf("Posterior distribution of apparent survival between consecutive periods (constant over time; Uniform(0,1) prior); dashed lines: %d%% credible interval.", round(100 * ci)))
  # density surface per period
  surf <- ob$surfaces
  if (is.null(surf)) surf <- tryCatch(so_surfaces(prep, so_pool(ob$chains, 100), prep$M), error = function(e) NULL)
  if (!is.null(surf)) {
    cellA <- prep$pix^2 * so_area_conv(prep)
    gd <- do.call(rbind, lapply(seq_len(TT), function(t) data.frame(x = prep$mask$x, y = prep$mask$y, value = surf$realised[, t] / cellA,
                                                                      period = factor(lab[t], levels = lab))))
    p <- fig_map(gd, sprintf("Animals per %s", al), points = data.frame(x = prep$traps$x, y = prep$traps$y)) +
      ggplot2::facet_wrap(~period, labeller = ggplot2::labeller(period = function(v) paste(pl, v))) +
      ggplot2::scale_x_continuous(n.breaks = 3, labels = kmlab) + ggplot2::scale_y_continuous(n.breaks = 4, labels = kmlab) +
      ggplot2::labs(x = fig_label("easting", sprintf("Easting (%s)", kmu), opts), y = fig_label("northing", sprintf("Northing (%s)", kmu), opts))
    figs$density_surface <- with_caption(p, sprintf(
      "Posterior mean realised density (animals per %s) in each period: alive animals placed at their estimated activity centres, including animals never detected. White circles: traps (%d posterior draws).", al, surf$R))
    attr(figs$density_surface, "fig_size") <- list(width = "double", aspect = if (TT > 3) 0.75 else 0.45)
  }
  # density covariate effects (response scale, per period)
  for (k in seq_along(prep$dcov)) {
    v <- prep$dcov[k]; sc <- prep$dscale[[v]]
    xs <- seq(min(prep$mask[[v]]), max(prep$mask[[v]]), length.out = 50)
    PD <- length(prep$dcov)
    B <- dr[, paste0("density_beta_", prep$dcov), drop = FALSE]
    lnorm <- apply(B, 1, function(b) so_lse(as.vector(prep$XD %*% b)))
    sel <- unique(c(1, TT))
    rows <- lapply(sel, function(t) {
      Dm <- sapply(xs, function(x) dr[, nm$N[t]] * exp(B[, k] * (x - sc$center) / sc$scale - lnorm) / (prep$pix^2 * so_area_conv(prep)))
      r <- data.frame(x = xs, estimate = colMeans(Dm), lower = apply(Dm, 2, stats::quantile, a), upper = apply(Dm, 2, stats::quantile, 1 - a))
      r$group <- paste(pl, lab[t]); r
    })
    pred <- do.call(rbind, rows); pred$group <- factor(pred$group, levels = unique(pred$group))
    obsx <- if (isTRUE(opts$show_data)) data.frame(x = prep$mask[[v]][so_thin_idx(prep$G, 300)]) else NULL
    figs[[paste0("density_", v)]] <- with_caption(
      fig_cond_effect(pred, so_lab(v, data, opts), fig_label("density", sprintf("Density (animals per %s)", al), opts), obs = obsx) +
        ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL) + ggplot2::scale_y_continuous(limits = c(0, NA)),
      sprintf("Expected density versus %s in the first and last period (other covariates at their state-space mean), with %d%% credible bands. Rug: covariate values across the state space.",
              so_lab(v, data, opts), round(100 * ci)))
  }
  # detection function
  sgd <- dr[, "sigma"]; p0d <- dr[, "baseline_detection"]
  dd <- seq(0, 4 * stats::median(sgd), length.out = 80)
  P <- sapply(dd, function(x) p0d * exp(-x^2 / (2 * sgd^2)))
  pred <- data.frame(x = dd, estimate = colMeans(P), lower = apply(P, 2, stats::quantile, a), upper = apply(P, 2, stats::quantile, 1 - a),
                     lower50 = apply(P, 2, stats::quantile, ai), upper50 = apply(P, 2, stats::quantile, 1 - ai))
  sp <- prep$trap_spacing
  figs$detection_function <- with_caption(
    fig_cond_effect(pred, fig_label("distance", sprintf("Distance from activity centre to trap (%s)", ul), opts),
                    fig_label("detection_prob", "Detection probability", opts)) +
      ggplot2::geom_vline(xintercept = sp, linetype = 2, colour = "grey45") +
      ggplot2::annotate("text", x = sp, y = Inf, label = " trap spacing", hjust = 0, vjust = 1.3, size = 2.6, colour = "grey30") +
      ggplot2::scale_y_continuous(limits = c(0, NA)),
    sprintf("Probability of detecting an animal at a trap on one occasion versus distance from its activity centre: half-normal p0 exp(-d^2 / 2 sigma^2), constant across periods, with %d%% and %d%% credible bands. Dashed line: trap spacing.", round(100 * ci), round(100 * ici)))
  # detection by period + gof
  ppc <- ob$ppc
  if (!is.null(ppc)) {
    pn <- paste0("detected_period_", seq_len(TT))
    R <- ppc$rep[, pn, drop = FALSE]
    dfp <- data.frame(time = tnum, obs = as.numeric(ppc$obs[pn]), est = colMeans(R),
                      lo = apply(R, 2, stats::quantile, a), hi = apply(R, 2, stats::quantile, 1 - a))
    figs$detection_by_period <- with_caption(
      ggplot2::ggplot(dfp, ggplot2::aes(time)) + ggplot2::geom_col(ggplot2::aes(y = obs), fill = "grey75", width = 0.6) +
        ggplot2::geom_pointrange(ggplot2::aes(y = est, ymin = lo, ymax = hi), colour = eco_pal$ink, size = 0.3) +
        ggplot2::scale_x_continuous(breaks = tnum, labels = lab) +
        ggplot2::labs(x = pl, y = "Individuals detected"),
      sprintf("Bars: individuals detected in each period. Points and bars: posterior predictive mean and %d%% interval under the fitted model (constant p0 and sigma). Bars outside the intervals suggest detection varied between periods.", round(100 * ci)))
    st <- c(n_detected = "Individuals detected (all periods)", detected_in_2plus_periods = "Detected in 2+ periods",
            total_detections = "Total detections", traps_per_individual_period = "Traps per individual-period",
            gaps_in_histories = "Missed periods within histories")
    pls <- lapply(names(st), function(s) fig_ppc(ppc$rep[, s], ppc$obs[[s]], st[[s]]))
    figs$gof <- with_caption(so_compose(pls, opts, ncol = 2),
      "Posterior predictive checks: histograms of statistics computed from data simulated under the fitted model; line = observed value. 'Detected in 2+ periods' and 'missed periods' react to transients and temporary emigration; 'traps per individual-period' to the detection function.")
    attr(figs$gof, "fig_size") <- list(width = "double", aspect = 1)
  }
  ec <- fit$estimates[grepl("^density_beta_", fit$estimates$param), ]
  if (nrow(ec)) figs$coefficients <- with_caption(fig_coef(ec, null_value = 0, xlab = "Log rate ratio of density per SD (95% CrI)", opts = opts),
    "Density covariate coefficients on the log scale per 1 SD of the covariate across the state space; dashed line: no association.")
  if (!is.null(ob$comparison)) {
    cm <- ob$comparison
    keep <- intersect(c(nm$N, "survival", nm$f, "sigma", "baseline_detection", "superpopulation_size", paste0("density_beta_", prep$dcov)), unique(cm$param))
    cm <- cm[cm$param %in% keep, ]
    cm$param <- factor(pretty(as.character(cm$param)), levels = pretty(keep))
    figs$engine_comparison <- with_caption(
      ggplot2::ggplot(cm, ggplot2::aes(estimate, engine, colour = engine)) +
        ggplot2::geom_point(size = 1.2) +
        ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), na.rm = TRUE) +
        ggplot2::facet_wrap(~param, scales = "free_x", ncol = 2, labeller = ggplot2::label_wrap_gen(32)) + ggplot2::labs(x = "Estimate (95% interval)", y = NULL) +
        ggplot2::scale_x_continuous(breaks = so_inner_breaks, expand = ggplot2::expansion(mult = 0.12), guide = ggplot2::guide_axis(check.overlap = TRUE)) + ggplot2::guides(colour = "none") +
        ggplot2::scale_y_discrete(labels = function(x) vapply(strwrap(x, 24, simplify = FALSE), paste, "", collapse = "\n")) +
        ggplot2::theme(panel.spacing.x = ggplot2::unit(3, "lines"), plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5)),
      paste("The same data fitted with each engine. openCR is maximum likelihood with a hazard half-normal detection function (its abundance is a point estimate from derived(); intervals not shown).",
            if (length(prep$dcov) && any(grepl("openCR", cm$engine))) "openCR is fitted with homogeneous density, whereas the Bayesian fits include the habitat covariate: its abundance differs whenever the traps sample habitat of unrepresentative density (see the design notes); survival, recruitment and sigma remain comparable." else ""))
    attr(figs$engine_comparison, "fig_size") <- list(width = "double", aspect = 1.25)
  }
  core <- c("survival", "baseline_detection", "sigma", nm$N[1], nm$N[TT], "superpopulation_size", paste0("density_beta_", prep$dcov))
  relab <- function(m) { m <- m[, core, drop = FALSE]; colnames(m) <- pretty(core); m }
  figs$trace <- with_caption(fig_trace(lapply(ob$rep_chains, relab)),
                             "MCMC traces of key quantities, one colour per chain (thinned draws used for derived quantities).")
  attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.7)
  figs$posterior <- with_caption(fig_posterior(relab(dr)) +
                                   ggplot2::facet_wrap(~param, scales = "free", ncol = 3, labeller = ggplot2::label_wrap_gen(28)) +
                                   ggplot2::scale_x_continuous(n.breaks = 4, guide = ggplot2::guide_axis(check.overlap = TRUE)),
                                 "Posterior densities of key quantities.")
  attr(figs$posterior, "fig_size") <- list(width = "double", aspect = 0.7)
  figs
}

# Short readable labels for facet strips (engine comparison, trace, posterior).
so_pretty_param <- function(x, prep, data, opts = NULL) {
  lab <- prep$period_labels; ul <- so_unit_len(prep)
  out <- x
  for (t in seq_along(lab)) {
    out[x == paste0("abundance_period_", lab[t])] <- paste("Abundance", lab[t])
    out[x == paste0("per_capita_recruitment_period_", lab[t])] <- paste("Per-capita recruitment", lab[t])
  }
  out[x == "survival"] <- "Apparent survival"
  out[x == "baseline_detection"] <- "Baseline detection (p0)"
  out[x == "sigma"] <- sprintf("Sigma (%s)", ul)
  out[x == "superpopulation_size"] <- "Superpopulation"
  for (v in prep$dcov) out[x == paste0("density_beta_", v)] <- paste("Density effect:", so_lab(v, data, opts))
  out
}

# ---------------------------------------------------------------- methods
so_methods_text <- function(fit, data, settings) {
  ob <- fit$object; prep <- ob$prep; st <- ob$settings
  al <- so_area_label(prep$units)
  dens <- if (length(prep$dcov)) sprintf("with pixel probabilities proportional to exp(beta x), x = %s standardised to mean 0 and SD 1 over the state space", paste(prep$dcov, collapse = ", ")) else "uniformly over the state space"
  eng <- switch(ob$engine,
    "nimble marginal" = sprintf("Activity centres and alive states were integrated out of the likelihood (alive states with the forward algorithm, centres by summing over mask pixels; Turek et al. 2016; Milleret et al. 2019) and the %d remaining parameters were sampled jointly with an adaptive block random-walk sampler in NIMBLE %s (de Valpine et al. 2017): %d chains of %d iterations after %d burn-in, thinned by %d.",
                                3 + prep$TT + length(prep$dcov), as.character(utils::packageVersion("nimble")), st$nchains, st$niter - st$nburnin, st$nburnin, st$thin),
    "nimble hmm" = sprintf("Alive states were integrated out with the forward algorithm per individual and activity centres sampled with block random-walk samplers in NIMBLE %s: %d chains of %d iterations after %d burn-in, thinned by %d.",
                           as.character(utils::packageVersion("nimble")), st$nchains, st$niter_latent - st$nburnin_latent, st$nburnin_latent, st$thin_latent),
    "nimble latent" = sprintf("Alive states and activity centres were sampled explicitly in NIMBLE %s: %d chains of %d iterations after %d burn-in, thinned by %d.",
                              as.character(utils::packageVersion("nimble")), st$nchains, st$niter_latent - st$nburnin_latent, st$nburnin_latent, st$thin_latent),
    sprintf("Activity centres and alive states were integrated out and the parameters sampled with the No-U-Turn sampler in NumPyro (Phan et al. 2019): %d chains of %d draws after %d warm-up.",
            st$nchains, st$numpyro_samples, st$numpyro_warmup))
  sprintf(paste(
    "We fitted an open-population (Jolly-Seber) spatial capture-recapture model (Gardner et al. 2010; Royle et al. 2014; Efford & Schofield 2020) to detections of %d individuals over %d primary periods of %d secondary occasions at %d traps.",
    "The population was assumed closed within periods and open between them. Using data augmentation (Royle & Dorazio 2012; M = %d potential individuals), each individual is not yet entered, alive or dead in each period: it enters in period 1 with probability gamma_1 and in later periods with probability gamma_t if not yet entered, and an alive individual survives to the next period with apparent survival phi (constant).",
    "Each individual has a static activity centre in a state space of %d mask pixels (%s %s spacing, %s %s, extending %s %s beyond the traps), distributed %s.",
    "While alive, the number of occasions on which it was detected at trap j in period t was Binomial(K_jt, p0 exp(-d^2 / (2 sigma^2))), where K_jt is the number of occasions the trap operated.",
    "Priors: p0 and phi Uniform(0,1); entry in the superpopulation parameterisation (Kery & Schaub 2012): the share psi of the M rows ever alive Uniform(0,1), so the superpopulation has a uniform prior on 0..M, and the entry probabilities of those animals over periods Dirichlet(1,...,1), from which gamma_t is derived; sigma Uniform(0, %s %s); density coefficients Normal(0, SD 2). %s",
    "Realised abundance N_t, recruits B_t, the superpopulation, per-capita recruitment B_t / N_t-1 and growth N_t / N_t-1 were obtained by drawing each individual's activity-centre pixel and state path from their conditional posterior (forward-filtering backward-sampling) for %d thinned posterior draws per chain; density is N_t divided by the state-space area.",
    "Convergence was assessed with rank-normalised R-hat and bulk effective sample size (posterior package); fit was checked with posterior predictive tests of the numbers of individuals detected per period and in two or more periods, total detections, traps per individual-period and missed periods within capture histories."),
    prep$n, prep$TT, prep$K, prep$J, prep$M, prep$G, format(signif(prep$spacing, 3)), so_unit_len(prep),
    so_fmt(prep$area * so_area_conv(prep), 3), al, format(signif(prep$buffer, 3)), so_unit_len(prep), dens,
    format(signif(prep$sigma_max * prep$cscale, 3)), so_unit_len(prep), eng, st$max_draws_per_chain)
}
