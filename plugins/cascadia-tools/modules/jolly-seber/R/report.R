# jolly-seber: interpretation, figures, methods text, references.

js_references <- c(
  "Jolly GM (1965) Explicit estimates from capture-recapture data with both death and immigration-stochastic model. Biometrika 52:225-247. doi:10.1093/biomet/52.1-2.225",
  "Seber GAF (1965) A note on the multiple-recapture census. Biometrika 52:249-259. doi:10.2307/2333827",
  "Crosbie SF, Manly BFJ (1985) Parsimonious modelling of capture-mark-recapture studies. Biometrics 41:385-398. doi:10.2307/2530864",
  "Schwarz CJ, Arnason AN (1996) A general methodology for the analysis of capture-recapture experiments in open populations. Biometrics 52:860-873. doi:10.2307/2533048",
  "Pradel R (1996) Utilization of capture-mark-recapture for the study of recruitment and population growth rate. Biometrics 52:703-709. doi:10.2307/2532908",
  "Kendall WL, Nichols JD, Hines JE (1997) Estimating temporary emigration using capture-recapture data with Pollock's robust design. Ecology 78:563-578. doi:10.2307/2266030",
  "Royle JA, Dorazio RM (2008) Hierarchical Modeling and Inference in Ecology. Academic Press (ch. 10, open populations / Jolly-Seber as restricted dynamic occupancy).",
  "Royle JA, Dorazio RM (2012) Parameter-expanded data augmentation for Bayesian analysis of capture-recapture models. Journal of Ornithology 152(S2):521-537. doi:10.1007/s10336-010-0619-4",
  "Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS: a hierarchical perspective. Academic Press (ch. 10, Jolly-Seber model).",
  "Dorazio RM (2020) Objective prior distributions for Jolly-Seber models of zero-augmented data. Biometrics 76:1285-1296. doi:10.1111/biom.13216",
  "Yackulic CB, Dodrill M, Dzul M, Sanderlin JS, Reid JA (2020) A need for speed in Bayesian population models: a practical guide to marginalizing and recovering discrete latent states. Ecological Applications 30:e02112. doi:10.1002/eap.2112",
  "Saracco JF, Yackulic CB (2023) Fast and flexible Bayesian Jolly Seber models and application to populations with transients. bioRxiv doi:10.1101/2023.05.25.542303",
  "Laake JL, Johnson DS, Conn PB (2013) marked: an R package for maximum likelihood and Markov chain Monte Carlo analysis of capture-recapture data. Methods in Ecology and Evolution 4:885-890. doi:10.1111/2041-210X.12065",
  "Gimenez O, Lebreton J-D, Choquet R, Pradel R (2018) R2ucare: an R package to perform goodness-of-fit tests for capture-recapture models. Methods in Ecology and Evolution 9:1749-1754. doi:10.1111/2041-210X.13014",
  "de Valpine P et al. (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413. doi:10.1080/10618600.2016.1172487")

js_row <- function(est, p) est[est$param == p, , drop = FALSE]
js_fmt <- function(r, d = 2, pct = FALSE) {
  if (!nrow(r)) return("NA")
  f <- if (pct) function(x) sprintf(paste0("%.", d, "f%%"), 100 * x) else function(x) formatC(x, digits = d, format = "f")
  sprintf("%s (95%% CI %s to %s)", f(r$estimate), f(r$lower), f(r$upper))
}
js_prob <- function(x) if (x > 0.995) "> 0.99" else if (x < 0.005) "< 0.01" else sprintf("= %.2f", x)
js_int <- function(r) if (!nrow(r)) "NA" else sprintf("%.0f (95%% CI %.0f to %.0f)", r$estimate, r$lower, r$upper)

js_draw_matrix <- function(fit) {
  o <- fit$object
  if (!is.null(o$rep_chains)) do.call(rbind, o$rep_chains) else NULL
}

js_interpret <- function(fit, data) {
  est <- fit$estimates; o <- fit$object; prep <- o$prep
  T <- prep$T; lab <- prep$labels; tu <- prep$time_unit
  sp <- if (!is.null(data$species)) data$species else "animals"
  ch <- prep$ch
  out <- character()
  out <- c(out, sprintf("%d %s were caught at least once over %d occasions (%s to %s); %d (%.0f%%) were caught on more than one occasion, and %.1f were caught per occasion on average.",
                        prep$n, sp, T, lab[1], lab[T], sum(rowSums(ch) > 1), 100 * mean(rowSums(ch) > 1), mean(colSums(ch))))
  N1 <- js_row(est, "abundance_1"); NT <- js_row(est, paste0("abundance_", T))
  Ns <- js_row(est, "superpopulation_size"); gm <- js_row(est, "growth_rate_mean")
  D <- js_draw_matrix(fit)
  pdec <- if (!is.null(D)) mean(D[, "growth_rate_mean"] < 1, na.rm = TRUE) else NA
  out <- c(out, sprintf("Estimated abundance was %s at the first occasion (%s) and %s at the last (%s). Across the study the average growth rate was %s per %s, i.e. a change of %+.1f%% per %s%s.",
                        js_int(N1), lab[1], js_int(NT), lab[T], js_fmt(gm, 3), tu, 100 * (gm$estimate - 1), tu,
                        if (is.finite(pdec)) sprintf("; the posterior probability that the population declined (growth rate < 1) %s", sub("^= ", "is ", js_prob(pdec))) else ""))
  out <- c(out, sprintf("In total %s %s are estimated to have used the study population at some time during the study (the superpopulation), of which %d (%.0f%%) were ever caught.",
                        js_int(Ns), sp, prep$n, 100 * prep$n / Ns$estimate))
  if (prep$G == 2) {
    l2 <- prep$group_levels[2]; pr2 <- js_row(est, paste0("prop_", l2))
    out <- c(out, sprintf("Animals with %s = %s made up %s of that superpopulation (%s: %d of the %d animals caught, %d of unknown %s).",
                          prep$group_name, l2, js_fmt(pr2, 1, pct = TRUE), "captured animals",
                          sum(prep$gcode == 2), prep$n, sum(prep$gcode == 0), prep$group_name))
  }
  # survival
  if (prep$Kphi_lev == 1) {
    s1 <- js_row(est, "survival")
    if (prep$G == 2) {
      l <- prep$group_levels; s2 <- js_row(est, paste0("survival_", l[2])); bg <- js_row(est, paste0("survival_beta_", l[2]))
      pneg <- if (!is.null(D)) mean(D[, paste0("survival_beta_", l[2])] < 0) else NA
      out <- c(out, sprintf("Apparent survival per %s was %s for %s and %s for %s. The odds of surviving were %.2f times as high for %s as for %s (95%% CI %.2f to %.2f)%s. 'Apparent' survival means that death and permanent emigration cannot be told apart.",
                            tu, js_fmt(s1), l[1], js_fmt(s2), l[2], exp(bg$estimate), l[2], l[1], exp(bg$lower), exp(bg$upper),
                            if (is.finite(pneg)) sprintf("; posterior probability that %s animals survive less well than %s animals %s", l[2], l[1], js_prob(pneg)) else ""))
    } else out <- c(out, sprintf("Apparent survival per %s was %s; death and permanent emigration cannot be told apart.", tu, js_fmt(s1)))
  } else {
    sv <- est[grepl("^survival_[0-9]+$", est$param), ]
    out <- c(out, sprintf("Apparent survival per %s varied between intervals from %.2f to %.2f (lowest %s-%s, highest %s-%s)%s.", tu, min(sv$estimate), max(sv$estimate),
                          lab[which.min(sv$estimate)], lab[which.min(sv$estimate) + 1], lab[which.max(sv$estimate)], lab[which.max(sv$estimate) + 1],
                          if (prep$Kp_lev == T) "; the last interval's survival is confounded with the last capture probability" else ""))
  }
  # recruitment
  fr <- est[grepl("^recruitment_rate_", est$param), ]
  if (!is.null(D)) {
    fm <- rowMeans(D[, grepl("^recruitment_rate_", colnames(D)), drop = FALSE], na.rm = TRUE)
    out <- c(out, sprintf("Per-capita recruitment (new animals entering per animal present, births plus immigration) averaged %.2f (95%% CI %.2f to %.2f) per interval, ranging from %.2f (%s-%s) to %.2f (%s-%s). Entrants at the first occasion are not recruitment: the first entry probability (%.2f) is the share of the superpopulation already present when sampling started.",
                          mean(fm), stats::quantile(fm, .025), stats::quantile(fm, .975), min(fr$estimate), lab[which.min(fr$estimate)], lab[which.min(fr$estimate) + 1],
                          max(fr$estimate), lab[which.max(fr$estimate)], lab[which.max(fr$estimate) + 1], js_row(est, "entry_prob_1")$estimate))
  }
  # detection
  if (prep$Kp_lev == T) {
    pv <- est[grepl("^detection_[0-9]+$", est$param), ]
    out <- c(out, sprintf("Capture probability varied by occasion from %.2f to %.2f.", min(pv$estimate), max(pv$estimate)))
  } else {
    p0 <- js_row(est, "detection")
    s <- sprintf("Capture probability of an animal present was %s%s.", js_fmt(p0), if (length(prep$scaling)) " at average sampling effort" else "")
    for (v in names(prep$scaling)) {
      bb <- js_row(est, paste0("detection_beta_", v)); sc <- prep$scaling[[v]]
      lo <- stats::plogis(stats::qlogis(p0$estimate) - bb$estimate); hi <- stats::plogis(stats::qlogis(p0$estimate) + bb$estimate)
      vl <- if (!is.null(data$covariate_labels[[v]])) tolower(data$covariate_labels[[v]]) else v
      s <- c(s, sprintf("Capture probability %s from %.2f to %.2f between %s = %.3g and %.3g (one SD below to one SD above the mean; odds ratio per SD %.2f, 95%% CI %.2f to %.2f). This is a sampling (detection) effect, not a change in the population.",
                        if (hi >= lo) "rose" else "fell", lo, hi, vl, sc$center - sc$scale, sc$center + sc$scale, exp(bb$estimate), exp(bb$lower), exp(bb$upper)))
    }
    out <- c(out, paste(s, collapse = " "))
  }
  if (!is.null(o$ppc)) out <- c(out, sprintf("Posterior predictive checks: Bayesian p-values %.2f (captures per animal) and %.2f (animals caught per occasion); values near 0 or 1 (below 0.025 or above 0.975) would indicate lack of fit, e.g. capture heterogeneity.",
                                             o$ppc$p[["individual"]], o$ppc$p[["occasion"]]))
  out <- c(out, "Assumptions to keep in mind: all animals present at an occasion (marked or not) share the same capture probability; marks are not lost or misread; animals that leave do so permanently (temporary emigration with a tendency to stay away biases survival and abundance; see the module's validation); and the population is closed within each occasion. Abundance refers to the animals that use the sampled area, not a density.")
  if (!isTRUE(fit$diagnostics$converged)) out <- c("WARNING: the model did not converge; do not interpret these numbers.", out)
  out
}

js_methods_text <- function(fit, data, settings) {
  o <- fit$object; prep <- o$prep; st <- o$settings
  T <- prep$T
  det <- if (prep$Kp_lev == T) "a separate capture probability for each occasion" else if (length(prep$scaling))
    sprintf("capture probability logit-linear in %s (standardised)", paste(names(prep$scaling), collapse = ", ")) else "constant capture probability"
  sv <- if (prep$Kphi_lev > 1) "apparent survival varying by interval" else "constant apparent survival per time unit"
  if (prep$G == 2) sv <- paste0(sv, sprintf(" with a logit-scale effect of %s (%s vs %s)", prep$group_name, prep$group_levels[2], prep$group_levels[1]))
  if (grepl("marked", o$engine)) {
    return(sprintf("We fitted a Jolly-Seber model in the POPAN parameterisation of Schwarz & Arnason (1996) by maximum likelihood with the R package marked (Laake et al. 2013), with %s, %s, occasion-specific entry probabilities and a superpopulation size; survival was raised to the power of the interval length. Abundance per occasion was derived as its expectation and 95%% intervals were obtained from 2000 draws of the coefficients from their asymptotic normal distribution.", sv, det))
  }
  sprintf(paste(
    "We estimated abundance, survival and recruitment with a Bayesian Jolly-Seber model (Jolly 1965; Seber 1965; Schwarz & Arnason 1996) fitted by data augmentation in the restricted dynamic occupancy formulation (Royle & Dorazio 2008, 2012; Kery & Schaub 2012).",
    "The %d animals caught were augmented with %d all-zero capture histories (M = %d); each potential animal enters the population at occasion t with hazard gamma_t, survives each interval with probability phi^dt (dt = interval length in %ss) and is caught with probability p_t when alive.",
    "Entry hazards were derived from the superpopulation inclusion probability psi ~ Uniform(0,1) and the entry probabilities b ~ Dirichlet(1,...,1), which makes the prior on the superpopulation size uniform on 0..M and the prior on the entry probabilities uniform.",
    "The model had %s and %s. Intercepts had Logistic(0,1) priors (uniform on the probability scale) and slopes Normal(0, sd = %g) priors on the logit scale%s.",
    "%s",
    "We ran %d chains of %d iterations (burn-in %d, thinning %d) in NIMBLE %s (de Valpine et al. 2017) and checked convergence with R-hat and bulk effective sample size (posterior package).",
    "We confirmed that the posterior of the superpopulation size was not truncated by M (99.5%% quantile below 0.9 M).",
    "Realised abundance N_t, entrants B_t, per-capita recruitment B_(t+1)/N_t and realised growth rates N_(t+1)/N_t are derived quantities. Goodness of fit was assessed with posterior predictive checks (Freeman-Tukey discrepancies of captures per animal and per occasion) and the R2ucare tests 3.SR and 2.CT (Gimenez et al. 2018)."),
    prep$n, prep$M - prep$n, prep$M, prep$time_unit, sv, det, st$prior_sd_slope,
    if (prep$G == 2) sprintf("; the proportion of %s had a Uniform(0,1) prior, animals of unknown %s were treated as a mixture, and the reported share of %s is the realised share among the animals in the superpopulation", prep$group_levels[2], prep$group_name, prep$group_levels[2]) else "",
    if (grepl("marginal", o$engine)) "The latent alive states were summed out of the likelihood with the forward algorithm over three states (not yet entered, alive, dead) for each unique capture history (Saracco & Yackulic 2023), and the latent states were drawn exactly from their conditional posterior for every retained draw." else "The latent alive states z[i,t] were sampled directly.",
    st$nchains, if (grepl("latent", o$engine)) st$niter_latent else st$niter, if (grepl("latent", o$engine)) st$nburnin_latent else st$nburnin,
    if (grepl("latent", o$engine)) st$thin_latent else st$thin, as.character(utils::packageVersion("nimble")))
}

# ---- figures --------------------------------------------------------------

js_time_axis <- function(prep) {
  num <- suppressWarnings(as.numeric(prep$labels))
  if (!anyNA(num)) list(x = num, numeric = TRUE) else list(x = seq_len(prep$T), numeric = FALSE)
}

js_qband <- function(D, cols, ci, inner) {
  a <- (1 - ci) / 2; b <- (1 - inner) / 2
  data.frame(estimate = colMeans(D[, cols, drop = FALSE], na.rm = TRUE),
             lower = apply(D[, cols, drop = FALSE], 2, stats::quantile, a, na.rm = TRUE),
             upper = apply(D[, cols, drop = FALSE], 2, stats::quantile, 1 - a, na.rm = TRUE),
             lower50 = apply(D[, cols, drop = FALSE], 2, stats::quantile, b, na.rm = TRUE),
             upper50 = apply(D[, cols, drop = FALSE], 2, stats::quantile, 1 - b, na.rm = TRUE))
}

js_series <- function(fit, D, cols, ci, inner) {
  if (!is.null(D)) return(js_qband(D, cols, ci, inner))
  est <- fit$estimates
  r <- est[match(cols, est$param), ]
  data.frame(estimate = r$estimate, lower = r$lower, upper = r$upper)
}

js_pointrange <- function(d, xlab, ylab, ax, null = NULL) {
  p <- ggplot2::ggplot(d, ggplot2::aes(x, estimate))
  if (!is.null(null)) p <- p + ggplot2::geom_hline(yintercept = null, linetype = 2, colour = eco_pal$null, linewidth = 0.4)
  p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$ink, linewidth = 0.5)
  if ("lower50" %in% names(d)) p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = eco_pal$ink, linewidth = 1.3)
  p <- p + ggplot2::geom_point(size = 1.8, colour = eco_pal$ink) + ggplot2::labs(x = xlab, y = ylab)
  if (!ax$numeric) p <- p + ggplot2::scale_x_continuous(breaks = d$x, labels = d$label)
  p
}

# Theme every panel before composing: the core adds the theme with `+`,
# which on a patchwork only reaches the last panel.
js_compose <- function(plots, opts, ncol = NULL) {
  th <- theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic", opts$legend %||% "right")
  compose_figures(lapply(plots, function(q) q + th), ncol = ncol)
}

js_figures <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2", "patchwork"))
  o <- fit$object; prep <- o$prep; est <- fit$estimates
  T <- prep$T; ci <- if (is.null(opts$ci)) 0.95 else opts$ci; inner <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci
  bayes <- !is.null(o$rep_chains)
  D <- js_draw_matrix(fit)
  ax <- js_time_axis(prep)
  cil <- sprintf("%d%%", round(100 * ci)); itype <- if (bayes) "credible" else "confidence"
  inner_txt <- if (bayes) sprintf(" (inner band %d%%)", round(100 * inner)) else ""
  tu <- prep$time_unit
  xlab <- fig_label("occasion", if (ax$numeric) tools::toTitleCase(tu) else "Occasion", opts)
  sp <- fig_label("animals", "animals", opts)
  figs <- list()
  # abundance trajectory
  ab <- js_series(fit, D, paste0("abundance_", seq_len(T)), ci, inner)
  ab$time <- ax$x
  obs <- if (isTRUE(opts$show_data)) data.frame(time = ax$x, y = colSums(prep$ch)) else NULL
  p <- fig_timeseries(ab, fig_label("abundance", "Abundance (N)", opts), obs, xlab = xlab)
  if (bayes) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25)
  p <- p + ggplot2::expand_limits(y = 0) + ggplot2::geom_point(size = 1.3, colour = eco_pal$band)
  if (!ax$numeric) p <- p + ggplot2::scale_x_continuous(breaks = ax$x, labels = prep$labels)
  figs$abundance_trajectory <- with_caption(p, sprintf(
    "Estimated abundance at each occasion (line and points: %s; band: %s %s interval%s). Black points: number of %s actually caught at each occasion (a minimum count).",
    if (bayes) "posterior mean" else "expected abundance (MLE)", cil, itype, inner_txt, sp))
  # recruitment: entrants and per-capita rate
  en <- js_series(fit, D, paste0("entrants_", 2:T), ci, inner); en$x <- ax$x[-1]; en$label <- prep$labels[-1]
  fr <- js_series(fit, D, paste0("recruitment_rate_", seq_len(T - 1)), ci, inner)
  fr$x <- (ax$x[-1] + ax$x[-T]) / 2; fr$label <- paste(prep$labels[-T], prep$labels[-1], sep = "-")
  pa <- js_pointrange(en, xlab, fig_label("entrants", "Entrants (animals)", opts), ax)
  pb <- js_pointrange(fr, xlab, fig_label("recruitment_rate", "Per-capita recruitment rate", opts), ax) + ggplot2::expand_limits(y = 0)
  if (!ax$numeric) pb <- pb + ggplot2::scale_x_continuous(breaks = fr$x, labels = fr$label)
  figs$recruitment <- with_caption(js_compose(list(pa, pb), opts, ncol = 2),
    sprintf("(A) Number of animals entering the population (births plus immigrants) that were first present at each occasion; (B) per-capita recruitment, entrants at t+1 divided by abundance at t, plotted at the midpoint of each interval. Points: %s; bars: %s %s interval%s.",
            if (bayes) "posterior means" else "MLE", cil, itype, if (bayes) " (thick bars 50%)" else ""))
  attr(figs$recruitment, "fig_size") <- list(width = "double", aspect = 0.42)
  # entry probabilities
  eb <- js_series(fit, D, paste0("entry_prob_", seq_len(T)), ci, inner); eb$x <- ax$x; eb$label <- prep$labels
  figs$entry_probabilities <- with_caption(js_pointrange(eb, xlab, fig_label("entry_prob", "Entry probability", opts), ax) + ggplot2::expand_limits(y = 0),
    sprintf("Entry probabilities b_t: the share of all animals that ever used the population (the superpopulation) that entered at each occasion; the first value is the share already present at the start. Bars: %s %s intervals.", cil, itype))
  # growth rates
  gr <- js_series(fit, D, paste0("growth_rate_", seq_len(T - 1)), ci, inner); gr$x <- fr$x; gr$label <- fr$label
  pg <- js_pointrange(gr, xlab, fig_label("growth_rate", "Growth rate N(t+1) / N(t)", opts), ax, null = 1)
  gmr <- js_row(est, "growth_rate_mean")
  pg <- pg + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = gmr$lower, ymax = gmr$upper, alpha = 0.12, fill = eco_pal$band) +
    ggplot2::geom_hline(yintercept = gmr$estimate, colour = eco_pal$band, linewidth = 0.6) +
    # redraw the no-change line on top so the mean line cannot hide it
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null, linewidth = 0.4)
  if (!ax$numeric) pg <- pg + ggplot2::scale_x_continuous(breaks = gr$x, labels = gr$label)
  figs$growth_rate <- with_caption(pg, sprintf(
    "Realised population growth rate between consecutive occasions (points, %s %s intervals); dashed line = no change (1). Green line and shading: mean growth rate per %s over the whole study with its %s interval.%s",
    cil, itype, tu, cil, if (any(prep$dt != 1)) " Rates between occasions further apart than one time unit cover the whole gap." else ""))
  # survival
  if (prep$Kphi_lev > 1) {
    sv <- js_series(fit, D, paste0("survival_", seq_len(T - 1)), ci, inner); sv$x <- fr$x; sv$label <- fr$label
    p <- js_pointrange(sv, xlab, fig_label("survival", sprintf("Apparent survival per %s", tu), opts), ax) + ggplot2::coord_cartesian(ylim = c(0, 1))
    if (!ax$numeric) p <- p + ggplot2::scale_x_continuous(breaks = sv$x, labels = sv$label)
    figs$survival_time <- with_caption(p, sprintf("Apparent survival per %s for each interval (%s %s intervals).%s", tu, cil, itype,
                                                  if (prep$Kp_lev == T) " The last interval is confounded with the last capture probability." else ""))
  } else if (prep$G == 2) {
    l <- prep$group_levels
    cols <- c("survival", paste0("survival_", l[2]))
    sv <- js_series(fit, D, cols, ci, inner)
    sv$level <- factor(vapply(l, fig_label, "", opts = opts), levels = vapply(l, fig_label, "", opts = opts))
    p <- fig_cat_effect(sv, fig_label(prep$group_name, tools::toTitleCase(prep$group_name), opts),
                        fig_label("survival", sprintf("Apparent survival per %s", tu), opts))
    if (bayes) p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = eco_pal$ink, linewidth = 1.3)
    figs[[paste0("survival_", prep$group_name)]] <- with_caption(p + ggplot2::expand_limits(y = c(min(sv$lower) - 0.05, 1)),
      sprintf("Apparent survival per %s by %s (points: %s; bars: %s %s intervals%s).", tu, prep$group_name, if (bayes) "posterior means" else "MLE", cil, itype, if (bayes) ", thick 50%" else ""))
  } else {
    sv <- js_series(fit, D, "survival", ci, inner); sv$level <- factor("All animals")
    figs$survival <- with_caption(fig_cat_effect(sv, NULL, fig_label("survival", sprintf("Apparent survival per %s", tu), opts)),
                                  sprintf("Apparent survival per %s (%s %s interval).", tu, cil, itype))
  }
  # detection by occasion (and covariate curves)
  if (bayes) {
    pr <- js_param_draws(prep, js_pool(o$chains, 3000))
    pt <- js_qband(`colnames<-`(pr$p, paste0("p", seq_len(T))), paste0("p", seq_len(T)), ci, inner)
  } else {
    pt <- if (prep$Kp_lev == T) js_series(fit, NULL, paste0("detection_", seq_len(T)), ci, inner) else NULL
    if (is.null(pt)) { # MLE with covariates: plug-in curve
      p0 <- js_row(est, "detection"); lp <- stats::qlogis(p0$estimate)
      for (v in names(prep$scaling)) lp <- lp + js_row(est, paste0("detection_beta_", v))$estimate * prep$Xp[, paste0("detection_beta_", v)]
      pt <- data.frame(estimate = stats::plogis(lp), lower = NA, upper = NA)
    }
  }
  pt$x <- ax$x; pt$label <- prep$labels
  pd <- js_pointrange(pt, xlab, fig_label("detection", "Capture probability", opts), ax) + ggplot2::coord_cartesian(ylim = c(0, 1))
  figs$detection_time <- with_caption(pd, sprintf("Capture probability of an animal present at each occasion (%s %s intervals)%s.", cil, itype,
                                                  if (length(prep$scaling)) paste0(", as predicted from ", paste(tolower(vapply(names(prep$scaling), function(v) fig_label(v, v, opts), "")), collapse = ", ")) else ""))
  if (bayes) for (v in names(prep$scaling)) {
    sc <- prep$scaling[[v]]
    xr <- range(data$occasions[[v]]); xs <- seq(xr[1], xr[2], length.out = 100)
    zs <- (xs - sc$center) / sc$scale
    S <- js_pool(o$chains, 3000)
    k <- which(colnames(prep$Xp) == paste0("detection_beta_", v))
    lp <- S[, "beta_p[1]"] + outer(S[, sprintf("beta_p[%d]", k)], zs)
    pm <- stats::plogis(lp)
    a <- (1 - ci) / 2; bq <- (1 - inner) / 2
    pred <- data.frame(x = xs, estimate = colMeans(pm), lower = apply(pm, 2, stats::quantile, a), upper = apply(pm, 2, stats::quantile, 1 - a),
                       lower50 = apply(pm, 2, stats::quantile, bq), upper50 = apply(pm, 2, stats::quantile, 1 - bq))
    obsv <- if (isTRUE(opts$show_data)) data.frame(x = data$occasions[[v]], y = colSums(prep$ch) / pmax(1, ab$estimate)) else NULL
    p <- fig_cond_effect(pred, fig_label(v, v, opts), fig_label("detection", "Capture probability", opts), obs = obsv, ylim01 = TRUE)
    figs[[paste0("detection_", v)]] <- with_caption(p, sprintf(
      "Capture probability against %s (line: posterior mean; bands: %s and %d%% credible intervals). Points: animals caught / estimated abundance at each occasion (a rough empirical check), plotted at that occasion's %s.%s",
      tolower(fig_label(v, v, opts)), cil, round(100 * inner), tolower(fig_label(v, v, opts)), ""))
  }
  # coefficients (logit-scale effects as odds ratios)
  cp <- est[grepl("^(survival_beta_|detection_beta_)", est$param), c("param", "estimate", "lower", "upper")]
  if (nrow(cp)) {
    o2 <- opts
    if (prep$G == 2) { k <- paste0("survival_beta_", prep$group_levels[2])
      if (is.null(o2$labels[[k]])) o2$labels[[k]] <- sprintf("Survival: %s vs %s", fig_label(prep$group_levels[2], prep$group_levels[2], opts), fig_label(prep$group_levels[1], prep$group_levels[1], opts)) }
    for (v in names(prep$scaling)) { k <- paste0("detection_beta_", v)
      if (is.null(o2$labels[[k]])) o2$labels[[k]] <- sprintf("Capture: +1 SD %s", gsub(" *\\(.*\\)", "", tolower(fig_label(v, v, opts)))) }
    figs$coefficients <- with_caption(fig_coef(cp, null_value = 1, transform = exp, log_x = TRUE,
                                               xlab = sprintf("Odds ratio (%s %s int.)", cil, if (bayes) "cred." else "conf."), opts = o2),
      sprintf("Effects on the odds of apparent survival (survival_beta_*: %s level vs reference) and of capture (detection_beta_*: per 1 SD of the covariate). Dashed line: no effect.",
              if (prep$G == 2) prep$group_levels[2] else "group"))
    attr(figs$coefficients, "fig_size")$width <- "onehalf"
  }
  # goodness of fit
  if (!is.null(o$ppc)) {
    d <- o$ppc$draws
    p1 <- fig_ppc(d[, "ind_rep"] - d[, "ind_obs"], 0, "Replicate minus observed discrepancy\n(captures per animal)")
    p2 <- fig_ppc(d[, "occ_rep"] - d[, "occ_obs"], 0, "Replicate minus observed discrepancy\n(animals caught per occasion)")
    figs$gof <- with_caption(js_compose(list(p1, p2), opts, ncol = 2),
      sprintf("Posterior predictive checks (Freeman-Tukey discrepancies, %d posterior draws). Each histogram shows replicated minus observed discrepancy; the line at 0 should lie well inside the histogram (Bayesian p-values %.2f and %.2f; extreme if < 0.025 or > 0.975).%s",
              nrow(d), o$ppc$p[["individual"]], o$ppc$p[["occasion"]],
              if (!is.null(o$gof)) paste0(" R2ucare tests: ", paste(sprintf("%s [%s] p = %.2f", o$gof$test, ifelse(o$gof$group == "unknown", paste("unknown", prep$group_name %||% "group"), o$gof$group), o$gof$p), collapse = "; "), ".") else ""))
    attr(figs$gof, "fig_size") <- list(width = "double", aspect = 0.4)
  }
  # superpopulation posterior
  if (bayes) {
    ns <- D[, "superpopulation_size"]
    p <- ggplot2::ggplot(data.frame(v = ns), ggplot2::aes(v)) +
      ggplot2::geom_histogram(binwidth = max(1, round(diff(range(ns)) / 40)), fill = eco_pal$band, alpha = 0.5, colour = NA) +
      ggplot2::geom_vline(xintercept = prep$n, linetype = 2, colour = eco_pal$ink) +
      ggplot2::labs(x = fig_label("superpopulation_size", "Superpopulation size (animals ever present)", opts), y = "Posterior draws")
    figs$superpopulation <- with_caption(p, sprintf(
      "Posterior distribution of the superpopulation size; dashed line = number of animals caught (%d). Data-augmentation limit M = %d; 99.5%% posterior quantile %.0f, %s.",
      prep$n, prep$M, stats::quantile(ns, 0.995),
      if (stats::quantile(ns, 0.995) < 0.9 * prep$M) "so M does not constrain the estimate" else "TOO CLOSE to M: the estimate may be truncated; refit with a larger settings$M"))
    keep <- intersect(c("superpopulation_size", "abundance_1", paste0("abundance_", T), "survival", "detection",
                        if (prep$G == 2) paste0("survival_beta_", prep$group_levels[2]), paste0("detection_beta_", names(prep$scaling))), colnames(D))
    short <- c(superpopulation_size = "N super", survival = sprintf("survival%s", if (prep$G == 2) paste0(" (", prep$group_levels[1], ")") else ""),
               detection = "capture prob.")
    short[paste0("abundance_", c(1, T))] <- paste("N", prep$labels[c(1, T)])
    if (prep$G == 2) short[paste0("survival_beta_", prep$group_levels[2])] <- paste("survival effect,", prep$group_levels[2])
    for (v in names(prep$scaling)) short[paste0("detection_beta_", v)] <- "capture effect /SD"
    sl <- unname(ifelse(keep %in% names(short), short[keep], keep))
    rn <- function(m) { m <- m[, keep, drop = FALSE]; colnames(m) <- sl; m }
    figs$trace <- with_caption(fig_trace(lapply(o$rep_chains, rn)) + ggplot2::scale_x_continuous(n.breaks = 3) +
                                 ggplot2::scale_y_continuous(n.breaks = 3),
                               "MCMC traces of key quantities (thinned), one colour per chain; chains should overlap.")
    attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.55)
    figs$posterior <- with_caption(fig_posterior(rn(D), sl) + ggplot2::scale_x_continuous(n.breaks = 3) +
                                     ggplot2::scale_y_continuous(n.breaks = 3),
                                   "Posterior densities of key quantities (N = abundance; survival per time unit; effects on the logit scale).")
    attr(figs$posterior, "fig_size") <- list(width = "double", aspect = 0.55)
  }
  if (!is.null(o$comparison)) {
    cm <- o$comparison
    keep <- intersect(c("superpopulation_size", "abundance_1", paste0("abundance_", ceiling(T / 2)), paste0("abundance_", T),
                        "survival", "detection", paste0("detection_beta_", names(prep$scaling))), cm$param)
    cm <- cm[cm$param %in% keep, ]; cm$param <- factor(cm$param, levels = keep)
    cmp_lab <- function(x) vapply(x, function(k) {
      if (k == "superpopulation_size") return("Superpopulation size")
      if (grepl("^abundance_", k)) return(paste("Abundance,", prep$labels[as.integer(sub("abundance_", "", k))]))
      if (k == "survival") return(if (prep$G == 2) sprintf("Survival per %s (%s)", prep$time_unit %||% "interval", prep$group_levels[1]) else "Survival")
      if (k == "detection") return("Capture probability")
      if (grepl("^detection_beta_", k)) return(sprintf("%s (logit/SD)", (gsub(" *\\(.*\\)", "", fig_label(sub("^detection_beta_", "", k), sub("^detection_beta_", "", k), opts)))))
      gsub("_", " ", k) }, "")
    # degenerate ML intervals (entry probability on a boundary) would blow up
    # the axis: clip each panel to the Bayesian interval range +- 60%
    ref <- cm[grepl("nimble", cm$engine), ]
    for (k in keep) { r <- ref[ref$param == k, ]; if (!nrow(r)) next
      w <- r$upper[1] - r$lower[1]; lo <- r$lower[1] - 0.6 * w; hi <- r$upper[1] + 0.6 * w
      i <- cm$param == k; cm$lower[i] <- pmax(cm$lower[i], lo); cm$upper[i] <- pmin(cm$upper[i], hi) }
    p <- ggplot2::ggplot(cm, ggplot2::aes(estimate, engine, colour = engine, shape = engine)) +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = lower, xmax = upper), size = 0.25) +
      ggplot2::facet_wrap(~param, scales = "free_x", ncol = 3,
                          labeller = ggplot2::as_labeller(cmp_lab)) +
      # drop ticks near panel edges so labels of adjacent facets cannot collide
      ggplot2::scale_x_continuous(breaks = function(lim) { b <- pretty(lim, 4); r <- diff(lim)
                                    b[b > lim[1] + 0.1 * r & b < lim[2] - 0.1 * r] },
                                  expand = ggplot2::expansion(mult = 0.2)) +
      ggplot2::guides(colour = "none", shape = "none") +
      ggplot2::labs(x = "Estimate (95% interval)", y = NULL, colour = NULL, shape = NULL) +
      ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    figs$engine_comparison <- with_caption(p, paste("Estimates of the same quantities from each engine fitted to the same data. marked gives expected (not realised) abundance;",
      "its intervals come from simulating the coefficients' asymptotic normal distribution and are clipped to the panel when an entry probability lies on its boundary.",
      if (prep$G == 2 && any(prep$gcode == 0)) "Animals of unknown group were present, so the marked fit has no group effect and its survival is pooled over groups." else ""))
    attr(figs$engine_comparison, "fig_size") <- list(width = "double", aspect = 0.55)
  }
  figs
}
