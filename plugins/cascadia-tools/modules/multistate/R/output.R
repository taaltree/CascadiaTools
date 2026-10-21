# multistate: interpretation, figures and methods text.
#
# All curves are computed from the posterior draws (Nimble) or parametric
# draws from the ML fit (engine "ml"), both stored in fit$object$draws with
# Nimble parameter names. Covariates were standardised in ms_prepare(); x axes
# are drawn in original units (des$scaling or the user's figures: units:).

ms_q <- function(M, ci = 0.95, inner = 0.5) {
  M <- as.matrix(M); a <- (1 - ci) / 2; b <- (1 - inner) / 2
  data.frame(estimate = colMeans(M, na.rm = TRUE),
             lower = apply(M, 2, stats::quantile, a, na.rm = TRUE),
             upper = apply(M, 2, stats::quantile, 1 - a, na.rm = TRUE),
             lower50 = apply(M, 2, stats::quantile, b, na.rm = TRUE),
             upper50 = apply(M, 2, stats::quantile, 1 - b, na.rm = TRUE), row.names = NULL)
}

ms_unit <- function(fit, data) {
  u <- data$interval_unit
  if (is.null(u) || !nzchar(u)) "interval between occasions" else u
}

ms_est <- function(fit, p) {
  e <- fit$estimates[fit$estimates$param == p, ]
  if (!nrow(e)) return(NULL)
  e
}

# Grid of standardised covariate values spanning the observed range, with
# original-unit values for the axis.
ms_cov_grid <- function(des, k, n = 60, opts = NULL) {
  xs <- des$raw_covs[[k]]$x
  if (des$raw_covs[[k]]$type == "occasion") {
    in_d <- k %in% des$covs$detection; in_g <- k %in% c(des$covs$survival, des$covs$transition)
    xs <- if (in_d && !in_g) xs[-1] else if (in_g && !in_d) xs[-length(xs)] else xs
  }
  g <- seq(min(xs, na.rm = TRUE), max(xs, na.rm = TRUE), length.out = n)
  sc <- des$scaling[[k]]
  u <- if (!is.null(opts$units[[k]])) opts$units[[k]] else list(center = sc[["center"]], scale = sc[["scale"]])
  list(z = g, x = g * u$scale + u$center, obs = xs * u$scale + u$center, zobs = xs)
}

ms_state_lab <- function(des, opts) unname(vapply(des$states, function(s) fig_label(s, s, opts), ""))

# ---- conditional-effect predictions ---------------------------------------

ms_pred_survival <- function(fit, k, opts = fig_opts()) {
  des <- fit$object$des; D <- as.matrix(fit$object$draws)
  g <- ms_cov_grid(des, k, opts = opts); j <- match(k, des$covs$survival)
  lab <- ms_state_lab(des, opts)
  do.call(rbind, lapply(seq_len(des$S), function(s) {
    lp <- outer(ms_get(D, sprintf("a_phi[%d]", s)), rep(1, length(g$z))) +
      outer(ms_get(D, sprintf("bphi[%d]", j)), g$z)
    cbind(data.frame(x = g$x, group = lab[s]), ms_q(stats::plogis(lp), opts$ci, opts$inner_ci))
  }))
}

ms_pred_detection <- function(fit, k, opts = fig_opts()) {
  des <- fit$object$des; D <- as.matrix(fit$object$draws)
  g <- ms_cov_grid(des, k, opts = opts); j <- match(k, des$covs$detection)
  lab <- ms_state_lab(des, opts)
  do.call(rbind, lapply(seq_len(des$S), function(s) {
    lp <- outer(ms_get(D, sprintf("a_p[%d]", s)), rep(1, length(g$z))) +
      outer(ms_get(D, sprintf("bp[%d]", j)), g$z)
    cbind(data.frame(x = g$x, group = lab[s]), ms_q(stats::plogis(lp), opts$ci, opts$inner_ci))
  }))
}

# Transition probabilities (each off-diagonal pair) along covariate k, other
# transition covariates at their mean (0 on the standardised scale).
ms_pred_transition <- function(fit, k, opts = fig_opts()) {
  des <- fit$object$des; D <- as.matrix(fit$object$draws)
  g <- ms_cov_grid(des, k, n = 40, opts = opts); j <- match(k, des$covs$transition)
  lab <- ms_state_lab(des, opts)
  arr <- lapply(g$z, function(z) { xs <- rep(0, des$Kpsi); xs[j] <- z; ms_psi_draws(D, des, xs) })
  do.call(rbind, lapply(seq_len(nrow(des$pairs)), function(m) {
    r <- des$pairs[m, 1]; s <- des$pairs[m, 2]
    M <- vapply(arr, function(P) P[, r, s], numeric(nrow(D)))
    cbind(data.frame(x = g$x, group = sprintf("%s to %s", lab[r], lab[s])), ms_q(matrix(M, nrow(D)), opts$ci, opts$inner_ci))
  }))
}

# Detection by occasion (occasion covariates on detection)
ms_pred_det_occasion <- function(fit, opts = fig_opts()) {
  des <- fit$object$des; D <- as.matrix(fit$object$draws); lab <- ms_state_lab(des, opts)
  occ <- if (!is.null(des$occasion_labels)) des$occasion_labels else as.character(seq_len(des$T))
  do.call(rbind, lapply(seq_len(des$S), function(s) {
    M <- vapply(2:des$T, function(t) {
      lp <- ms_get(D, sprintf("a_p[%d]", s))
      for (k in seq_len(des$Kp)) lp <- lp + ms_get(D, sprintf("bp[%d]", k)) * des$Xp[t, k]
      stats::plogis(lp)
    }, numeric(nrow(D)))
    cbind(data.frame(time = 2:des$T, occasion = occ[2:des$T], group = lab[s]),
          ms_q(matrix(M, nrow(D)), opts$ci, opts$inner_ci))
  }))
}

# ---- interpretation -------------------------------------------------------

ms_interpret <- function(fit, data) {
  des <- fit$object$des; st <- des$states; S <- des$S; unit <- ms_unit(fit, data)
  bayes <- identical(fit$object$engine, "nimble")
  ivl <- if (bayes) "credible interval" else "confidence interval"
  fmt <- function(e, d = 2) sprintf(paste0("%.", d, "f (95%% ", ivl, " %.", d, "f-%.", d, "f)"), e$estimate, e$lower, e$upper)
  out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge (see Diagnostics). The numbers below are not reliable.")
  covtxt <- if (des$Kphi) " for an animal with average covariate values" else ""
  sv <- vapply(st, function(s) fmt(ms_est(fit, paste0("survival_", s))), "")
  out <- c(out, sprintf("Survival over one %s%s was %s.", unit, covtxt,
                        paste(sprintf("%s in state '%s'", sv, st), collapse = "; ")))
  if (S == 2) {
    P <- ms_derive(fit$object$draws, des)
    diff <- P[, paste0("survival_", st[1])] - P[, paste0("survival_", st[2])]
    out <- c(out, sprintf("The difference in survival (%s minus %s) was %.2f (95%% interval %.2f to %.2f)%s.",
                          st[1], st[2], mean(diff), stats::quantile(diff, 0.025), stats::quantile(diff, 0.975),
                          if (bayes) sprintf("; posterior probability that survival is higher in '%s' = %.2f", st[1], mean(diff > 0)) else ""))
  }
  tr <- vapply(seq_len(nrow(des$pairs)), function(m) {
    r <- st[des$pairs[m, 1]]; s <- st[des$pairs[m, 2]]
    sprintf("from '%s' to '%s' %s", r, s, fmt(ms_est(fit, sprintf("transition_%s_to_%s", r, s))))
  }, "")
  out <- c(out, sprintf(paste("Transitions are probabilities of being in a new state at the next occasion, for animals that survive the %s",
                              "(at average transition covariates): %s."), unit, paste(tr, collapse = "; ")))
  # residence times
  P <- ms_psi_draws(as.matrix(fit$object$draws), des)
  res <- vapply(seq_len(S), function(r) {
    stay <- P[, r, r]; rt <- 1 / (1 - stay)
    sprintf("'%s' %.1f (%.1f-%.1f)", st[r], stats::median(rt), stats::quantile(rt, 0.025), stats::quantile(rt, 0.975))
  }, "")
  out <- c(out, sprintf(paste("Expected number of %ss a surviving animal stays in a state before switching (1 / switching probability, median and 95%% interval):",
                              "%s. Death is not counted here; with survival included, stays are shorter."), unit, paste(res, collapse = "; ")))
  stn <- vapply(st, function(s) fmt(ms_est(fit, paste0("stationary_", s))), "")
  out <- c(out, sprintf(paste("If the current transition probabilities persisted, the long-run (stationary) share of surviving animals in each state would be %s.",
                              "This describes the movement/switching process only: it ignores differences in survival and recruitment between states,",
                              "so it is not a prediction of the population's state structure."),
                        paste(sprintf("%s in '%s'", stn, st), collapse = "; ")))
  dt <- vapply(st, function(s) fmt(ms_est(fit, paste0("detection_", s))), "")
  out <- c(out, sprintf("Detection (the chance that a live animal is caught or resighted at an occasion%s) was %s. Low detection in a state means fewer data on that state, not fewer animals.",
                        if (des$Kp) ", at average detection covariates" else "", paste(sprintf("%s in '%s'", dt, st), collapse = "; ")))
  if (des$model == "multievent") {
    at <- vapply(st, function(s) fmt(ms_est(fit, paste0("assignment_", s))), "")
    it <- vapply(st, function(s) fmt(ms_est(fit, paste0("initial_", s))), "")
    out <- c(out, sprintf(paste("State uncertainty (multievent model): when a detected animal was in a given state, the probability that its state was recorded",
                                "(rather than coded unknown) was %s. Among newly marked animals, the estimated share in each state was %s.",
                                "Animals recorded with an unknown state were assigned to states probabilistically by the model, not discarded."),
                          paste(sprintf("%s in '%s'", at, st), collapse = "; "), paste(sprintf("%s in '%s'", it, st), collapse = "; ")))
  }
  # covariate effects
  D <- as.matrix(fit$object$draws)
  span <- function(k) { g <- ms_cov_grid(des, k); c(lo = stats::quantile(g$obs, 0.1, names = FALSE), hi = stats::quantile(g$obs, 0.9, names = FALSE),
                                                      zlo = stats::quantile(g$zobs, 0.1, names = FALSE), zhi = stats::quantile(g$zobs, 0.9, names = FALSE)) }
  pr_pos <- function(v) if (bayes) sprintf("; P(effect > 0) = %.2f", mean(v > 0)) else ""
  for (j in seq_len(des$Kphi)) {
    k <- des$covs$survival[j]; e <- ms_est(fit, paste0("survival_beta_", k)); sp <- span(k)
    b <- ms_get(D, sprintf("bphi[%d]", j))
    chg <- vapply(seq_len(S), function(s) {
      a <- ms_get(D, sprintf("a_phi[%d]", s))
      lo <- stats::plogis(a + b * sp[["zlo"]]); hi <- stats::plogis(a + b * sp[["zhi"]])
      sprintf("in '%s' from %.2f to %.2f", st[s], mean(lo), mean(hi))
    }, "")
    out <- c(out, sprintf(paste("Survival and %s: odds ratio per 1 SD (= %.3g) %.2f (95%% interval %.2f-%.2f)%s.",
                                "Going from the 10th to the 90th percentile of %s (%.3g to %.3g), predicted survival changes %s.",
                                "The same slope is assumed in every state. This is an association, not proof of cause."),
                          k, des$scaling[[k]][["scale"]], exp(e$estimate), exp(e$lower), exp(e$upper), pr_pos(b),
                          k, sp[["lo"]], sp[["hi"]], paste(chg, collapse = ", ")))
  }
  for (j in seq_len(des$Kpsi)) {
    k <- des$covs$transition[j]; sp <- span(k)
    xs_lo <- rep(0, des$Kpsi); xs_lo[j] <- sp[["zlo"]]; xs_hi <- xs_lo; xs_hi[j] <- sp[["zhi"]]
    Plo <- ms_psi_draws(D, des, xs_lo); Phi <- ms_psi_draws(D, des, xs_hi)
    txt <- vapply(seq_len(nrow(des$pairs)), function(m) {
      r <- des$pairs[m, 1]; s <- des$pairs[m, 2]
      e <- ms_est(fit, sprintf("transition_beta_%s_%s_to_%s", k, st[r], st[s]))
      d <- Phi[, r, s] - Plo[, r, s]
      sprintf("'%s' to '%s': %.2f at the 10th percentile vs %.2f at the 90th (difference %.2f, 95%% interval %.2f to %.2f; link-scale slope %.2f%s)",
              st[r], st[s], mean(Plo[, r, s]), mean(Phi[, r, s]), mean(d), stats::quantile(d, 0.025), stats::quantile(d, 0.975),
              e$estimate, pr_pos(ms_get(D, sprintf("bpsi[%d, %d, %d]", r, s, j))))
    }, "")
    out <- c(out, sprintf("Transitions and %s (from %.3g to %.3g, 10th to 90th percentile): %s. Associations only; transitions are on the multinomial-logit scale, so a slope describes the odds of moving relative to staying.",
                          k, sp[["lo"]], sp[["hi"]], paste(txt, collapse = "; ")))
  }
  for (j in seq_len(des$Kp)) {
    k <- des$covs$detection[j]; e <- ms_est(fit, paste0("detection_beta_", k))
    out <- c(out, sprintf("Detection and %s: odds ratio per 1 SD (= %.3g) %.2f (95%% interval %.2f-%.2f). This is about how easily animals were caught or seen, not about their survival or movement.",
                          k, des$scaling[[k]][["scale"]], exp(e$estimate), exp(e$lower), exp(e$upper)))
  }
  g <- fit$object$gof
  if (!is.null(g)) {
    bad <- c(g$p_ft, g$p_switch, g$p_never)
    out <- c(out, if (any(bad < 0.05 | bad > 0.95))
      sprintf("Goodness of fit: at least one posterior predictive check is extreme (p = %.2f, %.2f, %.2f for sightings, state switches, never-seen-again). Treat the estimates with caution; see Diagnostics.", g$p_ft, g$p_switch, g$p_never)
      else sprintf("Goodness of fit: posterior predictive checks showed no clear misfit (Bayesian p = %.2f sightings per occasion, %.2f state switches, %.2f animals never seen again).", g$p_ft, g$p_switch, g$p_never))
  }
  out
}

# ---- figures ---------------------------------------------------------------

ms_figures <- function(fit, data, opts = fig_opts()) {
  if (is.null(opts$ci)) opts <- fig_opts(opts)
  des <- fit$object$des; st <- des$states; lab <- ms_state_lab(des, opts); S <- des$S
  bayes <- identical(fit$object$engine, "nimble")
  ivl <- if (bayes) "credible" else "confidence"
  cil <- round(100 * opts$ci); inl <- round(100 * opts$inner_ci)
  bandtxt <- if (bayes) sprintf("dark band = %d%% and light band = %d%% %s intervals", inl, cil, ivl) else sprintf("bands = %d%% (light) and %d%% (dark) %s intervals from parametric draws", cil, inl, ivl)
  unit <- ms_unit(fit, data)
  est <- fit$estimates
  figs <- list()
  stfac <- function(v) factor(v, levels = lab)
  # 1. survival / detection / assignment by state
  rows <- function(prefix, what) {
    d <- est[est$param %in% paste0(prefix, st), ]
    d$level <- factor(lab[match(sub(prefix, "", d$param), st)], levels = lab); d$what <- what; d
  }
  vr <- rbind(rows("survival_", fig_label("survival", sprintf("Survival per %s", unit), opts)),
              rows("detection_", fig_label("detection", "Detection probability", opts)))
  if (des$model == "multievent") vr <- rbind(vr, rows("assignment_", fig_label("assignment", "State recorded\n(if detected)", opts)),
                                             rows("initial_", fig_label("initial", "Share at\nfirst capture", opts)))
  vr$what <- factor(vr$what, levels = unique(vr$what))
  p <- ggplot2::ggplot(vr, ggplot2::aes(level, estimate, colour = level, shape = level)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.3, linewidth = 0.6) +
    ggplot2::facet_wrap(~what, nrow = 1) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = fig_label("state_axis", "State", opts), y = "Probability", colour = NULL, shape = NULL) +
    ggplot2::guides(colour = "none", shape = "none")
  attr(p, "fig_size") <- list(width = "double", aspect = 0.35)
  figs$vital_rates_by_state <- with_caption(p, sprintf(
    "State-specific survival (per %s), detection%s, with posterior means and %d%% %s intervals%s.", unit,
    if (des$model == "multievent") ", probability that a detected animal's state was recorded, and the state distribution at first capture" else "",
    cil, ivl, if (des$Kphi + des$Kp) " (at average covariate values)" else ""))
  # 2. transition matrix (heatmap) + forest of transitions
  P <- ms_psi_draws(as.matrix(fit$object$draws), des)
  Pm <- apply(P, c(2, 3), mean)
  hd <- expand.grid(from = seq_len(S), to = seq_len(S))
  hd$value <- Pm[cbind(hd$from, hd$to)]
  hd$lo <- apply(P, c(2, 3), stats::quantile, (1 - opts$ci) / 2)[cbind(hd$from, hd$to)]
  hd$hi <- apply(P, c(2, 3), stats::quantile, 1 - (1 - opts$ci) / 2)[cbind(hd$from, hd$to)]
  hd$from <- factor(lab[hd$from], levels = rev(lab)); hd$to <- factor(lab[hd$to], levels = lab)
  hd$txt <- sprintf("%.2f\n(%.2f-%.2f)", hd$value, hd$lo, hd$hi)
  ph <- ggplot2::ggplot(hd, ggplot2::aes(to, from, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = txt, colour = value > 0.5), size = 2.6, lineheight = 0.9, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
    scale_fill_eco_seq(limits = c(0, 1)) + ggplot2::coord_equal() +
    ggplot2::labs(x = "State at next occasion", y = "State now", fill = "Probability")
  figs$transitions <- with_caption(ph,
    sprintf("Transition probabilities between states over one %s for animals that survive it: posterior mean (%d%% %s interval). Rows sum to 1; off-diagonal cells are movement/switching, diagonal cells are staying.%s", unit, cil, ivl,
            if (des$Kpsi) " Evaluated at average transition covariates." else ""))
  trd <- est[grepl("^transition_[^b]", est$param) & !grepl("^transition_beta_", est$param), ]
  trd$param <- vapply(trd$param, function(p0) { r <- regmatches(p0, regexec("^transition_(.*)_to_(.*)$", p0))[[1]]
    sprintf("%s to %s", fig_label(r[2], r[2], opts), fig_label(r[3], r[3], opts)) }, "")
  trd$param <- factor(trd$param, levels = rev(trd$param))
  pt <- ggplot2::ggplot(trd, ggplot2::aes(estimate, param)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.6) +
    ggplot2::geom_point(size = 1.8) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = sprintf("Transition probability per %s\n(%d%% %s interval)", unit, cil, ivl), y = NULL)
  attr(pt, "fig_size") <- list(width = "onehalf", height = 28 + 8 * nrow(trd))
  figs$transition_estimates <- with_caption(pt,
    sprintf("Estimated transition probabilities (points) with %d%% %s intervals.", cil, ivl))
  # 3. stationary distribution
  sd0 <- est[est$param %in% paste0("stationary_", st), ]
  sd0$level <- factor(lab[match(sub("stationary_", "", sd0$param), st)], levels = lab)
  figs$stationary <- with_caption(
    ggplot2::ggplot(sd0, ggplot2::aes(level, estimate, colour = level, shape = level)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.3, linewidth = 0.6) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) + ggplot2::guides(colour = "none", shape = "none") +
      ggplot2::labs(x = fig_label("state_axis", "State", opts), y = "Long-run share\nof surviving animals"),
    sprintf("Stationary state distribution implied by the transition matrix (posterior mean, %d%% %s interval): the share of time a surviving animal would spend in each state if current transition probabilities persisted. It ignores state differences in survival and recruitment.", cil, ivl))
  # 4. covariate effects
  for (k in des$covs$survival) {
    pr <- ms_pred_survival(fit, k, opts); pr$group <- stfac(pr$group)
    g <- ms_cov_grid(des, k, opts = opts)
    figs[[paste0("survival_", k)]] <- with_caption(
      fig_cond_effect(pr, fig_label(k, k, opts), fig_label("survival_y", sprintf("Survival per %s", unit), opts),
                      obs = if (isTRUE(opts$show_data)) data.frame(x = g$obs) else NULL, ylim01 = TRUE) +
        ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL),
      sprintf("Survival per %s against %s, one curve per state (other covariates at their mean); %s. Rug = values of the marked animals.", unit, fig_label(k, k, opts), bandtxt))
  }
  for (k in des$covs$transition) {
    pr <- ms_pred_transition(fit, k, opts)
    g <- ms_cov_grid(des, k, opts = opts)
    figs[[paste0("transition_", k)]] <- with_caption(
      fig_cond_effect(pr, fig_label(k, k, opts), fig_label("transition_y", "Transition probability", opts),
                      obs = if (isTRUE(opts$show_data)) data.frame(x = g$obs) else NULL, ylim01 = TRUE) +
        ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL),
      sprintf("Probability of switching state over one %s (for survivors) against %s, one curve per transition; %s. Rug = observed values.", unit, fig_label(k, k, opts), bandtxt))
  }
  for (k in des$covs$detection) {
    pr <- ms_pred_detection(fit, k, opts); pr$group <- stfac(pr$group)
    g <- ms_cov_grid(des, k, opts = opts)
    figs[[paste0("detection_", k)]] <- with_caption(
      fig_cond_effect(pr, fig_label(k, k, opts), fig_label("detection_y", "Detection probability", opts),
                      obs = if (isTRUE(opts$show_data)) data.frame(x = g$obs) else NULL, ylim01 = TRUE) +
        ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL),
      sprintf("Detection probability against %s, one curve per state; %s. Rug = observed values.", fig_label(k, k, opts), bandtxt))
  }
  if (des$Kp > 0) {
    pd <- ms_pred_det_occasion(fit, opts); pd$group <- stfac(pd$group)
    pdg <- ggplot2::position_dodge(width = 0.4)
    p <- ggplot2::ggplot(pd, ggplot2::aes(time, estimate, colour = group, shape = group)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), position = pdg, linewidth = 0.5) +
      ggplot2::geom_point(position = pdg, size = 1.4) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_x_continuous(breaks = pd$time, labels = pd$occasion, guide = ggplot2::guide_axis(angle = if (des$T > 5) 45 else 0)) +
      ggplot2::labs(x = fig_label("occasion", "Occasion", opts), y = "Detection probability", colour = NULL, shape = NULL)
    figs$detection_by_occasion <- with_caption(p, sprintf("Detection probability at each recapture occasion and state (posterior mean, %d%% %s interval), driven by the occasion covariates.", cil, ivl))
  }
  for (n in grep("^(survival|transition|detection)_", names(figs), value = TRUE))
    if (n %in% paste0(c("survival_", "transition_", "detection_"), rep(unique(unlist(des$covs)), each = 3)))
      attr(figs[[n]], "fig_size") <- list(width = "onehalf", aspect = 0.65)
  # 5. coefficients (link scale, odds ratios)
  co <- est[grepl("_beta_", est$param), ]
  if (nrow(co)) {
    co$group <- ifelse(grepl("^survival", co$param), "Survival", ifelse(grepl("^transition", co$param), "Transition", "Detection"))
    nice <- function(p0) {
      for (k in unique(unlist(des$covs))) {
        kl <- fig_label(k, k, opts)
        if (p0 == paste0("survival_beta_", k)) return(sprintf("Survival: %s", kl))
        if (p0 == paste0("detection_beta_", k)) return(sprintf("Detection: %s", kl))
        for (m in seq_len(nrow(des$pairs))) {
          r <- des$pairs[m, 1]; s2 <- des$pairs[m, 2]
          if (p0 == sprintf("transition_beta_%s_%s_to_%s", k, st[r], st[s2])) return(sprintf("%s to %s: %s", lab[r], lab[s2], kl))
        }
      }
      fig_label(p0, p0, opts)
    }
    co$label <- vapply(co$param, function(p0) paste(strwrap(fig_label(p0, nice(p0), opts), 30), collapse = "\n"), "")
    co$label <- factor(co$label, levels = rev(co$label))
    for (v in c("estimate", "lower", "upper")) co[[v]] <- exp(co[[v]])
    pc <- ggplot2::ggplot(co, ggplot2::aes(estimate, label, colour = group, shape = group)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.6) +
      ggplot2::geom_point(size = 1.8) + ggplot2::scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10), labels = function(x) sprintf("%g", x)) +
      ggplot2::labs(x = sprintf("Odds ratio per 1 SD\n(%d%% %s interval)", cil, ivl), y = NULL, colour = NULL, shape = NULL) +
      ggplot2::guides(colour = "none", shape = "none")   # the process is already named at the start of each label
    attr(pc, "fig_size") <- list(width = "onehalf", height = 30 + 11 * nrow(co))
    figs$coefficients <- with_caption(pc,
      "Covariate effects as odds ratios per 1 SD of the covariate (transitions: odds of moving relative to staying), log scale. Dashed line = no effect.")
  }
  # 6. observed data summary
  y <- des$y
  ev <- do.call(rbind, lapply(seq_len(des$T), function(t) {
    v <- y[, t]; tb <- tabulate(v[v > 0], S + 1)
    data.frame(time = t, event = c(lab, "State unknown"), n = tb)
  }))
  if (all(ev$n[ev$event == "State unknown"] == 0)) ev <- ev[ev$event != "State unknown", ]
  ev$event <- factor(ev$event, levels = c(lab, "State unknown"))
  occ <- if (!is.null(des$occasion_labels)) des$occasion_labels else as.character(seq_len(des$T))
  figs$captures_by_occasion <- with_caption(
    ggplot2::ggplot(ev, ggplot2::aes(time, n, fill = event)) + ggplot2::geom_col(width = 0.75) +
      ggplot2::scale_x_continuous(breaks = seq_len(des$T), labels = occ, guide = ggplot2::guide_axis(angle = if (des$T > 5) 45 else 0)) +
      ggplot2::labs(x = fig_label("occasion", "Occasion", opts), y = "Animals caught or seen", fill = NULL),
    "Number of animals caught or resighted at each occasion, by recorded state (raw data).")
  # 7. goodness of fit
  gf <- fit$object$gof
  if (!is.null(gf)) {
    figs$gof <- with_caption(compose_figures(list(
      fig_ppc(gf$ft_rep - gf$ft_obs, 0, "FT discrepancy,\nreplicate - observed") + ggplot2::labs(caption = NULL),
      fig_ppc(gf$switch_rep, gf$switch_obs, "Recorded state\nswitches") + ggplot2::labs(caption = NULL),
      fig_ppc(gf$never_rep, gf$never_obs, "Animals never\nseen again") + ggplot2::labs(caption = NULL)), ncol = 3) &
        theme_eco(opts$base_size, opts$font, opts$style, opts$legend),
      sprintf("Posterior predictive checks (%d draws). Histograms: replicated statistic; line: observed. A: Freeman-Tukey discrepancy of sightings per occasion and event, replicated minus observed (line at 0; p = %.2f). B: state switches between consecutive known-state sightings (p = %.2f). C: animals never seen after first capture (p = %.2f). p near 0 or 1 indicates misfit.",
              gf$n_draws, gf$p_ft, gf$p_switch, gf$p_never))
    attr(figs$gof, "fig_size") <- list(width = "double", aspect = 0.33)
  }
  # 8. MCMC trace and posterior (Bayesian)
  if (bayes && !is.null(fit$object$chains)) {
    keep <- c(paste0("survival_", st), sprintf("transition_%s_to_%s", st[des$pairs[, 1]], st[des$pairs[, 2]]))
    wrapn <- function(v) vapply(v, function(z) paste(strwrap(gsub("_", " ", fig_label(z, z, opts)), 24), collapse = "\n"), "")
    ch <- lapply(fit$object$chains, function(m) { x <- ms_derive(m, des)[, keep, drop = FALSE]; colnames(x) <- wrapn(colnames(x)); x })
    figs$trace <- with_caption(fig_trace(ch, wrapn(keep)), "MCMC traces of survival and transition probabilities, one colour per chain (after burn-in). Chains should overlap like a 'hairy caterpillar'.")
    attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.6)
    D2 <- ms_derive(fit$object$draws, des)
    D2 <- D2[, setdiff(colnames(D2), paste0("stationary_", st)), drop = FALSE]; colnames(D2) <- wrapn(colnames(D2))
    figs$posterior <- with_caption(fig_posterior(D2, colnames(D2)),
                                   "Posterior densities of all reported quantities.")
    attr(figs$posterior, "fig_size") <- list(width = "double", aspect = 0.8)
  }
  figs
}

ms_plots <- function(fit, data) {
  f <- tryCatch(ms_figures(fit, data, fig_opts()), error = function(e) list())
  f[intersect(names(f), c("vital_rates_by_state", "transitions", grep("^(survival|transition|detection)_[^b]", names(f), value = TRUE),
                          "stationary", "gof"))]
}

ms_methods_text <- function(fit, data, settings) {
  des <- fit$object$des; st <- des$states; s <- fit$object$settings; unit <- ms_unit(fit, data)
  me <- des$model == "multievent"
  base <- sprintf(paste(
    "We analysed capture histories of %d individually marked animals over %d occasions with a %s model (%s) with %d states (%s) plus dead.",
    "Each %s, an animal in state r survived with probability phi_r and, if it survived, moved to state s with probability psi_rs",
    "(multinomial logit with staying as the reference category); a live animal in state s was detected with probability p_s.%s",
    "The likelihood was conditioned on each animal's first capture and the latent state sequence was summed out with the forward algorithm",
    "(hidden Markov model formulation; Pradel 2005; Turek et al. 2016)."),
    des$N, des$T, if (me) "multievent" else "multistate Arnason-Schwarz",
    if (me) "Pradel 2005; Conn & Cooch 2009" else "Arnason 1973; Schwarz et al. 1993; Lebreton et al. 2009",
    des$S, paste(st, collapse = ", "), unit,
    if (me) paste(" When detected, an animal's state was recorded with probability delta_s and otherwise coded 'unknown' (no misclassification);",
                  "the state distribution at first capture (pi) was estimated from the first events.") else "")
  cv <- character()
  if (des$Kphi) cv <- c(cv, sprintf("survival on %s (logit link, slope shared among states)", paste(des$covs$survival, collapse = ", ")))
  if (des$Kpsi) cv <- c(cv, sprintf("transitions on %s (one slope per transition)", paste(des$covs$transition, collapse = ", ")))
  if (des$Kp) cv <- c(cv, sprintf("detection on %s (logit link, slope shared among states)", paste(des$covs$detection, collapse = ", ")))
  cvt <- if (length(cv)) sprintf(" Covariates (standardised to mean 0, SD 1) were modelled on %s.", paste(cv, collapse = "; ")) else ""
  if (identical(fit$object$engine, "nimble")) {
    eng <- sprintf(paste(
      "Models were fitted in NIMBLE (de Valpine et al. 2017) using nimbleEcology's dDHMMo distribution, with logistic(0,1) priors on all logit-scale intercepts",
      "(uniform on the probability scale for survival, detection%s), logistic(0,1) priors on multinomial-logit transition intercepts, and Normal(0, SD = %g) priors on covariate slopes.",
      "We ran %d chains of %d iterations (%d burn-in) with slice samplers and checked convergence with R-hat < 1.01 and bulk effective sample size >= 400 (posterior package).",
      "Derived quantities (transition probabilities, stationary state distribution as the left eigenvector of the transition matrix) were computed per posterior draw."),
      if (me) " and state assignment" else "", s$prior_sd, s$nchains, s$niter, s$nburnin)
  } else {
    eng <- "Parameters were estimated by maximum likelihood (BFGS then Nelder-Mead, R optim) using a hand-written forward algorithm; 95% intervals for derived quantities come from parametric draws from the asymptotic normal distribution of the link-scale estimates."
  }
  gof <- paste("Goodness of fit was assessed with posterior predictive checks (Freeman-Tukey discrepancy of sightings per occasion and state, number of recorded state switches,",
               "number of animals never recaptured) and with the U-CARE tests for the Arnason-Schwarz model (Pradel et al. 2003; R2ucare, Gimenez et al. 2018).")
  val <- "The model was validated by simulation (planted-truth recovery of all parameters at a comparable design) before being applied to these data."
  paste(base, cvt, eng, gof, val)
}
