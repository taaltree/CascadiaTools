# occu-rota: publication figure suite (docs/FIGURES.md).
# Required ids: conditional_occupancy, marginal_occupancy_<cov>,
# interaction_<cov>, co_occurrence_matrix, detection, coefficients.
# Extras: conditional_<A>_given_<B>_<cov>, conditional_occupancy_average,
# detection_<cov>, gof, and trace / posterior for the Nimble engine.
# Covariates are drawn in original units (predictions are made on raw-unit
# grids and pass through the same formula transforms and scaling as the fit).

rota_fig_state <- function(fit, newdata, opts) {
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  pr <- rota_predict_state(fit, newdata, ci)
  if (fit$object$engine == "nimble" && !is.null(opts$inner_ci)) {
    p2 <- rota_predict_state(fit, newdata, opts$inner_ci)
    pr$lower50 <- p2$lower; pr$upper50 <- p2$upper
  }
  pr
}

rota_fig_det <- function(fit, s, newdata, opts) {
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  d <- rota_predict_det(fit, s, newdata, ci)
  if (fit$object$engine == "nimble" && !is.null(opts$inner_ci)) {
    d2 <- rota_predict_det(fit, s, newdata, opts$inner_ci)
    d$lower50 <- d2$lower; d$upper50 <- d2$upper
  }
  d
}

rota_figures <- function(fit, data, opts = fig_opts()) {
  o <- fit$object; des <- o$des
  if (o$engine == "naive") return(list())
  sp <- des$species
  spl <- function(s) fig_label(s, rota_lab(s), opts)
  cvl <- function(v) fig_label(v, rota_lab(v), opts)
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  bayes <- o$engine == "nimble"
  band <- sprintf("%d%% %s", round(100 * ci), if (bayes) "credible band (inner band 50%)" else
    "confidence band (from the asymptotic normal distribution of the maximum-likelihood estimates)")
  rug <- if (isTRUE(opts$show_data)) " Rug: covariate values at the surveyed sites." else ""
  out <- list()
  numeric_state <- Filter(function(v) is.null(des$levels[[v]]), rota_state_vars(des))
  factor_state <- Filter(function(v) !is.null(des$levels[[v]]), rota_state_vars(des))
  pairs <- utils::combn(des$S, 2, simplify = FALSE)
  has_pair <- function(a, b) paste0(a, ":", b) %in% names(des$terms) || paste0(b, ":", a) %in% names(des$terms)
  obs_for <- function(v) if (isTRUE(opts$show_data)) data.frame(x = fig_x(des$raw[[v]], v, opts)) else NULL
  # user-facing words, all relabellable via figures: labels:
  L <- function(key, default) fig_label(key, default, opts)
  det_ax <- L("detection_axis", "Detection probability per survey")
  w_present <- L("present", "present"); w_absent <- L("absent", "absent")
  # Colour follows the entity: species keep Okabe-Ito 1-3 (the core default
  # order); the state of the OTHER species is not a species, so it gets neutral
  # ink (present) / grey (absent) and a black-grey marginal, never a species hue.
  # Species pairs use Okabe-Ito 7, 5, 6 (reddish purple first: it is furthest from the species hues). Greyscale: the grey palette.
  grey_mode <- identical(opts$palette, "grey") || isTRUE(opts$greyscale)
  state_cols <- function(lv) {
    v <- ifelse(grepl(paste0(w_absent, "$"), lv) | lv == L("other_absent", "Other species absent"), "#9a9a9a",
         ifelse(lv == L("marginal", "Marginal"), "#5d6762", eco_pal$ink))
    stats::setNames(v, lv)
  }
  state_scale <- function(lv) list(ggplot2::scale_colour_manual(values = state_cols(lv), aesthetics = c("colour", "fill")))
  pair_scale <- function(lv) {
    pal <- if (grey_mode) eco_pal$grey else eco_pal$okabe_ito[c(7, 5, 6)]
    list(ggplot2::scale_colour_manual(values = stats::setNames(rep_len(pal, length(lv)), lv), aesthetics = c("colour", "fill")))
  }
  sp_factor <- function(x) factor(x, levels = vapply(sp, spl, ""))
  # covariates that enter a formula as log(v): draw them on a log x axis
  logx <- function(v) any(grepl(paste0("log\\(\\s*", v, "\\s*\\)"), c(unlist(des$forms), unlist(des$dforms))))
  # plain-number labels (no 1e+01) and at most ~4 log breaks so labels do not collide in narrow panels
  plain <- function(b) format(b, scientific = FALSE, drop0trailing = TRUE, trim = TRUE, big.mark = ",")
  xs <- function(p, v) if (logx(v) && all(des$raw[[v]] > 0))
    p + ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 4), labels = plain) else p
  # prediction grid: evenly spaced on the axis that is drawn (log-spaced for log covariates,
  # otherwise a linear grid on a log axis gives visibly kinked curves at the low end)
  grid_for <- function(v, n = 80) {
    x <- des$raw[[v]]
    if (logx(v) && all(x > 0, na.rm = TRUE)) {
      r <- range(log(x), na.rm = TRUE)
      return(stats::setNames(data.frame(exp(seq(r[1], r[2], length.out = n))), v))
    }
    rota_grid(des, v, n)
  }
  # what the curves hold constant, in the user's labels and display units
  held <- function(vars, prefix = "other covariates held at ") {
    vars <- setdiff(vars, NULL)
    if (!length(vars)) return("")
    tx <- vapply(vars, function(v) {
      if (!is.null(des$levels[[v]])) return(sprintf("%s = %s (reference level)", cvl(v), des$levels[[v]][1]))
      h <- if (!is.null(des$hold[[v]])) des$hold[[v]] else des$means[[v]]
      hv <- fig_x(h, v, opts)
      sprintf("%s = %s (%s)", cvl(v), format(signif(hv, 3), big.mark = ","),
              if (logx(v)) "geometric mean, as the model uses its logarithm" else "mean")
    }, "")
    paste0(prefix, paste(tx, collapse = "; "))
  }
  ivl <- if (bayes) "credible" else "confidence"
  # size hint that yields to the user: a per_figure aspect (without height) replaces the module's height
  size_hint <- function(name, width, height = NULL, aspect = NULL) {
    u <- opts$per_figure[[name]]
    if (!is.null(u$aspect) && is.null(u$height)) height <- NULL
    list(width = width, height = height, aspect = aspect)
  }
  # Multi-panel composite that stays readable at the width the user asks for.
  # Panels are styled with the composite's own per_figure options; below ~120 mm
  # they stack in one column with legends underneath. The result is wrapped as a
  # single element so the export-time theme cannot re-style only the last panel.
  compose_fit <- function(raw, name, default_width, panel_h, legend_last_only = FALSE) {
    u <- opts$per_figure[[name]]
    wmm <- fig_width_mm(if (!is.null(u$width)) u$width else default_width, opts)
    narrow <- length(raw) > 1 && wmm < 120
    po <- opts
    if (narrow && is.null(u$legend)) po$legend <- "bottom"
    pp <- lapply(raw, function(p) {
      q <- apply_fig_opts(p, name, po)
      if (narrow) {
        yl <- q$labels$y
        if (is.character(yl) && nchar(yl) > 26) q <- q + ggplot2::labs(y = paste(strwrap(yl, 26), collapse = "\n"))
        # legend under the panel, aligned to the figure's left edge so long names use the full width
        q <- q + ggplot2::theme(legend.direction = "vertical", legend.justification = "left", legend.location = "plot")
      }
      q
    })
    if (legend_last_only && length(pp) > 1) for (k in seq_len(length(pp) - 1))
      pp[[k]] <- pp[[k]] + ggplot2::theme(legend.position = "none")
    nc <- if (narrow) 1 else min(2, length(pp))
    pc <- compose_figures(pp, ncol = nc, guides = "keep")
    pc <- patchwork::wrap_elements(full = pc)
    # stacked panels need their height: a per_figure aspect only applies to the side-by-side layout
    attr(pc, "fig_size") <- if (narrow) list(width = default_width, height = (panel_h + 18) * length(pp)) else
      size_hint(name, default_width, panel_h * ceiling(length(pp) / nc))
    pc
  }

  # ---- conditional occupancy at average conditions (dot plot)
  pr0 <- rota_fig_state(fit, data.frame(), opts)
  cond <- pr0[grepl("_given_", pr0$quantity), ]
  cond$target <- sub("^prob_(.*)_given_.*$", "\\1", cond$quantity)
  cond$other <- sub("^prob_.*_given_(.*)_(present|absent)$", "\\1", cond$quantity)
  cond$status <- ifelse(grepl("_present$", cond$quantity), L("other_present", "Other species present"), L("other_absent", "Other species absent"))
  cond$lab <- paste0(vapply(cond$target, spl, ""), " | ", vapply(cond$other, spl, ""))
  marg <- pr0[grepl("^occupancy_prob_", pr0$quantity), ]
  marg$target <- sub("^occupancy_prob_", "", marg$quantity); marg$other <- ""
  marg$status <- L("marginal", "Marginal"); marg$lab <- paste0(vapply(marg$target, spl, ""), " (", tolower(L("marginal", "Marginal")), ")")
  cc <- rbind(cond, marg)
  cc$status <- factor(cc$status, levels = c(L("marginal", "Marginal"), L("other_absent", "Other species absent"), L("other_present", "Other species present")))
  cc$lab <- factor(cc$lab, levels = rev(unique(c(marg$lab, sort(cond$lab)))))
  pavg <- ggplot2::ggplot(cc, ggplot2::aes(estimate, lab, colour = status, shape = status)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.5,
                            position = ggplot2::position_dodge(width = 0.5))
  if (bayes) pavg <- pavg + ggplot2::geom_linerange(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.3,
                                                    position = ggplot2::position_dodge(width = 0.5))
  pavg <- pavg + ggplot2::geom_point(size = 1.8, position = ggplot2::position_dodge(width = 0.5)) +
    state_scale(levels(cc$status)) + ggplot2::scale_shape_manual(values = c(1, 17, 16)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    ggplot2::labs(x = fig_label("occupancy_axis", "Occupancy probability", opts), y = NULL, colour = NULL, shape = NULL) +
    ggplot2::scale_y_discrete(labels = scales::label_wrap(24))
  attr(pavg, "fig_size") <- size_hint("conditional_occupancy_average", "onehalf", 30 + 7 * nlevels(cc$lab))
  out$conditional_occupancy_average <- with_caption(pavg, sprintf(
    "Probability that each species occupies a site at average covariate values (%s): marginal (any state of the other species) and conditional on another species being absent or present ('A | B' = A given B). Points: estimates; bars: %d%% %s intervals.%s",
    held(rota_state_vars(des), ""), round(100 * ci), if (bayes) "credible (thick: 50%)" else "confidence",
    if (des$S == 3) " With three species, conditional probabilities average over the third species." else ""))

  # ---- conditional occupancy vs each numeric covariate
  cond_panels <- list()
  for (v in numeric_state) {
    g <- grid_for(v)
    hs <- held(setdiff(rota_state_vars(des), v))
    pr <- rota_fig_state(fit, g, opts); pr$x <- fig_x(g[[v]][pr$row], v, opts)
    # marginal occupancy: one curve per species
    d <- pr[grepl("^occupancy_prob_", pr$quantity), ]
    d$group <- sp_factor(vapply(sub("^occupancy_prob_", "", d$quantity), spl, ""))
    p <- fig_cond_effect(d[, intersect(c("x", "estimate", "lower", "upper", "lower50", "upper50", "group"), names(d))],
                         cvl(v), fig_label("marginal_axis", "Occupancy probability (marginal)", opts),
                         obs = obs_for(v), ylim01 = TRUE) + ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL)
    p <- xs(p, v)
    out[[paste0("marginal_occupancy_", rota_clean(v))]] <- with_caption(p, sprintf(
      "Marginal occupancy probability of each species (line) against %s (%s), with %s.%s",
      cvl(v), hs, band, rug))
    for (a in sp) for (b in sp) if (a != b && has_pair(a, b)) {
      q <- paste0("prob_", a, "_given_", b, c("_present", "_absent"))
      d <- pr[pr$quantity %in% q, ]
      d$group <- factor(paste(spl(b), ifelse(grepl("_present$", d$quantity), w_present, w_absent)),
                        levels = paste(spl(b), c(w_absent, w_present)))
      p <- fig_cond_effect(d[, intersect(c("x", "estimate", "lower", "upper", "lower50", "upper50", "group"), names(d))],
                           cvl(v), sprintf(L("conditional_axis", "P(%s present)"), spl(a)), obs = obs_for(v), ylim01 = TRUE) +
        ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL) + state_scale(levels(d$group)) +
        ggplot2::scale_linetype_manual(values = c(2, 1))
      p <- xs(p, v)
      nm <- paste0("conditional_", a, "_given_", b, "_", rota_clean(v))
      out[[nm]] <- with_caption(p, sprintf(
        "Probability that %s occupies a site where %s is present vs absent, against %s (%s), with %s. Diverging curves mean the co-occurrence association changes along the gradient; parallel curves on the logit scale mean it does not. Co-occurrence is not evidence of an ecological interaction.%s",
        spl(a), spl(b), cvl(v), hs, band, rug))
      cond_panels[[v]][[nm]] <- p
    }
    # interaction natural parameter (as odds ratio) vs v, for pairs whose formula uses v
    dl <- list()
    for (cc2 in pairs) {
      a <- sp[cc2[1]]; b <- sp[cc2[2]]; key <- paste0(a, ":", b)
      k <- match(key, names(des$terms))
      if (is.na(k) || !v %in% all.vars(stats::as.formula(des$forms[[key]]))) next
      B <- o$draws[, seq_len(des$P), drop = FALSE]
      Xn <- rota_Xnew(des, g)
      Fl <- rota_F_from(des, B, Xn)
      fk <- vapply(Fl, function(m) matrix(m, nrow(B))[, k], numeric(nrow(B)))   # draws x grid
      est <- if (!is.null(o$point)) vapply(rota_F_from(des, matrix(o$point[seq_len(des$P)], 1), Xn), function(m) m[k], 0) else colMeans(fk)
      qs <- function(pp) apply(fk, 2, stats::quantile, pp)
      dd <- data.frame(x = fig_x(g[[v]], v, opts), estimate = exp(est), lower = exp(qs((1 - ci) / 2)), upper = exp(qs(1 - (1 - ci) / 2)),
                       group = paste(spl(a), "&", spl(b)))
      if (bayes && !is.null(opts$inner_ci)) { dd$lower50 <- exp(qs((1 - opts$inner_ci) / 2)); dd$upper50 <- exp(qs(1 - (1 - opts$inner_ci) / 2)) }
      dl[[key]] <- dd
    }
    if (length(dl)) {
      dd <- do.call(rbind, dl)
      p <- fig_cond_effect(dd, cvl(v), fig_label("interaction_axis", "Co-occurrence odds ratio", opts),
                           obs = obs_for(v)) +
        ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
        ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 6), labels = plain)
      p <- xs(p, v)
      p <- p + pair_scale(unique(dd$group))
      if (length(dl) > 1) p <- p + ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL) else p <- p + ggplot2::guides(colour = "none", fill = "none")
      out[[paste0("interaction_", rota_clean(v))]] <- with_caption(p, sprintf(
        "Second-order natural parameter f_AB (%s) expressed as a co-occurrence odds ratio exp(f_AB) (log axis) against %s, with %s. Dashed line: odds ratio 1 (species occur independently). Above 1: the pair occurs together more often than expected under independence; below 1: less often.%s Co-occurrence is not evidence of an ecological interaction.%s",
        paste(unique(dd$group), collapse = "; "), cvl(v), band, if (des$S == 3) " With three species this is the odds ratio where the third species is absent." else "", rug))
    }
  }
  # factor site covariates: marginal and conditional occupancy by level
  for (v in factor_state) {
    nd <- stats::setNames(data.frame(des$levels[[v]]), v)
    pr <- rota_fig_state(fit, nd, opts); pr$level <- factor(nd[[v]][pr$row], levels = des$levels[[v]])
    d <- pr[grepl("^occupancy_prob_", pr$quantity), ]
    d$group <- sp_factor(vapply(sub("^occupancy_prob_", "", d$quantity), spl, ""))
    p <- ggplot2::ggplot(d, ggplot2::aes(level, estimate, colour = group, shape = group)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(0.4), size = 0.35) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = cvl(v), y = L("marginal_axis", "Occupancy probability (marginal)"), colour = NULL, shape = NULL)
    out[[paste0("marginal_occupancy_", rota_clean(v))]] <- with_caption(p, sprintf(
      "Marginal occupancy probability of each species for each level of %s (points; %s) with %d%% %s intervals.",
      cvl(v), held(setdiff(rota_state_vars(des), v)), round(100 * ci), ivl))
  }

  # conditional_occupancy: the headline figure. Panels for every ordered pair
  # along the covariate on the interaction (else the first numeric covariate);
  # with no numeric covariate, the average-conditions dot plot.
  int_vars <- unique(unlist(lapply(names(des$terms)[lengths(des$terms) >= 2], function(k) all.vars(stats::as.formula(des$forms[[k]])))))
  main_v <- intersect(c(int_vars, numeric_state), names(cond_panels))[1]
  if (!is.na(main_v)) {
    pp <- cond_panels[[main_v]]
    pc <- compose_fit(pp, "conditional_occupancy", "double", 70)
    out$conditional_occupancy <- with_caption(pc, sprintf(
      "Conditional occupancy against %s: each panel shows the probability that one species occupies a site where the other is present vs absent (%s), with %s. Co-occurrence is not evidence of an ecological interaction.",
      cvl(main_v), held(setdiff(rota_state_vars(des), main_v)), band))
  } else out$conditional_occupancy <- out$conditional_occupancy_average

  # ---- co-occurrence matrix: P(row | column present) - P(row | column absent)
  dm <- do.call(rbind, lapply(sp, function(a) do.call(rbind, lapply(sp, function(b) {
    if (a == b) return(data.frame(row = spl(a), col = spl(b), value = NA_real_, lab = ""))
    if (!has_pair(a, b)) return(data.frame(row = spl(a), col = spl(b), value = 0, lab = "fixed 0"))
    D <- rota_contrast_draws(fit, a, b)
    est <- D$est; lo <- stats::quantile(D$draws, (1 - ci) / 2); hi <- stats::quantile(D$draws, 1 - (1 - ci) / 2)
    data.frame(row = spl(a), col = spl(b), value = est,
               lab = sprintf("%+.2f\n(%+.2f,\n%+.2f)", est, lo, hi))
  }))))
  dm$row <- factor(dm$row, levels = vapply(sp, spl, "")); dm$col <- factor(dm$col, levels = vapply(sp, spl, ""))
  lim <- max(0.05, abs(dm$value), na.rm = TRUE)
  pm <- ggplot2::ggplot(dm, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = if (des$S == 3) 2.1 else 2.4, colour = eco_pal$ink, lineheight = 0.9) +
    scale_fill_eco_div(0, limits = c(-lim, lim), na.value = "grey92") +
    ggplot2::scale_y_discrete(limits = rev) + ggplot2::coord_equal() +
    ggplot2::labs(x = L("matrix_x", "Other species (present vs absent)"), y = L("matrix_y", "Occupancy of"),
                  fill = L("matrix_fill", "Difference in\noccupancy\nprobability"))
  attr(pm, "fig_size") <- list(width = if (des$S == 3) "onehalf" else "single", aspect = if (des$S == 3) 0.7 else 0.6)
  out$co_occurrence_matrix <- with_caption(pm, sprintf(
    "Change in each row species' occupancy probability where the column species is present vs absent, at average covariate values (%s): P(row | column present) - P(row | column absent), with %d%% %s intervals. Red: the species is more likely to occupy sites where the other is present; blue: less likely. Co-occurrence is not evidence of an ecological interaction.",
    held(rota_state_vars(des), ""), round(100 * ci), if (bayes) "credible" else "confidence"))

  # ---- detection
  det_panels <- list()
  for (v in rota_det_vars(des)) {
    sps <- Filter(function(s) rota_dforms_has(des, s, v), sp)
    if (!length(sps)) next
    if (!is.null(des$levels[[v]])) {
      nd <- stats::setNames(data.frame(des$levels[[v]]), v)
      d <- do.call(rbind, lapply(sps, function(s) cbind(level = factor(nd[[v]], levels = des$levels[[v]]),
                                                         rota_fig_det(fit, s, nd, opts), group = spl(s))))
      d$group <- sp_factor(d$group)
      p <- ggplot2::ggplot(d, ggplot2::aes(level, estimate, colour = group, shape = group, linetype = group)) +
        ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(0.4), size = 0.35) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = cvl(v), y = det_ax, colour = NULL, shape = NULL, linetype = NULL)
    } else {
      if (is.null(des$raw[[v]])) next
      g <- grid_for(v)
      d <- do.call(rbind, lapply(sps, function(s) cbind(x = fig_x(g[[v]], v, opts), rota_fig_det(fit, s, g, opts), group = spl(s))))
      d$group <- sp_factor(d$group)
      p <- fig_cond_effect(d, cvl(v), det_ax, obs = obs_for(v), ylim01 = TRUE) +
        ggplot2::aes(linetype = group) + ggplot2::labs(linetype = NULL)
      p <- xs(p, v)
    }
    out[[paste0("detection_", rota_clean(v))]] <- with_caption(p, sprintf(
      "Probability of detecting each species on one survey where it is present, against %s (%s), with %d%% %s intervals.",
      cvl(v), sub("other covariates", "other detection covariates", held(setdiff(rota_det_vars(des), v))), round(100 * ci), ivl))
    det_panels[[v]] <- p
  }
  if (length(det_panels)) {
    # one shared legend: keep it on the last panel only (line and point keys differ, so patchwork cannot merge them)
    pd <- if (length(det_panels) == 1) det_panels[[1]] else compose_fit(det_panels, "detection", "double", 75, legend_last_only = TRUE)
    out$detection <- with_caption(pd, sprintf("Detection probability per survey for each species against each detection covariate (in each panel, %s), with %d%% %s intervals.",
      sub("other covariates", "the other detection covariates", held(rota_det_vars(des))), round(100 * ci), ivl))
  } else {
    e <- fit$estimates[grepl("^detection_prob_", fit$estimates$param), ]
    e$group <- sp_factor(vapply(sub("^detection_prob_", "", e$param), spl, ""))
    p <- ggplot2::ggplot(e, ggplot2::aes(group, estimate, colour = group)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.35) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) + ggplot2::guides(colour = "none") +
      ggplot2::labs(x = NULL, y = det_ax)
    out$detection <- with_caption(p, sprintf("Detection probability per survey for each species with %d%% %s intervals.", round(100 * ci), ivl))
  }

  # ---- coefficients (odds ratios, log axis)
  est <- fit$estimates
  keep <- grepl("^(occ_.*_beta_|interaction_|p_.*_beta_)", est$param)
  ce <- est[keep, c("param", "estimate", "lower", "upper")]
  if (nrow(ce)) {
    pl <- c(L("process_occupancy", "Occupancy"), L("process_interaction", "Interaction"), L("process_detection", "Detection"))
    ce$group <- ifelse(grepl("^occ_", ce$param), pl[1], ifelse(grepl("^interaction_", ce$param), pl[2], pl[3]))
    ce$group <- factor(ce$group, levels = pl)
    ce$param <- vapply(ce$param, function(p) fig_label(p, rota_coef_label(p, sp, spl, des, cvl), opts), "")
    ce$process <- ce$group; ce$group <- NULL
    ce$param <- factor(ce$param, levels = rev(unique(ce$param)))
    pc <- fig_coef(ce, null_value = 1, transform = exp, log_x = TRUE, order = FALSE,
                   xlab = sprintf(L("coef_axis", "Odds ratio (%d%% interval, log scale)"), round(100 * ci))) +
      ggplot2::facet_grid(process ~ ., scales = "free_y", space = "free_y",
                          labeller = ggplot2::label_wrap_gen(14)) +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0)) +
      ggplot2::scale_y_discrete(labels = scales::label_wrap(34))
    attr(pc, "fig_size") <- size_hint("coefficients", "onehalf", 30 + 7 * nrow(ce))
    out$coefficients <- with_caption(pc, paste(
      "Link-scale coefficients shown as odds ratios (log axis); dashed line = 1 (no effect).",
      "Covariate effects are per 1 SD of the (transformed) covariate, or a factor level vs the reference level.",
      "Interaction intercepts are co-occurrence odds ratios at average covariate values; interaction slopes multiply that odds ratio per 1 SD.",
      sprintf("Bars: %d%% %s intervals.", round(100 * ci), if (bayes) "credible" else "Wald confidence")))
  }

  # ---- goodness of fit
  gof <- fit$diagnostics$gof
  if (!is.null(gof) && length(gof$sims) > 1)
    out$gof <- with_caption(fig_ppc(gof$sims, gof$stat, L("gof_axis", "Fit statistic (chi-square)")),
      sprintf("Goodness of fit (%s): histogram of the fit statistic (chi-square on the number of sites with 0, 1, 2, ... detections of each species + chi-square on the number of sites where each pair was detected together) for data simulated from the fitted model; line = observed data%s. p = %.2f, c-hat = %.2f; p < 0.05 indicates lack of fit.",
              gof$method, if (bayes) " (posterior mean)" else "", gof$p_value, gof$c_hat))

  # ---- MCMC diagnostics (Nimble)
  if (bayes && !is.null(o$chains)) {
    nm <- colnames(o$chains[[1]])
    out$trace <- with_caption(fig_trace(o$chains, nm),
      "MCMC traces of all coefficients (link scale), one colour per chain; well-mixed chains overlap without trends.")
    attr(out$trace, "fig_size") <- list(width = "double", height = 30 * ceiling(length(nm) / 3) + 20)
    out$posterior <- with_caption(fig_posterior(o$draws, nm, ncol = 3),
      "Posterior densities of all coefficients (link scale), chains pooled.")
    attr(out$posterior, "fig_size") <- list(width = "double", height = 30 * ceiling(length(nm) / 3) + 20)
  }
  out
}

# Readable coefficient label built from the user's labels:
# "Coyote: Forest cover (%)", "Coyote & Gray fox: Housing density (log)",
# "Coyote detection: Camera on trail = yes".
rota_coef_label <- function(p, sp, spl, des = NULL, cvl = function(v) gsub("_", " ", v)) {
  pair_of <- function(x) { for (a in sp) for (b in sp) if (a != b && x == paste0(a, "_", b)) return(paste(spl(a), "&", spl(b))); NULL }
  term_lab <- function(t) {
    vars <- unique(c(names(des$raw), names(des$levels)))
    vars <- vars[order(-nchar(vars))]
    for (v in vars) {
      cv <- rota_clean(v)
      if (t == cv) return(cvl(v))
      if (t %in% paste0(c("log_", "log"), cv)) return(paste0(cvl(v), " (log)"))
      if (t %in% paste0(c("I_", "poly_"), cv)) return(cvl(v))
      if (!is.null(des$levels[[v]]) && startsWith(t, cv)) {
        lev <- substring(t, nchar(cv) + 1)
        if (nzchar(lev)) return(paste0(cvl(v), " = ", lev))
      }
    }
    gsub("_", " ", t)
  }
  for (s in sp) {
    if (startsWith(p, paste0("occ_", s, "_beta_"))) return(paste0(spl(s), ": ", term_lab(substring(p, nchar(s) + 11))))
    if (startsWith(p, paste0("p_", s, "_beta_"))) return(paste0(spl(s), " detection: ", term_lab(substring(p, nchar(s) + 9))))
  }
  if (grepl("^interaction_", p)) {
    q <- sub("^interaction_", "", p)
    if (grepl("_beta_", q)) {
      pr <- sub("_beta_.*$", "", q); t <- sub("^.*?_beta_", "", q, perl = TRUE)
      pl <- pair_of(pr); if (is.null(pl)) pl <- gsub("_", " ", pr)
      return(paste0(pl, ": ", term_lab(t)))
    }
    pl <- pair_of(q); if (is.null(pl)) pl <- gsub("_", " ", q)
    return(paste0(pl, " (intercept)"))
  }
  gsub("_", " ", p)
}

# Draws (and point estimate) of P(a | b present) - P(a | b absent) at average covariates.
rota_contrast_draws <- function(fit, a, b) {
  o <- fit$object; des <- o$des
  qa <- paste0("prob_", a, "_given_", b, "_present"); qb <- paste0("prob_", a, "_given_", b, "_absent")
  B <- o$draws[, seq_len(des$P), drop = FALSE]
  Xn <- rota_Xnew(des, data.frame())
  D <- rota_derived(rota_psi(matrix(rota_F_from(des, B, Xn)[[1]], nrow(B)), des$A), des$Z)
  dr <- D[, qa] - D[, qb]
  est <- if (!is.null(o$point)) {
    Dp <- rota_derived(rota_psi(matrix(rota_F_from(des, matrix(o$point[seq_len(des$P)], 1), Xn)[[1]], 1), des$A), des$Z)
    Dp[1, qa] - Dp[1, qb]
  } else mean(dr)
  list(est = est, draws = dr)
}
