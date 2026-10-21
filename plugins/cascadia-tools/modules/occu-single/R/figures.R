# occu-single: publication figure suite (docs/FIGURES.md).
#
# occ_figures(fit, data, opts) returns a named list of ggplots built with the
# core constructors (core/R/figures.R). Figure ids:
#   occupancy_<cov>      conditional effect on (initial) occupancy
#   detection_<cov>      conditional effect on per-survey detection
#   col_ext_<cov>        dynamic: colonisation and/or extinction vs covariate
#   coefficients         odds ratios of every slope, grouped by process
#   detection_by_survey  per-survey detection probability by visit number
#   occupancy_map        only when site coordinates are available
#   trajectory           dynamic: proportion of surveyed sites occupied by season
#   gof                  MacKenzie-Bailey bootstrap (unmarked) / PPC (nimble)
#   trace, posterior     nimble only
#
# Units: occ_predict() takes covariates in the units of the researcher's
# files (the module standardises inside the design matrices only), so x axes
# are already in original units. fig_x() is still applied so a user `units:`
# entry can back-transform covariates that were standardised before intake.

occ_fig_procs <- c("psi", "p", "col", "ext")

occ_fig_proc_name <- function(s, model, opts) {
  def <- c(psi = if (model == "dynamic") "Initial occupancy" else "Occupancy", p = "Detection",
           col = "Colonisation", ext = "Extinction")[[s]]
  fig_label(paste0("process_", s), def, opts)
}

occ_fig_ylab <- function(s, model, opts) {
  key <- c(psi = "occupancy", p = "detection", col = "colonisation", ext = "extinction")[[s]]
  def <- c(psi = if (model == "dynamic") "Initial occupancy probability" else "Occupancy probability",
           p = "Detection probability per survey", col = "Colonisation probability",
           ext = "Extinction probability")[[s]]
  fig_label(key, def, opts)
}

# Fixed colour per process, so a process keeps its colour in every figure.
# Okabe-Ito positions 1, 2, 3, 6 (orange, sky blue, green, vermillion): the
# yellow at position 4 is too light for lines on white.
occ_fig_proc_cols <- function(opts) {
  pal <- if (length(opts$palette) > 1) unlist(opts$palette) else eco_pal[[opts$palette %||% "okabe_ito"]]
  if (is.null(pal)) pal <- eco_pal$okabe_ito
  if (identical(opts$palette, "grey"))   # the 4th grey (#d9d9d9) is invisible as a line/band on white;
    return(stats::setNames(c("#1a1a1a", "#666666", "#1a1a1a", "#808080"), occ_fig_procs))  # col/ext also differ by linetype
  idx <- if (identical(opts$palette, "okabe_ito") || is.null(opts$palette)) c(1, 2, 3, 6) else 1:4
  stats::setNames(rep_len(pal, max(idx))[idx], occ_fig_procs)
}

# Readable default name for a covariate column, with the unit taken from a
# common suffix: forest_cover_pct -> "Forest cover (%)".
occ_fig_pretty_var <- function(v) {
  units <- c(pct = "%", percent = "%", m = "m", km = "km", ha = "ha", km2 = "km²", mm = "mm", cm = "cm",
             c = "°C", degc = "°C", min = "min", h = "h", hr = "h", d = "days", days = "days",
             yr = "years", years = "years")
  parts <- strsplit(v, "_", fixed = TRUE)[[1]]
  u <- NULL
  if (length(parts) > 1 && tolower(parts[length(parts)]) %in% names(units)) {
    u <- units[[tolower(parts[length(parts)])]]; parts <- parts[-length(parts)]
  }
  s <- paste(parts, collapse = " ")
  s <- paste0(toupper(substr(s, 1, 1)), substring(s, 2))
  if (!is.null(u)) sprintf("%s (%s)", s, u) else s
}

occ_fig_xlab <- function(v, opts) fig_label(v, occ_fig_pretty_var(v), opts)

occ_fig_q <- function(P, prob) {
  a <- (1 - prob) / 2
  cbind(apply(P, 1, stats::quantile, a, na.rm = TRUE), apply(P, 1, stats::quantile, 1 - a, na.rm = TRUE))
}

# Prediction frame for one submodel over a covariate (raw units), with the
# outer band at opts$ci and, for Bayesian fits, the inner band at opts$inner_ci.
occ_fig_pred <- function(fit, s, nd, ci, inner) {
  pd <- occ_predict(fit, s, nd, prob = ci)
  out <- pd[, c("estimate", "lower", "upper")]
  if (!is.null(inner)) {
    pi <- occ_predict(fit, s, nd, prob = inner)
    out$lower50 <- pi$lower; out$upper50 <- pi$upper
  }
  out
}

occ_fig_itype <- function(o) if (o$engine == "nimble") "credible" else "confidence"

occ_fig_band_text <- function(o, ci, inner) {
  if (o$engine == "nimble")
    sprintf("bands: %d%% (light) and %d%% (dark) posterior credible intervals", round(100 * ci), round(100 * inner))
  else sprintf("band: %d%% confidence interval from %d draws of the asymptotic normal distribution of the maximum-likelihood estimates",
               round(100 * ci), nrow(o$draws))
}

# Draws (rows) x coefficient columns -> reported quantities with public names
# (intercepts as probabilities, slopes on the logit scale).
occ_fig_public_draws <- function(des, D, opts = NULL) {
  inames <- occ_intercept_names(des$model)
  out <- list()
  for (s in des$subs) {
    nm <- paste0("beta_", s, "[", seq_along(des$colnames[[s]]), "]")
    pn <- occ_fig_proc_name(s, des$model, opts)
    out[[fig_label(inames[[s]], paste(pn, "probability"), opts)]] <- stats::plogis(D[, nm[1]])
    if (length(nm) > 1) for (k in 2:length(nm)) {
      cc <- des$colnames[[s]][k]
      out[[fig_label(occ_coef_public(s, cc), sprintf("%s slope: %s", pn, occ_fig_term_label(des, cc, opts)), opts)]] <- D[, nm[k]]
    }
  }
  as.data.frame(out, check.names = FALSE)
}

# Label for one design-matrix column of submodel s (slopes only).
occ_fig_term_label <- function(des, col, opts) {
  for (v in names(des$levels)) {
    lv <- des$levels[[v]]
    if (col %in% paste0(v, lv)) {
      l <- substring(col, nchar(v) + 1)   # level labels follow the same <var>_<level> keys as the effect plots
      return(sprintf("%s: %s vs %s", occ_fig_xlab(v, opts), fig_label(paste0(v, "_", l), l, opts),
                     fig_label(paste0(v, "_", lv[1]), lv[1], opts)))
    }
  }
  if (col %in% names(des$raw)) return(occ_fig_xlab(col, opts))
  # transformed term or interaction: put readable names inside, log(area_ha) -> log(Area (ha))
  for (v in names(des$raw)) col <- gsub(paste0("\\b", v, "\\b"), occ_fig_xlab(v, opts), col, perl = TRUE)
  gsub(":", " \u00d7 ", col, fixed = TRUE)
}

occ_coord_cols <- function(sc, opts) {
  if (!is.null(opts$coords)) {
    cc <- unlist(opts$coords)
    if (length(cc) == 2 && all(cc %in% names(sc))) return(cc)
  }
  cand <- list(c("x", "y"), c("X", "Y"), c("easting", "northing"), c("utm_x", "utm_y"), c("utm_e", "utm_n"),
               c("longitude", "latitude"), c("lon", "lat"), c("long", "lat"))
  for (cc in cand) {
    m <- match(tolower(cc), tolower(names(sc)))
    if (!anyNA(m) && all(vapply(sc[m], is.numeric, TRUE))) return(names(sc)[m])
  }
  NULL
}

occ_figures <- function(fit, data, opts = fig_opts()) {
  o <- fit$object; des <- o$des; est <- fit$estimates
  figs <- list()
  if (o$engine == "naive") return(figs)   # negative control: no publication figures
  ci <- opts$ci %||% 0.95
  bayes <- o$engine == "nimble"
  inner <- if (bayes) (opts$inner_ci %||% 0.5) else NULL
  itype <- occ_fig_itype(o)
  band_txt <- occ_fig_band_text(o, ci, inner)
  cols <- occ_fig_proc_cols(opts)
  pname <- function(s) occ_fig_proc_name(s, des$model, opts)
  proc_scales <- function(ss, linetype = FALSE) {
    v <- stats::setNames(cols[ss], vapply(ss, pname, ""))
    c(list(ggplot2::scale_colour_manual(values = v, aesthetics = c("colour", "fill"), name = NULL)),
      if (linetype) list(ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22", "solid", "42")[match(ss, occ_fig_procs)],
                                                                 vapply(ss, pname, "")), name = NULL)))
  }
  show <- isTRUE(opts$show_data)
  # What the other covariates of the plotted submodel(s) are held at, with values,
  # e.g. ", holding Elevation (m) at its mean (834) and Observer at its reference level (AB)".
  held_text <- function(ss, v) {
    others <- setdiff(unique(unlist(lapply(ss, function(s) all.vars(stats::as.formula(des$forms[[s]]))))), v)
    if (!length(others)) return("")
    parts <- vapply(others, function(u) {
      if (!is.null(des$levels[[u]])) {
        l1 <- des$levels[[u]][1]
        sprintf("%s at its reference level (%s)", occ_fig_xlab(u, opts), fig_label(paste0(u, "_", l1), l1, opts))
      } else {
        m <- des$scaling[[u]][["mean"]] %||% mean(des$raw[[u]], na.rm = TRUE)
        sprintf("%s at its mean (%s)", occ_fig_xlab(u, opts), format(signif(fig_x(m, u, opts), 3), big.mark = ","))
      }
    }, "")
    paste0(", holding ", paste(parts, collapse = " and "))
  }
  # Captions name the axes as drawn, including a per_figure axis override.
  cap_lab <- function(id, k, def) opts$per_figure[[id]]$axis[[k]] %||% def

  # ---- conditional effects ------------------------------------------------
  one_effect <- function(s, v, id) {
    xlab <- occ_fig_xlab(v, opts); ylab <- occ_fig_ylab(s, des$model, opts)
    cx <- cap_lab(id, "x", xlab); cy <- cap_lab(id, "y", ylab)
    if (!is.null(des$levels[[v]])) {
      lv <- des$levels[[v]]
      pd <- occ_fig_pred(fit, s, stats::setNames(data.frame(lv), v), ci, inner)
      pd$level <- factor(vapply(lv, function(l) fig_label(paste0(v, "_", l), l, opts), ""),
                         levels = vapply(lv, function(l) fig_label(paste0(v, "_", l), l, opts), ""))
      p <- fig_cat_effect(pd, xlab, ylab)
      # bars in the process colour and points in ink on top, as in detection_by_survey
      p$layers[[length(p$layers)]]$aes_params$colour <- cols[[s]]
      if (!is.null(inner)) p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50),
                                                            colour = cols[[s]], linewidth = 1.4)
      p <- p + ggplot2::geom_point(size = 1.8, colour = eco_pal$ink)
      p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02)))
      cap <- sprintf("%s for each level of %s (points: estimates; bars: %d%%%s %s intervals)%s.",
                     cy, cx, round(100 * ci), if (bayes) sprintf(" and %d%%", round(100 * inner)) else "", itype, held_text(s, v))
      return(list(p = p, cap = cap))
    }
    x <- des$raw[[v]]
    grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
    pd <- occ_fig_pred(fit, s, stats::setNames(data.frame(grid), v), ci, inner)
    pd$x <- fig_x(grid, v, opts)
    pd$group <- pname(s)
    obs <- if (show) data.frame(x = fig_x(x[is.finite(x)], v, opts)) else NULL
    p <- fig_cond_effect(pd, xlab, ylab, obs = obs, ylim01 = TRUE) + proc_scales(s) +
      ggplot2::guides(colour = "none", fill = "none")
    cap <- sprintf("%s versus %s (line: estimate; %s)%s.%s Curves span the observed range of %s.%s",
                   cy, cx, band_txt, held_text(s, v), if (show) sprintf(" Rug: observed %s values%s.", cx,
                                                                 if (v %in% des$sources$obs) " for surveys done" else "") else "",
                   cx, if (s == "p") " This is an effect on detecting the species where it is present, not on where it occurs." else "")
    list(p = p, cap = cap)
  }
  for (s in intersect(c("psi", "p"), des$subs)) {
    pre <- c(psi = "occupancy", p = "detection")[[s]]
    for (v in all.vars(stats::as.formula(des$forms[[s]]))) {
      id <- paste0(pre, "_", occ_clean(v))
      r <- one_effect(s, v, id)
      figs[[id]] <- with_caption(r$p, r$cap)
    }
  }
  if (des$model == "dynamic") {
    vc <- all.vars(stats::as.formula(des$forms$col)); ve <- all.vars(stats::as.formula(des$forms$ext))
    for (v in union(vc, ve)) {
      ss <- c(if (v %in% vc) "col", if (v %in% ve) "ext")
      xlab <- occ_fig_xlab(v, opts)
      ylab <- fig_label("col_ext", if (length(ss) == 2) "Probability between seasons" else occ_fig_ylab(ss, des$model, opts), opts)
      if (!is.null(des$levels[[v]])) {
        lv <- des$levels[[v]]
        pd <- do.call(rbind, lapply(ss, function(s) {
          d <- occ_fig_pred(fit, s, stats::setNames(data.frame(lv), v), ci, inner)
          d$level <- factor(lv, levels = lv); d$group <- pname(s); d }))
        pdg <- ggplot2::position_dodge(width = 0.4)
        p <- ggplot2::ggplot(pd, ggplot2::aes(level, estimate, colour = group, shape = group)) +
          ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), linewidth = 0.6, position = pdg)
        if (!is.null(inner)) p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), linewidth = 1.3, position = pdg)
        p <- p + ggplot2::geom_point(size = 1.8, position = pdg) + proc_scales(ss) +
          ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17)[seq_along(ss)], vapply(ss, pname, "")), name = NULL) +
          ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
          ggplot2::labs(x = xlab, y = ylab)
        if (length(ss) == 1) p <- p + ggplot2::guides(colour = "none", shape = "none")
      } else {
        x <- des$raw[[v]]
        grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
        pd <- do.call(rbind, lapply(ss, function(s) {
          d <- occ_fig_pred(fit, s, stats::setNames(data.frame(grid), v), ci, inner)
          d$x <- fig_x(grid, v, opts); d$group <- pname(s); d }))
        pd$group <- factor(pd$group, levels = vapply(ss, pname, ""))
        obs <- if (show) data.frame(x = fig_x(x[is.finite(x)], v, opts)) else NULL
        p <- fig_cond_effect(pd, xlab, ylab, obs = obs, ylim01 = TRUE) +
          ggplot2::aes(linetype = group) + proc_scales(ss, linetype = TRUE)
        if (length(ss) == 1) p <- p + ggplot2::guides(colour = "none", fill = "none", linetype = "none")

      }
      what <- paste(vapply(ss, function(s) tolower(occ_fig_ylab(s, des$model, opts)), ""), collapse = " and ")
      cap <- sprintf("%s between consecutive seasons versus %s (%s; %s)%s. Colonisation is the probability that an unoccupied site becomes occupied by the next season; extinction is the probability that an occupied site becomes unoccupied. %s",
                     occ_cap(what), cap_lab(paste0("col_ext_", occ_clean(v)), "x", xlab), if (is.null(des$levels[[v]])) "lines: estimates" else "points: estimates", band_txt, held_text(ss, v),
                     if (v %in% des$sources$yearly) sprintf("%s in season t applies to the change from season t to t+1.", xlab) else "")
      figs[[paste0("col_ext_", occ_clean(v))]] <- with_caption(p, trimws(cap))
    }
  }

  # ---- coefficients (odds ratios, grouped by process) ----------------------
  z <- stats::qnorm(1 - (1 - ci) / 2)
  rows <- list()
  for (s in des$subs) {
    cn <- des$colnames[[s]]; if (length(cn) < 2) next
    nm <- paste0("beta_", s, "[", seq_along(cn), "]")[-1]
    if (bayes) {
      D <- o$draws[, nm, drop = FALSE]
      e <- colMeans(D); q <- t(occ_fig_q(t(D), ci)); qi <- t(occ_fig_q(t(D), inner))
      r <- data.frame(estimate = e, lower = q[1, ], upper = q[2, ], lower50 = qi[1, ], upper50 = qi[2, ])
    } else {
      e <- o$point[nm]; se <- if (!is.null(o$vcov)) sqrt(pmax(diag(o$vcov)[nm], 0)) else rep(NA_real_, length(nm))
      r <- data.frame(estimate = e, lower = e - z * se, upper = e + z * se)
    }
    pub <- occ_coef_public(s, cn[-1])
    r$param <- vapply(seq_along(pub), function(k) fig_label(pub[k], occ_fig_term_label(des, cn[k + 1], opts), opts), "")
    r$group <- pname(s); r$sub <- s
    rows[[s]] <- r
  }
  if (length(rows)) {
    cd <- do.call(rbind, rows); rownames(cd) <- NULL
    ss <- unique(cd$sub)
    cd$group <- factor(cd$group, levels = vapply(ss, pname, ""))
    # unique y positions per process so the same covariate in two processes stays separate
    cd$param <- paste0(cd$param, strrep("​", match(cd$sub, occ_fig_procs) - 1))
    cd$param <- factor(cd$param, levels = unique(cd$param))   # formula order, first term on top
    has_fac <- any(vapply(ss, function(s) any(des$colnames[[s]][-1] %in% unlist(lapply(names(des$levels), function(v) paste0(v, des$levels[[v]])))), TRUE))
    # factor rows already read "<level> vs <reference>", so the axis only needs the numeric unit
    # (a longer axis title is clipped at single-column width; the caption gives the per-SD detail)
    xdef <- if (has_fac) "Odds ratio" else sprintf("Odds ratio (%s)", if (isTRUE(o$settings$standardize)) "per SD" else "per unit")
    p <- fig_coef(cd[, setdiff(names(cd), "sub")], null_value = 1, transform = exp, log_x = TRUE, order = FALSE,
                  xlab = fig_label("coefficients_x", xdef, opts), opts = NULL) +
      ggplot2::facet_wrap(~group, ncol = 1, scales = "free_y") + proc_scales(ss) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::guides(colour = "none")
    sdv <- unlist(lapply(ss, function(s) {
      cs <- des$colnames[[s]][-1]; sdx <- des$colscale[[s]]["sd", cs]
      keep <- sdx != 1 | des$colscale[[s]]["mean", cs] != 0
      if (!any(keep)) return(NULL)
      sprintf("%s %s", vapply(cs[keep], function(cc) occ_fig_term_label(des, cc, opts), ""), signif(sdx[keep], 3))
    }))
    attr(p, "fig_size") <- list(height = 22 + 8 * nrow(cd) + 8 * length(ss))
    figs$coefficients <- with_caption(p, sprintf(
      "Odds ratios for every covariate effect, grouped by process; factor levels are relative to the reference level (points: estimates; %s %d%% %s intervals%s; log scale). Dashed line: no effect (odds ratio 1). %s%s",
      if (bayes) "thin bars:" else "bars:", round(100 * ci), itype,
      if (bayes) sprintf("; thick bars: %d%%", round(100 * inner)) else "",
      if (isTRUE(o$settings$standardize)) "Numeric covariates were standardised, so each odds ratio is per 1 SD of the covariate" else "Odds ratios are per unit of the covariate",
      if (length(sdv) && isTRUE(o$settings$standardize)) sprintf(" (1 SD: %s).", paste(unique(sdv), collapse = "; ")) else "."))
  }

  # ---- detection by survey ---------------------------------------------------
  nd <- min(nrow(o$draws), 1000)
  idx <- unique(round(seq(1, nrow(o$draws), length.out = nd)))
  nmp <- paste0("beta_p[", seq_along(des$colnames$p), "]")
  Xp <- matrix(des$Xp, ncol = length(nmp))
  Pd <- stats::plogis(Xp %*% t(o$draws[idx, nmp, drop = FALSE]))   # cells x draws
  wv <- as.vector(des$w) == 1
  cell <- expand.grid(i = seq_len(des$M), j = seq_len(des$J), t = seq_len(des$T))
  key <- paste(cell$j, cell$t)
  ukeys <- unique(key[wv])
  pbar <- t(vapply(ukeys, function(k) colMeans(Pd[wv & key == k, , drop = FALSE]), numeric(length(idx))))
  pt_est <- if (!is.null(o$point)) {
    pp <- stats::plogis(Xp %*% o$point[nmp])
    vapply(ukeys, function(k) mean(pp[wv & key == k]), 1)
  } else rowMeans(pbar)
  q <- occ_fig_q(pbar, ci)
  ds <- data.frame(survey = as.integer(sub(" .*", "", ukeys)), t = as.integer(sub(".* ", "", ukeys)),
                   estimate = pt_est, lower = q[, 1], upper = q[, 2],
                   n = vapply(ukeys, function(k) sum(wv & key == k), 1))
  if (bayes) { qi <- occ_fig_q(pbar, inner); ds$lower50 <- qi[, 1]; ds$upper50 <- qi[, 2] }
  ds <- ds[ds$n > 0, ]
  ds$season <- factor(des$season_labels[ds$t], levels = des$season_labels)
  p <- ggplot2::ggplot(ds, ggplot2::aes(survey, estimate)) +
    ggplot2::geom_line(colour = cols[["p"]], linewidth = 0.6) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = cols[["p"]], linewidth = 0.6)
  if (bayes) p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = cols[["p"]], linewidth = 1.4)
  p <- p + ggplot2::geom_point(size = 1.8, colour = eco_pal$ink) +
    ggplot2::scale_x_continuous(breaks = seq_len(des$J), minor_breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = fig_label("survey", "Survey (visit number)", opts), y = occ_fig_ylab("p", des$model, opts))
  if (des$model == "dynamic") p <- p + ggplot2::facet_wrap(~season, nrow = if (des$T > 4) 2 else 1) +
    ggplot2::labs(x = fig_label("survey", "Survey (visit number within season)", opts))
  figs$detection_by_survey <- with_caption(p, sprintf(
    "Per-survey detection probability by visit number%s, averaged over the sites surveyed on that visit using each survey's own covariate values (points: estimates; %s %d%% %s intervals%s). Surveys per visit number: %s.",
    if (des$model == "dynamic") " within each season" else "", if (bayes) "thin bars:" else "bars:", round(100 * ci), itype,
    if (bayes) sprintf("; thick bars: %d%%", round(100 * inner)) else "",
    paste(if (des$model == "dynamic") sprintf("%s/%d %d", ds$season, ds$survey, ds$n) else sprintf("visit %d: %d", ds$survey, ds$n), collapse = ", ")))

  # ---- occupancy map (only with coordinates) ---------------------------------
  cc <- occ_coord_cols(data$site_covs, opts)
  if (!is.null(cc)) {
    nmap <- min(nrow(o$draws), 300)
    im <- unique(round(seq(1, nrow(o$draws), length.out = nmap)))
    post <- Reduce(`+`, lapply(im, function(k)
      occ_forward_backward(des, occ_probs(des, occ_split_coef(o$draws[k, ], des)))$post)) / length(im)
    y0 <- des$y; y0[des$w == 0] <- 0
    det <- apply(y0, c(1, 3), max) > 0
    md <- do.call(rbind, lapply(seq_len(des$T), function(t) data.frame(
      x = data$site_covs[[cc[1]]], y = data$site_covs[[cc[2]]], value = post[, t],
      detected = factor(ifelse(det[, t], fig_label("map_detected", "Detected", opts), fig_label("map_not_detected", "Not detected", opts)),
                        levels = c(fig_label("map_detected", "Detected", opts), fig_label("map_not_detected", "Not detected", opts))),
      season = des$season_labels[t])))
    md$season <- factor(md$season, levels = des$season_labels)
    p <- ggplot2::ggplot(md, ggplot2::aes(x, y)) +
      ggplot2::geom_point(ggplot2::aes(fill = value, shape = detected), size = 1.3, stroke = 0.25, colour = eco_pal$ink) +
      # interior breaks only, so tick labels of neighbouring facets cannot touch
      ggplot2::scale_x_continuous(breaks = function(l) { b <- pretty(l, 4); b[b > l[1] + 0.08 * diff(l) & b < l[2] - 0.08 * diff(l)] }) +
      ggplot2::scale_y_continuous(n.breaks = 4) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(21, 24), levels(md$detected)), name = NULL) +
      scale_fill_eco_seq(limits = c(0, 1)) + ggplot2::coord_equal() +
      ggplot2::labs(x = fig_label(cc[1], occ_fig_pretty_var(cc[1]), opts), y = fig_label(cc[2], occ_fig_pretty_var(cc[2]), opts),
                    fill = fig_label("occupancy_map_fill", "P(occupied)", opts))
    if (des$model == "dynamic") {
      p <- p + ggplot2::facet_wrap(~season, ncol = min(des$T, 3))
      attr(p, "fig_size") <- list(width = "double", height = 20 + 55 * ceiling(des$T / 3))
    }
    figs$occupancy_map <- with_caption(p, sprintf(
      "Estimated probability that each surveyed site was occupied%s, given its detection history (colour; averaged over %d %s). Circles: species detected at least once (probability 1); triangles: never detected.",
      if (des$model == "dynamic") " in each season" else "", length(im),
      if (bayes) "posterior draws" else "draws of the parameter estimates"))
  }

  # ---- trajectory (dynamic) ---------------------------------------------------
  if (des$model == "dynamic") {
    nm <- paste0("prop_sites_occupied_t", seq_len(des$T))
    pr <- o$prop
    q <- occ_fig_q(t(pr), ci)
    tr <- data.frame(time = seq_len(des$T), estimate = est$estimate[match(nm, est$param)], lower = q[, 1], upper = q[, 2])
    y0 <- des$y; y0[des$w == 0] <- 0
    obs <- data.frame(time = seq_len(des$T), y = colMeans(apply(y0, c(1, 3), max)))
    p <- fig_timeseries(tr, fig_label("trajectory", "Proportion of sites occupied", opts),
                        obs = if (show) obs else NULL, xlab = fig_label("season", "Season", opts))
    if (bayes) {
      qi <- occ_fig_q(t(pr), inner); tr$lower50 <- qi[, 1]; tr$upper50 <- qi[, 2]
      p <- p + ggplot2::geom_ribbon(data = tr, ggplot2::aes(time, ymin = lower50, ymax = upper50), inherit.aes = FALSE,
                                    fill = eco_pal$band, alpha = 0.25)
    }
    p <- p + ggplot2::geom_point(data = tr, ggplot2::aes(time, estimate), colour = eco_pal$band, size = 1.4) +
      ggplot2::scale_x_continuous(breaks = seq_len(des$T), labels = des$season_labels, minor_breaks = NULL) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02)))
    figs$trajectory <- with_caption(p, sprintf(
      "Proportion of the %d surveyed sites occupied in each season (line: estimate; %s %d%% %s interval%s)%s. This finite-sample proportion applies to the surveyed sites only.",
      des$M, if (bayes) "light band:" else "band:", round(100 * ci), itype,
      if (bayes) sprintf("; dark band: %d%%", round(100 * inner)) else "",
      if (show) "; black points: naive proportion of sites where the species was detected at least once, which is biased low because occupied sites can be missed" else ""))
  }

  # ---- goodness of fit --------------------------------------------------------
  g <- o$gof
  if (!is.null(g)) {
    if (!is.null(g$boot) && o$engine == "unmarked") {
      p <- fig_ppc(g$boot, g$stat, fig_label("gof_x", "MacKenzie-Bailey chi-square", opts)) + ggplot2::labs(caption = NULL)
      if (identical(opts$palette, "grey")) p$layers[[1]]$aes_params$fill <- "#bdbdbd"   # fig_ppc's tinted fill -> true grey
      figs$gof <- with_caption(p, sprintf(
        "Goodness of fit (MacKenzie & Bailey 2004): histogram of the chi-square statistic for %d data sets simulated from the fitted model and refitted (parametric bootstrap); vertical line: observed data (chi-square %.1f, c-hat %.2f, p = %.3f). A small p (< 0.05) or c-hat well above 1 means the detection histories vary more than the model allows.",
        g$nboot, g$stat, g$c_hat, g$p_value))
    } else if (!is.null(g$rep) && !is.null(g$obs)) {
      p <- fig_ppc(g$rep - g$obs, 0, fig_label("gof_x", "Replicated minus observed chi-square", opts)) + ggplot2::labs(caption = NULL)
      if (identical(opts$palette, "grey")) p$layers[[1]]$aes_params$fill <- "#bdbdbd"
      figs$gof <- with_caption(p, sprintf(
        "Posterior predictive check of the MacKenzie-Bailey chi-square statistic: histogram of the difference between the statistic for data replicated from the model and for the observed data, over %d posterior draws; vertical line at 0. Bayesian p = %.3f is the proportion of draws to the right of the line; values near 0 indicate lack of fit.",
        g$nboot, g$p_value))
    }
  }

  # ---- MCMC trace and posterior (nimble) ---------------------------------------
  if (bayes) {
    nch <- o$settings$nchains %||% 1
    pdw <- occ_fig_public_draws(des, o$draws, opts)
    n <- nrow(pdw) %/% nch
    ch <- lapply(seq_len(nch), function(k) as.matrix(pdw[(k - 1) * n + seq_len(n), , drop = FALSE]))
    np <- ncol(pdw)
    p <- fig_trace(ch) + ggplot2::facet_wrap(~ind, scales = "free_y", ncol = 2, labeller = ggplot2::label_wrap_gen(40)) +
      ggplot2::scale_colour_manual(values = rep_len(eco_categorical(min(nch, 8), opts$palette), nch), name = "Chain")
    attr(p, "fig_size") <- list(width = "double", height = 18 + 22 * ceiling(np / 2))
    figs$trace <- with_caption(p, sprintf(
      "MCMC traces of the reported parameters after burn-in (%d chains, one colour per chain): intercepts as probabilities at average covariate values, slopes on the logit scale per 1 SD. Chains should overlap with no trend.", nch))
    p <- fig_posterior(as.matrix(pdw), ncol = 3) +
      ggplot2::facet_wrap(~param, scales = "free", ncol = 3, labeller = ggplot2::label_wrap_gen(28))
    attr(p, "fig_size") <- list(width = "double", height = 18 + 34 * ceiling(np / 3))
    figs$posterior <- with_caption(p, "Posterior densities of the reported parameters: intercepts as probabilities at average covariate values (factors at reference level), slopes on the logit scale per 1 SD.")
  }
  figs
}
