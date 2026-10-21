# occu-waddle: predictions, figures, interpretation, methods text.

wd_fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
wd_ci <- function(e, l, u, d = 2) sprintf("%s (95%% %s %s to %s)", wd_fmt(e, d), "%INT%", wd_fmt(l, d), wd_fmt(u, d))

# Interval word depends on the engine.
wd_int <- function(fit) if (fit$object$engine == "nimble") "CrI" else "CI"
wd_ci2 <- function(fit, e, l, u, d = 2) sub("%INT%", wd_int(fit), wd_ci(e, l, u, d), fixed = TRUE)
wd_row <- function(fit, n) fit$estimates[fit$estimates$param == n, , drop = FALSE]

# One-row template of covariates: numeric at the mean (of the standardised
# model terms, see below), factors at the reference level.
wd_template <- function(df, mms) {
  out <- df[1, , drop = FALSE]
  for (v in names(df)) {
    x <- df[[v]]
    if (is.numeric(x)) out[[v]] <- mean(x, na.rm = TRUE)
    else {
      lv <- NULL
      for (m in mms) if (!is.null(m$xlev[[v]])) lv <- m$xlev[[v]]
      out[[v]] <- if (!is.null(lv)) lv[1] else sort(unique(stats::na.omit(as.character(x))))[1]
    }
  }
  # "Average conditions" must match the estimates table (intercepts = all
  # standardised model terms at 0). For a transformed term such as log(x) the
  # raw mean is not that point (log(mean x) != mean log x), so each numeric
  # covariate is moved to the value at which the standardised columns that
  # depend on it are closest to 0 (the geometric mean for a log term).
  for (v in names(df)) {
    x <- df[[v]]
    if (!is.numeric(x) || all(is.na(x))) next
    rng <- range(x, na.rm = TRUE)
    if (!(rng[2] > rng[1])) next
    mmat <- function(m, val) { nd <- out; nd[[v]] <- val; tryCatch(wd_mm_new(m, nd), error = function(e) NULL) }
    # standardised columns of each model matrix that change with v
    cols <- lapply(mms, function(m) {
      if (!(v %in% m$vars)) return(integer())
      a <- mmat(m, rng[1]); b <- mmat(m, rng[2])
      if (is.null(a) || is.null(b)) return(integer())
      which(is.finite(a[1, ]) & is.finite(b[1, ]) & abs(a[1, ] - b[1, ]) > 1e-12)
    })
    if (!length(unlist(cols))) next
    h <- function(val) {
      s <- 0
      for (k in seq_along(mms)) if (length(cols[[k]])) s <- s + sum(mmat(mms[[k]], val)[1, cols[[k]]]^2)
      if (is.finite(s)) s else 1e10
    }
    h0 <- h(out[[v]])
    if (h0 > 1e-12) {
      o <- tryCatch(stats::optimize(h, rng, tol = 1e-9 * diff(rng)), error = function(e) NULL)
      if (!is.null(o) && o$objective < h0) out[[v]] <- o$minimum
    }
  }
  out
}

# Coefficient draws (rows) and an optional point estimate.
wd_B <- function(fit) {
  o <- fit$object
  list(D = o$draws, pt = o$point)
}

# Bayesian point estimate of a draws matrix (columns), matching the estimates
# table: posterior median unless settings$point = "mean".
wd_ctr <- function(fit, m) {
  if (identical(fit$object$settings$point, "mean")) colMeans(m) else apply(m, 2, stats::median)
}

# Probabilities for parameter set `pars` at newdata nd (site-level for
# occupancy, survey-level for detection). Returns list per param: matrix
# n_draws x nrow(nd), plus point vector.
wd_prob_draws <- function(fit, nd, pars) {
  des <- fit$object$des; B <- wd_B(fit)
  out <- list()
  for (p in pars) {
    f <- des$map[[p]]; X <- wd_mm_new(des$mms[[f]], nd); ii <- des$idx[[f]]
    out[[p]] <- stats::plogis(B$D[, ii, drop = FALSE] %*% t(X))
    attr(out[[p]], "point") <- if (!is.null(B$pt)) as.vector(stats::plogis(X %*% B$pt[ii])) else wd_ctr(fit, out[[p]])
  }
  out
}

wd_summ <- function(m, point = attr(m, "point"), ci = 0.95, inner = NULL) {
  a <- (1 - ci) / 2
  d <- data.frame(estimate = if (is.null(point)) colMeans(m) else point,
                  lower = apply(m, 2, stats::quantile, a, na.rm = TRUE),
                  upper = apply(m, 2, stats::quantile, 1 - a, na.rm = TRUE))
  if (!is.null(inner)) {
    b <- (1 - inner) / 2
    d$lower50 <- apply(m, 2, stats::quantile, b, na.rm = TRUE); d$upper50 <- apply(m, 2, stats::quantile, 1 - b, na.rm = TRUE)
  }
  d
}

# Occupancy quantities (psiA, psiBA, psiBa, psiB, SIF) at nd.
wd_occ_draws <- function(fit, nd) {
  q <- wd_prob_draws(fit, nd, c("psiA", "psiBA", "psiBa"))
  psiB <- q$psiA * q$psiBA + (1 - q$psiA) * q$psiBa
  sif <- q$psiBA / psiB
  pt <- function(p) attr(q[[p]], "point")
  pB <- pt("psiA") * pt("psiBA") + (1 - pt("psiA")) * pt("psiBa")
  if (is.null(fit$object$point)) {          # Bayesian: summarise the derived draws directly
    attr(psiB, "point") <- wd_ctr(fit, psiB); attr(sif, "point") <- wd_ctr(fit, sif)
  } else { attr(psiB, "point") <- pB; attr(sif, "point") <- pt("psiBA") / pB }
  c(q, list(psiB = psiB, sif = sif))
}

wd_vars <- function(des, pars) {
  v <- unique(unlist(lapply(unique(des$map[pars]), function(p) des$mms[[p]]$vars)))
  setdiff(v, ".dummy")
}

wd_grid <- function(x, n = 60) {
  if (is.numeric(x)) seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n)
  else sort(unique(stats::na.omit(as.character(x))))
}

wd_species <- function(fit) {
  sp <- fit$object$des$species
  if (is.null(sp)) sp <- c(dominant = "dominant", subordinate = "subordinate")
  sp
}

wd_cov_df <- function(des, v, occ) {
  if (occ) des$site_covs else des$surv_df[des$obs_v, , drop = FALSE]
}

# Data frame of curves for a covariate v over parameters `what`.
wd_curve <- function(fit, v, occ = TRUE, opts = fig_opts()) {
  des <- fit$object$des
  df <- wd_cov_df(des, v, occ)
  tmpl <- wd_template(df, des$mms)
  g <- wd_grid(df[[v]])
  # log-entered covariates: grid evenly spaced on the log scale (drawn on a log axis)
  if (is.numeric(g) && wd_is_log(des, v) && min(g) > 0) g <- exp(seq(log(min(g)), log(max(g)), length.out = length(g)))
  nd <- tmpl[rep(1, length(g)), , drop = FALSE]; nd[[v]] <- g
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  inner <- if (fit$object$engine == "nimble") opts$inner_ci else NULL
  sp <- wd_species(fit)
  if (occ) {
    q <- wd_occ_draws(fit, nd)
    labs <- wd_occ_labels(sp, opts)[c("psiA", "psiBA", "psiBa", "psiB", "sif")]
  } else {
    pars <- c("pA", "pB", "rA", "rBA", "rBa")
    pars <- pars[!duplicated(des$map[pars])]
    q <- wd_prob_draws(fit, nd, pars)
    labs <- wd_det_labels(sp, opts)[pars]
  }
  do.call(rbind, lapply(names(labs), function(p) {
    s <- wd_summ(q[[p]], ci = ci, inner = inner)
    data.frame(x = if (is.numeric(g)) g else factor(g, levels = g), level = g, s, group = wd_wrap(labs[[p]], 28), par = p)
  }))
}

# Display names of every quantity. Each goes through fig_label(), so a config
# `figures: labels:` entry keyed by the parameter name (psiA, psiBA, psiBa,
# psiB, sif, pA, pB, rA, rBA, rBa) renames it in every figure.
wd_occ_labels <- function(sp, opts = NULL) {
  A <- sp[["dominant"]]; B <- sp[["subordinate"]]
  c(psiA = fig_label("psiA", A, opts),
    psiB = fig_label("psiB", sprintf("%s (overall)", B), opts),
    psiBa = fig_label("psiBa", sprintf("%s | %s absent", B, A), opts),
    psiBA = fig_label("psiBA", sprintf("%s | %s present", B, A), opts),
    sif = fig_label("sif", "Species interaction factor", opts))
}
wd_det_labels <- function(sp, opts = NULL) {
  A <- sp[["dominant"]]; B <- sp[["subordinate"]]
  c(pA = fig_label("pA", sprintf("%s, %s absent (pA)", A, B), opts),
    pB = fig_label("pB", sprintf("%s, %s absent (pB)", B, A), opts),
    rA = fig_label("rA", sprintf("%s, both present (rA)", A), opts),
    rBA = fig_label("rBA", sprintf("%s, both present, %s detected (rBA)", B, A), opts),
    rBa = fig_label("rBa", sprintf("%s, both present, %s missed (rBa)", B, A), opts))
}

# Fixed colour and linetype per quantity, so colour follows the entity across
# figures: dominant species = Okabe-Ito 1 (orange), subordinate overall = 2
# (sky blue), subordinate where the dominant is absent / was missed = 3
# (bluish green), where it is present / was detected = 6 (vermillion);
# rA = 7 (reddish purple). Linetype is a second, greyscale-safe encoding.
wd_col_index <- c(psiA = 1, psiB = 2, psiBa = 3, psiBA = 6, sif = 5, pA = 1, rA = 7, pB = 2, rBa = 3, rBA = 6)
wd_lty <- c(psiA = "solid", psiB = "dotdash", psiBa = "22", psiBA = "42", sif = "solid",
            pA = "solid", rA = "13", pB = "dotdash", rBa = "22", rBA = "42")
wd_cols <- function(pars, opts) {
  pal <- if (is.null(opts$palette)) "okabe_ito" else opts$palette
  p <- if (length(pal) > 1) unlist(pal) else eco_pal[[pal]]
  if (is.null(p)) p <- eco_pal$okabe_ito
  idx <- wd_col_index[pars]
  # short palette (greyscale or a user list): use its darkest shades in turn;
  # linetype / shape still separate the quantities
  if (max(idx) > length(p)) idx <- ((match(pars, unique(pars)) - 1) %% min(3, length(p))) + 1
  stats::setNames(p[idx], pars)
}
wd_legend_guide <- function(position = NULL) {   # spacing keeps two-line labels apart
  a <- list(theme = ggplot2::theme(legend.key.width = ggplot2::unit(1.8, "lines"), legend.key.spacing.y = ggplot2::unit(3, "pt"),
                                   legend.key.spacing.x = ggplot2::unit(8, "pt")))
  if (!is.null(position)) { a$position <- position; a$ncol <- 1 }   # guide-level position wins over the core theme (ggplot2 >= 3.5); one column never clips long names
  do.call(ggplot2::guide_legend, a)
}
# Move a curve plot's legend under the panel (stacked / narrow layouts keep the full width for the panels).
wd_bottom_legend <- function(p, aes = c("colour", "fill", "linetype")) {
  p + do.call(ggplot2::guides, stats::setNames(rep(list(wd_legend_guide("bottom")), length(aes)), aes))
}
wd_scales <- function(d, opts, linetype = TRUE) {
  u <- unique(d[, c("par", "group")])
  vals <- stats::setNames(wd_cols(u$par, opts), u$group)
  out <- list(ggplot2::scale_colour_manual(values = vals, breaks = u$group),
              ggplot2::scale_fill_manual(values = vals, breaks = u$group))
  if (linetype) out <- c(out, list(
    ggplot2::scale_linetype_manual(values = stats::setNames(unname(wd_lty[u$par]), u$group), breaks = u$group),
    ggplot2::guides(colour = wd_legend_guide(), fill = wd_legend_guide(), linetype = wd_legend_guide()),
    ggplot2::labs(linetype = NULL)))
  out
}

# Covariates entered as log(x) are drawn on a log axis labelled in original units.
wd_is_log <- function(des, v) {
  pat <- paste0("log(10|2|1p)?\\(\\Q", v, "\\E\\)")
  any(vapply(des$mms, function(m) any(grepl(pat, m$cols, perl = TRUE)), TRUE))
}
wd_xscale <- function(des, v, x) {
  if (is.numeric(x) && wd_is_log(des, v) && isTRUE(min(x, na.rm = TRUE) > 0))
    ggplot2::scale_x_log10(breaks = scales::breaks_log(6), labels = scales::label_comma())
  else NULL
}
wd_logx_note <- function(des, v, x) if (!is.null(wd_xscale(des, v, x))) " (log-scaled axis, as the covariate enters the model as a log)" else ""

wd_eff_plot <- function(d, v, ylab, opts, ylim01 = TRUE, des = NULL) {
  if (is.numeric(d$level)) {
    fig_cond_effect(d[, c("x", "estimate", "lower", "upper", intersect(c("lower50", "upper50"), names(d)), "group")],
                    fig_label(v, v, opts), ylab, ylim01 = ylim01) +
      ggplot2::aes(linetype = group) + wd_scales(d, opts) +
      (if (!is.null(des)) wd_xscale(des, v, d$x) else NULL)
  } else {
    ggplot2::ggplot(d, ggplot2::aes(level, estimate, colour = group, shape = group)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(0.5), size = 0.35) +
      ggplot2::labs(x = fig_label(v, v, opts), y = ylab, colour = NULL, shape = NULL) +
      (if (ylim01) ggplot2::scale_y_continuous(limits = c(0, 1)) else NULL) + wd_scales(d, opts, linetype = FALSE)
  }
}

# One row per coefficient: descriptive labels (renamable with fig_label keyed
# by the estimates-table name), the species it belongs to, and the SD used to
# standardise its column (original units).
wd_coef_table <- function(des, sp, opts) {
  A <- sp[["dominant"]]; B <- sp[["subordinate"]]
  plab <- c(psiA = sprintf("%s occupancy", A), psiBA = sprintf("%s occupancy, %s present", B, A),
            psiBa = sprintf("%s occupancy, %s absent", B, A), pA = sprintf("%s detection, %s absent", A, B),
            pB = sprintf("%s detection, %s absent", B, A), rA = sprintf("%s detection, both present", A),
            rBA = sprintf("%s detection, %s detected", B, A), rBa = sprintf("%s detection, %s missed", B, A))
  rows <- list()
  for (p in des$free) {
    m <- des$mms[[p]]; ii <- des$idx[[p]]
    for (k in seq_along(ii)) {
      col <- m$cols[k]
      term <- if (k == 1) "intercept" else {
        lab <- col; vv <- setdiff(m$vars, ".dummy")
        vv <- vv[order(-nchar(vv))]
        for (v in vv) {
          if (!is.null(m$xlev[[v]]) && startsWith(col, v)) {
            lab <- sprintf("%s: %s vs %s", fig_label(v, v, opts), substring(col, nchar(v) + 1), m$xlev[[v]][1]); break
          }
          lab <- if (identical(col, v)) fig_label(v, v, opts) else gsub(v, fig_label(v, v, opts), lab, fixed = TRUE)
        }
        lab
      }
      nm <- des$coef_names[ii[k]]
      rows[[nm]] <- data.frame(name = nm, key = p, col = col, term = term,
                               long = fig_label(nm, sprintf("%s: %s", plab[[p]], term), opts),
                               short = fig_label(nm, sprintf("%s: %s", p, term), opts),
                               species = if (p %in% c("psiA", "pA", "rA")) A else B,
                               scale = m$scale[k], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

wd_wrap <- function(x, w = 24) vapply(x, function(s) paste(strwrap(s, w), collapse = "\n"), "", USE.NAMES = FALSE)

# Pre-theme panels of a composed figure (the core styles only the last panel).
wd_theme_panels <- function(pl, opts) {
  lapply(pl, function(p) p + theme_eco(opts$base_size, opts$font, opts$style, "bottom"))
}
wd_size <- function(p, width = "onehalf", aspect = 0.62) { attr(p, "fig_size") <- list(width = width, aspect = aspect); p }

# Text naming the value every other covariate in `pars`' formulas is held at
# (template tmpl): mean, geometric mean for log-entered terms, reference level.
wd_held <- function(des, v, pars, tmpl, L) {
  vv <- setdiff(wd_vars(des, pars), c(v, ".survey"))
  if (!length(vv)) return("no other covariates in these formulas")
  paste0(if (nzchar(v)) "other covariates" else "covariates", " held at ", paste(vapply(vv, function(u) {
    x <- tmpl[[u]]
    if (is.numeric(x)) sprintf("%s = %s (%s)", L(u, u), trimws(formatC(signif(x, 3), digits = 3, format = "fg", big.mark = ",")),
                               if (wd_is_log(des, u)) "geometric mean" else "mean")
    else sprintf("%s = %s (reference level)", L(u, u), x)
  }, ""), collapse = "; "))
}

# Panel whose axis title should sit by its own tick labels, not aligned with
# a neighbour's wide labels (avoids tag / title collisions in stacked layouts).
wd_free <- function(p) if (utils::packageVersion("patchwork") >= "1.2.0") patchwork::free(p, type = "label", side = "l") else p

wd_add_rug <- function(p, x) {
  if (!is.numeric(x)) return(p)
  p + ggplot2::geom_rug(data = data.frame(x = x), ggplot2::aes(x = x), inherit.aes = FALSE, alpha = 0.3,
                        sides = "b", length = ggplot2::unit(0.02, "npc"))
}

# ------------------------------------------------------------------ figures
# Publication figure suite (docs/FIGURES.md). Covariates stay in raw units in
# the data (standardisation happens only inside the model matrices, wd_mm),
# so curves are computed on a raw-unit grid; data$scaling is honoured if a
# caller supplies pre-standardised covariates (fig_x). Every label goes
# through fig_label(); theme, palette and sizes are applied by the core.
wd_figures <- function(fit, data, opts = fig_opts()) {
  o <- fit$object; des <- o$des
  if (o$engine == "naive") return(list())
  eco_require(c("ggplot2", "patchwork", "scales"))
  sp <- wd_species(fit)
  sp <- c(dominant = fig_label("dominant", sp[["dominant"]], opts), subordinate = fig_label("subordinate", sp[["subordinate"]], opts))
  fit$object$des$species <- sp; des$species <- sp
  scaling <- data$scaling
  fx <- function(x, v) if (is.numeric(x)) fig_x(x, v, opts, scaling) else x
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci
  iw <- if (o$engine == "nimble") "credible" else "confidence"
  inner <- if (o$engine == "nimble") opts$inner_ci else NULL
  inner_txt <- if (!is.null(inner)) sprintf("; inner band / thick line = %d%%", round(100 * inner)) else ""
  pt_word <- if (o$engine == "nimble") (if (identical(o$settings$point, "mean")) "posterior mean" else "posterior median") else "maximum-likelihood estimate"
  L <- function(key, default) fig_label(key, default, opts)
  ylab_occ <- L("occupancy_axis", "Occupancy probability")
  ylab_sub <- L("sub_occupancy_axis", sprintf("Occupancy of %s", sp[["subordinate"]]))
  ylab_det <- L("detection_axis", "Detection probability per survey")
  sif_lab <- L("sif_axis", "Species interaction factor")
  ol <- wd_occ_labels(sp, opts); dl <- wd_det_labels(sp, opts)
  figs <- list()
  occ_int <- des$map[["psiBA"]] != des$map[["psiBa"]]
  # prepare a covariate curve: x in original units, rug of observed values
  curve <- function(v, occ, keep) {
    d <- wd_curve(fit, v, occ, opts); d <- d[d$par %in% keep, ]
    if (is.numeric(d$level)) { d$x <- fx(d$level, v); d$level <- d$x }
    d
  }
  obs_x <- function(v, occ) { x <- wd_cov_df(des, v, occ)[[v]]; if (is.numeric(x)) fx(x, v) else NULL }

  # ---- conditional_occupancy: (A) estimates at average conditions, (B) along a covariate
  nd0 <- wd_template(des$site_covs, des$mms)
  q0 <- wd_occ_draws(fit, nd0)
  p0 <- if (occ_int) c("psiA", "psiB", "psiBa", "psiBA") else c("psiA", "psiB")
  pd <- do.call(rbind, lapply(p0, function(p) data.frame(par = p, param = ol[[p]], wd_summ(q0[[p]], ci = ci, inner = inner))))
  pd$param <- factor(wd_wrap(pd$param), levels = rev(wd_wrap(ol[p0])))
  pd$par <- factor(pd$par, levels = p0)
  cols0 <- wd_cols(p0, opts)
  pA <- ggplot2::ggplot(pd, ggplot2::aes(estimate, param, colour = par, shape = par)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.5)
  if (!is.null(inner)) pA <- pA + ggplot2::geom_linerange(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.3)
  pA <- pA + ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_colour_manual(values = cols0, guide = "none") +
    ggplot2::scale_shape_manual(values = c(psiA = 16, psiB = 17, psiBa = 15, psiBA = 18)[p0], guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = ylab_occ, y = NULL)
  occ_pars <- c("psiA", "psiBA", "psiBa")
  held0 <- wd_held(des, "", occ_pars, nd0, L)
  occ_vars <- if (occ_int) wd_vars(des, c("psiBA", "psiBa")) else character()
  if (length(occ_vars)) {
    v <- occ_vars[1]
    d <- curve(v, TRUE, c("psiBA", "psiBa"))
    pB <- wd_add_rug(wd_eff_plot(d, v, ylab_sub, opts, des = des), obs_x(v, TRUE))
    tp <- wd_theme_panels(list(pA, wd_bottom_legend(pB)), opts)
    figs$conditional_occupancy <- with_caption(compose_figures(list(tp[[1]], wd_free(tp[[2]])), ncol = 1, guides = "keep"), sprintf(
      "(A) Occupancy at average covariate values (%s): %s, %s overall, and %s conditional on %s being absent or present (point = %s; lines = %d%% %s interval%s). (B) Occupancy of %s where %s is present (long dashes) vs absent (short dashes) along %s%s, %s; bands = %d%% %s intervals%s; rug = observed sites.",
      held0, sp[["dominant"]], sp[["subordinate"]], sp[["subordinate"]], sp[["dominant"]],
      pt_word,
      round(100 * ci), iw, inner_txt, sp[["subordinate"]], sp[["dominant"]], L(v, v), wd_logx_note(des, v, d$x),
      wd_held(des, v, c("psiBA", "psiBa"), nd0, L), round(100 * ci), iw, inner_txt))
    attr(figs$conditional_occupancy, "fig_size") <- list(width = "onehalf", aspect = 1.05)
    for (v2 in occ_vars[-1]) {
      d <- curve(v2, TRUE, c("psiBA", "psiBa"))
      figs[[paste0("conditional_occupancy_", wd_clean(v2))]] <- with_caption(
        wd_size(wd_add_rug(wd_eff_plot(d, v2, ylab_sub, opts, des = des), obs_x(v2, TRUE))),
        sprintf("Occupancy of %s where %s is present vs absent along %s%s (%s); bands = %d%% %s intervals%s; rug = observed sites.",
                sp[["subordinate"]], sp[["dominant"]], L(v2, v2), wd_logx_note(des, v2, d$x), wd_held(des, v2, c("psiBA", "psiBa"), nd0, L), round(100 * ci), iw, inner_txt))
    }
  } else {
    figs$conditional_occupancy <- with_caption(wd_size(pA, "single", 0.6), sprintf(
      "Occupancy at average covariate values (%s): %s, %s overall%s (point = %s; lines = %d%% %s interval%s).%s",
      held0, sp[["dominant"]], sp[["subordinate"]], if (occ_int) sprintf(", and %s conditional on %s present/absent", sp[["subordinate"]], sp[["dominant"]]) else "",
      pt_word, round(100 * ci), iw, inner_txt,
      if (occ_int) "" else " Occupancy of the subordinate species was constrained equal with the dominant present and absent (no psiBA formula), so no conditional contrast is estimated."))
  }

  # ---- sif: species interaction factor, null line at 1
  if (occ_int) {
    d0 <- data.frame(param = L("avg_conditions", "Average\nconditions"), wd_summ(q0$sif, ci = ci, inner = inner))
    ps <- ggplot2::ggplot(d0, ggplot2::aes(estimate, param)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.5, colour = eco_pal$ink)
    if (!is.null(inner)) ps <- ps + ggplot2::geom_linerange(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.3, colour = eco_pal$ink)
    ps <- ps + ggplot2::geom_point(size = 2, colour = eco_pal$ink) +
      ggplot2::labs(x = L("sif_axis_null", paste0(sif_lab, "\n(1 = independent occurrence)")), y = NULL)
    sv <- wd_vars(des, c("psiA", "psiBA", "psiBa"))
    sv <- sv[vapply(sv, function(v) is.numeric(des$site_covs[[v]]), TRUE)]
    if (length(sv)) {
      v <- sv[1]
      d <- curve(v, TRUE, "sif"); d$group <- NULL
      pc <- fig_cond_effect(d[, c("x", "estimate", "lower", "upper", intersect(c("lower50", "upper50"), names(d)))],
                            L(v, v), L("sif_short_axis", "SIF")) +
        ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null) + wd_xscale(des, v, d$x)
      pc <- wd_add_rug(pc, obs_x(v, TRUE))
      tp <- wd_theme_panels(list(ps, pc), opts)
      figs$sif <- with_caption(compose_figures(list(tp[[1]], wd_free(tp[[2]])), ncol = 1, guides = "keep") +
                                 patchwork::plot_layout(heights = c(1, 2.2)), sprintf(   # panel A is a single row
        "Species interaction factor SIF = P(both)/(P(%s) x P(%s)) (MacKenzie et al. 2004); dashed line = 1 (the species occur independently); < 1 = they co-occur less often than expected by chance. (A) At average covariate values (%s), point = %s with %d%% %s interval%s. (B) Along %s%s, %s; line = %s; band = %d%% %s interval%s; rug = observed sites. SIF depends on occupancy of %s, so report it with the log-odds ratio of psiBA vs psiBa.",
        sp[["dominant"]], sp[["subordinate"]], held0, pt_word, round(100 * ci), iw, inner_txt, L(v, v), wd_logx_note(des, v, d$x),
        wd_held(des, v, occ_pars, nd0, L), pt_word, round(100 * ci), iw, inner_txt, sp[["dominant"]]))
      attr(figs$sif, "fig_size") <- list(width = "single", aspect = 1.25)
    } else figs$sif <- with_caption(wd_size(ps, "single", 0.45), sprintf(
      "Species interaction factor at average covariate values (%s; point = %s, %d%% %s interval%s); dashed line = 1 (independence).", held0, pt_word, round(100 * ci), iw, inner_txt))
  }

  # ---- detection_interaction: pA, rA, pB, rBa, rBA
  tmpl_s <- wd_template(des$surv_df[des$obs_v, , drop = FALSE], des$mms)
  pars <- c("pA", "rA", "pB", "rBa", "rBA"); pars <- pars[!duplicated(des$map[pars])]
  qd <- wd_prob_draws(fit, tmpl_s, pars)
  spc <- stats::setNames(unname(wd_cols(c("psiA", "psiB"), opts)), sp)
  dd <- do.call(rbind, lapply(pars, function(p) data.frame(param = dl[[p]], wd_summ(qd[[p]], ci = ci, inner = inner),
                                                          species = if (p %in% c("pA", "rA")) sp[["dominant"]] else sp[["subordinate"]])))
  dd$param <- factor(wd_wrap(dd$param, 30), levels = rev(wd_wrap(dl[pars], 30)))
  dd$species <- factor(dd$species, levels = unname(sp))
  pdd <- ggplot2::ggplot(dd, ggplot2::aes(estimate, param, colour = species, shape = species)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.5)
  if (!is.null(inner)) pdd <- pdd + ggplot2::geom_linerange(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.3)
  pdd <- pdd + ggplot2::geom_point(size = 1.8) + ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = spc) + ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17), sp)) +
    ggplot2::labs(x = ylab_det, y = NULL, colour = NULL, shape = NULL)
  det_int <- des$map[["rBA"]] != des$map[["rBa"]]
  figs$detection_interaction <- with_caption(pdd, sprintf(
    "Per-survey detection probabilities at average survey conditions (%s; point = %s; lines = %d%% %s intervals%s).%s pA vs rA and pB vs rBa compare detection where the other species is absent vs present. Parameters without their own formula are constrained equal and shown once.",
    wd_held(des, "", pars, tmpl_s, L), pt_word, round(100 * ci), iw, inner_txt,
    if (det_int) sprintf(" rBA vs rBa compares detection of %s, where both species occur, on surveys that did vs did not detect %s.", sp[["subordinate"]], sp[["dominant"]]) else
      " rBA was constrained equal to rBa (no detection interaction estimated)."))
  attr(figs$detection_interaction, "fig_size") <- list(width = "onehalf", height = 30 + 11 * length(pars))

  # ---- coefficients: odds ratios per 1 SD (numeric) or vs reference level (factor)
  ct <- wd_coef_table(des, sp, opts)
  est <- fit$estimates
  sl <- est[est$param %in% ct$name[ct$term != "intercept"], ]
  if (nrow(sl)) {
    ct2 <- ct[match(sl$param, ct$name), ]
    if (!is.null(inner) && all(sl$param %in% colnames(o$draws))) {   # inner (50%) interval from the draws
      b <- (1 - inner) / 2
      sl$lower50 <- apply(o$draws[, sl$param, drop = FALSE], 2, stats::quantile, b)
      sl$upper50 <- apply(o$draws[, sl$param, drop = FALSE], 2, stats::quantile, 1 - b)
    }
    # width the figure will be exported at: narrow (single-column) layouts get
    # the short psiA/pB-style labels and a legend under the panel
    cw <- opts$per_figure$coefficients$width; if (is.null(cw)) cw <- "onehalf"
    narrow <- isTRUE(fig_width_mm(cw, opts) < 110)
    lab_col <- if (narrow) "short" else "long"; ww <- if (narrow) 26 else 40
    sl$param <- wd_wrap(ct2[[lab_col]], ww)
    sl$group <- factor(ct2$species, levels = unname(sp))
    lev <- rev(unique(wd_wrap(ct[[lab_col]][ct$name %in% ct2$name], ww)))   # model order: psiA, psiBA, psiBa, pA, ...
    pcf <- fig_coef(sl[, intersect(c("param", "estimate", "lower", "upper", "lower50", "upper50", "group"), names(sl))], null_value = 1, transform = exp, log_x = TRUE,
                    order = FALSE, xlab = L("odds_ratio_axis", sprintf("Odds ratio per 1 SD (%d%% %s interval)", round(100 * ci), iw)))
    pcf$data$param <- factor(as.character(pcf$data$param), levels = lev[lev %in% pcf$data$param])
    pcf <- pcf + ggplot2::aes(shape = group) + ggplot2::scale_colour_manual(values = spc) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17), sp)) + ggplot2::labs(shape = NULL)
    if (narrow) pcf <- pcf + ggplot2::guides(colour = ggplot2::guide_legend(position = "bottom", direction = "vertical"),
                                             shape = ggplot2::guide_legend(position = "bottom", direction = "vertical"))
    num <- ct2[ct2$scale != 1, ]
    num <- num[!duplicated(num$term), ]
    sd_txt <- if (nrow(num)) paste0(" 1 SD: ", paste(sprintf("%s = %s", num$term, formatC(num$scale, digits = 3, format = "g")), collapse = "; "), ".") else ""
    figs$coefficients <- with_caption(pcf, sprintf(
      "Covariate effects as odds ratios per 1 SD of each numeric covariate (or vs the reference level for factors): point = %s, thin line = %d%% %s interval%s; dashed line = no effect; log-scaled axis. Colour and shape = species whose probability the coefficient acts on.%s Parameter names follow Richmond et al. (2010): psiA (%s occupancy), psiBA / psiBa (%s occupancy with %s present / absent), pA, pB, rA, rBA, rBa (detection).",
      pt_word, round(100 * ci), iw, sub("inner band / thick line", "thick line", inner_txt), sd_txt, sp[["dominant"]], sp[["subordinate"]], sp[["dominant"]]))
    attr(figs$coefficients, "fig_size") <- list(width = "onehalf", height = 24 + 9 * nrow(sl) + if (narrow) 14 else 0)
  }

  # ---- occupancy_<cov> and detection_<cov>: conditional effects
  for (v in wd_vars(des, c("psiA", "psiBA", "psiBa"))) {
    d <- curve(v, TRUE, c(if (v %in% des$mms[[des$map[["psiA"]]]]$vars) "psiA", "psiB", if (occ_int) c("psiBA", "psiBa")))
    figs[[paste0("occupancy_", wd_clean(v))]] <- with_caption(
      wd_size(wd_add_rug(wd_eff_plot(d, v, ylab_occ, opts, des = des), obs_x(v, TRUE))),
      sprintf("Occupancy along %s%s (%s): %s%s overall (averaged over %s presence)%s. Bands = %d%% %s intervals%s; rug = observed sites. Curves that are flat do not depend on this covariate in the fitted model.",
              L(v, v), wd_logx_note(des, v, d$x), wd_held(des, v, occ_pars, nd0, L), if ("psiA" %in% d$par) paste0(sp[["dominant"]], ", ") else "", sp[["subordinate"]], sp[["dominant"]],
              if (occ_int) sprintf(", and %s with %s present / absent", sp[["subordinate"]], sp[["dominant"]]) else "",
              round(100 * ci), iw, inner_txt))
  }
  for (v in setdiff(wd_vars(des, c("pA", "pB", "rA", "rBA", "rBa")), ".survey")) {
    d <- wd_curve(fit, v, FALSE, opts)
    up <- unique(d$par); keepp <- up[vapply(up, function(p) v %in% des$mms[[des$map[[p]]]]$vars, TRUE)]
    d <- curve(v, FALSE, keepp)
    figs[[paste0("detection_", wd_clean(v))]] <- with_caption(
      wd_size(wd_add_rug(wd_eff_plot(d, v, ylab_det, opts, des = des), obs_x(v, FALSE))),
      sprintf("Per-survey detection probability along %s%s (%s); bands = %d%% %s intervals%s; rug = surveys done. Only detection parameters whose formula includes this covariate are drawn.",
              L(v, v), wd_logx_note(des, v, d$x), wd_held(des, v, keepp, tmpl_s, L), round(100 * ci), iw, inner_txt))
  }

  # ---- gof: predictive check of co-detection classes
  if (!is.null(o$gof)) {
    g <- o$gof
    figs$gof <- with_caption(fig_ppc(g$rep_pattern - g$obs_pattern, 0, L("gof_axis", "Replicated minus observed chi-square")) +
                               ggplot2::labs(y = L("replicates_axis", "Replicates")),
      sprintf("%s of the number of sites in each co-detection class (neither, %s only, %s only, both detected). Values right of the line mean replicated data fit worse than the observed data; p = %.2f (extreme values near 0 or 1 indicate lack of fit).",
              tools::toTitleCase(g$method), sp[["dominant"]], sp[["subordinate"]], g$p_pattern))
  }

  # ---- trace and posterior (Bayesian engine)
  if (o$engine == "nimble" && !is.null(o$chains)) {
    sn <- stats::setNames(wd_wrap(ct$short, 18), ct$name)
    ch <- lapply(o$chains, function(m) { m <- as.matrix(m)[, des$coef_names, drop = FALSE]; colnames(m) <- sn[colnames(m)]; m })
    figs$trace <- with_caption(fig_trace(ch) + ggplot2::scale_x_continuous(breaks = scales::breaks_extended(3)) +
                                 ggplot2::labs(x = L("iteration_axis", "Iteration (after burn-in)"), y = L("logit_axis", "Coefficient (logit scale)")),
      "MCMC traces of all coefficients on the logit scale (slopes per 1 SD of the covariate), one colour per chain; well-mixed chains overlap with no trend.")
    attr(figs$trace, "fig_size") <- list(width = "double", aspect = 0.8)
    D <- do.call(rbind, lapply(o$chains, function(m) wd_report_cols(des, m)))
    pp <- intersect(c("occ_sub_given_dom_present", "occ_sub_given_dom_absent", "species_interaction_factor",
                      "occ_log_odds_ratio_sub", "det_sub_with_dom_detected", "det_sub_with_dom_missed", "det_log_odds_ratio_sub"), colnames(D))
    plab <- c(occ_sub_given_dom_present = ol[["psiBA"]], occ_sub_given_dom_absent = ol[["psiBa"]], species_interaction_factor = sif_lab,
              occ_log_odds_ratio_sub = sprintf("Occupancy log-odds ratio, %s present vs absent", sp[["dominant"]]),
              det_sub_with_dom_detected = dl[["rBA"]], det_sub_with_dom_missed = dl[["rBa"]],
              det_log_odds_ratio_sub = "Detection log-odds ratio, rBA vs rBa")
    plab <- vapply(names(plab), function(k) wd_wrap(L(k, plab[[k]]), 34), "")
    Dp <- D[, pp, drop = FALSE]; colnames(Dp) <- plab[pp]
    tr <- data$truth
    tr <- if (!is.null(tr)) { tr <- tr[intersect(names(tr), pp)]; stats::setNames(unname(tr), plab[names(tr)]) } else NULL
    nulls <- data.frame(param = plab[intersect(pp, c("species_interaction_factor", "occ_log_odds_ratio_sub", "det_log_odds_ratio_sub"))])
    nulls$t <- ifelse(grepl("interaction", names(plab)[match(nulls$param, plab)]), 1, 0)
    pp_fig <- fig_posterior(Dp, truth = if (length(tr)) tr else NULL, ncol = 2)
    pp_fig$data$param <- factor(as.character(pp_fig$data$param), levels = unname(plab[pp]))   # logical, not alphabetical, order
    pp_fig <- pp_fig + ggplot2::scale_y_continuous(breaks = NULL)   # density height carries no units; keep panels legible
    if (nrow(nulls)) pp_fig <- pp_fig + ggplot2::geom_vline(data = nulls, ggplot2::aes(xintercept = t), linetype = 2, colour = eco_pal$null, linewidth = 0.4)
    pp_fig <- pp_fig + ggplot2::labs(x = L("value_axis", "Value (probability, ratio or log-odds ratio)"))
    for (k in seq_along(pp_fig$layers)) {          # same facet order in every layer (truth / null lines)
      ld <- pp_fig$layers[[k]]$data
      if (is.data.frame(ld) && "param" %in% names(ld)) pp_fig$layers[[k]]$data$param <- factor(as.character(ld$param), levels = unname(plab[pp]))
    }
    figs$posterior <- with_caption(pp_fig, sprintf(
      "Posterior densities of the interaction quantities at average covariate values (all chains pooled). Dashed lines = no interaction (SIF = 1, log-odds ratio = 0).%s",
      if (length(tr)) " Solid red lines = true values used to simulate the data." else ""))
    attr(figs$posterior, "fig_size") <- list(width = "onehalf", aspect = 0.3 * ceiling(length(pp) / 2) + 0.1)
  }
  figs
}

# Legacy plots for the HTML report (same figures, default styling).
wd_plots <- function(fit, data) {
  f <- tryCatch(wd_figures(fit, data, fig_opts()), error = function(e) { message("figures failed: ", conditionMessage(e)); list() })
  f$trace <- NULL
  nm <- c(conditional_occupancy = "Conditional occupancy of the subordinate species",
          sif = "Species interaction factor", detection_interaction = "Detection probabilities and the detection interaction",
          coefficients = "Covariate effects (odds ratios)", gof = "Goodness of fit")
  names(f) <- ifelse(names(f) %in% names(nm), nm[names(f)], gsub("_", " ", names(f)))
  lapply(f, function(p) if (inherits(p, "ggplot")) p + theme_eco(11) else p)
}

# ------------------------------------------------------------ interpretation
wd_cap <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
wd_plain_sp <- function(key, A, B) {
  x <- wd_plain[[key]]
  x <- gsub("the dominant species|the dominant", A, x)
  gsub("the subordinate species|the subordinate", B, x)
}

wd_interpret <- function(fit, data) {
  o <- fit$object; des <- o$des; e <- fit$estimates
  sp <- wd_species(fit); A <- sp[["dominant"]]; B <- sp[["subordinate"]]
  g <- function(n) wd_row(fit, n)
  ci <- function(n, d = 2) { r <- g(n); if (!nrow(r)) return("NA"); wd_ci2(fit, r$estimate, r$lower, r$upper, d) }
  it <- wd_int(fit)
  out <- character()
  if (o$engine == "naive") return("Naive model ignoring imperfect detection: a negative control only. Its estimates are biased whenever detection is imperfect; do not interpret them.")
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge (or estimates lie on a boundary). The numbers below are not reliable; see Diagnostics.")
  occ_int <- des$map[["psiBA"]] != des$map[["psiBa"]]
  out <- c(out, sprintf("At average covariate values, %s (treated as the dominant species) occupies an estimated %s of sites. %s occupies %s of sites overall.",
                        A, ci("occ_dominant"), wd_cap(B), ci("occ_subordinate")))
  if (occ_int) {
    lor <- g("occ_log_odds_ratio_sub")
    dir <- if (lor$upper < 0) sprintf("lower where %s is present", A) else if (lor$lower > 0) sprintf("higher where %s is present", A)
           else sprintf("not clearly different where %s is present (the interval includes no difference)", A)
    out <- c(out, sprintf(paste(
      "Where %s is present, %s occupies %s of sites; where %s is absent, %s. The odds of %s occupying a site are %s times as high where %s is present",
      "(odds ratio = exp of the log-odds ratio %s), so %s occupancy is %s."),
      A, B, ci("occ_sub_given_dom_present"), A, ci("occ_sub_given_dom_absent"), B,
      wd_ci2(fit, exp(lor$estimate), exp(lor$lower), exp(lor$upper)), A, ci("occ_log_odds_ratio_sub"), B, dir))
    s <- g("species_interaction_factor")
    sdir <- if (s$upper < 1) "less often than expected if they were distributed independently (a segregation pattern)"
            else if (s$lower > 1) "more often than expected if they were distributed independently (an aggregation pattern)"
            else "about as often as expected under independence; the data cannot distinguish segregation from aggregation"
    out <- c(out, sprintf("The species interaction factor (SIF = P(both) / [P(%s) x P(%s)]; MacKenzie et al. 2004) is %s: at average conditions the two species co-occur %s. Probability both occupy a site: %s.",
                          A, B, ci("species_interaction_factor"), sdir, ci("occ_both")))
  } else out <- c(out, sprintf("This model assumes %s occupancy does not depend on %s (psiBA = psiBa), so it cannot detect an occupancy interaction; add a psiBA formula to estimate one.", B, A))
  # covariate effects on occupancy
  sl <- e[grepl("_beta_", e$param), ]
  for (k in seq_len(nrow(sl))) {
    n <- sl$param[k]; proc <- sub("_beta_.*$", "", n); key <- names(wd_names)[match(proc, wd_names)]
    term <- sub("^.*_beta_", "", n)
    f <- des$map[[key]]; mm <- des$mms[[f]]; col <- mm$cols[wd_clean(mm$cols) == term][1]
    sdu <- mm$scale[match(col, mm$cols)]
    isbin <- !is.na(sdu) && sdu == 1 && mm$center[match(col, mm$cols)] == 0
    unit <- if (isbin) sprintf("for %s (vs the reference level)", col) else sprintf("per 1 SD of %s (1 SD = %s)", col, signif(sdu, 3))
    or <- wd_ci2(fit, exp(sl$estimate[k]), exp(sl$lower[k]), exp(sl$upper[k]))
    dirw <- if (sl$lower[k] > 0) "increases" else if (sl$upper[k] < 0) "decreases" else "shows no clear change"
    out <- c(out, sprintf("%s (%s): the odds are multiplied by %s %s; i.e. it %s with %s%s. %s",
                          wd_cap(wd_plain_sp(key, A, B)), key, or, unit, dirw, col,
                          if (!isbin) wd_range_text(fit, key, col) else "",
                          if (key %in% wd_occ_pars) "" else "(This is a detection effect, not a habitat effect.)"))
  }
  # detection
  det_int <- des$map[["rBA"]] != des$map[["rBa"]]
  out <- c(out, sprintf("Detection per survey at average conditions: %s where %s is absent %s; %s where %s is absent %s.",
                        A, B, ci("det_dominant_alone"), B, A, ci("det_sub_alone")))
  if (des$map[["rA"]] == "rA") out <- c(out, sprintf("Where both occur, %s is detected with probability %s.", A, ci("det_dominant_with_sub")))
  if (det_int) {
    d <- g("det_log_odds_ratio_sub")
    ddir <- if (d$upper < 0) sprintf("less likely to be detected on surveys that also detected %s", A)
            else if (d$lower > 0) sprintf("more likely to be detected on surveys that also detected %s", A)
            else sprintf("about equally likely to be detected whether or not %s was detected on the same survey", A)
    out <- c(out, sprintf("Where both occur, %s is detected on %s of surveys that detected %s (rBA) and on %s of surveys that missed %s (rBa); log-odds ratio %s. %s is %s. A detection interaction can reflect behaviour (e.g. the subordinate calling less or avoiding the station when the dominant is active) or shared survey conditions; it is a detection effect, not an occupancy effect.",
                          B, ci("det_sub_with_dom_detected"), A, ci("det_sub_with_dom_missed"), A, ci("det_log_odds_ratio_sub"),
                          wd_cap(B), ddir))
  }
  out <- c(out, sprintf(paste(
    "Caveats. (1) The model assumes %s is dominant: %s occupancy may depend on %s but not the reverse. That direction is an assumption chosen before fitting, not a result; with roles swapped the fit describes the same joint distribution differently.",
    "(2) A difference between psiBA and psiBa, or an SIF away from 1, is a co-occurrence pattern. It is consistent with competition or predation but equally with both species responding to habitat features not in the model (Blanchet et al. 2020); it does not by itself show that %s excludes %s.",
    "(3) Covariate associations are correlational."), A, B, A, A, B))
  if (!is.null(o$gof)) out <- c(out, sprintf("Model check: %s p-values were %.2f (co-detection pattern) and %.2f (per-site detection counts); values between about 0.05 and 0.95 mean the model reproduces these features of the data.",
                                             o$gof$method, o$gof$p_pattern, o$gof$p_counts))
  out <- c(out, sprintf("Validated limits (simulation recovery, see VALIDATION.md): with 200 sites x 4 surveys every estimate passed the bias check (mean error within about 10%% of the truth) and 95%% intervals covered it 91-98%% of the time. With 60 sites, occupancy of %s where %s is present and the probability that both occur were over-estimated by about 17-21%%, so the SIF was pulled toward 1 (+12%%): the interaction is understated, not exaggerated, and power to detect it was only about 0.4. This data set has %d sites%s.",
    sp[["subordinate"]], sp[["dominant"]], des$M, if (des$M < 100) " -- treat the size of the interaction as a conservative (too small) estimate" else ""))
  out
}

# "rises from 0.31 at 2 to 0.62 at 40" over the observed range of a numeric covariate.
wd_range_text <- function(fit, key, col) {
  des <- fit$object$des
  occ <- key %in% wd_occ_pars
  df <- wd_cov_df(des, col, occ)
  if (!col %in% names(df)) {               # transformed term, e.g. log(dist): vary the raw variable
    v <- tryCatch(all.vars(parse(text = col)[[1]]), error = function(e) character())
    if (length(v) != 1) return("")
    col <- v
  }
  if (!col %in% names(df) || !is.numeric(df[[col]])) return("")
  r <- range(df[[col]], na.rm = TRUE)
  nd <- wd_template(df, des$mms)[c(1, 1), , drop = FALSE]; nd[[col]] <- r
  q <- wd_prob_draws(fit, nd, key)[[key]]
  pt <- attr(q, "point")
  sprintf(" (from %s at %s = %s to %s at %s = %s, other covariates at their means)",
          wd_fmt(pt[1]), col, signif(r[1], 3), wd_fmt(pt[2]), col, signif(r[2], 3))
}

# ------------------------------------------------------------- methods text
wd_methods_text <- function(fit, data, settings) {
  o <- fit$object; des <- o$des; st <- o$settings
  sp <- wd_species(fit)
  fstr <- paste(sprintf("%s %s", names(des$forms), unlist(des$forms)), collapse = "; ")
  tied <- wd_pars[des$map[wd_pars] != wd_pars]
  eng <- switch(o$engine,
    nimble = sprintf(paste(
      "Models were fitted in a Bayesian framework with NIMBLE (de Valpine et al. 2017), using the marginal likelihood (latent occupancy states summed out) coded as a custom distribution.",
      "Intercepts had Logistic(0, 1) priors (uniform on the probability scale at average covariate values) and slopes Normal(0, SD = %s) priors on standardised covariates.",
      "We ran %d chains of %d iterations (%d burn-in) with slice samplers and assessed convergence with rank-normalised R-hat (< 1.01) and bulk effective sample size (>= 400)."),
      format(st$prior_sd), st$nchains, st$niter, st$nburnin),
    mle = "Models were fitted by maximum likelihood (hand-coded likelihood maximised with BFGS from several starting values); 95% confidence intervals are Wald intervals on the logit scale, and intervals for derived quantities come from 2000 draws from the asymptotic normal distribution of the estimates.",
    naive = "A naive model ignoring imperfect detection was fitted (negative control).")
  paste(
    sprintf("We modelled co-occurrence of %s (dominant species, A) and %s (subordinate species, B) with the conditional two-species occupancy model of Richmond et al. (2010) and Waddle et al. (2010).", sp[["dominant"]], sp[["subordinate"]]),
    "Occupancy of A (psiA) and of B conditional on A being present (psiBA) or absent (psiBa) were estimated together with per-survey detection probabilities of A where B is absent (pA) or present (rA), and of B where A is absent (pB) or, where both occur, on surveys that did (rBA) or did not (rBa) detect A.",
    sprintf("Each probability was modelled on the logit scale with covariates standardised to mean 0 and SD 1 (%s)%s.", fstr,
            if (length(tied)) sprintf("; %s were constrained equal to %s", paste(tied, collapse = ", "), paste(des$map[tied], collapse = ", ")) else ""),
    sprintf("Data comprised %d sites and up to %d surveys per site (%d surveys in total).", des$M, des$J, sum(des$obs)),
    "Derived quantities were the species interaction factor SIF = psiA psiBA / (psiA psiB) with psiB = psiA psiBA + (1 - psiA) psiBa (MacKenzie et al. 2004), and the log-odds ratio logit(psiBA) - logit(psiBa).",
    eng,
    if (!is.null(o$gof)) "Model fit was assessed with predictive checks of the number of sites in each co-detection class and of per-site detection counts." else "",
    "The estimation procedure was validated by simulation (ecotools occu-waddle VALIDATION.md).")
}
