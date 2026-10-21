# Output for movement-hmm: response-scale curves, figures, interpretation,
# methods text. Sourced by module.R.

state_label <- function(obj, k, opts = NULL) {
  K <- obj$K
  def <- if (!is.null(obj$state_names) && length(obj$state_names) >= k) obj$state_names[[k]] else
    sprintf("State %d (%s steps)", k, if (k == 1) "shortest" else if (k == K) "longest" else "medium")
  fig_label(paste0("state", k), def, opts)
}

# Draw matrix of beta (one per row of draws) -> list of K x ... computations.
beta_draws <- function(obj, n = 400) {
  dr <- obj$draws
  if (!dr$ok) return(list())
  W <- dr$W[seq_len(min(n, nrow(dr$W))), , drop = FALSE]
  lapply(seq_len(nrow(W)), function(i) w2nat(W[i, ], obj$K, obj$zm, obj$bdim)$beta)
}

beta_hat <- function(obj) w2nat(obj$draws$est, obj$K, obj$zm, obj$bdim)$beta

# Response-scale curves over hour of day or one covariate (original units).
# type = "stationary" (long-run share of time per state) or "transition"
# (probability of switching i -> j in one time step).
hmm_curve <- function(obj, var = "hour", type = c("stationary", "transition"), ngrid = 60, ci = 0.95, opts = NULL, inner = 0.5) {
  type <- match.arg(type)
  K <- obj$K; terms <- obj$terms
  dt <- obj$prep$interval_h; ph0 <- obj$prep$phase0 %||% 0
  per_hour <- var == "hour" && type == "stationary" && periodic_ok(dt)
  if (per_hour) {
    # state probabilities at each clock hour of the daily cycle (periodically
    # stationary distribution), drawn at the fix times and closed at 24 h
    hrs <- (ph0 + (seq_len(round(24 / dt)) - 1) * dt) %% 24
    o <- order(hrs); xs <- c(hrs[o], hrs[o[1]] + 24)
    rows <- NULL
  } else if (var == "hour") {
    xs <- seq(0, 24, length.out = ngrid); hours_avg <- NULL
    rows <- lapply(xs, function(h) list(design_row(terms, hour = h)))
  } else {
    raw <- obj$prep$grid[[var]]
    xs <- seq(stats::quantile(raw, 0.01, na.rm = TRUE), stats::quantile(raw, 0.99, na.rm = TRUE), length.out = ngrid)
    sc <- obj$prep$scaling[[var]]
    hh <- if ("tod_cos" %in% terms) 0:23 else NULL
    rows <- lapply(xs, function(x) {
      z <- stats::setNames(list((x - sc$center) / sc$scale), var)
      if (is.null(hh)) list(design_row(terms, z = z)) else lapply(hh, function(h) design_row(terms, hour = h, z = z))
    })
  }
  val <- function(B) {
    if (per_hour) { P <- periodic_probs(B, terms, K, dt, ph0)$P; return(rbind(P[o, , drop = FALSE], P[o[1], ])) }
    if (type == "stationary" && var != "hour") {
      sc <- obj$prep$scaling[[var]]
      return(t(vapply(xs, function(x) activity_budget(B, terms, K, dt, ph0,
                                                      stats::setNames(list((x - sc$center) / sc$scale), var)), numeric(K))))
    }
    do.call(rbind, lapply(rows, function(rr) {
      v <- lapply(rr, function(xr) {
        G <- tpm_from_eta(as.vector(xr %*% B), K)
        if (type == "stationary") stationary_dist(G) else unlist(lapply(seq_len(K), function(i) G[i, -i]))
      })
      colMeans(do.call(rbind, v))
    }))
  }
  labs <- if (type == "stationary") vapply(seq_len(K), function(k) state_label(obj, k, opts), "")
          else {
            # transition legend uses the same state names as the other figures
            # (independent verification 2026-10-20: it said "State 1 → 2" next
            # to figures labelled "Encamped"/"Travelling"); parenthetical
            # qualifiers dropped to keep the legend short
            sl <- vapply(seq_len(K), function(k) sub("\\s*\\(.*\\)$", "", state_label(obj, k, opts)), "")
            unlist(lapply(seq_len(K), function(i) sprintf("%s → %s", sl[i], sl[setdiff(seq_len(K), i)])))
          }
  hat <- val(beta_hat(obj))
  bd <- beta_draws(obj)
  a <- (1 - ci) / 2
  if (length(bd)) {
    arr <- simplify2array(lapply(bd, val))   # ngrid x nlab x ndraw
    lo <- apply(arr, c(1, 2), stats::quantile, a); hi <- apply(arr, c(1, 2), stats::quantile, 1 - a)
    a50 <- (1 - inner) / 2
    lo50 <- apply(arr, c(1, 2), stats::quantile, a50); hi50 <- apply(arr, c(1, 2), stats::quantile, 1 - a50)
  } else lo <- hi <- lo50 <- hi50 <- hat * NA
  out <- do.call(rbind, lapply(seq_along(labs), function(j)
    data.frame(x = xs, group = labs[j], estimate = hat[, j], lower = lo[, j], upper = hi[, j],
               lower50 = lo50[, j], upper50 = hi50[, j])))
  out$group <- factor(out$group, levels = labs)
  # origin / destination state index of each curve (figures colour by origin)
  if (type == "transition") {
    fr <- unlist(lapply(seq_len(K), function(i) rep(i, K - 1))); to <- unlist(lapply(seq_len(K), function(i) setdiff(seq_len(K), i)))
    out$from <- fr[as.integer(out$group)]; out$to <- to[as.integer(out$group)]
  } else out$state <- as.integer(out$group)
  out
}

dvm <- function(x, kappa, mu = 0) exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, 0))

eco_panel <- function(p, opts) p + theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic", opts$legend %||% "right")

# Legend labels wrapped so long state names don't eat half of a
# single-column figure (render check 2026-10-20).
# Transition labels ("A → B") break after the arrow.
wrap_lab <- function(x, width = 16) vapply(as.character(x), function(s) {
  parts <- if (nchar(s) > width && grepl(" → ", s)) strsplit(sub(" → ", " →\r", s), "\r")[[1]] else s
  paste(unlist(lapply(parts, strwrap, width)), collapse = "\n")
}, "", USE.NAMES = FALSE)
state_lty <- c("solid", "42", "13", "1343", "73", "22")
state_shape <- c(16, 17, 15, 18, 3, 4)
# guide with a gap between (wrapped, two-line) legend entries
leg_guide <- function(...) ggplot2::guide_legend(theme = ggplot2::theme(legend.key.spacing.y = ggplot2::unit(3, "pt")), ...)

# Figure suite (docs/FIGURES.md, movement-hmm row): state_distributions,
# transition_<tod|cov>, stationary_<tod|cov>, decoded_track, pseudo_residuals,
# plus coefficients. Colour follows the STATE everywhere (Okabe-Ito order =
# state order, states ordered by mean step); transitions take the colour of
# the state they leave. Linetype / shape duplicate the encoding for greyscale.
figures <- function(fit, data, opts = fig_opts()) {
  eco_require(c("ggplot2", "patchwork"))
  obj <- fit$object; K <- obj$K
  ci <- opts$ci %||% 0.95; inner <- opts$inner_ci %||% 0.5
  a <- (1 - ci) / 2; a50 <- (1 - inner) / 2
  cip <- round(100 * ci); inp <- round(100 * inner)
  dt <- obj$prep$interval_h; ph0 <- obj$prep$phase0 %||% 0
  pal <- eco_categorical(max(K, 2), opts$palette %||% "okabe_ito")[seq_len(K)]
  lty <- rep(state_lty, length.out = K); shp <- rep(state_shape, length.out = K)
  slabs <- vapply(seq_len(K), function(k) state_label(obj, k, opts), "")
  short <- sub("\\s*\\(.*\\)$", "", slabs)
  sv <- function(v) stats::setNames(v, slabs)
  sscales <- list(
    ggplot2::scale_colour_manual(values = sv(pal), labels = wrap_lab, name = NULL, guide = leg_guide()),
    ggplot2::scale_fill_manual(values = sv(pal), labels = wrap_lab, name = NULL, guide = leg_guide()),
    ggplot2::scale_linetype_manual(values = sv(lty), labels = wrap_lab, name = NULL, guide = leg_guide()))
  dens_lab <- fig_label("density", "Probability density", opts)
  # Wrap a long y-axis title onto two lines when the figure's export height
  # (global, module hint or per_figure override) is too short for it
  # (review 2026-10-20: 'Long-run share of time in state' was clipped at
  # per_figure {width: single, aspect: 0.6}).
  fit_ytitle <- function(lab, name, hint = NULL) {
    pf <- opts$per_figure[[name]]
    w <- pf$width %||% hint$width %||% opts$width %||% "single"
    asp <- pf$aspect %||% hint$aspect %||% opts$aspect %||% 0.75
    wmm <- tryCatch(fig_width_mm(w, opts), error = function(e) NULL)
    if (is.null(wmm) || is.null(opts$.journal)) return(lab)
    h <- min(pf$height %||% (wmm * asp), opts$.journal$max_h %||% Inf)
    nmax <- floor((h - 10) / (0.16 * (pf$base_size %||% opts$base_size %||% 9)))
    if (nchar(lab) <= nmax || grepl("\n", lab)) return(lab)
    paste(strwrap(lab, max(10, ceiling(nchar(lab) / 2) + 2)), collapse = "\n")
  }
  grey <- isTRUE(opts$greyscale)
  figs <- list()

  # ---- state-dependent distributions (with bands from the parameter draws)
  dr <- obj$draws
  nat <- w2nat(dr$est, K, obj$zm, obj$bdim)
  step <- obj$prep$pd$step; ang <- obj$prep$pd$angle
  st <- step[is.finite(step) & step > 0]
  xs <- seq(0.5, stats::quantile(st, 0.99), length.out = 300)
  as <- seq(-pi, pi, length.out = 200)
  # weight = long-run share of time in the state among non-zero steps
  wts <- function(n) { w <- activity_budget(n$beta, obj$terms, K, dt, ph0) * (1 - n$zm); w / sum(w) }
  fstep <- function(n) { w <- wts(n); vapply(seq_len(K), function(k) w[k] * stats::dgamma(xs, n$mu[k]^2 / n$sd[k]^2, n$mu[k] / n$sd[k]^2), numeric(length(xs))) }
  fang <- function(n) { w <- wts(n); vapply(seq_len(K), function(k) w[k] * dvm(as, n$kappa[k]), numeric(length(as))) }
  Wd <- if (isTRUE(dr$ok)) dr$W[seq_len(min(400, nrow(dr$W))), , drop = FALSE] else NULL
  band <- function(f, xg) {
    hat <- f(nat)
    if (!is.null(Wd)) {
      arr <- simplify2array(lapply(seq_len(nrow(Wd)), function(i) f(w2nat(Wd[i, ], K, obj$zm, obj$bdim))))
      q <- function(p) apply(arr, c(1, 2), stats::quantile, p, na.rm = TRUE)
      lo <- q(a); hi <- q(1 - a); lo50 <- q(a50); hi50 <- q(1 - a50)
    } else lo <- hi <- lo50 <- hi50 <- hat * NA
    d <- do.call(rbind, lapply(seq_len(K), function(k) data.frame(x = xg, state = slabs[k], estimate = hat[, k],
      lower = lo[, k], upper = hi[, k], lower50 = lo50[, k], upper50 = hi50[, k])))
    d$state <- factor(d$state, slabs)
    list(d = d, mix = data.frame(x = xg, estimate = rowSums(hat)))
  }
  dpanel <- function(b, hist, bins, boundary = NULL, xlab) {
    eco_panel(ggplot2::ggplot(b$d, ggplot2::aes(x, estimate)) +
      ggplot2::geom_histogram(data = data.frame(v = hist), ggplot2::aes(v, ggplot2::after_stat(density)), inherit.aes = FALSE,
                              bins = bins, boundary = boundary, fill = "grey88", colour = NA) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = state), alpha = 0.18, colour = NA, na.rm = TRUE) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50, fill = state), alpha = 0.3, colour = NA, na.rm = TRUE) +
      ggplot2::geom_line(data = b$mix, colour = eco_pal$ink, linewidth = 0.3) +
      ggplot2::geom_line(ggplot2::aes(colour = state, linetype = state), linewidth = 0.7) +
      sscales + ggplot2::labs(x = xlab, y = dens_lab), opts)
  }
  p1 <- dpanel(band(fstep, xs), st[st <= max(xs)], 60, NULL,
               fig_label("step", sprintf("Step length (m per %.3g h)", dt), opts))
  p2 <- dpanel(band(fang, as), ang[is.finite(ang)], 36, -pi, fig_label("angle", "Turning angle (radians)", opts)) +
    ggplot2::scale_x_continuous(breaks = c(-pi, -pi / 2, 0, pi / 2, pi), labels = c("-π", "-π/2", "0", "π/2", "π"))
  pd <- compose_figures(list(p1, p2), ncol = 2)
  attr(pd, "fig_size") <- list(width = "double", aspect = 0.4)
  figs$state_distributions <- with_caption(pd, sprintf(paste(
    "State-dependent distributions, each weighted by the long-run share of time in that state (activity budget):",
    "(A) gamma step-length densities, (B) von Mises turning-angle densities. Lines: maximum-likelihood estimates;",
    "bands: %d%% and %d%% confidence intervals from %d draws of the asymptotic sampling distribution.",
    "Grey histograms: observed steps and turns; thin black line: the mixture of the states."), cip, inp, if (is.null(Wd)) 0L else nrow(Wd)))

  # ---- transition and stationary curves
  vars <- c(if ("tod_cos" %in% obj$terms) "hour", obj$prep$covs)
  tz <- data$tz %||% "UTC"
  for (v in vars) {
    xl <- if (v == "hour") fig_label("hour", "Hour of day (local time)", opts) else fig_label(v, v, opts)
    key <- if (v == "hour") "tod" else v
    rug <- if (v == "hour") NULL else {
      rx <- stats::na.omit(obj$prep$grid[[v]][obj$prep$grid$fix_ok])
      data.frame(x = fig_x(if (length(rx) > 400) stats::quantile(rx, seq(0, 1, length.out = 400)) else rx, v, opts))
    }
    xsc <- if (v == "hour") ggplot2::scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) else NULL
    avg <- if (v != "hour" && "tod_cos" %in% obj$terms) ", averaged over the 24-h cycle" else ""
    # transitions: colour = state left, linetype distinguishes the curves
    tc <- hmm_curve(obj, v, "transition", ci = ci, opts = opts, inner = inner)
    if (v != "hour") tc$x <- fig_x(tc$x, v, opts)
    tl <- levels(tc$group)
    fr <- tc$from[match(tl, tc$group)]; to <- tc$to[match(tl, tc$group)]
    tlty <- if (K == 2) lty[fr] else lty[to]
    if (K > 2) tc$from_lab <- factor(paste("From", short[tc$from]), paste("From", short))
    ptc <- fig_cond_effect(tc, xl, fit_ytitle(fig_label("transition_prob", sprintf("Probability of switching per %.3g h", dt), opts),
                                              paste0("transition_", key), if (K > 2) list(width = "double", aspect = 0.4)), obs = rug) +
      ggplot2::aes(linetype = group) +
      ggplot2::scale_colour_manual(values = stats::setNames(pal[fr], tl), labels = wrap_lab, name = NULL, guide = leg_guide()) +
      ggplot2::scale_fill_manual(values = stats::setNames(pal[fr], tl), labels = wrap_lab, name = NULL, guide = leg_guide()) +
      ggplot2::scale_linetype_manual(values = stats::setNames(tlty, tl), labels = wrap_lab, name = NULL, guide = leg_guide()) +
      ggplot2::scale_y_continuous(limits = c(0, NA)) + xsc
    if (K > 2) {
      ptc <- ptc + ggplot2::facet_wrap(~from_lab, nrow = 1)
      attr(ptc, "fig_size") <- list(width = "double", aspect = 0.4)
    }
    figs[[paste0("transition_", key)]] <- with_caption(ptc,
      sprintf("Probability of switching from one behavioural state to another in one time step against %s, other covariates at their mean%s. Lines: maximum-likelihood estimates; bands: %d%% and %d%% confidence intervals from the asymptotic sampling distribution. Colour shows the state the animal leaves.%s",
              if (v == "hour") sprintf("time of day (clock time zone %s)", tz) else xl, avg, cip, inp,
              if (!is.null(rug)) " Rug: covariate values at the fixes." else ""))
    sc <- hmm_curve(obj, v, "stationary", ci = ci, opts = opts, inner = inner)
    if (v != "hour") sc$x <- fig_x(sc$x, v, opts)
    per <- v == "hour" && periodic_ok(dt)
    psc <- fig_cond_effect(sc, xl, fit_ytitle(if (per) fig_label("state_prob_hour", "Probability of being in state", opts)
                                   else fig_label("stationary_prob", "Long-run share of time in state", opts), paste0("stationary_", key)),
                           obs = rug, ylim01 = TRUE) +
      ggplot2::aes(linetype = group) + sscales + xsc
    if (per) psc <- psc + ggplot2::geom_point(ggplot2::aes(shape = group), size = 0.9, show.legend = FALSE) +
      ggplot2::scale_shape_manual(values = sv(shp), guide = "none")
    figs[[paste0("stationary_", key)]] <- with_caption(psc,
      if (per) sprintf("Probability of being in each state at each clock hour over the daily cycle (periodically stationary distribution of the fitted time-of-day model; Koslik et al. 2025), habitat covariates at their mean; clock time zone %s. Points: the fix times, joined by lines (maximum-likelihood estimates); bands: %d%% and %d%% confidence intervals.", tz, cip, inp)
      else sprintf("Long-run share of time in each state if %s were held at the value on the x-axis, other covariates at their mean%s. Lines: maximum-likelihood estimates; bands: %d%% and %d%% confidence intervals.%s",
                   if (v == "hour") "the hour of day" else xl,
                   if (v != "hour" && "tod_cos" %in% obj$terms) ", over the daily cycle of the time-of-day effect" else "", cip, inp,
                   if (!is.null(rug)) " Rug: covariate values at the fixes." else ""))
  }

  # ---- decoded tracks
  if (!is.null(obj$viterbi)) {
    g <- obj$prep$grid
    g$state <- factor(slabs[obj$viterbi], slabs)
    ids <- unlist(opts$ids %||% utils::head(unique(g$ID), 6))
    ids <- intersect(as.character(ids), unique(as.character(g$ID)))
    nc <- min(3, length(ids)); nr <- ceiling(length(ids) / nc)
    ratio <- numeric()
    pl <- lapply(ids, function(id) {
      d <- g[g$ID == id, ]
      d$xk <- (d$x - d$x[which(!is.na(d$x))[1]]) / 1000; d$yk <- (d$y - d$y[which(!is.na(d$y))[1]]) / 1000
      # segments only between consecutive grid times (no line across missed fixes)
      d$xend <- c(d$xk[-1], NA); d$yend <- c(d$yk[-1], NA)
      d <- d[!is.na(d$x), ]
      ratio[id] <<- min(2.5, max(0.3, diff(range(d$yk)) / max(diff(range(d$xk)), 1e-6)))
      eco_panel(ggplot2::ggplot(d, ggplot2::aes(xk, yk)) +
        ggplot2::geom_segment(ggplot2::aes(xend = xend, yend = yend, colour = state), linewidth = 0.3, na.rm = TRUE, show.legend = FALSE) +
        # points drawn longest-step state first, so the clustered short-step
        # fixes are not hidden under travelling segments (verification 2026-10-19)
        ggplot2::geom_point(data = d[order(-as.integer(d$state)), ], ggplot2::aes(colour = state, shape = state), size = 0.6) +
        ggplot2::coord_equal() + sscales[1] +
        ggplot2::scale_shape_manual(values = sv(shp), labels = wrap_lab, name = NULL, guide = leg_guide(override.aes = list(size = 1.8))) +
        ggplot2::facet_wrap(~ID) +
        ggplot2::scale_x_continuous(n.breaks = 4) + ggplot2::scale_y_continuous(n.breaks = 4) +
        ggplot2::labs(x = fig_label("track_x", "E (km)", opts), y = fig_label("track_y", "N (km)", opts)), opts)
    })
    # row heights follow the tallest track in each row (coord_equal panels)
    rmax <- vapply(seq_len(nr), function(r) max(ratio[ids[((r - 1) * nc + 1):min(r * nc, length(ids))]]), 0)
    pt <- patchwork::wrap_plots(pl, ncol = nc, guides = "collect", heights = rmax)
    attr(pt, "fig_size") <- list(width = "double", aspect = min(1.3, max(0.35, 0.185 * sum(rmax) + 0.075 * nr)))
    figs$decoded_track <- with_caption(pt,
      "Tracks coloured by the most likely state sequence (Viterbi algorithm) for each animal; coordinates in km relative to the first fix; each segment takes the colour of the state at its start. Decoded states are point estimates that ignore state uncertainty.")
  }

  # ---- pseudo-residuals
  if (!is.null(obj$pres)) {
    sr <- obj$pres$stepRes; ar <- obj$pres$angleRes
    qx <- fig_label("normal_quantile", "Normal quantile", opts)
    qq <- function(r, key, lab) {
      r <- r[is.finite(r)]
      eco_panel(ggplot2::ggplot(data.frame(r = r), ggplot2::aes(sample = r)) +
        ggplot2::stat_qq(size = 0.4, alpha = 0.5, colour = eco_pal$ink) +
        # greyscale: the reference line must not rely on red (review 2026-10-20)
        ggplot2::stat_qq_line(colour = if (grey) eco_pal$null else eco_pal$truth, linetype = if (grey) "22" else "solid") +
        ggplot2::labs(x = qx, y = fig_label(key, lab, opts)), opts)
    }
    ac <- stats::acf(sr, na.action = stats::na.pass, plot = FALSE, lag.max = 24)
    acd <- data.frame(lag = ac$lag[-1], acf = ac$acf[-1])
    bnd <- stats::qnorm(0.975) / sqrt(sum(is.finite(sr)))
    p3 <- eco_panel(ggplot2::ggplot(acd, ggplot2::aes(lag, acf)) +
      ggplot2::geom_hline(yintercept = 0, colour = eco_pal$muted) +
      ggplot2::geom_hline(yintercept = c(-bnd, bnd), linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0), colour = eco_pal$ink) +
      ggplot2::labs(x = fig_label("acf_lag", sprintf("Lag (time steps of %.3g h)", dt), opts),
                    y = fig_label("acf", "ACF, step residuals", opts)), opts)
    pr <- compose_figures(list(qq(sr, "step_residual", "Step pseudo-residual"), qq(ar, "angle_residual", "Angle pseudo-residual"), p3), ncol = 3)
    attr(pr, "fig_size") <- list(width = "double", aspect = 0.36)
    figs$pseudo_residuals <- with_caption(pr, sprintf(
      "Pseudo-residual checks (Zucchini et al. 2016): (A, B) normal QQ plots of step-length and turning-angle pseudo-residuals (points) against a reference line through the quartiles (%s), which they follow if the model is correct; (C) autocorrelation of step pseudo-residuals with 95%% bounds for white noise (dashed). Systematic curvature or strong autocorrelation means the model misses structure (more states, or covariates on the step distributions).",
      if (grey) "grey, dashed" else "red"))
  }

  # ---- coefficients (link scale, per SD for continuous covariates), at opts$ci
  terms <- obj$terms
  if (length(terms)) {
    nmT <- offdiag_names(K)
    Bh <- beta_hat(obj)[-1, , drop = FALSE]
    bd <- beta_draws(obj, n = if (isTRUE(dr$ok)) nrow(dr$W) else 0)
    q <- function(p) if (length(bd)) apply(simplify2array(lapply(bd, function(B) B[-1, , drop = FALSE])), c(1, 2), stats::quantile, p) else Bh * NA
    lo <- q(a); hi <- q(1 - a); lo50 <- q(a50); hi50 <- q(1 - a50)
    vname <- sub("_z$", "", terms)
    plab <- vapply(seq_along(terms), function(i) {
      t <- terms[i]; v <- vname[i]
      if (t == "tod_cos") return(fig_label("tod_cos", "Time of day (cos)", opts))
      if (t == "tod_sin") return(fig_label("tod_sin", "Time of day (sin)", opts))
      s <- obj$prep$scaling[[v]]
      if (grepl("_z$", t) && !is.null(s)) sprintf("%s, per SD = %.3g", fig_label(v, v, opts), s$scale) else fig_label(v, v, opts)
    }, "")
    fr <- as.integer(sub("to.*", "", nmT)); to <- as.integer(sub(".*to", "", nmT))
    glab <- sprintf("%s → %s", short[fr], short[to])
    cf <- do.call(rbind, lapply(seq_along(nmT), function(j) data.frame(param = plab, group = glab[j],
      estimate = Bh[, j], lower = lo[, j], upper = hi[, j], lower50 = lo50[, j], upper50 = hi50[, j])))
    cf$param <- factor(cf$param, rev(plab)); cf$group <- factor(cf$group, glab)
    pcf <- fig_coef(cf, 0, xlab = fig_label("coef_axis", sprintf("Effect on log-odds of switching (%d%% CI)", cip), opts), order = FALSE) +
      ggplot2::aes(shape = group) +
      ggplot2::scale_colour_manual(values = stats::setNames(pal[fr], glab), labels = function(x) wrap_lab(x, 13), name = NULL, guide = leg_guide()) +
      ggplot2::scale_shape_manual(values = stats::setNames(if (K == 2) shp[fr] else shp[to], glab), labels = function(x) wrap_lab(x, 13), name = NULL, guide = leg_guide()) +
      # wrapped row labels so the panel keeps most of a single-column width
      ggplot2::scale_y_discrete(labels = function(x) wrap_lab(x, 18)) +
      ggplot2::scale_x_continuous(n.breaks = 4, guide = ggplot2::guide_axis(check.overlap = TRUE))
    figs$coefficients <- with_caption(pcf,
      sprintf("Transition-probability coefficients on the multinomial-logit scale (log-odds of switching versus staying), per SD for continuous covariates (the SD in original units is given in the label): points are maximum-likelihood estimates, thin and thick bars %d%% and %d%% confidence intervals; colour shows the state the animal leaves. 'Time of day (cos)' and '(sin)' are the coefficients of cos(2π h/24) and sin(2π h/24), h = local clock hour. Dashed line: no effect.", cip, inp))
  }
  figs
}

plots <- function(fit, data) {
  f <- figures(fit, data, fig_opts())
  f[intersect(names(f), c("state_distributions", grep("^stationary_", names(f), value = TRUE),
                          grep("^transition_", names(f), value = TRUE), "decoded_track", "pseudo_residuals"))]
}

fmt_ci <- function(e, d = 2, pct = FALSE) {
  if (pct) sprintf("%.0f%% (95%% CI %.0f-%.0f%%)", 100 * e$estimate, 100 * e$lower, 100 * e$upper)
  else sprintf(paste0("%.", d, "f (95%% CI %.", d, "f to %.", d, "f)"), e$estimate, e$lower, e$upper)
}

interpret <- function(fit, data) {
  obj <- fit$object; K <- obj$K; e <- fit$estimates
  g <- function(p) e[e$param == p, ]
  dt <- obj$prep$interval_h
  out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge cleanly (see Diagnostics). The numbers below should not be interpreted until this is fixed.")
  desc <- vapply(seq_len(K), function(k) {
    mu <- g(paste0("step_mean_m_state", k)); kp <- g(paste0("turn_concentration_state", k))
    pers <- if (kp$estimate < 1) "little directional persistence (frequent sharp turns and reversals)"
            else if (kp$estimate < 3) "moderate directional persistence" else "strong directional persistence (mostly straight-line movement)"
    sprintf("%s: mean step %.0f m per %.3g h (95%% CI %.0f-%.0f m), turning concentration %.2f (%.2f-%.2f), i.e. %s",
            state_label(obj, k), mu$estimate, dt, mu$lower, mu$upper, kp$estimate, kp$lower, kp$upper, pers)
  }, "")
  out <- c(out, paste0("The model identifies ", K, " movement states, ordered by step length. ", paste(desc, collapse = ". "),
                       ". States are statistical summaries of movement; labels such as 'resting/foraging' or 'travelling' are interpretations that should be checked against independent observations."))
  ab <- vapply(seq_len(K), function(k) sprintf("%s %s", state_label(obj, k), fmt_ci(g(paste0("activity_budget_state", k)), pct = TRUE)), "")
  out <- c(out, paste0("Activity budget (expected long-run share of time in each state",
                       if ("tod_cos" %in% obj$terms) ", averaged over the 24-h cycle" else "",
                       if (length(obj$prep$covs)) ", at the average value of the habitat covariates" else "", "): ",
                       paste(ab, collapse = "; "), "."))
  if ("tod_cos" %in% obj$terms) {
    sc <- hmm_curve(obj, "hour", "stationary", ngrid = 25)
    fast <- sc[sc$group == levels(sc$group)[K], ]
    i1 <- which.max(fast$estimate); i0 <- which.min(fast$estimate)
    out <- c(out, sprintf("Time of day: the probability that an animal is in %s is highest around %02d:00 at %.0f%% (95%% CI %.0f-%.0f%%) and lowest around %02d:00 at %.0f%% (%.0f-%.0f%%), with habitat at its mean (clock time zone %s).",
                          state_label(obj, K), round(fast$x[i1]) %% 24, 100 * fast$estimate[i1], 100 * fast$lower[i1], 100 * fast$upper[i1],
                          round(fast$x[i0]) %% 24, 100 * fast$estimate[i0], 100 * fast$lower[i0], 100 * fast$upper[i0], data$tz %||% "UTC"))
  }
  for (v in obj$prep$covs) {
    sc <- obj$prep$scaling[[v]]
    raw <- obj$prep$grid[[v]][obj$prep$grid$fix_ok]   # values at actual fixes, not carried-forward gaps
    q <- stats::quantile(raw, c(0.1, 0.9), na.rm = TRUE)
    st <- stationary_at(obj, v, q)
    lines <- vapply(seq_len(K - 1), function(i) {
      b <- g(sprintf("tpm_%dto%d_%s", i, i + 1, v))
      sprintf("each 1-SD increase in %s (%.3g units) multiplies the odds of switching from %s to %s (versus staying) per time step by %.2f (95%% CI %.2f-%.2f)",
              v, sc$scale, state_label(obj, i), state_label(obj, i + 1), exp(b$estimate * sc$scale), exp(b$lower * sc$scale), exp(b$upper * sc$scale))
    }, "")
    back <- g(sprintf("tpm_%dto%d_%s", K, K - 1, v))
    out <- c(out, sprintf("%s: %s; switching back from %s to %s changes by a factor %.2f (%.2f-%.2f) per SD. Taken together, the model-predicted long-run share of time in %s would be %.0f%% (95%% CI %.0f-%.0f%%) if %s were held at %.3g (10th percentile of the fixes) and %.0f%% (%.0f-%.0f%%) if it were held at %.3g (90th percentile). This is an association with the covariate measured at the animal's location, not proof that %s causes the behaviour.",
                          v, paste(lines, collapse = "; "), state_label(obj, K), state_label(obj, K - 1), exp(back$estimate * sc$scale), exp(back$lower * sc$scale), exp(back$upper * sc$scale),
                          state_label(obj, 1), 100 * st$est[1], 100 * st$lo[1], 100 * st$hi[1], v, q[1],
                          100 * st$est[2], 100 * st$lo[2], 100 * st$hi[2], q[2], v))
  }
  chk <- sprintf("Fitting: %d of %d starting values reached the best log-likelihood (%.1f), so a local optimum is %s.",
                 obj$n_best, obj$nstarts, obj$llk, if (obj$n_best >= 2) "unlikely" else "possible: rerun with more starting values (settings: nstarts)")
  if (!is.null(obj$gof)) chk <- paste(chk, sprintf("Pseudo-residuals: step-length KS test p = %.3f, turning-angle KS test p = %.3f, lag-1 autocorrelation of step residuals %.2f%s.",
                                           obj$gof$step_ks_p, obj$gof$angle_ks_p, obj$gof$step_acf1,
                                           if (abs(obj$gof$step_acf1) > 0.1) " (residual autocorrelation suggests structure the model misses, e.g. an additional state or individual differences)" else ""))
  if (isTRUE(obj$prep$n_missing > 0))
    out <- c(out, sprintf("Missed fixes: %d of %d regular time steps (%.1f%%) had no usable location. The model treats them as missing at random, i.e. unrelated to behaviour. GPS fix success is often lower under dense canopy or while an animal is bedded, so if missed fixes cluster in one habitat or at certain hours, that state's share of time and its habitat association can be misstated. Check whether fix success varies with %s before interpreting small effects.",
                          obj$prep$n_missing, nrow(obj$prep$grid), 100 * obj$prep$n_missing / nrow(obj$prep$grid),
                          if (length(c(obj$prep$covs, obj$terms))) paste(c(obj$prep$covs, if ("tod_cos" %in% obj$terms) "time of day"), collapse = " or ") else "habitat or time of day"))
  out <- c(out, chk,
    "Limits (from simulation validation, see VALIDATION.md): the number of states is an analyst's choice, and information criteria tend to favour too many states (Pohle et al. 2017); a model with too few states merges behaviours and biases every state's parameters. With short tracks (a few hundred fixes in total) or states whose step lengths overlap strongly, covariate effects on switching are imprecise and intervals can under-cover. Even at the validated design, 95% intervals for switching coefficients were slightly too narrow: with no real effect they excluded zero in about 7% of simulated data sets instead of 5%, so treat an effect whose interval only just excludes zero as weak evidence.")
  out
}

stationary_at <- function(obj, v, xs) {
  sc <- obj$prep$scaling[[v]]; K <- obj$K; terms <- obj$terms
  f <- function(B) vapply(xs, function(x) {
    z <- stats::setNames(list((x - sc$center) / sc$scale), v)
    activity_budget(B, terms, K, obj$prep$interval_h, obj$prep$phase0 %||% 0, z)[1]
  }, 0)
  est <- f(beta_hat(obj))
  bd <- beta_draws(obj)
  if (!length(bd)) return(list(est = est, lo = est * NA, hi = est * NA))
  D <- sapply(bd, f)
  D <- matrix(D, nrow = length(xs))
  list(est = est, lo = apply(D, 1, stats::quantile, 0.025), hi = apply(D, 1, stats::quantile, 0.975))
}

methods_text <- function(fit, data, settings) {
  obj <- fit$object; K <- obj$K
  covs <- obj$prep$covs
  sprintf(paste(
    "We fitted a %d-state hidden Markov model to %s movement tracks (%d animals, %d regular %.3g-h time steps) with the R package momentuHMM %s (McClintock & Michelot 2018),",
    "following Langrock et al. (2012). Locations were projected to metres (%s) and %s",
    "Step lengths followed state-specific gamma distributions%s and turning angles von Mises distributions with mean zero.",
    "Transition probabilities were modelled with a multinomial logit link (staying in the current state as reference) as functions of %s;",
    "continuous covariates were standardised before fitting and coefficients are reported per original unit.",
    "The likelihood was maximised by numerical optimisation (nlm) from %d starting values (one data-driven, the rest randomly perturbed) to avoid local optima, and states were ordered by mean step length.",
    "Confidence intervals for derived quantities (state-specific parameters, state probabilities over the daily cycle from the periodically stationary distribution of the time-varying chain (Koslik et al. 2025), and the activity budget, i.e. the long-run share of time in each state, its average over the day) were obtained by drawing %d parameter vectors from the asymptotic normal distribution of the working-scale maximum-likelihood estimates.",
    "The most likely state sequence was obtained with the Viterbi algorithm and model fit was assessed with pseudo-residuals (Zucchini et al. 2016)."),
    K, if (obj$prep$regularise == "crawl") "crawl-regularised" else "GPS", length(unique(obj$prep$grid$ID)), nrow(obj$prep$grid), obj$prep$interval_h,
    as.character(utils::packageVersion("momentuHMM")), data$crs %||% "projected",
    if (obj$prep$regularise == "crawl")
      sprintf("locations at regular %.3g-h times were predicted from a continuous-time correlated random walk fitted to each track with crawl via momentuHMM::crawlWrap (Johnson et al. 2008; single best prediction, so location uncertainty is not propagated); covariates were taken from the nearest fix in time.", obj$prep$interval_h)
    else sprintf("placed on a regular %.3g-h grid; fixes more than %.0f%% of an interval from a grid time were dropped and grid times without a fix (%d) were treated as missing.",
                 obj$prep$interval_h, 100 * (settings$tolerance %||% 0.25), obj$prep$n_missing),
    if (obj$zm) " with a state-specific probability of zero-length steps" else "",
    { tt <- c(if ("tod_cos" %in% obj$terms) "time of day (cosine and sine of 2*pi*hour/24)", covs); if (length(tt)) paste(tt, collapse = " and ") else "no covariates (constant)" },
    settings$nstarts %||% 5, settings$ndraws %||% 2000)
}
