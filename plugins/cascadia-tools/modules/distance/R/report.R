# Interpretation, figures, methods text and references for the distance module.

dist_references <- c(
  "Buckland ST, Anderson DR, Burnham KP, Laake JL, Borchers DL, Thomas L (2001) Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford University Press.",
  "Buckland ST, Rexstad EA, Marques TA, Oedekoven CS (2015) Distance Sampling: Methods and Applications. Springer. doi:10.1007/978-3-319-19219-2",
  "Miller DL, Rexstad E, Thomas L, Marshall L, Laake JL (2019) Distance sampling in R. Journal of Statistical Software 89(1):1-28. doi:10.18637/jss.v089.i01",
  "Marques TA, Thomas L, Fancy SG, Buckland ST (2007) Improving estimates of bird density using multiple-covariate distance sampling. The Auk 124:1229-1243. doi:10.1642/0004-8038(2007)124[1229:IEOBDU]2.0.CO;2",
  "Fewster RM, Buckland ST, Burnham KP, Borchers DL, Jupp PE, Laake JL, Thomas L (2009) Estimating the encounter rate variance in distance sampling. Biometrics 65:225-236. doi:10.1111/j.1541-0420.2008.01018.x",
  "Innes S, Heide-Jorgensen MP, Laake JL, Laidre KL, Cleator HJ, Richard P, Stewart REA (2002) Surveys of belugas and narwhals in the Canadian High Arctic in 1996. NAMMCO Scientific Publications 4:169-190.",
  "Royle JA, Dawson DK, Bates S (2004) Modeling abundance effects in distance sampling. Ecology 85:1591-1597. doi:10.1890/03-3127",
  "Chandler RB, Royle JA, King DI (2011) Inference about density and temporary emigration in unmarked populations. Ecology 92:1429-1435. doi:10.1890/10-2433.1",
  "Fiske I, Chandler R (2011) unmarked: an R package for fitting hierarchical models of wildlife occurrence and abundance. Journal of Statistical Software 43(10):1-23. doi:10.18637/jss.v043.i10",
  "Buckland ST, Burnham KP, Augustin NH (1997) Model selection: an integral part of inference. Biometrics 53:603-618. doi:10.2307/2533961",
  "Hayes RJ, Buckland ST (1983) Radial-distance models for the line-transect method. Biometrics 39:29-42. doi:10.2307/2530804",
  "Kery M, Royle JA (2016) Applied Hierarchical Modeling in Ecology, Vol. 1, chapter 8 (distance sampling). Academic Press.")

# ---------------------------------------------------------------- helpers
dist_row <- function(est, p) { r <- est[est$param == p, ]; if (nrow(r)) r[1, ] else NULL }
dist_fmt <- function(x, d = NULL) {
  if (!is.finite(x)) return("NA")
  if (is.null(d)) d <- if (abs(x) >= 100) 0 else if (abs(x) >= 10) 1 else if (abs(x) >= 1) 2 else 3
  formatC(x, format = "f", digits = d, big.mark = ",")
}
dist_ci <- function(r, d = NULL) sprintf("%s (95%% CI %s to %s)", dist_fmt(r$estimate, d), dist_fmt(r$lower, d), dist_fmt(r$upper, d))
dist_keyname <- function(k) c(hn = "half-normal", hr = "hazard-rate", unif = "uniform")[[k]]
dist_lab <- function(v, opts) fig_label(v, gsub("_", " ", v), opts)

dist_main_part <- function(fit) { p <- fit$object$parts; p[[1]] }
dist_hds_part <- function(fit) { p <- fit$object$parts; if (!is.null(p$hds)) p$hds else if (!is.null(p$gds)) p$gds else NULL }

# ---------------------------------------------------------------- interpret
dist_interpret <- function(fit, data) {
  est <- fit$estimates
  R <- function(p) dist_row(est, p)
  pr <- fit$prepared
  sp <- data$species %||% "animals"
  gw <- data$group_word %||% "animal"
  eng <- fit$object$engine
  main <- dist_main_part(fit)
  w <- data$truncation_m
  out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge cleanly (see Diagnostics). Treat the numbers below as unreliable until this is fixed.")
  nobs <- nrow(pr$obs_w)
  ntr <- nrow(pr$tz)
  eff <- sum(data$transects$effort) * data$n_visits
  out <- c(out, sprintf("Data: %d detections of %s within the %s m truncation distance from %d %s (%s).",
                        nobs, if (gw == "animal") sp else paste(sp, paste0(gw, "s")), dist_fmt(w, 0), ntr,
                        if (data$survey == "line") "line transects" else "points",
                        if (data$survey == "line") paste(dist_fmt(eff, 1), "km of line walked") else paste(dist_fmt(eff, 0), "point visits")))
  # detection
  pa <- R("detection_prob")
  sel <- main$selection
  kd <- if (!is.null(main$key)) dist_keyname(main$key) else "half-normal"
  adjtxt <- if (!is.null(main$nadj) && main$nadj > 0) sprintf(" with %d cosine adjustment term(s)", main$nadj) else ""
  ma_txt <- if (isTRUE(fit$object$settings$model_average) && !is.null(sel) && nrow(sel) > 1)
    sprintf(" Density, abundance and detection probability in this report are averaged over the %d candidate models with AIC weights (weight of the best model %.2f), so the uncertainty about the shape of the detection function is included.",
            nrow(sel), max(exp(-sel$dAIC / 2) / sum(exp(-sel$dAIC / 2)))) else ""
  out <- c(out, sprintf("Detection: AIC selected a %s detection function%s%s. %s %s within %s m of the %s had probability %s of being detected. %s",
                        kd, adjtxt, if (!is.null(sel) && nrow(sel) > 1) sprintf(" (next best model %.1f AIC units worse)", sort(sel$dAIC)[2]) else "",
                        if (grepl("^[aeiou]", gw)) "An" else "A", gw, dist_fmt(w, 0), if (data$survey == "line") "line" else "point", dist_ci(pa, 2),
                        if (data$survey == "line" && !is.null(R("esw_m"))) sprintf("Equivalently, the effective strip half-width is %s m: the survey counted as many %ss as if every %s within %s m of the line had been seen and none beyond.",
                                                                                 dist_fmt(R("esw_m")$estimate, 0), gw, gw, dist_fmt(R("esw_m")$estimate, 0))
                        else if (!is.null(R("edr_m"))) sprintf("Equivalently, the effective detection radius is %s m: as many %ss were counted as if every %s within %s m of the point had been seen and none beyond.",
                                                               dist_fmt(R("edr_m")$estimate, 0), gw, gw, dist_fmt(R("edr_m")$estimate, 0)) else ""))
  out[length(out)] <- paste0(out[length(out)], ma_txt)
  if (!is.null(R("availability_prob")))
    out <- c(out, sprintf("Availability: on a given visit, a(n) %s using the area around a transect was present and available for detection with probability %s (temporary emigration; Chandler et al. 2011). Density below refers to the whole (super)population using the transect areas.",
                          gw, dist_ci(R("availability_prob"), 2)))
  # detection covariates
  db <- est[grepl("^detection_beta_", est$param), ]
  for (i in seq_len(nrow(db))) {
    nm <- sub("^detection_beta_", "", db$param[i])
    out <- c(out, sprintf("Detection covariate %s: the scale of the detection function is multiplied by %s (95%% CI %s to %s) %s; this is an effect on how far away %ss are detected, not on how many there are.",
                          nm, dist_fmt(exp(db$estimate[i]), 2), dist_fmt(exp(db$lower[i]), 2), dist_fmt(exp(db$upper[i]), 2),
                          if (nm %in% names(pr$scaling)) "per 1 SD of the covariate" else "relative to the reference level", gw))
  }
  # density and abundance
  D <- R("density_per_km2"); N <- R("abundance")
  unitw <- if (isTRUE(data$clusters)) "individuals" else paste0(sp)
  if (!is.null(D)) out <- c(out, sprintf("Density: %s %s per km2%s.", dist_ci(D), unitw,
                                         if (!is.null(N)) sprintf(", which over the %s km2 study area is %s %s", dist_fmt(ds_total_area(data), 0), dist_ci(N, 0), unitw) else ""))
  if (!is.null(R("group_density_per_km2")) && !is.null(R("mean_group_size")))
    out <- c(out, sprintf("This is %s %ss per km2 times a mean %s size of %s.", dist_ci(R("group_density_per_km2")), gw, gw, dist_ci(R("mean_group_size"), 2)))
  st <- est[grepl("^density_per_km2_", est$param) & est$param != "density_per_km2_hds", ]
  if (nrow(st)) out <- c(out, paste0("By stratum: ", paste(sprintf("%s %s per km2", sub("^density_per_km2_", "", st$param), vapply(seq_len(nrow(st)), function(i) dist_ci(st[i, ]), "")), collapse = "; "), "."))
  if (!is.null(D) && is.finite(D$estimate)) {
    cv <- D$sd / D$estimate
    er <- R("encounter_rate_per_km") %||% R("encounter_rate_per_point")
    out <- c(out, sprintf("Precision: the coefficient of variation of density is %.0f%%.%s", 100 * cv,
                          if (!is.null(er)) sprintf(" Most of this usually comes from variation in the number of detections between transects (encounter rate %s per %s, CV %.0f%%).",
                                                     dist_fmt(er$estimate, 2), if (data$survey == "line") "km" else "point visit", 100 * er$sd / er$estimate) else ""))
  }
  # HDS covariate effects
  bb <- est[grepl("^density_beta_", est$param), ]
  labs <- data$labels %||% list()
  for (i in seq_len(nrow(bb))) {
    v <- sub("^density_beta_", "", bb$param[i])
    vl <- if (!is.null(labs[[v]])) labs[[v]] else v
    vn <- sub("^([A-Z])([a-z])", "\\L\\1\\E\\2", sub(" \\(.*\\)$", "", vl), perl = TRUE); vu <- if (grepl("\\(.*\\)$", vl)) sub("^.*\\((.*)\\)$", "\\1", vl) else ""
    sc <- pr$scaling[[v]]
    if (!is.null(sc)) {
      step <- signif(sc$scale, 1)
      rr <- exp(bb$estimate[i] * step / sc$scale); rl <- exp(bb$lower[i] * step / sc$scale); ru <- exp(bb$upper[i] * step / sc$scale)
      xr <- range(data$transects[[v]])
      h <- dist_hds_part(fit)
      lo <- hi <- NA
      if (!is.null(h)) {
        zr <- (xr - sc$center) / sc$scale
        jb <- which(colnames(stats::model.matrix(stats::as.formula(ds_rhs(data$formulas$density)), pr$tz)) == v)
        if (length(jb)) { lo <- exp(h$theta[1] + h$theta[jb] * zr[1]) * mean(h$size_draws); hi <- exp(h$theta[1] + h$theta[jb] * zr[2]) * mean(h$size_draws) }
      }
      dirn <- if (bb$lower[i] > 0) "increases with" else if (bb$upper[i] < 0) "decreases with" else "is not clearly related to"
      out <- c(out, sprintf("Density and %s: expected density %s %s. Each increase of %s%s in %s multiplies expected density by %s (95%% CI %s to %s), i.e. %+.0f%% (1 SD = %s; the table's rate_ratio_%s is per 1 SD).%s This is an association from an observational survey, not proof that %s causes the difference.",
                            vn, dirn, vn, format(step), if (nzchar(vu)) paste0(" ", if (vu == "%") "percentage points" else vu) else "", vn,
                            dist_fmt(rr, 2), dist_fmt(rl, 2), dist_fmt(ru, 2), 100 * (rr - 1), format(signif(sc$scale, 3)), v,
                            if (is.finite(lo)) sprintf(" Across the surveyed range (%s to %s%s), with other covariates at their mean, predicted density goes from %s to %s per km2.",
                                                       format(signif(xr[1], 3)), format(signif(xr[2], 3)), if (nzchar(vu)) paste0(" ", vu) else "", dist_fmt(lo), dist_fmt(hi)) else "", vn))
    } else out <- c(out, sprintf("Density effect %s (factor level vs reference): rate ratio %s (95%% CI %s to %s).", v,
                                 dist_fmt(exp(bb$estimate[i]), 2), dist_fmt(exp(bb$lower[i]), 2), dist_fmt(exp(bb$upper[i]), 2)))
  }
  if (!is.null(R("density_per_km2_hds")) && !is.null(D))
    out <- c(out, sprintf("The hierarchical model's density over the study area (%s per km2) %s the design-based estimate above.",
                          dist_ci(R("density_per_km2_hds")),
                          if (abs(R("density_per_km2_hds")$estimate / D$estimate - 1) < 0.15) "agrees with" else "differs from (check that transects represent the study area and that the habitat model is adequate)"))
  # GOF
  g <- main$gof
  if (!is.null(g) && is.finite(g$p %||% NA))
    out <- c(out, sprintf("Goodness of fit: %s p = %.2f. %s", switch(g$type, cvm = "Cramer-von Mises test of the distances,", chisq = "chi-square test of binned distances,", freeman_tukey = "parametric-bootstrap Freeman-Tukey test,"), g$p,
                          if (g$p < 0.05) "The detection function does not fit the distances well: look at the Q-Q plot and histogram for heaping (many 0 or rounded distances) or a dip near the line (animals moving away), and consider binning." else "No evidence of lack of fit."))
  out <- c(out, "Assumptions behind these numbers: every animal on the line or point is detected (g(0) = 1), animals are detected at their initial location (no movement towards or away from the observer before detection), distances are measured accurately, and transects were placed independently of where the animals are. Simulations with this module show that evasive movement or heaping of distances at zero produce biased densities that the confidence intervals do not cover.")
  if (nobs < 60) out <- c(out, sprintf("With only %d detections the detection function is poorly determined; simulations with this module show that density intervals still cover well with ~50 detections but are wide, and AIC model choice becomes unstable.", nobs))
  c(out, unlist(lapply(data$warnings %||% dist_data_warnings(data), identity)))
}

# ---------------------------------------------------------------- detection curve
# Fitted detection probability g(x) over 0..w, with a band from parameter draws.
dist_detcurve <- function(fit, data, n = 101, ndraw = 300) {
  main <- dist_main_part(fit)
  pr <- fit$prepared
  w <- data$truncation_m
  xs <- seq(0, w, length.out = n)
  if (identical(main$engine, "Distance::ds")) {
    ddf <- main$fit$ddf
    par <- ddf$par
    V <- tryCatch(solve(ddf$hessian), error = function(e) NULL)
    drw <- if (is.null(V)) matrix(par, 1) else ds_mvn(par, V, ndraw)
    covs <- length(pr$dt_terms) > 0
    ks <- summary(ddf)$coeff$key.scale
    if (!covs && main$nadj > 0) {
      g1 <- function(p) {
        obj <- mrds:::assign.par(ddf$ds$aux$ddfobj, p)
        as.numeric(mrds:::detfct(xs, obj, width = w, standardize = TRUE))
      }
      G <- tryCatch(t(apply(drw, 1, g1)), error = function(e) NULL)
      if (is.null(G)) G <- t(apply(drw, 1, function(p) ds_g(xs, main$key, exp(p[length(p)]), if (main$key == "hr") exp(p[1]))))
    } else {
      # key only (with or without covariates): sigma_i = exp(Z_i alpha); average over observed covariates
      Z <- stats::model.matrix(stats::as.formula(ds_rhs(data$formulas$detection)), pr$obs_w)
      nsc <- ncol(Z)
      ishape <- if (main$key == "hr") 1 else NULL
      iscale <- if (main$key == "hr") 1 + seq_len(nsc) else seq_len(nsc)
      g1 <- function(p) {
        sig <- exp(Z %*% p[iscale]); b <- if (!is.null(ishape)) exp(p[ishape]) else NULL
        rowMeans(vapply(unique(as.numeric(sig)), function(s) ds_g(xs, main$key, s, b), numeric(length(xs))) |> as.matrix())
      }
      G <- t(apply(drw, 1, g1))
      G[1, ] <- g1(par)
    }
    pa <- summary(ddf)$average.p
  } else {
    h <- dist_hds_part(fit) %||% NULL
    if (is.null(h)) return(NULL)
    dr <- h$draws
    ia <- h$nb + 1
    b <- if (h$key == "hr") exp(dr[, h$nb + h$na + 1]) else rep(NA, nrow(dr))
    G <- t(vapply(seq_len(min(nrow(dr), ndraw)), function(i) ds_g(xs, h$key, exp(dr[i, ia]), b[i]), numeric(length(xs))))
    G[1, ] <- ds_g(xs, h$key, exp(h$theta[ia]), if (h$key == "hr") exp(h$theta[h$nb + h$na + 1]))
    pa <- fit$estimates$estimate[fit$estimates$param == "detection_prob"]
  }
  data.frame(x = xs, estimate = G[1, ], lower = apply(G, 2, stats::quantile, 0.025),
             upper = apply(G, 2, stats::quantile, 0.975),
             lower50 = apply(G, 2, stats::quantile, 0.25), upper50 = apply(G, 2, stats::quantile, 0.75),
             pa = pa)
}

dist_hist <- function(fit, data, pa) {
  pr <- fit$prepared
  w <- data$truncation_m
  o <- pr$obs_w
  if (data$distance_type == "binned") cut <- data$cutpoints[data$cutpoints <= w + 1e-9]
  else { nb <- max(6, min(20, round(sqrt(nrow(o))))); cut <- seq(0, w, length.out = nb + 1) }
  x <- if (data$distance_type == "binned") (o$distbegin + o$distend) / 2 else o$distance
  nj <- as.numeric(table(factor(findInterval(x, cut, rightmost.closed = TRUE, all.inside = TRUE), levels = seq_len(length(cut) - 1))))
  frac <- if (data$survey == "line") diff(cut) / w else diff(cut^2) / w^2
  data.frame(xmin = utils::head(cut, -1), xmax = cut[-1], n = nj, height = (nj / sum(nj)) / frac * pa)
}

# ---------------------------------------------------------------- figures
dist_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  figs <- list()
  pr <- fit$prepared
  est <- fit$estimates
  w <- data$truncation_m
  gw <- data$group_word %||% "animal"
  ci <- opts$ci %||% 0.95
  dlab <- if (data$survey == "line") fig_label("distance", "Perpendicular distance (m)", opts) else fig_label("distance", "Radial distance from point (m)", opts)
  # 1. detection function
  dc <- tryCatch(dist_detcurve(fit, data), error = function(e) { message("detection curve: ", conditionMessage(e)); NULL })
  if (!is.null(dc)) {
    hb <- dist_hist(fit, data, dc$pa[1])
    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = hb, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = height),
                         fill = "#d9e3dd", colour = "white", linewidth = 0.3) +
      ggplot2::geom_ribbon(data = dc, ggplot2::aes(x, ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.18) +
      ggplot2::geom_ribbon(data = dc, ggplot2::aes(x, ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25) +
      ggplot2::geom_line(data = dc, ggplot2::aes(x, estimate), colour = eco_pal$band, linewidth = 0.8) +
      ggplot2::scale_y_continuous(limits = c(0, max(1.05, max(hb$height) * 1.05)), expand = ggplot2::expansion(mult = c(0, 0.01))) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = dlab, y = fig_label("detection_probability", "Detection probability", opts))
    if (data$distance_type == "exact" && isTRUE(opts$show_data))
      p <- p + ggplot2::geom_rug(data = data.frame(x = pr$obs_w$distance), ggplot2::aes(x = x), inherit.aes = FALSE,
                                 alpha = 0.3, sides = "b", length = ggplot2::unit(0.02, "npc"))
    main <- dist_main_part(fit)
    figs$detection_function <- with_caption(p, sprintf(
      "Fitted detection function (%s%s; line) with %d%% and 50%% confidence bands from the parameter covariance, over a histogram of the %d observed %s distances scaled so that bar heights estimate detection probability in each distance band (a bar above the curve means more detections there than the model expects). %s",
      dist_keyname(main$key %||% "hn"), if (isTRUE(main$nadj > 0)) " with cosine adjustments" else "", round(100 * ci),
      nrow(pr$obs_w), if (data$survey == "line") "perpendicular" else "radial",
      if (length(pr$dt_terms)) "With detection covariates the curve is averaged over the observed covariate values." else ""))
  }
  # 2. goodness of fit
  g <- dist_main_part(fit)$gof
  if (!is.null(g) && identical(g$type, "cvm") && !is.null(g$edf)) {
    q <- data.frame(fitted = as.numeric(g$cdf), empirical = (g$edf[, 1] + g$edf[, 2]) / 2)
    figs$qq_gof <- with_caption(
      ggplot2::ggplot(q, ggplot2::aes(fitted, empirical)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = eco_pal$null) +
        ggplot2::geom_point(size = 0.9, colour = eco_pal$ink, alpha = 0.7) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "Fitted cumulative distribution", y = "Empirical cumulative distribution"),
      sprintf("Q-Q plot of the detection function: each point is one detection; points on the dashed line mean the fitted distribution of distances matches the observed one. Cramer-von Mises W = %.3f, p = %.2f.", g$stat, g$p))
  } else if (!is.null(fit$object$parts$hds) || !is.null(fit$object$parts$gds) || data$distance_type == "binned") {
    h <- dist_hds_part(fit)
    if (!is.null(h)) {
      yo <- colSums(matrix(unmarked::getY(h$fit@data), nrow = nrow(pr$tz)))
      ye <- colSums(matrix(unmarked::fitted(h$fit), nrow = nrow(pr$tz)), na.rm = TRUE)
      J <- length(h$cut) - 1
      yo <- tapply(yo, rep(seq_len(J), length.out = length(yo)), sum); ye <- tapply(ye, rep(seq_len(J), length.out = length(ye)), sum)
      bd <- data.frame(bin = factor(sprintf("%g-%g", utils::head(h$cut, -1), h$cut[-1]), levels = sprintf("%g-%g", utils::head(h$cut, -1), h$cut[-1])),
                       observed = as.numeric(yo), expected = as.numeric(ye))
      figs$gof_bins <- with_caption(
        ggplot2::ggplot(bd, ggplot2::aes(bin)) + ggplot2::geom_col(ggplot2::aes(y = observed), fill = "#d9e3dd") +
          ggplot2::geom_point(ggplot2::aes(y = expected), colour = eco_pal$ink, size = 1.6) +
          ggplot2::labs(x = paste(sub(" \\(m\\)", "", dlab), "band (m)"), y = paste("Number of", gw, "detections")) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
        sprintf("Observed (bars) and expected (points) detections per distance band from the hierarchical model%s.",
                if (!is.null(h$gof)) sprintf("; parametric-bootstrap Freeman-Tukey p = %.2f", h$gof$p) else ""))
      if (!is.null(h$gof)) figs$gof_ppc <- with_caption(fig_ppc(h$gof$rep, h$gof$stat, "Freeman-Tukey statistic"),
        "Parametric bootstrap goodness of fit of the hierarchical model: histogram of the Freeman-Tukey statistic for data simulated from the fitted model; line = observed data.")
    }
  }
  # 3. conditional effects on density (HDS)
  h <- dist_hds_part(fit)
  if (!is.null(h) && length(pr$st_terms)) {
    Xs <- stats::model.matrix(stats::as.formula(ds_rhs(data$formulas$density)), pr$tz)
    dr <- h$draws[, seq_len(h$nb), drop = FALSE]
    sz <- h$size_draws
    unit <- if (isTRUE(data$clusters)) "Density (individuals per km²)" else fig_label("density_unit", "Density (per km²)", opts)
    # per-transect naive density for show_data: n_i * mean size / (covered area_i * P_a)
    a_i <- ds_cov_area(data, max(h$cut))
    pa <- est$estimate[est$param == "detection_prob"][1]
    n_i <- as.numeric(table(factor(pr$obs_w$transect, levels = data$transects$transect)))
    dn <- n_i * mean(sz) / (a_i * pa)
    effs <- list()
    for (v in pr$st_terms) {
      xv <- pr$tz[[v]]
      if (is.factor(xv)) {
        lev <- levels(xv)
        nd <- pr$tz[rep(1, length(lev)), , drop = FALSE]
        for (u in pr$st_terms) if (u != v) nd[[u]] <- if (is.factor(pr$tz[[u]])) factor(levels(pr$tz[[u]])[1], levels(pr$tz[[u]])) else 0
        nd[[v]] <- factor(lev, levels = lev)
        Xn <- stats::model.matrix(stats::delete.response(stats::terms(stats::as.formula(ds_rhs(data$formulas$density)))), nd)
        Dd <- exp(dr %*% t(Xn)) * sz
        pdf <- data.frame(level = lev, estimate = apply(Dd, 2, stats::median), lower = apply(Dd, 2, stats::quantile, (1 - ci) / 2),
                          upper = apply(Dd, 2, stats::quantile, 1 - (1 - ci) / 2))
        p <- fig_cat_effect(pdf, dist_lab(v, opts), unit,
                            obs = if (isTRUE(opts$show_data)) data.frame(level = as.character(xv), y = dn) else NULL)
      } else {
        zg <- seq(min(xv), max(xv), length.out = 100)
        nd <- pr$tz[rep(1, 100), , drop = FALSE]
        for (u in pr$st_terms) if (u != v) nd[[u]] <- if (is.factor(pr$tz[[u]])) factor(levels(pr$tz[[u]])[1], levels(pr$tz[[u]])) else 0
        nd[[v]] <- zg
        Xn <- stats::model.matrix(stats::delete.response(stats::terms(stats::as.formula(ds_rhs(data$formulas$density)))), nd)
        Dd <- exp(dr %*% t(Xn)) * sz
        a <- (1 - ci) / 2; ia <- (1 - (opts$inner_ci %||% 0.5)) / 2
        pdf <- data.frame(x = fig_x(zg, v, opts, pr$scaling), estimate = apply(Dd, 2, stats::median),
                          lower = apply(Dd, 2, stats::quantile, a), upper = apply(Dd, 2, stats::quantile, 1 - a),
                          lower50 = apply(Dd, 2, stats::quantile, ia), upper50 = apply(Dd, 2, stats::quantile, 1 - ia))
        xo <- fig_x(xv, v, opts, pr$scaling)
        p <- fig_cond_effect(pdf, dist_lab(v, opts), unit,
                             obs = if (isTRUE(opts$show_data)) data.frame(x = xo, y = dn) else data.frame(x = xo))
      }
      effs[[v]] <- p
      figs[[paste0("density_", v)]] <- with_caption(p, sprintf(
        "Expected density versus %s from the hierarchical distance-sampling model (line: median; bands: %d%% and 50%% confidence bands from the parameter covariance)%s; other covariates at their mean or reference level.%s",
        dist_lab(v, opts), round(100 * ci), if (isTRUE(data$clusters)) ", multiplied by mean group size" else "",
        if (isTRUE(opts$show_data)) " Points are per-transect density estimates (detections x mean group size / (covered area x average detection probability)); they are noisy and shown only to display the data." else " Rug: observed covariate values."))
    }
    if (length(effs) > 1) {
      pe <- compose_figures(effs, ncol = min(3, length(effs))) & theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic")
      attr(pe, "fig_size") <- list(width = "double", height = fig_width_mm("double", opts) * 0.42 * ceiling(length(effs) / 3))
      figs$effects_panel <- with_caption(pe, "Conditional effects of each covariate on density (panels as in the individual figures).")
    }
    cf <- est[grepl("^density_beta_", est$param), c("param", "estimate", "lower", "upper")]
    if (nrow(cf)) {
      cf$param <- sub("^density_beta_", "", cf$param)
      cf$param <- vapply(cf$param, function(v) dist_lab(v, opts), "")
      figs$coefficients <- with_caption(fig_coef(cf, null_value = 1, transform = exp, log_x = TRUE,
                                                 xlab = fig_label("rate_ratio_axis", "Density ratio per 1 SD", opts), opts = NULL),
        "Multiplicative change in expected density per 1 SD increase of each covariate (factor levels: relative to the reference level), with 95% Wald confidence intervals; dashed line = no effect.")
    }
  }
  # 4. density by stratum
  st <- est[grepl("^density_per_km2_", est$param) & est$param != "density_per_km2_hds", ]
  if (nrow(st) > 1) {
    tot <- est[est$param == "density_per_km2", ]
    pdf <- data.frame(level = c(sub("^density_per_km2_", "", st$param), "Whole study area"),
                      estimate = c(st$estimate, tot$estimate), lower = c(st$lower, tot$lower), upper = c(st$upper, tot$upper))
    pdf$level <- factor(vapply(pdf$level, function(v) fig_label(v, v, opts), ""), levels = vapply(pdf$level, function(v) fig_label(v, v, opts), ""))
    figs$density_by_stratum <- with_caption(
      fig_cat_effect(pdf, fig_label("stratum", "Stratum", opts), fig_label("density_unit", "Density (per km²)", opts)) +
        ggplot2::expand_limits(y = 0),
      "Design-based density per stratum and for the whole study area (points) with 95% log-normal confidence intervals (encounter-rate variance by the R2/P2 estimator of Fewster et al. 2009, detection-function and group-size variance combined).")
  }
  # 5. density map from the prediction grid
  if (!is.null(h) && !is.null(data$grid) && all(c("x", "y") %in% names(data$grid))) {
    G <- data$grid
    gz <- G
    for (v in pr$st_terms) gz[[v]] <- if (is.factor(pr$tz[[v]])) factor(G[[v]], levels(pr$tz[[v]])) else (G[[v]] - pr$scaling[[v]]$center) / pr$scaling[[v]]$scale
    Xg <- stats::model.matrix(stats::delete.response(stats::terms(stats::as.formula(ds_rhs(data$formulas$density)))), gz)
    Dg <- exp(h$draws[, seq_len(h$nb), drop = FALSE] %*% t(Xg)) * h$size_draws
    grid <- data.frame(x = G$x, y = G$y, value = colMeans(Dg))
    sdg <- data.frame(x = G$x, y = G$y, value = apply(Dg, 2, stats::sd))
    tx <- intersect(c("x", "x_km", "easting"), names(data$transects))[1]
    ty <- intersect(c("y", "y_km", "northing"), names(data$transects))[1]
    pts <- if (!is.na(tx) && !is.na(ty)) data.frame(x = data$transects[[tx]], y = data$transects[[ty]]) else NULL
    m1 <- fig_map(grid, "Density\n(per km²)", points = pts)
    m2 <- fig_map(sdg, "SD", points = NULL, option = "magma")
    m1 <- m1 + ggplot2::labs(x = fig_label("map_x", "Easting (km)", opts), y = fig_label("map_y", "Northing (km)", opts))
    m2 <- m2 + ggplot2::labs(x = fig_label("map_x", "Easting (km)", opts), y = NULL)
    pm <- compose_figures(list(m1, m2), ncol = 2, guides = "keep") & theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic")
    asp <- diff(range(G$y)) / diff(range(G$x))
    attr(pm, "fig_size") <- list(width = "double", height = fig_width_mm("double", opts) * min(0.9, 0.42 * asp + 0.06))
    figs$density_map <- with_caption(pm,
      "A: predicted density per km² over the prediction grid from the hierarchical model (mean of parameter draws); white points are surveyed transects/points. B: its standard deviation. Predictions outside the surveyed covariate range are extrapolations.")
  }
  figs
}

dist_plots <- function(fit, data) {
  o <- fig_opts()
  f <- dist_figures(fit, data, o)
  lapply(stats::setNames(names(f), names(f)), function(n) if (inherits(f[[n]], "ggplot")) apply_fig_opts(f[[n]], n, o) else f[[n]])
}

# ---------------------------------------------------------------- methods
dist_methods_text <- function(fit, data, settings) {
  st <- fit$object$settings
  eng <- fit$object$engine
  main <- dist_main_part(fit)
  w <- data$truncation_m
  txt <- sprintf("We estimated density from %s-transect distance-sampling data (%d %s, %d detections within a truncation distance of %s m%s).",
                 data$survey, nrow(data$transects), if (data$survey == "line") "lines" else "points", nrow(fit$prepared$obs_w), dist_fmt(w, 0),
                 if (data$distance_type == "binned") sprintf(", distances recorded in %d bands", length(data$cutpoints) - 1) else "")
  if (eng %in% c("ds", "ds+hds")) txt <- paste(txt, sprintf(
    "Detection functions were fitted with the R package Distance (Miller et al. 2019): %s key functions%s, with the model chosen by AIC (Buckland et al. 2001, 2015)%s. Density and abundance were estimated with the encounter-rate variance estimator %s of Fewster et al. (2009), combined with the detection-function%s variance by the delta method, with log-normal confidence intervals using Satterthwaite degrees of freedom. Fit was assessed with the Cramer-von Mises test and a Q-Q plot.",
    paste(vapply(st$keys, dist_keyname, ""), collapse = " and "),
    if (length(fit$prepared$dt_terms)) sprintf(" with detection covariates on the scale parameter (%s; multiple-covariate distance sampling, Marques et al. 2007)", ds_rhs(data$formulas$detection))
    else if (!identical(st$adjustment, "none")) sprintf(" each with 0-%d %s adjustment terms added by forward selection", st$max_adjustments, if (st$adjustment == "cos") "cosine" else "Hermite polynomial") else "",
    paste0(if (isTRUE(st$model_average) && length(st$keys) > 1) "; density, abundance and detection probability were averaged over the best model of each key family with Akaike weights, with the unconditional standard error of Buckland, Burnham & Augustin (1997)" else "",
           if (!is.null(data$strata) && nrow(data$strata) > 1) ", and density was estimated separately for each stratum" else ""),
    if (data$survey == "line") "R2" else "P2", if (isTRUE(data$clusters)) " and mean group size" else ""))
  if (eng %in% c("ds", "ds+hds") && identical(st$variance, "bootstrap")) txt <- paste(txt, sprintf(
    "Standard errors and 95%% percentile intervals for density, abundance, detection probability and encounter rate came instead from %d nonparametric bootstrap resamples of transects within strata, with the detection-function selection repeated in every resample so that model-selection uncertainty is included (Buckland et al. 2001).", st$nboot))
  if (eng %in% c("hds", "ds+hds")) txt <- paste(txt, sprintf(
    "To relate density to transect covariates we fitted a hierarchical distance-sampling model (Royle et al. 2004) with unmarked::distsamp (Fiske & Chandler 2011): counts in %d distance bands per transect were modelled as Poisson with log(density) = %s (numeric covariates standardised) and a %s detection function (%s)%s. Derived quantities (density at an average site, predicted density over the study area, average detection probability) were computed from %d draws of the multivariate-normal approximation to the MLEs.",
    length(dist_hds_part(fit)$cut) - 1, sub("^~\\s*", "", ds_rhs(data$formulas$density)), dist_keyname(dist_hds_part(fit)$key),
    "chosen by AIC between half-normal and hazard-rate among keys whose detection probability was identified, CV < 50%",
    if (length(fit$prepared$dt_terms)) paste0(" with log(scale) = ", sub("^~\\s*", "", ds_rhs(data$formulas$detection))) else "", st$n_draws))
  if (eng == "gds") txt <- paste(txt, sprintf(
    "Because each transect was surveyed %d times, we used the generalised distance-sampling model of Chandler et al. (2011) (unmarked::gdistsamp): superpopulation size per transect M_i ~ Poisson(density x area), animals available on a visit ~ Binomial(M_i, phi), and detections per distance band multinomial given availability.", data$n_visits))
  if (eng == "rlik") txt <- paste(txt, "Estimates come from hand-written maximum-likelihood code (conditional distance likelihood and Poisson hierarchical likelihood).")
  paste(txt, "Simulation-recovery tests of this workflow (ecotools distance module, VALIDATION.md) are reported alongside.")
}
