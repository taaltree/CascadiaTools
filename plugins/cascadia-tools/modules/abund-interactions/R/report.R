# abund-interactions: interpretation, methods text and figures.

abi_row <- function(est, p) est[est$param == p, , drop = FALSE]
# Is the habitat modifier a factor? Then gamma0 is the interaction in its reference level.
abi_mod_ref <- function(des) {
  if (!length(des$colnames$mod)) return(NULL)
  raw <- des$site_covs_raw[[des$modifier]]
  if (is.numeric(raw)) NULL else levels(factor(raw))[1]
}
abi_hab_txt <- function(des) {
  if (!length(des$colnames$mod)) return("")
  ref <- abi_mod_ref(des)
  if (is.null(ref)) " at average habitat" else sprintf(" where %s = %s (reference level)", des$modifier, ref)
}
# Describe a 0/1 dummy column of a factor site covariate, e.g. "habitat = conifer (vs broadleaf)"
abi_dummy_desc <- function(des, col) {
  for (v in names(des$site_covs_std)) {
    x <- des$site_covs_std[[v]]
    if (!is.factor(x)) next
    lv <- levels(x)
    hit <- which(abi_clean(paste0(v, lv)) == col)
    if (length(hit)) return(sprintf("%s = %s (compared with %s)", v, lv[hit[1]], lv[1]))
  }
  NULL
}
abi_prob <- function(x) if (x > 0.99) "> 0.99" else if (x < 0.01) "< 0.01" else sprintf("%.2f", x)
# effect of a 10-unit change of a raw covariate, from a per-SD coefficient row
abi_per10 <- function(r, sc, what) {
  if (is.null(sc) || !is.finite(sc$scale) || sc$scale <= 0) return("")
  f <- function(x) exp(x * 10 / sc$scale)
  sprintf(" (per 10 units: %s %s)", what, abi_ci(r, 2, f))
}
abi_fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
abi_ci <- function(r, d = 2, f = identity) sprintf("%s (95%% CI %s to %s)", abi_fmt(f(r$estimate), d), abi_fmt(f(r$lower), d), abi_fmt(f(r$upper), d))

# Covariate row of the psi design matrix at "average" conditions, with the
# habitat modifier set to the standardised value(s) `modval` when given.
abi_psi_rows <- function(des, modval = NULL) {
  x <- colMeans(des$Xpsi)
  for (k in colnames(des$Xpsi)) if (all(des$Xpsi[, k] %in% c(0, 1)) && k != "_Intercept_" && k != "Intercept")
    x[k] <- 0                                   # factor dummies: reference level
  m <- if (length(des$colnames$mod)) {
    mv <- if (is.null(modval)) {
      raw <- des$site_covs_raw[[des$modifier]]
      if (is.numeric(raw)) 0 else rep(0, length(des$colnames$mod))
    } else modval
    x[des$colnames$mod] <- mv
    mv
  } else numeric()
  list(x = x, mod = m)
}

# Draws (n x length(Nseq)) of subordinate occupancy vs dominant abundance.
abi_psi_curve <- function(des, draws, Nseq, modval = NULL, xrow = NULL) {
  pr <- abi_psi_rows(des, modval)
  x <- if (is.null(xrow)) pr$x else xrow
  bpsi <- draws[, paste0("beta_psi[", seq_along(des$colnames$psi), "]"), drop = FALSE]
  eta <- drop(bpsi %*% x)
  slope <- draws[, "gamma0"]
  if (length(des$colnames$mod)) slope <- slope + drop(draws[, paste0("gamma_mod[", seq_along(des$colnames$mod), "]"), drop = FALSE] %*% pr$mod)
  stats::plogis(eta + outer(slope, Nseq))
}

abi_band <- function(mat, x, ci = 0.95, inner = 0.5) {
  q <- function(p) apply(mat, 2, stats::quantile, p, na.rm = TRUE)
  data.frame(x = x, estimate = apply(mat, 2, stats::median), lower = q((1 - ci) / 2), upper = q(1 - (1 - ci) / 2),
             lower50 = q((1 - inner) / 2), upper50 = q(1 - (1 - inner) / 2))
}

abi_post_N <- function(fit) {
  o <- fit$object
  if (!is.null(o$N)) colMeans(o$N) else if (!is.null(o$bup)) o$bup else rep(NA_real_, o$des$M)
}

# habitat levels at which interaction curves are drawn
abi_hab_levels <- function(des) {
  if (!length(des$colnames$mod)) return(NULL)
  raw <- des$site_covs_raw[[des$modifier]]
  if (is.numeric(raw)) {
    q <- stats::quantile(raw, c(0.1, 0.5, 0.9), na.rm = TRUE)
    sc <- des$scaling[[des$modifier]]
    z <- if (is.null(sc)) q else (q - sc$center) / sc$scale
    list(label = sprintf("%s (%s pct.)", signif(q, 3), c("10th", "50th", "90th")), mod = lapply(z, function(v) v), raw = q)
  } else {
    lv <- levels(factor(raw))
    X <- stats::model.matrix(stats::as.formula(paste("~", des$modifier)),
                             stats::setNames(data.frame(factor(lv, levels = lv)), des$modifier))[, -1, drop = FALSE]
    list(label = lv, mod = lapply(seq_along(lv), function(k) X[k, ]), raw = lv)
  }
}

abi_newrow <- function(des, f, var, vals) {
  sc <- des$site_covs_std
  base <- lapply(sc, function(v) if (is.numeric(v)) rep(0, length(vals)) else factor(rep(names(which.max(table(v))), length(vals)), levels = levels(v)))
  base <- as.data.frame(base); base[[var]] <- vals
  X <- stats::model.matrix(stats::as.formula(f), base)
  X
}

# ------------------------------------------------------------- interpret
abi_interpret <- function(fit, data) {
  est <- fit$estimates; o <- fit$object; des <- o$des; dg <- fit$diagnostics
  A <- des$species[["dominant"]]; S <- des$species[["subordinate"]]
  out <- character()
  if (!isTRUE(dg$converged)) out <- c(out, "WARNING: the model did not converge; the numbers below should not be interpreted.")
  lam <- abi_row(est, "dominant_mean_abundance"); tot <- abi_row(est, "dominant_total_abundance")
  occA <- abi_row(est, "dominant_prop_occupied")
  out <- c(out, sprintf(paste(
    "%s abundance: at average covariate values a site is used by %s individuals on average (Royle-Nichols abundance: the number of",
    "individuals exposed to the camera, not a density). %s"),
    A, abi_ci(lam), if (nrow(tot)) sprintf("Summed over the %d surveyed sites this is %s individual-site associations (an individual can use several sites), and %s of sites are used by at least one %s.",
                                             des$M, abi_ci(tot, 0), abi_ci(occA, 2, function(x) x), A) else ""))
  for (k in grep("^lambda_beta_", est$param, value = TRUE)) {
    r <- abi_row(est, k); cov <- sub("^lambda_beta_", "", k); sc <- des$scaling[[cov]]
    unit <- if (!is.null(sc)) sprintf(" (1 SD = %s units)", signif(sc$scale, 3)) else ""
    out <- c(out, sprintf("Each 1-SD increase in %s%s is associated with %s times the expected %s abundance (a %+.0f%% change)%s.",
                          cov, unit, abi_ci(r, 2, exp), A, 100 * (exp(r$estimate) - 1), abi_per10(r, sc, "rate ratio")))
  }
  g <- abi_row(est, "psi_beta_dominant_abundance")
  dr <- o$draws
  pneg <- mean(dr[, "gamma0"] < 0)
  hab_txt <- abi_hab_txt(des)
  out <- c(out, sprintf(paste(
    "Interaction%s: each additional %s individual at a site is associated with %s times the odds that %s occupies the site",
    "(odds ratio per individual); on the logit scale the coefficient is %s per individual. The posterior probability that the association is negative is %s."),
    hab_txt, A, abi_ci(g, 2, exp), S, abi_ci(g), abi_prob(pneg)))
  cur <- abi_psi_curve(des, dr, 0:3)
  cq <- function(k) sprintf("%.2f (%.2f-%.2f)", stats::median(cur[, k]), stats::quantile(cur[, k], 0.025), stats::quantile(cur[, k], 0.975))
  out <- c(out, sprintf("Predicted %s occupancy%s: %s at a site with no %s, %s with 1, %s with 2 and %s with 3 individuals (posterior medians, 95%% CI).",
                        S, hab_txt, cq(1), A, cq(2), cq(3), cq(4)))
  if (length(des$colnames$mod)) {
    for (k in grep("^psi_beta_dominant_abundance_x_", est$param, value = TRUE)) {
      r <- abi_row(est, k)
      sc <- des$scaling[[des$modifier]]
      dd <- abi_dummy_desc(des, sub("^psi_beta_dominant_abundance_x_", "", k))
      out <- c(out, sprintf("Habitat dependence: the interaction changes by %s logit units per individual %s; %s.",
                            abi_ci(r), if (!is.null(sc)) sprintf("for each 1-SD (%s units) increase in %s", signif(sc$scale, 3), des$modifier)
                                       else if (!is.null(dd)) paste("where", dd) else paste("per level change in", des$modifier),
                            if (r$lower > 0 || r$upper < 0) "the 95% interval excludes zero, so the data support an interaction whose strength depends on habitat"
                            else "the 95% interval includes zero, so the data do not show clearly that the interaction depends on habitat"))
    }
    sp <- abi_slope_points(des)
    for (k in seq_len(nrow(sp))) {
      r <- abi_row(est, sp$label[k])
      if (!nrow(r)) next
      dirn <- if (r$upper < 0) "negative" else if (r$lower > 0) "positive" else "not clearly different from zero"
      out <- c(out, sprintf("Where %s: odds ratio per %s individual %s (%s).", sp$desc[k], A, abi_ci(r, 2, exp), dirn))
    }
  }
  for (k in setdiff(grep("^psi_beta_", est$param, value = TRUE), grep("dominant_abundance", est$param, value = TRUE))) {
    r <- abi_row(est, k); cov <- sub("^psi_beta_", "", k)
    sc <- des$scaling[[cov]]
    out <- c(out, sprintf("%s occupancy (holding %s abundance fixed): odds ratio %s for %s%s%s.", S, A, abi_ci(r, 2, exp),
                          if (!is.null(sc)) sprintf("each 1-SD (%s units) increase in %s", signif(sc$scale, 3), cov)
                          else { dd <- abi_dummy_desc(des, cov); if (is.null(dd)) paste("unit change in", cov) else dd },
                          "", abi_per10(r, sc, "odds ratio")))
  }
  rr <- abi_row(est, "dominant_detection_per_individual"); pp <- abi_row(est, "subordinate_detection")
  unit <- if (!is.null(data$occasion_days) || any(des$K > 1)) "camera-day" else "survey occasion"
  out <- c(out, sprintf("Detection (a detection process, not ecology): each %s individual is detected with probability %s per %s; an occupied site yields a %s detection with probability %s per %s.",
                        A, abi_ci(rr, 3), unit, S, abi_ci(pp, 3), unit))
  if (!is.null(o$twostep)) {
    ts <- abi_row(o$twostep, "psi_beta_dominant_abundance")
    if (nrow(ts)) out <- c(out, sprintf(paste(
      "Comparison with the common two-step shortcut (estimate %s abundance first, then use it as a fixed covariate): the shortcut gives %s,",
      "with an interval %.0f%% as wide as the joint model's. The joint model carries the uncertainty in each site's abundance into the",
      "interaction; the shortcut treats estimated abundance as if it were known, and in simulations it is biased towards zero with too-narrow intervals."),
      A, abi_ci(ts), 100 * (ts$upper - ts$lower) / (g$upper - g$lower)))
  }
  if (!is.null(o$gof)) out <- c(out, sprintf("Goodness of fit: posterior predictive p = %.2f (%s) and %.2f (%s); %s",
                                             o$gof$p[["dominant"]], A, o$gof$p[["subordinate"]], S,
                                             if (all(o$gof$p > 0.05 & o$gof$p < 0.95)) "no evidence of lack of fit." else "LACK OF FIT detected: treat the estimates with caution (unmodelled heterogeneity, overdispersed abundance or non-independent detections are likely causes)."))
  out <- c(out, paste(
    "Caution: this is an association between the dominant's estimated abundance and the subordinate's occupancy, not proof of predation",
    "or competition. Unmodelled habitat preferences shared (or opposed) by the two species produce the same signal, so include the main",
    "habitat drivers of both species. Royle-Nichols abundance assumes individuals are detected independently and abundance is Poisson",
    "among sites; absolute abundance is less reliable than relative differences. Simulation checks for this module (VALIDATION.md)",
    "recovered the interaction without bias at 600 sites (4 weekly occasions of daily records), including when the dominant's per-individual daily detection was as low as 0.03;",
    "in an independent check with 700 sites and a weaker interaction its size was overstated by about 12%, and with 150-300 sites by roughly 10-20%",
    "(credible intervals kept near-nominal coverage). With fewer than about 1000 sites, read the interval rather than the point estimate."))
  out
}

# ------------------------------------------------------------ methods text
# render a one-sided formula as a linear predictor, e.g. "~ a + b" -> "beta_lambda0 + beta_lambda1 a + beta_lambda2 b"
fx <- function(f, nm) {
  tl <- attr(stats::terms(stats::as.formula(f)), "term.labels")
  paste(c(paste0(nm, "0"), if (length(tl)) paste0(nm, "_", seq_along(tl), " x ", tl)), collapse = " + ")
}

abi_methods_text <- function(fit, data, settings) {
  o <- fit$object; des <- o$des; s <- o$settings
  A <- des$species[["dominant"]]; S <- des$species[["subordinate"]]
  eng <- switch(o$engine,
    nimble = sprintf("The model was fitted in NIMBLE (de Valpine et al. 2017) with %d chains of %d iterations (%d burn-in, thinning %d); regression coefficients were updated with automated factor slice samplers (abundance and detection jointly; occupancy and interaction jointly), and latent abundances and occupancy states with NIMBLE's default samplers. Convergence was assessed with rank-normalised R-hat (< 1.01) and bulk effective sample size (>= 400).",
                     s$nchains, s$niter, s$nburnin, s$thin),
    ml = sprintf("The model was fitted by maximising the exact marginal likelihood, summing latent abundance over 0-%d and occupancy over {0,1}; standard errors come from the inverse Hessian.", s$Nmax),
    twostep = "Two-step plug-in: unmarked::occuRN for the dominant, then unmarked::occu for the subordinate with the posterior mean abundance as a fixed covariate.")
  paste0(
    "We modelled the occurrence of ", S, " as a function of the latent abundance of ", A, " following the abundance-mediated interaction model of Twining et al. (2025). ",
    "Abundance of ", A, " at site i was N_i ~ Poisson(lambda_i), log(lambda_i) = ", fx(des$formulas$lambda, "lambda"), ", and its detections were modelled with the Royle-Nichols (2003) model: ",
    "the number of days with a detection in occasion j, y_ij ~ Binomial(K_ij, 1 - (1 - r_ij)^N_i), where K_ij is the number of camera-days and r_ij the per-individual daily detection probability (logit(r_ij) = ", fx(des$formulas$r, "r"), "). ",
    "Occupancy of ", S, " was z_i ~ Bernoulli(psi_i), logit(psi_i) = ", fx(des$formulas$psi, "psi"), " + N_i (gamma0",
    if (length(des$colnames$mod)) paste0(" + gamma1 x ", des$modifier) else "", "), and its detections y_ij ~ Binomial(K_ij, z_i p_ij), logit(p_ij) = ", fx(des$formulas$p, "p"), ". ",
    "N_i therefore enters as a latent covariate whose uncertainty propagates into the interaction coefficient", if (length(des$colnames$mod)) ", and gamma1 lets the strength of the interaction vary with habitat (Twining et al. 2022)" else "", ". ",
    "Continuous covariates were standardised (mean 0, SD 1). Priors: intercepts of detection and occupancy Logistic(0, 1) (uniform on the probability scale); log-abundance intercept, slopes and interaction coefficients Normal(0, SD ", s$prior_sd, "). ",
    eng, " Goodness of fit was checked with posterior predictive Freeman-Tukey statistics on site detection totals for each species. ",
    "Before fitting, the model was checked by simulation: it recovered known parameter values at comparable designs (ecotools abund-interactions VALIDATION.md).")
}

# Short display labels for coefficient names (users can override any with figures: labels:)
abi_coef_label <- function(param, des, A, opts, per_sd = TRUE) {
  vapply(param, function(p) {
    if (!is.null(opts$labels[[p]])) return(opts$labels[[p]])
    cov <- sub("^(lambda|r|psi|p)_beta_", "", p)
    if (p == "psi_beta_dominant_abundance") {
      ref <- abi_mod_ref(des)
      return(if (length(des$colnames$mod) && !is.null(ref)) sprintf("Per %s (%s)", A, ref) else if (length(des$colnames$mod)) sprintf("Per %s (average habitat)", A) else sprintf("Per %s", A))
    }
    sdtag <- function(k) { sc <- des$scaling[[k]]; if (isTRUE(per_sd) && !is.null(sc) && is.finite(sc$scale)) sprintf(" [per SD = %s]", signif(sc$scale, 3)) else "" }
    if (grepl("^psi_beta_dominant_abundance_x_", p)) { k <- sub("^psi_beta_dominant_abundance_x_", "", p); return(sprintf("Per %s x %s%s", A, fig_label(k, k, opts), sdtag(k))) }
    if (grepl("^interaction_at_(low|high)_", p)) {
      k <- sub("^interaction_at_(low|high)_", "", p); hi <- grepl("^interaction_at_high_", p)
      raw <- des$site_covs_raw[[des$modifier]]
      q <- if (is.numeric(raw)) signif(stats::quantile(raw, if (hi) 0.9 else 0.1, na.rm = TRUE), 3) else NA
      return(if (is.na(q)) sprintf("Per %s, %s %s", A, if (hi) "high" else "low", fig_label(k, k, opts))
             else sprintf("Per %s, %s = %s (%s pct.)", A, fig_label(k, k, opts), q, if (hi) "90th" else "10th"))
    }
    if (grepl("^interaction_in_", p)) return(sprintf("Per %s in %s", A, sub("^interaction_in_", "", p)))
    paste0(fig_label(cov, cov, opts), sdtag(cov))
  }, "", USE.NAMES = FALSE)
}

# ------------------------------------------------------------------ figures
abi_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  o <- fit$object; des <- o$des; est <- fit$estimates; dr <- o$draws
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; inner <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci
  A <- fig_label("dominant", des$species[["dominant"]], opts); S <- fig_label("subordinate", des$species[["subordinate"]], opts)
  cil <- sprintf("%d%%", round(100 * ci))
  ityp <- if (identical(o$engine, "nimble")) "credible" else "approximate confidence (draws from the asymptotic normal of the MLE)"
  pmN <- abi_post_N(fit)
  Nhi <- max(4, min(15, ceiling(stats::quantile(pmN, 0.99, na.rm = TRUE) + 1)))
  Nseq <- seq(0, Nhi, by = 0.1)
  out <- list()
  # 1. subordinate occupancy vs dominant abundance
  pr <- abi_band(abi_psi_curve(des, dr, Nseq), Nseq, ci, inner)
  out$subordinate_vs_dominant <- with_caption(
    fig_cond_effect(pr, sprintf("%s abundance at site (individuals)", A), sprintf("%s occupancy probability", S),
                    obs = if (isTRUE(opts$show_data)) data.frame(x = pmN) else NULL, ylim01 = TRUE),
    sprintf("Predicted probability that %s occupies a site against the number of %s individuals using it (posterior median, %s and %d%% %s bands; full posterior uncertainty of all coefficients), other covariates at their means%s. Rug: posterior mean %s abundance at each surveyed site.",
            S, A, cil, round(100 * inner), ityp,
            if (!length(des$colnames$mod)) "" else if (is.null(abi_mod_ref(des))) " and habitat at its average"
            else sprintf(" and %s at its reference level (%s)", des$modifier, abi_mod_ref(des)), A))
  # 2. interaction by habitat
  hl <- abi_hab_levels(des)
  if (!is.null(hl)) {
    pd <- do.call(rbind, lapply(seq_along(hl$label), function(k) {
      b <- abi_band(abi_psi_curve(des, dr, Nseq, hl$mod[[k]]), Nseq, ci, inner); b$group <- hl$label[k]; b
    }))
    pd$group <- factor(pd$group, levels = hl$label)
    mlab <- fig_label(des$modifier, des$modifier, opts)
    pd$group <- factor(paste(mlab, "=", pd$group), levels = paste(mlab, "=", hl$label))
    # short wide figure: wrap the y title and strips so they are not clipped at narrow widths;
    # strip wrap width follows the exported width (about 1.9 mm per bold character in each of 3 panels)
    wsel <- opts$per_figure$interaction_by_habitat$width; if (is.null(wsel)) wsel <- "double"
    wmm <- tryCatch(fig_width_mm(wsel, opts), error = function(e) NULL); if (!is.numeric(wmm) || !length(wmm)) wmm <- 170
    stripw <- max(14, floor((wmm / 3 - 6) / 1.9))
    p <- fig_cond_effect(pd, sprintf("%s per site (individuals)", A),
                         paste(strwrap(sprintf("%s occupancy probability", S), 24), collapse = "\n"), ylim01 = TRUE) +
      ggplot2::facet_wrap(~group, nrow = 1, labeller = ggplot2::label_wrap_gen(stripw)) +
      ggplot2::guides(colour = "none", fill = "none")
    attr(p, "fig_size") <- list(width = "double", aspect = 0.42)
    out$interaction_by_habitat <- with_caption(p, sprintf(
      "%s occupancy against %s abundance at the 10th, 50th and 90th percentiles of %s (posterior median, %s and %d%% %s bands). Different slopes across panels mean the association between the species depends on habitat.",
      S, A, mlab, cil, round(100 * inner), ityp))
    raw <- des$site_covs_raw[[des$modifier]]
    if (is.numeric(raw)) {
      sc <- des$scaling[[des$modifier]]
      xr <- seq(min(raw), max(raw), length.out = 80); xz <- if (is.null(sc)) xr else (xr - sc$center) / sc$scale
      sl <- dr[, "gamma0"] + outer(dr[, "gamma_mod[1]"], xz)
      b <- abi_band(sl, xr, ci, inner)
      out[[paste0("interaction_strength_", abi_clean(des$modifier))]] <- with_caption(
        fig_cond_effect(b, mlab, sprintf("Interaction (logit per %s)", A),
                        obs = if (isTRUE(opts$show_data)) data.frame(x = raw) else NULL) +
          ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = eco_pal$null),
        sprintf("Strength of the abundance-mediated interaction across %s (logit scale; below the dashed line = fewer %s where there are more %s). Band: %s %s interval; rug: surveyed sites.",
                mlab, S, A, cil, ityp))
    }
  }
  # 3. dominant abundance vs each numeric site covariate
  bl <- dr[, paste0("beta_lam[", seq_along(des$colnames$lambda), "]"), drop = FALSE]
  for (v in intersect(all.vars(stats::as.formula(des$formulas$lambda)), names(des$site_covs_std))) {
    z <- des$site_covs_std[[v]]
    if (!is.numeric(z)) next
    zs <- seq(min(z), max(z), length.out = 80)
    X <- abi_newrow(des, des$formulas$lambda, v, zs)
    b <- abi_band(exp(bl %*% t(X)), fig_x(zs, v, opts, des$scaling), ci, inner)
    out[[paste0("dominant_abundance_", abi_clean(v))]] <- with_caption(
      fig_cond_effect(b, fig_label(v, v, opts), sprintf("%s per site (individuals)", A),
                      obs = if (isTRUE(opts$show_data)) data.frame(x = fig_x(z, v, opts, des$scaling)) else NULL),
      sprintf("Expected number of %s individuals using a site against %s (posterior median with %s and %d%% %s bands), other covariates at their means. Rug: surveyed sites.",
              A, fig_label(v, v, opts), cil, round(100 * inner), ityp))
  }
  # 4. subordinate occupancy vs its covariates (no dominant vs typical high dominant abundance)
  Nq <- max(1, round(stats::quantile(pmN, 0.9, na.rm = TRUE)))
  for (v in setdiff(intersect(all.vars(stats::as.formula(des$formulas$psi)), names(des$site_covs_std)), des$modifier)) {
    z <- des$site_covs_std[[v]]
    if (!is.numeric(z)) next
    zs <- seq(min(z), max(z), length.out = 60)
    X <- abi_newrow(des, des$formulas$psi, v, zs)
    pd <- do.call(rbind, lapply(c(0, Nq), function(n) {
      m <- sapply(seq_along(zs), function(k) abi_psi_curve(des, dr, n, xrow = X[k, ])[, 1])
      b <- abi_band(m, fig_x(zs, v, opts, des$scaling), ci, inner)
      b$group <- sprintf("%d %s", n, A); b[, c("x", "estimate", "lower", "upper", "group")]
    }))
    pd$group <- factor(pd$group, levels = unique(pd$group))
    out[[paste0("subordinate_occupancy_", abi_clean(v))]] <- with_caption(
      fig_cond_effect(pd, fig_label(v, v, opts), sprintf("%s occupancy probability", S),
                      obs = if (isTRUE(opts$show_data)) data.frame(x = fig_x(z, v, opts, des$scaling)) else NULL, ylim01 = TRUE) +
        ggplot2::aes(linetype = group) + ggplot2::labs(colour = NULL, fill = NULL, linetype = NULL),

      sprintf("%s occupancy against %s at sites with no %s and with %d (90th percentile of estimated site abundance); %s %s bands; other covariates at their means.",
              S, fig_label(v, v, opts), A, Nq, cil, ityp))
    # the legend sits beside the panel; at single-column width it squeezed the panel to half the figure
    attr(out[[paste0("subordinate_occupancy_", abi_clean(v))]], "fig_size") <- list(width = "onehalf", aspect = 0.6)
  }
  # 5. map of dominant abundance
  if (!is.null(des$coords)) {
    md <- data.frame(des$coords, value = pmN)
    # projected coordinates in metres (6-digit eastings) are shown in km so the axis labels stay short
    km <- max(diff(range(md$x)), diff(range(md$y))) > 5000
    if (km) { md$x <- md$x / 1000; md$y <- md$y / 1000 }
    asp <- diff(range(md$y)) / max(diff(range(md$x)), 1e-9)
    # greyscale: core only re-maps discrete scales, so the continuous colour scale is switched here
    grey <- isTRUE(opts$greyscale)
    cscale <- if (grey) ggplot2::scale_colour_gradient(name = sprintf("%s\nper site", A), low = "#cccccc", high = "#000000")
              else ggplot2::scale_colour_viridis_c(name = sprintf("%s\nper site", A))
    out$dominant_abundance_map <- with_caption(
      ggplot2::ggplot(md, ggplot2::aes(x, y, colour = value)) + ggplot2::geom_point(size = 1.2, alpha = 0.9) +
        cscale + ggplot2::coord_equal() +
        # long projected coordinates (e.g. 6-digit eastings) collide at single-column width
        ggplot2::scale_x_continuous(n.breaks = 3, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::labs(x = fig_label("x", if (km) "Easting (km)" else "Easting", opts),
                      y = fig_label("y", if (km) "Northing (km)" else "Northing", opts)),
      sprintf("Posterior mean number of %s individuals using each surveyed site (%s).%s", A,
              if (grey) "grey scale, darker = more" else "colour, viridis scale",
              if (km) " Coordinates are the projected site coordinates divided by 1000 (assumed metres, shown as km)." else ""))
    attr(out$dominant_abundance_map, "fig_size") <- list(aspect = min(1.2, max(0.45, asp * 0.85)))
  }
  # 6. coefficients (link scale) by process
  cf <- est[grepl("_beta_", est$param) & !is.na(est$estimate), c("param", "estimate", "lower", "upper")]
  if (nrow(cf)) {
    cf$group <- ifelse(grepl("^lambda_", cf$param), paste(A, "abundance (log)"),
                ifelse(grepl("^r_", cf$param), paste(A, "detection (logit)"),
                ifelse(grepl("dominant_abundance", cf$param), "Interaction (logit)",
                ifelse(grepl("^psi_", cf$param), paste(S, "occupancy (logit)"), paste(S, "detection (logit)")))))
    gl <- c(paste(A, "abundance (log)"), paste(A, "detection (logit)"), paste(S, "occupancy (logit)"), "Interaction (logit)", paste(S, "detection (logit)"))
    cf$group <- factor(cf$group, levels = intersect(gl, cf$group))
    cf$param <- abi_coef_label(cf$param, des, A, opts)
    # rows grouped by process, first process at the top
    cf <- cf[order(cf$group), ]
    cf$param <- factor(cf$param, levels = rev(unique(cf$param)))
    pc0 <- fig_coef(cf, null_value = 0, xlab = sprintf("Coefficient (%s interval, link scale)", cil), order = FALSE)
    # long relabelled rows (habitat names, SD tags) squeezed the panel: wrap them
    wrapy <- ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(z) paste(strwrap(z, 34), collapse = "\n"), ""))
    pc <- pc0 + ggplot2::aes(shape = group) + ggplot2::labs(shape = NULL) + wrapy  # shape = second encoding for greyscale
    attr(pc, "fig_size") <- attr(pc0, "fig_size")
    attr(pc, "fig_size")$width <- "double"
    out$coefficients <- with_caption(pc,
                                     sprintf("Regression coefficients on the link scale (log for abundance, logit for occupancy and detection) with %s %s intervals; dashed line = no effect. Continuous covariates were standardised, so each slope is the change per 1 SD of that covariate (SD in original units in square brackets); factor levels are contrasts with the reference level. Interaction terms: change in %s log-odds of occupancy per %s individual.",
                                             cil, if (o$engine == "nimble") "credible" else "confidence", S, A))
  }
  # 7. joint vs two-step
  if (!is.null(o$twostep)) {
    keys <- c("psi_beta_dominant_abundance", grep("^psi_beta_dominant_abundance_x_|^interaction_", est$param, value = TRUE))
    a <- est[est$param %in% keys, c("param", "estimate", "lower", "upper")]; a$group <- "Joint model (this analysis)"
    b <- o$twostep[o$twostep$param %in% keys, c("param", "estimate", "lower", "upper")]; b$group <- "Two-step plug-in"
    ab <- rbind(a, b); ab$param <- abi_coef_label(ab$param, des, A, opts)
    ab$group <- factor(ab$group, levels = c("Joint model (this analysis)", "Two-step plug-in"))  # dodge draws level 1 lowest; legend reversed below to match
    ab$param <- factor(ab$param, levels = rev(unique(abi_coef_label(a$param, des, A, opts))))  # average, x, 10th, 90th top-down
    pj0 <- fig_coef(ab, null_value = 0, xlab = sprintf("Interaction coefficient (logit scale, %s interval)", cil), order = FALSE)
    pj <- pj0 + ggplot2::aes(shape = group) + ggplot2::labs(shape = NULL) +
      ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(z) paste(strwrap(z, 42), collapse = "\n"), "")) +
      ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE), shape = ggplot2::guide_legend(reverse = TRUE))
    attr(pj, "fig_size") <- attr(pj0, "fig_size")
    attr(pj, "fig_size")$width <- "double"
    out$joint_vs_twostep <- with_caption(pj,
                                         sprintf("Interaction estimates from the joint model and from the two-step shortcut that plugs in estimated %s abundance as if known (points: estimates; bars: %s intervals). The joint model propagates abundance uncertainty.%s Interaction slopes are per %s individual on the logit scale of %s occupancy.",
                                                 A, cil, if (length(des$colnames$mod) && is.numeric(des$site_covs_raw[[des$modifier]])) " Rows at the 10th and 90th percentiles give the slope at those habitat values; the 'x' row is the change in slope per 1 SD of the habitat variable." else "", A, S))
  }
  # 8. goodness of fit
  if (!is.null(o$gof)) {
    st <- o$gof$stats
    th <- theme_eco(if (is.null(opts$base_size)) 9 else opts$base_size, opts$font, opts$style)
    pp <- function(d, lab, pv) fig_ppc(d, 0, sprintf("%s: T(rep) - T(obs)\np = %.2f", lab, pv)) + ggplot2::labs(caption = NULL) + th
    pg <- compose_figures(list(pp(st[, "repA"] - st[, "obsA"], A, o$gof$p[["dominant"]]),
                               pp(st[, "repS"] - st[, "obsS"], S, o$gof$p[["subordinate"]])), ncol = 2)
    attr(pg, "fig_size") <- list(width = "double", aspect = 0.4)
    out$gof <- with_caption(pg,
      "Posterior predictive check: distribution of the Freeman-Tukey discrepancy of replicated minus observed site detection totals; values centred on 0 (line) indicate adequate fit.")
  }
  # 9. traces (Bayesian)
  if (!is.null(o$chains)) {
    kp <- c("beta_lam[1]", "beta_r[1]", "beta_psi[1]", "gamma0", if (length(des$colnames$mod)) "gamma_mod[1]", "beta_p[1]")
    kp <- intersect(kp, colnames(o$chains[[1]]))
    mlab0 <- if (length(des$colnames$mod)) fig_label(des$modifier, des$modifier, opts) else ""
    nice <- c("beta_lam[1]" = sprintf("%s log abundance, intercept", A),
              "beta_r[1]" = sprintf("%s detection (logit), intercept", A),
              "beta_psi[1]" = sprintf("%s occupancy (logit), intercept", S),
              "gamma0" = sprintf("Interaction per %s (logit)", A),
              "gamma_mod[1]" = sprintf("Interaction x %s (per SD)", mlab0),
              "beta_p[1]" = sprintf("%s detection (logit), intercept", S))
    nm <- vapply(kp, function(k) fig_label(k, nice[[k]], opts), "")
    ch <- lapply(o$chains, function(m) { m <- as.matrix(m)[, kp, drop = FALSE]; colnames(m) <- nm; m })
    pt <- fig_trace(ch) + ggplot2::facet_wrap(~ind, scales = "free_y", ncol = 2, labeller = ggplot2::label_wrap_gen(30)) +
      ggplot2::scale_x_continuous(n.breaks = 3) + ggplot2::labs(x = "Iteration (after burn-in and thinning)")
    attr(pt, "fig_size") <- list(width = "double", aspect = 0.6)
    out$trace <- with_caption(pt, "MCMC traces for the intercepts and interaction coefficients (link scale; one colour per chain). Well-mixed chains overlap with no trend.")
    pall <- do.call(rbind, ch)
    pp <- fig_posterior(pall, ncol = 2) + ggplot2::facet_wrap(~param, scales = "free", ncol = 2, labeller = ggplot2::label_wrap_gen(30)) +
      ggplot2::geom_vline(data = data.frame(param = factor(nm[grepl("^gamma", kp)], levels = levels(utils::stack(as.data.frame(pall))$ind)), z = 0),
                          ggplot2::aes(xintercept = z), linetype = 2, colour = eco_pal$null, linewidth = 0.3) +
      ggplot2::labs(x = "Value (link scale)")
    attr(pp, "fig_size") <- list(width = "double", aspect = 0.6)
    out$posterior <- with_caption(pp, "Posterior densities (all chains pooled) of the intercepts and interaction coefficients on the link scale; dashed line = 0 (no interaction), drawn on the interaction panels only.")
  }
  out
}

abi_plots <- function(fit, data) {
  f <- abi_figures(fit, data, fig_opts())
  lapply(f, function(p) if (inherits(p, "ggplot")) p + theme_eco(11) else p)
}
