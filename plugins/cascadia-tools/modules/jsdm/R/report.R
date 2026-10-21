# jsdm: interpretation, methods text and the publication figure suite.

jsdm_est <- function(fit, p) fit$estimates[match(p, fit$estimates$param), ]

jsdm_support <- function(fit, cols) colMeans(fit$draws[, cols, drop = FALSE] > 0)

# Names of the interpreted response / community / trait parameters. For
# presence-absence (probit) these are the population-averaged quantities
# (latent site factors integrated out), which are identified; Hmsc's raw Beta
# and Gamma are conditional on the latent factors and their scale depends on
# the latent-factor variance (VALIDATION.md).
jsdm_presence <- function(fit) fit$prep$response == "presence"
jsdm_resp_names <- function(fit, k, sp = fit$prep$sp)
  paste0(if (jsdm_presence(fit)) "popavg_" else "beta_", k, "_", sp)
jsdm_comm_name <- function(fit, k) paste0(if (jsdm_presence(fit)) "comm_popavg_" else "comm_mean_", k)
jsdm_trait_re <- function(fit) if (jsdm_presence(fit)) "^popavg_trait_(.*)_on_(.*)$" else "^trait_(.*)_on_(.*)$"

# symmetric species x species matrix from pair columns
jsdm_pair_matrix <- function(vals, sp) {
  ns <- length(sp)
  M <- matrix(NA_real_, ns, ns, dimnames = list(sp, sp))
  M[upper.tri(M)] <- vals
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

strwrap_lab <- function(x, width) vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), "", USE.NAMES = FALSE)

jsdm_cov_label <- function(k, fit, opts = NULL) {
  sc <- fit$prep$X$scaling[[k]]
  fig_label(k, if (!is.null(sc)) sc$label else k, opts)
}

jsdm_sp_label <- function(sp, fit, opts = NULL)
  vapply(sp, function(s) fig_label(s, as.character(fit$prep$sp_labels[[s]] %||% s), opts), "")

# one row per model term that is a numeric covariate (for gradients)
jsdm_gradient_covs <- function(fit, data) {
  vars <- all.vars(stats::as.formula(fit$prep$x_formula))
  vars[vapply(vars, function(v) is.numeric(data$site_data[[v]]), logical(1))]
}

interpret_jsdm <- function(fit, data) {
  pr <- fit$prep; sp <- pr$sp; ns <- length(sp); n <- nrow(pr$Y)
  presence <- pr$response == "presence"
  out <- character()
  dg <- fit$diagnostics
  if (!isTRUE(dg$converged) || dg$rhat_max >= 1.01 || dg$ess_min < 400)
    out <- c(out, sprintf("CAUTION: the MCMC has not reached publication-grade convergence (max R-hat %.3f, min bulk ESS %.0f; the targets are < 1.01 and >= 400). Treat the numbers below as provisional and rerun with more samples/thinning.",
                          dg$rhat_max, dg$ess_min))
  prev <- colMeans(pr$Y > 0)
  out <- c(out, sprintf("The model describes %d species at %d sampling units (%s). Species were recorded at %.0f%% to %.0f%% of sites (median %.0f%%).",
                        ns, n, if (presence) "presence/absence, probit link" else "counts, lognormal Poisson",
                        100 * min(prev), 100 * max(prev), 100 * stats::median(prev)))
  out <- c(out, sprintf("All responses below describe where species were %s, not true occurrence: HMSC assumes detection is perfect or constant within a species. A covariate that also changes detectability (for example canopy cover or noise in point counts, or observer and effort) shows up here as a 'response', so ecological and detection effects cannot be separated with this model (use a multi-species occupancy model with repeat visits for that).",
                        if (presence) "recorded" else "counted"))
  covs <- setdiff(pr$X$names, "intercept")
  bins <- vapply(covs, function(k) isTRUE(pr$X$scaling[[k]]$binary), logical(1))
  avg_txt <- if (any(bins)) sprintf("at average conditions, i.e. continuous covariates at their mean and %s = 0", paste(covs[bins], collapse = ", ")) else "at average conditions"
  p0 <- if (presence) mean(jsdm_est(fit, paste0("prev_", sp))$estimate) else NA
  for (k in covs) {
    bn <- jsdm_resp_names(fit, k)
    sup <- jsdm_support(fit, bn)
    npos <- sum(sup >= 0.95); nneg <- sum(sup <= 0.05)
    g <- jsdm_est(fit, jsdm_comm_name(fit, k))
    sc <- pr$X$scaling[[k]]
    lab <- data$labels[[sc$label]] %||% data$labels[[k]] %||% sc$label
    unit <- if (isTRUE(sc$binary)) sprintf("'%s' (0/1)", lab) else sprintf("1 SD of %s (1 SD = %.3g in the units of %s)", lab, sc$scale, sc$label)
    resp_txt <- if (presence)
      sprintf(" For a species with the average occurrence probability (%.2f %s), that is a change to %.2f (95%% CrI %.2f to %.2f) per %s.",
              p0, avg_txt, stats::pnorm(stats::qnorm(p0) + g$estimate), stats::pnorm(stats::qnorm(p0) + g$lower),
              stats::pnorm(stats::qnorm(p0) + g$upper), if (isTRUE(sc$binary)) "change from 0 to 1" else "SD increase")
    else sprintf(" That is a %+.0f%% change in expected count (95%% CrI %+.0f%% to %+.0f%%).",
                 100 * (exp(g$estimate) - 1), 100 * (exp(g$lower) - 1), 100 * (exp(g$upper) - 1))
    gm <- jsdm_est(fit, paste0("comm_mean_", k))
    pool_txt <- if (presence && !is.na(gm$estimate))
      sprintf(" This is the average over these %d species; for the wider species pool they represent, the hierarchical community mean (Gamma; on Hmsc's conditional probit scale, which is larger than the population-averaged scale) is %.2f (95%% CrI %.2f to %.2f), %s.",
              ns, gm$estimate, gm$lower, gm$upper,
              if (gm$lower > 0 || gm$upper < 0) "which also excludes zero" else "which does not exclude zero: a new species from the same pool could respond either way")
    else ""
    out <- c(out, sprintf("%s: %d of %d species respond positively and %d negatively with at least 95%% posterior support. The average response across these species is %.2f per %s on the %s scale (95%% CrI %.2f to %.2f).%s%s",
                          lab, npos, ns, nneg, g$estimate, unit,
                          if (presence) "probit (population-averaged)" else "log", g$lower, g$upper, resp_txt, pool_txt))
  }
  gc <- jsdm_gradient_covs(fit, data)
  for (v in gc) {
    q <- stats::quantile(data$site_data[[v]], c(0.1, 0.9), na.rm = TRUE)
    gr <- tryCatch(jsdm_gradient(fit, data, v, xs = as.numeric(q)), error = function(e) NULL)
    if (is.null(gr)) next
    b <- summ_band(gr$community)
    what <- if (presence) "Expected species richness (number of these species present)" else "Expected total abundance (summed over species)"
    out <- c(out, sprintf("%s is %.1f (95%% CrI %.1f-%.1f) at %s = %.3g (10th percentile of sites) and %.1f (%.1f-%.1f) at %.3g (90th percentile), other covariates at their mean. This is an association, not evidence that %s causes the change.",
                          what, b$estimate[1], b$lower[1], b$upper[1], data$labels[[v]] %||% v, q[1], b$estimate[2], b$lower[2], b$upper[2], q[2], data$labels[[v]] %||% v))
  }
  trn <- grep(jsdm_trait_re(fit), fit$estimates$param, value = TRUE)
  if (length(trn)) {
    tt <- fit$estimates[fit$estimates$param %in% trn, ]
    sup <- jsdm_support(fit, tt$param)
    strong <- tt[sup >= 0.95 | sup <= 0.05, ]
    if (nrow(strong)) {
      txt <- vapply(seq_len(nrow(strong)), function(i) {
        p <- strong$param[i]; parts <- regmatches(p, regexec(jsdm_trait_re(fit), p))[[1]]
        tsc <- pr$Tr$scaling[[parts[2]]]
        tl <- if (!is.null(tsc)) tsc$label else parts[2]
        tl <- data$labels[[parts[2]]] %||% tl
        clab <- function(k) data$labels[[pr$X$scaling[[k]]$label]] %||% pr$X$scaling[[k]]$label
        what <- if (parts[3] == "intercept")
          paste(if (strong$estimate[i] > 0) "higher" else "lower", if (presence) "baseline prevalence" else "baseline abundance")
        else paste(if (strong$estimate[i] > 0) "more positive" else "more negative", "response to", clab(parts[3]))
        # pool-level check (Hmsc Gamma): the population-averaged slope above
        # describes these particular species and is "supported" in roughly half of
        # simulated datasets whose true trait effect (Gamma) is 0 (VALIDATION.md,
        # Independent verification); Gamma carries the generalisable test.
        gm <- if (presence) jsdm_est(fit, paste0("trait_", parts[2], "_on_", parts[3])) else NULL
        pool <- if (!is.null(gm) && !is.na(gm$estimate))
          sprintf("; pool-level trait effect Gamma %.2f, 95%% CrI %.2f to %.2f, which %s", gm$estimate, gm$lower, gm$upper,
                  if (gm$lower > 0 || gm$upper < 0) "also excludes zero" else "does not exclude zero, so this pattern may not hold beyond these species")
          else ""
        sprintf("species with higher %s have a %s (%.2f per 1 SD of the trait on the %s scale, 95%% CrI %.2f to %.2f%s)",
                tl, what, strong$estimate[i], if (presence) "population-averaged probit" else "log", strong$lower[i], strong$upper[i], pool)
      }, "")
      out <- c(out, paste0("Traits (across the modelled species): ", paste(txt, collapse = "; "), ".",
                           if (presence) " The first number describes the species modelled here; the Gamma interval is the test of whether the trait-response relationship extends to the wider species pool (when Gamma is truly 0 it excludes zero in about 5% of simulated datasets, while the species-level description does so in about half)." else ""))
    } else out <- c(out, "Traits: no trait effect on species' responses had at least 95% posterior support.")
    if (!is.null(fit$vp$R2T$Beta)) {
      r2 <- fit$vp$R2T$Beta
      r2 <- r2[names(r2) != "intercept" & !grepl("Intercept", names(r2))]
      if (length(r2)) out <- c(out, sprintf("Traits explain %s of the among-species variation in responses (R2 of Beta on traits).",
                                            paste(sprintf("%.0f%% (%s)", 100 * r2, vapply(names(r2), function(k) data$labels[[k]] %||% k, "")), collapse = ", ")))
    }
  }
  if (!is.null(fit$vp)) {
    m <- rowMeans(fit$vp$vals)
    out <- c(out, sprintf("Variance partitioning (average over species, of the explained variance on the linear-predictor scale): %s.",
                          paste(sprintf("%s %.0f%%", vapply(names(m), function(k) data$labels[[k]] %||% k, ""), 100 * m), collapse = "; ")))
  }
  pn <- paste0("latcor_", pair_names(sp))
  sup <- jsdm_support(fit, pn)
  est <- jsdm_est(fit, pn)
  npos <- sum(sup >= 0.95); nneg <- sum(sup <= 0.05)
  o <- order(-abs(est$estimate))[seq_len(min(3, length(pn)))]
  o <- o[sup[o] >= 0.95 | sup[o] <= 0.05]
  pretty_pair <- function(p) { a <- strsplit(sub("^latcor_", "", p), "__")[[1]]; paste(jsdm_sp_label(a, fit), collapse = " & ") }
  out <- c(out, sprintf("Residual associations: of %d species pairs, %d co-occur more and %d less often than expected from the covariates, with at least 95%% posterior support.%s These residual correlations (latent-scale) are what remains after the covariates; they can arise from unmeasured habitat variables, shared history or sampling, as well as from biotic interactions, so they are not evidence of interactions (Blanchet et al. 2020).",
                        length(pn), npos, nneg,
                        if (length(o)) paste0(" Strongest: ", paste(sprintf("%s (r = %.2f, 95%% CrI %.2f to %.2f)",
                                                                          vapply(est$param[o], pretty_pair, ""), est$estimate[o], est$lower[o], est$upper[o]), collapse = "; "), ".") else ""))
  if ("phylo_signal_rho" %in% fit$estimates$param) {
    r <- jsdm_est(fit, "phylo_signal_rho")
    out <- c(out, sprintf("Phylogenetic signal in species' responses (rho, 0 = none, 1 = fully phylogenetic, beyond what traits explain): %.2f (95%% CrI %.2f to %.2f).",
                          r$estimate, r$lower, r$upper))
  }
  if (!is.null(fit$model_fit)) {
    mf <- fit$model_fit
    if (presence && "TjurR2" %in% names(mf))
      out <- c(out, sprintf("Explanatory power (in-sample, including the latent site factors): mean Tjur R2 %.2f (range %.2f-%.2f); mean AUC %.2f. In-sample fit is optimistic; use cross-validation for predictive performance.",
                            mean(mf$TjurR2, na.rm = TRUE), min(mf$TjurR2, na.rm = TRUE), max(mf$TjurR2, na.rm = TRUE), mean(mf$AUC, na.rm = TRUE)))
    if (!presence && "SR2" %in% names(mf))
      out <- c(out, sprintf("Explanatory power (in-sample): mean pseudo-R2 (squared Spearman correlation of observed and predicted counts) %.2f.",
                            mean(mf$SR2, na.rm = TRUE)))
  }
  if (!is.null(fit$ppc)) {
    pp <- fit$ppc
    out <- c(out, sprintf("Posterior predictive check (new sites): the observed variance of site richness (%.2f) lies at P = %.2f of the replicated distribution, and the mean pairwise co-occurrence at P = %.2f; values near 0 or 1 (< 0.05 or > 0.95) indicate the model misses community structure.",
                          pp$observed[1], pp$p_value[1], pp$p_value[2]))
  }
  if (length(data$read_notes)) out <- c(out, paste("Data notes:", paste(data$read_notes, collapse = " ")))
  out
}

methods_jsdm <- function(fit, data, settings) {
  pr <- fit$prep; st <- fit$settings
  presence <- pr$response == "presence"
  paste0(
    "We fitted a joint species distribution model with the Hierarchical Modelling of Species Communities framework (HMSC; Ovaskainen et al. 2017; Ovaskainen & Abrego 2020) using the R package Hmsc ",
    as.character(utils::packageVersion("Hmsc")), " (Tikhonov et al. 2020). ",
    if (presence) "Presence/absence of each species at each sampling unit was modelled with a probit link. "
    else "Counts of each species were modelled with a lognormal Poisson distribution (Poisson with a species-specific normally distributed overdispersion term on the log scale). ",
    sprintf("The linear predictor included the site covariates %s (continuous covariates standardised to mean 0 and SD 1; coefficients are per SD). ",
            paste(vapply(setdiff(pr$X$names, "intercept"), function(k) pr$X$scaling[[k]]$label, ""), collapse = ", ")),
    if (!is.null(pr$Tr)) sprintf("Species' responses (Beta) were modelled hierarchically as a function of species traits (%s; standardised), so the trait effects (Gamma) describe how traits shape responses to the environment. ",
                                  paste(vapply(setdiff(pr$Tr$names, "intercept"), function(k) pr$Tr$scaling[[k]]$label, ""), collapse = ", "))
    else "Species' responses were drawn from a common community distribution (no traits). ",
    if (!is.null(pr$C)) "Residual variation in responses was allowed a phylogenetic correlation structure, with phylogenetic signal parameter rho estimated (Hmsc default prior: 50% mass at rho = 0, the rest uniform on a grid over (0, 1]). " else "",
    sprintf("Residual co-occurrence was modelled with %s latent factors at the level of the sampling unit (multiplicative gamma shrinkage prior of Bhattacharya & Dunson 2011; number of factors adapted during burn-in, maximum %d)%s. ",
            fit$extra_nf %||% "site-level", st$max_factors %||% 5,
            if (identical(fit$spatial, "none")) "" else paste0(", with spatially structured factors (", fit$spatial, " Gaussian process)")),
    "Residual species-to-species association matrices were derived from the factor loadings: Hmsc's conventional correlation of the latent-factor covariance Omega = Lambda'Lambda, and the latent-scale residual correlation that also includes the link-level residual variance (cov2cor(Omega + I) for probit), which we used for inference; associations were reported when the posterior probability of a positive (or negative) sign was at least 0.95. ",
    if (presence) "With a probit link, Hmsc's species responses (Beta) and trait effects (Gamma) are conditional on the latent site factors, and their scale depends on each species' latent-factor variance Omega_jj, which binary data identify only weakly; we therefore report and interpret species responses on the population-averaged probit scale, beta_kj / sqrt(1 + Omega_jj) (latent factors integrated out), and trait effects as the least-squares regression of these population-averaged responses on the standardised traits, computed for every posterior draw. " else "",
    sprintf("Priors were the Hmsc defaults except for the latent-factor shrinkage prior (a1 = %s, b1 = %s, a2 = %s, b2 = %s; Bhattacharya & Dunson 2011), which shrinks residual associations less than the Hmsc default (a1 = a2 = 50). ",
            st$factor_prior$a1 %||% 50, st$factor_prior$b1 %||% 1, st$factor_prior$a2 %||% 50, st$factor_prior$b2 %||% 1),
    sprintf("We ran %d MCMC chains of %d samples each, thinning %d, after a transient (burn-in) of %d iterations, and assessed convergence with the rank-normalised R-hat and bulk effective sample size (Vehtari et al. 2021) of all community-level parameters, species responses, latent-scale residual correlations and derived prevalences. ",
            st$n_chains %||% 4, st$samples %||% 250, st$thin %||% 10, st$transient %||% 2500),
    "Variance in the linear predictor was partitioned among covariates and the latent factors (computeVariancePartitioning), explanatory power was summarised by Tjur's R2 (Tjur 2009) and AUC (presence/absence) or pseudo-R2 (counts), and fit was checked with a posterior predictive check on the variance of site richness and pairwise co-occurrence. ",
    "The module was validated by simulation (recovery of Beta, Gamma and residual correlations from data simulated from the HMSC model; see the module VALIDATION.md).")
}

# ---- figures -----------------------------------------------------------------

figures_jsdm <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  pr <- fit$prep; sp <- pr$sp; ns <- length(sp)
  presence <- pr$response == "presence"
  ci <- opts$ci %||% 0.95; inner <- opts$inner_ci %||% 0.5
  splab <- jsdm_sp_label(sp, fit, opts)
  covs <- setdiff(pr$X$names, "intercept")
  figs <- list()

  # species_responses: Beta heatmap, cells without 95% support blank
  if (length(covs)) {
    bn <- outer(sp, covs, function(s, k) jsdm_resp_names(fit, k, s))
    B <- matrix(jsdm_est(fit, as.vector(bn))$estimate, ns, length(covs))
    S <- matrix(jsdm_support(fit, as.vector(bn)), ns, length(covs))
    Bm <- B; Bm[S < 0.95 & S > 0.05] <- NA
    dimnames(Bm) <- list(splab, vapply(covs, jsdm_cov_label, "", fit = fit, opts = opts))
    lim <- max(abs(B), na.rm = TRUE)
    p <- fig_matrix(Bm, fig_label("beta_fill", "Response\n(per SD)", opts), diverging = TRUE, label = ns <= 30) +
      ggplot2::scale_fill_gradient2(low = eco_pal$diverging[["low"]], mid = eco_pal$diverging[["mid"]],
                                    high = eco_pal$diverging[["high"]], midpoint = 0, limits = c(-lim, lim),
                                    na.value = "#eeeeee", name = fig_label("beta_fill", "Response\n(per SD)", opts)) +
      ggplot2::coord_cartesian() + ggplot2::theme(aspect.ratio = NULL) +
      ggplot2::scale_x_discrete(labels = function(x) vapply(strwrap_lab(x, 14), identity, ""))
    attr(p, "fig_size") <- list(height = 30 + 4.2 * ns, width = if (length(covs) > 4) "double" else "onehalf")
    figs$species_responses <- with_caption(p, sprintf(
      "Species-specific responses (%s, per 1 SD of each covariate). Coloured cells have at least 95%% posterior support for a positive (red) or negative (blue) response; grey cells lack that support. Numbers are posterior means.",
      if (presence) "probit scale, population-averaged over the latent site factors" else "Beta, log scale"))
  }

  # residual_associations: latent-scale residual correlations with support mask
  pn <- paste0("latcor_", pair_names(sp))
  R <- jsdm_pair_matrix(jsdm_est(fit, pn)$estimate, sp)
  Sx <- jsdm_pair_matrix(jsdm_support(fit, pn), sp)
  Rm <- R; Rm[Sx < 0.95 & Sx > 0.05] <- NA
  dimnames(Rm) <- list(splab, splab)
  p <- fig_matrix(Rm, fig_label("assoc_fill", "Residual\ncorrelation", opts), diverging = TRUE, label = ns <= 15) +
    ggplot2::scale_fill_gradient2(low = eco_pal$diverging[["low"]], mid = eco_pal$diverging[["mid"]],
                                  high = eco_pal$diverging[["high"]], midpoint = 0, limits = c(-1, 1),
                                  na.value = "#eeeeee", name = fig_label("assoc_fill", "Residual\ncorrelation", opts)) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))
  attr(p, "fig_size") <- list(width = if (ns > 12) "double" else "onehalf", aspect = 0.9)
  figs$residual_associations <- with_caption(p,
    "Residual species associations: posterior mean latent-scale residual correlations (after the covariates). Coloured cells have at least 95% posterior support for a positive (red) or negative (blue) association; grey = not supported, blank diagonal. Associations are not evidence of biotic interactions (Blanchet et al. 2020): unmeasured shared habitat preferences produce the same pattern.")

  # variance_partitioning
  if (!is.null(fit$vp)) {
    V <- fit$vp$vals
    d <- data.frame(species = rep(splab, each = nrow(V)), component = rep(rownames(V), ns), value = as.vector(V))
    comp_lab <- vapply(rownames(V), function(r) fig_label(r, if (r %in% names(pr$X$scaling)) pr$X$scaling[[r]]$label else r, opts), "")
    d$component <- factor(comp_lab[d$component], levels = comp_lab)
    m <- rowMeans(V); names(m) <- comp_lab
    d$species <- factor(d$species, levels = splab[order(V[nrow(V), ])])
    p <- ggplot2::ggplot(d, ggplot2::aes(value, species, fill = component)) +
      ggplot2::geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
      ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = function(x) paste0(100 * x, "%"),
                                  expand = ggplot2::expansion(mult = c(0, 0.04))) +
      ggplot2::scale_fill_manual(values = stats::setNames(c(eco_categorical(nrow(V) - 1), "#9e9e9e"), comp_lab),
                                 labels = function(x) strwrap_lab(x, 28)) +
      ggplot2::labs(x = fig_label("vp_x", "Share of explained variance", opts), y = NULL, fill = NULL) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, position = "bottom"))
    attr(p, "fig_size") <- list(height = 55 + 4.2 * ns, width = "onehalf")
    figs$variance_partitioning <- with_caption(p, sprintf(
      "Variance partitioning of the linear predictor for each species among covariates and residual associations (latent site factors). Community averages: %s.",
      paste(sprintf("%s %.0f%%", names(m), 100 * m), collapse = "; ")))
  }

  # trait_effects: Gamma with 95% and 50% intervals
  gn <- jsdm_comm_params(fit$estimates$param, pr$response)
  pretty <- function(p) {
    if (grepl("^comm_(mean|popavg)_", p)) {
      k <- sub("^comm_(mean|popavg)_", "", p); return(paste0("Community mean: ", jsdm_cov_label(k, fit, opts)))
    }
    parts <- regmatches(p, regexec(jsdm_trait_re(fit), p))[[1]]
    tl <- if (!is.null(pr$Tr$scaling[[parts[2]]])) pr$Tr$scaling[[parts[2]]]$label else parts[2]
    tgt <- if (parts[3] == "intercept") "baseline" else jsdm_cov_label(parts[3], fit, opts)
    paste0(fig_label(parts[2], tl, opts), " → ", tgt)
  }
  if (length(gn)) {
    D <- fit$draws[, gn, drop = FALSE]
    a <- (1 - ci) / 2; b <- (1 - inner) / 2
    e <- data.frame(param = gn, estimate = colMeans(D),
                    lower = apply(D, 2, stats::quantile, a), upper = apply(D, 2, stats::quantile, 1 - a),
                    lower50 = apply(D, 2, stats::quantile, b), upper50 = apply(D, 2, stats::quantile, 1 - b))
    e$param <- vapply(gn, function(p) fig_label(p, pretty(p), opts), "")
    e$group <- ifelse(grepl("^comm_", gn), "Community mean response", "Trait effect on response")
    p <- fig_coef(e, null_value = 0, xlab = sprintf("Effect on %s scale\n(per SD; %d%% and %d%% credible intervals)",
                                                    if (presence) "population-averaged probit" else "log", round(100 * ci), round(100 * inner)),
                  order = FALSE) +
      ggplot2::scale_colour_manual(values = c("Community mean response" = "#0072B2", "Trait effect on response" = "#D55E00")) +
      ggplot2::guides(colour = ggplot2::guide_legend(position = "bottom"))
    attr(p, "fig_size") <- list(height = 35 + 8 * length(gn), width = "double")
    figs$trait_effects <- with_caption(p, sprintf(
      "Community-level parameters%s. Blue: average response across species (per 1 SD of the covariate); orange: change in species' responses per 1 SD of the trait (trait-environment interaction; 'baseline' = effect on prevalence/abundance at average conditions). Points are posterior means, thick and thin bars %d%% and %d%% credible intervals; dashed line = no effect.",
      if (presence) " (population-averaged probit scale: species responses with latent site factors integrated out, regressed on traits; these describe the modelled species, and whether a trait effect extends to the wider species pool is judged from Hmsc's Gamma, reported in the text)" else " (Gamma, log scale)",
      round(100 * inner), round(100 * ci)))
  }

  # richness_<cov> and occurrence_<cov>
  obs_rich <- rowSums(if (presence) pr$Y > 0 else pr$Y)
  for (v in jsdm_gradient_covs(fit, data)) {
    g <- jsdm_gradient(fit, data, v)
    xx <- fig_x(g$x, v, opts)
    b <- summ_band(g$community, ci, inner); b$x <- xx
    ylab <- if (presence) fig_label("richness", "Expected species richness", opts) else fig_label("total_abundance", "Expected total abundance", opts)
    obs <- if (isTRUE(opts$show_data)) data.frame(x = fig_x(data$site_data[[v]], v, opts), y = obs_rich) else data.frame(x = fig_x(data$site_data[[v]], v, opts))
    xl <- fig_label(v, v, opts)
    p <- fig_cond_effect(b, xl, ylab, obs = obs)
    figs[[paste0("richness_", v)]] <- with_caption(p, sprintf(
      "%s along %s, other covariates at their mean (line: posterior mean; bands: %d%% and %d%% credible intervals; latent site factors averaged out). Points: observed %s per site.",
      ylab, xl, round(100 * ci), round(100 * inner), if (presence) "richness" else "total count"))
    # species curves, one panel per species
    d <- do.call(rbind, lapply(seq_len(ns), function(j) {
      bb <- summ_band(g$sp[, j, ], ci, inner); bb$x <- xx; bb$species <- splab[j]; bb
    }))
    d$species <- factor(d$species, levels = splab)
    p <- ggplot2::ggplot(d, ggplot2::aes(x, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.18) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25) +
      ggplot2::geom_line(colour = eco_pal$band, linewidth = 0.6) +
      ggplot2::facet_wrap(~species, ncol = if (ns > 12) 5 else 4, scales = if (presence) "fixed" else "free_y",
                          labeller = ggplot2::as_labeller(function(x) strwrap_lab(x, 20))) +
      ggplot2::scale_x_continuous(breaks = function(l) base::pretty(l, n = 3)) +
      ggplot2::labs(x = xl, y = if (presence) fig_label("occ_prob", "Probability of occurrence", opts) else fig_label("exp_count", "Expected count", opts)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = ggplot2::rel(0.75)),
                     panel.spacing.x = ggplot2::unit(2.4, "lines"))
    if (presence) p <- p + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1))
    attr(p, "fig_size") <- list(width = "double", height = 18 + 30 * ceiling(ns / (if (ns > 12) 5 else 4)))
    figs[[paste0("occurrence_", v)]] <- with_caption(p, sprintf(
      "Species-specific %s along %s (other covariates at their mean, latent factors averaged out); bands are %d%% and %d%% credible intervals.",
      if (presence) "probability of occurrence" else "expected count", xl, round(100 * ci), round(100 * inner)))
  }

  # model fit
  if (!is.null(fit$model_fit)) {
    mf <- fit$model_fit
    yv <- if (presence) "TjurR2" else "SR2"
    if (yv %in% names(mf)) {
      d <- data.frame(prev = mf$prevalence, r2 = mf[[yv]], species = splab)
      p <- ggplot2::ggplot(d, ggplot2::aes(prev, r2)) + ggplot2::geom_point(colour = eco_pal$ink, size = 1.6) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = fig_label("prevalence", "Observed prevalence (share of sites where recorded)", opts),
                      y = if (presence) "Tjur R²" else "Pseudo-R²")
      if (requireNamespace("ggrepel", quietly = TRUE) && ns <= 20)
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = species), size = 2, colour = eco_pal$muted, max.overlaps = 30,
                                          box.padding = 0.4, min.segment.length = 0, seed = 1)
      attr(p, "fig_size") <- list(width = "onehalf", aspect = 0.75)
      figs$model_fit <- with_caption(p, sprintf(
        "In-sample explanatory power per species (%s, including latent site factors) against prevalence. Mean %.2f.",
        if (presence) "Tjur R²" else "pseudo-R²", mean(d$r2, na.rm = TRUE)))
    }
  }
  if (!is.null(fit$ppc)) {
    pp <- fit$ppc
    figs$ppc <- with_caption(fig_ppc(pp$replicated_new_sites[, 1], pp$observed[1], "Variance of site species richness"),
      "Posterior predictive check: variance of site richness in data replicated from the posterior at new sites (latent factors redrawn) vs the observed value (line).")
  }
  # trace and posterior of the community-level parameters
  if (!is.null(fit$chains_gamma) && length(gn)) {
    nice <- strwrap_lab(vapply(gn, function(p) fig_label(p, pretty(p), opts), ""), 34)
    ren <- function(M) { M <- M[, gn, drop = FALSE]; colnames(M) <- nice; M }
    pt <- fig_trace(lapply(fit$chains_gamma, ren)) +
      ggplot2::facet_wrap(~ind, scales = "free_y", ncol = 2) +
      ggplot2::scale_x_continuous(breaks = function(l) base::pretty(l, n = 3)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = ggplot2::rel(0.8)))
    attr(pt, "fig_size") <- list(width = "double", height = 20 + 32 * ceiling(length(gn) / 2))
    figs$trace <- with_caption(pt,
      "MCMC traces of the community-level parameters, one colour per chain (post burn-in, thinned).")
    pp <- fig_posterior(ren(fit$draws), nice, ncol = 2) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = ggplot2::rel(0.8)))
    attr(pp, "fig_size") <- list(width = "double", height = 20 + 32 * ceiling(length(gn) / 2))
    figs$posterior <- with_caption(pp, "Posterior densities of the community-level parameters.")
  }
  figs
}

plots_jsdm <- function(fit, data) {
  f <- figures_jsdm(fit, data, fig_opts())
  keep <- c("species_responses", "residual_associations", "variance_partitioning", "trait_effects",
            grep("^richness_", names(f), value = TRUE), "model_fit", "ppc")
  lapply(f[intersect(keep, names(f))], function(p) p + theme_eco(10))
}
