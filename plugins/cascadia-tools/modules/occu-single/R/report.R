# occu-single: conditional effects, interpretation and methods text.

occ_sub_label <- function(s) c(psi = "occupancy", p = "detection probability per survey",
                               col = "colonisation probability", ext = "local extinction probability")[[s]]
occ_axis_label <- function(s) c(psi = "occupancy", p = "detection prob. per survey",
                                col = "colonisation prob.", ext = "local extinction prob.")[[s]]
occ_sub_short <- function(s, model) c(psi = if (model == "dynamic") "Initial occupancy" else "Occupancy",
                                      p = "Detection per survey", col = "Colonisation", ext = "Extinction")[[s]]

# Response-scale predictions for one submodel over a newdata frame given in
# RAW covariate units (numeric) or level names (factors). Variables not in
# newdata are held at their mean (numeric) or reference level (factor).
occ_predict <- function(fit, sub, newdata, prob = 0.95) {
  o <- fit$object; des <- o$des
  f <- stats::as.formula(des$forms[[sub]])
  vars <- all.vars(f)
  n <- max(1, nrow(newdata))
  nd <- data.frame(row.names = seq_len(n))
  for (v in vars) {
    if (!is.null(des$levels[[v]])) {
      val <- if (v %in% names(newdata)) as.character(newdata[[v]]) else rep(des$levels[[v]][1], n)
      nd[[v]] <- factor(val, levels = des$levels[[v]])
    } else nd[[v]] <- if (v %in% names(newdata)) newdata[[v]] else rep(des$scaling[[v]][["mean"]], n)
  }
  X <- occ_apply_scale(stats::model.matrix(f, nd), des$colscale[[sub]])
  nm <- paste0("beta_", sub, "[", seq_along(des$colnames[[sub]]), "]")
  B <- o$draws[, nm, drop = FALSE]
  P <- stats::plogis(X %*% t(B))
  a <- (1 - prob) / 2
  est <- if (!is.null(o$point)) as.vector(stats::plogis(X %*% o$point[nm])) else rowMeans(P)
  data.frame(estimate = est, lower = apply(P, 1, stats::quantile, a), upper = apply(P, 1, stats::quantile, 1 - a),
             draws_mean = rowMeans(P))
}

occ_cap <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))
occ_fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
occ_ci <- function(e, l, u, d = 2) sprintf("%s (95%% CI %s to %s)", occ_fmt(e, d), occ_fmt(l, d), occ_fmt(u, d))

occ_pretty_unit <- function(s) 10^floor(log10(s) + 1e-6)

# " (with observer = AB, the reference level)" when a submodel has factor covariates:
# intercept-scale probabilities refer to the reference level, not an average over levels.
occ_ref_text <- function(des, subs) {
  v <- unique(unlist(lapply(subs, function(s) all.vars(stats::as.formula(des$forms[[s]])))))
  v <- v[v %in% names(des$levels)]
  if (!length(v)) return("")
  sprintf(" (with %s, the reference level%s)", paste(sprintf("%s = %s", v, vapply(v, function(k) des$levels[[k]][1], "")), collapse = ", "),
          if (length(v) > 1) "s" else "")
}

occ_interpret <- function(fit, data) {
  o <- fit$object; des <- o$des; est <- fit$estimates
  g <- function(p) est[est$param == p, , drop = FALSE]
  out <- character()
  if (o$engine == "naive")
    out <- c(out, "WARNING: this is the naive control analysis that ignores imperfect detection. Its occupancy estimates are biased low and must not be reported as occupancy.")
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the model did not converge cleanly or some parameters are at a boundary (see Diagnostics). The numbers below should not be interpreted until this is resolved.")
  y0 <- des$y; y0[des$w == 0] <- 0
  naive <- colMeans(apply(y0, c(1, 3), max))
  nsurv <- stats::median(apply(des$w, c(1, 3), sum))
  inames <- occ_intercept_names(des$model)
  # occupancy
  r <- g(inames[["psi"]])
  if (des$model == "static") {
    out <- c(out, sprintf(
      "The species was detected at least once at %d of %d sites (naive occupancy %.2f). Allowing for missed detections, the estimated probability that a site with average covariate values%s is occupied is %s.",
      sum(apply(y0, 1, max)), des$M, naive[1], occ_ref_text(des, "psi"), occ_ci(r$estimate, r$lower, r$upper)))
    pr <- g("prop_sites_occupied")
    out <- c(out, sprintf(
      "Of the %d surveyed sites, an estimated %.0f (95%% CI %.0f to %.0f) were actually occupied (proportion %s); this finite-sample estimate applies to these sites only, not to the wider landscape.",
      des$M, pr$estimate * des$M, pr$lower * des$M, pr$upper * des$M, occ_ci(pr$estimate, pr$lower, pr$upper)))
  } else {
    out <- c(out, sprintf("In season %s the probability that a site with average covariate values%s was occupied is estimated at %s (naive: %.2f).",
                          des$season_labels[1], occ_ref_text(des, "psi"), occ_ci(r$estimate, r$lower, r$upper), naive[1]))
    rc <- g(inames[["col"]]); re <- g(inames[["ext"]])
    out <- c(out, sprintf(
      "Between consecutive seasons, an unoccupied site with average covariates%s became occupied with probability %s (colonisation), and an occupied site became unoccupied with probability %s (local extinction). If these rates persisted, occupancy would settle near %.2f (colonisation / (colonisation + extinction)).",
      occ_ref_text(des, c("col", "ext")), occ_ci(rc$estimate, rc$lower, rc$upper), occ_ci(re$estimate, re$lower, re$upper), rc$estimate / (rc$estimate + re$estimate)))
    p1 <- g("prop_sites_occupied_t1"); pT <- g(paste0("prop_sites_occupied_t", des$T)); ch <- g("prop_occupied_change")
    if (nrow(ch)) out <- c(out, sprintf(
      "Among the %d surveyed sites, the estimated proportion occupied went from %s in season %s to %s in season %s, a change of %s.",
      des$M, occ_ci(p1$estimate, p1$lower, p1$upper), des$season_labels[1], occ_ci(pT$estimate, pT$lower, pT$upper),
      des$season_labels[des$T], occ_ci(ch$estimate, ch$lower, ch$upper)))
  }
  # detection
  if (o$engine != "naive") {
    r <- g(inames[["p"]])
    pstar <- 1 - (1 - c(r$estimate, r$lower, r$upper))^nsurv
    out <- c(out, sprintf(
      "Per-survey detection probability under average survey conditions%s is %s. With the typical %d surveys per site%s, an occupied site is detected at least once%s with probability %s, so about %.0f%% of occupied sites go undetected with this survey effort.",
      occ_ref_text(des, "p"), occ_ci(r$estimate, r$lower, r$upper), nsurv, if (des$model == "dynamic") " per season" else "", if (des$model == "dynamic") " within a season" else "",
      occ_ci(pstar[1], pstar[2], pstar[3]), 100 * (1 - pstar[1])))
    if (pstar[1] < 0.5 || r$estimate < 0.2)
      out <- c(out, sprintf("Caution: cumulative detection is low (%.2f). In this module's simulation tests, designs with 2 surveys and p = 0.2 (cumulative 0.36) gave unreliable occupancy estimates (see VALIDATION.md); treat occupancy as imprecise and consider more surveys.", pstar[1]))
  }
  # covariate effects
  subs <- if (o$engine == "naive") character() else des$subs
  for (s in subs) {
    f <- stats::as.formula(des$forms[[s]])
    for (v in all.vars(f)) {
      if (!is.null(des$levels[[v]])) {
        lv <- des$levels[[v]]
        pd <- occ_predict(fit, s, stats::setNames(data.frame(lv), v))
        out <- c(out, sprintf("%s by %s: %s.", occ_sub_short(s, des$model), v,
                              paste(sprintf("%s %s", lv, occ_ci(pd$estimate, pd$lower, pd$upper)), collapse = "; ")))
        next
      }
      x <- des$raw[[v]]; q <- stats::quantile(x, c(0.1, 0.9), na.rm = TRUE)
      pd <- occ_predict(fit, s, stats::setNames(data.frame(q), v))
      sen <- sprintf("%s and %s: going from the 10th to the 90th percentile of %s observed in the data (%s to %s), %s changes from %s to %s.",
                     occ_sub_short(s, des$model), v, v, signif(q[1], 3), signif(q[2], 3), occ_sub_label(s),
                     occ_ci(pd$estimate[1], pd$lower[1], pd$upper[1]), occ_ci(pd$estimate[2], pd$lower[2], pd$upper[2]))
      b <- g(occ_coef_public(s, v))
      if (nrow(b)) {
        sdv <- des$colscale[[s]]["sd", v]; u <- occ_pretty_unit(sdv)
        dir <- if (b$lower > 0) "positively" else if (b$upper < 0) "negatively" else "not clearly"
        sen <- paste(sen, sprintf(
          "On the logit scale the coefficient is %s per 1 SD of %s (1 SD = %s): each 1 SD increase multiplies the odds by %s, and each %s-unit increase by %.2f. %s is %s associated with %s%s.",
          occ_ci(b$estimate, b$lower, b$upper), v, signif(sdv, 3),
          occ_ci(exp(b$estimate), exp(b$lower), exp(b$upper)), format(u), exp(b$estimate * u / sdv),
          occ_cap(occ_sub_label(s)), dir, v,
          if (dir == "not clearly") " (the 95% interval includes no effect)" else ""))
      }
      if (s == "p") sen <- paste(sen, "This is an effect on detection (how easily the species is found when present), not on where it occurs.")
      if (s %in% c("col", "ext") && v %in% des$sources$yearly)
        sen <- paste(sen, sprintf("(%s is taken from season t and applies to the change between season t and season t+1.)", v))
      out <- c(out, sen)
    }
    # transformed terms / interactions: coefficient per 1 SD of the term itself
    cols <- des$colnames[[s]][-1]
    fac <- unlist(lapply(names(des$levels), function(v) paste0(v, des$levels[[v]])))
    for (cc in setdiff(cols, c(all.vars(f), fac))) {
      b <- g(occ_coef_public(s, cc))
      if (nrow(b)) out <- c(out, sprintf("%s model term %s: logit-scale coefficient %s per 1 SD of the term (1 SD = %s); odds ratio %s.",
                                         occ_sub_short(s, des$model), cc, occ_ci(b$estimate, b$lower, b$upper),
                                         signif(des$colscale[[s]]["sd", cc], 3), occ_ci(exp(b$estimate), exp(b$lower), exp(b$upper))))
    }
  }
  if (o$engine != "naive") out <- c(out,
    "These are associations from an observational design; they do not by themselves show that a covariate causes the change.")
  gof <- o$gof
  if (!is.null(gof)) {
    lowpow <- "Note that this test has low power against unmodelled differences in detection among sites (in validation it flagged only about 1 in 10 data sets where such differences biased occupancy by about -20%), so a pass rules out gross misfit, not all bias."
    if (o$engine == "nimble") {
      # Posterior predictive check: only the p-value is interpretable. The
      # observed/replicated ratio is not a bootstrap c-hat (VALIDATION.md known
      # limit 4), so it is neither used for the verdict nor for SE inflation.
      msg <- sprintf("Goodness of fit (posterior predictive check of the MacKenzie-Bailey statistic, %d posterior draws): Bayesian p = %.3f (ratio of observed to replicated statistic %.2f; this ratio is not a c-hat and should not be used to inflate uncertainty). ",
                     gof$nboot, gof$p_value, gof$c_hat)
      msg <- paste0(msg, if (gof$p_value < 0.05)
        "The detection histories are less consistent with the model than data simulated from it (often unmodelled heterogeneity in detection or occupancy). Consider adding detection or occupancy covariates, and check the same model with the unmarked engine, whose bootstrap c-hat can be used to inflate standard errors."
        else paste("No evidence of lack of fit.", lowpow))
    } else {
      msg <- sprintf("Goodness of fit (MacKenzie-Bailey test, %d bootstrap replicates): c-hat = %.2f, p = %.3f. ", gof$nboot, gof$c_hat, gof$p_value)
      msg <- paste0(msg, if (gof$p_value < 0.05 || gof$c_hat > 1.5)
        sprintf("The model does not fully fit the detection histories (overdispersion, often unmodelled heterogeneity in detection or occupancy). Consider adding covariates; if kept, inflate standard errors by sqrt(c-hat) = %.2f.", sqrt(max(gof$c_hat, 1)))
        else paste("No evidence of lack of fit.", lowpow))
    }
    out <- c(out, msg)
  } else if (o$engine == "unmarked") out <- c(out, "Goodness of fit was not run (settings gof_nboot = 0).")
  out
}

occ_plots <- function(fit, data) {
  o <- fit$object; des <- o$des; pl <- list()
  if (o$engine == "naive") return(pl)
  for (s in des$subs) {
    for (v in all.vars(stats::as.formula(des$forms[[s]]))) {
      ttl <- sprintf("%s vs %s", occ_sub_short(s, des$model), v)
      if (!is.null(des$levels[[v]])) {
        lv <- des$levels[[v]]
        pd <- occ_predict(fit, s, stats::setNames(data.frame(lv), v)); pd$x <- factor(lv, levels = lv)
        pl[[ttl]] <- ggplot2::ggplot(pd, ggplot2::aes(x, estimate)) +
          ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), colour = "#2f6f4f") +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(x = v, y = occ_axis_label(s), title = ttl,
                        caption = "Other covariates at their mean / reference level\n95% intervals") +
          ggplot2::theme_minimal(base_size = 12)
        next
      }
      x <- des$raw[[v]]
      grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
      pd <- occ_predict(fit, s, stats::setNames(data.frame(grid), v)); pd$x <- grid
      pl[[ttl]] <- cond_effect_plot(pd, xlab = v, ylab = occ_axis_label(s), observed_x = x, title = ttl) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(caption = "Other covariates at their mean (factors at reference level)\n95% band; rug = observed values")
    }
  }
  if (des$model == "dynamic") {
    e <- fit$estimates; nm <- paste0("prop_sites_occupied_t", seq_len(des$T))
    d <- e[match(nm, e$param), ]; d$season <- factor(des$season_labels, levels = des$season_labels)
    y0 <- des$y; y0[des$w == 0] <- 0
    d$naive <- colMeans(apply(y0, c(1, 3), max))
    pl[["Proportion of surveyed sites occupied by season"]] <- ggplot2::ggplot(d, ggplot2::aes(season, estimate, group = 1)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#2f6f4f", alpha = 0.22) +
      ggplot2::geom_line(colour = "#2f6f4f", linewidth = 1) + ggplot2::geom_point(colour = "#2f6f4f") +
      ggplot2::geom_point(ggplot2::aes(y = naive), shape = 4, size = 2.5, colour = "grey30") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "Season", y = "Proportion of sites occupied",
                    caption = "Line and band: estimate with 95% interval\nCrosses: naive (detected at least once)") +
      ggplot2::theme_minimal(base_size = 12)
  }
  pl
}

occ_methods_text <- function(fit, data, settings) {
  o <- fit$object; des <- o$des; st <- o$settings
  fm <- function(s) des$forms[[s]]
  model <- if (des$model == "static") sprintf(
    "We fitted a single-season occupancy model (MacKenzie et al. 2002) to detection/non-detection data from %d sites with up to %d repeat surveys each (%d surveys in total). Occupancy probability (psi) was modelled as logit(psi) %s and per-survey detection probability (p) as logit(p) %s.",
    des$M, des$J, sum(des$w), fm("psi"), fm("p"))
  else sprintf(
    "We fitted a dynamic (multi-season) occupancy model (MacKenzie et al. 2003) to detection/non-detection data from %d sites over %d seasons with up to %d surveys per season (%d surveys in total). Sites were assumed closed to changes in occupancy within a season. Initial occupancy was modelled as logit(psi1) %s, colonisation as logit(gamma) %s, local extinction as logit(epsilon) %s (covariates for the transition from season t to t+1 taken from season t) and per-survey detection as logit(p) %s.",
    des$M, des$T, des$J, sum(des$w), fm("psi"), fm("col"), fm("ext"), fm("p"))
  stdz <- if (isTRUE(st$standardize)) " Each numeric model term was standardised (mean 0, SD 1 over the sites, transitions or surveys it applies to), so intercepts are probabilities at average covariate values (factors at their reference level) and slopes are per 1 SD." else ""
  eng <- switch(o$engine,
    unmarked = sprintf("The model was fitted by maximum likelihood with unmarked %s (Fiske & Chandler 2011); 95%% confidence intervals are Wald intervals on the logit scale. The finite-sample proportion of surveyed sites occupied was computed from the conditional (empirical Bayes) probability that each site was occupied given its detection history, propagating parameter uncertainty by drawing %d parameter sets from the asymptotic normal distribution of the estimates. Goodness of fit was assessed with the MacKenzie & Bailey (2004) chi-square test on detection histories (cohorts defined by missing-survey pattern%s) with a %d-replicate parametric bootstrap; c-hat is the observed statistic divided by the bootstrap mean.",
                       as.character(utils::packageVersion("unmarked")), st$n_prop_draws,
                       if (des$model == "dynamic") "; computed within each season using model-implied occupancy and summed" else "", st$gof_nboot),
    nimble = sprintf("The model was fitted in a Bayesian framework with NIMBLE %s (de Valpine et al. 2017), with the latent occupancy state of every site sampled explicitly (Royle & Kery 2007). Priors: Logistic(0,1) on each intercept (uniform on the probability scale) and Normal(0, SD = %g) on each slope. We ran %d chains of %d iterations (%d burn-in)%s and checked convergence with R-hat and bulk effective sample size. The proportion of surveyed sites occupied is the posterior of the mean latent state. Fit was checked with a posterior predictive test of the MacKenzie & Bailey (2004) chi-square statistic (%d draws).",
                     as.character(utils::packageVersion("nimble")), st$prior_sd, st$nchains, st$niter, st$nburnin,
                     if (isTRUE(st$slice)) ", with slice samplers on the regression coefficients" else "", st$gof_ndraws),
    naive = "NEGATIVE CONTROL: logistic regression of whether the species was ever detected at a site, ignoring imperfect detection.")
  paste0(model, stdz, " ", eng)
}
