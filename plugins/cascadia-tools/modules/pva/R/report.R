# Interpretation, methods text and publication figures for the pva module.

pva_pct <- function(p) ifelse(p < 0.001, "<0.1%", ifelse(p > 0.999, ">99.9%", ifelse(p < 0.01 | p > 0.99, sprintf("%.1f%%", 100 * p), sprintf("%.0f%%", 100 * p))))
pva_f <- function(x, d = 2) formatC(x, digits = d, format = "f")
pva_cap <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))
pva_rate_labels <- function(rates) {
  d <- if (!is.null(rates$description)) as.character(rates$description) else rep(NA_character_, nrow(rates))
  stats::setNames(ifelse(is.na(d) | !nzchar(d), rates$rate, paste0(gsub("\\s*\\(([^)]*)\\)", ", \\1", d), " [", rates$rate, "]")), rates$rate)
}
pva_row <- function(est, p) est[est$param == p, , drop = FALSE]
pva_risk_word <- function(lo, hi) {
  if (hi - lo > 0.6) "the data cannot pin this risk down: the interval spans most of the 0-1 range"
  else if (hi < 0.05) "the risk is low under the model's assumptions"
  else if (lo > 0.5) "the risk is high under the model's assumptions"
  else "the risk is uncertain: the interval includes both low and substantial values"
}

# ---- interpretation ----------------------------------------------------------

count_interpret <- function(fit, data) {
  o <- fit$object; e <- fit$estimates
  mu <- pva_row(e, "growth_rate_mu"); lam <- pva_row(e, "lambda_stochastic"); s2 <- pva_row(e, "process_variance_sigma2")
  u <- data$meta$units; sp <- data$meta$species
  yrs <- range(o$df$year)
  meth <- if (o$method == "diffusion") "the diffusion approximation (Dennis et al. 1991), which treats every count as exact"
          else "an exponential-growth state-space model that separates counting error from real year-to-year change (Bayesian, exact grid integration over the REML likelihood; Staples et al. 2004, Humbert et al. 2009)"
  dirn <- if (lam$upper < 1) "declining" else if (lam$lower > 1) "increasing" else "not clearly increasing or decreasing"
  out <- c(
    sprintf("Using %d counts of %s from %g to %g and %s, the population was %s: the long-run (median) growth rate is lambda = %s per year (95%% CI %s to %s), i.e. a typical change of %+.1f%% per year (95%% CI %+.1f%% to %+.1f%%).",
            nrow(o$df), sp, yrs[1], yrs[2], meth, dirn, pva_f(lam$estimate, 3), pva_f(lam$lower, 3), pva_f(lam$upper, 3),
            100 * (lam$estimate - 1), 100 * (lam$lower - 1), 100 * (lam$upper - 1)),
    sprintf("Year-to-year variability in the log growth rate (process variance sigma^2) is %s (95%% CI %s to %s); this means that in roughly one year in three the annual change departs from the trend by more than about %.0f%% (one standard deviation, sigma = %s on the log scale).",
            pva_f(s2$estimate, 4), pva_f(s2$lower, 4), pva_f(s2$upper, 4), 100 * (exp(sqrt(s2$estimate)) - 1), pva_f(sqrt(s2$estimate), 3)))
  if (o$method == "state_space") {
    t2 <- pva_row(e, "obs_variance_tau2"); nc <- pva_row(e, "N_current")
    out <- c(out, sprintf("Counting (observation) error variance on the log scale is %s (95%% CI %s to %s), equivalent to a coefficient of variation of roughly %.0f%% in a single count. The model's estimate of the true current abundance is %s %s (95%% CI %s to %s).",
                          pva_f(t2$estimate, 4), pva_f(t2$lower, 4), pva_f(t2$upper, 4), 100 * sqrt(exp(t2$estimate) - 1),
                          format(round(nc$estimate)), u, format(round(nc$lower)), format(round(nc$upper))))
  }
  rk <- vapply(o$horizons, function(h) {
    r <- pva_row(e, paste0("p_qext_", h, "yr"))
    sprintf("%d years: %s (95%% CI %s to %s)", h, pva_pct(r$estimate), pva_pct(r$lower), pva_pct(r$upper))
  }, "")
  hl <- pva_row(e, paste0("p_qext_", max(o$horizons), "yr"))
  out <- c(out,
    sprintf("Probability that the population falls to or below the quasi-extinction threshold of %g %s in at least one annual census within %s. At the longest horizon %s.",
            o$threshold, u, paste(rk, collapse = "; "), pva_risk_word(hl$lower, hl$upper)),
    sprintf("How to read these numbers: quasi-extinction probabilities assume the past %d years of growth and variability continue unchanged, with no density dependence, catastrophes outside the observed range, or management change. Their %s come from the uncertainty in the estimated growth rate and variance%s (%s, %d draws). Forecasts beyond about 20%% of the monitoring span (here %.0f years) are known to be very imprecise (Fieberg & Ellner 2000); use them to compare options, not as absolute predictions. In ecotools simulations the point estimate was as often above as below the true risk, but it tended to be more extreme than the truth (too high when high, too low when low), whereas the 95%% interval contained the true risk in 88-99%% of simulated data sets (the lowest values, 0.88-0.89, at a 5-year horizon where almost all misses were at true risks below 0.3%%; VALIDATION.md).",
            round(o$span), if (o$method == "diffusion") "confidence intervals" else "credible intervals",
            if (o$method == "diffusion") "" else " and in the current true abundance",
            if (o$method == "diffusion") "parametric bootstrap" else "posterior", o$nboot, 0.2 * o$span))
  if (o$correction == "discrete")
    out <- c(out, sprintf("The classic continuous-time formula (which assumes the population is watched continuously) gives slightly higher values (%s at %d years); the values above are corrected for once-a-year censuses, which matches how quasi-extinction is actually detected.",
                          pva_pct(o$classic[length(o$classic)]), max(o$horizons)))
  out
}

matrix_interpret <- function(fit, data) {
  o <- fit$object; e <- fit$estimates
  u <- data$meta$units; sp <- data$meta$species
  ld <- pva_row(e, "lambda_deterministic"); ls <- pva_row(e, "lambda_stochastic"); lt <- pva_row(e, "lambda_stochastic_tuljapurkar")
  out <- c(
    sprintf("With the mean vital rates, %s would change by a factor of %s per year (deterministic lambda; 95%% CI %s to %s from the uncertainty in the vital-rate estimates).",
            sp, pva_f(ld$estimate, 3), pva_f(ld$lower, 3), pva_f(ld$upper, 3)),
    sprintf("Because vital rates vary from year to year, the long-run growth rate %s: stochastic lambda = %s (95%% CI %s to %s) by simulation, %s by Tuljapurkar's small-variance approximation. %s",
            if (ls$estimate < ld$estimate - 0.005) "is lower" else if (ls$estimate < ld$estimate - 0.0005) sprintf("is slightly lower (by %s)", pva_f(ld$estimate - ls$estimate, 4)) else "is about the same (environmental variation is small relative to growth)",
            pva_f(ls$estimate, 3), pva_f(ls$lower, 3), pva_f(ls$upper, 3), pva_f(lt$estimate, 3),
            if (ls$upper < 1) "The population is projected to decline in the long run."
            else if (ls$lower > 1) "The population is projected to grow in the long run (before any ceiling)."
            else "The data cannot tell whether the population will grow or decline in the long run."))
  rk <- vapply(o$horizons, function(h) {
    r <- pva_row(e, paste0("p_qext_", h, "yr"))
    sprintf("%d years: %s (95%% CI %s to %s)", h, pva_pct(r$estimate), pva_pct(r$lower), pva_pct(r$upper))
  }, "")
  hl <- pva_row(e, paste0("p_qext_", max(o$horizons), "yr"))
  cs <- if (is.null(data$count_stages)) "all stages" else paste(data$count_stages, collapse = " + ")
  out <- c(out, sprintf("Probability that the number of %s (%s) falls to or below %g in any year within %s. At the longest horizon %s.",
                        u, cs, o$threshold, paste(rk, collapse = "; "), pva_risk_word(hl$lower, hl$upper)))
  el <- sort(o$vr_elast, decreasing = TRUE)
  rl <- pva_rate_labels(o$rates); names(el) <- ifelse(names(el) %in% names(rl), rl[names(el)], names(el))
  out <- c(out, sprintf("Elasticities: a 1%% increase in %s raises lambda by %.2f%%; the corresponding values for the other rates are %s. Rates with high elasticity are where a given proportional improvement does most for growth, but they are often the hardest to change (e.g. adult survival); weigh elasticity against what management can realistically achieve.",
                        names(el)[1], el[1], paste(sprintf("%s %.2f%%", names(el)[-1], el[-1]), collapse = ", ")))
  if (length(o$scenarios)) {
    h <- max(o$horizons)
    for (nm in names(o$scenarios)) {
      slug <- gsub("[^a-z0-9]+", "_", tolower(nm))
      a <- pva_row(e, paste0("scenario_", slug, "_p_qext_", h, "yr"))
      d <- pva_row(e, paste0("scenario_", slug, "_change_p_qext_", h, "yr"))
      l <- pva_row(e, paste0("scenario_", slug, "_lambda_stochastic"))
      out <- c(out, sprintf("Scenario '%s': the model projects a %d-year quasi-extinction risk of %s (baseline %s), a change of %+.1f percentage points (95%% CI %+.1f to %+.1f), and stochastic lambda %s.%s This is a projection of what the model implies if the change were achieved and sustained, not evidence that the action will achieve it.",
                            nm, h, pva_pct(a$estimate), pva_pct(hl$estimate), 100 * d$estimate, 100 * d$lower, 100 * d$upper, pva_f(l$estimate, 3),
                            if (abs(d$upper) < 0.0005 || abs(d$lower) < 0.0005) " An interval limit of 0 means that for some plausible vital-rate values the risk is near zero with or without the action, so the action would make no difference there." else ""))
    }
  }
  cm <- o$mx$corr; rho <- if (is.matrix(cm)) cm[upper.tri(cm)] else 0
  vary <- if (any(rho != 0)) sprintf("vital rates vary between years around constant means (no trends), with within-year correlation among rates (%s)",
                                     if (length(unique(rho)) == 1) sprintf("r = %.2f", rho[1]) else sprintf("r from %.2f to %.2f", min(rho), max(rho)))
          else "vital rates vary independently between years around constant means (no trends)"
  out <- c(out, sprintf("Assumptions: %s, %s, %s; parameter uncertainty is propagated from the standard errors %s. Extinction risks over long horizons are imprecise (Fieberg & Ellner 2000) and sensitive to the threshold and to the ceiling; compare scenarios with the same assumptions rather than quoting absolute risk. In ecotools simulations the point estimate of risk was as often above as below the truth but more extreme than it, while the 95%% interval contained the true risk in 96-100%% of simulated data sets (VALIDATION.md).",
                        vary, if (is.finite(data$ceiling %||% Inf)) sprintf("abundance is capped at a ceiling of %g", data$ceiling) else "there is no density dependence",
                        if (isTRUE(o$settings$demographic)) "demographic stochasticity is included (binomial survival, Poisson births)" else "demographic stochasticity is ignored",
                        if (is.null(data$rates) && !is.null(data$demography)) "estimated from the monitoring data" else "supplied"))
  out
}

# ---- methods ------------------------------------------------------------------

count_methods <- function(fit, data, settings) {
  o <- fit$object
  core <- if (o$method == "diffusion")
    sprintf("We estimated the parameters of the stochastic exponential growth (diffusion) model, log N(t+s) = log N(t) + mu*s + e, e ~ Normal(0, sigma^2*s), from %d counts (%g-%g) following Dennis et al. (1991) and Morris & Doak (2002): mu and sigma^2 were estimated by regressing log(N(t+s)/N(t))/sqrt(s) on sqrt(s) through the origin, with a t-based interval for mu and a chi-square interval for sigma^2. Counts were treated as exact.",
            nrow(o$df), min(o$df$year), max(o$df$year))
  else sprintf("We fitted the exponential-growth state-space model (EGSS: log N(t+s) = log N(t) + mu*s + process noise with variance sigma^2*s; log count = log N + observation error with variance tau^2; Dennis et al. 2006; Humbert et al. 2009) to %d counts (%g-%g); unsurveyed years were left out. Inference was Bayesian with flat priors on the initial log abundance and mu and uniform priors on the standard deviations sigma and tau (Gelman 2006). Because flat priors on the two fixed effects make the marginal likelihood of (sigma^2, tau^2) the restricted (REML) likelihood (Harville 1974; Staples et al. 2004), the joint posterior of (sigma, tau) was computed exactly on a %d x %d grid, and mu, the initial state and the true log abundances were then drawn from their exact Gaussian conditional posteriors (%d draws; no MCMC). Point estimates of mu, lambda, current abundance and quasi-extinction risk are posterior medians; for sigma^2 and tau^2 the point estimates are the REML estimates (the posterior median under uniform priors on the SDs is biased upward in short series); all intervals are 95%% equal-tailed credible intervals. Current abundance is the posterior of the true abundance in the last count year. REML estimates: sigma^2 = %s, tau^2 = %s, mu = %s.",
               nrow(o$df), min(o$df$year), max(o$df$year), length(o$ss$grid$sigma), length(o$ss$grid$tau), o$nboot,
               if (is.null(o$ss$reml)) "NA" else pva_f(o$ss$reml$s2, 4), if (is.null(o$ss$reml)) "NA" else pva_f(o$ss$reml$tau2, 4),
               if (is.null(o$ss$reml)) "NA" else pva_f(o$ss$reml$mu, 4))
  c(core, sprintf("The probability of falling to or below a quasi-extinction threshold of %g within t years was computed from the inverse-Gaussian first-passage distribution (Dennis et al. 1991; Holmes 2004)%s. %s Analyses used the ecotools pva module (R %s).",
                  o$threshold, if (o$correction == "discrete") ", with the barrier shifted by 0.5826*sigma to correct for annual rather than continuous censuses (Broadie, Glasserman & Kou 1997)" else "",
                  if (o$method == "diffusion") sprintf("Confidence intervals for these probabilities are percentiles of the probabilities recomputed for %d parameter sets drawn from the exact sampling distributions of the estimates (scaled chi-square for sigma^2, normal for mu given sigma^2; parametric bootstrap as in Dennis et al. 1991 and Holmes 2004). The point estimate is the plug-in value at the estimated parameters.", o$nboot)
                  else sprintf("The probability was computed for each of the %d posterior draws of (mu, sigma^2, current log abundance); we report the posterior median and 95%% credible interval.", o$nboot),
                  getRversion()))
}

matrix_methods <- function(fit, data, settings) {
  o <- fit$object; r <- o$rates
  src <- if (!is.null(data$demography) && is.null(data$rates))
    sprintf("Vital rates were estimated from %d years of monitoring data as the mean of annual estimates (standard error = SD of annual estimates / sqrt(years), multiplied by t(0.975, years-1)/1.96 to give small-sample t-type intervals); their environmental (temporal) standard deviations were estimated by subtracting the mean binomial (survival) or Poisson (fecundity) sampling variance from the variance of the annual estimates (cf. Kendall 1998; Morris & Doak 2002).", max(r$n_years, na.rm = TRUE))
  else "Vital-rate means, standard errors and environmental standard deviations were taken from the supplied estimates."
  c(src, sprintf("We projected a %d-stage matrix model (stages: %s) for %d years. Each year, survival-type rates were drawn from beta distributions and fecundities from %s distributions with the estimated means and environmental SDs%s (Morris & Doak 2002). %s%s Parameter uncertainty was propagated by drawing %d sets of mean vital rates from their sampling distributions (beta for probabilities, lognormal for fecundity%s) and running %d trajectories for each; point estimates use the estimated rates with %d trajectories. Stochastic lambda was computed by simulation (mean annual log growth of an environment-only projection) and by Tuljapurkar's approximation; vital-rate elasticities of the deterministic lambda by finite differences, checked against popbio::vitalsens (Stubben & Milligan 2007). %s",
                 length(o$mx$stages), paste(o$mx$stages, collapse = ", "), o$tmax,
                 if (any(r$dist == "stretched_beta")) "lognormal or stretched-beta" else "lognormal",
                 if (is.matrix(o$mx$corr) && any(o$mx$corr[upper.tri(o$mx$corr)] != 0)) ", correlated within years through a Gaussian copula" else ", independently",
                 if (isTRUE(o$settings$demographic)) "Demographic stochasticity was included (binomial survival and stage transitions, Poisson births). " else "",
                 if (is.finite(data$ceiling %||% Inf)) sprintf("Abundance was capped at a ceiling of %g. ", data$ceiling) else "",
                 o$settings$ndraw, if (any(is.finite(r$n_years))) "; environmental variances from scaled inverse chi-square distributions with years-1 degrees of freedom (the chi-square pivot)" else "",
                 o$settings$ntraj, o$settings$npt,
                 if (length(o$scenarios)) sprintf("Management scenarios (%s) were run with the same random numbers as the baseline so that differences are paired.", paste(names(o$scenarios), collapse = "; ")) else ""))
}

# ---- figures ------------------------------------------------------------------

pva_qext_fig <- function(cdf, opts, units, thr, extra = NULL, itype = "confidence interval from parameter uncertainty (parametric bootstrap)") {
  ci <- round(100 * (opts$ci %||% 0.95))
  pred <- data.frame(x = cdf$year, estimate = cdf$estimate, lower = cdf$lower, upper = cdf$upper,
                     lower50 = cdf$lower50, upper50 = cdf$upper50)
  top <- max(c(pred$upper, pred$estimate), na.rm = TRUE)
  zoom <- top < 0.5
  p <- fig_cond_effect(pred, fig_label("years_ahead", "Years from now", opts),
                       fig_label("p_qext", sprintf("P(N ≤ %g %s)", thr, units), opts), ylim01 = !zoom)
  if (zoom) p <- p + ggplot2::scale_y_continuous(limits = c(0, max(0.02, 1.08 * top)), expand = ggplot2::expansion(mult = c(0, 0.02)))
  with_caption(p, sprintf("Cumulative probability that abundance falls to or below the quasi-extinction threshold (%g %s) by each year. Line: estimate; dark band: 50%% and light band: %d%% %s.%s", thr, units, ci, itype,
                          if (zoom) sprintf(" The y axis is truncated at %.2f because all risks are low.", max(0.02, 1.08 * top)) else ""))
}

count_figures <- function(fit, data, opts = fig_opts()) {
  o <- fit$object; u <- data$meta$units
  figs <- list()
  figs$quasi_extinction <- pva_qext_fig(o$cdf, opts, u, o$threshold,
    itype = if (o$method == "diffusion") "confidence interval from parameter uncertainty (parametric bootstrap)" else "credible interval (posterior of growth rate, variance and current abundance)")
  # Fan chart: observed counts, smoothed states (state-space), projections
  set.seed(11)
  nb <- min(nrow(o$boot), 2000); bi <- sample.int(nrow(o$boot), nb, replace = nrow(o$boot) < nb)
  H <- o$tmax; last <- max(o$df$year)
  x <- if (!is.null(o$boot$xT)) o$boot$xT[bi] else o$xT + if (o$PT > 0) sqrt(o$PT) * stats::rnorm(nb) else 0
  qs <- matrix(NA, H, 5)
  for (h in seq_len(H)) {
    x <- x + o$boot$mu[bi] + sqrt(pmax(o$boot$s2[bi], 0)) * stats::rnorm(nb)
    qs[h, ] <- stats::quantile(exp(x), c(.025, .25, .5, .75, .975))
  }
  proj <- data.frame(time = last + seq_len(H), estimate = qs[, 3], lower = qs[, 1], upper = qs[, 5], lower50 = qs[, 2], upper50 = qs[, 4])
  p <- ggplot2::ggplot(proj, ggplot2::aes(time, estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.15) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25) +
    ggplot2::geom_line(colour = eco_pal$band, linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = o$threshold, linetype = 2, colour = eco_pal$truth, linewidth = 0.4)
  if (!is.null(o$ss)) {
    sm <- data.frame(time = o$df$year, estimate = exp(o$ss$xs), lower = exp(o$ss$xs_lo), upper = exp(o$ss$xs_hi))
    p <- p + ggplot2::geom_ribbon(data = sm, ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$muted, alpha = 0.2) +
      ggplot2::geom_line(data = sm, colour = eco_pal$muted, linewidth = 0.6)
  }
  p <- p + ggplot2::geom_point(data = data.frame(time = o$df$year, estimate = o$df$count), size = 1, colour = eco_pal$ink) +
    ggplot2::geom_line(data = data.frame(time = o$df$year, estimate = o$df$count), colour = eco_pal$ink, linewidth = 0.3) +
    ggplot2::scale_y_log10() + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.04)),
                                breaks = function(l) { b <- pretty(l, 5); b[b >= min(o$df$year) & b <= last + H] }) +
    ggplot2::annotate("text", x = min(o$df$year), y = o$threshold, label = "quasi-extinction threshold", hjust = 0, vjust = -0.5, size = 2.4, colour = eco_pal$truth) +
    ggplot2::labs(x = fig_label("year", "Year", opts), y = fig_label("count", sprintf("%s (log scale)", pva_cap(u)), opts))
  figs$trajectories <- with_caption(p, sprintf("Observed counts (points)%s and projected abundance: median (line), 50%% and 95%% of simulated futures (bands) including uncertainty in growth rate and variance. Dashed line: quasi-extinction threshold.",
                                               if (!is.null(o$ss)) ", state-space estimate of true abundance with 95% band (grey)" else ""))
  lam <- exp(o$boot$mu)
  p <- ggplot2::ggplot(data.frame(v = lam), ggplot2::aes(v)) +
    ggplot2::geom_histogram(bins = 40, fill = eco_pal$band, alpha = 0.6, colour = "white") +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = eco_pal$null) +
    ggplot2::geom_vline(xintercept = exp(o$mu), colour = eco_pal$ink, linewidth = 0.6) +
    ggplot2::labs(x = fig_label("lambda", "Stochastic growth rate, lambda = exp(mu)", opts), y = if (o$method == "diffusion") "Bootstrap replicates" else "Posterior draws")
  figs$lambda_distribution <- with_caption(p, sprintf("%s of the long-run growth rate. Solid line: estimate; dashed line: lambda = 1 (stable). Share of %s below 1: %s.",
    if (o$method == "diffusion") "Sampling distribution (parametric bootstrap)" else "Posterior distribution",
    if (o$method == "diffusion") "replicates" else "posterior draws", pva_pct(mean(lam < 1))))
  inc <- o$checks$inc
  p1 <- ggplot2::ggplot(inc, ggplot2::aes(year, increment)) + ggplot2::geom_hline(yintercept = o$mu, colour = eco_pal$band) +
    ggplot2::geom_point(size = 1, colour = eco_pal$ink) + ggplot2::geom_line(linewidth = 0.3, colour = eco_pal$muted) +
    ggplot2::labs(x = fig_label("year", "Year", opts), y = fig_label("log_growth", "Annual log growth rate", opts))
  p2 <- ggplot2::ggplot(inc, ggplot2::aes(logN, increment)) + ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = eco_pal$band, fill = eco_pal$band, alpha = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(size = 1, colour = eco_pal$ink) +
    ggplot2::labs(x = fig_label("logN", if (!is.null(o$ss)) "Log abundance (state-space estimate)" else "Log count at start of interval", opts),
                  y = fig_label("log_growth", "Annual log growth rate", opts))
  th <- theme_eco(opts$base_size %||% 9, opts$font %||% "", opts$style %||% "classic")
  figs$growth_checks <- with_caption(compose_figures(list(p1, p2), ncol = 2) & th,
    sprintf("Model checks. A: annual log growth rates (per year of interval) with the estimated mu; no trend or strong autocorrelation expected (lag-1 r = %s). B: growth rate against log abundance; a clear negative slope (here %s, p = %s) suggests density dependence, which this model ignores.",
            pva_f(o$checks$lag1, 2), pva_f(o$checks$dd_slope, 3), pva_f(o$checks$dd_p, 3)))
  attr(figs$growth_checks, "fig_size") <- list(width = "double", aspect = 0.4)
  figs
}

matrix_figures <- function(fit, data, opts = fig_opts()) {
  o <- fit$object; u <- data$meta$units; e <- fit$estimates
  figs <- list()
  figs$quasi_extinction <- pva_qext_fig(o$cdf, opts, u, o$threshold,
    itype = "interval from parameter uncertainty (vital-rate means and environmental SDs drawn from their sampling distributions)")
  qN <- o$base$point$qN
  pr <- data.frame(time = 0:o$tmax, estimate = qN[, 3], lower = qN[, 1], upper = qN[, 5], lower50 = qN[, 2], upper50 = qN[, 4])
  sm <- o$base$point$samp[, 1:min(15, ncol(o$base$point$samp)), drop = FALSE]
  sl <- data.frame(time = rep(0:o$tmax, ncol(sm)), n = as.vector(sm), id = rep(seq_len(ncol(sm)), each = nrow(sm)))
  p <- ggplot2::ggplot(pr, ggplot2::aes(time, estimate)) +
    ggplot2::geom_line(data = sl, ggplot2::aes(time, n, group = id), colour = eco_pal$muted, alpha = 0.35, linewidth = 0.25) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = eco_pal$band, alpha = 0.15) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25) +
    ggplot2::geom_line(colour = eco_pal$band, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = o$threshold, linetype = 2, colour = eco_pal$truth, linewidth = 0.4) +
    ggplot2::annotate("text", x = 0, y = o$threshold, label = "quasi-extinction threshold", hjust = 0, vjust = -0.5, size = 2.4, colour = eco_pal$truth) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = fig_label("years_ahead", "Years from now", opts), y = fig_label("abundance", sprintf("Abundance (%s)", u), opts))
  if (is.finite(data$ceiling %||% Inf)) p <- p + ggplot2::geom_hline(yintercept = data$ceiling, linetype = 3, colour = eco_pal$null, linewidth = 0.4)
  figs$trajectories <- with_caption(p, sprintf("Projected abundance at the estimated vital rates: median (line), 50%% and 95%% of %d simulated trajectories (bands), and 15 example trajectories (grey). Dashed: quasi-extinction threshold%s.",
                                               o$settings$npt, if (is.finite(data$ceiling %||% Inf)) "; dotted: ceiling" else ""))
  el <- e[startsWith(e$param, "elasticity_"), ]
  el$param <- sub("^elasticity_", "", el$param)
  # readable labels: the rate's description (text before any parenthesis or comma), else a default name
  dsc <- if (!is.null(o$rates$description)) stats::setNames(as.character(o$rates$description), o$rates$rate) else character()
  defaults <- c(s_juv = "First-year survival", s_sub = "Second-year survival", s_ad = "Adult survival", fec = "Fecundity")
  el$param <- vapply(el$param, function(r) {
    d <- if (r %in% names(dsc) && !is.na(dsc[[r]]) && nzchar(dsc[[r]])) sub("\\s*[(,].*$", "", dsc[[r]]) else if (r %in% names(defaults)) defaults[[r]] else r
    fig_label(r, pva_cap(d), opts)
  }, "")
  pe <- fig_coef(el, null_value = 0, xlab = fig_label("elasticity", "Elasticity of lambda", opts), opts = opts)
  # long vital-rate descriptions need the 1.5-column width or the x axis is squeezed / clipped
  attr(pe, "fig_size") <- utils::modifyList(attr(pe, "fig_size") %||% list(), list(width = "onehalf"))
  figs$elasticity <- with_caption(pe,
    "Vital-rate elasticities (proportional change in lambda per proportional change in the rate) of the deterministic growth rate at the estimated rates, with 95% intervals from parameter uncertainty. Values sum to more than 1 when a rate appears in several matrix entries.")
  em <- o$elast_matrix; st_lab <- vapply(o$mx$stages, function(s) fig_label(s, s, opts), "")
  dimnames(em) <- list(st_lab, st_lab)
  dm <- as.data.frame(as.table(em)); names(dm) <- c("row", "col", "value")
  pm <- fig_matrix(em, "Elasticity", diverging = FALSE, label = FALSE) +
    ggplot2::geom_text(data = dm, ggplot2::aes(col, row, label = formatC(value, digits = 2, format = "f")),
                       inherit.aes = FALSE, size = 2.4, colour = ifelse(dm$value > 0.6 * max(dm$value), eco_pal$ink, "white")) +
    ggplot2::labs(x = fig_label("stage_from", "Stage this year", opts), y = fig_label("stage_to", "Stage next year", opts))
  figs$elasticity_matrix <- with_caption(pm,
    "Elasticities of lambda to each matrix entry (they sum to 1). Rows: stage in the next year; columns: stage this year.")
  # Scenarios: risk at each horizon with intervals
  if (length(o$scenarios)) {
    sc <- c("Baseline", names(o$scenarios))
    d <- do.call(rbind, lapply(o$horizons, function(h) {
      b <- pva_row(e, paste0("p_qext_", h, "yr"))
      rows <- data.frame(scenario = "Baseline", horizon = h, estimate = b$estimate, lower = b$lower, upper = b$upper)
      for (nm in names(o$scenarios)) {
        r <- pva_row(e, paste0("scenario_", gsub("[^a-z0-9]+", "_", tolower(nm)), "_p_qext_", h, "yr"))
        rows <- rbind(rows, data.frame(scenario = nm, horizon = h, estimate = r$estimate, lower = r$lower, upper = r$upper))
      }
      rows
    }))
    d$scenario <- factor(d$scenario, levels = rev(sc))
    d$horizon <- factor(paste(d$horizon, "years"), levels = paste(o$horizons, "years"))
    p <- ggplot2::ggplot(d, ggplot2::aes(estimate, scenario)) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), colour = eco_pal$ink, linewidth = 0.5) +
      ggplot2::geom_point(size = 1.8, colour = eco_pal$ink) +
      ggplot2::facet_wrap(~horizon, nrow = 1) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      ggplot2::theme(panel.spacing.x = ggplot2::unit(5, "mm")) +
      ggplot2::labs(x = fig_label("p_qext", "Quasi-extinction probability", opts), y = NULL)
    attr(p, "fig_size") <- list(width = "double", height = 25 + 9 * length(sc))
    figs$scenarios <- with_caption(p, "Quasi-extinction probability by management scenario and time horizon (point: estimate; line: 95% interval from parameter uncertainty). Scenarios share random numbers with the baseline, so differences are paired; see the estimates table for the change and its interval.")
  }
  ld <- data.frame(v = c(exp(o$loglams_draws), o$lam_det_draws),
                   what = rep(c("Stochastic", "Deterministic"), c(length(o$loglams_draws), length(o$lam_det_draws))))
  dd <- lapply(split(ld$v, ld$what), stats::density)
  ymax <- max(vapply(dd, function(z) max(z$y), 0))
  lab <- data.frame(what = c("Deterministic (mean matrix)", "Stochastic (simulated)"), x = min(ld$v), y = ymax * c(1.12, 1.04))
  ld$what <- ifelse(ld$what == "Stochastic", "Stochastic (simulated)", "Deterministic (mean matrix)")
  p <- ggplot2::ggplot(ld, ggplot2::aes(v, colour = what, fill = what, linetype = what)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.5, show.legend = FALSE) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x, y, label = what, colour = what), inherit.aes = FALSE, size = 2.4, hjust = 0, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = eco_pal$null) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(4)) +
    ggplot2::labs(x = fig_label("lambda", "Population growth rate, lambda", opts), y = "Density across parameter draws")
  figs$lambda_distribution <- with_caption(p, sprintf("Distribution of the deterministic (mean-matrix) and stochastic (simulated long-run) growth rates across %d vital-rate sets drawn from their sampling uncertainty. Dashed: lambda = 1. P(stochastic lambda < 1) = %s.",
                                                      length(o$loglams_draws), pva_pct(mean(o$loglams_draws < 0))))
  figs
}
