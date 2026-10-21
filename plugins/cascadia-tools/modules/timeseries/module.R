# ecotools module `timeseries`: population trends and dynamics from annual
# counts or indices, separating real population change (process noise) from
# counting error (observation error).
#
# Models (settings$model; design$model in simulation):
#   exponential  EGSS: exponential growth with drift + process + observation
#                error, Kalman-filter ML (Dennis et al. 2006; Humbert et al. 2009)
#   gompertz     GSS: stochastic Gompertz density dependence + observation
#                error, Kalman-filter ML (Dennis et al. 2006; Knape 2008)
#   marss        several sub-populations, shared or separate trends (MARSS;
#                Holmes et al. 2012)
#   bayes        Nimble state-space model: covariates on growth rate, optional
#                density dependence, missing years (Kery & Schaub 2012 ch. 5)
#   loglinear    log-linear trend GLMM with AR(1) residuals (glmmTMB) for
#                comparison with conventional trend analyses
#
# Data structure: see R/common.R.

for (f in c("common.R", "kalman.R", "reml.R", "uni.R", "marss.R", "bayes.R", "loglinear.R", "figures.R"))
  sys.source(file.path(MODULE_DIR, "R", f), envir = environment())

ts_simulate <- function(params, design, seed) {
  switch(design$model %||% "exponential",
         exponential = uni_simulate(params, design, seed, "exponential"),
         gompertz = uni_simulate(params, design, seed, "gompertz"),
         marss = marss_simulate(params, design, seed),
         bayes = bayes_simulate(params, design, seed),
         loglinear = ll_simulate(params, design, seed),
         stop("design$model must be one of ", paste(TS_MODELS, collapse = ", ")))
}

ts_fit <- function(data, settings) {
  st <- merge_list(TS_SETTINGS, settings)
  m <- ts_model(data, st)
  f <- switch(m,
              exponential = uni_fit(data, st, "exponential"),
              gompertz = uni_fit(data, st, "gompertz"),
              marss = marss_fit(data, st),
              bayes = bayes_fit(data, st),
              loglinear = ll_fit(data, st))
  f$object$settings <- st
  f
}

TS_SETTINGS <- list(model = NULL, estimation = "reml", observation_error = TRUE, zero_offset = 0, n_draws = 1000, seed = 1,
                    trend = "shared", states = "separate", process_cov = "diagonal and equal",
                    density_dependence = NULL, prior_sd_process = c(0, 1), prior_sd_observation = c(0, 1),
                    niter = NULL, nburnin = NULL, thin = NULL, nchains = NULL, fast = FALSE,
                    family = NULL, ar1 = TRUE, dd_test_boot = 199)

# ---- interpretation ---------------------------------------------------------

ts_ci <- function(r, d = 1, pct = FALSE) {
  f <- if (pct) function(v) sprintf(paste0("%+.", d, "f%%"), v) else function(v) ts_fmt(v, d)
  sprintf("%s (95%% %s %s to %s)", f(r$estimate), "CI", f(r$lower), f(r$upper))
}

ts_interpret <- function(fit, data) {
  obj <- fit$object; est <- fit$estimates; out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the fit did not converge or a variance is at its boundary (see Diagnostics). Do not interpret the numbers below without fixing this.")
  cl <- tolower(data$meta$count_label %||% "count")
  ci_word <- if (obj$model == "bayes") "credible interval" else "confidence interval"
  yrs <- range(obj$trajectory$year)
  # trend(s)
  tr <- est[grepl("^trend_pct_per_year", est$param), ]
  for (i in seq_len(nrow(tr))) {
    r <- tr[i, ]
    who <- if (r$param == "trend_pct_per_year") "The population" else paste0("Site ", sub("trend_pct_per_year_", "", r$param))
    dir <- if (r$lower > 0) "increased by" else if (r$upper < 0) "declined by" else "showed no clear trend: its estimated change was"
    amt <- if (r$lower > 0 || r$upper < 0) sprintf("%.1f%%", abs(r$estimate)) else sprintf("%+.1f%%", r$estimate)
    out <- c(out, sprintf("%s %s %s per year on average over %d-%d (95%% %s %s to %s). Over the %d years this compounds to a %+.0f%% change.",
                          who, dir, amt, yrs[1], yrs[2], ci_word,
                          sprintf("%+.1f%%", r$lower), sprintf("%+.1f%%", r$upper), diff(yrs),
                          100 * ((1 + r$estimate / 100)^diff(yrs) - 1)))
    if (r$lower <= 0 && r$upper >= 0)
      out <- c(out, sprintf("%s interval runs from %+.1f%% to %+.1f%% per year, i.e. compounded over %d years from a %+.0f%% to a %+.0f%% change. \"No clear trend\" therefore means the counts cannot tell these apart, not that the population was stable. How large a decline this design could detect is a power question (core/run.R power; VALIDATION.md: the default REML intervals are conservative, especially when counting error is large relative to real year-to-year change).",
                            if (r$param == "trend_pct_per_year") "The trend" else paste0("Site ", sub("trend_pct_per_year_", "", r$param), "'s trend"),
                            r$lower, r$upper, diff(yrs), 100 * ((1 + r$lower / 100)^diff(yrs) - 1), 100 * ((1 + r$upper / 100)^diff(yrs) - 1)))
  }
  rt <- ts_est(fit, "realised_trend_pct_per_year")
  if (nrow(rt)) out <- c(out, sprintf("The estimated true population changed by %+.1f%% per year on average between %d and %d (95%% credible interval %+.1f%% to %+.1f%%). This realised trend describes what happened; with density dependence it is not a constant growth rate to project forward.",
                                      rt$estimate, yrs[1], yrs[2], rt$lower, rt$upper))
  # abundance
  nf <- ts_est(fit, "N_final")
  if (nrow(nf)) out <- c(out, sprintf("Estimated true abundance (on the scale of the %s) in %d: %s (95%% %s %s to %s). If counts miss a constant fraction of animals this is an index of abundance, not a census; trends are unaffected by a constant detection rate.",
                                      cl, yrs[2], ts_sig(nf$estimate), ci_word, ts_sig(nf$lower), ts_sig(nf$upper)))
  # variances
  sp <- ts_est(fit, "sd_process"); so <- ts_est(fit, "sd_observation")
  if (nrow(sp) == 1 && nrow(so) == 1) {
    out <- c(out, sprintf("Real year-to-year variation in growth (process SD) is %.2f (95%% %s %.2f-%.2f) on the log scale, i.e. growth in a typical year departs from the average by about %.0f%%; counting error (observation SD) is %.2f (%.2f-%.2f), i.e. a single count is typically off by about %.0f%%. %s",
                          sp$estimate, ci_word, sp$lower, sp$upper, 100 * (exp(sp$estimate) - 1), so$estimate, so$lower, so$upper,
                          100 * (exp(so$estimate) - 1),
                          if (so$estimate > sp$estimate) "More of the bounce in the raw counts is attributed to counting error than to real change, so single year-to-year changes in raw counts should not be read as real population changes." else "More of the bounce in the raw counts is attributed to real population change than to counting error."))
    zero <- c(if (sp$estimate < 1e-3) "process", if (so$estimate < 1e-3) "observation")
    if (length(zero)) out <- c(out, sprintf("The %s SD was estimated at zero (a boundary). This does not mean there is none: the series is too short or too smooth for the model to split the noise between real change and counting error. %s",
                                            paste(zero, collapse = " and "),
                                            if (identical(obj$method, "REML")) "The trend and abundance intervals reported here average over every split of the noise the data allow, so they remain valid (VALIDATION.md); only the SD point estimates are unreliable."
                                            else if (obj$model == "exponential") "Abundance intervals average over the possible splits, but the profile interval for the trend under ML is then too narrow (simulation coverage 0.68 at 25 years; VALIDATION.md): use the default estimation: reml."
                                            else "Abundance intervals average over the possible splits of the noise; the SD and density-dependence estimates are unreliable (VALIDATION.md)."))
  }
  else if (nrow(sp) == 1 && isFALSE(obj$obs_error))
    out <- c(out, sprintf("Process SD %.2f: this model assumes counts have NO observation error, so counting noise is attributed to the population itself.", sp$estimate))
  # density dependence
  b <- ts_est(fit, "density_dependence_b")
  if (nrow(b)) {
    hl <- if (b$estimate < 0 && b$estimate > -1) log(0.5) / log(1 + b$estimate) else NA
    K <- ts_est(fit, "equilibrium_N")
    at_lo <- obj$model == "gompertz" && b$lower <= -1.98
    out <- c(out, sprintf("Density dependence b = %.2f (95%% %s %.2f to %.2f; 0 = none, -1 = complete return to equilibrium each year).%s%s%s",
                          b$estimate, ci_word, b$lower, b$upper,
                          if (at_lo) " The lower limit is the edge of the range searched (b = -2), i.e. the data do not bound b from below." else "",
                          if (is.finite(hl) && b$upper < 0) sprintf(" Taken at face value, a population displaced from equilibrium would close half the (log-scale) gap in about %.1f years; because b is overestimated from single count series (see below), the true return is likely slower.", hl)
                          else if (is.finite(hl)) " The interval for b includes 0 (no density dependence), so no return time to equilibrium can be given." else "",
                          if (nrow(K) && is.finite(K$upper) && K$lower > 0 && K$upper / K$lower < 20)
                            sprintf(" Equilibrium abundance (where growth averages zero%s): %s (95%% %s %s to %s).",
                                    if (length(obj$covariates)) ", at average covariates" else "",
                                    ts_sig(K$estimate), ci_word, ts_sig(K$lower), ts_sig(K$upper))
                          else if (nrow(K)) " Equilibrium abundance cannot be estimated from these data: its interval spans more than a 20-fold range, because the data are consistent with almost no density dependence, and then the equilibrium is undefined." else ""))
    if (obj$model == "bayes")
      out <- c(out, "Because the prior allows only b between -2 and 0, the credible interval for b can never include 0, so it is not evidence that density dependence exists. In simulations with true b = -0.3 and 25 annual counts this model's posterior median for b averaged -0.49, i.e. it overstates density dependence (VALIDATION.md). To test for density dependence use model: gompertz, which reports a bootstrap likelihood-ratio test.")
    out <- c(out, sprintf("Caution: density dependence estimated from a single count series is weakly identified when observation error is also estimated (Knape 2008; Dennis et al. 2006): the data often cannot tell strong density dependence plus little counting error from weak density dependence plus more counting error. Our simulation checks (VALIDATION.md) show how reliable b is at this series length; an independent estimate of counting error (replicate counts, sightability trials) sharpens it considerably; supply it as settings: prior_sd_observation%s. Ignoring observation error altogether makes density dependence look stronger than it is.",
      if (obj$model == "bayes") "" else " with model: bayes (this maximum-likelihood model cannot use a prior)"))
  }
  # covariates
  for (cv in obj$covariates %||% character()) {
    lab <- data$meta$covariate_labels[[cv]] %||% gsub("_", " ", cv)
    r <- ts_est(fit, paste0("growth_effect_pct_", cv)); bb <- ts_est(fit, paste0("growth_beta_", cv))
    if (nrow(r)) {
      sc <- obj$scaling[[cv]]
      pos <- mean(unlist(lapply(obj$derived, function(d) d[, paste0("growth_beta_", cv)])) > 0)
      q <- stats::quantile(obj$cov_raw[-1, cv], c(0.1, 0.9), na.rm = TRUE, names = FALSE)
      g <- 100 * (exp(stats::median(obj$pdraws$r0) + bb$estimate * (q - sc$center) / sc$scale) - 1)
      out <- c(out, sprintf("An increase of one standard deviation (%.3g) in %s is associated with a %+.1f%% change in the annual growth multiplier (95%% credible interval %+.1f%% to %+.1f%%; P(effect > 0) = %.2f). At average abundance the population is expected to change by %+.1f%% in a year with %s at its 10th percentile (%.3g) vs %+.1f%% at its 90th percentile (%.3g). This is an association, not proof that %s causes the change.",
                            sc$scale, lab, r$estimate, r$lower, r$upper, pos,
                            g[1], lab, q[1], g[2], q[2], lab))
    }
    r <- ts_est(fit, paste0("count_effect_pct_", cv))
    if (nrow(r)) out <- c(out, sprintf("An increase of one standard deviation (%.3g) in %s is associated with a %+.1f%% change in the expected %s in that year (95%% CI %+.1f%% to %+.1f%%). In this log-linear model the covariate shifts the level of counts in a year; it is not an effect on population growth.",
                                       obj$scaling[[cv]]$scale, lab, r$estimate, cl, r$lower, r$upper))
  }
  g <- obj$gof
  if (!is.null(g) && is.finite(g$ljung_box_p)) {
    bad <- c(if (g$ljung_box_p < 0.05) sprintf("autocorrelated (Ljung-Box p = %.3f: a pattern the model misses, e.g. a change in trend, cycles or density dependence)", g$ljung_box_p),
             if (is.finite(g$shapiro_p) && g$shapiro_p < 0.05) sprintf("not normal (Shapiro-Wilk p = %.3f: outlying counts or a crash/boom year)", g$shapiro_p))
    out <- c(out, if (length(bad)) sprintf("Model check: the standardised prediction errors are %s. See figure residuals and treat the estimates with caution.", paste(bad, collapse = " and "))
                  else sprintf("Model check: the standardised one-step-ahead prediction errors show no autocorrelation (Ljung-Box p = %.2f) and no departure from normality (Shapiro-Wilk p = %.2f), so there is no sign of misfit. With short series these tests have little power.",
                               g$ljung_box_p, g$shapiro_p))
  }
  pp <- obj$ppc
  if (!is.null(pp)) {
    ext <- pp$p < 0.05 | pp$p > 0.95
    out <- c(out, sprintf("Posterior predictive check: the variability (Bayesian p = %.2f) and lag-1 autocorrelation (p = %.2f) of year-to-year changes in the replicated counts %s.",
                          pp$p[[1]], pp$p[[2]], if (any(ext, na.rm = TRUE)) "differ from the observed counts: the model misrepresents the noise. See figure ppc" else "match the observed counts (p between 0.05 and 0.95), so there is no sign of misfit"))
  }
  rho <- ts_est(fit, "ar1_rho")
  if (obj$model == "loglinear" && nrow(rho) && is.finite(rho$estimate) && rho$estimate > 0.9)
    out <- c(out, sprintf("Warning: the year-to-year deviations are almost a random walk (AR(1) correlation %.2f). The counts then behave like a population whose real changes accumulate, and the straight-line trend interval of this regression is likely too narrow. Use model: marss (several sites) or exponential (one series) for the trend.", rho$estimate))
  if (obj$model == "loglinear")
    out <- c(out, "The log-linear regression treats year-to-year deviations as autocorrelated noise around a straight trend; it does not separate real population change from counting error. If the population itself wanders (real year-to-year changes accumulate, as in the state-space models), its trend interval is too narrow: in simulations of 4 sites x 25 years with process and counting noise the 95% interval covered the true trend in only 82.5% of series (VALIDATION.md, independent verification). Compare with model: exponential (single series) or marss (several series) when the trend is the key result.")
  if (obj$model == "marss" && identical(obj$trend, "shared"))
    out <- c(out, "All sites were assumed to share one long-term growth rate (settings: trend: shared). Fit trend: separate as well and compare the per-site trends before concluding that sites behave alike.")
  tst <- obj$kf$dd_test
  if (!is.null(tst) && is.finite(tst$p))
    out <- c(out, sprintf("Formal test for density dependence (parametric-bootstrap likelihood ratio, %d simulated series): p = %.3f. %s Note that the interval for b is NOT a valid test of b = 0: in simulations with no density dependence it excluded 0 in about half of 30-year series, whereas this test held its 5%% error rate (VALIDATION.md).",
                          tst$B, tst$p, if (tst$p < 0.05) "The counts are more consistent with density-dependent growth than with density-independent growth with drift." else "The counts are consistent with density-independent growth (exponential growth with drift); this is not evidence that density dependence is absent, only that this series cannot show it."))
  out
}

ts_methods <- function(fit, data, settings) {
  obj <- fit$object; yrs <- range(obj$trajectory$year)
  base <- switch(obj$model,
    exponential = if (identical(obj$method, "REML"))
      "We fitted an exponential-growth state-space model (EGSS; Dennis et al. 2006; Humbert et al. 2009) to log-transformed counts: x_t = x_{t-1} + mu + e_t with process noise e_t ~ N(0, sigma_p^2), and log count_t = x_t + f_t with observation error f_t ~ N(0, sigma_o^2). The two variances were estimated by restricted maximum likelihood (REML) from the first differences of the log counts, whose covariance is tridiagonal (Staples et al. 2004; differences spanning unsurveyed years have variance k sigma_p^2 + 2 sigma_o^2); profile-likelihood intervals were used for the standard deviations. Uncertainty about the two variances was then integrated out: the restricted likelihood was evaluated on a 60 x 60 grid of the two SDs (flat prior), and mu given each pair of SDs is normal with its generalised-least-squares mean and variance, so the reported trend is the mean and its 95% interval the 2.5% and 97.5% quantiles of this mixture of normals. The true trajectory was estimated with a diffuse-start Kalman smoother, propagating parameter uncertainty with 1000 draws. Trend = 100(exp(mu) - 1)% per year."
    else "We fitted an exponential-growth state-space model (EGSS; Dennis et al. 2006; Humbert et al. 2009) to log-transformed counts: x_t = x_{t-1} + mu + e_t with process noise e_t ~ N(0, sigma_p^2), and log count_t = x_t + f_t with observation error f_t ~ N(0, sigma_o^2), by exact maximum likelihood using the Kalman filter (initial state one year before the first count estimated; years without counts skipped), with profile-likelihood intervals. Abundance intervals come from the Kalman smoother with parameter uncertainty propagated by adaptive sampling-importance-resampling from the likelihood (flat priors on the parameters and on the SD scale). Trend = 100(exp(mu) - 1)% per year.",
    gompertz = paste0("We fitted a stochastic Gompertz state-space model (GSS; Dennis et al. 2006): x_t = x_{t-1} + b (x_{t-1} - log K) + e_t, e_t ~ N(0, sigma_p^2), with log count_t = x_t + f_t, f_t ~ N(0, sigma_o^2), where b measures density dependence (b = 0: none) and K is the equilibrium abundance. Parameters were estimated by exact maximum likelihood using the Kalman filter (initial state one year before the first count estimated); intervals for the parameters are profile-likelihood intervals; abundance intervals come from the Kalman smoother with parameter uncertainty propagated by adaptive sampling-importance-resampling from the likelihood (flat priors). Density dependence in this model is known to be weakly identifiable when observation error is estimated (Knape 2008; Auger-Methe et al. 2016).",
                      if (!is.null(obj$kf$dd_test)) sprintf(" Density dependence was tested with a parametric-bootstrap likelihood-ratio test of the exponential-growth state-space model (H0) against the Gompertz model (%d bootstrap series simulated from the fitted H0 model; Dennis & Taper 1994; Dennis et al. 2006).", obj$kf$dd_test$B) else ""),
    marss = if (identical(obj$method, "REML")) sprintf("We fitted a multivariate state-space model to log counts from %d sites (one state per site; MARSS form with B = I, Z = I, Q and R diagonal and equal; Holmes et al. 2012): x_{s,t} = x_{s,t-1} + %s + w_{s,t}, w ~ N(0, sigma_p^2); log count_{s,t} = x_{s,t} + v_{s,t}, v ~ N(0, sigma_o^2). Because sites are independent given the two variances, the variances were estimated by restricted maximum likelihood pooled over the first differences of each site's log counts (Staples et al. 2004; Humbert et al. 2009), with profile intervals; trend intervals integrate the restricted likelihood over the two variances (flat prior). True trajectories: diffuse-start Kalman smoother per site with parameter uncertainty propagated by 1000 draws. Trend = 100(exp(u) - 1)%% per year.",
                    length(obj$sites), if (identical(obj$trend, "shared")) "u (shared by all sites)" else "u_s (one per site)")
    else sprintf("We fitted a multivariate autoregressive state-space model (MARSS; Holmes et al. 2012) with %s to log counts from %d sites: x_t = x_{t-1} + u + w_t, w_t ~ MVN(0, Q) with Q '%s'; y_t = Z x_t + a + v_t, v_t ~ MVN(0, R) with R diagonal and equal. Parameters were estimated by maximum likelihood (MARSS package; BFGS, or EM for equalvarcov) with Hessian-based intervals.",
                    if (identical(obj$states, "single")) "all sites treated as replicate surveys of one population (Z = column of ones)" else "one state per site",
                    length(obj$sites), obj$process_cov),
    bayes = sprintf("We fitted a Bayesian state-space population model (Kery & Schaub 2012) in NIMBLE (de Valpine et al. 2017): x_t = x_{t-1} + r0 %s%s + e_t, e_t ~ N(0, sigma_p^2), log count_t ~ N(x_t, sigma_o^2) for surveyed years. %sPriors: r0 ~ N(0, 1)%s%s, sigma_p ~ half-N(%g, %g), sigma_o ~ half-N(%g, %g). MCMC: %d chains of %d iterations (%d burn-in, thinning %d); convergence assessed with R-hat and bulk ESS (posterior package).",
                    if (isTRUE(obj$dd)) "+ b (x_{t-1} - mean log count)" else "",
                    if (length(obj$covariates)) " + sum_k beta_k z_{k,t} (standardised covariates acting on growth into year t)" else "",
                    "Missing years were handled as unobserved latent states. ",
                    if (isTRUE(obj$dd)) ", b ~ U(-2, 0)" else "", if (length(obj$covariates)) ", beta_k ~ N(0, 1)" else "",
                    obj$priors$sd_process[1], obj$priors$sd_process[2], obj$priors$sd_observation[1], obj$priors$sd_observation[2],
                    obj$mcmc$nchains, obj$mcmc$niter, obj$mcmc$nburnin, obj$mcmc$thin),
    loglinear = sprintf("We fitted a log-linear trend model with %s family and AR(1) year effects within %s using glmmTMB (Brooks et al. 2017): log E[count] = alpha%s + mu (year - mean year)%s + AR(1) deviation. Trend = 100(exp(mu) - 1)%% per year with Wald intervals.",
                        obj$family, if (length(obj$sites) > 1) "each site, with random site intercepts" else "the series",
                        if (length(obj$sites) > 1) " + alpha_site" else "", if (length(obj$covariates)) " + covariates" else ""))
  paste(base, sprintf("Data: %d-%d, %d site(s). Simulation-based validation of the method at comparable designs is reported in ecotools modules/timeseries/VALIDATION.md.",
                      yrs[1], yrs[2], length(obj$sites)))
}

ts_plots <- function(fit, data) {
  f <- ts_figures(fit, data, fig_opts())
  f[setdiff(names(f), c("trace"))]
}

module <- ecomodule(
  id = "timeseries",
  title = "Population time series: trends, density dependence and drivers of growth",
  description = paste(
    "Is the population increasing or declining, by how much per year, how large is it really once counting error is removed,",
    "is growth density dependent, and do environmental covariates (e.g. winter severity) explain good and bad years?",
    "State-space models separate real population change (process noise) from counting error (observation error)."),
  defaults = list(
    params = list(growth_rate_mu = -0.03, sd_process = 0.1, sd_observation = 0.15,
                  density_dependence_b = -0.3, equilibrium_N = 500, growth_beta = -0.15,
                  count_beta = 0.15, sd_site = 0.5, ar1_rho = 0.5, sd_ar1 = 0.2, mean_count = 40, nb_size = 10,
                  process_corr = 0),
    design = list(model = "exponential", n_years = 25, start_N = 500, frac_missing = 0, n_sites = 4,
                  density_dependence = TRUE, n_covariates = 1, first_year = 2000)),
  simulate = ts_simulate, fit = ts_fit, read_data = ts_read_data,
  check_data = ts_check_data, interpret = ts_interpret, figures = ts_figures, plots = ts_plots,
  methods_text = ts_methods,
  settings = list(),
  fast_settings = list(fast = TRUE, n_draws = 300, dd_test_boot = 0),
  # density_dependence_b is deliberately not a power parameter: its interval
  # is not a valid test of b = 0 (ML: 47% false positives; bayes: the U(-2, 0)
  # prior never lets the interval reach 0). The valid test is the bootstrap
  # LRT reported for model: gompertz (VALIDATION.md).
  power_params = c("growth_rate_mu", "growth_beta_winter_severity", "count_beta_spring_rain_mm"),
  engine = "Kalman-filter ML (hand-written), MARSS, nimble, glmmTMB",
  references = c(
    "Dennis B, Ponciano JM, Lele SR, Taper ML, Staples DF (2006) Estimating density dependence, process noise, and observation error. Ecological Monographs 76:323-341.",
    "Dennis B, Taper ML (1994) Density dependence in time series observations of natural populations: estimation and testing. Ecological Monographs 64:205-224.",
    "Shenk TM, White GC, Burnham KP (1998) Sampling-variance effects on detecting density dependence from temporal trends in natural populations. Ecological Monographs 68:445-463.",
    "Knape J (2008) Estimability of density dependence in models of time series data. Ecology 89:2994-3000.",
    "Auger-Methe M, Field C, Albertsen CM, Derocher AE, Lewis MA, Jonsen ID, Mills Flemming J (2016) State-space models' dirty little secrets: even simple linear Gaussian models can have estimation problems. Scientific Reports 6:26677.",
    "Humbert J-Y, Mills LS, Horne JS, Dennis B (2009) A better way to estimate population trends. Oikos 118:1940-1946.",
    "Staples DF, Taper ML, Dennis B (2004) Estimating population trend and process variation for PVA in the presence of sampling error. Ecology 85:923-929.",
    "Holmes EE, Ward EJ, Wills K (2012) MARSS: multivariate autoregressive state-space models for analyzing time-series data. The R Journal 4:11-19.",
    "Holmes EE, Scheuerell MD, Ward EJ (2024) Analysis of multivariate time series using the MARSS package. NOAA Fisheries, Seattle (MARSS User Guide).",
    "Kery M, Schaub M (2012) Bayesian Population Analysis using WinBUGS. Academic Press, chapter 5 (state-space models).",
    "Auger-Methe M, Newman K, Cole D, et al. (2021) A guide to state-space modeling of ecological time series. Ecological Monographs 91:e01470.",
    "Brooks ME, Kristensen K, van Benthem KJ, et al. (2017) glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. The R Journal 9:378-400.",
    "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26:403-413.",
    "Harvey AC (1989) Forecasting, Structural Time Series Models and the Kalman Filter. Cambridge University Press.")
)
