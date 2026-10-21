# ipm: interpretation, methods text, references and publication figures.

ipm_references <- c(
  "Besbeas P, Freeman SN, Morgan BJT, Catchpole EA (2002) Integrating mark-recapture-recovery and census data to estimate animal abundance and demographic parameters. Biometrics 58: 540-547.",
  "Schaub M, Abadi F (2011) Integrated population models: a novel analysis framework for deeper insights into population dynamics. Journal of Ornithology 152 (Suppl 1): S227-S237.",
  "Schaub M, Kery M (2022) Integrated Population Models: Theory and Ecological Applications with R and JAGS. Academic Press, London. (chapters 5, 6, 9 and 11; code: github.com/mikemeredith/IPM_code; R package IPMbook)",
  "Abadi F, Gimenez O, Ullrich B, Arlettaz R, Schaub M (2010) Estimation of immigration rate using integrated population models. Journal of Applied Ecology 47: 393-400.",
  "Schaub M, Reichlin TS, Abadi F, Kery M, Jenni L, Arlettaz R (2012) The demographic drivers of local population dynamics in two rare migratory birds. Oecologia 168: 97-108.",
  "Koons DN, Iles DT, Schaub M, Caswell H (2016) A life history perspective on the demographic drivers of structured population dynamics in changing environments. Ecology Letters 19: 1023-1031.",
  "Koons DN, Arnold TW, Schaub M (2017) Understanding the demographic drivers of realized population growth rates. Ecological Applications 27: 2102-2115.",
  "Burnham KP, Anderson DR, White GC, Brownie C, Pollock KH (1987) Design and analysis methods for fish survival experiments based on release-recapture. American Fisheries Society Monograph 5.",
  "de Valpine P, Turek D, Paciorek CJ, Anderson-Bergman C, Temple Lang D, Bodik R (2017) Programming with models: writing statistical algorithms for general model structures with NIMBLE. Journal of Computational and Graphical Statistics 26: 403-413.",
  "Abadi F, Gimenez O, Arlettaz R, Schaub M (2010) An assessment of integrated population models: bias, accuracy, and violation of the assumption of independence. Ecology 91: 7-14.",
  "Schaub M, Fletcher D (2015) Estimating immigration using a Bayesian integrated population model: choice of parametrization and priors. Environmental and Ecological Statistics 22: 535-549.",
  "Riecke TV, Williams PJ, Behnke TL, Gibson D, Leach AG, Sedinger BS, Street PA, Sedinger JS (2019) Integrated population models: model assumptions and inference. Methods in Ecology and Evolution 10.",
  "Paquet M, Knape J, Arlt D, et al. (2021) Integrated population models poorly estimate the demographic contribution of immigration. Methods in Ecology and Evolution 12: 1899-1910. doi:10.1111/2041-210X.13667",
  "Schaub M, Ullrich B (2021) A drop in immigration results in the extinction of a local woodchat shrike population. Animal Conservation 24: 335-345.")

# covariate label for text: the config figure label if given (figures: labels:), else the column name
ipm_cov_label <- function(cv, data) {
  l <- data$figure_labels[[cv]]
  if (is.null(l)) cv else as.character(l)
}

ipm_row <- function(est, p) est[est$param == p, , drop = FALSE]
ipm_fmt <- function(r, d = 2, pct = FALSE) {
  if (!nrow(r)) return("NA")
  if (pct) return(sprintf("%.0f%% (95%% CrI %.0f-%.0f%%)", 100 * r$estimate, 100 * r$lower, 100 * r$upper))
  sprintf(paste0("%.", d, "f (95%% CrI %.", d, "f-%.", d, "f)"), r$estimate, r$lower, r$upper)
}

# quantities computed from annual draws (immigrant share, local growth)
ipm_extra_draws <- function(fit) {
  o <- fit$object; dr <- o$draws; yrs <- o$years; T <- length(yrs)
  out <- list()
  nt <- paste0("N_total[", yrs, "]")
  if (all(nt %in% colnames(dr)) && o$flags$IMM == 1) {
    im <- dr[, paste0("immigrants[", yrs[-1], "]"), drop = FALSE]
    out$immigrant_share <- rowMeans(im / pmax(dr[, nt[-1], drop = FALSE], 1e-9))
  }
  if (o$flags$IPM == 1) {
    rho <- o$settings$sex_ratio
    f1 <- dr[, if (o$flags$AGEF == 1) "fecundity_yearling" else "fecundity"]
    f2 <- if (o$flags$AGEF == 1) dr[, "fecundity_adult"] else f1
    N2 <- dr[, paste0("surviving_adults[", yrs, "]"), drop = FALSE]
    share2 <- rowMeans(N2 / pmax(dr[, nt, drop = FALSE], 1e-9))
    out$lambda_without_immigration <- dr[, "adult_survival"] + rho * dr[, "juvenile_survival"] * ((1 - share2) * f1 + share2 * f2)
  }
  out
}

ipm_q <- function(x) c(estimate = mean(x), lower = unname(stats::quantile(x, 0.025)), upper = unname(stats::quantile(x, 0.975)))

ipm_interpret <- function(fit, data) {
  est <- fit$estimates; o <- fit$object; fl <- o$flags
  sp <- data$species; unit <- data$count_unit; yrs <- data$years; T <- length(yrs)
  out <- character()
  if (!isTRUE(fit$diagnostics$converged))
    out <- c(out, "WARNING: the MCMC did not converge (R-hat >= 1.05 for at least one quantity). Do not interpret the numbers below; run longer.")
  if (fl$IPM == 1) {
    lg <- ipm_row(est, "lambda_geomean")
    out <- c(out, sprintf(
      "Population size. The model estimates %s %s in %d and %s in %d (counts were %s and %s). Over the %d years the population changed on average by a factor of %s per year, i.e. %+.1f%% per year (95%% CrI %+.1f%% to %+.1f%%).",
      ipm_fmt(ipm_row(est, "N_first"), 0), unit, yrs[1], ipm_fmt(ipm_row(est, "N_final"), 0), yrs[T],
      if (is.na(data$counts[1])) "not taken" else format(data$counts[1]), if (is.na(data$counts[T])) "not taken" else format(data$counts[T]), T, ipm_fmt(lg, 3),
      100 * (lg$estimate - 1), 100 * (lg$lower - 1), 100 * (lg$upper - 1)))
  }
  out <- c(out, sprintf(
    "Survival. Apparent survival from one breeding season to the next (at average conditions) was %s for first-year birds (from fledging to their first breeding season) and %s for adults. 'Apparent' means that animals which emigrate permanently are counted as dead, so first-year survival in particular is a local-recruitment probability rather than true survival.",
    ipm_fmt(ipm_row(est, "juvenile_survival"), 2), ipm_fmt(ipm_row(est, "adult_survival"), 2)))
  if (fl$AGEF == 1) {
    out <- c(out, sprintf("Fecundity. Females breeding for the first time (yearlings, and immigrants, which are assumed to breed like yearlings) raised %s fledglings (both sexes) per year; older females raised %s.",
                          ipm_fmt(ipm_row(est, "fecundity_yearling"), 2), ipm_fmt(ipm_row(est, "fecundity_adult"), 2)))
  } else {
    out <- c(out, sprintf("Fecundity. Each breeding female raised %s fledglings (both sexes) per year on average.", ipm_fmt(ipm_row(est, "fecundity"), 2)))
  }
  ex <- tryCatch(ipm_extra_draws(fit), error = function(e) list())
  if (fl$IMM == 1) {
    sh <- if (!is.null(ex$immigrant_share)) ipm_q(ex$immigrant_share) else NULL
    out <- c(out, sprintf(
      "Immigration. There are no direct data on immigration: it is estimated from the mismatch between the observed counts and the number of local recruits and surviving adults that the survival and fecundity data predict. The estimated immigration rate, averaged over the study years, is %s new female breeders per female already present, i.e. on average %s immigrant %s per year%s.%s",
      ipm_fmt(ipm_row(est, "immigration_rate"), 2), ipm_fmt(ipm_row(est, "immigrants_mean"), 1), unit,
      if (!is.null(sh)) sprintf(", making up %.0f%% (95%% CrI %.0f-%.0f%%) of the breeders each year", 100 * sh[1], 100 * sh[2], 100 * sh[3]) else "",
      " Immigration is the least certain part of the model: in simulations (VALIDATION.md) the study-period immigration rate was recovered well when about 12% of breeders were immigrants (bias -5%), but was overestimated by about 19% (mean number of immigrants by 14%) when immigration was low (about 8%), and its year-to-year SD was overestimated severalfold. Treat a low immigration estimate as an upper value."))
  }
  if (!is.null(ex$lambda_without_immigration) && fl$IMM == 1) {
    q <- ipm_q(ex$lambda_without_immigration)
    out <- c(out, sprintf(
      "Local recruitment alone. With the estimated survival and fecundity but no immigrants, the population would change by a factor of %.2f per year (95%% CrI %.2f-%.2f; %s). %s",
      q[1], q[2], q[3], if (q[3] < 1) "a decline: the local population behaves as a sink that is maintained by immigration"
                        else if (q[2] > 1) "growth even without immigrants"
                        else "the data cannot tell whether it would grow or decline",
      "This projection assumes the vital rates would stay at their average values."))
  }
  # covariate effects
  for (rate in names(data$effects)) for (cv in data$effects[[rate]]) {
    r <- ipm_row(est, paste0(rate, "_beta_", cv))
    if (!nrow(r)) next
    link <- if (rate %in% c("juvenile_survival", "adult_survival")) "logit" else "log"
    sc <- o$scaling[[cv]]
    pr <- mean(fit$object$draws[, paste0(rate, "_beta_", cv)] > 0)
    rng <- range(data$covariates[[cv]][if (rate == "fecundity") seq_len(T) else seq_len(T - 1)], na.rm = TRUE)
    base <- ipm_row(est, if (rate == "fecundity" && fl$AGEF == 1) "fecundity_adult" else if (rate == "fecundity") "fecundity" else if (rate == "immigration") "immigration_rate_median" else rate)
    zlo <- (rng[1] - sc$center) / sc$scale; zhi <- (rng[2] - sc$center) / sc$scale
    inv <- if (link == "logit") stats::plogis else exp
    lk <- if (link == "logit") stats::qlogis else log
    vlo <- inv(lk(base$estimate) + r$estimate * zlo); vhi <- inv(lk(base$estimate) + r$estimate * zhi)
    cvl <- ipm_cov_label(cv, data)
    eff <- if (link == "logit") sprintf("odds ratio %.2f (95%% CrI %.2f-%.2f) per 1 SD of %s (1 SD = %.3g in the units of the covariate file)", exp(r$estimate), exp(r$lower), exp(r$upper), cvl, sc$scale)
           else sprintf("rate ratio %.2f (95%% CrI %.2f-%.2f) per 1 SD of %s (1 SD = %.3g in the units of the covariate file)", exp(r$estimate), exp(r$lower), exp(r$upper), cvl, sc$scale)
    out <- c(out, sprintf(
      "Covariate. %s was associated with %s: %s; P(effect > 0) = %.2f. Across the observed range of %s (%.3g to %.3g), the expected %s goes from about %.2f to %.2f (other conditions average). This is an association across %d years, not proof of cause.",
      cvl, gsub("_", " ", rate), eff, pr, cvl, rng[1], rng[2], gsub("_", " ", rate), vlo, vhi, T - 1))
  }
  # LTRE
  lt <- est[grepl("^ltre_", est$param), ]
  if (nrow(lt)) {
    tot <- ipm_row(est, "var_lambda_explained")
    lt <- lt[order(-lt$estimate), ]
    lab <- ipm_ltre_label(lt$param)
    out <- c(out, sprintf(
      "What drove the ups and downs of population growth (LTRE). The year-to-year variance in the growth rate attributable to the vital rates was %s. The largest contribution came from %s (%s, %.0f%% of the total), followed by %s (%s). Contributions are products of how much a rate varied over time and how sensitive growth is to it; negative values arise from negative covariation between rates. Read the ranking, not the absolute sizes: in simulations (VALIDATION.md) the posterior means of the contributions were biased upward (about +25%% for the leading rate with 15 years, more for small contributions), while the leading contributor was identified correctly in 96%% of simulated 15-year data sets and 77%% of 6-year data sets.",
      ipm_fmt(tot, 4), lab[1], ipm_fmt(lt[1, ], 4), 100 * lt$estimate[1] / tot$estimate,
      lab[2], ipm_fmt(lt[2, ], 4)))
    if (any(lt$param == "ltre_immigration")) {
      ld <- lt[!lt$param %in% c("ltre_immigration", "ltre_population_structure"), ]
      out <- c(out, sprintf(paste(
        "Caution on immigration in the LTRE. The immigration contribution (%s) and therefore the total explained variance are NOT reliable:",
        "there are no direct data on immigration, its year-to-year SD (estimated %s on the log scale) is largely set by the prior, and",
        "simulations show that IPMs overestimate the contribution of immigration, often severalfold (Paquet et al. 2021; this module's VALIDATION.md).",
        "Treat it as an upper bound. Among the rates with direct data, the largest contribution came from %s (%s)."),
        ipm_fmt(ipm_row(est, "ltre_immigration"), 4), ipm_fmt(ipm_row(est, "immigration_sd"), 2),
        ipm_ltre_label(ld$param[1]), ipm_fmt(ld[1, ], 4)))
    }
  }
  if (!is.null(o$comparison)) {
    cm <- o$comparison[!grepl("_beta_", o$comparison$param), ]
    out <- c(out, paste0(
      "Benefit of the joint analysis. Compared with analysing the capture-recapture and productivity data separately, the IPM changed the posterior SD by: ",
      paste(sprintf("%s %+.0f%%", gsub("_", " ", cm$param), 100 * (cm$sd_ratio - 1)), collapse = "; "),
      ". (Negative = more precise in the IPM.) The counts add information mostly about the rate with the weakest direct data, typically juvenile survival."))
  }
  out <- c(out, paste(
    "Caveats. The IPM assumes the three data sets are independent samples from the same population and that the population model is right.",
    "Immigration is estimated as the residual that reconciles counts with local demography, so any error elsewhere (count error that is not random, permanent emigration of adults, breeding dispersal of marked birds) is absorbed by it.",
    "Simulation (VALIDATION.md) shows that leaving immigration out when it occurs biases juvenile survival (and, less, adult survival) upwards; fecundity is held in place by the productivity data."))
  out
}

ipm_ltre_label <- function(p) {
  m <- c(ltre_juvenile_survival = "juvenile survival", ltre_adult_survival = "adult survival",
         ltre_immigration = "immigration", ltre_fecundity_yearling = "fecundity of yearlings",
         ltre_fecundity_adult = "fecundity of adults", ltre_fecundity = "fecundity",
         ltre_population_structure = "population (age) structure")
  unname(ifelse(p %in% names(m), m[p], p))
}

ipm_methods_text <- function(fit, data, settings) {
  st <- merge_list(ipm_settings, settings); fl <- fit$object$flags
  T <- length(data$years)
  sprintf(paste(
    "We fitted a female-based integrated population model (Besbeas et al. 2002; Schaub & Abadi 2011; Schaub & Kery 2022) to %d years (%d-%d) of %s.",
    "The population was described by a pre-breeding-census model with %s: the number of local recruits in year t+1 was Poisson with mean %.2f x phi_j,t x (sum of fecundity x females), surviving adults were binomial with adult apparent survival phi_a,t%s.",
    "Annual counts were modelled as %s around the total number of breeding females (state-space model).",
    "Capture-recapture data of %s animals marked as juveniles and %s marked as adults were summarised in age-dependent m-arrays and modelled with a Cormack-Jolly-Seber multinomial likelihood with first-year and adult survival and a constant recapture probability.",
    "Productivity (number of fledglings from monitored broods) was modelled as Poisson with %s.",
    "%s%s",
    "Priors were uniform: mean survival and recapture U(0,1), mean fecundity U(0,%g)%s, temporal SDs U(0,3); covariate coefficients N(0, SD 2) on standardised covariates; initial stage sizes discrete uniform.",
    "We used Nimble %s (de Valpine et al. 2017), %d chains of %d iterations (%d burn-in, thinning %d); convergence was judged by R-hat < 1.01 and bulk ESS >= 400.",
    "%s",
    "Model fit was assessed with posterior predictive checks (Freeman-Tukey discrepancies) for each data set. The model recovered known parameters in simulations at a comparable design (ecotools VALIDATION.md)."),
    T, data$years[1], data$years[T], data$species,
    if (fl$AGEF == 1) "two age classes (1-year-olds and older adults) with age-specific fecundity" else "two age classes (1-year-olds and older adults) and the same fecundity for both",
    st$sex_ratio,
    if (fl$IMM == 1) ", and immigrants were Poisson with mean omega_t x total breeding females (immigrants assumed to breed like yearlings)" else "; immigration was not modelled",
    if (fl$CNORM == 1) "normal with estimated SD" else "Poisson",
    format(data$n_marked[["juveniles"]]), if (is.na(data$n_marked[["adults"]])) "an unknown number of" else format(data$n_marked[["adults"]]),
    "mean fecundity x number of broods monitored",
    if (fl$RE == 1) "All vital rates varied among years: logit (survival) or log (fecundity, immigration) normal year effects around the mean. " else "Vital rates were constant over time. ",
    if (length(data$effects)) paste0("Year covariates entered the linear predictors: ",
                                     paste(sprintf("%s ~ %s", names(data$effects), sapply(data$effects, paste, collapse = " + ")), collapse = "; "), ". ") else "",
    st$fmax, if (fl$IMM == 1) ", mean immigration rate U(0,3)" else "",
    as.character(utils::packageVersion("nimble")), st$nchains, st$niter, st$nburnin, st$thin,
    if (fl$RE == 1 && fl$IPM == 1) "Contributions of the vital rates and of population structure to the temporal variance of the population growth rate were computed for every posterior draw with a fixed-design (variance decomposition) life-table response experiment using analytic sensitivities at the mean rates (Koons et al. 2016, 2017; Schaub & Kery 2022, section 9.4)." else "")
}

# ---- figures (docs/FIGURES.md: population_trajectory, vital_rates, ltre, ipm_vs_separate) ----

ipm_annual <- function(dr, base, yrs, ci = 0.95, inner = 0.5) {
  cn <- paste0(base, "[", yrs, "]"); ok <- cn %in% colnames(dr)
  if (!any(ok)) return(NULL)
  m <- dr[, cn[ok], drop = FALSE]
  a <- (1 - ci) / 2; b <- (1 - inner) / 2
  data.frame(time = yrs[ok], estimate = colMeans(m),
             lower = apply(m, 2, stats::quantile, a), upper = apply(m, 2, stats::quantile, 1 - a),
             lower50 = apply(m, 2, stats::quantile, b), upper50 = apply(m, 2, stats::quantile, 1 - b),
             row.names = NULL)
}

ipm_figures <- function(fit, data, opts = fig_opts()) {
  eco_require("ggplot2")
  o <- fit$object; fl <- o$flags; dr <- o$draws; est <- fit$estimates
  yrs <- data$years; T <- length(yrs)
  ci <- if (is.null(opts$ci)) 0.95 else opts$ci; inner <- if (is.null(opts$inner_ci)) 0.5 else opts$inner_ci
  cil <- round(100 * ci); inl <- round(100 * inner)
  unit <- data$count_unit
  figs <- list()
  xl <- fig_label("year", "Year", opts)
  # population trajectory
  if (fl$IPM == 1) {
    nt <- ipm_annual(dr, "N_total", yrs, ci, inner)
    obs <- data.frame(time = yrs, y = data$counts); obs <- obs[!is.na(obs$y), ]
    p <- fig_timeseries(nt, fig_label("N_total", paste0("Number of ", unit), opts),
                        obs = if (isTRUE(opts$show_data)) obs else NULL, xlab = xl) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = eco_pal$band, alpha = 0.25, colour = NA)
    figs$population_trajectory <- with_caption(p, sprintf(
      "Estimated number of %s per year (line: posterior mean; bands: %d%% and %d%% credible intervals) with the observed counts (points).",
      unit, cil, inl))
    # composition
    comp <- do.call(rbind, lapply(c("local_recruits", "surviving_adults", if (fl$IMM == 1) "immigrants"), function(b) {
      a <- ipm_annual(dr, b, yrs, ci, inner)
      data.frame(time = a$time, value = a$estimate,
                 component = fig_label(b, c(local_recruits = "Local recruits (1-year-olds)", surviving_adults = "Surviving adults",
                                            immigrants = "Immigrants")[[b]], opts))
    }))
    comp$component <- factor(comp$component, levels = unique(comp$component))
    p <- ggplot2::ggplot(comp, ggplot2::aes(time, value, fill = component)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::labs(x = xl, y = fig_label("N_total", paste0("Number of ", unit), opts), fill = NULL) +
      ggplot2::theme(legend.position = "bottom")
    figs$population_composition <- with_caption(p, sprintf(
      "Composition of the breeding population each year (posterior means): females recruited locally (born in the study area the year before), surviving adults%s.",
      if (fl$IMM == 1) " and immigrants" else ""))
    # growth rate
    lam <- ipm_annual(dr, "lambda", yrs[-T], ci, inner)
    lam$time <- lam$time + 0.5
    p <- ggplot2::ggplot(lam, ggplot2::aes(time, estimate)) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$ink, linewidth = 0.4) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = eco_pal$ink, linewidth = 1.1) +
      ggplot2::geom_point(colour = eco_pal$ink, size = 1.4) +
      ggplot2::labs(x = xl, y = fig_label("lambda", "Population growth rate (N[t+1] / N[t])", opts))
    figs$growth_rate <- with_caption(p, sprintf(
      "Realised annual population growth rate between successive breeding seasons (points: posterior means plotted between the two years; thick and thin lines: %d%% and %d%% credible intervals); dashed line = no change.", inl, cil))
  }
  # vital rates through time
  fec1 <- if (fl$AGEF == 1) "fecundity_yearling" else "fecundity"
  vr <- list(
    list(b = "juvenile_survival", panel = fig_label("juvenile_survival", "Juvenile (first-year) survival", opts), grp = "Juvenile survival", shift = 0.5, mean = "juvenile_survival"),
    list(b = "adult_survival", panel = fig_label("adult_survival", "Adult survival", opts), grp = "Adult survival", shift = 0.5, mean = "adult_survival"),
    list(b = fec1, panel = fig_label("fecundity", "Fecundity (fledglings per female)", opts),
         grp = if (fl$AGEF == 1) fig_label("fecundity_yearling", "Yearlings", opts) else "Fecundity", shift = 0, mean = fec1))
  if (fl$AGEF == 1) vr[[3]]$panel <- fig_label("fecundity_yearling", "Fecundity of yearlings (fledglings per female)", opts)
  if (fl$AGEF == 1) vr[[4]] <- list(b = "fecundity_adult", panel = fig_label("fecundity_adult", "Fecundity of older adults (fledglings per female)", opts),
                                    grp = fig_label("fecundity_adult", "Older adults", opts), shift = 0, mean = "fecundity_adult")
  if (fl$IMM == 1) vr[[length(vr) + 1]] <- list(b = "immigration_rate", panel = fig_label("immigration_rate", "Immigration rate (per female)", opts),
                                                grp = "Immigration", shift = 0.5, mean = "immigration_rate")
  vd <- do.call(rbind, lapply(vr, function(v) {
    a <- ipm_annual(dr, v$b, if (v$shift > 0) yrs[-T] else yrs, ci, inner)
    if (is.null(a)) return(NULL)
    a$time <- a$time + v$shift; a$panel <- v$panel; a$series <- v$grp
    a$mean <- ipm_row(est, v$mean)$estimate
    a
  }))
  vd$panel <- factor(vd$panel, levels = unique(vd$panel))
  vd$series <- factor(vd$series, levels = unique(vd$series))
  pd <- ggplot2::position_dodge(width = 0.35)
  p <- ggplot2::ggplot(vd, ggplot2::aes(time, estimate, colour = series, shape = series)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = mean, colour = series), linetype = 2, linewidth = 0.3, show.legend = FALSE) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), linewidth = 0.35, position = pd) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), linewidth = 0.9, position = pd) +
    ggplot2::geom_point(size = 1.2, position = pd) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 1) +
    ggplot2::guides(colour = "none", shape = "none") +
    ggplot2::labs(x = xl, y = NULL, colour = NULL, shape = NULL)
  attr(p, "fig_size") <- list(width = "double", height = 34 * length(levels(vd$panel)) + 15)
  figs$vital_rates <- with_caption(p, sprintf(
    "Annual vital rates (posterior means with %d%% [thick] and %d%% [thin] credible intervals); dashed lines = mean rate. Survival and immigration refer to the interval from one breeding season to the next and are plotted between the two years.%s",
    inl, cil, if (fl$AGEF == 1) " Fecundity: fledglings (both sexes) per breeding female, yearlings vs older adults." else ""))
  # LTRE
  lt <- est[grepl("^ltre_", est$param), ]
  if (nrow(lt)) {
    a <- (1 - ci) / 2
    lq <- t(apply(dr[, lt$param, drop = FALSE], 2, stats::quantile, c(a, 1 - a, (1 - inner) / 2, 1 - (1 - inner) / 2)))
    ld <- data.frame(param = vapply(lt$param, function(k) fig_label(k, tools::toTitleCase(ipm_ltre_label(k)), opts), ""),
                     estimate = lt$estimate, lower = lq[, 1], upper = lq[, 2], lower50 = lq[, 3], upper50 = lq[, 4],
                     kind = ifelse(lt$param == "ltre_population_structure", "Population structure", "Vital rate"))
    ld$param <- stats::reorder(factor(ld$param), ld$estimate)
    p <- ggplot2::ggplot(ld, ggplot2::aes(estimate, param)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = eco_pal$null) +
      ggplot2::geom_col(fill = eco_pal$accent %||% "#4E9A8A", width = 0.6, alpha = 0.6) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.4, colour = eco_pal$ink) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower50, xmax = upper50), linewidth = 1.1, colour = eco_pal$ink) +
      ggplot2::scale_x_continuous(n.breaks = 4, labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
      ggplot2::labs(x = fig_label("ltre_x", "Contribution to variance in population growth rate", opts), y = NULL, fill = NULL)
    attr(p, "fig_size") <- list(width = "double", height = 30 + 9 * nrow(ld))
    figs$ltre <- with_caption(p, sprintf(
      "Life-table response experiment: contribution of each vital rate and of population structure to the temporal variance of the population growth rate (bars: posterior means; thick and thin lines: %d%% and %d%% credible intervals). Read the ranking rather than the sizes: posterior means of contributions are biased upward, and the immigration contribution is overestimated, often severalfold, because immigration has no direct data (Paquet et al. 2021; VALIDATION.md).", inl, cil))
  }
  if (!is.null(o$transient)) {
    tc <- o$transient
    keys <- intersect(c("phij", "phia", "f1", "f2", "om", "n1", "n2", "n3"), colnames(tc))
    labs <- c(phij = "Juvenile survival", phia = "Adult survival", om = "Immigration", f1 = if (fl$AGEF == 1) "Fecundity yearlings" else "Fecundity",
              f2 = "Fecundity adults", n1 = "Structure: local recruits", n2 = "Structure: adults", n3 = "Structure: immigrants")
    td <- do.call(rbind, lapply(keys, function(k) data.frame(time = yrs[2:(T - 1)], value = tc[, k],
                                                             component = fig_label(paste0("ltre_", k), labs[[k]], opts))))
    td$component <- factor(td$component, levels = unique(td$component))
    dl <- data.frame(time = yrs[2:(T - 1)], value = o$lambda_diff)
    p <- ggplot2::ggplot(td, ggplot2::aes(time, value, fill = component)) +
      ggplot2::geom_hline(yintercept = 0, colour = eco_pal$muted, linewidth = 0.3) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_point(data = dl, ggplot2::aes(time, value), inherit.aes = FALSE, shape = 21, fill = "white", colour = "black", size = 1.6) +
      ggplot2::labs(x = fig_label("year_change", "Year (change from previous growth rate)", opts),
                    y = fig_label("transient_y", "Contribution to change in growth rate", opts), fill = NULL) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
      ggplot2::theme(legend.position = "bottom")
    attr(p, "fig_size") <- list(width = "double", aspect = 0.55)
    figs$transient_ltre <- with_caption(p, "Transient LTRE: contributions of each vital rate and of population structure to the change in population growth rate from one year to the next (bars, posterior means); white points = estimated change in growth rate.")
  }
  # IPM vs separate analyses
  if (!is.null(o$comparison)) {
    cm <- o$comparison
    cd <- rbind(data.frame(param = cm$param, analysis = "Integrated (IPM)", estimate = cm$ipm_estimate, lower = cm$ipm_lower, upper = cm$ipm_upper, cv = cm$ipm_cv),
                data.frame(param = cm$param, analysis = "Separate analyses", estimate = cm$separate_estimate, lower = cm$separate_lower, upper = cm$separate_upper, cv = cm$separate_cv))
    cd$param <- vapply(cd$param, function(k) fig_label(k, ipm_param_label(k, opts), opts), "")
    cd$label <- sprintf("CV %.0f%%", 100 * cd$cv)
    cd$analysis <- factor(cd$analysis, levels = c("Separate analyses", "Integrated (IPM)"))
    p <- ggplot2::ggplot(cd, ggplot2::aes(estimate, analysis, colour = analysis, shape = analysis)) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.6) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::geom_text(ggplot2::aes(x = upper, label = label), hjust = -0.15, size = 2.3, colour = eco_pal$muted, show.legend = FALSE) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.45))) +
      ggplot2::facet_wrap(~param, scales = "free_x", ncol = 2, labeller = ggplot2::label_wrap_gen(width = 38)) +
      ggplot2::labs(x = sprintf("Posterior mean and %d%% credible interval", cil), y = NULL, colour = NULL, shape = NULL) +
      ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    attr(p, "fig_size") <- list(width = "double", height = 25 + 22 * ceiling(length(unique(cd$param)) / 2))
    figs$ipm_vs_separate <- with_caption(p, sprintf(
      "Precision of the integrated population model versus separate analyses of the capture-recapture and productivity data (same priors and year effects). Lines: %d%% credible intervals; labels: coefficient of variation (posterior SD / mean).", cil))
  }
  # covariate effects (response scale, original units)
  for (rate in names(data$effects)) for (cv in data$effects[[rate]]) {
    bn <- paste0(rate, "_beta_", cv)
    if (!bn %in% colnames(dr)) next
    sc <- o$scaling[[cv]]
    meanp <- switch(rate, juvenile_survival = "juvenile_survival", adult_survival = "adult_survival",
                    fecundity = if (fl$AGEF == 1) "fecundity_adult" else "fecundity", immigration = "immigration_rate_median")
    link <- if (rate %in% c("juvenile_survival", "adult_survival")) "logit" else "log"
    inv <- if (link == "logit") stats::plogis else exp; lk <- if (link == "logit") stats::qlogis else log
    used <- if (rate == "fecundity") seq_len(T) else seq_len(T - 1)  # covariate years that enter the model
    xr <- range(data$covariates[[cv]][used], na.rm = TRUE)
    xs <- seq(xr[1], xr[2], length.out = 80); z <- (xs - sc$center) / sc$scale
    idx <- ipm_thin_idx(nrow(dr), 2000)
    pm <- sapply(z, function(zz) inv(lk(dr[idx, meanp]) + dr[idx, bn] * zz))
    a <- (1 - ci) / 2; b <- (1 - inner) / 2
    pred <- data.frame(x = fig_x(xs, cv, opts, NULL), estimate = colMeans(pm),
                       lower = apply(pm, 2, stats::quantile, a), upper = apply(pm, 2, stats::quantile, 1 - a),
                       lower50 = apply(pm, 2, stats::quantile, b), upper50 = apply(pm, 2, stats::quantile, 1 - b))
    ann_b <- switch(rate, juvenile_survival = "juvenile_survival", adult_survival = "adult_survival",
                    fecundity = if (fl$AGEF == 1) "fecundity_adult" else "fecundity", immigration = "immigration_rate")
    ay <- if (rate == "fecundity") yrs else yrs[-T]
    an <- ipm_annual(dr, ann_b, ay, ci, inner)
    obs <- if (isTRUE(opts$show_data) && !is.null(an))
      data.frame(x = fig_x(data$covariates[[cv]][match(an$time, yrs)], cv, opts, NULL), y = an$estimate) else NULL
    ylab <- fig_label(rate, c(juvenile_survival = "Juvenile (first-year) survival", adult_survival = "Adult survival",
                              fecundity = "Fecundity (fledglings per female)", immigration = "Immigration rate (per female)")[[rate]], opts)
    p <- fig_cond_effect(pred, fig_label(cv, cv, opts), ylab, obs = obs, ylim01 = FALSE)
    if (link == "logit") p <- p + ggplot2::coord_cartesian(ylim = c(0, max(pred$upper, if (!is.null(obs)) obs$y) * 1.05))
    figs[[paste0("effect_", cv)]] <- with_caption(p, sprintf(
      "Expected %s in relation to %s (line: posterior mean; bands: %d%% and %d%% credible intervals; year random effects at zero%s). Points: annual estimates from the IPM.",
      gsub("_", " ", rate), fig_label(cv, cv, opts), cil, inl, if (rate == "fecundity" && fl$AGEF == 1) ", older adults" else ""))
  }
  # goodness of fit
  if (!is.null(o$ppc)) {
    pp <- o$ppc
    nm <- c(counts = "Population counts", productivity = "Productivity", capture_recapture = "Capture-recapture")
    gd <- do.call(rbind, lapply(colnames(pp$obs), function(k) data.frame(obs = pp$obs[, k], rep = pp$rep[, k],
      panel = sprintf("%s (P = %.2f)", fig_label(paste0("gof_", k), nm[[k]], opts), pp$p[[k]]))))
    p <- ggplot2::ggplot(gd, ggplot2::aes(obs, rep)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = eco_pal$muted) +
      ggplot2::geom_point(alpha = 0.3, size = 0.6, colour = eco_pal$ink) +
      ggplot2::facet_wrap(~panel, scales = "free", nrow = 1) +
      ggplot2::labs(x = "Discrepancy, observed data", y = "Discrepancy, replicated data")
    attr(p, "fig_size") <- list(width = "double", aspect = 0.38)
    figs$gof <- with_caption(p, "Posterior predictive checks (Schaub & Kery 2022, ch. 11): Freeman-Tukey discrepancy of the observed vs replicated data for each data set, one point per posterior draw. P = Bayesian p-value (share of points above the 1:1 line); values near 0 or 1 indicate lack of fit.")
  }
  # MCMC diagnostics
  main <- intersect(c("juvenile_survival", "adult_survival", fec1, "fecundity_adult", "immigration_rate", "recapture_prob",
                      "N_mean", "lambda_geomean"), colnames(dr))
  if (length(o$chains) > 1) {
    pt <- fig_trace(lapply(o$chains, function(m) m[, main, drop = FALSE]))
    attr(pt, "fig_size") <- list(width = "double", height = 30 * ceiling(length(main) / 3) + 20)
    figs$trace <- with_caption(pt, "MCMC traces (thinned draws, one colour per chain) of the main estimands; chains should overlap without trends.")
  }
  pp <- fig_posterior(dr[, main, drop = FALSE], main, ncol = 3)
  attr(pp, "fig_size") <- list(width = "double", height = 30 * ceiling(length(main) / 3) + 15)
  figs$posterior <- with_caption(pp, "Posterior densities of the main estimands.")
  figs
}


# readable names for summary parameters (figure labels)
ipm_param_label <- function(k, opts = NULL) {
  m <- c(juvenile_survival = "Juvenile (first-year) survival", adult_survival = "Adult survival",
         fecundity_yearling = "Fecundity of yearlings", fecundity_adult = "Fecundity of older adults",
         fecundity = "Fecundity", recapture_prob = "Recapture probability",
         immigration_rate = "Immigration rate")
  if (k %in% names(m)) return(m[[k]])
  b <- regmatches(k, regexec("^(.*)_beta_(.*)$", k))[[1]]
  if (length(b) == 3) {
    rl <- if (b[2] %in% names(m)) m[[b[2]]] else tools::toTitleCase(gsub("_", " ", b[2]))
    cv <- if (!is.null(opts)) fig_label(b[3], gsub("_", " ", b[3]), opts) else gsub("_", " ", b[3])
    return(sprintf("%s: effect of %s (per SD, %s scale)", rl, cv, if (grepl("survival", b[2])) "logit" else "log"))
  }
  tools::toTitleCase(gsub("_", " ", k))
}
