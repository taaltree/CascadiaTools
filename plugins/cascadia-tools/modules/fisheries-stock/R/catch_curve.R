# Catch curves: total instantaneous mortality Z from catch-at-age.
#
# Estimators (Smith et al. 2012 recommendations):
#   * Chapman-Robson (1960) survival S = T / (n + T - 1) on recoded ages,
#     Z = -log(S) with the Hoenig et al. (1983) bias correction, SE of Z
#     (Smith et al. 2012, eq. 2 as in FSA) inflated by an overdispersion
#     factor c-hat = Pearson chi-square / df (floored at 1), first age =
#     peak age + 1 ("peak plus").
#   * Weighted linear regression of log(catch) on age (Maceina & Bettoli
#     1998; FSA::catchCurve(weighted = TRUE)), first age = peak age.
#
# data$df: data.frame(age, catch) with one row per age (zeros filled).

cc_selectivity <- function(age, a50, slope) stats::plogis(slope * (age - a50))

cc_simulate <- function(params, design, seed) {
  set.seed(seed)
  amax <- design$max_age %||% 20
  ages <- 0:amax
  Z <- params$Z
  sR <- params$recruit_sd_cc %||% 0
  phiR <- params$recruit_ar1_cc %||% 0
  # Year-class strength R_{year-a}: lognormal with optional AR(1) in time.
  e <- numeric(length(ages)); e[1] <- stats::rnorm(1, 0, sR)
  for (i in seq_along(ages)[-1]) e[i] <- phiR * e[i - 1] + stats::rnorm(1, 0, sR * sqrt(1 - phiR^2))
  trend <- params$recruit_trend_cc %||% 0          # log change in recruitment per year class (neg. control)
  R <- exp(e - sR^2 / 2 - trend * ages)            # age a was born a years ago
  Za <- Z + (params$Z_age_slope %||% 0) * pmax(ages - (params$sel_age50 %||% 2.5), 0)
  surv <- exp(-c(0, cumsum(Za[-length(Za)])))
  N <- R * surv
  sel <- cc_selectivity(ages, params$sel_age50 %||% 2.5, params$sel_slope %||% 4)
  if (!is.null(params$sel_dome_age)) sel <- sel * ifelse(ages > params$sel_dome_age, exp(-params$sel_dome_rate * (ages - params$sel_dome_age)), 1)
  p <- N * sel; p <- p / sum(p)
  catch <- as.vector(stats::rmultinom(1, design$n_fish, p))
  df <- data.frame(age = ages, catch = catch)
  # trim trailing zero ages beyond the oldest fish (as a real table would be)
  last <- max(df$age[df$catch > 0]); first <- min(df$age[df$catch > 0])
  df <- df[df$age >= first & df$age <= last, ]
  list(data = list(model = "catch_curve", df = df, units = list(),
                   meta = list(source = "simulated")),
       truth = c(total_mortality_Z = Z, annual_survival = exp(-Z), annual_mortality = 1 - exp(-Z),
                 Z_chapman_robson = Z, Z_weighted_regression = Z, Z_glmm = Z))
}

# Choose first age for an estimator.
cc_first_age <- function(df, rule) {
  peak <- df$age[which.max(df$catch)]
  if (is.numeric(rule)) return(as.numeric(rule))
  switch(rule, peak_plus = peak + 1, peak = peak,
         stop("age_full must be 'peak_plus', 'peak' or an age (number)"))
}

# Chapman-Robson with Hoenig/Smith correction and Pearson c-hat.
cc_chapman_robson <- function(age, catch, overdispersion = TRUE) {
  x <- age - min(age)
  n <- sum(catch); Tt <- sum(x * catch)
  S <- Tt / (n + Tt - 1)
  Z <- -log(S) - ((n - 1) * (n - 2)) / (n * (Tt + 1) * (n + Tt - 1))
  se0 <- (1 - exp(-Z)) / sqrt(n * exp(-Z))
  # Goodness of fit of the geometric age distribution: cells x = 0..k-1 with
  # the last cell holding the tail (x >= k-1); old cells pooled until the
  # expected count is >= 1. df = cells - 2 (1 constraint, 1 parameter).
  k <- length(x)
  pr <- (1 - S) * S^x; pr[k] <- S^x[k]
  obs <- catch; ex <- n * pr
  while (length(ex) > 2 && ex[length(ex)] < 1) {
    j <- length(ex); ex[j - 1] <- ex[j - 1] + ex[j]; obs[j - 1] <- obs[j - 1] + obs[j]
    ex <- ex[-j]; obs <- obs[-j]
  }
  chi <- sum((obs - ex)^2 / ex)
  df <- length(ex) - 2
  chat <- if (df >= 1) max(1, chi / df) else 1
  se <- if (overdispersion) se0 * sqrt(chat) else se0
  list(S = S, Z = Z, se = se, se_raw = se0, chat = chat, chisq = chi, df = df, n = n, T = Tt,
       p_gof = if (df >= 1) stats::pchisq(chi, df, lower.tail = FALSE) else NA_real_)
}

# Millar (2015) random-intercept Poisson GLMM: catch_a ~ Poisson(exp(b0 - Z a + e_a)),
# e_a ~ N(0, s^2) (year-class strength variability), fitted with glmmTMB.
cc_glmm <- function(age, catch) {
  eco_require("glmmTMB")
  d <- data.frame(age = age, catch = catch, yc = factor(age))
  m <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(catch ~ age + (1 | yc), family = stats::poisson, data = d)),
                error = function(e) e)
  if (inherits(m, "error")) return(list(Z = NA_real_, se = NA_real_, ok = FALSE, msg = conditionMessage(m), sd_re = NA_real_))
  co <- summary(m)$coefficients$cond
  sdre <- sqrt(as.numeric(glmmTMB::VarCorr(m)$cond$yc))
  ok <- isTRUE(m$sdr$pdHess) && is.finite(co["age", 2])
  list(Z = -co["age", 1], se = co["age", 2], ok = ok, msg = if (ok) "" else "Hessian not positive definite", sd_re = sdre, model = m)
}

cc_fit <- function(data, settings) {
  eco_require("FSA")
  df <- data$df
  df <- df[order(df$age), ]
  first_cr <- cc_first_age(df, settings$age_full %||% "peak_plus")
  first_rg <- cc_first_age(df, settings$age_full_regression %||% "peak")
  last <- settings$age_last %||% max(df$age)
  use_cr <- df$age >= first_cr & df$age <= last
  use_rg <- df$age >= first_rg & df$age <= last & df$catch > 0
  notes <- character()
  if (sum(use_cr) < 3) stop(sprintf("Only %d ages from age %g on: too few for a catch curve.", sum(use_cr), first_cr))
  cr <- cc_chapman_robson(df$age[use_cr], df$catch[use_cr], isTRUE(settings$overdispersion %||% TRUE))
  z <- stats::qnorm(0.975)
  # quasi-likelihood: with c-hat estimated, use t on the chi-square df
  zc <- if (isTRUE(settings$overdispersion %||% TRUE) && cr$chat > 1 && cr$df >= 1) stats::qt(0.975, cr$df) else z
  # Weighted regression (FSA engine)
  wn <- character()
  cc <- withCallingHandlers(FSA::catchCurve(df$age, df$catch, ages2use = df$age[use_rg], weighted = TRUE),
                            warning = function(w) { wn <<- c(wn, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (length(wn)) notes <- c(notes, "Weighted regression: the oldest age(s) had a predicted log catch <= 0 and were given weight 0 (Maceina & Bettoli 1998 weighting).")
  co <- summary(cc)
  Zr <- -unname(stats::coef(cc$lm)[2]); seR <- co["Z", "Std. Error"]
  ci_r <- stats::confint(cc$lm, level = 0.95)[2, ]
  gm <- cc_glmm(df$age[use_cr], df$catch[use_cr])
  zt <- rbind(
    chapman_robson = c(cr$Z, cr$se, cr$Z - zc * cr$se, cr$Z + zc * cr$se),
    weighted_regression = c(Zr, seR, -ci_r[2], -ci_r[1]),
    glmm = c(gm$Z, gm$se, gm$Z - z * gm$se, gm$Z + z * gm$se))
  prim <- match.arg(settings$primary %||% "weighted_regression", rownames(zt))
  P <- zt[prim, ]
  if (!all(is.finite(P))) { notes <- c(notes, sprintf("Primary estimator '%s' failed; Chapman-Robson used instead.", prim)); prim <- "chapman_robson"; P <- zt[prim, ] }
  est <- rbind(
    fs_rows("total_mortality_Z", P[1], P[2], P[3], P[4]),
    fs_rows("annual_survival", exp(-P[1]), exp(-P[1]) * P[2], exp(-P[4]), exp(-P[3])),
    fs_rows("annual_mortality", 1 - exp(-P[1]), exp(-P[1]) * P[2], 1 - exp(-P[3]), 1 - exp(-P[4])),
    fs_rows("Z_chapman_robson", zt[1, 1], zt[1, 2], zt[1, 3], zt[1, 4]),
    fs_rows("Z_weighted_regression", zt[2, 1], zt[2, 2], zt[2, 3], zt[2, 4]),
    fs_rows("Z_glmm", zt[3, 1], zt[3, 2], zt[3, 3], zt[3, 4]))
  if (!gm$ok) notes <- c(notes, paste("Poisson-lognormal GLMM catch curve did not converge cleanly:", gm$msg))
  if (cr$n < 50) notes <- c(notes, sprintf("Only %d fish at or above the age of full recruitment: intervals will be wide and the age choice matters a lot.", cr$n))
  notes <- c(notes, sprintf("Chapman-Robson used ages %g-%g (n = %d fish; peak age %g + 1); regression used ages %g-%g with nonzero catch (%d ages).",
                            first_cr, last, cr$n, df$age[which.max(df$catch)], first_rg, last, sum(use_rg)))
  notes <- c(notes, sprintf("Overdispersion c-hat = %.2f (Pearson chi-square %.1f on %d df, P = %.3f)%s.", cr$chat, cr$chisq, cr$df, cr$p_gof,
                            if (cr$chat > 1 && isTRUE(settings$overdispersion %||% TRUE)) "; the Chapman-Robson SE was inflated by sqrt(c-hat)" else ""))
  notes <- c(notes, sprintf("Primary estimator (total_mortality_Z, annual survival and mortality): %s.", gsub("_", " ", prim)))
  if (abs(Zr - cr$Z) > 1.96 * sqrt(cr$se^2 + seR^2))
    notes <- c(notes, sprintf("Chapman-Robson and regression disagree (Z = %.3f vs %.3f): look for curvature in the catch curve (dome-shaped selectivity, age-varying mortality or a recruitment trend).", cr$Z, Zr))
  list(estimates = est,
       diagnostics = list(converged = is.finite(cr$Z) && is.finite(cr$se), rhat_max = NA_real_, ess_min = NA_real_,
                          notes = notes, chat = cr$chat, p_gof = cr$p_gof),
       object = list(model = "catch_curve", df = df, cr = cr, reg = cc, glmm = gm, zt = zt, primary = prim, first_cr = first_cr, first_rg = first_rg,
                     last = last, use_cr = use_cr, use_rg = use_rg),
       settings = settings)
}

cc_read_data <- function(config, dc) {
  f <- file.path(config$.dir, dc$file)
  raw <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  ages <- fs_num(fs_col(raw, dc$age_col %||% "age", "age"), dc$age_col %||% "age", "age")
  if (all(is.na(ages))) stop("The age column has no ages.", call. = FALSE)
  if (sum(is.na(ages)) > 0) message(sprintf("%d of %d rows have no age and were left out.", sum(is.na(ages)), length(ages)))
  if (!is.null(dc$count_col)) {
    cnt <- fs_num(fs_col(raw, dc$count_col, "count"), dc$count_col, "count")
    tab <- stats::aggregate(cnt, list(age = ages), sum, na.rm = TRUE); names(tab)[2] <- "catch"
  } else {
    ages <- ages[!is.na(ages)]
    tab <- as.data.frame(table(age = ages), stringsAsFactors = FALSE); names(tab)[2] <- "catch"
    tab$age <- as.numeric(tab$age)
  }
  full <- data.frame(age = seq(min(tab$age), max(tab$age)))
  full$catch <- 0
  full$catch[match(tab$age, full$age)] <- tab$catch
  list(model = "catch_curve", df = full, raw_n = nrow(raw),
       n_missing_age = sum(is.na(ages)),
       meta = list(file = dc$file, species = dc$species %||% "fish", gear = dc$gear %||% ""))
}

cc_check_data <- function(data) {
  pr <- character(); df <- data$df
  if (any(df$age != round(df$age))) pr <- c(pr, "Ages must be whole years (integers). Round or check the age column.")
  if (any(df$catch < 0, na.rm = TRUE)) pr <- c(pr, "Negative catches found in the catch-at-age table.")
  if (anyNA(df$catch)) pr <- c(pr, "Missing catches in the catch-at-age table: use 0 for ages with no fish.")
  if (sum(df$catch) < 30) pr <- c(pr, sprintf("Only %d aged fish: a catch curve needs at least ~50-100 aged fish on the descending limb.", sum(df$catch)))
  peak <- df$age[which.max(df$catch)]
  if (sum(df$age > peak + 1) < 2) pr <- c(pr, sprintf("Fewer than 3 ages older than the peak age + 1 (peak = age %g): no descending limb to fit.", peak))
  pr
}

CC_LABELS <- c(chapman_robson = "Chapman-Robson", weighted_regression = "weighted regression", glmm = "Poisson-lognormal GLMM")

cc_interpret <- function(fit, data) {
  e <- fit$estimates; g <- function(p) e[e$param == p, ]
  Z <- g("total_mortality_Z"); S <- g("annual_survival"); A <- g("annual_mortality")
  o <- fit$object; zt <- o$zt
  sp <- data$meta$species %||% "fish"
  first <- if (o$primary == "weighted_regression") o$first_rg else o$first_cr
  others <- setdiff(rownames(zt), o$primary)
  agree <- vapply(others, function(k) is.finite(zt[k, 1]) && abs(zt[k, 1] - zt[o$primary, 1]) <= 1.96 * sqrt(zt[k, 2]^2 + zt[o$primary, 2]^2), logical(1))
  c(sprintf("Total instantaneous mortality of fully recruited %s (ages %g and older) is Z = %.2f per year (95%% CI %.2f to %.2f; %s). Equivalently, about %.0f%% (%.0f%% to %.0f%%) of fully recruited fish survive each year and %.0f%% (%.0f%% to %.0f%%) die from all causes combined (fishing plus natural mortality).",
            sp, first, Z$estimate, Z$lower, Z$upper, CC_LABELS[[o$primary]], 100 * S$estimate, 100 * S$lower, 100 * S$upper, 100 * A$estimate, 100 * A$lower, 100 * A$upper),
    sprintf("The other estimators give %s. %s",
            paste(sprintf("Z = %.2f (%.2f to %.2f; %s)", zt[others, 1], zt[others, 3], zt[others, 4], CC_LABELS[others]), collapse = " and "),
            if (all(agree)) "All estimators agree within their uncertainty." else "They do not all agree: the catch curve is probably not a straight line, so check the constant-mortality, constant-recruitment and full-selectivity assumptions."),
    sprintf("The age distribution is %s the smooth decline expected from constant mortality and recruitment (Chapman-Robson c-hat = %.2f). %s",
            if (o$cr$chat > 1.5) "much more variable than" else if (o$cr$chat > 1.1) "somewhat more variable than" else "consistent with",
            o$cr$chat,
            if (o$cr$chat > 1.1) "This usually reflects variable year-class strength. In this module's simulations with variable year classes, all three 95% intervals were too narrow: the weighted-regression interval covered the true Z in about 84-93% of simulated samples, Chapman-Robson about 77-92% and the GLMM 79-91% (see VALIDATION.md), so treat the interval as closer to a 90% interval. A strong or weak year class at either end of the age range also shifts Z." else ""),
    "A catch curve assumes constant recruitment (or no trend in it), constant mortality and full, equal vulnerability of every age from the first age used onwards. Z mixes fishing and natural mortality; separating them needs an independent estimate of natural mortality M (F = Z - M).")
}

cc_methods <- function(fit, data, settings) {
  o <- fit$object
  sprintf(paste(
    "Total instantaneous mortality (Z) was estimated from catch-at-age data (%d aged fish, ages %g-%g) with three catch-curve estimators following Smith et al. (2012) and Millar (2015).",
    "(1) Weighted linear regression of log catch on age from the modal age (age %g) onwards, excluding ages with zero catch, with weights equal to the predicted log catch of an unweighted fit (Maceina & Bettoli 1998), using FSA::catchCurve (Ogle et al. 2025); 95%% t intervals.",
    "(2) Chapman-Robson (Chapman & Robson 1960) from one year older than the modal age (age %g; the 'peak plus' criterion of Smith et al. 2012), with the Hoenig et al. (1983) bias correction for Z and the standard error of Z (Smith et al. 2012) inflated by sqrt(c-hat), where c-hat (= %.2f) is the Pearson chi-square goodness-of-fit statistic of the geometric age distribution divided by its degrees of freedom (floored at 1), with t intervals on those degrees of freedom.",
    "(3) A Poisson generalized linear mixed model of catch on age with an age-class (year-class) random intercept (Millar 2015), fitted with glmmTMB on the same ages as Chapman-Robson; Wald intervals.",
    "The %s estimate is reported as the primary estimate; annual survival and mortality are exp(-Z) and 1 - exp(-Z).",
    "Bias and interval coverage of each estimator were checked by simulation from an age-structured population with lognormal year-class variability (ecotools fisheries-stock module, VALIDATION.md)."),
    sum(o$df$catch), min(o$df$age), max(o$df$age), o$first_rg, o$first_cr, o$cr$chat, CC_LABELS[[o$primary]])
}

cc_figures <- function(fit, data, opts) {
  o <- fit$object; df <- o$df
  ci <- opts$ci %||% 0.95; z <- stats::qnorm(1 - (1 - ci) / 2)
  # Estimators are entities: fixed order -> same colour in every figure.
  est_lab <- c(chapman_robson = fig_label("chapman_robson", "Chapman-Robson", opts),
               weighted_regression = fig_label("weighted_regression", "Weighted regression", opts),
               glmm = fig_label("glmm", "Poisson-lognormal GLMM", opts))
  est_col <- stats::setNames(fs_pal(3, opts), est_lab)
  est_lty <- stats::setNames(c(1, 2, 4), est_lab)
  pt_lab <- c(cr = fig_label("used_cr", "Used by both", opts),
              rg = fig_label("used_rg", "Regression only", opts),
              no = fig_label("not_used", "Not fully recruited", opts))
  d <- df[df$catch > 0, ]
  d$used <- factor(ifelse(o$use_cr[df$catch > 0], pt_lab[["cr"]],
                          ifelse(o$use_rg[df$catch > 0], pt_lab[["rg"]], pt_lab[["no"]])), levels = pt_lab)
  # Chapman-Robson expected log catch: n (1-S) S^x on the used ages
  cr <- o$cr; a0 <- o$first_cr
  ages_cr <- seq(a0, o$last, length.out = 100)
  S <- exp(-cr$Z)
  # drawn on the response scale (number of fish) with a log axis
  lcr <- data.frame(x = ages_cr, y = cr$n * (1 - S) * exp(-cr$Z * (ages_cr - a0)),
                    est = factor(est_lab[["chapman_robson"]], levels = est_lab))
  # weighted-regression line and confidence band
  ages_rg <- seq(o$first_rg, o$last, length.out = 100)
  pr <- stats::predict(o$reg$lm, newdata = data.frame(age.e = ages_rg), se.fit = TRUE)
  # t quantile on the residual df, matching the reported confint() interval for Z
  tq <- stats::qt(1 - (1 - ci) / 2, stats::df.residual(o$reg$lm))
  rg <- data.frame(x = ages_rg, y = exp(pr$fit), lower = exp(pr$fit - tq * pr$se.fit), upper = exp(pr$fit + tq * pr$se.fit),
                   est = factor(est_lab[["weighted_regression"]], levels = est_lab))
  # Default: estimator key above the panel, age classes below, one row each
  # (two side-by-side keys do not fit a single column). A user-chosen legend
  # position puts both keys there.
  lpos <- fs_legend_pos(opts); stacked <- identical(opts$legend %||% "right", "right")
  lpos_est <- if (stacked) "top" else lpos
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = rg, ggplot2::aes(x, ymin = lower, ymax = upper), fill = est_col[[2]], alpha = 0.18) +
    ggplot2::geom_line(data = rbind(lcr, rg[, c("x", "y", "est")]), ggplot2::aes(x, y, colour = est, linetype = est), linewidth = 0.7) +
    ggplot2::geom_point(data = d, ggplot2::aes(age, catch, shape = used), size = 1.8, colour = eco_pal$ink, stroke = 0.6) +
    ggplot2::scale_colour_manual(values = est_col, drop = TRUE) +
    ggplot2::scale_linetype_manual(values = est_lty, drop = TRUE) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 0, 4), pt_lab), drop = TRUE) +
    ggplot2::scale_x_continuous(breaks = function(l) { b <- seq(ceiling(l[1]), floor(l[2]), by = 1); if (length(b) > 10) b[b %% 2 == 0] else b }) +
    ggplot2::scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000),
                           labels = function(b) format(b, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)) +
    ggplot2::labs(x = fig_label("age", "Age (years)", opts), y = fig_label("log_catch", "Number caught (log scale)", opts),
                  colour = fig_label("estimator_title", NULL, opts), linetype = fig_label("estimator_title", NULL, opts),
                  shape = fig_label("ages_title", NULL, opts)) +
    ggplot2::guides(colour = ggplot2::guide_legend(position = lpos_est, order = 1, nrow = if (stacked) 1 else NULL, keywidth = ggplot2::unit(1.8, "lines")),
                    linetype = ggplot2::guide_legend(position = lpos_est, order = 1, nrow = if (stacked) 1 else NULL, keywidth = ggplot2::unit(1.8, "lines")),
                    shape = ggplot2::guide_legend(position = lpos, order = 2, nrow = if (stacked) 2 else NULL, byrow = TRUE))
  zt <- o$zt
  p <- p + ggplot2::labs(title = sprintf("Z = %.2f (Chapman-Robson), %.2f (regression), %.2f (GLMM)", zt[1, 1], zt[2, 1], zt[3, 1]))
  p <- fs_size(p, "catch_curve", opts, aspect = 0.95)
  # Intervals at the figure's ci level, built exactly as cc_fit builds the 95%
  # ones: Chapman-Robson t on the c-hat df (when c-hat > 1), regression t on
  # its residual df, GLMM Wald. At ci = 0.95 these equal the estimates table.
  hw <- c(chapman_robson = if (isTRUE(fit$settings$overdispersion %||% TRUE) && o$cr$chat > 1 && o$cr$df >= 1)
                             stats::qt(1 - (1 - ci) / 2, o$cr$df) else z,
          weighted_regression = tq, glmm = z) * zt[, 2]
  zz <- data.frame(param = est_lab, estimate = zt[, 1], lower = zt[, 1] - hw, upper = zt[, 1] + hw,
                   group = factor(est_lab, levels = est_lab))
  zz <- zz[is.finite(zz$estimate), ]
  zz$param <- factor(zz$param, levels = rev(est_lab))
  pz <- fig_coef(zz, null_value = zt[o$primary, 1], order = FALSE,
                 xlab = fig_label("Z", "Total mortality, Z (per year)", opts)) +
    ggplot2::scale_colour_manual(values = est_col, guide = "none")
  # fig_coef dodges grouped rows; one row per estimator, so no dodge needed
  pz <- fs_size(pz, "mortality_estimates", opts, aspect = 0.5)
  list(
    catch_curve = with_caption(p, sprintf(paste(
      "Catch curve: number of fish caught at each age on a log scale (points; filled circles = ages used by Chapman-Robson and the regression,",
      "open squares = used by the regression only, crosses = not fully recruited and excluded). Solid line: Chapman-Robson expected catch (Z = %.2f);",
      "dashed line and band: weighted regression with its %d%% confidence band for expected catch (Z = %.2f). The GLMM uses the same ages as Chapman-Robson. Ages with zero catch are not shown."),
      zt[1, 1], round(100 * ci), zt[2, 1])),
    mortality_estimates = with_caption(pz, sprintf(paste(
      "Total instantaneous mortality Z by estimator (points) with %d%% confidence intervals (Chapman-Robson: t interval including the overdispersion correction;",
      "weighted regression: t interval; GLMM: Wald interval). Dashed line: the primary (%s) estimate. Annual survival is exp(-Z)."), round(100 * ci), CC_LABELS[[o$primary]]))
  )
}
