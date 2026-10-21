# Catch curves: total instantaneous mortality Z from catch-at-age.
#
# Single-year estimators (one age sample):
#   * Chapman-Robson (1960) survival S = T / (n + T - 1) on recoded ages,
#     Z = -log(S) with the Hoenig et al. (1983) bias correction, SE of Z
#     (Smith et al. 2012, eq. 2 as in FSA) inflated by an overdispersion
#     factor c-hat = Pearson chi-square / df (floored at 1), first age =
#     peak age + 1 ("peak plus").
#   * Weighted linear regression of log(catch) on age (Maceina & Bettoli
#     1998; FSA::catchCurve(weighted = TRUE)), first age = peak age.
#   * Poisson-lognormal GLMM with a year-class random intercept (Millar 2015),
#     glmmTMB, with a choice of interval (settings$glmm_interval): Wald (ML),
#     REML + t, profile likelihood, or parametric bootstrap (percentile or
#     studentized). Zero catches are added for `zero_tail` ages beyond the
#     oldest fish (the sample contains no older fish: that is information).
#   * Bayesian Poisson-lognormal catch curve in Nimble (settings$bayes).
#
# Multi-year (cohort) catch curve, used when the data hold age samples from
# two or more years (data$long):
#   log E[catch(a, t)] = alpha_t + c_(t-a) - Z a   [ - dZ (t - tbar) a ]
#   alpha_t: fixed effect per sampling year (sample size / effort),
#   c: year-class (cohort) random effect, shared by the same cohort in every
#   year it is sampled. Fitted with glmmTMB and (optionally) in Nimble.
#
# data$df: data.frame(age, catch) with one row per age (zeros filled); for
# multi-year data this is the age composition pooled over years, and
# data$long = data.frame(year, age, catch) (all year x age cells, zeros filled).

cc_selectivity <- function(age, a50, slope) stats::plogis(slope * (age - a50))

CC_ROWS <- c("chapman_robson", "weighted_regression", "glmm", "bayes", "cohort_glmm", "cohort_bayes")
CC_LABELS <- c(chapman_robson = "Chapman-Robson", weighted_regression = "weighted regression",
               glmm = "Poisson-lognormal GLMM", bayes = "Bayesian Poisson-lognormal catch curve",
               cohort_glmm = "multi-year cohort GLMM", cohort_bayes = "multi-year Bayesian cohort catch curve")
CC_SHORT <- c(chapman_robson = "Chapman-Robson", weighted_regression = "Weighted regression", glmm = "GLMM (year-class RE)",
              bayes = "Bayesian (Nimble)", cohort_glmm = "Cohort GLMM (multi-year)", cohort_bayes = "Cohort Bayesian (multi-year)")

# Year-class strength deviations (log scale): lognormal, optional AR(1) over cohorts.
cc_yc_dev <- function(n, sR, phiR) {
  e <- numeric(n); e[1] <- stats::rnorm(1, 0, sR)
  for (i in seq_len(n)[-1]) e[i] <- phiR * e[i - 1] + stats::rnorm(1, 0, sR * sqrt(1 - phiR^2))
  e
}

cc_simulate <- function(params, design, seed) {
  if ((design$cc_years %||% 1) > 1) return(cc_simulate_multi(params, design, seed))
  set.seed(seed)
  amax <- design$max_age %||% 20
  ages <- 0:amax
  Z <- params$Z
  sR <- params$recruit_sd_cc %||% 0
  phiR <- params$recruit_ar1_cc %||% 0
  # Year-class strength R_{year-a}: lognormal with optional AR(1) in time.
  e <- cc_yc_dev(length(ages), sR, phiR)
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
                 Z_chapman_robson = Z, Z_weighted_regression = Z, Z_glmm = Z, Z_bayes = Z))
}

# Multi-year survey: cc_years annual (design key cc_years, because n_years is
# the stock-recruitment series length) age samples of n_fish each from a
# population with lognormal cohort strength. Z may change linearly over the
# sampling years (param Z_trend, per year; Z before the first sampling year =
# Z in year 1); each cohort's survival accumulates the Z of the years it lived
# through (Z at year y = Z + Z_trend (y - ybar), ybar = middle sampling year),
# so the truth for Z is the value at the middle year.
cc_simulate_multi <- function(params, design, seed) {
  set.seed(seed)
  amax <- design$max_age %||% 20; TT <- design$cc_years
  ages <- 0:amax
  Z <- params$Z; dZ <- params$Z_trend %||% 0
  sR <- params$recruit_sd_cc %||% 0; phiR <- params$recruit_ar1_cc %||% 0
  trend <- params$recruit_trend_cc %||% 0
  # cohorts born in years (1 - amax) .. TT; index k = year_born + amax
  nc <- TT + amax
  e <- cc_yc_dev(nc, sR, phiR)
  born <- seq_len(nc) - amax
  R <- exp(e - sR^2 / 2 + trend * (born - mean(born)))
  ybar <- (TT + 1) / 2
  # Z changes linearly during the sampling years only; before year 1 it is Z(1)
  Zy <- function(y) pmax(Z + dZ * (pmin(pmax(y, 1), TT) - ybar), 0)
  sel <- cc_selectivity(ages, params$sel_age50 %||% 2.5, params$sel_slope %||% 4)
  n_fish <- design$n_fish
  long <- do.call(rbind, lapply(seq_len(TT), function(t) {
    N <- vapply(ages, function(a) {
      b <- t - a
      if (a == 0) return(R[b + amax])
      R[b + amax] * exp(-sum(Zy(b:(t - 1))))   # survival through years b .. t-1
    }, numeric(1))
    p <- N * sel; p <- p / sum(p)
    data.frame(year = (design$first_year %||% 2015) + t - 1, age = ages, catch = as.vector(stats::rmultinom(1, n_fish, p)))
  }))
  pos <- long$age[long$catch > 0]
  long <- long[long$age >= min(pos) & long$age <= max(pos), ]
  pooled <- stats::aggregate(catch ~ age, long, sum)
  list(data = list(model = "catch_curve", df = pooled, long = long, units = list(), meta = list(source = "simulated")),
       truth = c(total_mortality_Z = Z, annual_survival = exp(-Z), annual_mortality = 1 - exp(-Z),
                 Z_chapman_robson = Z, Z_weighted_regression = Z, Z_cohort_glmm = Z, Z_cohort_bayes = Z, Z_trend = dZ))
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

# ---- GLMM catch curves (single year and multi-year cohort) ----------------

# d: data.frame(catch, age, yc = factor (year class), year = factor or NULL).
# trend (multi-year only): Z_y = Z + dZ (y - ybar) during the sampling years and
# Z_(first year) before them; a fish of age a sampled in year t has lived through
# years t-a .. t-1, so log survival = -Z a - dZ tc with the exposure covariate
# tc = sum over those years of (clamp(y, first, last) - ybar). Z is then the
# value in the middle sampling year.
cc_glmm_formula <- function(multi, trend) {
  if (!multi) return(catch ~ age + (1 | yc))
  if (trend) catch ~ 0 + year + age + tc + (1 | yc) else catch ~ 0 + year + age + (1 | yc)
}

cc_glmm_once <- function(d, fo, REML = FALSE) {
  m <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(fo, family = stats::poisson, data = d, REML = REML)),
                error = function(e) e)
  if (inherits(m, "error")) return(NULL)
  co <- summary(m)$coefficients$cond
  if (!"age" %in% rownames(co) || !is.finite(co["age", 2])) return(NULL)
  list(m = m, co = co, Z = -co["age", 1], se = co["age", 2],
       dZ = if ("tc" %in% rownames(co)) -co["tc", 1] else NA_real_,
       se_dZ = if ("tc" %in% rownames(co)) co["tc", 2] else NA_real_,
       sd_re = sqrt(as.numeric(glmmTMB::VarCorr(m)$cond$yc)), pd = isTRUE(m$sdr$pdHess))
}

# Interval for Z (and the Z trend) by the chosen method.
#   wald      ML estimate, Z +- 1.96 SE (Millar 2015; the module's original GLMM)
#   reml_t    REML fit, Z +- t(df) SE with df = number of year classes - number of fixed effects
#   profile   ML estimate, profile-likelihood interval (glmmTMB::confint method = "profile")
#   boot_pct  ML estimate, parametric-bootstrap percentile interval
#   boot_t    ML estimate, parametric-bootstrap studentized (bootstrap-t) interval
cc_glmm <- function(d, multi = FALSE, trend = FALSE, interval = "wald", n_boot = 999, seed = 1) {
  eco_require("glmmTMB")
  fo <- cc_glmm_formula(multi, trend)
  g <- cc_glmm_once(d, fo)
  fail <- function(msg) list(Z = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_, ok = FALSE, msg = msg, sd_re = NA_real_,
                             dZ = NA_real_, se_dZ = NA_real_, lower_dZ = NA_real_, upper_dZ = NA_real_, interval = interval)
  if (is.null(g)) return(fail("glmmTMB fit failed"))
  nfix <- nrow(g$co); nyc <- nlevels(droplevels(d$yc)); dft <- max(nyc - nfix, 1)
  z <- stats::qnorm(0.975)
  Z <- g$Z; se <- g$se; dZ <- g$dZ; sedZ <- g$se_dZ
  lo <- Z - z * se; hi <- Z + z * se; lodZ <- dZ - z * sedZ; hidZ <- dZ + z * sedZ
  extra <- list()
  if (interval == "reml_t") {
    r <- cc_glmm_once(d, fo, REML = TRUE)
    if (is.null(r)) return(fail("REML fit failed"))
    tq <- stats::qt(0.975, dft)
    Z <- r$Z; se <- r$se; lo <- Z - tq * se; hi <- Z + tq * se
    dZ <- r$dZ; sedZ <- r$se_dZ; lodZ <- dZ - tq * sedZ; hidZ <- dZ + tq * sedZ
    g$sd_re <- r$sd_re; g$pd <- r$pd
  } else if (interval == "profile") {
    pn <- c("age", if (trend) "tc")
    ci <- tryCatch(suppressWarnings(stats::confint(g$m, parm = pn, method = "profile")), error = function(e) NULL)
    if (is.null(ci) || !all(is.finite(ci[1, 1:2]))) return(fail("profile likelihood failed"))
    lo <- -ci[1, 2]; hi <- -ci[1, 1]
    if (trend) { lodZ <- -ci[2, 2]; hidZ <- -ci[2, 1] }
  } else if (interval %in% c("boot_pct", "boot_t")) {
    # Parametric bootstrap: new year-class effects ~ N(0, sd_re^2) and Poisson
    # catches from the fitted fixed effects, refit by ML, B times.
    X <- stats::model.matrix(stats::update(fo, . ~ . - (1 | yc)), d)
    eta0 <- as.vector(X %*% glmmTMB::fixef(g$m)$cond[colnames(X)])
    yc <- as.integer(droplevels(d$yc)); nyc <- max(yc)
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    bs <- t(vapply(seq_len(n_boot), function(b) {
      db <- d; db$catch <- stats::rpois(nrow(d), exp(eta0 + stats::rnorm(nyc, 0, g$sd_re)[yc]))
      gb <- cc_glmm_once(db, fo)
      if (is.null(gb)) rep(NA_real_, 4) else c(gb$Z, gb$se, gb$dZ, gb$se_dZ)
    }, numeric(4)))
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
    ok <- is.finite(bs[, 1]) & is.finite(bs[, 2])
    if (sum(ok) < 0.8 * n_boot) return(fail(sprintf("only %d of %d bootstrap refits succeeded", sum(ok), n_boot)))
    q <- function(v, p) stats::quantile(v, p, names = FALSE, na.rm = TRUE)
    if (interval == "boot_pct") {
      lo <- q(bs[ok, 1], 0.025); hi <- q(bs[ok, 1], 0.975)
      if (trend) { lodZ <- q(bs[ok, 3], 0.025); hidZ <- q(bs[ok, 3], 0.975) }
    } else {
      ts <- (bs[ok, 1] - g$Z) / bs[ok, 2]
      lo <- g$Z - q(ts, 0.975) * se; hi <- g$Z - q(ts, 0.025) * se
      if (trend) { td <- (bs[ok, 3] - g$dZ) / bs[ok, 4]; lodZ <- g$dZ - q(td, 0.975) * sedZ; hidZ <- g$dZ - q(td, 0.025) * sedZ }
    }
    extra <- list(n_boot_ok = sum(ok), boot_sd = stats::sd(bs[ok, 1]), boot_mean = mean(bs[ok, 1]))
  }
  ok <- isTRUE(g$pd)
  c(list(Z = Z, se = se, lower = lo, upper = hi, ok = ok, msg = if (ok) "" else "Hessian not positive definite",
         sd_re = g$sd_re, dZ = dZ, se_dZ = sedZ, lower_dZ = lodZ, upper_dZ = hidZ, model = g$m, interval = interval,
         df_t = dft, n_yc = nyc), extra)
}

# ---- Bayesian catch curves (Nimble) ---------------------------------------
# log E[catch_i] = a0[year_i] - Z x_i [- dZ tc_i] + sigma u[yc_i], u ~ N(0, 1)
# x = age - mean age used; single-year data have one year and one year class per age.

cc_nimble_code <- function(trend) {
  if (trend) nimble::nimbleCode({
    Z ~ dunif(0, 5)
    dZ ~ dnorm(0, sd = 1)
    sigma ~ T(dnorm(0, sd = sigma_scale), 0, )
    for (j in 1:ny) a0[j] ~ dnorm(0, sd = 10)
    for (c in 1:nc) u[c] ~ dnorm(0, 1)
    for (i in 1:n) {
      log(mu[i]) <- a0[yr[i]] - Z * x[i] - dZ * tc[i] + sigma * u[yc[i]]
      y[i] ~ dpois(mu[i])
    }
  }) else nimble::nimbleCode({
    Z ~ dunif(0, 5)
    sigma ~ T(dnorm(0, sd = sigma_scale), 0, )
    for (j in 1:ny) a0[j] ~ dnorm(0, sd = 10)
    for (c in 1:nc) u[c] ~ dnorm(0, 1)
    for (i in 1:n) {
      log(mu[i]) <- a0[yr[i]] - Z * x[i] + sigma * u[yc[i]]
      y[i] ~ dpois(mu[i])
    }
  })
}

cc_mcmc <- function(settings) {
  fast <- isTRUE(settings$fast)
  list(niter = settings$niter %||% if (fast) 10000 else 40000, nburnin = settings$nburnin %||% if (fast) 2000 else 10000,
       thin = settings$thin %||% if (fast) 4 else 10)
}

cc_bayes <- function(d, trend = FALSE, settings = list()) {
  suppressPackageStartupMessages(library(nimble))
  yr <- as.integer(factor(d$year)); yc <- as.integer(droplevels(d$yc))
  ny <- max(yr); nc <- max(yc); n <- nrow(d)
  x <- d$age - mean(d$age); tc <- if (trend) d$tc else rep(0, n)
  mc <- cc_mcmc(settings)
  y <- d$catch
  a0i <- vapply(seq_len(ny), function(j) log(mean(y[yr == j]) + 0.5), numeric(1))
  inits <- function() c(list(Z = stats::runif(1, 0.2, 1), sigma = stats::runif(1, 0.2, 0.6), a0 = a0i + stats::rnorm(ny, 0, 0.1), u = rep(0, nc)),
                        if (trend) list(dZ = 0))
  pars <- c("Z", if (trend) "dZ", "sigma")
  res <- nimble_run(cc_nimble_code(trend), constants = list(n = n, ny = ny, nc = nc, yr = yr, yc = yc),
                    data = list(y = y, x = x, tc = tc, sigma_scale = settings$prior_yc_sd %||% 1),
                    inits = inits, monitors = c(pars, "a0"),
                    niter = mc$niter, nburnin = mc$nburnin, nchains = settings$nchains %||% 3, thin = mc$thin,
                    # yr / yc are fixed by (n, ny, nc) for a full year x age grid, so they belong in the key
                    cache_key = paste("fs_cc", trend, n, ny, nc, paste(yr, collapse = "."), paste(yc, collapse = "."), sep = "_"),
                    seed = settings$seed %||% 1,
                    configure = function(conf, m) {
                      conf$addSampler(c("Z", "sigma", if (trend) "dZ", "a0"), "AF_slice")
                    })
  all <- stack_draws(res$samples)
  e <- res$estimates
  row <- function(p) {
    v <- all[, p]
    c(estimate = stats::median(v), sd = stats::sd(v), stats::quantile(v, c(0.025, 0.975), names = FALSE),
      rhat = e$rhat[e$param == p], ess = e$ess_bulk[e$param == p])
  }
  rz <- row("Z"); rs <- row("sigma"); rd <- if (trend) row("dZ") else rep(NA_real_, 6)
  rh <- max(e$rhat[e$param %in% pars], na.rm = TRUE); es <- min(e$ess_bulk[e$param %in% pars], na.rm = TRUE)
  list(Z = rz[[1]], se = rz[[2]], lower = rz[[3]], upper = rz[[4]], sd_re = rs[[1]], sd_re_sd = rs[[2]], sd_re_lower = rs[[3]], sd_re_upper = rs[[4]],
       dZ = rd[[1]], se_dZ = rd[[2]], lower_dZ = rd[[3]], upper_dZ = rd[[4]],
       rhat = rh, ess = es, ok = is.finite(rh) && rh < 1.05, samples = res$samples, mcmc = mc,
       msg = sprintf("R-hat %.3f, min bulk ESS %.0f", rh, es))
}

# Data frame for the GLMM / Bayesian fits from first age a0 onwards, with
# zero catches for `tail` ages beyond the oldest fish.
cc_model_frame <- function(age, catch, a0, last, tail) {
  # age_last below the oldest fish: older fish exist, so no zero tail (just truncate)
  if (last < max(age[catch > 0])) tail <- 0
  u <- age >= a0 & age <= last
  a <- age[u]; y <- catch[u]
  if (tail > 0) { a <- c(a, max(a) + seq_len(tail)); y <- c(y, rep(0, tail)) }
  data.frame(catch = y, age = a, yc = factor(a), year = factor(1))
}

cc_model_frame_multi <- function(long, a0, last, tail, trend) {
  yrs <- sort(unique(long$year))
  amax <- max(long$age[long$catch > 0]) + tail
  # age_last below the oldest fish: truncate there, without zero-tail ages (older fish exist)
  if (!is.null(last) && last < max(long$age[long$catch > 0])) amax <- last
  g <- expand.grid(age = seq(a0, amax), year = yrs)
  g$catch <- 0
  key <- paste(long$year, long$age); k2 <- paste(g$year, g$age)
  g$catch[k2 %in% key] <- long$catch[match(k2[k2 %in% key], key)]
  g$cohort <- g$year - g$age
  # exposure to the Z trend: sum over the years lived, t-a .. t-1, of (clamped year - mid year)
  y1 <- min(yrs); yT <- max(yrs); ybar <- (y1 + yT) / 2
  tc <- mapply(function(a, t) if (a == 0) 0 else sum(pmin(pmax((t - a):(t - 1), y1), yT) - ybar), g$age, g$year)
  data.frame(catch = g$catch, age = g$age, yc = factor(g$cohort), year = factor(g$year),
             tc = tc, year_num = g$year)
}

cc_interval_label <- function(m) switch(m, wald = "Wald (ML)", reml_t = "REML with t quantile", profile = "profile likelihood",
                                        boot_pct = "parametric-bootstrap percentile", boot_t = "parametric-bootstrap studentized (bootstrap-t)", m)

# Defaults chosen by planted-truth recovery (VALIDATION.md section 5): the
# Bayesian catch curve was the only single-year estimator that passed the gate
# with margin at the default, worked-example and harder designs; REML + t is the
# best GLMM interval (Wald under-covers, profile and bootstrap did not help),
# but it covered only 0.89 in the independent simulator, so it is not the default.
CC_DEFAULT_PRIMARY <- "bayes"
CC_DEFAULT_INTERVAL <- "reml_t"

cc_fit <- function(data, settings) {
  eco_require("FSA")
  df <- data$df
  df <- df[order(df$age), ]
  multi <- !is.null(data$long) && length(unique(data$long$year)) > 1
  first_cr <- cc_first_age(df, settings$age_full %||% "peak_plus")
  first_rg <- cc_first_age(df, settings$age_full_regression %||% "peak")
  last <- settings$age_last %||% max(df$age)
  use_cr <- df$age >= first_cr & df$age <= last
  use_rg <- df$age >= first_rg & df$age <= last & df$catch > 0
  tail <- settings$zero_tail %||% 5
  notes <- character()
  if (last < max(df$age[df$catch > 0])) {   # older fish exist beyond age_last: truncate, no zero tail
    if (tail > 0) notes <- c(notes, sprintf("age_last = %g is below the oldest fish (age %g): ages above %g were left out and no zero-catch ages were added.", last, max(df$age[df$catch > 0]), last))
    tail <- 0
  }
  interval <- settings$glmm_interval %||% CC_DEFAULT_INTERVAL
  if (!interval %in% c("wald", "reml_t", "profile", "boot_pct", "boot_t"))
    stop("settings: glmm_interval must be one of wald, reml_t, profile, boot_pct, boot_t", call. = FALSE)
  n_boot <- settings$n_boot %||% if (isTRUE(settings$fast)) 299 else 999
  trend <- multi && isTRUE(settings$z_trend)
  prim <- settings$primary %||% if (multi) "cohort_glmm" else CC_DEFAULT_PRIMARY
  do_bayes <- isTRUE(settings$bayes) || prim %in% c("bayes", "cohort_bayes")
  if (prim %in% c("cohort_glmm", "cohort_bayes") && !multi)
    stop(sprintf("settings: primary '%s' needs age samples from two or more years (data: year_col).", prim), call. = FALSE)
  if (prim %in% c("glmm", "bayes") && multi) prim <- sub("^", "cohort_", prim)
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
  zt <- rbind(chapman_robson = c(cr$Z, cr$se, cr$Z - zc * cr$se, cr$Z + zc * cr$se),
              weighted_regression = c(Zr, seR, -ci_r[2], -ci_r[1]))
  gm <- bay <- NULL; conv_bayes <- TRUE; rh <- NA_real_; es <- NA_real_
  if (!multi) {
    d <- cc_model_frame(df$age, df$catch, first_cr, last, tail)
    gm <- cc_glmm(d, FALSE, FALSE, interval, n_boot, settings$seed %||% 1)
    zt <- rbind(zt, glmm = c(gm$Z, gm$se, gm$lower, gm$upper))
    if (do_bayes) { bay <- cc_bayes(d, FALSE, settings); zt <- rbind(zt, bayes = c(bay$Z, bay$se, bay$lower, bay$upper)) }
  } else {
    d <- cc_model_frame_multi(data$long, first_cr, settings$age_last, tail, trend)
    gm <- cc_glmm(d, TRUE, trend, interval, n_boot, settings$seed %||% 1)
    zt <- rbind(zt, cohort_glmm = c(gm$Z, gm$se, gm$lower, gm$upper))
    if (do_bayes) { bay <- cc_bayes(d, trend, settings); zt <- rbind(zt, cohort_bayes = c(bay$Z, bay$se, bay$lower, bay$upper)) }
  }
  colnames(zt) <- c("Z", "se", "lower", "upper")
  if (!is.null(bay)) {
    rh <- bay$rhat; es <- bay$ess; conv_bayes <- bay$ok
    notes <- c(notes, sprintf("Bayesian catch curve (Nimble, %d chains x %d iterations, %d burn-in, thin %d): %s. Priors: Z ~ Uniform(0, 5); year-class SD ~ half-Normal(0, %g); log-scale intercepts ~ Normal(0, 10)%s. Estimate = posterior median, interval = 95%% equal-tailed credible interval.",
                              settings$nchains %||% 3, bay$mcmc$niter, bay$mcmc$nburnin, bay$mcmc$thin, bay$msg, settings$prior_yc_sd %||% 1,
                              if (trend) "; Z trend ~ Normal(0, 1)" else ""))
    if (!bay$ok) notes <- c(notes, "WARNING: the Bayesian catch curve has R-hat >= 1.05; rerun with more iterations (settings niter, nburnin) before using it.")
  }
  if (!prim %in% rownames(zt)) stop(sprintf("settings: primary must be one of %s", paste(rownames(zt), collapse = ", ")), call. = FALSE)
  P <- zt[prim, ]
  if (!all(is.finite(P))) { notes <- c(notes, sprintf("Primary estimator '%s' failed; Chapman-Robson used instead.", prim)); prim <- "chapman_robson"; P <- zt[prim, ] }
  est <- rbind(
    fs_rows("total_mortality_Z", P[1], P[2], P[3], P[4]),
    fs_rows("annual_survival", exp(-P[1]), exp(-P[1]) * P[2], exp(-P[4]), exp(-P[3])),
    fs_rows("annual_mortality", 1 - exp(-P[1]), exp(-P[1]) * P[2], 1 - exp(-P[3]), 1 - exp(-P[4])),
    do.call(rbind, lapply(rownames(zt), function(k) fs_rows(paste0("Z_", k), zt[k, 1], zt[k, 2], zt[k, 3], zt[k, 4]))))
  if (!is.null(bay)) {
    est$rhat[est$param %in% c(paste0("Z_", c("bayes", "cohort_bayes")))] <- rh
    est$ess_bulk[est$param %in% c(paste0("Z_", c("bayes", "cohort_bayes")))] <- es
    if (prim %in% c("bayes", "cohort_bayes")) { est$rhat[1:3] <- rh; est$ess_bulk[1:3] <- es }
  }
  if (trend) {
    src <- if (prim == "cohort_bayes" && !is.null(bay)) bay else gm
    est <- rbind(est, fs_rows("Z_trend", src$dZ, src$se_dZ, src$lower_dZ, src$upper_dZ))
  }
  # year-class SD (log scale): Bayesian posterior when that model was fitted, else the GLMM (no interval)
  est <- rbind(est, if (!is.null(bay)) fs_rows("yearclass_sd", bay$sd_re, bay$sd_re_sd, bay$sd_re_lower, bay$sd_re_upper)
                    else fs_rows("yearclass_sd", gm$sd_re, NA, NA, NA))
  if (!is.null(gm)) {
    if (!gm$ok) notes <- c(notes, paste("Poisson-lognormal GLMM catch curve did not converge cleanly:", gm$msg))
    notes <- c(notes, sprintf("%s: year-class SD (log scale) %.2f; interval for Z: %s%s; %d zero-catch ages added beyond the oldest fish.",
                              if (multi) "Multi-year cohort GLMM" else "GLMM", gm$sd_re, cc_interval_label(gm$interval),
                              if (gm$interval == "reml_t") sprintf(" (t on %d df)", gm$df_t) else if (grepl("^boot", gm$interval)) sprintf(" (%d of %d refits succeeded)", gm$n_boot_ok %||% NA, n_boot) else "",
                              tail))
  }
  if (cr$n < 50) notes <- c(notes, sprintf("Only %d fish at or above the age of full recruitment: intervals will be wide and the age choice matters a lot.", cr$n))
  notes <- c(notes, sprintf("Chapman-Robson used ages %g-%g (n = %d fish; peak age %g + 1); regression used ages %g-%g with nonzero catch (%d ages)%s.",
                            first_cr, last, cr$n, df$age[which.max(df$catch)], first_rg, last, sum(use_rg),
                            if (multi) sprintf("; both on the age composition pooled over %d years", length(unique(data$long$year))) else ""))
  notes <- c(notes, sprintf("Overdispersion c-hat = %.2f (Pearson chi-square %.1f on %d df, P = %.3f)%s.", cr$chat, cr$chisq, cr$df, cr$p_gof,
                            if (cr$chat > 1 && isTRUE(settings$overdispersion %||% TRUE)) "; the Chapman-Robson SE was inflated by sqrt(c-hat)" else ""))
  notes <- c(notes, sprintf("Primary estimator (total_mortality_Z, annual survival and mortality): %s.", CC_LABELS[[prim]]))
  if (abs(Zr - cr$Z) > 1.96 * sqrt(cr$se^2 + seR^2))
    notes <- c(notes, sprintf("Chapman-Robson and regression disagree (Z = %.3f vs %.3f): look for curvature in the catch curve (dome-shaped selectivity, age-varying mortality or a recruitment trend).", cr$Z, Zr))
  conv <- is.finite(P[1]) && is.finite(P[2]) && (conv_bayes || !prim %in% c("bayes", "cohort_bayes"))
  list(estimates = est,
       diagnostics = list(converged = conv, rhat_max = rh, ess_min = es,
                          notes = notes, chat = cr$chat, p_gof = cr$p_gof),
       object = list(model = "catch_curve", df = df, long = data$long, multi = multi, mframe = d, cr = cr, reg = cc,
                     glmm = gm[setdiff(names(gm), "model")], glmm_model = gm$model, bayes = bay, zt = zt, primary = prim,
                     first_cr = first_cr, first_rg = first_rg, last = last, use_cr = use_cr, use_rg = use_rg, trend = trend, tail = tail),
       settings = settings)
}

# ---- data -----------------------------------------------------------------

cc_read_data <- function(config, dc) {
  f <- file.path(config$.dir, dc$file)
  raw <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  ages <- fs_num(fs_col(raw, dc$age_col %||% "age", "age"), dc$age_col %||% "age", "age")
  if (all(is.na(ages))) stop("The age column has no ages.", call. = FALSE)
  nw <- which(!is.na(ages) & ages != round(ages))
  if (length(nw)) stop(sprintf("Column '%s' (age) has %d age(s) that are not whole years, e.g. %s (data rows %s). Catch curves need ages in whole years: check the ageing or round in a copy of the file.",
                               dc$age_col %||% "age", length(nw), paste(utils::head(unique(ages[nw]), 4), collapse = ", "),
                               paste(utils::head(nw, 6), collapse = ", ")), call. = FALSE)
  n_missing_age <- sum(is.na(ages))
  years <- if (!is.null(dc$year_col)) fs_num(fs_col(raw, dc$year_col, "year"), dc$year_col, "year") else NULL
  cnt <- if (!is.null(dc$count_col)) fs_num(fs_col(raw, dc$count_col, "count"), dc$count_col, "count") else rep(1, length(ages))
  drop <- is.na(ages) | (!is.null(years) & is.na(years %||% 0))
  if (!is.null(years) && any(is.na(years))) message(sprintf("%d of %d rows have no year and were left out.", sum(is.na(years) & !is.na(ages)), length(ages)))
  if (sum(is.na(ages)) > 0) message(sprintf("%d of %d rows have no age and were left out.", sum(is.na(ages)), length(ages)))
  if (!is.null(dc$count_col) && any(is.na(cnt[!drop]))) {
    message(sprintf("%d rows have an empty count and were treated as 0.", sum(is.na(cnt[!drop])))); cnt[is.na(cnt)] <- 0
  }
  ages <- ages[!drop]; cnt <- cnt[!drop]; if (!is.null(years)) years <- years[!drop]
  tab <- stats::aggregate(cnt, list(age = ages), sum); names(tab)[2] <- "catch"
  full <- data.frame(age = seq(min(tab$age), max(tab$age)))
  full$catch <- 0
  full$catch[match(tab$age, full$age)] <- tab$catch
  long <- NULL
  if (!is.null(years)) {
    lt <- stats::aggregate(cnt, list(year = years, age = ages), sum); names(lt)[3] <- "catch"
    g <- expand.grid(age = full$age, year = sort(unique(years)))
    g$catch <- 0
    key <- paste(lt$year, lt$age); k2 <- paste(g$year, g$age)
    g$catch[k2 %in% key] <- lt$catch[match(k2[k2 %in% key], key)]
    long <- g[, c("year", "age", "catch")]
    # warnings, not errors: small samples and missing years are allowed
    nyr <- tapply(long$catch, long$year, sum); small <- names(nyr)[nyr < 30]
    if (length(small)) message(sprintf("Year(s) %s have fewer than 30 aged fish (%s): their age compositions are very noisy.",
                                       paste(small, collapse = ", "), paste(nyr[small], collapse = ", ")))
    yrs <- sort(unique(years)); gaps <- setdiff(seq(min(yrs), max(yrs)), yrs)
    if (length(gaps)) message(sprintf("No age sample in year(s) %s. That is allowed (year classes are matched by year - age).", paste(gaps, collapse = ", ")))
    if (length(yrs) == 1) long <- NULL
  }
  list(model = "catch_curve", df = full, long = long, raw_n = nrow(raw),
       n_missing_age = n_missing_age,
       meta = list(file = dc$file, species = dc$species %||% "fish", gear = dc$gear %||% ""))
}

cc_check_data <- function(data) {
  pr <- character(); df <- data$df
  if (any(df$age != round(df$age))) pr <- c(pr, "Ages must be whole years (integers). Round or check the age column.")
  if (any(df$catch < 0, na.rm = TRUE)) pr <- c(pr, "Negative catches found in the catch-at-age table.")
  if (anyNA(df$catch)) pr <- c(pr, "Missing catches in the catch-at-age table: use 0 for ages with no fish.")
  cnt <- if (!is.null(data$long)) data$long$catch else df$catch
  if (any(abs(cnt - round(cnt)) > 1e-8, na.rm = TRUE))
    pr <- c(pr, sprintf("Counts must be whole numbers of fish that were actually aged (found e.g. %s). Catch-at-age expanded with an age-length key or scaled up to the total catch is not a count of aged fish: the Poisson-based estimators would treat every expanded fish as independently aged, and intervals would be far too narrow. Give the number of fish aged at each age (or one row per aged fish).",
                        paste(utils::head(unique(signif(cnt[abs(cnt - round(cnt)) > 1e-8], 4)), 3), collapse = ", ")))
  if (sum(df$catch) < 30) pr <- c(pr, sprintf("Only %d aged fish: a catch curve needs at least ~50-100 aged fish on the descending limb.", sum(df$catch)))
  peak <- df$age[which.max(df$catch)]
  if (sum(df$age > peak + 1) < 2) pr <- c(pr, sprintf("Fewer than 3 ages older than the peak age + 1 (peak = age %g): no descending limb to fit.", peak))
  if (!is.null(data$long)) {
    L <- data$long
    if (any(L$year != round(L$year))) pr <- c(pr, "Years must be whole numbers (e.g. 2021). Check the year column.")
    if (any(tapply(L$catch, L$year, sum) == 0)) pr <- c(pr, "At least one year has no aged fish: remove that year from the file.")
  }
  pr
}

# ---- interpretation, methods, figures --------------------------------------

# Measured 95%-interval coverage under lognormal year-class variability
# (VALIDATION.md section 5): ranges over the single-year designs tried.
CC_COVERAGE <- c(chapman_robson = "about 77-92%", weighted_regression = "about 84-93%", glmm = "89-96%", glmm_wald = "about 80-90%",
                 bayes = "95-98%", cohort_glmm = "94-97%", cohort_glmm_wald = "91-94%", cohort_bayes = "93-97%",
                 pooled = "58-85%")
# Did this estimator's 95% interval pass the coverage gate at every design tried?
# (single-year GLMM + REML-t passed the builder's designs but covered 0.89 in the
# independent simulator, so it is not counted as honest)
cc_honest <- function(prim, interval) prim %in% c("bayes", "cohort_bayes") || (prim == "cohort_glmm" && identical(interval, "reml_t"))
cc_cov_key <- function(prim, interval, multi) {
  if (multi && prim %in% c("chapman_robson", "weighted_regression")) return("pooled")
  if (prim %in% c("glmm", "cohort_glmm") && !identical(interval, "reml_t")) return(paste0(prim, "_wald"))   # Wald / profile / bootstrap: all under-covered
  prim
}

cc_and <- function(x) if (length(x) <= 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])

cc_interpret <- function(fit, data) {
  e <- fit$estimates; g <- function(p) e[e$param == p, ]
  Z <- g("total_mortality_Z"); S <- g("annual_survival"); A <- g("annual_mortality")
  o <- fit$object; zt <- o$zt
  sp <- data$meta$species %||% "fish"
  first <- if (o$primary == "weighted_regression") o$first_rg else o$first_cr
  others <- setdiff(rownames(zt), o$primary)
  agree <- vapply(others, function(k) is.finite(zt[k, 1]) && abs(zt[k, 1] - zt[o$primary, 1]) <= 1.96 * sqrt(zt[k, 2]^2 + zt[o$primary, 2]^2), logical(1))
  nyr <- if (isTRUE(o$multi)) length(unique(o$long$year)) else 1
  out <- c(sprintf("Total instantaneous mortality of fully recruited %s (ages %g and older) is Z = %.2f per year (95%% %s %.2f to %.2f; %s%s). Equivalently, about %.0f%% (%.0f%% to %.0f%%) of fully recruited fish survive each year and %.0f%% (%.0f%% to %.0f%%) die from all causes combined (fishing plus natural mortality).",
            sp, first, Z$estimate, if (o$primary %in% c("bayes", "cohort_bayes")) "credible interval" else "CI", Z$lower, Z$upper, CC_LABELS[[o$primary]],
            if (nyr > 1) sprintf(", %d years of age samples", nyr) else "",
            100 * S$estimate, 100 * S$lower, 100 * S$upper, 100 * A$estimate, 100 * A$lower, 100 * A$upper),
    sprintf("For comparison, the other estimators give %s.%s %s",
            paste(sprintf("Z = %.2f (%.2f to %.2f; %s)", zt[others, 1], zt[others, 3], zt[others, 4], CC_LABELS[others]), collapse = "; "),
            if (any(!vapply(others, cc_honest, logical(1), interval = o$glmm$interval)))
              sprintf(" The intervals of %s were too narrow (below the 95%% coverage limit) in some of this module's simulations with variable year classes (VALIDATION.md), so they are shown only as a cross-check of the estimate.",
                      cc_and(CC_LABELS[others[!vapply(others, cc_honest, logical(1), interval = o$glmm$interval)]])) else "",
            if (all(agree)) "All estimators agree within their uncertainty." else "They do not all agree: the catch curve is probably not a straight line, so check the constant-mortality, constant-recruitment and full-selectivity assumptions."))
  sdyc <- if (!is.null(o$bayes)) o$bayes$sd_re else if (!is.null(o$glmm)) o$glmm$sd_re else NA
  scatter <- if (isTRUE(o$multi))
    sprintf("Following each year class through the %d sampling years, the estimated year-class SD is %.2f on the log scale, i.e. a typical year class is about %.0f%% stronger or %.0f%% weaker than average (pooling the years hides much of this: the pooled age composition has Chapman-Robson c-hat = %.2f)",
            nyr, sdyc, 100 * (exp(sdyc) - 1), 100 * (1 - exp(-sdyc)), o$cr$chat)
  else sprintf("Here the scatter is %s what constant recruitment would give (Chapman-Robson c-hat = %.2f)%s",
            if (o$cr$chat > 1.5) "much larger than" else if (o$cr$chat > 1.1) "somewhat larger than" else "about",
            o$cr$chat,
            if (is.finite(sdyc)) sprintf("; the estimated year-class SD is %.2f on the log scale, i.e. a typical year class is about %.0f%% stronger or %.0f%% weaker than average", sdyc, 100 * (exp(sdyc) - 1), 100 * (1 - exp(-sdyc))) else "")
  out <- c(out, sprintf("Real fish populations have strong and weak year classes, so the counts at age scatter around the straight line. %s. %s",
            scatter,
            if (cc_honest(o$primary, o$glmm$interval))
              sprintf("The %s treats each year class as a random draw and carries that variability into its interval; in this module's simulations with variable year classes its 95%% intervals covered the true Z in %s of samples, as they should.", CC_LABELS[[o$primary]], CC_COVERAGE[[cc_cov_key(o$primary, o$glmm$interval, isTRUE(o$multi))]])
            else sprintf("The %s interval does not fully allow for this: in this module's simulations with variable year classes its 95%% interval covered the true Z in only %s of samples, so read it as closer to a 90%% interval or less. %s gives intervals that held their nominal coverage.",
                         CC_LABELS[[o$primary]], CC_COVERAGE[[cc_cov_key(o$primary, o$glmm$interval, isTRUE(o$multi))]],
                         if (isTRUE(o$multi)) "The multi-year cohort GLMM with glmm_interval: reml_t (the default for several years of samples)" else "The Bayesian catch curve (settings: primary: bayes, the default)")))
  if (isTRUE(o$trend)) {
    tr <- g("Z_trend")
    out <- c(out, sprintf("Z changed by %+.3f per year across the sampled years (95%% interval %.3f to %.3f)%s. The model lets Z change linearly during the sampling years, assumes Z before the first sample equalled Z in the first year, and follows each year class through the years it lived; the Z above is the value in the middle sampling year.",
                          tr$estimate, tr$lower, tr$upper, if (tr$lower > 0 || tr$upper < 0) "; the interval excludes zero" else "; the interval includes zero, so the data do not show a clear change"))
  }
  c(out, paste("A catch curve assumes no trend in recruitment over the year classes sampled, constant mortality over the ages used and full, equal vulnerability of every age from the first age used onwards.",
               "It also assumes the aged fish are a random sample of the population's fully recruited fish.",
               if (!isTRUE(o$multi)) "If the fish came from a few net sets, hauls or sites that each caught a different mix of ages (cluster sampling), the age composition is noisier than the model allows: in this module's verification with 20 hauls the default Bayesian Z was biased +8% (intra-haul correlation 0.1) to +15% (correlation 0.2, which fails the module's 10% bias limit), the GLMM +2% to +8% and Chapman-Robson -2% to -3%, so with such data treat the Bayesian Z as possibly too high and compare it with the GLMM and Chapman-Robson estimates (Nelson 2019)." else "Fish aged from a few net sets or hauls each year (cluster sampling) add noise; in this module's verification (20 hauls a year, intra-haul correlation 0.2) the multi-year cohort GLMM still passed (Z biased +2.5%, interval coverage 93%).",
               "Z mixes fishing and natural mortality; separating them needs an independent estimate of natural mortality M (F = Z - M)."))
}

cc_methods <- function(fit, data, settings) {
  o <- fit$object
  nyr <- if (isTRUE(o$multi)) length(unique(o$long$year)) else 1
  gl <- o$glmm
  txt <- c(sprintf("Total instantaneous mortality (Z) was estimated from catch-at-age data (%d aged fish, ages %g-%g%s) following Smith et al. (2012) and Millar (2015).",
                   sum(o$df$catch), min(o$df$age), max(o$df$age), if (nyr > 1) sprintf(", %d annual samples", nyr) else ""),
           sprintf("(1) Weighted linear regression of log catch on age from the modal age (age %g) onwards%s, excluding ages with zero catch, with weights equal to the predicted log catch of an unweighted fit (Maceina & Bettoli 1998), using FSA::catchCurve (Ogle et al. 2025); 95%% t intervals.",
                   o$first_rg, if (nyr > 1) " on the age composition pooled over years" else ""),
           sprintf("(2) Chapman-Robson (Chapman & Robson 1960) from one year older than the modal age (age %g; the 'peak plus' criterion of Smith et al. 2012), with the Hoenig et al. (1983) bias correction and the standard error of Z (Smith et al. 2012) inflated by sqrt(c-hat), where c-hat (= %.2f) is the Pearson chi-square of the geometric age distribution divided by its degrees of freedom (floored at 1), with t intervals on those degrees of freedom.",
                   o$first_cr, o$cr$chat))
  if (nyr == 1)
    txt <- c(txt, sprintf("(3) A Poisson generalized linear mixed model of catch on age with a year-class random intercept (Millar 2015), fitted with glmmTMB (Brooks et al. 2017) on the same ages as Chapman-Robson plus %d zero-catch ages beyond the oldest fish; %s 95%% interval.",
                          o$tail, cc_interval_label(gl$interval)))
  else
    txt <- c(txt, sprintf("(3) A multi-year (cohort) catch curve: a Poisson generalized linear mixed model with log E[catch at age a in year t] = alpha_t + c_(t-a) - Z a%s, where alpha_t is a fixed effect for each sampling year (sample size) and c_(t-a) ~ Normal(0, s^2) is a random effect for the year class, shared by that year class in every year it was sampled; fitted with glmmTMB (Brooks et al. 2017) to ages %g and older, with zero catches filled in for unobserved age-year cells and %d ages beyond the oldest fish; %s 95%% interval.",
                          if (isTRUE(o$trend)) " - dZ e_(a,t) (e_(a,t) = sum over the years the fish lived of (year - middle sampling year), so Z can change linearly over the sampling years)" else "", o$first_cr, o$tail, cc_interval_label(gl$interval)))
  if (!is.null(o$bayes))
    txt <- c(txt, sprintf("(4) The same %s model fitted in a Bayesian framework with Nimble (de Valpine et al. 2017): Z ~ Uniform(0, 5), year-class SD ~ half-Normal(0, %g), %s ~ Normal(0, 10)%s; %d chains of %d iterations (%d burn-in, thinning %d); posterior median and 95%% equal-tailed credible interval; R-hat %.3f, minimum bulk ESS %.0f.",
                          if (nyr > 1) "cohort" else "Poisson-lognormal", settings$prior_yc_sd %||% 1, if (nyr > 1) "year intercepts" else "intercept (log scale)", if (isTRUE(o$trend)) ", dZ ~ Normal(0, 1)" else "",
                          settings$nchains %||% 3, o$bayes$mcmc$niter, o$bayes$mcmc$nburnin, o$bayes$mcmc$thin, o$bayes$rhat, o$bayes$ess))
  c(txt, sprintf("The %s estimate is reported as the primary estimate; annual survival and mortality are exp(-Z) and 1 - exp(-Z). Bias and interval coverage of each estimator were checked by simulation from an age-structured population with lognormal year-class variability at this study's design (ecotools fisheries-stock module, VALIDATION.md).",
                 CC_LABELS[[o$primary]])) |> paste(collapse = " ")
}

cc_est_labels <- function(opts) c(chapman_robson = fig_label("chapman_robson", "Chapman-Robson", opts),
                                  weighted_regression = fig_label("weighted_regression", "Weighted regression", opts),
                                  glmm = fig_label("glmm", "GLMM (year-class RE)", opts),
                                  bayes = fig_label("bayes", "Bayesian (Nimble)", opts),
                                  cohort_glmm = fig_label("cohort_glmm", "Cohort GLMM", opts),
                                  cohort_bayes = fig_label("cohort_bayes", "Cohort Bayesian", opts))

cc_figures <- function(fit, data, opts) {
  o <- fit$object; df <- o$df
  ci <- opts$ci %||% 0.95; z <- stats::qnorm(1 - (1 - ci) / 2)
  # Estimators are entities: fixed order -> same colour in every figure.
  est_all <- cc_est_labels(opts)
  col_all <- stats::setNames(fs_pal(length(est_all), opts), est_all)
  lty_all <- stats::setNames(c(1, 2, 4, 3, 5, 6), est_all)
  zt <- o$zt
  est_lab <- est_all[rownames(zt)]
  est_col <- col_all[est_lab]; est_lty <- lty_all[est_lab]
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
  lines <- data.frame(x = ages_cr, y = cr$n * (1 - S) * exp(-cr$Z * (ages_cr - a0)), est = est_all[["chapman_robson"]])
  # weighted-regression line and confidence band
  ages_rg <- seq(o$first_rg, o$last, length.out = 100)
  pr <- stats::predict(o$reg$lm, newdata = data.frame(age.e = ages_rg), se.fit = TRUE)
  tq <- stats::qt(1 - (1 - ci) / 2, stats::df.residual(o$reg$lm))
  rg <- data.frame(x = ages_rg, y = exp(pr$fit), lower = exp(pr$fit - tq * pr$se.fit), upper = exp(pr$fit + tq * pr$se.fit),
                   est = est_all[["weighted_regression"]])
  lines <- rbind(lines, rg[, c("x", "y", "est")])
  # Primary GLMM / Bayesian single-year curve: expected catch for an average
  # year class, scaled to the pooled numbers (same Z as reported).
  if (o$primary %in% c("glmm", "bayes")) {
    k <- o$primary; u <- o$use_cr
    Zk <- zt[k, 1]; b <- log(sum(df$catch[u])) - log(sum(exp(-Zk * (df$age[u] - a0))))
    lines <- rbind(lines, data.frame(x = ages_cr, y = exp(b - Zk * (ages_cr - a0)), est = est_all[[k]]))
  }
  lines$est <- factor(lines$est, levels = est_lab)
  lpos <- fs_legend_pos(opts); stacked <- identical(opts$legend %||% "right", "right")
  lpos_est <- if (stacked) "top" else lpos
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = rg, ggplot2::aes(x, ymin = lower, ymax = upper), fill = col_all[[2]], alpha = 0.18) +
    ggplot2::geom_line(data = lines, ggplot2::aes(x, y, colour = est, linetype = est), linewidth = 0.7) +
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
    ggplot2::guides(colour = ggplot2::guide_legend(position = lpos_est, order = 1, nrow = if (stacked) ceiling(nlevels(droplevels(lines$est)) / 2) else NULL, byrow = TRUE, keywidth = ggplot2::unit(1.8, "lines")),
                    linetype = ggplot2::guide_legend(position = lpos_est, order = 1, nrow = if (stacked) ceiling(nlevels(droplevels(lines$est)) / 2) else NULL, byrow = TRUE, keywidth = ggplot2::unit(1.8, "lines")),
                    shape = ggplot2::guide_legend(position = lpos, order = 2, nrow = if (stacked) 2 else NULL, byrow = TRUE))
  abbr <- c(chapman_robson = "CR", weighted_regression = "WR", glmm = "GLMM", bayes = "Bayes", cohort_glmm = "cohort GLMM", cohort_bayes = "cohort Bayes")
  p <- p + ggplot2::labs(title = paste0("Z: ", paste(sprintf("%s %.2f", abbr[rownames(zt)], zt[, 1]), collapse = ", ")))
  p <- fs_size(p, "catch_curve", opts, aspect = 0.95)
  # Intervals at the figure's ci level, built exactly as cc_fit builds the 95%
  # ones. At ci = 0.95 these equal the estimates table; at other levels the
  # bootstrap / profile / Bayesian intervals are drawn from the 95% table
  # (their quantiles at other levels are not stored), which the caption says.
  hw <- c(chapman_robson = if (isTRUE(fit$settings$overdispersion %||% TRUE) && o$cr$chat > 1 && o$cr$df >= 1)
                             stats::qt(1 - (1 - ci) / 2, o$cr$df) else z,
          weighted_regression = tq)
  zz <- data.frame(param = est_lab, estimate = zt[, 1], lower = zt[, 3], upper = zt[, 4], group = factor(est_lab, levels = est_lab))
  for (k in names(hw)) { zz$lower[rownames(zt) == k] <- zt[k, 1] - hw[[k]] * zt[k, 2]; zz$upper[rownames(zt) == k] <- zt[k, 1] + hw[[k]] * zt[k, 2] }
  zz <- zz[is.finite(zz$estimate), ]
  zz$param <- factor(zz$param, levels = rev(est_lab))
  pz <- fig_coef(zz, null_value = zt[o$primary, 1], order = FALSE,
                 xlab = fig_label("Z", "Total mortality, Z (per year)", opts)) +
    ggplot2::scale_colour_manual(values = est_col, guide = "none")
  pz <- fs_size(pz, "mortality_estimates", opts, aspect = 0.5)
  out <- list(
    catch_curve = with_caption(p, sprintf(paste(
      "Catch curve: number of fish caught at each age on a log scale%s (points; filled circles = ages used by Chapman-Robson and the regression,",
      "open squares = used by the regression only, crosses = not fully recruited and excluded). Solid line: Chapman-Robson expected catch (Z = %.2f);",
      "dashed line and band: weighted regression with its %d%% confidence band for expected catch (Z = %.2f).%s The GLMM and Bayesian fits use the same ages as Chapman-Robson. Ages with zero catch are not shown."),
      if (isTRUE(o$multi)) ", pooled over years" else "", zt[1, 1], round(100 * ci), zt[2, 1],
      if (o$primary %in% c("glmm", "bayes")) sprintf(" Line '%s': the primary estimate (Z = %.2f) for an average year class.", est_all[[o$primary]], zt[o$primary, 1]) else "")),
    mortality_estimates = with_caption(pz, sprintf(paste(
      "Total instantaneous mortality Z by estimator (points) with intervals (Chapman-Robson: t interval including the overdispersion correction;",
      "weighted regression: t interval; GLMM: %s; Bayesian: equal-tailed credible interval).%s Dashed line: the primary (%s) estimate. Annual survival is exp(-Z)."),
      cc_interval_label(o$glmm$interval %||% "wald"),
      if (abs(ci - 0.95) > 1e-9) sprintf(" Chapman-Robson and regression intervals are drawn at %d%%; GLMM and Bayesian intervals at 95%%.", round(100 * ci)) else " All intervals are 95%.",
      CC_LABELS[[o$primary]])))
  if (isTRUE(o$multi)) out <- c(out, cc_figures_multi(fit, data, opts, col_all))
  if (!is.null(o$bayes)) {
    tp <- c("Z", if (isTRUE(o$trend)) "dZ", "sigma")
    out$trace <- with_caption(fs_mcmc_fig(fig_trace(lapply(o$bayes$samples, function(m) m[, tp, drop = FALSE])), "trace", opts), sprintf("MCMC traces of the Bayesian catch curve, one colour per chain: Z (total mortality), sigma (year-class SD, log scale)%s. Chains should overlap and show no trend.",
                                                 if (isTRUE(o$trend)) ", dZ (change in Z per year)" else ""))
  }
  out
}

# Multi-year figures: catch at age by sampling year with the fitted cohort
# model, and estimated relative year-class strength.
cc_figures_multi <- function(fit, data, opts, col_all) {
  o <- fit$object; mf <- o$mframe; m <- o$glmm_model
  if (is.null(m)) return(list())
  eco_require("glmmTMB"); suppressPackageStartupMessages(requireNamespace("glmmTMB"))  # predict()/ranef() methods when redrawing from a saved fit
  mf$fit_re <- stats::predict(m, newdata = mf, type = "response", allow.new.levels = TRUE)
  nd <- mf; nd$yc <- NA
  mf$fit_fix <- stats::predict(m, newdata = nd, type = "response", re.form = NA)
  obs <- o$long[o$long$catch > 0, ]
  yrs <- sort(unique(mf$year_num))
  lab_fit <- fig_label("fit_with_yc", "Fitted, with year-class effects", opts)
  lab_fix <- fig_label("fit_average_yc", "Fitted, average year class", opts)
  ld <- rbind(data.frame(age = mf$age, year = mf$year_num, y = mf$fit_re, what = lab_fit),
              data.frame(age = mf$age, year = mf$year_num, y = mf$fit_fix, what = lab_fix))
  ld <- ld[ld$age <= max(obs$age), ]   # not the zero-catch tail ages
  ld$what <- factor(ld$what, levels = c(lab_fit, lab_fix))
  ccol <- col_all[[5]]
  nc_f <- if (length(yrs) <= 4) length(yrs) else if (length(yrs) <= 6) 3 else if (length(yrs) <= 8) 4 else 5
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ld, ggplot2::aes(age, y, linetype = what), colour = ccol, linewidth = 0.6) +
    ggplot2::geom_point(data = obs, ggplot2::aes(age, catch, shape = age >= o$first_cr), colour = eco_pal$ink, size = 1.3) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), labels = c(`TRUE` = fig_label("used_cr", "Used", opts), `FALSE` = fig_label("not_used", "Not fully recruited", opts)), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(1, 2), name = NULL) +
    ggplot2::scale_y_log10(labels = function(b) format(b, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)) +
    ggplot2::facet_wrap(~year, ncol = nc_f) +
    ggplot2::labs(x = fig_label("age", "Age (years)", opts), y = fig_label("log_catch", "Number caught (log scale)", opts)) +
    ggplot2::guides(linetype = ggplot2::guide_legend(position = fs_legend_pos(opts), order = 1),
                    shape = ggplot2::guide_legend(position = fs_legend_pos(opts), order = 2))
  p1 <- fs_size(p1, "catch_by_year", opts, aspect = 0.33 * ceiling(length(yrs) / nc_f) * 5 / nc_f + 0.15, width = "double")
  # year-class effects: conditional modes and SDs
  re <- glmmTMB::ranef(m, condVar = TRUE)$cond$yc
  pv <- attr(re, "condVar")
  sdv <- if (is.list(pv)) vapply(pv, function(v) sqrt(v[1]), numeric(1)) else sqrt(as.vector(pv))
  if (length(sdv) != nrow(re)) sdv <- rep(NA_real_, nrow(re))
  z <- stats::qnorm(1 - (1 - (opts$ci %||% 0.95)) / 2)
  yc <- data.frame(cohort = as.numeric(rownames(re)), est = exp(re[, 1]), lower = exp(re[, 1] - z * sdv), upper = exp(re[, 1] + z * sdv))
  nobs <- tapply(mf$catch, as.character(mf$yc), sum)
  nyrs <- tapply(mf$catch > 0, as.character(mf$yc), sum)
  yc$n_years <- as.vector(nyrs[as.character(yc$cohort)])
  yc <- yc[yc$n_years > 0, ]
  p2 <- ggplot2::ggplot(yc, ggplot2::aes(cohort, est)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = ccol, linewidth = 0.5) +
    ggplot2::geom_point(colour = ccol, size = 1.6) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = fig_label("year_class", "Year class (year of birth)", opts), y = fig_label("yc_strength", "Year-class strength\n(1 = average)", opts))
  p2 <- fs_size(p2, "yearclass_strength", opts, aspect = 0.75)
  list(catch_by_year = with_caption(p1, sprintf(paste(
         "Catch at age in each sampling year (points, log scale; crosses = younger than the first fully recruited age %g). Solid line: the multi-year cohort GLMM's",
         "expected catch including each year class's estimated strength; dashed line: expected catch for an average year class (slope = -Z = %.2f). A strong year class shows up",
         "as points above the dashed line, moving one age to the right each year."), o$first_cr, -o$zt["cohort_glmm", 1])),
       yearclass_strength = with_caption(p2, sprintf(paste(
         "Relative strength of each year class (cohort) from the multi-year cohort GLMM: exp(random effect), 1 = average, with %d%% intervals from the conditional",
         "variance. Year classes seen only as old fish or in few years have wide intervals and are shrunk towards 1. Estimated year-class SD (log scale) = %.2f."),
         round(100 * (opts$ci %||% 0.95)), o$glmm$sd_re)))
}
