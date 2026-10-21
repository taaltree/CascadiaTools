# Stock-recruitment: Ricker and Beverton-Holt with lognormal error.
#
#   Ricker:          R_t = alpha S_t exp(-beta S_t) exp(eps_t)
#   Beverton-Holt:   R_t = alpha S_t / (1 + beta S_t) exp(eps_t)
#   eps_t ~ N(0, sigma^2) (engine "ls"), or AR(1) with marginal SD sigma and
#   lag-1 autocorrelation phi (engine "nimble", ar1 = TRUE).
#
# Reference points use the "returns" convention (recruits are adult
# returns in the same units as spawners, as in salmon brood tables, or are
# converted to spawner equivalents): replacement line R = S, equilibrium S_eq,
# S_MSY maximises R(S) - S, U_MSY = 1 - S_MSY / R(S_MSY). Ricker S_MSY and
# U_MSY use the exact Lambert-W solution (Scheuerell 2016); Beverton-Holt
# has the closed form S_MSY = (sqrt(alpha) - 1)/beta, U_MSY = 1 - 1/sqrt(alpha).
# Steepness h = R(0.2 S0)/R0 needs unfished spawners per recruit phi0
# (settings$spr0; 1 under the returns convention):
#   Ricker h = 0.2 (alpha phi0)^0.8,  BH h = alpha phi0 / (4 + alpha phi0).
# Reference points are on the median scale (no lognormal bias correction).

sr_mean <- function(form, S, ln_alpha, beta) {
  if (form == "ricker") exp(ln_alpha) * S * exp(-beta * S) else exp(ln_alpha) * S / (1 + beta * S)
}

# Vectorised reference points from parameter draws.
sr_refpoints <- function(form, ln_alpha, beta, spr0 = 1) {
  a <- exp(ln_alpha)
  if (form == "ricker") {
    # alpha <= 1: no surplus, equilibrium and S_MSY are 0 (stock declines to
    # extinction even unfished); beta <= 0 (possible only in LS draws): no
    # density dependence, S_eq and S_MSY are infinite.
    w <- lambert_w0(exp(1 - pmax(ln_alpha, 0)))
    smsy <- ifelse(ln_alpha <= 0, 0, ifelse(beta > 0, (1 - w) / beta, Inf))
    umsy <- ifelse(ln_alpha > 0, 1 - w, 0)
    seq_ <- ifelse(ln_alpha <= 0, 0, ifelse(beta > 0, ln_alpha / beta, Inf))
    h <- 0.2 * (a * spr0)^0.8
    smax <- ifelse(beta > 0, 1 / beta, Inf)
    out <- cbind(spawners_at_MSY = smsy, harvest_rate_at_MSY = umsy, spawners_at_equilibrium = seq_,
                 spawners_at_max_recruitment = smax, steepness = h)
  } else {
    smsy <- ifelse(a <= 1, 0, ifelse(beta > 0, (sqrt(a) - 1) / beta, Inf))
    umsy <- ifelse(a > 1, 1 - 1 / sqrt(a), 0)
    seq_ <- ifelse(a <= 1, 0, ifelse(beta > 0, (a - 1) / beta, Inf))
    h <- a * spr0 / (4 + a * spr0)
    rmax <- ifelse(beta > 0, a / beta, Inf)
    out <- cbind(spawners_at_MSY = smsy, harvest_rate_at_MSY = umsy, spawners_at_equilibrium = seq_,
                 max_recruitment = rmax, steepness = h)
  }
  out
}

sr_simulate <- function(params, design, seed) {
  set.seed(seed)
  n <- design$n_years %||% 30
  form <- design$sr_form %||% "ricker"
  la <- params$sr_ln_alpha; b <- params$sr_beta
  sig <- params$sr_sigma; phi <- params$sr_phi %||% 0
  # Spawners: exogenous, lognormal AR(1) around a mean (escapement varies with
  # harvest policy and implementation error), independent of recruitment
  # residuals.
  Sm <- design$S_mean %||% 60; Scv <- design$S_cv %||% 0.6; Sphi <- design$S_phi %||% 0.5
  ssd <- sqrt(log(1 + Scv^2))
  ls <- numeric(n); ls[1] <- stats::rnorm(1, 0, ssd)
  for (t in 2:n) ls[t] <- Sphi * ls[t - 1] + stats::rnorm(1, 0, ssd * sqrt(1 - Sphi^2))
  S <- Sm * exp(ls - ssd^2 / 2)
  eps <- numeric(n); eps[1] <- stats::rnorm(1, 0, sig)
  for (t in 2:n) eps[t] <- phi * eps[t - 1] + stats::rnorm(1, 0, sig * sqrt(1 - phi^2))
  R <- sr_mean(form, S, la, b) * exp(eps)
  df <- data.frame(year = seq_len(n) + 1990, spawners = S, recruits = R)
  tr <- c(productivity_ln_alpha = la, density_dependence_beta = b, recruitment_sd = sig)
  if (b > 0 && la > 0) {
    rp <- sr_refpoints(form, la, b, design$spr0 %||% 1)[1, ]
    tr <- c(tr, rp[c("spawners_at_MSY", "harvest_rate_at_MSY", "spawners_at_equilibrium", "steepness")])
  }
  if (isTRUE(design$truth_ar1)) tr <- c(tr, recruitment_ar1 = phi)
  list(data = list(model = "stock_recruit", sr_form = form, df = df, spr0 = design$spr0 %||% 1,
                   meta = list(source = "simulated", spawner_units = "thousand fish", recruit_units = "thousand fish")),
       truth = tr)
}

sr_fit <- function(data, settings) {
  form <- settings$sr_form %||% data$sr_form %||% "ricker"
  form <- match.arg(form, c("ricker", "beverton_holt"))
  engine <- settings$sr_engine %||% "ls"
  spr0 <- settings$spr0 %||% data$spr0 %||% 1
  df <- data$df
  if (engine == "nimble") return(sr_fit_nimble(df, form, spr0, settings))
  S <- df$spawners; lR <- log(df$recruits)
  n <- nrow(df)
  notes <- character()
  if (form == "ricker") {
    m <- stats::lm(I(lR - log(S)) ~ S)
    co <- stats::coef(m); V <- stats::vcov(m)
    theta <- c(ln_alpha = unname(co[1]), beta = -unname(co[2]))
    Vt <- V; Vt[1, 2] <- Vt[2, 1] <- -V[1, 2]
    sig <- summary(m)$sigma; dfr <- m$df.residual
    resid <- stats::residuals(m)
    conv <- TRUE
  } else {
    st <- tryCatch({
      mr <- stats::lm(I(lR - log(S)) ~ S); c(la = unname(stats::coef(mr)[1]), lb = log(max(-stats::coef(mr)[2], 1e-3 / mean(S))))
    }, error = function(e) c(la = 1, lb = log(1 / mean(S))))
    fobj <- function(p) sum((lR - (p[1] + log(S) - log1p(exp(p[2]) * S)))^2)
    o <- stats::optim(st, fobj, method = "BFGS", hessian = FALSE, control = list(reltol = 1e-12, maxit = 1000))
    o2 <- stats::optim(o$par, fobj, method = "Nelder-Mead", control = list(reltol = 1e-14, maxit = 5000))
    if (o2$value < o$value) o <- o2
    H <- stats::optimHess(o$par, fobj)
    dfr <- n - 2
    sig <- sqrt(o$value / dfr)
    Vlog <- tryCatch(solve(H / (2 * sig^2)), error = function(e) matrix(NA, 2, 2))
    conv <- o$convergence == 0 && all(is.finite(Vlog)) && all(diag(Vlog) > 0) && sqrt(Vlog[2, 2]) < 5
    if (!conv) notes <- c(notes, "Beverton-Holt fit is poorly determined (density-dependence parameter not identified: the spawner data may not reach the flat part of the curve). Do not use the reference points.")
    theta <- c(ln_alpha = unname(o$par[1]), beta = exp(unname(o$par[2])))
    resid <- lR - (o$par[1] + log(S) - log1p(exp(o$par[2]) * S))
    lb_hat <- o$par[2]
  }
  set.seed(settings$seed %||% 1)
  nd <- settings$n_draws %||% 4000
  if (form == "ricker") {
    # t-scaled multivariate draws (exact for the linear Ricker regression)
    L <- tryCatch(chol(Vt), error = function(e) NULL)
    zz <- matrix(stats::rnorm(nd * 2), nd) %*% (if (is.null(L)) diag(sqrt(pmax(diag(Vt), 0))) else L)
    w <- sqrt(dfr / stats::rchisq(nd, dfr))
    dr <- sweep(zz * w, 2, theta, "+")
    colnames(dr) <- c("ln_alpha", "beta")
  } else {
    if (conv) {
      L <- chol(Vlog)
      zz <- matrix(stats::rnorm(nd * 2), nd) %*% L
      w <- sqrt(dfr / stats::rchisq(nd, dfr))
      dl <- sweep(zz * w, 2, c(theta[1], lb_hat), "+")
      dr <- cbind(ln_alpha = dl[, 1], beta = exp(dl[, 2]))
    } else dr <- cbind(ln_alpha = rep(NA_real_, nd), beta = rep(NA_real_, nd))
  }
  rp_draws <- sr_refpoints(form, dr[, 1], dr[, 2], spr0)
  rp_hat <- sr_refpoints(form, theta[1], theta[2], spr0)[1, ]
  D <- cbind(productivity_ln_alpha = dr[, 1], density_dependence_beta = dr[, 2], rp_draws)
  est <- fs_draw_rows(D, est = c(productivity_ln_alpha = theta[[1]], density_dependence_beta = theta[[2]], rp_hat))
  if (form == "ricker") {   # exact t intervals for the regression coefficients
    tq <- stats::qt(0.975, dfr); se <- sqrt(diag(Vt))
    est$sd[1:2] <- se; est$lower[1:2] <- theta - tq * se; est$upper[1:2] <- theta + tq * se
  }
  sig_ci <- sig * sqrt(dfr / stats::qchisq(c(0.975, 0.025), dfr))
  est <- rbind(est, data.frame(param = "recruitment_sd", estimate = sig, sd = sig / sqrt(2 * dfr),
                               lower = sig_ci[1], upper = sig_ci[2]))
  est$rhat <- NA_real_; est$ess_bulk <- NA_real_
  ac <- if (n > 3) stats::cor(resid[-1], resid[-n]) else NA_real_
  dw <- sum(diff(resid)^2) / sum(resid^2)
  notes <- c(notes, sprintf("%s fitted by least squares on the log scale to %d brood years (residual SD %.2f, %d df).",
                            if (form == "ricker") "Ricker" else "Beverton-Holt", n, sig, dfr))
  if (is.finite(ac) && abs(ac) > 2 / sqrt(n))
    notes <- c(notes, sprintf("Recruitment residuals are autocorrelated (lag-1 r = %.2f; Durbin-Watson %.2f). Least-squares intervals ignore this and are too narrow: refit with settings sr_engine: nimble, ar1: true.", ac, dw))
  if (theta[2] <= 0) notes <- c(notes, "No evidence of density dependence (beta <= 0): equilibrium and MSY reference points are undefined.")
  smax <- if (form == "ricker") 1 / theta[2] else NA
  if (form == "ricker" && is.finite(smax) && max(S) < 0.5 * smax)
    notes <- c(notes, sprintf("Observed spawners (max %.3g) stay well below the spawner level giving maximum recruitment (%.3g): density dependence and S_MSY are extrapolated.", max(S), smax))
  list(estimates = est,
       diagnostics = list(converged = conv, rhat_max = NA_real_, ess_min = NA_real_, notes = notes,
                          resid_ar1 = ac, durbin_watson = dw),
       object = list(model = "stock_recruit", form = form, engine = "ls", df = df, theta = theta, draws = dr,
                     resid = resid, sigma = sig, spr0 = spr0, ar1 = FALSE),
       settings = settings)
}

sr_nimble_code <- function(form) {
  if (form == "ricker") {
    nimble::nimbleCode({
      ln_alpha ~ dnorm(0, sd = 5)
      beta ~ dunif(0, beta_max)
      sigma_nu ~ dunif(0, 3)
      phi ~ dunif(phi_lo, phi_hi)
      for (t in 1:n) mu[t] <- ln_alpha + logS[t] - beta * S[t]
      sigma <- sigma_nu / sqrt(1 - phi^2)          # marginal SD of the residuals
      lR[1] ~ dnorm(mu[1], sd = sigma)
      for (t in 2:n) lR[t] ~ dnorm(mu[t] + phi * (lR[t - 1] - mu[t - 1]), sd = sigma_nu)
    })
  } else {
    nimble::nimbleCode({
      ln_alpha ~ dnorm(0, sd = 5)
      beta ~ dunif(0, beta_max)
      sigma_nu ~ dunif(0, 3)
      phi ~ dunif(phi_lo, phi_hi)
      for (t in 1:n) mu[t] <- ln_alpha + logS[t] - log(1 + beta * S[t])
      sigma <- sigma_nu / sqrt(1 - phi^2)          # marginal SD of the residuals
      lR[1] ~ dnorm(mu[1], sd = sigma)
      for (t in 2:n) lR[t] ~ dnorm(mu[t] + phi * (lR[t - 1] - mu[t - 1]), sd = sigma_nu)
    })
  }
}

sr_mcmc <- function(settings) {
  fast <- isTRUE(settings$fast)
  list(niter = settings$niter %||% if (fast) 10000 else 60000, nburnin = settings$nburnin %||% if (fast) 3000 else 10000,
       thin = settings$thin %||% if (fast) 2 else 10)
}

sr_fit_nimble <- function(df, form, spr0, settings) {
  suppressPackageStartupMessages(library(nimble))  # model building needs nimble attached
  ar1 <- isTRUE(settings$ar1 %||% TRUE)
  n <- nrow(df); S <- df$spawners; lR <- log(df$recruits)
  bmax <- settings$beta_max %||% (20 / stats::median(S))
  inits <- function() {
    mr <- stats::lm(I(lR - log(S)) ~ S)
    list(ln_alpha = unname(stats::coef(mr)[1]) + stats::rnorm(1, 0, 0.1),
         beta = min(max(-unname(stats::coef(mr)[2]), bmax / 100), bmax * 0.9) * exp(stats::rnorm(1, 0, 0.1)),
         sigma_nu = summary(mr)$sigma * exp(stats::rnorm(1, 0, 0.1)), phi = if (ar1) stats::runif(1, -0.2, 0.2) else 0)
  }
  mc <- sr_mcmc(settings)
  set.seed(settings$seed %||% 1)
  res <- nimble_run(sr_nimble_code(form), constants = list(n = n),
                    data = list(lR = lR, S = S, logS = log(S), beta_max = bmax,
                                phi_lo = if (ar1) -0.99 else 0, phi_hi = if (ar1) 0.99 else 1e-8),
                    inits = inits, monitors = c("ln_alpha", "beta", "sigma", "sigma_nu", "phi"),
                    niter = mc$niter, nburnin = mc$nburnin,
                    nchains = settings$nchains %||% 3, thin = mc$thin,
                    cache_key = paste("fs_sr", form, ar1, n, sep = "_"),
                    configure = function(conf, m) {
                      conf$removeSamplers(c("ln_alpha", "beta", "sigma_nu", "phi"))
                      conf$addSampler(c("ln_alpha", "beta"), "AF_slice")
                      conf$addSampler(c("sigma_nu", "phi"), "AF_slice")
                      conf$addSampler(c("ln_alpha", "beta", "sigma_nu", "phi"), "AF_slice")
                    })
  chains <- lapply(res$samples, function(m) {
    rp <- sr_refpoints(form, m[, "ln_alpha"], m[, "beta"], spr0)
    out <- cbind(productivity_ln_alpha = m[, "ln_alpha"], density_dependence_beta = m[, "beta"], rp,
                 recruitment_sd = m[, "sigma"])
    if (ar1) out <- cbind(out, recruitment_ar1 = m[, "phi"], recruitment_innovation_sd = m[, "sigma_nu"])
    out
  })
  # drop columns that are NA in all draws (e.g. no positive draws)
  est <- fs_use_medians(fs_summarise_chains(chains), chains)
  core <- est[est$param %in% c("productivity_ln_alpha", "density_dependence_beta", "recruitment_sd", "recruitment_ar1", "recruitment_innovation_sd"), ]
  rh <- max(core$rhat, na.rm = TRUE); es <- min(core$ess_bulk, na.rm = TRUE)
  allm <- stack_draws(res$samples)
  med <- c(ln_alpha = stats::median(allm[, "ln_alpha"]), beta = stats::median(allm[, "beta"]))
  mu <- if (form == "ricker") med[1] + log(S) - med[2] * S else med[1] + log(S) - log1p(med[2] * S)
  resid <- lR - mu
  notes <- sprintf("Bayesian %s model in Nimble (%d chains x %d iterations, %d burn-in)%s. Priors: ln(alpha) ~ Normal(0, 5); beta ~ Uniform(0, %.3g); innovation SD ~ Uniform(0, 3)%s.",
                   if (form == "ricker") "Ricker" else "Beverton-Holt", settings$nchains %||% 3, mc$niter, mc$nburnin,
                   if (ar1) " with AR(1) recruitment residuals" else "", bmax, if (ar1) "; phi ~ Uniform(-0.99, 0.99)" else "")
  if (bmax < 5 * max(med[2], 1e-12)) notes <- c(notes, "The posterior of beta approaches its prior upper bound: increase settings beta_max.")
  list(estimates = est,
       diagnostics = list(converged = is.finite(rh) && rh < 1.05, rhat_max = rh, ess_min = es, notes = notes,
                          t_total = res$diagnostics$t_total),
       object = list(model = "stock_recruit", form = form, engine = "nimble", df = df, theta = med,
                     draws = cbind(ln_alpha = allm[, "ln_alpha"], beta = allm[, "beta"]), samples = res$samples,
                     resid = resid, sigma = stats::median(allm[, "sigma"]), spr0 = spr0, ar1 = ar1,
                     phi = if (ar1) stats::median(allm[, "phi"]) else 0),
       settings = settings)
}

sr_read_data <- function(config, dc) {
  raw <- utils::read.csv(file.path(config$.dir, dc$file), check.names = FALSE, stringsAsFactors = FALSE)
  yc <- dc$year_col %||% "brood_year"; sc <- dc$spawners_col %||% "spawners"; rc <- dc$recruits_col %||% "recruits"
  df <- data.frame(year = fs_num(fs_col(raw, yc, "year"), yc, "year"),
                   spawners = fs_num(fs_col(raw, sc, "spawners"), sc, "spawners"),
                   recruits = fs_num(fs_col(raw, rc, "recruits"), rc, "recruits"))
  nmiss <- sum(!stats::complete.cases(df))
  if (nmiss) {
    bad <- which(!stats::complete.cases(df))
    lab <- ifelse(is.na(df$year[bad]), paste("data row", bad), as.character(df$year[bad]))
    message(sprintf("%d brood year(s) with a missing year, spawner or recruit value were left out: %s.", nmiss, paste(lab, collapse = ", ")))
  }
  df <- df[stats::complete.cases(df), ]
  df <- df[order(df$year), ]
  list(model = "stock_recruit", sr_form = config$settings$sr_form %||% dc$sr_form %||% "ricker", df = df,
       spr0 = config$settings$spr0 %||% 1, n_incomplete = nmiss,
       meta = list(file = dc$file, species = dc$species %||% "stock", spawner_units = dc$spawner_units %||% "fish",
                   recruit_units = dc$recruit_units %||% dc$spawner_units %||% "fish"))
}

sr_check_data <- function(data) {
  pr <- character(); df <- data$df
  if (nrow(df) < 10) pr <- c(pr, sprintf("Only %d complete brood years: stock-recruitment parameters need at least ~15-20 years with contrast in spawner abundance.", nrow(df)))
  if (any(df$spawners <= 0)) pr <- c(pr, sprintf("Spawners must be > 0 (years: %s).", paste(df$year[df$spawners <= 0], collapse = ", ")))
  if (any(df$recruits <= 0)) pr <- c(pr, sprintf("Recruits must be > 0 for a log-scale fit (years: %s). Zero recruitment needs a different error model.", paste(df$year[df$recruits <= 0], collapse = ", ")))
  if (any(duplicated(df$year))) pr <- c(pr, sprintf("Brood year(s) %s appear more than once: each brood year should be one row (check for a typo in the year column).", paste(unique(df$year[duplicated(df$year)]), collapse = ", ")))
  if (nrow(df) >= 3 && all(df$spawners > 0) && max(df$spawners) / min(df$spawners) < 3)
    pr <- c(pr, sprintf("Spawner abundance only ranges %.1f-fold: without contrast the stock-recruitment curve cannot be estimated.", max(df$spawners) / min(df$spawners)))
  pr
}

sr_interpret <- function(fit, data) {
  e <- fit$estimates; g <- function(p) e[e$param == p, ]
  o <- fit$object; u <- data$meta$spawner_units %||% "fish"
  la <- g("productivity_ln_alpha")
  out <- sprintf("At low spawner abundance each spawner produces about %.1f recruits (alpha = exp(ln alpha) = %.2f, 95%% interval %.2f to %.2f). %s",
                 exp(la$estimate), exp(la$estimate), exp(la$lower), exp(la$upper),
                 if (la$lower > 0) "Productivity is above replacement (alpha > 1), so the stock can sustain some harvest."
                 else "The interval for productivity includes replacement (alpha = 1): sustainable harvest cannot be established.")
  b <- g("density_dependence_beta")
  if (b$lower > 0) {
    smsy <- g("spawners_at_MSY"); umsy <- g("harvest_rate_at_MSY"); seq_ <- g("spawners_at_equilibrium")
    out <- c(out, sprintf("Recruitment is density dependent (beta = %.3g per %s of spawners, 95%% interval %.3g to %.3g). The spawner abundance expected to give maximum sustainable yield is S_MSY = %.3g %s (%.3g to %.3g), harvested at U_MSY = %.0f%% (%.0f%% to %.0f%%); the unfished equilibrium is about %.3g %s.",
                          b$estimate, u, b$lower, b$upper, smsy$estimate, u, smsy$lower, smsy$upper,
                          100 * umsy$estimate, 100 * umsy$lower, 100 * umsy$upper, seq_$estimate, u))
  } else out <- c(out, "The interval for density dependence includes zero, so the data cannot pin down where recruitment levels off; S_MSY and equilibrium abundance are very uncertain or undefined.")
  h <- g("steepness")
  out <- c(out, sprintf("Steepness (recruitment at 20%% of unfished spawners relative to unfished recruitment) is %.2f (%.2f to %.2f)%s.", h$estimate, h$lower, h$upper,
                        if (abs(o$spr0 - 1) < 1e-9) " under the convention that recruits are measured as adult returns comparable to spawners (unfished spawners per recruit = 1)" else sprintf(" with unfished spawners per recruit = %g", o$spr0)))
  sd <- g("recruitment_sd")
  out <- c(out, sprintf("Recruitment varies around the curve with a log-scale SD of %.2f: a typical year is within a factor of %.1f of the curve (2 SD), so single-year recruitment is hard to predict from spawners alone.", sd$estimate, exp(2 * sd$estimate)))
  ar <- g("recruitment_ar1")
  if (nrow(ar)) out <- c(out, sprintf("Recruitment residuals are autocorrelated with lag-1 correlation %.2f (%.2f to %.2f): good and poor years come in runs, consistent with persistent environmental (e.g. ocean) conditions.", ar$estimate, ar$lower, ar$upper))
  c(out, "Reference points assume the stock-recruitment relationship is stationary over the years analysed and that spawner counts are measured without large error; errors in spawner counts bias productivity and density dependence (errors-in-variables) and short series overestimate productivity (time-series bias).")
}

sr_methods <- function(fit, data, settings) {
  o <- fit$object
  sprintf(paste(
    "%s stock-recruitment models (%s) were fitted to %d brood years of spawner and recruit abundance on the log scale with lognormal errors (Hilborn & Walters 1992).",
    "%s",
    "Reference points assume recruits are expressed in spawner units: S_MSY and U_MSY for the Ricker model use the exact Lambert-W solution of Scheuerell (2016); for Beverton-Holt, S_MSY = (sqrt(alpha) - 1)/beta and U_MSY = 1 - 1/sqrt(alpha).",
    "Steepness was computed as 0.2(alpha phi0)^0.8 (Ricker) or alpha phi0/(4 + alpha phi0) (Beverton-Holt) with unfished spawners per recruit phi0 = %g. All reference points are on the median (not mean) recruitment scale."),
    if (o$form == "ricker") "Ricker (1954)" else "Beverton-Holt (1957)",
    if (o$form == "ricker") "R = alpha S exp(-beta S) exp(e)" else "R = alpha S/(1 + beta S) exp(e)",
    nrow(o$df),
    if (o$engine == "ls") "Parameters were estimated by least squares (linear regression of log(R/S) on S for the Ricker model; nonlinear least squares for Beverton-Holt), and intervals for derived quantities were obtained by propagating the multivariate-t sampling distribution of the parameters (4000 draws)."
    else sprintf("The model was fitted in a Bayesian framework with Nimble (de Valpine et al. 2017)%s, priors ln(alpha) ~ N(0, 5^2), beta ~ U(0, %.3g), innovation SD sigma_nu ~ U(0, 3)%s; %d chains of %d iterations after %d burn-in; point estimates are posterior medians with 95%% equal-tailed credible intervals; convergence was checked with R-hat and bulk effective sample size.",
                 if (o$ar1) ", with first-order autoregressive recruitment residuals (e_t = phi e_{t-1} + nu_t, nu_t ~ N(0, sigma_nu^2), stationary start; marginal SD sigma = sigma_nu/sqrt(1 - phi^2))" else "",
                 settings$beta_max %||% (20 / stats::median(o$df$spawners)), if (o$ar1) ", phi ~ U(-0.99, 0.99)" else "",
                 settings$nchains %||% 3, sr_mcmc(settings)$niter - sr_mcmc(settings)$nburnin, sr_mcmc(settings)$nburnin),
    o$spr0)
}

sr_figures <- function(fit, data, opts) {
  o <- fit$object; df <- o$df
  bayes <- o$engine == "nimble"
  ci <- opts$ci %||% 0.95; a <- (1 - ci) / 2; ia <- (1 - fs_inner(opts)) / 2
  iv <- if (bayes) "credible" else "confidence"
  su <- data$meta$spawner_units %||% "fish"; ru <- data$meta$recruit_units %||% su
  col <- fs_pal(2, opts)
  Sg <- seq(0, max(df$spawners) * 1.15, length.out = 120)
  dr <- o$draws[stats::complete.cases(o$draws), , drop = FALSE]
  if (nrow(dr) > 2000) dr <- dr[round(seq(1, nrow(dr), length.out = 2000)), ]
  Rm <- vapply(Sg, function(s) sr_mean(o$form, s, dr[, 1], dr[, 2]), numeric(nrow(dr)))
  qf <- function(M, p) apply(M, 2, stats::quantile, p, na.rm = TRUE)
  pred <- data.frame(x = Sg, estimate = if (bayes) qf(Rm, 0.5) else sr_mean(o$form, Sg, o$theta[1], o$theta[2]),
                     lower = qf(Rm, a), upper = qf(Rm, 1 - a))
  if (bayes) { pred$lower50 <- qf(Rm, ia); pred$upper50 <- qf(Rm, 1 - ia) }
  lab_s <- fig_label("spawners", sprintf("Spawners (%s)", su), opts)
  lab_r <- fig_label("recruits", sprintf("Recruits (%s)", ru), opts)
  ytop <- max(df$recruits, stats::quantile(pred$upper, 0.9)) * 1.08
  p <- ggplot2::ggplot(pred, ggplot2::aes(x, estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = col[1], alpha = 0.15)
  if (bayes) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), fill = col[1], alpha = 0.25)
  p <- p + ggplot2::geom_line(colour = col[1], linewidth = 0.8)
  repl <- abs(o$spr0 - 1) < 1e-9
  if (repl) p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    ggplot2::annotate("text", x = min(max(Sg), ytop) * 0.97, y = min(max(Sg), ytop) * 0.97, label = fig_label("replacement", "R = S", opts),
                      hjust = 1, vjust = -0.4, size = 2.3, colour = eco_pal$muted, angle = 0)
  smsy <- fit$estimates$estimate[fit$estimates$param == "spawners_at_MSY"]
  has_smsy <- length(smsy) && is.finite(smsy) && smsy < max(Sg)
  if (has_smsy) p <- p + ggplot2::geom_vline(xintercept = smsy, linetype = 3, colour = col[2], linewidth = 0.6) +
    ggplot2::annotate("text", x = smsy, y = ytop, label = fig_label("S_MSY_lab", "S[MSY]", opts), parse = TRUE,
                      hjust = -0.15, vjust = 1, size = 2.5, colour = col[2])
  p <- p + ggplot2::geom_point(data = df, ggplot2::aes(spawners, recruits), inherit.aes = FALSE, size = 1.3, colour = eco_pal$ink)
  if (isTRUE(opts$show_data %||% TRUE)) {
    yl <- data.frame(spawners = df$spawners, recruits = df$recruits, lab = sprintf("'%02d", df$year %% 100))
    p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = yl, ggplot2::aes(spawners, recruits, label = lab), inherit.aes = FALSE, size = 2.15,
                               colour = eco_pal$muted, max.overlaps = Inf, box.padding = 0.12, point.padding = 0.05,
                               min.segment.length = 0.3, segment.size = 0.2, segment.colour = "#b3b3b3", seed = 1)
    else ggplot2::geom_text(data = yl, ggplot2::aes(spawners, recruits, label = lab), inherit.aes = FALSE,
                            size = 2.15, vjust = -0.7, colour = eco_pal$muted, check_overlap = TRUE)
  }
  p <- p + ggplot2::labs(x = lab_s, y = lab_r) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::coord_cartesian(ylim = c(0, ytop))
  p <- fs_size(p, "stock_recruitment", opts, aspect = 0.85)
  # residual time series with the interval from parameter uncertainty
  Ry <- vapply(df$spawners, function(s) log(sr_mean(o$form, s, dr[, 1], dr[, 2])), numeric(nrow(dr)))
  res_d <- sweep(-Ry, 2, log(df$recruits), "+")   # draws x years: log R_t - log E[R_t]
  rts <- data.frame(time = df$year,
                    estimate = if (bayes) qf(res_d, 0.5) else log(df$recruits) - log(sr_mean(o$form, df$spawners, o$theta[1], o$theta[2])),
                    lower = qf(res_d, a), upper = qf(res_d, 1 - a), lower50 = qf(res_d, ia), upper50 = qf(res_d, 1 - ia))
  rr <- ggplot2::ggplot(rts, ggplot2::aes(time, estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = eco_pal$null, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), colour = eco_pal$muted, linewidth = 0.35)
  if (bayes) rr <- rr + ggplot2::geom_linerange(ggplot2::aes(ymin = lower50, ymax = upper50), colour = eco_pal$muted, linewidth = 1)
  rr <- rr + ggplot2::geom_line(colour = col[1], linewidth = 0.5) + ggplot2::geom_point(size = 1.2, colour = eco_pal$ink) +
    ggplot2::labs(x = fig_label("year", "Brood year", opts), y = fig_label("resid", "Recruitment residual (log)", opts)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.08)))
  rr <- fs_size(rr, "recruitment_residuals", opts, aspect = 0.7)
  ac <- if (nrow(rts) > 3) stats::cor(rts$estimate[-1], rts$estimate[-nrow(rts)]) else NA
  figs <- list(
    stock_recruitment = with_caption(p, sprintf(paste(
      "Fitted %s stock-recruitment curve (line: %s) with %s; points are brood years labelled by the last two digits of the year.%s%s"),
      if (o$form == "ricker") "Ricker" else "Beverton-Holt",
      if (bayes) "posterior median" else "least-squares estimate",
      if (bayes) sprintf("%d%% and %d%% credible bands", round(100 * ci), round(100 * fs_inner(opts))) else sprintf("%d%% confidence band", round(100 * ci)),
      if (repl) " Dashed line: replacement (recruits = spawners)." else "",
      if (has_smsy) sprintf(" Dotted line: spawners at maximum sustainable yield, S_MSY = %s %s.", fs_sig(smsy), su) else "")),
    recruitment_residuals = with_caption(rr, sprintf(paste(
      "Recruitment residuals log(R / expected R) by brood year (points: %s) with %s from parameter uncertainty;",
      "dashed line: no deviation. Runs of positive or negative values indicate autocorrelation (lag-1 r = %.2f)."),
      if (bayes) "posterior medians" else "estimates",
      if (bayes) sprintf("%d%% (thin) and %d%% (thick) credible intervals", round(100 * ci), round(100 * fs_inner(opts))) else sprintf("%d%% confidence intervals", round(100 * ci)),
      ac)))
  if (bayes) {
    ch <- fs_relabel_chains(o$samples, opts)
    figs$trace <- with_caption(fs_mcmc_fig(fig_trace(ch), "trace", opts), "MCMC traces of the stock-recruitment parameters (one colour per chain); chains should overlap with no trend.")
    figs$posterior <- with_caption(fs_mcmc_fig(fig_posterior(stack_draws(ch)), "posterior", opts), "Posterior densities of the stock-recruitment parameters (natural scale of each parameter).")
  }
  figs
}
