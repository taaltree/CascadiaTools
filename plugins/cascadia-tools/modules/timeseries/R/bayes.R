# Bayesian state-space population model in Nimble with covariates on the
# annual growth rate, optional Gompertz density dependence and missing years
# (Kery & Schaub 2012 ch. 5; Dennis et al. 2006; Knape 2008).
#
#   x_t = log N_t
#   x_t = x_{t-1} + r0 + b (x_{t-1} - xbar) + sum_k beta_k z_{k,t} + e_t,  e_t ~ N(0, sd_process^2)
#   y_t = log(count_t) ~ N(x_t, sd_observation^2)     (only years with a count)
#
# xbar = mean log count (centring, so r0 is the growth rate at the average
# abundance and is only weakly correlated with b). Covariates z are
# standardised; the value in the row for year t acts on growth from t-1 to t.
# b = 0 (density_dependence = FALSE) gives exponential growth with
# covariates. Equilibrium (at average covariates): log K = xbar - r0 / b.
#
# Priors: r0 ~ N(0, 1); b ~ U(-2, 0); beta_k ~ N(0, 1) (per SD of covariate);
#   sd_process ~ half-N(prior_sd_process) and sd_observation ~ half-N(prior_sd_observation),
#   each given as c(mean, sd) of a normal truncated at 0 (default c(0, 1)). An
#   informative sd_observation prior (e.g. from replicate counts or a
#   sightability study) is the standard cure for the Knape (2008) problem.
#   x_1 ~ N(first log count, 2).

bayes_code <- function() nimble::nimbleCode({
  r0 ~ dnorm(0, sd = 1)
  if (HAS_DD) {
    b ~ dunif(-2, 0)
  }
  if (HAS_COV) {
    for (k in 1:nc) {
      beta[k] ~ dnorm(0, sd = 1)
    }
    # padding element (Z has an extra zero column) so inprod never sees a
    # length-1 slice, which nimble's C++ generator cannot compile
    beta[ncp] <- 0
  }
  sigma_p ~ T(dnorm(sp_mu, sd = sp_sd), 0, 10)
  sigma_o ~ T(dnorm(so_mu, sd = so_sd), 0, 10)
  x[1] ~ dnorm(x1_mu, sd = 2)
  for (t in 2:nT) {
    if (HAS_DD) {
      if (HAS_COV) {
        mx[t] <- x[t - 1] + r0 + b * (x[t - 1] - xbar) + inprod(beta[1:ncp], Z[t, 1:ncp])
      } else {
        mx[t] <- x[t - 1] + r0 + b * (x[t - 1] - xbar)
      }
    } else {
      if (HAS_COV) {
        mx[t] <- x[t - 1] + r0 + inprod(beta[1:ncp], Z[t, 1:ncp])
      } else {
        mx[t] <- x[t - 1] + r0
      }
    }
    x[t] ~ dnorm(mx[t], sd = sigma_p)
  }
  for (i in 1:nobs) {
    y[i] ~ dnorm(x[tobs[i]], sd = sigma_o)
  }
})

bayes_simulate <- function(params, design, seed) {
  set.seed(seed)
  n <- design$n_years %||% 25
  sp <- params$sd_process; so <- params$sd_observation
  dd <- isTRUE(design$density_dependence %||% TRUE)
  b <- if (dd) params$density_dependence_b else 0
  covs <- design$covariate_names %||% "winter_severity"
  nc <- design$n_covariates %||% length(covs)
  covs <- covs[seq_len(nc)]
  beta <- rep_len(unlist(params$growth_beta %||% 0), nc)
  # raw covariates in "original units" (e.g. snow depth, cm): mean 60, sd 20
  raw <- matrix(stats::rnorm(n * nc, 60, 20), n, nc, dimnames = list(NULL, covs))
  zs <- if (nc) scale(raw) else matrix(0, n, 0)
  lK <- log(params$equilibrium_N %||% 500)
  x <- numeric(n)
  if (dd) {
    sd0 <- sp / sqrt(1 - (1 + b)^2)
    x[1] <- lK + stats::rnorm(1, 0, sd0)
  } else x[1] <- log(design$start_N %||% 500)
  for (t in 2:n) {
    cz <- if (nc) sum(beta * zs[t, ]) else 0
    x[t] <- x[t - 1] + (if (dd) b * (x[t - 1] - lK) else params$growth_rate_mu) + cz + stats::rnorm(1, 0, sp)
  }
  y <- exp(x + stats::rnorm(n, 0, so))
  if (isTRUE(design$round_counts %||% TRUE)) y <- pmax(round(y), 1)
  y[ts_drop_years(n, design$frac_missing)] <- NA
  df <- data.frame(year = ts_years0(design) + seq_len(n) - 1, site = "all", count = y)
  for (k in covs) df[[k]] <- raw[, k]
  truth <- c(if (dd) c(density_dependence_b = b, equilibrium_N = exp(lK)) else c(growth_rate_r0 = params$growth_rate_mu),
             sd_process = sp, sd_observation = so, N_final = exp(x[n]),
             realised_trend_pct_per_year = 100 * (exp((x[n] - x[1]) / (n - 1)) - 1))
  if (nc) {
    # truth per 1 SD of the covariate *as observed in this dataset* (fit standardises with sample SD)
    bsd <- beta * 1  # zs already standardised with sample SD
    truth <- c(truth, stats::setNames(bsd, paste0("growth_beta_", covs)),
               stats::setNames(100 * (exp(bsd) - 1), paste0("growth_effect_pct_", covs)))
  }
  list(data = list(model = "bayes", df = df, covariates = covs, response = "count",
                   meta = list(source = "simulated", count_label = "Count", x_true = x,
                               density_dependence = dd)),
       truth = truth)
}

bayes_mcmc <- function(settings) {
  fast <- isTRUE(settings$fast)
  list(niter = settings$niter %||% if (fast) 12000 else 40000,
       nburnin = settings$nburnin %||% if (fast) 4000 else 10000,
       thin = settings$thin %||% if (fast) 4 else 10,
       nchains = settings$nchains %||% if (fast) 3 else 4)
}

bayes_fit <- function(data, settings) {
  suppressPackageStartupMessages(library(nimble))
  tm <- ts_matrix(data$df)
  y_all <- ts_logY(tm$Y[1, ], settings)
  # trim leading years without counts (growth before the first count is unobserved);
  # trailing years without counts are kept only if forecasting is not intended
  ok <- which(!is.na(y_all)); keep <- min(ok):max(ok)
  yrs <- tm$years[keep]; y_all <- y_all[keep]
  covs <- data$covariates %||% character()
  cv <- ts_cov_matrix(data$df, covs, yrs)
  Z <- cv$Z
  if (length(covs)) { Z[1, ] <- 0; if (anyNA(Z[-1, , drop = FALSE])) stop("covariate missing in a year inside the series; see check_data") }
  dd <- if (!is.null(settings$density_dependence)) isTRUE(settings$density_dependence) else !isFALSE(data$meta$density_dependence)
  nT <- length(yrs); tobs <- which(!is.na(y_all)); nobs <- length(tobs)
  xbar <- mean(y_all, na.rm = TRUE)
  psp <- unlist(settings$prior_sd_process %||% c(0, 1)); pso <- unlist(settings$prior_sd_observation %||% c(0, 1))
  nc <- length(covs)
  const <- list(nT = nT, nobs = nobs, nc = max(nc, 1), ncp = max(nc, 1) + 1, HAS_DD = dd, HAS_COV = nc > 0)
  dat <- list(y = y_all[tobs], tobs = tobs, xbar = xbar, x1_mu = y_all[tobs[1]],
              sp_mu = psp[1], sp_sd = psp[2], so_mu = pso[1], so_sd = pso[2])
  if (nc) dat$Z <- cbind(Z, 0)
  # tobs (which years have counts) is a constant (keeps the latent states
  # conjugate); it is part of the cache key, so a new gap pattern recompiles.
  const$tobs <- tobs; dat$tobs <- NULL
  mc <- bayes_mcmc(settings)
  set.seed(settings$seed %||% 1)
  inits <- function() {
    yi <- stats::approx(tobs, y_all[tobs], xout = seq_len(nT), rule = 2)$y
    l <- list(r0 = stats::rnorm(1, 0, 0.05), sigma_p = stats::runif(1, 0.05, 0.3), sigma_o = stats::runif(1, 0.05, 0.3),
              x = yi + stats::rnorm(nT, 0, 0.05))
    if (dd) l$b <- stats::runif(1, -0.6, -0.1)
    if (nc) l$beta <- c(stats::rnorm(nc, 0, 0.1), 0)
    l
  }
  mon <- c("r0", "sigma_p", "sigma_o", "x", if (dd) "b", if (nc) "beta")
  key <- paste("timeseries_bayes", nT, nobs, nc, dd, paste(tobs, collapse = "."), sep = "_")
  conf_fun <- function(conf, m) {
    nodes <- c("r0", "sigma_p", "sigma_o", if (dd) "b", if (nc) paste0("beta[", seq_len(nc), "]"))
    conf$removeSamplers(nodes)
    conf$addSampler(target = nodes, type = "AF_slice")
    conf$addSampler(target = c("sigma_p", "sigma_o"), type = "RW_block")
  }
  r <- nimble_run(bayes_code(), constants = const, data = dat, inits = inits, monitors = mon,
                  niter = mc$niter, nburnin = mc$nburnin, nchains = mc$nchains, thin = mc$thin,
                  cache_key = key, configure = conf_fun, summary_params = NULL, seed = settings$seed %||% 1)
  # derived quantities per chain
  der <- lapply(r$samples, function(s) {
    s <- as.matrix(s)
    X <- s[, paste0("x[", seq_len(nT), "]"), drop = FALSE]
    out <- cbind(
      if (dd) cbind(density_dependence_b = s[, "b"], equilibrium_N = exp(xbar - s[, "r0"] / s[, "b"]))
      else cbind(growth_rate_r0 = s[, "r0"]),
      sd_process = s[, "sigma_p"], sd_observation = s[, "sigma_o"],
      N_final = exp(X[, nT]),
      realised_trend_pct_per_year = 100 * (exp((X[, nT] - X[, 1]) / (nT - 1)) - 1))
    if (nc) {
        B <- s[, paste0("beta[", seq_len(nc), "]"), drop = FALSE]
      out <- cbind(out, stats::setNames(as.data.frame(B), paste0("growth_beta_", covs)) |> as.matrix(),
                   stats::setNames(as.data.frame(100 * (exp(B) - 1)), paste0("growth_effect_pct_", covs)) |> as.matrix())
    }
    out
  })
  sm <- summarise_draws_eco(der)
  est <- sm$estimates
  # Point estimates are posterior medians: several reported quantities are
  # skewed (SDs near 0, abundance, and equilibrium N, which has a heavy right
  # tail when b is near 0). sd = posterior SD; interval = 2.5-97.5% quantiles.
  allD <- do.call(rbind, der)
  est$estimate <- apply(allD[, est$param, drop = FALSE], 2, stats::median)
  # Convergence over everything reported (incl. all latent states)
  sx <- summarise_draws_eco(r$samples)$estimates
  rh <- max(c(est$rhat, sx$rhat), na.rm = TRUE); ess <- min(c(est$ess_bulk, sx$ess_bulk), na.rm = TRUE)
  S <- stack_draws(r$samples)
  X <- S[, paste0("x[", seq_len(nT), "]"), drop = FALSE]
  q <- function(p) exp(apply(X, 2, stats::quantile, p))
  traj <- data.frame(year = yrs, site = tm$sites[1], estimate = q(0.5), lower = q(0.025), upper = q(0.975),
                     lower50 = q(0.25), upper50 = q(0.75))
  pd <- data.frame(r0 = S[, "r0"], b = if (dd) S[, "b"] else 0, xc = xbar, sd_process = S[, "sigma_p"], xT = X[, nT])
  if (nc) pd$cov_beta <- I(S[, paste0("beta[", seq_len(nc), "]"), drop = FALSE])
  # posterior predictive check: replicate log counts in surveyed years and
  # compare the SD and lag-1 autocorrelation of their year-to-year changes with
  # the data (too little process noise -> observed SD too large; observation
  # error mis-split -> lag-1 autocorrelation off; Auger-Methe et al. 2021).
  stat <- function(v) { d <- diff(v); c(sd_change = stats::sd(d), acf1_change = if (length(d) > 3) stats::cor(d[-1], d[-length(d)]) else NA_real_) }
  set.seed(settings$seed %||% 1)
  ii <- sample.int(nrow(S), min(1000, nrow(S)))
  Trep <- t(vapply(ii, function(i) stat(X[i, tobs] + stats::rnorm(nobs, 0, S[i, "sigma_o"])), numeric(2)))
  Tobs <- stat(y_all[tobs])
  ppc <- list(T_obs = Tobs, T_rep = Trep,
              p = c(sd_change = mean(Trep[, 1] >= Tobs[1]), acf1_change = mean(Trep[, 2] >= Tobs[2], na.rm = TRUE)))
  notes <- character()
  if (rh >= 1.01 || ess < 400) notes <- c(notes, sprintf("Real-data standard (R-hat < 1.01 and bulk ESS >= 400 for every reported quantity and state) not met: max R-hat %.3f, min ESS %.0f. Increase settings: niter.", rh, ess))
  list(estimates = est[, c("param", "estimate", "sd", "lower", "upper", "rhat", "ess_bulk")],
       diagnostics = list(converged = is.finite(rh) && rh < 1.05, rhat_max = rh, ess_min = ess, notes = notes,
                          ppc_p_sd_change = ppc$p[["sd_change"]], ppc_p_acf1_change = ppc$p[["acf1_change"]],
                          t_total = r$diagnostics$t_total),
       object = list(model = "bayes", engine = "Nimble MCMC", samples = r$samples, derived = der, trajectory = traj,
                     years = yrs, sites = tm$sites, Y = tm$Y[, keep, drop = FALSE], pdraws = pd, xbar = xbar,
                     dd = dd, covariates = covs, scaling = cv$scaling, cov_raw = cv$raw, mcmc = mc,
                     priors = list(sd_process = psp, sd_observation = pso), ppc = ppc))
}
