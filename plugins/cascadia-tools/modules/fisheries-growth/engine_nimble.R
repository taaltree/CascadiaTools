# Bayesian engine (Nimble).
#
# group_effects = "random" (hierarchical; Helser & Lai 2004 style):
#   log Linf_g ~ N(mu_logLinf + Linf_beta * z_g, sd_log_Linf)
#   log K_g    ~ N(mu_logK    + K_beta    * z_g, sd_log_K)
#   t0 common to all groups. z_g is the standardised group covariate (optional).
# group_effects = "fixed": independent priors per group, all three parameters
#   group-specific (the Bayesian version of Kimura's general model). Groups
#   are compared through posterior differences.
#
# Likelihood: log L_i ~ N(log VB_i, sigma) (multiplicative) or
#             L_i ~ N(VB_i, sigma) (additive).
# Priors (weakly informative, on scales set by the data):
#   mu_logLinf, logLinf_g (fixed) ~ N(log(max observed length), sd 1)
#   mu_logK, logK_g (fixed) ~ N(log 0.3, sd 1.5)
#   t0 ~ N(0, sd 2) truncated below the youngest age (keeps lengths > 0)
#   sd_log_Linf, sd_log_K ~ half-normal(0, 0.5); betas ~ N(0, sd 1)
#   sigma ~ U(0, 1) (log scale) or U(0, max length) (additive)

fg_nimble_code <- function() nimble::nimbleCode({
  if (RANDOM == 1) {
    mu_lL ~ dnorm(prior_lL, sd = 1)
    mu_lK ~ dnorm(prior_lK, sd = 1.5)
    sd_lL ~ T(dnorm(0, sd = 0.5), 0, )
    sd_lK ~ T(dnorm(0, sd = 0.5), 0, )
    t0c ~ T(dnorm(0, sd = 2), , t0max)
    if (HASCOV == 1) {
      bL ~ dnorm(0, sd = 1)
      bK ~ dnorm(0, sd = 1)
      for (g in 1:G) {
        lL[g] ~ dnorm(mu_lL + bL * x[g], sd = sd_lL)
        lK[g] ~ dnorm(mu_lK + bK * x[g], sd = sd_lK)
      }
    } else {
      for (g in 1:G) {
        lL[g] ~ dnorm(mu_lL, sd = sd_lL)
        lK[g] ~ dnorm(mu_lK, sd = sd_lK)
      }
    }
    for (g in 1:G) { t0[g] <- t0c }
  } else {
    for (g in 1:G) {
      lL[g] ~ dnorm(prior_lL, sd = 1)
      lK[g] ~ dnorm(prior_lK, sd = 1.5)
      t0[g] ~ T(dnorm(0, sd = 2), , t0max_g[g])
    }
  }
  sigma ~ dunif(0, sig_up)
  for (i in 1:N) {
    pred[i] <- exp(lL[gi[i]]) * (1 - exp(-exp(lK[gi[i]]) * (age[i] - t0[gi[i]])))
    if (MULT == 1) {
      y[i] ~ dnorm(log(pred[i]), sd = sigma)
    } else {
      y[i] ~ dnorm(pred[i], sd = sigma)
    }
  }
})

fit_nimble_fg <- function(data, settings) {
  # Definition-time if() blocks in the model code are evaluated by nimble with
  # unqualified calls (getNimbleOption), so nimble must be attached.
  suppressPackageStartupMessages(library(nimble))
  df <- data$df; G <- nlevels(df$group); N <- nrow(df)
  random <- settings$group_effects == "random"
  if (random && G < 3) stop("random group effects need at least 3 groups (5+ recommended); use group_effects: fixed")
  mult <- settings$error == "multiplicative"
  hascov <- random && !is.null(data$cov_name)
  x <- rep(0, G); sc <- NULL
  if (hascov) {
    raw <- data$group_cov[[data$cov_name]][match(levels(df$group), data$group_cov$group)]
    sc <- list(center = mean(raw), scale = stats::sd(raw))
    x <- (raw - sc$center) / sc$scale
  }
  mina <- vapply(split(df$age, df$group), min, 0)
  consts <- list(N = N, G = G, RANDOM = as.integer(random), HASCOV = as.integer(hascov), MULT = as.integer(mult))
  dat <- list(y = if (mult) log(df$length) else df$length, age = df$age, gi = as.integer(df$group),
              x = x, prior_lL = log(max(df$length)), prior_lK = log(0.3),
              t0max = min(df$age) - 0.01, t0max_g = unname(mina) - 0.01,
              sig_up = if (mult) 1 else max(df$length))
  # gi is an index: it must be a constant for Nimble to build the graph
  consts$gi <- dat$gi; dat$gi <- NULL
  if (!hascov) dat$x <- NULL
  if (random) dat$t0max_g <- NULL else dat$t0max <- NULL
  st <- t(vapply(split(df, df$group), fg_starts_one, numeric(3)))
  inits <- function() {
    lL <- log(st[, "Linf"]) + stats::rnorm(G, 0, 0.05); lK <- log(st[, "K"]) + stats::rnorm(G, 0, 0.1)
    t0g <- pmin(st[, "t0"], mina - 0.3)
    ini <- list(lL = unname(lL), lK = unname(lK), sigma = if (mult) 0.1 else 0.05 * max(df$length))
    if (random) {
      ini$mu_lL <- mean(lL); ini$mu_lK <- mean(lK); ini$sd_lL <- 0.1; ini$sd_lK <- 0.2
      ini$t0c <- min(min(t0g), min(df$age) - 0.3)
      if (hascov) { ini$bL <- 0; ini$bK <- 0 }
    } else ini$t0 <- unname(t0g)
    ini
  }
  mons <- c("lL", "lK", "t0", "sigma", if (random) c("mu_lL", "mu_lK", "sd_lL", "sd_lK", "t0c"),
            if (hascov) c("bL", "bK"))
  configure <- function(conf, m) {
    if (random) {
      conf$removeSamplers(c("lL", "lK", "mu_lL", "mu_lK", if (hascov) c("bL", "bK")))
      for (g in seq_len(G)) conf$addSampler(c(sprintf("lL[%d]", g), sprintf("lK[%d]", g)), "AF_slice")
      conf$addSampler(c("mu_lL", "mu_lK", if (hascov) c("bL", "bK")), "RW_block")
      conf$addSampler(c("mu_lL", "mu_lK"), "AF_slice")
      # t0c is tied to every K: also sample it jointly with the mean log K
      conf$addSampler(c("t0c", "mu_lK"), "AF_slice")
      # A shift in t0 is compensated by all K_g together: joint block move
      conf$addSampler(c("t0c", "mu_lK", sprintf("lK[%d]", seq_len(G))), "RW_block")
      if (hascov) { conf$addSampler(c("mu_lK", "bK"), "AF_slice"); conf$addSampler(c("mu_lL", "bL"), "AF_slice") }
    } else {
      conf$removeSamplers(c("lL", "lK", "t0"))
      for (g in seq_len(G))
        conf$addSampler(c(sprintf("lL[%d]", g), sprintf("lK[%d]", g), sprintf("t0[%d]", g)), "AF_slice")
    }
  }
  key <- paste("fisheries-growth", random, hascov, mult, G, N, sep = "_")
  r <- nimble_run(fg_nimble_code(), consts, dat, inits, mons,
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, configure = configure,
                  seed = settings$seed)
  gn <- levels(df$group); sg <- fg_safe(gn)
  derive <- function(M) {
    lL <- M[, sprintf("lL[%d]", seq_len(G)), drop = FALSE]; lK <- M[, sprintf("lK[%d]", seq_len(G)), drop = FALSE]
    t0 <- M[, sprintf("t0[%d]", seq_len(G)), drop = FALSE]
    out <- list()
    for (g in seq_len(G)) {
      out[[paste0("Linf_", sg[g])]] <- exp(lL[, g]); out[[paste0("K_", sg[g])]] <- exp(lK[, g])
      if (!random) out[[paste0("t0_", sg[g])]] <- t0[, g]
      out[[paste0("age_at_90pct_Linf_", sg[g])]] <- t0[, g] + log(10) / exp(lK[, g])
    }
    for (g in seq_len(G)[-1]) {
      out[[paste0("Linf_diff_", sg[g], "_vs_", sg[1])]] <- exp(lL[, g]) - exp(lL[, 1])
      out[[paste0("K_diff_", sg[g], "_vs_", sg[1])]] <- exp(lK[, g]) - exp(lK[, 1])
      if (!random) out[[paste0("t0_diff_", sg[g], "_vs_", sg[1])]] <- t0[, g] - t0[, 1]
    }
    if (random) {
      out$Linf_typical <- exp(M[, "mu_lL"]); out$K_typical <- exp(M[, "mu_lK"]); out$t0 <- M[, "t0c"]
      out$sd_log_Linf <- M[, "sd_lL"]; out$sd_log_K <- M[, "sd_lK"]
      if (hascov) {
        out[[paste0("Linf_beta_", data$cov_name)]] <- M[, "bL"]
        out[[paste0("K_beta_", data$cov_name)]] <- M[, "bK"]
      }
    }
    out$residual_sd <- M[, "sigma"]
    do.call(cbind, out)
  }
  dsamp <- lapply(r$samples, derive)
  s <- summarise_draws_eco(dsamp)
  # Among-group SDs have right-skewed posteriors when there are few groups, so
  # their posterior mean sits above the bulk; report the posterior median
  # (the 95% interval is unchanged). See VALIDATION.md.
  for (p in intersect(c("sd_log_Linf", "sd_log_K"), s$estimates$param))
    s$estimates$estimate[s$estimates$param == p] <- stats::median(unlist(lapply(dsamp, function(d) d[, p])))
  s$diagnostics$t_setup <- r$diagnostics$t_setup; s$diagnostics$t_total <- r$diagnostics$t_total
  s$diagnostics$notes <- character()
  list(estimates = s$estimates, diagnostics = s$diagnostics, samples = r$samples,
       derived = dsamp, cov_scaling = sc, x = x, random = random, hascov = hascov)
}

# Posterior curves per group: median-of-posterior curve with credible band
# and a posterior predictive interval for single fish.
fg_curve_nimble <- function(fit, df, error, ci = 0.95, inner = 0.5, ndraw = 600) {
  M <- stack_draws(fit$samples); G <- nlevels(df$group); gn <- levels(df$group)
  set.seed(1); M <- M[sample.int(nrow(M), min(ndraw, nrow(M))), , drop = FALSE]
  a2 <- (1 - ci) / 2; i2 <- (1 - inner) / 2
  do.call(rbind, lapply(seq_len(G), function(g) {
    a <- seq(max(0, min(df$age) - 0.5), max(df$age), length.out = 80)
    t0 <- M[, sprintf("t0[%d]", g)]; a <- a[a > max(t0) + 0.02]
    P <- outer(exp(M[, sprintf("lL[%d]", g)]), rep(1, length(a))) *
      (1 - exp(-outer(exp(M[, sprintf("lK[%d]", g)]), rep(1, length(a))) * outer(-t0, a, "+")))
    Y <- if (error == "multiplicative") P * exp(matrix(stats::rnorm(length(P)), nrow(P)) * M[, "sigma"])
         else P + matrix(stats::rnorm(length(P)), nrow(P)) * M[, "sigma"]
    q <- function(X, p) apply(X, 2, stats::quantile, p)
    data.frame(group = gn[g], x = a, estimate = q(P, 0.5), lower = q(P, a2), upper = q(P, 1 - a2),
               lower50 = q(P, i2), upper50 = q(P, 1 - i2), pi_lower = q(Y, a2), pi_upper = q(Y, 1 - a2))
  }))
}
