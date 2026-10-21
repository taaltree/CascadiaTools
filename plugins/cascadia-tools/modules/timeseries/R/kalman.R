# Univariate state-space models fitted by exact maximum likelihood with a
# hand-written Kalman filter (Dennis et al. 2006; Humbert et al. 2009).
#
#   state (log abundance):  x_t = x_{t-1} + mu + e_t                  (exponential, EGSS)
#                           x_t = x_{t-1} + b (x_{t-1} - log K) + e_t (Gompertz, GSS)
#                           e_t ~ N(0, sd_process^2)
#   observation:            y_t = log(count_t) = x_t + f_t,  f_t ~ N(0, sd_observation^2)
#
# Years without a count are skipped by the filter (state predicted forward),
# so gaps are handled exactly. The initial state x_0 (one year BEFORE the
# first count) is a free parameter (MARSS x0 = "unequal", tinitx = 0), so the
# population need not start at equilibrium. Putting the free state at t = 0
# rather than t = 1 matters: with x_1 free and fixed (variance 0) the first
# count can be fitted exactly and the likelihood is unbounded as
# sd_observation -> 0, which drags the MLE to a spurious sd_observation = 0
# (found in the cross-check against MARSS, 2026-10-19). observation_error = FALSE fits the process-noise-only model
# (the classic AR(1)/diffusion regression on log counts), used as the
# negative control: it attributes observation noise to dynamics.

# Kalman filter + RTS smoother for x_t = cc * x_{t-1} + dd + e_t. y may have NA.
# Initial state x_0 = `x0` (variance 0) one step before y[1], so
# x_1 ~ N(cc x0 + dd, sp^2). Returns log-likelihood and smoothed states.
kf_uni <- function(y, cc, dd, sp, so, x0, smooth = FALSE) {
  n <- length(y)
  ap <- Pp <- af <- Pf <- numeric(n)
  ll <- 0
  for (t in seq_len(n)) {
    if (t == 1) { a <- cc * x0 + dd; P <- sp^2 } else { a <- cc * af[t - 1] + dd; P <- cc^2 * Pf[t - 1] + sp^2 }
    ap[t] <- a; Pp[t] <- P
    if (!is.na(y[t])) {
      Fv <- P + so^2
      v <- y[t] - a
      ll <- ll - 0.5 * (log(2 * pi * Fv) + v^2 / Fv)
      K <- P / Fv
      af[t] <- a + K * v; Pf[t] <- P * (1 - K)
    } else { af[t] <- a; Pf[t] <- P }
  }
  out <- list(ll = ll)
  if (smooth) {
    xs <- af; Ps <- Pf
    if (n > 1) for (t in (n - 1):1) {
      J <- if (Pp[t + 1] > 0) Pf[t] * cc / Pp[t + 1] else 0
      xs[t] <- af[t] + J * (xs[t + 1] - ap[t + 1])
      Ps[t] <- Pf[t] + J^2 * (Ps[t + 1] - Pp[t + 1])
    }
    out$xs <- xs; out$Ps <- pmax(Ps, 0)
  }
  out
}

# Unpack a parameter vector into the filter's cc, dd, sp, so, x0.
kf_unpack <- function(th, form, obs_error) {
  sp <- exp(th[["log_sd_process"]])
  so <- if (obs_error) exp(th[["log_sd_observation"]]) else 0
  if (form == "exponential") list(cc = 1, dd = th[["mu"]], sp = sp, so = so, x0 = th[["x0"]])
  else list(cc = 1 + th[["b"]], dd = -th[["b"]] * th[["logK"]], sp = sp, so = so, x0 = th[["x0"]])
}

kf_negll <- function(th, y, form, obs_error, nm) {
  names(th) <- nm
  p <- kf_unpack(th, form, obs_error)
  -kf_uni(y, p$cc, p$dd, p$sp, p$so, p$x0)$ll
}

# Likelihood for obs_error = FALSE (process noise only): condition on the
# first count (x_1 = y_1) and use exact multi-year transitions across gaps.
kf_negll_noobs <- function(th, y, form, nm) {
  names(th) <- nm
  sp <- exp(th[["log_sd_process"]])
  cc <- if (form == "exponential") 1 else 1 + th[["b"]]
  dd <- if (form == "exponential") th[["mu"]] else -th[["b"]] * th[["logK"]]
  n <- length(y); ll <- 0; a <- y[1]; P <- 0
  for (t in 2:n) {
    a <- cc * a + dd; P <- cc^2 * P + sp^2
    if (!is.na(y[t])) { ll <- ll + stats::dnorm(y[t], a, sqrt(P), log = TRUE); a <- y[t]; P <- 0 }
  }
  -ll
}

kf_fit_uni <- function(y, form = c("exponential", "gompertz"), obs_error = TRUE, n_draws = 1000,
                       seed = 1, prob = 0.95, profile = TRUE) {
  form <- match.arg(form)
  # trim leading/trailing years without counts
  ok <- which(!is.na(y)); y <- y[min(ok):max(ok)]
  n <- length(y)
  d <- diff(y); vd <- stats::var(d, na.rm = TRUE); md <- mean(d, na.rm = TRUE)
  if (!is.finite(vd) || vd <= 0) vd <- 0.05
  ly <- mean(y, na.rm = TRUE)
  nm <- c(if (form == "exponential") "mu" else c("b", "logK"), "log_sd_process",
          if (obs_error) "log_sd_observation", if (obs_error) "x0")
  lo <- c(if (form == "exponential") -2 else c(-1.99, ly - 15), -8, if (obs_error) -8, if (obs_error) min(y, na.rm = TRUE) - 5)
  hi <- c(if (form == "exponential") 2 else c(0.99, ly + 15), 3, if (obs_error) 3, if (obs_error) max(y, na.rm = TRUE) + 5)
  starts <- list()
  for (f in c(0.2, 0.5, 0.8)) for (b0 in if (form == "gompertz") c(-0.1, -0.5, -1) else NA) {
    sp0 <- sqrt(f * vd); so0 <- sqrt((1 - f) * vd / 2)
    s <- c(if (form == "exponential") md else c(b0, ly), log(sp0), if (obs_error) log(so0),
           if (obs_error) (if (form == "exponential") y[1] - md else y[1]))
    if (!obs_error && f != 0.5) next
    starts[[length(starts) + 1]] <- s
  }
  obj <- if (obs_error) function(th) kf_negll(th, y, form, TRUE, nm) else function(th) kf_negll_noobs(th, y, form, nm)
  if (is.na(y[1]) && !obs_error) stop("first year must have a count")
  best <- NULL
  for (s in starts) {
    s <- pmin(pmax(s, lo + 1e-3), hi - 1e-3)
    o <- tryCatch(stats::nlminb(s, obj, lower = lo, upper = hi, control = list(iter.max = 500, eval.max = 1000)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$objective) && (is.null(best) || o$objective < best$objective)) best <- o
  }
  if (is.null(best)) stop("likelihood optimisation failed from every start")
  th <- stats::setNames(best$par, nm)
  notes <- character()
  at_bound <- nm[abs(th - lo) < 1e-3 | abs(th - hi) < 1e-3]
  free <- setdiff(nm, at_bound)
  # Hessian of the free parameters (boundary ones held fixed)
  objf <- function(v) { t2 <- th; t2[free] <- v; obj(t2) }
  H <- tryCatch(stats::optimHess(th[free], objf), error = function(e) NULL)
  Vf <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
  pd <- !is.null(Vf) && all(is.finite(Vf)) && all(diag(Vf) > 0)
  V <- NULL
  if (pd) { V <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm)); V[free, free] <- Vf }
  if (length(at_bound)) notes <- c(notes, sprintf("Maximum-likelihood estimate on the boundary for %s (a standard deviation estimated at ~0): the data cannot separate process from observation error well; the profile-likelihood interval for it starts at 0.",
                                                 paste(at_bound, collapse = ", ")))
  if (!pd) notes <- c(notes, "The Hessian is not invertible: the likelihood is flat in some direction (non-identifiable parameters).")
  converged <- best$convergence == 0 && pd
  kf <- list(theta = th, V = V, y = y, n = n, form = form, obs_error = obs_error,
             loglik = -best$objective, converged = converged, notes = notes, nm = nm, at_bound = at_bound,
             aic = 2 * best$objective + 2 * length(nm), obj = obj, lo = lo, hi = hi)
  if (!profile) return(kf)
  kf$profile <- kf_profile_all(kf, c(if (form == "exponential") "mu" else c("b", "logK"), "log_sd_process",
                                     if (obs_error) "log_sd_observation"), prob = prob)
  kf
}

# Profile-likelihood interval for one parameter: the set of values whose
# profile log-likelihood is within qchisq(prob, 1)/2 of the maximum. Works at
# boundaries (a variance estimated at 0 gets a lower limit at the bound).
kf_profile_ci <- function(kf, p, prob = 0.95) {
  th <- kf$theta; nm <- kf$nm; j <- match(p, nm)
  crit <- stats::qchisq(prob, 1) / 2
  target <- kf$loglik - crit
  others <- setdiff(seq_along(nm), j)
  warm <- th[others]
  prof <- function(v) {
    f <- function(o) { t2 <- th; t2[others] <- o; t2[j] <- v; kf$obj(t2) }
    best <- Inf; bp <- warm
    for (st in list(warm, th[others])) {
      o <- tryCatch(stats::nlminb(pmin(pmax(st, kf$lo[others] + 1e-6), kf$hi[others] - 1e-6), f,
                                  lower = kf$lo[others], upper = kf$hi[others], control = list(iter.max = 300, eval.max = 600)),
                    error = function(e) NULL)
      if (!is.null(o) && is.finite(o$objective) && o$objective < best) { best <- o$objective; bp <- o$par }
    }
    warm <<- bp
    -best
  }
  se0 <- if (!is.null(kf$V) && kf$V[j, j] > 0) sqrt(kf$V[j, j]) else 0.1 * (kf$hi[j] - kf$lo[j]) / 10
  side <- function(dir) {
    bound <- if (dir < 0) kf$lo[j] else kf$hi[j]
    warm <<- th[others]
    inside <- th[[j]]; step <- max(se0, 1e-3)
    for (k in 1:30) {
      v <- th[[j]] + dir * step
      if ((dir < 0 && v <= bound) || (dir > 0 && v >= bound)) {
        v <- bound
        if (prof(v) >= target) return(bound)
      }
      lp <- prof(v)
      if (lp < target) {
        r <- tryCatch(stats::uniroot(function(u) prof(u) - target, sort(c(inside, v)), tol = 1e-4)$root,
                      error = function(e) (inside + v) / 2)
        return(r)
      }
      if (v == bound) return(bound)
      inside <- v; step <- step * 1.8
    }
    bound
  }
  c(lower = side(-1), upper = side(1))
}

kf_profile_all <- function(kf, params, prob = 0.95) {
  out <- lapply(params, function(p) kf_profile_ci(kf, p, prob))
  names(out) <- params
  out
}

# Derived quantities + trajectory, propagating parameter uncertainty by
# drawing parameters from their asymptotic normal distribution and then
# states from the smoother given each draw.
kf_summaries <- function(kf, n_draws = 1000, seed = 1, prob = 0.95) {
  set.seed(seed)
  th <- kf$theta; z <- stats::qnorm(1 - (1 - prob) / 2)
  se <- if (!is.null(kf$V)) sqrt(diag(kf$V)) else stats::setNames(rep(NA_real_, length(th)), names(th))
  rows <- list()
  pr <- kf$profile
  row <- function(name, p, tf = identity, dtf = function(v) 1) {
    ci <- pr[[p]]; atb <- p %in% kf$at_bound
    lo <- tf(ci[["lower"]]); hi <- tf(ci[["upper"]])
    # a log-SD profiled down to its bound means the interval includes 0
    if (grepl("^log_sd", p) && abs(ci[["lower"]] - kf$lo[match(p, kf$nm)]) < 1e-6) lo <- 0
    if (atb && grepl("^log_sd", p)) est <- 0 else est <- tf(th[[p]])
    ts_rows(name, est, abs(dtf(th[[p]])) * se[[p]], min(lo, hi), max(lo, hi))
  }
  if (kf$form == "exponential") {
    rows$mu <- row("growth_rate_mu", "mu")
    rows$tr <- row("trend_pct_per_year", "mu", function(v) 100 * (exp(v) - 1), function(v) 100 * exp(v))
  } else {
    rows$b <- row("density_dependence_b", "b")
    rows$K <- row("equilibrium_N", "logK", exp, exp)
  }
  rows$sp <- row("sd_process", "log_sd_process", exp, exp)
  if (kf$obs_error) rows$so <- row("sd_observation", "log_sd_observation", exp, exp)
  # states
  sm0 <- if (kf$obs_error) {
    p0 <- kf_unpack(th, kf$form, kf$obs_error)
    kf_uni(kf$y, p0$cc, p0$dd, p0$sp, p0$so, p0$x0, smooth = TRUE)
  } else list(xs = kf$y, Ps = rep(0, kf$n))
  X <- matrix(NA_real_, n_draws, kf$n)
  TH <- matrix(th, n_draws, length(th), byrow = TRUE, dimnames = list(NULL, names(th)))
  if (kf$obs_error) {
    # Parameter uncertainty by sampling-importance-resampling from the
    # likelihood (flat priors on the SDs, b, log K / mu and x0; kf_sir). The
    # earlier asymptotic-normal draws collapsed whenever an SD sat on its
    # boundary (sd_observation = 0 -> state = count, zero width) and exploded
    # when b was drawn near -2: Gompertz N_final coverage was 0.57 (VALIDATION.md).
    s <- kf_sir(kf, n_draws, seed = seed)
    TH <- s$draws; kf$sir_ess <- s$ess
    for (i in seq_len(n_draws)) {
      pd <- kf_unpack(TH[i, ], kf$form, TRUE)
      sm <- kf_uni(kf$y, pd$cc, pd$dd, pd$sp, pd$so, pd$x0, smooth = TRUE)
      X[i, ] <- sm$xs + sqrt(sm$Ps) * stats::rnorm(kf$n)
    }
  } else if (!is.null(kf$V)) {
    E <- eigen(kf$V, symmetric = TRUE)
    L <- E$vectors %*% diag(sqrt(pmax(E$values, 0)), length(th))
    for (i in seq_len(n_draws)) {
      thd <- pmin(pmax(th + as.vector(L %*% stats::rnorm(length(th))), kf$lo), kf$hi)
      names(thd) <- names(th); TH[i, ] <- thd
      X[i, ] <- sm0$xs
    }
  } else {
    for (i in seq_len(n_draws)) X[i, ] <- sm0$xs + sqrt(sm0$Ps) * stats::rnorm(kf$n)
  }
  a <- (1 - prob) / 2
  traj <- data.frame(t = seq_len(kf$n), estimate = exp(sm0$xs),
                     lower = exp(apply(X, 2, stats::quantile, a, na.rm = TRUE)),
                     upper = exp(apply(X, 2, stats::quantile, 1 - a, na.rm = TRUE)),
                     lower50 = exp(apply(X, 2, stats::quantile, 0.25, na.rm = TRUE)),
                     upper50 = exp(apply(X, 2, stats::quantile, 0.75, na.rm = TRUE)))
  nf <- exp(X[, kf$n])
  rows$nf <- ts_rows("N_final", exp(sm0$xs[kf$n]), stats::sd(nf, na.rm = TRUE),
                     stats::quantile(nf, a, na.rm = TRUE, names = FALSE), stats::quantile(nf, 1 - a, na.rm = TRUE, names = FALSE))
  est <- do.call(rbind, rows); rownames(est) <- NULL
  list(estimates = est, trajectory = traj, state_draws = X, theta_draws = TH, sir_ess = kf$sir_ess)
}

# Kalman log-likelihood for many parameter vectors at once (vectors cc, dd,
# sp, so, x0 of equal length); same recursion as kf_uni.
kf_uni_vec <- function(y, cc, dd, sp, so, x0) {
  a <- cc * x0 + dd; P <- sp^2; ll <- 0
  for (t in seq_along(y)) {
    if (t > 1) { a <- cc * a + dd; P <- cc^2 * P + sp^2 }
    if (!is.na(y[t])) {
      Fv <- P + so^2; v <- y[t] - a
      ll <- ll - 0.5 * (log(2 * pi * Fv) + v^2 / Fv)
      K <- P / Fv; a <- a + K * v; P <- P * (1 - K)
    }
  }
  ll
}

# Sampling-importance-resampling from the likelihood x flat priors on
# (mu or b, log K, sd_process, sd_observation, x0) within the optimiser box.
# Proposal: multivariate t (5 df) centred at the MLE with the Hessian
# correlation and SDs inflated 1.5x; SDs on a boundary get a wide proposal
# reaching up to their profile upper limit. Flat on the SD scale = Jacobian
# exp(log_sd) on the log scale. Returns n_draws resampled parameter vectors
# and the importance-sampling effective sample size.
kf_sir <- function(kf, n_draws, n_prop = 20000, df = 5, seed = 1) {
  set.seed(seed + 7)
  th <- kf$theta; nm <- kf$nm; p <- length(th)
  sdv <- if (!is.null(kf$V)) sqrt(pmax(diag(kf$V), 0)) else rep(0, p); names(sdv) <- nm
  ctr <- th
  C <- diag(p)
  if (!is.null(kf$V)) { s0 <- ifelse(sdv > 0, sdv, 1); C <- kf$V / outer(s0, s0); diag(C) <- 1 }
  for (j in seq_len(p)) {
    pr <- kf$profile[[nm[j]]]
    w <- if (!is.null(pr)) (pr[["upper"]] - pr[["lower"]]) / 3.92 else NA_real_
    if (!is.finite(sdv[j]) || sdv[j] <= 0) { sdv[j] <- if (is.finite(w) && w > 0) w else 0.5; C[j, ] <- C[, j] <- 0; C[j, j] <- 1 }
    if (nm[j] %in% kf$at_bound && grepl("^log_sd", nm[j])) {
      ctr[j] <- (if (!is.null(pr)) pr[["upper"]] else th[[j]]) - 1.5
      sdv[j] <- max(sdv[j], 1.5); C[j, ] <- C[, j] <- 0; C[j, j] <- 1
    }
  }
  sdv <- 1.5 * sdv
  Sig <- C * outer(sdv, sdv)
  get <- function(TH, k) if (k %in% nm) TH[, k] else rep(NA_real_, nrow(TH))
  one_round <- function(ctr, Sig, np) {
    R <- tryCatch(chol(Sig), error = function(e) chol(diag(diag(Sig)) + 1e-10))
    Z <- matrix(stats::rnorm(np * p), np, p)
    u <- sqrt(stats::rchisq(np, df) / df)
    TH <- sweep((Z %*% R) / u, 2, ctr, "+"); colnames(TH) <- nm
    inb <- rowSums(sweep(TH, 2, kf$lo, ">") & sweep(TH, 2, kf$hi, "<")) == p
    lq <- -(df + p) / 2 * log1p(rowSums((Z / u)^2) / df)   # t density up to a constant
    if (kf$form == "exponential") { cc <- rep(1, np); dd <- get(TH, "mu") }
    else { cc <- 1 + get(TH, "b"); dd <- -get(TH, "b") * get(TH, "logK") }
    ll <- kf_uni_vec(kf$y, cc, dd, exp(get(TH, "log_sd_process")), exp(get(TH, "log_sd_observation")), get(TH, "x0"))
    lt <- ll + get(TH, "log_sd_process") + get(TH, "log_sd_observation")
    lw <- ifelse(inb & is.finite(lt), lt - lq, -Inf)
    w <- exp(lw - max(lw)); w <- w / sum(w)
    list(TH = TH, w = w, ess = 1 / sum(w^2))
  }
  # adaptive importance sampling: three pilot rounds re-centre and re-scale the
  # proposal on the weighted posterior mean and covariance (inflated 1.3x SD)
  for (it in 1:3) {
    r <- one_round(ctr, Sig, 4000)
    if (r$ess < 10) break
    ctr <- colSums(r$TH * r$w)
    D <- sweep(r$TH, 2, ctr); Sig <- 1.69 * crossprod(D * sqrt(r$w)) / (1 - sum(r$w^2))
    Sig <- Sig + diag(1e-8, p)
  }
  r <- one_round(ctr, Sig, n_prop)
  ess <- r$ess
  idx <- sample.int(n_prop, n_draws, replace = TRUE, prob = r$w)
  TH <- r$TH
  list(draws = TH[idx, , drop = FALSE], ess = ess)
}

# Parametric-bootstrap likelihood-ratio test of H0: no density dependence
# (EGSS, random walk with drift) against H1: Gompertz (GSS), the test Dennis
# & Taper (1994) and Dennis et al. (2006) recommend. Why not the interval for
# b: under H0 the ML estimate of b is biased towards negative values (a
# unit-root problem made worse by observation error), so "the 95% interval
# for b excludes 0" is a false positive in ~47% of series of 30 years
# (VALIDATION.md). The bootstrap reference distribution of the LR statistic is
# generated under the fitted H0 model with the same years missing, so it
# carries that bias with it.
gss_dd_test <- function(y, obs_error = TRUE, B = 199, seed = 1) {
  ok <- which(!is.na(y)); y <- y[min(ok):max(ok)]
  lr_of <- function(yy) {
    fE <- tryCatch(kf_fit_uni(yy, "exponential", obs_error, profile = FALSE), error = function(e) NULL)
    fG <- tryCatch(kf_fit_uni(yy, "gompertz", obs_error, profile = FALSE), error = function(e) NULL)
    if (is.null(fE) || is.null(fG)) return(list(lr = NA_real_, fE = fE))
    list(lr = max(0, 2 * (fG$loglik - fE$loglik)), fE = fE)
  }
  o <- lr_of(y)
  if (is.null(o$fE) || !is.finite(o$lr)) return(list(lr = NA_real_, p = NA_real_, B = 0L))
  th <- o$fE$theta; n <- length(y); miss <- is.na(y)
  sp <- exp(th[["log_sd_process"]]); so <- if (obs_error) exp(th[["log_sd_observation"]]) else 0
  x0 <- if (obs_error) th[["x0"]] else y[1] - th[["mu"]]
  set.seed(seed)
  lrb <- vapply(seq_len(B), function(i) {
    yb <- if (obs_error) x0 + cumsum(th[["mu"]] + stats::rnorm(n, 0, sp)) + stats::rnorm(n, 0, so)
          else y[1] + c(0, cumsum(th[["mu"]] + stats::rnorm(n - 1, 0, sp)))
    yb[miss] <- NA
    lr_of(yb)$lr
  }, numeric(1))
  lrb <- lrb[is.finite(lrb)]
  list(lr = o$lr, p = (1 + sum(lrb >= o$lr - 1e-8)) / (length(lrb) + 1), B = length(lrb),
       crit95 = if (length(lrb)) stats::quantile(lrb, 0.95, names = FALSE) else NA_real_)
}

# ---- goodness of fit: standardised one-step-ahead innovations ---------------
# If the model is right, v_t / sqrt(F_t) (count minus its prediction from all
# earlier counts, divided by its predicted SD) are independent N(0, 1)
# (Harvey 1989 ch. 5; Auger-Methe et al. 2021). Autocorrelated innovations
# point to a missing process (trend change, cycles, density dependence);
# heavy tails to catastrophes or outlying counts.
# diffuse = TRUE: first state has a flat prior (REML fits), first innovation dropped.
kf_innov <- function(y, cc, dd, sp, so, x0 = NA, diffuse = FALSE) {
  n <- length(y); e <- rep(NA_real_, n); a <- P <- NA
  for (t in seq_len(n)) {
    if (t == 1) {
      if (diffuse) { a <- 0; P <- 1e7 } else { a <- cc * x0 + dd; P <- sp^2 }
    } else { a <- cc * a + dd; P <- cc^2 * P + sp^2 }
    if (!is.na(y[t])) {
      Fv <- P + so^2; v <- y[t] - a
      if (!(diffuse && P > 1e6)) e[t] <- v / sqrt(Fv)
      K <- P / Fv; a <- a + K * v; P <- P * (1 - K)
    }
  }
  e
}

ts_resid_check <- function(e) {
  e <- e[is.finite(e)]; n <- length(e)
  lb <- if (n >= 8) stats::Box.test(e, lag = min(5, floor(n / 4)), type = "Ljung-Box")$p.value else NA_real_
  sw <- if (n >= 8 && stats::sd(e) > 0) stats::shapiro.test(e)$p.value else NA_real_
  list(n = n, ljung_box_p = lb, shapiro_p = sw, acf1 = if (n > 3) stats::acf(e, plot = FALSE, lag.max = 1)$acf[2] else NA_real_)
}
