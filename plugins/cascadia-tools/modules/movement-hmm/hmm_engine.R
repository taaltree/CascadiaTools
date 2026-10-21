# Engine for the movement-hmm module: track preparation, simulation, fitting
# (momentuHMM, multiple starts, state ordering) and parameter-uncertainty draws.
#
# Notation (matches momentuHMM / moveHMM):
#   row t of a track holds location x_t, the step from x_t to x_{t+1} and the
#   turning angle at x_t; state s_t governs both.
#   Transition s_{t-1} -> s_t uses the covariates of row t:
#     gamma_ij(t) = exp(eta_ij) / sum_k exp(eta_ik), eta_ii = 0,
#     eta_ij = beta0_ij + sum_l beta_lij * w_l(t)
#   beta columns are the off-diagonal cells filled row by row:
#   (1,2), (1,3), (2,1), (2,3), (3,1), (3,2).

`%||%` <- function(a, b) if (is.null(a)) b else a

offdiag_names <- function(K) {
  out <- character()
  for (i in seq_len(K)) for (j in seq_len(K)) if (i != j) out <- c(out, paste0(i, "to", j))
  out
}

# K x K transition matrix from one linear-predictor vector (length K(K-1)).
tpm_from_eta <- function(eta, K) {
  tG <- matrix(0, K, K)
  tG[!diag(K)] <- eta          # column-major fill of t(G) = row-wise fill of G
  G <- exp(t(tG))
  G / rowSums(G)
}

stationary_dist <- function(G) {
  K <- nrow(G)
  A <- t(diag(K) - G)
  A <- rbind(A, rep(1, K))
  as.vector(qr.solve(A, c(rep(0, K), 1)))
}

# ---- covariates -----------------------------------------------------------

# Time-of-day harmonics from the local clock hour (0-24).
tod_cols <- function(hour) cbind(tod_cos = cos(2 * pi * hour / 24), tod_sin = sin(2 * pi * hour / 24))

# Standardise transition covariates. Centre/scale are computed over the rows
# the HMM sees (regular grid, after filling covariates at missed fixes).
scale_covs <- function(df, covs) {
  sc <- list()
  for (v in covs) {
    x <- df[[v]]
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- 1
    sc[[v]] <- list(center = mean(x, na.rm = TRUE), scale = s)
    df[[paste0(v, "_z")]] <- (x - sc[[v]]$center) / sc[[v]]$scale
  }
  list(df = df, scaling = sc)
}

# Fill covariates at rows without a location (missed fixes): last value
# carried forward, then first value carried backward, within each animal.
fill_locf <- function(x) {
  if (all(is.na(x))) return(x)
  idx <- cumsum(!is.na(x)); idx[idx == 0] <- NA
  y <- x[!is.na(x)][idx]
  y[is.na(y)] <- x[!is.na(x)][1]
  y
}

# ---- regularisation ------------------------------------------------------

# Put each animal's fixes on a regular time grid. A grid time takes the
# nearest fix within `tol_frac` x interval; grid times without such a fix get
# NA coordinates (the HMM treats them as missing). Fixes that match no grid
# time are dropped and counted.
regularise_grid <- function(tr, interval_h, tol_frac = 0.25) {
  out <- list(); n_drop <- 0; n_missing <- 0
  dt <- interval_h * 3600
  for (id in unique(tr$ID)) {
    d <- tr[tr$ID == id, ]
    d <- d[order(d$time), ]
    t0 <- as.numeric(d$time[1]); t1 <- as.numeric(d$time[nrow(d)])
    grid <- seq(t0, t1 + dt * tol_frac, by = dt)
    k <- round((as.numeric(d$time) - t0) / dt) + 1
    off <- abs(as.numeric(d$time) - grid[pmin(k, length(grid))])
    keep <- off <= tol_frac * dt & k >= 1 & k <= length(grid)
    # if two fixes map to the same grid time keep the closer one
    o <- order(k, off); keep2 <- keep
    keep2[o] <- keep[o] & !duplicated(ifelse(keep[o], k[o], -seq_along(o)))
    n_drop <- n_drop + sum(!keep2)
    g <- data.frame(ID = id, time = as.POSIXct(grid, origin = "1970-01-01", tz = attr(d$time, "tzone") %||% "UTC"))
    oth <- setdiff(names(d), c("ID", "time"))
    for (v in oth) g[[v]] <- if (is.numeric(d[[v]])) NA_real_ else NA
    g[k[keep2], oth] <- d[keep2, oth]
    g$fix_ok <- FALSE; g$fix_ok[k[keep2]] <- TRUE
    n_missing <- n_missing + sum(!g$fix_ok)
    out[[length(out) + 1]] <- g
  }
  res <- do.call(rbind, out)
  attr(res, "n_dropped") <- n_drop
  attr(res, "n_missing") <- n_missing
  res
}

# Optional: continuous-time correlated random walk (crawl) to predict
# locations on the regular grid (single best prediction; see README).
regularise_crawl <- function(tr, interval_h, covs) {
  eco_require(c("momentuHMM", "crawl", "sf"))
  d <- tr[, c("ID", "time", "x", "y")]
  names(d)[1] <- "ID"
  cw <- momentuHMM::crawlWrap(d, timeStep = paste(interval_h * 60, "min"),
                              coord = c("x", "y"), Time.name = "time", attempts = 5, retryFits = 2)
  p <- as.data.frame(cw$crwPredict)
  p <- p[p$locType == "p", c("ID", "time", "mu.x", "mu.y")]
  names(p) <- c("ID", "time", "x", "y")
  # covariates: nearest original fix in time (covariates measured at fixes)
  for (v in covs) {
    p[[v]] <- NA_real_
    for (id in unique(p$ID)) {
      a <- tr[tr$ID == id, ]; b <- which(p$ID == id)
      j <- findInterval(as.numeric(p$time[b]), as.numeric(a$time), all.inside = TRUE)
      jn <- ifelse(abs(as.numeric(a$time[j]) - as.numeric(p$time[b])) <= abs(as.numeric(a$time[pmin(j + 1, nrow(a))]) - as.numeric(p$time[b])), j, pmin(j + 1, nrow(a)))
      p[[v]][b] <- a[[v]][jn]
    }
  }
  p$fix_ok <- TRUE
  attr(p, "n_dropped") <- 0; attr(p, "n_missing") <- 0
  p
}

# Nominal fix interval (hours) from the within-animal gaps. The median gap is
# inflated when many fixes are missed (gaps of 2, 3 ... intervals) and by
# fix-time jitter, so use the MODE of the gap distribution (the scheduled
# interval is the most common gap) and snap it to a standard collar schedule
# when within 10%; otherwise round to whole minutes. Found with an irregular
# test set (30% fixes lost, +/-25 min jitter): the median gave 77 min for a
# 60-min schedule, inflating step lengths.
guess_interval <- function(gaps_h) {
  g <- gaps_h[is.finite(gaps_h) & gaps_h > 0] * 60
  if (!length(g)) stop("cannot work out the fix interval: fewer than two fixes per animal")
  g <- g[g <= 5 * stats::quantile(g, 0.25)]          # drop very long outages before the density
  md <- if (length(unique(g)) > 1) { d <- stats::density(g); d$x[which.max(d$y)] } else g[1]
  nice <- c(1, 2, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 240, 300, 360, 480, 720, 1440)
  j <- which.min(abs(log(nice / md)))
  m <- if (abs(nice[j] / md - 1) <= 0.10) nice[j] else max(1, round(md))
  m / 60
}

# From the module data structure to a momentuHMMData object + bookkeeping.
prep_track <- function(data, settings) {
  eco_require("momentuHMM")
  tr <- data$track
  covs <- settings$covariates %||% data$covariates %||% character()
  use_tod <- isTRUE(settings$tod %||% data$tod %||% TRUE) && "time" %in% names(tr)
  tr <- tr[order(tr$ID, tr$time), ]
  # interval: median within-animal gap (hours) unless given
  gaps <- unlist(lapply(split(as.numeric(tr$time), tr$ID), function(t) diff(t))) / 3600
  # The median gap is inflated slightly by missed fixes and fix-time jitter
  # (e.g. 60.1 min); round it to whole minutes so the grid does not drift
  # away from the fix schedule over hundreds of steps.
  interval <- settings$interval_h %||% data$interval_h %||% guess_interval(gaps)
  reg <- settings$regularise %||% "grid"
  g <- if (identical(reg, "crawl")) regularise_crawl(tr, interval, covs)
       else regularise_grid(tr, interval, settings$tolerance %||% 0.25)
  for (v in covs) g[[v]] <- stats::ave(g[[v]], g$ID, FUN = fill_locf)
  hour <- as.numeric(format(g$time, "%H", tz = data$tz %||% "UTC")) +
          as.numeric(format(g$time, "%M", tz = data$tz %||% "UTC")) / 60
  g$hour <- hour
  if (use_tod) { tc <- tod_cols(hour); g$tod_cos <- tc[, 1]; g$tod_sin <- tc[, 2] }
  sc <- scale_covs(g, covs)
  g <- sc$df
  terms <- c(if (use_tod) c("tod_cos", "tod_sin"), if (length(covs)) paste0(covs, "_z"))
  df <- data.frame(ID = g$ID, x = g$x, y = g$y, g[, c(terms, "hour"), drop = FALSE])
  pd <- momentuHMM::prepData(df, type = "UTM", coordNames = c("x", "y"),
                             covNames = c(terms, "hour"))
  pd$time <- g$time
  list(pd = pd, grid = g, terms = terms, covs = covs, tod = use_tod, scaling = sc$scaling,
       interval_h = interval, phase0 = grid_phase(hour, interval), n_dropped = attr(g, "n_dropped"), n_missing = attr(g, "n_missing"),
       regularise = reg)
}

# ---- simulation ----------------------------------------------------------

# Smooth habitat surface (e.g. % forest cover) over a square landscape (m).
forest_field <- function(L, seed) {
  set.seed(seed + 7)
  ph <- stats::runif(4, 0, 2 * pi); wl <- stats::runif(4, 3000, 7000)
  function(x, y) {
    v <- 50 + 22 * sin(2 * pi * x / wl[1] + ph[1]) * cos(2 * pi * y / wl[2] + ph[2]) +
              18 * sin(2 * pi * (x + y) / wl[3] + ph[3]) * sin(2 * pi * (x - y) / wl[4] + ph[4])
    pmin(pmax(v, 0), 100)
  }
}

rvonmises <- function(n, mu, kappa) {
  # Best & Fisher (1979) rejection sampler
  out <- numeric(n)
  for (i in seq_len(n)) {
    if (kappa < 1e-8) { out[i] <- stats::runif(1, -pi, pi); next }
    tau <- 1 + sqrt(1 + 4 * kappa^2); rho <- (tau - sqrt(2 * tau)) / (2 * kappa); r <- (1 + rho^2) / (2 * rho)
    repeat {
      u <- stats::runif(3); z <- cos(pi * u[1]); f <- (1 + r * z) / (r + z); c <- kappa * (r - f)
      if (c * (2 - c) - u[2] > 0 || log(c / u[2]) + 1 - c >= 0) break
    }
    th <- sign(u[3] - 0.5) * acos(f) + mu
    out[i] <- (th + pi) %% (2 * pi) - pi
  }
  out
}

# eta (length K(K-1)) for raw covariate values: hour and forest (%).
sim_eta <- function(params, hour, forest) {
  tc <- tod_cols(hour)
  params$tpm_intercept + params$tpm_tod_cos * tc[1] + params$tpm_tod_sin * tc[2] +
    params$tpm_forest * (forest - (params$forest_ref %||% 50))
}

simulate_tracks <- function(params, design, seed) {
  set.seed(seed)
  K <- length(params$step_mean)
  stopifnot(length(params$tpm_intercept) == K * (K - 1))
  L <- (design$landscape_km %||% 20) * 1000
  fld <- forest_field(L, seed)
  dt <- design$fix_interval_h %||% 1
  zm <- params$zero_mass %||% rep(0, K)
  shp <- params$step_mean^2 / params$step_sd^2; rte <- params$step_mean / params$step_sd^2
  t_start <- as.POSIXct(design$start %||% "2024-05-01 00:00:00", tz = "UTC")
  rows <- list()
  for (a in seq_len(design$n_animals)) {
    n <- design$n_fixes
    x <- y <- forest <- numeric(n); s <- integer(n); step <- angle <- numeric(n)
    x[1] <- stats::runif(1, 0.25 * L, 0.75 * L); y[1] <- stats::runif(1, 0.25 * L, 0.75 * L)
    hour0 <- (as.numeric(format(t_start, "%H")) + stats::runif(1, 0, 24)) %/% dt * dt
    hours <- (hour0 + (seq_len(n) - 1) * dt) %% 24
    head <- stats::runif(1, -pi, pi)
    for (t in seq_len(n)) {
      forest[t] <- fld(x[t], y[t])
      if (t == 1) {
        s[t] <- sample.int(K, 1, prob = stationary_dist(tpm_from_eta(sim_eta(params, hours[t], forest[t]), K)))
      } else {
        G <- tpm_from_eta(sim_eta(params, hours[t], forest[t]), K)
        s[t] <- sample.int(K, 1, prob = G[s[t - 1], ])
      }
      if (t == n) break
      step[t] <- if (stats::runif(1) < zm[s[t]]) 0 else stats::rgamma(1, shp[s[t]], rte[s[t]])
      angle[t] <- rvonmises(1, 0, params$angle_kappa[s[t]])
      head <- head + angle[t]
      x[t + 1] <- x[t] + step[t] * cos(head); y[t + 1] <- y[t] + step[t] * sin(head)
    }
    tm <- t_start + (hour0 + (seq_len(n) - 1) * dt) * 3600
    d <- data.frame(ID = sprintf("A%02d", a), time = tm, x = x, y = y, forest = forest, state_true = s)
    rows[[a]] <- d
  }
  tr <- do.call(rbind, rows)
  # GPS realism: failed fixes are absent rows; fix times jitter by a few minutes
  miss <- stats::runif(nrow(tr)) < (design$missing_prop %||% 0)
  first <- !duplicated(tr$ID)
  miss[first] <- FALSE
  jit <- design$jitter_min %||% 0
  tr$time <- tr$time + round(stats::runif(nrow(tr), -jit, jit) * 60)
  list(track = tr[!miss, ], true_states = tr$state_true, field = fld)
}

# True values of the estimands, in exactly the parameterisation fit() reports:
# slopes per original covariate unit; intercepts at the covariate centre
# (sample mean, as fit() standardises) and with tod harmonics at 0.
truth_from_params <- function(params, center_forest, dt = 1) {
  if (!periodic_ok(dt)) stop("simulation needs a fix interval that divides 24 h")
  K <- length(params$step_mean)
  nm <- offdiag_names(K)
  fr <- params$forest_ref %||% 50
  tr <- c(stats::setNames(params$step_mean, paste0("step_mean_m_state", 1:K)),
          stats::setNames(params$step_sd, paste0("step_sd_m_state", 1:K)),
          stats::setNames(params$angle_kappa, paste0("turn_concentration_state", 1:K)),
          stats::setNames(params$tpm_intercept + params$tpm_forest * (center_forest - fr), paste0("tpm_", nm, "_intercept")),
          stats::setNames(params$tpm_tod_cos, paste0("tpm_", nm, "_tod_cos")),
          stats::setNames(params$tpm_tod_sin, paste0("tpm_", nm, "_tod_sin")),
          stats::setNames(params$tpm_forest, paste0("tpm_", nm, "_forest")))
  # long-run share of time: periodically stationary distribution of the
  # simulated chain (fix times are at multiples of dt from midnight), averaged
  # over the steps of one day; written directly from sim_eta, not via the
  # fit-side activity_budget().
  hrs <- seq(0, 24 - dt, by = dt)
  Gs <- lapply(hrs, function(h) tpm_from_eta(sim_eta(params, h, center_forest), K))
  P <- diag(K); for (G in Gs) P <- P %*% G
  ev <- eigen(t(P)); d <- Re(ev$vectors[, which.min(abs(ev$values - 1))]); d <- d / sum(d)
  acc <- 0; for (G in Gs) { d <- as.vector(d %*% G); acc <- acc + d }
  ab <- acc / length(Gs)
  c(tr, stats::setNames(ab, paste0("activity_budget_state", 1:K)))
}

# ---- fitting -------------------------------------------------------------

# Data-driven starting values: state means spread over step-length
# quantiles (state 1 = slowest), sd = mean, concentration rising with speed.
start_values <- function(step, K, jitter = FALSE, zm = FALSE) {
  s <- step[is.finite(step) & step > 0]
  q <- switch(as.character(K), "2" = c(0.25, 0.8), "3" = c(0.15, 0.55, 0.9),
              seq(0.1, 0.9, length.out = K))
  mu <- as.numeric(stats::quantile(s, q))
  kap <- seq(0.3, 3, length.out = K)
  if (jitter) {
    mu <- sort(mu * exp(stats::rnorm(K, 0, 0.5)))
    kap <- sort(exp(stats::runif(K, log(0.1), log(6))))
  }
  sdv <- mu * (if (jitter) exp(stats::rnorm(K, 0, 0.3)) else 1)
  stepPar <- c(mu, sdv, if (zm) rep(0.05, K))
  list(step = stepPar, angle = kap)
}

hmm_fit_once <- function(prep, K, Par0, beta0 = NULL, delta0 = NULL, retryFits = 0) {
  form <- if (length(prep$terms)) stats::as.formula(paste("~", paste(prep$terms, collapse = " + "))) else ~1
  suppressMessages(suppressWarnings(
    momentuHMM::fitHMM(prep$pd, nbStates = K, dist = list(step = "gamma", angle = "vm"),
                       Par0 = Par0, beta0 = beta0, delta0 = delta0, formula = form,
                       estAngleMean = list(angle = FALSE), retryFits = retryFits,
                       nlmPar = list(iterlim = 1000))))
}

# Relabel states so that state 1 has the shortest mean step. Refits from the
# permuted MLE, which reproduces the same optimum with ordered labels.
order_states <- function(m, prep, K) {
  mu <- m$mle$step["mean", ]
  p <- order(mu)
  if (identical(p, seq_len(K))) return(m)
  st <- m$mle$step; an <- m$mle$angle
  Par0 <- list(step = as.vector(t(st[, p, drop = FALSE])), angle = an["concentration", p])
  B <- m$mle$beta; nm <- offdiag_names(K)
  B2 <- B
  for (a in seq_len(K)) for (b in seq_len(K)) if (a != b)
    B2[, match(paste0(a, "to", b), nm)] <- B[, match(paste0(p[a], "to", p[b]), nm)]
  d0 <- as.vector(m$mle$delta[1, p])
  m2 <- hmm_fit_once(prep, K, Par0, beta0 = B2, delta0 = d0)
  attr(m2, "relabel_llk_diff") <- m2$mod$minimum - m$mod$minimum
  m2
}

# Parameter draws from the asymptotic normal distribution of the working
# parameters (MLE, inverse Hessian). Natural-scale quantities are transformed
# draw by draw, so intervals respect bounds (log scale for means/sd/kappa).
hmm_draws <- function(m, K, ndraws, zm) {
  est <- m$mod$estimate
  H <- m$mod$hessian
  nst <- 2 * K + (if (zm) K else 0)
  nb <- length(m$mle$beta)
  idx <- seq_len(nst + K + nb)
  # Positive definite with a RELATIVE tolerance and a working Cholesky factor.
  # (Independent verification 2026-10-19: with delta on a boundary the
  # numerical Hessian has H[delta, delta] = 0 and tiny off-diagonals, solve(H)
  # "succeeds", but the step/angle/beta block of that inverse is rank-deficient
  # (smallest eigenvalue ~1e-19 > 0), chol() then failed and ~8% of fits were
  # wrongly reported as having no standard errors instead of using the
  # conditional fallback below.)
  # Zero-step probabilities estimated at ~0 (e.g. a single zero step in the
  # data turns the zero-mass branch on for every state) sit on a boundary
  # (logit -> -Inf) and make the Hessian singular; condition on them too.
  zm_bd <- if (zm) (2 * K + seq_len(K))[est[2 * K + seq_len(K)] < stats::qlogis(1e-6)] else integer()
  free <- setdiff(idx, zm_bd)
  pd_ok <- function(V) {
    if (is.null(V) || !all(is.finite(V[free, free])) || !all(diag(V)[free] > 0)) return(FALSE)
    S <- (V[free, free] + t(V[free, free])) / 2
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    min(ev) > 1e-10 * max(ev) && !is.null(tryCatch(chol(S), error = function(e) NULL))
  }
  # delta on a boundary (an animal-start state with probability ~0): its
  # Hessian row is ~0 and contaminates the full inverse, so condition on it.
  dl <- m$mle$delta
  delta_bd <- !is.null(dl) && any(dl < 1e-6)
  V <- if (delta_bd || length(zm_bd)) NULL else tryCatch(solve(H), error = function(e) NULL)
  ok <- pd_ok(V)
  conditional <- FALSE
  if (!ok && ncol(H) > length(free)) {
    # The initial-state distribution (delta) is often on a boundary (all
    # animals start in the same state), which makes the full Hessian
    # singular. It is a nuisance parameter: condition on its estimate (and on
    # boundary zero-step probabilities) and invert the block for the
    # remaining step, angle and transition parameters.
    V <- matrix(NA_real_, ncol(H), ncol(H))
    Vs <- tryCatch(solve(H[free, free]), error = function(e) NULL)
    if (!is.null(Vs)) { V[free, free] <- Vs; ok <- pd_ok(V); conditional <- ok }
  }
  W <- if (ok) {
    L <- tryCatch(chol((V[free, free] + t(V[free, free])) / 2), error = function(e) NULL)
    if (is.null(L)) { ok <- FALSE; NULL } else {
      Wf <- sweep(matrix(stats::rnorm(ndraws * length(free)), ndraws) %*% L, 2, est[free], "+")
      W <- matrix(est[idx], ndraws, length(idx), byrow = TRUE); W[, match(free, idx)] <- Wf; W
    }
  } else NULL
  list(ok = ok, est = est[idx], W = W, nst = nst, nb = nb, conditional = conditional && (delta_bd || !length(zm_bd)),
       zm_boundary = zm_bd - 2 * K)
}

# Working vector -> natural parameters (one draw).
w2nat <- function(w, K, zm, beta_dim) {
  mu <- exp(w[1:K]); sdv <- exp(w[K + 1:K])
  z <- if (zm) stats::plogis(w[2 * K + 1:K]) else rep(0, K)
  o <- 2 * K + (if (zm) K else 0)
  kap <- exp(w[o + 1:K])
  B <- matrix(w[o + K + seq_len(prod(beta_dim))], beta_dim[1], beta_dim[2])
  list(mu = mu, sd = sdv, zm = z, kappa = kap, beta = B)
}

# With time of day on the transitions the state process is periodically
# inhomogeneous (period 24 h). Its long-run behaviour is described by the
# periodically stationary distribution (Koslik et al. 2025, Ann. Appl. Stat.):
# delta_h = probability of each state at a time step with clock hour h,
# satisfying delta_{h} = delta_{h-dt} Gamma(h). The hour-specific stationary
# distributions delta(Gamma(h)) are NOT this: they describe a hypothetical
# chain frozen at hour h (independent verification, 2026-10-19: averaging
# them misstated the 3-state activity budget by ~2 percentage points).
periodic_ok <- function(dt) is.numeric(dt) && length(dt) == 1 && dt > 0 && dt <= 24 &&
  abs(24 / dt - round(24 / dt)) < 1e-6
periodic_probs <- function(B, terms, K, dt = 1, phase0 = 0, z = list()) {
  nT <- round(24 / dt)
  hrs <- (phase0 + (seq_len(nT) - 1) * dt) %% 24
  Gs <- lapply(hrs, function(h) tpm_from_eta(as.vector(design_row(terms, hour = h, z = z) %*% B), K))
  d <- stationary_dist(Reduce(`%*%`, Gs))   # state distribution at a step with hour hrs[nT]
  P <- matrix(NA_real_, nT, K)
  for (i in seq_len(nT)) { d <- as.vector(d %*% Gs[[i]]); P[i, ] <- d }
  list(hours = hrs, P = P)
}

# Activity budget: long-run share of time in each state at the mean (or given
# z) of the habitat covariates. With time of day: the periodically stationary
# distribution averaged over the time steps of one day (exact when 24 h is a
# whole number of fix intervals; otherwise the hour-specific stationary
# distributions averaged over the day, an approximation).
activity_budget <- function(B, terms, K, dt = 1, phase0 = 0, z = list()) {
  if (!"tod_cos" %in% terms)
    return(stationary_dist(tpm_from_eta(as.vector(design_row(terms, z = z) %*% B), K)))
  if (periodic_ok(dt)) return(colMeans(periodic_probs(B, terms, K, dt, phase0, z)$P))
  hh <- seq(0, 24, length.out = 97)[-97]
  rowMeans(matrix(sapply(hh, function(h)
    stationary_dist(tpm_from_eta(as.vector(design_row(terms, hour = h, z = z) %*% B), K))), nrow = K))
}

# Clock-hour phase of the regular grid within one fix interval (quarter-hour
# resolution, most common value over all animals), e.g. 1 for 2-h fixes at
# 01:00, 03:00, ...
grid_phase <- function(hour, dt) {
  if (!periodic_ok(dt)) return(0)
  ph <- (round(hour * 4) / 4) %% dt
  ph <- ph[is.finite(ph)]
  if (!length(ph)) return(0)
  as.numeric(names(which.max(table(ph))))
}

# Linear-predictor row for given hour / standardised covariates.
design_row <- function(terms, hour = NULL, z = list()) {
  xr <- c(1, stats::setNames(rep(0, length(terms)), terms))
  if (!is.null(hour) && "tod_cos" %in% terms) { xr["tod_cos"] <- cos(2 * pi * hour / 24); xr["tod_sin"] <- sin(2 * pi * hour / 24) }
  for (v in names(z)) xr[paste0(v, "_z")] <- z[[v]]
  xr
}
