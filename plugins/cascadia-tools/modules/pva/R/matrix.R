# Matrix-model PVA (settings$model = "matrix").
#
# A stage/age-structured projection matrix whose entries are expressions in
# named vital rates (e.g. "fec * s_juv"). Each year every vital rate is drawn
# from its environmental distribution:
#   survival / probabilities: beta with the given mean and temporal SD
#   fecundity: lognormal (default) or stretched beta on [min, max]
# with optional within-year correlation among rates via a Gaussian copula
# (correlated normal deviates pushed through each rate's inverse CDF; Morris
# & Doak 2002 ch. 8). Demographic stochasticity (binomial survival and stage
# transitions, Poisson births), an optional population ceiling, and
# parameter uncertainty (means drawn from their sampling distribution; the
# environmental variance from a scaled inverse chi-square when the number of years
# behind it is known) are layered on top. Quasi-extinction = the counted
# stages fall to or below the threshold in any year up to the horizon.
#
# Engine: hand-written vectorised projection (all trajectories at once).
# Cross-checked against popbio::vitalsim / stoch.growth.rate (VALIDATION.md).

PVA_ZGRID <- seq(-6, 6, length.out = 1201)

# --- distributions -----------------------------------------------------------

beta_ab <- function(m, v) {
  v <- pmin(v, 0.95 * m * (1 - m))
  k <- m * (1 - m) / v - 1
  list(a = m * k, b = (1 - m) * k)
}
lnorm_ms <- function(m, s) { sl2 <- log1p((s / m)^2); list(mu = log(m) - sl2 / 2, sd = sqrt(sl2)) }

# Draw parameter means from their sampling distributions (uncertainty about
# the estimates). Returns ndraw x nrates matrices of means and env SDs. u: a
# list of stored uniforms so scenarios can reuse the same draws.
mx_param_draws <- function(rates, ndraw, u = NULL, corr = NULL) {
  nr <- nrow(rates)
  if (is.null(u)) {
    u <- list(m = matrix(stats::runif(ndraw * nr), ndraw, nr),
              v = matrix(stats::rchisq(ndraw * nr, df = pmax(rep(rates$n_years, each = ndraw) - 1, 1)), ndraw, nr))
    # Means estimated from the SAME years (raw monitoring, n_years and sampling_var
    # known) share the year effects, so with correlated environments their errors
    # are correlated: corr(mean_i, mean_j) = rho_ij sd_i sd_j / sqrt(V_i V_j), V =
    # sd_env^2 + sampling variance. Drawing them independently understated the
    # uncertainty of lambda (coverage 0.895 at rho = 0.4; independent verification,
    # VALIDATION.md). Rates without sampling_var (vital-rate table) stay independent.
    if (!is.null(corr) && is.matrix(corr) && ndraw > 1) {
      sd <- rates$sd_env; sp <- rates$sampling_var
      ok <- is.finite(rates$n_years) & rates$n_years >= 2 & is.finite(sp) & is.finite(sd) & rates$se > 0
      if (sum(ok) >= 2) {
        V <- sqrt(sd^2 + sp)
        R <- diag(nr)
        for (i in which(ok)) for (j in which(ok)) if (i != j) R[i, j] <- corr[i, j] * sd[i] * sd[j] / (V[i] * V[j])
        L <- tryCatch(chol(R), error = function(e) NULL)
        if (!is.null(L) && any(R[upper.tri(R)] != 0)) u$m <- stats::pnorm(stats::qnorm(u$m) %*% L)
      }
    }
  }
  M <- S <- matrix(NA_real_, ndraw, nr, dimnames = list(NULL, rates$rate))
  for (r in seq_len(nr)) {
    m <- rates$mean[r]; se <- rates$se[r]
    if (ndraw == 1 || !is.finite(se) || se <= 0) M[, r] <- m
    else if (rates$type[r] == "survival") {
      ab <- beta_ab(m, se^2); M[, r] <- stats::qbeta(u$m[, r], ab$a, ab$b)
    } else if (rates$dist[r] == "stretched_beta") {
      lo <- rates$min[r]; hi <- rates$max[r]; ab <- beta_ab((m - lo) / (hi - lo), (se / (hi - lo))^2)
      M[, r] <- lo + (hi - lo) * stats::qbeta(u$m[, r], ab$a, ab$b)
    } else {
      ln <- lnorm_ms(m, se); M[, r] <- stats::qlnorm(u$m[, r], ln$mu, ln$sd)
    }
    sv <- rates$sd_env[r]^2; sp <- rates$sampling_var[r]; ny <- rates$n_years[r]
    if (ndraw == 1 || !is.finite(ny) || ny < 2) S[, r] <- sqrt(sv)
    else {
      sp <- if (is.finite(sp)) sp else 0
      # chi-square pivot: (ny-1) s^2 / sigma^2 ~ chi2(ny-1), so draws of the
      # total variance are s^2 (ny-1) / chi2 (scaled inverse chi-square; the
      # posterior under the Jeffreys prior), not s^2 chi2 / (ny-1), which is
      # skewed the wrong way and gave too-low upper limits.
      tot <- (sv + sp) * (ny - 1) / u$v[, r]
      S[, r] <- sqrt(pmax(tot - sp, 0))
    }
  }
  list(M = M, S = S, u = u)
}

# Inverse-CDF lookup tables (one row per parameter draw) for beta-type rates,
# so correlated normal deviates can be transformed quickly.
mx_tables <- function(rates, M, S) {
  tabs <- vector("list", nrow(rates))
  pz <- stats::pnorm(PVA_ZGRID)
  for (r in seq_len(nrow(rates))) {
    if (rates$type[r] == "survival" || rates$dist[r] == "stretched_beta") {
      lo <- if (rates$type[r] == "survival") 0 else rates$min[r]
      hi <- if (rates$type[r] == "survival") 1 else rates$max[r]
      m <- (M[, r] - lo) / (hi - lo); v <- (S[, r] / (hi - lo))^2
      tab <- matrix(NA_real_, nrow(M), length(PVA_ZGRID))
      for (b in seq_len(nrow(M))) {
        if (v[b] <= 1e-12) tab[b, ] <- m[b]
        else { ab <- beta_ab(m[b], v[b]); tab[b, ] <- stats::qbeta(pz, ab$a, ab$b) }
      }
      tabs[[r]] <- lo + (hi - lo) * tab
    }
  }
  tabs
}

mx_interp <- function(tab, didx, z) {
  g <- PVA_ZGRID; dz <- g[2] - g[1]
  pos <- (pmin(pmax(z, g[1]), g[length(g)]) - g[1]) / dz + 1
  i <- pmin(floor(pos), length(g) - 1); w <- pos - i
  tab[cbind(didx, i)] * (1 - w) + tab[cbind(didx, i + 1)] * w
}

# Annual vital rates for every trajectory: list of vectors.
mx_vr_year <- function(rates, M, S, tabs, didx, L) {
  nr <- nrow(rates); n <- length(didx)
  z <- matrix(stats::rnorm(n * nr), n, nr)
  if (!is.null(L)) z <- z %*% L
  out <- vector("list", nr); names(out) <- rates$rate
  for (r in seq_len(nr)) {
    out[[r]] <- if (!is.null(tabs[[r]])) mx_interp(tabs[[r]], didx, z[, r])
    else {
      m <- M[didx, r]; s <- S[didx, r]
      ln <- lnorm_ms(m, pmax(s, 1e-12)); exp(ln$mu + ln$sd * z[, r])
    }
  }
  out
}

# Parse the matrix definition: k x k character -> list of entries with
# expression, row, col and whether it is reproduction (uses a fecundity rate).
mx_parse <- function(A, rates) {
  k <- nrow(A); ent <- list()
  fec <- rates$rate[rates$type == "fecundity"]
  for (i in seq_len(k)) for (j in seq_len(k)) {
    s <- trimws(as.character(A[i, j]))
    if (s %in% c("", "0", "0.0")) next
    e <- parse(text = s)[[1]]
    v <- all.vars(e)
    ent[[length(ent) + 1]] <- list(i = i, j = j, expr = e, vars = v, txt = s,
                                   repro = any(v %in% fec))
  }
  ent
}

mx_eval <- function(ent, vr, n) {
  lapply(ent, function(e) { x <- eval(e$expr, vr); if (length(x) == 1) rep(x, n) else x })
}

mx_mean_matrix <- function(ent, vals, k) {
  A <- matrix(0, k, k)
  for (e in ent) A[e$i, e$j] <- eval(e$expr, as.list(vals))
  A
}

mx_lambda1 <- function(A) max(Re(eigen(A, only.values = TRUE)$values))

mx_sens <- function(A) {
  ev <- eigen(A); i <- which.max(Re(ev$values))
  w <- abs(Re(ev$vectors[, i])); w <- w / sum(w)
  evl <- eigen(t(A)); il <- which.max(Re(evl$values))
  v <- abs(Re(evl$vectors[, il]))
  S <- outer(v, w) / sum(v * w)
  list(lambda = Re(ev$values[i]), S = S, w = w, v = v)
}

# Vital-rate elasticities of lambda1 by central finite differences.
mx_vr_elast <- function(ent, vals, k) {
  l0 <- mx_lambda1(mx_mean_matrix(ent, vals, k))
  vapply(names(vals), function(r) {
    h <- 1e-5 * max(abs(vals[[r]]), 1e-3)
    up <- vals; dn <- vals; up[[r]] <- up[[r]] + h; dn[[r]] <- dn[[r]] - h
    (mx_lambda1(mx_mean_matrix(ent, up, k)) - mx_lambda1(mx_mean_matrix(ent, dn, k))) / (2 * h) * vals[[r]] / l0
  }, numeric(1))
}

# Core projection. ndraw parameter sets x ntraj trajectories each.
# Returns per-draw quasi-extinction CDF (ndraw x tmax), per-draw log lambda_s
# (environment-only companion projection), Tuljapurkar approximation per draw,
# and (optionally) quantiles of counted N over time.
mx_project <- function(mx, M, S, ntraj, tmax, seed, keep_paths = FALSE, burn = 10, tul_years = 20) {
  rates <- mx$rates; k <- length(mx$stages); ent <- mx$ent
  ndraw <- nrow(M); ntot <- ndraw * ntraj
  didx <- rep(seq_len(ndraw), each = ntraj)
  tabs <- mx_tables(rates, M, S)
  L <- mx$cholcor
  cnt <- mx$count_idx
  N <- matrix(rep(mx$n0, each = ntot), ntot, k)
  demog <- isTRUE(mx$demographic)
  if (demog) N <- round(N)
  # stable-stage start for the environment-only lambda_s companion
  A0 <- mx_mean_matrix(ent, rates$mean |> stats::setNames(rates$rate), k)
  w0 <- mx_sens(A0)$w
  W <- matrix(rep(w0, each = ntot), ntot, k)
  dead <- rowSums(N[, cnt, drop = FALSE]) <= mx$threshold
  cdf <- matrix(0, ndraw, tmax)
  lgr <- numeric(ndraw); nlg <- 0
  qN <- if (keep_paths) matrix(NA_real_, tmax + 1, 5) else NULL
  if (keep_paths) qN[1, ] <- stats::quantile(rowSums(N[, cnt, drop = FALSE]), c(.025, .25, .5, .75, .975))
  samp <- if (keep_paths) matrix(NA_real_, tmax + 1, min(ntot, 60)) else NULL
  if (keep_paths) samp[1, ] <- rowSums(N[seq_len(ncol(samp)), cnt, drop = FALSE])
  # Tuljapurkar accumulators: sums and cross-products of nonzero entries
  ne <- length(ent)
  s1 <- matrix(0, ndraw, ne); s2 <- array(0, c(ndraw, ne, ne)); nt <- 0
  for (t in seq_len(tmax)) {
    set.seed(seed * 1000 + t)   # common random numbers across scenarios
    vr <- mx_vr_year(rates, M, S, tabs, didx, L)
    a <- mx_eval(ent, vr, ntot)
    # --- environment-only projection for lambda_s
    Wn <- matrix(0, ntot, k)
    for (e in seq_along(ent)) Wn[, ent[[e]]$i] <- Wn[, ent[[e]]$i] + a[[e]] * W[, ent[[e]]$j]
    g <- rowSums(Wn)
    if (t > burn) { lgr <- lgr + as.vector(rowsum(log(g), didx)) / ntraj; nlg <- nlg + 1 }
    W <- Wn / g
    if (t <= tul_years) {
      am <- do.call(cbind, a)
      s1 <- s1 + rowsum(am, didx) / ntraj
      for (e1 in seq_len(ne)) s2[, e1, ] <- s2[, e1, ] + rowsum(am * am[, e1], didx) / ntraj
      nt <- nt + 1
    }
    # --- population projection
    Nn <- matrix(0, ntot, k)
    if (demog) {
      for (j in seq_len(k)) {
        nj <- N[, j]
        if (!any(nj > 0)) next
        rem <- nj; prem <- rep(1, ntot)
        for (e in seq_along(ent)) {
          en <- ent[[e]]
          if (en$j != j) next
          if (en$repro) {
            Nn[, en$i] <- Nn[, en$i] + stats::rpois(ntot, a[[e]] * nj)
          } else {
            p <- pmin(pmax(a[[e]] / pmax(prem, 1e-12), 0), 1)
            x <- stats::rbinom(ntot, rem, p)
            Nn[, en$i] <- Nn[, en$i] + x
            rem <- rem - x; prem <- prem - a[[e]]
          }
        }
      }
    } else {
      for (e in seq_along(ent)) Nn[, ent[[e]]$i] <- Nn[, ent[[e]]$i] + a[[e]] * N[, ent[[e]]$j]
    }
    if (is.finite(mx$ceiling %||% Inf)) {
      tot <- rowSums(Nn[, cnt, drop = FALSE])
      over <- tot > mx$ceiling
      if (any(over)) {
        f <- mx$ceiling / tot[over]
        Nn[over, ] <- if (demog) matrix(stats::rbinom(sum(over) * k, Nn[over, ], f), sum(over), k) else Nn[over, ] * f
      }
    }
    N <- Nn
    nc <- rowSums(N[, cnt, drop = FALSE])
    dead <- dead | nc <= mx$threshold
    cdf[, t] <- as.vector(rowsum(as.numeric(dead), didx)) / ntraj
    if (keep_paths) {
      qN[t + 1, ] <- stats::quantile(nc, c(.025, .25, .5, .75, .975))
      samp[t + 1, ] <- nc[seq_len(ncol(samp))]
    }
  }
  # Tuljapurkar: log lambda_s ~ log lambda1(Abar) - tau2 / (2 lambda1^2)
  tul <- if (nt == 0) rep(NA_real_, ndraw) else vapply(seq_len(ndraw), function(b) {
    m1 <- s1[b, ] / nt
    C <- s2[b, , ] / nt - outer(m1, m1)
    Ab <- matrix(0, k, k)
    for (e in seq_along(ent)) Ab[ent[[e]]$i, ent[[e]]$j] <- m1[e]
    sn <- mx_sens(Ab)
    sv <- vapply(ent, function(en) sn$S[en$i, en$j], numeric(1))
    log(sn$lambda) - as.numeric(t(sv) %*% C %*% sv) / (2 * sn$lambda^2)
  }, numeric(1))
  list(cdf = cdf, loglams = lgr / max(nlg, 1), tul = tul, qN = qN, samp = samp)
}

# --- estimating rates from raw annual monitoring data --------------------------

# demog: data.frame(year, rate, n, events). For survival-type rates events =
# survivors of n at risk; for fecundity events = offspring produced by n
# breeders. Returns mean of annual estimates, its SE, the environmental SD
# with mean sampling variance removed (moment estimator; cf. Kendall 1998,
# Morris & Doak 2002 ch. 8), number of years and mean sampling variance.
mx_rates_from_demography <- function(demog, types) {
  out <- lapply(split(demog, demog$rate), function(d) {
    r <- d$rate[1]; ty <- types[[r]] %||% "survival"
    est <- d$events / d$n
    sv <- if (ty == "survival") est * (1 - est) / d$n else est / d$n
    ny <- length(est)
    vt <- if (ny > 1) stats::var(est) else NA
    # SE of a mean of ny annual estimates; inflated by t(ny-1)/z so that the
    # parameter draws (beta/lognormal with this SE) give t-type intervals in
    # short series (Student small-sample correction; checked by recovery).
    tf <- if (ny > 1) stats::qt(0.975, ny - 1) / stats::qnorm(0.975) else 1
    data.frame(rate = r, type = ty, mean = mean(est), se = sqrt(vt / ny) * tf,
               sd_env = sqrt(max(vt - mean(sv), 0)), n_years = ny, sampling_var = mean(sv),
               dist = if (ty == "survival") "beta" else "lognormal", min = NA_real_, max = NA_real_)
  })
  do.call(rbind, out)
}

# --- simulation (planted truth) -----------------------------------------------

.pva_truth_cache <- new.env(parent = emptyenv())

mx_default_definition <- function(params, design) {
  stages <- c("yearling", "adult")
  A <- matrix(c("0", "fec * s_juv", "s_sub", "s_ad"), 2, 2, byrow = TRUE)
  rates <- data.frame(rate = c("s_juv", "s_sub", "s_ad", "fec"),
                      type = c("survival", "survival", "survival", "fecundity"),
                      mean = c(params$s_juv, params$s_sub, params$s_ad, params$fec),
                      se = 0,
                      sd_env = c(params$sd_s_juv, params$sd_s_sub, params$sd_s_ad, params$sd_fec),
                      n_years = NA_real_, sampling_var = NA_real_,
                      dist = c("beta", "beta", "beta", "lognormal"), min = NA_real_, max = NA_real_)
  list(stages = stages, A = A, rates = rates)
}

mx_setup <- function(def, n0, threshold, count_stages, ceiling, corr, demographic) {
  k <- length(def$stages)
  rates <- def$rates
  for (cn in c("se", "n_years", "sampling_var", "min", "max")) if (is.null(rates[[cn]])) rates[[cn]] <- NA_real_
  if (is.null(rates$dist)) rates$dist <- ifelse(rates$type == "survival", "beta", "lognormal")
  rates$dist[is.na(rates$dist) | rates$dist == ""] <- ifelse(rates$type[is.na(rates$dist) | rates$dist == ""] == "survival", "beta", "lognormal")
  rates$se[is.na(rates$se)] <- 0
  rates$sd_env[is.na(rates$sd_env)] <- 0
  nr <- nrow(rates)
  C <- if (is.null(corr)) diag(nr) else if (length(corr) == 1) { m <- matrix(corr, nr, nr); diag(m) <- 1; m } else corr
  if (!is.null(dimnames(C)) && !is.null(rownames(C))) C <- C[rates$rate, rates$rate]
  L <- if (all(C == diag(nr))) NULL else chol(C)
  cs <- if (is.null(count_stages)) seq_len(k) else match(count_stages, def$stages)
  list(stages = def$stages, A = def$A, rates = rates, ent = mx_parse(def$A, rates),
       n0 = as.numeric(n0), threshold = threshold, count_idx = cs, ceiling = ceiling %||% Inf,
       cholcor = L, corr = C, demographic = demographic)
}

mx_truth <- function(mx, horizons, npaths, seed = 99) {
  key <- paste(deparse(list(mx$rates$mean, mx$rates$sd_env, mx$n0, mx$threshold, mx$ceiling, mx$corr,
                            mx$demographic, horizons, npaths, mx$A)), collapse = "")
  if (!is.null(.pva_truth_cache[[key]])) return(.pva_truth_cache[[key]])
  M <- matrix(mx$rates$mean, 1, dimnames = list(NULL, mx$rates$rate))
  S <- matrix(mx$rates$sd_env, 1, dimnames = list(NULL, mx$rates$rate))
  # in chunks to bound memory
  nch <- ceiling(npaths / 20000); cdf <- 0; ll <- 0
  for (ch in seq_len(nch)) {
    p <- mx_project(mx, M, S, ntraj = 20000, tmax = max(horizons), seed = seed + ch * 7, tul_years = 0)
    cdf <- cdf + p$cdf / nch; ll <- ll + p$loglams / nch
  }
  # long-run lambda_s: environment-only projection, 2000 traj x 400 yrs
  mxe <- mx; mxe$demographic <- FALSE; mxe$threshold <- -1; mxe$ceiling <- Inf
  pl <- mx_project(mxe, M, S, ntraj = 2000, tmax = 400, seed = seed + 1, burn = 20, tul_years = 0)
  vals <- stats::setNames(mx$rates$mean, mx$rates$rate)
  A <- mx_mean_matrix(mx$ent, vals, length(mx$stages))
  el <- mx_vr_elast(mx$ent, vals, length(mx$stages))
  tr <- c(stats::setNames(mx$rates$mean, mx$rates$rate),
          stats::setNames(mx$rates$sd_env, paste0("sd_env_", mx$rates$rate)),
          lambda_deterministic = mx_lambda1(A), lambda_stochastic = exp(pl$loglams),
          stats::setNames(cdf[1, horizons], paste0("p_qext_", horizons, "yr")),
          stats::setNames(el, paste0("elasticity_", names(el))))
  .pva_truth_cache[[key]] <- tr
  tr
}

matrix_simulate <- function(params, design, seed) {
  def <- mx_default_definition(params, design)
  mxt <- mx_setup(def, design$n0, design$threshold, design$count_stages, design$ceiling,
                  design$env_correlation %||% 0, demographic = TRUE)
  truth <- mx_truth(mxt, design$horizons, design$truth_paths %||% 1e5)
  # monitoring data: year-specific rates from the environmental distributions
  set.seed(seed)
  ny <- design$n_mon_years
  M <- matrix(def$rates$mean, 1); S <- matrix(def$rates$sd_env, 1)
  tabs <- mx_tables(mxt$rates, M, S)
  yr <- lapply(seq_len(ny), function(t) mx_vr_year(mxt$rates, M, S, tabs, 1L, mxt$cholcor))
  nm <- c(s_juv = design$n_marked_juv, s_sub = design$n_marked_sub, s_ad = design$n_marked_ad, fec = design$n_breeders)
  demog <- do.call(rbind, lapply(seq_len(ny), function(t) {
    v <- yr[[t]]
    data.frame(year = design$first_year + t - 1, rate = names(nm), n = unname(nm),
               events = c(stats::rbinom(3, nm[1:3], unlist(v[c("s_juv", "s_sub", "s_ad")])),
                          stats::rpois(1, nm[4] * v$fec)))
  }))
  data <- list(model = "matrix", stages = def$stages, A = def$A,
               rate_types = list(s_juv = "survival", s_sub = "survival", s_ad = "survival", fec = "fecundity"),
               rates = NULL, demography = demog, n0 = design$n0, threshold = design$threshold,
               count_stages = design$count_stages, ceiling = design$ceiling,
               env_correlation = if (isTRUE(design$fit_ignores_correlation)) 0 else design$env_correlation %||% 0,
               horizons = design$horizons, scenarios = NULL,
               meta = list(species = "simulated population", units = "females"))
  truth <- truth[c(def$rates$rate, paste0("sd_env_", def$rates$rate), "lambda_deterministic",
                   "lambda_stochastic", paste0("p_qext_", design$horizons, "yr"),
                   paste0("elasticity_", def$rates$rate))]
  list(data = data, truth = truth)
}

# --- fit -------------------------------------------------------------------

matrix_fit <- function(data, st) {
  rates <- data$rates
  if (is.null(rates)) {
    if (is.null(data$demography)) stop("Matrix PVA needs either a vital-rate table or annual demographic monitoring data.", call. = FALSE)
    rates <- mx_rates_from_demography(data$demography, data$rate_types)
  }
  rates <- rates[match(unique(rates$rate), rates$rate), ]
  hs <- sort(unique(as.integer(st$horizons %||% data$horizons)))
  thr <- st$threshold %||% data$threshold
  tmax <- max(c(hs, st$tmax %||% 0))
  if (identical(st$ignore_parameter_uncertainty, TRUE)) { rates$se <- 0; rates$n_years <- NA }
  demographic <- if (is.null(st$demographic)) TRUE else isTRUE(st$demographic)
  mx <- mx_setup(list(stages = data$stages, A = data$A, rates = rates), data$n0, thr, data$count_stages,
                 data$ceiling, data$env_correlation, demographic)
  k <- length(mx$stages)
  ndraw <- st$n_param_draws %||% 400; ntraj <- st$n_traj_per_draw %||% 250; npt <- st$n_traj_point %||% 10000
  seed <- st$seed %||% 1
  set.seed(seed)
  pd <- mx_param_draws(mx$rates, ndraw, corr = mx$corr)
  vals <- stats::setNames(mx$rates$mean, mx$rates$rate)
  run_scen <- function(mult = NULL, add = NULL) {
    r2 <- mx$rates
    for (nm in names(mult)) r2$mean[r2$rate == nm] <- r2$mean[r2$rate == nm] * mult[[nm]]
    for (nm in names(add)) r2$mean[r2$rate == nm] <- r2$mean[r2$rate == nm] + add[[nm]]
    r2$mean[r2$type == "survival"] <- pmin(r2$mean[r2$type == "survival"], 0.999)
    mx2 <- mx; mx2$rates <- r2
    d2 <- mx_param_draws(r2, ndraw, pd$u)
    p0 <- mx_project(mx2, matrix(r2$mean, 1, dimnames = list(NULL, r2$rate)),
                     matrix(r2$sd_env, 1, dimnames = list(NULL, r2$rate)), ntraj = npt, tmax = tmax,
                     seed = seed + 1, keep_paths = TRUE)
    pu <- mx_project(mx2, d2$M, d2$S, ntraj = ntraj, tmax = tmax, seed = seed + 2)
    list(point = p0, unc = pu, draws = d2, rates = r2)
  }
  base <- run_scen()
  qs <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  rows <- list()
  add_row <- function(p, est, draws) {
    ci <- qs(draws)
    rows[[length(rows) + 1]] <<- data.frame(param = p, estimate = est, sd = stats::sd(draws), lower = ci[1], upper = ci[2])
  }
  # vital rates and their environmental SDs
  for (r in seq_len(nrow(mx$rates))) {
    nm <- mx$rates$rate[r]
    add_row(nm, mx$rates$mean[r], base$draws$M[, r])
    if (is.finite(mx$rates$n_years[r]) && mx$rates$n_years[r] >= 2)
      add_row(paste0("sd_env_", nm), mx$rates$sd_env[r], base$draws$S[, r])
    else rows[[length(rows) + 1]] <- data.frame(param = paste0("sd_env_", nm), estimate = mx$rates$sd_env[r], sd = NA, lower = NA, upper = NA)
  }
  lam_d <- apply(base$draws$M, 1, function(v) mx_lambda1(mx_mean_matrix(mx$ent, stats::setNames(v, mx$rates$rate), k)))
  A0 <- mx_mean_matrix(mx$ent, vals, k)
  add_row("lambda_deterministic", mx_lambda1(A0), lam_d)
  add_row("lambda_stochastic", exp(base$point$loglams), exp(base$unc$loglams))
  add_row("lambda_stochastic_tuljapurkar", exp(base$point$tul), exp(base$unc$tul))
  for (h in hs) add_row(paste0("p_qext_", h, "yr"), base$point$cdf[1, h], base$unc$cdf[, h])
  el0 <- mx_vr_elast(mx$ent, vals, k)
  eld <- t(apply(base$draws$M, 1, function(v) mx_vr_elast(mx$ent, stats::setNames(as.list(v), mx$rates$rate), k)))
  for (r in names(el0)) add_row(paste0("elasticity_", r), el0[[r]], eld[, r])
  # management scenarios (common random numbers -> paired differences)
  scen <- list()
  for (sc in st$scenarios %||% data$scenarios %||% list()) {
    nm <- sc$name %||% "scenario"; slug <- gsub("[^a-z0-9]+", "_", tolower(nm))
    s <- run_scen(sc$multiply, sc$add)
    scen[[nm]] <- s
    add_row(paste0("scenario_", slug, "_lambda_stochastic"), exp(s$point$loglams), exp(s$unc$loglams))
    for (h in hs) {
      add_row(paste0("scenario_", slug, "_p_qext_", h, "yr"), s$point$cdf[1, h], s$unc$cdf[, h])
      add_row(paste0("scenario_", slug, "_change_p_qext_", h, "yr"), s$point$cdf[1, h] - base$point$cdf[1, h],
              s$unc$cdf[, h] - base$unc$cdf[, h])
    }
  }
  est <- do.call(rbind, rows)
  est$rhat <- NA_real_; est$ess_bulk <- NA_real_
  notes <- character()
  mcse <- sqrt(base$point$cdf[1, hs] * (1 - base$point$cdf[1, hs]) / npt)
  notes <- c(notes, sprintf("Monte Carlo SE of the point quasi-extinction probabilities: %s (%d trajectories).",
                            paste(sprintf("%.3f", mcse), collapse = ", "), npt))
  csum <- vapply(seq_len(k), function(j) sum(vapply(mx$ent, function(e) if (e$j == j && !e$repro) eval(e$expr, as.list(vals)) else 0, 0)), 0)
  if (any(csum > 1 + 1e-8)) notes <- c(notes, sprintf("WARNING: survival/transition entries out of stage(s) %s sum to more than 1 at the mean rates; demographic stochasticity then truncates them.",
                                                    paste(mx$stages[csum > 1], collapse = ", ")))
  if (any(mx$rates$type == "survival" & mx$rates$sd_env^2 >= 0.95 * mx$rates$mean * (1 - mx$rates$mean)))
    notes <- c(notes, "Some survival rates have environmental SD near the maximum possible for their mean (sqrt(m(1-m))); variance was capped.")
  span <- if (!is.null(data$demography)) diff(range(data$demography$year)) + 1 else max(mx$rates$n_years, na.rm = TRUE)
  if (is.finite(span) && any(hs > 2 * span))
    notes <- c(notes, sprintf("Horizons %s yr are long relative to the %g years of demographic data behind the rates: long-horizon extinction probabilities are imprecise (Fieberg & Ellner 2000); compare scenarios rather than reading absolute risk.",
                              paste(hs[hs > 2 * span], collapse = ", "), span))
  tab_cdf <- data.frame(year = seq_len(tmax), estimate = base$point$cdf[1, ],
                        lower = apply(base$unc$cdf, 2, stats::quantile, 0.025, names = FALSE),
                        upper = apply(base$unc$cdf, 2, stats::quantile, 0.975, names = FALSE),
                        lower50 = apply(base$unc$cdf, 2, stats::quantile, 0.25, names = FALSE),
                        upper50 = apply(base$unc$cdf, 2, stats::quantile, 0.75, names = FALSE))
  sens <- mx_sens(A0)
  list(estimates = est,
       diagnostics = list(converged = TRUE, notes = notes),
       object = list(model = "matrix", mx = mx, rates = mx$rates, A0 = A0, sens = sens,
                     elast_matrix = sens$S * A0 / sens$lambda, vr_elast = el0, vr_elast_draws = eld,
                     cdf = tab_cdf, base = base[c("point")], loglams_draws = base$unc$loglams,
                     tul_draws = base$unc$tul, lam_det_draws = lam_d,
                     scenarios = lapply(scen, function(s) list(point_cdf = s$point$cdf[1, ], unc_cdf = s$unc$cdf,
                                                               loglams = s$unc$loglams, point_loglam = s$point$loglams,
                                                               rates = s$rates)),
                     base_unc_cdf = base$unc$cdf, horizons = hs, tmax = tmax, threshold = thr,
                     settings = list(ndraw = ndraw, ntraj = ntraj, npt = npt, demographic = demographic)))
}

# --- read / check ------------------------------------------------------------

matrix_read_data <- function(config, dc) {
  st <- config$settings %||% list()
  stages <- unlist(dc$stages)
  if (is.null(stages)) stop("config data: stages: [names of the stages, in matrix order] is required.", call. = FALSE)
  k <- length(stages)
  A <- do.call(rbind, lapply(dc$matrix, function(r) as.character(unlist(r))))
  if (is.null(A) || any(dim(A) != c(k, k)))
    stop(sprintf("config data: matrix must have %d rows of %d entries (one per stage).", k, k), call. = FALSE)
  rates <- NULL; demog <- NULL; types <- NULL
  if (!is.null(dc$vital_rates)) {
    f <- pva_path(config, dc$vital_rates)
    if (!file.exists(f)) stop("Vital-rate file not found: ", f, call. = FALSE)
    rates <- utils::read.csv(f, strip.white = TRUE, stringsAsFactors = FALSE)
    need <- c("rate", "type", "mean")
    miss <- setdiff(need, names(rates))
    if (length(miss)) stop(sprintf("%s is missing column(s) %s. Columns present: %s.", dc$vital_rates,
                                   paste(miss, collapse = ", "), paste(names(rates), collapse = ", ")), call. = FALSE)
    rates$type <- tolower(rates$type)
    for (cn in c("se", "sd_env", "n_years", "sampling_var", "min", "max")) if (is.null(rates[[cn]])) rates[[cn]] <- NA_real_
    if (is.null(rates$dist)) rates$dist <- NA_character_
  }
  if (!is.null(dc$demography)) {
    f <- pva_path(config, dc$demography)
    if (!file.exists(f)) stop("Demography file not found: ", f, call. = FALSE)
    demog <- utils::read.csv(f, strip.white = TRUE, stringsAsFactors = FALSE)
    miss <- setdiff(c("year", "rate", "n", "events"), names(demog))
    if (length(miss)) stop(sprintf("%s is missing column(s) %s (need year, rate, n, events).", dc$demography, paste(miss, collapse = ", ")), call. = FALSE)
    types <- dc$rate_types %||% list()
  }
  n0 <- unlist(dc$initial_abundance)
  if (!is.null(names(n0)) && all(stages %in% names(n0))) n0 <- n0[stages]
  else if (!is.null(names(n0)) && any(nzchar(names(n0))))
    stop(sprintf("config data: initial_abundance names (%s) must match the stage names (%s); not found: %s.",
                 paste(names(n0), collapse = ", "), paste(stages, collapse = ", "), paste(setdiff(stages, names(n0)), collapse = ", ")), call. = FALSE)
  corr <- dc$env_correlation %||% 0
  if (is.character(corr)) {
    cm <- as.matrix(utils::read.csv(pva_path(config, corr), row.names = 1, check.names = FALSE))
    corr <- cm
  }
  list(model = "matrix", stages = stages, A = A, rates = rates, demography = demog, rate_types = types,
       n0 = as.numeric(n0), threshold = st$threshold, count_stages = unlist(dc$count_stages),
       ceiling = dc$ceiling, env_correlation = corr, horizons = st$horizons %||% c(10, 25, 50),
       scenarios = st$scenarios,
       meta = list(species = dc$species %||% "the population", units = dc$units %||% "individuals",
                   file = dc$vital_rates %||% dc$demography))
}

matrix_check_data <- function(data) {
  pr <- character()
  k <- length(data$stages)
  if (length(data$n0) != k || any(!is.finite(data$n0)) || any(data$n0 < 0))
    pr <- c(pr, sprintf("initial_abundance must give %d non-negative numbers (one per stage: %s).", k, paste(data$stages, collapse = ", ")))
  rates <- data$rates
  if (is.null(rates) && !is.null(data$demography)) {
    d <- data$demography
    # one type per row (unlist() would drop rates without a type and misalign the rows)
    rt <- vapply(as.character(d$rate), function(r) tolower(as.character(data$rate_types[[r]] %||% "survival"))[1], "")
    badn <- which(!is.finite(d$n) | !is.finite(d$events) | d$n <= 0 | d$events < 0)
    if (length(badn)) pr <- c(pr, sprintf("Missing, zero or negative n / events in row(s) %s of the monitoring file (%s); fix or delete those rows.",
                                         paste(head(badn, 10), collapse = ", "), paste(head(sprintf("%s %s", d$year[badn], d$rate[badn]), 10), collapse = "; ")))
    over <- which(d$events > d$n & rt == "survival")
    if (length(over)) pr <- c(pr, sprintf("Survival rows with more survivors (events) than animals at risk (n): %s.",
                                         paste(head(sprintf("year %s, %s: %s of %s", d$year[over], d$rate[over], d$events[over], d$n[over]), 10), collapse = "; ")))
    used <- unique(unlist(lapply(data$A, function(s) if (s %in% c("", "0")) character() else all.vars(parse(text = s)[[1]]))))
    unused <- setdiff(unique(d$rate), used)
    for (u in unused) {
      near <- used[utils::adist(u, used) <= 2]
      pr <- c(pr, sprintf("Rate '%s' (%d rows) in the monitoring file is not used in the matrix%s", u, sum(d$rate == u),
                          if (length(near)) sprintf("; is it a misspelling of '%s'?", near[1]) else "; check the spelling or remove those rows."))
    }
    miss <- setdiff(intersect(unique(d$rate), used), names(data$rate_types))
    if (length(miss)) pr <- c(pr, sprintf("Give data: rate_types for %s (survival or fecundity).", paste(miss, collapse = ", ")))
    okrow <- setdiff(seq_len(nrow(d)), badn)
    rates <- tryCatch(mx_rates_from_demography(d[okrow, , drop = FALSE], data$rate_types), error = function(e) NULL)
    yrs <- table(d$rate)
    if (any(yrs < 3)) pr <- c(pr, sprintf("Rate(s) %s have fewer than 3 years of data: environmental variation cannot be estimated.", paste(names(yrs)[yrs < 3], collapse = ", ")))
  }
  if (is.null(rates)) return(c(pr, "No vital rates: give data: vital_rates (a CSV of rate estimates) or data: demography (annual monitoring)."))
  if (any(!rates$type %in% c("survival", "fecundity"), na.rm = TRUE))
    pr <- c(pr, sprintf("Vital-rate type must be 'survival' (any probability) or 'fecundity'; found: %s.", paste(unique(rates$type[!rates$type %in% c("survival", "fecundity")]), collapse = ", ")))
  sv <- rates$type == "survival"
  badm <- which(!is.finite(rates$mean))
  if (length(badm)) pr <- c(pr, sprintf("Missing or non-numeric mean for rate(s) %s.", paste(rates$rate[badm], collapse = ", ")))
  bads <- which(sv & is.finite(rates$mean) & (rates$mean <= 0 | rates$mean >= 1))
  if (length(bads)) pr <- c(pr, sprintf("Survival/probability rates must lie strictly between 0 and 1: %s.", paste(sprintf("%s = %g", rates$rate[bads], rates$mean[bads]), collapse = ", ")))
  if (any(!sv & is.finite(rates$mean) & rates$mean <= 0)) pr <- c(pr, "Fecundity rates must be positive.")
  vars <- unique(unlist(lapply(data$A, function(s) if (s %in% c("", "0")) character() else all.vars(parse(text = s)[[1]]))))
  miss <- setdiff(vars, rates$rate)
  if (length(miss)) pr <- c(pr, sprintf("The matrix uses rate(s) %s that are not in the vital-rate data.", paste(miss, collapse = ", ")))
  if (is.null(data$threshold) || !is.numeric(data$threshold)) pr <- c(pr, "settings: threshold (quasi-extinction threshold, number of individuals in the counted stages) is required.")
  if (!is.null(data$count_stages) && !all(data$count_stages %in% data$stages))
    pr <- c(pr, sprintf("count_stages %s are not stage names (%s).", paste(setdiff(data$count_stages, data$stages), collapse = ", "), paste(data$stages, collapse = ", ")))
  C <- data$env_correlation
  if (is.matrix(C)) {
    if (!isTRUE(all.equal(C, t(C), check.attributes = FALSE))) pr <- c(pr, "The environmental correlation matrix is not symmetric.")
    else if (min(eigen(C, only.values = TRUE)$values) <= 0) pr <- c(pr, "The environmental correlation matrix is not positive definite (some correlations are impossible together; Morris & Doak 2002 discuss this for Hudsonia).")
    if (!is.null(rownames(C)) && !all(rates$rate %in% rownames(C))) pr <- c(pr, "Row/column names of the correlation matrix must include every vital rate.")
  } else if (is.numeric(C) && (C <= -1 / (nrow(rates) - 1) || C >= 1)) pr <- c(pr, "env_correlation must be between -1/(n_rates-1) and 1.")
  pr
}
