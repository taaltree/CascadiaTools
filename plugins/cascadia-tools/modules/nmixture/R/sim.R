# nmixture: simulation of repeated counts (planted truth) and goodness of fit.

nm_rmix <- function(mu, mixture, size = NA, zi = 0) {
  n <- length(mu)
  switch(mixture,
         P = stats::rpois(n, mu),
         NB = stats::rnbinom(n, size = size, mu = mu),
         ZIP = stats::rbinom(n, 1, 1 - zi) * stats::rpois(n, mu / (1 - zi)))
}

nm_simulate <- function(params, design, seed) {
  set.seed(seed)
  model <- design$model %||% "static"
  M <- design$n_sites; J <- design$n_visits
  T <- if (model == "dynamic") design$n_seasons else 1
  dist <- design$abundance_dist %||% "P"
  zstd <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(as.vector(x), na.rm = TRUE)
  w <- array(stats::rbinom(M * J * T, 1, 1 - (design$prop_missing %||% 0)), c(M, J, T))
  for (i in seq_len(M)) for (t in seq_len(T)) if (sum(w[i, , t]) == 0) w[i, sample.int(J, 1), t] <- 1
  forest <- zstd(stats::rnorm(M))
  wind <- array(stats::rnorm(M * J * T), c(M, J, T)); wind[w == 0] <- NA; wind <- zstd(wind)
  mean_ab <- if (model == "static") params$abundance_at_average_site else params$initial_abundance_per_site
  mu <- mean_ab * exp(params$lambda_beta_forest * forest)
  zi <- if (dist == "ZIP") params$zero_inflation_prob else 0
  N <- matrix(0L, M, T)
  N[, 1] <- nm_rmix(mu, dist, params$nb_size, zi)
  if (T > 1) for (t in 2:T)
    N[, t] <- stats::rbinom(M, N[, t - 1], params$apparent_survival) + stats::rpois(M, params$recruitment_per_site)
  # detection: covariate + optional unmodelled heterogeneity (negative controls)
  site_re <- if ((design$p_site_sd %||% 0) > 0) stats::rnorm(M, 0, design$p_site_sd) else rep(0, M)
  visit_re <- if ((design$p_visit_sd %||% 0) > 0) array(stats::rnorm(M * J * T, 0, design$p_visit_sd), c(M, J, T)) else 0
  eta <- stats::qlogis(params$detection_prob) + params$p_beta_wind * ifelse(is.na(wind), 0, wind) + site_re + visit_re
  p <- stats::plogis(eta)
  NN <- aperm(array(N, c(M, T, J)), c(1, 3, 2))
  y <- array(stats::rbinom(M * J * T, NN, p), c(M, J, T))
  # optional double counting: each detected animal counted twice with this probability
  if ((design$double_count_prob %||% 0) > 0) y <- y + array(stats::rbinom(M * J * T, y, design$double_count_prob), c(M, J, T))
  y[w == 0] <- NA
  ids <- sprintf("S%03d", seq_len(M))
  rr <- c(lambda_beta_forest = params$lambda_beta_forest, rate_ratio_forest = exp(params$lambda_beta_forest),
          p_beta_wind = params$p_beta_wind, odds_ratio_detection_wind = exp(params$p_beta_wind))
  extra <- c(if (dist == "NB") c(nb_size = params$nb_size), if (dist == "ZIP") c(zero_inflation_prob = zi))
  if (model == "static") {
    data <- list(model = "static", y = matrix(y, M, J), site_covs = data.frame(forest = forest),
                 obs_covs = list(wind = matrix(wind, M, J)),
                 formulas = list(lambda = "~ forest", p = "~ wind"), site_ids = ids)
    truth <- c(abundance_at_average_site = mean_ab, detection_prob = params$detection_prob, rr, extra,
               total_abundance = sum(N), mean_abundance_per_site = mean(N))
  } else {
    data <- list(model = "dynamic", y = y, site_covs = data.frame(forest = forest),
                 obs_covs = list(wind = wind), yearly_covs = list(),
                 formulas = list(lambda = "~ forest", p = "~ wind", gamma = "~ 1", omega = "~ 1"),
                 site_ids = ids, season_labels = as.character(seq_len(T)))
    tot <- colSums(N)
    truth <- c(initial_abundance_per_site = mean_ab, recruitment_per_site = params$recruitment_per_site,
               apparent_survival = params$apparent_survival, detection_prob = params$detection_prob, rr, extra,
               stats::setNames(tot, paste0("total_abundance_t", seq_len(T))),
               abundance_ratio_last_first = tot[T] / max(tot[1], 1))
  }
  list(data = data, truth = truth)
}

# Simulate a new count array from fitted parameters on the observed design.
nm_sim_from <- function(u, des, mixture) {
  M <- des$M; J <- des$J; T <- des$T
  N <- matrix(0L, M, T)
  N[, 1] <- switch(mixture, P = stats::rpois(M, u$lam), NB = stats::rnbinom(M, size = u$size, mu = u$lam),
                   ZIP = stats::rbinom(M, 1, 1 - u$zi) * stats::rpois(M, u$lam))
  if (T > 1) for (t in 2:T) N[, t] <- stats::rbinom(M, N[, t - 1], u$om[, t - 1]) + stats::rpois(M, u$gam[, t - 1])
  NN <- aperm(array(N, c(M, T, J)), c(1, 3, 2))
  y <- array(stats::rbinom(M * J * T, NN, u$p), c(M, J, T))
  y[!des$w] <- NA
  list(y = y, N = N)
}

# Marginal expected counts E[y] (as unmarked's fitted()).
nm_expected <- function(u, des, mixture) {
  M <- des$M; J <- des$J; T <- des$T
  EN <- matrix(0, M, T)
  EN[, 1] <- u$lam * (1 - u$zi)
  if (T > 1) for (t in 2:T) EN[, t] <- EN[, t - 1] * u$om[, t - 1] + u$gam[, t - 1]
  aperm(array(EN, c(M, T, J)), c(1, 3, 2)) * u$p
}

nm_fitstats <- function(y, E) {
  ok <- !is.na(y) & E > 0
  c(chisq = sum((y[ok] - E[ok])^2 / E[ok]), ft = sum((sqrt(y[ok]) - sqrt(E[ok]))^2))
}

# Parametric bootstrap (likelihood engines) or posterior predictive check (nimble).
nm_gof <- function(eng, des, settings) {
  mixture <- settings$mixture; dyn <- settings$dynamics
  nb <- settings$gof_nboot
  if (eng$engine == "nimble") {
    idx <- sample.int(nrow(eng$draws), min(nb, nrow(eng$draws)))
    st <- t(vapply(idx, function(d) {
      th <- eng$draws[d, ]; u <- nm_unpack(th, des, mixture, dyn)
      Nd <- nm_N_draws(matrix(th, 1, dimnames = list(NULL, names(th))), des, eng$K, mixture, dyn)[1, , ]
      NN <- aperm(array(Nd, c(des$M, des$T, des$J)), c(1, 3, 2))
      E <- NN * u$p
      yr <- array(stats::rbinom(length(E), NN, u$p), dim(E)); yr[!des$w] <- NA
      c(nm_fitstats(des$y, E), nm_fitstats(yr, E))
    }, numeric(4)))
    return(list(method = "posterior predictive check (conditional on N)", nboot = nrow(st),
                chisq_obs = mean(st[, 1]), chisq_rep = st[, 3], ft_obs = mean(st[, 2]), ft_rep = st[, 4],
                c_hat = mean(st[, 1]) / mean(st[, 3]),
                p_chisq = mean(st[, 3] >= st[, 1]), p_ft = mean(st[, 4] >= st[, 2])))
  }
  u0 <- nm_unpack(eng$theta, des, mixture, dyn)
  obs <- nm_fitstats(des$y, nm_expected(u0, des, mixture))
  ncores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  st0 <- merge_list(settings, list(n_draws = 2, K = eng$K))
  one <- function(b) {
    set.seed(settings$seed * 1000 + b)
    d2 <- des; d2$y <- nm_sim_from(u0, des, mixture)$y
    e <- tryCatch(if (eng$engine == "unmarked") nm_engine_unmarked(d2, st0, K = eng$K) else nm_engine_rlik(d2, st0, K = eng$K),
                  error = function(e) NULL)
    if (is.null(e)) return(c(NA, NA))
    nm_fitstats(d2$y, nm_expected(nm_unpack(e$theta, d2, mixture, dyn), d2, mixture))
  }
  # Forked workers hang in unmarked::pcountOpen (threaded C++ code initialised in the parent;
  # diagnosed with dev/nmixture/test_gof_dyn.R), so the dynamic model uses socket workers.
  bs <- if (ncores <= 1) lapply(seq_len(nb), one)
        else if (des$model == "dynamic") {
          cl <- parallel::makeCluster(min(ncores, nb)); on.exit(parallel::stopCluster(cl), add = TRUE)
          parallel::clusterCall(cl, function(lp, root) {   # core helpers live in globalenv; the module env travels with `one`
            .libPaths(lp); Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", ECOTOOLS_ROOT = root)
            source(file.path(root, "core", "load.R")); NULL }, .libPaths(), eco_root())
          parallel::parLapply(cl, seq_len(nb), one)
        } else parallel::mclapply(seq_len(nb), one, mc.cores = ncores)
  bs <- do.call(rbind, bs); bs <- bs[stats::complete.cases(bs), , drop = FALSE]
  list(method = "parametric bootstrap", nboot = nrow(bs), chisq_obs = obs[["chisq"]], chisq_rep = bs[, 1],
       ft_obs = obs[["ft"]], ft_rep = bs[, 2], c_hat = obs[["chisq"]] / mean(bs[, 1]),
       p_chisq = mean(bs[, 1] >= obs[["chisq"]]), p_ft = mean(bs[, 2] >= obs[["ft"]]))
}
