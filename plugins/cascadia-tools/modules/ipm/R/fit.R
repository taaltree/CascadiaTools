# ipm: model fitting with Nimble, derived quantities, checks.

ipm_settings <- list(
  components = "ipm",          # "ipm" (joint) or "separate" (CJS + productivity only; no counts)
  immigration = TRUE,          # estimate immigration as a latent rate
  temporal_variation = TRUE,   # random year effects on every rate
  count_model = "poisson",     # "poisson" or "normal" (normal: sd estimated)
  sex_ratio = 0.5,             # proportion of females among fledglings
  fmax = 20,                   # upper limit of the uniform prior on mean fecundity
  niter = 150000, nburnin = 30000, thin = 30, nchains = 3, seed = 1,
  max_draws_per_chain = 4000,
  slice_hyper = TRUE,          # slice samplers for mean rates and temporal SDs (better mixing)
  compare_separate = TRUE,     # also fit the separate analyses and compare precision
  ppc_draws = 600)

ipm_fast <- merge_list(ipm_settings, list(niter = 35000, nburnin = 5000, thin = 6,
                                          max_draws_per_chain = 5000,
                                          compare_separate = FALSE, ppc_draws = 0))

ipm_flags <- function(data, settings) {
  ipm <- !identical(settings$components, "separate")
  eff <- data$effects
  K <- function(r) length(eff[[r]])
  list(IPM = as.integer(ipm), IMM = as.integer(ipm && isTRUE(settings$immigration)),
       AGEF = as.integer(isTRUE(data$age_fecundity)), RE = as.integer(isTRUE(settings$temporal_variation)),
       CNORM = as.integer(identical(settings$count_model, "normal")),
       K_phij = K("juvenile_survival"), K_phia = K("adult_survival"), K_f = K("fecundity"),
       K_om = if (ipm && isTRUE(settings$immigration)) K("immigration") else 0L)
}

# standardised covariate design matrices (original-unit scaling kept for plots)
ipm_X <- function(data, rate, rows) {
  cv <- data$effects[[rate]]
  if (!length(cv)) return(NULL)
  X <- sapply(cv, function(v) data$covariates[[v]])
  X <- matrix(X, nrow = length(data$years))
  sc <- lapply(seq_along(cv), function(k) list(center = mean(X[, k], na.rm = TRUE), scale = stats::sd(X[, k], na.rm = TRUE)))
  for (k in seq_along(cv)) X[, k] <- (X[, k] - sc[[k]]$center) / sc[[k]]$scale
  X[is.na(X)] <- 0
  list(X = X[rows, , drop = FALSE], scaling = stats::setNames(sc, cv))
}

ipm_init_N <- function(counts, phia0, om0, imm, U) {
  T <- length(counts)
  C <- counts
  if (all(is.na(C))) C <- rep(20, T)
  ok <- which(!is.na(C))
  C <- stats::approx(ok, C[ok], xout = seq_len(T), rule = 2)$y
  C <- pmax(round(C), 3)
  N <- matrix(0, 3, T)
  N[2, 1] <- max(1, round(0.45 * C[1])); N[3, 1] <- if (imm) max(1, round(0.15 * C[1])) else 0
  N[1, 1] <- max(1, C[1] - N[2, 1] - N[3, 1])
  N[, 1] <- pmin(N[, 1], U); if (!imm) N[3, 1] <- 0
  for (t in 1:(T - 1)) {
    Nt <- sum(N[, t])
    N[2, t + 1] <- min(Nt, round(phia0 * Nt))
    N[3, t + 1] <- if (imm) round(om0 * Nt) else 0
    N[1, t + 1] <- max(0, C[t + 1] - N[2, t + 1] - N[3, t + 1])
  }
  N
}

ipm_nimble_inputs <- function(data, settings) {
  T <- length(data$years)
  fl <- ipm_flags(data, settings)
  cst <- list(T = T)
  dat <- list(marr_j = data$marr_j, marr_a = data$marr_a,
              rel_j = rowSums(data$marr_j), rel_a = rowSums(data$marr_a),
              J1 = data$prod$J1, B1 = data$prod$B1, fmax = settings$fmax)
  if (fl$AGEF == 1) { dat$J2 <- data$prod$J2; dat$B2 <- data$prod$B2 }
  scaling <- list()
  for (r in list(c("juvenile_survival", "X_phij", "T1"), c("adult_survival", "X_phia", "T1"),
                 c("fecundity", "X_f", "T"), c("immigration", "X_om", "T1"))) {
    if (r[1] == "immigration" && fl$K_om == 0) next
    x <- ipm_X(data, r[1], if (r[3] == "T") seq_len(T) else seq_len(T - 1))
    if (!is.null(x)) { dat[[r[2]]] <- x$X; scaling <- c(scaling, x$scaling) }
  }
  C1 <- data$counts[!is.na(data$counts)][1]
  upper <- max(20, ceiling(2 * C1))
  U <- as.integer(ceiling(upper / 50) * 50)
  if (fl$IPM == 1) {
    dat$C <- data$counts
    dat$pinit <- c(rep(1 / upper, upper), rep(0, U - upper))
    dat$rho <- settings$sex_ratio
    if (fl$CNORM == 1) dat$sig_max <- max(10, max(data$counts, na.rm = TRUE))
    cst$U <- U
  }
  fbar <- function(J, B) if (sum(B) > 0) sum(J) / sum(B) else 3
  f1b <- fbar(data$prod$J1, data$prod$B1); f2b <- fbar(data$prod$J2, data$prod$B2)
  # data-based log-fecundity deviations for initial values (jittered; years
  # without broods get a random draw)
  eps_f <- function(J, B, mf) {
    e <- ifelse(B > 0, log((J + 0.5) / (B + 0.5 / mf)) - log(mf), NA)
    e <- e + stats::rnorm(length(e), 0, 0.05)
    e[is.na(e)] <- stats::rnorm(sum(is.na(e)), 0, 0.1)
    pmin(pmax(e, -2), 2)
  }
  inits <- function() {
    phij0 <- stats::runif(1, 0.15, 0.4); phia0 <- stats::runif(1, 0.35, 0.65)
    om0 <- stats::runif(1, 0.05, 0.25); p0 <- stats::runif(1, 0.3, 0.8)
    ini <- list(mean_phij = phij0, mean_phia = phia0, p = p0,
                mean_f1 = min(settings$fmax * 0.9, max(0.2, f1b * stats::runif(1, 0.9, 1.1))))
    if (fl$AGEF == 1) ini$mean_f2 <- min(settings$fmax * 0.9, max(0.2, f2b * stats::runif(1, 0.9, 1.1)))
    if (fl$RE == 1) {
      ini$sd_phij <- stats::runif(1, 0.1, 0.5); ini$sd_phia <- stats::runif(1, 0.1, 0.5)
      ini$sd_f1 <- stats::runif(1, 0.05, 0.3)
      # Year effects must NOT all start at exactly 0: with the centred
      # parameterisation, eps = 0 makes sd -> 0 an absorbing point for the
      # sampler (funnel collapse; one chain stuck at sd ~ 1e-35 was the cause
      # of R-hat ~1.7 in about a third of replicates). Start survival effects
      # jittered around 0 and fecundity effects at their data-based values.
      ini$eps_phij <- stats::rnorm(T - 1, 0, ini$sd_phij); ini$eps_phia <- stats::rnorm(T - 1, 0, ini$sd_phia)
      ini$eps_f1 <- eps_f(data$prod$J1, data$prod$B1, ini$mean_f1)
      ini$sd_f1 <- max(0.05, stats::sd(ini$eps_f1))
      if (fl$AGEF == 1) {
        ini$eps_f2 <- eps_f(data$prod$J2, data$prod$B2, ini$mean_f2)
        ini$sd_f2 <- max(0.05, stats::sd(ini$eps_f2))
      }
    }
    if (fl$K_phij > 0) ini$beta_phij <- rep(0, fl$K_phij)
    if (fl$K_phia > 0) ini$beta_phia <- rep(0, fl$K_phia)
    if (fl$K_f > 0) ini$beta_f <- rep(0, fl$K_f)
    if (fl$IPM == 1) {
      N <- ipm_init_N(data$counts, phia0, om0, fl$IMM == 1, upper)
      ini$N1 <- N[1, ]; ini$N2 <- N[2, ]
      if (fl$IMM == 1) {
        ini$N3 <- N[3, ]; ini$mean_omega <- om0
        if (fl$RE == 1) { ini$sd_omega <- stats::runif(1, 0.1, 0.5); ini$eps_om <- stats::rnorm(T - 1, 0, ini$sd_omega) }
        if (fl$K_om > 0) ini$beta_om <- rep(0, fl$K_om)
      }
      if (fl$CNORM == 1) ini$sigma_count <- stats::runif(1, 1, 0.5 * dat$sig_max)
      if (anyNA(data$counts)) {
        ci <- data$counts; ci[is.na(ci)] <- colSums(N)[is.na(ci)]; ini$C <- ci
      }
    }
    ini
  }
  mon <- c("mean_phij", "mean_phia", "mean_f1", "p", "phij", "phia", "f1")
  if (fl$AGEF == 1) mon <- c(mon, "mean_f2", "f2")
  if (fl$RE == 1) mon <- c(mon, "sd_phij", "sd_phia", "sd_f1", if (fl$AGEF == 1) "sd_f2")
  if (fl$K_phij > 0) mon <- c(mon, "beta_phij")
  if (fl$K_phia > 0) mon <- c(mon, "beta_phia")
  if (fl$K_f > 0) mon <- c(mon, "beta_f")
  if (fl$IPM == 1) {
    mon <- c(mon, "N1", "N2")
    if (fl$IMM == 1) mon <- c(mon, "N3", "mean_omega", "omega", if (fl$RE == 1) "sd_omega", if (fl$K_om > 0) "beta_om")
    if (fl$CNORM == 1) mon <- c(mon, "sigma_count")
  }
  key <- paste("ipm", T, paste(names(fl), unlist(fl), sep = "=", collapse = ","), if (fl$IPM == 1) U,
               paste(which(is.na(data$counts)), collapse = "."), sep = "|")
  list(constants = cst, data = dat, inits = inits, monitors = mon, flags = fl, key = key,
       scaling = scaling, U = U)
}

# Build per-chain matrices of reported quantities (biological names + derived)
ipm_report_chain <- function(m, data, fl, settings) {
  T <- length(data$years); yrs <- data$years
  g <- function(b, n) { cn <- paste0(b, "[", seq_len(n), "]"); m[, cn, drop = FALSE] }
  fec1 <- if (fl$AGEF == 1) "fecundity_yearling" else "fecundity"
  out <- cbind(juvenile_survival = m[, "mean_phij"], adult_survival = m[, "mean_phia"])
  out <- cbind(out, stats::setNames(data.frame(m[, "mean_f1"]), fec1))
  if (fl$AGEF == 1) out <- cbind(out, fecundity_adult = m[, "mean_f2"])
  out <- cbind(out, recapture_prob = m[, "p"])
  if (fl$RE == 1) {
    out <- cbind(out, juvenile_survival_sd = m[, "sd_phij"], adult_survival_sd = m[, "sd_phia"])
    out <- cbind(out, stats::setNames(data.frame(m[, "sd_f1"]), paste0(fec1, "_sd")))
    if (fl$AGEF == 1) out <- cbind(out, fecundity_adult_sd = m[, "sd_f2"])
  }
  addb <- function(out, node, K, rate) {
    if (K == 0) return(out)
    b <- m[, paste0(node, "[", seq_len(K), "]"), drop = FALSE]
    colnames(b) <- paste0(rate, "_beta_", data$effects[[rate]])
    cbind(out, b)
  }
  out <- addb(out, "beta_phij", fl$K_phij, "juvenile_survival")
  out <- addb(out, "beta_phia", fl$K_phia, "adult_survival")
  out <- addb(out, "beta_f", fl$K_f, "fecundity")
  out <- as.matrix(out)
  ann <- list(phij = g("phij", T - 1), phia = g("phia", T - 1), f1 = g("f1", T))
  ann$f2 <- if (fl$AGEF == 1) g("f2", T) else ann$f1
  if (fl$IPM == 1) {
    if (fl$IMM == 1) {
      # immigration_rate = mean of the annual rates over the study years (the
      # quantity the counts inform); immigration_rate_median = exp(mu), the
      # rate in a typical year. The median depends on the poorly identified
      # immigration SD (VALIDATION.md) and is not interpreted on its own.
      out <- cbind(out, immigration_rate = rowMeans(g("omega", T - 1)),
                   immigration_rate_median = m[, "mean_omega"])
      if (fl$RE == 1) out <- cbind(out, immigration_sd = m[, "sd_omega"])
      out <- addb(out, "beta_om", fl$K_om, "immigration")
      ann$omega <- g("omega", T - 1)
    }
    if (fl$CNORM == 1) out <- cbind(out, count_sd = m[, "sigma_count"])
    N1 <- g("N1", T); N2 <- g("N2", T); N3 <- if (fl$IMM == 1) g("N3", T) else N1 * 0
    Nt <- N1 + N2 + N3
    lam <- Nt[, -1, drop = FALSE] / pmax(Nt[, -T, drop = FALSE], 1e-9)
    out <- cbind(out, N_first = Nt[, 1], N_final = Nt[, T], N_mean = rowMeans(Nt),
                 lambda_geomean = (Nt[, T] / pmax(Nt[, 1], 1e-9))^(1 / (T - 1)))
    if (fl$IMM == 1) out <- cbind(out, immigrants_mean = rowMeans(N3[, -1, drop = FALSE]))
    if (fl$RE == 1) {
      lt <- ipm_ltre_draws(c(ann, list(N1 = N1, N2 = N2, N3 = N3)),
                           means = list(phij = m[, "mean_phij"], phia = m[, "mean_phia"], f1 = m[, "mean_f1"],
                                        f2 = if (fl$AGEF == 1) m[, "mean_f2"] else m[, "mean_f1"],
                                        omega = if (fl$IMM == 1) m[, "mean_omega"] else NULL),
                           flags = fl, rho = settings$sex_ratio)
      out <- cbind(out, as.matrix(lt$contrib))
    }
    colnames(Nt) <- paste0("N_total[", yrs, "]")
    colnames(lam) <- paste0("lambda[", yrs[-T], "]")
    out <- cbind(out, Nt, lam)
    a1 <- N1; colnames(a1) <- paste0("local_recruits[", yrs, "]")
    a2 <- N2; colnames(a2) <- paste0("surviving_adults[", yrs, "]")
    out <- cbind(out, a1, a2)
    if (fl$IMM == 1) {
      a3 <- N3; colnames(a3) <- paste0("immigrants[", yrs, "]")
      om <- ann$omega; colnames(om) <- paste0("immigration_rate[", yrs[-T], "]")
      out <- cbind(out, a3, om)
    }
  }
  pj <- ann$phij; colnames(pj) <- paste0("juvenile_survival[", yrs[-T], "]")
  pa <- ann$phia; colnames(pa) <- paste0("adult_survival[", yrs[-T], "]")
  f1 <- ann$f1; colnames(f1) <- paste0(fec1, "[", yrs, "]")
  out <- cbind(out, pj, pa, f1)
  if (fl$AGEF == 1) { f2 <- ann$f2; colnames(f2) <- paste0("fecundity_adult[", yrs, "]"); out <- cbind(out, f2) }
  out
}

ipm_thin_idx <- function(n, k) if (n <= k) seq_len(n) else unique(round(seq(1, n, length.out = k)))

ipm_run <- function(data, settings) {
  ipm_define_nimble()
  inp <- ipm_nimble_inputs(data, settings)
  code <- ipm_nimble_code(inp$flags)
  # Slice samplers for the hyperparameters (mean rates, temporal SDs): the
  # temporal SDs (immigration SD especially) mix slowly with the default
  # adaptive random-walk samplers.
  conf_fun <- if (isTRUE(settings$slice_hyper)) function(conf, m) {
    hyp <- intersect(c("mean_phij", "mean_phia", "mean_f1", "mean_f2", "mean_omega", "p",
                       "sd_phij", "sd_phia", "sd_f1", "sd_f2", "sd_omega", "sigma_count"),
                     m$getNodeNames(stochOnly = TRUE, includeData = FALSE))
    for (n in hyp) { conf$removeSamplers(n); conf$addSampler(n, type = "slice") }
  } else NULL
  res <- nimble_run(code, inp$constants, inp$data, inp$inits, inp$monitors,
                    niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                    thin = settings$thin, cache_key = paste0(inp$key, "|slice=", isTRUE(settings$slice_hyper)),
                    configure = conf_fun, seed = settings$seed,
                    summary_params = c("mean_phij", "mean_phia"))
  chains <- lapply(res$samples, function(m) {
    m <- m[ipm_thin_idx(nrow(m), settings$max_draws_per_chain), , drop = FALSE]
    ipm_report_chain(m, data, inp$flags, settings)
  })
  list(chains = chains, raw = res$samples, inp = inp, t_total = res$diagnostics$t_total,
       t_setup = res$diagnostics$t_setup)
}

ipm_fit <- function(data, settings) {
  settings <- merge_list(ipm_settings, settings)
  set.seed(settings$seed)
  r <- ipm_run(data, settings)
  fl <- r$inp$flags
  s <- summarise_draws_eco(r$chains)
  est <- s$estimates
  annual <- grepl("\\[", est$param)
  key_rows <- !annual
  dg <- list(converged = s$diagnostics$converged, rhat_max = s$diagnostics$rhat_max,
             ess_min = s$diagnostics$ess_min, t_total = r$t_total,
             rhat_max_summary = max(est$rhat[key_rows], na.rm = TRUE),
             ess_min_summary = min(est$ess_bulk[key_rows], na.rm = TRUE))
  notes <- data$notes
  if (is.finite(dg$rhat_max) && (dg$rhat_max >= 1.01 || dg$ess_min < 400))
    notes <- c(notes, sprintf("Real-data standard not met (max R-hat %.3f, min bulk ESS %.0f over all %d reported quantities; need < 1.01 and >= 400): run longer (settings niter, nburnin, thin).",
                              dg$rhat_max, dg$ess_min, nrow(est)))
  draws <- stack_draws(r$chains)
  obj <- list(chains = r$chains, draws = draws, flags = fl, scaling = r$inp$scaling, settings = settings,
              years = data$years)
  if (fl$IPM == 1 && fl$RE == 1) {
    # transient LTRE (mean over draws) for the figure
    idx <- ipm_thin_idx(nrow(draws), 1000)
    rm <- stack_draws(r$raw)[ipm_thin_idx(nrow(stack_draws(r$raw)), 1000), , drop = FALSE]
    T <- length(data$years)
    g <- function(b, n) rm[, paste0(b, "[", seq_len(n), "]"), drop = FALSE]
    d <- list(phij = g("phij", T - 1), phia = g("phia", T - 1), f1 = g("f1", T),
              f2 = if (fl$AGEF == 1) g("f2", T) else NULL, omega = if (fl$IMM == 1) g("omega", T - 1) else NULL,
              N1 = g("N1", T), N2 = g("N2", T), N3 = if (fl$IMM == 1) g("N3", T) else NULL)
    lt <- ipm_ltre_draws(d, list(phij = rm[, "mean_phij"], phia = rm[, "mean_phia"], f1 = rm[, "mean_f1"],
                                 f2 = if (fl$AGEF == 1) rm[, "mean_f2"] else NULL,
                                 omega = if (fl$IMM == 1) rm[, "mean_omega"] else NULL),
                         fl, settings$sex_ratio, transient = TRUE)
    if (!is.null(lt$transient)) obj$transient <- apply(lt$transient, 2:3, mean)
    Nt <- g("N1", T) + g("N2", T) + (if (fl$IMM == 1) g("N3", T) else 0)
    lam <- Nt[, -1] / Nt[, -T]
    obj$lambda_diff <- colMeans(lam[, -1, drop = FALSE] - lam[, -(T - 1), drop = FALSE])
  }
  if (settings$ppc_draws > 0) {
    obj$ppc <- tryCatch(ipm_ppc(stack_draws(r$raw), data, fl, settings), error = function(e) NULL)
    if (!is.null(obj$ppc)) notes <- c(notes, sprintf(
      "Posterior predictive checks (Freeman-Tukey discrepancy, %d draws; Bayesian p-value, extreme if < 0.025 or > 0.975): %s.",
      nrow(obj$ppc$obs), paste(sprintf("%s %.2f", names(obj$ppc$p), obj$ppc$p), collapse = ", ")))
  }
  if (isTRUE(settings$compare_separate) && fl$IPM == 1) {
    s2 <- settings; s2$components <- "separate"; s2$compare_separate <- FALSE; s2$ppc_draws <- 0
    sep <- tryCatch(ipm_fit(data, s2), error = function(e) { message("separate fit failed: ", conditionMessage(e)); NULL })
    if (!is.null(sep)) {
      obj$separate <- sep$estimates
      obj$comparison <- ipm_compare(est, sep$estimates)
      if (!isTRUE(sep$diagnostics$converged)) notes <- c(notes, "The separate-analysis comparison fit did not converge (R-hat >= 1.05); its precision numbers are unreliable.")
    }
  }
  dg$notes <- c(notes, ipm_design_notes(data, est, fl))
  list(estimates = est, diagnostics = dg, object = obj)
}

ipm_compare <- function(ipm, sep) {
  p <- intersect(sep$param[!grepl("\\[", sep$param)], ipm$param)
  p <- p[!grepl("_sd$|_beta_", p) | grepl("_beta_", p)]
  a <- ipm[match(p, ipm$param), ]; b <- sep[match(p, sep$param), ]
  data.frame(param = p, ipm_estimate = a$estimate, ipm_sd = a$sd, ipm_lower = a$lower, ipm_upper = a$upper,
             separate_estimate = b$estimate, separate_sd = b$sd, separate_lower = b$lower, separate_upper = b$upper,
             ipm_cv = a$sd / abs(a$estimate), separate_cv = b$sd / abs(b$estimate),
             sd_ratio = a$sd / b$sd, row.names = NULL)
}

ipm_design_notes <- function(data, est, fl) {
  n <- character()
  T <- length(data$years)
  if (T < 10 && fl$RE == 1)
    n <- c(n, sprintf("Only %d years: temporal SDs of the vital rates and the LTRE contributions rest on few years and are imprecise (see VALIDATION.md, short-series design).", T))
  rj <- sum(data$marr_j[, -ncol(data$marr_j)])
  if (rj < 30) n <- c(n, sprintf("Only %d recaptures of animals marked as juveniles: juvenile survival is informed mainly by the counts through the IPM (it confounds survival with permanent emigration of young).", rj))
  n
}

# m-array cell probabilities in R (for predictive checks)
ipm_marr_probs <- function(phi1, phia, p) {
  T <- length(phi1) + 1
  pr <- matrix(0, T - 1, T)
  for (t in 1:(T - 1)) {
    s <- phi1[t]; nr <- 1
    for (j in t:(T - 1)) {
      if (j > t) { s <- s * phia[j]; nr <- nr * (1 - p) }
      pr[t, j] <- s * nr * p
    }
    pr[t, T] <- max(0, 1 - sum(pr[t, 1:(T - 1)]))
  }
  pr
}

ipm_ppc <- function(raw, data, fl, settings) {
  T <- length(data$years)
  idx <- ipm_thin_idx(nrow(raw), settings$ppc_draws)
  ft <- function(o, e) sum((sqrt(o) - sqrt(e))^2)
  g <- function(m, b, n) m[paste0(b, "[", seq_len(n), "]")]
  rel_j <- rowSums(data$marr_j); rel_a <- rowSums(data$marr_a)
  obs <- rep_ <- matrix(NA_real_, length(idx), 3, dimnames = list(NULL, c("counts", "productivity", "capture_recapture")))
  for (k in seq_along(idx)) {
    m <- raw[idx[k], ]
    phij <- g(m, "phij", T - 1); phia <- g(m, "phia", T - 1); p <- m[["p"]]
    f1 <- g(m, "f1", T); f2 <- if (fl$AGEF == 1) g(m, "f2", T) else f1
    # productivity
    e1 <- f1 * data$prod$B1; e2 <- f2 * data$prod$B2
    jo <- c(data$prod$J1, if (fl$AGEF == 1) data$prod$J2); je <- c(e1, if (fl$AGEF == 1) e2)
    jr <- stats::rpois(length(je), je)
    obs[k, 2] <- ft(jo, je); rep_[k, 2] <- ft(jr, je)
    # m-arrays
    pj <- ipm_marr_probs(phij, phia, p); pa <- ipm_marr_probs(phia, phia, p)
    ej <- pj * rel_j; ea <- pa * rel_a
    rj <- t(sapply(seq_len(T - 1), function(t) if (rel_j[t] > 0) stats::rmultinom(1, rel_j[t], pj[t, ]) else rep(0, T)))
    ra <- t(sapply(seq_len(T - 1), function(t) if (rel_a[t] > 0) stats::rmultinom(1, rel_a[t], pa[t, ]) else rep(0, T)))
    obs[k, 3] <- ft(data$marr_j, ej) + ft(data$marr_a, ea); rep_[k, 3] <- ft(rj, ej) + ft(ra, ea)
    if (fl$IPM == 1) {
      Nt <- g(m, "N1", T) + g(m, "N2", T) + (if (fl$IMM == 1) g(m, "N3", T) else 0)
      ok <- !is.na(data$counts)
      if (fl$CNORM == 1) {
        sg <- m[["sigma_count"]]
        cr <- stats::rnorm(sum(ok), Nt[ok], sg)
        obs[k, 1] <- sum((data$counts[ok] - Nt[ok])^2) / sg^2; rep_[k, 1] <- sum((cr - Nt[ok])^2) / sg^2
      } else {
        cr <- stats::rpois(sum(ok), Nt[ok])
        obs[k, 1] <- ft(data$counts[ok], Nt[ok]); rep_[k, 1] <- ft(cr, Nt[ok])
      }
    }
  }
  keep <- colSums(!is.na(obs)) > 0
  obs <- obs[, keep, drop = FALSE]; rep_ <- rep_[, keep, drop = FALSE]
  list(obs = obs, rep = rep_, p = colMeans(rep_ >= obs))
}
