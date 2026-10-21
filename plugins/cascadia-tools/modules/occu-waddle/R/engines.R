# occu-waddle engines. Each returns list(engine, point (MLE or NULL),
# draws (n x P coefficient draws, columns = des$coef_names), chains (Bayesian),
# diag (converged, rhat_max, ess_min, notes, loglik)).

# ------------------------------------------------------------------ Nimble
# Marginal likelihood (latent states summed out) as a custom distribution.
# x: (M*J) x 2 matrix of detections (col 1 dominant, col 2 subordinate),
# row r = (i-1)*J + j; obs[r] = 1 if survey j at site i was done.
wd_nimble_dist <- function() {
  if (exists("dwaddle", envir = .GlobalEnv, inherits = FALSE)) return(invisible(TRUE))
  dwaddle <- nimble::nimbleFunction(
    run = function(x = double(2), psiA = double(1), psiBA = double(1), psiBa = double(1),
                   pA = double(1), pB = double(1), rA = double(1), rBA = double(1), rBa = double(1),
                   obs = double(1), J = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      M <- length(psiA)
      ll <- 0
      for (i in 1:M) {
        l1 <- 0; l2 <- 0; l3 <- 0; nA <- 0; nB <- 0
        for (j in 1:J) {
          r <- (i - 1) * J + j
          if (obs[r] > 0) {
            ya <- x[r, 1]; yb <- x[r, 2]
            nA <- nA + ya; nB <- nB + yb
            if (ya > 0) { l1 <- l1 + log(rA[r]); q <- rBA[r]; l2 <- l2 + log(pA[r]) }
            else { l1 <- l1 + log(1 - rA[r]); q <- rBa[r]; l2 <- l2 + log(1 - pA[r]) }
            if (yb > 0) { l1 <- l1 + log(q); l3 <- l3 + log(pB[r]) }
            else { l1 <- l1 + log(1 - q); l3 <- l3 + log(1 - pB[r]) }
          }
        }
        lik <- psiA[i] * psiBA[i] * exp(l1)
        if (nB == 0) lik <- lik + psiA[i] * (1 - psiBA[i]) * exp(l2)
        if (nA == 0) lik <- lik + (1 - psiA[i]) * psiBa[i] * exp(l3)
        if (nA + nB == 0) lik <- lik + (1 - psiA[i]) * (1 - psiBa[i])
        ll <- ll + log(lik)
      }
      if (log) return(ll)
      return(exp(ll))
    })
  rwaddle <- nimble::nimbleFunction(
    run = function(n = integer(0), psiA = double(1), psiBA = double(1), psiBa = double(1),
                   pA = double(1), pB = double(1), rA = double(1), rBA = double(1), rBa = double(1),
                   obs = double(1), J = double(0)) {
      returnType(double(2))
      out <- matrix(0, nrow = length(pA), ncol = 2)
      return(out)
    })
  assign("dwaddle", dwaddle, envir = .GlobalEnv)
  assign("rwaddle", rwaddle, envir = .GlobalEnv)
  suppressMessages(nimble::registerDistributions(list(dwaddle = list(
    BUGSdist = "dwaddle(psiA, psiBA, psiBa, pA, pB, rA, rBA, rBa, obs, J)",
    types = c("value = double(2)", "psiA = double(1)", "psiBA = double(1)", "psiBa = double(1)",
              "pA = double(1)", "pB = double(1)", "rA = double(1)", "rBA = double(1)", "rBa = double(1)",
              "obs = double(1)", "J = double(0)"),
    discrete = TRUE, mixedSizes = TRUE)), verbose = FALSE))
  invisible(TRUE)
}

# Model code generated for the free parameters (ties pass the same vector).
wd_nimble_code <- function(des, prior_sd) {
  L <- character()
  for (p in des$free) {
    N <- if (p %in% wd_occ_pars) "M" else "NJ"
    K <- des$K[[p]]
    L <- c(L, sprintf("b_%s[1] ~ dlogis(0, 1)", p))
    if (K > 1) L <- c(L, sprintf("for (k in 2:%d) { b_%s[k] ~ dnorm(0, sd = %s) }", K, p, format(prior_sd)))
    L <- c(L, if (K > 1) sprintf("lp_%s[1:%s] <- (X_%s[1:%s, 1:%d] %%*%% b_%s[1:%d])[1:%s, 1]", p, N, p, N, K, p, K, N)
              else sprintf("lp_%s[1:%s] <- X_%s[1:%s, 1] * b_%s[1]", p, N, p, N, p),
           sprintf("q_%s[1:%s] <- ilogit(lp_%s[1:%s])", p, N, p, N))
  }
  arg <- function(p) { f <- des$map[[p]]; sprintf("q_%s[1:%s]", f, if (p %in% wd_occ_pars) "M" else "NJ") }
  L <- c(L, sprintf("y[1:NJ, 1:2] ~ dwaddle(%s, obs[1:NJ], J)", paste(vapply(wd_pars, arg, ""), collapse = ", ")))
  eval(parse(text = paste0("nimble::nimbleCode({\n", paste(L, collapse = "\n"), "\n})")))
}

wd_engine_nimble <- function(des, settings) {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))  # model building needs nimble attached
  wd_nimble_dist()
  code <- wd_nimble_code(des, settings$prior_sd)
  NJ <- des$M * des$J
  data <- list(y = cbind(des$ya, des$yb), obs = as.numeric(des$obs_v))
  for (p in des$free) data[[paste0("X_", p)]] <- des$mms[[p]]$X
  constants <- list(M = des$M, NJ = NJ, J = des$J)
  # inits: start from the maximum-likelihood fit when available (jittered), else near 0
  st <- tryCatch(wd_optim(des)$par, error = function(e) rep(0, length(des$coef_names)))
  if (any(!is.finite(st)) || any(abs(st) > 8)) st <- rep(0, length(des$coef_names))
  inits <- function() {
    v <- st + stats::rnorm(length(st), 0, 0.3)
    stats::setNames(lapply(des$free, function(p) v[des$idx[[p]]]), paste0("b_", des$free))
  }
  key <- paste("occu-waddle", des$M, des$J, paste(des$map, collapse = ","),
               paste(des$K, collapse = ","), settings$prior_sd, settings$sampler, sep = "|")
  sampler <- settings$sampler
  conf_fun <- function(conf, m) {
    if (sampler == "af_slice") {
      conf$removeSamplers(); conf$addSampler(target = paste0("b_", des$free), type = "AF_slice")
    } else if (sampler == "slice") {
      conf$removeSamplers()
      for (nd in m$expandNodeNames(paste0("b_", des$free))) conf$addSampler(target = nd, type = "slice")
    }
  }
  r <- nimble_run(code, constants, data, inits, monitors = paste0("b_", des$free),
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, configure = conf_fun, seed = settings$seed)
  # rename columns b_psiA[1] -> descriptive names
  ren <- function(m) {
    m <- as.matrix(m)
    nm <- unlist(lapply(des$free, function(p) sprintf("b_%s[%d]", p, seq_len(des$K[[p]]))))
    m <- m[, nm, drop = FALSE]; colnames(m) <- des$coef_names; m
  }
  chains <- lapply(r$samples, ren)
  list(engine = "nimble", point = NULL, draws = stack_draws(chains), chains = chains,
       diag = list(converged = r$diagnostics$converged, rhat_max = r$diagnostics$rhat_max,
                   ess_min = r$diagnostics$ess_min, t_total = r$diagnostics$t_total, notes = character()))
}

# ------------------------------------------------------- maximum likelihood
wd_optim <- function(des, starts = 3) {
  P <- length(des$coef_names)
  nll <- function(b) { v <- -wd_loglik(des, b); if (!is.finite(v)) 1e10 else v }
  best <- NULL
  for (s in seq_len(starts)) {
    b0 <- if (s == 1) rep(0, P) else stats::rnorm(P, 0, 0.5)
    o <- tryCatch(stats::optim(b0, nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value - 1e-8)) best <- o
  }
  if (is.null(best)) stop("optimisation failed from every start")
  best$hessian <- stats::optimHess(best$par, nll)
  best
}

wd_engine_mle <- function(des, settings) {
  o <- wd_optim(des, settings$mle_starts)
  P <- length(o$par)
  V <- tryCatch(solve(o$hessian), error = function(e) NULL)
  notes <- character()
  ok <- o$convergence == 0 && !is.null(V) && all(is.finite(diag(V))) && all(diag(V) > 0)
  if (!ok) notes <- c(notes, "Maximum-likelihood fit did not converge or the Hessian is singular: some parameters are not identifiable from these data.")
  if (any(abs(o$par) > 7)) {
    ok <- FALSE
    notes <- c(notes, sprintf("Estimate(s) on a boundary (probability near 0 or 1): %s. Intervals for these are not reliable.",
                              paste(des$coef_names[abs(o$par) > 7], collapse = ", ")))
  }
  se <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else rep(NA_real_, P)
  draws <- if (!is.null(V) && all(is.finite(V))) {
    R <- tryCatch(chol((V + t(V)) / 2), error = function(e) NULL)
    if (is.null(R)) matrix(o$par, settings$n_draws, P, byrow = TRUE)
    else sweep(matrix(stats::rnorm(settings$n_draws * P), ncol = P) %*% R, 2, o$par, "+")
  } else matrix(o$par, settings$n_draws, P, byrow = TRUE)
  colnames(draws) <- des$coef_names
  list(engine = "mle", point = stats::setNames(o$par, des$coef_names), se = stats::setNames(se, des$coef_names),
       vcov = V, draws = draws, chains = NULL,
       diag = list(converged = ok, notes = notes, loglik = -o$value, aic = 2 * o$value + 2 * P))
}

# --------------------------------------------- naive (negative control only)
# Ignores imperfect detection: a site counts as occupied if the species was
# ever detected. Occupancy parameters only, by logistic regression.
wd_engine_naive <- function(des, settings) {
  eA <- as.numeric(rowSums(des$yA, na.rm = TRUE) > 0); eB <- as.numeric(rowSums(des$yB, na.rm = TRUE) > 0)
  par <- rep(NA_real_, length(des$coef_names)); se <- par
  fitp <- function(p, y, rows) {
    X <- des$mms[[p]]$X[rows, , drop = FALSE]
    g <- suppressWarnings(stats::glm.fit(X, y[rows], family = stats::binomial()))
    V <- tryCatch(solve(crossprod(X * sqrt(g$weights))), error = function(e) matrix(NA, ncol(X), ncol(X)))
    par[des$idx[[p]]] <<- g$coefficients; se[des$idx[[p]]] <<- sqrt(diag(V))
  }
  fitp("psiA", eA, rep(TRUE, des$M))
  if ("psiBA" %in% des$free) { fitp("psiBA", eB, eA == 1); fitp("psiBa", eB, eA == 0) }
  else fitp("psiBa", eB, rep(TRUE, des$M))
  occ <- unlist(des$idx[intersect(des$free, wd_occ_pars)])
  P <- length(occ)
  D <- matrix(stats::rnorm(settings$n_draws * P), ncol = P) %*% diag(se[occ], P)
  D <- sweep(D, 2, par[occ], "+")
  draws <- matrix(0, settings$n_draws, length(par)); draws[, occ] <- D; colnames(draws) <- des$coef_names
  list(engine = "naive", point = stats::setNames(par, des$coef_names), se = stats::setNames(se, des$coef_names),
       draws = draws, chains = NULL, diag = list(converged = TRUE, notes = "Naive model ignoring detection (negative control)."))
}
