# occu-single engines. Each returns a common structure:
#   draws   matrix (n x coefficients, nimble-style names beta_psi[1], ...) for
#           derived quantities: posterior draws (nimble) or draws from the
#           asymptotic normal of the MLE (unmarked)
#   chains  list of per-chain matrices of reported quantities (nimble only)
#   point   named coefficient vector (MLE) or NULL; vcov (unmarked)
#   prop    matrix of draws of proportion of sites occupied (n x T)
#   diag    list(converged, notes, ...)

# ---------------------------------------------------------------- unmarked
occ_umf <- function(des, y = des$y) {
  eco_require("unmarked")
  M <- des$M; J <- des$J; T <- des$T
  cols <- function(X, pre) {  # drop intercept, generic syntactic names
    K <- dim(X)[length(dim(X))]
    if (K < 2) return(NULL)
    stats::setNames(seq(2, K), paste0(pre, "_c", seq(2, K)))
  }
  rhs <- function(cc) if (is.null(cc)) "~ 1" else paste("~", paste(names(cc), collapse = " + "))
  cpsi <- cols(des$Xpsi, "psi"); cp <- cols(des$Xp, "p")
  sc <- if (is.null(cpsi)) NULL else as.data.frame(stats::setNames(lapply(cpsi, function(k) des$Xpsi[, k]), names(cpsi)))
  # obsCovs rows: site-major, then occasions in y-column order (season-major, survey within season)
  oc <- if (is.null(cp)) NULL else as.data.frame(stats::setNames(lapply(cp, function(k)
    as.vector(aperm(des$Xp[, , , k, drop = TRUE] |> array(c(M, J, T)), c(2, 3, 1)))), names(cp)))
  ymat <- matrix(aperm(array(y, c(M, J, T)), c(1, 2, 3)), M, J * T)
  if (des$model == "static") {
    umf <- unmarked::unmarkedFrameOccu(y = ymat, siteCovs = sc, obsCovs = oc)
    return(list(umf = umf, f = list(psi = rhs(cpsi), p = rhs(cp))))
  }
  ccol <- cols(des$Xcol, "col"); cext <- cols(des$Xext, "ext")
  pad <- function(X, k) as.vector(t(cbind(X[, , k], 0)))  # M x (T-1) -> rows site-major, season fastest
  ys <- c(if (!is.null(ccol)) stats::setNames(lapply(ccol, function(k) pad(des$Xcol, k)), names(ccol)),
          if (!is.null(cext)) stats::setNames(lapply(cext, function(k) pad(des$Xext, k)), names(cext)))
  ys <- if (length(ys)) as.data.frame(ys) else NULL
  umf <- unmarked::unmarkedMultFrame(y = ymat, siteCovs = sc, obsCovs = oc, yearlySiteCovs = ys, numPrimary = T)
  list(umf = umf, f = list(psi = rhs(cpsi), col = rhs(ccol), ext = rhs(cext), p = rhs(cp)))
}

occ_fit_unmarked <- function(des, settings, y = des$y, se = TRUE) {
  u <- occ_umf(des, y)
  fm <- suppressWarnings(tryCatch({
    if (des$model == "static")
      unmarked::occu(stats::as.formula(paste(u$f$p, u$f$psi)), data = u$umf, se = se)
    else unmarked::colext(psiformula = stats::as.formula(u$f$psi), gammaformula = stats::as.formula(u$f$col),
                          epsilonformula = stats::as.formula(u$f$ext), pformula = stats::as.formula(u$f$p),
                          data = u$umf, se = se, method = "BFGS")
  }, error = function(e) e))
  if (inherits(fm, "error")) stop("unmarked fit failed: ", conditionMessage(fm))
  # unmarked orders coefficients state/psi, (col, ext), det: same as des$subs with p last
  ord <- if (des$model == "static") c("psi", "p") else c("psi", "col", "ext", "p")
  cf <- unmarked::coef(fm)
  nm <- unlist(lapply(ord, function(s) paste0("beta_", s, "[", seq_along(des$colnames[[s]]), "]")))
  if (length(cf) != length(nm)) stop("unexpected number of unmarked coefficients")
  names(cf) <- nm
  V <- if (se) tryCatch(unmarked::vcov(fm), error = function(e) NULL) else NULL
  if (!is.null(V)) dimnames(V) <- list(nm, nm)
  list(fm = fm, coef = cf, vcov = V)
}

occ_mvn_draws <- function(mu, V, n) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(mu))
  Z <- matrix(stats::rnorm(n * length(mu)), n)
  out <- sweep(Z %*% t(L), 2, mu, "+")
  colnames(out) <- names(mu)
  out
}

occ_engine_unmarked <- function(des, settings) {
  u <- occ_fit_unmarked(des, settings)
  cf <- u$coef; V <- u$vcov
  notes <- character()
  se <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else rep(NA_real_, length(cf))
  opt_ok <- isTRUE(u$fm@opt$convergence == 0)
  if (!opt_ok) notes <- c(notes, "Optimiser did not report convergence.")
  hess_ok <- !is.null(V) && all(is.finite(V)) && all(diag(V) > 0)
  if (!hess_ok) notes <- c(notes, "Hessian not invertible: standard errors unavailable (parameters not identified).")
  boundary <- any(abs(cf) > 8) || any(!is.finite(se)) || any(se > 8)
  if (boundary) notes <- c(notes, sprintf(
    "Some estimates are at or near a boundary (|logit coefficient| up to %.1f, SE up to %.1f). Probabilities near 0 or 1 with huge uncertainty mean the data cannot estimate them; do not interpret.",
    max(abs(cf)), max(se[is.finite(se)], 0)))
  converged <- opt_ok && hess_ok && !boundary
  draws <- if (hess_ok) occ_mvn_draws(cf, V, settings$n_draws) else matrix(cf, 1, dimnames = list(NULL, names(cf)))
  list(engine = "unmarked", point = cf, vcov = V, se = stats::setNames(se, names(cf)), draws = draws,
       fm = u$fm, diag = list(converged = converged, rhat_max = NA_real_, ess_min = NA_real_, notes = notes,
                               loglik = -u$fm@negLogLike))
}

# ------------------------------------------------------------------ nimble
occ_nimble_code <- function(model) {
  if (model == "static") nimble::nimbleCode({
    beta_psi[1] ~ dlogis(0, 1)
    beta_p[1] ~ dlogis(0, 1)
    if (Kpsi > 1) { for (k in 2:Kpsi) { beta_psi[k] ~ dnorm(0, sd = prior_sd) } }
    if (Kp > 1) { for (k in 2:Kp) { beta_p[k] ~ dnorm(0, sd = prior_sd) } }
    for (i in 1:M) {
      logit(psi[i]) <- inprod(beta_psi[1:Kpsi], Xpsi[i, 1:Kpsi])
      z[i] ~ dbern(psi[i])
      for (j in 1:J) {
        logit(p[i, j]) <- inprod(beta_p[1:Kp], Xp[i, j, 1:Kp])
        y[i, j] ~ dbern(z[i] * p[i, j] * w[i, j])
      }
    }
    n_occ <- sum(z[1:M])
  }) else nimble::nimbleCode({
    beta_psi[1] ~ dlogis(0, 1)
    beta_p[1] ~ dlogis(0, 1)
    beta_col[1] ~ dlogis(0, 1)
    beta_ext[1] ~ dlogis(0, 1)
    if (Kpsi > 1) { for (k in 2:Kpsi) { beta_psi[k] ~ dnorm(0, sd = prior_sd) } }
    if (Kp > 1) { for (k in 2:Kp) { beta_p[k] ~ dnorm(0, sd = prior_sd) } }
    if (Kcol > 1) { for (k in 2:Kcol) { beta_col[k] ~ dnorm(0, sd = prior_sd) } }
    if (Kext > 1) { for (k in 2:Kext) { beta_ext[k] ~ dnorm(0, sd = prior_sd) } }
    for (i in 1:M) {
      logit(psi[i]) <- inprod(beta_psi[1:Kpsi], Xpsi[i, 1:Kpsi])
      z[i, 1] ~ dbern(psi[i])
      for (t in 2:T) {
        logit(col[i, t - 1]) <- inprod(beta_col[1:Kcol], Xcol[i, t - 1, 1:Kcol])
        logit(ext[i, t - 1]) <- inprod(beta_ext[1:Kext], Xext[i, t - 1, 1:Kext])
        z[i, t] ~ dbern(z[i, t - 1] * (1 - ext[i, t - 1]) + (1 - z[i, t - 1]) * col[i, t - 1])
      }
      for (t in 1:T) {
        for (j in 1:J) {
          logit(p[i, j, t]) <- inprod(beta_p[1:Kp], Xp[i, j, t, 1:Kp])
          y[i, j, t] ~ dbern(z[i, t] * p[i, j, t] * w[i, j, t])
        }
      }
    }
    for (t in 1:T) { n_occ[t] <- sum(z[1:M, t]) }
  })
}

occ_engine_nimble <- function(des, settings) {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))  # nimble model building needs the package attached
  K <- lapply(des$colnames, length)
  y0 <- des$y; y0[des$w == 0] <- 0
  if (des$model == "static") {
    constants <- list(M = des$M, J = des$J, Kpsi = K$psi, Kp = K$p, prior_sd = settings$prior_sd)
    data <- list(y = matrix(y0, des$M, des$J), w = matrix(des$w, des$M, des$J), Xpsi = des$Xpsi,
                 Xp = array(des$Xp, c(des$M, des$J, K$p)))
    zinit <- function() rep(1, des$M)
  } else {
    constants <- list(M = des$M, J = des$J, T = des$T, Kpsi = K$psi, Kp = K$p, Kcol = K$col, Kext = K$ext,
                      prior_sd = settings$prior_sd)
    data <- list(y = y0, w = des$w, Xpsi = des$Xpsi, Xp = des$Xp, Xcol = des$Xcol, Xext = des$Xext)
    zinit <- function() matrix(1, des$M, des$T)
  }
  dimnames(data$Xpsi) <- NULL; dimnames(data$Xp) <- NULL
  if (!is.null(data$Xcol)) { dimnames(data$Xcol) <- NULL; dimnames(data$Xext) <- NULL }
  inits <- function() {
    l <- list(z = zinit())
    for (s in des$subs) l[[paste0("beta_", s)]] <- c(stats::rnorm(1, 0, 0.5), rep(0, K[[s]] - 1))
    l
  }
  mons <- c(paste0("beta_", des$subs), "n_occ")
  key <- paste("occu-single", des$model, des$M, des$J, des$T, paste(unlist(K), collapse = "-"),
               settings$prior_sd, sep = "_")
  r <- nimble_run(occ_nimble_code(des$model), constants, data, inits, mons,
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    if (isTRUE(settings$slice)) for (s in des$subs) {
                      nd <- m$expandNodeNames(paste0("beta_", s))
                      conf$removeSamplers(nd); for (n in nd) conf$addSampler(n, "slice")
                    }
                  })
  draws <- stack_draws(r$samples)
  prop <- draws[, grep("^n_occ", colnames(draws)), drop = FALSE] / des$M
  list(engine = "nimble", point = NULL, draws = draws[, occ_coef_internal(des), drop = FALSE],
       chains = r$samples, prop = prop,
       diag = list(converged = r$diagnostics$converged, notes = character(),
                   t_total = r$diagnostics$t_total))
}

# ------------------------------------------- naive control (ignores detection)
# Logistic regression of "ever detected" on the occupancy covariates. Exists
# only as a negative control: it treats non-detection as absence, so
# occupancy is biased low whenever p* = 1 - (1 - p)^J < 1.
occ_engine_naive <- function(des, settings) {
  y0 <- des$y; y0[des$w == 0] <- 0
  ever <- apply(y0, c(1, 3), max)  # M x T
  fitglm <- function(yv, X) {
    g <- suppressWarnings(stats::glm.fit(X, yv, family = stats::binomial()))
    V <- tryCatch(solve(crossprod(X * sqrt(g$weights))), error = function(e) matrix(NA, ncol(X), ncol(X)))
    list(coef = g$coefficients, se = sqrt(diag(V)), conv = g$converged)
  }
  res <- list(psi = fitglm(ever[, 1], des$Xpsi))
  if (des$model == "dynamic") {
    T1 <- des$T - 1
    prev <- as.vector(ever[, 1:T1]); nxt <- as.vector(ever[, 2:des$T])
    Xc <- matrix(des$Xcol, ncol = dim(des$Xcol)[3]); Xe <- matrix(des$Xext, ncol = dim(des$Xext)[3])
    res$col <- fitglm(nxt[prev == 0], Xc[prev == 0, , drop = FALSE])
    res$ext <- fitglm(1 - nxt[prev == 1], Xe[prev == 1, , drop = FALSE])
  }
  list(engine = "naive", glm = res, ever = ever,
       diag = list(converged = all(vapply(res, function(r) isTRUE(r$conv), TRUE)), rhat_max = NA_real_,
                   notes = "NAIVE CONTROL: detection ignored; occupancy is the proportion of sites where the species was seen (biased low)."))
}

# ------------------------------------------------- finite-sample occupancy
# Draws of the proportion of the surveyed sites occupied in each season:
# for each coefficient draw, P(z = 1 | detections) from forward-backward,
# then z ~ Bernoulli. Sites are independent given the parameters, so drawing
# each site's z from its marginal gives the exact distribution of the
# per-season total.
occ_prop_draws <- function(des, draws) {
  n <- nrow(draws)
  out <- matrix(NA_real_, n, des$T)
  for (d in seq_len(n)) {
    pr <- occ_probs(des, occ_split_coef(draws[d, ], des))
    post <- occ_forward_backward(des, pr)$post
    z <- matrix(stats::rbinom(length(post), 1, post), des$M)
    out[d, ] <- colMeans(z)
  }
  out
}
