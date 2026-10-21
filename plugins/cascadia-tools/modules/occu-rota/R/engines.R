# occu-rota engines. Each returns a list with
#   engine, coef_names (P + Q link-scale names), point (MLE or NULL),
#   draws (n x (P+Q) matrix of beta then alpha), chains (Bayesian: list of
#   per-chain matrices), diag (converged, notes, ...), fm (raw fit).

# ------------------------------------------------------------------ unmarked
rota_engine_unmarked <- function(des, settings) {
  eco_require("unmarked")
  suppressPackageStartupMessages(requireNamespace("unmarked"))
  S <- des$S; M <- des$M; J <- des$J
  ylist <- stats::setNames(lapply(seq_len(S), function(s) des$y[, , s]), des$species)
  sc <- as.data.frame(des$X); names(sc) <- paste0("x", seq_len(des$P))
  oc <- as.data.frame(apply(des$W, 3, function(m) as.vector(t(m))))
  names(oc) <- paste0("d", seq_len(des$Q))
  umf <- suppressWarnings(unmarked::unmarkedFrameOccuMulti(y = ylist, siteCovs = sc, obsCovs = oc,
                                                           maxOrder = des$max_order))
  fcols <- colnames(umf@fDesign)
  fterm <- sub("^f[0-9]+\\[(.*)\\]$", "\\1", fcols)
  sform <- character(); perm <- integer()
  for (ft in fterm) {
    k <- match(ft, names(des$terms))
    if (is.na(k)) { sform <- c(sform, "0"); next }
    cols <- which(des$term_of == k)
    nonint <- cols[-1]
    sform <- c(sform, if (length(nonint)) paste("~", paste0("x", nonint, collapse = " + ")) else "~ 1")
    perm <- c(perm, cols)
  }
  dform <- vapply(seq_len(S), function(s) {
    cols <- which(des$det_of == s)[-1]
    if (length(cols)) paste("~", paste0("d", cols, collapse = " + ")) else "~ 1"
  }, "")
  dperm <- unlist(lapply(seq_len(S), function(s) which(des$det_of == s)))
  t0 <- proc.time()[["elapsed"]]
  # do.call so the stored call holds values (unmarked's vcov() re-reads call$penalty)
  fm <- tryCatch(suppressWarnings(do.call(unmarked::occuMulti, list(
                   detformulas = dform, stateformulas = sform, data = umf,
                   maxOrder = des$max_order, penalty = settings$penalty,
                   boot = settings$penalty_boot, silent = TRUE, control = list(maxit = 2000)))),
                 error = function(e) e)
  if (inherits(fm, "error")) stop("occuMulti failed: ", conditionMessage(fm))
  cf <- unmarked::coef(fm)
  verr <- NULL
  V <- tryCatch(unmarked::vcov(fm),
                error = function(e) { verr <<- conditionMessage(e); matrix(NA_real_, length(cf), length(cf)) })
  if (length(cf) != des$P + des$Q) stop("internal: occuMulti returned ", length(cf), " coefficients, expected ", des$P + des$Q)
  # reorder unmarked's (state by fDesign order, then detection by species) into ours
  ord <- c(perm, des$P + dperm)
  point <- numeric(des$P + des$Q); point[ord] <- cf
  Vm <- matrix(NA_real_, length(cf), length(cf)); Vm[ord, ord] <- V
  nm <- c(colnames(des$X), rota_det_names(des))
  names(point) <- nm; dimnames(Vm) <- list(nm, nm)
  se <- sqrt(pmax(diag(Vm), 0))
  notes <- if (!is.null(verr)) paste("Variance-covariance matrix unavailable:", verr) else character()
  ok <- fm@opt$convergence == 0 && all(is.finite(se)) && all(se > 0)
  if (fm@opt$convergence != 0) notes <- c(notes, "The optimiser did not report convergence.")
  if (!all(is.finite(se)) || any(se <= 0))
    notes <- c(notes, "Some standard errors could not be computed (singular Hessian): at least one parameter is not identifiable from these data.")
  big <- nm[abs(point) > 8 | se > 10]
  if (length(big)) {
    ok <- FALSE
    notes <- c(notes, paste0("Parameter(s) at or near a boundary (|estimate| > 8 or SE > 10 on the link scale): ",
                             paste(big, collapse = ", "),
                             ". Usually a species pair never or always co-occurs in the detections; consider settings penalty (Clipp et al. 2021) or a simpler interaction formula."))
  }
  draws <- if (all(is.finite(Vm))) {
    set.seed(settings$seed)
    Vs <- (Vm + t(Vm)) / 2
    ev <- eigen(Vs, symmetric = TRUE)
    L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(point))
    d <- matrix(stats::rnorm(settings$n_draws * length(point)), settings$n_draws) %*% t(L)
    sweep(d, 2, point, "+")
  } else matrix(point, 2, length(point), byrow = TRUE)
  colnames(draws) <- nm
  list(engine = "unmarked", point = point, vcov = Vm, se = se, draws = draws, fm = fm,
       diag = list(converged = ok, notes = notes, loglik = -fm@negLogLike, aic = fm@AIC,
                   t_fit = proc.time()[["elapsed"]] - t0))
}

# link-scale names of the detection coefficients (intercepts named by the
# probability they become after plogis)
rota_det_names <- function(des) unlist(lapply(seq_along(des$wlist), function(s) {
  cn <- colnames(des$wlist[[s]])
  ifelse(cn == "(Intercept)", paste0("detection_prob_", des$species[s]),
         paste0("p_", des$species[s], "_beta_", rota_public_col(cn)))
}))

# ------------------------------------------------------------------- nimble
rota_nimble_code <- function() nimble::nimbleCode({
  for (h in 1:P) beta[h] ~ dnorm(0, sd = prior_sd)
  for (h in 1:Q) alpha[h] ~ dnorm(0, sd = prior_sd)
  for (i in 1:M) {
    # natural parameters f[i, k] = X_ik' beta_k (XB masks columns not in term k)
    for (k in 1:K) f[i, k] <- inprod(XB[i, 1:P, k], beta[1:P])
    # multivariate Bernoulli state probabilities over the C = 2^S latent states
    for (c in 1:C) ew[i, c] <- exp(inprod(A[c, 1:K], f[i, 1:K]))
    psi[i, 1:C] <- ew[i, 1:C] / sum(ew[i, 1:C])
    z[i] ~ dcat(psi[i, 1:C])
    for (s in 1:S) {
      for (j in 1:J) {
        logit(p[i, j, s]) <- inprod(WB[i, j, 1:Q, s], alpha[1:Q])
        y[i, j, s] ~ dbern(Zm[z[i], s] * p[i, j, s] * w[i, j])
      }
    }
  }
})

rota_engine_nimble <- function(des, settings) {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))  # nimble model building needs the package attached
  M <- des$M; J <- des$J; S <- des$S; K <- des$K; P <- des$P; Q <- des$Q; C <- des$C
  XB <- array(0, c(M, P, K))
  for (h in seq_len(P)) XB[, h, des$term_of[h]] <- des$X[, h]
  WB <- array(0, c(M, J, Q, S))
  for (q in seq_len(Q)) WB[, , q, des$det_of[q]] <- des$W[, , q]
  y <- des$y; y[is.na(y)] <- 0
  # initial state: species detected -> present; undetected species present with prob 0.5
  det <- apply(y, c(1, 3), max)
  code_state <- function(zz) 1 + as.vector(zz %*% (2^(seq_len(S) - 1)))   # matches rota_states ordering
  stopifnot(all(code_state(des$Z) == seq_len(C)))
  inits <- function() {
    zz <- det
    zz[det == 0] <- stats::rbinom(sum(det == 0), 1, 0.5)
    list(beta = stats::rnorm(P, 0, 0.3), alpha = stats::rnorm(Q, 0, 0.3), z = code_state(zz))
  }
  key <- paste("occu-rota", M, J, S, K, P, Q, C, paste(des$term_of, collapse = "."),
               paste(des$det_of, collapse = "."), paste(as.vector(des$A), collapse = ""),
               settings$prior_sd, settings$sampler, sep = "_")
  conf_fun <- function(conf, m) {
    if (identical(settings$sampler, "block")) {
      conf$removeSamplers(c("beta", "alpha"))
      conf$addSampler(target = "beta", type = "RW_block")
      conf$addSampler(target = "alpha", type = "RW_block")
    } else if (identical(settings$sampler, "slice")) {
      conf$removeSamplers(c("beta", "alpha"))
      for (h in seq_len(P)) conf$addSampler(target = paste0("beta[", h, "]"), type = "slice")
      for (h in seq_len(Q)) conf$addSampler(target = paste0("alpha[", h, "]"), type = "slice")
    }
  }
  r <- nimble_run(rota_nimble_code(),
                  constants = list(M = M, J = J, S = S, K = K, P = P, Q = Q, C = C, prior_sd = settings$prior_sd),
                  data = list(y = y, w = des$w, XB = XB, WB = WB, A = des$A, Zm = des$Z),
                  inits = inits, monitors = c("beta", "alpha"),
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, configure = conf_fun, seed = settings$seed)
  nm <- c(colnames(des$X), rota_det_names(des))
  ch <- lapply(r$samples, function(m) {
    m <- as.matrix(m)
    out <- cbind(m[, paste0("beta[", seq_len(P), "]"), drop = FALSE], m[, paste0("alpha[", seq_len(Q), "]"), drop = FALSE])
    colnames(out) <- nm; out
  })
  list(engine = "nimble", point = NULL, draws = stack_draws(ch), chains = ch,
       diag = list(converged = TRUE, notes = character(), t_setup = r$diagnostics$t_setup,
                   t_fit = r$diagnostics$t_total))
}

# -------------------------------------------------- naive (negative control)
# Treats "detected at least once" as the true occupancy state (detection = 1)
# and fits the same multivariate-Bernoulli state model by maximum likelihood.
rota_engine_naive <- function(des, settings) {
  y <- des$y; y[is.na(y)] <- 0
  det <- apply(y, c(1, 3), max)
  obs_state <- 1 + as.vector(det %*% (2^(seq_len(des$S) - 1)))
  nll <- function(b) {
    Psi <- rota_psi(rota_F_sites(des, b), des$A)
    -sum(log(pmax(Psi[cbind(seq_len(des$M), obs_state)], 1e-300)))
  }
  o <- stats::optim(rep(0, des$P), nll, method = "BFGS", hessian = TRUE, control = list(maxit = 2000))
  V <- tryCatch(solve(o$hessian), error = function(e) matrix(NA_real_, des$P, des$P))
  point <- o$par; names(point) <- colnames(des$X)
  ok <- o$convergence == 0 && all(is.finite(V))
  set.seed(settings$seed)
  draws <- if (all(is.finite(V))) {
    Vs <- (V + t(V)) / 2; ev <- eigen(Vs, symmetric = TRUE)
    L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), des$P)
    sweep(matrix(stats::rnorm(settings$n_draws * des$P), settings$n_draws) %*% t(L), 2, point, "+")
  } else matrix(point, 2, des$P, byrow = TRUE)
  colnames(draws) <- names(point)
  list(engine = "naive", point = point, vcov = V, se = sqrt(pmax(diag(V), 0)), draws = draws,
       diag = list(converged = ok, notes = "NAIVE CONTROL: detection ignored (detected-at-least-once treated as presence).",
                   loglik = -o$value))
}
