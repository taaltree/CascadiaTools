# abund-interactions engines. Each returns
#   chains   list of per-chain matrices of internal coefficients (nimble), or NULL
#   draws    matrix of coefficient draws (posterior, or MVN around the MLE)
#   N, z     matrices (draws x sites) of latent dominant abundance and
#            subordinate occupancy (posterior or conditional draws)
#   point, vcov  (likelihood engines)
#   diag     list(converged, rhat_max, ess_min, notes)

# ------------------------------------------------------------------ Nimble
# Twining et al. (2025, Ecology 106:e4468) abundance-occupancy model:
#   N_i ~ Poisson(lambda_i),           log(lambda_i) = X_i beta_lam
#   yA_ij ~ Binomial(K_ij, 1 - (1 - r_ij)^N_i)        (Royle & Nichols 2003)
#   z_i ~ Bernoulli(psi_i), logit(psi_i) = X_i beta_psi + N_i (gamma0 + Mod_i gamma_mod)
#   yS_ij ~ Binomial(K_ij, z_i p_ij)
abi_nimble_code <- function(P) {
  # Linear predictors are written out as text: inprod() over a length-1
  # index range (intercept-only formulas) does not compile in nimble 1.4.
  lin <- function(b, X, n, idx) if (n == 1) sprintf("%s[1] * %s[%s, 1]", b, X, idx)
                                else sprintf("inprod(%s[1:%d], %s[%s, 1:%d])", b, n, X, idx, n)
  pri <- function(b, n) if (n > 1) sprintf("for (k in 2:%d) { %s[k] ~ dnorm(0, sd = prior_sd) }", n, b) else ""
  slope <- if (P$mod == 0) "gamma0" else sprintf("(gamma0 + %s)", lin("gamma_mod", "Mod", P$mod, "i"))
  noint <- isTRUE(P$noint)   # reduced model: no interaction (subordinate = single-species occupancy)
  txt <- paste(
    "{",
    "beta_lam[1] ~ dnorm(0, sd = prior_sd)", "beta_r[1] ~ dlogis(0, 1)", "beta_psi[1] ~ dlogis(0, 1)", "beta_p[1] ~ dlogis(0, 1)",
    pri("beta_lam", P$lambda), pri("beta_r", P$r), pri("beta_psi", P$psi), pri("beta_p", P$p),
    if (noint) "" else "gamma0 ~ dnorm(0, sd = prior_sd)",
    if (P$mod > 0) sprintf("for (k in 1:%d) { gamma_mod[k] ~ dnorm(0, sd = prior_sd) }", P$mod) else "",
    "for (i in 1:M) {",
    sprintf("  log(lambda[i]) <- %s", lin("beta_lam", "Xl", P$lambda, "i")),
    "  N[i] ~ dpois(lambda[i])",
    if (noint) sprintf("  logit(psi[i]) <- %s", lin("beta_psi", "Xpsi", P$psi, "i"))
    else sprintf("  logit(psi[i]) <- %s + N[i] * %s", lin("beta_psi", "Xpsi", P$psi, "i"), slope),
    "  z[i] ~ dbern(psi[i])",
    "  for (j in 1:J) {",
    sprintf("    logit(r[i, j]) <- %s", lin("beta_r", "Xr", P$r, "i, j")),
    sprintf("    logit(p[i, j]) <- %s", lin("beta_p", "Xp", P$p, "i, j")),
    "    yA[i, j] ~ dbin(1 - pow(1 - r[i, j], N[i]), K[i, j])",
    "    yS[i, j] ~ dbin(z[i] * p[i, j], K[i, j])",
    "  }",
    "  occA[i] <- step(N[i] - 1)",
    "}",
    "N_total <- sum(N[1:M])", "n_occ_dominant <- sum(occA[1:M])", "n_occ_subordinate <- sum(z[1:M])",
    "}", sep = "\n")
  parse(text = txt)[[1]]
}

abi_engine_nimble <- function(des, settings) {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))
  P <- lapply(des$colnames, length)
  noint <- identical(settings$interaction, FALSE)
  P$noint <- noint
  has_mod <- P$mod > 0
  constants <- list(M = des$M, J = des$J, Pl = P$lambda, Pr = P$r, Ppsi = P$psi, Pp = P$p,
                    prior_sd = settings$prior_sd)
  if (has_mod) constants$Pm <- P$mod
  strip <- function(x) { dimnames(x) <- NULL; x }
  data <- list(yA = strip(des$yA), yS = strip(des$yS), K = strip(des$K), Xl = strip(des$Xl),
               Xr = strip(des$Xr), Xpsi = strip(des$Xpsi), Xp = strip(des$Xp))
  if (has_mod) data$Mod <- strip(des$Mod)
  detA <- rowSums(des$yA) > 0; detS <- rowSums(des$yS) > 0
  inits <- function() {
    l <- list(beta_lam = c(stats::rnorm(1, 0, 0.3), rep(0, P$lambda - 1)),
              beta_r = c(stats::rnorm(1, -1, 0.3), rep(0, P$r - 1)),
              beta_psi = c(stats::rnorm(1, 0, 0.3), rep(0, P$psi - 1)),
              beta_p = c(stats::rnorm(1, -1, 0.3), rep(0, P$p - 1)),
              N = ifelse(detA, 2L, 1L), z = rep(1L, des$M))
    if (!noint) l$gamma0 <- 0
    if (has_mod) l$gamma_mod <- rep(0, P$mod)
    l
  }
  coef_nodes <- c("beta_lam", "beta_r", "beta_psi", "beta_p", if (!noint) "gamma0", if (has_mod) "gamma_mod")
  mons <- c(coef_nodes, "N_total", "n_occ_dominant", "n_occ_subordinate", "N", "z")
  key <- paste("abund-interactions", des$M, des$J, P$lambda, P$r, P$psi, P$p, P$mod, settings$prior_sd,
               settings$block, noint, sep = "_")
  r <- nimble_run(abi_nimble_code(P), constants, data, inits, mons,
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed,
                  summary_params = c("N_total"),
                  configure = function(conf, m) {
                    nd <- m$expandNodeNames(coef_nodes)
                    conf$removeSamplers(nd)
                    if (isTRUE(settings$block)) {
                      # the abundance-detection trade-off (lambda vs r) is the slow direction in RN models
                      conf$addSampler(m$expandNodeNames(c("beta_lam", "beta_r")), "AF_slice")
                      conf$addSampler(m$expandNodeNames(c("beta_psi", if (!noint) "gamma0", if (has_mod) "gamma_mod")), "AF_slice")
                      for (n in m$expandNodeNames("beta_p")) conf$addSampler(n, "slice")
                    } else for (n in nd) conf$addSampler(n, "slice")
                  })
  int <- abi_internal(des)
  if (noint) r$samples <- lapply(r$samples, function(s) cbind(s, gamma0 = 0))
  chains <- lapply(r$samples, function(s) s[, c(int, "N_total", "n_occ_dominant", "n_occ_subordinate"), drop = FALSE])
  all <- stack_draws(r$samples)
  Ncols <- paste0("N[", seq_len(des$M), "]"); zcols <- paste0("z[", seq_len(des$M), "]")
  list(engine = "nimble", chains = chains, draws = all[, int, drop = FALSE],
       N = all[, Ncols, drop = FALSE], z = all[, zcols, drop = FALSE],
       diag = list(notes = character(), t_total = r$diagnostics$t_total))
}

# -------------------------------------------- marginal likelihood (ML) engine
# Independent route: N_i (0..Nmax) and z_i are summed out exactly and the
# likelihood is maximised with optim(). mediator = "abundance" is the same
# model as the Nimble engine; mediator = "occupancy" makes subordinate
# occupancy depend only on whether the dominant is present (1[N_i > 0]),
# the Waddle et al. (2010)-style occupancy-mediated model that Twining et
# al. (2025) show is biased when interactions act through abundance.
abi_site_loglik <- function(des, th, Nmax, mediator = "abundance") {
  lp <- abi_linpred(des, th)
  Ns <- 0:Nmax
  M <- des$M
  lPN <- outer(lp$lam, Ns, function(l, n) stats::dpois(n, l, log = TRUE))
  lA <- matrix(0, M, length(Ns))
  for (j in seq_len(des$J)) {
    Kj <- des$K[, j]; yj <- des$yA[, j]; rj <- lp$r[, j]
    q <- outer(log1p(-rj), Ns)                 # log (1 - r)^N
    pA <- -expm1(q)                            # 1 - (1 - r)^N
    ll <- yj * log(pA) + (Kj - yj) * q
    ll[yj == 0 & matrix(Ns, M, length(Ns), byrow = TRUE) == 0] <- 0
    ll[Kj == 0, ] <- 0
    lA <- lA + ll
  }
  lS1 <- rowSums(des$yS * log(lp$p) + (des$K - des$yS) * log1p(-lp$p))
  any_S <- rowSums(des$yS) > 0
  f <- if (mediator == "abundance") Ns else as.numeric(Ns > 0)
  eta <- lp$eta + outer(lp$slope, f)
  psi <- stats::plogis(eta)
  a <- log(psi) + lS1                          # occupied
  b <- log1p(-psi)                             # unoccupied: only possible if never detected
  b[any_S, ] <- -Inf
  lsub <- pmax(a, b) + log1p(exp(-abs(a - b)))
  lj <- lPN + lA + lsub
  list(lj = lj, lsub_occ = log(psi) + lS1, lsub = lsub, Ns = Ns)
}

abi_negll <- function(theta, des, Nmax, mediator) {
  names(theta) <- abi_internal(des)
  lj <- abi_site_loglik(des, theta, Nmax, mediator)$lj
  mx <- apply(lj, 1, max)
  v <- -sum(mx + log(rowSums(exp(lj - mx))))
  if (!is.finite(v)) 1e10 else v
}

abi_engine_ml <- function(des, settings) {
  int <- abi_internal(des)
  Nmax <- settings$Nmax
  med <- if (is.null(settings$mediator)) "abundance" else settings$mediator
  st <- stats::setNames(rep(0, length(int)), int)
  st["beta_r[1]"] <- -1; st["beta_p[1]"] <- -1
  o <- stats::optim(st, abi_negll, des = des, Nmax = Nmax, mediator = med, method = "BFGS",
                    control = list(maxit = 1000, reltol = 1e-10))
  o2 <- stats::optim(o$par, abi_negll, des = des, Nmax = Nmax, mediator = med, method = "Nelder-Mead",
                     control = list(maxit = 4000, reltol = 1e-12))
  if (o2$value < o$value - 1e-6) o <- stats::optim(o2$par, abi_negll, des = des, Nmax = Nmax, mediator = med,
                                                   method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))
  H <- tryCatch(stats::optimHess(o$par, abi_negll, des = des, Nmax = Nmax, mediator = med), error = function(e) NULL)
  V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
  notes <- character()
  if (o$convergence != 0) notes <- c(notes, "Optimiser did not report convergence.")
  hess_ok <- !is.null(V) && all(is.finite(V)) && all(diag(V) > 0)
  if (!hess_ok) notes <- c(notes, "Hessian not invertible: parameters not identified at this design.")
  se <- if (hess_ok) sqrt(diag(V)) else rep(NA_real_, length(int))
  boundary <- any(abs(o$par) > 10) || any(!is.finite(se)) || any(se > 10)
  if (boundary) notes <- c(notes, "Some estimates are at or near a boundary (|coefficient| or SE > 10): do not interpret.")
  names(o$par) <- int
  # truncation check: posterior mass of N at Nmax
  sl <- abi_site_loglik(des, o$par, Nmax, med)
  post <- exp(sl$lj - apply(sl$lj, 1, max)); post <- post / rowSums(post)
  if (max(post[, Nmax + 1]) > 1e-4) notes <- c(notes, sprintf("Abundance truncation Nmax = %d may be too low (posterior mass %.1e at Nmax): raise settings$Nmax.", Nmax, max(post[, Nmax + 1])))
  # draws: MVN for coefficients, then latent N and z from their conditional distribution
  nd <- settings$n_draws
  draws <- if (hess_ok) abi_mvn(o$par, V, nd) else matrix(o$par, 1, dimnames = list(NULL, int))
  latent <- abi_latent_draws(des, draws, Nmax, med)
  if (nrow(draws) == 1) draws <- draws[c(1, 1), , drop = FALSE]
  list(engine = "ml", mediator = med, point = o$par, vcov = V, se = stats::setNames(se, int), draws = draws,
       N = latent$N, z = latent$z, post_N = post,
       diag = list(converged = o$convergence == 0 && hess_ok && !boundary, rhat_max = NA_real_, ess_min = NA_real_,
                   notes = notes, loglik = -o$value))
}

abi_mvn <- function(mu, V, n) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(mu))
  out <- sweep(matrix(stats::rnorm(n * length(mu)), n) %*% t(L), 2, mu, "+")
  colnames(out) <- names(mu); out
}

# Sample N_i and z_i from p(N_i, z_i | y_i, theta) for each coefficient draw.
abi_latent_draws <- function(des, draws, Nmax, mediator = "abundance") {
  n <- nrow(draws); M <- des$M
  Nd <- matrix(0L, n, M); zd <- matrix(0L, n, M)
  for (d in seq_len(n)) {
    th <- draws[d, ]
    sl <- abi_site_loglik(des, th, Nmax, mediator)
    pr <- exp(sl$lj - apply(sl$lj, 1, max)); pr <- pr / rowSums(pr)
    u <- stats::runif(M)
    Ni <- rowSums(pr[, -ncol(pr), drop = FALSE] |> apply(1, cumsum) |> t() < u)
    Nd[d, ] <- Ni
    a <- sl$lsub_occ[cbind(seq_len(M), Ni + 1)]; b <- sl$lsub[cbind(seq_len(M), Ni + 1)]
    zd[d, ] <- stats::rbinom(M, 1, pmin(1, exp(a - b)))
  }
  list(N = Nd, z = zd)
}

# ------------------------------------------------- two-step plug-in engine
# Step 1: Royle-Nichols model for the dominant alone (unmarked::occuRN on
# per-occasion detection/non-detection), and the best unbiased predictor
# E[N_i | y_i] (unmarked::ranef). Step 2: single-species occupancy for the
# subordinate (unmarked::occu) with E[N_i | y] as a fixed, error-free
# covariate. This is the common "plug-in" shortcut; it ignores uncertainty
# in N_i (intervals too narrow) and suffers regression dilution.
abi_engine_twostep <- function(des, settings) {
  eco_require("unmarked")
  bin <- function(y) { y <- (y > 0) * 1; y[des$K == 0] <- NA; y }
  P <- lapply(des$colnames, length)
  oc <- list()
  add_obs <- function(X, pre) {
    if (dim(X)[3] < 2) return(NULL)
    nm <- paste0(pre, "_c", 2:dim(X)[3])
    for (k in 2:dim(X)[3]) oc[[nm[k - 1]]] <<- as.vector(t(X[, , k]))
    nm
  }
  rn <- add_obs(des$Xr, "r"); pn <- add_obs(des$Xp, "p")
  Kvar <- length(unique(des$K[des$K > 0])) > 1
  if (Kvar) { oc$logK <- as.vector(t(log(pmax(des$K, 1)))); rn <- c(rn, "logK"); pn <- c(pn, "logK") }
  ocdf <- if (length(oc)) as.data.frame(oc) else NULL
  sc <- data.frame(des$Xl[, -1, drop = FALSE], check.names = FALSE)
  names(sc) <- if (ncol(sc)) paste0("l_c", seq_len(ncol(sc))) else character()
  rhs <- function(v) if (length(v)) paste("~", paste(v, collapse = " + ")) else "~ 1"
  umfA <- unmarked::unmarkedFrameOccu(y = bin(des$yA), siteCovs = if (ncol(sc)) sc else NULL, obsCovs = ocdf)
  fA <- stats::as.formula(paste(rhs(rn), rhs(names(sc))))
  fmA <- suppressWarnings(unmarked::occuRN(fA, data = umfA, K = settings$Nmax))
  bup <- unmarked::bup(unmarked::ranef(fmA, K = settings$Nmax), stat = "mean")
  # step 2
  scS <- data.frame(des$Xpsi[, -1, drop = FALSE], check.names = FALSE)
  names(scS) <- if (ncol(scS)) paste0("s_c", seq_len(ncol(scS))) else character()
  scS$Nhat <- bup
  mterms <- character()
  if (P$mod > 0) for (k in seq_len(P$mod)) { nm <- paste0("NxM", k); scS[[nm]] <- bup * des$Mod[, k]; mterms <- c(mterms, nm) }
  umfS <- unmarked::unmarkedFrameOccu(y = bin(des$yS), siteCovs = scS, obsCovs = ocdf)
  fS <- stats::as.formula(paste(rhs(pn), rhs(c(setdiff(names(scS), c("Nhat", mterms)), "Nhat", mterms))))
  fmS <- suppressWarnings(unmarked::occu(fS, data = umfS))
  cA <- unmarked::coef(fmA); VA <- unmarked::vcov(fmA)
  cS <- unmarked::coef(fmS); VS <- unmarked::vcov(fmS)
  # map to internal names: occuRN order = state (lambda), det (r [+ logK]); occu = state (psi..., Nhat, NxM), det (p [+logK])
  nl <- P$lambda; nr <- length(cA) - nl
  npsi <- P$psi; ns <- length(cS)
  idx_r <- nl + seq_len(P$r)
  idx_p <- npsi + 1 + P$mod + seq_len(P$p)
  th <- c(cA[seq_len(nl)], cA[idx_r], cS[seq_len(npsi)], cS[npsi + 1],
          if (P$mod > 0) cS[npsi + 1 + seq_len(P$mod)], cS[idx_p])
  names(th) <- abi_internal(des)
  se <- c(sqrt(diag(VA))[c(seq_len(nl), idx_r)], sqrt(diag(VS))[c(seq_len(npsi + 1 + P$mod), idx_p)])
  names(se) <- names(th)
  Kc <- if (Kvar) NA else unique(des$K[des$K > 0])
  notes <- c("TWO-STEP PLUG-IN (comparison only): dominant abundance from occuRN, then E[N|y] used as a fixed covariate for the subordinate; ignores uncertainty in N.",
             if (!Kvar && Kc > 1) sprintf("Detection probabilities are per occasion of %d sub-occasions (not per day).", Kc))
  iA <- c(seq_len(nl), idx_r); iS <- c(seq_len(npsi + 1 + P$mod), idx_p)
  V <- matrix(0, length(th), length(th), dimnames = list(names(th), names(th)))
  V[seq_along(iA), seq_along(iA)] <- VA[iA, iA]
  V[length(iA) + seq_along(iS), length(iA) + seq_along(iS)] <- VS[iS, iS]
  ok <- all(is.finite(se)) && max(abs(th)) < 10
  list(engine = "twostep", point = th, se = se, vcov = V, bup = bup, fmA = fmA, fmS = fmS,
       draws = abi_mvn(th, V, settings$n_draws), N = NULL, z = NULL,
       diag = list(converged = ok, rhat_max = NA_real_, ess_min = NA_real_, notes = notes))
}
