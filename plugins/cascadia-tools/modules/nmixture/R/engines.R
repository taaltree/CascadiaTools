# nmixture engines. Each returns a common structure:
#   list(engine, theta (named, link scale, internal order), vcov, loglik, aic,
#        K, converged, notes, fm (raw fit), draws (likelihood engines: MVN draws of theta)
#        chains (nimble: list of per-chain theta matrices))

nm_mvn_draws <- function(theta, V, n) {
  ok <- all(is.finite(V)) && all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > -1e-10)
  if (!ok) return(matrix(theta, n, length(theta), byrow = TRUE, dimnames = list(NULL, names(theta))))
  ev <- eigen(V, symmetric = TRUE)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(theta))
  Z <- matrix(stats::rnorm(n * length(theta)), n)
  out <- sweep(Z %*% t(A), 2, theta, "+")
  colnames(out) <- names(theta)
  out
}

nm_rhs <- function(f) { f <- trimws(f); if (!startsWith(f, "~")) f <- paste("~", f); f }

# ---------------------------------------------------------------- unmarked
nm_umf <- function(des) {
  if (des$model == "static") {
    y <- matrix(des$y, des$M, des$J)
    oc <- des$odf; if (!ncol(oc)) oc <- NULL
    sc <- des$sc; if (!ncol(sc)) sc <- NULL
    unmarked::unmarkedFramePCount(y = y, siteCovs = sc, obsCovs = oc)
  } else {
    y <- matrix(aperm(des$y, c(1, 2, 3)), des$M, des$J * des$T)   # columns: season1 visits, season2 visits, ...
    oc <- des$odf; if (!ncol(oc)) oc <- NULL
    sc <- des$sc; if (!ncol(sc)) sc <- NULL
    yc <- des$ydf; if (!ncol(yc)) yc <- NULL
    unmarked::unmarkedFramePCO(y = y, siteCovs = sc, obsCovs = oc, yearlySiteCovs = yc, numPrimary = des$T)
  }
}

nm_fit_unmarked_once <- function(des, umf, K, mixture, dynamics, starts = NULL) {
  f <- des$formulas
  if (des$model == "static") {
    form <- stats::as.formula(paste(nm_rhs(f$p), nm_rhs(f$lambda)))
    args <- list(formula = form, data = umf, K = K, mixture = mixture, se = TRUE)
    if (!is.null(starts)) args$starts <- starts
    fm <- suppressWarnings(do.call(unmarked::pcount, args))
  } else {
    args <- list(lambdaformula = stats::as.formula(nm_rhs(f$lambda)),
                 gammaformula = stats::as.formula(nm_rhs(f$gamma)),
                 omegaformula = stats::as.formula(nm_rhs(f$omega)),
                 pformula = stats::as.formula(nm_rhs(f$p)),
                 data = umf, mixture = mixture, K = K, dynamics = dynamics, se = TRUE)
    if (!is.null(starts)) args$starts <- starts
    fm <- suppressWarnings(do.call(unmarked::pcountOpen, args))
  }
  fm
}

# Map unmarked's coefficient vector to the internal order/names.
nm_from_unmarked <- function(fm, des, mixture, dynamics) {
  est <- fm@estimates@estimates
  lens <- vapply(est, function(e) length(e@estimates), 1L)
  idx <- split(seq_len(sum(lens)), rep(names(est), lens))
  cf <- unlist(lapply(est, function(e) e@estimates), use.names = FALSE)
  V <- tryCatch(unmarked::vcov(fm), error = function(e) matrix(NA_real_, length(cf), length(cf)))
  pick <- c(if (des$model == "static") "state" else "lambda", "det",
            if (des$model == "dynamic" && dynamics != "notrend") "gamma",
            if (des$model == "dynamic") "omega",
            if (mixture == "NB") "alpha", if (mixture == "ZIP") "psi")
  ord <- unlist(lapply(pick, function(k) idx[[k]]))
  th <- cf[ord]; names(th) <- nm_theta_names(des, mixture, dynamics)
  V <- V[ord, ord, drop = FALSE]; dimnames(V) <- list(names(th), names(th))
  list(theta = th, vcov = V)
}

nm_to_unmarked_starts <- function(th, des, mixture, dynamics) {
  # internal order -> unmarked order (state/lambda, [gamma], [omega], det, [alpha/psi])
  g <- function(pre) th[grep(paste0("^", pre, "\\["), names(th))]
  if (des$model == "static") c(g("lam"), g("p"), if (mixture == "NB") th["log_size"], if (mixture == "ZIP") th["logit_zi"])
  else c(g("lam"), if (dynamics != "notrend") g("gam"), g("om"), g("p"),
         if (mixture == "NB") th["log_size"], if (mixture == "ZIP") th["logit_zi"])
}

nm_engine_unmarked <- function(des, settings, K = NULL) {
  eco_require("unmarked")
  mixture <- settings$mixture; dyn <- settings$dynamics
  umf <- nm_umf(des)
  ymax <- max(des$y, na.rm = TRUE)
  auto <- is.null(K) && is.null(settings[["K"]])
  K <- K %||% settings[["K"]] %||% (ymax + max(30, 2 * ymax))
  notes <- character()
  for (attempt in 1:4) {
    fm <- tryCatch(nm_fit_unmarked_once(des, umf, K, mixture, dyn), error = function(e) e)
    if (inherits(fm, "error")) stop("unmarked fit failed: ", conditionMessage(fm))
    m <- nm_from_unmarked(fm, des, mixture, dyn)
    if (!auto) break
    u <- nm_unpack(m$theta, des, mixture, dyn)
    tp <- nm_tail_prob(u, K, mixture, des)
    if (!is.finite(tp) || tp < settings$K_tail_tol) break
    if (attempt == 4) { notes <- c(notes, sprintf("K = %d still truncates the abundance distribution (tail probability %.2g): abundance may be unbounded (see Known limits).", K, tp)); break }
    K <- min(ceiling(K * 2), 2000)
  }
  conv <- fm@opt$convergence == 0
  se <- sqrt(pmax(diag(m$vcov), 0))
  if (any(!is.finite(diag(m$vcov))) || any(diag(m$vcov) < 0)) { conv <- FALSE; notes <- c(notes, "Hessian not invertible: some parameters are not identifiable from these data.") }
  if (any(se > 5, na.rm = TRUE)) {
    conv <- FALSE
    notes <- c(notes, sprintf("Very large standard errors (%s): estimates are poorly identified (often abundance and detection trading off).",
                              paste(names(se)[se > 5], collapse = ", ")))
  }
  if (fm@opt$convergence != 0) notes <- c(notes, "optimiser did not report convergence")
  list(engine = "unmarked", theta = m$theta, vcov = m$vcov, loglik = -fm@negLogLike, aic = fm@AIC,
       K = K, converged = conv, notes = notes, fm = fm,
       draws = nm_mvn_draws(m$theta, m$vcov, settings$n_draws))
}

# ------------------------------------------------- hand-written ML (rlik)
nm_engine_rlik <- function(des, settings, K = NULL) {
  mixture <- settings$mixture; dyn <- settings$dynamics
  nms <- nm_theta_names(des, mixture, dyn)
  ymax <- max(des$y, na.rm = TRUE)
  K <- K %||% settings[["K"]] %||% (ymax + max(30, 2 * ymax))
  th0 <- stats::setNames(rep(0, length(nms)), nms)
  th0["lam[1]"] <- log(mean(apply(matrix(des$y, des$M), 1, max, na.rm = TRUE)) + 1)
  nll <- function(th) { names(th) <- nms; v <- -nm_loglik(th, des, K, mixture, dyn); if (!is.finite(v)) 1e10 else v }
  o <- stats::optim(th0, nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  H <- stats::optimHess(o$par, nll)
  V <- tryCatch(solve(H), error = function(e) matrix(NA_real_, length(nms), length(nms)))
  th <- stats::setNames(o$par, nms); dimnames(V) <- list(nms, nms)
  conv <- o$convergence == 0 && all(is.finite(diag(V))) && all(diag(V) > 0)
  list(engine = "rlik", theta = th, vcov = V, loglik = -o$value, aic = 2 * o$value + 2 * length(nms),
       K = K, converged = conv, notes = character(), fm = o,
       draws = nm_mvn_draws(th, V, settings$n_draws))
}

# ---------------------------------------------------------------- nimble
nm_register_zip <- function() {
  if (exists("dZIPNmix_v", envir = globalenv())) return(invisible())
  suppressPackageStartupMessages(library(nimbleEcology))
  assign("dZIPNmix_v", nimble::nimbleFunction(
    run = function(x = double(1), lambda = double(), zi = double(), prob = double(1),
                   len = integer(), log = integer(0, default = 0)) {
      returnType(double())
      lp <- dNmixture_v(x, lambda, prob, -1, -1, len, 1)
      if (sum(x) == 0) ll <- log(zi + (1 - zi) * exp(lp)) else ll <- log(1 - zi) + lp
      if (log) return(ll)
      return(exp(ll))
    }), envir = globalenv())
  assign("rZIPNmix_v", nimble::nimbleFunction(
    run = function(n = integer(), lambda = double(), zi = double(), prob = double(1), len = integer()) {
      returnType(double(1))
      N <- rbinom(1, 1, 1 - zi) * rpois(1, lambda)
      out <- numeric(len)
      for (j in 1:len) out[j] <- rbinom(1, N, prob[j])
      return(out)
    }), envir = globalenv())
  invisible()
}

nm_nimble_code <- function(mixture) {
  obs <- switch(mixture,
    P = quote(y[i, 1:J] ~ dNmixture_v(lambda[i], pp[i, 1:J], -1, -1, J)),
    NB = quote(y[i, 1:J] ~ dNmixture_BNB_v(lambda[i], theta_nb, pp[i, 1:J], -1, -1, J)),
    ZIP = quote(y[i, 1:J] ~ dZIPNmix_v(lambda[i], zi, pp[i, 1:J], J)))
  extra <- switch(mixture,
    P = quote({}),
    NB = quote({ log_size ~ dnorm(0, sd = sd_size); theta_nb <- exp(-log_size) }),
    ZIP = quote({ zi ~ dunif(0, 1) }))
  lam_line <- if (mixture == "ZIP") quote(log(lambda[i]) <- inprod(beta_lam[1:nl], Xlam[i, 1:nl]) - log(1 - zi))
              else quote(log(lambda[i]) <- inprod(beta_lam[1:nl], Xlam[i, 1:nl]))
  # In the nimble ZIP model beta_lam[1] is log E[N] (unmarked: log of the Poisson part); converted in nm_engine_nimble.
  eval(substitute(nimble::nimbleCode({
    for (k in 1:nl) beta_lam[k] ~ dnorm(0, sd = sd_lam[k])
    for (k in 1:np) beta_p[k] ~ dnorm(0, sd = sd_p[k])
    EXTRA
    for (i in 1:M) {
      LAM
      for (j in 1:J) {
        logit(p[i, j]) <- inprod(beta_p[1:np], Xp[i, j, 1:np])
        pp[i, j] <- p[i, j] * w[i, j]
      }
      OBS
    }
  }), list(EXTRA = extra, LAM = lam_line, OBS = obs)))
}

nm_engine_nimble <- function(des, settings) {
  if (des$model != "static") stop("the nimble engine supports the static (closed) model only; use engine: unmarked for model: dynamic")
  eco_require(c("nimble", "nimbleEcology"))
  suppressPackageStartupMessages(library(nimbleEcology))
  mixture <- settings$mixture
  if (mixture == "ZIP") nm_register_zip()
  M <- des$M; J <- des$J; nl <- ncol(des$Xlam); np <- dim(des$Xp)[4]
  y <- matrix(des$y, M, J); w <- matrix(des$w, M, J) * 1
  y[is.na(y)] <- 0
  Xp <- array(des$Xp, c(M, J, np))
  sd_lam <- c(settings$prior_sd_lambda_int, rep(settings$prior_sd_beta, nl - 1))
  sd_p <- c(settings$prior_sd_p_int, rep(settings$prior_sd_beta, np - 1))
  data <- list(y = y, w = w, Xlam = des$Xlam, Xp = Xp, sd_lam = sd_lam, sd_p = sd_p)
  if (mixture == "NB") data$sd_size <- settings$prior_sd_log_size
  const <- list(M = M, J = J, nl = nl, np = np)
  ymax <- apply(y, 1, max)
  inits <- function() {
    l <- list(beta_lam = c(log(mean(ymax) + 1) + 0.5 + stats::rnorm(1, 0, 0.2), stats::rnorm(nl - 1, 0, 0.1)),
              beta_p = c(stats::rnorm(1, 0, 0.3), stats::rnorm(np - 1, 0, 0.1)))
    if (mixture == "NB") l$log_size <- stats::rnorm(1, 1, 0.3)
    if (mixture == "ZIP") l$zi <- stats::runif(1, 0.05, 0.3)
    l
  }
  tops <- c("beta_lam", "beta_p", if (mixture == "NB") "log_size", if (mixture == "ZIP") "zi")
  cfg <- function(conf, m) {
    conf$removeSamplers(tops)
    conf$addSampler(target = m$expandNodeNames(tops), type = "AF_slice")
  }
  r <- nimble_run(nm_nimble_code(mixture), const, data, inits, monitors = tops,
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = paste("nmixture", mixture, M, J, nl, np, sep = "_"),
                  configure = cfg, seed = settings$seed)
  nms <- nm_theta_names(des, mixture)
  conv_chain <- function(s) {
    s <- as.matrix(s)
    th <- cbind(s[, paste0("beta_lam[", seq_len(nl), "]"), drop = FALSE], s[, paste0("beta_p[", seq_len(np), "]"), drop = FALSE])
    if (mixture == "NB") th <- cbind(th, s[, "log_size"])
    if (mixture == "ZIP") {
      zi <- s[, "zi"]
      th[, 1] <- th[, 1] - log(1 - zi)        # to unmarked convention: log of the Poisson part
      th <- cbind(th, stats::qlogis(zi))
    }
    colnames(th) <- nms
    th
  }
  chains <- lapply(r$samples, conv_chain)
  all <- do.call(rbind, chains)
  list(engine = "nimble", theta = colMeans(all), vcov = stats::cov(all), loglik = NA_real_, aic = NA_real_,
       K = max(des$y, na.rm = TRUE) + max(30, 2 * max(des$y, na.rm = TRUE)),
       converged = isTRUE(r$diagnostics$converged), notes = character(), chains = chains,
       draws = all, nimble = r)
}
