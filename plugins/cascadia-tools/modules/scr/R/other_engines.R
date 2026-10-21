# scr: NumPyro (marginal, NUTS; subprocess with JSON in/out) and secr
# (maximum likelihood; Borchers & Efford 2008; Efford secr package) engines.

scr_engine_numpyro <- function(prep, settings) {
  if (!scr_is_basic(prep))
    stop("the NumPyro engine supports binomial/Bernoulli/Poisson models with density covariates and detector-level detection covariates; ",
         "occasion covariates, sex-specific detection, sigma covariates and behavioural response need engine: nimble (or secr)")
  eco_require("jsonlite")
  py <- eco_python()
  if (!file.exists(py)) stop("NumPyro engine needs the Python env at ", py)
  td <- tempfile("scr_np_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ini <- scr_inits_top(prep, prep$M)
  pois <- prep$obs == "poisson"
  PD <- ncol(prep$XD); PP <- ncol(prep$XP)
  init <- list(psi = ini$psi, det0 = if (pois) ini$lam0 else ini$p0, sigma = ini$sigma)
  if (PD) init$bD <- I(as.numeric(ini$beta[seq_len(PD)]))
  if (PP) init$bP <- I(as.numeric(ini$beta[PD + seq_len(PP)]))
  inp <- list(d2 = prep$d2, use = I(as.numeric(prep$use[, 1])), XD = prep$XD, XP = prep$XP,
              PD = PD, PP = PP, poisson = pois, lam_max = prep$lam_max,
              nzi = I(prep$nzi - 1L), nzj = I(prep$nzj - 1L), nzc = I(as.numeric(prep$nzc)), n = prep$n, M = prep$M,
              sigma_max = prep$sigma_max, chains = settings$nchains,
              warmup = settings$numpyro_warmup, samples = settings$numpyro_samples,
              seed = settings$seed, x64 = !isTRUE(settings$numpyro_float32), init = init)
  fi <- file.path(td, "in.json"); fo <- file.path(td, "out.json")
  jsonlite::write_json(inp, fi, digits = NA, auto_unbox = TRUE, matrix = "rowmajor")
  env <- c("XLA_PYTHON_CLIENT_PREALLOCATE=false",
           if (!isTRUE(settings$numpyro_gpu)) "JAX_PLATFORMS=cpu",
           "XLA_FLAGS=--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
  t0 <- proc.time()[["elapsed"]]
  msg <- system2(py, c(shQuote(file.path(MODULE_DIR, "python", "scr_numpyro.py")), shQuote(fi), shQuote(fo)),
                 env = env, stdout = TRUE, stderr = TRUE)
  if (!file.exists(fo)) stop("NumPyro run failed:\n", paste(utils::tail(msg, 20), collapse = "\n"))
  res <- jsonlite::read_json(fo, simplifyVector = TRUE)
  dr <- res$draws
  C <- settings$nchains
  chains <- lapply(seq_len(C), function(k) {
    m <- sapply(names(dr), function(v) as.numeric(matrix(dr[[v]], nrow = C)[k, ]))
    if (!is.matrix(m)) m <- matrix(m, nrow = 1, dimnames = list(NULL, names(dr)))
    # bD1.., bP1.. -> beta layout (density first, then detection)
    cn <- colnames(m)
    bcols <- c(grep("^bD", cn, value = TRUE), grep("^bP", cn, value = TRUE))
    out <- m[, c("psi", "det0", "sigma"), drop = FALSE]
    if (length(bcols)) { bb <- m[, bcols, drop = FALSE]; colnames(bb) <- paste0("beta", seq_along(bcols)); out <- cbind(out, bb) }
    out
  })
  list(chains = chains, t_total = proc.time()[["elapsed"]] - t0, engine = "numpyro_marginal",
       timing = res$timing)
}

# ---------------------------------------------------------------- secr
# Builds a secr capthist + mask identical to the Bayesian engines' (same mask
# pixels, same standardised covariates, coordinates in internal units x 1000
# so secr sees metre-like numbers). Detector type: proximity (binomial /
# Bernoulli, half-normal g0) or count with binomN = 0 (Poisson, hazard
# half-normal lambda0); usage = effort; occasion covariates via timecov; sex
# via the hybrid mixture hcov (finite mixture over sex with known class for sexed animals; Pledger 2000); trap-specific
# behavioural response via bk.
scr_secr_build <- function(prep) {
  eco_require("secr")
  suppressPackageStartupMessages(requireNamespace("secr"))
  f <- 1000
  pois <- prep$obs == "poisson"
  tdf <- data.frame(x = prep$txy[, 1] * f, y = prep$txy[, 2] * f)
  rownames(tdf) <- prep$traps$id
  tr <- secr::read.traps(data = tdf, detector = if (pois) "count" else "proximity")
  secr::usage(tr) <- prep$usage
  tc <- unique(c(prep$pcov, prep$scov))
  if (length(tc)) {
    cv <- data.frame(row.names = prep$traps$id)
    for (v in prep$pcov) cv[[v]] <- prep$XP[, v]
    for (v in prep$scov) cv[[v]] <- prep$XS[, v]
    secr::covariates(tr) <- cv
  }
  # one capture row per detection (per sample for Poisson)
  w <- which(prep$y3 > 0, arr.ind = TRUE)
  reps <- if (pois) prep$y3[w] else rep(1, nrow(w))
  cap <- data.frame(Session = 1, ID = prep$ids[rep(w[, 1], reps)], Occasion = rep(w[, 3], reps),
                    TrapID = prep$traps$id[rep(w[, 2], reps)])
  ch <- secr::make.capthist(cap, tr, fmt = "trapID", noccasions = prep$K)
  if (prep$S == 2) {
    sx <- factor(ifelse(prep$sexi == 0, "F", ifelse(prep$sexi == 1, "M", NA)), levels = c("F", "M"))
    secr::covariates(ch) <- data.frame(sex = sx[match(rownames(ch), prep$ids)])
  }
  mdf <- data.frame(x = prep$mxy[, 1] * f, y = prep$mxy[, 2] * f)
  msk <- secr::read.mask(data = mdf, spacing = prep$pix * f)
  if (length(prep$dcov)) secr::covariates(msk) <- as.data.frame(prep$XD)
  dp <- if (pois) "lambda0" else "g0"
  rhs <- function(v) if (length(v)) paste(v, collapse = " + ") else "1"
  fD <- stats::as.formula(paste("D ~", rhs(prep$dcov)))
  fg <- stats::as.formula(paste(dp, "~", rhs(c(prep$pcov, prep$ocov, if ("p0" %in% prep$sexeff) "h2", if (prep$behav == "trap") "bk"))))
  fs <- stats::as.formula(paste("sigma ~", rhs(c(prep$scov, if ("sigma" %in% prep$sexeff) "h2"))))
  timecov <- if (length(prep$ocov)) as.data.frame(prep$XO) else NULL
  list(ch = ch, mask = msk, model = list(fD, fg, fs), detectfn = if (pois) "HHN" else "HN",
       binomN = if (pois) 0 else NULL, timecov = timecov, hcov = if (prep$S == 2) "sex" else NULL, dp = dp, f = f)
}

# secr link-scale coefficient vector -> internal draw (psi, det0, sigma, beta*, phi)
scr_secr_to_th <- function(prep, b, Mbig, f = 1000) {
  g <- function(nm) {
    for (cand in c(nm, paste0(nm, "TRUE"), paste0(nm, "1"), sub("h22$", "h2M", nm), paste0(nm, ".h2M"))) if (cand %in% names(b)) return(b[[cand]])
    stop("secr coefficient '", nm, "' not found (have: ", paste(names(b), collapse = ", "), ")")
  }
  pois <- prep$obs == "poisson"; dp <- if (pois) "lambda0" else "g0"
  bet <- numeric(prep$nb)
  PD <- length(prep$dcov)
  if (PD) bet[seq_len(PD)] <- vapply(paste0("D.", prep$dcov), g, 1)
  pc <- c(prep$pcov, prep$ocov)
  if (length(pc)) bet[prep$bidx[3] + seq_along(pc)] <- vapply(paste0(dp, ".", pc), g, 1)
  if (length(prep$scov)) bet[prep$bidx[5] + seq_along(prep$scov)] <- vapply(paste0("sigma.", prep$scov), g, 1)
  if (prep$bidx[7] > 0) bet[prep$bidx[7]] <- g(paste0(dp, ".h22"))
  if (prep$bidx[8] > 0) bet[prep$bidx[8]] <- g("sigma.h22")
  if (prep$bidx[9] > 0) bet[prep$bidx[9]] <- g(paste0(dp, ".bk"))
  lw <- if (PD) as.vector(prep$XD %*% bet[seq_len(PD)]) else rep(0, prep$G)
  cell_ha <- (prep$pix * f)^2 / 1e4
  EN <- sum(exp(g("D") + lw)) * cell_ha
  th <- c(psi = min(EN / Mbig, 0.999), det0 = if (pois) exp(g(dp)) else stats::plogis(g(dp)), sigma = exp(g("sigma")) / f)
  if (prep$nb) th <- c(th, stats::setNames(bet, paste0("beta", seq_len(prep$nb))))
  # secr pmix (logit link) is the proportion in the SECOND class (M); proportion female = 1 - pmix
  if (prep$S == 2) th <- c(th, phi = 1 - stats::plogis(g("pmix")))
  th
}

scr_engine_secr <- function(prep, settings) {
  sb <- scr_secr_build(prep)
  t0 <- proc.time()[["elapsed"]]
  start <- list(D = max(prep$n, 1) / (prep$area * sb$f^2 / 1e4),
                sigma = if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale * sb$f else prep$trap_spacing / prep$cscale * sb$f)
  start[[sb$dp]] <- if (sb$dp == "g0") 0.2 else 0.3
  # secr's default 'fastproximity' collapses occasions into one count; with
  # non-binary usage (partial effort) that shortcut does not reproduce the exact
  # per-occasion likelihood (dev/scr/round2/secr_fast.R: the exact path,
  # fastproximity = FALSE, gives the same MLE as this module's likelihood), so
  # it is switched off whenever effort is not 0/1.
  fit <- suppressWarnings(secr::secr.fit(sb$ch, mask = sb$mask, model = sb$model, detectfn = sb$detectfn,
                                         binomN = sb$binomN, timecov = sb$timecov, hcov = sb$hcov,
                                         details = if (prep$binary_usage) list() else list(fastproximity = FALSE),
                                         trace = FALSE, start = start, ncores = 1))
  b <- fit$fit$par; names(b) <- fit$betanames
  V <- fit$beta.vcv
  Mbig <- 1000 * max(prep$n, 1)
  B <- MASS_mvrnorm(4000, b, V); colnames(B) <- names(b)
  th_hat <- scr_secr_to_th(prep, b, Mbig, sb$f)
  dr <- t(apply(B, 1, function(bb) scr_secr_to_th(prep, stats::setNames(bb, names(b)), Mbig, sb$f)))
  m_hat <- scr_report_chain(prep, cbind(rbind(th_hat), N = NA), 1)
  m_dr <- scr_report_chain(prep, cbind(dr, N = NA), seq_len(nrow(dr)))
  rows <- data.frame(param = colnames(m_hat), estimate = as.numeric(m_hat[1, ]),
                     sd = apply(m_dr, 2, stats::sd), lower = apply(m_dr, 2, stats::quantile, 0.025, na.rm = TRUE),
                     upper = apply(m_dr, 2, stats::quantile, 0.975, na.rm = TRUE))
  # density and abundance from region.N (expected and realised N in the mask)
  rn <- secr::region.N(fit, se.N = TRUE)
  dn <- scr_density_name(prep$units); aconv <- scr_area_conv(prep)
  EN <- rn["E.N", ]; RN <- rn["R.N", ]
  rows[rows$param == dn, c("estimate", "sd", "lower", "upper")] <- c(EN$estimate, EN$SE.estimate, EN$lcl, EN$ucl) / prep$area * aconv
  rows[rows$param == "abundance", c("estimate", "sd", "lower", "upper")] <- c(RN$estimate, RN$SE.estimate, RN$lcl, RN$ucl)
  rownames(rows) <- NULL
  ok <- !is.null(V) && all(is.finite(rows$sd)) && all(rows$sd > 0) && all(abs(b) < 20) && all(is.finite(diag(V)))
  list(estimates = rows, engine = "secr", t_total = proc.time()[["elapsed"]] - t0, fit = fit,
       converged = isTRUE(ok) && fit$fit$code %in% c(1, 2), coef = stats::coef(fit), b = b, V = V, Mbig = Mbig)
}

MASS_mvrnorm <- function(n, mu, S) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(mu))
  sweep(matrix(stats::rnorm(n * length(mu)), n) %*% t(L), 2, mu, "+")
}

# Conversion factors from internal units to reported units
scr_area_conv <- function(prep) {
  # internal area unit -> reported density denominator (km2 when units known, else map unit2)
  if (!is.na(prep$units) && prep$units %in% c("m", "km")) 1 else 1 / prep$cscale^2
}
scr_hr_conv <- function(prep) {
  # internal area -> km2 (units known) or map unit2
  if (!is.na(prep$units) && prep$units %in% c("m", "km")) 1 else prep$cscale^2
}
