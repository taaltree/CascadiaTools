# cjs engines. Each returns
#   list(engine, theta (named, internal order: phi coefs, p coefs[, phi:sd_time]),
#        vcov (likelihood engines), chains (Nimble: list of per-chain theta matrices),
#        eps_chains (Nimble with random time effect), loglik, aic, converged, notes)

cjs_marked_reserved <- c("age", "Age", "group", "time", "Time", "cohort", "Cohort", "occ", "id", "ind",
                         "ch", "freq", "fix", "order", "initial.age")

cjs_norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# ---------------------------------------------------------------- marked
cjs_engine_marked <- function(des, settings) {
  eco_require("marked")
  suppressPackageStartupMessages(requireNamespace("marked"))
  if (des$random_time) stop("random time effects ('(1 | time)') need engine: nimble")
  data <- des$data; ch <- data$ch; n <- des$n; T <- des$T
  df <- data.frame(ch = apply(ch, 1, paste, collapse = ""), ind = seq_len(n), stringsAsFactors = FALSE)
  bt <- suppressWarnings(as.numeric(data$occasions[1]))
  numeric_labels <- all(is.finite(suppressWarnings(as.numeric(data$occasions))))
  if (!numeric_labels || !is.finite(bt)) bt <- 1
  pr <- marked::process.data(df, model = "CJS", begin.time = bt, time.intervals = data$intervals, accumulate = FALSE)
  ddl <- marked::make.design.data(pr)
  add <- function(dl, my, mmo, Tcol) {
    idx <- (dl$ind - 1) * Tcol + dl$occ
    for (v in setdiff(mmo$vars, c("time", "Time"))) {
      val <- my[[v]][idx]
      if (is.factor(val)) val <- factor(sub("^_", "", as.character(val)), levels = sub("^_", "", levels(val)))
      dl[[v]] <- val
    }
    dl
  }
  ddl$Phi <- add(ddl$Phi, des$dd$phi, des$mphi, T - 1)
  ddl$p <- add(ddl$p, des$dd$p, des$mp, T)
  pform <- des$formulas$p
  if (isTRUE(des$fixed_last)) {
    # last-occasion p fixed to 1 (phi(t) p(t) confounding); our time factor has the last
    # occasion merged into the reference level, passed to marked as 'ptime'
    idx <- (ddl$p$ind - 1) * T + ddl$p$occ
    v <- des$dd$p$time[idx]
    ddl$p$ptime <- factor(sub("^_", "", as.character(v)), levels = sub("^_", "", levels(v)))
    ddl$p$fix <- ifelse(des$dd$p$.fixT[idx] == 1, 1, NA)
    pform <- gsub("\\btime\\b", "ptime", pform)
  }
  mp <- list(Phi = list(formula = stats::as.formula(des$formulas$phi)), p = list(formula = stats::as.formula(pform)))
  t0 <- proc.time()[["elapsed"]]
  # marked's default (non-TMB) optimiser returned wrong MLEs whenever covariates were in the model
  # (marked 1.2.8 here; -2lnL worse than the nested intercept-only model, see VALIDATION.md), so the
  # TMB likelihood is used. marked compiles/loads cjsre_tmb in getwd(): use the module's tmb/ folder.
  tdir <- file.path(MODULE_DIR, "tmb")
  owd <- setwd(tdir); on.exit(setwd(owd), add = TRUE)
  out <- utils::capture.output(fm <- suppressWarnings(suppressMessages(
    marked::crm(pr, ddl, model.parameters = mp, hessian = TRUE, accumulate = FALSE, use.tmb = TRUE))))
  setwd(owd)
  b_phi <- fm$results$beta$Phi; b_p <- fm$results$beta$p
  chk <- function(mine, theirs, what) {
    if (length(mine) != length(theirs))
      stop(sprintf("marked %s design has %d columns, ecotools expects %d", what, length(theirs), length(mine)))
    nt <- !grepl("^time|^tsm_first_time", mine)
    if (any(cjs_norm(mine[nt]) != cjs_norm(theirs[nt])))
      stop("marked ", what, " coefficient order differs: ", paste(theirs, collapse = ","), " vs ", paste(mine, collapse = ","))
  }
  chk(colnames(des$Xphi), names(b_phi), "survival"); chk(colnames(des$Xp), names(b_p), "detection")
  th <- c(b_phi, b_p); names(th) <- des$theta_names
  V <- fm$results$beta.vcv
  if (is.null(V)) V <- matrix(NA_real_, length(th), length(th))
  V <- as.matrix(V); dimnames(V) <- list(names(th), names(th))
  conv <- isTRUE(fm$results$convergence == 0) && all(is.finite(V)) &&
    all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0)
  notes <- character()
  if (!conv) notes <- c(notes, "marked: optimiser did not converge or the Hessian is not invertible (some parameters not estimable).")
  list(engine = "marked", theta = th, vcov = V, loglik = -fm$results$neg2lnl / 2, aic = fm$results$AIC,
       converged = conv, notes = notes, fm = fm, time = proc.time()[["elapsed"]] - t0)
}

# ---------------------------------------------------------------- Nimble
cjs_nimble_code <- function(Kphi, Kp, int_phi, int_p, re, pow, latent) {
  L <- character()
  pri <- function(b, K, int) {
    s <- character()
    if (int) s <- c(s, sprintf("%s[1] ~ dlogis(0, 1)", b))
    k0 <- if (int) 2 else 1
    if (K >= k0) s <- c(s, sprintf("for (k in %d:%d) { %s[k] ~ dnorm(0, sd = %s) }", k0, K, b, "prior_sd"))
    s
  }
  L <- c(L, pri("bphi", Kphi, int_phi), pri("bp", Kp, int_p))
  xb <- function(X, b, R, K, Kn) if (K == 1) sprintf("%s[1:%s, 1] * %s[1]", X, R, b) else sprintf("(%s[1:%s, 1:%s] %%*%% %s[1:%s])[1:%s, 1]", X, R, Kn, b, Kn, R)
  lp <- paste("lpphi[1:Rphi] <-", xb("Xphi", "bphi", "Rphi", Kphi, "Kphi"))
  if (re) {
    # non-centred random interval effects: eps[t] = sd_time * zeps[t] (mixes far better than
    # eps[t] ~ N(0, sd_time) when there are few intervals)
    L <- c(L, "sd_time ~ T(dnorm(0, sd = 1), 0, )",
           "for (t in 1:Tm1) { zeps[t] ~ dnorm(0, 1)", "  eps[t] <- sd_time * zeps[t] }")
    lp <- paste("lpphi[1:Rphi] <-", xb("Xphi", "bphi", "Rphi", Kphi, "Kphi"), "+ (Z[1:Rphi, 1:Tm1] %*% eps[1:Tm1])[1:Rphi, 1]")
  }
  L <- c(L, lp,
         # phi^dt for unequal intervals (vectorised pow() needs a scalar exponent in NIMBLE)
         if (pow) "phiv[1:Rphi] <- prephi[1:Rphi] + (1 - prephi[1:Rphi]) * exp(dt[1:Rphi] * log(expit(lpphi[1:Rphi])))"
         else "phiv[1:Rphi] <- prephi[1:Rphi] + (1 - prephi[1:Rphi]) * expit(lpphi[1:Rphi])",
         paste("lpp[1:Rp] <-", xb("Xp", "bp", "Rp", Kp, "Kp")),
         "pv[1:Rp] <- upto[1:Rp] + (1 - upto[1:Rp]) * expit(lpp[1:Rp])")
  if (!latent) {
    L <- c(L, "for (i in 1:n) { y[i, 1:T] ~ dCJS_vv(probSurvive = phiv[(Tm1 * (i - 1) + 1):(Tm1 * i)], probCapture = pv[(T * (i - 1) + 1):(T * i)], len = T) }")
  } else {
    L <- c(L, "for (i in 1:n) {",
           "  z[i, 1] ~ dbern(1)",
           "  for (t in 2:T) { z[i, t] ~ dbern(z[i, t - 1] * phiv[Tm1 * (i - 1) + t - 1]) }",
           "  for (t in 1:T) { y[i, t] ~ dbern(z[i, t] * pv[T * (i - 1) + t]) }",
           "}")
  }
  str2lang(paste0("{\n", paste(L, collapse = "\n"), "\n}"))
}

cjs_engine_nimble <- function(des, settings) {
  eco_require(c("nimble", "nimbleEcology"))
  suppressPackageStartupMessages({ requireNamespace("nimble"); requireNamespace("nimbleEcology") })
  # nimbleEcology's distributions must be visible to nimbleModel()
  if (!"package:nimbleEcology" %in% search()) suppressPackageStartupMessages(attachNamespace("nimbleEcology"))
  if (!"package:nimble" %in% search()) suppressPackageStartupMessages(attachNamespace("nimble"))
  latent <- identical(settings$formulation, "latent")
  n <- des$n; T <- des$T; Tm1 <- T - 1
  Xphi <- des$Xphi; Xp <- des$Xp
  Kphi <- ncol(Xphi); Kp <- ncol(Xp)
  int_phi <- colnames(Xphi)[1] == "(Intercept)"; int_p <- colnames(Xp)[1] == "(Intercept)"
  re <- des$random_time
  dt <- des$dd$phi$.dt
  pow <- any(abs(dt - 1) > 1e-12)
  ch <- des$data$ch
  first <- des$first
  ypad <- ch
  for (i in seq_len(n)) if (first[i] > 1) ypad[i, seq_len(first[i] - 1)] <- 1
  code <- cjs_nimble_code(Kphi, Kp, int_phi, int_p, re, pow, latent)
  constants <- list(n = n, T = T, Tm1 = Tm1, Rphi = n * Tm1, Rp = n * T, Kphi = Kphi, Kp = Kp)
  dat <- list(y = ypad, Xphi = Xphi, Xp = Xp, prephi = des$dd$phi$.pre, upto = pmax(des$dd$p$.upto, des$dd$p$.fixT),
              prior_sd = settings$prior_sd %||% 2)
  if (pow) dat$dt <- dt
  if (re) dat$Z <- outer(des$dd$phi$.t, seq_len(Tm1), "==") * 1
  last <- cjs_last(ch)
  inits <- function() {
    b1 <- c(stats::qlogis(0.6) + stats::rnorm(1, 0, 0.3), stats::rnorm(Kphi - 1, 0, 0.2))
    b2 <- c(stats::qlogis(0.5) + stats::rnorm(1, 0, 0.3), stats::rnorm(Kp - 1, 0, 0.2))
    if (!int_phi) b1 <- stats::rnorm(Kphi, 0, 0.3)
    if (!int_p) b2 <- stats::rnorm(Kp, 0, 0.3)
    l <- list(bphi = b1[seq_len(Kphi)], bp = b2[seq_len(Kp)])
    if (re) { l$sd_time <- stats::runif(1, 0.1, 0.5); l$zeps <- stats::rnorm(Tm1, 0, 0.5) }
    if (latent) l$z <- outer(seq_len(n), seq_len(T), function(i, t) as.numeric(t <= last[i]))
    l
  }
  configure <- function(conf, m) {
    conf$removeSamplers(c("bphi", "bp", if (re) c("sd_time", "zeps")))
    if (re) {
      # survival intercept, random effects and their SD are strongly correlated: one block
      conf$addSampler(c(sprintf("bphi[1:%d]", Kphi), sprintf("zeps[1:%d]", Tm1), "sd_time"), "AF_slice")
    } else if (Kphi > 1) conf$addSampler(sprintf("bphi[1:%d]", Kphi), "AF_slice") else conf$addSampler("bphi[1]", "slice")
    if (Kp > 1) conf$addSampler(sprintf("bp[1:%d]", Kp), "AF_slice") else conf$addSampler("bp[1]", "slice")
  }
  key <- paste("cjs", if (latent) "latent" else "marg", n, T, Kphi, Kp, int_phi, int_p, re, pow, sep = "_")
  t0 <- proc.time()[["elapsed"]]
  r <- nimble_run(code, constants, dat, inits, monitors = c("bphi", "bp", if (re) c("sd_time", "eps")),
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin %||% 1, cache_key = if (isTRUE(settings$cache)) key else NULL,
                  configure = configure, seed = settings$seed)
  thn <- des$theta_names
  map <- c(sprintf("bphi[%d]", seq_len(Kphi)), sprintf("bp[%d]", seq_len(Kp)), if (re) "sd_time")
  chains <- lapply(r$samples, function(s) { m <- as.matrix(s)[, map, drop = FALSE]; colnames(m) <- thn; m })
  eps_chains <- if (re) lapply(r$samples, function(s) as.matrix(s)[, sprintf("eps[%d]", seq_len(Tm1)), drop = FALSE]) else NULL
  draws <- do.call(rbind, chains)
  list(engine = if (latent) "nimble_latent" else "nimble", theta = colMeans(draws), chains = chains, draws = draws,
       eps_chains = eps_chains, converged = TRUE, notes = character(), time = proc.time()[["elapsed"]] - t0,
       t_setup = r$diagnostics$t_setup)
}
