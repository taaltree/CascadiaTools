# jolly-seber: maximum-likelihood POPAN cross-check with marked::crm(model = "JS")
# (Schwarz & Arnason 1996 parameterisation, Laake et al. 2013).
#
# marked parameters: Phi (logit), p (logit), pent (multinomial logit over
# occasions 2..T; occasion 1 gets the remainder), N = number of animals in
# the superpopulation never caught (log link, one per group). Derived
# quantities (Nsuper, expected N_t, entrants, growth) and their intervals come
# from 2000 draws of the betas from their asymptotic normal distribution.
# N_t here is the EXPECTED abundance E[N_t] (POPAN), not the realised N_t of
# the Bayesian fit; the two agree closely in moderate to large populations.

js_engine_marked <- function(prep, settings) {
  eco_require("marked")
  suppressPackageStartupMessages(requireNamespace("marked"))
  T <- prep$T
  use_g <- prep$G == 2 && all(prep$gcode > 0)
  notes <- character()
  if (prep$G == 2 && !use_g) notes <- "marked cross-check: animals of unknown group present, so the group effect was dropped."
  df <- data.frame(ch = apply(prep$ch, 1, paste, collapse = ""), stringsAsFactors = FALSE)
  if (use_g) df$grp <- factor(prep$group_levels[prep$gcode], levels = prep$group_levels)
  t0 <- proc.time()[["elapsed"]]
  pr <- marked::process.data(df, model = "JS", groups = if (use_g) "grp" else NULL, time.intervals = prep$dt)
  ddl <- marked::make.design.data(pr)
  pf <- "~1"
  if (prep$Kp_lev == T) pf <- "~time" else {
    for (v in names(prep$scaling)) {
      # raw (unstandardised) values: marked's starting-value step fails with
      # the standardised covariate (NA/NaN in the initial evaluation); the
      # slope is rescaled to per-SD units in js_marked_derived()
      sc <- prep$scaling[[v]]
      x <- unname(prep$Xp[, paste0("detection_beta_", v)]) * sc$scale + sc$center
      ddl$p[[v]] <- x[as.integer(as.character(ddl$p$time))]
      pf <- paste(pf, "+", v)
    }
  }
  phif <- if (prep$Kphi_lev > 1) "~time" else "~1"
  if (use_g) phif <- paste(phif, "+ grp")
  mp <- list(Phi = list(formula = stats::as.formula(phif)), p = list(formula = stats::as.formula(pf)),
             pent = list(formula = ~time), N = list(formula = if (use_g) ~grp else ~1))
  fit <- suppressMessages(suppressWarnings(utils::capture.output(
    m <- marked::crm(pr, ddl, model.parameters = mp, hessian = TRUE, accumulate = FALSE,
                     method = "nlminb", itnmax = 5000))))
  beta <- m$results$beta
  vc <- m$results$beta.vcv
  b <- unlist(lapply(names(beta), function(k) stats::setNames(beta[[k]], paste0(k, ".", names(beta[[k]])))))
  ok <- all(is.finite(vc)) && all(diag(vc) > 0)
  set.seed(settings$seed)
  nd <- 2000
  D <- if (ok) MASS_mvrnorm(nd, b, vc) else matrix(b, 1, dimnames = list(NULL, names(b)))
  der <- function(bb) js_marked_derived(bb, prep, use_g, names(beta$p), names(beta$Phi))
  point <- der(b)
  draws <- t(apply(D, 1, der))
  est <- data.frame(param = names(point), estimate = point,
                    sd = if (ok) apply(draws, 2, stats::sd) else NA_real_,
                    lower = if (ok) apply(draws, 2, stats::quantile, 0.025, na.rm = TRUE) else NA_real_,
                    upper = if (ok) apply(draws, 2, stats::quantile, 0.975, na.rm = TRUE) else NA_real_, row.names = NULL)
  conv <- isTRUE(m$results$convergence == 0) && ok
  list(estimates = est, converged = conv, fit = m, t_total = proc.time()[["elapsed"]] - t0,
       notes = c(notes, if (!ok) "marked: Hessian not invertible; intervals unavailable."))
}

MASS_mvrnorm <- function(n, mu, S) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  ev <- pmax(e$values, 0)
  z <- matrix(stats::rnorm(n * length(mu)), n)
  sweep(z %*% t(e$vectors %*% diag(sqrt(ev), length(ev))), 2, mu, "+") |> `colnames<-`(names(mu))
}

js_marked_derived <- function(b, prep, use_g, pn, phin) {
  T <- prep$T
  g <- function(k) if (k %in% names(b)) b[[k]] else 0
  # detection
  if (prep$Kp_lev == T) {
    lp <- g("p.(Intercept)") + c(0, vapply(2:T, function(t) g(paste0("p.time", t)), 0))
  } else {
    # marked slopes are per raw unit: convert to the standardised scale
    a0 <- g("p.(Intercept)")
    for (v in names(prep$scaling)) a0 <- a0 + g(paste0("p.", v)) * prep$scaling[[v]]$center
    lp <- rep(a0, T)
    for (v in names(prep$scaling)) lp <- lp + g(paste0("p.", v)) * prep$scaling[[v]]$scale * prep$Xp[, paste0("detection_beta_", v)]
  }
  # survival per time unit (marked Phi is per unit interval; time.intervals exponentiate)
  lphi <- if (prep$Kphi_lev > 1) g("Phi.(Intercept)") + c(0, vapply(2:(T - 1), function(t) g(paste0("Phi.time", t)), 0))
          else rep(g("Phi.(Intercept)"), T - 1)
  bgr <- if (use_g) g(paste0("Phi.grp", prep$group_levels[2])) else 0
  eta <- g("pent.(Intercept)") + c(0, vapply(3:T, function(t) g(paste0("pent.time", t)), 0))
  pent <- c(1, exp(eta)) / (1 + sum(exp(eta)))
  n_g <- if (use_g) tabulate(prep$gcode, 2) else prep$n
  f0 <- if (use_g) exp(g("N.(Intercept)") + c(0, g(paste0("N.grp", prep$group_levels[2])))) else exp(g("N.(Intercept)"))
  Ng <- n_g + f0
  N <- numeric(T); B <- numeric(T - 1)
  for (k in seq_along(Ng)) {
    ph <- stats::plogis(lphi + if (k == 2) bgr else 0)^prep$dt
    a <- Ng[k] * pent[1]; N[1] <- N[1] + a
    for (t in 2:T) { a <- a * ph[t - 1] + Ng[k] * pent[t]; N[t] <- N[t] + a; B[t - 1] <- B[t - 1] + Ng[k] * pent[t] }
  }
  out <- c(superpopulation_size = sum(Ng), stats::setNames(N, paste0("abundance_", seq_len(T))),
           stats::setNames(B, paste0("entrants_", 2:T)), stats::setNames(pent, paste0("entry_prob_", seq_len(T))),
           stats::setNames(B / N[-T], paste0("recruitment_rate_", seq_len(T - 1))),
           stats::setNames(N[-1] / N[-T], paste0("growth_rate_", seq_len(T - 1))),
           growth_rate_mean = (N[T] / N[1])^(1 / sum(prep$dt)))
  if (prep$Kphi_lev > 1) out <- c(out, stats::setNames(stats::plogis(lphi), paste0("survival_", seq_len(T - 1))))
  else out <- c(out, survival = stats::plogis(lphi[1]))
  if (use_g) {
    l2 <- prep$group_levels[2]
    out <- c(out, stats::setNames(bgr, paste0("survival_beta_", l2)), stats::setNames(Ng[2] / sum(Ng), paste0("prop_", l2)))
    if (prep$Kphi_lev == 1) out <- c(out, stats::setNames(stats::plogis(lphi[1] + bgr), paste0("survival_", l2)))
  }
  if (prep$Kp_lev == T) out <- c(out, stats::setNames(stats::plogis(lp), paste0("detection_", seq_len(T))))
  else {
    out <- c(out, detection = stats::plogis(a0))
    for (v in names(prep$scaling)) out <- c(out, stats::setNames(g(paste0("p.", v)) * prep$scaling[[v]]$scale, paste0("detection_beta_", v)))
  }
  out
}
