# scr-smr: secr engine (Efford & Hunter 2018 spatial capture-mark-resight;
# secr 5.x, experimental mark-resight support). Used as an independent
# cross-check where the models coincide:
#   SMR : all occasions are sighting-only (markocc 0), marked animals known
#         beforehand (details$knownmarks = TRUE, all-zero rows for marked
#         animals never resighted), count detectors with Poisson counts
#         (binomN = 0), Tu = unmarked counts, Tm = marked-unidentified counts.
#   gSMR: marking occasions (markocc 1, proximity detectors at the marking
#         traps) followed by sighting occasions (count detectors at the
#         cameras); detection g0 differs by occasion type (~ ts).
# secr maximises a pseudo-likelihood in which the unmarked counts enter as
# independent Poisson counts and adjusts variances for their overdispersion
# (chat, by simulation). The half-normal detection function (HN) is used for
# both phases; for Poisson count detectors g0 is the expected number of
# sightings per occasion at distance 0 (lambda0 here).
# Bernoulli resighting and telemetry are not passed to secr.

smr_secr_inputs <- function(prep) {
  eco_require("secr")
  suppressPackageStartupMessages(requireNamespace("secr"))
  f <- 1000                       # internal km -> metre-like numbers for secr
  J <- prep$J; K <- prep$K; m <- prep$m
  cams <- data.frame(x = prep$txy[, 1] * f, y = prep$txy[, 2] * f)
  rownames(cams) <- prep$cameras$id
  si <- prep$sightings
  gsmr <- prep$has_marking
  if (!gsmr) {
    tr <- secr::read.traps(data = cams, detector = "count")
    secr::usage(tr) <- prep$usage
    if (length(prep$pcov)) secr::covariates(tr) <- as.data.frame(prep$XP)
    secr::markocc(tr) <- rep(0, K)
    ch <- array(0L, c(m, K, J), dimnames = list(prep$ids, seq_len(K), prep$cameras$id))
    sid <- si[si$status == "marked", ]
    for (r in seq_len(nrow(sid))) {
      i <- match(sid$animal_id[r], prep$ids); j <- match(sid$camera[r], prep$cameras$id)
      ch[i, sid$occasion[r], j] <- ch[i, sid$occasion[r], j] + 1L
    }
    occ_s <- seq_len(K); Kall <- K; Jall <- J; jcam <- seq_len(J)
  } else {
    mk <- prep$mk; L <- mk$L; Km <- mk$Km
    mt <- data.frame(x = mk$mtxy[, 1] * f, y = mk$mtxy[, 2] * f); rownames(mt) <- mk$traps$id
    all <- rbind(mt, cams)
    Kall <- Km + K; Jall <- L + J; jcam <- L + seq_len(J); occ_s <- Km + seq_len(K)
    tr <- secr::read.traps(data = all, detector = "count")
    secr::detector(tr) <- c(rep("proximity", Km), rep("count", K))
    u <- matrix(0L, Jall, Kall)
    u[seq_len(L), seq_len(Km)] <- mk$usage
    u[jcam, occ_s] <- prep$usage
    secr::usage(tr) <- u
    if (length(prep$pcov)) secr::covariates(tr) <- as.data.frame(rbind(matrix(0, L, ncol(prep$XP), dimnames = list(NULL, prep$pcov)), prep$XP))
    secr::markocc(tr) <- c(rep(1, Km), rep(0, K))
    ch <- array(0L, c(m, Kall, Jall), dimnames = list(prep$ids, seq_len(Kall), rownames(all)))
    mc <- mk$captures
    for (r in seq_len(nrow(mc))) ch[match(mc$individual[r], prep$ids), mc$occasion[r], match(mc$trap[r], mk$traps$id)] <- 1L
    sid <- si[si$status == "marked", ]
    for (r in seq_len(nrow(sid))) {
      i <- match(sid$animal_id[r], prep$ids); j <- L + match(sid$camera[r], prep$cameras$id); k <- Km + sid$occasion[r]
      ch[i, k, j] <- ch[i, k, j] + 1L
    }
  }
  class(ch) <- "capthist"
  secr::traps(ch) <- tr
  cntmat <- function(st) {
    mm <- matrix(0L, Jall, Kall)
    w <- si[si$status == st, ]
    for (r in seq_len(nrow(w))) { j <- jcam[match(w$camera[r], prep$cameras$id)]; k <- occ_s[w$occasion[r]]; mm[j, k] <- mm[j, k] + 1L }
    mm
  }
  Tu <- cntmat("unmarked")
  Tm <- if (prep$has_tm) cntmat("marked_unidentified") else NULL
  ch <- secr::addSightings(ch, unmarked = Tu, nonID = Tm, verify = FALSE)
  msk <- secr::read.mask(data = data.frame(x = prep$mxy[, 1] * f, y = prep$mxy[, 2] * f), spacing = prep$pix * f)
  if (length(prep$dcov)) secr::covariates(msk) <- as.data.frame(prep$XD)
  list(ch = ch, mask = msk, f = f)
}

smr_engine_secr <- function(prep, settings) {
  if (!prep$pois) stop("the secr cross-check supports the Poisson resight model only")
  inp <- smr_secr_inputs(prep)
  f <- inp$f
  notes <- character()
  if (prep$has_tel) notes <- c(notes, "secr fit ignores the telemetry fixes (secr's telemetry 'marking' type is not used here).")
  fD <- stats::as.formula(paste("D ~", if (length(prep$dcov)) paste(prep$dcov, collapse = " + ") else "1"))
  rhs <- c(if (prep$has_marking) "ts", prep$pcov)
  fg <- stats::as.formula(paste("g0 ~", if (length(rhs)) paste(rhs, collapse = " + ") else "1"))
  mods <- list(fD, fg)
  if (prep$has_tm) mods[[3]] <- pID ~ 1
  sg0 <- if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale * f else prep$trap_spacing / prep$cscale * f
  t0 <- proc.time()[["elapsed"]]
  # Two steps (secr ?details): fit with chat = 1, then estimate the
  # overdispersion of the unmarked (Tu) and unidentified (Tm) counts by
  # simulation at the fitted values and recompute the variances only
  # (method = "none"), which widens the intervals for overdispersed counts.
  fit <- suppressWarnings(secr::secr.fit(inp$ch, mask = inp$mask, model = mods, detectfn = "HN", binomN = 0,
                                         trace = FALSE, details = list(knownmarks = TRUE, nsim = 0), ncores = 1,
                                         start = list(D = max(prep$m, 1) * 2 / (prep$area * f^2 / 1e4), g0 = 0.2, sigma = sg0)))
  nsim <- if (is.null(settings$secr_nsim)) 0 else settings$secr_nsim
  if (nsim > 0) {
    fit2 <- tryCatch(suppressWarnings(secr::secr.fit(inp$ch, mask = inp$mask, model = mods, detectfn = "HN", binomN = 0,
                                                     trace = FALSE, details = list(knownmarks = TRUE, nsim = nsim), ncores = 1,
                                                     start = fit, method = "none")),
                     error = function(e) { notes <<- c(notes, paste("secr overdispersion step failed:", conditionMessage(e))); NULL })
    if (!is.null(fit2)) fit <- fit2
  } else notes <- c(notes, "secr variances not adjusted for overdispersion of unmarked counts (secr_nsim = 0); intervals may be too narrow.")
  z <- stats::qnorm(0.975)
  cf <- stats::coef(fit)
  # density (per ha in secr) -> E[N] over the mask; realised-N style CI from log-normal
  A_ha <- prep$area * f^2 / 1e4
  dD <- secr::derived(fit, se.D = TRUE)
  pr <- stats::predict(fit)
  if (is.data.frame(pr)) pr <- list(pr)
  aconv <- smr_area_conv(prep)
  rows <- list()
  dn <- smr_density_name(prep$units)
  if (length(prep$dcov)) {
    rn <- secr::region.N(fit)
    EN <- rn["E.N", ]
    D <- EN$estimate / prep$area / aconv; Dse <- EN$SE.estimate / prep$area / aconv
    Nl <- EN$lcl; Nu <- EN$ucl; Ne <- EN$estimate; Nse <- EN$SE.estimate
  } else {
    dr <- cf["D", ]
    Ne <- exp(dr$beta) * A_ha; Nse <- Ne * dr$SE.beta
    Nl <- exp(dr$beta - z * dr$SE.beta) * A_ha; Nu <- exp(dr$beta + z * dr$SE.beta) * A_ha
  }
  rows[[1]] <- data.frame(param = dn, estimate = Ne / (prep$area * aconv), sd = Nse / (prep$area * aconv),
                          lower = Nl / (prep$area * aconv), upper = Nu / (prep$area * aconv))
  rows[[2]] <- data.frame(param = "abundance", estimate = Ne, sd = Nse, lower = Nl, upper = Nu)
  rows[[3]] <- data.frame(param = "marked_fraction", estimate = prep$m / Ne, sd = prep$m / Ne^2 * Nse,
                          lower = prep$m / Nu, upper = prep$m / Nl)
  s <- cf["sigma", ]
  sg <- function(v) exp(v) / f * prep$cscale
  rows[[4]] <- data.frame(param = "sigma", estimate = sg(s$beta), sd = sg(s$beta) * s$SE.beta,
                          lower = sg(s$beta - z * s$SE.beta), upper = sg(s$beta + z * s$SE.beta))
  hr <- function(v) pi * (sqrt(stats::qchisq(0.95, 2)) * exp(v) / f)^2 * aconv
  rows[[5]] <- data.frame(param = "home_range_area_95", estimate = hr(s$beta), sd = 2 * hr(s$beta) * s$SE.beta,
                          lower = hr(s$beta - z * s$SE.beta), upper = hr(s$beta + z * s$SE.beta))
  lg <- function(b, se, nm) data.frame(param = nm, estimate = stats::plogis(b), sd = stats::plogis(b) * (1 - stats::plogis(b)) * se,
                                       lower = stats::plogis(b - z * se), upper = stats::plogis(b + z * se))
  g0r <- cf["g0", ]
  if (prep$has_marking) {
    # ts: level "marking" is the reference in secr; sighting = g0 + g0.tssighting
    rows[[length(rows) + 1]] <- lg(g0r$beta, g0r$SE.beta, "marking_detection")
    V <- fit$beta.vcv; bn <- rownames(cf)
    ii <- match(c("g0", "g0.tssighting"), bn)
    b <- sum(cf$beta[ii]); se <- sqrt(sum(V[ii, ii]))
    rows[[length(rows) + 1]] <- lg(b, se, smr_rate_name(prep))
  } else rows[[length(rows) + 1]] <- lg(g0r$beta, g0r$SE.beta, smr_rate_name(prep))
  if (prep$has_tm && "pID" %in% rownames(cf)) rows[[length(rows) + 1]] <- lg(cf["pID", "beta"], cf["pID", "SE.beta"], "id_probability")
  for (v in prep$dcov) { r <- cf[paste0("D.", v), ]
    rows[[length(rows) + 1]] <- data.frame(param = paste0("density_beta_", v), estimate = r$beta, sd = r$SE.beta,
                                           lower = r$beta - z * r$SE.beta, upper = r$beta + z * r$SE.beta) }
  for (v in prep$pcov) { r <- cf[paste0("g0.", v), ]
    rows[[length(rows) + 1]] <- data.frame(param = paste0("detection_beta_", v), estimate = r$beta, sd = r$SE.beta,
                                           lower = r$beta - z * r$SE.beta, upper = r$beta + z * r$SE.beta) }
  est <- do.call(rbind, rows); rownames(est) <- NULL
  ch <- fit$details$chat
  if (!is.null(ch)) notes <- c(notes, sprintf("secr overdispersion of sighting counts (chat, Tu and Tm): %s.", paste(format(signif(unlist(ch), 3)), collapse = ", ")))
  notes <- c(notes, "secr detection is logit-link g0 (bounded at 1), the Bayesian Poisson model uses lambda0 with a log link; detection-covariate effects are therefore on different scales (logit vs log) and not directly comparable when lambda0 is not small.")
  ok <- all(is.finite(est$sd)) && (is.null(fit$fit$code) || fit$fit$code %in% c(0, 1, 2))
  list(estimates = est, fit = fit, converged = isTRUE(ok), t_total = proc.time()[["elapsed"]] - t0, notes = notes)
}
