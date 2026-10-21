# scr-search: Nimble engine. Marginal (activity centres and z summed over
# the habitat mask) data-augmentation likelihood for search-encounter data
# discretised into detector cells with effort (Royle et al. 2011; Thompson
# et al. 2012; Milleret et al. 2019; Moqanaki et al. 2021), with LOCAL
# EVALUATION: the expected number of samples an animal centred at pixel g
# leaves in cell c is only computed for pixel-cell pairs closer than a radius
# R, and a detected animal's centre is only evaluated on pixels within R of
# every cell where it was detected (Milleret et al. 2019's "local evaluation
# of the individual state space", applied to the marginal likelihood).
#
#   lambda[g, u] = lambda0 * E_u * exp(w_u' alpha) * exp(-d(g, c(u))^2 / (2 sigma^2))
#   u = (cell, track) search unit, E_u = length of that track inside the cell
#   Poisson:   y[i, u] ~ Poisson(lambda)          (number of samples)
#   Bernoulli: y[i, u] ~ Bernoulli(1 - exp(-lambda))   (detected in unit or not)
#   Lambda[g] = sum_u lambda[g, u]; P(never detected | g) = exp(-Lambda[g])
#
# Arrays are padded to geometric buckets so recovery replicates reuse one
# compiled model (see nimble_run()).

ss_define_nimble <- function() {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))
  if (exists("dSearchMarg", envir = .GlobalEnv, inherits = FALSE)) return(invisible())
  nf <- nimble::nimbleFunction
  dSearchMarg <- nf(
    run = function(x = double(0), psi = double(0), lam0 = double(0), sigma = double(0),
                   lpi = double(1), lEu = double(1), ucell = double(1), mxy = double(2), cxy = double(2),
                   pg = double(1), pj = double(1), pd2 = double(1),
                   rcell = double(1), runit = double(1), rcnt = double(1), rstart = double(1), rend = double(1),
                   sstart = double(1), send = double(1), spix = double(1),
                   nobs = double(0), M = double(0), bern = double(0),
                   G = integer(0), C = integer(0), U = integer(0), P = integer(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      a <- -1 / (2 * sigma^2)
      ll0 <- log(lam0)
      Et <- numeric(C)
      for (u in 1:U) Et[ucell[u]] <- Et[ucell[u]] + exp(lEu[u])
      Lam <- numeric(G)
      for (p in 1:P) Lam[pg[p]] <- Lam[pg[p]] + Et[pj[p]] * exp(a * pd2[p])
      base <- numeric(G)
      for (g in 1:G) base[g] <- lpi[g] - lam0 * Lam[g]
      tot <- 0
      nob <- nobs
      if (nob > 0) {
        for (i in 1:nob) {
          m <- -Inf; s <- 0
          for (k in sstart[i]:send[i]) {
            g <- spix[k]
            v <- base[g]
            for (r in rstart[i]:rend[i]) {
              j <- rcell[r]
              d2 <- (mxy[g, 1] - cxy[j, 1])^2 + (mxy[g, 2] - cxy[j, 2])^2
              lr <- ll0 + lEu[runit[r]] + a * d2
              if (bern > 0) {
                lam <- exp(lr)
                if (lam > 30) {
                  v <- v + lam
                } else if (lam < 1e-6) {
                  v <- v + lr + 0.5 * lam
                } else {
                  v <- v + log(exp(lam) - 1)
                }
              } else {
                v <- v + rcnt[r] * lr
              }
            }
            if (v > m) {
              s <- s * exp(m - v) + 1
              m <- v
            } else {
              s <- s + exp(v - m)
            }
          }
          tot <- tot + m + log(s)
        }
      }
      mx0 <- max(base)
      lf0 <- mx0 + log(sum(exp(base - mx0)))
      tot <- tot + nobs * log(psi) + (M - nobs) * log(psi * exp(lf0) + 1 - psi)
      if (log) return(tot) else return(exp(tot))
    })
  rSearchMarg <- nf(
    run = function(n = integer(0), psi = double(0), lam0 = double(0), sigma = double(0),
                   lpi = double(1), lEu = double(1), ucell = double(1), mxy = double(2), cxy = double(2),
                   pg = double(1), pj = double(1), pd2 = double(1),
                   rcell = double(1), runit = double(1), rcnt = double(1), rstart = double(1), rend = double(1),
                   sstart = double(1), send = double(1), spix = double(1),
                   nobs = double(0), M = double(0), bern = double(0),
                   G = integer(0), C = integer(0), U = integer(0), P = integer(0)) {
      returnType(double(0)); return(1)
    })
  assign("dSearchMarg", dSearchMarg, envir = .GlobalEnv)
  assign("rSearchMarg", rSearchMarg, envir = .GlobalEnv)
  invisible()
}

ss_code <- function() nimble::nimbleCode({
  psi ~ dunif(0, 1)
  lam0 ~ dunif(0, lam_max)
  sigma ~ dunif(0, sigma_max)
  if (HASD) {
    for (k in 1:PD) bD[k] ~ dnorm(0, sd = 2)
    if (PD == 1) {
      lw[1:G] <- XD[1:G, 1] * bD[1]
    } else {
      lw[1:G] <- (XD[1:G, 1:PD] %*% bD[1:PD])[1:G, 1]
    }
    lpi[1:G] <- lw[1:G] - log(sum(exp(lw[1:G])))
  } else {
    lpi[1:G] <- lpi0[1:G]
  }
  if (HASW) {
    for (k in 1:PW) bW[k] ~ dnorm(0, sd = 2)
    if (PW == 1) {
      lEu[1:U] <- logE[1:U] + W[1:U, 1] * bW[1]
    } else {
      lEu[1:U] <- logE[1:U] + (W[1:U, 1:PW] %*% bW[1:PW])[1:U, 1]
    }
  } else {
    lEu[1:U] <- logE[1:U]
  }
  ones ~ dSearchMarg(psi = psi, lam0 = lam0, sigma = sigma, lpi = lpi[1:G], lEu = lEu[1:U], ucell = ucell[1:U],
                     mxy = mxy[1:G, 1:2], cxy = cxy[1:C, 1:2], pg = pg[1:P], pj = pj[1:P], pd2 = pd2[1:P],
                     rcell = rcell[1:RB], runit = runit[1:RB], rcnt = rcnt[1:RB], rstart = rstart[1:NMAX], rend = rend[1:NMAX],
                     sstart = sstart[1:NMAX], send = send[1:NMAX], spix = spix[1:SB],
                     nobs = nobs, M = M, bern = bern, G = G, C = C, U = U, P = P)
})

ss_tops <- function(prep) c("psi", "lam0", "sigma",
                            if (ncol(prep$XD)) paste0("bD[", seq_len(ncol(prep$XD)), "]"),
                            if (ncol(prep$W)) paste0("bW[", seq_len(ncol(prep$W)), "]"))

# Rough starting values: sigma from the sample spread, lambda0 from the
# number of samples per detected animal and the effort around it.
ss_inits <- function(prep, M) {
  sg <- if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else 2 * prep$csz
  sg <- min(sg, prep$sigma_max / 2)
  Etot <- sum(exp(prep$logE))
  dens <- 2 * max(prep$n, 1) / prep$area
  lam <- prep$n_samples / max(Etot * 2 * pi * sg^2 * dens, 1e-6)
  lam <- min(max(lam, 1e-3), 0.2 * ss_lam_max(prep))
  ini <- list(psi = min(0.9, 2 * max(prep$n, 1) / M) * stats::runif(1, 0.8, 1.2),
              lam0 = lam * stats::runif(1, 0.7, 1.3), sigma = sg * stats::runif(1, 0.8, 1.25))
  if (ncol(prep$XD)) ini$bD <- stats::rnorm(ncol(prep$XD), 0, 0.2)
  if (ncol(prep$W)) ini$bW <- stats::rnorm(ncol(prep$W), 0, 0.2)
  ini
}

# Upper bound of the flat prior on lambda0: 50 samples per internal length unit
# (per km) at the activity centre, far above any search-encounter study.
ss_lam_max <- function(prep) 50

ss_engine_nimble <- function(prep, settings) {
  ss_define_nimble()
  PD <- ncol(prep$XD); PW <- ncol(prep$W)
  n <- prep$n
  P <- ss_bucket(nrow(prep$pairs), 2000); CB <- ss_bucket(prep$C, 50); UB <- ss_bucket(prep$U, 100)
  rec <- prep$rec
  RB <- ss_bucket(max(nrow(rec), 1), 100); NMAX <- ss_bucket(max(n, 1), 50)
  setlen <- vapply(prep$sets, length, 1L)
  SB <- ss_bucket(max(sum(setlen), 1), 2000)
  pad <- function(v, L, val = 1) c(v, rep(val, L - length(v)))
  rstart <- match(seq_len(n), rec$ind); rend <- vapply(seq_len(n), function(i) max(which(rec$ind == i)), 1L)
  send <- cumsum(setlen); sstart <- send - setlen + 1
  cxy <- rbind(prep$cxy, matrix(1e4, CB - prep$C, 2))
  W <- if (PW) rbind(prep$W, matrix(0, UB - prep$U, PW)) else matrix(0, UB, 1)
  constants <- list(G = prep$G, C = CB, U = UB, P = P, RB = RB, NMAX = NMAX, SB = SB,
                    HASD = PD > 0, HASW = PW > 0, PD = max(PD, 1), PW = max(PW, 1))
  data <- list(ones = 1, mxy = prep$mxy, cxy = cxy,
               logE = pad(prep$logE, UB, -1e3), ucell = pad(prep$ucell, UB, 1), W = W,
               pg = pad(prep$pairs[, 1], P, 1), pj = pad(prep$pairs[, 2], P, 1), pd2 = pad(prep$pd2, P, 1e8),
               rcell = pad(rec$cell, RB), runit = pad(rec$unit, RB), rcnt = pad(as.numeric(rec$count), RB),
               rstart = pad(rstart, NMAX), rend = pad(rend, NMAX), sstart = pad(sstart, NMAX), send = pad(send, NMAX),
               spix = pad(unlist(prep$sets), SB), nobs = n, M = prep$M, bern = as.numeric(prep$observation == "bernoulli"),
               lam_max = ss_lam_max(prep), sigma_max = prep$sigma_max, lpi0 = rep(-log(prep$G), prep$G),
               XD = if (PD) prep$XD else matrix(0, prep$G, 1))
  tops <- ss_tops(prep)
  key <- paste("scrsearch", settings$sampler, prep$G, CB, UB, P, RB, NMAX, SB, PD, PW, sep = "_")
  r <- nimble_run(ss_code(), constants, data,
                  inits = function() ss_inits(prep, prep$M),
                  monitors = c("psi", "lam0", "sigma", if (PD) "bD", if (PW) "bW"),
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    conf$removeSamplers(tops)
                    if (identical(settings$sampler, "AF_slice")) conf$addSampler(target = tops, type = "AF_slice")
                    else conf$addSampler(target = tops, type = "RW_block", control = list(adaptInterval = 100))
                  })
  list(chains = lapply(r$samples, ss_rename_draws), t_total = r$diagnostics$t_total,
       t_setup = r$diagnostics$t_setup, engine = "nimble_marginal_local")
}

ss_rename_draws <- function(m) {
  m <- as.matrix(m); cn <- colnames(m)
  cn <- sub("^bD\\[(\\d+)\\]$", "bD\\1", cn); cn <- sub("^bW\\[(\\d+)\\]$", "bW\\1", cn)
  colnames(m) <- cn; m
}
