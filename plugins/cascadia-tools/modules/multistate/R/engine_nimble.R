# multistate: Nimble engine. The latent state sequence is summed out with the
# forward algorithm of nimbleEcology::dDHMMo (one HMM per individual, from the
# occasion after first capture). The first capture is handled exactly as in
# Pradel (2005): the event at first capture is categorical given detection
# (pi_s * delta_s for "seen in s", sum_s pi_s (1 - delta_s) for "state unknown"),
# and the state distribution at first capture given that event, propagated one
# step through the transition matrix, is the HMM's initial distribution.
#
# Occasions are padded with two dummy occasions (T+1, T+2) where nothing can be
# seen and the transition is the identity, so every individual's HMM has
# length >= 3 and all slices keep their dimensions.

ms_define_nimble_functions <- function() {
  if (exists("ms_gamma_seq", envir = globalenv()) && exists("ms_obs_seq", envir = globalenv()) &&
      exists("ms_init_vec", envir = globalenv())) return(invisible())
  suppressPackageStartupMessages({ library(nimble); library(nimbleEcology) })
  assign("ms_gamma_seq", nimble::nimbleFunction(
    run = function(aphi = double(1), bphi = double(1), Xphi = double(2),
                   apsi = double(2), bpsi = double(3), Xpsi = double(2),
                   S = integer(0), T1 = integer(0), Kphi = integer(0), Kpsi = integer(0)) {
      S1 <- S + 1
      G <- nimArray(0, dim = c(S1, S1, T1 + 2))
      e <- nimNumeric(S)
      for (t in 1:T1) {
        for (s in 1:S) {
          lp <- aphi[s]
          if (Kphi > 0) for (k in 1:Kphi) lp <- lp + bphi[k] * Xphi[t, k]
          phi <- 1 / (1 + exp(-lp))
          den <- 0
          for (r in 1:S) {
            eta <- 0
            if (r != s) {
              eta <- apsi[s, r]
              if (Kpsi > 0) for (k in 1:Kpsi) eta <- eta + bpsi[s, r, k] * Xpsi[t, k]
            }
            e[r] <- exp(eta)
            den <- den + e[r]
          }
          for (r in 1:S) G[s, r, t] <- phi * e[r] / den
          G[s, S1, t] <- 1 - phi
        }
        G[S1, S1, t] <- 1
      }
      for (t in (T1 + 1):(T1 + 2)) for (s in 1:S1) G[s, s, t] <- 1
      returnType(double(3))
      return(G)
    }), envir = globalenv())
  assign("ms_obs_seq", nimble::nimbleFunction(
    run = function(ap = double(1), bp = double(1), Xp = double(2), delta = double(1),
                   S = integer(0), T = integer(0), Kp = integer(0), O = integer(0)) {
      S1 <- S + 1
      P <- nimArray(0, dim = c(S1, O, T + 2))
      for (t in 1:T) {
        for (s in 1:S) {
          lp <- ap[s]
          if (Kp > 0) for (k in 1:Kp) lp <- lp + bp[k] * Xp[t, k]
          p <- 1 / (1 + exp(-lp))
          P[s, 1, t] <- 1 - p
          P[s, 1 + s, t] <- p * delta[s]
          if (O > S + 1) P[s, O, t] <- p * (1 - delta[s])
        }
        P[S1, 1, t] <- 1
      }
      for (t in (T + 1):(T + 2)) for (s in 1:S1) P[s, 1, t] <- 1
      returnType(double(3))
      return(P)
    }), envir = globalenv())
  # State distribution at first capture given the first event e1 (1..S known,
  # S+1 unknown), pushed one step through G (S1 x S1).
  assign("ms_init_vec", nimble::nimbleFunction(
    run = function(e1 = double(0), pi = double(1), delta = double(1), G = double(2), S = integer(0)) {
      S1 <- S + 1
      post <- nimNumeric(S)
      if (e1 <= S) {
        post[e1] <- 1
      } else {
        for (s in 1:S) post[s] <- pi[s] * (1 - delta[s])
        post <- post / sum(post)
      }
      out <- nimNumeric(S1)
      for (r in 1:S1) for (s in 1:S) out[r] <- out[r] + post[s] * G[s, r]
      returnType(double(1))
      return(out)
    }), envir = globalenv())
  assign("ms_softmax", nimble::nimbleFunction(
    run = function(a = double(1)) {
      e <- exp(a - max(a))
      returnType(double(1))
      return(e / sum(e))
    }), envir = globalenv())
  invisible()
}

ms_nimble_code <- function() nimble::nimbleCode({
  for (s in 1:S) {
    a_phi[s] ~ dlogis(0, 1)
    a_p[s] ~ dlogis(0, 1)
  }
  for (m in 1:M) {
    apsi[from[m], to[m]] ~ dlogis(0, 1)
  }
  for (s in 1:S) {
    apsi[s, s] <- 0
  }
  # Coefficient vectors carry at least one trailing dummy element fixed at 0 (length >= 2) (with a
  # matching all-zero covariate column): nimble cannot pass a length-1 slice
  # such as bp[1:1] to a nimbleFunction expecting a vector.
  if (Kphi > 0) {
    for (k in 1:Kphi) { bphi[k] ~ dnorm(0, sd = prior_sd) }
  }
  for (k in (Kphi + 1):Kphid) { bphi[k] <- 0 }
  for (k in 1:Kpsid) {
    for (s in 1:S) { bpsi[s, s, k] <- 0 }
  }
  if (Kpsi > 0) {
    for (k in 1:Kpsi) {
      for (m in 1:M) { bpsi[from[m], to[m], k] ~ dnorm(0, sd = prior_sd) }
    }
  }
  for (k in (Kpsi + 1):Kpsid) {
    for (m in 1:M) { bpsi[from[m], to[m], k] <- 0 }
  }
  if (Kp > 0) {
    for (k in 1:Kp) { bp[k] ~ dnorm(0, sd = prior_sd) }
  }
  for (k in (Kp + 1):Kpd) { bp[k] <- 0 }
  if (ME == 1) {
    for (s in 1:S) {
      a_delta[s] ~ dlogis(0, 1)
      delta[s] <- ilogit(a_delta[s])
    }
    a_pi[1] <- 0
    for (s in 2:S) { a_pi[s] ~ dlogis(0, 1) }
    pi[1:S] <- ms_softmax(a_pi[1:S])
    for (s in 1:S) {
      pfirst[s] <- pi[s] * delta[s]
      qfirst[s] <- pi[s] * (1 - delta[s])
    }
    pfirst[S + 1] <- sum(qfirst[1:S])
    nfirst[1:(S + 1)] ~ dmulti(pfirst[1:(S + 1)], Nall)
  } else {
    for (s in 1:S) {
      delta[s] <- 1
      pi[s] <- 1 / S
    }
  }
  Pobs[1:S1, 1:O, 1:Tp] <- ms_obs_seq(a_p[1:S], bp[1:Kpd], Xp[1:T, 1:Kpd], delta[1:S], S, T, Kp, O)
  if (BYIND == 1) {
    for (i in 1:Nh) {
      G[i, 1:S1, 1:S1, 1:Tm] <- ms_gamma_seq(a_phi[1:S], bphi[1:Kphid], Xphi[i, 1:T1, 1:Kphid],
                                             apsi[1:S, 1:S], bpsi[1:S, 1:S, 1:Kpsid], Xpsi[i, 1:T1, 1:Kpsid],
                                             S, T1, Kphi, Kpsi)
      init[i, 1:S1] <- ms_init_vec(e1[i], pi[1:S], delta[1:S], G[i, 1:S1, 1:S1, f[i]], S)
      y[i, 1:L[i]] ~ dDHMMo(init = init[i, 1:S1], probObs = Pobs[1:S1, 1:O, (f[i] + 1):Tp],
                            probTrans = G[i, 1:S1, 1:S1, (f[i] + 1):Tm], len = L[i], checkRowSums = 0)
    }
  } else {
    G[1:S1, 1:S1, 1:Tm] <- ms_gamma_seq(a_phi[1:S], bphi[1:Kphid], Xphi[1:T1, 1:Kphid],
                                        apsi[1:S, 1:S], bpsi[1:S, 1:S, 1:Kpsid], Xpsi[1:T1, 1:Kpsid],
                                        S, T1, Kphi, Kpsi)
    for (i in 1:Nh) {
      init[i, 1:S1] <- ms_init_vec(e1[i], pi[1:S], delta[1:S], G[1:S1, 1:S1, f[i]], S)
      y[i, 1:L[i]] ~ dDHMMo(init = init[i, 1:S1], probObs = Pobs[1:S1, 1:O, (f[i] + 1):Tp],
                            probTrans = G[1:S1, 1:S1, (f[i] + 1):Tm], len = L[i], checkRowSums = 0)
    }
  }
})

# Build constants / data for the model from ms_prepare() output.
ms_nimble_inputs <- function(des, settings) {
  S <- des$S; T <- des$T; ME <- as.integer(des$model == "multievent")
  O <- S + 1 + ME
  h <- which(des$f <= T - 1)                       # individuals with an HMM part
  Nh <- length(h)
  Tp <- T + 2
  # event codes for dDHMMo: 1 = not seen, 1+s = seen in s, S+2 = unknown
  yy <- matrix(1L, Nh, T + 1)
  for (j in seq_len(Nh)) {
    i <- h[j]; f <- des$f[i]
    ev <- des$y[i, (f + 1):T] + 1L
    yy[j, seq_along(ev)] <- ev
  }
  L <- Tp - des$f[h]
  nfirst <- tabulate(des$e1, S + 1)
  # pad covariate arrays to max(K + 1, 2) columns (extra ones all zero, see model code)
  padX <- function(X, K) {
    Kd <- max(K + 1, 2)
    if (length(dim(X)) == 3) {
      out <- array(0, c(dim(X)[1:2], Kd))
      if (K > 0) out[, , 1:K] <- X[, , 1:K, drop = FALSE]
    } else {
      out <- matrix(0, nrow(X), Kd)
      if (K > 0) out[, 1:K] <- X[, 1:K, drop = FALSE]
    }
    out
  }
  Xphi <- padX(des$Xphi, des$Kphi); Xpsi <- padX(des$Xpsi, des$Kpsi); Xp <- padX(des$Xp, des$Kp)
  if (des$gam_by_ind) { Xphi <- Xphi[h, , , drop = FALSE]; Xpsi <- Xpsi[h, , , drop = FALSE] }
  constants <- list(S = S, S1 = S + 1, T = T, T1 = T - 1, Tp = Tp, Tm = Tp - 1, O = O, M = nrow(des$pairs),
                    from = des$pairs[, 1], to = des$pairs[, 2], Nh = Nh, f = des$f[h], L = L,
                    Kphi = des$Kphi, Kpsi = des$Kpsi, Kp = des$Kp, Kphid = max(des$Kphi + 1, 2),
                    Kpsid = max(des$Kpsi + 1, 2), Kpd = max(des$Kp + 1, 2), ME = ME,
                    BYIND = as.integer(des$gam_by_ind), prior_sd = settings$prior_sd, Nall = des$N)
  data <- list(y = yy, e1 = as.numeric(des$e1[h]), Xphi = Xphi, Xpsi = Xpsi, Xp = Xp)
  if (ME == 1) data$nfirst <- nfirst
  list(constants = constants, data = data, h = h)
}

ms_engine_nimble <- function(des, settings) {
  eco_require(c("nimble", "nimbleEcology"))
  ms_define_nimble_functions()
  inp <- ms_nimble_inputs(des, settings)
  cs <- inp$constants
  S <- des$S
  inits <- function() {
    l <- list(a_phi = stats::rnorm(S, 1, 0.3), a_p = stats::rnorm(S, 0, 0.3),
              apsi = matrix(NA_real_, S, S))
    for (m in seq_len(cs$M)) l$apsi[cs$from[m], cs$to[m]] <- stats::rnorm(1, -1.5, 0.3)
    l$bphi <- c(stats::rnorm(cs$Kphi, 0, 0.1), rep(0, cs$Kphid - cs$Kphi))
    b <- array(0, c(S, S, cs$Kpsid))
    if (cs$Kpsi > 0) for (m in seq_len(cs$M)) b[cs$from[m], cs$to[m], 1:cs$Kpsi] <- stats::rnorm(cs$Kpsi, 0, 0.1)
    l$bpsi <- b
    l$bp <- c(stats::rnorm(cs$Kp, 0, 0.1), rep(0, cs$Kpd - cs$Kp))
    if (cs$ME == 1) { l$a_delta <- stats::rnorm(S, 1, 0.3); l$a_pi <- c(NA, stats::rnorm(S - 1, 0, 0.3)) }
    l
  }
  mons <- c("a_phi", "a_p", "apsi", if (cs$Kphi > 0) "bphi", if (cs$Kpsi > 0) "bpsi", if (cs$Kp > 0) "bp",
            if (cs$ME == 1) c("a_delta", "a_pi"))
  key <- paste("multistate", S, des$T, cs$Nh, cs$Kphi, cs$Kpsi, cs$Kp, cs$ME, cs$BYIND,
               paste(cs$f, collapse = "."), settings$prior_sd, settings$sampler, sep = "_")
  key <- paste0("multistate_", digest_key(key))
  r <- nimble_run(ms_nimble_code(), cs, inp$data, inits, mons,
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    nd <- m$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
                    nd <- nd[!grepl("^(y|nfirst)\\[", nd)]
                    if (identical(settings$sampler, "slice")) {
                      conf$removeSamplers(nd); for (n in nd) conf$addSampler(n, "slice")
                    } else if (identical(settings$sampler, "AF_slice")) {
                      conf$removeSamplers(nd); conf$addSampler(nd, "AF_slice")
                    } else if (identical(settings$sampler, "RW_block")) {
                      conf$removeSamplers(nd); conf$addSampler(nd, "RW_block", control = list(adaptInterval = 100))
                    }
                  })
  list(engine = "nimble", chains = r$samples, t_total = r$diagnostics$t_total, t_setup = r$diagnostics$t_setup)
}

# short stable hash for cache keys (first-capture vectors can be long)
digest_key <- function(x) {
  v <- utf8ToInt(x); h <- 0
  for (k in v) h <- (h * 31 + k) %% 2147483647
  sprintf("%s_%d", substr(gsub("[^A-Za-z0-9]", "", x), 1, 40), h)
}
