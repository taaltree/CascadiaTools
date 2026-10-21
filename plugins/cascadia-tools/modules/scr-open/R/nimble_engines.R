# scr-open: Nimble engines. Ported from ../../R/fit_nimble_open.R (the SCR GPU
# benchmark's open-population code) and extended with trap usage per period,
# irregular habitat masks and named density covariates.
#
# States: 1 = not yet entered, 2 = alive, 3 = dead (Royle & Dorazio 2008 /
# Kery & Schaub 2012 restricted-occupancy Jolly-Seber with data augmentation;
# Gardner et al. 2010 spatial version).
#   z[i,1] ~ Bern(gamma_1);  z[i,t] ~ Bern(phi z[i,t-1] + gamma_t * recruitable[i,t])
#   y[i,j,t] ~ Binomial(K_jt, z[i,t] p0 exp(-d(s_i, x_j)^2 / (2 sigma^2)))
#
#  marginal (default): activity centres AND alive states summed out for every
#            augmented individual: the alive-state chain with a forward
#            algorithm, the centre over habitat-mask pixels. One custom
#            likelihood node; RW_block on the 3 + T (+ covariate) parameters.
#            Pixels whose seen-period likelihood is < exp(-30) of the best
#            pixel are skipped (local evaluation; Milleret et al. 2019).
#  hmm     : alive-state chain summed out per individual (dOSCRhmm, as the
#            benchmark's dOpenSCR); activity centres sampled (RW_block).
#  latent  : z[i,t] and s[i,] sampled explicitly (the textbook formulation).

so_bucket <- function(x, b = 50) as.integer(max(b, ceiling(x / b) * b))

so_define_nimble <- function() {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))
  if (exists("dlogitbeta", envir = .GlobalEnv, inherits = FALSE)) return(invisible())
  nf <- nimble::nimbleFunction
  # logit of a Beta(1, a) variable (stick-breaking pieces of a Dirichlet(1,..,1))
  dlogitbeta <- nf(
    run = function(x = double(0), a = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      ld <- log(a) - log1p(exp(-x)) - a * log1p(exp(x))
      if (log) return(ld) else return(exp(ld))
    })
  rlogitbeta <- nf(
    run = function(n = integer(0), a = double(0)) {
      returnType(double(0)); u <- rbeta(1, 1, a); return(log(u / (1 - u)))
    })
  dOSCRmarg <- nf(
    run = function(x = double(0), phi = double(0), gam = double(1), sigma = double(0), lp0 = double(0),
                   lpi = double(1), d2 = double(2), Kjt = double(2), nzj = double(1), nzt = double(1),
                   nzc = double(1), nzstart = double(1), nzend = double(1), seen = double(2),
                   nobs = double(0), M = double(0), G = integer(0), J = integer(0), TT = integer(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      a <- -1 / (2 * sigma^2)
      lg <- matrix(0, G, J)
      base <- matrix(0, G, TT)
      for (g in 1:G) {
        for (j in 1:J) {
          lp <- lp0 + a * d2[g, j]
          l1 <- log1p(-exp(lp))
          lg[g, j] <- lp - l1
          for (t in 1:TT) base[g, t] <- base[g, t] + Kjt[j, t] * l1
        }
      }
      pg <- exp(lpi)
      # never-detected history, per pixel
      L0 <- 0
      for (g in 1:G) {
        a1 <- 1 - gam[1]; a2 <- gam[1] * exp(base[g, 1]); a3 <- 0
        for (t in 2:TT) {
          n1 <- a1 * (1 - gam[t])
          n2 <- (a1 * gam[t] + a2 * phi) * exp(base[g, t])
          n3 <- a2 * (1 - phi) + a3
          a1 <- n1; a2 <- n2; a3 <- n3
        }
        L0 <- L0 + pg[g] * (a1 + a2 + a3)
      }
      if (L0 <= 0) { if (log) return(-Inf) else return(0) }
      tot <- (M - nobs) * log(L0)
      e <- matrix(0, G, TT)
      mt <- numeric(TT)
      sc <- numeric(G)
      nob <- nobs
      if (nob > 0) {
        for (i in 1:nob) {
          e <- base
          for (k in nzstart[i]:nzend[i]) {
            j <- nzj[k]; tt <- nzt[k]; cc <- nzc[k]
            for (g in 1:G) e[g, tt] <- e[g, tt] + cc * lg[g, j]
          }
          shift <- 0
          for (t in 1:TT) {
            mt[t] <- 0
            if (seen[i, t] > 0) { mt[t] <- max(e[1:G, t]); shift <- shift + mt[t] }
          }
          for (g in 1:G) {
            sc[g] <- 0
            for (t in 1:TT) if (seen[i, t] > 0) sc[g] <- sc[g] + e[g, t] - mt[t]
          }
          scmax <- max(sc)
          Li <- 0
          for (g in 1:G) {
            if (sc[g] > scmax - 30) {
              u <- 1 - seen[i, 1]
              a1 <- (1 - gam[1]) * u; a2 <- gam[1] * exp(e[g, 1] - mt[1]); a3 <- 0
              for (t in 2:TT) {
                u <- 1 - seen[i, t]
                et <- exp(e[g, t] - mt[t])
                n1 <- a1 * (1 - gam[t]) * u
                n2 <- (a1 * gam[t] + a2 * phi) * et
                n3 <- (a2 * (1 - phi) + a3) * u
                a1 <- n1; a2 <- n2; a3 <- n3
              }
              Li <- Li + pg[g] * (a1 + a2 + a3)
            }
          }
          if (Li <= 0) { if (log) return(-Inf) else return(0) }
          tot <- tot + log(Li) + shift
        }
      }
      if (log) return(tot) else return(exp(tot))
    })
  rOSCRmarg <- nf(
    run = function(n = integer(0), phi = double(0), gam = double(1), sigma = double(0), lp0 = double(0),
                   lpi = double(1), d2 = double(2), Kjt = double(2), nzj = double(1), nzt = double(1),
                   nzc = double(1), nzstart = double(1), nzend = double(1), seen = double(2),
                   nobs = double(0), M = double(0), G = integer(0), J = integer(0), TT = integer(0)) {
      returnType(double(0)); return(1)
    })
  # detection of one individual in one period (skips z = 0)
  dOSCRbin <- nf(
    run = function(x = double(1), s = double(1), lp0 = double(0), sigma = double(0),
                   traps = double(2), Kj = double(1), J = integer(0), z = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      if (z == 0) {
        if (sum(x) > 0) { if (log) return(-Inf) else return(0) }
        if (log) return(0) else return(1)
      }
      ll <- 0
      a <- -1 / (2 * sigma^2)
      for (j in 1:J) if (Kj[j] > 0) {
        p <- exp(lp0 + a * ((traps[j, 1] - s[1])^2 + (traps[j, 2] - s[2])^2))
        ll <- ll + dbinom(x[j], Kj[j], p, log = TRUE)
      }
      if (log) return(ll) else return(exp(ll))
    })
  rOSCRbin <- nf(
    run = function(n = integer(0), s = double(1), lp0 = double(0), sigma = double(0),
                   traps = double(2), Kj = double(1), J = integer(0), z = double(0)) {
      returnType(double(1)); out <- numeric(J); return(out)
    })
  # whole history of one individual, alive-state chain summed out
  dOSCRhmm <- nf(
    run = function(x = double(2), s = double(1), lp0 = double(0), sigma = double(0),
                   phi = double(0), gam = double(1), traps = double(2), Kjt = double(2),
                   J = integer(0), TT = integer(0), log = integer(0, default = 0)) {
      returnType(double(0))
      a <- -1 / (2 * sigma^2)
      lp <- numeric(J); l1m <- numeric(J)
      for (j in 1:J) {
        lpj <- lp0 + a * ((traps[j, 1] - s[1])^2 + (traps[j, 2] - s[2])^2)
        lp[j] <- lpj; l1m[j] <- log1p(-exp(lpj))
      }
      la <- numeric(TT); seen <- numeric(TT)
      for (t in 1:TT) {
        ll <- 0; tot <- 0
        for (j in 1:J) {
          ll <- ll + x[j, t] * lp[j] + (Kjt[j, t] - x[j, t]) * l1m[j]
          tot <- tot + x[j, t]
        }
        la[t] <- ll; seen[t] <- tot
      }
      neg <- -1e30
      zt <- neg; if (seen[1] == 0) zt <- 0
      f1 <- log(1 - gam[1]) + zt
      f2 <- log(gam[1]) + la[1]
      f3 <- neg
      for (t in 2:TT) {
        zt <- neg; if (seen[t] == 0) zt <- 0
        g1 <- f1 + log(1 - gam[t]) + zt
        mx <- f1 + log(gam[t]); m2 <- f2 + log(phi)
        if (m2 > mx) { tmp <- mx; mx <- m2; m2 <- tmp }
        g2 <- mx + log(1 + exp(m2 - mx)) + la[t]
        mx3 <- f2 + log(1 - phi); m3 <- f3
        if (m3 > mx3) { tmp <- mx3; mx3 <- m3; m3 <- tmp }
        g3 <- mx3 + log(1 + exp(m3 - mx3)) + zt
        f1 <- g1; f2 <- g2; f3 <- g3
      }
      mx <- f1
      if (f2 > mx) mx <- f2
      if (f3 > mx) mx <- f3
      out <- mx + log(exp(f1 - mx) + exp(f2 - mx) + exp(f3 - mx))
      if (log) return(out) else return(exp(out))
    })
  rOSCRhmm <- nf(
    run = function(n = integer(0), s = double(1), lp0 = double(0), sigma = double(0),
                   phi = double(0), gam = double(1), traps = double(2), Kjt = double(2),
                   J = integer(0), TT = integer(0)) {
      returnType(double(2)); out <- matrix(0, J, TT); return(out)
    })
  # activity-centre point process on the bounding rectangle of the mask
  dOSCRhab <- nf(
    run = function(x = double(1), lpi = double(1), hab = double(1), pix = double(0), nx = integer(0),
                   ny = integer(0), log = integer(0, default = 0)) {
      returnType(double(0))
      ix <- floor(x[1] / pix); iy <- floor(x[2] / pix)
      if (ix < 0 | iy < 0 | ix >= nx | iy >= ny) { if (log) return(-Inf) else return(0) }
      g <- ix + iy * nx + 1
      if (hab[g] == 0) { if (log) return(-Inf) else return(0) }
      ld <- lpi[g] - 2 * log(pix)
      if (log) return(ld) else return(exp(ld))
    })
  rOSCRhab <- nf(
    run = function(n = integer(0), lpi = double(1), hab = double(1), pix = double(0), nx = integer(0),
                   ny = integer(0)) {
      returnType(double(1))
      g <- rcat(1, exp(lpi) * hab) - 1
      out <- numeric(2)
      out[1] <- (g %% nx + runif(1)) * pix
      out[2] <- (floor(g / nx) + runif(1)) * pix
      return(out)
    })
  for (nm in c("dlogitbeta", "rlogitbeta", "dOSCRmarg", "rOSCRmarg", "dOSCRbin", "rOSCRbin", "dOSCRhmm", "rOSCRhmm", "dOSCRhab", "rOSCRhab"))
    assign(nm, get(nm), envir = .GlobalEnv)
  invisible()
}

# Probabilities are sampled on the logit scale with a standard logistic prior,
# which is exactly the Uniform(0, 1) prior on the probability scale; this
# lets the block random walk move freely near 0/1.
#
# Entry: superpopulation parameterisation (Kery & Schaub 2012 sec. 10.3.3;
# Schwarz & Arnason 1996): psi ~ Uniform(0,1) is the share of the M augmented
# rows that are ever alive (so Nsuper has a uniform prior on 0..M), and
# b ~ Dirichlet(1,..,1) the entry probabilities of those animals (built by
# stick-breaking, v_t ~ Beta(1, T - t)). The per-period entry hazards gamma_t
# used by the likelihood are derived: gamma_1 = psi b_1,
# gamma_t = psi b_t / (1 - psi sum_{k<t} b_k). (Independent Uniform priors on
# gamma_1..T, used before 2026-10-20, imply a prior on Nsuper/M with median
# 0.97 at T = 4 and biased abundance upward; see VALIDATION.md.)
so_code_entry <- function(logit = TRUE) {
  if (logit) quote({
    lpsi ~ dlogis(0, 1)
    psi <- ilogit(lpsi)
    for (t in 1:(TT - 1)) {
      lv[t] ~ dlogitbeta(a = TT - t)
      v[t] <- ilogit(lv[t])
    }
  }) else quote({
    psi ~ dunif(0, 1)
    for (t in 1:(TT - 1)) {
      v[t] ~ dbeta(1, TT - t)
    }
  })
}
so_code_entry_derived <- quote({
  bb[1] <- v[1]
  rem[1] <- 1 - v[1]
  if (TT > 2) {
    for (t in 2:(TT - 1)) {
      bb[t] <- v[t] * rem[t - 1]
      rem[t] <- rem[t - 1] * (1 - v[t])
    }
  }
  bb[TT] <- rem[TT - 1]
  gam[1] <- psi * bb[1]
  cb[1] <- bb[1]
  for (t in 2:TT) {
    gam[t] <- psi * bb[t] / (1 - psi * cb[t - 1])
    cb[t] <- cb[t - 1] + bb[t]
  }
})
# splice the entry block into a model body
so_with_entry <- function(code, logit) {
  body <- as.list(code)[-1]
  i <- which(vapply(body, function(e) identical(e, quote(ENTRY)), logical(1)))
  ent <- c(as.list(so_code_entry(logit))[-1], as.list(so_code_entry_derived)[-1])
  as.call(c(as.name("{"), body[seq_len(i - 1)], ent, body[-seq_len(i)]))
}

so_code_marginal <- function() so_with_entry(quote({
  lp0l ~ dlogis(0, 1)
  p0 <- ilogit(lp0l)
  sigma ~ dunif(0, sigma_max)
  lphi ~ dlogis(0, 1)
  phi <- ilogit(lphi)
  ENTRY
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
  lp0 <- log(p0)
  ones ~ dOSCRmarg(phi = phi, gam = gam[1:TT], sigma = sigma, lp0 = lp0, lpi = lpi[1:G], d2 = d2[1:G, 1:J],
                   Kjt = Kjt[1:J, 1:TT], nzj = nzj[1:NNZMAX], nzt = nzt[1:NNZMAX], nzc = nzc[1:NNZMAX],
                   nzstart = nzstart[1:NMAX], nzend = nzend[1:NMAX], seen = seen[1:NMAX, 1:TT],
                   nobs = nobs, M = M, G = G, J = J, TT = TT)
}), logit = TRUE)

# Shared top of the hmm / latent codes (activity centres sampled on the
# bounding rectangle of the mask).
so_code_latent <- function() so_with_entry(quote({
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, sigma_max)
  phi ~ dunif(0, 1)
  ENTRY
  if (HASD) {
    for (k in 1:PD) bD[k] ~ dnorm(0, sd = 2)
    if (PD == 1) {
      lw[1:GR] <- XR[1:GR, 1] * bD[1] + offR[1:GR]
    } else {
      lw[1:GR] <- (XR[1:GR, 1:PD] %*% bD[1:PD])[1:GR, 1] + offR[1:GR]
    }
    lpi[1:GR] <- lw[1:GR] - log(sum(exp(lw[1:GR])))
  } else {
    lpi[1:GR] <- lpi0[1:GR]
  }
  lp0 <- log(p0)
  for (i in 1:M) {
    s[i, 1:2] ~ dOSCRhab(lpi = lpi[1:GR], hab = hab[1:GR], pix = pix, nx = nx, ny = ny)
    z[i, 1] ~ dbern(gam[1])
    rec[i, 1] <- 1
    ent[i, 1] <- z[i, 1]
    y[i, 1:J, 1] ~ dOSCRbin(s = s[i, 1:2], lp0 = lp0, sigma = sigma, traps = traps[1:J, 1:2],
                            Kj = Kjt[1:J, 1], J = J, z = z[i, 1])
    for (t in 2:TT) {
      rec[i, t] <- rec[i, t - 1] * (1 - z[i, t - 1])
      mu[i, t] <- phi * z[i, t - 1] + gam[t] * rec[i, t]
      z[i, t] ~ dbern(mu[i, t])
      ent[i, t] <- z[i, t] * rec[i, t]
      y[i, 1:J, t] ~ dOSCRbin(s = s[i, 1:2], lp0 = lp0, sigma = sigma, traps = traps[1:J, 1:2],
                              Kj = Kjt[1:J, t], J = J, z = z[i, t])
    }
    ever[i] <- 1 - prod(1 - z[i, 1:TT])
  }
  for (t in 1:TT) {
    N[t] <- sum(z[1:M, t])
    B[t] <- sum(ent[1:M, t])
  }
  Nsuper <- sum(ever[1:M])
}), logit = FALSE)

so_code_hmm <- function() so_with_entry(quote({
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, sigma_max)
  phi ~ dunif(0, 1)
  ENTRY
  if (HASD) {
    for (k in 1:PD) bD[k] ~ dnorm(0, sd = 2)
    if (PD == 1) {
      lw[1:GR] <- XR[1:GR, 1] * bD[1] + offR[1:GR]
    } else {
      lw[1:GR] <- (XR[1:GR, 1:PD] %*% bD[1:PD])[1:GR, 1] + offR[1:GR]
    }
    lpi[1:GR] <- lw[1:GR] - log(sum(exp(lw[1:GR])))
  } else {
    lpi[1:GR] <- lpi0[1:GR]
  }
  lp0 <- log(p0)
  for (i in 1:M) {
    s[i, 1:2] ~ dOSCRhab(lpi = lpi[1:GR], hab = hab[1:GR], pix = pix, nx = nx, ny = ny)
    y[i, 1:J, 1:TT] ~ dOSCRhmm(s = s[i, 1:2], lp0 = lp0, sigma = sigma, phi = phi, gam = gam[1:TT],
                               traps = traps[1:J, 1:2], Kjt = Kjt[1:J, 1:TT], J = J, TT = TT)
  }
}), logit = FALSE)

so_tops <- function(prep) c("lp0l", "sigma", "lphi", "lpsi", paste0("lv[", seq_len(prep$TT - 1), "]"),
                            if (ncol(prep$XD)) paste0("bD[", seq_len(ncol(prep$XD)), "]"))

# initial entry hazards gamma_1..T -> superpopulation parameterisation (psi, stick-breaking v)
so_gam_to_entry <- function(gam) {
  TT <- length(gam)
  ent <- 1 - cumprod(1 - gam)                   # P(entered by t)
  psi <- ent[TT]; cb <- ent / psi
  b <- diff(c(0, cb))
  v <- b[-TT] / (1 - c(0, cb[-TT])[-TT])
  list(psi = psi, v = pmin(pmax(v, 1e-4), 1 - 1e-4))
}

so_inits_top <- function(prep, M) {
  sg <- if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else prep$trap_spacing / prep$cscale
  gam <- numeric(prep$TT)
  gam[1] <- min(0.9, 1.6 * max(prep$n_by_period[1], 1) / M)
  for (t in 2:prep$TT) gam[t] <- min(0.5, 0.6 * max(prep$n_by_period[t], 1) / M)
  ini <- list(p0 = stats::runif(1, 0.1, 0.3), sigma = min(prep$sigma_max / 2, sg * stats::runif(1, 0.8, 1.25)),
              phi = stats::runif(1, 0.5, 0.85), gam = gam * stats::runif(prep$TT, 0.8, 1.2))
  if (ncol(prep$XD)) ini$bD <- stats::rnorm(ncol(prep$XD), 0, 0.2)
  ini
}

so_engine_nimble_marginal <- function(prep, settings) {
  so_define_nimble()
  PD <- ncol(prep$XD); TT <- prep$TT
  NNZMAX <- so_bucket(max(length(prep$nzc), 1), 100); NMAX <- so_bucket(max(prep$n, 1), 50)
  pad <- function(v, L) c(v, rep(1, L - length(v)))
  seen <- matrix(0, NMAX, TT); if (prep$n) seen[seq_len(prep$n), ] <- prep$seen * 1
  constants <- list(G = prep$G, J = prep$J, TT = TT, NNZMAX = NNZMAX, NMAX = NMAX, HASD = PD > 0, PD = max(PD, 1))
  data <- list(ones = 1, d2 = prep$d2, Kjt = prep$usage * 1,
               nzj = pad(prep$nzj, NNZMAX), nzt = pad(prep$nzt, NNZMAX), nzc = pad(prep$nzc, NNZMAX),
               nzstart = pad(prep$nzstart, NMAX), nzend = pad(prep$nzend, NMAX), seen = seen,
               nobs = prep$n, M = prep$M, sigma_max = prep$sigma_max,
               lpi0 = rep(-log(prep$G), prep$G), XD = if (PD) prep$XD else matrix(0, prep$G, 1))
  tops <- so_tops(prep)
  sg0 <- if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else prep$trap_spacing / prep$cscale
  pc <- diag(c(0.1, 0.03 * sg0, 0.2, rep(0.3, TT), rep(0.1, PD))^2, length(tops))
  key <- paste("scro_marg3", if (is.null(settings$top_sampler)) "RW_block" else settings$top_sampler, prep$G, prep$J, TT, NNZMAX, NMAX, PD, sep = "_")
  r <- nimble_run(so_code_marginal(), constants, data,
                  inits = function() {
                    i0 <- so_inits_top(prep, prep$M)
                    out <- list(lp0l = stats::qlogis(i0$p0), sigma = i0$sigma, lphi = stats::qlogis(i0$phi),
                                lpsi = stats::qlogis(so_gam_to_entry(i0$gam)$psi),
                                lv = stats::qlogis(so_gam_to_entry(i0$gam)$v))
                    if (PD) out$bD <- i0$bD
                    out
                  },
                  monitors = c("p0", "sigma", "phi", "gam", if (PD) "bD"),
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    conf$removeSamplers(tops)
                    if (identical(settings$top_sampler, "AF_slice")) {
                      conf$addSampler(target = tops, type = "AF_slice")
                    } else {
                      conf$addSampler(target = tops, type = "RW_block",
                                      control = list(adaptInterval = 100, adaptFactorExponent = 0.6, propCov = pc))
                    }
                  })
  list(chains = lapply(r$samples, so_rename_draws), t_total = r$diagnostics$t_total,
       t_setup = r$diagnostics$t_setup, engine = "nimble marginal")
}

so_rect <- function(prep) {
  ix <- floor(prep$mxy[, 1] / prep$pix + 1e-6); iy <- floor(prep$mxy[, 2] / prep$pix + 1e-6)
  nx <- max(ix) + 1L; ny <- max(iy) + 1L
  list(nx = as.integer(nx), ny = as.integer(ny), GR = as.integer(nx * ny), g = ix + iy * nx + 1)
}

# hmm and latent engines share data/inits; approach picks the code
so_engine_nimble_s <- function(prep, settings, approach = c("latent", "hmm")) {
  approach <- match.arg(approach)
  so_define_nimble()
  PD <- ncol(prep$XD); TT <- prep$TT; J <- prep$J
  rc <- so_rect(prep)
  M <- prep$M
  hab <- numeric(rc$GR); hab[rc$g] <- 1
  XR <- matrix(0, rc$GR, max(PD, 1)); if (PD) XR[rc$g, ] <- prep$XD
  offR <- rep(-1e10, rc$GR); offR[rc$g] <- 0
  lpi0 <- rep(-1e10, rc$GR); lpi0[rc$g] <- -log(prep$G)
  yaug <- array(0, c(M, J, TT)); if (prep$n) yaug[seq_len(prep$n), , ] <- prep$y
  constants <- list(M = M, J = J, TT = TT, GR = rc$GR, nx = rc$nx, ny = rc$ny, HASD = PD > 0, PD = max(PD, 1))
  data <- list(y = yaug, traps = prep$txy, Kjt = prep$usage * 1, XR = XR, offR = offR, lpi0 = lpi0,
               hab = hab, pix = prep$pix, sigma_max = prep$sigma_max)
  inits <- function() {
    ini <- so_inits_top(prep, M)
    s0 <- prep$mxy[sample.int(prep$G, M, replace = TRUE), , drop = FALSE]
    for (i in seq_len(prep$n)) {
      w <- rowSums(prep$y[i, , , drop = FALSE][1, , , drop = TRUE])
      c0 <- colSums(prep$txy * w) / sum(w)
      gi <- which.min((prep$mxy[, 1] - c0[1])^2 + (prep$mxy[, 2] - c0[2])^2)
      s0[i, ] <- prep$mxy[gi, ]
    }
    ini$s <- s0 + matrix(stats::runif(2 * M, -0.45, 0.45) * prep$pix, M, 2)
    if (approach == "latent") {
      # alive from first to last detected period (anything else starts at -Inf)
      z0 <- matrix(0, M, TT)
      for (i in seq_len(prep$n)) { per <- which(prep$seen[i, ]); z0[i, min(per):max(per)] <- 1 }
      ini$z <- z0
      ini$gam[1] <- max(ini$gam[1], (sum(z0[, 1]) + 1) / M)
    }
    en <- so_gam_to_entry(ini$gam); ini$psi <- en$psi; ini$v <- en$v; ini$gam <- NULL
    ini
  }
  sc <- 0.5 * (if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else prep$trap_spacing / prep$cscale)
  key <- paste("scro3", approach, M, J, TT, rc$GR, PD, sep = "_")
  mons <- c("p0", "sigma", "phi", "gam", if (PD) "bD", if (approach == "latent") c("N", "B", "Nsuper"))
  r <- nimble_run(if (approach == "latent") so_code_latent() else so_code_hmm(), constants, data, inits = inits,
                  monitors = mons, niter = settings$niter_latent, nburnin = settings$nburnin_latent,
                  nchains = settings$nchains, thin = settings$thin_latent, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    conf$removeSamplers("s")
                    for (i in seq_len(M)) conf$addSampler(target = sprintf("s[%d, 1:2]", i), type = "RW_block",
                                                          control = list(adaptScaleOnly = TRUE, scale = sc))
                  })
  list(chains = lapply(r$samples, so_rename_draws), t_total = r$diagnostics$t_total,
       t_setup = r$diagnostics$t_setup, engine = paste("nimble", approach))
}

# Nimble names -> internal names p0, sigma, phi, gam1.., bD1.., N1.., B1.., Nsuper
so_rename_draws <- function(m) {
  m <- as.matrix(m)
  cn <- colnames(m)
  cn <- sub("^(bD|gam|N|B)\\[(\\d+)\\]$", "\\1\\2", cn)
  colnames(m) <- cn
  m
}
