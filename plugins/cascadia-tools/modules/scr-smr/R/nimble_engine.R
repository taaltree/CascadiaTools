# scr-smr: Nimble engine (Bayesian, data augmentation, latent activity
# centres), following Chandler & Royle (2013), Sollmann et al. (2013b),
# Whittington et al. (2018) and Jimenez et al. (2019).
#
# Rows of the augmented population (M rows):
#   1..m        marked animals (identities known; z = 1)
#   m+1..MK     padding rows reserved for marked animals (unmarked/augmented)
#   MK+1..M     augmented rows: unmarked animals or "not in the population"
# For each row: z[i] ~ Bern(pz[i]), s[i] ~ pixelwise point process.
#   SMR : pz = 1 for marked rows (they are a given random sample; no psi term),
#         psi for the rest, so N = m + sum(unmarked z) with N_u ~ Bin(M - m, psi).
#   gSMR: pz = psi for every row; the marking captures x[i, l] ~ Bin(Km_l,
#         pm0 exp(-d^2 / 2 sigma^2)) are modelled for every row (zero rows for
#         animals never caught), so marked animals are NOT a random sample.
# Resighting (camera j, K_j operating occasions), lam[i, j] = lambda0_j
# exp(-d_ij^2 / 2 sigma^2) (Poisson rate) or p (Bernoulli):
#   identified sightings of marked i:  y[i, j] ~ Pois(K_j delta lam_ij)
#                                      (Bernoulli: Bin(K_j, delta p_ij))
#   unmarked sightings at camera j:    nU[j] ~ Pois(K_j sum_i z_i (1 - mk_i) lam_ij)
#   marked, not identified:            nTm[j] ~ Pois(K_j (1 - delta) sum_i mk_i lam_ij)
# (exact for the Poisson model by Poisson thinning/superposition; for the
# Bernoulli model the two pooled counts use the Poisson approximation of
# Chandler & Royle 2013.) Telemetry fixes of marked animals ~ N(s_i, sigma^2 I),
# entered through their sufficient statistics.

smr_bucket <- function(x, b = 50) as.integer(max(b, ceiling(x / b) * b))

smr_define_nimble <- function() {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))
  nimble::nimbleOptions(verbose = FALSE, MCMCprogressBar = FALSE)
  if (exists("dSMRPool", envir = .GlobalEnv, inherits = FALSE)) return(invisible())
  nf <- nimble::nimbleFunction
  # Activity-centre prior. Rows in the population (z = 1) follow the
  # pixelwise inhomogeneous point process pi_g / pixel area; rows not in the
  # population (z = 0) get a uniform pseudo-prior over habitat (lunif =
  # -log(habitat area)). s does not enter the likelihood when z = 0, so any
  # proper pseudo-prior leaves the posterior of everything else unchanged,
  # and the density coefficients are then informed only by real animals
  # (much better mixing than letting all M rows follow pi).
  dSMRHab <- nf(
    run = function(x = double(1), lpi = double(1), hab = double(1), pix = double(0), nx = integer(0),
                   ny = integer(0), z = double(0), lunif = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      ix <- floor(x[1] / pix); iy <- floor(x[2] / pix)
      if (ix < 0 | iy < 0 | ix >= nx | iy >= ny) { if (log) return(-Inf) else return(0) }
      g <- ix + iy * nx + 1
      if (hab[g] == 0) { if (log) return(-Inf) else return(0) }
      if (z == 0) ld <- lunif else ld <- lpi[g] - 2 * log(pix)
      if (log) return(ld) else return(exp(ld))
    })
  rSMRHab <- nf(
    run = function(n = integer(0), lpi = double(1), hab = double(1), pix = double(0), nx = integer(0), ny = integer(0),
                   z = double(0), lunif = double(0)) {
      returnType(double(1))
      if (z == 0) g <- rcat(1, hab) - 1 else g <- rcat(1, exp(lpi) * hab) - 1
      out <- numeric(2)
      out[1] <- (g %% nx + runif(1)) * pix
      out[2] <- (floor(g / nx) + runif(1)) * pix
      return(out)
    })
  # identified resightings of one marked animal (rows with mk = 0 carry zeros)
  dSMRSight <- nf(
    run = function(x = double(1), lam = double(1), Kj = double(1), delta = double(0), mk = double(0),
                   z = double(0), pois = double(0), J = integer(0), log = integer(0, default = 0)) {
      returnType(double(0))
      if (mk == 0) { if (log) return(0) else return(1) }
      if (z == 0) { if (log) return(-Inf) else return(0) }
      ll <- 0
      for (j in 1:J) if (Kj[j] > 0) {
        if (pois == 1) ll <- ll + dpois(x[j], Kj[j] * delta * lam[j], log = TRUE)
        else ll <- ll + dbinom(x[j], Kj[j], delta * lam[j], log = TRUE)
      }
      if (log) return(ll) else return(exp(ll))
    })
  rSMRSight <- nf(
    run = function(n = integer(0), lam = double(1), Kj = double(1), delta = double(0), mk = double(0),
                   z = double(0), pois = double(0), J = integer(0)) {
      returnType(double(1)); out <- numeric(J); return(out)
    })
  # pooled Poisson counts per camera. type 0: unmarked (w = z (1 - mk));
  # type 1: marked but unidentified (w = mk (1 - delta)).
  dSMRPool <- nf(
    run = function(x = double(1), lam = double(2), z = double(1), mk = double(1), Kj = double(1),
                   delta = double(0), type = double(0), R = integer(0), J = integer(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      Lam <- numeric(J)
      for (i in 1:R) {
        if (type == 0) w <- z[i] * (1 - mk[i]) else w <- mk[i] * (1 - delta)
        if (w > 0) for (j in 1:J) Lam[j] <- Lam[j] + w * lam[i, j]
      }
      ll <- 0
      for (j in 1:J) if (Kj[j] > 0) ll <- ll + dpois(x[j], Kj[j] * Lam[j], log = TRUE)
      if (log) return(ll) else return(exp(ll))
    })
  rSMRPool <- nf(
    run = function(n = integer(0), lam = double(2), z = double(1), mk = double(1), Kj = double(1),
                   delta = double(0), type = double(0), R = integer(0), J = integer(0)) {
      returnType(double(1)); out <- numeric(J); return(out)
    })
  # marking-phase captures (binomial, proximity-type traps)
  dSMRMark <- nf(
    run = function(x = double(1), s = double(1), lp0 = double(0), sigma = double(0), traps = double(2),
                   Kj = double(1), L = integer(0), z = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      if (z == 0) {
        if (sum(x) > 0) { if (log) return(-Inf) else return(0) }
        if (log) return(0) else return(1)
      }
      ll <- 0
      a <- -1 / (2 * sigma^2)
      for (l in 1:L) if (Kj[l] > 0) {
        p <- exp(lp0 + a * ((traps[l, 1] - s[1])^2 + (traps[l, 2] - s[2])^2))
        ll <- ll + dbinom(x[l], Kj[l], p, log = TRUE)
      }
      if (log) return(ll) else return(exp(ll))
    })
  rSMRMark <- nf(
    run = function(n = integer(0), s = double(1), lp0 = double(0), sigma = double(0), traps = double(2),
                   Kj = double(1), L = integer(0), z = double(0)) {
      returnType(double(1)); out <- numeric(L); return(out)
    })
  # telemetry sufficient statistics x = (n fixes, mean x, mean y, SS about mean)
  dSMRTel <- nf(
    run = function(x = double(1), s = double(1), sigma = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      if (x[1] == 0) { if (log) return(0) else return(1) }
      ll <- -x[1] * log(2 * pi * sigma^2) - (x[4] + x[1] * ((x[2] - s[1])^2 + (x[3] - s[2])^2)) / (2 * sigma^2)
      if (log) return(ll) else return(exp(ll))
    })
  rSMRTel <- nf(
    run = function(n = integer(0), s = double(1), sigma = double(0)) {
      returnType(double(1)); out <- numeric(4); return(out)
    })
  # Joint sweep over all activity centres s[i] and inclusion indicators z[i]
  # with cached camera-level pooled rates, so one sweep costs O(M J) instead
  # of O(M^2 J) for the graph-based samplers (the pooled unmarked-count node
  # depends on every row). Each s[i] gets a random-walk Metropolis step
  # (adaptive scale for z = 1 rows, wider fixed multiple for z = 0 rows);
  # each augmented, unmarked row gets a Metropolised flip of z[i].
  sampler_SMRlatent <- nf(
    name = "sampler_SMRlatent",
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      M <- control$M; MK <- control$MK; J <- control$J; L <- control$L
      nx <- control$nx; ny <- control$ny
      pois <- control$pois; hasTm <- control$hasTm; hasMark <- control$hasMark
      scale <- control$scale; scale0 <- control$scale
      calcNodes <- model$getDependencies(c("s", "z"))
      lamC <- matrix(0, M, J)
      lamNew <- numeric(J); LU <- numeric(J); LT <- numeric(J); LUn <- numeric(J); LTn <- numeric(J)
      nAcc <- 0; nTry <- 0; nRuns <- 0
    },
    run = function() {
      g <- 1; go <- 1; gs <- 1; ix <- 0; iy <- 0; p1 <- 0; p2 <- 0; lr <- 0; pm <- 0; pn <- 0; po <- 0
      lm <- 0; s1 <- 0; s2 <- 0; sc <- 0; zi <- 0; mki <- 0; dz <- 0; ok <- 0
      pix <- 0; sigma <- 0; psi <- 0; dl <- 0; lpm <- 0; a <- 0
      pix <- model$pix[1]
      sigma <- model$sigma[1]; psi <- model$psi[1]; dl <- model$dID[1]; lpm <- model$lpm[1]
      a <- -1 / (2 * sigma^2)
      lamC <<- model$lam
      LU <<- numeric(J); LT <<- numeric(J)
      for (i in 1:M) {
        if (model$mk[i] == 1) { for (j in 1:J) LT[j] <<- LT[j] + (1 - dl) * lamC[i, j] }
        else if (model$z[i] == 1) { for (j in 1:J) LU[j] <<- LU[j] + lamC[i, j] }
      }
      dU <- 0; dT <- 0
      for (i in 1:M) {
        zi <- model$z[i]; mki <- model$mk[i]
        # ---- activity centre
        s1 <- model$s[i, 1]; s2 <- model$s[i, 2]
        sc <- scale
        if (zi == 0) sc <- 3 * scale
        p1 <- s1 + rnorm(1, 0, sc); p2 <- s2 + rnorm(1, 0, sc)
        ix <- floor(p1 / pix); iy <- floor(p2 / pix)
        ok <- 1
        if (ix < 0 | iy < 0 | ix >= nx | iy >= ny) ok <- 0
        if (ok == 1) {
          g <- ix + iy * nx + 1
          if (model$hab[g] == 0) ok <- 0
        }
        if (zi == 1) nTry <<- nTry + 1
        if (ok == 1) {
          go <- floor(s1 / pix) + floor(s2 / pix) * nx + 1
          lr <- model$lpi[g] - model$lpi[go]
          if (zi == 0) lr <- lr * 0
          if (zi == 1) {
            for (j in 1:J) lamNew[j] <<- exp(model$lr0[j] + a * ((model$traps[j, 1] - p1)^2 + (model$traps[j, 2] - p2)^2))
            if (mki == 1) {
              for (j in 1:J) if (model$Kj[j] > 0) {
                if (pois == 1) lr <- lr + dpois(model$y[i, j], model$Kj[j] * dl * lamNew[j], log = TRUE) -
                                     dpois(model$y[i, j], model$Kj[j] * dl * lamC[i, j], log = TRUE)
                else lr <- lr + dbinom(model$y[i, j], model$Kj[j], dl * lamNew[j], log = TRUE) -
                              dbinom(model$y[i, j], model$Kj[j], dl * lamC[i, j], log = TRUE)
              }
              if (model$tel[i, 1] > 0) {
                lr <- lr - model$tel[i, 1] * ((model$tel[i, 2] - p1)^2 + (model$tel[i, 3] - p2)^2 -
                                               (model$tel[i, 2] - s1)^2 - (model$tel[i, 3] - s2)^2) / (2 * sigma^2)
              }
            }
            if (hasMark == 1) {
              for (l in 1:L) if (model$Kmj[l] > 0) {
                pn <- exp(lpm + a * ((model$mtraps[l, 1] - p1)^2 + (model$mtraps[l, 2] - p2)^2))
                po <- exp(lpm + a * ((model$mtraps[l, 1] - s1)^2 + (model$mtraps[l, 2] - s2)^2))
                lr <- lr + dbinom(model$xm[i, l], model$Kmj[l], pn, log = TRUE) - dbinom(model$xm[i, l], model$Kmj[l], po, log = TRUE)
              }
            }
            # change in the pooled Poisson log-likelihoods, without lgamma:
            # sum_j n_j (log L'_j - log L_j) - K_j (L'_j - L_j)
            dU <- 0; dT <- 0
            if (mki == 0) {
              for (j in 1:J) {
                LUn[j] <<- LU[j] - lamC[i, j] + lamNew[j]
                if (model$Kj[j] > 0) {
                  dU <- dU - model$Kj[j] * (LUn[j] - LU[j])
                  if (model$nU[j] > 0 & LUn[j] != LU[j]) dU <- dU + model$nU[j] * (log(LUn[j]) - log(LU[j]))
                }
              }
            } else if (hasTm == 1) {
              for (j in 1:J) {
                LTn[j] <<- LT[j] + (1 - dl) * (lamNew[j] - lamC[i, j])
                if (model$Kj[j] > 0) {
                  dT <- dT - model$Kj[j] * (LTn[j] - LT[j])
                  if (model$nTm[j] > 0 & LTn[j] != LT[j]) dT <- dT + model$nTm[j] * (log(LTn[j]) - log(LT[j]))
                }
              }
            }
            lr <- lr + dU + dT
            if (log(runif(1)) < lr) {
              model$s[i, 1] <<- p1; model$s[i, 2] <<- p2
              for (j in 1:J) lamC[i, j] <<- lamNew[j]
              if (mki == 0) { for (j in 1:J) LU[j] <<- LUn[j] }
              else if (hasTm == 1) { for (j in 1:J) LT[j] <<- LTn[j] }
              nAcc <<- nAcc + 1
            }
          } else {
            # z = 0: only the point-process prior
            if (log(runif(1)) < lr) {
              model$s[i, 1] <<- p1; model$s[i, 2] <<- p2
              s1 <- p1; s2 <- p2
            }
          }
        }
        # ---- inclusion indicator (augmented rows without marks only)
        if (model$aug[i] == 1 & mki == 0) {
          zi <- model$z[i]
          s1 <- model$s[i, 1]; s2 <- model$s[i, 2]
          if (zi == 0) {
            for (j in 1:J) lamC[i, j] <<- exp(model$lr0[j] + a * ((model$traps[j, 1] - s1)^2 + (model$traps[j, 2] - s2)^2))
          }
          dz <- 1 - 2 * zi
          gs <- floor(s1 / pix) + floor(s2 / pix) * nx + 1
          lr <- dz * (log(psi) - log(1 - psi) + model$lpi[gs] - 2 * log(pix) - model$lunif[1])
          if (hasMark == 1) {
            lm <- 0
            for (l in 1:L) if (model$Kmj[l] > 0) {
              pm <- exp(lpm + a * ((model$mtraps[l, 1] - s1)^2 + (model$mtraps[l, 2] - s2)^2))
              lm <- lm + model$Kmj[l] * log(1 - pm)
            }
            lr <- lr + dz * lm
          }
          dU <- 0
          for (j in 1:J) {
            LUn[j] <<- LU[j] + dz * lamC[i, j]
            if (LUn[j] < 0) LUn[j] <<- 0
            if (model$Kj[j] > 0) {
              dU <- dU - model$Kj[j] * (LUn[j] - LU[j])
              if (model$nU[j] > 0 & LUn[j] != LU[j]) dU <- dU + model$nU[j] * (log(LUn[j]) - log(LU[j]))
            }
          }
          lr <- lr + dU
          if (log(runif(1)) < lr) {
            model$z[i] <<- 1 - zi
            for (j in 1:J) LU[j] <<- LUn[j]
          }
        }
      }
      model$calculate(calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      nRuns <<- nRuns + 1
      if (nRuns %% 50 == 0) {
        if (nTry > 0) {
          ar <- nAcc / nTry
          scale <<- scale * exp(0.5 * (ar - 0.3) / sqrt(nRuns / 50))
        }
        nAcc <<- 0; nTry <<- 0
      }
    },
    methods = list(reset = function() {
      scale <<- scale0; nAcc <<- 0; nTry <<- 0; nRuns <<- 0
    })
  )
  for (nm in c("dSMRHab", "rSMRHab", "dSMRSight", "rSMRSight", "dSMRPool", "rSMRPool",
               "dSMRMark", "rSMRMark", "dSMRTel", "rSMRTel", "sampler_SMRlatent"))
    assign(nm, get(nm), envir = .GlobalEnv)
  invisible()
}

smr_code <- function() nimble::nimbleCode({
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, sigma_max)
  lam0 ~ dunif(0, lam0_max)
  if (HASTM) {
    delta ~ dunif(0, 1)
    dID <- delta
  } else {
    dID <- one
  }
  # marking-phase baseline (its sampler is removed for SMR, where Kmj = 0
  # makes every marking node contribute exactly 0)
  pm0 ~ dunif(0, 1)
  lpm <- log(pm0)
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
  if (HASP) {
    for (k in 1:PP) bP[k] ~ dnorm(0, sd = 2)
    if (POIS) {
      if (PP == 1) {
        lr0[1:J] <- log(lam0) + XP[1:J, 1] * bP[1]
      } else {
        lr0[1:J] <- log(lam0) + (XP[1:J, 1:PP] %*% bP[1:PP])[1:J, 1]
      }
    } else {
      if (PP == 1) {
        lr0[1:J] <- log(expit(logit(lam0) + XP[1:J, 1] * bP[1]))
      } else {
        lr0[1:J] <- log(expit(logit(lam0) + (XP[1:J, 1:PP] %*% bP[1:PP])[1:J, 1]))
      }
    }
  } else {
    lr0[1:J] <- log(lam0) + zeroJ[1:J]
  }
  for (i in 1:M) {
    pz[i] <- aug[i] * psi + 1 - aug[i]
    z[i] ~ dbern(pz[i])
    s[i, 1:2] ~ dSMRHab(lpi = lpi[1:GR], hab = hab[1:GR], pix = pix, nx = nx, ny = ny, z = z[i], lunif = lunif)
    lam[i, 1:J] <- exp(lr0[1:J] - ((traps[1:J, 1] - s[i, 1])^2 + (traps[1:J, 2] - s[i, 2])^2) / (2 * sigma^2))
    xm[i, 1:L] ~ dSMRMark(s = s[i, 1:2], lp0 = lpm, sigma = sigma, traps = mtraps[1:L, 1:2],
                          Kj = Kmj[1:L], L = L, z = z[i])
  }
  for (i in 1:MK) {
    y[i, 1:J] ~ dSMRSight(lam = lam[i, 1:J], Kj = Kj[1:J], delta = dID, mk = mk[i], z = z[i], pois = POISD, J = J)
    tel[i, 1:4] ~ dSMRTel(s = s[i, 1:2], sigma = sigma)
  }
  nU[1:J] ~ dSMRPool(lam = lam[1:M, 1:J], z = z[1:M], mk = mk[1:M], Kj = Kj[1:J], delta = one, type = 0, R = M, J = J)
  nTm[1:J] ~ dSMRPool(lam = lam[1:MK, 1:J], z = z[1:MK], mk = mk[1:MK], Kj = Kj[1:J], delta = dID, type = 1, R = MK, J = J)
  N <- sum(z[1:M])
})

# Bounding rectangle of the mask (pixel index g = ix + iy * nx), hab flags.
smr_rect <- function(prep) {
  ix <- floor(prep$mxy[, 1] / prep$pix + 1e-6); iy <- floor(prep$mxy[, 2] / prep$pix + 1e-6)
  nx <- max(ix) + 1L; ny <- max(iy) + 1L
  list(nx = as.integer(nx), ny = as.integer(ny), GR = as.integer(nx * ny), g = ix + iy * nx + 1)
}

# Rough moment estimate of the number of unmarked animals (inits and M only)
smr_nu_guess <- function(prep) {
  per <- if (prep$m > 0) sum(prep$yID) / prep$m else NA
  dl <- if (prep$has_tm && (sum(prep$yID) + sum(prep$nTm)) > 0) sum(prep$yID) / (sum(prep$yID) + sum(prep$nTm)) else 1
  per <- per / max(dl, 0.05)
  if (!is.finite(per) || per <= 0) return(max(sum(prep$nU), 20))
  sum(prep$nU) / per
}

smr_tops <- function(prep) c("psi", "sigma", "lam0", if (prep$has_tm) "delta", if (prep$has_marking) "pm0",
                             if (ncol(prep$XD)) paste0("bD[", seq_len(ncol(prep$XD)), "]"),
                             if (ncol(prep$XP)) paste0("bP[", seq_len(ncol(prep$XP)), "]"))

smr_engine_nimble <- function(prep, settings) {
  smr_define_nimble()
  PD <- ncol(prep$XD); PP <- ncol(prep$XP)
  rc <- smr_rect(prep)
  m <- prep$m; J <- prep$J
  MK <- smr_bucket(max(m, 1), 25)
  nug <- smr_nu_guess(prep)
  naug <- if (!is.null(settings$M_aug)) settings$M_aug else smr_bucket(max(3 * nug, 120), 40)
  M <- MK + naug
  gsmr <- prep$has_marking
  hab <- numeric(rc$GR); hab[rc$g] <- 1
  XR <- matrix(0, rc$GR, max(PD, 1)); if (PD) XR[rc$g, ] <- prep$XD
  offR <- rep(-1e10, rc$GR); offR[rc$g] <- 0
  lpi0 <- rep(-1e10, rc$GR); lpi0[rc$g] <- -log(prep$G)
  mk <- c(rep(1, m), rep(0, M - m))
  aug <- if (gsmr) rep(1, M) else c(rep(0, m), rep(1, M - m))
  y <- matrix(0, MK, J); if (m) y[seq_len(m), ] <- prep$yID
  tel <- matrix(0, MK, 4); if (m) tel[seq_len(m), ] <- prep$telS
  pois <- prep$pois
  lam0_max <- if (pois) (if (is.null(settings$lam0_max)) 10 else settings$lam0_max) else 1
  constants <- list(M = M, MK = MK, J = J, GR = rc$GR, nx = rc$nx, ny = rc$ny, HASD = PD > 0, HASP = PP > 0,
                    PD = max(PD, 1), PP = max(PP, 1), HASTM = prep$has_tm,
                    POIS = pois, L = if (gsmr) prep$mk$L else 2L)
  data <- list(y = y, nU = prep$nU, traps = prep$txy, Kj = as.numeric(prep$Kj), XR = XR, offR = offR, lpi0 = lpi0,
               hab = hab, pix = prep$pix, lunif = -log(prep$G * prep$pix^2), sigma_max = prep$sigma_max, lam0_max = lam0_max, zeroJ = rep(0, J),
               XP = if (PP) prep$XP else matrix(0, J, 1), mk = mk, aug = aug, one = 1, POISD = as.numeric(pois))
  data$nTm <- if (prep$has_tm) prep$nTm else rep(0, J)
  data$tel <- if (prep$has_tel) tel else matrix(0, MK, 4)
  if (gsmr) {
    xm <- matrix(0, M, prep$mk$L); xm[seq_len(m), ] <- prep$mk$xm
    data$xm <- xm; data$mtraps <- prep$mk$mtxy; data$Kmj <- as.numeric(prep$mk$Kmj)
  } else {
    data$xm <- matrix(0, M, 2); data$mtraps <- matrix(0, 2, 2); data$Kmj <- c(0, 0)
  }
  # starting values
  sg0 <- if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else prep$trap_spacing / prep$cscale
  sg0 <- min(sg0, prep$sigma_max / 2)
  inits <- function() {
    psi0 <- min(0.8, max(0.1, 1.3 * nug / (M - (if (gsmr) 0 else m))))
    if (gsmr) psi0 <- min(0.8, (m + 1.3 * nug) / M)
    sg <- sg0 * stats::runif(1, 0.85, 1.2)
    # lambda0 guess: identified sightings per marked animal / (delta * sum_j K_j exp(-d^2/2 sg^2)) at the camera-grid centre
    cen <- colMeans(prep$txy)
    ex <- sum(prep$Kj * exp(-((prep$txy[, 1] - cen[1])^2 + (prep$txy[, 2] - cen[2])^2) / (2 * sg^2)))
    per <- if (m) max(sum(prep$yID) / m, 0.2) else 1
    l0 <- min(per / max(ex, 1e-3), 0.9 * lam0_max) * stats::runif(1, 0.8, 1.2)
    if (!pois) l0 <- min(l0, 0.5)
    ini <- list(psi = psi0, sigma = sg, lam0 = l0)
    if (prep$has_tm) ini$delta <- stats::runif(1, 0.6, 0.9)
    ini$pm0 <- if (gsmr) stats::runif(1, 0.1, 0.3) else 0.5
    if (PD) ini$bD <- stats::rnorm(PD, 0, 0.2)
    if (PP) ini$bP <- stats::rnorm(PP, 0, 0.2)
    s0 <- prep$mxy[sample.int(prep$G, M, replace = TRUE), , drop = FALSE]
    for (i in seq_len(m)) {
      wy <- prep$yID[i, ]; xyl <- prep$txy[rep(seq_len(J), wy), , drop = FALSE]
      if (gsmr) { wx <- prep$mk$xm[i, ]; xyl <- rbind(xyl, prep$mk$mtxy[rep(seq_len(prep$mk$L), wx), , drop = FALSE]) }
      if (prep$has_tel && prep$telS[i, 1] > 0) xyl <- rbind(xyl, matrix(prep$telS[i, 2:3], 1))
      if (nrow(xyl)) {
        c0 <- colMeans(xyl)
        gi <- which.min((prep$mxy[, 1] - c0[1])^2 + (prep$mxy[, 2] - c0[2])^2)
        s0[i, ] <- prep$mxy[gi, ]
      }
    }
    s0 <- s0 + matrix(stats::runif(2 * M, -0.45, 0.45) * prep$pix, M, 2)
    ini$s <- s0
    z0 <- c(rep(1, m), stats::rbinom(M - m, 1, psi0))
    if (gsmr) z0[seq_len(m)] <- 1
    ini$z <- z0
    ini
  }
  tops <- smr_tops(prep)
  mon_s <- isTRUE(settings$.monitor_latent)
  key <- paste("smr", M, MK, J, rc$GR, PD, PP, prep$has_tm, gsmr, prep$has_tel, pois,
               if (gsmr) prep$mk$L, mon_s, settings$sampler_scale, sep = "_")
  sc <- 0.5 * sg0
  mons <- c("psi", "sigma", "lam0", "N", if (prep$has_tm) "delta", if (gsmr) "pm0", if (PD) "bD", if (PP) "bP",
            if (mon_s) c("s", "z"))
  sumP <- c("psi", "sigma", "lam0", "N", if (prep$has_tm) "delta", if (gsmr) "pm0",
            if (PD) paste0("bD[", seq_len(PD), "]"), if (PP) paste0("bP[", seq_len(PP), "]"))
  r <- nimble_run(smr_code(), constants, data, inits = inits, monitors = mons,
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed, summary_params = sumP,
                  configure = function(conf, mdl) {
                    conf$removeSamplers(c("s", "z"))
                    if (!gsmr) conf$removeSamplers("pm0")
                    conf$addSampler(target = c("s", "z"), type = "sampler_SMRlatent",
                                    control = list(M = M, MK = MK, J = J, L = constants$L, nx = rc$nx, ny = rc$ny,
                                                   pois = as.numeric(pois), hasTm = as.numeric(prep$has_tm),
                                                   hasMark = as.numeric(gsmr), scale = sc))
                    # extra joint update of the correlated top-level parameters
                    conf$addSampler(target = tops, type = "RW_block", control = list(adaptInterval = 100))
                  })
  list(chains = r$samples, t_total = r$diagnostics$t_total, t_setup = r$diagnostics$t_setup,
       M = M, MK = MK, naug = naug)
}
