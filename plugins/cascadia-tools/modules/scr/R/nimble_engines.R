# scr: Nimble engines, ported from ../../R/fit_nimble.R (SCR GPU benchmark)
# and generalised.
#
#  marginal: activity centres, sex (when unknown) and z integrated over the
#            habitat mask (Borchers & Efford 2008 likelihood under data
#            augmentation, Royle et al. 2014 ch. 5/9/11); one compiled
#            distribution dSCRgen covers
#              obs 1 binomial  y_ij ~ Bin(K_j, p_ij) (occasions aggregated)
#              obs 2 Bernoulli y_ijk ~ Bern(1 - (1 - p_ijk)^u_jk)  (per occasion)
#              obs 3 Poisson   y_ijk ~ Pois(u_jk lambda_ijk)        (per occasion or aggregated)
#            with half-normal p (or hazard lambda) = det0_jk(s) exp(-d^2 / 2 sigma_j(s)^2),
#            link(det0_jk(s)) = link(det0) + X_jk beta_P + male * beta_sexP
#            + trap-specific behavioural response, log sigma_j(s) = log sigma
#            + X_j beta_S + male * beta_sexS, sex ratio phi (= proportion female)
#            with unknown sex summed out (sex as a latent individual covariate, Royle et al. 2014; secr hcov).
#            RW_block on all top-level parameters.
#  latent  : data augmentation with z[i], s[i, ] sampled (Royle & Young 2008;
#            Royle et al. 2014); custom dBinomSCR skips z = 0 individuals,
#            dHabPP is the pixelwise inhomogeneous point process. Binomial
#            model with detector-level p0 covariates only (cross-check).
#
# Arrays are padded to bucketed sizes so that simulation replicates of the
# same design reuse one compiled model (see nimble_run()).

scr_bucket <- function(x, b = 50) as.integer(max(b, ceiling(x / b) * b))

scr_define_nimble <- function() {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))  # model building needs nimble attached
  if (exists("dSCRgen", envir = .GlobalEnv, inherits = FALSE)) return(invisible())
  nf <- nimble::nimbleFunction
  dSCRgen <- nf(
    run = function(x = double(0), psi = double(0), ldet0 = double(0), sigma = double(0), beta = double(1),
                   phi = double(0), bidx = double(1), XD = double(2), XP = double(2), XS = double(2),
                   d2 = double(2), use = double(2), pj = double(1), yp = double(2), bp = double(2),
                   pstart = double(1), pend = double(1), sexi = double(1), nobs = double(0), M = double(0),
                   G = integer(0), J = integer(0), Kq = integer(0), obs = integer(0), S = integer(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      PD <- bidx[2]; oP <- bidx[3]; PP <- bidx[4]; oS <- bidx[5]; PS <- bidx[6]
      iSP <- bidx[7]; iSS <- bidx[8]; iB <- bidx[9]
      # density: pi_g = softmax(XD beta_D)
      lpi <- numeric(G)
      if (PD > 0) {
        for (g in 1:G) { tv <- 0; for (q in 1:PD) tv <- tv + XD[g, q] * beta[q]; lpi[g] <- tv }
        mx <- max(lpi); lz <- mx + log(sum(exp(lpi - mx)))
        for (g in 1:G) lpi[g] <- lpi[g] - lz
      } else {
        for (g in 1:G) lpi[g] <- -log(G)
      }
      JK <- J * Kq
      bB <- 0
      if (iB > 0) bB <- beta[iB]
      emB <- exp(bB) - 1
      nob <- nobs
      NI <- length(pstart)
      LI <- matrix(0, NI, 2)
      # Work per sex s: kernel Km[g, j] = exp(-d2 / (2 sigma_js^2)) is computed once
      # (G x J exps); p (or lambda) at (g, j, k) = P0[r] Km[g, j] with r = (k-1) J + j.
      # Binomial / Bernoulli: Lm[g, r] = log(1 - p) is stored and reused by the
      # detected individuals; Poisson needs no transcendental per (g, j, k).
      a <- numeric(J); eta <- numeric(JK); lq <- numeric(JK); P0 <- numeric(JK)
      Km <- matrix(0, G, J)
      nL <- JK
      if (obs == 3) nL <- 1
      Lm <- matrix(0, G, nL)
      base <- numeric(G)
      v <- numeric(G)
      f0 <- 0
      for (s in 1:S) {
        for (j in 1:J) {
          ls <- log(sigma)
          if (PS > 0) for (q in 1:PS) ls <- ls + XS[j, q] * beta[oS + q]
          if (s == 2 & iSS > 0) ls <- ls + beta[iSS]
          a[j] <- -0.5 * exp(-2 * ls)
        }
        for (r in 1:JK) {
          tv <- ldet0
          if (PP > 0) for (q in 1:PP) tv <- tv + XP[r, q] * beta[oP + q]
          if (s == 2 & iSP > 0) tv <- tv + beta[iSP]
          eta[r] <- tv
          if (obs == 3) lq[r] <- tv
          else {
            if (tv > 0) lq[r] <- -log1p(exp(-tv)) else lq[r] <- tv - log1p(exp(tv))
          }
          P0[r] <- exp(lq[r])
        }
        for (j in 1:J) for (g in 1:G) Km[g, j] <- exp(a[j] * d2[g, j])
        # base[g] = log pi_g + log P(no detection anywhere | centre at g, sex s)
        # (pixel g is always the innermost loop: column-major arrays are read contiguously)
        for (g in 1:G) base[g] <- lpi[g]
        if (obs == 3) {
          for (j in 1:J) {
            tv <- 0; for (k in 1:Kq) tv <- tv + use[j, k] * P0[(k - 1) * J + j]
            for (g in 1:G) base[g] <- base[g] - tv * Km[g, j]
          }
        } else {
          for (k in 1:Kq) {
            for (j in 1:J) {
              u <- use[j, k]
              if (u > 0) {
                r <- (k - 1) * J + j
                pr <- P0[r]
                for (g in 1:G) {
                  xx <- pr * Km[g, j]
                  # log(1 - x): series exact to double precision for x < 1e-4 (most far pixels)
                  if (xx < 1e-4) L <- -xx * (1 + xx * (0.5 + xx / 3)) else L <- log1p(-xx)
                  Lm[g, r] <- L
                  base[g] <- base[g] + u * L
                }
              }
            }
          }
        }
        w <- 1
        if (S == 2) { if (s == 1) w <- phi else w <- 1 - phi }
        mx <- max(base)
        f0 <- f0 + w * exp(mx) * sum(exp(base - mx))
        if (nob > 0) {
          for (i in 1:nob) {
            doit <- 1
            if (S == 2) { if (sexi[i] >= 0 & sexi[i] != (s - 1)) doit <- 0 }
            if (doit == 0) LI[i, s] <- -Inf
            else {
              for (g in 1:G) v[g] <- base[g]
              for (p in pstart[i]:pend[i]) {
                j <- pj[p]
                for (k in 1:Kq) {
                  yy <- yp[p, k]; bk <- bp[p, k]
                  if (yy > 0 | bk > 0) {
                    u <- use[j, k]; r <- (k - 1) * J + j
                    if (obs == 1) {
                      for (g in 1:G) v[g] <- v[g] + yy * (lq[r] + a[j] * d2[g, j] - Lm[g, r])
                    } else if (obs == 2) {
                      if (bk > 0) {
                        e <- eta[r] + bB
                        pb <- 1 / (1 + exp(-e))
                        for (g in 1:G) {
                          xx <- pb * Km[g, j]
                          if (xx < 1e-4) L1 <- -xx * (1 + xx * (0.5 + xx / 3)) else L1 <- log1p(-xx)
                          if (yy > 0) v[g] <- v[g] + log1p(-exp(u * L1)) - u * Lm[g, r]
                          else v[g] <- v[g] + u * (L1 - Lm[g, r])
                        }
                      } else {
                        for (g in 1:G) v[g] <- v[g] + log1p(-exp(u * Lm[g, r])) - u * Lm[g, r]
                      }
                    } else {
                      lqb <- lq[r]; cm <- 0
                      if (bk > 0) { lqb <- lq[r] + bB; cm <- u * P0[r] * emB }
                      cy <- yy * (log(u) + lqb)
                      for (g in 1:G) v[g] <- v[g] + cy + yy * a[j] * d2[g, j] - cm * Km[g, j]
                    }
                  }
                }
              }
              mx <- max(v)
              LI[i, s] <- mx + log(sum(exp(v - mx)))
            }
          }
        }
      }
      tot <- 0
      if (nob > 0) {
        if (S == 1) { for (i in 1:nob) tot <- tot + LI[i, 1] }
        else {
          lphi <- log(phi); lphm <- log(1 - phi)
          for (i in 1:nob) {
            c1 <- lphi + LI[i, 1]; c2 <- lphm + LI[i, 2]
            mx <- max(c1, c2)
            tot <- tot + mx + log(exp(c1 - mx) + exp(c2 - mx))
          }
        }
      }
      tot <- tot + nobs * log(psi) + (M - nobs) * log(psi * f0 + 1 - psi)
      if (log) return(tot) else return(exp(tot))
    })
  rSCRgen <- nf(
    run = function(n = integer(0), psi = double(0), ldet0 = double(0), sigma = double(0), beta = double(1),
                   phi = double(0), bidx = double(1), XD = double(2), XP = double(2), XS = double(2),
                   d2 = double(2), use = double(2), pj = double(1), yp = double(2), bp = double(2),
                   pstart = double(1), pend = double(1), sexi = double(1), nobs = double(0), M = double(0),
                   G = integer(0), J = integer(0), Kq = integer(0), obs = integer(0), S = integer(0)) {
      returnType(double(0)); return(1)
    })
  dBinomSCR <- nf(
    run = function(x = double(1), s = double(1), lp0 = double(1), sigma = double(0),
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
        p <- exp(lp0[j] + a * ((traps[j, 1] - s[1])^2 + (traps[j, 2] - s[2])^2))
        ll <- ll + dbinom(x[j], Kj[j], p, log = TRUE)
      }
      if (log) return(ll) else return(exp(ll))
    })
  rBinomSCR <- nf(
    run = function(n = integer(0), s = double(1), lp0 = double(1), sigma = double(0),
                   traps = double(2), Kj = double(1), J = integer(0), z = double(0)) {
      returnType(double(1))
      out <- numeric(J)
      if (z == 0) return(out)
      for (j in 1:J) out[j] <- rbinom(1, Kj[j], exp(lp0[j] - ((traps[j, 1] - s[1])^2 + (traps[j, 2] - s[2])^2) / (2 * sigma^2)))
      return(out)
    })
  # Pixelwise point process for s on the bounding rectangle of the mask;
  # hab[g] = 0 marks non-habitat pixels.
  dHabPP <- nf(
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
  rHabPP <- nf(
    run = function(n = integer(0), lpi = double(1), hab = double(1), pix = double(0), nx = integer(0),
                   ny = integer(0)) {
      returnType(double(1))
      g <- rcat(1, exp(lpi) * hab) - 1
      out <- numeric(2)
      out[1] <- (g %% nx + runif(1)) * pix
      out[2] <- (floor(g / nx) + runif(1)) * pix
      return(out)
    })
  for (nm in c("dSCRgen", "rSCRgen", "dBinomSCR", "rBinomSCR", "dHabPP", "rHabPP"))
    assign(nm, get(nm), envir = .GlobalEnv)
  invisible()
}

# General marginal model. POIS, HASB, SEX are build-time constants.
scr_code_marginal <- function() nimble::nimbleCode({
  psi ~ dunif(0, 1)
  if (POIS) {
    lam0 ~ dunif(0, lam_max)
    ldet0 <- log(lam0)
  } else {
    p0 ~ dunif(0, 1)
    ldet0 <- logit(p0)
  }
  sigma ~ dunif(0, sigma_max)
  if (HASB) {
    for (k in 1:NB) beta[k] ~ dnorm(0, sd = 2)
  }
  if (SEX) {
    phi ~ dunif(0, 1)
  }
  ones ~ dSCRgen(psi = psi, ldet0 = ldet0, sigma = sigma, beta = beta[1:NB], phi = phi, bidx = bidx[1:9],
                 XD = XD[1:G, 1:PDX], XP = XP[1:JK, 1:PPX], XS = XS[1:J, 1:PSX], d2 = d2[1:G, 1:J],
                 use = use[1:J, 1:KqX], pj = pj[1:NPMAX], yp = yp[1:NPMAX, 1:KqX], bp = bp[1:NPMAX, 1:KqX],
                 pstart = pstart[1:NMAX], pend = pend[1:NMAX], sexi = sexi[1:NMAX], nobs = nobs, M = M,
                 G = G, J = J, Kq = Kq, obs = OBS, S = S)
})

scr_code_latent <- function() nimble::nimbleCode({
  psi ~ dunif(0, 1)
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, sigma_max)
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
    if (PP == 1) {
      lp0[1:J] <- log(expit(logit(p0) + XP[1:J, 1] * bP[1]))
    } else {
      lp0[1:J] <- log(expit(logit(p0) + (XP[1:J, 1:PP] %*% bP[1:PP])[1:J, 1]))
    }
  } else {
    lp0[1:J] <- log(p0) + zeroJ[1:J]
  }
  for (i in 1:M) {
    z[i] ~ dbern(psi)
    s[i, 1:2] ~ dHabPP(lpi = lpi[1:GR], hab = hab[1:GR], pix = pix, nx = nx, ny = ny)
    y[i, 1:J] ~ dBinomSCR(s = s[i, 1:2], lp0 = lp0[1:J], sigma = sigma, traps = traps[1:J, 1:2],
                          Kj = Kj[1:J], J = J, z = z[i])
  }
  N <- sum(z[1:M])
})

# Internal draw columns: psi, det0 (p0 or lambda0, natural scale), sigma,
# beta1..betaNB (layout prep$bidx / prep$bnames), phi (proportion female).
scr_tops <- function(prep) c("psi", if (prep$obs == "poisson") "lam0" else "p0", "sigma",
                             if (prep$nb) paste0("beta[", seq_len(prep$nb), "]"), if (prep$S == 2) "phi")

scr_inits_top <- function(prep, M) {
  sg <- min(prep$sigma_max / 2, if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else prep$trap_spacing / prep$cscale)
  ini <- list(psi = min(0.9, 1.5 * max(prep$n, 1) / M) * stats::runif(1, 0.8, 1.2),
              sigma = sg * stats::runif(1, 0.8, 1.25))
  if (prep$obs == "poisson") {
    # rough lambda0: mean samples per detected individual-trap pair / effort
    ini$lam0 <- min(prep$lam_max / 2, max(0.05, 0.3 * stats::runif(1, 0.7, 1.4)))
  } else ini$p0 <- stats::runif(1, 0.1, 0.3)
  if (prep$nb) ini$beta <- stats::rnorm(max(prep$nb, 2), 0, 0.2)   # model vector has >= 2 entries (Nimble dimension rule)
  if (prep$S == 2) ini$phi <- stats::runif(1, 0.4, 0.6)
  ini
}

scr_nimble_inputs <- function(prep) {
  NP <- length(prep$pj)
  NPMAX <- scr_bucket(max(NP, 1), 100); NMAX <- scr_bucket(max(prep$n, 1), 50)
  padv <- function(v, L, fill = 1) c(v, rep(fill, L - length(v)))
  # Nimble drops length-1 dimensions, so every matrix gets >= 2 columns
  # (extra zero columns are never read: the loops use PD/PP/PS/Kq)
  pm <- function(m, nr = nrow(m)) { m <- if (ncol(m)) m else matrix(0, nr, 0); if (ncol(m) < 2) cbind(m, matrix(0, nr, 2 - ncol(m))) else m }
  padm <- function(m, L) pm(rbind(m, matrix(0, L - nrow(m), ncol(m))))
  pois <- prep$obs == "poisson"
  JK <- prep$J * prep$Kq
  constants <- list(G = prep$G, J = prep$J, Kq = prep$Kq, KqX = max(prep$Kq, 2), JK = JK, NPMAX = NPMAX, NMAX = NMAX,
                    NB = max(prep$nb, 2), HASB = prep$nb > 0, SEX = prep$S == 2, POIS = pois,
                    OBS = prep$obs_code, S = prep$S,
                    PDX = max(ncol(prep$XD), 2), PPX = max(ncol(prep$XPq), 2), PSX = max(ncol(prep$XS), 2))
  data <- list(ones = 1, bidx = as.numeric(prep$bidx), XD = pm(prep$XD, prep$G), XP = pm(prep$XPq, JK),
               XS = pm(prep$XS, prep$J), d2 = prep$d2, use = pm(prep$use),
               pj = padv(prep$pj, NPMAX), yp = padm(prep$yp, NPMAX), bp = padm(prep$bp, NPMAX),
               pstart = padv(prep$pstart, NMAX), pend = padv(prep$pend, NMAX), sexi = padv(prep$sexi, NMAX, -1),
               nobs = prep$n, M = prep$M, sigma_max = prep$sigma_max)
  if (pois) data$lam_max <- prep$lam_max
  if (prep$nb == 0) data$beta <- c(0, 0)
  if (prep$S == 1) data$phi <- 1
  list(constants = constants, data = data, NPMAX = NPMAX, NMAX = NMAX, pois = pois)
}

scr_engine_nimble_marginal <- function(prep, settings) {
  scr_define_nimble()
  ni <- scr_nimble_inputs(prep)
  constants <- ni$constants; data <- ni$data; pois <- ni$pois; NPMAX <- ni$NPMAX; NMAX <- ni$NMAX
  tops <- scr_tops(prep)
  key <- paste("scr_gen", settings$sampler, prep$G, prep$J, prep$Kq, NPMAX, NMAX, prep$obs_code, prep$S,
               prep$nb, ncol(prep$XD), ncol(prep$XPq), ncol(prep$XS), sep = "_")
  r <- nimble_run(scr_code_marginal(), constants, data,
                  inits = function() scr_inits_top(prep, prep$M),
                  monitors = c("psi", if (pois) "lam0" else "p0", "sigma", if (prep$nb) "beta", if (prep$S == 2) "phi"),
                  niter = settings$niter, nburnin = settings$nburnin, nchains = settings$nchains,
                  thin = settings$thin, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    conf$removeSamplers(tops)
                    if (identical(settings$sampler, "AF_slice")) conf$addSampler(target = tops, type = "AF_slice")
                    else conf$addSampler(target = tops, type = "RW_block", control = list(adaptInterval = 100))
                  })
  list(chains = lapply(r$samples, scr_rename_draws), t_total = r$diagnostics$t_total,
       t_setup = r$diagnostics$t_setup, engine = "nimble_marginal")
}

# Latent engine: bounding rectangle of the mask, hab = 1 on mask pixels.
scr_rect <- function(prep) {
  ix <- floor(prep$mxy[, 1] / prep$pix + 1e-6); iy <- floor(prep$mxy[, 2] / prep$pix + 1e-6)
  nx <- max(ix) + 1L; ny <- max(iy) + 1L
  g <- ix + iy * nx + 1
  list(nx = as.integer(nx), ny = as.integer(ny), GR = as.integer(nx * ny), g = g)
}

scr_engine_nimble_latent <- function(prep, settings) {
  if (!scr_is_basic(prep) || prep$obs == "poisson")
    stop("the latent Nimble engine supports the binomial model with detector-level p0 covariates only; use approach: marginal")
  scr_define_nimble()
  PD <- ncol(prep$XD); PP <- ncol(prep$XP)
  rc <- scr_rect(prep)
  M <- scr_bucket(if (!is.null(settings$M)) settings$M else max(3 * prep$n, prep$n + 80), 50)
  hab <- numeric(rc$GR); hab[rc$g] <- 1
  XR <- matrix(0, rc$GR, max(PD, 1)); if (PD) XR[rc$g, ] <- prep$XD
  offR <- rep(-1e10, rc$GR); offR[rc$g] <- 0
  lpi0 <- rep(-1e10, rc$GR); lpi0[rc$g] <- -log(prep$G)
  yaug <- matrix(0, M, prep$J); yaug[seq_len(prep$n), ] <- prep$y
  constants <- list(M = M, J = prep$J, GR = rc$GR, nx = rc$nx, ny = rc$ny, HASD = PD > 0, HASP = PP > 0,
                    PD = max(PD, 1), PP = max(PP, 1))
  data <- list(y = yaug, traps = prep$txy, Kj = as.numeric(prep$Kj), XR = XR, offR = offR, lpi0 = lpi0,
               hab = hab, pix = prep$pix, sigma_max = prep$sigma_max, zeroJ = rep(0, prep$J),
               XP = if (PP) prep$XP else matrix(0, prep$J, 1))
  # initial s: detected individuals at the habitat pixel nearest their mean
  # detection location; others at random habitat pixels
  inits <- function() {
    i0 <- scr_inits_top(prep, M)
    ini <- list(psi = i0$psi, p0 = i0$p0, sigma = i0$sigma)
    if (PD) ini$bD <- stats::rnorm(PD, 0, 0.2)
    if (PP) ini$bP <- stats::rnorm(PP, 0, 0.2)
    s0 <- prep$mxy[sample.int(prep$G, M, replace = TRUE), , drop = FALSE]
    for (i in seq_len(prep$n)) {
      w <- prep$y[i, ]
      c0 <- colSums(prep$txy * w) / sum(w)
      gi <- which.min((prep$mxy[, 1] - c0[1])^2 + (prep$mxy[, 2] - c0[2])^2)
      s0[i, ] <- prep$mxy[gi, ]
    }
    s0 <- s0 + matrix(stats::runif(2 * M, -0.45, 0.45) * prep$pix, M, 2)
    ini$s <- s0
    ini$z <- c(rep(1, prep$n), stats::rbinom(M - prep$n, 1, ini$psi))
    ini$psi <- max(ini$psi, (prep$n + 1) / M)
    ini
  }
  key <- paste("scr_lat", M, prep$J, rc$GR, PD, PP, sep = "_")
  sc <- 0.5 * (if (is.finite(prep$sigma_guess)) prep$sigma_guess / prep$cscale else prep$trap_spacing / prep$cscale)
  r <- nimble_run(scr_code_latent(), constants, data, inits = inits,
                  monitors = c("psi", "p0", "sigma", "N", if (PD) "bD", if (PP) "bP"),
                  niter = settings$niter_latent, nburnin = settings$nburnin_latent, nchains = settings$nchains,
                  thin = settings$thin_latent, cache_key = key, seed = settings$seed,
                  configure = function(conf, m) {
                    conf$removeSamplers("s")
                    for (i in seq_len(M)) conf$addSampler(target = sprintf("s[%d, 1:2]", i), type = "RW_block",
                                                          control = list(adaptScaleOnly = TRUE, scale = sc))
                  })
  # bD then bP is exactly the beta layout of a basic model
  ren <- function(m) {
    m <- as.matrix(m); cn <- colnames(m)
    b <- c(grep("^bD\\[", cn, value = TRUE), grep("^bP\\[", cn, value = TRUE))
    out <- m[, c("psi", "p0", "sigma", "N"), drop = FALSE]; colnames(out)[2] <- "det0"
    if (length(b)) { bb <- m[, b, drop = FALSE]; colnames(bb) <- paste0("beta", seq_along(b)); out <- cbind(out, bb) }
    out
  }
  list(chains = lapply(r$samples, ren), t_total = r$diagnostics$t_total,
       t_setup = r$diagnostics$t_setup, engine = "nimble_latent", M = M)
}

# Nimble column names -> internal names psi, det0, sigma, beta1.., phi, N
scr_rename_draws <- function(m) {
  m <- as.matrix(m)
  cn <- colnames(m)
  cn[cn %in% c("p0", "lam0")] <- "det0"
  cn <- sub("^beta\\[(\\d+)\\]$", "beta\\1", cn)
  colnames(m) <- cn
  m
}
