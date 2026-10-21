# ipm: Nimble model code.
#
# Female-based, pre-breeding-census integrated population model with two age
# classes at breeding (1-year-olds and older adults) plus immigrants, as in
# Schaub & Kery (2022) section 9.2 (hoopoe) and chapter 11 (woodchat shrike):
#
#   N1[t+1] ~ Pois( rho * phij[t] * ( f1[t] * (N1[t] + N3[t]) + f2[t] * N2[t] ) )  local recruits
#   N2[t+1] ~ Bin ( phia[t], Ntot[t] )                                             surviving adults
#   N3[t+1] ~ Pois( omega[t] * Ntot[t] )                                           immigrants
#   Ntot[t] = N1[t] + N2[t] + N3[t]
#
# (rho = proportion of females among fledglings, 0.5 by default; immigrants
# breed like 1-year-olds, as in the book's hoopoe model.) Data likelihoods:
#   counts       C[t] ~ Pois(Ntot[t])  or  Normal(Ntot[t], sigma_count)
#   productivity J1[t] ~ Pois(B1[t] f1[t]);  J2[t] ~ Pois(B2[t] f2[t])
#   capture-recapture: age-dependent CJS m-arrays (multinomial; custom dMarr)
# Annual rates have logit/log-normal year effects plus optional covariates.
#
# The model structure is switched with constants (evaluated when the model
# is built): IPM (0 = separate analyses: CJS + productivity only), IMM
# (immigration), AGEF (age-specific fecundity), RE (random year effects),
# CNORM (normal count error), K_* (number of covariates per rate).

ipm_define_nimble <- function() {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))
  if (exists("dMarr", envir = .GlobalEnv, inherits = FALSE)) return(invisible())
  # Multinomial likelihood of one CJS m-array (rows = release occasion
  # t = 1..T-1, columns = first recapture at occasion 2..T, last column =
  # never recaptured). phi1[t] = survival over the first interval after
  # release (juvenile survival for animals marked as young), phia = adult
  # survival afterwards, p[t] = recapture probability at occasion t+1.
  # Loops are written so they never run backwards (nimble `a:b` with a > b
  # counts down, unlike JAGS).
  dMarr <- nimble::nimbleFunction(
    run = function(x = double(2), phi1 = double(1), phia = double(1), p = double(1),
                   rel = double(1), log = integer(0, default = 0)) {
      returnType(double(0))
      nT <- dim(x)[2]
      ll <- 0
      for (t in 1:(nT - 1)) {
        if (rel[t] > 0) {
          s <- phi1[t]
          nr <- 1
          tot <- 0
          ll <- ll + lgamma(rel[t] + 1)
          for (j in t:(nT - 1)) {
            if (j > t) {
              s <- s * phia[j]
              nr <- nr * (1 - p[j - 1])
            }
            pr <- s * nr * p[j]
            tot <- tot + pr
            if (x[t, j] > 0) ll <- ll + x[t, j] * log(pr) - lgamma(x[t, j] + 1)
          }
          pn <- 1 - tot
          if (pn < 1e-300) pn <- 1e-300
          if (x[t, nT] > 0) ll <- ll + x[t, nT] * log(pn) - lgamma(x[t, nT] + 1)
        }
      }
      if (log) return(ll) else return(exp(ll))
    })
  rMarr <- nimble::nimbleFunction(
    run = function(n = integer(0), phi1 = double(1), phia = double(1), p = double(1),
                   rel = double(1)) {
      returnType(double(2))
      nT <- length(rel) + 1
      x <- matrix(0, nT - 1, nT)
      for (t in 1:(nT - 1)) {
        pr <- numeric(nT)
        s <- phi1[t]
        nr <- 1
        tot <- 0
        for (j in t:(nT - 1)) {
          if (j > t) {
            s <- s * phia[j]
            nr <- nr * (1 - p[j - 1])
          }
          pr[j] <- s * nr * p[j]
          tot <- tot + pr[j]
        }
        pr[nT] <- max(0, 1 - tot)
        if (rel[t] > 0) x[t, 1:nT] <- rmulti(1, size = rel[t], prob = pr[1:nT])
      }
      return(x)
    })
  assign("dMarr", dMarr, envir = .GlobalEnv)
  assign("rMarr", rMarr, envir = .GlobalEnv)
  suppressMessages(nimble::registerDistributions(list(dMarr = list(
    BUGSdist = "dMarr(phi1, phia, p, rel)",
    types = c("value = double(2)", "phi1 = double(1)", "phia = double(1)", "p = double(1)", "rel = double(1)"),
    discrete = TRUE, mixedSizes = TRUE)), verbose = FALSE))
  invisible()
}

ipm_code <- quote({
  # ---- priors (Schaub & Kery 2022, ch. 9.2) ----
  mean_phij ~ dunif(0, 1)
  mean_phia ~ dunif(0, 1)
  mean_f1 ~ dunif(0, fmax)
  mu_phij <- logit(mean_phij)
  mu_phia <- logit(mean_phia)
  mu_f1 <- log(mean_f1)
  if (AGEF == 1) {
    mean_f2 ~ dunif(0, fmax)
    mu_f2 <- log(mean_f2)
  }
  p ~ dunif(0, 1)
  if (RE == 1) {
    sd_phij ~ dunif(0, 3)
    sd_phia ~ dunif(0, 3)
    sd_f1 ~ dunif(0, 3)
    if (AGEF == 1) {
      sd_f2 ~ dunif(0, 3)
    }
  }
  if (K_phij > 0) { for (k in 1:K_phij) { beta_phij[k] ~ dnorm(0, sd = 2) } }
  if (K_phia > 0) { for (k in 1:K_phia) { beta_phia[k] ~ dnorm(0, sd = 2) } }
  if (K_f > 0) { for (k in 1:K_f) { beta_f[k] ~ dnorm(0, sd = 2) } }

  # ---- annual rates ----
  for (t in 1:(T - 1)) {
    if (RE == 1) {
      eps_phij[t] ~ dnorm(0, sd = sd_phij)
      eps_phia[t] ~ dnorm(0, sd = sd_phia)
    }
    if (RE == 0) {
      eps_phij[t] <- 0
      eps_phia[t] <- 0
    }
    if (K_phij > 1) { xb_phij[t] <- inprod(beta_phij[1:K_phij], X_phij[t, 1:K_phij]) }
    if (K_phij == 1) { xb_phij[t] <- beta_phij[1] * X_phij[t, 1] }
    if (K_phij == 0) { xb_phij[t] <- 0 }
    if (K_phia > 1) { xb_phia[t] <- inprod(beta_phia[1:K_phia], X_phia[t, 1:K_phia]) }
    if (K_phia == 1) { xb_phia[t] <- beta_phia[1] * X_phia[t, 1] }
    if (K_phia == 0) { xb_phia[t] <- 0 }
    logit(phij[t]) <- mu_phij + xb_phij[t] + eps_phij[t]
    logit(phia[t]) <- mu_phia + xb_phia[t] + eps_phia[t]
    pvec[t] <- p
  }
  for (t in 1:T) {
    if (RE == 1) {
      eps_f1[t] ~ dnorm(0, sd = sd_f1)
    }
    if (RE == 0) {
      eps_f1[t] <- 0
    }
    if (K_f > 1) { xb_f[t] <- inprod(beta_f[1:K_f], X_f[t, 1:K_f]) }
    if (K_f == 1) { xb_f[t] <- beta_f[1] * X_f[t, 1] }
    if (K_f == 0) { xb_f[t] <- 0 }
    log(f1[t]) <- mu_f1 + xb_f[t] + eps_f1[t]
    if (AGEF == 1) {
      if (RE == 1) {
        eps_f2[t] ~ dnorm(0, sd = sd_f2)
      }
      if (RE == 0) {
        eps_f2[t] <- 0
      }
      log(f2[t]) <- mu_f2 + xb_f[t] + eps_f2[t]
    }
    if (AGEF == 0) {
      f2[t] <- f1[t]
    }
  }

  # ---- productivity data (Poisson; J = fledglings, B = broods monitored) ----
  for (t in 1:T) {
    J1[t] ~ dpois(f1[t] * B1[t])
    if (AGEF == 1) {
      J2[t] ~ dpois(f2[t] * B2[t])
    }
  }

  # ---- capture-recapture data (age-dependent CJS m-arrays) ----
  marr_j[1:(T - 1), 1:T] ~ dMarr(phi1 = phij[1:(T - 1)], phia = phia[1:(T - 1)], p = pvec[1:(T - 1)], rel = rel_j[1:(T - 1)])
  marr_a[1:(T - 1), 1:T] ~ dMarr(phi1 = phia[1:(T - 1)], phia = phia[1:(T - 1)], p = pvec[1:(T - 1)], rel = rel_a[1:(T - 1)])

  # ---- population model and counts (IPM only) ----
  if (IPM == 1) {
    if (IMM == 1) {
      mean_omega ~ dunif(0, 3)
      mu_omega <- log(mean_omega)
      if (RE == 1) {
        sd_omega ~ dunif(0, 3)
      }
      if (K_om > 0) { for (k in 1:K_om) { beta_om[k] ~ dnorm(0, sd = 2) } }
      for (t in 1:(T - 1)) {
        if (RE == 1) {
          eps_om[t] ~ dnorm(0, sd = sd_omega)
        }
        if (RE == 0) {
          eps_om[t] <- 0
        }
        if (K_om > 1) { xb_om[t] <- inprod(beta_om[1:K_om], X_om[t, 1:K_om]) }
    if (K_om == 1) { xb_om[t] <- beta_om[1] * X_om[t, 1] }
        if (K_om == 0) { xb_om[t] <- 0 }
        log(omega[t]) <- mu_omega + xb_om[t] + eps_om[t]
      }
    }
    N1[1] ~ dcat(pinit[1:U])
    N2[1] ~ dcat(pinit[1:U])
    if (IMM == 1) {
      N3[1] ~ dcat(pinit[1:U])
    }
    if (IMM == 0) {
      N3[1] <- 0
    }
    for (t in 1:(T - 1)) {
      N1[t + 1] ~ dpois(rho * phij[t] * (f1[t] * (N1[t] + N3[t]) + f2[t] * N2[t]))
      N2[t + 1] ~ dbin(phia[t], Ntot[t])
      if (IMM == 1) {
        N3[t + 1] ~ dpois(omega[t] * Ntot[t])
      }
      if (IMM == 0) {
        N3[t + 1] <- 0
      }
    }
    for (t in 1:T) {
      Ntot[t] <- N1[t] + N2[t] + N3[t]
    }
    if (CNORM == 1) {
      sigma_count ~ dunif(0.5, sig_max)
      for (t in 1:T) {
        C[t] ~ dnorm(Ntot[t], sd = sigma_count)
      }
    }
    if (CNORM == 0) {
      for (t in 1:T) {
        C[t] ~ dpois(Ntot[t])
      }
    }
  }
})

# Build the nimbleCode object with the structural switches substituted in
# (nimble evaluates `if` on constants at build time; substituting literal
# values makes that explicit and keeps one code path).
ipm_nimble_code <- function(flags) do.call(substitute, list(ipm_code, flags))
