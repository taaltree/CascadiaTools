# jolly-seber: Nimble engines.
#
# Both engines fit the same model: the restricted dynamic occupancy (RDO)
# formulation of the Jolly-Seber model with data augmentation (Royle &
# Dorazio 2008, 2012; Kery & Schaub 2012 ch. 10):
#   z[i,1] ~ Bern(gamma[1]);  z[i,t] ~ Bern(phi * z[i,t-1] + gamma[t] * A[i,t-1])
#   y[i,t] ~ Bern(p[t] * z[i,t]),   A[i,t-1] = 1 if i has not yet entered.
# The entry hazards gamma[t] are not given independent uniform priors (that
# puts most prior mass on Nsuper near M); they are derived from the
# superpopulation parameters psi ~ U(0,1) and b ~ Dirichlet(1,...,1):
#   gamma[1] = psi b[1],  gamma[t] = psi b[t] / (1 - psi sum_{k<t} b[k]),
# so that Nsuper ~ Binomial(M, psi) is a priori uniform on 0..M and the entry
# probabilities b are uniform on the simplex (Schwarz & Arnason 1996 POPAN
# parameters; the "superpopulation" formulation of Kery & Schaub 2012 10.3).
#
#  marginal (default): z summed out by the forward algorithm over the three
#    states (not yet entered, alive, dead) for each unique capture history,
#    with the M - n all-zero augmented histories handled as one weighted
#    pattern (Saracco & Yackulic 2023 marginalisation). Latent N_t, entrants
#    and Nsuper are then drawn exactly from their conditional posterior per
#    draw (post.R).
#  latent: the literal RDO model with every z[i,t] sampled (binary samplers);
#    slower, used as an in-module cross-check.

js_define_nimble <- function() {
  eco_require("nimble")
  suppressPackageStartupMessages(library(nimble))
  if (exists("dJSmarg_eco", envir = .GlobalEnv, inherits = FALSE)) return(invisible())
  nf <- nimble::nimbleFunction
  dJSmarg_eco <- nf(
    run = function(x = double(1), ch = double(2), grp = double(1), gamma = double(1), phi = double(2),
                   p = double(1), pig = double(1), nT = integer(0), U = integer(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      ll <- 0
      for (u in 1:U) {
        if (x[u] > 0) {
          lik <- 0
          for (g in 1:2) {
            if (pig[g] > 0 & (grp[u] == 0 | grp[u] == g)) {
              a1 <- 1 - gamma[1]; a2 <- gamma[1]; a3 <- 0
              if (ch[u, 1] == 1) { a1 <- 0; a2 <- a2 * p[1] } else { a2 <- a2 * (1 - p[1]) }
              for (t in 2:nT) {
                n1 <- a1 * (1 - gamma[t])
                n2 <- a1 * gamma[t] + a2 * phi[g, t - 1]
                n3 <- a3 + a2 * (1 - phi[g, t - 1])
                if (ch[u, t] == 1) { a1 <- 0; a2 <- n2 * p[t]; a3 <- 0 }
                else { a1 <- n1; a2 <- n2 * (1 - p[t]); a3 <- n3 }
              }
              lik <- lik + pig[g] * (a1 + a2 + a3)
            }
          }
          if (lik <= 0) return(-Inf)
          ll <- ll + x[u] * log(lik)
        }
      }
      if (log) return(ll) else return(exp(ll))
    })
  rJSmarg_eco <- nf(
    run = function(n = integer(0), ch = double(2), grp = double(1), gamma = double(1), phi = double(2),
                   p = double(1), pig = double(1), nT = integer(0), U = integer(0)) {
      returnType(double(1))
      return(numeric(U))
    })
  assign("dJSmarg_eco", dJSmarg_eco, envir = .GlobalEnv)
  assign("rJSmarg_eco", rJSmarg_eco, envir = .GlobalEnv)
  invisible()
}

js_priors_code <- quote({
  psi ~ dunif(0, 1)
  b[1:T] ~ ddirch(alpha[1:T])
  cb[1] <- 0
  for (t in 2:T) cb[t] <- cb[t - 1] + b[t - 1]
  for (t in 1:T) gamma[t] <- psi * b[t] / (1 - psi * cb[t])
  for (k in 1:Kp_lev) beta_p[k] ~ dlogis(0, 1)
  if (Kp > Kp_lev) {
    for (k in (Kp_lev + 1):Kp) beta_p[k] ~ dnorm(0, sd = slope_sd)
  }
  for (k in 1:Kphi_lev) beta_phi[k] ~ dlogis(0, 1)
  # design matrices carry one extra all-zero column with coefficient fixed at
  # 0, so inprod never runs over a length-1 range (nimble C++ codegen issue)
  beta_p[Kp + 1] <- 0
  beta_phi[Kphi + 1] <- 0
  for (t in 1:T) logit(p[t]) <- inprod(Xp[t, 1:(Kp + 1)], beta_p[1:(Kp + 1)])
  for (t in 1:(T - 1)) lphi[t] <- inprod(Xphi[t, 1:(Kphi + 1)], beta_phi[1:(Kphi + 1)])
  if (G == 1) {
    pig[1] <- 1
    pig[2] <- 0
    for (t in 1:(T - 1)) {
      phi[1, t] <- ilogit(lphi[t])^dt[t]
      phi[2, t] <- phi[1, t]
    }
  } else {
    beta_grp ~ dnorm(0, sd = slope_sd)
    pi_grp ~ dunif(0, 1)
    pig[1] <- 1 - pi_grp
    pig[2] <- pi_grp
    for (t in 1:(T - 1)) {
      phi[1, t] <- ilogit(lphi[t])^dt[t]
      phi[2, t] <- ilogit(lphi[t] + beta_grp)^dt[t]
    }
  }
})

js_code <- function(approach) {
  lik <- if (approach == "marginal") quote({
    w[1:U] ~ dJSmarg_eco(ch = ch[1:U, 1:T], grp = grp[1:U], gamma = gamma[1:T], phi = phi[1:2, 1:(T - 1)],
                         p = p[1:T], pig = pig[1:2], nT = T, U = U)
  }) else quote({
    for (i in 1:M) {
      gi[i] ~ dbern(pig[2])
      z[i, 1] ~ dbern(gamma[1])
      rec[i, 1] <- 1 - z[i, 1]
      for (t in 2:T) {
        z[i, t] ~ dbern(phi[gi[i] + 1, t - 1] * z[i, t - 1] + gamma[t] * rec[i, t - 1])
        rec[i, t] <- rec[i, t - 1] * (1 - z[i, t])
      }
      for (t in 1:T) y[i, t] ~ dbern(p[t] * z[i, t])
    }
    for (t in 1:T) N[t] <- sum(z[1:M, t])
    for (t in 2:T) B[t - 1] <- inprod(z[1:M, t], rec[1:M, t - 1])
    Nsuper <- M - sum(rec[1:M, T])
    for (i in 1:M) ent2[i] <- gi[i] * (1 - rec[i, T])
    Nsuper2 <- sum(ent2[1:M])
  })
  body <- c(as.list(js_priors_code)[-1], as.list(lik)[-1])
  as.call(c(as.name("{"), body))
}

js_inits_fun <- function(prep, latent = FALSE) {
  force(prep)
  function() {
    ini <- list(psi = stats::runif(1, 0.3, 0.8), b = rep(1 / prep$T, prep$T),
                beta_p = c(stats::rnorm(prep$Kp_lev, 0, 0.3), rep(0, ncol(prep$Xp) - prep$Kp_lev), 0),
                beta_phi = c(stats::rnorm(ncol(prep$Xphi), 1, 0.3), 0))
    if (prep$G == 2) { ini$beta_grp <- 0; ini$pi_grp <- 0.5 }
    if (latent) {
      M <- prep$M; T <- prep$T
      z <- matrix(0L, M, T)
      f <- apply(prep$ch, 1, function(r) min(which(r > 0)))
      l <- apply(prep$ch, 1, function(r) max(which(r > 0)))
      for (i in seq_len(prep$n)) z[i, f[i]:l[i]] <- 1L
      ini$z <- z
      # full vector: nimble_run copies every init into the model, so the
      # known (data) entries must carry their observed values
      gi <- rep(0, M)
      if (prep$G == 2) {
        gi[seq_len(prep$n)] <- pmax(prep$gcode - 1, 0)
        unk <- c(prep$gcode == 0, rep(TRUE, M - prep$n))
        gi[unk] <- stats::rbinom(sum(unk), 1, 0.5)
      }
      ini$gi <- gi
      ini$psi <- min(0.95, max(0.2, 2 * prep$n / M))
    }
    ini
  }
}

js_engine_nimble <- function(prep, settings) {
  js_define_nimble()
  approach <- if (identical(settings$approach, "latent")) "latent" else "marginal"
  T <- prep$T
  const <- list(T = T, Kp = ncol(prep$Xp), Kp_lev = prep$Kp_lev, Kphi = ncol(prep$Xphi),
                Kphi_lev = prep$Kphi_lev, G = prep$G, slope_sd = settings$prior_sd_slope)
  data <- list(alpha = rep(1, T), Xp = cbind(prep$Xp, 0), Xphi = cbind(prep$Xphi, 0), dt = prep$dt)
  mon <- c("psi", "b", "beta_p", "beta_phi", if (prep$G == 2) c("beta_grp", "pi_grp"))
  if (approach == "marginal") {
    const$U <- prep$Upad
    data$w <- prep$wp; data$ch <- prep$chp; data$grp <- prep$gp
    key <- paste("js_marg", T, prep$Upad, const$Kp, const$Kp_lev, const$Kphi, const$Kphi_lev, prep$G, const$slope_sd, sep = "_")
    niter <- settings$niter; nburn <- settings$nburnin; thin <- settings$thin
  } else {
    M <- prep$M
    const$M <- M
    y <- matrix(0L, M, T); y[seq_len(prep$n), ] <- prep$ch
    gi <- rep(NA_real_, M)
    if (prep$G == 2) gi[seq_len(prep$n)][prep$gcode > 0] <- prep$gcode[prep$gcode > 0] - 1
    else gi[] <- 0
    data$y <- y; data$gi <- gi
    mon <- c(mon, "N", "B", "Nsuper", "Nsuper2")
    key <- paste("js_lat", T, M, const$Kp, const$Kp_lev, const$Kphi, const$Kphi_lev, prep$G, const$slope_sd, sep = "_")
    niter <- settings$niter_latent; nburn <- settings$nburnin_latent; thin <- settings$thin_latent
  }
  conf_fun <- function(conf, m) {
    # an extra joint block update for the correlated scalar parameters
    tops <- c("psi", sprintf("beta_p[%d]", seq_len(const$Kp)), sprintf("beta_phi[%d]", seq_len(const$Kphi)),
              if (prep$G == 2) c("beta_grp", "pi_grp"))
    conf$addSampler(target = tops, type = "RW_block", silent = TRUE)
  }
  r <- nimble_run(js_code(approach), const, data, js_inits_fun(prep, approach == "latent"), mon,
                  niter = niter, nburnin = nburn, nchains = settings$nchains, thin = thin,
                  cache_key = key, configure = conf_fun, seed = settings$seed, calculate = FALSE)
  list(chains = r$samples, t_total = r$diagnostics$t_total, approach = approach)
}
