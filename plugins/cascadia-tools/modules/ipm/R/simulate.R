# ipm: data simulation with known truth.
#
# The population (counts and productivity) and the marked sample (capture-
# recapture) are simulated from the same annual vital rates but as separate
# draws, which is the independence assumption of the IPM joint likelihood
# (Schaub & Kery 2022, section 5.5; Besbeas et al. 2002).

ipm_sim_rates <- function(params, design, X) {
  T <- design$n_years
  re <- isTRUE(design$temporal_variation %||% TRUE)
  rn <- function(n, sd) if (re && sd > 0) stats::rnorm(n, 0, sd) else rep(0, n)
  b <- params$juvenile_survival_beta %||% 0
  ba <- params$adult_survival_beta %||% 0
  x <- if (is.null(X)) rep(0, T) else X
  list(
    phij = stats::plogis(stats::qlogis(params$juvenile_survival) + b * x[1:(T - 1)] + rn(T - 1, params$juvenile_survival_sd)),
    phia = stats::plogis(stats::qlogis(params$adult_survival) + ba * x[1:(T - 1)] + rn(T - 1, params$adult_survival_sd)),
    f1 = exp(log(params$fecundity_yearling) + rn(T, params$fecundity_yearling_sd)),
    f2 = exp(log(params$fecundity_adult) + rn(T, params$fecundity_adult_sd)),
    omega = if (isTRUE(design$immigration)) exp(log(params$immigration_rate) + rn(T - 1, params$immigration_sd)) else rep(0, T - 1))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Individual capture histories of the marked sample -> m-arrays.
ipm_sim_marked <- function(rates, design, p) {
  T <- design$n_years
  nj <- rep_len(design$marked_juveniles, T); na <- rep_len(design$marked_adults, T)
  rows <- list()
  for (t in seq_len(T - 1)) {
    for (age in 1:2) {
      n <- if (age == 1) nj[t] else na[t]
      if (n <= 0) next
      ch <- matrix(0L, n, T); ch[, t] <- 1L
      alive <- rep(TRUE, n)
      for (k in t:(T - 1)) {
        s <- if (k == t && age == 1) rates$phij[k] else rates$phia[k]
        alive <- alive & (stats::runif(n) < s)
        ch[alive & stats::runif(n) < p, k + 1] <- 1L
      }
      rows[[length(rows) + 1]] <- list(ch = ch, age = rep(age, n))
    }
  }
  ch <- do.call(rbind, lapply(rows, `[[`, "ch"))
  age <- unlist(lapply(rows, `[[`, "age"))
  list(ch = ch, age = age)
}

ipm_simulate <- function(params, design, seed) {
  set.seed(seed)
  T <- design$n_years
  cov_name <- design$covariate_name %||% "winter_severity"
  use_cov <- isTRUE(design$covariate)
  X <- if (use_cov) as.numeric(scale(stats::rnorm(T))) else NULL
  if (!use_cov) params$juvenile_survival_beta <- 0
  agef <- isTRUE(design$age_specific_fecundity %||% TRUE)
  rates <- ipm_sim_rates(params, design, X)
  if (!agef) rates$f2 <- rates$f1
  rho <- design$sex_ratio %||% 0.5
  # population process
  N <- matrix(0, 3, T)
  ni <- design$N_initial * c(0.4, 0.45, 0.15)
  if (!isTRUE(design$immigration)) ni <- design$N_initial * c(0.45, 0.55, 0)
  N[, 1] <- stats::rpois(3, ni)
  for (t in 1:(T - 1)) {
    Nt <- sum(N[, t])
    N[1, t + 1] <- stats::rpois(1, rho * rates$phij[t] * (rates$f1[t] * (N[1, t] + N[3, t]) + rates$f2[t] * N[2, t]))
    N[2, t + 1] <- stats::rbinom(1, Nt, rates$phia[t])
    N[3, t + 1] <- stats::rpois(1, rates$omega[t] * Nt)
  }
  Ntot <- colSums(N)
  counts <- if (identical(design$count_model, "normal"))
    pmax(0, round(stats::rnorm(T, Ntot, design$count_sd))) else stats::rpois(T, Ntot)
  # productivity: a share of broods monitored, by mother age class
  pm <- design$prop_broods_monitored
  B1 <- stats::rbinom(T, N[1, ] + N[3, ], pm); B2 <- stats::rbinom(T, N[2, ], pm)
  J1 <- stats::rpois(T, B1 * rates$f1); J2 <- stats::rpois(T, B2 * rates$f2)
  prod <- if (agef) data.frame(year = seq_len(T), J1 = J1, B1 = B1, J2 = J2, B2 = B2)
          else data.frame(year = seq_len(T), J1 = J1 + J2, B1 = B1 + B2, J2 = 0, B2 = 0)
  # marked sample
  mk <- ipm_sim_marked(rates, design, params$recapture_prob)
  ma <- ipm_marray(mk$ch, mk$age)
  years <- seq_len(T) + (design$first_year %||% 2001) - 1
  prod$year <- years
  cov <- if (use_cov) stats::setNames(data.frame(years, X), c("year", cov_name)) else NULL
  data <- list(species = "Simulated population", years = years, counts = counts,
               marr_j = ma$marr_j, marr_a = ma$marr_a, prod = prod, age_fecundity = agef,
               covariates = cov,
               effects = if (use_cov) list(juvenile_survival = cov_name) else list(),
               n_marked = c(juveniles = sum(mk$age == 1), adults = sum(mk$age == 2)),
               count_unit = "breeding females")
  # truth
  tr <- c(juvenile_survival = params$juvenile_survival, adult_survival = params$adult_survival)
  tr <- c(tr, if (agef) c(fecundity_yearling = params$fecundity_yearling, fecundity_adult = params$fecundity_adult)
              else c(fecundity = params$fecundity_yearling))
  tr <- c(tr, recapture_prob = params$recapture_prob)
  re <- isTRUE(design$temporal_variation %||% TRUE)
  if (re) {
    tr <- c(tr, juvenile_survival_sd = params$juvenile_survival_sd, adult_survival_sd = params$adult_survival_sd)
    tr <- c(tr, if (agef) c(fecundity_yearling_sd = params$fecundity_yearling_sd, fecundity_adult_sd = params$fecundity_adult_sd)
                else c(fecundity_sd = params$fecundity_yearling_sd))
  }
  if (use_cov) tr[paste0("juvenile_survival_beta_", cov_name)] <- params$juvenile_survival_beta
  if (isTRUE(design$immigration)) {
    tr <- c(tr, immigration_rate = mean(rates$omega), immigration_rate_median = params$immigration_rate)
    if (re) tr <- c(tr, immigration_sd = params$immigration_sd)
  }
  tr <- c(tr, N_first = Ntot[1], N_final = Ntot[T], N_mean = mean(Ntot),
          lambda_geomean = (Ntot[T] / Ntot[1])^(1 / (T - 1)))
  if (isTRUE(design$immigration)) tr <- c(tr, immigrants_mean = mean(N[3, -1]))
  if (re) {
    flags <- list(IMM = as.integer(isTRUE(design$immigration)), AGEF = as.integer(agef))
    lt <- ipm_ltre_draws(
      list(phij = matrix(rates$phij, 1), phia = matrix(rates$phia, 1), f1 = matrix(rates$f1, 1),
           f2 = matrix(if (agef) rates$f2 else rates$f1, 1), omega = matrix(rates$omega, 1),
           N1 = matrix(N[1, ], 1), N2 = matrix(N[2, ], 1), N3 = matrix(N[3, ], 1)),
      means = list(phij = params$juvenile_survival, phia = params$adult_survival,
                   f1 = params$fecundity_yearling, f2 = if (agef) params$fecundity_adult else params$fecundity_yearling,
                   omega = if (isTRUE(design$immigration)) params$immigration_rate else 0),
      flags = flags, rho = rho)
    tr <- c(tr, unlist(lt$contrib[1, ]))
  }
  list(data = data, truth = tr,
       latent = list(N = N, rates = rates, ch = mk$ch, age = mk$age, X = X))
}
