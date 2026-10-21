# jsdm: simulate community data from the HMSC generative model
# (Ovaskainen et al. 2017; Tikhonov et al. 2020).
#
#   L[i, j]  = sum_k X[i, k] Beta[k, j] + sum_h Eta[i, h] Lambda[h, j]
#   presence: Y[i, j] = 1(L[i, j] + e[i, j] > 0),       e ~ N(0, 1)   (probit)
#   count:    Y[i, j] ~ Poisson(exp(L[i, j] + e[i, j])), e ~ N(0, sigma^2) (lognormal Poisson)
#   Beta[, j] ~ N(Gamma t_j, V)  (rows of Beta also phylogenetically correlated
#              with Cov = V[k,k] (rho C + (1 - rho) I) when phylo_rho > 0)
#   Eta[, h] ~ N(0, 1) independent across sites (or an exponential GP over
#              site coordinates when spatial_range_km > 0)
#
# Species traits (and the phylogeny) are drawn once from
# design$community_seed, so species keep their identity across replicates.
# design$community = "random" (default): each replicate draws new species
#   deviations from the trait model (Beta - Gamma T') and new loadings, so
#   Gamma is the hyperparameter the hierarchical model estimates and its
#   coverage is judged over communities (a finite set of 12 species has a
#   realised trait regression that differs from Gamma by sampling error).
# design$community = "fixed": Beta and loadings are also fixed by
#   community_seed (same true species responses in every replicate). Covariates and the trait are drawn
# so that their sample mean/SD are exactly the centre/scale used by fit(),
# so Beta in the truth is per 1 SD of each covariate, as estimated.

`%||%` <- function(a, b) if (is.null(a)) b else a

std_exact <- function(z) as.numeric((z - mean(z)) / stats::sd(z))

jsdm_species_names <- function(ns) sprintf("sp%02d", seq_len(ns))

# Build the fixed true community from params/design.
jsdm_community <- function(params, design, seed = NULL) {
  ns <- design$n_species
  nf <- design$n_factors %||% 2
  covs <- c("intercept", "temperature", "forest")
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(design$community_seed %||% 2024)
  sp <- jsdm_species_names(ns)
  tz <- std_exact(stats::rnorm(ns))                      # standardised log body mass
  Tr <- cbind(intercept = 1, log_body_mass = tz)
  G <- rbind(intercept   = c(params$comm_mean_intercept,   params$trait_log_body_mass_on_intercept),
             temperature = c(params$comm_mean_temperature, params$trait_log_body_mass_on_temperature),
             forest      = c(params$comm_mean_forest,      params$trait_log_body_mass_on_forest))
  colnames(G) <- colnames(Tr)
  sdv <- params$species_sd %||% 0.5
  rho <- design$phylo_rho %||% 0
  tree <- NULL; Q <- diag(ns)
  if (rho > 0) {
    eco_require("ape")
    tree <- ape::rcoal(ns, tip.label = sp)
    C <- ape::vcv(tree, corr = TRUE)[sp, sp]
    Q <- rho * C + (1 - rho) * diag(ns)
  }
  if (!identical(design$community %||% "random", "fixed") && !is.null(seed)) set.seed(seed + 7919)
  E <- t(vapply(seq_along(covs), function(k) as.numeric(t(chol(Q)) %*% stats::rnorm(ns)) * sdv, numeric(ns)))
  B <- G %*% t(Tr) + E
  dimnames(B) <- list(covs, sp)
  lsd <- c(params$loading_sd1 %||% 1, params$loading_sd2 %||% 0.7, params$loading_sd3 %||% 0.5)
  lsd <- rep(lsd, length.out = nf)
  Lambda <- matrix(stats::rnorm(nf * ns), nf, ns) * lsd   # nf = 0: no residual associations (null)
  colnames(Lambda) <- sp
  sig <- rep(params$residual_sd %||% 0.5, ns)
  names(sig) <- sp
  list(species = sp, covs = covs, Tr = Tr, Gamma = G, Beta = B, Lambda = Lambda,
       sigma = sig, tree = tree, rho = rho, body_mass_g = round(exp(3 + 1.2 * tz), 1))
}

# Derived "truth" quantities shared by simulate() and fit()
jsdm_assoc <- function(Omega, extra_var) {
  ns <- ncol(Omega)
  rc <- suppressWarnings(stats::cov2cor(Omega))
  Sig <- Omega + diag(extra_var, ns)
  lc <- stats::cov2cor(Sig)
  list(rescor = rc[upper.tri(rc)], latcor = lc[upper.tri(lc)])
}

pair_names <- function(sp) {
  ns <- length(sp)
  idx <- which(upper.tri(matrix(0, ns, ns)), arr.ind = TRUE)
  paste0(sp[idx[, 1]], "__", sp[idx[, 2]])
}

simulate_jsdm <- function(params, design, seed) {
  comm <- jsdm_community(params, design, seed)
  set.seed(seed)
  n <- design$n_sites; ns <- design$n_species; nf <- nrow(comm$Lambda)
  resp <- design$response %||% "presence"
  ext <- design$extent_km %||% 20
  xy <- cbind(x = stats::runif(n, 0, ext), y = stats::runif(n, 0, ext))
  zt <- std_exact(stats::rnorm(n))
  zf <- std_exact((design$cov_correlation %||% 0) * zt + stats::rnorm(n))
  X <- cbind(intercept = 1, temperature = zt, forest = zf)
  rng <- design$spatial_range_km %||% 0
  if (rng > 0) {
    D <- as.matrix(stats::dist(xy))
    R <- t(chol(exp(-D / rng) + diag(1e-8, n)))
    Eta <- R %*% matrix(stats::rnorm(n * nf), n, nf)
  } else Eta <- matrix(stats::rnorm(n * nf), n, nf)
  L <- X %*% comm$Beta + Eta %*% comm$Lambda
  if (resp == "presence") {
    Y <- (L + matrix(stats::rnorm(n * ns), n, ns) > 0) * 1
  } else {
    Y <- matrix(stats::rpois(n * ns, exp(L + matrix(stats::rnorm(n * ns), n, ns) %*% diag(comm$sigma, ns))), n, ns)
  }
  sites <- sprintf("S%03d", seq_len(n))
  dimnames(Y) <- list(sites, comm$species)
  site_data <- data.frame(site = sites, x = xy[, 1], y = xy[, 2],
                          temperature = 12 + 3 * zt, forest = 55 + 20 * zf)
  trait_data <- data.frame(species = comm$species, body_mass = comm$body_mass_g)
  data <- list(Y = Y, site_data = site_data, trait_data = trait_data,
               x_formula = "~ temperature + forest",
               tr_formula = "~ log(body_mass)",
               response = resp,
               coords = if (isTRUE(design$fit_spatial)) xy else NULL,
               tree = comm$tree,
               species_labels = stats::setNames(comm$species, comm$species))
  # truth
  B <- comm$Beta; G <- comm$Gamma
  Omega <- crossprod(comm$Lambda)
  extra <- if (resp == "presence") rep(1, ns) else comm$sigma^2
  as <- jsdm_assoc(Omega, extra)
  pn <- pair_names(comm$species)
  b0 <- B["intercept", ]
  tv <- c(
    stats::setNames(G[, "intercept"], paste0("comm_mean_", rownames(G))),
    stats::setNames(G[, "log_body_mass"], paste0("trait_log_body_mass_on_", rownames(G))),
    stats::setNames(as.numeric(t(B)), paste0("beta_", rep(rownames(B), each = ns), "_", comm$species)),
    if (nf > 0) stats::setNames(as$rescor, paste0("rescor_", pn)),
    stats::setNames(as$latcor, paste0("latcor_", pn)))
  if (resp == "presence") {
    pr <- stats::pnorm(b0 / sqrt(1 + diag(Omega)))
    # population-averaged (marginal) probit slopes: latent factors integrated out
    PAall <- sweep(B, 2, sqrt(1 + diag(Omega)), "/")
    PA <- PAall[-1, , drop = FALSE]
    Tr <- comm$Tr
    Gpa <- PAall %*% Tr %*% solve(crossprod(Tr))   # population-averaged trait regression (realised community)
    tv <- c(tv, stats::setNames(as.numeric(t(PA)), paste0("popavg_", rep(rownames(PA), each = ns), "_", comm$species)),
            stats::setNames(rowMeans(PA), paste0("comm_popavg_", rownames(PA))),
            stats::setNames(Gpa[, "log_body_mass"], paste0("popavg_trait_log_body_mass_on_", rownames(Gpa))),
            stats::setNames(pr, paste0("prev_", comm$species)), richness_mean_env = sum(pr))
  } else {
    ab <- exp(b0 + (diag(Omega) + comm$sigma^2) / 2)
    tv <- c(tv, stats::setNames(ab, paste0("abund_", comm$species)), total_abund_mean_env = sum(ab))
  }
  if (comm$rho > 0) tv <- c(tv, phylo_signal_rho = comm$rho)
  list(data = data, truth = tv)
}
