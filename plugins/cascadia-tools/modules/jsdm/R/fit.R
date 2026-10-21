# jsdm: fit HMSC with the Hmsc R package and turn the posterior into
# contract estimates on biological scales.

clean_name <- function(x) {
  x <- ifelse(x == "(Intercept)", "intercept", x)
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", x)
}

# Model matrix with continuous columns standardised (mean 0, SD 1); 0/1
# columns are kept as they are. Returns the scaling so predictions and
# figures can be expressed in original units.
jsdm_model_matrix <- function(formula, df, scaling = NULL, xlev = NULL, what = "sites") {
  f <- stats::as.formula(formula)
  miss <- setdiff(all.vars(f), names(df))
  if (length(miss))
    stop(sprintf("the formula '%s' uses %s, but the %s table has no such column. Its columns are: %s.",
                 paste(deparse(f), collapse = ""), paste0("'", miss, "'", collapse = ", "), what,
                 paste(names(df), collapse = ", ")), call. = FALSE)
  mf <- stats::model.frame(f, df, na.action = stats::na.pass, xlev = xlev)
  mm <- stats::model.matrix(f, mf)
  raw <- colnames(mm); nm <- clean_name(raw)
  colnames(mm) <- nm
  new <- is.null(scaling)
  if (new) scaling <- list()
  for (k in seq_along(nm)) {
    if (nm[k] == "intercept") next
    v <- mm[, k]
    if (new) {
      bin <- all(v[!is.na(v)] %in% c(0, 1))
      scaling[[nm[k]]] <- if (bin) list(center = 0, scale = 1, label = raw[k], binary = TRUE)
                          else list(center = mean(v, na.rm = TRUE), scale = stats::sd(v, na.rm = TRUE),
                                    label = raw[k], binary = FALSE)
    }
    mm[, k] <- (v - scaling[[nm[k]]]$center) / scaling[[nm[k]]]$scale
  }
  xl <- lapply(mf, function(c) if (is.factor(c) || is.character(c)) levels(factor(c)))
  list(M = mm, names = nm, raw = raw, scaling = scaling, assign = attr(mm, "assign"),
       term_labels = attr(stats::terms(f), "term.labels"), xlev = Filter(Negate(is.null), xl))
}

jsdm_prepare <- function(data, settings = list()) {
  Y <- as.matrix(data$Y)
  sp_raw <- colnames(Y)
  sp <- clean_name(sp_raw)
  if (anyDuplicated(sp)) sp <- make.unique(sp, sep = "_")
  colnames(Y) <- sp
  sd <- data$site_data
  omit <- settings$omit_covariates
  xf <- data$x_formula
  if (length(omit)) {
    # negative-control / sensitivity: drop named covariates from the formula
    tl <- attr(stats::terms(stats::as.formula(xf)), "term.labels")
    keep <- tl[!tl %in% omit]
    xf <- if (length(keep)) paste("~", paste(keep, collapse = " + ")) else "~ 1"
  }
  Xm <- jsdm_model_matrix(xf, sd, what = "site covariates")
  Trm <- NULL
  if (!is.null(data$tr_formula) && !is.null(data$trait_data) && !isTRUE(settings$no_traits)) {
    td <- data$trait_data
    Trm <- jsdm_model_matrix(data$tr_formula, td, what = "species traits")
    rownames(Trm$M) <- sp
  }
  C <- NULL
  if (!is.null(data$tree) && !isTRUE(settings$no_phylogeny)) {
    eco_require("ape")
    tr <- data$tree
    tr$tip.label <- clean_name(tr$tip.label)
    C <- ape::vcv(tr, corr = TRUE)[sp, sp]
  }
  lab <- data$species_labels
  if (is.null(lab)) lab <- stats::setNames(sp_raw, sp) else names(lab) <- clean_name(names(lab))
  list(Y = Y, sp = sp, sp_labels = lab, X = Xm, Tr = Trm, C = C, x_formula = xf,
       response = data$response %||% "presence")
}

fit_jsdm <- function(data, settings) {
  eco_require(c("Hmsc", "posterior"))
  suppressPackageStartupMessages(requireNamespace("Hmsc"))
  st <- settings
  pr <- jsdm_prepare(data, st)
  Y <- pr$Y; n <- nrow(Y); ns <- ncol(Y)
  distr <- if (pr$response == "presence") "probit" else "lognormal poisson"
  sites <- rownames(Y) %||% sprintf("S%03d", seq_len(n))
  sdes <- data.frame(site = factor(sites))
  spatial <- st$spatial %||% "none"
  if (!identical(spatial, "none") && !is.null(data$coords)) {
    xy <- as.matrix(data$coords); rownames(xy) <- sites
    rL <- Hmsc::HmscRandomLevel(sData = xy, sMethod = if (spatial == "full") "Full" else "NNGP",
                                nNeighbours = st$n_neighbours %||% 10)
  } else {
    spatial <- "none"
    rL <- Hmsc::HmscRandomLevel(units = sites)
  }
  rL <- Hmsc::setPriors(rL, nfMax = st$max_factors %||% 5, nfMin = st$min_factors %||% 2)
  fp <- st$factor_prior   # list(a1, b1, a2, b2); NULL = Hmsc defaults (a1 = a2 = 50, b1 = b2 = 1)
  if (!is.null(fp)) rL <- do.call(Hmsc::setPriors, c(list(rL), fp[intersect(names(fp), c("nu", "a1", "b1", "a2", "b2"))]))
  args <- list(Y = Y, X = pr$X$M, XScale = FALSE, distr = distr, studyDesign = sdes,
               ranLevels = list(site = rL))
  if (!is.null(pr$Tr)) { args$Tr <- pr$Tr$M; args$TrScale <- FALSE }
  if (!is.null(pr$C)) args$C <- pr$C
  m <- do.call(Hmsc::Hmsc, args)
  set.seed(st$seed %||% 1)
  nch <- st$n_chains %||% 4
  t0 <- proc.time()[["elapsed"]]
  # parallel chains run in fresh PSOCK workers: give them the ecotools library
  # path (the site Renviron would otherwise hide Rlib)
  if (nch > 1 && (st$n_parallel %||% 1) > 1) Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"))
  m <- Hmsc::sampleMcmc(m, samples = st$samples %||% 250, thin = st$thin %||% 10,
                        transient = st$transient %||% 2500, nChains = nch,
                        nParallel = min(nch, st$n_parallel %||% 1), verbose = 0,
                        updater = list(Gamma2 = FALSE, GammaEta = isTRUE(st$gamma_eta_updater)))
  t_mcmc <- proc.time()[["elapsed"]] - t0
  dr <- jsdm_draws(m, pr)
  conv_on <- st$converge_on %||% "core"
  core_cols <- grep(if (conv_on == "all") "." else "^(comm_mean_|comm_popavg_|trait_|beta_|popavg_|latcor_|prev_|abund_|richness_|total_abund_)",
                    colnames(dr$chains[[1]]), value = TRUE)
  s <- summarise_draws_eco(dr$chains)
  est <- s$estimates
  ec <- est[est$param %in% core_cols, ]
  rhat_max <- max(ec$rhat, na.rm = TRUE); ess_min <- min(ec$ess_bulk, na.rm = TRUE)
  notes <- sprintf("Hmsc %s, %d chains x %d samples (thin %d, transient %d); MCMC time %.0f s; %s random level; %d latent factors retained (max %d).",
                   as.character(utils::packageVersion("Hmsc")), nch, st$samples %||% 250, st$thin %||% 10,
                   st$transient %||% 2500, t_mcmc,
                   if (spatial == "none") "non-spatial site-level" else paste0("spatial (", spatial, ") site-level"),
                   dr$nf, st$max_factors %||% 5)
  if (rhat_max >= 1.01 || ess_min < 400)
    notes <- c(notes, sprintf("Convergence for real-data reporting needs R-hat < 1.01 and bulk ESS >= 400 for every reported quantity; here max R-hat = %.3f and min ESS = %.0f. Increase samples/thin (settings) before interpreting.",
                              rhat_max, ess_min))
  if (pr$response != "presence")
    notes <- c(notes, "Counts (lognormal Poisson) are NOT validated in ecotools: in simulation, species-level intercepts, residual SDs and residual correlations mixed very slowly (R-hat 1.3-2.0 even at 2.5x the recovery budget; VALIDATION.md section 5). Use thin >= 500 and transient >= 100000, and do not interpret species-level results unless R-hat < 1.01; or analyse presence/absence instead.")
  wr <- est[grepl("^rescor_", est$param), ]
  if (nrow(wr) && max(wr$rhat, na.rm = TRUE) >= 1.05)
    notes <- c(notes, sprintf("Hmsc-style residual correlations (rescor_*) mix slowly for species with small loadings (max R-hat %.2f); use latcor_* (latent-scale residual correlations) for inference.",
                              max(wr$rhat, na.rm = TRUE)))
  out <- list(estimates = est,
              diagnostics = list(converged = is.finite(rhat_max) && rhat_max < 1.05,
                                 rhat_max = rhat_max, ess_min = ess_min, notes = notes,
                                 t_mcmc = t_mcmc),
              prep = pr, draws = stack_draws(dr$chains), extra = dr$extra, settings = st,
              spatial = spatial, extra_nf = dr$nf,
              chains_gamma = lapply(dr$chains, function(M) M[, jsdm_comm_params(colnames(M), pr$response), drop = FALSE]))
  if (isTRUE(st$full_output %||% TRUE)) {
    out$object <- m
    out$vp <- tryCatch(jsdm_vp(m, pr), error = function(e) { message("variance partitioning failed: ", conditionMessage(e)); NULL })
    out$model_fit <- tryCatch(jsdm_model_fit(m, pr), error = function(e) { message("model fit failed: ", conditionMessage(e)); NULL })
    out$ppc <- tryCatch(jsdm_ppc(dr, pr, seed = st$seed %||% 1), error = function(e) NULL)
  }
  out
}

# Community-level parameters that are interpreted: population-averaged
# (identified) ones for probit, Gamma on the log scale for counts.
jsdm_comm_params <- function(nms, response) {
  if (response == "presence") grep("^(comm_popavg_|popavg_trait_)", nms, value = TRUE)
  else setdiff(grep("^(comm_mean_|trait_)", nms, value = TRUE), "comm_mean_intercept")
}

# Per-chain matrices of named estimands from the Hmsc posterior.
jsdm_draws <- function(m, pr) {
  sp <- pr$sp; ns <- length(sp); xn <- pr$X$names
  trn <- if (is.null(pr$Tr)) "intercept" else pr$Tr$names
  pn <- pair_names(sp)
  ut <- upper.tri(matrix(0, ns, ns))
  presence <- pr$response == "presence"
  nf_used <- 0
  Trm <- if (is.null(pr$Tr)) matrix(1, ns, 1, dimnames = list(sp, "intercept")) else pr$Tr$M
  trp <- setdiff(trn, "intercept")
  # OLS projection of species responses on traits (population-averaged trait effects)
  Pt <- Trm %*% solve(crossprod(Trm))
  one <- function(s) {
    B <- s$Beta; G <- s$Gamma
    Lam <- s$Lambda[[1]]; nf_used <<- max(nf_used, nrow(Lam))
    Om <- crossprod(Lam)
    extra <- if (presence) rep(1, ns) else s$sigma
    rc <- suppressWarnings(stats::cov2cor(Om))[ut]
    lc <- stats::cov2cor(Om + diag(extra, ns))[ut]
    b0 <- B[1, ]
    v <- c(as.numeric(G), as.numeric(t(B)), rc, lc)
    if (presence) {
      sc <- sqrt(1 + diag(Om)); p <- stats::pnorm(b0 / sc)
      PAall <- sweep(B, 2, sc, "/")
      PA <- PAall[-1, , drop = FALSE]
      Gpa <- PAall %*% Pt                     # covariates x traits
      v <- c(v, as.numeric(t(PA)), rowMeans(PA),
             if (length(trp)) as.numeric(Gpa[, trp, drop = FALSE]), p, sum(p))
    }
    else { a <- exp(b0 + (diag(Om) + s$sigma) / 2); v <- c(v, a, sum(a), sqrt(s$sigma)) }
    if (!is.null(pr$C)) v <- c(v, s$rho)
    v
  }
  gnames <- as.vector(outer(xn, trn, function(k, t) ifelse(t == "intercept", paste0("comm_mean_", k),
                                                          paste0("trait_", t, "_on_", k))))
  nms <- c(gnames, paste0("beta_", rep(xn, each = ns), "_", sp), paste0("rescor_", pn), paste0("latcor_", pn),
           if (presence) c(paste0("popavg_", rep(xn[-1], each = ns), "_", sp), paste0("comm_popavg_", xn[-1]),
                           if (length(trp)) as.vector(outer(xn, trp, function(k, t) paste0("popavg_trait_", t, "_on_", k))),
                           paste0("prev_", sp), "richness_mean_env")
           else c(paste0("abund_", sp), "total_abund_mean_env", paste0("residual_sd_", sp)),
           if (!is.null(pr$C)) "phylo_signal_rho")
  chains <- lapply(m$postList, function(pl) {
    M <- t(vapply(pl, one, numeric(length(nms))))
    colnames(M) <- nms
    M
  })
  # Raw arrays needed for predictions (Beta, Omega diag, sigma) and PPC
  extra <- lapply(m$postList, function(pl) lapply(pl, function(s)
    list(Beta = s$Beta, Lambda = s$Lambda[[1]], Eta = s$Eta[[1]], sigma = s$sigma)))
  list(chains = chains, extra = unlist(extra, recursive = FALSE), nf = nf_used)
}

# Variance partitioning (Ovaskainen et al. 2017): one group per model term.
jsdm_vp <- function(m, pr) {
  a <- pr$X$assign
  tl <- pr$X$term_labels
  if (!length(tl)) return(NULL)
  grp <- ifelse(a == 0, 1, a)
  vp <- Hmsc::computeVariancePartitioning(m, group = grp, groupnames = tl)
  vals <- vp$vals
  colnames(vals) <- pr$sp
  rownames(vals) <- c(tl, "Residual associations (site latent factors)")
  list(vals = vals, R2T = vp$R2T)
}

# In-sample explanatory power: Tjur R2 and AUC (presence) or pseudo-R2 (counts)
jsdm_model_fit <- function(m, pr) {
  pred <- Hmsc::computePredictedValues(m, expected = TRUE)
  ev <- Hmsc::evaluateModelFit(m, pred)
  df <- data.frame(species = pr$sp, prevalence = colMeans(pr$Y > 0),
                   do.call(cbind, lapply(ev, function(v) as.numeric(v))))
  # fixed-effects-only explanatory power (without the latent site scores)
  df
}

# Posterior predictive check: site species richness variance and pairwise
# co-occurrence. Replicates are drawn conditional on each posterior draw.
jsdm_ppc <- function(dr, pr, seed = 1, ndraw = 300) {
  set.seed(seed)
  X <- pr$X$M; Y <- pr$Y; n <- nrow(Y); ns <- ncol(Y)
  ex <- dr$extra
  ii <- unique(round(seq(1, length(ex), length.out = min(ndraw, length(ex)))))
  presence <- pr$response == "presence"
  stat <- function(Yb) c(var_richness = stats::var(rowSums(Yb > 0)),
                         mean_cooccur = mean(crossprod(Yb > 0)[upper.tri(diag(ns))]))
  rep <- t(vapply(ii, function(k) {
    s <- ex[[k]]
    L <- X %*% s$Beta + s$Eta %*% s$Lambda
    Yr <- if (presence) (L + matrix(stats::rnorm(n * ns), n, ns) > 0) * 1
          else matrix(stats::rpois(n * ns, exp(L + matrix(stats::rnorm(n * ns), n, ns) %*% diag(sqrt(s$sigma), ns))), n, ns)
    stat(Yr)
  }, numeric(2)))
  # Marginal replicate (new sites): latent site scores redrawn, tests
  # whether the factors reproduce the richness variance in new sites
  rep_new <- t(vapply(ii, function(k) {
    s <- ex[[k]]
    Eta <- matrix(stats::rnorm(n * nrow(s$Lambda)), n, nrow(s$Lambda))
    L <- X %*% s$Beta + Eta %*% s$Lambda
    Yr <- if (presence) (L + matrix(stats::rnorm(n * ns), n, ns) > 0) * 1
          else matrix(stats::rpois(n * ns, exp(L + matrix(stats::rnorm(n * ns), n, ns) %*% diag(sqrt(s$sigma), ns))), n, ns)
    stat(Yr)
  }, numeric(2)))
  list(observed = stat(Y), replicated = rep, replicated_new_sites = rep_new,
       p_value = colMeans(sweep(rep_new, 2, stat(Y), ">=")))
}

# Response-scale predictions along one covariate (original units), other
# covariates at their mean (factors at their most common level).
# Returns per-species and community (richness / total abundance) draws.
jsdm_gradient <- function(fit, data, cov, n_grid = 60, ndraw = 400, xs = NULL) {
  pr <- fit$prep
  sd <- data$site_data
  vars <- all.vars(stats::as.formula(pr$x_formula))
  base <- lapply(vars, function(v) {
    x <- sd[[v]]
    if (is.numeric(x)) mean(x, na.rm = TRUE) else names(sort(table(x), decreasing = TRUE))[1]
  })
  names(base) <- vars
  if (is.null(xs)) xs <- seq(min(sd[[cov]], na.rm = TRUE), max(sd[[cov]], na.rm = TRUE), length.out = n_grid)
  n_grid <- length(xs)
  nd <- as.data.frame(lapply(base, function(v) rep(v, n_grid)), stringsAsFactors = FALSE)
  nd[[cov]] <- xs
  for (v in vars) if (!is.numeric(sd[[v]])) nd[[v]] <- factor(nd[[v]], levels = levels(factor(sd[[v]])))
  Xn <- jsdm_model_matrix(pr$x_formula, nd, scaling = pr$X$scaling, xlev = pr$X$xlev)$M
  ex <- fit$extra
  ii <- unique(round(seq(1, length(ex), length.out = min(ndraw, length(ex)))))
  presence <- pr$response == "presence"
  ns <- length(pr$sp)
  arr <- vapply(ii, function(k) {
    s <- ex[[k]]
    om <- colSums(s$Lambda^2)
    L <- Xn %*% s$Beta
    if (presence) stats::pnorm(sweep(L, 2, sqrt(1 + om), "/"))
    else exp(sweep(L, 2, (om + s$sigma) / 2, "+"))
  }, matrix(0, n_grid, ns))
  list(x = xs, sp = arr, community = apply(arr, c(1, 3), sum))
}

summ_band <- function(M, ci = 0.95, inner = 0.5) {
  a <- (1 - ci) / 2; b <- (1 - inner) / 2
  data.frame(estimate = apply(M, 1, mean),
             lower = apply(M, 1, stats::quantile, a), upper = apply(M, 1, stats::quantile, 1 - a),
             lower50 = apply(M, 1, stats::quantile, b), upper50 = apply(M, 1, stats::quantile, 1 - b))
}
