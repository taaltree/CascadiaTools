# Generic GLM/GLMM simulator, driven by `design` (sample size, groups, family,
# covariates, offset, trials) and `params` (coefficients per 1 SD of each
# continuous covariate, random-intercept SD, dispersion, zero-inflation).
#
# Continuous covariates are standardised with the SAMPLE mean and SD before
# the linear predictor is formed, exactly as fit() does, so beta_<x> is the
# true per-SD coefficient of the fitted model and beta_intercept is the true
# value at the covariate means and reference factor levels.

glm_simulate <- function(params, design, seed) {
  set.seed(seed)
  fam <- canon_family(design$family)
  link <- design$link %||% default_link(fam)
  ng <- as.integer(design$n_groups %||% 0)
  gname <- design$group_name %||% "site"
  if (ng > 0) {
    npg <- as.integer(design$n_per_group %||% 5)
    n <- ng * npg; grp <- rep(seq_len(ng), each = npg)
  } else { n <- as.integer(design$n %||% 200); grp <- NULL }
  covs <- design$covariates %||% list()
  df <- data.frame(row.names = seq_len(n))
  for (v in names(covs)) {
    cs <- covs[[v]]
    at_group <- identical(cs$level, "group") && ng > 0
    m <- if (at_group) ng else n
    val <- switch(cs$type %||% "normal",
      factor = factor(sample(rep(unlist(cs$levels), length.out = m)), levels = unlist(cs$levels)),
      uniform = stats::runif(m, cs$min, cs$max),
      normal = stats::rnorm(m, cs$mean %||% 0, cs$sd %||% 1),
      stop("covariate type must be normal, uniform or factor"))
    if (at_group) val <- val[grp]
    df[[v]] <- val
  }
  if (ng > 0) df[[gname]] <- factor(sprintf(paste0("%s%0", nchar(ng), "d"), toupper(substr(gname, 1, 1)), grp))
  terms <- c(names(covs), unlist(design$interactions))
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
  Z <- df
  for (v in names(covs)) if (is.numeric(Z[[v]])) Z[[v]] <- (Z[[v]] - mean(Z[[v]])) / stats::sd(Z[[v]])
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), Z)
  xlev <- lapply(df[vapply(df, is.factor, logical(1))], levels)
  labs <- vapply(colnames(X), coef_label, "", xlev = xlev)
  pnames <- paste0("beta_", labs)
  unknown <- setdiff(grep("^beta_", names(params), value = TRUE), pnames)
  if (length(unknown)) stop("params not in the simulated design: ", paste(unknown, collapse = ", "),
                            ". Coefficients available: ", paste(pnames, collapse = ", "))
  b <- vapply(pnames, function(p) as.numeric(params[[p]] %||% 0), 0)
  eta <- drop(X %*% b)
  re_sd <- as.numeric(params$re_sd %||% 0)
  if (ng > 0) eta <- eta + stats::rnorm(ng, 0, re_sd)[grp]
  oname <- NULL
  if (!is.null(design$offset)) {
    oname <- design$offset$name %||% "effort"
    eff <- round(stats::runif(n, design$offset$min %||% 0.5, design$offset$max %||% 2), 2)
    df[[oname]] <- eff
    eta <- eta + log(eff)
  }
  mu <- link_inv(link)(eta)
  disp <- as.numeric(params$dispersion %||% default_dispersion(fam))
  tr <- design$trials %||% 1
  trials <- if (length(tr) == 2) sample(seq(tr[1], tr[2]), n, replace = TRUE) else rep(as.integer(tr), n)
  y <- switch(fam,
    gaussian = stats::rnorm(n, mu, disp),
    binomial = stats::rbinom(n, trials, mu),
    poisson = , zipoisson = stats::rpois(n, mu),
    nbinom2 = , zinb = stats::rnbinom(n, size = disp, mu = mu),
    beta = pmin(pmax(stats::rbeta(n, mu * disp, (1 - mu) * disp), 1e-6), 1 - 1e-6),
    Gamma = stats::rgamma(n, shape = 1 / disp^2, rate = 1 / (disp^2 * mu)))
  zp <- 0
  if (is_zi(fam)) { zp <- as.numeric(params$zi_prob %||% 0.2); y[stats::runif(n) < zp] <- 0 }
  resp <- design$response %||% "y"
  df[[resp]] <- y
  tcol <- NULL
  if (fam == "binomial" && any(trials > 1)) { tcol <- "trials"; df$trials <- trials }
  spec <- list(formula = paste(resp, "~", rhs, if (ng > 0) sprintf("+ (1|%s)", gname)),
               family = fam, link = link, offset = oname, offset_per = 1, trials = tcol,
               scale = "all", units = design$units, labels = design$labels, response_label = resp,
               zi = if (is_zi(fam)) "~1")
  truth <- stats::setNames(b, pnames)
  rp <- ratio_prefix(fam, link)
  if (!is.na(rp) && length(b) > 1) truth <- c(truth, stats::setNames(exp(b[-1]), paste0(rp, "_", labs[-1])))
  truth["mean_response_at_reference"] <- link_inv(link)(b[["beta_intercept"]]) * (1 - zp)
  if (!is.na(disp_name(fam))) truth[disp_name(fam)] <- disp
  if (is_zi(fam)) truth["zero_inflation_prob"] <- zp
  if (ng > 0) truth[paste0("re_sd_", clean_name(gname))] <- re_sd
  list(data = list(df = df, spec = spec), truth = truth)
}
